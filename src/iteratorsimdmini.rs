//! SIMD-accelerated syncmer-based superkmer iterator.
//!
//! Uses `simd-minimizers` to find canonical closed syncmer positions via SIMD,
//! then runs a rightmost-wins MSP sliding window to produce superkmers.
//! Requires odd l (e.g. l=9). Takes ASCII `&[u8]` input directly.
//! Mint is canonical by default; use `SuperkmersIterator::non_canonical()` for
//! forward-strand mint.
//!
//! ```ignore
//! use rust_superkmers::iteratorsimdmini::SuperkmersIterator;
//! let iter = SuperkmersIterator::new(seq, 31, 9);           // canonical
//! let iter = SuperkmersIterator::non_canonical(seq, 31, 9); // forward-strand
//! ```
use crate::{Superkmer, SplitMode};

const SMER_SIZE: usize = 2; // syncmer's s parameter

/// Lookup table: ASCII byte -> 2-bit encoding (A=0, C=1, G=2, T=3)
const ASCII_TO_2BIT: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

/// Encode forward l-mer from ASCII bytes and look up canonical form in table.
#[inline]
fn canonical_lmer_index(ascii: &[u8], pos: usize, l: usize) -> (usize, bool) {
    let mut fwd = 0usize;
    for i in 0..l {
        let base = ASCII_TO_2BIT[ascii[pos + i] as usize] as usize;
        fwd = (fwd << 2) | base;
    }
    let canonical_table = match l {
        8 => &*crate::CANONICAL_8,
        9 => &*crate::CANONICAL_9,
        10 => &*crate::CANONICAL_10,
        12 => &*crate::CANONICAL_12,
        _ => panic!("Unsupported l={} for canonical lookup", l),
    };
    let (canonical, is_rc) = canonical_table[fwd];
    (canonical as usize, is_rc)
}

/// Encode forward l-mer from ASCII bytes (no canonicalization).
#[inline]
fn forward_lmer_index(ascii: &[u8], pos: usize, l: usize) -> usize {
    let mut fwd = 0usize;
    for i in 0..l {
        let base = ASCII_TO_2BIT[ascii[pos + i] as usize] as usize;
        fwd = (fwd << 2) | base;
    }
    fwd
}

/// Collect superkmers from a single N-free fragment using SIMD syncmer detection.
fn superkmers_from_fragment(
    ascii_slice: &[u8],
    k: usize,
    l: usize,
    offset: usize,
    canonical: bool,
    mode: SplitMode,
    results: &mut Vec<Superkmer>,
    syncmer_pos: &mut Vec<u32>,
) {
    use simd_minimizers::packed_seq::AsciiSeq;

    let seq_len = ascii_slice.len();
    if seq_len < k {
        return;
    }

    assert!(
        l % 2 == 1,
        "simd-mini canonical syncmers require odd l (try l=9), got l={}",
        l
    );

    // simd-minimizers handles both uppercase and lowercase ASCII (ACTGactg)
    let w_sync = l - SMER_SIZE + 1;
    syncmer_pos.clear();
    simd_minimizers::canonical_closed_syncmers(SMER_SIZE, w_sync)
        .run(AsciiSeq(ascii_slice), syncmer_pos);

    // Check if the l-mer at a syncmer position is all-A or all-T (demoted).
    let is_demoted = |pos: u32| -> bool {
        let p = pos as usize;
        let b0 = ascii_slice[p] & 0xDF; // case-insensitive: clear bit 5
        (b0 == b'A' || b0 == b'T') && ascii_slice[p..p + l].iter().all(|&b| b & 0xDF == b0)
    };

    let num_kmers = seq_len - k + 1;
    let max_lmer_offset = k - l;
    results.reserve(syncmer_pos.len());

    // Helper: emit a superkmer from sk_start to the k-mer at last_kmer_idx.
    let emit = |sk_start: usize, last_kmer_idx: usize, min_pos: u32,
                results: &mut Vec<Superkmer>| {
        let start = sk_start;
        let end = last_kmer_idx + k;
        let size = end - start;
        let mpos_fwd = min_pos as usize - start;
        let (mint, mint_is_rc) = if canonical {
            canonical_lmer_index(ascii_slice, min_pos as usize, l)
        } else {
            (forward_lmer_index(ascii_slice, min_pos as usize, l), false)
        };
        results.push(Superkmer {
            start: start + offset,
            mint: mint as u32,
            size: size as u16,
            mpos: mpos_fwd as u16,
            mint_is_rc,
        });
    };

    if mode != SplitMode::Sticky {
        // Classical: rightmost non-demoted syncmer in window, split on any change.
        // Msp/MspXor: lowest-score non-demoted syncmer in window, split on any change.
        // All are context-independent.
        if num_kmers == 0 { return; }

        const XOR_CONSTANT: usize = 0xACE5_ACE5;

        // Precompute scores for all syncmer positions (avoids per-k-mer canonical_lmer_index calls)
        let syncmer_scores: Vec<usize> = if mode == SplitMode::Classical {
            Vec::new() // Classical doesn't use scores
        } else {
            syncmer_pos.iter().map(|&p| {
                let (canon, _) = canonical_lmer_index(ascii_slice, p as usize, l);
                let tiebreak = if mode == SplitMode::MspXor { canon ^ XOR_CONSTANT } else { canon };
                if is_demoted(p) { (1usize << 32) | tiebreak } else { tiebreak }
            }).collect()
        };

        let mut lo = 0usize; // first syncmer index >= window left
        let mut hi = 0usize; // first syncmer index > window right

        // Rescan the window [lo..hi) to find the best syncmer.
        let find_best = |lo: usize, hi: usize| -> (u32, usize) {
            if lo >= hi { return (u32::MAX, usize::MAX); }
            if mode == SplitMode::Classical {
                // Rightmost non-demoted, fallback to rightmost demoted
                if !is_demoted(syncmer_pos[hi - 1]) {
                    return (syncmer_pos[hi - 1], 0);
                }
                let mut j = hi - 1;
                while j > lo {
                    j -= 1;
                    if !is_demoted(syncmer_pos[j]) { return (syncmer_pos[j], 0); }
                }
                (syncmer_pos[hi - 1], 0)
            } else {
                // Msp/MspXor: rightmost with lowest precomputed score
                let mut best_pos = u32::MAX;
                let mut best_score = usize::MAX;
                for j in lo..hi {
                    if syncmer_scores[j] <= best_score {
                        best_score = syncmer_scores[j];
                        best_pos = syncmer_pos[j];
                    }
                }
                (best_pos, best_score)
            }
        };

        // Initialize for first k-mer
        while hi < syncmer_pos.len() && (syncmer_pos[hi] as usize) <= max_lmer_offset {
            hi += 1;
        }
        let (mut curr_min, mut curr_score) = find_best(lo, hi);
        let mut sk_start = 0usize;
        // Track which syncmer index curr_min corresponds to, for Classical falloff detection
        let mut prev_hi = hi;

        for i in 1..num_kmers {
            let old_lo = lo;
            // Advance lo: remove syncmers that fell off left edge
            while lo < syncmer_pos.len() && (syncmer_pos[lo] as usize) < i {
                lo += 1;
            }
            let old_hi = hi;
            // Advance hi: add syncmers entering from right edge
            while hi < syncmer_pos.len() && (syncmer_pos[hi] as usize) <= i + max_lmer_offset {
                hi += 1;
            }

            // Did the current minimum fall off the left edge?
            let fell_off = curr_min != u32::MAX && (curr_min as usize) < i;

            if fell_off {
                // Rescan to find new best
                let (new_min, new_score) = find_best(lo, hi);
                if new_min != curr_min {
                    if curr_min != u32::MAX {
                        emit(sk_start, i - 1, curr_min, results);
                    }
                    sk_start = i;
                }
                curr_min = new_min;
                curr_score = new_score;
            } else if mode == SplitMode::Classical {
                // Classical: rightmost syncmer wins. Check if a new syncmer entered.
                if hi > old_hi {
                    // New syncmer(s) entered from right — rightmost is at hi-1
                    let new_pos = syncmer_pos[hi - 1];
                    let new_demoted = is_demoted(new_pos);
                    let curr_demoted = curr_min != u32::MAX && is_demoted(curr_min);
                    // New non-demoted beats current, or new rightmost replaces current if same demotion status
                    if (!new_demoted && curr_demoted) || (!new_demoted || curr_demoted) {
                        if new_pos != curr_min {
                            if curr_min != u32::MAX {
                                emit(sk_start, i - 1, curr_min, results);
                            }
                            sk_start = i;
                            curr_min = new_pos;
                        }
                    }
                }
            } else {
                // Msp/MspXor: check new syncmers entering from right
                for j in old_hi..hi {
                    let new_score_j = syncmer_scores[j];
                    if new_score_j <= curr_score {
                        let new_pos = syncmer_pos[j];
                        if new_pos != curr_min {
                            if curr_min != u32::MAX {
                                emit(sk_start, i - 1, curr_min, results);
                            }
                            sk_start = i;
                            curr_min = new_pos;
                            curr_score = new_score_j;
                        }
                    }
                }
            }
        }
        // Flush final superkmer
        if curr_min != u32::MAX {
            emit(sk_start, num_kmers - 1, curr_min, results);
        }
    } else {
        // Sticky MSP: track RIGHTMOST syncmer as minimizer, only re-evaluate
        // when current minimizer falls off the left edge.
        let mut ptr = 0usize;
        let mut curr_min_pos: u32 = u32::MAX;
        let mut sk_start: usize = 0;

        // After taking rightmost in syncmer_pos[from..j_end], demote if needed.
        macro_rules! demote_check {
            ($from:expr, $j_end:expr) => {
                if curr_min_pos != u32::MAX && is_demoted(curr_min_pos) {
                    let mut best = u32::MAX;
                    let mut di = $from;
                    while di < $j_end {
                        if !is_demoted(syncmer_pos[di]) { best = syncmer_pos[di]; }
                        di += 1;
                    }
                    if best != u32::MAX { curr_min_pos = best; }
                }
            };
        }

        // Initialize: find rightmost syncmer in first k-mer window [0, k-l]
        if num_kmers > 0 {
            let start_ptr = ptr;
            while ptr < syncmer_pos.len() && (syncmer_pos[ptr] as usize) <= max_lmer_offset {
                curr_min_pos = syncmer_pos[ptr];
                ptr += 1;
            }
            if ptr > start_ptr { demote_check!(start_ptr, ptr); }
        }

        let mut i = 1usize;
        while i <= num_kmers {
            let need_new_min = if i < num_kmers {
                (curr_min_pos as usize) < i
            } else {
                true // sentinel: flush last superkmer
            };

            if need_new_min {
                if curr_min_pos != u32::MAX {
                    emit(sk_start, i - 1, curr_min_pos, results);
                }

                sk_start = i;
                curr_min_pos = u32::MAX;
                if i < num_kmers {
                    while ptr < syncmer_pos.len() && (syncmer_pos[ptr] as usize) < i {
                        ptr += 1;
                    }
                    let scan_from = ptr;
                    let mut j = ptr;
                    while j < syncmer_pos.len() && (syncmer_pos[j] as usize) <= i + max_lmer_offset {
                        curr_min_pos = syncmer_pos[j];
                        j += 1;
                    }
                    if j > scan_from { demote_check!(scan_from, j); }
                    // No syncmer in this window — jump to next syncmer position
                    if curr_min_pos == u32::MAX && ptr < syncmer_pos.len() {
                        let next_sync = syncmer_pos[ptr] as usize;
                        let jump_to = if next_sync > max_lmer_offset {
                            next_sync - max_lmer_offset
                        } else {
                            1
                        };
                        if jump_to > i {
                            sk_start = jump_to;
                            i = jump_to;
                            let scan_from2 = ptr;
                            let mut j2 = ptr;
                            while j2 < syncmer_pos.len() && (syncmer_pos[j2] as usize) <= i + max_lmer_offset {
                                curr_min_pos = syncmer_pos[j2];
                                j2 += 1;
                            }
                            if j2 > scan_from2 { demote_check!(scan_from2, j2); }
                        }
                    }
                }
            }
            // Jump directly to when the current minimizer falls off the left edge.
            if curr_min_pos != u32::MAX {
                let jump = curr_min_pos as usize + 1;
                i = if jump > i + 1 { std::cmp::min(jump, num_kmers) } else { i + 1 };
            } else {
                i += 1;
            }
        }
    }
}

pub struct SuperkmersIterator {
    superkmers: Vec<Superkmer>,
    storage: Vec<u64>,
    pos: usize,
}

/// Generate iterator constructors for each (mode, canonical) combination.
macro_rules! iter_constructors {
    ($($name:ident, $name_n:ident, $canonical:expr, $mode:expr;)*) => {
        $(
            pub fn $name(seq: &[u8], k: usize, l: usize) -> Self {
                Self::new_inner_full(seq, k, l, $canonical, $mode)
            }
            pub fn $name_n(seq: &[u8], k: usize, l: usize) -> Self {
                Self::new_with_n_inner_full(seq, k, l, $canonical, $mode)
            }
        )*
    };
}

impl SuperkmersIterator {
    iter_constructors! {
        new,                            new_with_n,                            true,  SplitMode::Sticky;
        non_canonical,                  non_canonical_with_n,                  false, SplitMode::Sticky;
        classical,                      classical_with_n,                      true,  SplitMode::Classical;
        classical_non_canonical,        classical_non_canonical_with_n,        false, SplitMode::Classical;
        msp,                            msp_with_n,                            true,  SplitMode::Msp;
        msp_non_canonical,              msp_non_canonical_with_n,              false, SplitMode::Msp;
        mspxor,                         mspxor_with_n,                         true,  SplitMode::MspXor;
        mspxor_non_canonical,           mspxor_non_canonical_with_n,           false, SplitMode::MspXor;
    }

    /// Access the 2-bit packed representation of the input sequence.
    pub fn storage(&self) -> &[u64] {
        &self.storage
    }

    fn new_inner_full(seq: &[u8], k: usize, l: usize, canonical: bool, mode: SplitMode) -> Self {
        let storage = crate::utils::bitpack_fragment(seq);
        let mut superkmers = Vec::new();
        let mut syncmer_pos = Vec::new();
        superkmers_from_fragment(seq, k, l, 0, canonical, mode, &mut superkmers, &mut syncmer_pos);
        SuperkmersIterator { superkmers, storage, pos: 0 }
    }

    fn new_with_n_inner_full(seq: &[u8], k: usize, l: usize, canonical: bool, mode: SplitMode) -> Self {
        let storage = crate::utils::bitpack_fragment(seq);
        let fragments = crate::utils::split_on_n(seq, k);
        let mut superkmers = Vec::new();
        let mut syncmer_pos = Vec::new();
        for (offset, fragment) in &fragments {
            superkmers_from_fragment(fragment, k, l, *offset, canonical, mode, &mut superkmers, &mut syncmer_pos);
        }
        SuperkmersIterator { superkmers, storage, pos: 0 }
    }
}

impl Iterator for SuperkmersIterator {
    type Item = Superkmer;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.superkmers.len() {
            return None;
        }
        let idx = self.pos;
        self.pos += 1;
        // Move out of vec to avoid clone
        Some(Superkmer {
            start: self.superkmers[idx].start,
            mint: self.superkmers[idx].mint,
            size: self.superkmers[idx].size,
            mpos: self.superkmers[idx].mpos,
            mint_is_rc: self.superkmers[idx].mint_is_rc,
        })
    }
}

/// Reusable superkmer extractor that avoids per-read allocations.
/// Create once, call `process()` or `process_with_n()` for each read.
pub struct SuperkmerExtractor {
    superkmers: Vec<Superkmer>,
    syncmer_pos: Vec<u32>,
    storage: Vec<u64>,
    k: usize,
    l: usize,
    canonical: bool,
    mode: SplitMode,
}

/// Generate extractor constructors for each (mode, canonical) combination.
macro_rules! extractor_constructors {
    ($($name:ident, $canonical:expr, $mode:expr;)*) => {
        $(
            pub fn $name(k: usize, l: usize) -> Self {
                Self::new_inner_full(k, l, $canonical, $mode)
            }
        )*
    };
}

impl SuperkmerExtractor {
    extractor_constructors! {
        new,                       true,  SplitMode::Sticky;
        non_canonical,             false, SplitMode::Sticky;
        classical,                 true,  SplitMode::Classical;
        classical_non_canonical,   false, SplitMode::Classical;
        msp,                       true,  SplitMode::Msp;
        msp_non_canonical,         false, SplitMode::Msp;
        mspxor,                    true,  SplitMode::MspXor;
        mspxor_non_canonical,      false, SplitMode::MspXor;
    }

    fn new_inner_full(k: usize, l: usize, canonical: bool, mode: SplitMode) -> Self {
        SuperkmerExtractor {
            superkmers: Vec::new(),
            syncmer_pos: Vec::new(),
            storage: Vec::new(),
            k,
            l,
            canonical,
            mode,
        }
    }

    /// Process a sequence with no N characters. Returns the superkmers slice.
    pub fn process(&mut self, seq: &[u8]) -> &[Superkmer] {
        self.superkmers.clear();
        let num_words = (seq.len() + 31) / 32;
        self.storage.resize(num_words, 0);
        crate::utils::bitpack_fragment_into(seq, &mut self.storage);
        superkmers_from_fragment(seq, self.k, self.l, 0, self.canonical, self.mode, &mut self.superkmers, &mut self.syncmer_pos);
        &self.superkmers
    }

    /// Process a sequence that may contain N/n characters. Returns the superkmers slice.
    pub fn process_with_n(&mut self, seq: &[u8]) -> &[Superkmer] {
        self.superkmers.clear();
        let num_words = (seq.len() + 31) / 32;
        self.storage.resize(num_words, 0);
        crate::utils::bitpack_fragment_into(seq, &mut self.storage);
        let fragments = crate::utils::split_on_n(seq, self.k);
        for (offset, fragment) in &fragments {
            superkmers_from_fragment(fragment, self.k, self.l, *offset, self.canonical, self.mode, &mut self.superkmers, &mut self.syncmer_pos);
        }
        &self.superkmers
    }

    /// Access the 2-bit packed representation of the last processed sequence.
    pub fn storage(&self) -> &[u64] {
        &self.storage
    }
}

