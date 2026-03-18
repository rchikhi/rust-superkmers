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

/// Encode ASCII base to 2-bit (A=0, C=1, G=2, T=3) via bit tricks.
/// Works for both uppercase and lowercase (case-insensitive).
#[inline(always)]
fn encode_base(b: u8) -> usize {
    (((b >> 1) ^ (b >> 2)) & 3) as usize
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
    deque: &mut Vec<usize>,
    score_buf: &mut Vec<usize>,
    mint_buf: &mut Vec<u32>,
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

    let num_kmers = seq_len - k + 1;
    let max_lmer_offset = k - l;
    results.reserve(syncmer_pos.len());

    if mode != SplitMode::Sticky {
        // Non-sticky modes: monotonic deque over syncmer positions.
        // No lookup tables — scores and mints computed inline from ASCII.
        if num_kmers == 0 { return; }

        let xor_mask = if mode == SplitMode::MspXor { 0xACE5_ACE5usize } else { 0 };
        let n_sync = syncmer_pos.len();

        // Precompute score and mint for each syncmer position.
        // Compute forward l-mer from ASCII (bytes hot in cache from SIMD scan),
        // then canonical via bit ops. No bitpacking, no tables.
        score_buf.clear();
        score_buf.reserve(n_sync);
        mint_buf.clear();
        mint_buf.reserve(n_sync);

        for &p in syncmer_pos.iter() {
            // Encode forward l-mer directly from ASCII (cache-hot from SIMD scan)
            let pos = p as usize;
            let mut fwd = 0usize;
            for j in 0..l {
                fwd = (fwd << 2) | encode_base(ascii_slice[pos + j]);
            }

            // Compute canonical = min(fwd, rc(fwd)) via bit manipulation
            let mut rc = 0usize;
            let mut v = fwd;
            for _ in 0..l {
                rc = (rc << 2) | (3 - (v & 3));
                v >>= 2;
            }
            let (canon_val, is_rc) = if rc < fwd { (rc as u32, true) } else { (fwd as u32, false) };

            // Score: for Classical all syncmers equal (score=0), for Msp/MspXor
            // use composite (demoted_priority, tiebreak).
            let score = if mode == SplitMode::Classical {
                0usize
            } else {
                let tiebreak = canon_val as usize ^ xor_mask;
                // Demoted = all-A or all-T. canon(AAA...A) = 0, canon(TTT...T) = 0.
                if canon_val == 0 { (1usize << 32) | tiebreak } else { tiebreak }
            };
            score_buf.push(score);

            // Cache mint with is_rc packed in bit 31
            if canonical {
                mint_buf.push(canon_val | if is_rc { 1u32 << 31 } else { 0 });
            } else {
                mint_buf.push(fwd as u32);
            }
        }
        let syncmer_score = &score_buf[..];
        let syncmer_mint = &mint_buf[..];

        // Helper: emit a superkmer. All lookups are from cached per-syncmer arrays.
        let emit_cached = |sk_start: usize, last_kmer_idx: usize, sync_idx: usize,
                    results: &mut Vec<Superkmer>| {
            let min_pos = syncmer_pos[sync_idx] as usize;
            let end = last_kmer_idx + k;
            let packed_mint = syncmer_mint[sync_idx];
            let mint = packed_mint & 0x7FFF_FFFF;
            let mint_is_rc = canonical && (packed_mint & (1u32 << 31)) != 0;
            results.push(Superkmer {
                start: sk_start + offset,
                mint,
                size: (end - sk_start) as u16,
                mpos: (min_pos - sk_start) as u16,
                mint_is_rc,
            });
        };

        // Monotonic deque (Vec + head index): stores syncmer indices. Front = current
        // best (lowest score, rightmost on tie). Scores are non-decreasing from front
        // to back. For "rightmost wins on tie" (<=), pop from back while back >= new.
        deque.clear();
        let mut dq_head = 0usize;

        // Inline helpers for the Vec-based deque (head..len).
        // front = deque[dq_head], back = deque[len-1], empty = dq_head >= len
        macro_rules! dq_front {
            () => { if dq_head < deque.len() { Some(deque[dq_head]) } else { None } };
        }
        macro_rules! dq_push {
            ($idx:expr) => {
                let idx = $idx;
                while dq_head < deque.len() && syncmer_score[*deque.last().unwrap()] >= syncmer_score[idx] {
                    deque.pop();
                }
                deque.push(idx);
            };
        }
        macro_rules! dq_pop_front_expired {
            ($min_pos:expr) => {
                while dq_head < deque.len() && (syncmer_pos[deque[dq_head]] as usize) < $min_pos {
                    dq_head += 1;
                }
            };
        }
        macro_rules! dq_best {
            () => { if let Some(f) = dq_front!() { f } else { usize::MAX } };
        }

        let mut hi = 0usize; // first syncmer index with pos > window right

        // Initialize for first k-mer (window [0, max_lmer_offset])
        while hi < n_sync && (syncmer_pos[hi] as usize) <= max_lmer_offset {
            dq_push!(hi);
            hi += 1;
        }

        let mut curr_best_idx = dq_best!();
        let mut sk_start = 0usize;

        // Event-driven loop with monotonic deque
        let mut i = 1usize;
        while i < num_kmers {
            if curr_best_idx == usize::MAX {
                // No minimizer — jump to when next syncmer enters
                if hi >= n_sync { break; }
                let next_i = (syncmer_pos[hi] as usize).saturating_sub(max_lmer_offset).max(i);

                dq_pop_front_expired!(next_i);

                while hi < n_sync && (syncmer_pos[hi] as usize) <= next_i + max_lmer_offset {
                    dq_push!(hi);
                    hi += 1;
                }

                curr_best_idx = dq_best!();
                sk_start = next_i;
                i = next_i + 1;
                continue;
            }

            let curr_min_pos = syncmer_pos[curr_best_idx] as usize;
            let falloff_at = curr_min_pos + 1; // k-mer index where current min leaves window

            let next_entry_at = if hi < n_sync {
                (syncmer_pos[hi] as usize).saturating_sub(max_lmer_offset).max(1)
            } else {
                num_kmers
            };

            if next_entry_at < falloff_at && next_entry_at < num_kmers {
                // A syncmer enters before the current minimum falls off.
                let entry_i = next_entry_at;

                dq_pop_front_expired!(entry_i);

                while hi < n_sync && (syncmer_pos[hi] as usize) <= entry_i + max_lmer_offset {
                    dq_push!(hi);
                    hi += 1;
                }

                let new_best_idx = dq_best!();
                if new_best_idx != curr_best_idx {
                    emit_cached(sk_start, entry_i - 1, curr_best_idx, results);
                    sk_start = entry_i;
                    curr_best_idx = new_best_idx;
                }
                i = entry_i + 1;
            } else {
                // Current minimum falls off first (or no more entries).
                let emit_end = (falloff_at - 1).min(num_kmers - 1);
                emit_cached(sk_start, emit_end, curr_best_idx, results);

                i = falloff_at;
                if i >= num_kmers { break; }

                dq_pop_front_expired!(i);

                while hi < n_sync && (syncmer_pos[hi] as usize) <= i + max_lmer_offset {
                    dq_push!(hi);
                    hi += 1;
                }

                curr_best_idx = dq_best!();
                sk_start = i;
                i += 1;
            }
        }
        // Flush final superkmer
        if curr_best_idx != usize::MAX && sk_start < num_kmers {
            emit_cached(sk_start, num_kmers - 1, curr_best_idx, results);
        }
    } else {
        // Sticky MSP: track RIGHTMOST syncmer as minimizer, only re-evaluate
        // when current minimizer falls off the left edge.

        // Check if the l-mer at a syncmer position is all-A or all-T (demoted).
        let is_demoted = |pos: u32| -> bool {
            let p = pos as usize;
            let b0 = ascii_slice[p] & 0xDF; // case-insensitive: clear bit 5
            (b0 == b'A' || b0 == b'T') && ascii_slice[p..p + l].iter().all(|&b| b & 0xDF == b0)
        };

        // Helper: emit a superkmer from sk_start to the k-mer at last_kmer_idx.
        // Computes canonical l-mer inline via bit manipulation (no table needed).
        let emit = |sk_start: usize, last_kmer_idx: usize, min_pos: u32,
                    results: &mut Vec<Superkmer>| {
            let start = sk_start;
            let end = last_kmer_idx + k;
            let size = end - start;
            let mpos_fwd = min_pos as usize - start;
            let p = min_pos as usize;
            let mut fwd = 0usize;
            for j in 0..l {
                fwd = (fwd << 2) | encode_base(ascii_slice[p + j]);
            }
            let (mint, mint_is_rc) = if canonical {
                let mut rc = 0usize;
                let mut v = fwd;
                for _ in 0..l {
                    rc = (rc << 2) | (3 - (v & 3));
                    v >>= 2;
                }
                if rc < fwd { (rc, true) } else { (fwd, false) }
            } else {
                (fwd, false)
            };
            results.push(Superkmer {
                start: start + offset,
                mint: mint as u32,
                size: size as u16,
                mpos: mpos_fwd as u16,
                mint_is_rc,
            });
        };

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
        let mut deque = Vec::new();
        let mut score_buf = Vec::new();
        let mut mint_buf = Vec::new();
        superkmers_from_fragment(seq, k, l, 0, canonical, mode, &mut superkmers, &mut syncmer_pos, &mut deque, &mut score_buf, &mut mint_buf);
        SuperkmersIterator { superkmers, storage, pos: 0 }
    }

    fn new_with_n_inner_full(seq: &[u8], k: usize, l: usize, canonical: bool, mode: SplitMode) -> Self {
        let storage = crate::utils::bitpack_fragment(seq);
        let fragments = crate::utils::split_on_n(seq, k);
        let mut superkmers = Vec::new();
        let mut syncmer_pos = Vec::new();
        let mut deque = Vec::new();
        let mut score_buf = Vec::new();
        let mut mint_buf = Vec::new();
        for (offset, fragment) in &fragments {
            superkmers_from_fragment(fragment, k, l, *offset, canonical, mode, &mut superkmers, &mut syncmer_pos, &mut deque, &mut score_buf, &mut mint_buf);
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
    deque: Vec<usize>,
    score_buf: Vec<usize>,
    mint_buf: Vec<u32>,
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
            deque: Vec::new(),
            score_buf: Vec::new(),
            mint_buf: Vec::new(),
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
        superkmers_from_fragment(seq, self.k, self.l, 0, self.canonical, self.mode, &mut self.superkmers, &mut self.syncmer_pos, &mut self.deque, &mut self.score_buf, &mut self.mint_buf);
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
            superkmers_from_fragment(fragment, self.k, self.l, *offset, self.canonical, self.mode, &mut self.superkmers, &mut self.syncmer_pos, &mut self.deque, &mut self.score_buf, &mut self.mint_buf);
        }
        &self.superkmers
    }

    /// Access the 2-bit packed representation of the last processed sequence.
    pub fn storage(&self) -> &[u64] {
        &self.storage
    }
}

