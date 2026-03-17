//! AVX2-accelerated syncmer-based superkmer iterator.
//!
//! Uses 2-bit packed storage and a precomputed syncmer score table for fast MSP.
//! Supports l=8 and l=9. Mint is canonical by default; call `.non_canonical()`
//! on the iterator to get forward-strand mint instead.
//!
//! ```
//! let seq = b"AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCAC";
//! let iter = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(seq, 21, 8);
//! for sk in iter { /* canonical mint */ }
//!
//! let iter = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::non_canonical(seq, 21, 8);
//! for sk in iter { /* forward-strand mint */ }
//! ```
use crate::{Superkmer, SplitMode};
use crate::minimizer_core::{canonical_table, minimizer_positions_deque, materialize_superkmers};
use lazy_static::lazy_static;

// Re-export for public API compatibility (bench_iterators uses these)
pub use crate::minimizer_core::{get_kmer_value, get_base};

const S: usize = 2; //syncmer's s parameter

lazy_static! {
static ref SYNCMER_SCORES_8: Vec<usize> = generate_syncmer_scores::<8>();
static ref SYNCMER_SCORES_9: Vec<usize> = generate_syncmer_scores::<9>();
static ref MSP_SYNCMER_SCORES_8: Vec<usize> = generate_msp_syncmer_scores::<8>();
static ref MSP_SYNCMER_SCORES_9: Vec<usize> = generate_msp_syncmer_scores::<9>();
static ref MSPXOR_SYNCMER_SCORES_8: Vec<usize> = generate_mspxor_syncmer_scores::<8>();
static ref MSPXOR_SYNCMER_SCORES_9: Vec<usize> = generate_mspxor_syncmer_scores::<9>();
}

pub(crate) fn syncmer_scores(l: usize) -> &'static [usize] {
    match l {
        8 => &SYNCMER_SCORES_8[..],
        9 => &SYNCMER_SCORES_9[..],
        _ => panic!("Unsupported l={} for syncmer scores", l),
    }
}

pub(crate) fn msp_syncmer_scores(l: usize) -> &'static [usize] {
    match l {
        8 => &MSP_SYNCMER_SCORES_8[..],
        9 => &MSP_SYNCMER_SCORES_9[..],
        _ => panic!("Unsupported l={} for MSP syncmer scores", l),
    }
}

pub fn mspxor_syncmer_scores(l: usize) -> &'static [usize] {
    match l {
        8 => &MSPXOR_SYNCMER_SCORES_8[..],
        9 => &MSPXOR_SYNCMER_SCORES_9[..],
        _ => panic!("Unsupported l={} for MSP-xor syncmer scores", l),
    }
}



fn generate_syncmer_scores<const K: usize>() -> Vec<usize> {
    generate_syncmer_scores_with_s(K, S)
}

/// Generate syncmer scores for arbitrary (l, s) parameters.
pub fn generate_syncmer_scores_with_s(l: usize, s: usize) -> Vec<usize> {
    let num_lmers = 1 << (2 * l);
    let mut scores = vec![0usize; num_lmers];
    let mut kmer_bytes = vec![0u8; l];
    for kmer_int in 0..num_lmers {
        for (i, byte) in kmer_bytes.iter_mut().enumerate() {
            *byte = match (kmer_int >> (2 * (l - 1 - i))) & 3 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => unreachable!(),
            };
        }
        let syncmer = crate::syncmers::find_syncmers(l, s, &[0, l - s], None, &kmer_bytes);
        scores[kmer_int] = syncmer.is_empty() as usize;
    }
    scores[0] = 1; // Demote all-A l-mer
    scores[num_lmers - 1] = 1; // Demote all-T l-mer
    scores
}

/// MSP scores: composite (syncmer_priority << 32 | canonical_value).
/// Syncmers sort before non-syncmers; within each group, lower canonical value wins.
/// Scores are unique per canonical l-mer, so ties only occur between forward/RC
/// pairs (same bucket), making the sticky loop context-independent.
fn generate_msp_syncmer_scores<const K: usize>() -> Vec<usize> {
    let base = generate_syncmer_scores::<K>();
    let canon_table = canonical_table(K);
    let mut scores = vec![0usize; 1 << (2 * K)];
    for fwd in 0..(1 << (2 * K)) {
        let (canon_val, _) = canon_table[fwd];
        scores[fwd] = (base[canon_val as usize] << 32) | canon_val as usize;
    }
    scores
}

/// Random constant for XOR tiebreaker. Breaks lexicographic ordering without a full hash.
const XOR_CONSTANT: usize = 0xACE5_ACE5;

/// MSP-xor scores: composite (syncmer_priority << 32 | (canonical_value ^ XOR_CONSTANT)).
/// Like MSP but XORs a constant to break A-rich lexicographic bias.
fn generate_mspxor_syncmer_scores<const K: usize>() -> Vec<usize> {
    let base = generate_syncmer_scores::<K>();
    let canon_table = canonical_table(K);
    let mut scores = vec![0usize; 1 << (2 * K)];
    for fwd in 0..(1 << (2 * K)) {
        let (canon_val, _) = canon_table[fwd];
        scores[fwd] = (base[canon_val as usize] << 32) | (canon_val as usize ^ XOR_CONSTANT);
    }
    scores
}

/// Generate mspxor scores with arbitrary (l, s) parameters.
pub fn generate_mspxor_syncmer_scores_with_s(l: usize, s: usize) -> Vec<usize> {
    let base = generate_syncmer_scores_with_s(l, s);
    let canon_table = canonical_table(l);
    let num_lmers = 1 << (2 * l);
    let mut scores = vec![0usize; num_lmers];
    for fwd in 0..num_lmers {
        let (canon_val, _) = canon_table[fwd];
        scores[fwd] = (base[canon_val as usize] << 32) | (canon_val as usize ^ XOR_CONSTANT);
    }
    scores
}


use crate::utils::bitpack_fragment;

/// MinPos tracks the best minimizer position within a k-mer window.
/// Ordering: lower val wins; on tie, higher pos (rightmost) wins.
/// This matches debruijn's msp.rs MinPos exactly.
#[derive(Clone, Copy, PartialEq, Eq)]
struct MinPos {
    val: usize,
    pos: usize,
    kmer: usize,
}

impl Ord for MinPos {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let val_cmp = self.val.cmp(&other.val);
        if val_cmp != std::cmp::Ordering::Equal {
            return val_cmp;
        }
        // Reverse position: higher pos is "less" (preferred by min)
        match self.pos.cmp(&other.pos) {
            std::cmp::Ordering::Equal => std::cmp::Ordering::Equal,
            std::cmp::Ordering::Less => std::cmp::Ordering::Greater,
            std::cmp::Ordering::Greater => std::cmp::Ordering::Less,
        }
    }
}

impl PartialOrd for MinPos {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}


/// Run MSP sliding window on a single fragment, returning (kmer_start, minimizer_pos, minimizer_kmer, fragment_end).
/// All positions are absolute (offset already added).
/// `scores_buf` and `deque` are reusable scratch buffers (unused for Sticky mode).
fn msp_syncmer_positions_into(storage: &[u64], frag_len: usize, k: usize, l: usize, offset: usize, mode: SplitMode, min_positions: &mut Vec<(usize, usize, usize, usize)>, scores_buf: &mut Vec<usize>, deque: &mut Vec<usize>) {
    match mode {
        SplitMode::Sticky => msp_syncmer_positions_sticky(storage, frag_len, k, l, offset, min_positions),
        SplitMode::Classical => {
            let scores = syncmer_scores(l);
            minimizer_positions_deque::<true>(storage, frag_len, k, l, offset, scores, min_positions, scores_buf, deque);
        }
        SplitMode::Msp => {
            let scores = msp_syncmer_scores(l);
            minimizer_positions_deque::<false>(storage, frag_len, k, l, offset, scores, min_positions, scores_buf, deque);
        }
        SplitMode::MspXor => {
            let scores = mspxor_syncmer_scores(l);
            minimizer_positions_deque::<false>(storage, frag_len, k, l, offset, scores, min_positions, scores_buf, deque);
        }
    }
}

/// Sticky mode: single-pass with rescan on falloff. Rescans are rare because
/// equal-score l-mers (all syncmers have score 0) never replace the current minimizer.
#[inline(always)]
fn msp_syncmer_positions_sticky(storage: &[u64], frag_len: usize, k: usize, l: usize, offset: usize, min_positions: &mut Vec<(usize, usize, usize, usize)>) {
    let scores = syncmer_scores(l);
    let mp = |pos: usize| -> MinPos {
        let kmer = get_kmer_value(storage, pos, l);
        let val = scores[kmer];
        MinPos { val, pos, kmer }
    };

    let find_min = |start: usize, stop: usize| -> MinPos {
        let mut min_pos = mp(start);
        for pos in (start + 1)..=stop {
            let current = mp(pos);
            min_pos = std::cmp::min(min_pos, current);
        }
        min_pos
    };

    let frag_end = offset + frag_len;

    if frag_len >= k {
        let mut min_pos = find_min(0, k - l);
        min_positions.push((offset, min_pos.pos + offset, min_pos.kmer, frag_end));

        let mask = (1usize << (l * 2)) - 1;
        let mut rolling_kmer = get_kmer_value(storage, k - l, l);

        for i in 1..(frag_len - k + 1) {
            let new_base = get_base(storage, i + k - 1);
            rolling_kmer = ((rolling_kmer << 2) | new_base) & mask;
            let end_val = scores[rolling_kmer];
            let end_pos = MinPos { val: end_val, pos: i + k - l, kmer: rolling_kmer };

            if i > min_pos.pos {
                min_pos = find_min(i, i + k - l);
                min_positions.push((i + offset, min_pos.pos + offset, min_pos.kmer, frag_end));
            } else if end_pos.val < min_pos.val {
                min_pos = end_pos;
                min_positions.push((i + offset, min_pos.pos + offset, min_pos.kmer, frag_end));
            }
        }
    }
}

pub struct SuperkmersIterator {
    min_positions: Vec<(usize, usize, usize, usize)>, // (kmer_start, minimizer_pos, minimizer_kmer, fragment_end)
    storage: Vec<u64>,
    p: usize,
    k: usize,
    l: usize,
    canonical: bool,
}

/// Generate iterator constructors for each (mode, canonical) combination.
macro_rules! iter_constructors {
    ($($name:ident, $name_n:ident, $canonical:expr, $mode:expr;)*) => {
        $(
            pub fn $name(seq_str: &[u8], k: usize, l: usize) -> Self {
                Self::new_inner_full(seq_str, k, l, $canonical, $mode)
            }
            pub fn $name_n(seq_str: &[u8], k: usize, l: usize) -> Self {
                Self::new_with_n_inner_full(seq_str, k, l, $canonical, $mode)
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

    fn new_inner_full(seq_str: &[u8], k: usize, l: usize, canonical: bool, mode: SplitMode) -> Self {
        let storage = bitpack_fragment(seq_str);
        let num_lmers = seq_str.len().saturating_sub(l - 1);
        let mut min_positions = Vec::new();
        let mut scores_buf = Vec::with_capacity(num_lmers);
        let mut deque = Vec::with_capacity(num_lmers);
        msp_syncmer_positions_into(&storage, seq_str.len(), k, l, 0, mode, &mut min_positions, &mut scores_buf, &mut deque);

        SuperkmersIterator {
            min_positions,
            storage,
            p: 0,
            k,
            l,
            canonical,
        }
    }

    fn new_with_n_inner_full(seq_str: &[u8], k: usize, l: usize, canonical: bool, mode: SplitMode) -> Self {
        let fragments = crate::utils::split_on_n(seq_str, k);
        let full_storage = bitpack_fragment(seq_str);

        let mut all_min_positions: Vec<(usize, usize, usize, usize)> = Vec::new();
        let num_lmers = seq_str.len().saturating_sub(l - 1);
        let mut scores_buf = Vec::with_capacity(num_lmers);
        let mut deque = Vec::with_capacity(num_lmers);

        for (offset, fragment) in &fragments {
            let frag_storage = bitpack_fragment(fragment);
            msp_syncmer_positions_into(&frag_storage, fragment.len(), k, l, *offset, mode, &mut all_min_positions, &mut scores_buf, &mut deque);
        }

        SuperkmersIterator {
            min_positions: all_min_positions,
            storage: full_storage,
            p: 0,
            k,
            l,
            canonical,
        }
    }
}

impl Iterator for SuperkmersIterator {
    type Item = Superkmer;

    fn next(&mut self) -> Option<Self::Item> {
        if self.p >= self.min_positions.len() {
            return None;
        }

        let (start_pos, min_abs_pos, min_kmer, frag_end) = self.min_positions[self.p];

        let size = if self.p < self.min_positions.len() - 1 {
            let (next_pos, _, _, next_frag_end) = self.min_positions[self.p + 1];
            if next_frag_end == frag_end {
                // Same fragment: superkmer extends to overlap with next k-mer
                next_pos + self.k - 1 - start_pos
            } else {
                // Last superkmer in this fragment
                frag_end - start_pos
            }
        } else {
            frag_end - start_pos
        };

        self.p += 1;

        let (mint, mint_is_rc) = if self.canonical {
            canonical_table(self.l)[min_kmer]
        } else {
            (min_kmer as u32, false)
        };
        Some(Superkmer {
            start: start_pos,
            mint,
            size: size as u16,
            mpos: (min_abs_pos - start_pos) as u16,
            mint_is_rc,
        })
    }
}

/// Reusable superkmer extractor that avoids per-read allocations.
/// Create once, call `process()` or `process_with_n()` for each read.
pub struct SuperkmerExtractor {
    superkmers: Vec<Superkmer>,
    min_positions: Vec<(usize, usize, usize, usize)>,
    storage: Vec<u64>,
    frag_storage: Vec<u64>,
    scores_buf: Vec<usize>,
    deque: Vec<usize>,
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
            min_positions: Vec::new(),
            storage: Vec::new(),
            frag_storage: Vec::new(),
            scores_buf: Vec::new(),
            deque: Vec::new(),
            k,
            l,
            canonical,
            mode,
        }
    }

    /// Process a sequence with no N characters. Returns the superkmers slice.
    pub fn process(&mut self, seq: &[u8]) -> &[Superkmer] {
        self.superkmers.clear();
        self.min_positions.clear();
        let num_words = (seq.len() + 31) / 32;
        self.storage.resize(num_words, 0);
        crate::utils::bitpack_fragment_into(seq, &mut self.storage);
        msp_syncmer_positions_into(&self.storage, seq.len(), self.k, self.l, 0, self.mode, &mut self.min_positions, &mut self.scores_buf, &mut self.deque);
        materialize_superkmers(&self.min_positions, self.k, self.l, self.canonical, &mut self.superkmers);
        &self.superkmers
    }

    /// Process a sequence that may contain N/n characters. Returns the superkmers slice.
    pub fn process_with_n(&mut self, seq: &[u8]) -> &[Superkmer] {
        self.superkmers.clear();
        self.min_positions.clear();
        let num_words = (seq.len() + 31) / 32;
        self.storage.resize(num_words, 0);
        crate::utils::bitpack_fragment_into(seq, &mut self.storage);
        let fragments = crate::utils::split_on_n(seq, self.k);
        for (offset, fragment) in &fragments {
            let frag_words = (fragment.len() + 31) / 32;
            self.frag_storage.resize(frag_words, 0);
            crate::utils::bitpack_fragment_into(fragment, &mut self.frag_storage);
            msp_syncmer_positions_into(&self.frag_storage, fragment.len(), self.k, self.l, *offset, self.mode, &mut self.min_positions, &mut self.scores_buf, &mut self.deque);
        }
        materialize_superkmers(&self.min_positions, self.k, self.l, self.canonical, &mut self.superkmers);
        &self.superkmers
    }


    /// Access the 2-bit packed representation of the last processed sequence.
    pub fn storage(&self) -> &[u64] {
        &self.storage
    }
}
