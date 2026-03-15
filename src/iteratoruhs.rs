//! UHS (Universal Hitting Set) minimizer-based superkmer iterator.
//!
//! Uses purine/pyrimidine (ry) patterns from Martin Frith's UHS construction.
//! Each l-mer is classified by its ry pattern (r=purine A/G, y=pyrimidine C/T).
//! L-mers whose ry pattern matches a UHS pattern are valid minimizers (score 0);
//! non-matching l-mers get score 1. The UHS guarantees that every k=31 k-mer
//! contains at least one matching l-mer.
//!
//! Supports l=7 (21 patterns), l=8 (39 patterns), l=9 (71 patterns).
//! Uses the generic deque sliding window from `minimizer_core`.
//!
//! ```
//! let seq = b"AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCAC";
//! let iter = rust_superkmers::iteratoruhs::SuperkmersIterator::new(seq, 31, 8);
//! for sk in iter { /* canonical mint */ }
//! ```
use crate::{Superkmer, SplitMode};
use crate::minimizer_core::{canonical_table, minimizer_positions_deque, minimizer_positions_sticky, materialize_superkmers};
use crate::utils::bitpack_fragment;
use lazy_static::lazy_static;

/// UHS ry patterns for l=7 (21 patterns, density 0.1640625, sparsity 6.095).
/// Each pattern is l bits: r=0, y=1, MSB-first (leftmost letter = highest bit).
const UHS_PATTERNS_7: &[u8] = &[
    0, 7, 8, 12, 13, 18, 19,
    20, 21, 22, 23, 51, 63, 71,
    83, 85, 107, 109, 115, 119, 127,
];

/// UHS ry patterns for l=8 (39 patterns, density 0.15234375, sparsity 6.564).
const UHS_PATTERNS_8: &[u16] = &[
    0, 8, 12, 13, 18, 22, 23, 28,
    29, 30, 31, 36, 38, 39, 40, 41,
    42, 43, 68, 70, 102, 103, 150, 151,
    169, 170, 171, 198, 211, 214, 215, 219,
    221, 222, 223, 230, 231, 251, 255,
];

/// UHS ry patterns for l=9 (71 patterns, density 0.138671875, sparsity 7.211).
const UHS_PATTERNS_9: &[u16] = &[
    0, 8, 9, 12, 20, 72, 73,
    84, 86, 88, 89, 90, 92, 94,
    104, 105, 106, 108, 109, 112, 113,
    114, 115, 116, 117, 118, 120, 121,
    122, 124, 125, 126, 128, 132, 136,
    140, 148, 168, 170, 186, 187, 190,
    191, 220, 246, 268, 276, 304, 305,
    306, 307, 328, 342, 348, 350, 364,
    365, 396, 428, 429, 442, 443, 444,
    445, 446, 447, 476, 502, 508, 509,
    511,
];

/// Convert an l-mer integer to its ry pattern.
/// For 2-bit encoding (A=0, C=1, G=2, T=3): bit 0 gives ry class (A,G → 0=r; C,T → 1=y).
#[inline(always)]
fn lmer_to_ry(lmer: usize, l: usize) -> usize {
    let mut ry = 0usize;
    for i in 0..l {
        let base = (lmer >> (2 * (l - 1 - i))) & 3;
        ry = (ry << 1) | (base & 1);
    }
    ry
}

/// Check if an ry pattern is in the UHS set.
fn is_uhs_pattern(ry: usize, l: usize) -> bool {
    match l {
        7 => UHS_PATTERNS_7.contains(&(ry as u8)),
        8 => UHS_PATTERNS_8.contains(&(ry as u16)),
        9 => UHS_PATTERNS_9.contains(&(ry as u16)),
        _ => panic!("UHS not defined for l={}", l),
    }
}

/// Generate UHS binary scores: 0 = UHS member (good minimizer), 1 = non-member.
/// An l-mer is valid if either its forward OR reverse-complement ry pattern is
/// in the UHS. This ensures fwd/RC pairs get the same priority, which is required
/// for canonical minimizer selection.
fn generate_uhs_scores(l: usize) -> Vec<usize> {
    let num_lmers = 1 << (2 * l);
    let canon_table = canonical_table(l);
    let mut scores = vec![1usize; num_lmers];

    for lmer in 0..num_lmers {
        let ry_fwd = lmer_to_ry(lmer, l);
        let (canon_val, is_rc) = canon_table[lmer];
        // Get the RC l-mer integer
        let rc_lmer = if is_rc { canon_val as usize } else {
            // lmer is already canonical, find its RC
            // RC is the other one: if canon == lmer, then RC is the complement
            let mut rc = 0usize;
            let mut v = lmer;
            for _ in 0..l {
                rc = (rc << 2) | (3 - (v & 3));
                v >>= 2;
            }
            rc
        };
        let ry_rc = lmer_to_ry(rc_lmer, l);
        if is_uhs_pattern(ry_fwd, l) || is_uhs_pattern(ry_rc, l) {
            scores[lmer] = 0;
        }
    }

    // Demote homopolymers (all-A and all-T)
    scores[0] = 1;
    scores[num_lmers - 1] = 1;
    scores
}

/// Generate UHS scores with XOR tiebreaker for context-independent splitting.
/// Composite: (uhs_priority << 32) | (canonical_value ^ XOR_CONSTANT).
fn generate_uhs_mspxor_scores(l: usize) -> Vec<usize> {
    let base = generate_uhs_scores(l);
    let canon_table = canonical_table(l);
    let num_lmers = 1 << (2 * l);
    let mut scores = vec![0usize; num_lmers];
    const XOR_CONSTANT: usize = 0xACE5_ACE5;
    for fwd in 0..num_lmers {
        let (canon_val, _) = canon_table[fwd];
        scores[fwd] = (base[canon_val as usize] << 32) | (canon_val as usize ^ XOR_CONSTANT);
    }
    scores
}

lazy_static! {
    static ref UHS_SCORES_7: Vec<usize> = generate_uhs_scores(7);
    static ref UHS_SCORES_8: Vec<usize> = generate_uhs_scores(8);
    static ref UHS_SCORES_9: Vec<usize> = generate_uhs_scores(9);
    static ref UHS_MSPXOR_SCORES_7: Vec<usize> = generate_uhs_mspxor_scores(7);
    static ref UHS_MSPXOR_SCORES_8: Vec<usize> = generate_uhs_mspxor_scores(8);
    static ref UHS_MSPXOR_SCORES_9: Vec<usize> = generate_uhs_mspxor_scores(9);
}

fn uhs_scores(l: usize) -> &'static [usize] {
    match l {
        7 => &UHS_SCORES_7[..],
        8 => &UHS_SCORES_8[..],
        9 => &UHS_SCORES_9[..],
        _ => panic!("UHS not defined for l={}", l),
    }
}

fn uhs_mspxor_scores(l: usize) -> &'static [usize] {
    match l {
        7 => &UHS_MSPXOR_SCORES_7[..],
        8 => &UHS_MSPXOR_SCORES_8[..],
        9 => &UHS_MSPXOR_SCORES_9[..],
        _ => panic!("UHS not defined for l={}", l),
    }
}

/// Sliding window using scores selected by mode.
fn uhs_positions_into(
    storage: &[u64], frag_len: usize, k: usize, l: usize, offset: usize,
    mode: SplitMode, min_positions: &mut Vec<(usize, usize, usize, usize)>,
    scores_buf: &mut Vec<usize>, deque: &mut Vec<usize>,
) {
    match mode {
        SplitMode::Sticky => {
            let scores = uhs_scores(l);
            minimizer_positions_sticky(storage, frag_len, k, l, offset, scores, min_positions);
        }
        SplitMode::Classical => {
            let scores = uhs_scores(l);
            minimizer_positions_deque::<true>(storage, frag_len, k, l, offset, scores, min_positions, scores_buf, deque);
        }
        SplitMode::MspXor => {
            let scores = uhs_mspxor_scores(l);
            minimizer_positions_deque::<false>(storage, frag_len, k, l, offset, scores, min_positions, scores_buf, deque);
        }
        SplitMode::Msp => {
            let scores = uhs_scores(l);
            minimizer_positions_deque::<false>(storage, frag_len, k, l, offset, scores, min_positions, scores_buf, deque);
        }
    }
}

pub struct SuperkmersIterator {
    min_positions: Vec<(usize, usize, usize, usize)>,
    storage: Vec<u64>,
    p: usize,
    k: usize,
    l: usize,
    canonical: bool,
}

macro_rules! iter_constructors {
    ($($name:ident, $name_n:ident, $canonical:expr, $mode:expr;)*) => {
        $(
            pub fn $name(seq_str: &[u8], k: usize, l: usize) -> Self {
                Self::new_inner(seq_str, k, l, $canonical, $mode)
            }
            pub fn $name_n(seq_str: &[u8], k: usize, l: usize) -> Self {
                Self::new_with_n_inner(seq_str, k, l, $canonical, $mode)
            }
        )*
    };
}

impl SuperkmersIterator {
    iter_constructors! {
        new,                            new_with_n,                            true,  SplitMode::Sticky;
        non_canonical,                  non_canonical_with_n,                  false, SplitMode::Sticky;
        mspxor,                         mspxor_with_n,                         true,  SplitMode::MspXor;
        mspxor_non_canonical,           mspxor_non_canonical_with_n,           false, SplitMode::MspXor;
    }

    pub fn storage(&self) -> &[u64] {
        &self.storage
    }

    fn new_inner(seq_str: &[u8], k: usize, l: usize, canonical: bool, mode: SplitMode) -> Self {
        let storage = bitpack_fragment(seq_str);
        let num_lmers = seq_str.len().saturating_sub(l - 1);
        let mut min_positions = Vec::new();
        let mut scores_buf = Vec::with_capacity(num_lmers);
        let mut deque = Vec::with_capacity(num_lmers);
        uhs_positions_into(&storage, seq_str.len(), k, l, 0, mode, &mut min_positions, &mut scores_buf, &mut deque);
        SuperkmersIterator { min_positions, storage, p: 0, k, l, canonical }
    }

    fn new_with_n_inner(seq_str: &[u8], k: usize, l: usize, canonical: bool, mode: SplitMode) -> Self {
        let fragments = crate::utils::split_on_n(seq_str, k);
        let full_storage = bitpack_fragment(seq_str);
        let mut all_min_positions = Vec::new();
        let num_lmers = seq_str.len().saturating_sub(l - 1);
        let mut scores_buf = Vec::with_capacity(num_lmers);
        let mut deque = Vec::with_capacity(num_lmers);
        for (offset, fragment) in &fragments {
            let frag_storage = bitpack_fragment(fragment);
            uhs_positions_into(&frag_storage, fragment.len(), k, l, *offset, mode, &mut all_min_positions, &mut scores_buf, &mut deque);
        }
        SuperkmersIterator { min_positions: all_min_positions, storage: full_storage, p: 0, k, l, canonical }
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
                next_pos + self.k - 1 - start_pos
            } else {
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

macro_rules! extractor_constructors {
    ($($name:ident, $canonical:expr, $mode:expr;)*) => {
        $(
            pub fn $name(k: usize, l: usize) -> Self {
                Self::new_inner(k, l, $canonical, $mode)
            }
        )*
    };
}

impl SuperkmerExtractor {
    extractor_constructors! {
        new,                       true,  SplitMode::Sticky;
        non_canonical,             false, SplitMode::Sticky;
        mspxor,                    true,  SplitMode::MspXor;
        mspxor_non_canonical,      false, SplitMode::MspXor;
    }

    fn new_inner(k: usize, l: usize, canonical: bool, mode: SplitMode) -> Self {
        SuperkmerExtractor {
            superkmers: Vec::new(),
            min_positions: Vec::new(),
            storage: Vec::new(),
            frag_storage: Vec::new(),
            scores_buf: Vec::new(),
            deque: Vec::new(),
            k, l, canonical, mode,
        }
    }

    pub fn process(&mut self, seq: &[u8]) -> &[Superkmer] {
        self.superkmers.clear();
        self.min_positions.clear();
        let num_words = (seq.len() + 31) / 32;
        self.storage.resize(num_words, 0);
        crate::utils::bitpack_fragment_into(seq, &mut self.storage);
        uhs_positions_into(&self.storage, seq.len(), self.k, self.l, 0, self.mode, &mut self.min_positions, &mut self.scores_buf, &mut self.deque);
        materialize_superkmers(&self.min_positions, self.k, self.l, self.canonical, &mut self.superkmers);
        &self.superkmers
    }

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
            uhs_positions_into(&self.frag_storage, fragment.len(), self.k, self.l, *offset, self.mode, &mut self.min_positions, &mut self.scores_buf, &mut self.deque);
        }
        materialize_superkmers(&self.min_positions, self.k, self.l, self.canonical, &mut self.superkmers);
        &self.superkmers
    }

    pub fn storage(&self) -> &[u64] {
        &self.storage
    }
}
