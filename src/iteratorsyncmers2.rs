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
use crate::{Superkmer, SuperkmerParts, SplitMode};
use lazy_static::lazy_static;

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

pub(crate) fn canonical_table(l: usize) -> &'static [(u32, bool)] {
    match l {
        8 => &*crate::CANONICAL_8,
        9 => &*crate::CANONICAL_9,
        10 => &*crate::CANONICAL_10,
        12 => &*crate::CANONICAL_12,
        _ => panic!("Unsupported l={} for canonical lookup", l),
    }
}


fn generate_syncmer_scores<const K: usize>() -> Vec<usize> {
    let mut scores = vec![0usize; 1 << (2 * K)];
    for kmer_int in 0..(1 << (2 * K)) {
        let kmer_bytes = {
            let mut bytes = [0u8; K];
            for (i, byte) in bytes.iter_mut().enumerate() {
                *byte = match (kmer_int >> (2 * (K - 1 - i))) & 3 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    3 => b'T',
                    _ => unreachable!(),
                };
            }
            bytes
        };
        let syncmer = crate::syncmers::find_syncmers(K as usize, S, &[0, K - S], None, &kmer_bytes);
        scores[kmer_int] = syncmer.is_empty() as usize;
    }
    scores[0] = 1; // Demote all-A l-mer: valid syncmer but causes hot buckets
    scores[(1 << (2 * K)) - 1] = 1; // Demote all-T l-mer (RC of all-A, same canonical bucket)
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
        // Use canonical syncmer priority so fwd and rc(fwd) get the same score
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
        // Use canonical syncmer priority so fwd and rc(fwd) get the same score
        scores[fwd] = (base[canon_val as usize] << 32) | (canon_val as usize ^ XOR_CONSTANT);
    }
    scores
}

use crate::utils::bitpack_fragment;

/// Extract a single base (2 bits) from bitpacked storage at the given position.
#[inline(always)]
fn get_base(data: &[u64], pos: usize) -> usize {
    let word = pos / 32;
    let bit_offset = (31 - (pos % 32)) * 2;
    ((data[word] >> bit_offset) & 3) as usize
}

/// Extract an l-mer value from the packed storage, MSB-first encoding
/// (matching debruijn's Kmer8::to_u64() convention).
/// Each u64 in `data` holds 32 bases, with base 0 at bits 63-62 (MSB).
#[inline(always)]
pub fn get_kmer_value(data: &[u64], base_pos: usize, l: usize) -> usize {
    let word = base_pos / 32;
    let offset = base_pos % 32; // base offset within the word
    let bit_len = l * 2;

    if offset + l <= 32 {
        // All bases fit in one word
        let shift = 64 - (offset + l) * 2;
        ((data[word] >> shift) & ((1u64 << bit_len) - 1)) as usize
    } else {
        // Spans two words
        let bases_in_first = 32 - offset;
        let bases_in_second = l - bases_in_first;

        // High bits: last bases_in_first bases of the first word (lowest bits)
        let first_bits = (data[word] & ((1u64 << (bases_in_first * 2)) - 1)) as usize;

        // Low bits: first bases_in_second bases of the second word (highest bits)
        let shift = 64 - bases_in_second * 2;
        let second_bits = ((data[word + 1] >> shift) & ((1u64 << (bases_in_second * 2)) - 1)) as usize;

        (first_bits << (bases_in_second * 2)) | second_bits
    }
}

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
        SplitMode::Classical => msp_syncmer_positions_deque::<true>(storage, frag_len, k, l, offset, mode, min_positions, scores_buf, deque),
        _ => msp_syncmer_positions_deque::<false>(storage, frag_len, k, l, offset, mode, min_positions, scores_buf, deque),
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

/// Non-sticky modes (Classical, Msp, MspXor): two-pass with monotonic deque.
/// Pass 1: precompute all l-mer scores into `scores_buf`.
/// Pass 2: sliding window minimum via deque (amortized O(1) per step, no rescans).
/// CLASSICAL=true: rightmost wins on tie (pop back on >=).
/// CLASSICAL=false: leftmost wins on tie (pop back on >).
#[inline(always)]
fn msp_syncmer_positions_deque<const CLASSICAL: bool>(
    storage: &[u64], frag_len: usize, k: usize, l: usize, offset: usize,
    mode: SplitMode, min_positions: &mut Vec<(usize, usize, usize, usize)>,
    scores_buf: &mut Vec<usize>, deque: &mut Vec<usize>,
) {
    let scores = match mode {
        SplitMode::Msp => msp_syncmer_scores(l),
        SplitMode::MspXor => mspxor_syncmer_scores(l),
        _ => syncmer_scores(l),
    };

    let frag_end = offset + frag_len;
    if frag_len < k { return; }

    let num_lmers = frag_len - l + 1;
    let w = k - l + 1; // number of l-mers per k-mer window

    // Ensure buffers have enough capacity; we index by position directly.
    if scores_buf.len() < num_lmers { scores_buf.resize(num_lmers, 0); }
    if deque.len() < num_lmers { deque.resize(num_lmers, 0); }
    let sc = &mut scores_buf[..num_lmers];
    let dq = &mut deque[..num_lmers];

    // Pass 1: precompute all l-mer scores via rolling kmer
    let mask = (1usize << (l * 2)) - 1;
    let mut rolling = get_kmer_value(storage, 0, l);
    sc[0] = scores[rolling];
    for pos in 1..num_lmers {
        let new_base = get_base(storage, pos + l - 1);
        rolling = ((rolling << 2) | new_base) & mask;
        sc[pos] = scores[rolling];
    }

    // Pass 2: monotonic deque sliding window minimum
    // dq[head..tail] stores l-mer positions in increasing order, with scores
    // monotonically non-decreasing from front to back.
    // Safety: head <= tail <= num_lmers, and each element is pushed/popped at most once.
    let mut head = 0usize;
    let mut tail = 0usize;

    // Fill first window (l-mer positions 0..w-1)
    for pos in 0..w {
        while tail > head {
            let back = dq[tail - 1];
            let dominated = if CLASSICAL {
                sc[back] >= sc[pos]
            } else {
                sc[back] > sc[pos]
            };
            if dominated { tail -= 1; } else { break; }
        }
        dq[tail] = pos;
        tail += 1;
    }

    // Emit first superkmer
    let mut prev_min_pos = dq[head];
    let min_kmer = get_kmer_value(storage, prev_min_pos, l);
    min_positions.push((offset, prev_min_pos + offset, min_kmer, frag_end));

    // Slide: k-mer at position i covers l-mer positions i..i+w-1
    for i in 1..(frag_len - k + 1) {
        let new_lmer_pos = i + w - 1;

        // Pop front if it fell off the left edge
        while dq[head] < i {
            head += 1;
        }

        // Push new element, popping dominated elements from back
        while tail > head {
            let back = dq[tail - 1];
            let dominated = if CLASSICAL {
                sc[back] >= sc[new_lmer_pos]
            } else {
                sc[back] > sc[new_lmer_pos]
            };
            if dominated { tail -= 1; } else { break; }
        }
        dq[tail] = new_lmer_pos;
        tail += 1;

        let cur_min_pos = dq[head];
        if cur_min_pos != prev_min_pos {
            let min_kmer = get_kmer_value(storage, cur_min_pos, l);
            min_positions.push((i + offset, cur_min_pos + offset, min_kmer, frag_end));
            prev_min_pos = cur_min_pos;
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

/// Convert min_positions to superkmers.
fn materialize_superkmers(min_positions: &[(usize, usize, usize, usize)], k: usize, l: usize, canonical: bool, out: &mut Vec<Superkmer>) {
    let canon_table = if canonical { Some(canonical_table(l)) } else { None };
    for p in 0..min_positions.len() {
        let (start_pos, min_abs_pos, min_kmer, frag_end) = min_positions[p];
        let size = if p < min_positions.len() - 1 {
            let (next_pos, _, _, next_frag_end) = min_positions[p + 1];
            if next_frag_end == frag_end {
                next_pos + k - 1 - start_pos
            } else {
                frag_end - start_pos
            }
        } else {
            frag_end - start_pos
        };
        let (mint, mint_is_rc) = if let Some(table) = canon_table {
            table[min_kmer]
        } else {
            (min_kmer as u32, false)
        };
        out.push(Superkmer {
            start: start_pos,
            mint,
            size: size as u16,
            mpos: (min_abs_pos - start_pos) as u16,
            mint_is_rc,
        });
    }
}

/// Reverse-complement a right-aligned packed sequence of `len` bases in a u64.
/// A=0↔T=3, C=1↔G=2. Result is right-aligned.
#[inline]
fn rc_packed(val: u64, len: usize) -> u64 {
    let mut rc = 0u64;
    let mut v = val;
    for _ in 0..len {
        rc = (rc << 2) | (3 - (v & 3));
        v >>= 2;
    }
    rc
}

/// Convert min_positions to self-contained SuperkmerParts (left + canonical_mint + right).
///
/// When `canonical` is true and `mint_is_rc`, parts are pre-oriented to canonical orientation:
/// left and right are swapped and reverse-complemented so the consumer always sees
/// (left_context, canonical_minimizer, right_context) in the canonical strand.
///
/// Values are right-aligned (LSB) in u64. To MSB-align for byte-oriented output:
/// `word << (64 - len * 2)` then `.to_be_bytes()[..ceil(len * 2, 8)]`.
fn materialize_parts(storage: &[u64], min_positions: &[(usize, usize, usize, usize)], k: usize, l: usize, canonical: bool, out: &mut Vec<SuperkmerParts>) {
    let canon_table = if canonical { Some(canonical_table(l)) } else { None };
    for p in 0..min_positions.len() {
        let (start_pos, min_abs_pos, min_kmer, frag_end) = min_positions[p];
        let size = if p < min_positions.len() - 1 {
            let (next_pos, _, _, next_frag_end) = min_positions[p + 1];
            if next_frag_end == frag_end {
                next_pos + k - 1 - start_pos
            } else {
                frag_end - start_pos
            }
        } else {
            frag_end - start_pos
        };
        let mpos = min_abs_pos - start_pos;
        let (canonical_mint, mint_is_rc) = if let Some(table) = canon_table {
            table[min_kmer]
        } else {
            (min_kmer as u32, false)
        };
        let fwd_left_len = mpos;
        let fwd_right_len = size - mpos - l;
        let fwd_left = if fwd_left_len > 0 { get_kmer_value(storage, start_pos, fwd_left_len) as u64 } else { 0 };
        let fwd_right = if fwd_right_len > 0 { get_kmer_value(storage, min_abs_pos + l, fwd_right_len) as u64 } else { 0 };

        // Pre-orient: when mint_is_rc, the canonical strand is the RC of the forward strand.
        // RC reverses the whole superkmer, so left↔right swap and each gets RC'd.
        let (left, right, left_len, right_len) = if mint_is_rc {
            (rc_packed(fwd_right, fwd_right_len), rc_packed(fwd_left, fwd_left_len), fwd_right_len, fwd_left_len)
        } else {
            (fwd_left, fwd_right, fwd_left_len, fwd_right_len)
        };
        out.push(SuperkmerParts {
            canonical_mint,
            left,
            right,
            left_len: left_len as u8,
            right_len: right_len as u8,
            mint_is_rc,
        });
    }
}

/// Reusable superkmer extractor that avoids per-read allocations.
/// Create once, call `process()` or `process_with_n()` for each read.
pub struct SuperkmerExtractor {
    superkmers: Vec<Superkmer>,
    parts: Vec<SuperkmerParts>,
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
            parts: Vec::new(),
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

    /// Process a sequence (no Ns) and return self-contained SuperkmerParts.
    /// Each part contains (left, canonical_mint, right) — no external storage needed.
    pub fn process_parts(&mut self, seq: &[u8]) -> &[SuperkmerParts] {
        self.parts.clear();
        self.min_positions.clear();
        let num_words = (seq.len() + 31) / 32;
        self.storage.resize(num_words, 0);
        crate::utils::bitpack_fragment_into(seq, &mut self.storage);
        msp_syncmer_positions_into(&self.storage, seq.len(), self.k, self.l, 0, self.mode, &mut self.min_positions, &mut self.scores_buf, &mut self.deque);
        materialize_parts(&self.storage, &self.min_positions, self.k, self.l, self.canonical, &mut self.parts);
        &self.parts
    }

    /// Process a sequence (may contain Ns) and return self-contained SuperkmerParts.
    pub fn process_parts_with_n(&mut self, seq: &[u8]) -> &[SuperkmerParts] {
        self.parts.clear();
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
        materialize_parts(&self.storage, &self.min_positions, self.k, self.l, self.canonical, &mut self.parts);
        &self.parts
    }

    /// Access the 2-bit packed representation of the last processed sequence.
    pub fn storage(&self) -> &[u64] {
        &self.storage
    }
}
