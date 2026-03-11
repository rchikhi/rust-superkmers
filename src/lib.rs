#![feature(generic_const_exprs)]
#[allow(incomplete_features)]
pub mod utils;
pub mod naive;
pub mod iterator1;
pub mod iteratormsp;
pub mod syncmers;
pub mod iteratorsyncmersmsp;
pub mod iteratorsyncmers2;
pub mod iteratorkmc2;
#[cfg(feature = "simd-mini")]
pub mod iteratorsimdmini;
#[cfg(feature = "multi-mini")]
pub mod iteratormultiminimizers;
pub mod naivesyncmers;
use std::cmp::Ordering;
use lazy_static::lazy_static;

/// Controls how superkmer boundaries are determined.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SplitMode {
    /// Default sticky MSP: ties keep current minimizer. Longest superkmers, context-dependent.
    Sticky,
    /// Classical: split on every minimizer change (rightmost wins on tie). Context-independent.
    Classical,
    /// Syncmer-priority lexicographic: composite score (syncmer_priority, canonical_value).
    /// No ties, context-independent, fewer splits than Classical. A-rich bias.
    Msp,
    /// Syncmer-priority with XOR tiebreaker: composite score (syncmer_priority, canonical_value ^ constant).
    /// No ties, context-independent, better bucket balance than Msp.
    MspXor,
}

pub fn generate_canonical_table<const K: usize>() -> Vec<(u32, bool)> {
    let mut table = vec![(0u32, false); 1 << (2 * K)];
    for fwd in 0..(1 << (2 * K)) {
        let mut rc = 0usize;
        let mut v = fwd;
        for _ in 0..K {
            rc = (rc << 2) | (3 - (v & 3));
            v >>= 2;
        }
        table[fwd] = if rc < fwd { (rc as u32, true) } else { (fwd as u32, false) };
    }
    table
}

lazy_static! {
    pub static ref CANONICAL_8: Vec<(u32, bool)> = generate_canonical_table::<8>();
    pub static ref CANONICAL_9: Vec<(u32, bool)> = generate_canonical_table::<9>();
    pub static ref CANONICAL_10: Vec<(u32, bool)> = generate_canonical_table::<10>();
    pub static ref CANONICAL_12: Vec<(u32, bool)> = generate_canonical_table::<12>();
}

/// A superkmer: a maximal run of consecutive k-mers sharing the same minimizer.
///
/// # Fields
/// - `start` — position of the first base in the original sequence.
/// - `mint` — 2-bit packed minimizer value (canonical by default, forward-strand if
///   `.non_canonical()` is used). Encoding: A=0, C=1, G=2, T=3, MSB-first.
/// - `size` — length of the superkmer in bases (≥ k).
/// - `mpos` — relative position of the minimizer within the superkmer (0-based offset from `start`).
/// - `mint_is_rc` — `true` if the canonical minimizer is the reverse complement of the forward-strand l-mer
///   at position `start + mpos`. Always `false` when using `.non_canonical()`.
#[derive(PartialEq, Eq, Hash, Debug)]
pub struct Superkmer {
    pub start: usize,
    pub mint: u32,
    pub size: u16,
    pub mpos: u16,
    pub mint_is_rc: bool
}


impl PartialOrd for Superkmer {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Superkmer {
    fn cmp(&self, other: &Self) -> Ordering {
        self.start.cmp(&other.start)
    }
}





/// Example: extract superkmers with the syncmers2 iterator.
///
/// ```
/// let seq = b"AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCACACTGCACTGCACTGACTGCACTGCACTGCACTGCACTGCCTGC";
///
/// // Canonical mint (default) — RC-equivalent l-mers share the same bucket
/// let iter = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(seq, 21, 8);
/// for sk in iter {
///     println!("start={} mint={} size={} mpos={} mint_is_rc={}", sk.start, sk.mint, sk.size, sk.mpos, sk.mint_is_rc);
/// }
///
/// // Non-canonical — forward-strand mint, mint_is_rc always false
/// let iter = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::non_canonical(seq, 21, 8);
/// for sk in iter {
///     assert!(!sk.mint_is_rc);
/// }
/// ```
fn noop() {}
