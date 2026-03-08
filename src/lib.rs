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
pub mod naivesyncmers;
use std::cmp::Ordering;
use lazy_static::lazy_static;

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
/// - `rc` — `true` if the canonical minimizer is the reverse complement of the forward-strand l-mer
///   at position `start + mpos`. Always `false` when using `.non_canonical()`.
#[derive(PartialEq, Eq, Hash, Debug)]
pub struct Superkmer {
    pub start: usize,
    pub mint: u32,
    pub size: u8,
    pub mpos: u8,
    pub rc: bool
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

#[derive(PartialEq, Eq, Hash, Debug)]
pub struct SuperkmerVerbose {
    pub sequence: String,
    pub minimizer: String,
    pub mpos: usize,
}


impl PartialOrd for SuperkmerVerbose {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for SuperkmerVerbose {
    fn cmp(&self, other: &Self) -> Ordering {
        self.mpos.cmp(&other.mpos)
    }
}




/// Example: extract superkmers with the syncmers2 iterator.
///
/// ```
/// let seq = b"AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCACACTGCACTGCACTGACTGCACTGCACTGCACTGCACTGCCTGC";
///
/// // Canonical mint (default) — RC-equivalent l-mers share the same bucket
/// let (_, iter) = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(seq, 21, 8);
/// for sk in iter {
///     println!("start={} mint={} size={} mpos={} rc={}", sk.start, sk.mint, sk.size, sk.mpos, sk.rc);
/// }
///
/// // Non-canonical — forward-strand mint, rc always false
/// let (_, iter) = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(seq, 21, 8);
/// for sk in iter.non_canonical() {
///     assert!(!sk.rc);
/// }
/// ```
fn noop() {}
