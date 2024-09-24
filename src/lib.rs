#![feature(generic_const_exprs)]
#[allow(incomplete_features)]
pub mod utils;
pub mod naive;
pub mod iterator1;
pub mod iteratormsp;
pub mod syncmers;
pub mod iteratorsyncmersmsp;
pub mod iteratorsyncmers2;
pub mod naivesyncmers;
use std::cmp::Ordering;

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




// An iterator for getting superkmers out of a DNA sequence
///
/// Parameters:
///     * l: minimizer length
///     * k: k-mer length
///
/// Hashing is performed by NtHash
///
/// Example usage:
/// ```
///     use debruijn::dna_string::DnaString; 
///     fn main() {
///        let seq = b"AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCACACTGCACTGCACTGACTGCACTGCACTGCACTGCACTGCCTGC";
///        let iter = rust_superkmers::iterator1::SuperkmersIterator::new(seq, 10, 5);
///        for superkmer in iter
///        {
///            println!("superkmer: {:?}",superkmer);
///        }
///        let dnastring = DnaString::from_acgt_bytes(seq).to_bytes();
///        let iter = rust_superkmers::iteratormsp::SuperkmersIterator::new(&dnastring,21, 8);
///        for superkmer in iter
///        {
///            println!("superkmer: {:?}",superkmer);
///        }

///     }
/// ```
fn noop() {}
