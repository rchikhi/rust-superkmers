#![feature(generic_const_exprs)]
pub mod utils;
pub mod naive;
pub mod iterator1;
pub mod iteratormsp;
pub mod syncmers;
pub mod iteratorsyncmers;

#[derive(PartialEq, Eq, Hash, Debug)]
pub struct Superkmer {
    pub start: usize,
    pub mint: u32,
    pub size: u8,
    pub mpos: u8,
    pub rc: bool
}

#[derive(PartialEq, Eq, Hash, Debug)]
pub struct SuperkmerVerbose {
    pub sequence: String,
    pub minimizer: String,
    pub mpos: usize,
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
///        let dnastring = DnaString::from_acgt_bytes(seq);
///        let iter = rust_superkmers::iteratormsp::SuperkmersIterator::new(&dnastring, 21, 8);
///        for superkmer in iter
///        {
///            println!("superkmer: {:?}",superkmer);
///        }

///     }
/// ```
fn noop() {}
