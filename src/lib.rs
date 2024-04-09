pub mod utils;
pub mod naive;
pub mod iterator1;
pub mod iteratormsp;

#[derive(PartialEq, Eq, Hash, Debug)]
pub struct Superkmer {
    pub start: usize,
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
///     fn main() {
///        let seq = b"AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCACACTGCACTGCACTGACTGCACTGCACTGCACTGCACTGCCTGC";
///        let iter = rust_superkmers::iterator1::SuperkmersIterator::new(seq, 10, 5);
///        for superkmer in iter
///        {
///            println!("superkmer: {:?}",superkmer);
///        }
///        let iter = rust_superkmers::iteratormsp::SuperkmersIterator::new(seq, 21, 8);
///        for superkmer in iter
///        {
///            println!("superkmer: {:?}",superkmer);
///        }

///     }
/// ```
fn test() { }


