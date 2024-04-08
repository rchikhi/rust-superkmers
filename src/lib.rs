mod utils;
pub mod naive;
pub mod iterator1;

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
///        let iter = SuperkmersIterator::new(seq, 10, 5, ).unwrap();
///        for superkmer in iter
///        {
///            println!("superkmer: {:?}",superkmer);
///        }
///     }
/// ```
fn test() { }


