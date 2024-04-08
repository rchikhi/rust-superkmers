#![feature(core_intrinsics)]
#![feature(stdarch_x86_avx512)]

use std::io::{Result};

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


