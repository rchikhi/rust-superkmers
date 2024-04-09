# Rust-superkmers

A little investigation into constructing super-kmers, also known as minimum substring partitioning, in Rust.

## Speed considerations

* The Rust Nthash crate from Luiz Irber runs at around 180MB/s on my machine (it doesn't compute superkmers).
* The lexicographic minimizer version of `rust-debruijn` MSP from 10XGenomics reaches 200MB/s.
* From this crate:
  * `naive`, which uses nthash, but recomputes minimizers for each kmer, runs at 16 MB/s
  * `iterator1`, which 100% matches the result of naive, runs at 40MB/s.
  
