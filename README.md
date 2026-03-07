# Rust-superkmers

A little investigation into constructing super-kmers using syncmers for minimum substring partitioning in Rust.

## Speed considerations

Benchmarked on 150bp random sequences, k=31, l=8:

* The Rust Nthash iterator from Luiz Irber runs at ~380 MB/s (it doesn't compute superkmers).
* The lexicographic minimizer version of `rust-debruijn` MSP from 10XGenomics reaches ~285 MB/s.
* From this crate:
  * `iteratorsyncmers2`, which uses AVX2 for sequence bit-packing and reimplements the MSP sliding window with syncmer lookup tables, runs at ~153 MB/s. Produces identical output to `iteratorsyncmersmsp`.
  * `iteratorsyncmersmsp`, which wraps the debruijn Scanner with syncmer scoring, runs at ~145 MB/s
  * `iterator1`, which 100% matches the result of naive, runs at ~102 MB/s
  * `naive`, which uses nthash but recomputes minimizers for each kmer, runs at ~32 MB/s
