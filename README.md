# Rust-superkmers

A little investigation into constructing super-kmers using syncmers for minimum substring partitioning in Rust.

## Speed considerations

Benchmarked on 150bp random sequences, k=31:

* The Rust Nthash iterator from Luiz Irber runs at ~380 MB/s (it doesn't compute superkmers).
* The lexicographic minimizer version of `rust-debruijn` MSP from 10XGenomics reaches ~285 MB/s.
* From this crate:
  * `iteratorsimdmini` (l=9), which uses the `simd-minimizers` crate for SIMD-accelerated canonical closed syncmer detection with a sticky MSP sliding window, runs at ~142 MB/s. Requires odd l.
  * `iteratorsyncmers2` (l=8), which uses AVX2 for sequence bit-packing and reimplements the MSP sliding window with syncmer lookup tables, runs at ~134 MB/s. Produces identical output to `iteratorsyncmersmsp`.
  * `iteratorsyncmersmsp` (l=8), which wraps the debruijn Scanner with syncmer scoring, runs at ~85 MB/s
  * `iterator1`, which 100% matches the result of naive, runs at ~102 MB/s
  * `naive`, which uses nthash but recomputes minimizers for each kmer, runs at ~32 MB/s

Note: `iteratorsimdmini` scales significantly better on longer sequences (~211 MB/s at 1Mb vs ~107 MB/s for `iteratorsyncmers2`).
