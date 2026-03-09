# Rust-superkmers

A little investigation into constructing super-kmers using syncmers for minimum substring partitioning in Rust.

## Speed considerations

Benchmarked on random sequences, k=31:

| Method | 150bp | 10Kbp | 1Mbp |
|--------|------:|------:|-----:|
| `iteratorsimdmini` (l=9) | 213 MB/s | 367 MB/s | 253 MB/s |
| `iteratorsyncmersmsp` (l=8) | 202 MB/s | 186 MB/s | 186 MB/s |
| `iteratorsyncmers2` (l=8) | 179 MB/s | 155 MB/s | 114 MB/s |
| `msp` (l=8) | 175 MB/s | 168 MB/s | 157 MB/s |
| `kmc2` (l=8) | 168 MB/s | 136 MB/s | 109 MB/s |

- `iteratorsimdmini` uses the `simd-minimizers` crate for SIMD-accelerated canonical closed syncmer detection with a sticky MSP sliding window. Requires odd l.
- `iteratorsyncmers2` uses AVX2 for sequence bit-packing and reimplements the MSP sliding window with syncmer lookup tables.
- `iteratorsyncmersmsp` wraps the debruijn Scanner with syncmer scoring. Produces identical output to `iteratorsyncmers2`.
- `msp` wraps the debruijn Scanner with lexicographic scoring.
- `kmc2` uses KMC2 disqualification rules for minimizer selection.
