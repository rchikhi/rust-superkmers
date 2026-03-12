# Rust-superkmers

A little investigation into constructing super-kmers using syncmers for minimum substring partitioning in Rust.

## Speed considerations

Benchmarked on 150bp random sequences, k=31:

* The Rust Nthash iterator from Luiz Irber runs at ~380 MB/s (it doesn't compute superkmers).
* The lexicographic minimizer version of `rust-debruijn` MSP from 10XGenomics reaches ~285 MB/s.
* From this crate (sticky mode, context-dependent):
  * `iteratorsimdmini` (l=9), SIMD-accelerated canonical closed syncmer detection, runs at ~165 MB/s. Requires odd l.
  * `iteratorsyncmersmsp` (l=8), wraps the debruijn Scanner with syncmer scoring, runs at ~150 MB/s
  * `iteratorsyncmers2` (l=8), AVX2 bit-packing with syncmer lookup tables, runs at ~174 MB/s. Produces identical output to `iteratorsyncmersmsp`.
  * `iterator1`, which 100% matches the result of naive, runs at ~95 MB/s
  * `naive`, which uses nthash but recomputes minimizers for each kmer, runs at ~28 MB/s
* MspXor mode (context-independent, same k-mer always maps to same bucket):
  * `iteratorsyncmers2:mspxor` (l=8) runs at ~148 MB/s.
  * `iteratorsimdmini:mspxor` (l=9) runs at ~94 MB/s.

Note: `iteratorsimdmini` scales significantly better on longer sequences (~307 MB/s at 1Mb vs ~136 MB/s for `iteratorsyncmers2`).

## API

### One-shot iteration (SuperkmersIterator)

```rust
use rust_superkmers::iteratorsyncmers2::SuperkmersIterator;

// Default: sticky MSP, canonical minimizers
for sk in SuperkmersIterator::new(seq, 21, 8) {
    println!("mint={} size={}", sk.mint, sk.size);
}

// Context-independent MspXor mode
for sk in SuperkmersIterator::mspxor(seq, 21, 8) {
    println!("mint={} size={}", sk.mint, sk.size);
}

// Non-canonical (forward-strand) minimizers
for sk in SuperkmersIterator::mspxor_non_canonical(seq, 21, 8) {
    assert!(!sk.mint_is_rc);
}
```

### Reusable extractor (SuperkmerExtractor) — allocation-free repeated processing

```rust
use rust_superkmers::iteratorsyncmers2::SuperkmerExtractor;

// Create once, reuse across many sequences
let mut extractor = SuperkmerExtractor::mspxor(31, 8);

for seq in sequences {
    let superkmers = extractor.process(seq);
    for sk in superkmers {
        println!("mint={} size={}", sk.mint, sk.size);
    }
}
```

Available constructors for both types: `new`, `non_canonical`, `classical`, `classical_non_canonical`, `msp`, `msp_non_canonical`, `mspxor`, `mspxor_non_canonical`. Add `_with_n` suffix (iterator only) for sequences containing N characters.
