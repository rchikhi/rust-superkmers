# Rust-superkmers

A little investigation into constructing super-kmers using syncmers for minimum substring partitioning in Rust.

## Speed considerations

Benchmarked on 150bp random sequences, k=31:

* The Rust Nthash iterator from Luiz Irber runs at ~380 MB/s (it doesn't compute superkmers).
* The lexicographic minimizer version of `rust-debruijn` MSP from 10XGenomics reaches ~285 MB/s.
* From this crate (sticky mode, context-dependent):
  * `iteratorsimdmini` (l=9), SIMD-accelerated canonical closed syncmer detection, runs at ~239 MB/s. Requires odd l.
  * `iteratorsyncmersmsp` (l=8), wraps the debruijn Scanner with syncmer scoring, runs at ~186 MB/s
  * `iteratorsyncmers2` (l=8), AVX2 bit-packing with syncmer lookup tables, runs at ~172 MB/s. Produces identical output to `iteratorsyncmersmsp`.
  * `iterator1`, which 100% matches the result of naive, runs at ~95 MB/s
  * `naive`, which uses nthash but recomputes minimizers for each kmer, runs at ~28 MB/s
* MspXor mode (context-independent, same k-mer always maps to same bucket):
  * `iteratorsyncmers2:mspxor` (l=8) runs at ~186 MB/s — faster than sticky thanks to unique scores reducing window rescans.
  * `iteratorsimdmini:mspxor` (l=9) runs at ~81 MB/s — slower due to linear scan in `find_best` (needs monotonic deque optimization).

Note: `iteratorsimdmini` scales significantly better on longer sequences (~333 MB/s at 1Mb vs ~106 MB/s for `iteratorsyncmers2`). `iteratorsyncmers2:mspxor` also scales well (~138 MB/s at 1Mb, faster than sticky's 106 MB/s).

## Bucket Distribution (CHM13 human genome, k=31)

Measured with `bucket_stats` on the full CHM13v2.0 T2T human genome (~3.1B k-mers).

### Split modes

Each iterator supports multiple **split modes** that control how superkmer boundaries are determined:

- **Sticky** (default): ties keep the current minimizer. Longest superkmers, but context-dependent (same k-mer may get different minimizers depending on neighbors).
- **Classical**: rightmost wins on ties. Context-independent but produces many more superkmers because all syncmers have equal score.
- **Msp**: composite score `(syncmer_priority, canonical_value)`. No ties, context-independent. Has A-rich lexicographic bias.
- **MspXor**: composite score `(syncmer_priority, canonical_value ^ 0xACE5ACE5)`. No ties, context-independent, best bucket balance.

### Comparison

| Method | l | Superkmers | Distinct min. | Max bucket | Max/Mean | Context-indep? |
|--------|---|-----------|--------------|-----------|----------|----------------|
| Syncmer sticky | 8 | **140M** | 20,483 | 28.0M | 166x | No |
| SIMD-mini sticky | 9 | 148M | 35,960 | 18.1M | 210x | No |
| Multi-mini N=4 | 9 | 179M | 115,040 | 5.5M | 171x | No |
| **Syncmer:mspxor** | **8** | **246M** | **16,920** | **11.5M** | **63x** | **Yes** |
| SIMD-mini:mspxor | 9 | 258M | 32,166 | 16.4M | 170x | Yes |
| SIMD-mini:msp | 9 | 254M | 34,428 | 28.3M | 303x | Yes |
| Syncmer:msp | 8 | 274M | 19,484 | 35.8M | 179x | Yes |
| MSP/lex | 8 | 294M | 51,145 | 18.2M | 298x | Yes |
| SIMD-mini:classical | 9 | 954M | 35,960 | 18.2M | 210x | Yes |
| Syncmer:classical | 8 | 1,021M | 20,483 | 28.3M | 167x | Yes |

**Syncmer:mspxor** is the clear winner for context-independent minimizers: 63x max/mean (best bucket balance of any method), 11.5M max bucket (less than half of sticky's 28M), and 246M superkmers (reasonable 1.75x overhead vs sticky).

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
