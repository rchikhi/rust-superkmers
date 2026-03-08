# Rust Superkmers Project Guide

## Build & Test Commands
- Build: `cargo build`
- Test all: `cargo test`
- Run single test: `cargo test main`
- Benchmark: `./bench.sh` or `RUSTFLAGS='-C target-cpu=native' cargo +nightly bench`
- Format code: `cargo fmt`
- Check code: `cargo clippy`

## Code Style Guidelines
- **Naming**: Use snake_case for variables/functions, CamelCase for types/structs
- **Imports**: Group by external crates first, then internal modules
- **Error handling**: Use `expect()` with descriptive messages for file ops, `unwrap()` for trusted operations
- **Types**: Include explicit type annotations, especially for public interfaces
- **Comments**: Use doc comments (`///`) for public API, regular comments (`//`) for implementation details
- **Module structure**: Separate modules by implementation approach, main library in lib.rs
- **Function design**: Prefer small, focused functions with clear inputs/outputs
- **Performance**: For this bioinformatics project, performance is critical - benchmark different approaches

Follow standard Rust idioms and patterns from the official Rust style guide.

## Bucket Distribution Analysis (CHM13 human genome, k=31)

Measured with `bucket_stats` binary on the full CHM13v2.0 T2T human genome (3.1B kmers).

| Metric | Syncmer (l=8) | SIMD-mini (l=9) | KMC2 (l=8) | MSP (l=8) |
|--------|--------------|-----------------|------------|-----------|
| Distinct minimizers | 20,483 | 35,960 | 12,196 | 51,145 |
| Superkmers | **140M** | 148M | 268M | 294M |
| Max bucket | **16.6M** | 18.1M | 34.7M | 18.2M |
| Median bucket | 79,750 | 54,013 | 23,201 | 57 |
| Mean bucket | 152,189 | 86,683 | 255,600 | 60,950 |
| Max/Mean | **109x** | 209x | 136x | 298x |
| Max/Median | **208x** | 336x | 1,497x | 318,858x |

- **Syncmers** (`iteratorsyncmers2`, l=8): best overall — fewest superkmers, most uniform buckets.
- **SIMD-mini** (`iteratorsimdmini`, l=9): uses `simd-minimizers` crate for SIMD-accelerated
  canonical closed syncmer detection. Competitive superkmer count (148M vs 140M), more distinct
  minimizers (36K vs 20K), but worse max/mean (209x vs 109x). Requires odd l.
- **KMC2** (`iteratorkmc2`, l=8): disqualification rules backfire on repeat-rich genomes;
  concentrates k-mers into `CAG*`/`CCA*` Alu-related signatures.
- **MSP/lexicographic** (`iteratormsp`, l=8): most distinct minimizers but extreme skew.

Top offenders per method:
- Syncmer: AATGGAAT (16.6M), ATTCCATT (11.4M) — Alu-related RC pair
- SIMD-mini: AATGGAATG (18.1M), AAAAAAAAA (9.6M) — Alu + homopolymer
- KMC2: CCATTCCA (34.7M), CAGCCTGG (26.4M) — Alu consensus
- MSP: AAAAAAAT (18.2M), AAAAAAAG (14.6M) — homopolymer-adjacent

## SIMD-mini Iterator (`iteratorsimdmini`)

Ported from rust-notbcalm3's `simd-mini` feature. Uses `simd-minimizers` crate (optional dep
behind `simd-mini` feature, enabled by default).

Architecture:
1. `simd_minimizers::canonical_closed_syncmers(s=2, w=l-1)` finds all syncmer positions via SIMD
2. MSP sliding window tracks the rightmost syncmer as minimizer (sticky: only re-evaluates when
   current minimizer falls off left edge of k-mer window)

This sticky behavior matches `iteratorsyncmers2` / debruijn's `Scanner`: a new l-mer entering the
right edge only replaces the current minimizer if it has a strictly lower score. Since all syncmers
have score 0 (equal), a new syncmer never displaces an existing one — only falling off the left
triggers a change. This is correct MSP behavior, not a bug.

Constraints:
- **Requires odd l** (simd-minimizers constraint). Use l=9 instead of l=8.
- **Requires uppercase ASCII input** — the iterator handles this internally.
- `bucket_stats` auto-selects l per method (8 for syncmer/kmc2/msp, 9 for simdmini).

Bug fixed during porting: when a k-mer window has no syncmer, `curr_min_pos` was set to
`u32::MAX`. The condition `(u32::MAX as usize) < i` is never true, so the loop got stuck and
skipped the rest of the sequence. Fixed by jumping forward to the next syncmer position.
This bug also exists in rust-notbcalm3 but doesn't manifest there because it processes short
reads (150-300bp), not whole chromosomes.

## Dead Ends

- **SyncKMC2** (syncmer + KMC2 disqualification hybrid): tested with both lexicographic and
  hash-based tiebreaking. Both worse than plain syncmers. In A-rich regions all l-mers are
  KMC2-disqualified, so the filter can't redirect load. Max bucket 40.7M (vs 16.6M syncmers).
  KMC2 disqualification is counterproductive on top of syncmers — it removes valid candidates
  without providing better alternatives.

Run: `cargo +nightly run --release --bin bucket_stats <genome.fa> [k] [l] [syncmer|kmc2|msp|simdmini]`