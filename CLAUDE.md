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

## Bucket Distribution Analysis (CHM13 human genome, k=31, l=8)

Measured with `bucket_stats` binary on the full CHM13v2.0 T2T human genome (3.1B kmers).

| Metric | Syncmer | KMC2 | MSP (lex) |
|--------|---------|------|-----------|
| Distinct minimizers | 20,483 | 12,196 | 51,145 |
| Superkmers | **140M** | 268M | 294M |
| Max bucket | **16.6M** | 34.7M | 18.2M |
| Median bucket | 79,750 | 23,201 | 57 |
| Mean bucket | 152,189 | 255,600 | 60,950 |
| Max/Mean | **109x** | 136x | 298x |
| Max/Median | **208x** | 1,497x | 318,858x |

- **Syncmers** (`iteratorsyncmers2`): best overall — fewest superkmers, most uniform buckets.
- **KMC2** (`iteratorkmc2`): disqualification rules backfire on repeat-rich genomes; concentrates k-mers into `CAG*`/`CCA*` Alu-related signatures. 2x worse max bucket than syncmers.
- **MSP/lexicographic** (`iteratormsp`): most distinct minimizers but extreme skew; half are nearly empty while `AAAA*` signatures are massive.

Top offenders per method:
- Syncmer: AATGGAAT (16.6M), ATTCCATT (11.4M) — Alu-related RC pair
- KMC2: CCATTCCA (34.7M), CAGCCTGG (26.4M) — Alu consensus
- MSP: AAAAAAAT (18.2M), AAAAAAAG (14.6M) — homopolymer-adjacent

Run: `cargo +nightly run --release --bin bucket_stats <genome.fa> [k] [l] [syncmer|kmc2|msp]`