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

## Iterator Implementations

Each module implements a different minimizer scheme for superkmer extraction.

| Module | Scheme | Score table | Split modes | l constraints | Feature flag | API |
|--------|--------|------------|-------------|---------------|-------------|-----|
| `iteratorsyncmers2` | Closed syncmers (s=2) | `Vec<u32>`, l=8,9 | Sticky, Classical, Msp, MspXor | 8, 9 | — | Iterator + Extractor |
| `iteratorsimdmini` | SIMD closed syncmers (s=2) | None (inline) | Sticky, Classical, Msp, MspXor | odd l, any size | `simd-mini` (default) | Iterator + Extractor |
| `iteratorsimdmini_cminim` | ntHash random minimizers | None (SIMD) | single (internal) | odd k | `simd-mini` (default) | Iterator |
| `iteratoruhs` | UHS ry-alphabet patterns | `Vec<u32>`, l=7,8,9 | Sticky, Classical, Msp, MspXor | 7, 8, 9 | — | Iterator + Extractor |
| `iteratorsyncmersmsp` | Syncmers via debruijn Scanner | `Vec<u32>`, l=8,10,12 | single (debruijn) | 8, 10, 12 | — | Iterator |
| `iteratormsp` | Lexicographic (canonical value) | via debruijn | single (debruijn) | 8, 10, 12 | — | Iterator |
| `iteratorkmc2` | KMC2 disqualification | precomputed u64 | Sticky only | 8 | — | Iterator |
| `iterator1` | ntHash | None (VecDeque) | Sticky only | any | — | Iterator |
| `iteratormultiminimizers` | N independent random hashes | None (multiminimizers) | single (internal) | k-l even | `multi-mini` | Iterator |

**Extractor** (`SuperkmerExtractor`): reusable across sequences, zero allocation on repeated `.process()` calls. Preferred for multi-read pipelines.

**Iterator** (`SuperkmersIterator`): one-shot, loads entire sequence. Simpler API.

Key modules:
- `minimizer_core`: shared sliding window (block-decomposition + sticky), generic over `Score` trait (u16/u32)
- `syncmers`: syncmer detection utility (`find_syncmers`)
- `utils`: AVX2/scalar bitpacking, `split_on_n`

## Bucket Distribution Analysis (CHM13 human genome, k=31)

Measured with `bucket_stats` binary on the full CHM13v2.0 T2T human genome (3.1B kmers).

### Sticky mode (default)

| Metric | Syncmer (l=8) | SIMD-mini (l=9) | UHS (l=8) | UHS (l=9) | KMC2 (l=8) | MSP (l=8) |
|--------|--------------|-----------------|-----------|-----------|------------|-----------|
| Distinct minimizers | 18,420 | 35,960 | 11,092 | 27,456 | 12,196 | 51,145 |
| Superkmers | **140M** | 148M | 144M | 158M | 268M | 294M |
| Max bucket | 28.0M | 18.1M | **11.7M** | 20.5M | 34.7M | 18.2M |
| Median bucket | 84,998 | 53,982 | 113,558 | 32,523 | 23,201 | 57 |
| Mean bucket | 169,234 | 86,683 | 281,040 | 113,538 | 255,600 | 60,950 |
| Max/Mean | 165x | 209x | **42x** | 181x | 136x | 298x |
| Max/Median | 329x | 336x | **103x** | 631x | 1,497x | 318,858x |

### MspXor mode (context-independent)

| Metric | Syncmer (l=8) | SIMD-mini (l=9) | UHS (l=8) | UHS (l=9) |
|--------|--------------|-----------------|-----------|-----------|
| Distinct minimizers | 13,309 | 32,166 | 9,729 | 24,537 |
| Superkmers | 251M | 258M | 246M | **245M** |
| Max bucket | 22.1M | 16.4M | 48.3M | 20.5M |
| Median bucket | 4,157 | 9,509 | 13,361 | 20,615 |
| Mean bucket | 234,224 | 96,907 | 320,412 | 127,045 |
| Max/Mean | **94x** | 170x | 151x | 161x |
| Max/Median | 5,312x | **1,729x** | 3,614x | 995x |

Notes:
- **Syncmers** (`iteratorsyncmers2`, l=8): best max/mean in mspxor mode (94x). Sticky mode
  has best superkmer count (140M) but high max bucket (28M) due to Alu-related AATGGAAT.
- **UHS** (`iteratoruhs`): RC-closed ry-alphabet patterns (Martin Frith). l=8 sticky has
  best max/mean (42x) and max/median (103x) of any method. l=9 uses 84 RC-closed patterns.
- **SIMD-mini** (`iteratorsimdmini`, l=9): SIMD-accelerated canonical closed syncmers.
  Competitive across modes, best max bucket in mspxor (16.4M). Requires odd l.
- **KMC2** (`iteratorkmc2`, l=8): disqualification rules backfire on repeat-rich genomes.
- **MSP/lexicographic** (`iteratormsp`, l=8): most distinct minimizers but extreme skew.

Top offenders (sticky):
- Syncmer: AATGGAAT (28.0M), AAAAAAAG (3.8M) — Alu-related
- SIMD-mini: AATGGAATG (18.1M) — Alu
- UHS l=8: CCATTCCA (11.7M), GGAATGGA (10.9M) — Alu RC pair
- UHS l=9: AATGGAATG (20.5M), ATGGAATCA (7.3M) — Alu
- KMC2: CCATTCCA (34.7M), CAGCCTGG (26.4M) — Alu consensus
- MSP: AAAAAAAT (18.2M), AAAAAAAG (14.6M) — homopolymer-adjacent

Top offenders (mspxor):
- Syncmer: GGGATTAC (22.1M), GCCTCCCA (12.3M)
- SIMD-mini: CGAATGGAA (16.4M), CAGCCTCCC (12.3M)
- UHS l=8: GGAATGGA (48.3M), GGGATTAC (22.6M) — Alu hotspot
- UHS l=9: CTGGGATTA (20.5M), CTGAGGCAG (10.8M)

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

## Historical Artifacts

- **`src/iteratorsyncmers2.claudemadness.rs`** — A failed 1107-line attempt by an early Claude Code
  version (2025) to write an AVX2-accelerated syncmers2 iterator. Kept as a historical artifact.
  Features hardcoded test outputs for specific sequence lengths, fake AVX2 batching (loads registers
  then discards results), debug prints everywhere, and a fundamentally wrong MSP algorithm. Do not
  use, do not integrate, do not update.

## Dead Ends

- **SyncKMC2** (syncmer + KMC2 disqualification hybrid): tested with both lexicographic and
  hash-based tiebreaking. Both worse than plain syncmers. In A-rich regions all l-mers are
  KMC2-disqualified, so the filter can't redirect load. Max bucket 40.7M (vs 16.6M syncmers).
  KMC2 disqualification is counterproductive on top of syncmers — it removes valid candidates
  without providing better alternatives.

## Homopolymer Minimizer Demotion

All-A and all-T l-mers are valid closed syncmers but cause hot buckets on repeat-rich genomes.
They are demoted (score 0→1) so they're only selected when no other syncmer exists in the window.
Both forward-strand encodings must be demoted (`scores[0]` for AAA...A, `scores[4^l-1]` for
TTT...T) because scoring happens before canonicalization. Applied to `iteratorsyncmers2`,
`iteratorsyncmersmsp`, and `iteratorsimdmini`. The `iteratormsp` uses `usize::MAX` (original
fix from 2024, commit 3a7b239). The SIMD iterator demotes lazily during rescan only (when the
current minimizer falls off the left edge), not mid-slide — correct and simple but adds ~16%
overhead on 1M random DNA (295→249 MB/s) due to closure capture and macro expansion in the
hot loop.

Run: `cargo +nightly run --release --bin bucket_stats -- <genome.fa> [k] [l] [method[:mode]]`
Methods: `syncmer`, `simdmini`, `cminim`, `uhs`, `kmc2`, `msp`, `multimini[:N]`. Modes: `mspxor`, `classical`, `msp`.

## Profiling

On hybrid CPUs (Intel Alder Lake / Raptor Lake), `perf record` captures two event types:
`cpu_atom/cycles/` (efficiency cores) and `cpu_core/cycles/` (performance cores). The atom
profile appears first in `perf report` and is misleading — focus on the `cpu_core/cycles/`
section (the one with ~100x more samples). Always `rm -f perf.data perf.data.old` between runs.

```bash
# Profile:
perf record -o /tmp/perf_superkmers.data --call-graph dwarf -- target/release/bench_iterators
# Report (flat, ≥1% overhead):
perf report -i /tmp/perf_superkmers.data --stdio --no-children -g none --percent-limit 1
```