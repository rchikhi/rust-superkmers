# Autoresearch: syncmers_simd_l8k40max optimization

SIMD 8-read batch closed syncmer superkmer extractor (AVX2, l=8, k=31). Goal: 1.8 GB/s at 150bp (currently ~780 MB/s). Reference: csyncmer_fast achieves 1.8 GB/s for syncmer detection alone (no superkmers). Must produce identical results to scalar syncmers2-ext:mspxor.

## Setup

- **File to modify:** `src/syncmers_simd_l8k40max.rs` (only file in scope)
- **Build:** `RUSTFLAGS='-C target-cpu=native' cargo +nightly build --release`
- **Benchmark:** `RUSTFLAGS='-C target-cpu=native' cargo +nightly run --release --bin bench_iterators -- "simd-batch-ext" 2>&1 | grep "simd-batch-ext (l=8)" | head -1`
- **Test:** `cargo +nightly test simd_batch 2>&1 | grep "passed"`

### Permissions (auto-approve for uninterrupted loop)

Add these to `~/.claude/settings.json` under `permissions.allow`:
```json
"Bash(RUSTFLAGS=*)",
"Bash(cargo +nightly build *)",
"Bash(cargo +nightly test *)",
"Bash(cargo +nightly run *)",
"Bash(cp src/syncmers_simd_l8k40max.rs /tmp/simd_backup.rs)",
"Bash(cp /tmp/simd_backup.rs src/syncmers_simd_l8k40max.rs)",
"Edit(src/syncmers_simd_l8k40max.rs)",
"Edit(AUTORESEARCH.md)",
"Edit(results.tsv)",
"Write(results.tsv)"
```
**Important:** Permissions must be EXACT command prefixes, not broad wildcards.
`Bash(for i in *)` is too broad — use the full benchmark command instead:
`Bash(for i in 1 2 3; do RUSTFLAGS=... done)`. Any unexpected permission prompt
breaks the loop. Stick to commands already in `.claude/settings.local.json`.

- **Metric:** ns/read at 150bp (lower is better). Target: ≤ 83 ns (1.8 GB/s). Baseline: ~192 ns.
- **Log:** `results.tsv` (untracked)
- **Constraint:** Results must match scalar syncmers2-ext:mspxor (run tests after each change).

## Current architecture (what you're optimizing)

```
Pack (17 ns): AVX2 32-base chunks → 2-bit LSB-first bytes
Kernel (140 ns): 32-op hot loop × 120 iters @ 9.4 cycles/iter, 3.4 IPC
  - load_base: gather/16 + srlv + and (3 ops amortized)
  - update_rolling: fwd slli+or+and (3) + RC sub+slli+srli+or+and (5) = 8 ops
  - issue_gather: srli+and+gather+min = 3 ops + 1 L1 gather (12 cycles, hidden by pipeline)
  - build_elem: srlv+and+cmpeq+xor+slli+and+or+andnot+or = 9 ops
  - window_push: store+min + amortized suffix recomp = 4.9 ops
  - query: load+min = 2 ops
  - store min_elem: 1 op
  - counter: add = 1 op
Materialization (19 ns): canonical value extracted from window element, ~6 superkmers
Post-scan (5 ns): SIMD cmpeq scan for change detection
```

8KB syncmer bit table (L1-resident). Software-pipelined gather (issue at iter N, use at N+1). Two-stack block-decomposition sliding window on stack array (ring[32]). Deferred change detection.

## Results TSV format

```
experiment	ns_per_read	mb_per_sec	status	description
```

Status: `keep` (improved), `discard` (equal or worse), `crash` (broke tests).

## Experiment loop

LOOP FOREVER:

1. Read current `src/syncmers_simd_l8k40max.rs` and this file for context
2. Pick an idea from the list below (or invent one)
3. **Snapshot:** `cp src/syncmers_simd_l8k40max.rs /tmp/simd_backup.rs`
4. Edit `src/syncmers_simd_l8k40max.rs`
5. **Build:** `RUSTFLAGS='-C target-cpu=native' cargo +nightly build --release`
6. If compile fails, fix and retry (max 3 attempts, then discard)
7. **Test:** `cargo +nightly test simd_batch` — must pass all 6 tests
8. **Benchmark 3 times**, take the median:
   ```
   for i in 1 2 3; do RUSTFLAGS='-C target-cpu=native' cargo +nightly run --release --bin bench_iterators -- "simd-batch-ext" 2>&1 | grep "simd-batch-ext (l=8)" | head -1; done
   ```
9. Log to `results.tsv`
10. If faster (median < previous best): **commit** (`git add src/syncmers_simd_l8k40max.rs && git commit -m "autoresearch: <description>"`) — do NOT push yet
11. If equal or slower: **revert:** `git checkout -- src/syncmers_simd_l8k40max.rs` (reverts to last commit)
12. Cross off the idea below, note the result
13. Pick next idea

NEVER STOP. NEVER ASK. Run until interrupted.

## Key constraint: the compiler is smart

Previous attempt to "remove ops" (RC mask, compact syncmer→priority, pos mask removal) made things SLOWER despite fewer instructions. The compiler's register allocator and instruction scheduler are near-optimal. Changes that look like wins on paper may regress in practice. ALWAYS measure.

## Ideas to try

Roughly ordered by expected impact. Check off as tested.

### Hot loop restructuring
- [ ] Unroll by 2: process 2 l-mers per loop iteration (share gather, reduce loop overhead)
- [ ] Unroll by 4: process 4 l-mers (load full byte from packed data with fixed shifts)
- [ ] Move suffix recomputation to a separate loop after main loop (separate hot path from cold path)
- [ ] Try ring buffer size = next power of 2 (32 already) with bitwise masking instead of comparison

### Reducing ops
- [ ] Combine `srli(fwd,3)` and `and(fwd,7)` into single `srli` + reuse fwd for and (check if compiler already does this)
- [ ] Use `vpxor` instead of `vpsubd` for complement (same result, more port options)
- [ ] Eliminate `vpand(pos, pos_mask)` in build_elem (pos < 0x7FFF for 150bp reads)
- [ ] Pre-shift the XOR constant: compute `xor_const_shifted = xor_const << 15` once, then `slli(canon,15) ^ xor_const_shifted` instead of `xor(canon,xor) + slli(tb,15)`
- [ ] Merge the `and(shifted, one_vec) + cmpeq(_, one_vec)` into `slli(shifted,31) + srai(result,31)`

### Memory layout
- [ ] Align ring buffer to cache line boundary (64 bytes)
- [ ] Use `_mm256_store_si256` (aligned) instead of indexing for ring writes
- [ ] Pre-allocate min_pos_buf on stack instead of heap (`[__m256i; 128]`)
- [ ] Reduce min_pos_buf to store only position (not full element) — recover canon during phase 5

### Packing optimization
- [ ] Skip packing entirely: load ASCII bytes directly in the kernel, compute 2-bit inline
- [ ] Process packing and first few kernel iterations in parallel (overlap pack with warmup)

### Gather optimization
- [ ] Try 2-iteration pipeline delay (issue gather at N, use at N+2) to fully hide L1 latency
- [ ] Prefetch next gather address with `_mm_prefetch`
- [ ] Try byte table (64KB) instead of bit table (8KB) — eliminates srlv+and bit extraction (saves 2 ops) but 64KB may spill L1

### Sliding window
- [ ] Try smaller window via reduced k (k=21 instead of 31 → w=14, faster suffix recomp)
- [ ] Lazy suffix recomputation: only recompute when the window minimum might have changed
- [ ] Try Lemire's "ascending minima" deque instead of two-stack (different amortization pattern)

### Phase 5 (post-loop)
- [ ] Process phase 5 change detection with SIMD comparison on 4 elements at a time
- [ ] Fuse phase 5 into materialization (single pass instead of two)

### Precompute fwd values + 256KB combined table with prefetch
The BIG idea. Currently: 8KB bit table (L1) + inline RC + canonical + XOR + pack = 32 ops.
Alternative: precompute all rolling_fwd values into a buffer BEFORE the kernel (143 × 3 ops = cheap).
Then use a 256KB combined table: `table[fwd] = (priority << 31) | (tiebreaker << 15) | (is_rc_bit)`.
One gather gives EVERYTHING — score + is_rc. Eliminates 19 ops → 13 ops/iter.
The precomputed fwd_buf lets us PREFETCH future table entries:
```
prefetch(table + fwd_buf[q + 8])  // 8 iters ahead → data in L1 when gather executes
```
With prefetch: gather hits L1 (~4 cycles) instead of L2 (~14 cycles).
13 ops at 3.4 IPC = 3.8 cycles + 4 cycle gather = ~8 cycles/iter.
120 × 8 = 960 cycles / 8 = 120 ns kernel → ~83 ns total → 1.8 GB/s.
Challenge: need 8 separate _mm_prefetch per iter (1 per lane). Use storeu + 8 scalar prefetches.
- [ ] Phase 0: precompute fwd_buf[0..143] as __m256i array (rolling fwd, no RC)
- [ ] Combined table: generate at startup (256KB, lazy_static)
- [ ] Kernel: gather from combined table + prefetch 8 iters ahead
- [ ] Measure prefetch benefit vs overhead

### Eliminate packing: ASCII-direct loading
- [ ] Remove pack_seq_2bit entirely (save 17 ns)
- [ ] In kernel, load ASCII bytes via gather: `gather(seq_ptrs, offset, 1)` → 4 bytes/lane
- [ ] Extract base inline: `((b >> 1) ^ (b >> 2)) & 3` in SIMD
- [ ] Load 4 bases per gather (unroll by 4) to amortize gather cost

### Interleave two groups of 8 reads (16 total)
- [ ] Process Groups A and B in alternating half-iterations
- [ ] Group A's gather executes while Group B's compute runs (and vice versa)
- [ ] Doubles compute between gather issue and use → hides L2 latency
- [ ] Doubles register usage (2 × window state) — exactly 16 YMM, very tight
- [ ] Net: same ops/l-mer but better resource utilization

### Move is_rc out of hot loop
- [ ] Remove cmpeq + xor + movemask + store from hot loop (5 ops saved)
- [ ] During Phase 5, recompute fwd from packed data at change points (~6 per read)
- [ ] `is_rc = (fwd_from_packed != canon_val)` — only ~48 byte accesses total

### Radical ideas
- [ ] Inline the entire kernel (no function call) — `#[inline(always)]` on simd_kernel
- [ ] Use `core::arch::asm!` for the hot loop (hand-written assembly)
- [ ] Process 16 reads with AVX-512 (if available on the CPU)
- [ ] Split kernel into syncmer-detection pass + minimizer pass (two simpler loops)
- [ ] Fuse packing with fwd_buf precomputation (pack + rolling fwd in one pass)
- [ ] Use SIMD comparison in the hot loop instead of deferred (the old approach was 8 ops but some were free on idle ports)
