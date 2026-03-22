//! SIMD 8-read batch closed syncmer extractor (AVX2, l=8, k≤40, s=2, mspxor).
//!
//! Processes 8 short reads simultaneously using AVX2 SIMD. Syncmer status
//! is determined by an 8KB bit table lookup via SIMD gather (L1-resident,
//! latency hidden by software pipelining). A single two-stack sliding window
//! performs minimizer selection with mspxor tiebreaking.
//!
//! Constraints: l=8, k≤40, s=2, reads ≤32KB. Results match syncmers2-ext:mspxor.
//! ~800 MB/s at 150bp (2× scalar).
//!
//! Requires AVX2.

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

use crate::Superkmer;
use crate::minimizer_core::{canonical_table, base_from_ascii};
use lazy_static::lazy_static;

// ---------------------------------------------------------------------------
// 8KB syncmer bit table
// ---------------------------------------------------------------------------

/// 8KB syncmer bit table: bit N set if l-mer N is a closed syncmer.
fn generate_syncmer_bit_table(l: usize) -> Vec<u8> {
    let scores = crate::iteratorsyncmers2::generate_syncmer_scores_with_s(l, 2);
    let num_lmers = 1 << (2 * l);
    let mut bits = vec![0u8; num_lmers / 8];
    for lmer in 0..num_lmers {
        if scores[lmer] == 0 {
            bits[lmer / 8] |= 1 << (lmer % 8);
        }
    }
    bits
}

lazy_static! {
    static ref SYNCMER_BITS_8: Vec<u8> = generate_syncmer_bit_table(8);
}

// ---------------------------------------------------------------------------
// Packing helpers
// ---------------------------------------------------------------------------

/// Pack ASCII base to 2-bit: A=0, C=1, G=2, T=3 (standard encoding, matches get_base).
#[inline(always)]
fn pack_base(b: u8) -> u8 {
    (((b >> 1) ^ (b >> 2)) & 3) as u8
}

/// Pack a sequence into 2-bit (4 bases per byte, LSB-first within byte).
/// Uses AVX2 when available for 32-byte chunks.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn pack_seq_2bit(seq: &[u8], out: &mut [u8]) {
    let len = seq.len();
    let mut i = 0;

    // AVX2 path: process 32 ASCII bases → 8 packed bytes
    // Encoding: ((b >> 1) ^ (b >> 2)) & 3, packed 4 bases per byte LSB-first
    while i + 32 <= len {
        let input = _mm256_loadu_si256(seq.as_ptr().add(i) as *const __m256i);
        let shifted1 = _mm256_srli_epi16(input, 1);
        let shifted2 = _mm256_srli_epi16(input, 2);
        let xored = _mm256_xor_si256(shifted1, shifted2);
        let bases = _mm256_and_si256(xored, _mm256_set1_epi8(3));

        // Pack 4 bases per byte: byte[n] = b[4n] | (b[4n+1]<<2) | (b[4n+2]<<4) | (b[4n+3]<<6)
        // Use pmaddubsw to combine pairs, then pmaddwd for final packing
        // Step 1: multiply odd bytes by 4 (shift left 2), add to even bytes
        let mul_pairs = _mm256_set1_epi16(0x0401); // even=1, odd=4
        let paired = _mm256_maddubs_epi16(bases, mul_pairs); // 2 bases per u16

        // Step 2: combine u16 pairs: low_u16 | (high_u16 << 4)
        // Use pmaddwd: low_u16 * 1 + high_u16 * 16 → u32
        let mul_quads = _mm256_set1_epi32(0x00100001); // low=1, high=16
        let quads = _mm256_madd_epi16(paired, mul_quads); // 4 bases per u32

        // Step 3: pack u32 → u8 via shuffle + permute
        // Each u32 has its packed byte in the lowest byte
        let shuffle_mask = _mm256_set_epi8(
            -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, 12,8,4,0,
            -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, 12,8,4,0);
        let shuffled = _mm256_shuffle_epi8(quads, shuffle_mask);
        // Low 4 bytes of each 128-bit lane contain the packed data
        // Lane 0: bytes 0-3, Lane 1: bytes 16-19
        let lo = _mm256_extracti128_si256(shuffled, 0);
        let hi = _mm256_extracti128_si256(shuffled, 1);
        // Combine: lo has 4 bytes at position 0-3, hi has 4 bytes at position 0-3
        let combined = _mm_unpacklo_epi32(lo, hi); // 8 bytes at positions 0-7
        _mm_storel_epi64(out.as_mut_ptr().add(i / 4) as *mut __m128i, combined);

        i += 32;
    }

    // Scalar remainder
    while i + 4 <= len {
        out[i / 4] = pack_base(seq[i])
            | (pack_base(seq[i + 1]) << 2)
            | (pack_base(seq[i + 2]) << 4)
            | (pack_base(seq[i + 3]) << 6);
        i += 4;
    }
    if i < len {
        let mut byte = 0u8;
        for j in i..len {
            byte |= pack_base(seq[j]) << (((j % 4) * 2) as u8);
        }
        out[i / 4] = byte;
    }
}

// ---------------------------------------------------------------------------
// SIMD helpers
// ---------------------------------------------------------------------------

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline(always)]
unsafe fn simd_min_u32(a: __m256i, b: __m256i) -> __m256i {
    _mm256_min_epu32(a, b)
}

// ---------------------------------------------------------------------------
// Superkmer boundary tracking
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct MinChange {
    kmer_pos: u32,
    min_lmer_pos: u32,
    canon_val: u32,
    is_rc: bool,
}

// ---------------------------------------------------------------------------
// SIMD kernel: single window + table gather
// ---------------------------------------------------------------------------

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn simd_kernel(
    packed_seqs: &[u8],
    packed_offsets: &[usize; 8],
    seq_lengths: &[usize; 8],
    syncmer_bits: &[u8],       // 8KB bit table (L1 resident)
    k: usize,
    l: usize,
    min_pos_buf: &mut [__m256i],
    is_rc_buf: &mut [u8],     // 1 byte per l-mer position (8 lane bits)
    changes: &mut [Vec<MinChange>; 8],
) {
    let w = k - l + 1;
    let lmer_mask = ((1u32 << (2 * l)) - 1) as i32;
    let xor_const = 0xACE5i32;

    // Per-lane k-mer counts
    let mut per_lane_kmers = [0usize; 8];
    let mut max_num_kmers = 0usize;
    for i in 0..8 {
        if seq_lengths[i] >= k {
            per_lane_kmers[i] = seq_lengths[i] - k + 1;
            if per_lane_kmers[i] > max_num_kmers {
                max_num_kmers = per_lane_kmers[i];
            }
        }
    }
    if max_num_kmers == 0 { return; }

    // Total base positions = max_seq_len = max_num_kmers + k - 1
    let max_bases = max_num_kmers + k - 1;
    // Total l-mer positions
    let max_lmers = max_bases - l + 1; // = max_num_kmers + k - l = max_num_kmers + w - 1

    for v in changes.iter_mut() { v.clear(); }

    // Gather offsets (byte offset of each lane's packed data)
    let mut gather_arr = [0i32; 8];
    for i in 0..8 { gather_arr[i] = packed_offsets[i] as i32; }
    let gather_base_vec = _mm256_loadu_si256(gather_arr.as_ptr() as *const __m256i);

    // Lane limits
    let mut lane_limit_arr = [0i32; 8];
    for i in 0..8 { lane_limit_arr[i] = per_lane_kmers[i] as i32; }

    // Rolling l-mer state: fwd + RC for canonical computation
    let mut rolling_fwd = _mm256_setzero_si256();
    let mut rolling_rc = _mm256_setzero_si256();
    let lmer_mask_vec = _mm256_set1_epi32(lmer_mask);
    let three_vec = _mm256_set1_epi32(3);
    debug_assert!(l == 8, "SIMD kernel currently only supports l=8");
    const RC_SHIFT: i32 = 14;

    // Pre-hoisted constants (avoid per-iteration broadcasts)
    let one_vec = _mm256_set1_epi32(1);
    let seven_vec = _mm256_set1_epi32(7);
    let pos_mask = _mm256_set1_epi32(0x7FFF);
    let sign_bit = _mm256_set1_epi32(0x80000000u32 as i32);
    let all_ones = _mm256_set1_epi32(-1i32);
    let xor_const_vec = _mm256_set1_epi32(xor_const);

    // Software-pipelined gather: buffer prev iteration's results
    let mut prev_gather = _mm256_setzero_si256();
    let mut prev_bit_idx = _mm256_setzero_si256();
    let mut prev_canon = _mm256_setzero_si256();
    let mut prev_pos = _mm256_setzero_si256();

    // Minimizer sliding window (two-stack) — STACK array, no heap indirection
    assert!(w <= 32, "window size must be <= 32");
    let mut ring = [all_ones; 32]; // fixed-size stack array (max w=32)
    let mut prefix_min = all_ones;
    let mut ring_idx = 0usize;

    // min_pos_buf passed in (pre-allocated)

    let mut pos_vec = _mm256_setzero_si256();
    let mut cur_data = _mm256_setzero_si256();
    let mut seq_pos: usize = 0;

    macro_rules! load_base {
        () => {{
            if seq_pos % 16 == 0 {
                let bv_idx = _mm256_add_epi32(gather_base_vec,
                    _mm256_set1_epi32((seq_pos / 4) as i32));
                cur_data = _mm256_i32gather_epi32::<1>(
                    packed_seqs.as_ptr() as *const i32, bv_idx);
            }
            let b = _mm256_and_si256(
                _mm256_srlv_epi32(cur_data, _mm256_set1_epi32(((seq_pos % 16) * 2) as i32)),
                three_vec);
            seq_pos += 1;
            b
        }};
    }

    // Rolling fwd + RC (8 ops)
    macro_rules! update_rolling {
        ($base:expr) => {{
            rolling_fwd = _mm256_and_si256(
                _mm256_or_si256(_mm256_slli_epi32(rolling_fwd, 2), $base),
                lmer_mask_vec);
            let comp = _mm256_sub_epi32(three_vec, $base);
            rolling_rc = _mm256_and_si256(
                _mm256_or_si256(
                    _mm256_srli_epi32(rolling_rc, 2),
                    _mm256_slli_epi32(comp, RC_SHIFT)),
                lmer_mask_vec);
        }};
    }

    // Issue L1 gather + canonical + is_rc (3 ops + 1 gather + 3 ops off critical path)
    let mut lmer_idx: usize = 0; // l-mer position counter for is_rc_buf
    macro_rules! issue_gather {
        () => {{
            let byte_idx = _mm256_srli_epi32(rolling_fwd, 3);
            let bit_idx = _mm256_and_si256(rolling_fwd, seven_vec);
            let gathered = _mm256_i32gather_epi32::<1>(
                syncmer_bits.as_ptr() as *const i32, byte_idx);
            let canon = simd_min_u32(rolling_fwd, rolling_rc);
            // is_rc: off critical path (canon and fwd ready early, result only needed in materialization)
            let not_rc = _mm256_cmpeq_epi32(canon, rolling_fwd);
            let is_rc_mask = _mm256_xor_si256(not_rc, all_ones);
            *is_rc_buf.get_unchecked_mut(lmer_idx) =
                _mm256_movemask_ps(_mm256_castsi256_ps(is_rc_mask)) as u8;
            lmer_idx += 1;
            (gathered, bit_idx, canon)
        }};
    }

    // Build score from pipelined gather result (9 ops)
    macro_rules! build_elem {
        ($gather:expr, $bit_idx:expr, $canon:expr, $pos:expr) => {{
            let shifted = _mm256_srlv_epi32($gather, $bit_idx);
            let is_syncmer = _mm256_cmpeq_epi32(
                _mm256_and_si256(shifted, one_vec), one_vec);
            let tb = _mm256_xor_si256($canon, xor_const_vec);
            let score = _mm256_or_si256(
                _mm256_slli_epi32(tb, 15),
                _mm256_and_si256($pos, pos_mask));
            _mm256_or_si256(score, _mm256_andnot_si256(is_syncmer, sign_bit))
        }};
    }

    macro_rules! window_push {
        ($elem:expr) => {{
            ring[ring_idx] = $elem;
            prefix_min = simd_min_u32(prefix_min, $elem);
            ring_idx += 1;
            if ring_idx == w {
                ring_idx = 0;
                let mut smin = ring[w - 1];
                for j in (0..w - 1).rev() {
                    smin = simd_min_u32(smin, ring[j]);
                    ring[j] = smin;
                }
                prefix_min = all_ones;
            }
        }};
    }

    // ===== Phase 1: Warmup =====
    for _ in 0..l - 1 {
        let base = load_base!();
        update_rolling!(base);
    }

    // ===== Phase 2: Prime pipeline =====
    {
        let base = load_base!();
        update_rolling!(base);
        let (g, bi, c) = issue_gather!();
        prev_gather = g; prev_bit_idx = bi; prev_canon = c;
        prev_pos = pos_vec;
        pos_vec = _mm256_add_epi32(pos_vec, one_vec);
    }

    // ===== Phase 3: Fill window =====
    for _ in 1..w {
        let base = load_base!();
        update_rolling!(base);
        let (g, bi, c) = issue_gather!();
        let elem = build_elem!(prev_gather, prev_bit_idx, prev_canon, prev_pos);
        window_push!(elem);
        prev_gather = g; prev_bit_idx = bi; prev_canon = c;
        prev_pos = pos_vec;
        pos_vec = _mm256_add_epi32(pos_vec, one_vec);
    }

    // ===== Phase 4a: Hot loop — no branches, no flush check =====
    let num_queries = max_num_kmers;
    let main_iters = max_lmers - w; // iterations with new bases available

    for q in 0..main_iters {
        let base = load_base!();
        update_rolling!(base);
        let (g, bi, c) = issue_gather!();

        let elem = build_elem!(prev_gather, prev_bit_idx, prev_canon, prev_pos);
        window_push!(elem);

        let suf_min = ring[ring_idx];
        let min_elem = simd_min_u32(prefix_min, suf_min);
        // Store FULL min_elem (not just position) — contains tiebreaker for materialization
        *min_pos_buf.get_unchecked_mut(q) = min_elem;

        prev_gather = g; prev_bit_idx = bi; prev_canon = c;
        prev_pos = pos_vec;
        pos_vec = _mm256_add_epi32(pos_vec, one_vec);
    }

    // ===== Phase 4b: Flush =====
    for q in main_iters..num_queries {
        let elem = build_elem!(prev_gather, prev_bit_idx, prev_canon, prev_pos);
        window_push!(elem);

        let suf_min = ring[ring_idx];
        let min_elem = simd_min_u32(prefix_min, suf_min);
        *min_pos_buf.get_unchecked_mut(q) = min_elem;

        prev_gather = _mm256_setzero_si256();
        prev_bit_idx = _mm256_setzero_si256();
        prev_canon = _mm256_setzero_si256();
        prev_pos = pos_vec;
        pos_vec = _mm256_add_epi32(pos_vec, one_vec);
    }

    // ===== Phase 5: Post-loop change detection =====
    // min_pos_buf stores full window min elements: (priority|tiebreaker|position).
    // Compare on position (low 15 bits) to detect changes.
    // Extract canonical value from tiebreaker: canon = ((elem >> 15) & 0xFFFF) ^ 0xACE5.
    {
        let xor_recover = 0xACE5u32;

        // Helper: extract position, canon_val, and is_rc from element + is_rc_buf
        let extract = |elem_u32: u32, rc_byte: u8, lane: usize| -> (u32, u32, bool) {
            let pos = elem_u32 & 0x7FFF;
            let canon = ((elem_u32 >> 15) & 0xFFFF) ^ xor_recover;
            let is_rc = (rc_byte >> lane) & 1 != 0;
            (pos, canon, is_rc)
        };

        // First k-mer: always emit
        let first_elem = *min_pos_buf.get_unchecked(0);
        let mut prev_pos_simd = _mm256_and_si256(first_elem, pos_mask);
        {
            let mut arr = [0u32; 8];
            _mm256_storeu_si256(arr.as_mut_ptr() as *mut __m256i, first_elem);
            for lane in 0..8 {
                if per_lane_kmers[lane] > 0 {
                    // The minimizer's l-mer position tells us which is_rc_buf entry to look up
                    let (pos, canon, _) = extract(arr[lane], 0, lane);
                    let rc_byte = *is_rc_buf.get_unchecked(pos as usize);
                    let is_rc = (rc_byte >> lane) & 1 != 0;
                    changes[lane].push(MinChange {
                        kmer_pos: 0, min_lmer_pos: pos, canon_val: canon, is_rc,
                    });
                }
            }
        }

        // SIMD scan: compare positions, extract canon only on changes (~6 per read)
        for q in 1..num_queries {
            let cur_elem = *min_pos_buf.get_unchecked(q);
            let cur_pos_simd = _mm256_and_si256(cur_elem, pos_mask);
            let diff = _mm256_xor_si256(
                _mm256_cmpeq_epi32(cur_pos_simd, prev_pos_simd), all_ones);
            let mask = _mm256_movemask_ps(_mm256_castsi256_ps(diff)) as u32;
            if mask != 0 {
                let mut arr = [0u32; 8];
                _mm256_storeu_si256(arr.as_mut_ptr() as *mut __m256i, cur_elem);
                let mut km = mask;
                while km != 0 {
                    let lane = km.trailing_zeros() as usize;
                    if (q as usize) < per_lane_kmers[lane] {
                        let pos = arr[lane] & 0x7FFF;
                        let canon = ((arr[lane] >> 15) & 0xFFFF) ^ xor_recover;
                        let rc_byte = *is_rc_buf.get_unchecked(pos as usize);
                        let is_rc = (rc_byte >> lane) & 1 != 0;
                        changes[lane].push(MinChange {
                            kmer_pos: q as u32, min_lmer_pos: pos, canon_val: canon, is_rc,
                        });
                    }
                    km &= km - 1;
                }
            }
            prev_pos_simd = cur_pos_simd;
        }
    }
}

// ---------------------------------------------------------------------------
// Superkmer materialization
// ---------------------------------------------------------------------------

fn materialize_from_changes(
    seq_len: usize, changes: &[MinChange], k: usize, _l: usize,
    out: &mut Vec<Superkmer>,
) {
    if changes.is_empty() { return; }

    for p in 0..changes.len() {
        let start_pos = changes[p].kmer_pos as usize;
        let min_lmer_pos = changes[p].min_lmer_pos as usize;
        let size = if p < changes.len() - 1 {
            let next_pos = changes[p + 1].kmer_pos as usize;
            next_pos + k - 1 - start_pos
        } else {
            seq_len - start_pos
        };

        let mint = changes[p].canon_val;
        let mint_is_rc = changes[p].is_rc;

        out.push(Superkmer {
            start: start_pos, mint,
            size: size as u16,
            mpos: (min_lmer_pos - start_pos) as u16,
            mint_is_rc,
        });
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

pub struct SimdBatchExtractor {
    k: usize,
    l: usize,
    packed_buf: Vec<u8>,
    packed_stride: usize, // bytes between lanes in packed_buf
    min_pos_buf: Vec<__m256i>,
    is_rc_buf: Vec<u8>,
    changes: [Vec<MinChange>; 8],
    superkmers: [Vec<Superkmer>; 8],
}

impl SimdBatchExtractor {
    pub fn new(k: usize, l: usize) -> Self {
        // Pre-allocate for 8 × 300bp reads
        let max_read = 300;
        let ps_uniform = ((max_read + 3) / 4 + 32 + 31) & !31;
        let max_kmers = max_read - k + 1;
        let max_lmers = max_read - l + 1;
        SimdBatchExtractor {
            k, l,
            packed_buf: vec![0u8; 8 * ps_uniform + 32],
            packed_stride: ps_uniform,
            min_pos_buf: vec![unsafe { _mm256_setzero_si256() }; max_kmers],
            is_rc_buf: vec![0u8; max_lmers],
            changes: Default::default(),
            superkmers: Default::default(),
        }
    }

    /// Pack only (no kernel). For benchmarking packing overhead.
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    pub unsafe fn bench_pack_only(&mut self, seqs: &[&[u8]; 8]) {
        self.prepare_packed(seqs);
    }

    pub unsafe fn bench_kernel_only(&mut self, seqs: &[&[u8]; 8]) -> usize {
        self.prepare_and_run_kernel(seqs);
        self.changes.iter().map(|c| c.len()).sum()
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    pub unsafe fn process_batch(&mut self, seqs: &[&[u8]; 8]) -> &[Vec<Superkmer>; 8] {
        let seq_lengths = self.prepare_and_run_kernel(seqs);

        for i in 0..8 {
            self.superkmers[i].clear();
            if seq_lengths[i] >= self.k {
                materialize_from_changes(
                    seq_lengths[i], &self.changes[i],
                    self.k, self.l, &mut self.superkmers[i]);
            }
        }
        &self.superkmers
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    unsafe fn prepare_and_run_kernel(&mut self, seqs: &[&[u8]; 8]) -> [usize; 8] {
        self.prepare_packed(seqs);
        let (packed_offsets, seq_lengths, _) = Self::batch_params_static(seqs);
        let max_len = seq_lengths.iter().copied().max().unwrap_or(0);
        let max_kmers = if max_len >= self.k { max_len - self.k + 1 } else { 0 };
        if self.min_pos_buf.len() < max_kmers {
            self.min_pos_buf.resize(max_kmers, _mm256_setzero_si256());
        }
        let max_lmers = if max_len >= self.l { max_len - self.l + 1 } else { 0 };
        if self.is_rc_buf.len() < max_lmers {
            self.is_rc_buf.resize(max_lmers, 0);
        }
        simd_kernel(&self.packed_buf, &packed_offsets, &seq_lengths,
            &SYNCMER_BITS_8, self.k, self.l, &mut self.min_pos_buf,
            &mut self.is_rc_buf, &mut self.changes);
        seq_lengths
    }

    /// Access the 2-bit packed representation of read `lane` from the last batch.
    /// Format: 4 bases per byte, LSB-first (base 0 at bits 0-1, base 1 at bits 2-3, etc.).
    /// Encoding: A=0, C=1, G=2, T=3 (same as standard, just LSB-first byte packing).
    pub fn packed_storage(&self, lane: usize) -> &[u8] {
        let off = lane * self.packed_stride;
        &self.packed_buf[off..off + self.packed_stride]
    }

    fn batch_params_static(seqs: &[&[u8]; 8]) -> ([usize; 8], [usize; 8], usize) {
        let mut max_len = 0;
        let mut seq_lengths = [0usize; 8];
        for i in 0..8 {
            seq_lengths[i] = seqs[i].len();
            if seqs[i].len() > max_len { max_len = seqs[i].len(); }
        }
        let ps_uniform = ((max_len + 3) / 4 + 32 + 31) & !31;
        let mut packed_offsets = [0usize; 8];
        for i in 0..8 { packed_offsets[i] = i * ps_uniform; }
        (packed_offsets, seq_lengths, ps_uniform)
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    unsafe fn prepare_packed(&mut self, seqs: &[&[u8]; 8]) {
        let mut max_len = 0;
        for s in seqs { if s.len() > max_len { max_len = s.len(); } }
        let ps_uniform = ((max_len + 3) / 4 + 32 + 31) & !31;
        self.packed_stride = ps_uniform;
        let total = 8 * ps_uniform + 32;
        if self.packed_buf.len() < total {
            self.packed_buf.resize(total, 0);
        }
        // No zeroing needed: pack_seq_2bit writes complete bytes (=, not |=),
        // and stale padding bytes beyond the read are masked out by the kernel
        // (invalid lanes have kmer_pos >= per_lane_kmers → changes ignored).
        for i in 0..8 {
            let off = i * ps_uniform;
            if seqs[i].len() > 0 {
                pack_seq_2bit(seqs[i], &mut self.packed_buf[off..]);
            }
        }
    }
}
