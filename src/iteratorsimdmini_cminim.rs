//! Thin wrapper over simd-minimizers' `canonical_minimizers` with built-in
//! super-kmer extraction. Uses ntHash-based random minimizers (not syncmers).
//!
//! The library handles the sliding window internally via SIMD, so no manual
//! deque or scoring is needed here.
//!
//! Requires odd k (e.g. k=31) for canonical mode (simd-minimizers constraint:
//! l = w + k_min - 1 must be odd, and their l = our k).

use crate::Superkmer;
use simd_minimizers::packed_seq::AsciiSeq;

/// Encode ASCII base to 2-bit (A=0, C=1, G=2, T=3).
#[inline(always)]
fn encode_base(b: u8) -> u64 {
    (((b >> 1) ^ (b >> 2)) & 3) as u64
}

/// Collect superkmers from a single N-free fragment.
fn superkmers_from_fragment(
    ascii_slice: &[u8],
    k: usize,
    l: usize,
    offset: usize,
    canonical: bool,
    results: &mut Vec<Superkmer>,
    min_pos_buf: &mut Vec<u32>,
    sk_pos_buf: &mut Vec<u32>,
) {
    let seq_len = ascii_slice.len();
    if seq_len < k {
        return;
    }

    assert!(
        k % 2 == 1,
        "canonical_minimizers requires odd k (their l = our k), got k={}",
        k
    );

    let w = k - l + 1; // window size in simd-minimizers notation

    min_pos_buf.clear();
    sk_pos_buf.clear();

    // Collect canonical values eagerly to release the borrow on min_pos_buf.
    let vals: Vec<u64> = simd_minimizers::canonical_minimizers(l, w)
        .super_kmers(sk_pos_buf)
        .run(AsciiSeq(ascii_slice), min_pos_buf)
        .values_u64()
        .collect();

    let n = min_pos_buf.len();
    if n == 0 {
        return;
    }

    results.reserve(n);
    for i in 0..n {
        let sk_start = sk_pos_buf[i] as usize;
        let last_kmer = if i + 1 < n {
            sk_pos_buf[i + 1] as usize - 1
        } else {
            seq_len - k
        };
        let sk_end = last_kmer + k;
        let size = sk_end - sk_start;
        let mpos = min_pos_buf[i] as usize - sk_start;
        let mint = vals[i] as u32;

        let mint_is_rc = if canonical {
            // Compare canonical value with forward l-mer to determine RC
            let pos = min_pos_buf[i] as usize;
            let mut fwd = 0u64;
            for j in 0..l {
                fwd = (fwd << 2) | encode_base(ascii_slice[pos + j]);
            }
            mint as u64 != fwd
        } else {
            false
        };

        results.push(Superkmer {
            start: sk_start + offset,
            mint,
            size: size as u16,
            mpos: mpos as u16,
            mint_is_rc,
        });
    }
}

pub struct SuperkmersIterator {
    superkmers: Vec<Superkmer>,
    pos: usize,
}

impl SuperkmersIterator {
    pub fn new(seq: &[u8], k: usize, l: usize) -> Self {
        let mut superkmers = Vec::new();
        let mut min_pos_buf = Vec::new();
        let mut sk_pos_buf = Vec::new();
        superkmers_from_fragment(seq, k, l, 0, true, &mut superkmers, &mut min_pos_buf, &mut sk_pos_buf);
        SuperkmersIterator { superkmers, pos: 0 }
    }

    pub fn new_with_n(seq: &[u8], k: usize, l: usize) -> Self {
        let fragments = crate::utils::split_on_n(seq, k);
        let mut superkmers = Vec::new();
        let mut min_pos_buf = Vec::new();
        let mut sk_pos_buf = Vec::new();
        for (offset, fragment) in &fragments {
            superkmers_from_fragment(fragment, k, l, *offset, true, &mut superkmers, &mut min_pos_buf, &mut sk_pos_buf);
        }
        SuperkmersIterator { superkmers, pos: 0 }
    }
}

impl Iterator for SuperkmersIterator {
    type Item = Superkmer;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.superkmers.len() {
            return None;
        }
        let idx = self.pos;
        self.pos += 1;
        Some(Superkmer {
            start: self.superkmers[idx].start,
            mint: self.superkmers[idx].mint,
            size: self.superkmers[idx].size,
            mpos: self.superkmers[idx].mpos,
            mint_is_rc: self.superkmers[idx].mint_is_rc,
        })
    }
}
