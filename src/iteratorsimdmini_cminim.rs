//! Thin wrapper over simd-minimizers' `canonical_minimizers` with built-in
//! super-kmer extraction. Uses ntHash-based random minimizers (not syncmers).
//!
//! The library handles the sliding window internally via SIMD, so no manual
//! deque or scoring is needed here. Mint (canonical l-mer value) is computed
//! inline at superkmer boundaries only.
//!
//! Requires odd k (e.g. k=31) for canonical mode (simd-minimizers constraint:
//! l = w + k_min - 1 must be odd, and their l = our k).

use crate::Superkmer;
use simd_minimizers::packed_seq::AsciiSeq;

/// Encode ASCII base to 2-bit (A=0, C=1, G=2, T=3).
#[inline(always)]
fn encode_base(b: u8) -> usize {
    (((b >> 1) ^ (b >> 2)) & 3) as usize
}

/// Compute canonical l-mer value and RC flag from ASCII at a given position.
#[inline(always)]
fn canonical_lmer(ascii: &[u8], pos: usize, l: usize) -> (u32, bool) {
    let mut fwd = 0usize;
    for j in 0..l {
        fwd = (fwd << 2) | encode_base(ascii[pos + j]);
    }
    let mut rc = 0usize;
    let mut v = fwd;
    for _ in 0..l {
        rc = (rc << 2) | (3 - (v & 3));
        v >>= 2;
    }
    if rc < fwd { (rc as u32, true) } else { (fwd as u32, false) }
}

/// Convert SIMD output buffers into Superkmer structs.
#[inline]
fn emit_superkmers(
    ascii_slice: &[u8],
    k: usize,
    l: usize,
    offset: usize,
    canonical: bool,
    min_pos: &[u32],
    sk_pos: &[u32],
    results: &mut Vec<Superkmer>,
) {
    let seq_len = ascii_slice.len();
    let n = min_pos.len();
    if n == 0 { return; }

    results.reserve(n);
    for i in 0..n {
        let sk_start = sk_pos[i] as usize;
        let last_kmer = if i + 1 < n {
            sk_pos[i + 1] as usize - 1
        } else {
            seq_len - k
        };
        let size = last_kmer + k - sk_start;
        let mpos = min_pos[i] as usize - sk_start;

        let (mint, mint_is_rc) = if canonical {
            canonical_lmer(ascii_slice, min_pos[i] as usize, l)
        } else {
            let mut fwd = 0usize;
            let pos = min_pos[i] as usize;
            for j in 0..l {
                fwd = (fwd << 2) | encode_base(ascii_slice[pos + j]);
            }
            (fwd as u32, false)
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

/// Reusable superkmer extractor. Create once, call `process()` per read.
/// Pre-creates the NtHasher and SIMD cache to avoid per-call overhead.
pub struct SuperkmerExtractor {
    superkmers: Vec<Superkmer>,
    storage: Vec<u64>,
    min_pos: Vec<u32>,
    sk_pos: Vec<u32>,
    hasher: simd_minimizers::seq_hash::NtHasher<true>,
    cache: simd_minimizers::Cache,
    k: usize,
    l: usize,
    w: usize,
    canonical: bool,
}

impl SuperkmerExtractor {
    pub fn new(k: usize, l: usize) -> Self {
        let w = k - l + 1;
        Self { superkmers: Vec::new(), storage: Vec::new(), min_pos: Vec::new(), sk_pos: Vec::new(), hasher: simd_minimizers::seq_hash::NtHasher::new(l), cache: simd_minimizers::Cache::default(), k, l, w, canonical: true }
    }

    pub fn non_canonical(k: usize, l: usize) -> Self {
        let w = k - l + 1;
        Self { superkmers: Vec::new(), storage: Vec::new(), min_pos: Vec::new(), sk_pos: Vec::new(), hasher: simd_minimizers::seq_hash::NtHasher::new(l), cache: simd_minimizers::Cache::default(), k, l, w, canonical: false }
    }

    pub fn process(&mut self, seq: &[u8]) -> &[Superkmer] {
        self.superkmers.clear();
        if seq.len() < self.k { return &self.superkmers; }
        let num_words = (seq.len() + 31) / 32;
        self.storage.resize(num_words, 0);
        crate::utils::bitpack_fragment_into(seq, &mut self.storage);
        self.min_pos.clear();
        self.sk_pos.clear();
        simd_minimizers::canonical_minimizers(self.l, self.w)
            .hasher(&self.hasher)
            .super_kmers(&mut self.sk_pos)
            .run_with_buf(AsciiSeq(seq), &mut self.min_pos, &mut self.cache);
        emit_superkmers(seq, self.k, self.l, 0, self.canonical, &self.min_pos, &self.sk_pos, &mut self.superkmers);
        &self.superkmers
    }

    pub fn process_with_n(&mut self, seq: &[u8]) -> &[Superkmer] {
        self.superkmers.clear();
        if seq.len() < self.k { return &self.superkmers; }
        let num_words = (seq.len() + 31) / 32;
        self.storage.resize(num_words, 0);
        crate::utils::bitpack_fragment_into(seq, &mut self.storage);
        let fragments = crate::utils::split_on_n(seq, self.k);
        for (offset, fragment) in &fragments {
            self.min_pos.clear();
            self.sk_pos.clear();
            simd_minimizers::canonical_minimizers(self.l, self.w)
                .hasher(&self.hasher)
                .super_kmers(&mut self.sk_pos)
                .run_with_buf(AsciiSeq(fragment), &mut self.min_pos, &mut self.cache);
            emit_superkmers(fragment, self.k, self.l, *offset, self.canonical, &self.min_pos, &self.sk_pos, &mut self.superkmers);
        }
        &self.superkmers
    }

    /// Access the 2-bit packed representation of the last processed sequence.
    pub fn storage(&self) -> &[u64] {
        &self.storage
    }
}

pub struct SuperkmersIterator {
    superkmers: Vec<Superkmer>,
    pos: usize,
}

impl SuperkmersIterator {
    pub fn new(seq: &[u8], k: usize, l: usize) -> Self {
        let mut ext = SuperkmerExtractor::new(k, l);
        let sks = ext.process(seq);
        SuperkmersIterator { superkmers: sks.to_vec(), pos: 0 }
    }

    pub fn new_with_n(seq: &[u8], k: usize, l: usize) -> Self {
        let mut ext = SuperkmerExtractor::new(k, l);
        let sks = ext.process_with_n(seq);
        SuperkmersIterator { superkmers: sks.to_vec(), pos: 0 }
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
        Some(self.superkmers[idx])
    }
}
