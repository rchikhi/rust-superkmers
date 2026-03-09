//! SIMD-accelerated syncmer-based superkmer iterator.
//!
//! Uses `simd-minimizers` to find canonical closed syncmer positions via SIMD,
//! then runs a rightmost-wins MSP sliding window to produce superkmers.
//! Requires odd l (e.g. l=9). Takes ASCII `&[u8]` input directly.
//! Mint is canonical by default; use `SuperkmersIterator::non_canonical()` for
//! forward-strand mint.
//!
//! ```ignore
//! use rust_superkmers::iteratorsimdmini::SuperkmersIterator;
//! let iter = SuperkmersIterator::new(seq, 31, 9);           // canonical
//! let iter = SuperkmersIterator::non_canonical(seq, 31, 9); // forward-strand
//! ```
use crate::Superkmer;

const SMER_SIZE: usize = 2; // syncmer's s parameter

/// Lookup table: ASCII byte -> 2-bit encoding (A=0, C=1, G=2, T=3)
const ASCII_TO_2BIT: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

/// Encode forward l-mer from ASCII bytes and look up canonical form in table.
#[inline]
fn canonical_lmer_index(ascii: &[u8], pos: usize, l: usize) -> (usize, bool) {
    let mut fwd = 0usize;
    for i in 0..l {
        let base = ASCII_TO_2BIT[ascii[pos + i] as usize] as usize;
        fwd = (fwd << 2) | base;
    }
    let canonical_table = match l {
        8 => &*crate::CANONICAL_8,
        9 => &*crate::CANONICAL_9,
        10 => &*crate::CANONICAL_10,
        12 => &*crate::CANONICAL_12,
        _ => panic!("Unsupported l={} for canonical lookup", l),
    };
    let (canonical, is_rc) = canonical_table[fwd];
    (canonical as usize, is_rc)
}

/// Encode forward l-mer from ASCII bytes (no canonicalization).
#[inline]
fn forward_lmer_index(ascii: &[u8], pos: usize, l: usize) -> usize {
    let mut fwd = 0usize;
    for i in 0..l {
        let base = ASCII_TO_2BIT[ascii[pos + i] as usize] as usize;
        fwd = (fwd << 2) | base;
    }
    fwd
}

/// Collect superkmers from a single N-free fragment using SIMD syncmer detection.
fn superkmers_from_fragment(
    ascii_slice: &[u8],
    k: usize,
    l: usize,
    offset: usize,
    canonical: bool,
    results: &mut Vec<Superkmer>,
) {
    use simd_minimizers::packed_seq::AsciiSeq;

    let seq_len = ascii_slice.len();
    if seq_len < k {
        return;
    }

    assert!(
        l % 2 == 1,
        "simd-mini canonical syncmers require odd l (try l=9), got l={}",
        l
    );

    // simd-minimizers requires uppercase ASCII
    let upper: Vec<u8> = ascii_slice.iter().map(|&b| b.to_ascii_uppercase()).collect();

    let w_sync = l - SMER_SIZE + 1;
    let mut syncmer_pos: Vec<u32> = Vec::new();
    simd_minimizers::canonical_closed_syncmers(SMER_SIZE, w_sync)
        .run(AsciiSeq(&upper), &mut syncmer_pos);

    // Mark all-A and all-T l-mer positions as demoted (valid syncmers but cause
    // hot buckets). Only used as minimizer when no other syncmer exists in window.
    let demoted: Vec<bool> = syncmer_pos.iter().map(|&pos| {
        let p = pos as usize;
        let slice = &upper[p..p + l];
        slice.iter().all(|&b| b == b'A') || slice.iter().all(|&b| b == b'T')
    }).collect();

    // Scan syncmer_pos[from..] within [lo, hi] and return the best position:
    // rightmost non-demoted, or rightmost demoted if no non-demoted exists.
    let scan_best = |from: usize, lo: usize, hi: usize| -> u32 {
        let mut best: u32 = u32::MAX;
        let mut best_demoted = true;
        let mut j = from;
        while j < syncmer_pos.len() && (syncmer_pos[j] as usize) <= hi {
            let p = syncmer_pos[j] as usize;
            if p >= lo {
                let d = demoted[j];
                if best == u32::MAX || (!d && best_demoted) || (d == best_demoted) {
                    best = syncmer_pos[j];
                    best_demoted = d;
                }
            }
            j += 1;
        }
        best
    };

    // MSP sliding window: track RIGHTMOST syncmer as minimizer.
    // Two-tier: prefer non-demoted during rescan; sticky between rescans.
    let num_kmers = seq_len - k + 1;
    let max_lmer_offset = k - l;
    let mut ptr = 0usize;
    let mut curr_min_pos: u32 = u32::MAX;
    let mut sk_start: usize = 0;

    // Initialize: find best syncmer in first k-mer window [0, k-l]
    if num_kmers > 0 {
        curr_min_pos = scan_best(0, 0, max_lmer_offset);
        while ptr < syncmer_pos.len() && (syncmer_pos[ptr] as usize) <= max_lmer_offset {
            ptr += 1;
        }
    }

    let mut i = 1usize;
    while i <= num_kmers {
        let need_new_min = if i < num_kmers {
            (curr_min_pos as usize) < i
        } else {
            true // sentinel: flush last superkmer
        };

        if need_new_min {
            if curr_min_pos != u32::MAX {
                let start = sk_start;
                let end = i - 1 + k;
                let size = end - start;
                let mpos_fwd = curr_min_pos as usize - start;

                let (mint, mint_is_rc) = if canonical {
                    canonical_lmer_index(ascii_slice, curr_min_pos as usize, l)
                } else {
                    (forward_lmer_index(ascii_slice, curr_min_pos as usize, l), false)
                };

                results.push(Superkmer {
                    start: start + offset,
                    mint: mint as u32,
                    size: size as u16,
                    mpos: mpos_fwd as u16,
                    mint_is_rc,
                });
            }

            sk_start = i;
            curr_min_pos = u32::MAX;
            if i < num_kmers {
                while ptr < syncmer_pos.len() && (syncmer_pos[ptr] as usize) < i {
                    ptr += 1;
                }
                curr_min_pos = scan_best(ptr, i, i + max_lmer_offset);
                // No syncmer in this window — jump to next syncmer position
                if curr_min_pos == u32::MAX && ptr < syncmer_pos.len() {
                    let next_sync = syncmer_pos[ptr] as usize;
                    // Jump i so next_sync is within the window [i, i+max_lmer_offset]
                    let jump_to = if next_sync > max_lmer_offset {
                        next_sync - max_lmer_offset
                    } else {
                        1
                    };
                    if jump_to > i {
                        sk_start = jump_to;
                        i = jump_to;
                        curr_min_pos = scan_best(ptr, i, i + max_lmer_offset);
                    }
                }
            }
        }
        i += 1;
    }
}

pub struct SuperkmersIterator {
    superkmers: Vec<Superkmer>,
    storage: Vec<u64>,
    pos: usize,
}

impl SuperkmersIterator {
    /// Process a sequence assumed to contain no N characters. Mint is canonical.
    pub fn new(seq: &[u8], k: usize, l: usize) -> Self {
        Self::new_inner(seq, k, l, true)
    }

    /// Process a sequence that may contain N/n characters. Mint is canonical.
    pub fn new_with_n(seq: &[u8], k: usize, l: usize) -> Self {
        Self::new_with_n_inner(seq, k, l, true)
    }

    /// Non-canonical version (forward-strand mint, mint_is_rc=false).
    pub fn non_canonical(seq: &[u8], k: usize, l: usize) -> Self {
        Self::new_inner(seq, k, l, false)
    }

    /// Non-canonical version with N-splitting.
    pub fn non_canonical_with_n(seq: &[u8], k: usize, l: usize) -> Self {
        Self::new_with_n_inner(seq, k, l, false)
    }

    /// Access the 2-bit packed representation of the input sequence.
    pub fn storage(&self) -> &[u64] {
        &self.storage
    }

    fn new_inner(seq: &[u8], k: usize, l: usize, canonical: bool) -> Self {
        let storage = crate::utils::bitpack_fragment(seq);
        let mut superkmers = Vec::new();
        superkmers_from_fragment(seq, k, l, 0, canonical, &mut superkmers);
        SuperkmersIterator { superkmers, storage, pos: 0 }
    }

    fn new_with_n_inner(seq: &[u8], k: usize, l: usize, canonical: bool) -> Self {
        let storage = crate::utils::bitpack_fragment(seq);
        let fragments = crate::utils::split_on_n(seq, k);
        let mut superkmers = Vec::new();
        for (offset, fragment) in &fragments {
            superkmers_from_fragment(fragment, k, l, *offset, canonical, &mut superkmers);
        }
        SuperkmersIterator { superkmers, storage, pos: 0 }
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
        // Move out of vec to avoid clone
        Some(Superkmer {
            start: self.superkmers[idx].start,
            mint: self.superkmers[idx].mint,
            size: self.superkmers[idx].size,
            mpos: self.superkmers[idx].mpos,
            mint_is_rc: self.superkmers[idx].mint_is_rc,
        })
    }
}

