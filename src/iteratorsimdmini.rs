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

    // simd-minimizers handles both uppercase and lowercase ASCII (ACTGactg)
    let w_sync = l - SMER_SIZE + 1;
    let mut syncmer_pos: Vec<u32> = Vec::new();
    simd_minimizers::canonical_closed_syncmers(SMER_SIZE, w_sync)
        .run(AsciiSeq(ascii_slice), &mut syncmer_pos);

    // Check if the l-mer at a syncmer position is all-A or all-T (demoted).
    let is_demoted = |pos: u32| -> bool {
        let p = pos as usize;
        let b0 = ascii_slice[p] & 0xDF; // case-insensitive: clear bit 5
        (b0 == b'A' || b0 == b'T') && ascii_slice[p..p + l].iter().all(|&b| b & 0xDF == b0)
    };

    // MSP sliding window: track RIGHTMOST syncmer as minimizer.
    // After each rescan, if the rightmost is demoted (all-A/T), fall back to
    // the rightmost non-demoted in the same range. ~16% overhead on 1M random DNA.
    let num_kmers = seq_len - k + 1;
    let max_lmer_offset = k - l;
    let mut ptr = 0usize;
    let mut curr_min_pos: u32 = u32::MAX;
    let mut sk_start: usize = 0;

    // After taking rightmost in syncmer_pos[from..j_end], demote if needed.
    macro_rules! demote_check {
        ($from:expr, $j_end:expr) => {
            if curr_min_pos != u32::MAX && is_demoted(curr_min_pos) {
                let mut best = u32::MAX;
                let mut di = $from;
                while di < $j_end {
                    if !is_demoted(syncmer_pos[di]) { best = syncmer_pos[di]; }
                    di += 1;
                }
                if best != u32::MAX { curr_min_pos = best; }
            }
        };
    }

    // Initialize: find rightmost syncmer in first k-mer window [0, k-l]
    if num_kmers > 0 {
        let start_ptr = ptr;
        while ptr < syncmer_pos.len() && (syncmer_pos[ptr] as usize) <= max_lmer_offset {
            curr_min_pos = syncmer_pos[ptr];
            ptr += 1;
        }
        if ptr > start_ptr { demote_check!(start_ptr, ptr); }
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
                let scan_from = ptr;
                let mut j = ptr;
                while j < syncmer_pos.len() && (syncmer_pos[j] as usize) <= i + max_lmer_offset {
                    curr_min_pos = syncmer_pos[j];
                    j += 1;
                }
                if j > scan_from { demote_check!(scan_from, j); }
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
                        // Re-scan for syncmers in [jump_to, jump_to+max_lmer_offset]
                        let scan_from2 = ptr;
                        let mut j2 = ptr;
                        while j2 < syncmer_pos.len() && (syncmer_pos[j2] as usize) <= i + max_lmer_offset {
                            curr_min_pos = syncmer_pos[j2];
                            j2 += 1;
                        }
                        if j2 > scan_from2 { demote_check!(scan_from2, j2); }
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

