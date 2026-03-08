/// SIMD-accelerated syncmer-based superkmer iterator.
///
/// Uses `simd-minimizers` to find canonical closed syncmer positions via SIMD,
/// then runs a rightmost-wins MSP sliding window to produce superkmers.
/// Ported from rust-notbcalm3's simd-mini path.
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

/// Compute canonical l-mer index from ASCII bytes at given position.
/// Returns (canonical_index, is_rc) where is_rc means the RC was smaller.
#[inline]
fn canonical_lmer_index(ascii: &[u8], pos: usize, l: usize) -> (usize, bool) {
    let mut fwd = 0usize;
    for i in 0..l {
        let base = ASCII_TO_2BIT[ascii[pos + i] as usize] as usize;
        fwd = (fwd << 2) | base;
    }
    let mut rc = 0usize;
    for i in (0..l).rev() {
        let base = ASCII_TO_2BIT[ascii[pos + i] as usize] as usize;
        rc = (rc << 2) | (3 - base);
    }
    (fwd.min(rc), rc < fwd)
}

/// Collect superkmers from a single N-free fragment using SIMD syncmer detection.
fn superkmers_from_fragment(
    ascii_slice: &[u8],
    k: usize,
    l: usize,
    offset: usize,
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

    // MSP sliding window: track RIGHTMOST syncmer as minimizer.
    let num_kmers = seq_len - k + 1;
    let max_lmer_offset = k - l;
    let mut ptr = 0usize;
    let mut curr_min_pos: u32 = u32::MAX;
    let mut sk_start: usize = 0;

    // Initialize: find rightmost syncmer in first k-mer window [0, k-l]
    if num_kmers > 0 {
        while ptr < syncmer_pos.len() && (syncmer_pos[ptr] as usize) <= max_lmer_offset {
            curr_min_pos = syncmer_pos[ptr];
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

                let (mint, need_rc) = canonical_lmer_index(ascii_slice, curr_min_pos as usize, l);

                results.push(Superkmer {
                    start: start + offset,
                    mint: mint as u32,
                    size: size as u8,
                    mpos: mpos_fwd as u8,
                    rc: need_rc,
                });
            }

            sk_start = i;
            curr_min_pos = u32::MAX;
            if i < num_kmers {
                while ptr < syncmer_pos.len() && (syncmer_pos[ptr] as usize) < i {
                    ptr += 1;
                }
                let mut j = ptr;
                while j < syncmer_pos.len() && (syncmer_pos[j] as usize) <= i + max_lmer_offset {
                    curr_min_pos = syncmer_pos[j];
                    j += 1;
                }
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
                        let mut j2 = ptr;
                        while j2 < syncmer_pos.len() && (syncmer_pos[j2] as usize) <= i + max_lmer_offset {
                            curr_min_pos = syncmer_pos[j2];
                            j2 += 1;
                        }
                    }
                }
            }
        }
        i += 1;
    }
}

pub struct SuperkmersIterator {
    superkmers: Vec<Superkmer>,
    pos: usize,
}

impl SuperkmersIterator {
    /// Process a sequence assumed to contain no N characters.
    pub fn new(seq: &[u8], k: usize, l: usize) -> Self {
        let mut superkmers = Vec::new();
        superkmers_from_fragment(seq, k, l, 0, &mut superkmers);
        SuperkmersIterator { superkmers, pos: 0 }
    }

    /// Process a sequence that may contain N/n characters.
    /// Splits on N's so superkmers never span across them.
    pub fn new_with_n(seq: &[u8], k: usize, l: usize) -> Self {
        let fragments = crate::utils::split_on_n(seq, k);
        let mut superkmers = Vec::new();
        for (offset, fragment) in &fragments {
            superkmers_from_fragment(fragment, k, l, *offset, &mut superkmers);
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
        // Move out of vec to avoid clone
        Some(Superkmer {
            start: self.superkmers[idx].start,
            mint: self.superkmers[idx].mint,
            size: self.superkmers[idx].size,
            mpos: self.superkmers[idx].mpos,
            rc: self.superkmers[idx].rc,
        })
    }
}

