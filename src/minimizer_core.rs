//! Generic minimizer sliding window core.
//!
//! Score-table-agnostic sliding window minimum via monotonic deque.
//! Used by iteratorsyncmers2, iteratoruhs, and other minimizer-based iterators.

use crate::Superkmer;

/// Extract a single base (2 bits) from bitpacked storage at the given position.
#[inline(always)]
pub fn get_base(data: &[u64], pos: usize) -> usize {
    let word = pos / 32;
    let bit_offset = (31 - (pos % 32)) * 2;
    ((data[word] >> bit_offset) & 3) as usize
}

/// Extract an l-mer value from the packed storage, MSB-first encoding
/// (matching debruijn's Kmer8::to_u64() convention).
/// Each u64 in `data` holds 32 bases, with base 0 at bits 63-62 (MSB).
#[inline(always)]
pub fn get_kmer_value(data: &[u64], base_pos: usize, l: usize) -> usize {
    let word = base_pos / 32;
    let offset = base_pos % 32; // base offset within the word
    let bit_len = l * 2;

    if offset + l <= 32 {
        // All bases fit in one word
        let shift = 64 - (offset + l) * 2;
        ((data[word] >> shift) & ((1u64 << bit_len) - 1)) as usize
    } else {
        // Spans two words
        let bases_in_first = 32 - offset;
        let bases_in_second = l - bases_in_first;

        // High bits: last bases_in_first bases of the first word (lowest bits)
        let first_bits = (data[word] & ((1u64 << (bases_in_first * 2)) - 1)) as usize;

        // Low bits: first bases_in_second bases of the second word (highest bits)
        let shift = 64 - bases_in_second * 2;
        let second_bits = ((data[word + 1] >> shift) & ((1u64 << (bases_in_second * 2)) - 1)) as usize;

        (first_bits << (bases_in_second * 2)) | second_bits
    }
}

/// Look up the canonical table for a given l-mer length.
pub fn canonical_table(l: usize) -> &'static [(u32, bool)] {
    match l {
        7 => &*crate::CANONICAL_7,
        8 => &*crate::CANONICAL_8,
        9 => &*crate::CANONICAL_9,
        10 => &*crate::CANONICAL_10,
        12 => &*crate::CANONICAL_12,
        _ => panic!("Unsupported l={} for canonical lookup", l),
    }
}

/// Pack (score, position) into a single usize for block-decomposition comparison.
/// High 32 bits: compressed score. Low 32 bits: position (or complement for CLASSICAL).
/// Comparing packed values directly gives correct minimizer ordering.
#[inline(always)]
fn pack_score_pos<const CLASSICAL: bool>(score: usize, pos: usize) -> usize {
    // Compress score to 32 bits: syncmer bit (bit 32 of score) → bit 31, low 31 bits of hash.
    // For mspxor l=8: 16 bits of tiebreaker entropy (full resolution).
    // For mspxor l=9: 18 bits of tiebreaker entropy (full resolution).
    // For binary scores (0/1): maps to 0 and 1 in score32.
    let score32 = ((score >> 32) << 31) | (score & 0x7FFF_FFFF);
    if CLASSICAL {
        // Rightmost wins on tie: complement position so min() picks higher pos
        (score32 << 32) | (0xFFFF_FFFF - pos)
    } else {
        // Leftmost wins on tie: lower position = lower packed value
        (score32 << 32) | pos
    }
}

/// Extract position from a packed (score, position) value.
#[inline(always)]
fn unpack_pos<const CLASSICAL: bool>(packed: usize) -> usize {
    if CLASSICAL {
        0xFFFF_FFFF - (packed & 0xFFFF_FFFF)
    } else {
        packed & 0xFFFF_FFFF
    }
}

/// Block-decomposition sliding window minimum (no deque).
///
/// Divides positions into blocks of size w. Precomputes prefix-min and suffix-min
/// within each block. Window minimum = min(suffix[i], prefix[i+w-1]).
/// All passes are sequential reads/writes — no random access, no unpredictable branches.
///
/// CLASSICAL=true: rightmost wins on tie.
/// CLASSICAL=false: leftmost wins on tie.
///
/// Emits (kmer_start, minimizer_pos, minimizer_kmer, fragment_end) tuples.
/// All positions are absolute (offset already added).
#[inline(always)]
pub fn minimizer_positions_deque<const CLASSICAL: bool>(
    storage: &[u64], frag_len: usize, k: usize, l: usize, offset: usize,
    scores: &[usize], min_positions: &mut Vec<(usize, usize, usize, usize)>,
    scores_buf: &mut Vec<usize>, deque: &mut Vec<usize>,
) {
    let frag_end = offset + frag_len;
    if frag_len < k { return; }

    let num_lmers = frag_len - l + 1;
    let w = k - l + 1;
    let num_kmers = frag_len - k + 1;
    let mask = (1usize << (l * 2)) - 1;

    // Reuse buffers: scores_buf → packed values then prefix_min, deque → suffix_min
    if scores_buf.len() < num_lmers { scores_buf.resize(num_lmers, 0); }
    if deque.len() < num_lmers { deque.resize(num_lmers, 0); }

    // Pass 1: compute packed (score, position) values via rolling kmer
    unsafe {
        let mut rolling = get_kmer_value(storage, 0, l);
        *scores_buf.get_unchecked_mut(0) = pack_score_pos::<CLASSICAL>(*scores.get_unchecked(rolling), 0);
        for pos in 1..num_lmers {
            let new_base = get_base(storage, pos + l - 1);
            rolling = ((rolling << 2) | new_base) & mask;
            *scores_buf.get_unchecked_mut(pos) = pack_score_pos::<CLASSICAL>(*scores.get_unchecked(rolling), pos);
        }
    }

    // Pass 2: build suffix_min (right-to-left within each block of size w)
    {
        let mut block_start = 0;
        while block_start < num_lmers {
            let block_end = std::cmp::min(block_start + w, num_lmers);
            unsafe {
                *deque.get_unchecked_mut(block_end - 1) = *scores_buf.get_unchecked(block_end - 1);
                for i in (block_start..block_end - 1).rev() {
                    let cur = *scores_buf.get_unchecked(i);
                    let next = *deque.get_unchecked(i + 1);
                    *deque.get_unchecked_mut(i) = std::cmp::min(cur, next);
                }
            }
            block_start += w;
        }
    }

    // Pass 3: build prefix_min in-place over scores_buf (left-to-right within each block)
    {
        let mut block_start = 0;
        while block_start < num_lmers {
            let block_end = std::cmp::min(block_start + w, num_lmers);
            for i in block_start + 1..block_end {
                unsafe {
                    let prev = *scores_buf.get_unchecked(i - 1);
                    let cur = *scores_buf.get_unchecked(i);
                    *scores_buf.get_unchecked_mut(i) = std::cmp::min(prev, cur);
                }
            }
            block_start += w;
        }
    }

    // Pass 4: query — window min = min(suffix[i], prefix[i+w-1])
    unsafe {
        let first = std::cmp::min(*deque.get_unchecked(0), *scores_buf.get_unchecked(w - 1));
        let mut prev_min_pos = unpack_pos::<CLASSICAL>(first);
        let min_kmer = get_kmer_value(storage, prev_min_pos, l);
        min_positions.push((offset, prev_min_pos + offset, min_kmer, frag_end));

        for i in 1..num_kmers {
            let s = *deque.get_unchecked(i);
            let p = *scores_buf.get_unchecked(i + w - 1);
            let window_min = std::cmp::min(s, p);
            let cur_min_pos = unpack_pos::<CLASSICAL>(window_min);
            if cur_min_pos != prev_min_pos {
                let min_kmer = get_kmer_value(storage, cur_min_pos, l);
                min_positions.push((i + offset, cur_min_pos + offset, min_kmer, frag_end));
                prev_min_pos = cur_min_pos;
            }
        }
    }
}

/// Sticky sliding window minimum: single-pass with rescan on falloff.
/// On entry: new l-mer only replaces if strictly less (ties keep current).
/// On falloff: rescan picks rightmost among ties (stays in window longest).
/// This maximizes superkmer length for binary/low-cardinality score tables.
#[inline(always)]
pub fn minimizer_positions_sticky(
    storage: &[u64], frag_len: usize, k: usize, l: usize, offset: usize,
    scores: &[usize], min_positions: &mut Vec<(usize, usize, usize, usize)>,
) {
    let frag_end = offset + frag_len;
    if frag_len < k { return; }

    let w = k - l + 1;
    let mask = (1usize << (l * 2)) - 1;

    // Find rightmost minimum in a range
    let find_min = |start: usize, stop: usize| -> (usize, usize, usize) {
        let mut best_val = usize::MAX;
        let mut best_pos = start;
        let mut best_kmer = 0;
        for pos in start..=stop {
            let kmer = get_kmer_value(storage, pos, l);
            let val = scores[kmer];
            if val <= best_val {  // <= means rightmost wins on tie
                best_val = val;
                best_pos = pos;
                best_kmer = kmer;
            }
        }
        (best_val, best_pos, best_kmer)
    };

    let (mut min_val, mut min_pos, mut min_kmer) = find_min(0, w - 1);
    min_positions.push((offset, min_pos + offset, min_kmer, frag_end));

    let mut rolling_kmer = get_kmer_value(storage, w - 1, l);

    for i in 1..(frag_len - k + 1) {
        let new_base = get_base(storage, i + k - 1);
        rolling_kmer = ((rolling_kmer << 2) | new_base) & mask;
        let end_val = scores[rolling_kmer];
        let end_pos = i + w - 1;

        if i > min_pos {
            // Current minimizer fell off — rescan
            let (v, p, km) = find_min(i, end_pos);
            min_val = v;
            min_pos = p;
            min_kmer = km;
            min_positions.push((i + offset, min_pos + offset, min_kmer, frag_end));
        } else if end_val < min_val {
            // Strictly better l-mer entered from right
            min_val = end_val;
            min_pos = end_pos;
            min_kmer = rolling_kmer;
            min_positions.push((i + offset, min_pos + offset, min_kmer, frag_end));
        }
    }
}

/// Convert min_positions to superkmers.
pub fn materialize_superkmers(min_positions: &[(usize, usize, usize, usize)], k: usize, l: usize, canonical: bool, out: &mut Vec<Superkmer>) {
    let canon_table = if canonical { Some(canonical_table(l)) } else { None };
    for p in 0..min_positions.len() {
        let (start_pos, min_abs_pos, min_kmer, frag_end) = min_positions[p];
        let size = if p < min_positions.len() - 1 {
            let (next_pos, _, _, next_frag_end) = min_positions[p + 1];
            if next_frag_end == frag_end {
                next_pos + k - 1 - start_pos
            } else {
                frag_end - start_pos
            }
        } else {
            frag_end - start_pos
        };
        let (mint, mint_is_rc) = if let Some(table) = canon_table {
            table[min_kmer]
        } else {
            (min_kmer as u32, false)
        };
        out.push(Superkmer {
            start: start_pos,
            mint,
            size: size as u16,
            mpos: (min_abs_pos - start_pos) as u16,
            mint_is_rc,
        });
    }
}
