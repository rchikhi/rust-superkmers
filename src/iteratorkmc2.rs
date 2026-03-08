// KMC2-style signature-based superkmer iterator.
// Uses the lexicographically smallest canonical l-mer as minimizer,
// with KMC2 disqualification rules (no signatures starting with A
// or containing AA/AC) to improve bucket balance.
//
// Same MSP sliding-window approach as iteratorsyncmers2, just a different
// scoring function backed by a precomputed lookup table.

use crate::Superkmer;
use crate::utils::bitpack_fragment;
use lazy_static::lazy_static;

const K8: usize = 8;

lazy_static! {
    /// For each possible 8-mer value, stores the KMC2 score.
    /// Score = canonical_value for non-disqualified signatures,
    ///       = (1 << 16) | canonical_value for disqualified ones.
    /// This ensures non-disqualified signatures always sort before disqualified.
    /// The canonical value can be recovered as score & 0xFFFF.
    static ref KMC2_SCORES_8: Box<[usize; 1 << (2 * K8)]> = {
        let mut table = Box::new([0usize; 1 << (2 * K8)]);
        for val in 0..(1usize << (2 * K8)) {
            let rc = revcomp_lmer(val, K8);
            let canonical = std::cmp::min(val, rc);
            let disqualified = kmc2_is_disqualified(canonical, K8);
            table[val] = if disqualified { (1 << (2 * K8)) | canonical } else { canonical };
        }
        table
    };
}

#[inline]
fn revcomp_lmer(val: usize, l: usize) -> usize {
    let mut rc = 0usize;
    let mut v = val;
    for _ in 0..l {
        rc = (rc << 2) | (3 - (v & 3));
        v >>= 2;
    }
    rc
}

/// KMC2 disqualification: a canonical l-mer is disqualified if it
/// starts with A, or contains AA or AC as a dinucleotide anywhere.
#[inline]
fn kmc2_is_disqualified(canonical_val: usize, l: usize) -> bool {
    // Starts with A (top base = 0)
    if (canonical_val >> (2 * (l - 1))) & 3 == 0 {
        return true;
    }
    // Contains AA or AC dinucleotide
    for j in 0..l - 1 {
        // Extract dinucleotide (base_j, base_{j+1}) as 4 bits
        // base_j in bits 3:2, base_{j+1} in bits 1:0
        let dinuc = (canonical_val >> (2 * (l - 2 - j))) & 0xF;
        // AA = 0b0000 = 0, AC = 0b0001 = 1
        if dinuc <= 1 {
            return true;
        }
    }
    false
}

/// Extract an l-mer value from packed storage, MSB-first encoding.
/// Each u64 in `data` holds 32 bases, with base 0 at bits 63-62 (MSB).
fn get_kmer_value(data: &[u64], base_pos: usize, l: usize) -> usize {
    let word = base_pos / 32;
    let offset = base_pos % 32;
    let bit_len = l * 2;

    if offset + l <= 32 {
        let shift = 64 - (offset + l) * 2;
        ((data[word] >> shift) & ((1u64 << bit_len) - 1)) as usize
    } else {
        let bases_in_first = 32 - offset;
        let bases_in_second = l - bases_in_first;
        let first_bits = (data[word] & ((1u64 << (bases_in_first * 2)) - 1)) as usize;
        let shift = 64 - bases_in_second * 2;
        let second_bits = ((data[word + 1] >> shift) & ((1u64 << (bases_in_second * 2)) - 1)) as usize;
        (first_bits << (bases_in_second * 2)) | second_bits
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
struct MinPos {
    val: usize,
    pos: usize,
    kmer: usize,
}

impl Ord for MinPos {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let val_cmp = self.val.cmp(&other.val);
        if val_cmp != std::cmp::Ordering::Equal {
            return val_cmp;
        }
        // Tie: higher pos (rightmost) is preferred (sorts as "less")
        match self.pos.cmp(&other.pos) {
            std::cmp::Ordering::Equal => std::cmp::Ordering::Equal,
            std::cmp::Ordering::Less => std::cmp::Ordering::Greater,
            std::cmp::Ordering::Greater => std::cmp::Ordering::Less,
        }
    }
}

impl PartialOrd for MinPos {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

fn msp_minimizer_positions(storage: &[u64], frag_len: usize, k: usize, l: usize, offset: usize) -> Vec<(usize, usize, usize, usize)> {
    let l_mask = (1usize << (2 * l)) - 1;

    let mp = |pos: usize| -> MinPos {
        let lmer = get_kmer_value(storage, pos, l);
        let score = KMC2_SCORES_8[lmer];
        let canonical = score & l_mask;
        MinPos { val: score, pos, kmer: canonical }
    };

    let find_min = |start: usize, stop: usize| -> MinPos {
        let mut min_pos = mp(start);
        for pos in (start + 1)..=stop {
            let current = mp(pos);
            min_pos = std::cmp::min(min_pos, current);
        }
        min_pos
    };

    let frag_end = offset + frag_len;
    let mut min_positions: Vec<(usize, usize, usize, usize)> = Vec::new();

    if frag_len >= k {
        let mut min_pos = find_min(0, k - l);
        min_positions.push((offset, min_pos.pos + offset, min_pos.kmer, frag_end));

        for i in 1..(frag_len - k + 1) {
            let end_pos = mp(i + k - l);

            if i > min_pos.pos {
                min_pos = find_min(i, i + k - l);
                min_positions.push((i + offset, min_pos.pos + offset, min_pos.kmer, frag_end));
            } else if end_pos.val < min_pos.val {
                min_pos = end_pos;
                min_positions.push((i + offset, min_pos.pos + offset, min_pos.kmer, frag_end));
            }
        }
    }

    min_positions
}

pub struct SuperkmersIterator {
    min_positions: Vec<(usize, usize, usize, usize)>,
    p: usize,
    k: usize,
}

impl SuperkmersIterator {
    pub fn new(seq_str: &[u8], k: usize, l: usize) -> (Vec<u64>, Self) {
        let storage = bitpack_fragment(seq_str);
        let min_positions = msp_minimizer_positions(&storage, seq_str.len(), k, l, 0);

        (storage, SuperkmersIterator {
            min_positions,
            p: 0,
            k,
        })
    }

    pub fn new_with_n(seq_str: &[u8], k: usize, l: usize) -> (Vec<u64>, Self) {
        let fragments = crate::utils::split_on_n(seq_str, k);
        let full_storage = bitpack_fragment(seq_str);

        let mut all_min_positions: Vec<(usize, usize, usize, usize)> = Vec::new();

        for (offset, fragment) in &fragments {
            let frag_storage = bitpack_fragment(fragment);
            let positions = msp_minimizer_positions(&frag_storage, fragment.len(), k, l, *offset);
            all_min_positions.extend(positions);
        }

        (full_storage, SuperkmersIterator {
            min_positions: all_min_positions,
            p: 0,
            k,
        })
    }
}

impl Iterator for SuperkmersIterator {
    type Item = Superkmer;

    fn next(&mut self) -> Option<Self::Item> {
        if self.p >= self.min_positions.len() {
            return None;
        }

        let (start_pos, min_abs_pos, min_kmer, frag_end) = self.min_positions[self.p];

        let size = if self.p < self.min_positions.len() - 1 {
            let (next_pos, _, _, next_frag_end) = self.min_positions[self.p + 1];
            if next_frag_end == frag_end {
                next_pos + self.k - 1 - start_pos
            } else {
                frag_end - start_pos
            }
        } else {
            frag_end - start_pos
        };

        self.p += 1;

        Some(Superkmer {
            start: start_pos,
            mint: min_kmer as u32,
            size: size as u16,
            mpos: (min_abs_pos - start_pos) as u16,
            mint_is_rc: false,
        })
    }
}
