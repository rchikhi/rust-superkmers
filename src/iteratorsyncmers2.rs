// much taken from rust-debruijn
use crate::Superkmer;
use lazy_static::lazy_static;

const S: usize = 2; //syncmer's s parameter

lazy_static! {
static ref SYNCMER_SCORES_8: Vec<usize> = generate_syncmer_scores::<8>();
static ref SYNCMER_SCORES_9: Vec<usize> = generate_syncmer_scores::<9>();
}

fn syncmer_scores(l: usize) -> &'static [usize] {
    match l {
        8 => &SYNCMER_SCORES_8[..],
        9 => &SYNCMER_SCORES_9[..],
        _ => panic!("Unsupported l={} for syncmer scores", l),
    }
}

fn canonical_table(l: usize) -> &'static [(u32, bool)] {
    match l {
        8 => &*crate::CANONICAL_8,
        9 => &*crate::CANONICAL_9,
        10 => &*crate::CANONICAL_10,
        12 => &*crate::CANONICAL_12,
        _ => panic!("Unsupported l={} for canonical lookup", l),
    }
}


fn generate_syncmer_scores<const K: usize>() -> Vec<usize> {
    let mut scores = vec![0usize; 1 << (2 * K)];
    for kmer_int in 0..(1 << (2 * K)) {
        let kmer_bytes = {
            let mut bytes = [0u8; K];
            for (i, byte) in bytes.iter_mut().enumerate() {
                *byte = match (kmer_int >> (2 * (K - 1 - i))) & 3 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    3 => b'T',
                    _ => unreachable!(),
                };
            }
            bytes
        };
        let syncmer = crate::syncmers::find_syncmers(K as usize, S, &[0, K - S], None, &kmer_bytes);
        scores[kmer_int] = syncmer.is_empty() as usize;
    }
    scores
}

use crate::utils::bitpack_fragment;

/// Extract an l-mer value from the packed storage, MSB-first encoding
/// (matching debruijn's Kmer8::to_u64() convention).
/// Each u64 in `data` holds 32 bases, with base 0 at bits 63-62 (MSB).
fn get_kmer_value(data: &[u64], base_pos: usize, l: usize) -> usize {
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

/// MinPos tracks the best minimizer position within a k-mer window.
/// Ordering: lower val wins; on tie, higher pos (rightmost) wins.
/// This matches debruijn's msp.rs MinPos exactly.
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
        // Reverse position: higher pos is "less" (preferred by min)
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


/// Run MSP sliding window on a single fragment, returning (kmer_start, minimizer_pos, minimizer_kmer, fragment_end).
/// All positions are absolute (offset already added).
/// `scores` maps each l-mer integer value to its score (lower = better minimizer).
fn msp_minimizer_positions(storage: &[u64], frag_len: usize, k: usize, l: usize, offset: usize, scores: &[usize]) -> Vec<(usize, usize, usize, usize)> {
    let mp = |pos: usize| -> MinPos {
        let kmer = get_kmer_value(storage, pos, l);
        let val = scores[kmer];
        MinPos { val, pos, kmer }
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
    min_positions: Vec<(usize, usize, usize, usize)>, // (kmer_start, minimizer_pos, minimizer_kmer, fragment_end)
    p: usize,
    k: usize,
    l: usize,
    canonical: bool,
}

impl SuperkmersIterator {
    /// Process a sequence assumed to contain no N characters.
    pub fn new(seq_str: &[u8], k: usize, l: usize) -> (Vec<u64>, Self) {
        let storage = bitpack_fragment(seq_str);
        let scores = syncmer_scores(l);
        let min_positions = msp_minimizer_positions(&storage, seq_str.len(), k, l, 0, scores);

        (storage, SuperkmersIterator {
            min_positions,
            p: 0,
            k,
            l,
            canonical: true,
        })
    }

    /// Process a sequence that may contain N/n characters.
    /// Splits on N's so superkmers never span across them.
    pub fn new_with_n(seq_str: &[u8], k: usize, l: usize) -> (Vec<u64>, Self) {
        let fragments = crate::utils::split_on_n(seq_str, k);
        let full_storage = bitpack_fragment(seq_str);
        let scores = syncmer_scores(l);

        let mut all_min_positions: Vec<(usize, usize, usize, usize)> = Vec::new();

        for (offset, fragment) in &fragments {
            let frag_storage = bitpack_fragment(fragment);
            let positions = msp_minimizer_positions(&frag_storage, fragment.len(), k, l, *offset, scores);
            all_min_positions.extend(positions);
        }

        (full_storage, SuperkmersIterator {
            min_positions: all_min_positions,
            p: 0,
            k,
            l,
            canonical: true,
        })
    }

    /// Disable canonical mint (use forward-strand l-mer, rc=false).
    pub fn non_canonical(mut self) -> Self {
        self.canonical = false;
        self
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
                // Same fragment: superkmer extends to overlap with next k-mer
                next_pos + self.k - 1 - start_pos
            } else {
                // Last superkmer in this fragment
                frag_end - start_pos
            }
        } else {
            frag_end - start_pos
        };

        self.p += 1;

        let (mint, rc) = if self.canonical {
            canonical_table(self.l)[min_kmer]
        } else {
            (min_kmer as u32, false)
        };
        Some(Superkmer {
            start: start_pos,
            mint,
            size: size as u8,
            mpos: (min_abs_pos - start_pos) as u8,
            rc,
        })
    }
}
