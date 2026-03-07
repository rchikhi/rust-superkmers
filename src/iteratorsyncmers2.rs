// much taken from rust-debruijn
use crate::Superkmer;
use lazy_static::lazy_static;

const S: usize = 2; //syncmer's s parameter
const K8: usize = 8;

lazy_static! {
static ref SYNCMERS_8: Box<[bool; 1 << (2 * K8)]> = generate_syncmers::<K8>();
}

fn generate_syncmers<const K: usize>() -> Box<[bool; 1 << (2 * K)]> {
    let mut syncmers_arr = Box::new([false; 1 << (2 * K)]);
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
        syncmers_arr[kmer_int] = !syncmer.is_empty();
    }
    syncmers_arr
}

// bitops_avx2
#[allow(clippy::wildcard_imports)]
use std::arch::x86_64::*;

/// pack the lowest 2 bits of each byte of a 32 byte slice into a u64
/// all bytes must be have values in the range 0..3 incorrect results will be returned.
/// The first byte in the m256 will become the highest 2bits in the output, consistent
/// with the lexicographic sorting convention in this crate.
#[target_feature(enable = "avx2")]
pub(crate) unsafe fn pack_32_bases(bases: __m256i) -> u64 {
    // bases = d c b a

    let reverse_mask = _mm256_set_epi8(
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        12, 13, 14, 15,
    );

    // step 1: reverse all bytes within lanes
    // reversed = c d a b
    let reversed = _mm256_shuffle_epi8(bases, reverse_mask);

    // step 2: use a lane crossing permute to reverse lanes and
    // swap the middle two 64-bit chunks
    // permuted = a c b d
    let permuted = _mm256_permute4x64_epi64(reversed, 0b01_11_00_10);

    // step 3: interleave the bytes that contain the first and second bits
    let first_bits = _mm256_slli_epi16(permuted, 7);
    let second_bits = _mm256_slli_epi16(permuted, 6);

    // i(x) = interleave first and second bits of each byte of x
    // lo_half = i(c) i(d)
    let lo_half = _mm256_unpacklo_epi8(first_bits, second_bits);

    // hi_half = i(a) i(b)
    let hi_half = _mm256_unpackhi_epi8(first_bits, second_bits);

    // step 4: extract bits using movemask (zero extend)
    let packed_lo = (_mm256_movemask_epi8(lo_half) as u32) as u64;
    let packed_hi = (_mm256_movemask_epi8(hi_half) as u32) as u64;

    (packed_hi << 32) | packed_lo
}

/// Convert a slice of 32 ACGT bytes into a __m256i using 0-4 encoding.
/// Second element in tuple will be false if any of the bases are not in [aAcCgGtT].
/// Bases not in [aAcCgGtT] will be converted to A / 0.
#[target_feature(enable = "avx2")]
pub(crate) unsafe fn convert_bases(bytes: &[u8]) -> (__m256i, bool) {
    assert!(bytes.len() == 32);

    // a lookup table to map a number from 0..128 to a single bit
    let hi_lut = {
        let mut lut_hi = 0i64;
        lut_hi |= 1i64 << ((b'A' as i64) - 64i64);
        lut_hi |= 1i64 << ((b'C' as i64) - 64i64);
        lut_hi |= 1i64 << ((b'G' as i64) - 64i64);
        lut_hi |= 1i64 << ((b'T' as i64) - 64i64);
        lut_hi |= 1i64 << ((b'a' as i64) - 64i64);
        lut_hi |= 1i64 << ((b'c' as i64) - 64i64);
        lut_hi |= 1i64 << ((b'g' as i64) - 64i64);
        lut_hi |= 1i64 << ((b't' as i64) - 64i64);
        _mm256_set_epi64x(lut_hi, 0i64, lut_hi, 0i64)
    };

    let lo_lut = _mm256_set_epi8(
        1i8 << 7,
        1i8 << 6,
        1i8 << 5,
        1i8 << 4,
        1i8 << 3,
        1i8 << 2,
        1i8 << 1,
        1i8 << 0,
        1i8 << 7,
        1i8 << 6,
        1i8 << 5,
        1i8 << 4,
        1i8 << 3,
        1i8 << 2,
        1i8 << 1,
        1i8 << 0,
        1i8 << 7,
        1i8 << 6,
        1i8 << 5,
        1i8 << 4,
        1i8 << 3,
        1i8 << 2,
        1i8 << 1,
        1i8 << 0,
        1i8 << 7,
        1i8 << 6,
        1i8 << 5,
        1i8 << 4,
        1i8 << 3,
        1i8 << 2,
        1i8 << 1,
        1i8 << 0,
    );

    let lo_mask = _mm256_set1_epi8(0b00001111);

    // convert ACGT (case-insensitive) to a 2-bit representation by looking up low 4 bits
    let lut = _mm256_set_epi8(
        0, 0, 0, 0, 0, 0, 0, 0, /* G */ 2, 0, 0, /* T */ 3, /* C */ 1, 0,
        /* A */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* G */ 2, 0, 0, /* T */ 3,
        /* C */ 1, 0, /* A */ 0, 0,
    );

    let input = _mm256_loadu_si256(bytes.as_ptr() as *const __m256i);

    // step 1: lookup a byte using the high 4 bits of each character
    let hi = _mm256_and_si256(_mm256_srli_epi16(input, 3), lo_mask);
    let hi_lookup = _mm256_shuffle_epi8(hi_lut, hi);

    // step 2: convert the low 3 bits of each character to a mask
    // byte x is converted to (1 << x)
    let lo_lookup = _mm256_shuffle_epi8(lo_lut, input);

    // step 3: AND the byte derived from the high 4 bits and the mask derived from the low 3 bits
    let mask = _mm256_cmpeq_epi8(
        _mm256_and_si256(lo_lookup, hi_lookup),
        _mm256_setzero_si256(),
    );
    let valid = _mm256_testc_si256(_mm256_setzero_si256(), mask) != 0;

    // step 4: use lookup table to convert nucleotides to their 2-bit representation,
    // while zeroing out invalid characters (shuffle returns a zero byte if MSB of a byte is 1)
    let shuffled = _mm256_shuffle_epi8(lut, input);
    let res = _mm256_andnot_si256(mask, shuffled);

    (res, valid)
}

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

fn bitpack_fragment(fragment: &[u8]) -> Vec<u64> {
    let mut storage = Vec::new();
    let mut buffer = [b'A'; 32];
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") {
            for chunk in fragment.chunks(32) {
                if chunk.len() < 32 {
                    buffer[..chunk.len()].copy_from_slice(chunk);
                    let (conv_chunk, _) = unsafe { convert_bases(&buffer) };
                    let packed = unsafe { pack_32_bases(conv_chunk) };
                    storage.push(packed);
                } else {
                    let (conv_chunk, _) = unsafe { convert_bases(chunk) };
                    let packed = unsafe { pack_32_bases(conv_chunk) };
                    storage.push(packed);
                }
            }
        }
    }
    storage
}

/// Run MSP sliding window on a single fragment, returning (kmer_start, minimizer_pos, minimizer_kmer, fragment_end).
/// All positions are absolute (offset already added).
fn msp_minimizer_positions(storage: &[u64], frag_len: usize, k: usize, l: usize, offset: usize) -> Vec<(usize, usize, usize, usize)> {
    let score = |kmer_val: usize| -> usize { !SYNCMERS_8[kmer_val] as usize };

    let mp = |pos: usize| -> MinPos {
        let kmer = get_kmer_value(storage, pos, l);
        let val = score(kmer);
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
}

impl SuperkmersIterator {
    /// Process a sequence assumed to contain no N characters.
    pub fn new(seq_str: &[u8], k: usize, l: usize) -> (Vec<u64>, Self) {
        let storage = bitpack_fragment(seq_str);
        let min_positions = msp_minimizer_positions(&storage, seq_str.len(), k, l, 0);

        (storage, SuperkmersIterator {
            min_positions,
            p: 0,
            k,
        })
    }

    /// Process a sequence that may contain N/n characters.
    /// Splits on N's so superkmers never span across them.
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

        Some(Superkmer {
            start: start_pos,
            mint: min_kmer as u32,
            size: size as u8,
            mpos: (min_abs_pos - start_pos) as u8,
            rc: false,
        })
    }
}
