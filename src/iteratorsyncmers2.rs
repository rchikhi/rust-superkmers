// much taken from rust-debruijn
use crate::Superkmer;
//use static_init::{dynamic};
use lazy_static::lazy_static;
//use debruijn::dna_string::DnaString;

const K8: usize = 8;
const K10: usize = 10;
const K12: usize = 12;

lazy_static! {
static ref SYNCMERS_8: Box<[bool; 1 << (2 * K8)]> = generate_syncmers::<K8>();
}
/*#[dynamic]
static SYNCMERS_10: Box<[bool; 1 << (2 * K10)]> = generate_syncmers::<K10>();
#[dynamic]
static SYNCMERS_12: Box<[bool; 1 << (2 * K12)]> = generate_syncmers::<K12>();
*/

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
        let syncmer = crate::syncmers::find_syncmers(K as usize, 2, &[0, K - 2], None, &kmer_bytes);
        syncmers_arr[kmer_int] = !syncmer.is_empty();
    }
    syncmers_arr
}


#[inline(always)]
pub fn base_to_bits(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0u8,
        b'C' | b'c' => 1u8,
        b'G' | b'g' => 2u8,
        b'T' | b't' => 3u8,
        _ => 0u8,
    }
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
    // TODO: this test doesn't work

    // step 4: use lookup table to convert nucleotides to their 2-bit representation,
    // while zeroing out invalid characters (shuffle returns a zero byte if MSB of a byte is 1)
    let shuffled = _mm256_shuffle_epi8(lut, input);
    let res = _mm256_andnot_si256(mask, shuffled);

    (res, valid)
}


pub struct SuperkmersIterator {
    min_positions : Vec<(usize, u16)>,
	p : usize,
    k :usize,
    l :usize,
    m :usize,
}

fn get_storage_window(data: &[u64], pos: usize) -> usize {

    let idx = pos / 64;
    let offset = pos % 64;
    let mut window: usize;

    let window_size = 16;
    if offset + window_size > 64 {
        // Window spans across two u64 entries
        let bits_in_first = 64 - offset;
        let bits_in_second = window_size - bits_in_first;
        // Extract bits from the first u64
        window = ((data[idx] >> offset) & ((1 << bits_in_first) - 1)) as usize;
        // Shift left and extract bits from the second u64
        window |= ((data[idx + 1] & ((1 << bits_in_second) - 1)) as usize) << bits_in_first;
    } else {
        // All bits are within the current u64
        window = ((data[idx] >> offset) & ((1 << window_size) - 1)) as usize;
    }
    window
}


impl<'a> SuperkmersIterator {
    pub fn new(seq_str: &[u8], k: usize, l: usize) -> (Vec<u64>, Self) {
        let m = seq_str.len();
		let mut storage: Vec<u64> = Vec::new();
		let mut len = 0;
        let mut _has_N = false;
        let mut buffer = [b'A'; 32];
		// Accelerated avx2 mode. Should run on most machines made since 2013.
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                /*let seq_str_transmuted = unsafe {
                    // Cast raw pointer to a slice of size N+2
                    std::slice::from_raw_parts(seq_str.as_ptr(), seq_str.len()+(32-seq_str.len()%32))
                };*/
                for chunk in seq_str.chunks(32) {
                    if chunk.len() < 32 {
                        // If the chunk is less than 32, create a temporary vector and resize it
                        buffer[..chunk.len()].copy_from_slice(chunk);
                        // Perform operations on the adjusted chunk
                        let (conv_chunk, has_N_chunk) = unsafe { convert_bases(&buffer) };
                        _has_N |= has_N_chunk;
                        let packed = unsafe { pack_32_bases(conv_chunk) };
                        storage.push(packed);
                    } else {
                        // Perform operations directly on the original chunk
                        let (conv_chunk, has_N_chunk) = unsafe { convert_bases(chunk) };
                        _has_N |= has_N_chunk;
                        /*if has_N // TODO: fix the validity check in convert_bases, it's not right
                        {println!("N's in chunk {}: {}", i, std::str::from_utf8(chunk).unwrap());
                        }*/
                        let packed = unsafe { pack_32_bases(conv_chunk) };
                        storage.push(packed);
                    }
                }
                len = m;
            }
        }
        /*
		let mut Npos: Vec<usize> = Vec::new();
        if has_N {
            //println!("investigating N's in {}", std::str::from_utf8(seq_str).unwrap());
            for i in 0..seq_str.len() {
                if seq_str[i] == b'N' {
                    Npos.push(len);
                }
            }
            // TODO do something with Npos
            if Npos.len() > 0 {
                println!("ignoring seq containing N's: {:?} len {}",seq_str, len);
                //min_positions = Vec::new(); 
            }
        }*/

        // populate min_positions
        let mut min_positions = Vec::with_capacity(((1.3* (((m-k+1)*2)/(k-l+1)) as f32)) as usize); //expected closed syncmer density * 1.3

        let mut i = 0;
        while i < m/2 as usize {
            let window = get_storage_window(&storage, i*2);
            if SYNCMERS_8[window] {
                min_positions.push((i, window as u16));
                i += k-7;
            }
            i += 1
		}
        println!("its2 minimizer pos {:?}",min_positions);
        (storage, SuperkmersIterator { 
			min_positions,
			p: 0,
            k,
            l,
            m
		})
    }
}

impl<'a> Iterator for SuperkmersIterator {
    type Item = Superkmer;

    fn next(&mut self) -> Option<Self::Item> {
        if self.p >= self.min_positions.len() { return None; }
        if self.p == 0 {
            if self.m < self.k { return None; }
            let (start_pos, mint) = self.min_positions[0];
            let (next_pos, _) = self.min_positions[1];
            self.p += 1;
            return Some(Superkmer {
                start: 0 as usize,
                mint: mint as u32,
                size: (next_pos as usize + self.k - 1 - start_pos as usize) as u8,
                mpos: start_pos as u8,
                rc: false,
            });
        }
        else {
            if self.p < self.min_positions.len() - 1 {
                let (start_pos, mint) = self.min_positions[self.p];
                let (_next_pos, _) = self.min_positions[self.p+1];
                //let mut end : usize= next_pos as usize + self.k - 1;
                let mut end : usize = start_pos as usize + self.k as usize - self.l;
                if end >= self.m { end = (self.m-1) as usize; }
                let start = start_pos as usize - (self.k-self.l);
                self.p += 1;
                return Some(Superkmer {
                    start: start as usize,
                    mint: mint as u32,
                    size: (end as usize - start) as u8,
                    mpos: start_pos as u8,
                    rc: false,
                });
            }
            else {
                self.p += 1;
                return Some(Superkmer {
                    start: 0 as usize,
                    mint: 0 as u32,
                    size: 0 as u8,
                    mpos: 0 as u8,
                    rc: false,
                });
            }
        }
    }
}

