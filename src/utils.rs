pub fn revcomp(seq: &str) -> String
{
    std::str::from_utf8(&bio::alphabets::dna::revcomp(seq.as_bytes())).unwrap().to_string()
}

/// Split a byte sequence on N/n characters, returning (offset, fragment) pairs.
/// Fragments shorter than `min_len` are skipped.
pub fn split_on_n(seq: &[u8], min_len: usize) -> Vec<(usize, &[u8])> {
    let mut fragments = Vec::new();
    let mut start = 0;
    for (i, &b) in seq.iter().enumerate() {
        if b == b'N' || b == b'n' {
            if i - start >= min_len {
                fragments.push((start, &seq[start..i]));
            }
            start = i + 1;
        }
    }
    if seq.len() - start >= min_len {
        fragments.push((start, &seq[start..]));
    }
    fragments
}

pub fn superkmer_to_verbose(superkmer :crate::Superkmer, read: &String, l :usize) -> crate::SuperkmerVerbose
{
        let sequence = &read[superkmer.start..superkmer.start+(superkmer.size as usize)];
        let sequence = if superkmer.mint_is_rc { revcomp(sequence) } else { sequence.to_string() };

        let mpos = superkmer.mpos as usize;
        let minimizer = sequence[mpos..mpos+l].to_string();

        let superkmer_verbose = crate::SuperkmerVerbose {
            sequence,
            minimizer,
            mpos: mpos.try_into().unwrap()
        };
        superkmer_verbose
}

// --- AVX2 2-bit DNA encoding utilities ---

#[allow(clippy::wildcard_imports)]
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use std::arch::x86_64::*;

/// Pack the lowest 2 bits of each byte of a 32-byte __m256i into a u64.
/// All bytes must have values in 0..3; incorrect results otherwise.
/// The first byte becomes the highest 2 bits in the output (MSB-first).
///
/// # Safety
/// Requires AVX2 support.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
pub unsafe fn pack_32_bases(bases: __m256i) -> u64 {
    let reverse_mask = _mm256_set_epi8(
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        12, 13, 14, 15,
    );

    let reversed = _mm256_shuffle_epi8(bases, reverse_mask);
    let permuted = _mm256_permute4x64_epi64(reversed, 0b01_11_00_10);

    let first_bits = _mm256_slli_epi16(permuted, 7);
    let second_bits = _mm256_slli_epi16(permuted, 6);

    let lo_half = _mm256_unpacklo_epi8(first_bits, second_bits);
    let hi_half = _mm256_unpackhi_epi8(first_bits, second_bits);

    let packed_lo = (_mm256_movemask_epi8(lo_half) as u32) as u64;
    let packed_hi = (_mm256_movemask_epi8(hi_half) as u32) as u64;

    (packed_hi << 32) | packed_lo
}

/// Convert a slice of exactly 32 ASCII ACGT bytes into a __m256i of 2-bit encoded bases.
/// Returns `(encoded, valid)` where `valid` is false if any byte is not in [aAcCgGtT].
/// Invalid bytes are silently mapped to A (0).
///
/// Encoding: A=0, C=1, G=2, T=3 (case-insensitive).
///
/// # Safety
/// Requires AVX2 support. `bytes` must have length exactly 32.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
pub unsafe fn convert_bases(bytes: &[u8]) -> (__m256i, bool) {
    assert!(bytes.len() == 32);

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

    let lut = _mm256_set_epi8(
        0, 0, 0, 0, 0, 0, 0, 0, /* G */ 2, 0, 0, /* T */ 3, /* C */ 1, 0,
        /* A */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* G */ 2, 0, 0, /* T */ 3,
        /* C */ 1, 0, /* A */ 0, 0,
    );

    let input = _mm256_loadu_si256(bytes.as_ptr() as *const __m256i);

    let hi = _mm256_and_si256(_mm256_srli_epi16(input, 3), lo_mask);
    let hi_lookup = _mm256_shuffle_epi8(hi_lut, hi);

    let lo_lookup = _mm256_shuffle_epi8(lo_lut, input);

    let mask = _mm256_cmpeq_epi8(
        _mm256_and_si256(lo_lookup, hi_lookup),
        _mm256_setzero_si256(),
    );
    let valid = _mm256_testc_si256(_mm256_setzero_si256(), mask) != 0;

    let shuffled = _mm256_shuffle_epi8(lut, input);
    let res = _mm256_andnot_si256(mask, shuffled);

    (res, valid)
}

/// Bit-pack an ASCII DNA fragment into a Vec of u64 words.
/// Each u64 holds 32 bases in MSB-first 2-bit encoding (A=0, C=1, G=2, T=3).
/// The last word is padded with A's if the fragment length is not a multiple of 32.
/// Requires AVX2 support at runtime.
/// Encode a single ASCII base to its 2-bit representation.
/// Invalid characters (including N) map to 0 (same as A).
#[inline]
pub fn encode_base(b: u8) -> u8 {
    match b | 0x20 { // lowercase
        b'a' => 0,
        b'c' => 1,
        b'g' => 2,
        b't' => 3,
        _ => 0,
    }
}

/// Scalar packing: pack up to 32 ASCII bases into a u64, MSB-first.
/// Bases beyond the slice are padded with A (0).
pub fn scalar_pack_chunk(bases: &[u8]) -> u64 {
    let mut packed = 0u64;
    for i in 0..32 {
        let val = if i < bases.len() { encode_base(bases[i]) as u64 } else { 0 };
        packed = (packed << 2) | val;
    }
    packed
}

/// Pack a 32-byte buffer and complement the 2-bit values.
/// `count` is the number of real bases (rest is padding that stays 0 after complement).
#[inline]
fn pack_and_complement(buffer: &[u8; 32], count: usize) -> u64 {
    let complement_mask = if count >= 32 { !0u64 } else { !0u64 << ((32 - count) * 2) };

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") {
            let packed = unsafe {
                let (conv, _) = convert_bases(buffer);
                pack_32_bases(conv)
            };
            return packed ^ complement_mask;
        }
    }

    scalar_pack_chunk(buffer) ^ complement_mask
}

/// Scalar-only bit-packing. Same output as `bitpack_fragment` but never uses AVX2.
/// Useful for testing the scalar path on machines that have AVX2.
pub fn bitpack_fragment_scalar(fragment: &[u8]) -> Vec<u64> {
    fragment.chunks(32).map(|chunk| scalar_pack_chunk(chunk)).collect()
}

/// Scalar-only RC bit-packing. Same output as `bitpack_fragment_rc` but never uses AVX2.
pub fn bitpack_fragment_rc_scalar(fragment: &[u8]) -> Vec<u64> {
    let n = fragment.len();
    let num_words = (n + 31) / 32;
    let mut storage = Vec::with_capacity(num_words);
    for w in 0..num_words {
        let rc_start = w * 32;
        let count = std::cmp::min(32, n - rc_start);
        let mut buffer = [b'A'; 32];
        for j in 0..count {
            buffer[j] = fragment[n - 1 - rc_start - j];
        }
        let complement_mask = if count >= 32 { !0u64 } else { !0u64 << ((32 - count) * 2) };
        storage.push(scalar_pack_chunk(&buffer) ^ complement_mask);
    }
    storage
}

/// Bit-pack an ASCII DNA fragment into a caller-provided buffer.
/// Each u64 holds 32 bases in MSB-first 2-bit encoding (A=0, C=1, G=2, T=3).
/// The last word is padded with A's if the fragment length is not a multiple of 32.
/// Uses AVX2 when available, falls back to scalar otherwise.
///
/// `out` must have length >= `(fragment.len() + 31) / 32`.
pub fn bitpack_fragment_into(fragment: &[u8], out: &mut [u64]) {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") {
            let mut buffer = [b'A'; 32];
            for (i, chunk) in fragment.chunks(32).enumerate() {
                if chunk.len() < 32 {
                    buffer[..chunk.len()].copy_from_slice(chunk);
                    buffer[chunk.len()..].fill(b'A');
                    let (conv_chunk, _) = unsafe { convert_bases(&buffer) };
                    out[i] = unsafe { pack_32_bases(conv_chunk) };
                } else {
                    let (conv_chunk, _) = unsafe { convert_bases(chunk) };
                    out[i] = unsafe { pack_32_bases(conv_chunk) };
                }
            }
            return;
        }
    }

    // Scalar fallback
    for (i, chunk) in fragment.chunks(32).enumerate() {
        out[i] = scalar_pack_chunk(chunk);
    }
}

/// Bit-pack an ASCII DNA fragment into a Vec of u64 words.
/// Each u64 holds 32 bases in MSB-first 2-bit encoding (A=0, C=1, G=2, T=3).
/// The last word is padded with A's if the fragment length is not a multiple of 32.
/// Uses AVX2 when available, falls back to scalar otherwise.
pub fn bitpack_fragment(fragment: &[u8]) -> Vec<u64> {
    let num_words = (fragment.len() + 31) / 32;
    let mut storage = vec![0u64; num_words];
    bitpack_fragment_into(fragment, &mut storage);
    storage
}

/// Bit-pack the reverse complement of an ASCII DNA fragment into a caller-provided buffer.
/// Equivalent to `bitpack_fragment_into(rc(fragment), out)` but without materializing the RC string.
/// Each u64 holds 32 bases in MSB-first 2-bit encoding (A=0, C=1, G=2, T=3).
/// The last word is padded with A's if the fragment length is not a multiple of 32.
///
/// `out` must have length >= `(fragment.len() + 31) / 32`.
pub fn bitpack_fragment_rc_into(fragment: &[u8], out: &mut [u64]) {
    let n = fragment.len();
    let num_words = (n + 31) / 32;

    for w in 0..num_words {
        let rc_start = w * 32;
        let count = std::cmp::min(32, n - rc_start);
        let mut buffer = [b'A'; 32];
        for j in 0..count {
            buffer[j] = fragment[n - 1 - rc_start - j];
        }
        // Pack reversed bases, then XOR to complement the 2-bit values.
        // In 2-bit encoding (A=0,C=1,G=2,T=3), XOR 3 gives A↔T, C↔G.
        out[w] = pack_and_complement(&buffer, count);
    }
}

/// Bit-pack the reverse complement of an ASCII DNA fragment into a Vec of u64 words.
/// Equivalent to `bitpack_fragment(rc(fragment))` but without materializing the RC string.
/// Each u64 holds 32 bases in MSB-first 2-bit encoding (A=0, C=1, G=2, T=3).
/// The last word is padded with A's if the fragment length is not a multiple of 32.
pub fn bitpack_fragment_rc(fragment: &[u8]) -> Vec<u64> {
    let n = fragment.len();
    let num_words = (n + 31) / 32;
    let mut storage = vec![0u64; num_words];
    bitpack_fragment_rc_into(fragment, &mut storage);
    storage
}
