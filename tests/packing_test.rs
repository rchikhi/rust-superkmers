use rand::Rng;
use rust_superkmers::utils::{
    convert_bases, pack_32_bases, bitpack_fragment, bitpack_fragment_rc,
    bitpack_fragment_scalar, bitpack_fragment_rc_scalar, encode_base,
};

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T', b'T' => b'A',
        b'C' => b'G', b'G' => b'C',
        b'a' => b't', b't' => b'a',
        b'c' => b'g', b'g' => b'c',
        _ => b'A',
    }
}

fn rc(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Scalar reference implementation for verifying AVX2 packing.
fn scalar_pack_32(bases: &[u8]) -> u64 {
    assert!(bases.len() == 32);
    let mut packed = 0u64;
    for &b in bases {
        packed = (packed << 2) | encode_base(b) as u64;
    }
    packed
}

#[test]
fn test_convert_and_pack_known() {
    if !is_x86_feature_detected!("avx2") { return; }

    // All A's → 0
    let all_a = [b'A'; 32];
    let result = unsafe {
        let (conv, valid) = convert_bases(&all_a);
        assert!(valid);
        pack_32_bases(conv)
    };
    assert_eq!(result, 0u64);

    // All T's → all 11 bits
    let all_t = [b'T'; 32];
    let result = unsafe {
        let (conv, valid) = convert_bases(&all_t);
        assert!(valid);
        pack_32_bases(conv)
    };
    assert_eq!(result, 0xFFFF_FFFF_FFFF_FFFF);

    // All C's → all 01 bits
    let all_c = [b'C'; 32];
    let result = unsafe {
        let (conv, valid) = convert_bases(&all_c);
        assert!(valid);
        pack_32_bases(conv)
    };
    assert_eq!(result, 0x5555_5555_5555_5555);

    // All G's → all 10 bits
    let all_g = [b'G'; 32];
    let result = unsafe {
        let (conv, valid) = convert_bases(&all_g);
        assert!(valid);
        pack_32_bases(conv)
    };
    assert_eq!(result, 0xAAAA_AAAA_AAAA_AAAA);

    // ACGT repeated → 00_01_10_11 = 0x1B repeated
    let mut acgt_rep = [0u8; 32];
    for i in 0..8 {
        acgt_rep[i * 4] = b'A';
        acgt_rep[i * 4 + 1] = b'C';
        acgt_rep[i * 4 + 2] = b'G';
        acgt_rep[i * 4 + 3] = b'T';
    }
    let result = unsafe {
        let (conv, valid) = convert_bases(&acgt_rep);
        assert!(valid);
        pack_32_bases(conv)
    };
    assert_eq!(result, 0x1B1B_1B1B_1B1B_1B1B);
}

#[test]
fn test_convert_and_pack_random() {
    if !is_x86_feature_detected!("avx2") { return; }

    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::rng();

    for _ in 0..1000 {
        let mut seq = [0u8; 32];
        for b in seq.iter_mut() {
            *b = bases[rng.random_range(0..4)];
        }
        let expected = scalar_pack_32(&seq);
        let result = unsafe {
            let (conv, valid) = convert_bases(&seq);
            assert!(valid);
            pack_32_bases(conv)
        };
        assert_eq!(result, expected, "mismatch for {:?}", std::str::from_utf8(&seq).unwrap());
    }
}

#[test]
fn test_convert_bases_case_insensitive() {
    if !is_x86_feature_detected!("avx2") { return; }

    let upper = b"ACGTACGTACGTACGTacgtacgtacgtacgt";
    let all_upper = b"ACGTACGTACGTACGTACGTACGTACGTACGT";

    let result_mixed = unsafe {
        let (conv, valid) = convert_bases(upper);
        assert!(valid);
        pack_32_bases(conv)
    };
    let result_upper = unsafe {
        let (conv, valid) = convert_bases(all_upper);
        assert!(valid);
        pack_32_bases(conv)
    };
    assert_eq!(result_mixed, result_upper);
}

#[test]
fn test_convert_bases_invalid_chars() {
    if !is_x86_feature_detected!("avx2") { return; }

    // N should be flagged as invalid
    let mut with_n = [b'A'; 32];
    with_n[0] = b'N';
    let (_, valid) = unsafe { convert_bases(&with_n) };
    assert!(!valid, "N should be flagged as invalid");

    // lowercase n
    with_n[0] = b'n';
    let (_, valid) = unsafe { convert_bases(&with_n) };
    assert!(!valid, "n should be flagged as invalid");

    // Various invalid characters
    for &invalid in &[b'X', b'0', b' ', b'.', b'-'] {
        let mut seq = [b'A'; 32];
        seq[15] = invalid;
        let (_, valid) = unsafe { convert_bases(&seq) };
        assert!(!valid, "byte {:?} should be flagged as invalid", invalid as char);
    }

    // All valid → valid=true
    let valid_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
    let (_, valid) = unsafe { convert_bases(valid_seq) };
    assert!(valid);
}

#[test]
fn test_convert_bases_n_mapped_to_zero() {
    if !is_x86_feature_detected!("avx2") { return; }

    // N's should produce 0 (same as A) in the encoded output
    let mut with_n = [b'A'; 32];
    with_n[0] = b'N';
    let result_n = unsafe {
        let (conv, _) = convert_bases(&with_n);
        pack_32_bases(conv)
    };
    let all_a = [b'A'; 32];
    let result_a = unsafe {
        let (conv, _) = convert_bases(&all_a);
        pack_32_bases(conv)
    };
    assert_eq!(result_n, result_a, "N should encode the same as A (0)");
}

#[test]
fn test_bitpack_fragment_various_lengths() {
    if !is_x86_feature_detected!("avx2") { return; }

    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::rng();

    for len in [1, 15, 31, 32, 33, 50, 63, 64, 65, 100, 128] {
        let seq: Vec<u8> = (0..len).map(|_| bases[rng.random_range(0..4)]).collect();
        let packed = bitpack_fragment(&seq);

        assert_eq!(packed.len(), (len + 31) / 32, "wrong number of words for len={}", len);

        // Verify each base can be extracted correctly
        for (i, &b) in seq.iter().enumerate() {
            let word = i / 32;
            let offset = i % 32;
            let shift = 62 - offset * 2;
            let extracted = ((packed[word] >> shift) & 3) as u8;
            assert_eq!(extracted, encode_base(b),
                "base mismatch at position {} for len={}", i, len);
        }
    }
}

#[test]
fn test_bitpack_fragment_padding() {
    if !is_x86_feature_detected!("avx2") { return; }

    // Fragment of 5 bases: padded positions should be 0 (A)
    let seq = b"TTTTT";
    let packed = bitpack_fragment(seq);
    assert_eq!(packed.len(), 1);

    for i in 0..5 {
        let shift = 62 - i * 2;
        assert_eq!((packed[0] >> shift) & 3, 3, "base {} should be T", i);
    }
    for i in 5..32 {
        let shift = 62 - i * 2;
        assert_eq!((packed[0] >> shift) & 3, 0, "padding base {} should be A", i);
    }
}

#[test]
fn test_bitpack_fragment_with_n() {
    if !is_x86_feature_detected!("avx2") { return; }

    let mut seq = [b'C'; 32];
    seq[10] = b'N';
    seq[20] = b'N';
    let packed = bitpack_fragment(&seq);

    // Position 10 and 20 should be 0 (A), rest should be 1 (C)
    for i in 0..32 {
        let shift = 62 - i * 2;
        let expected = if i == 10 || i == 20 { 0 } else { 1 };
        assert_eq!((packed[0] >> shift) & 3, expected,
            "base {} should be {}", i, if expected == 0 { "A(from N)" } else { "C" });
    }
}

#[test]
fn test_bitpack_fragment_random_with_invalid() {
    if !is_x86_feature_detected!("avx2") { return; }

    let mut rng = rand::rng();
    let bases = [b'A', b'C', b'G', b'T'];
    let invalid = [b'N', b'n', b'X', b'.', b'0'];

    for _ in 0..100 {
        let len = rng.random_range(32..=128);
        let mut seq: Vec<u8> = (0..len).map(|_| bases[rng.random_range(0..4)]).collect();

        // Sprinkle in some invalid characters
        let n_invalid = rng.random_range(1..=5);
        for _ in 0..n_invalid {
            let pos = rng.random_range(0..len);
            seq[pos] = invalid[rng.random_range(0..invalid.len())];
        }

        let packed = bitpack_fragment(&seq);
        assert_eq!(packed.len(), (len + 31) / 32);

        // Invalid chars should map to 0 (A), valid chars should match
        for (i, &b) in seq.iter().enumerate() {
            let word = i / 32;
            let offset = i % 32;
            let shift = 62 - offset * 2;
            let extracted = ((packed[word] >> shift) & 3) as u8;
            assert_eq!(extracted, encode_base(b),
                "base mismatch at position {} (byte={:?})", i, b as char);
        }
    }
}

// --- bitpack_fragment_rc tests ---

#[test]
fn test_bitpack_fragment_rc_known() {
    // RC of "AAAA..." (32 A's) = "TTTT..." (32 T's)
    let all_a = [b'A'; 32];
    assert_eq!(bitpack_fragment_rc(&all_a), vec![0xFFFF_FFFF_FFFF_FFFF]);

    // RC of "TTTT..." = "AAAA..."
    let all_t = [b'T'; 32];
    assert_eq!(bitpack_fragment_rc(&all_t), vec![0x0000_0000_0000_0000]);

    // RC of "CCCC..." = "GGGG..."
    let all_c = [b'C'; 32];
    assert_eq!(bitpack_fragment_rc(&all_c), vec![0xAAAA_AAAA_AAAA_AAAA]);

    // RC of "GGGG..." = "CCCC..."
    let all_g = [b'G'; 32];
    assert_eq!(bitpack_fragment_rc(&all_g), vec![0x5555_5555_5555_5555]);
}

#[test]
fn test_bitpack_fragment_rc_matches_naive() {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::rng();

    for _ in 0..1000 {
        let len = rng.random_range(1..=200);
        let seq: Vec<u8> = (0..len).map(|_| bases[rng.random_range(0..4)]).collect();

        let rc_seq = rc(&seq);
        let expected = bitpack_fragment(&rc_seq);
        let result = bitpack_fragment_rc(&seq);

        assert_eq!(result, expected,
            "RC mismatch for len={}, seq={}", len, std::str::from_utf8(&seq).unwrap());
    }
}

#[test]
fn test_bitpack_fragment_rc_short() {
    // Single base
    assert_eq!(bitpack_fragment_rc(b"A"), bitpack_fragment(b"T"));
    assert_eq!(bitpack_fragment_rc(b"C"), bitpack_fragment(b"G"));

    // A few bases
    assert_eq!(bitpack_fragment_rc(b"ACGT"), bitpack_fragment(b"ACGT")); // ACGT is its own RC
    assert_eq!(bitpack_fragment_rc(b"AACCC"), bitpack_fragment(&rc(b"AACCC")));
}

#[test]
fn test_bitpack_fragment_rc_word_boundary() {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::rng();

    // Test lengths around 32-byte word boundaries
    for len in [31, 32, 33, 63, 64, 65, 95, 96, 97] {
        for _ in 0..50 {
            let seq: Vec<u8> = (0..len).map(|_| bases[rng.random_range(0..4)]).collect();
            let expected = bitpack_fragment(&rc(&seq));
            let result = bitpack_fragment_rc(&seq);
            assert_eq!(result, expected, "RC mismatch at boundary len={}", len);
        }
    }
}

#[test]
fn test_bitpack_fragment_rc_double_rc_is_identity() {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::rng();

    for _ in 0..200 {
        let len = rng.random_range(1..=150);
        let seq: Vec<u8> = (0..len).map(|_| bases[rng.random_range(0..4)]).collect();

        let forward = bitpack_fragment(&seq);
        let rc_of_rc = bitpack_fragment_rc(&rc(&seq));

        assert_eq!(forward, rc_of_rc, "double RC should be identity for len={}", len);
    }
}

// --- Scalar path tests (run even on AVX2 machines) ---

#[test]
fn test_scalar_matches_avx2() {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::rng();

    for _ in 0..500 {
        let len = rng.random_range(1..=200);
        let seq: Vec<u8> = (0..len).map(|_| bases[rng.random_range(0..4)]).collect();

        let avx2_result = bitpack_fragment(&seq);
        let scalar_result = bitpack_fragment_scalar(&seq);
        assert_eq!(avx2_result, scalar_result, "scalar vs avx2 mismatch for len={}", len);
    }
}

#[test]
fn test_scalar_rc_matches_avx2_rc() {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::rng();

    for _ in 0..500 {
        let len = rng.random_range(1..=200);
        let seq: Vec<u8> = (0..len).map(|_| bases[rng.random_range(0..4)]).collect();

        let avx2_result = bitpack_fragment_rc(&seq);
        let scalar_result = bitpack_fragment_rc_scalar(&seq);
        assert_eq!(avx2_result, scalar_result, "scalar RC vs avx2 RC mismatch for len={}", len);
    }
}

#[test]
fn test_scalar_rc_matches_naive() {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::rng();

    for _ in 0..500 {
        let len = rng.random_range(1..=200);
        let seq: Vec<u8> = (0..len).map(|_| bases[rng.random_range(0..4)]).collect();

        let rc_seq = rc(&seq);
        let expected = bitpack_fragment_scalar(&rc_seq);
        let result = bitpack_fragment_rc_scalar(&seq);
        assert_eq!(result, expected, "scalar RC mismatch for len={}", len);
    }
}

#[test]
fn test_scalar_known_values() {
    assert_eq!(bitpack_fragment_scalar(&[b'A'; 32]), vec![0u64]);
    assert_eq!(bitpack_fragment_scalar(&[b'T'; 32]), vec![0xFFFF_FFFF_FFFF_FFFF]);
    assert_eq!(bitpack_fragment_scalar(&[b'C'; 32]), vec![0x5555_5555_5555_5555]);
    assert_eq!(bitpack_fragment_scalar(&[b'G'; 32]), vec![0xAAAA_AAAA_AAAA_AAAA]);

    assert_eq!(bitpack_fragment_rc_scalar(&[b'A'; 32]), vec![0xFFFF_FFFF_FFFF_FFFF]);
    assert_eq!(bitpack_fragment_rc_scalar(&[b'T'; 32]), vec![0u64]);
}
