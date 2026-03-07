use std::fs;
use debruijn::dna_string::DnaString;
use rand::Rng;
use rust_superkmers::Superkmer;
use rust_superkmers::utils::split_on_n;

fn random_dna(len: usize) -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::rng();
    (0..len).map(|_| bases[rng.random_range(0..4)] as char).collect()
}

/// Compare iteratorsyncmers2 vs iteratorsyncmersmsp on a given sequence,
/// asserting that all fields match exactly.
fn assert_implementations_match(seq: &str, k: usize, l: usize) {
    let dnastring = DnaString::from_dna_string(seq).to_bytes();
    let iter_msp = rust_superkmers::iteratorsyncmersmsp::SuperkmersIterator::new(&dnastring, k, l);
    let msp: Vec<Superkmer> = iter_msp.collect();

    let (_, iter2) = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(seq.as_bytes(), k, l);
    let avx2: Vec<Superkmer> = iter2.collect();

    assert_eq!(msp.len(), avx2.len(),
        "Count mismatch for len={}, k={}, l={}: msp={}, avx2={}",
        seq.len(), k, l, msp.len(), avx2.len());

    for (i, (a, b)) in msp.iter().zip(avx2.iter()).enumerate() {
        assert_eq!(a.start, b.start,
            "start differs at superkmer {} for len={}, k={}, l={}", i, seq.len(), k, l);
        assert_eq!(a.size, b.size,
            "size differs at superkmer {} for len={}, k={}, l={}", i, seq.len(), k, l);
        assert_eq!(a.mint, b.mint,
            "mint differs at superkmer {} for len={}, k={}, l={}: msp={:#x}, avx2={:#x}",
            i, seq.len(), k, l, a.mint, b.mint);
        assert_eq!(a.mpos, b.mpos,
            "mpos differs at superkmer {} for len={}, k={}, l={}", i, seq.len(), k, l);
    }
}

/// Verify that superkmers tile the sequence: every k-mer is covered exactly once.
fn assert_tiling(superkmers: &[Superkmer], seq_len: usize, k: usize) {
    let num_kmers = seq_len - k + 1;
    let mut covered = vec![false; num_kmers];

    for sk in superkmers {
        let sk_start = sk.start;
        let sk_end = sk.start + sk.size as usize;
        assert!(sk_end <= seq_len,
            "superkmer extends past sequence: start={}, size={}, seq_len={}",
            sk.start, sk.size, seq_len);
        for j in sk_start..=(sk_end - k) {
            assert!(!covered[j],
                "k-mer at position {} covered by multiple superkmers", j);
            covered[j] = true;
        }
    }

    for (i, &c) in covered.iter().enumerate() {
        assert!(c, "k-mer at position {} not covered by any superkmer", i);
    }
}

#[test]
fn test_compare_ecoli_220() {
    let contents = fs::read_to_string("tests/ecoli.genome.220.fa")
        .expect("Failed to read test genome file")
        .split("\n").collect::<Vec<&str>>()[1].to_string();

    for k in [17, 21, 31, 41] {
        assert_implementations_match(&contents, k, 8);
    }
}

#[test]
fn test_compare_ecoli_76() {
    let contents = fs::read_to_string("tests/ecoli.genome.76.fa")
        .expect("Failed to read test genome file")
        .split("\n").collect::<Vec<&str>>()[1].to_string();

    for k in [17, 21, 31] {
        assert_implementations_match(&contents, k, 8);
    }
}

#[test]
fn test_compare_ecoli_100k() {
    let contents = fs::read_to_string("tests/ecoli.genome.100k.fa")
        .expect("Failed to read test genome file")
        .split("\n").collect::<Vec<&str>>()[1].to_string();

    assert_implementations_match(&contents, 21, 8);
    assert_implementations_match(&contents, 31, 8);
}

#[test]
fn test_tiling_ecoli_220() {
    let contents = fs::read_to_string("tests/ecoli.genome.220.fa")
        .expect("Failed to read test genome file")
        .split("\n").collect::<Vec<&str>>()[1].to_string();

    for k in [17, 21, 31] {
        let (_, iter) = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(
            contents.as_bytes(), k, 8);
        let superkmers: Vec<Superkmer> = iter.collect();
        assert_tiling(&superkmers, contents.len(), k);
    }
}

#[test]
fn test_random_sequences_various_lengths() {
    for k in [17, 21, 31] {
        let l = 8;
        // Test lengths near k, around word boundary (32), and longer
        for len in [k, k + 1, k + 5, 32, 33, 50, 63, 64, 65, 100, 128, 150, 200, 300, 500] {
            if len < k { continue; }
            for _ in 0..10 {
                let seq = random_dna(len);
                assert_implementations_match(&seq, k, l);
            }
        }
    }
}

#[test]
fn test_random_sequences_tiling() {
    for k in [17, 21, 31] {
        let l = 8;
        for len in [k + 5, 64, 150, 300] {
            for _ in 0..20 {
                let seq = random_dna(len);
                let (_, iter) = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(
                    seq.as_bytes(), k, l);
                let superkmers: Vec<Superkmer> = iter.collect();
                assert_tiling(&superkmers, len, k);
            }
        }
    }
}

#[test]
fn test_edge_case_exact_k_length() {
    let l = 8;
    for k in [17, 21, 31] {
        let seq = random_dna(k);
        assert_implementations_match(&seq, k, l);

        let (_, iter) = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(
            seq.as_bytes(), k, l);
        let superkmers: Vec<Superkmer> = iter.collect();
        assert_eq!(superkmers.len(), 1,
            "sequence of exactly k={} should produce exactly 1 superkmer", k);
        assert_eq!(superkmers[0].start, 0);
        assert_eq!(superkmers[0].size as usize, k);
    }
}

#[test]
fn test_homopolymer_sequences() {
    let l = 8;
    for k in [17, 21] {
        for base in ['A', 'C', 'G', 'T'] {
            let seq: String = std::iter::repeat(base).take(100).collect();
            assert_implementations_match(&seq, k, l);
        }
    }
}

#[test]
fn test_dinucleotide_repeats() {
    let l = 8;
    for k in [17, 21] {
        for pattern in ["AC", "GT", "AG", "CT"] {
            let seq: String = pattern.chars().cycle().take(150).collect();
            assert_implementations_match(&seq, k, l);
        }
    }
}

#[test]
fn test_word_boundary_crossing() {
    // Sequences where l-mers straddle the 32-base word boundary in AVX2 storage
    let l = 8;
    let k = 21;
    // Lengths 25..=40 force l-mers to cross the boundary at position 32
    for len in 50..=70 {
        for _ in 0..10 {
            let seq = random_dna(len);
            assert_implementations_match(&seq, k, l);
        }
    }
}

#[test]
fn test_long_random_sequences() {
    let l = 8;
    for k in [21, 31] {
        for len in [1000, 5000, 10000] {
            let seq = random_dna(len);
            assert_implementations_match(&seq, k, l);

            let (_, iter) = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(
                seq.as_bytes(), k, l);
            let superkmers: Vec<Superkmer> = iter.collect();
            assert_tiling(&superkmers, len, k);
        }
    }
}

#[test]
fn test_split_on_n_basic() {
    let seq = b"ACGTNNNTACG";
    let frags = split_on_n(seq, 1);
    assert_eq!(frags.len(), 2);
    assert_eq!(frags[0], (0, &b"ACGT"[..]));
    assert_eq!(frags[1], (7, &b"TACG"[..]));

    // min_len filtering
    let frags = split_on_n(seq, 5);
    assert_eq!(frags.len(), 0);

    // no N's
    let seq = b"ACGTACGT";
    let frags = split_on_n(seq, 1);
    assert_eq!(frags.len(), 1);
    assert_eq!(frags[0], (0, &b"ACGTACGT"[..]));
}

#[test]
fn test_n_splitting_iteratorsyncmers2() {
    let k = 21;
    let l = 8;

    // Build a sequence with N's in the middle
    let left = random_dna(50);
    let right = random_dna(50);
    let seq_with_n = format!("{}NNNNN{}", left, right);

    // new_with_n on the combined sequence
    let (_, iter) = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new_with_n(
        seq_with_n.as_bytes(), k, l);
    let superkmers_split: Vec<Superkmer> = iter.collect();

    // Process each fragment independently with new (no N)
    let (_, iter_left) = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(
        left.as_bytes(), k, l);
    let left_sks: Vec<Superkmer> = iter_left.collect();

    let (_, iter_right) = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(
        right.as_bytes(), k, l);
    let right_sks: Vec<Superkmer> = iter_right.collect();

    // Same number of superkmers
    assert_eq!(superkmers_split.len(), left_sks.len() + right_sks.len(),
        "N-split should produce same superkmers as processing fragments independently");

    // Left fragment superkmers should match exactly
    for (a, b) in superkmers_split.iter().zip(left_sks.iter()) {
        assert_eq!(a.start, b.start);
        assert_eq!(a.size, b.size);
        assert_eq!(a.mint, b.mint);
        assert_eq!(a.mpos, b.mpos);
    }

    // Right fragment superkmers should match with offset
    let right_offset = left.len() + 5; // 5 N's
    for (a, b) in superkmers_split[left_sks.len()..].iter().zip(right_sks.iter()) {
        assert_eq!(a.start, b.start + right_offset);
        assert_eq!(a.size, b.size);
        assert_eq!(a.mint, b.mint);
        assert_eq!(a.mpos, b.mpos);
    }
}

#[test]
fn test_n_at_edges() {
    let k = 21;
    let l = 8;
    let dna = random_dna(60);

    // N's at start
    let seq = format!("NNN{}", dna);
    let (_, iter) = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new_with_n(
        seq.as_bytes(), k, l);
    let sks: Vec<Superkmer> = iter.collect();
    assert!(sks.iter().all(|sk| sk.start >= 3), "no superkmer should start before the N's");

    // N's at end
    let seq = format!("{}NNN", dna);
    let (_, iter) = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new_with_n(
        seq.as_bytes(), k, l);
    let sks: Vec<Superkmer> = iter.collect();
    assert!(sks.iter().all(|sk| (sk.start + sk.size as usize) <= dna.len()),
        "no superkmer should extend into trailing N's");

    // Fragment too short to produce k-mers
    let seq = format!("ACGTNNNNN{}", dna);
    let (_, iter) = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new_with_n(
        seq.as_bytes(), k, l);
    let sks: Vec<Superkmer> = iter.collect();
    // "ACGT" is too short for k=21, should only get superkmers from the right fragment
    assert!(sks.iter().all(|sk| sk.start >= 9));
}
