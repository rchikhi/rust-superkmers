use std::fs;
use debruijn::dna_string::DnaString;
use rust_superkmers::Superkmer;
use rust_superkmers::utils::split_on_n;

fn random_dna(len: usize, seed: u64) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut x = seed;
    (0..len).map(|_| {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
        bases[((x >> 33) % 4) as usize]
    }).collect()
}

/// Encode l-mer at position `pos` in `seq` to a 2-bit integer.
/// Uses debruijn convention: A=0, C=1, G=2, T=3, MSB-first.
fn encode_lmer(seq: &[u8], pos: usize, l: usize) -> usize {
    let mut v = 0usize;
    for i in 0..l {
        let base = match seq[pos + i] {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 0,
        };
        v = (v << 2) | base;
    }
    v
}

fn rc_lmer(val: usize, l: usize) -> usize {
    let mut rc = 0usize;
    let mut v = val;
    for _ in 0..l {
        rc = (rc << 2) | (3 - (v & 3));
        v >>= 2;
    }
    rc
}

fn canonical_lmer(seq: &[u8], pos: usize, l: usize) -> (usize, bool) {
    let fwd = encode_lmer(seq, pos, l);
    let rc = rc_lmer(fwd, l);
    if rc < fwd { (rc, true) } else { (fwd, false) }
}

/// Read a single-sequence FASTA file and return the sequence (second line).
fn read_fasta_seq(path: &str) -> Vec<u8> {
    fs::read_to_string(path)
        .expect("Failed to read test genome file")
        .split("\n").collect::<Vec<&str>>()[1]
        .as_bytes().to_vec()
}

/// Collect superkmers from iteratorsyncmers2 and assert they tile the sequence.
fn assert_syncmers2_tiling(seq: &[u8], k: usize, l: usize) {
    let iter = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(seq, k, l);
    let superkmers: Vec<Superkmer> = iter.collect();
    assert_tiling(&superkmers, seq.len(), k);
}

/// Assert two superkmer slices have identical fields, with an optional position offset
/// applied to the `expected` side (for comparing N-split results).
fn assert_superkmers_eq(actual: &[Superkmer], expected: &[Superkmer], offset: usize) {
    assert_eq!(actual.len(), expected.len(), "superkmer count mismatch");
    for (i, (a, b)) in actual.iter().zip(expected.iter()).enumerate() {
        assert_eq!(a.start, b.start + offset,
            "start differs at superkmer {}", i);
        assert_eq!(a.size, b.size,
            "size differs at superkmer {}", i);
        assert_eq!(a.mint, b.mint,
            "mint differs at superkmer {}", i);
        assert_eq!(a.mpos, b.mpos,
            "mpos differs at superkmer {}", i);
    }
}

/// Compare iteratorsyncmers2 vs iteratorsyncmersmsp on a given sequence,
/// asserting that all fields match exactly.
fn assert_implementations_match(seq: &[u8], k: usize, l: usize) {
    let seq_str = std::str::from_utf8(seq).unwrap();
    let dnastring = DnaString::from_dna_string(seq_str).to_bytes();
    let iter_msp = rust_superkmers::iteratorsyncmersmsp::SuperkmersIterator::new(&dnastring, k, l);
    let msp: Vec<Superkmer> = iter_msp.collect();

    let iter2 = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(seq, k, l);
    let avx2: Vec<Superkmer> = iter2.collect();

    assert_superkmers_eq(&msp, &avx2, 0);
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
    let seq = read_fasta_seq("tests/ecoli.genome.220.fa");
    for k in [17, 21, 31, 41] {
        assert_implementations_match(&seq, k, 8);
    }
}

#[test]
fn test_compare_ecoli_76() {
    let seq = read_fasta_seq("tests/ecoli.genome.76.fa");
    for k in [17, 21, 31] {
        assert_implementations_match(&seq, k, 8);
    }
}

#[test]
fn test_compare_ecoli_100k() {
    let seq = read_fasta_seq("tests/ecoli.genome.100k.fa");
    assert_implementations_match(&seq, 21, 8);
    assert_implementations_match(&seq, 31, 8);
}

#[test]
fn test_tiling_ecoli_220() {
    let seq = read_fasta_seq("tests/ecoli.genome.220.fa");
    for k in [17, 21, 31] {
        assert_syncmers2_tiling(&seq, k, 8);
    }
}

#[test]
fn test_random_sequences_various_lengths() {
    for k in [17, 21, 31] {
        let l = 8;
        // Test lengths near k, around word boundary (32), and longer
        for len in [k, k + 1, k + 5, 32, 33, 50, 63, 64, 65, 100, 128, 150, 200, 300, 500] {
            if len < k { continue; }
            for i in 0..10 {
                let seq = random_dna(len, (k * 1000 + len * 10 + i) as u64);
                assert_implementations_match(&seq, k, l);
            }
        }
    }
}

#[test]
fn test_random_sequences_tiling() {
    for k in [17, 21, 31] {
        for len in [k + 5, 64, 150, 300] {
            for i in 0..20 {
                let seq = random_dna(len, (k * 1000 + len * 100 + i) as u64);
                assert_syncmers2_tiling(&seq, k, 8);
            }
        }
    }
}

#[test]
fn test_edge_case_exact_k_length() {
    for k in [17, 21, 31] {
        let seq = random_dna(k, k as u64 * 7 + 13);
        assert_implementations_match(&seq, k, 8);

        let iter = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(
            &seq, k, 8);
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
        for &base in &[b'A', b'C', b'G', b'T'] {
            let seq: Vec<u8> = std::iter::repeat(base).take(100).collect();
            assert_implementations_match(&seq, k, l);
        }
    }
}

#[test]
fn test_dinucleotide_repeats() {
    let l = 8;
    for k in [17, 21] {
        for pattern in [b"AC", b"GT", b"AG", b"CT"] {
            let seq: Vec<u8> = pattern.iter().cycle().take(150).cloned().collect();
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
        for i in 0..10 {
            let seq = random_dna(len, (len * 100 + i) as u64);
            assert_implementations_match(&seq, k, l);
        }
    }
}

#[test]
fn test_long_random_sequences() {
    for k in [21, 31] {
        for len in [1000, 5000, 10000] {
            let seq = random_dna(len, (k * 10000 + len) as u64);
            assert_implementations_match(&seq, k, 8);
            assert_syncmers2_tiling(&seq, k, 8);
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
    let left = random_dna(50, 42);
    let right = random_dna(50, 137);
    let mut seq_with_n = Vec::new();
    seq_with_n.extend_from_slice(&left);
    seq_with_n.extend_from_slice(b"NNNNN");
    seq_with_n.extend_from_slice(&right);

    // new_with_n on the combined sequence
    let iter = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new_with_n(
        &seq_with_n, k, l);
    let superkmers_split: Vec<Superkmer> = iter.collect();

    // Process each fragment independently with new (no N)
    let iter_left = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(
        &left, k, l);
    let left_sks: Vec<Superkmer> = iter_left.collect();

    let iter_right = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(
        &right, k, l);
    let right_sks: Vec<Superkmer> = iter_right.collect();

    // Left fragment superkmers should match exactly
    assert_superkmers_eq(&superkmers_split[..left_sks.len()], &left_sks, 0);

    // Right fragment superkmers should match with offset
    let right_offset = left.len() + 5; // 5 N's
    assert_superkmers_eq(&superkmers_split[left_sks.len()..], &right_sks, right_offset);
}

#[test]
fn test_n_at_edges() {
    let k = 21;
    let l = 8;
    let dna = random_dna(60, 77);

    // N's at start
    let mut seq = b"NNN".to_vec();
    seq.extend_from_slice(&dna);
    let iter = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new_with_n(
        &seq, k, l);
    let sks: Vec<Superkmer> = iter.collect();
    assert!(sks.iter().all(|sk| sk.start >= 3), "no superkmer should start before the N's");

    // N's at end
    let mut seq = dna.clone();
    seq.extend_from_slice(b"NNN");
    let iter = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new_with_n(
        &seq, k, l);
    let sks: Vec<Superkmer> = iter.collect();
    assert!(sks.iter().all(|sk| (sk.start + sk.size as usize) <= dna.len()),
        "no superkmer should extend into trailing N's");

    // Fragment too short to produce k-mers
    let mut seq = b"ACGTNNNNN".to_vec();
    seq.extend_from_slice(&dna);
    let iter = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new_with_n(
        &seq, k, l);
    let sks: Vec<Superkmer> = iter.collect();
    // "ACGT" is too short for k=21, should only get superkmers from the right fragment
    assert!(sks.iter().all(|sk| sk.start >= 9));
}

// ---- iteratorsyncmers2 naive correctness tests ----

mod syncmers2_correctness {
    use super::{encode_lmer, random_dna, read_fasta_seq, canonical_lmer};
    use rust_superkmers::Superkmer;

    const S: usize = 2; // syncmer s parameter

    /// Naive: is a given l-mer (as bytes) a closed syncmer?
    /// Closed syncmer: minimum s-mer (by lexicographic byte comparison, matching
    /// syncmers.rs `find_syncmers_pos`) is at position 0 or l-s.
    fn is_syncmer_naive(lmer: &[u8], l: usize, s: usize) -> bool {
        let num_smers = l - s + 1;

        // Find the position of the lexicographically smallest s-mer
        // (matching .min_by(|(_, a), (_, b)| a.cmp(b)) in syncmers.rs)
        let min_pos = (0..num_smers)
            .min_by(|&i, &j| lmer[i..i + s].cmp(&lmer[j..j + s]))
            .unwrap();

        // Closed syncmer: min s-mer at position 0 or l-s
        min_pos == 0 || min_pos == l - s
    }

    /// Build the syncmer score table naively for all 4^l l-mers.
    /// Score 0 = syncmer (good), 1 = non-syncmer.
    fn build_score_table(l: usize) -> Vec<usize> {
        let num_lmers = 1 << (2 * l);
        let mut scores = vec![0usize; num_lmers];
        for val in 0..num_lmers {
            // Decode to bytes
            let mut lmer = vec![0u8; l];
            for i in 0..l {
                lmer[i] = match (val >> (2 * (l - 1 - i))) & 3 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    3 => b'T',
                    _ => unreachable!(),
                };
            }
            scores[val] = if is_syncmer_naive(&lmer, l, S) { 0 } else { 1 };
        }
        scores
    }

    /// Naive MSP following the exact algorithm from debruijn's Scanner::scan():
    /// - Initialize: find best minimizer in first window (lowest score, rightmost on tie)
    /// - Slide: only change minimizer when it falls off the left edge (rescan full window)
    ///   or a strictly better (lower score) l-mer enters from the right
    /// - Group consecutive k-mers sharing the same minimizer into superkmers
    fn naive_superkmers(seq: &[u8], k: usize, l: usize) -> Vec<(usize, usize, usize, usize)> {
        // Returns: (start, size, mpos_relative, mint)
        let scores = build_score_table(l);
        let num_kmers = seq.len() - k + 1;
        let max_lmer_offset = k - l;

        // Find best minimizer in a window [start, stop] (lowest score, rightmost on tie)
        let find_min = |start: usize, stop: usize| -> (usize, usize) {
            // Returns (score, position)
            let mut best_score = usize::MAX;
            let mut best_pos = start;
            for pos in start..=stop {
                let lmer_val = encode_lmer(seq, pos, l);
                let score = scores[lmer_val];
                // Lower score wins; on tie, higher pos wins (rightmost)
                if score < best_score || (score == best_score && pos >= best_pos) {
                    best_score = score;
                    best_pos = pos;
                }
            }
            (best_score, best_pos)
        };

        // Initialize with first k-mer window
        let (mut min_score, mut min_pos) = find_min(0, max_lmer_offset);

        // Track minimizer changes (kmer_index, minimizer_position)
        let mut min_positions: Vec<(usize, usize)> = vec![(0, min_pos)];

        for i in 1..num_kmers {
            let end_lmer_pos = i + max_lmer_offset;
            let end_val = encode_lmer(seq, end_lmer_pos, l);
            let end_score = scores[end_val];

            if i > min_pos {
                // Current minimizer fell off left edge → rescan full window
                let (new_score, new_pos) = find_min(i, i + max_lmer_offset);
                min_score = new_score;
                min_pos = new_pos;
                min_positions.push((i, min_pos));
            } else if end_score < min_score {
                // Strictly better l-mer entered from the right
                min_score = end_score;
                min_pos = end_lmer_pos;
                min_positions.push((i, min_pos));
            }
            // Otherwise: keep current minimizer (sticky)
        }

        // Convert minimizer change points into superkmers
        let mut result = Vec::new();
        for p in 0..min_positions.len() {
            let (start_pos, mpos_abs) = min_positions[p];
            let end = if p + 1 < min_positions.len() {
                let (next_pos, _) = min_positions[p + 1];
                next_pos + k - 1 // last k-mer of this superkmer ends at next_pos-1+k-1
            } else {
                seq.len()
            };
            let size = end - start_pos;
            let mpos = mpos_abs - start_pos;
            let (mint, _) = canonical_lmer(seq, mpos_abs, l);
            result.push((start_pos, size, mpos, mint));
        }

        result
    }

    fn check_correctness(seq: &[u8], k: usize, l: usize) {
        let naive = naive_superkmers(seq, k, l);

        let iter = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(seq, k, l);
        let actual: Vec<Superkmer> = iter.collect();

        assert_eq!(naive.len(), actual.len(),
            "superkmer count mismatch: naive={} actual={} (seq_len={}, k={}, l={})",
            naive.len(), actual.len(), seq.len(), k, l);

        for (i, ((n_start, n_size, n_mpos, n_mint), a)) in naive.iter().zip(actual.iter()).enumerate() {
            assert_eq!(*n_start, a.start,
                "superkmer {} start: naive={} actual={}", i, n_start, a.start);
            assert_eq!(*n_size, a.size as usize,
                "superkmer {} size: naive={} actual={}", i, n_size, a.size);
            assert_eq!(*n_mpos, a.mpos as usize,
                "superkmer {} mpos: naive={} actual={}", i, n_mpos, a.mpos);
            assert_eq!(*n_mint, a.mint as usize,
                "superkmer {} mint: naive={:#x} actual={:#x}", i, n_mint, a.mint);
        }
    }

    #[test]
    fn test_syncmers2_naive_short() {
        let seq = b"AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCACACTGCACTGCACTGACTGCACTGCACTGCACTGCACTGCCTGC";
        check_correctness(seq, 21, 8);
        check_correctness(seq, 31, 8);
    }

    #[test]
    fn test_syncmers2_naive_random() {
        for k in [17, 21, 31] {
            for len in [k, k + 1, 50, 100, 300, 1000] {
                if len < k { continue; }
                for seed in 0..10 {
                    check_correctness(&random_dna(len, seed * 137 + 42), k, 8);
                }
            }
        }
    }

    #[test]
    fn test_syncmers2_naive_homopolymer() {
        for base in [b'A', b'C', b'G', b'T'] {
            let seq: Vec<u8> = std::iter::repeat(base).take(100).collect();
            check_correctness(&seq, 21, 8);
        }
    }

    #[test]
    fn test_syncmers2_naive_dinucleotide() {
        for pattern in [b"AC", b"GT", b"AG", b"CT"] {
            let seq: Vec<u8> = pattern.iter().cycle().take(200).cloned().collect();
            check_correctness(&seq, 21, 8);
            check_correctness(&seq, 31, 8);
        }
    }

    #[test]
    fn test_syncmers2_naive_large() {
        check_correctness(&random_dna(100_000, 999), 31, 8);
    }

    #[test]
    fn test_syncmers2_naive_ecoli() {
        let seq = read_fasta_seq("tests/ecoli.genome.220.fa");
        check_correctness(&seq, 21, 8);
        check_correctness(&seq, 31, 8);
    }
}

// ---- simdmini tests ----

#[cfg(feature = "simd-mini")]
mod simdmini {
    use super::{random_dna, canonical_lmer};
    use rust_superkmers::utils::split_on_n;
    use simd_minimizers::packed_seq::AsciiSeq;
    use std::collections::HashSet;

    /// Compute canonical l-mer using string-based reverse complement (independent of 2-bit encoding).
    fn canonical_lmer_string(seq: &[u8], pos: usize, l: usize) -> (usize, bool) {
        let lmer = &seq[pos..pos + l];
        // Build RC by reversing and complementing each base
        let rc_bytes: Vec<u8> = lmer.iter().rev().map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b,
        }).collect();
        // Encode both as 2-bit integers
        let encode = |s: &[u8]| -> usize {
            let mut v = 0usize;
            for &b in s {
                let bits = match b {
                    b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3,
                    _ => 0,
                };
                v = (v << 2) | bits;
            }
            v
        };
        let fwd_upper: Vec<u8> = lmer.iter().map(|&b| b.to_ascii_uppercase()).collect();
        let fwd_val = encode(&fwd_upper);
        let rc_val = encode(&rc_bytes);
        (fwd_val.min(rc_val), rc_val < fwd_val)
    }

    /// Get syncmer positions using the same simd-minimizers function as the iterator.
    fn simd_syncmer_positions(seq: &[u8], l: usize) -> HashSet<usize> {
        let s = 2usize;
        let w = l - s + 1;
        let upper: Vec<u8> = seq.iter().map(|&b| b.to_ascii_uppercase()).collect();
        let mut pos_vec: Vec<u32> = Vec::new();
        simd_minimizers::canonical_closed_syncmers(s, w)
            .run(AsciiSeq(&upper), &mut pos_vec);
        pos_vec.into_iter().map(|p| p as usize).collect()
    }

    /// Full correctness check:
    /// 1. Tiling: total k-mers match, no gaps/overlaps
    /// 2. Minimizer position is a syncmer (per simd-minimizers)
    /// 3. mint matches canonical l-mer encoding at that position
    /// 4. rc flag matches canonical orientation
    /// 5. Minimizer is within every k-mer's window bounds
    /// 6. mpos within bounds
    fn check_correctness(seq: &[u8], k: usize, l: usize) {
        let upper: Vec<u8> = seq.iter().map(|&b| b.to_ascii_uppercase()).collect();
        let iter = rust_superkmers::iteratorsimdmini::SuperkmersIterator::new(seq, k, l);
        let sks: Vec<_> = iter.collect();

        // 1. Tiling
        let total_kmers: usize = sks.iter().map(|s| s.size as usize - k + 1).sum();
        let expected = seq.len() - k + 1;
        assert_eq!(total_kmers, expected,
            "k-mer coverage: got {} expected {} (seq_len={})", total_kmers, expected, seq.len());

        for i in 1..sks.len() {
            let prev_end = sks[i-1].start + sks[i-1].size as usize - k + 1;
            assert_eq!(sks[i].start, prev_end,
                "gap/overlap at superkmer {}", i);
        }

        // Build syncmer set using the same SIMD function as the iterator
        let syncmer_set = simd_syncmer_positions(&upper, l);

        for (si, sk) in sks.iter().enumerate() {
            let mpos_abs = sk.start + sk.mpos as usize;

            // 2. Minimizer is a syncmer
            assert!(syncmer_set.contains(&mpos_abs),
                "superkmer {} (start={}): minimizer at pos {} is NOT a syncmer",
                si, sk.start, mpos_abs);

            // 3. mint matches canonical l-mer
            let (canon, is_rc) = canonical_lmer(&upper, mpos_abs, l);
            assert_eq!(sk.mint as usize, canon,
                "superkmer {} (start={}): mint={} but canonical at pos {} = {}",
                si, sk.start, sk.mint, mpos_abs, canon);

            // 4. rc flag
            assert_eq!(sk.mint_is_rc, is_rc,
                "superkmer {} (start={}): rc={} but expected {}", si, sk.start, sk.mint_is_rc, is_rc);

            // 5. Minimizer is within every k-mer's window
            let num_kmers = sk.size as usize - k + 1;
            for j in 0..num_kmers {
                let kmer_start = sk.start + j;
                let window_end = kmer_start + k - l;
                assert!(mpos_abs >= kmer_start && mpos_abs <= window_end,
                    "superkmer {} kmer {} (pos={}): minimizer at {} outside window [{}, {}]",
                    si, j, kmer_start, mpos_abs, kmer_start, window_end);
            }

            // 6. Boundary: at superkmer transitions, the previous minimizer
            //    must have fallen off the left edge of the new first k-mer's window
            if si > 0 {
                let prev_mpos_abs = sks[si-1].start + sks[si-1].mpos as usize;
                let first_kmer_start = sk.start;
                assert!(prev_mpos_abs < first_kmer_start,
                    "superkmer {}: previous minimizer at {} should have fallen off \
                     (first kmer starts at {})",
                    si, prev_mpos_abs, first_kmer_start);
            }

            // 7. mpos within bounds
            assert!((sk.mpos as usize) + l <= sk.size as usize,
                "superkmer {} mpos out of bounds", si);
        }
    }

    #[test]
    fn test_canonical_lmer_index_correctness() {
        // Verify the 2-bit canonical_lmer matches the string-based implementation
        // for every l-mer position across diverse sequences.
        let seqs: Vec<Vec<u8>> = vec![
            b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec(),
            b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_vec(),
            b"ATATATATATATATATATATATATATATATAT".to_vec(),
            b"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG".to_vec(),
            b"AACTTGGCCAATGGTTCCAAGGTTAACCGGTT".to_vec(),
            random_dna(1000, 42),
            random_dna(1000, 123),
            random_dna(1000, 999),
        ];

        for l in [3, 5, 7, 9, 11] {
            for seq in &seqs {
                if seq.len() < l { continue; }
                for pos in 0..=seq.len() - l {
                    let (val_2bit, rc_2bit) = canonical_lmer(seq, pos, l);
                    let (val_str, rc_str) = canonical_lmer_string(seq, pos, l);
                    assert_eq!(val_2bit, val_str,
                        "mint mismatch at pos={} l={} lmer={:?}",
                        pos, l, std::str::from_utf8(&seq[pos..pos+l]).unwrap_or("?"));
                    assert_eq!(rc_2bit, rc_str,
                        "rc mismatch at pos={} l={} lmer={:?}",
                        pos, l, std::str::from_utf8(&seq[pos..pos+l]).unwrap_or("?"));
                }
            }
        }
    }

    #[test]
    fn test_canonical_lmer_known_values() {
        // Hand-verified cases
        // ACG -> fwd=0b000110=6, rc=CGT=0b_01_10_11=27 -> canonical=6, rc=false
        assert_eq!(canonical_lmer(b"ACG", 0, 3), (6, false));
        assert_eq!(canonical_lmer_string(b"ACG", 0, 3), (6, false));

        // AAA -> fwd=0, rc=TTT=0b111111=63 -> canonical=0, rc=false
        assert_eq!(canonical_lmer(b"AAA", 0, 3), (0, false));

        // TTT -> fwd=0b111111=63, rc=AAA=0 -> canonical=0, rc=true
        assert_eq!(canonical_lmer(b"TTT", 0, 3), (0, true));

        // Palindromic: ACGT -> fwd=0b00011011=27, rc=ACGT=27 -> canonical=27, rc=false (rc==fwd, not rc<fwd)
        assert_eq!(canonical_lmer(b"ACGT", 0, 4), (27, false));
        assert_eq!(canonical_lmer_string(b"ACGT", 0, 4), (27, false));
    }

    #[test]
    fn test_simdmini_mint_rc_independent_verification() {
        // Run simdmini iterator and verify every superkmer's mint/rc
        // using the string-based canonical computation (fully independent of 2-bit encoding).
        for seed in [42u64, 77, 123, 456, 789] {
            let seq = random_dna(5000, seed);
            let k = 31;
            let l = 9;
            let upper: Vec<u8> = seq.iter().map(|&b| b.to_ascii_uppercase()).collect();
            let iter = rust_superkmers::iteratorsimdmini::SuperkmersIterator::new(&seq, k, l);
            for (si, sk) in iter.enumerate() {
                let mpos_abs = sk.start + sk.mpos as usize;
                let (expected_mint, expected_rc) = canonical_lmer_string(&upper, mpos_abs, l);
                assert_eq!(sk.mint as usize, expected_mint,
                    "seed={} superkmer {}: mint mismatch at pos {}", seed, si, mpos_abs);
                assert_eq!(sk.mint_is_rc, expected_rc,
                    "seed={} superkmer {}: rc mismatch at pos {}", seed, si, mpos_abs);
            }
        }
    }

    #[test]
    fn test_simdmini_correctness_short() {
        let seq = b"AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCACACTGCACTGCACTGACTGCACTGCACTGCACTGCACTGCCTGCAACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCACACTGCACTGCACTGACTGCACTGCACTGCACTGCACTGCCTGC";
        check_correctness(seq, 31, 9);
    }

    #[test]
    fn test_simdmini_correctness_random_10k() {
        check_correctness(&random_dna(10000, 42), 31, 9);
    }

    #[test]
    fn test_simdmini_correctness_homopolymer() {
        let mut seq = Vec::new();
        seq.extend_from_slice(b"ACGTACGTACGTACGT");
        for _ in 0..1000 { seq.push(b'A'); }
        seq.extend_from_slice(b"ACGTACGTACGTACGT");
        check_correctness(&seq, 31, 9);
    }

    #[test]
    fn test_simdmini_correctness_random_100k() {
        check_correctness(&random_dna(100_000, 99), 31, 9);
    }

    #[test]
    fn test_simdmini_correctness_with_n() {
        let mut seq = random_dna(500, 77);
        seq[200] = b'N';
        seq[201] = b'N';
        let iter = rust_superkmers::iteratorsimdmini::SuperkmersIterator::new_with_n(&seq, 31, 9);
        let sks: Vec<_> = iter.collect();
        let total_kmers: usize = sks.iter().map(|s| s.size as usize - 31 + 1).sum();
        let fragments = split_on_n(&seq, 31);
        let expected: usize = fragments.iter().filter(|(_, f)| f.len() >= 31).map(|(_, f)| f.len() - 30).sum();
        assert_eq!(total_kmers, expected, "N-split coverage mismatch");
    }
}

// ---- canonical vs non-canonical tests ----

mod canonical_toggle {
    use super::{random_dna, encode_lmer, rc_lmer, canonical_lmer};
    use rust_superkmers::Superkmer;
    use debruijn::dna_string::DnaString;

    /// Verify canonical vs non-canonical for iteratorsyncmers2:
    /// same start/size/mpos, canonical mint = min(fwd, rc), non-canonical rc=false.
    fn check_syncmers2(seq: &[u8], k: usize, l: usize) {
        let iter_c = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(seq, k, l);
        let canonical: Vec<Superkmer> = iter_c.collect();

        let iter_nc = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::non_canonical(seq, k, l);
        let non_canonical: Vec<Superkmer> = iter_nc.collect();

        assert_eq!(canonical.len(), non_canonical.len(), "superkmer count differs");

        for (i, (c, nc)) in canonical.iter().zip(non_canonical.iter()).enumerate() {
            // Same superkmer boundaries
            assert_eq!(c.start, nc.start, "start differs at {}", i);
            assert_eq!(c.size, nc.size, "size differs at {}", i);
            assert_eq!(c.mpos, nc.mpos, "mpos differs at {}", i);

            // Non-canonical must have rc=false
            assert!(!nc.mint_is_rc, "non-canonical rc should be false at {}", i);

            // Non-canonical mint is the forward-strand l-mer
            let fwd = encode_lmer(seq, c.start + c.mpos as usize, l);
            assert_eq!(nc.mint as usize, fwd,
                "non-canonical mint should be forward l-mer at {}", i);

            // Canonical mint is min(fwd, rc)
            let rc = rc_lmer(fwd, l);
            assert_eq!(c.mint as usize, fwd.min(rc),
                "canonical mint should be min(fwd, rc) at {}", i);

            // rc flag matches
            assert_eq!(c.mint_is_rc, rc < fwd,
                "canonical rc flag wrong at {}", i);
        }
    }

    #[test]
    fn test_syncmers2_canonical_toggle_random() {
        for k in [21, 31] {
            for len in [50, 150, 1000] {
                for seed in 0..5 {
                    check_syncmers2(&random_dna(len, seed * 137 + 42), k, 8);
                }
            }
        }
    }

    #[test]
    fn test_syncmers2_canonical_toggle_l9() {
        for seed in 0..5 {
            check_syncmers2(&random_dna(200, seed + 100), 31, 9);
        }
    }

    #[test]
    fn test_syncmers2_canonical_toggle_homopolymer() {
        for base in [b'A', b'C', b'G', b'T'] {
            let seq: Vec<u8> = std::iter::repeat(base).take(100).collect();
            check_syncmers2(&seq, 21, 8);
        }
    }

    #[test]
    fn test_syncmers2_canonical_toggle_with_n() {
        let mut seq = random_dna(300, 42);
        seq[100] = b'N';
        seq[200] = b'N';
        seq[201] = b'N';

        let iter_c = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new_with_n(&seq, 31, 8);
        let canonical: Vec<Superkmer> = iter_c.collect();

        let iter_nc = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::non_canonical_with_n(&seq, 31, 8);
        let non_canonical: Vec<Superkmer> = iter_nc.collect();

        assert_eq!(canonical.len(), non_canonical.len());
        for (c, nc) in canonical.iter().zip(non_canonical.iter()) {
            assert_eq!(c.start, nc.start);
            assert_eq!(c.size, nc.size);
            assert_eq!(c.mpos, nc.mpos);
            assert!(!nc.mint_is_rc);
        }
    }

    /// Verify canonical vs non-canonical for iteratorsyncmersmsp.
    fn check_syncmersmsp(seq: &[u8], k: usize, l: usize) {
        let dnastring = DnaString::from_acgt_bytes(seq).to_bytes();

        let iter_c = rust_superkmers::iteratorsyncmersmsp::SuperkmersIterator::new(&dnastring, k, l);
        let canonical: Vec<Superkmer> = iter_c.collect();

        let iter_nc = rust_superkmers::iteratorsyncmersmsp::SuperkmersIterator::non_canonical(&dnastring, k, l);
        let non_canonical: Vec<Superkmer> = iter_nc.collect();

        assert_eq!(canonical.len(), non_canonical.len(), "superkmer count differs");

        for (i, (c, nc)) in canonical.iter().zip(non_canonical.iter()).enumerate() {
            assert_eq!(c.start, nc.start, "start differs at {}", i);
            assert_eq!(c.size, nc.size, "size differs at {}", i);
            assert_eq!(c.mpos, nc.mpos, "mpos differs at {}", i);
            assert!(!nc.mint_is_rc, "non-canonical rc should be false at {}", i);

            let fwd = encode_lmer(seq, c.start + c.mpos as usize, l);
            assert_eq!(nc.mint as usize, fwd,
                "non-canonical mint should be forward l-mer at {}", i);

            let rc = rc_lmer(fwd, l);
            assert_eq!(c.mint as usize, fwd.min(rc),
                "canonical mint should be min(fwd, rc) at {}", i);
            assert_eq!(c.mint_is_rc, rc < fwd,
                "canonical rc flag wrong at {}", i);
        }
    }

    #[test]
    fn test_syncmersmsp_canonical_toggle_random() {
        for k in [21, 31] {
            for seed in 0..5 {
                check_syncmersmsp(&random_dna(200, seed * 73 + 11), k, 8);
            }
        }
    }

    #[cfg(feature = "simd-mini")]
    mod simdmini_toggle {
        use super::*;

        fn check_simdmini(seq: &[u8], k: usize, l: usize) {
            let iter_c = rust_superkmers::iteratorsimdmini::SuperkmersIterator::new(seq, k, l);
            let canonical: Vec<Superkmer> = iter_c.collect();

            let iter_nc = rust_superkmers::iteratorsimdmini::SuperkmersIterator::non_canonical(seq, k, l);
            let non_canonical: Vec<Superkmer> = iter_nc.collect();

            assert_eq!(canonical.len(), non_canonical.len(), "superkmer count differs");

            for (i, (c, nc)) in canonical.iter().zip(non_canonical.iter()).enumerate() {
                assert_eq!(c.start, nc.start, "start differs at {}", i);
                assert_eq!(c.size, nc.size, "size differs at {}", i);
                assert_eq!(c.mpos, nc.mpos, "mpos differs at {}", i);
                assert!(!nc.mint_is_rc, "non-canonical rc should be false at {}", i);

                let (can_val, _) = canonical_lmer(seq, c.start + c.mpos as usize, l);
                assert_eq!(c.mint as usize, can_val,
                    "canonical mint wrong at {}", i);

                let fwd = encode_lmer(seq, c.start + c.mpos as usize, l);
                assert_eq!(nc.mint as usize, fwd,
                    "non-canonical mint should be forward l-mer at {}", i);
            }
        }

        #[test]
        fn test_simdmini_canonical_toggle_random() {
            for seed in 0..10 {
                check_simdmini(&random_dna(300, seed * 53 + 7), 31, 9);
            }
        }

        #[test]
        fn test_simdmini_canonical_toggle_with_n() {
            let mut seq = random_dna(300, 42);
            seq[100] = b'N';
            seq[200] = b'N';

            let iter_c = rust_superkmers::iteratorsimdmini::SuperkmersIterator::new_with_n(&seq, 31, 9);
            let canonical: Vec<Superkmer> = iter_c.collect();

            let iter_nc = rust_superkmers::iteratorsimdmini::SuperkmersIterator::non_canonical_with_n(&seq, 31, 9);
            let non_canonical: Vec<Superkmer> = iter_nc.collect();

            assert_eq!(canonical.len(), non_canonical.len());
            for (c, nc) in canonical.iter().zip(non_canonical.iter()) {
                assert_eq!(c.start, nc.start);
                assert_eq!(c.size, nc.size);
                assert_eq!(c.mpos, nc.mpos);
                assert!(!nc.mint_is_rc);
            }
        }

        #[test]
        fn test_simdmini_canonical_toggle_homopolymer() {
            for base in [b'A', b'C', b'G', b'T'] {
                let seq: Vec<u8> = std::iter::repeat(base).take(100).collect();
                check_simdmini(&seq, 31, 9);
            }
        }
    }

    /// Encode a k-mer at position `pos` in `seq` as a canonical 2-bit integer.
    fn canonical_kmer(seq: &[u8], pos: usize, k: usize) -> u128 {
        let mut fwd = 0u128;
        let mut rc = 0u128;
        for i in 0..k {
            let base = match seq[pos + i] {
                b'A' | b'a' => 0u128,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 0,
            };
            fwd = (fwd << 2) | base;
            rc |= (3 - base) << (2 * i);
        }
        fwd.min(rc)
    }

    /// Test context-independence: the same canonical k-mer must always map to the
    /// same minimizer bucket, regardless of surrounding context.
    ///
    /// Strategy: take a long sequence, extract overlapping "reads" from different
    /// positions (and also their reverse complements), run the superkmer iterator
    /// on each read, and check that every canonical k-mer is assigned to the same mint.
    fn check_context_independence<F>(name: &str, seq: &[u8], k: usize, l: usize, make_iter: F)
    where
        F: Fn(&[u8]) -> Vec<Superkmer>,
    {
        use std::collections::HashMap;

        let rc_seq = |s: &[u8]| -> Vec<u8> {
            s.iter().rev().map(|&b| match b {
                b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C',
                b'a' => b't', b't' => b'a', b'c' => b'g', b'g' => b'c',
                _ => b,
            }).collect()
        };

        // Map: canonical k-mer -> (mint, first_read_description)
        let mut kmer_to_mint: HashMap<u128, (u32, String)> = HashMap::new();
        let mut violations = 0usize;

        let mut check_read = |read: &[u8], desc: &str| {
            let superkmers = make_iter(read);
            for sk in &superkmers {
                let sk_start = sk.start;
                let sk_end = sk.start + sk.size as usize;
                for pos in sk_start..=(sk_end - k) {
                    let canon = canonical_kmer(read, pos, k);
                    match kmer_to_mint.get(&canon) {
                        Some((prev_mint, prev_desc)) => {
                            if *prev_mint != sk.mint {
                                if violations < 10 {
                                    eprintln!("  VIOLATION in {}: k-mer at {} got mint={} but previously got mint={} from {}",
                                        desc, pos, sk.mint, prev_mint, prev_desc);
                                }
                                violations += 1;
                            }
                        }
                        None => {
                            kmer_to_mint.insert(canon, (sk.mint, desc.to_string()));
                        }
                    }
                }
            }
        };

        // Overlapping reads of various sizes from different offsets
        for read_len in [100, 150, 200, 500] {
            if read_len > seq.len() { continue; }
            let step = read_len / 3; // heavy overlap
            let mut offset = 0;
            while offset + read_len <= seq.len() {
                let read = &seq[offset..offset + read_len];
                check_read(read, &format!("fwd[{}..{}]", offset, offset + read_len));

                // Also check reverse complement of the same read
                let rc = rc_seq(read);
                check_read(&rc, &format!("rc[{}..{}]", offset, offset + read_len));

                offset += step;
            }
        }

        assert_eq!(violations, 0,
            "{}: {} k-mers assigned to different buckets depending on context ({} distinct k-mers tested)",
            name, violations, kmer_to_mint.len());
    }

    #[test]
    fn test_context_independence_syncmers2_mspxor() {
        for &len in &[1000, 10000] {
            for seed in 0..3 {
                let seq = random_dna(len, seed * 137 + 42);
                check_context_independence(
                    "syncmers2:mspxor", &seq, 31, 8,
                    |s| rust_superkmers::iteratorsyncmers2::SuperkmersIterator::mspxor(s, 31, 8).collect(),
                );
            }
        }
    }

    #[test]
    fn test_context_independence_syncmers2_msp() {
        for seed in 0..3 {
            let seq = random_dna(5000, seed * 73 + 11);
            check_context_independence(
                "syncmers2:msp", &seq, 31, 8,
                |s| rust_superkmers::iteratorsyncmers2::SuperkmersIterator::msp(s, 31, 8).collect(),
            );
        }
    }

    #[cfg(feature = "simd-mini")]
    #[test]
    fn test_context_independence_simdmini_mspxor() {
        for &len in &[1000, 10000] {
            for seed in 0..3 {
                let seq = random_dna(len, seed * 137 + 42);
                check_context_independence(
                    "simdmini:mspxor", &seq, 31, 9,
                    |s| rust_superkmers::iteratorsimdmini::SuperkmersIterator::mspxor(s, 31, 9).collect(),
                );
            }
        }
    }

    #[cfg(feature = "simd-mini")]
    #[test]
    fn test_context_independence_simdmini_msp() {
        for seed in 0..3 {
            let seq = random_dna(5000, seed * 73 + 11);
            check_context_independence(
                "simdmini:msp", &seq, 31, 9,
                |s| rust_superkmers::iteratorsimdmini::SuperkmersIterator::msp(s, 31, 9).collect(),
            );
        }
    }

    /// Verify that sticky mode is NOT context-independent (sanity check).
    /// This test expects violations — if sticky were context-independent,
    /// we wouldn't need the Msp/MspXor modes.
    #[test]
    fn test_sticky_is_context_dependent() {
        let seq = random_dna(5000, 42);
        let mut kmer_to_mint: std::collections::HashMap<u128, u32> = std::collections::HashMap::new();
        let mut violations = 0;

        let rc_seq = |s: &[u8]| -> Vec<u8> {
            s.iter().rev().map(|&b| match b {
                b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', _ => b,
            }).collect()
        };

        for read_len in [150, 200] {
            let mut offset = 0;
            while offset + read_len <= seq.len() {
                let read = &seq[offset..offset + read_len];
                let superkmers: Vec<Superkmer> = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(read, 31, 8).collect();
                for sk in &superkmers {
                    for pos in sk.start..=(sk.start + sk.size as usize - 31) {
                        let canon = canonical_kmer(read, pos, 31);
                        if let Some(&prev_mint) = kmer_to_mint.get(&canon) {
                            if prev_mint != sk.mint { violations += 1; }
                        } else {
                            kmer_to_mint.insert(canon, sk.mint);
                        }
                    }
                }
                let rc = rc_seq(read);
                let superkmers: Vec<Superkmer> = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(&rc, 31, 8).collect();
                for sk in &superkmers {
                    for pos in sk.start..=(sk.start + sk.size as usize - 31) {
                        let canon = canonical_kmer(&rc, pos, 31);
                        if let Some(&prev_mint) = kmer_to_mint.get(&canon) {
                            if prev_mint != sk.mint { violations += 1; }
                        } else {
                            kmer_to_mint.insert(canon, sk.mint);
                        }
                    }
                }
                offset += read_len / 3;
            }
        }
        assert!(violations > 0,
            "Expected sticky mode to have context-dependent violations, but found none");
    }

    /// Test palindromic l-mer: ACGT repeated (palindromic 4-mer). For l=8,
    /// ACGTACGT has rc=ACGTACGT — a self-palindrome. Both modes should agree on mint.
    #[test]
    fn test_palindromic_minimizer() {
        // ACGTACGT is palindromic (rc of ACGT is ACGT, so ACGTACGT rc = ACGTACGT)
        // Build a sequence where this is likely to be a minimizer
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let k = 21;
        let l = 8;

        let iter_c = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(seq, k, l);
        let canonical: Vec<Superkmer> = iter_c.collect();

        let iter_nc = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::non_canonical(seq, k, l);
        let non_canonical: Vec<Superkmer> = iter_nc.collect();

        for (c, nc) in canonical.iter().zip(non_canonical.iter()) {
            let fwd = encode_lmer(seq, c.start + c.mpos as usize, l);
            let rc = rc_lmer(fwd, l);
            if fwd == rc {
                // Palindromic: both modes should give same mint, rc=false
                assert_eq!(c.mint, nc.mint, "palindromic l-mer: mint should match");
                assert!(!c.mint_is_rc, "palindromic l-mer: canonical rc should be false");
            }
        }
    }

}

// =============================================================================
// SIMD batch extractor tests
// =============================================================================

mod simd_batch {
    use super::*;

    fn check_tiling(superkmers: &[Superkmer], seq_len: usize, k: usize) {
        if seq_len < k { return; }
        assert!(!superkmers.is_empty(), "should have at least one superkmer for len={}", seq_len);
        assert_eq!(superkmers[0].start, 0);

        for i in 0..superkmers.len() {
            let sk = &superkmers[i];
            let sk_end = sk.start + sk.size as usize;
            assert!(sk.size as usize >= k, "superkmer {} too small: size={}", i, sk.size);

            if i < superkmers.len() - 1 {
                assert_eq!(superkmers[i + 1].start, sk_end - k + 1,
                    "gap/overlap between superkmers {} and {}", i, i + 1);
            } else {
                assert_eq!(sk_end, seq_len,
                    "last superkmer doesn't reach end: sk_end={}, seq_len={}", sk_end, seq_len);
            }
        }
    }

    #[test]
    fn test_simd_batch_tiling_150bp() {
        let k = 31;
        let l = 8;
        let mut ext = rust_superkmers::syncmers_simd_l8k40max::SimdBatchExtractor::new(k, l);

        for seed in 0..100 {
            let seq = random_dna(150, 42 + seed);
            let batch: [&[u8]; 8] = [&seq; 8];
            let results = unsafe { ext.process_batch(&batch) };

            for lane in 0..8 {
                check_tiling(&results[lane], seq.len(), k);
            }
            // All lanes get same input → same output
            for lane in 1..8 {
                assert_eq!(results[0].len(), results[lane].len(),
                    "seed={} lane {} count mismatch", seed, lane);
                for (i, (a, b)) in results[0].iter().zip(results[lane].iter()).enumerate() {
                    assert_eq!((a.start, a.size, a.mpos), (b.start, b.size, b.mpos),
                        "seed={} sk={} mismatch lane 0 vs {}", seed, i, lane);
                }
            }
        }
    }

    #[test]
    fn test_simd_batch_various_lengths() {
        let k = 31;
        let l = 8;
        let mut ext = rust_superkmers::syncmers_simd_l8k40max::SimdBatchExtractor::new(k, l);

        for &len in &[50, 100, 150, 200, 300, 500] {
            let seq = random_dna(len, 123 + len as u64);
            let batch: [&[u8]; 8] = [&seq; 8];
            let results = unsafe { ext.process_batch(&batch) };
            for lane in 0..8 {
                check_tiling(&results[lane], seq.len(), k);
            }
        }
    }

    #[test]
    fn test_simd_batch_different_reads() {
        let k = 31;
        let l = 8;
        let mut ext = rust_superkmers::syncmers_simd_l8k40max::SimdBatchExtractor::new(k, l);

        let seqs: Vec<Vec<u8>> = (0..8).map(|i| random_dna(150, 100 + i)).collect();
        let batch: [&[u8]; 8] = std::array::from_fn(|i| seqs[i].as_slice());
        let results = unsafe { ext.process_batch(&batch) };
        for lane in 0..8 {
            check_tiling(&results[lane], seqs[lane].len(), k);
        }
    }

    #[test]
    fn test_simd_batch_mixed_lengths() {
        let k = 31;
        let l = 8;
        let mut ext = rust_superkmers::syncmers_simd_l8k40max::SimdBatchExtractor::new(k, l);

        let lens = [150, 100, 200, 50, 150, 300, 75, 150];
        let seqs: Vec<Vec<u8>> = lens.iter().enumerate()
            .map(|(i, &len)| random_dna(len, 200 + i as u64)).collect();
        let batch: [&[u8]; 8] = std::array::from_fn(|i| seqs[i].as_slice());
        let results = unsafe { ext.process_batch(&batch) };
        for lane in 0..8 {
            check_tiling(&results[lane], seqs[lane].len(), k);
        }
    }
}
