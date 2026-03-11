/// Cross-bucket test: does mspxor produce k-mer duplicates across buckets?
/// A k-mer is "cross-bucket" if it appears in superkmers with different mints
/// across different reads.
///
/// Uses tests/read50x_ref10K_e001.fa (5000 reads, ~100bp each, 50x coverage of a 10K reference, 1% error).
///
/// Run:  cargo +nightly test test_mspxor_crossbucket -- --nocapture

use std::collections::HashMap;

fn canonical(kmer: &[u8]) -> Vec<u8> {
    let rc: Vec<u8> = kmer.iter().rev().map(|&b| match b {
        b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', _ => b
    }).collect();
    if kmer <= &rc[..] { kmer.to_vec() } else { rc }
}

fn load_reads() -> Vec<Vec<u8>> {
    let content = std::fs::read_to_string("tests/read50x_ref10K_e001.fa")
        .expect("Need tests/read50x_ref10K_e001.fa");
    content.lines()
        .filter(|l| !l.starts_with('>'))
        .map(|l| l.trim().as_bytes().to_vec())
        .collect()
}

fn check_crossbucket<F>(name: &str, reads: &[Vec<u8>], k: usize, mut process: F)
where
    F: FnMut(&[u8]) -> Vec<(usize, u16, u32)>, // returns (start, size, mint) per superkmer
{
    let mut kmer_mints: HashMap<Vec<u8>, Vec<u32>> = HashMap::new();

    for read in reads {
        if read.len() < k { continue; }
        for (start, size, mint) in process(read) {
            let n_kmers = size as usize - k + 1;
            for i in 0..n_kmers {
                let kmer = &read[start + i..start + i + k];
                let can = canonical(kmer);
                let entry = kmer_mints.entry(can).or_default();
                if !entry.contains(&mint) {
                    entry.push(mint);
                }
            }
        }
    }

    let total = kmer_mints.len();
    let multi_bucket: Vec<_> = kmer_mints.iter()
        .filter(|(_, mints)| mints.len() > 1)
        .collect();

    println!("{}: {} distinct canonical k-mers, {} cross-bucket", name, total, multi_bucket.len());

    for (kmer, mints) in multi_bucket.iter().take(5) {
        println!("  {} -> mints {:?}", std::str::from_utf8(kmer).unwrap(), mints);
    }

    assert_eq!(multi_bucket.len(), 0,
        "{}: {} k-mers appear in multiple buckets (expected 0)",
        name, multi_bucket.len());
}

#[test]
fn test_mspxor_crossbucket_simdmini() {
    let reads = load_reads();
    println!("{} reads loaded", reads.len());

    let mut extractor = rust_superkmers::iteratorsimdmini::SuperkmerExtractor::mspxor(31, 9);
    check_crossbucket("simdmini:mspxor", &reads, 31, |read| {
        extractor.process(read).iter()
            .map(|sk| (sk.start, sk.size, sk.mint))
            .collect()
    });
}

#[test]
fn test_mspxor_crossbucket_syncmers2() {
    let reads = load_reads();
    println!("{} reads loaded", reads.len());

    let mut extractor = rust_superkmers::iteratorsyncmers2::SuperkmerExtractor::mspxor(31, 8);
    check_crossbucket("syncmers2:mspxor", &reads, 31, |read| {
        extractor.process(read).iter()
            .map(|sk| (sk.start, sk.size, sk.mint))
            .collect()
    });
}

#[test]
fn test_msp_crossbucket_simdmini() {
    let reads = load_reads();
    let mut extractor = rust_superkmers::iteratorsimdmini::SuperkmerExtractor::msp(31, 9);
    check_crossbucket("simdmini:msp", &reads, 31, |read| {
        extractor.process(read).iter()
            .map(|sk| (sk.start, sk.size, sk.mint))
            .collect()
    });
}

#[test]
fn test_msp_crossbucket_syncmers2() {
    let reads = load_reads();
    let mut extractor = rust_superkmers::iteratorsyncmers2::SuperkmerExtractor::msp(31, 8);
    check_crossbucket("syncmers2:msp", &reads, 31, |read| {
        extractor.process(read).iter()
            .map(|sk| (sk.start, sk.size, sk.mint))
            .collect()
    });
}
