use std::fs;
use std::collections::HashSet;
use std::hash::Hash;
use std::fmt::Debug;
use rust_superkmers::{Superkmer, naive, naive::superkmer_to_verbose, iterator1};
use debruijn::dna_string::DnaString;

fn compare_sets<T: Eq + Hash + Debug + std::cmp::Ord>(set1: HashSet<T>, set2: HashSet<T>) {
    if set1 == set2 {
        return;
    }
    let only_in_set1: Vec<_> = set1.difference(&set2).collect();
    let only_in_set2: Vec<_> = set2.difference(&set1).collect();

    panic!("Sets not equal: {}/{} only in set1, {}/{} only in set2",
        only_in_set1.len(), set1.len(), only_in_set2.len(), set2.len());
}

#[test]
fn main() {
    let genome_file = "tests/ecoli.genome.220.fa";
    let contents = fs::read_to_string(genome_file)
        .expect("Failed to read test genome file")
        .split("\n").collect::<Vec<&str>>()[1].to_string();

    // Test iteratorsyncmers2 vs iteratorsyncmersmsp (should be identical)
    let dnastring = DnaString::from_dna_string(&contents).to_bytes();
    let iter = rust_superkmers::iteratorsyncmersmsp::SuperkmersIterator::new(&dnastring, 31, 8);
    let syncmers_truth: Vec<Superkmer> = iter.collect();

    let iter = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(contents.as_bytes(), 31, 8);
    let syncmers_test: Vec<Superkmer> = iter.collect();

    let set1: HashSet<_> = syncmers_truth.into_iter().collect();
    let set2: HashSet<_> = syncmers_test.into_iter().collect();
    println!("syncmers k=31 l=8 testfile={}", genome_file);
    println!("iteratorsyncmersmsp: {} superkmers, iteratorsyncmers2: {} superkmers", set1.len(), set2.len());
    compare_sets(set1, set2);

    // Stress test iterator1 vs naive over various parameters
    for m in vec![4, 5, 7, 11, 14] {
        for k in vec![17, 21, 31, 41] {
            let iter = iterator1::SuperkmersIterator::new(contents.as_bytes(), k, m);
            let iterator1_superkmers_verbose: Vec<_> = iter.into_iter().map(|superkmer| {
                superkmer_to_verbose(superkmer, &contents, m)
            }).collect();

            let (_, iter_verbose) = naive::extract_superkmers(contents.as_bytes(), k, m);
            let naive_superkmers_verbose: Vec<_> = iter_verbose.into_iter().collect();

            let set1: HashSet<_> = iterator1_superkmers_verbose.into_iter().collect();
            let set2: HashSet<_> = naive_superkmers_verbose.into_iter().collect();
            println!("lexi k={} m={} testfile={}", k, m, genome_file);
            compare_sets(set1, set2);
        }
    }
}
