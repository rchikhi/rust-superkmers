use std::fs;
use rust_superkmers::{Superkmer, SuperkmerVerbose, naive, iterator1};
use std::collections::HashSet;
use debruijn::dna_string::DnaString;

#[test]
fn main() {

    //let genome_file = "tests/ecoli.genome.100k.fa";
    //let genome_file = "tests/ecoli.genome.76.fa";
    let genome_file = "tests/ecoli.genome.220.fa";
    let contents = fs::read_to_string(genome_file).expect("Failed to read test genome file").split("\n").collect::<Vec<&str>>()[1].to_string();

    // testing syncmers with l=10
    let dnastring = DnaString::from_dna_string(&contents).to_bytes();
    let iter = rust_superkmers::iteratorsyncmers::SuperkmersIterator::new(&dnastring, 31, 10);
    iter.into_iter().collect::<Vec<Superkmer>>(); 


    // a stresstest over some parameters
    for m in vec![4,5,7,11,14] {
        for k in vec![17,21,31,41] {

            let iter = iterator1::SuperkmersIterator::new(contents.as_bytes(), k, m);
            let iterator1_superkmers = iter.into_iter().collect::<Vec<Superkmer>>(); 
            let iterator1_superkmers_verbose: Vec<_> = iterator1_superkmers.into_iter().map(|superkmer| {
                rust_superkmers::utils::superkmer_to_verbose(superkmer, &contents, m)
            }).collect();

            let (iter, iter_verbose) = naive::extract_superkmers(contents.as_bytes(), k, m);
            let naive_superkmers = iter.into_iter().collect::<Vec<Superkmer>>();
            let naive_superkmers_verbose = iter_verbose.into_iter().collect::<Vec<SuperkmerVerbose>>();

            /*
            let set1: HashSet<_> = iterator1_superkmers.into_iter().collect();
            let set2: HashSet<_> = naive_superkmers.into_iter().collect();
            */
            let set1: HashSet<_> = iterator1_superkmers_verbose.into_iter().collect();
            let set2: HashSet<_> = naive_superkmers_verbose.into_iter().collect();

            println!("k={} m={} testfile={}",k,m, genome_file);
            if set1 == set2 {
                println!("Both sets of superkmers are equal.\n\n");
            } else {
                let only_in_set1: HashSet<_> = set1.difference(&set2).collect();
                let only_in_set2: HashSet<_> = set2.difference(&set1).collect();

                if !only_in_set1.is_empty() {
                    println!("Superkmers only in iterator1 set:");
                    for superkmer in only_in_set1 {
                        println!("{:?}", superkmer);
                    }
                }

                if !only_in_set2.is_empty() {
                    println!("Superkmers only in naive set:");
                    for superkmer in only_in_set2 {
                        println!("{:?}", superkmer);
                    }
                }
                assert!(false);
            }
        }
    }
}
