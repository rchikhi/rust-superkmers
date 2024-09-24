use std::fs;
use rust_superkmers::{naive, iterator1, iteratormsp, iteratorsyncmersmsp, iteratorsyncmers2, naivesyncmers};
use debruijn::dna_string::DnaString;

fn main() {
    //let genome_file = "../../tests/ecoli.genome.76.fa";
    let genome_file = "tests/ecoli.genome.220.fa";
    let contents = fs::read_to_string(genome_file).expect("Failed to read test genome file").split("\n").collect::<Vec<&str>>()[1].to_string();
    let k = 31;
    let mut m = 8;
    let dnastring = DnaString::from_dna_string(&contents).to_bytes();

    println!("Lexicographic minimizers:");
    println!("-------------------------");
    // classical minimizers, msp superkmers construction
    let result = iteratormsp::SuperkmersIterator::new(&dnastring, k,m);
    println!("Iteratormsp:");
    for superkmer in result {
        println!("{:?}",superkmer);
    }

    println!("NtHash, window minimizers");
    println!("-------------------------");
    // nthash, naive superkmer construction 
    let (result, osef) = naive::extract_superkmers(contents.as_bytes(), k, m);
    println!("Naive:");
    for superkmer  in result {
        println!("{:?}",superkmer);
    }

    // nthash, slightly more sophisticated superkmers construction
    let result = iterator1::SuperkmersIterator::new(contents.as_bytes(), k,m);
    println!("Iterator1:");
    for superkmer in result {
        println!("{:?}",superkmer);
    }

    //m=20; // need to increase minimizer length because otherwise too many syncmers per kmer
    println!("Syncmers (m={}):", m);
    println!("-------------------------");
    // syncmers, naive superkmer construction 
    let (result, osef) = naivesyncmers::extract_superkmers(contents.as_bytes(), k, m);
    println!("Naivesyncmers:");
    for superkmer  in result {
        println!("{:?}",superkmer);
    }

    // syncmers, msp
    // doesn't support large m
    if m <= 12 {
        let result = iteratorsyncmersmsp::SuperkmersIterator::new(&dnastring, k,m);
        println!("Iteratorsyncmersmsp:");
        for superkmer in result {
            println!("{:?}",superkmer);
        }
    }


    let result = iteratorsyncmers2::SuperkmersIterator::new(contents.as_bytes(), k,m);
    println!("Iteratorsyncmers2:");
    for superkmer in result.1 {
        println!("{:?}",superkmer);
    }

}

