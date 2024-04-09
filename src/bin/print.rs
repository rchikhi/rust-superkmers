use std::fs;
use rust_superkmers::{iterator1, iteratormsp};

fn main() {
    //let genome_file = "../../tests/ecoli.genome.76.fa";
    let genome_file = "tests/ecoli.genome.220.fa";
    let contents = fs::read_to_string(genome_file).expect("Failed to read test genome file").split("\n").collect::<Vec<&str>>()[1].to_string();
    let k = 31;
    let m = 8;

    let result = iteratormsp::SuperkmersIterator::new(contents.as_bytes(), k,m);
    println!("Iteratormsp:");
    for superkmer in result {
        println!("{:?}",superkmer);
    }

    let result = iterator1::SuperkmersIterator::new(contents.as_bytes(), k,m);
    println!("Iterator1:");
    for superkmer in result {
        println!("{:?}",superkmer);
    }
}

