use std::fs;
use rust_superkmers::{iteratorsyncmers2, iteratorsyncmersmsp};
use debruijn::dna_string::DnaString;

fn main() {
    let args: Vec<String> = std::env::args().collect();

    let genome_file = if args.len() > 1 {
        &args[1]
    } else {
        "tests/ecoli.genome.100k.fa"
    };

    let k = 21;
    let m = 8;

    let contents = fs::read_to_string(genome_file)
        .expect("Failed to read genome file")
        .split('\n')
        .collect::<Vec<&str>>()[1]
        .to_string();

    println!("File: {}, size: {} bytes", genome_file, contents.len());

    // Run using iteratorsyncmers2 (AVX2)
    let start = std::time::Instant::now();
    let (_, iter) = iteratorsyncmers2::SuperkmersIterator::new(contents.as_bytes(), k, m);
    let result2: Vec<_> = iter.collect();
    let elapsed2 = start.elapsed();

    println!("\niteratorsyncmers2 (AVX2) results:");
    println!("--------------------------------");
    println!("Total time: {:?}", elapsed2);
    println!("Superkmers formed: {}", result2.len());
    println!("Speed: {:.1} MiB/s", contents.len() as f64 / elapsed2.as_micros() as f64);

    // Run using iteratorsyncmersmsp (reference implementation)
    let dnastring = DnaString::from_dna_string(&contents).to_bytes();
    let start = std::time::Instant::now();
    let iter = iteratorsyncmersmsp::SuperkmersIterator::new(&dnastring, k, m);
    let result_msp: Vec<_> = iter.collect();
    let elapsed_msp = start.elapsed();

    println!("\niteratorsyncmersmsp (reference) results:");
    println!("----------------------------------------");
    println!("Total time: {:?}", elapsed_msp);
    println!("Superkmers formed: {}", result_msp.len());
    println!("Speed: {:.1} MiB/s", contents.len() as f64 / elapsed_msp.as_micros() as f64);

    // Compare performance
    println!("\nPerformance comparison:");
    println!("----------------------");
    println!("iteratorsyncmers2 is {:.1}x vs reference",
        elapsed2.as_nanos() as f64 / elapsed_msp.as_nanos() as f64);

    // Verify results
    if result2.len() != result_msp.len() {
        println!("\nWARNING: Result count mismatch!");
        println!("iteratorsyncmers2: {} superkmers", result2.len());
        println!("iteratorsyncmersmsp: {} superkmers", result_msp.len());
    } else {
        println!("\nResult count matches: {} superkmers", result2.len());

        let mut all_match = true;
        for (i, (sk1, sk2)) in result2.iter().zip(result_msp.iter()).enumerate() {
            if sk1.start != sk2.start || sk1.mint != sk2.mint || sk1.size != sk2.size || sk1.mpos != sk2.mpos {
                all_match = false;
                println!("Mismatch at position {}:", i);
                println!("  iteratorsyncmers2: {:?}", sk1);
                println!("  iteratorsyncmersmsp: {:?}", sk2);
                break;
            }
        }

        if all_match {
            println!("All superkmers match exactly between implementations");
        }
    }
}
