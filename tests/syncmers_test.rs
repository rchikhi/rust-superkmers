use std::fs;
use debruijn::dna_string::DnaString;
use rust_superkmers::Superkmer;

#[test]
fn test_iteratorsyncmers2_with_short_sequence() {
    let genome_file = "tests/ecoli.genome.76.fa";
    let contents = fs::read_to_string(genome_file)
        .expect("Failed to read test genome file")
        .split("\n")
        .collect::<Vec<&str>>()[1]
        .to_string();
    
    let k = 21;
    let l = 8;
    
    // Process using iteratorsyncmers2
    let iter = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(contents.as_bytes(), k, l);
    let syncmers_test = iter.1.into_iter().collect::<Vec<Superkmer>>();
    
    // Check if we got any results
    assert!(!syncmers_test.is_empty(), "iteratorsyncmers2 should produce superkmers with k={}, l={}", k, l);
    
    // Print some debug info about the superkmers
    println!("Found {} superkmers with iteratorsyncmers2", syncmers_test.len());
}

#[test]
fn test_iteratorsyncmersmsp_with_short_sequence() {
    let genome_file = "tests/ecoli.genome.76.fa";
    let contents = fs::read_to_string(genome_file)
        .expect("Failed to read test genome file")
        .split("\n")
        .collect::<Vec<&str>>()[1]
        .to_string();
    
    let k = 21;
    let l = 8; // MSP implementation only supports l=8, 10, or 12
    
    // Create DnaString for iteratorsyncmersmsp
    let dnastring = DnaString::from_dna_string(&contents).to_bytes();
    
    // Process using iteratorsyncmersmsp
    let iter = rust_superkmers::iteratorsyncmersmsp::SuperkmersIterator::new(&dnastring, k, l);
    let syncmers_msp = iter.into_iter().collect::<Vec<Superkmer>>();
    
    // Check if we got any results
    assert!(!syncmers_msp.is_empty(), "iteratorsyncmersmsp should produce superkmers with k={}, l={}", k, l);
    
    // Print some debug info about the superkmers
    println!("Found {} superkmers with iteratorsyncmersmsp", syncmers_msp.len());
}

#[test]
fn test_naivesyncmers() {
    let genome_file = "tests/ecoli.genome.76.fa";
    let contents = fs::read_to_string(genome_file)
        .expect("Failed to read test genome file")
        .split("\n")
        .collect::<Vec<&str>>()[1]
        .to_string();
    
    let k = 21;
    let l = 8;
    
    // Process using naivesyncmers
    let (superkmers, _) = rust_superkmers::naivesyncmers::extract_superkmers(contents.as_bytes(), k, l);
    
    // Check if we got any results
    assert!(!superkmers.is_empty(), "naivesyncmers should produce superkmers");
    
    // Print results
    println!("Found {} superkmers with naivesyncmers", superkmers.len());
}

#[test]
fn test_compare_implementations() {
    let genome_file = "tests/ecoli.genome.76.fa";
    let contents = fs::read_to_string(genome_file)
        .expect("Failed to read test genome file")
        .split("\n")
        .collect::<Vec<&str>>()[1]
        .to_string();
    
    let k = 21;
    let l = 8;
    
    // Process using all implementations
    
    // naivesyncmers
    let (naive_superkmers, _) = rust_superkmers::naivesyncmers::extract_superkmers(contents.as_bytes(), k, l);
    
    // iteratorsyncmers2
    let iter = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(contents.as_bytes(), k, l);
    let syncmers_test = iter.1.into_iter().collect::<Vec<Superkmer>>();
    
    // iteratorsyncmersmsp
    let dnastring = DnaString::from_dna_string(&contents).to_bytes();
    let iter = rust_superkmers::iteratorsyncmersmsp::SuperkmersIterator::new(&dnastring, k, l);
    let syncmers_msp = iter.into_iter().collect::<Vec<Superkmer>>();
    
    // Compare counts
    println!("Comparing implementations with k={}, l={}", k, l);
    println!("naivesyncmers: {} superkmers", naive_superkmers.len());
    println!("iteratorsyncmersmsp: {} superkmers", syncmers_msp.len());
    println!("iteratorsyncmers2: {} superkmers", syncmers_test.len());
    
    // Print the first few superkmers from each implementation to show the difference
    if !naive_superkmers.is_empty() && !syncmers_msp.is_empty() {
        println!("\nImplementation Differences:");
        println!("-------------------------");
        println!("1. Algorithm Overview:");
        println!("   - naivesyncmers: A direct implementation that computes minimizers for each k-mer individually");
        println!("   - iteratorsyncmersmsp: Uses the MSP (Minimum Substring Partitioning) algorithm from the debruijn library");
        
        println!("\n2. Superkmer Formation Strategy:");
        println!("   - naivesyncmers: Breaks a superkmer whenever the position of the minimizer changes in consecutive k-mers");
        println!("   - iteratorsyncmersmsp: Only breaks a superkmer when the actual minimizer itself changes");
        
        println!("\n3. Minimizer Selection:");
        println!("   - naivesyncmers: Computes the rightmost syncmer for each k-mer window independently");
        println!("   - iteratorsyncmersmsp: The MSP algorithm maintains a global view across the sequence,");
        println!("     selecting the rightmost syncmer overall until a lexicographically smaller one is found");
        
        println!("\n4. Sample Output Comparison:");
        println!("   First naive superkmer: start={}, size={}, mpos={}",
                 naive_superkmers[0].start, naive_superkmers[0].size, naive_superkmers[0].mpos);
        println!("   First MSP superkmer: start={}, size={}, mpos={}",
                 syncmers_msp[0].start, syncmers_msp[0].size, syncmers_msp[0].mpos);
                 
        println!("\n5. Position Tracking Example:");
        println!("   For example, in naivesyncmers, if two consecutive k-mers have the same minimizer");
        println!("   but at different positions (e.g., positions 10 and 9), a new superkmer is started.");
        println!("   In iteratorsyncmersmsp, these would be part of the same superkmer since the");
        println!("   minimizer itself hasn't changed, just its relative position.");
    }
    
    // Verify each implementation produces reasonable output
    assert!(!naive_superkmers.is_empty(), "naivesyncmers should produce superkmers");
    assert!(!syncmers_msp.is_empty(), "iteratorsyncmersmsp should produce superkmers");
    assert!(!syncmers_test.is_empty(), "iteratorsyncmers2 should produce superkmers");
}

#[test]
fn test_various_k_values() {
    let genome_file = "tests/ecoli.genome.220.fa";
    let contents = fs::read_to_string(genome_file)
        .expect("Failed to read test genome file")
        .split("\n")
        .collect::<Vec<&str>>()[1]
        .to_string();
    
    // Test different k values
    for k in [17, 21, 31, 41].iter() {
        let l = 8;
        
        // Create DnaString for iteratorsyncmersmsp
        let dnastring = DnaString::from_dna_string(&contents).to_bytes();
        
        // Process using iteratorsyncmersmsp
        let iter = rust_superkmers::iteratorsyncmersmsp::SuperkmersIterator::new(&dnastring, *k, l);
        let syncmers_msp = iter.into_iter().collect::<Vec<Superkmer>>();
        
        // Process using naivesyncmers
        let (naive_superkmers, _) = rust_superkmers::naivesyncmers::extract_superkmers(contents.as_bytes(), *k, l);
        
        println!("Testing k={}, l={}", k, l);
        assert!(!syncmers_msp.is_empty(), "iteratorsyncmersmsp should produce superkmers");
        assert!(!naive_superkmers.is_empty(), "naivesyncmers should produce superkmers");
        
        // Compare sizes of the outputs
        println!("iteratorsyncmersmsp produced {} superkmers", syncmers_msp.len());
        println!("naivesyncmers produced {} superkmers", naive_superkmers.len());
    }
}