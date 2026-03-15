use std::collections::HashMap;
use std::env;
use std::fs;
use std::io::{BufRead, BufReader};

use rust_superkmers::iteratorsyncmers2;
use rust_superkmers::iteratorkmc2;
use rust_superkmers::iteratormsp;
use rust_superkmers::iteratoruhs;
#[cfg(feature = "simd-mini")]
use rust_superkmers::iteratorsimdmini;
#[cfg(feature = "simd-mini")]
use rust_superkmers::iteratorsimdmini_cminim;
#[cfg(feature = "multi-mini")]
use rust_superkmers::iteratormultiminimizers;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <genome.fa> [k] [l] [method]", args[0]);
        eprintln!("  k: kmer length (default 31)");
        eprintln!("  l: minimizer length (default 8)");
        eprintln!("  method: 'syncmer' (default), 'kmc2', 'msp', 'simdmini', or 'multimini[:N]' (N=nb_hash, default 2)");
        std::process::exit(1);
    }

    let fasta_path = &args[1];
    let k: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(31);
    let l_arg: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(0); // 0 = auto
    let method = args.get(4).map(|s| s.as_str()).unwrap_or("syncmer");

    let base_method = method.split(':').next().unwrap();
    if !["syncmer", "kmc2", "msp", "simdmini", "cminim", "multimini", "uhs"].contains(&base_method) {
        eprintln!("Unknown method: {}. Use 'syncmer', 'kmc2', 'msp', 'simdmini', 'cminim', 'multimini[:N]', or 'uhs'.", method);
        eprintln!("Append ':classical' or ':msp' or ':mspxor' for context-independent minimizers (e.g. 'syncmer:msp', 'uhs:mspxor').");
        std::process::exit(1);
    }

    let suffix = method.split(':').nth(1).unwrap_or("");
    let split_mode = match suffix {
        "classical" => rust_superkmers::SplitMode::Classical,
        "msp" => rust_superkmers::SplitMode::Msp,
        "mspxor" => rust_superkmers::SplitMode::MspXor,
        _ => rust_superkmers::SplitMode::Sticky,
    };

    // Parse nb_hash for multimini (e.g. "multimini:4")
    // If just "multimini" with no :N, run all variants (1, 2, 4)
    let multimini_nb_hashes: Vec<usize> = if base_method == "multimini" {
        match method.split(':').nth(1).and_then(|s| s.parse().ok()) {
            Some(n) => vec![n],
            None => vec![2, 4, 8],
        }
    } else { vec![0] };

    // simdmini/multimini require odd l or even k-l; other methods use even l (8).
    let l = if base_method == "simdmini" || base_method == "cminim" {
        if l_arg == 0 { 9 } else if l_arg % 2 == 0 { eprintln!("Note: {} requires odd l, using l={}", base_method, l_arg + 1); l_arg + 1 } else { l_arg }
    } else if base_method == "multimini" {
        if l_arg == 0 { 9 } else if (k - l_arg) % 2 != 0 { eprintln!("Note: multimini requires k-l even, using l={}", l_arg + 1); l_arg + 1 } else { l_arg }
    } else {
        if l_arg == 0 { 8 } else { l_arg }
    };

    // Read FASTA once
    let sequences = read_fasta(fasta_path);

    for &nb_hash in &multimini_nb_hashes {
        if base_method == "multimini" {
            eprintln!("Running k={}  l={}  method=multimini  nb_hash={}", k, l, nb_hash);
        } else {
            eprintln!("Running k={}  l={}  method={}", k, l, method);
        }

        let mut bucket_counts: HashMap<u32, u64> = HashMap::new();
        let mut total_kmers: u64 = 0;
        let mut total_superkmers: u64 = 0;

        for (i, seq) in sequences.iter().enumerate() {
            process_seq(seq, k, l, method, nb_hash, split_mode, &mut bucket_counts, &mut total_kmers, &mut total_superkmers);
            if (i + 1) % 10 == 0 {
                eprintln!("  processed {} sequences, {} superkmers, {} kmers so far", i + 1, total_superkmers, total_kmers);
            }
        }

        eprintln!("Done. {} sequences, {} superkmers, {} total kmers, {} distinct minimizers",
            sequences.len(), total_superkmers, total_kmers, bucket_counts.len());

        if base_method == "multimini" && multimini_nb_hashes.len() > 1 {
            println!("--- multimini nb_hash={} ---", nb_hash);
        }
        print_stats(&bucket_counts, total_kmers, total_superkmers, l);
        if base_method == "multimini" && multimini_nb_hashes.len() > 1 {
            println!();
        }
    }
}

fn read_fasta(path: &str) -> Vec<Vec<u8>> {
    eprintln!("Reading {}", path);
    let file = fs::File::open(path).expect("Failed to open FASTA file");
    let reader = BufReader::new(file);
    let mut sequences = Vec::new();
    let mut current_seq = Vec::<u8>::new();
    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if line.starts_with('>') {
            if !current_seq.is_empty() {
                sequences.push(std::mem::take(&mut current_seq));
            }
        } else {
            current_seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if !current_seq.is_empty() {
        sequences.push(current_seq);
    }
    eprintln!("Read {} sequences", sequences.len());
    sequences
}

fn process_seq(seq: &[u8], k: usize, l: usize, method: &str, nb_hash: usize, split_mode: rust_superkmers::SplitMode, bucket_counts: &mut HashMap<u32, u64>, total_kmers: &mut u64, total_superkmers: &mut u64) {
    let base_method = method.split(':').next().unwrap();
    match base_method {
        "syncmer" => {
            let iter = match split_mode {
                rust_superkmers::SplitMode::Classical => iteratorsyncmers2::SuperkmersIterator::classical_with_n(seq, k, l),
                rust_superkmers::SplitMode::Msp => iteratorsyncmers2::SuperkmersIterator::msp_with_n(seq, k, l),
                rust_superkmers::SplitMode::MspXor => iteratorsyncmers2::SuperkmersIterator::mspxor_with_n(seq, k, l),
                rust_superkmers::SplitMode::Sticky => iteratorsyncmers2::SuperkmersIterator::new_with_n(seq, k, l),
            };
            count_superkmers(iter, k, bucket_counts, total_kmers, total_superkmers);
        }
        "kmc2" => {
            let (_storage, iter) = iteratorkmc2::SuperkmersIterator::new_with_n(seq, k, l);
            count_superkmers(iter, k, bucket_counts, total_kmers, total_superkmers);
        }
        "msp" => {
            let superkmers = iteratormsp::superkmers_with_n(seq, k, l);
            count_superkmers(superkmers.into_iter(), k, bucket_counts, total_kmers, total_superkmers);
        }
        #[cfg(feature = "simd-mini")]
        "simdmini" => {
            let iter = match split_mode {
                rust_superkmers::SplitMode::Classical => iteratorsimdmini::SuperkmersIterator::classical_with_n(seq, k, l),
                rust_superkmers::SplitMode::Msp => iteratorsimdmini::SuperkmersIterator::msp_with_n(seq, k, l),
                rust_superkmers::SplitMode::MspXor => iteratorsimdmini::SuperkmersIterator::mspxor_with_n(seq, k, l),
                rust_superkmers::SplitMode::Sticky => iteratorsimdmini::SuperkmersIterator::new_with_n(seq, k, l),
            };
            count_superkmers(iter, k, bucket_counts, total_kmers, total_superkmers);
        }
        #[cfg(feature = "simd-mini")]
        "cminim" => {
            let iter = iteratorsimdmini_cminim::SuperkmersIterator::new_with_n(seq, k, l);
            count_superkmers(iter, k, bucket_counts, total_kmers, total_superkmers);
        }
        #[cfg(feature = "multi-mini")]
        "multimini" => {
            let iter = iteratormultiminimizers::SuperkmersIterator::new_with_n(seq, k, l, nb_hash);
            count_superkmers(iter, k, bucket_counts, total_kmers, total_superkmers);
        }
        "uhs" => {
            let iter = match split_mode {
                rust_superkmers::SplitMode::MspXor => iteratoruhs::SuperkmersIterator::mspxor_with_n(seq, k, l),
                _ => iteratoruhs::SuperkmersIterator::new_with_n(seq, k, l),
            };
            count_superkmers(iter, k, bucket_counts, total_kmers, total_superkmers);
        }
        _ => unreachable!(),
    }
}

fn count_superkmers(iter: impl Iterator<Item = rust_superkmers::Superkmer>, k: usize, bucket_counts: &mut HashMap<u32, u64>, total_kmers: &mut u64, total_superkmers: &mut u64) {
    for superkmer in iter {
        let num_kmers = superkmer.size as u64 - k as u64 + 1;
        *bucket_counts.entry(superkmer.mint).or_insert(0) += num_kmers;
        *total_kmers += num_kmers;
        *total_superkmers += 1;
    }
}

fn print_stats(bucket_counts: &HashMap<u32, u64>, total_kmers: u64, total_superkmers: u64, l: usize) {
    let mut sizes: Vec<u64> = bucket_counts.values().cloned().collect();
    sizes.sort_unstable();

    let n = sizes.len();
    if n == 0 {
        println!("No data.");
        return;
    }

    let sum: u64 = sizes.iter().sum();
    let mean = sum as f64 / n as f64;
    let median = sizes[n / 2];
    let p1 = sizes[n / 100];
    let p5 = sizes[n * 5 / 100];
    let p25 = sizes[n * 25 / 100];
    let p75 = sizes[n * 75 / 100];
    let p95 = sizes[n * 95 / 100];
    let p99 = sizes[n * 99 / 100];
    let max = *sizes.last().unwrap();
    let min = sizes[0];

    println!("=== Bucket size distribution (kmers per minimizer) ===");
    println!("Distinct minimizers: {}", n);
    println!("Total kmers:        {}", total_kmers);
    println!("Total superkmers:   {}", total_superkmers);
    println!();
    println!("Min:    {}", min);
    println!("P1:     {}", p1);
    println!("P5:     {}", p5);
    println!("P25:    {}", p25);
    println!("Median: {}", median);
    println!("Mean:   {:.1}", mean);
    println!("P75:    {}", p75);
    println!("P95:    {}", p95);
    println!("P99:    {}", p99);
    println!("Max:    {}", max);
    println!("Max/Mean: {:.1}x", max as f64 / mean);
    println!("Max/Median: {:.1}x", max as f64 / median as f64);
    println!();

    // Log2 histogram
    let mut log_bins: HashMap<u32, (u64, u64)> = HashMap::new();
    for &size in &sizes {
        let bin = if size == 0 { 0 } else { 64 - size.leading_zeros() };
        let entry = log_bins.entry(bin).or_insert((0, 0));
        entry.0 += 1;
        entry.1 += size;
    }
    let mut bin_keys: Vec<u32> = log_bins.keys().cloned().collect();
    bin_keys.sort_unstable();

    println!("=== Histogram (bucket size -> count of minimizers) ===");
    println!("{:<20} {:>15} {:>15} {:>10}", "Bucket size range", "# minimizers", "# kmers", "% kmers");
    for bin in bin_keys {
        let (count, kmer_count) = log_bins[&bin];
        let lo = if bin == 0 { 0 } else { 1u64 << (bin - 1) };
        let hi = (1u64 << bin) - 1;
        let pct = 100.0 * kmer_count as f64 / total_kmers as f64;
        println!("{:<20} {:>15} {:>15} {:>9.2}%", format!("[{}, {}]", lo, hi), count, kmer_count, pct);
    }

    // Top 20
    println!();
    println!("=== Top 20 largest buckets ===");
    let mut by_size: Vec<(u32, u64)> = bucket_counts.iter().map(|(&k, &v)| (k, v)).collect();
    by_size.sort_by(|a, b| b.1.cmp(&a.1));
    for (i, (mint, count)) in by_size.iter().take(20).enumerate() {
        let minimizer_str = decode_kmer(*mint as usize, l);
        println!("  {:>2}. {} (mint={:>10})  {} kmers", i + 1, minimizer_str, mint, count);
    }
}

fn decode_kmer(val: usize, l: usize) -> String {
    let mut s = String::with_capacity(l);
    for i in (0..l).rev() {
        let base = (val >> (2 * i)) & 3;
        s.push(match base {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => unreachable!(),
        });
    }
    s
}
