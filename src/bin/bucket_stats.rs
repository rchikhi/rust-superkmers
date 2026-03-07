use std::collections::HashMap;
use std::env;
use std::fs;
use std::io::{BufRead, BufReader};

use rust_superkmers::iteratorsyncmers2;
use rust_superkmers::iteratorkmc2;
use rust_superkmers::iteratormsp;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <genome.fa> [k] [l] [method]", args[0]);
        eprintln!("  k: kmer length (default 31)");
        eprintln!("  l: minimizer length (default 8)");
        eprintln!("  method: 'syncmer' (default), 'kmc2', or 'msp'");
        std::process::exit(1);
    }

    let fasta_path = &args[1];
    let k: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(31);
    let l: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(8);
    let method = args.get(4).map(|s| s.as_str()).unwrap_or("syncmer");

    if !["syncmer", "kmc2", "msp"].contains(&method) {
        eprintln!("Unknown method: {}. Use 'syncmer', 'kmc2', or 'msp'.", method);
        std::process::exit(1);
    }

    eprintln!("Reading {}  k={}  l={}  method={}", fasta_path, k, l, method);

    let mut bucket_counts: HashMap<u32, u64> = HashMap::new();
    let mut total_kmers: u64 = 0;
    let mut total_superkmers: u64 = 0;
    let mut seq_count: u64 = 0;

    let file = fs::File::open(fasta_path).expect("Failed to open FASTA file");
    let reader = BufReader::new(file);

    let mut current_seq = Vec::<u8>::new();

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if line.starts_with('>') {
            if !current_seq.is_empty() {
                process_seq(&current_seq, k, l, method, &mut bucket_counts, &mut total_kmers, &mut total_superkmers);
                seq_count += 1;
                if seq_count % 10 == 0 {
                    eprintln!("  processed {} sequences, {} superkmers, {} kmers so far", seq_count, total_superkmers, total_kmers);
                }
                current_seq.clear();
            }
        } else {
            current_seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if !current_seq.is_empty() {
        process_seq(&current_seq, k, l, method, &mut bucket_counts, &mut total_kmers, &mut total_superkmers);
        seq_count += 1;
    }

    eprintln!("Done. {} sequences, {} superkmers, {} total kmers, {} distinct minimizers",
        seq_count, total_superkmers, total_kmers, bucket_counts.len());

    print_stats(&bucket_counts, total_kmers, total_superkmers, l);
}

fn process_seq(seq: &[u8], k: usize, l: usize, method: &str, bucket_counts: &mut HashMap<u32, u64>, total_kmers: &mut u64, total_superkmers: &mut u64) {
    match method {
        "syncmer" => {
            let (_storage, iter) = iteratorsyncmers2::SuperkmersIterator::new_with_n(seq, k, l);
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
