use debruijn::dna_string::DnaString;
use std::time::Instant;

fn random_dna(len: usize, seed: u64) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
            bases[((x >> 33) % 4) as usize]
        })
        .collect()
}

fn bench<F: FnMut(&[u8])>(name: &str, seq: &[u8], iters: u32, mut f: F) {
    // warmup
    for _ in 0..3 {
        f(seq);
    }
    let start = Instant::now();
    for _ in 0..iters {
        f(seq);
    }
    let elapsed = start.elapsed();
    let ns_per_iter = elapsed.as_nanos() as f64 / iters as f64;
    let mb_per_sec = seq.len() as f64 / ns_per_iter * 1000.0;
    println!("{:25} {:10.0} ns    {:8.1} MB/s", name, ns_per_iter, mb_per_sec);
}

fn main() {
    let k = 31;
    for &len in &[150, 10_000, 1_000_000] {
        let seq = random_dna(len, 42);
        let iters = if len <= 150 {
            100_000
        } else if len <= 10_000 {
            5_000
        } else {
            50
        };

        let dna_bytes = DnaString::from_acgt_bytes(&seq).to_bytes();

        println!("\n--- seq_len={} k={} ({} iters) ---", len, k, iters);

        bench("syncmers2 (l=8)", &seq, iters, |s| {
            let iter =
                rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(s, k, 8);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("syncmers2:classical", &seq, iters, |s| {
            let iter =
                rust_superkmers::iteratorsyncmers2::SuperkmersIterator::classical(s, k, 8);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("syncmers2:msp", &seq, iters, |s| {
            let iter =
                rust_superkmers::iteratorsyncmers2::SuperkmersIterator::msp(s, k, 8);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("syncmers2:mspxor", &seq, iters, |s| {
            let iter =
                rust_superkmers::iteratorsyncmers2::SuperkmersIterator::mspxor(s, k, 8);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("kmc2 (l=8)", &seq, iters, |s| {
            let (_, iter) =
                rust_superkmers::iteratorkmc2::SuperkmersIterator::new(s, k, 8);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        let dna_bytes_ref = &dna_bytes;
        bench("syncmersmsp (l=8)", &seq, iters, |_s| {
            let iter = rust_superkmers::iteratorsyncmersmsp::SuperkmersIterator::new(dna_bytes_ref, k, 8);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("msp (l=8)", &seq, iters, |_s| {
            let iter = rust_superkmers::iteratormsp::SuperkmersIterator::new(dna_bytes_ref, k, 8);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        #[cfg(feature = "simd-mini")]
        {
        bench("simdmini (l=9)", &seq, iters, |s| {
            let iter =
                rust_superkmers::iteratorsimdmini::SuperkmersIterator::new(s, k, 9);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("simdmini:classical", &seq, iters, |s| {
            let iter =
                rust_superkmers::iteratorsimdmini::SuperkmersIterator::classical(s, k, 9);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("simdmini:msp", &seq, iters, |s| {
            let iter =
                rust_superkmers::iteratorsimdmini::SuperkmersIterator::msp(s, k, 9);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("simdmini:mspxor", &seq, iters, |s| {
            let iter =
                rust_superkmers::iteratorsimdmini::SuperkmersIterator::mspxor(s, k, 9);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });
        }

        // Extractor benchmarks (reusable buffers, amortized across reads)
        {
            let mut ext2 = rust_superkmers::iteratorsyncmers2::SuperkmerExtractor::new(k, 8);
            bench("syncmers2-ext (l=8)", &seq, iters, |s| {
                let sks = ext2.process(s);
                std::hint::black_box(sks);
            });
        }

        {
            let mut ext2 = rust_superkmers::iteratorsyncmers2::SuperkmerExtractor::mspxor(k, 8);
            bench("syncmers2-ext:mspxor", &seq, iters, |s| {
                let sks = ext2.process(s);
                std::hint::black_box(sks);
            });
        }

        #[cfg(feature = "simd-mini")]
        {
            let mut ext_simd = rust_superkmers::iteratorsimdmini::SuperkmerExtractor::new(k, 9);
            bench("simdmini-ext (l=9)", &seq, iters, |s| {
                let sks = ext_simd.process(s);
                std::hint::black_box(sks);
            });
        }

        #[cfg(feature = "simd-mini")]
        {
            let mut ext_simd = rust_superkmers::iteratorsimdmini::SuperkmerExtractor::mspxor(k, 9);
            bench("simdmini-ext:mspxor", &seq, iters, |s| {
                let sks = ext_simd.process(s);
                std::hint::black_box(sks);
            });
        }
    }
}
