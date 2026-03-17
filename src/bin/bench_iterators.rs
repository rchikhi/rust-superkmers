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

fn bench<F: FnMut(&[u8])>(name: &str, seqs: &[Vec<u8>], iters: u32, filter: Option<&str>, mut f: F) {
    if let Some(pat) = filter {
        if !name.contains(pat) { return; }
    }
    let seq_len = seqs[0].len();
    // warmup
    for i in 0..3 {
        f(&seqs[i % seqs.len()]);
    }
    let start = Instant::now();
    for i in 0..iters {
        f(&seqs[i as usize % seqs.len()]);
    }
    let elapsed = start.elapsed();
    let ns_per_iter = elapsed.as_nanos() as f64 / iters as f64;
    let mb_per_sec = seq_len as f64 / ns_per_iter * 1000.0;
    println!("{:25} {:10.0} ns    {:8.1} MB/s", name, ns_per_iter, mb_per_sec);
}

fn main() {
    let k = 31;
    let filter: Option<String> = std::env::args().nth(1);
    let filter_ref = filter.as_deref();
    for &len in &[150, 10_000, 1_000_000] {
        // Generate a pool of distinct sequences to defeat branch predictor/cache learning
        let pool_size = if len <= 150 { 1024 } else if len <= 10_000 { 128 } else { 4 };
        let seqs: Vec<Vec<u8>> = (0..pool_size).map(|i| random_dna(len, 42 + i as u64)).collect();
        let iters = if len <= 150 {
            100_000
        } else if len <= 10_000 {
            5_000
        } else {
            50
        };

        let dna_bytes = DnaString::from_acgt_bytes(&seqs[0]).to_bytes();

        println!("\n--- seq_len={} k={} ({} iters, {} seqs) ---", len, k, iters, pool_size);

        bench("syncmers2 (l=8)", &seqs, iters, filter_ref, |s| {
            let iter =
                rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(s, k, 8);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("syncmers2:classical", &seqs, iters, filter_ref, |s| {
            let iter =
                rust_superkmers::iteratorsyncmers2::SuperkmersIterator::classical(s, k, 8);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("syncmers2:msp", &seqs, iters, filter_ref, |s| {
            let iter =
                rust_superkmers::iteratorsyncmers2::SuperkmersIterator::msp(s, k, 8);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("syncmers2:mspxor", &seqs, iters, filter_ref, |s| {
            let iter =
                rust_superkmers::iteratorsyncmers2::SuperkmersIterator::mspxor(s, k, 8);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("kmc2 (l=8)", &seqs, iters, filter_ref, |s| {
            let (_, iter) =
                rust_superkmers::iteratorkmc2::SuperkmersIterator::new(s, k, 8);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        let dna_bytes_ref = &dna_bytes;
        bench("syncmersmsp (l=8)", &seqs, iters, filter_ref, |_s| {
            let iter = rust_superkmers::iteratorsyncmersmsp::SuperkmersIterator::new(dna_bytes_ref, k, 8);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("msp (l=8)", &seqs, iters, filter_ref, |_s| {
            let iter = rust_superkmers::iteratormsp::SuperkmersIterator::new(dna_bytes_ref, k, 8);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        #[cfg(feature = "simd-mini")]
        {
        bench("simdmini (l=9)", &seqs, iters, filter_ref, |s| {
            let iter =
                rust_superkmers::iteratorsimdmini::SuperkmersIterator::new(s, k, 9);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("simdmini:classical", &seqs, iters, filter_ref, |s| {
            let iter =
                rust_superkmers::iteratorsimdmini::SuperkmersIterator::classical(s, k, 9);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("simdmini:msp", &seqs, iters, filter_ref, |s| {
            let iter =
                rust_superkmers::iteratorsimdmini::SuperkmersIterator::msp(s, k, 9);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("simdmini:mspxor", &seqs, iters, filter_ref, |s| {
            let iter =
                rust_superkmers::iteratorsimdmini::SuperkmersIterator::mspxor(s, k, 9);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("cminim (l=9)", &seqs, iters, filter_ref, |s| {
            let iter =
                rust_superkmers::iteratorsimdmini_cminim::SuperkmersIterator::new(s, k, 9);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });
        }

        // Syncmer detection only (no superkmer extraction)
        {
            let scores = rust_superkmers::iteratorsyncmers2::mspxor_syncmer_scores(8);
            bench("syncmers2 scores-only l=8", &seqs, iters, filter_ref, |s| {
                let storage = rust_superkmers::utils::bitpack_fragment(s);
                let mut total = 0usize;
                for i in 0..s.len().saturating_sub(7) {
                    total += scores[rust_superkmers::iteratorsyncmers2::get_kmer_value(&storage, i, 8)];
                }
                std::hint::black_box(total);
            });
        }

        {
            let scores = rust_superkmers::iteratorsyncmers2::mspxor_syncmer_scores(9);
            bench("syncmers2 scores-only l=9", &seqs, iters, filter_ref, |s| {
                let storage = rust_superkmers::utils::bitpack_fragment(s);
                let mut total = 0usize;
                for i in 0..s.len().saturating_sub(8) {
                    total += scores[rust_superkmers::iteratorsyncmers2::get_kmer_value(&storage, i, 9)];
                }
                std::hint::black_box(total);
            });
        }

        #[cfg(feature = "simd-mini")]
        {
            use simd_minimizers::packed_seq::AsciiSeq;
            let mut positions = Vec::new();
            bench("simdmini SIMD-detect l=9", &seqs, iters, filter_ref, |s| {
                positions.clear();
                simd_minimizers::canonical_closed_syncmers(2, 8)
                    .run(AsciiSeq(s), &mut positions);
                std::hint::black_box(positions.len());
            });
        }

        // Extractor benchmarks (reusable buffers, amortized across reads)
        {
            let mut ext2 = rust_superkmers::iteratorsyncmers2::SuperkmerExtractor::new(k, 8);
            bench("syncmers2-ext (l=8)", &seqs, iters, filter_ref, |s| {
                let sks = ext2.process(s);
                std::hint::black_box(sks);
            });
        }

        {
            let mut ext2 = rust_superkmers::iteratorsyncmers2::SuperkmerExtractor::mspxor(k, 8);
            bench("syncmers2-ext:mspxor", &seqs, iters, filter_ref, |s| {
                let sks = ext2.process(s);
                std::hint::black_box(sks);
            });
        }

        {
            let mut ext2 = rust_superkmers::iteratorsyncmers2::SuperkmerExtractor::mspxor(k, 9);
            bench("syncmers2-ext:mspxor l=9", &seqs, iters, filter_ref, |s| {
                let sks = ext2.process(s);
                std::hint::black_box(sks);
            });
        }

        #[cfg(feature = "simd-mini")]
        {
            let mut ext_simd = rust_superkmers::iteratorsimdmini::SuperkmerExtractor::new(k, 9);
            bench("simdmini-ext (l=9)", &seqs, iters, filter_ref, |s| {
                let sks = ext_simd.process(s);
                std::hint::black_box(sks);
            });
        }

        #[cfg(feature = "simd-mini")]
        {
            let mut ext_simd = rust_superkmers::iteratorsimdmini::SuperkmerExtractor::mspxor(k, 9);
            bench("simdmini-ext:mspxor", &seqs, iters, filter_ref, |s| {
                let sks = ext_simd.process(s);
                std::hint::black_box(sks);
            });
        }

        // UHS benchmarks
        bench("uhs (l=8)", &seqs, iters, filter_ref, |s| {
            let iter =
                rust_superkmers::iteratoruhs::SuperkmersIterator::new(s, k, 8);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        bench("uhs:mspxor (l=8)", &seqs, iters, filter_ref, |s| {
            let iter =
                rust_superkmers::iteratoruhs::SuperkmersIterator::mspxor(s, k, 8);
            let v: Vec<_> = iter.collect();
            std::hint::black_box(v);
        });

        {
            let mut ext_uhs = rust_superkmers::iteratoruhs::SuperkmerExtractor::new(k, 8);
            bench("uhs-ext (l=8)", &seqs, iters, filter_ref, |s| {
                let sks = ext_uhs.process(s);
                std::hint::black_box(sks);
            });
        }

        {
            let mut ext_uhs = rust_superkmers::iteratoruhs::SuperkmerExtractor::mspxor(k, 8);
            bench("uhs-ext:mspxor (l=8)", &seqs, iters, filter_ref, |s| {
                let sks = ext_uhs.process(s);
                std::hint::black_box(sks);
            });
        }
    }
}
