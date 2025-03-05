// based on rust-seq2kminmers 

#[macro_use]
extern crate criterion;

use criterion::{Bencher, Criterion, Throughput, BenchmarkId, black_box};
use rand::distributions::{Distribution, Uniform};

use nthash::{nthash, NtHashIterator};

//use cocktail::tokenizer::TokenizerMini;
use debruijn::dna_string::DnaString;
use debruijn::kmer::{Kmer8};
use debruijn::Kmer;
use debruijn::msp::Scanner;

fn superkmers_bench(c: &mut Criterion) {
    let k = 31;
    let m = 8;
    let range = Uniform::from(0..4);
    let mut rng = rand::thread_rng();
    let seq_len = 150;
    let seq = (0..seq_len)
        .map(|_| match range.sample(&mut rng) {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            //_ => 'N',
            _ => 'A',
        })
        .collect::<String>();
    let dnastring = DnaString::from_dna_string(&seq);

    let mut group = c.benchmark_group("BenchmarkGroup");
    group.throughput(Throughput::Bytes(seq.len() as u64));

    group.bench_with_input(BenchmarkId::new("superkmers_iteratorsyncmers2", seq_len), &seq,
	|b: &mut Bencher, i: &String| {
		b.iter(|| {
            let iter = rust_superkmers::iteratorsyncmers2::SuperkmersIterator::new(i.as_bytes(), k, m);
            let _res = iter.1.into_iter().collect::<Vec<rust_superkmers::Superkmer>>(); 
			black_box(_res);
		})});


    group.bench_with_input(BenchmarkId::new("superkmers_iteratorsyncmersmsp", seq_len), &seq,
	|b: &mut Bencher, i: &String| {
		b.iter(|| {
            let binding = dnastring.to_bytes();
            let iter = rust_superkmers::iteratorsyncmersmsp::SuperkmersIterator::new(&binding, k, m);
            let _res = iter.into_iter().collect::<Vec<rust_superkmers::Superkmer>>(); 
			black_box(_res);
		})});


    group.bench_with_input(BenchmarkId::new("syncmers", seq_len), &seq,
	|b: &mut Bencher, i: &String| {
		b.iter(|| {
            let syncmers = rust_superkmers::syncmers::find_syncmers(8, 2, &[0,6], None, &seq.as_bytes());
			black_box(syncmers);
		})});



    group.bench_with_input(BenchmarkId::new("rust_debruijn_msp", seq_len), &seq,
	|b: &mut Bencher, i: &String| {
		b.iter(|| {
            let score = |p: &Kmer8| p.to_u64() as usize;
            let scanner = Scanner::new(&dnastring, score, k);
            let _res = scanner.scan();
			black_box(_res);
		})});


    group.bench_with_input(BenchmarkId::new("superkmers_iterator1", seq_len), &seq,
    |b: &mut Bencher, i: &String| {
        b.iter(|| {
            let iter = rust_superkmers::iterator1::SuperkmersIterator::new(i.as_bytes(), k, m);
            let _res = iter.into_iter().collect::<Vec<rust_superkmers::Superkmer>>(); 
            black_box(_res);
    })});

    /*
     // it was about 68 MB/s for 150bp sequences, k=32, l=8
     // to enable: uncomment cocktail in Cargo.toml
	group.bench_with_input(BenchmarkId::new("cocktail_tokenizermini", seq_len), &seq,
	|b: &mut Bencher, i: &String| {
		b.iter(|| {
			let tokenizer = TokenizerMini::new(i.as_bytes(), k as u8, m as u8);
			let _res = tokenizer.into_iter().collect::<Vec<(u64,u64)>>(); 
			black_box(_res);
		})});
    */

    group.bench_with_input(BenchmarkId::new("superkmers_naive", seq_len), &seq,
    |b: &mut Bencher, i: &String| {
        b.iter(|| {
            let (iter, iter_verbose) = rust_superkmers::naive::extract_superkmers(i.as_bytes(), k, m);
            let _res = iter.into_iter().collect::<Vec<rust_superkmers::Superkmer>>();
            black_box(_res);
    })});

    group.bench_with_input(BenchmarkId::new("nthashiterator_luiz", seq_len), &seq,
    |b: &mut Bencher, i: &String| {
        b.iter(|| {
            let iter = NtHashIterator::new(i.as_bytes(), k).unwrap();
            //  iter.for_each(drop);
            let _res = iter.collect::<Vec<u64>>(); // original nthash iterator only has 64 bits
            black_box(_res);
        })});

     group.bench_with_input(BenchmarkId::new("nthash (64)", seq_len), &seq,
    |b: &mut Bencher, i: &String| {
        b.iter(|| {
            nthash(i.as_bytes(), k);
        })});



    group.finish();
}

criterion_group!(benches, superkmers_bench);
criterion_main!(benches);
