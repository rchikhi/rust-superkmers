// based on rust-seq2kminmers 

#[macro_use]
extern crate criterion;

use criterion::{Bencher, Criterion, Throughput, BenchmarkId, black_box};
use rand::distributions::{Distribution, Uniform};

use nthash::{nthash, NtHashIterator};
use nthash32;
use rust_seq2kminmers::{NtHashHPCIterator, NtHashSIMDIterator, HashMode};

fn superkmers_bench(c: &mut Criterion) {
    let range = Uniform::from(0..4);
    let mut rng = rand::thread_rng();
    let seq_len = 100000;
    let seq = (0..seq_len)
        .map(|_| match range.sample(&mut rng) {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => 'N',
        })
        .collect::<String>();

    let mut group = c.benchmark_group("BenchmarkGroup");
    group.throughput(Throughput::Bytes(seq.len() as u64));


    group.bench_with_input(BenchmarkId::new("nthash (64)", seq_len), &seq,
    |b: &mut Bencher, i: &String| {
        b.iter(|| {
            nthash(i.as_bytes(), 5);
        })});

    group.bench_with_input(BenchmarkId::new("nthashiterator_luiz", seq_len), &seq,
    |b: &mut Bencher, i: &String| {
        b.iter(|| {
            let iter = NtHashIterator::new(i.as_bytes(), 5).unwrap();
            //  iter.for_each(drop);
            let _res = iter.collect::<Vec<u64>>(); // original nthash iterator only has 64 bits
            black_box(_res);
        })});

    group.finish();
}

criterion_group!(benches, superkmers_bench);
criterion_main!(benches);
