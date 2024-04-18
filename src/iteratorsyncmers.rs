use crate::Superkmer;
use debruijn::dna_string::DnaString;
use debruijn::kmer::{Kmer8, Kmer10, Kmer12};
use debruijn::Kmer;
use debruijn::msp::Scanner;
use lazy_static::lazy_static;

pub struct SuperkmersIterator<'a> {
    iter: Box<dyn Iterator<Item = Superkmer> + 'a>,
}


const K8: usize = 8;
const K10: usize = 10;
const K12: usize = 12;

lazy_static! {
    static ref SYNCMERS_8: [bool; 1 << (2 * K8)] = generate_syncmers::<K8>();
    static ref SYNCMERS_10: [bool; 1 << (2 * K10)] = generate_syncmers::<K10>();
    static ref SYNCMERS_12: [bool; 1 << (2 * K12)] = generate_syncmers::<K12>();
}

fn generate_syncmers<const K: usize>() -> [bool; 1 << (2 * K)] {
    let mut syncmers_arr = [false; 1 << (2 * K)];
    for kmer_int in 0..(1 << (2 * K)) {
        let kmer_bytes = {
            let mut bytes = [0u8; K];
            for (i, byte) in bytes.iter_mut().enumerate() {
                *byte = match (kmer_int >> (2 * (K - 1 - i))) & 3 {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    3 => b'T',
                    _ => unreachable!(),
                };
            }
            bytes
        };
        let syncmer = crate::syncmers::find_syncmers(K as usize, 2, &[0, K - 2], None, &kmer_bytes);
        syncmers_arr[kmer_int] = !syncmer.is_empty();
    }
    syncmers_arr
}

fn score8(p: &Kmer8) -> usize {
    let kmer = p.to_u64() as usize;
    !SYNCMERS_8[kmer] as usize
}

fn score10(p: &Kmer10) -> usize {
    let kmer = p.to_u64() as usize;
    !SYNCMERS_10[kmer] as usize
}

fn score12(p: &Kmer12) -> usize {
    let kmer = p.to_u64() as usize;
    !SYNCMERS_12[kmer] as usize
}

 /* wrapper around rust-debruijn msp funtions
  */
impl<'a> SuperkmersIterator<'a> {
    pub fn new(dnastring: &'a [u8], k: usize, l: usize) -> Self {
        let dnastring = &DnaString::from_bytes(dnastring);
		let superkmer_iter: Box<dyn Iterator<Item = Superkmer>> = match l {
			8 => {
				let scanner8 = Scanner::new(dnastring, score8, k);
				Box::new(scanner8.scan().into_iter().map(|msp| Superkmer {
					start: msp.start as usize,
					mint: msp.minimizer.to_u64() as u32,
					size: msp.len as u8,
					mpos: (msp.minimizer_pos - msp.start) as u8,
					rc: false,
				}))
			}
			10 => {
				let scanner10 = Scanner::new(dnastring, score10, k);
				Box::new(scanner10.scan().into_iter().map(|msp| Superkmer {
					start: msp.start as usize,
					mint: msp.minimizer.to_u64() as u32,
					size: msp.len as u8,
					mpos: (msp.minimizer_pos - msp.start) as u8,
                    rc: false,
                }))
            }
            12 => {
                let scanner12 = Scanner::new(dnastring, score12, k);
                Box::new(scanner12.scan().into_iter().map(|msp| Superkmer {
                    start: msp.start as usize,
                    mint: msp.minimizer.to_u64() as u32,
                    size: msp.len as u8,
                    mpos: (msp.minimizer_pos - msp.start) as u8,
                    rc: false,
                }))
            }
            _ => panic!("Unsupported l size for MSP iteration"),
        };
        SuperkmersIterator { iter: superkmer_iter }
    }
}

impl<'a> Iterator for SuperkmersIterator<'a> {
    type Item = Superkmer;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

