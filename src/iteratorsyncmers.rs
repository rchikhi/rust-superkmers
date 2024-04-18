use crate::Superkmer;
use debruijn::dna_string::DnaString;
use debruijn::kmer::{Kmer8, Kmer10, Kmer12};
use debruijn::Kmer;
use debruijn::msp::Scanner;

pub struct SuperkmersIterator<'a> {
    iter: Box<dyn Iterator<Item = Superkmer> + 'a>,
}

use lazy_static::lazy_static;
use std::sync::Mutex;

lazy_static! {
    static ref SYNCMERS_8: [bool; 65536] = {
        let mut syncmers_arr = [false; 65536];
        for kmer_int in 0..65536 {
            let kmer_bytes = {
                let mut bytes = [0u8; 8];
                for (i, byte) in bytes.iter_mut().enumerate() {
                    *byte = match (kmer_int >> (2 * (7 - i))) & 3 {
                        0 => b'A',
                        1 => b'C',
                        2 => b'G',
                        3 => b'T',
                        _ => unreachable!(),
                    };
                }
                bytes
            };
			let syncmer = crate::syncmers::find_syncmers(8, 2, &[0, 6], None, &kmer_bytes);
            syncmers_arr[kmer_int] = !syncmer.is_empty();
        }
        syncmers_arr
    };
}


fn score8(p: &Kmer8) -> usize {
    let kmer = p.to_u64() as usize;
    if SYNCMERS_8[kmer] {
        0
    } else {
        std::usize::MAX
    }
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

