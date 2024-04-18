use crate::Superkmer;
use debruijn::dna_string::DnaString;
use debruijn::kmer::{Kmer8, Kmer10, Kmer12};
use debruijn::Kmer;
use debruijn::msp::Scanner;

pub struct SuperkmersIterator<'a> {
    iter: Box<dyn Iterator<Item = Superkmer> + 'a>,
}

fn score8(p: &Kmer8) -> usize {
    if p.to_u64() == 0 { std::usize::MAX } else { p.to_u64() as usize }
}

fn score10(p: &Kmer10) -> usize {
    if p.to_u64() == 0 { std::usize::MAX } else { p.to_u64() as usize }
}

fn score12(p: &Kmer12) -> usize {
    if p.to_u64() == 0 { std::usize::MAX } else { p.to_u64() as usize }
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

