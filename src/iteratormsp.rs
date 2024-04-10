use crate::Superkmer; 
use debruijn::dna_string::DnaString;
use debruijn::kmer::{Kmer8, Kmer12};
use debruijn::Kmer;
use debruijn::msp::{Scanner};

pub struct SuperkmersIterator<'a> {
    iter: Box<dyn Iterator<Item = Superkmer> + 'a>,
}



 /* wrapper around rust-debruijn msp funtions
  */
impl<'a> SuperkmersIterator<'a> {
    pub fn new(dnastring: &'a DnaString, k: usize, l: usize) -> Self {
            let superkmer_iter : Box<dyn Iterator<Item = Superkmer>>;
            if l == 8 {
                let score8 = |p: &Kmer8| p.to_u64() as usize;
                let scanner8 = Scanner::new(dnastring, score8, k);
				superkmer_iter = Box::new(scanner8.scan().into_iter().map(|msp| Superkmer {
					start: msp.start as usize,
					size: msp.len as u8,
					mpos: (msp.minimizer_pos - msp.start) as u8,
					rc: false,
					mseq: msp.minimizer.to_string(),
				}));
            }
            else {  if l == 12 {
                let score12 = |p: &Kmer12| p.to_u64() as usize;
                let scanner12 = Scanner::new(dnastring, score12, k);
				superkmer_iter = Box::new(scanner12.scan().into_iter().map(|msp| Superkmer {
					start: msp.start as usize,
					size: msp.len as u8,
					mpos: (msp.minimizer_pos - msp.start) as u8,
					rc: false,
					mseq: msp.minimizer.to_string(),
				}));
            }
            else {
                panic!("unsupported l size for MSP iteration");
            } }

            SuperkmersIterator {
                iter: superkmer_iter,
            }
    }
}


// Implementing Iterator for SuperkmersIterator
impl<'a> Iterator for SuperkmersIterator<'a> {
    type Item = Superkmer;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}
