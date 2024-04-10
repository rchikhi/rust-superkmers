use crate::Superkmer; 
use debruijn::dna_string::DnaString;
use debruijn::kmer::{Kmer8};
use debruijn::Kmer;
use debruijn::msp::{Scanner};

pub struct SuperkmersIterator<'a> {
    iter: Box<dyn Iterator<Item = Superkmer> + 'a>,
}

 /* wrapper around rust-debruijn msp funtions
  */
impl<'a> SuperkmersIterator<'a> {
    pub fn new(dnastring: &'a DnaString, k: usize, l: usize) -> Self {
            if l != 8 {
                panic!("unsupported l size for MSP iteration");
            }

            let score = |p: &Kmer8| p.to_u64() as usize;
            let scanner = Scanner::new(dnastring, score, k);

            let superkmer_iter = scanner.scan().into_iter().map(|msp| Superkmer {
                start: msp.start as usize,
                size: msp.len as u8,
                mpos: (msp.minimizer_pos - msp.start) as u8,
                rc: false,
            });

            SuperkmersIterator {
                iter: Box::new(superkmer_iter),
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
