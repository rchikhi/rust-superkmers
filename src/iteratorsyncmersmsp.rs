use crate::Superkmer;
use debruijn::dna_string::DnaString;
use debruijn::kmer::{Kmer8, Kmer10, Kmer12};
use debruijn::Kmer;
use debruijn::msp::Scanner;
use lazy_static::lazy_static;

/// Extract superkmers from ASCII sequence, splitting on N/n characters.
pub fn superkmers_with_n(seq: &[u8], k: usize, l: usize) -> Vec<Superkmer> {
    let fragments = crate::utils::split_on_n(seq, k);
    let mut all = Vec::new();
    for (offset, fragment) in fragments {
        let dnastring = DnaString::from_acgt_bytes(fragment).to_bytes();
        let iter = SuperkmersIterator::new(&dnastring, k, l);
        for mut sk in iter {
            sk.start += offset;
            all.push(sk);
        }
    }
    all
}

pub struct SuperkmersIterator<'a> {
    iter: Box<dyn Iterator<Item = Superkmer> + 'a>,
}

const S: usize = 2;
const K8: usize = 8;
const K10: usize = 10;
const K12: usize = 12;

lazy_static! {
    static ref SYNCMERS_8: Vec<bool> = generate_syncmers::<K8>();
    static ref SYNCMERS_10: Vec<bool> = generate_syncmers::<K10>();
    static ref SYNCMERS_12: Vec<bool> = generate_syncmers::<K12>();
}


fn generate_syncmers<const K: usize>() -> Vec<bool> {
    let mut syncmers_arr = vec![false; 1 << (2 * K)];
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
        let syncmer = crate::syncmers::find_syncmers(K as usize, S, &[0, K - S], None, &kmer_bytes);
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

/// Convert an MSP result to a Superkmer, optionally canonicalizing the mint.
#[inline]
fn make_superkmer<K: Kmer>(msp_start: usize, msp_len: u16, msp_minimizer_pos: usize, mint_fwd: u64, l: usize, canonical: bool) -> Superkmer {
    let (mint, rc) = if canonical {
        let table = match l {
            8 => &*crate::CANONICAL_8,
            10 => &*crate::CANONICAL_10,
            12 => &*crate::CANONICAL_12,
            _ => panic!("Unsupported l={}", l),
        };
        table[mint_fwd as usize]
    } else {
        (mint_fwd as u32, false)
    };
    Superkmer {
        start: msp_start, mint,
        size: msp_len as u8, mpos: (msp_minimizer_pos - msp_start) as u8, rc,
    }
}

/* wrapper around rust-debruijn msp functions */
impl<'a> SuperkmersIterator<'a> {
    pub fn new(dnastring: &'a [u8], k: usize, l: usize) -> Self {
        Self::new_inner(dnastring, k, l, true)
    }

    /// Non-canonical version (forward-strand mint, rc=false).
    pub fn non_canonical(dnastring: &'a [u8], k: usize, l: usize) -> Self {
        Self::new_inner(dnastring, k, l, false)
    }

    fn new_inner(dnastring: &'a [u8], k: usize, l: usize, canonical: bool) -> Self {
        let dnastring = &DnaString::from_bytes(dnastring);
        let superkmer_iter: Box<dyn Iterator<Item = Superkmer>> = match l {
            8 => {
                let scanner8 = Scanner::new(dnastring, score8, k);
                let msps = scanner8.scan();
                Box::new(msps.into_iter().map(move |msp| {
                    make_superkmer::<Kmer8>(msp.start as usize, msp.len, msp.minimizer_pos as usize, msp.minimizer.to_u64(), l, canonical)
                }))
            }
            10 => {
                let scanner10 = Scanner::new(dnastring, score10, k);
                Box::new(scanner10.scan().into_iter().map(move |msp| {
                    make_superkmer::<Kmer10>(msp.start as usize, msp.len, msp.minimizer_pos as usize, msp.minimizer.to_u64(), l, canonical)
                }))
            }
            12 => {
                let scanner12 = Scanner::new(dnastring, score12, k);
                Box::new(scanner12.scan().into_iter().map(move |msp| {
                    make_superkmer::<Kmer12>(msp.start as usize, msp.len, msp.minimizer_pos as usize, msp.minimizer.to_u64(), l, canonical)
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
