//! Multi-minimizer superkmer iterator using the `multiminimizers` crate.
//!
//! Runs N independent minimizer schemes in parallel and picks the one yielding
//! the longest superkmer at each position. More hashes = fewer, longer superkmers.
//!
//! Constraint: `k - l` must be even (i.e. k and l must have the same parity).
//!
//! ```ignore
//! use rust_superkmers::iteratormultiminimizers::SuperkmersIterator;
//! let iter = SuperkmersIterator::new(seq, 31, 9, 2);           // canonical, 2 hashes
//! let iter = SuperkmersIterator::non_canonical(seq, 31, 9, 4); // forward-strand, 4 hashes
//! ```
use crate::Superkmer;

/// Lookup table: ASCII byte -> 2-bit encoding (A=0, C=1, G=2, T=3)
const ASCII_TO_2BIT: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

#[inline]
fn canonical_lmer_index(ascii: &[u8], pos: usize, l: usize) -> (usize, bool) {
    let mut fwd = 0usize;
    for i in 0..l {
        let base = ASCII_TO_2BIT[ascii[pos + i] as usize] as usize;
        fwd = (fwd << 2) | base;
    }
    let canonical_table = match l {
        8 => &*crate::CANONICAL_8,
        9 => &*crate::CANONICAL_9,
        10 => &*crate::CANONICAL_10,
        12 => &*crate::CANONICAL_12,
        _ => panic!("Unsupported l={} for canonical lookup", l),
    };
    let (canonical, is_rc) = canonical_table[fwd];
    (canonical as usize, is_rc)
}

#[inline]
fn forward_lmer_index(ascii: &[u8], pos: usize, l: usize) -> usize {
    let mut fwd = 0usize;
    for i in 0..l {
        let base = ASCII_TO_2BIT[ascii[pos + i] as usize] as usize;
        fwd = (fwd << 2) | base;
    }
    fwd
}

/// Collect superkmers from a single N-free fragment using multiminimizers.
fn superkmers_from_fragment<const N: usize>(
    ascii_slice: &[u8],
    k: usize,
    l: usize,
    offset: usize,
    canonical: bool,
    results: &mut Vec<Superkmer>,
) {
    if ascii_slice.len() < k {
        return;
    }

    assert!(
        (k - l) % 2 == 0,
        "multiminimizers requires k - l to be even (k={}, l={})",
        k, l
    );

    // multiminimizers requires uppercase ASCII
    let upper: Vec<u8> = ascii_slice.iter().map(|&b| b.to_ascii_uppercase()).collect();

    let iter = multiminimizers::compute_superkmers_linear_streaming::<N, true>(&upper, k, l);
    if let Some(iter) = iter {
        for sk in iter {
            let sk_start = sk.superkmer.start();
            let sk_end = sk.superkmer.end();
            let size = sk_end - sk_start;
            let mini_pos = sk.start_of_minimizer();
            let mpos = mini_pos - sk_start;

            let (mint, mint_is_rc) = if canonical {
                canonical_lmer_index(ascii_slice, mini_pos, l)
            } else {
                (forward_lmer_index(ascii_slice, mini_pos, l), false)
            };

            results.push(Superkmer {
                start: sk_start + offset,
                mint: mint as u32,
                size: size as u16,
                mpos: mpos as u16,
                mint_is_rc,
            });
        }
    }
}

/// Dispatch to the correct const-generic N at runtime.
fn superkmers_from_fragment_dispatch(
    ascii_slice: &[u8],
    k: usize,
    l: usize,
    offset: usize,
    canonical: bool,
    nb_hash: usize,
    results: &mut Vec<Superkmer>,
) {
    match nb_hash {
        1 => superkmers_from_fragment::<1>(ascii_slice, k, l, offset, canonical, results),
        2 => superkmers_from_fragment::<2>(ascii_slice, k, l, offset, canonical, results),
        3 => superkmers_from_fragment::<3>(ascii_slice, k, l, offset, canonical, results),
        4 => superkmers_from_fragment::<4>(ascii_slice, k, l, offset, canonical, results),
        8 => superkmers_from_fragment::<8>(ascii_slice, k, l, offset, canonical, results),
        16 => superkmers_from_fragment::<16>(ascii_slice, k, l, offset, canonical, results),
        _ => panic!("Unsupported nb_hash={}. Supported values: 1, 2, 3, 4, 8, 16", nb_hash),
    }
}

pub struct SuperkmersIterator {
    superkmers: Vec<Superkmer>,
    storage: Vec<u64>,
    pos: usize,
}

impl SuperkmersIterator {
    /// Process a sequence assumed to contain no N characters. Mint is canonical.
    pub fn new(seq: &[u8], k: usize, l: usize, nb_hash: usize) -> Self {
        Self::new_inner(seq, k, l, nb_hash, true)
    }

    /// Process a sequence that may contain N/n characters. Mint is canonical.
    pub fn new_with_n(seq: &[u8], k: usize, l: usize, nb_hash: usize) -> Self {
        Self::new_with_n_inner(seq, k, l, nb_hash, true)
    }

    /// Non-canonical version (forward-strand mint, mint_is_rc=false).
    pub fn non_canonical(seq: &[u8], k: usize, l: usize, nb_hash: usize) -> Self {
        Self::new_inner(seq, k, l, nb_hash, false)
    }

    /// Non-canonical version with N-splitting.
    pub fn non_canonical_with_n(seq: &[u8], k: usize, l: usize, nb_hash: usize) -> Self {
        Self::new_with_n_inner(seq, k, l, nb_hash, false)
    }

    /// Access the 2-bit packed representation of the input sequence.
    pub fn storage(&self) -> &[u64] {
        &self.storage
    }

    fn new_inner(seq: &[u8], k: usize, l: usize, nb_hash: usize, canonical: bool) -> Self {
        let storage = crate::utils::bitpack_fragment(seq);
        let mut superkmers = Vec::new();
        superkmers_from_fragment_dispatch(seq, k, l, 0, canonical, nb_hash, &mut superkmers);
        SuperkmersIterator { superkmers, storage, pos: 0 }
    }

    fn new_with_n_inner(seq: &[u8], k: usize, l: usize, nb_hash: usize, canonical: bool) -> Self {
        let storage = crate::utils::bitpack_fragment(seq);
        let fragments = crate::utils::split_on_n(seq, k);
        let mut superkmers = Vec::new();
        for (offset, fragment) in &fragments {
            superkmers_from_fragment_dispatch(fragment, k, l, *offset, canonical, nb_hash, &mut superkmers);
        }
        SuperkmersIterator { superkmers, storage, pos: 0 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_output() {
        let seq = b"AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCACACTGCACTGCACTGACTGCACTGCACTGCACTGCACTGCCTGC";
        let sks: Vec<_> = SuperkmersIterator::new(seq, 31, 9, 2).collect();
        assert!(!sks.is_empty(), "should produce at least one superkmer");
        // Last superkmer should reach the end of the sequence
        let last = sks.last().unwrap();
        assert_eq!(last.start as usize + last.size as usize, seq.len(),
            "superkmers should cover to end of sequence");
        for sk in &sks {
            assert!(sk.size >= 31, "superkmer size must be >= k");
            assert!((sk.mpos as usize) + 9 <= sk.size as usize, "minimizer must fit within superkmer");
        }
    }

    #[test]
    fn test_nb_hash_values() {
        let seq = b"AACTGCACTGCACTGCACTGCACACTGCACTGCACTGCACTGCACACTGCACTGCACTGACTGCACTGCACTGCACTGCACTGCCTGC";
        // All supported nb_hash values should work
        for nb in [1, 2, 3, 4, 8] {
            let sks: Vec<_> = SuperkmersIterator::new(seq, 31, 9, nb).collect();
            assert!(!sks.is_empty(), "nb_hash={} should produce superkmers", nb);
        }
    }
}

impl Iterator for SuperkmersIterator {
    type Item = Superkmer;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.superkmers.len() {
            return None;
        }
        let idx = self.pos;
        self.pos += 1;
        Some(Superkmer {
            start: self.superkmers[idx].start,
            mint: self.superkmers[idx].mint,
            size: self.superkmers[idx].size,
            mpos: self.superkmers[idx].mpos,
            mint_is_rc: self.superkmers[idx].mint_is_rc,
        })
    }
}
