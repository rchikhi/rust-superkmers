// taken from syncmers crate, with a small fix

/// Syncmers as defined by Dutta et al. 2022, https://www.biorxiv.org/content/10.1101/2022.01.10.475696v2.full
/// Esp Fig 1b
/// Planning to implement other methods soon
///
/// TODO: Add Iterator impl's
// use std::iter::{FilterMap, Enumerate};
// use std::slice::Windows;
use std::cmp::Ordering;

//use pulp::Arch;
use twox_hash::XxHash64;
use std::hash::Hasher;

// TODO:Denote the reverse complement of x by Embedded Image. For a given order, the canonical form of a k-mer x, denoted by Canonical(x), is the smaller of x and Embedded Image. For example, under the lexicographic order, Canonical(CGGT) = ACCG.
// Canonical(x) = min(x, revcomp(x))

// Copied from ffforf. Really fast thanks to @sarah-ek
/// Complement a sequence, primarily used with revcomp function
pub fn complement(c: &mut u8) {
    let val = *c;
    let new_val = if val != b'N' {
        if val & 2 != 0 {
            val ^ 4
        } else {
            val ^ 21
        }
    } else {
        val
    };
    *c = new_val;
}

/*
/// Reverse complement a DNA Sequence. Case preserved
pub fn revcomp(sequence: &mut [u8]) {
    let arch = Arch::new();
    arch.dispatch(|| {
        sequence.reverse();
        sequence.iter_mut().for_each(complement);
    });
}*/

/// Test if the reverse complement is smaller, lexicographically, than the original sequence
/// Syncmers should be obtained from the Canonical strand (minimum strand).
/// Use this function for your own tests. It is not used in the Syncmer implementation.
pub fn is_revcomp_min(seq: &[u8]) -> bool {
    assert!(!seq.is_empty());
    for i in 0..seq.len() {
        let mut c = seq[seq.len() - i - 1];
        complement(&mut c);
        match seq[i].cmp(&c) {
            Ordering::Less => return true,
            Ordering::Greater => return false,
            Ordering::Equal => continue,
        }
    }

    false
}

// Best as determined by criterion benchmarks
// 303.62 MiB/s
/// Find syncmers from &[u8] and return Vec<&[u8]>
///
/// Parameterized syncmers as defined by Dutta et al. 2022, https://www.biorxiv.org/content/10.1101/2022.01.10.475696v2.full
/// Not all implemented yet (downsampling, windows, are not, for example).
///
/// # Arguments
/// k: kmer length
/// s: smer length
/// ts: Target positions, set at beginning or end for open/closed syncmers only.
///     Smallest smer must appear in one of these position of the kmer to be a valid syncmer
/// downsample fraction: None, or Some(float) between 0 and 1. If Some, only return some syncmers.
///
/// ```rust
/// # use rust_superkmers::syncmers::find_syncmers;
/// let sequence = b"CCAGTGTTTACGG";
/// let syncmers = find_syncmers(5, 2, &[2], None, sequence);
/// assert!(syncmers == vec![b"CCAGT", b"TTACG"]);
///
/// // You may also use multiple values for ts
/// let syncmers = find_syncmers(5, 2, &[2, 3], None, sequence);
/// ```
pub fn find_syncmers<'a, const N: usize>(
    k: usize,
    s: usize,
    ts: &[usize; N],
    downsample: Option<f64>,
    seq: &'a [u8],
) -> Vec<&'a [u8]> {
    assert!(seq.len() >= k);
    assert!(s < k);
    assert!(ts.iter().all(|&t| t <= k - s));
    assert!(N < 5);
    assert!(N == ts.len());

    let mut downsample_threshold = None;
    if let Some(downsample) = downsample {
        assert!(downsample > 0.0);
        assert!(downsample <= 1.0);
        downsample_threshold = Some((std::u64::MAX as f64 * downsample) as u64);
    }

    let syncmer_positions = find_syncmers_pos(k, s, ts, seq);

    if downsample.is_none() {
        syncmer_positions
            .iter()
            .map(|&pos| &seq[pos..pos + k])
            .collect()
    } else {
        syncmer_positions
            .iter()
            .map(|&pos| &seq[pos..pos + k])
            .filter(|&syncmer| {
                let mut hash = XxHash64::with_seed(42);
                hash.write(syncmer);
                let hash = hash.finish();
                hash < downsample_threshold.unwrap()
            })
            .collect()
    }
}

// Best as determined by criterion benchmarks
// 340.19 MiB/s
/// Find positions of syncmers
///
/// # Arguments
/// k: kmer length
/// s: smer length
/// ts: Target positions, set at beginning or end for open/closed syncmers only.
///    Smallest smer must appear in one of these position of the kmer to be a valid syncmer
///
/// # Returns
/// Vec<usize> of positions of syncmers (kmers meeting above critera) in the sequence
#[allow(clippy::if_same_then_else)]
pub fn find_syncmers_pos<const N: usize>(
    k: usize,
    s: usize,
    ts: &[usize; N],
    seq: &[u8],
) -> Vec<usize> {
    assert!(seq.len() >= k);
    assert!(s < k);
    assert!(ts.iter().all(|&t| t <= k - s));
    assert!(N < 5);
    assert!(N == ts.len());

    seq.windows(k)
        .enumerate()
        .filter_map(|(i, kmer)| {
            let min_pos = kmer
                .windows(s)
                .enumerate()
                .min_by(|(_, a), (_, b)| a.cmp(b));

            if N == 1 && ts[0] == min_pos.unwrap().0 {
                Some(i)
            } else if N != 1 && ts[0..N].contains(&min_pos.unwrap().0) {
                Some(i)
            } else {
                None
            }
        })
        .collect::<Vec<_>>()
}

#[cfg(test)]
mod test {
    #[test]
    pub fn test_syncmers_fig1b() {
        let sequence = b"CCAGTGTTTACGG";
        let syncmer_positions = crate::syncmers::find_syncmers_pos(5, 2, &[2], sequence);
        println!("{:?}", syncmer_positions);
        assert!(syncmer_positions == vec![0, 7]);

        let sequence = b"CCAGTGTTTACGG";
        let syncmers = crate::syncmers::find_syncmers(5, 2, &[2], None, sequence);
        assert!(syncmers == vec![b"CCAGT", b"TTACG"]);
        println!("{:?}", syncmers);

        let sequence = b"CCAGTGTTTACGG";
        let syncmer_positions = crate::syncmers::find_syncmers_pos(5, 2, &[2, 3], sequence);
        println!("{:?}", syncmer_positions);
        assert!(syncmer_positions == vec![0, 6, 7]);
    }
}
