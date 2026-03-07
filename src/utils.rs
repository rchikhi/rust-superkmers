pub fn revcomp(seq: &str) -> String
{
    std::str::from_utf8(&bio::alphabets::dna::revcomp(seq.as_bytes())).unwrap().to_string()
}

/// Split a byte sequence on N/n characters, returning (offset, fragment) pairs.
/// Fragments shorter than `min_len` are skipped.
pub fn split_on_n(seq: &[u8], min_len: usize) -> Vec<(usize, &[u8])> {
    let mut fragments = Vec::new();
    let mut start = 0;
    for (i, &b) in seq.iter().enumerate() {
        if b == b'N' || b == b'n' {
            if i - start >= min_len {
                fragments.push((start, &seq[start..i]));
            }
            start = i + 1;
        }
    }
    if seq.len() - start >= min_len {
        fragments.push((start, &seq[start..]));
    }
    fragments
}

pub fn superkmer_to_verbose(superkmer :crate::Superkmer, read: &String, l :usize) -> crate::SuperkmerVerbose
{
        let sequence = &read[superkmer.start..superkmer.start+(superkmer.size as usize)];
        let sequence = if superkmer.rc { revcomp(sequence) } else { sequence.to_string() };

        let mpos = superkmer.mpos as usize;
        let minimizer = sequence[mpos..mpos+l].to_string();

        let superkmer_verbose = crate::SuperkmerVerbose {
            sequence,
            minimizer,
            mpos: mpos.try_into().unwrap()
        };
        superkmer_verbose
}
