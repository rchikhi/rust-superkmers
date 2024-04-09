pub fn revcomp(seq: &str) -> String
{
    std::str::from_utf8(&bio::alphabets::dna::revcomp(seq.as_bytes())).unwrap().to_string()
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
