pub fn revcomp(dna: &str) -> String {
    dna.chars()
        .rev()
        .map(switch_base)
        .collect::<String>()
}

fn switch_base(c: char) -> char {
    match c {
        'a' => 't',
        'c' => 'g',
        't' => 'a',
        'g' => 'c',
        'u' => 'a',
        'A' => 'T',
        'C' => 'G',
        'T' => 'A',
        'G' => 'C',
        'U' => 'A',
        _ => 'N'
    }
}

pub fn superkmer_to_verbose(superkmer :crate::Superkmer, read: &String, l :usize) -> crate::SuperkmerVerbose
{
        let sequence = &read[superkmer.start..superkmer.start+(superkmer.size as usize)];
        let sequence = if superkmer.rc { revcomp(&sequence) } else { sequence.to_string() };

        let mpos = superkmer.mpos as usize;
        let minimizer = sequence[mpos..mpos+l].to_string();

        let superkmer_verbose = crate::SuperkmerVerbose {
            sequence,
            minimizer,
            mpos: mpos.try_into().unwrap()
        };
        superkmer_verbose
}
