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

pub fn normalize_mpos(sequence: String, sequence_rc: String, minimizer: String, minimizer_rc: String, mpos: usize, l: usize, mm: bool)
    -> (String, String, String, String, usize)
{
    let mut mpos = mpos;
    let mut sequence = sequence;
    let mut sequence_rc = sequence_rc;
    let mut minimizer = minimizer;
    let mut minimizer_rc = minimizer_rc;
    // if multiple minimizers, explore whole seq until mpos, otherwise just normalize by looking at revcomp
    if ! mm {
        if sequence.len()-(mpos+l) < mpos {
            (sequence, sequence_rc) = (sequence_rc, sequence);
            mpos = sequence.len()-(mpos+l);
        }
        if sequence.len()-(mpos+l) == mpos && sequence_rc < sequence {
            (sequence, sequence_rc) = (sequence_rc, sequence);
        }
    }
    else
    {
        for i in 0..mpos {
            let found_forward = sequence[i..i+l] == *minimizer  || sequence[i..i+l]  == *minimizer_rc;
            let found_rc      = sequence_rc[i..i+l] == *minimizer  || sequence_rc[i..i+l]  == *minimizer_rc;
            if found_forward || found_rc {
                if found_forward && found_rc { 
                    if sequence_rc < sequence { 
                        (sequence, sequence_rc) = (sequence_rc, sequence);
                    }
                } 
                #[allow(unused_assignments)]
                if found_rc {
                    (sequence, sequence_rc) = (sequence_rc, sequence);
                }
                mpos = i;
                break;
            }
        }
    }
    minimizer =  sequence[mpos..mpos+l].to_string();
    minimizer_rc = revcomp(&minimizer);
    (sequence, sequence_rc, minimizer, minimizer_rc, mpos)
}
