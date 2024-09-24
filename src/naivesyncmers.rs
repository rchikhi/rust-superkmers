use crate::utils::revcomp;
use crate::{Superkmer,SuperkmerVerbose};

// a naive O(nk) superkmer implementation that stores all superkmers of a sequence in a vector
// (problematic for long chromosomes)
#[allow(dead_code)]
pub fn extract_superkmers(read: &[u8], k: usize, l: usize) -> (Vec<Superkmer>, Vec<SuperkmerVerbose>) {
    let read_len = read.len();

    let S = 2;

    let read_str = std::str::from_utf8(&read[..]).unwrap().to_string();

    // Compute the hash sequence for each l-mer in the read
    let syncmers = crate::syncmers::find_syncmers(l, S, &[0, l-S], None, read);

    // Find the minimizer (syncmer) in each k-mer
    let mut kmer_minimizers_pos  = vec![usize::MAX; read_len - k + 1];
    let mut has_two_minimizer_occ = vec![false; read_len - k + 1];
    for i in 0..(read_len - k + 1) {
		// compute minimizer of that "window" (=kmer), i.e get the syncmer(s) inside
		let mut syncmer_count = 0;
		let mut first_syncmer_pos = None;

		for j in 0..(k - l + 1) {
			let lmer = &read[i + j..i + j + l];
			let lmer_str = &read_str[i + j..i + j + l];
			if syncmers.contains(&lmer) { // do we want this?: //|| syncmers.contains(&revcomp(&lmer_str).as_bytes())  {
				syncmer_count += 1;
				if first_syncmer_pos.is_none() {
					first_syncmer_pos = Some(i + j);
				}
			}
		}

		has_two_minimizer_occ[i] = syncmer_count > 1;
		if let Some(pos) = first_syncmer_pos {
			kmer_minimizers_pos[i] = pos;
		}
    }
                println!("yncmers {:?}",syncmers);
                println!("kmer minimizers pos {:?}",kmer_minimizers_pos);

    // Find super-kmers by scanning the read
    let mut super_kmers = vec![];
    let mut super_kmers_verbose = vec![];
    let mut start = 0;
    let mut end = k+1;
    let mut last_minimizer_pos = kmer_minimizers_pos[0];
    let debug = false;
    while end <= read_len {
        let kmer_minimizer_pos = kmer_minimizers_pos[end-k];
        let next_mm = has_two_minimizer_occ[end-k];
        if kmer_minimizer_pos != last_minimizer_pos { //|| next_mm {
            let sequence = std::str::from_utf8(&read[start..end-1]).unwrap().to_string();
            let sequence_rc = revcomp(&sequence);
            let mpos = last_minimizer_pos-start;
            let minimizer = sequence[mpos..mpos+l].to_string();
            let minimizer_rc = revcomp(&minimizer);
            let mm = has_two_minimizer_occ[end-k-1];
            let (new_sequence, _new_sequence_rc, new_minimizer, _new_minimizer_rc, new_mpos) = 
                normalize_mpos(&sequence, &sequence_rc, &minimizer, &minimizer_rc, mpos, l, mm);
            
            if debug 
            {
                println!("new superkmer: {} len {} minimizer {} pos {} mm {}",
                      new_sequence, 
                      sequence.len(),
                      new_minimizer,
                      new_mpos,
                      mm);
            }

            let superkmer = Superkmer {
                start,
                mint: 0, // don't record the minimizer here
                size: (end-1-start).try_into().unwrap(),
                mpos: new_mpos.try_into().unwrap(),
                rc: new_sequence == sequence_rc,
            };

            let superkmer_verbose = SuperkmerVerbose {
                sequence: new_sequence,
                minimizer: new_minimizer,
                mpos: new_mpos
            };

            super_kmers.push(superkmer);
            super_kmers_verbose.push(superkmer_verbose);
            start = end - k;
        }
        last_minimizer_pos = kmer_minimizer_pos;
        end += 1;
    }
    // output last superkmer
    if start <= read_len - k {
        let sequence = std::str::from_utf8(&read[start..read_len]).unwrap().to_string();
        let sequence_rc = revcomp(&sequence);
        let mpos = last_minimizer_pos-start;
        let minimizer = sequence[mpos..mpos+l].to_string();
        let minimizer_rc = revcomp(&minimizer);
        let mm = has_two_minimizer_occ[end-k-1];
        let (new_sequence, _new_sequence_rc, new_minimizer, _new_minimizer_rc, new_mpos) = 
            normalize_mpos(&sequence, &sequence_rc, &minimizer, &minimizer_rc, mpos, l, mm);
  
        if debug 
            {
                println!("new superkmer: {} len {} minimizer {} pos {} mm {}",
                      new_sequence, 
                      new_sequence.len(),
                      new_minimizer,
                      new_mpos,
                      mm);
            }

        let superkmer = Superkmer {
            start,
            mint: 0, // don't record the minimizer here
            size: (read_len-start).try_into().unwrap(),
            mpos: new_mpos.try_into().unwrap(),
            rc: new_sequence == sequence_rc
        };

        let superkmer_verbose = SuperkmerVerbose {
            sequence: new_sequence,
            minimizer: new_minimizer,
            mpos: new_mpos
        };

        super_kmers.push(superkmer);
        super_kmers_verbose.push(superkmer_verbose);
    }

    //println!("Maximum RSS at end of superkmers extraction: {:?}GB", (get_memory_rusage() as f32) / 1024.0 / 1024.0 / 1024.0);    
    (super_kmers, super_kmers_verbose)
}

// this function is actually agnostic to minimizer scheme
fn normalize_mpos(sequence: &String, sequence_rc: &String, minimizer: &String, minimizer_rc: &String, mpos: usize, l: usize, mm: bool)
    -> (String, String, String, String, usize)
{
    let mut mpos = mpos;
    let mut sequence = sequence;
    let mut sequence_rc = sequence_rc;
    let minimizer = minimizer;
    let minimizer_rc = minimizer_rc;
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
    let minimizer =  sequence[mpos..mpos+l].to_string();
    let minimizer_rc = revcomp(&minimizer);
    (sequence.to_string(), sequence_rc.to_string(), minimizer, minimizer_rc, mpos)
}
