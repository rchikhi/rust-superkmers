use crate::utils::{revcomp, normalize_mpos};
use nthash::NtHashIterator;
use crate::Superkmer;

// a naive O(nk) superkmer implementation that stores all superkmers of a sequence in a vector
// (problematic for long chromosomes)
#[allow(dead_code)]
pub fn extract_superkmers(read: &[u8], k: usize, l: usize) -> Vec<Superkmer> {
    let read_len = read.len();

    // Compute the hash sequence for each l-mer in the read
    let nthash_iter = NtHashIterator::new(read, l).unwrap();
    let hashes : Vec<u64> = nthash_iter.collect();

    // Find the minimizer in each k-mer
    let mut kmer_minimizers_pos  = vec![0; read_len - k + 1];
    let mut has_two_minimizer_occ = vec![false; read_len - k + 1];
    for i in 0..(read_len - k + 1) {
        let mut min_minimizer     = hashes[i];
        let mut min_minimizer_pos = i;
        let mut nb_minimizer_pos = 1;
        for j in 1..k-l+1 {
            let minimizer = hashes[i+j];
            if minimizer == min_minimizer {
                nb_minimizer_pos += 1;
            }
            if minimizer < min_minimizer {
                min_minimizer     = minimizer;
                min_minimizer_pos = i+j;
                nb_minimizer_pos  = 1;
            }
        }
        has_two_minimizer_occ[i] = nb_minimizer_pos > 1;
        kmer_minimizers_pos[i] = min_minimizer_pos;
        // mm == multiple minimizers in kmer
        //println!("kmer: {} minimizer pos {} mm {}",std::str::from_utf8(&read[i..i+k]).unwrap(),
        //                                        kmer_minimizers_pos[i],has_two_minimizer_occ[i]);
    }

    // Find super-kmers by scanning the read
    let mut super_kmers = vec![];
    let mut start = 0;
    let mut end = k+1;
    let mut last_minimizer_pos = kmer_minimizers_pos[0];
    let debug = false;
    while end <= read_len {
        let kmer_minimizer_pos = kmer_minimizers_pos[end-k];
        let next_mm = has_two_minimizer_occ[end-k];
        if kmer_minimizer_pos != last_minimizer_pos || next_mm {
            let mut sequence = std::str::from_utf8(&read[start..end-1]).unwrap().to_string();
            let mut sequence_rc = revcomp(&sequence);
            let mut mpos = last_minimizer_pos-start;
            let mut minimizer = sequence[mpos..mpos+l].to_string();
            let mut minimizer_rc = revcomp(&minimizer);
            let mm = has_two_minimizer_occ[end-k-1];
            (sequence, sequence_rc, minimizer, minimizer_rc, mpos) = 
                normalize_mpos(sequence, sequence_rc, minimizer, minimizer_rc, mpos, l, mm);
            
            if debug 
            {
                println!("new superkmer: {} len {} minimizer {} pos {} mm {}",
                      sequence, 
                      sequence.len(),
                      minimizer,
                      mpos,
                      mm);
            }

            let superkmer = Superkmer {
                sequence,
                minimizer,
                mpos,
            };
            super_kmers.push(superkmer);
            start = end - k;
        }
        last_minimizer_pos = kmer_minimizer_pos;
        end += 1;
    }
    // output last superkmer
    if start <= read_len - k {
        let mut sequence = std::str::from_utf8(&read[start..read_len]).unwrap().to_string();
        let mut sequence_rc = revcomp(&sequence);
        let mut mpos = last_minimizer_pos-start;
        let mut minimizer = sequence[mpos..mpos+l].to_string();
        let mut minimizer_rc = revcomp(&minimizer);
        let mm = has_two_minimizer_occ[end-k-1];
        (sequence, sequence_rc, minimizer, minimizer_rc, mpos) = 
            normalize_mpos(sequence, sequence_rc, minimizer, minimizer_rc, mpos, l, mm);
  
        if debug 
            {
                println!("new superkmer: {} len {} minimizer {} pos {} mm {}",
                      sequence, 
                      sequence.len(),
                      minimizer,
                      mpos,
                      mm);
            }

        let superkmer = Superkmer {
            sequence,
            minimizer,
            mpos,
        };

        super_kmers.push(superkmer);
    }

    //println!("Maximum RSS at end of superkmers extraction: {:?}GB", (get_memory_rusage() as f32) / 1024.0 / 1024.0 / 1024.0);    
    super_kmers
}

