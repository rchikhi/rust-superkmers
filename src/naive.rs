use crate::utils::revcomp;
use nthash::NtHashIterator;

#[derive(Debug)]
pub struct Superkmer {
    pub sequence: String,
    pub minimizer: String,
    pub mpos: usize,
}

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
    for i in 0..(read_len - k + 1) {
        let mut min_minimizer     = hashes[i];
        let mut min_minimizer_pos = i;
        for j in 1..k-l+1 {
            let minimizer = hashes[i+j];
            if minimizer < min_minimizer {
                min_minimizer     = minimizer;
                min_minimizer_pos = i+j;
            }
        }
        kmer_minimizers_pos[i] = min_minimizer_pos;
    }

    // Find super-kmers by scanning the read
    let mut super_kmers = vec![];
    let mut start = 0;
    let mut end = k+1;
    let mut last_minimizer_pos = kmer_minimizers_pos[0];
    while end < read_len {
        let kmer_minimizer_pos = kmer_minimizers_pos[end-k];
        if kmer_minimizer_pos != last_minimizer_pos {
            let sequence = std::str::from_utf8(&read[start..end-1]).unwrap();
            let mpos = last_minimizer_pos-start;
            let minimizer =  sequence[mpos..mpos+l].to_string();
            let minimizer_rc = revcomp(&minimizer);
            let rc = minimizer_rc < minimizer;
            let sequence = if rc { revcomp(sequence) } else { sequence.to_string() };
            let mpos = if rc { sequence.len()-(mpos+l) } else { mpos };
            let minimizer = if rc { minimizer_rc } else { minimizer };
            
            let debug = false;
            if debug 
            {
                println!("new superkmer: {} len {} minimizer {} pos {}",
                      sequence, 
                      sequence.len(),
                      minimizer,
                      mpos);
            }

            let superkmer = Superkmer {
                sequence: sequence.clone(),
                minimizer,
                mpos,
            };
            super_kmers.push(superkmer);
            start = end - k;
        }
        last_minimizer_pos = kmer_minimizer_pos;
        end += 1;
    }
    if start != end - k {
        let mpos = last_minimizer_pos-start;
        let sequence = std::str::from_utf8(&read[start..end]).unwrap();
        let minimizer =  sequence[mpos..mpos+l].to_string();
        let minimizer_rc = revcomp(&minimizer);
        let rc = minimizer_rc < minimizer;
        let sequence = if rc { revcomp(sequence) } else { sequence.to_string() };
        let mpos = if rc { sequence.len()-(mpos+l) } else { mpos };
        let minimizer = if rc { minimizer_rc } else { minimizer };

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

