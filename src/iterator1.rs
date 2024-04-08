use nthash::NtHashIterator;
use crate::utils::revcomp;
use std::collections::VecDeque;
use colored::Colorize;
 

#[derive(Debug)]
pub struct Superkmer {
    pub sequence: String,
    pub minimizer: String,
    pub mpos: usize,
}

pub struct SuperkmersIterator<'a> {
    read: &'a [u8],
    k: usize,
    l: usize,
    nthash_iterator: NtHashIterator<'a>,
    start: usize,
    end: usize,
    i: usize,
    dq: VecDeque<usize>,
    buffer: Vec<usize>,
    done: bool,
    mm: Option<usize>, // whether the previous k-mer had multiple occurrences of the minimizer
}
/* convention for k-mers having multiple minimizers, possibly in fw or rc orientation:
 * they're currently in their own superkmer, quarantined!
 * and the normalization is done as follows:
 * consider the canonical representation of the minimizer (ie mseq, where mseq <= rc(mseq))
 * then the k-mer is recorded as the orientation where mseq appears and the position
 * of the minimizer is the leftmost one
 * */
impl<'a> SuperkmersIterator<'a> {
    pub fn new(read: &'a [u8], k: usize, l: usize) -> SuperkmersIterator<'a> {
        //println!("getting superkmers for read {}",std::str::from_utf8(read).unwrap());
        let mut nthash_iterator = NtHashIterator::new(read, l).unwrap();
        let mut dq: VecDeque<usize> = VecDeque::new(); // holds position of the smallest hash in front
        let mut buffer = vec![usize::max_value(); k];
		let mut mm = None;
        let mut minim = usize::max_value();

        for j in 0..k-l+1 {
            let hash = nthash_iterator.next().unwrap() as usize;
			//println!("prelim dq {:?} buffer {:?} hash {:?}",dq,
            //buffer.clone().into_iter().map(|c| if c<usize::max_value() { c } else { 0 } ).collect::<Vec<usize>>(),hash);
            // Remove elements from the back of the deque that are greater than the current element
            while !dq.is_empty() && buffer[dq[dq.len() - 1] % k] >= hash {
                dq.pop_back();
            }

            // Add the current index to the back of the deque
            dq.push_back(j);

            // Update the buffer with the current element
            buffer[j] = hash;

            // update the multiple minimizer signal
            if hash == minim {
                mm = Some(hash);
            }
            if hash < minim { 
                minim = hash;
                mm = None;
            }
       }
 
        SuperkmersIterator {
            read,
            k,
            l,
            nthash_iterator,
            start: 0,
            end: k,
            i: k-l+1,
            dq,
            buffer,
            done: false,
            mm,
        }
    }
}

impl<'a> Iterator for SuperkmersIterator<'a> {
    type Item = Superkmer;
    /* invariant: we're at a k-mer that's either the first one in the sequence,
     * or one that has a minimizer different from the one previously.
     * The following hash is the last l-mer of the k-mer*/
    fn next(&mut self) -> Option<Superkmer> {
        let verbose = false;
        // FIXME: there is one final issue i haven't solved yet
        // RUST_BACKTRACE=1 bash run.sh test/minimizer_twice_in_kmer4.fa  --threads 16  -l 6 -k 9 --stats
        // it's the case of two consecutive k-mers inside a superkmer, where one is the rc of the
        // other: AGAGCTCTT - AAGAGCTCT. in this case we should break the superkmer and normalize,
        // but how to detect that without computing the kmers or the minimizers explicitly?
        // I think it can be fixed in post, when doing queries, check for both the forward and the
        // reverse. All it does is add redundancy to the structure
        if self.done { return None; }
		let mut mpos = self.dq[0]-self.start;
		assert!(mpos <= self.k-self.l);
        let mut next_mm = self.mm;
        loop {
            // get the rightmost l-mer of the current k-mer (position self.i)
            let hash = self.nthash_iterator.next();
            if verbose {
        		println!("dq {:?} buffer {:?} hash {:?}",self.dq,self.
                         buffer.clone().into_iter().map(|c| if c<usize::max_value() { c/10000000000000000 } else { 0 } )
                         .collect::<Vec<usize>>(), hash.unwrap_or(0)/10000000000000000);
            }
            if hash.is_none() { 
                self.done = true;
                self.end += 1;
                break;
            }
            
            let hash = hash.unwrap() as usize;

            // Remove minimizers that fall off of the window
            let minimizer_out_of_scope = !self.dq.is_empty() && 
                ( self.dq[0] - self.start > self.k-self.l ||
                  self.end - self.dq[0] >= self.k );
            if minimizer_out_of_scope {
                self.dq.pop_front();
            }

            // Remove elements from the back of the deque that are greater than the current element
            while !self.dq.is_empty() && self.buffer[self.dq[self.dq.len() - 1] % self.k] > hash {
                self.dq.pop_back();
            }

            // Testing if multiple minimizers are present
            if self.dq.len() > 0 && self.buffer[self.dq[0] % self.k] == hash {
                next_mm = Some(hash);
                if verbose { println!("{}: {}",
                                  "multiple minimizers detected".red(),std::str::from_utf8(&self.read[self.end-self.k..self.end]).unwrap()); } 
            }
        
            let new_minimizer = self.dq.len() == 0;

            // Add the current index to the back of the deque
            self.dq.push_back(self.i);

            // Maybe reset multiple minimizers if new hash
            if self.mm.is_some() && self.dq.len() > 0 && self.buffer[self.dq[0] % self.k] != self.mm.unwrap() {
                next_mm = None;
            }

            // Update the buffer with the current element
            self.buffer[self.i % self.k] = hash;

            self.end += 1;
            self.i += 1;

            if minimizer_out_of_scope || new_minimizer || self.mm.is_some() || next_mm.is_some() {break;}
        }
        
        /*let superkmer = Superkmer {
            sequence: String::new(),
            minimizer: String::new(),
            mpos: 0,
        };
        self.start = self.end - self.k;
        return Some(superkmer);*/

        let mut sequence = std::str::from_utf8(&self.read[self.start..self.end-1]).unwrap().to_string();
        let mut sequence_rc = revcomp(&sequence);
        if verbose { println!("new superkmer {} mpos {}",sequence,mpos); }
        let minimizer =  sequence[mpos..mpos+self.l].to_string();
        let minimizer_rc = revcomp(&minimizer);
        if self.mm.is_some()  || minimizer == minimizer_rc { // in the double minimizer case, or if minimizer is its own rc, try to minimize mpos
            if self.mm.is_some() {
                if verbose { println!("{}: {}, minimizer {} dq {:?} buffer {:?}",
                                      "multiple minimizers in kmer".red(),
                                      sequence,minimizer,self.dq,self.buffer.clone().into_iter()
                                      .map(|c| if c<usize::max_value() { c } else { 0 } ).collect::<Vec<usize>>());}
            }
            else {
                if verbose { println!("{}: {}, minimizer {} minimizer_rc {}",
                                      "minimizer is its own rc".red(),sequence,minimizer,minimizer_rc);
                }
            }
            for i in 0..mpos {
                // if multiple minimizers, normalize mpos
                let found_forward = sequence[i..i+self.l] == minimizer  || sequence[i..i+self.l]  == minimizer_rc;
                let found_rc      = sequence_rc[i..i+self.l] == minimizer  || sequence_rc[i..i+self.l]  == minimizer_rc;
                if found_forward || found_rc {
                    if found_forward && found_rc { } // do we do anything in this case?
                    #[allow(unused_assignments)]
                    if found_rc {
                        (sequence, sequence_rc) = (sequence_rc, sequence);
                    }
                    mpos = i;
                    break;
                }
            }
        }
        self.mm = next_mm;
        let minimizer =  sequence[mpos..mpos+self.l].to_string();
        let minimizer_rc = revcomp(&minimizer);
        let rc = minimizer_rc < minimizer;
        sequence = if rc { revcomp(&sequence) } else { sequence };
        let mpos = if rc { sequence.len()-(mpos+self.l) } else { mpos };
        let minimizer = if rc { minimizer_rc } else { minimizer };
        if verbose { println!("rc {} mpos {} sequence {} minimizer {}",rc, mpos,sequence,minimizer); }

        let debug = false;
        if debug 
        {
            println!("new superkmer: {} len {} minimizer {} pos {}",
                     sequence, 
                     sequence.len(),
                     minimizer,
                     mpos);
        
            let sequence_rc = revcomp(&sequence);
            let minimizer_rc = revcomp(&minimizer);
            if ! (sequence.contains(&minimizer) || sequence_rc.contains(&minimizer) || 
                  sequence.contains(&minimizer_rc) || sequence_rc.contains(&minimizer_rc)) {
               println!("huh?! neither kmer and revcomp contain minimizer {}: {} - {}, mpos {}", minimizer, sequence, sequence_rc, mpos) 
            };
        }

        let superkmer = Superkmer {
            sequence: sequence,
            minimizer,
            mpos,
        };
        self.start = self.end - self.k;
        Some(superkmer)
    }
}

