use nthash::NtHashIterator;
use crate::utils::{revcomp, normalize_mpos};
use std::collections::VecDeque;
use colored::Colorize;
use crate::Superkmer; 

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

        for j in 0..k-l+1 {
            let hash = nthash_iterator.next().unwrap() as usize;
			//println!("prelim dq {:?} buffer {:?} hash {:?}",dq,
            //buffer.clone().into_iter().map(|c| if c<usize::max_value() { c } else { 0 } ).collect::<Vec<usize>>(),hash);
            // Remove elements from the back of the deque that are greater than the current element
            while !dq.is_empty() && buffer[dq[dq.len() - 1] % k] > hash {
                dq.pop_back();
            }

            // Add the current index to the back of the deque
            dq.push_back(j);

            // Update the buffer with the current element
            buffer[j] = hash;
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
        if self.done { return None; }

        let mut mpos = self.dq[0]-self.start;
		assert!(mpos <= self.k-self.l);
        let mm = self.dq.len() > 1 && self.buffer[self.dq[0] % self.k ] == self.buffer[self.dq[1] % self.k];

        // given a starting k-mer, a minimizer pos mpos, and the fact that minimizer is maybe seen
        // twice in starting kmer (but if seen again later, it will start another superkmer),
        // extend superkmer as far as possible
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

            if verbose
            {
                println!("new hash {} mmer {}",hash/10000000000000000, 
                         std::str::from_utf8(&self.read[self.i..self.i+self.l]).unwrap().to_string());
            }
 
            // Remove elements from the back of the deque that are greater than the current element
            while !self.dq.is_empty() && self.buffer[self.dq[self.dq.len() - 1] % self.k] > hash {
                self.dq.pop_back();
            }

            let new_minimizer = self.dq.len() == 0;

            // Add the current index to the back of the deque
            self.dq.push_back(self.i);

            // Update the buffer with the current element
            self.buffer[self.i % self.k] = hash;

            self.end += 1;
            self.i += 1;

            let seen_minimizer_again = hash == self.buffer[self.dq[0] % self.k];
        
            if minimizer_out_of_scope || new_minimizer || mm || seen_minimizer_again {break;}
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
        let mut minimizer =  sequence[mpos..mpos+self.l].to_string();
        let mut minimizer_rc = revcomp(&minimizer);
        if verbose 
        {
            if mm {
                 println!("{}: {}, minimizer {} dq {:?} buffer {:?}",
                                      "multiple minimizers in kmer".red(),
                                      sequence,minimizer,self.dq,self.buffer.clone().into_iter()
                                      .map(|c| if c<usize::max_value() { c } else { 0 } ).collect::<Vec<usize>>());
            }
            if minimizer == minimizer_rc {
                 println!("{}: {}, minimizer {} minimizer_rc {}",
                                      "minimizer is its own rc".red(),sequence,minimizer,minimizer_rc);
            }
        }
        (sequence, sequence_rc, minimizer, minimizer_rc, mpos) = 
            normalize_mpos(sequence, sequence_rc, minimizer, minimizer_rc, mpos, self.l, mm);

        if verbose
        {
            println!("new superkmer: {} len {} minimizer {} pos {} mm {}",
                     sequence, 
                     sequence.len(),
                     minimizer,
                     mpos,
                     mm
                     );
        
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

