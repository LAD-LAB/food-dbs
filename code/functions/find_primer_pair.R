find_primer_pair <- function(seqs, fwd, rev, err=0.2, removePrimers = FALSE){
     ### Given seqs (a DNAStringSet) and a primer pair (DNAStrings), returns 
     ### a trimmed sequence list containing only the sequences that lie between
     ### and include the two primers and would be amplified by PCR. 
     
     ### Requires the following packages: BioStrings, tidyverse
     fwd_rc <- reverseComplement(fwd)
     rev_rc <- reverseComplement(rev)
     
     fwd_err <- err*length(fwd)
     rev_err <- err*length(rev)
     
     # Note parameter settings. 
     # fixed=FALSE accommodates degenerate nucleotides, but will match Ns
     # if they haven't been removed first
     
     # Could put in a check to make sure the length of the match isn't more than
     # error rate Ns.
     
     if (removePrimers){
          fwd_matches <- 
               vmatchPattern(fwd, 
                             seqs, 
                             max.mismatch = fwd_err,
                             fixed = TRUE) %>%
               as.data.frame() %>%
               mutate(trim_start = end + 1) %>% 
               select(group, trim_start)
          
          fwd_rc_matches <- 
               vmatchPattern(fwd_rc, 
                             seqs, 
                             max.mismatch = fwd_err,
                             fixed = TRUE) %>%
               as.data.frame() %>%
               mutate(trim_end = start - 1) %>%
               select(group, trim_end)
          
          rev_matches <- 
               vmatchPattern(rev, 
                             seqs, 
                             max.mismatch = rev_err,
                             fixed = TRUE) %>%
               as.data.frame() %>%
               mutate(trim_start = end + 1) %>%
               select(group, trim_start)
          
          rev_rc_matches <- 
               vmatchPattern(rev_rc, 
                             seqs, 
                             max.mismatch = rev_err,
                             fixed = TRUE) %>%
               as.data.frame() %>%
               mutate(trim_end = start - 1) %>%
               select(group, trim_end)
          
     } else {
          fwd_matches <- 
               vmatchPattern(fwd, 
                             seqs, 
                             max.mismatch = fwd_err,
                             fixed = TRUE) %>%
               as.data.frame() %>%
               mutate(trim_start = start) %>%
               select(group, trim_start)
          
          fwd_rc_matches <- 
               vmatchPattern(fwd_rc, 
                             seqs, 
                             max.mismatch = fwd_err,
                             fixed = TRUE) %>%
               as.data.frame() %>%
               mutate(trim_end = end) %>%
               select(group, trim_end)
          
          rev_matches <- 
               vmatchPattern(rev, 
                             seqs, 
                             max.mismatch = rev_err,
                             fixed = TRUE) %>%
               as.data.frame() %>%
               mutate(trim_start = start) %>%
               select(group, trim_start)
          
          rev_rc_matches <- 
               vmatchPattern(rev_rc, 
                             seqs, 
                             max.mismatch = rev_err,
                             fixed = TRUE) %>%
               as.data.frame() %>%
               mutate(trim_end = end) %>%
               select(group, trim_end)
     }

     
     # Take intersect of matches between primer pairs to find sequences
     # containing hits to both
     combined <-
          inner_join(fwd_matches,
                     rev_rc_matches,
                     by='group',
                     multiple = 'all',
                     relationship = 'many-to-many') %>%
          rbind(inner_join(rev_matches,
                           fwd_rc_matches,
                           by='group',
                           multiple = 'all',
                           relationship = 'many-to-many')) %>%
          # Confirm amplicon size is within reason
          mutate(size = trim_end-trim_start) %>%
          filter(size > 0 & size <1000)
     
     # This piece isn't necessary if primers are being trimmed
     # Errors can occur when there's a mismatch at the terminal bases of the
     # primer, which allow it to extend "off" the end of the sequence
     if(removePrimers == FALSE){
          combined <-
               combined %>%
               # Add length column so I can compare it to current 'end'
               mutate(len = width(seqs[combined$group])) %>%
               # Update end variable conditionally
               mutate(trim_end = case_when(trim_end < len ~ trim_end,
                                           trim_end > len ~ len)) %>%
               # Correct for 'starts' that extend off front of seq, giving negative
               # coordinates
               mutate(trim_start = case_when(trim_start > 1 ~ trim_start,
                                             trim_start < 1 ~ as.integer(1))) # RHS has to be same type of vector
     }
          
     # Trim input sequences based on paired primer coordinates
     # Need to figure out how to deal with duplicates??
     matches <- mapply(subseq, 
                       seqs[combined$group], 
                       combined$trim_start, combined$trim_end) 
     
     DNAStringSet(matches)
}

# # This is helpful for figuring out errors that I don't have a hypothesis about
# # what they are, so I can't find appropriate filter
# for (index in combined$group){
#      cat(index)
#      match <- subseq(seqs[index], 
#                      combined[combined$group==index, 'start'],
#                      combined[combined$group==index, 'end'])
# }