# Trim RefSeq mitochondrial database to 12SV5 region

args <- commandArgs(trailingOnly=TRUE)
print(args)

setwd(args[2]) 

library(Biostrings)
library(tidyverse)

# Load parsed RefSeq sequences -------------------------------------------------

seqs <- readDNAStringSet('~/Box/project_davidlab/LAD_LAB_Personnel/Brianna P/diet/reference/human-diets/dada2-compatible/refseq_mito.fasta')
     
# 12SV5
fwd <- DNAString('TAGAACAGGCTCCTCTAG')
rev <- DNAString('TTAGATACCCCACTATGC')
err <- 0.2

fwd_rc <- reverseComplement(fwd)
rev_rc <- reverseComplement(rev)

fwd_err <- err*length(fwd)
rev_err <- err*length(rev)

# Note parameter settings. 
# No generate nucleotides in 12SV5 primers, so leaving default fixed = TRUE.
# Being careful about choosing "start" and "end" columns here to exclude
# primer sequence from being retained.  This is important because primers
# are ultimately trimmed off the reads in DADA2 step
fwd_matches <- vmatchPattern(fwd, seqs, max.mismatch = fwd_err) %>%
     as.data.frame() %>%
     mutate(trim_start = end + 1) %>%
     select(group, trim_start)
fwd_rc_matches <- vmatchPattern(fwd_rc, seqs, max.mismatch = fwd_err) %>%
     as.data.frame() %>%
     mutate(trim_end = start - 1) %>%
     select(group, trim_end)
rev_matches <- vmatchPattern(rev, seqs, max.mismatch = rev_err) %>%
     as.data.frame() %>%
     mutate(trim_start = end + 1) %>%
     select(group, trim_start)
rev_rc_matches <- vmatchPattern(rev_rc, seqs, max.mismatch = rev_err) %>%
     as.data.frame() %>%
     mutate(trim_end = start - 1) %>%
     select(group, trim_end)

# TODO: Should I build in checks for multiple matches within the same sequence?

# Take intersect of matches between primer pairs to find sequences
# containing hits to both
combined <- 
     inner_join(fwd_matches, rev_rc_matches, by='group') %>%
     rbind(inner_join(rev_matches, fwd_rc_matches, by='group')) %>%
     # Confirm amplicon size is within reason
     mutate(size=trim_end-trim_start) %>%
     filter(size > 0 & size <1000)

head(combined)

# Trim input sequences based on paired primer coordinates
# Need to figure out how to deal with duplicates??
matches <- mapply(subseq, 
                  seqs[combined$group], 
                  combined$trim_start, combined$trim_end) 

writeXStringSet(DNAStringSet(matches), filepath = '~/Box/project_davidlab/LAD_LAB_Personnel/Brianna P/diet/reference/human-diets/dada2-compatible/12S/12SV5_RefSeqMito_trimmed.fasta')
