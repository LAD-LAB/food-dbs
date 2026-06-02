acc_to_dada <- function(ncbi_result){
     ### Formats the output of query_ncbi.R in a DADA2-compatible format.
     ### Handles processing of repeated sequences at the level of individual
     ### species. 
     
     # Strip name to accession number
     accs <- 
          ncbi_result %>%
          names() %>%
          str_split(pattern = ' ') %>%
          lapply(sub, pattern = '>', replacement = '') %>%
          lapply('[[', 1) %>%
          unlist()
     
     # Query NCBI to link accession number to taxid
     taxmap <- taxa::lookup_tax_data(accs, 
                                     type = 'seq_id')
     
     # Split into 50-accession number chunks
     accs.50 <- split(accs, ceiling(seq_along(accs)/50))
     # Use accession number to link to taxid with entrez_link
     acc.links <- lapply(accs.50, function(x) {
          entrez_link(id=x, dbfrom='nucleotide', db='taxonomy', by_id = T)
     })
     
     # Unless by_id = T, this gives a "collapsed" set of taxids rather than ones 
     # that are repeated for each entry
     
     # Un-chunk
     # Ordering is preserved, but names are concatenated
     # Item 1 of list 1 becomes item '11' in the unlisted version
     # Add underscore to eliminate risk of overlap
     names(acc.links) <- paste0(names(acc.links), '_')
     acc.links.unlist <- unlist(acc.links, recursive = FALSE) 
     
     # Extract taxids
     accs.retreived <- 
          lapply(acc.links.unlist, 
                 function(x) x$links$nuccore_taxonomy)
     
     taxids <- 
          lapply(acc.links.unlist, 
                 function(x) x$links$nuccore_taxonomy) %>% 
          unlist() %>% 
          data.frame()
     
     # Organize
     revised_names <- data.frame(taxid = taxids, 
                                 accession = accs, 
                                 stringsAsFactors = F)
     
     taxmap <- taxa::lookup_tax_data(taxtab.species, 
                                     type = 'taxon_name', 
                                     column = 'Species')
     
     taxonomy <- taxa::taxonomy_table(taxmap, 
                                      use_ranks = c('superkingdom', 'kingdom',
                                                    'phylum', 'order', 'family',
                                                    'genus', 'species'))
     
     taxonomy <- taxonomizr::getTaxonomy(revised_names$taxid, '/Volumes/brianna/david/taxonomy/prior/accessionTaxa.sql')
     
     rows <- trimws(row.names(taxonomy)) # Pull taxid information before converting to dataframe, which changes formatting
     # Put taxonomy into dada2-compatible format (single string with taxonomic entries separated by semicolons as described in https://benjjneb.github.io/dada2/training.html)
     taxonomy <- 
          taxonomy %>% 
          data.frame() %>% 
          mutate(taxid=rows) %>% 
          mutate(dadaname = apply(select(., superkingdom:species), 1, paste, collapse=';'))  
     
     ## TODO: figure out a strategy to preserve accession number
     # Rename sequences with taxonomic name
     names(ncbi_result) <- taxonomy$dadaname
     # Remove seqs with 'NA'
     remove <- names(ncbi_result)=='NA;NA;NA;NA;NA;NA;NA'
     ncbi_result <- ncbi_result[!remove]
     
     # Split by names-- check on multiple entries for a single species
     result_species <- split(ncbi_result, names(ncbi_result))
     
     only_one <- lapply(result_species, length) == 1
     add_ncbi <- unlist(result_species[only_one])
     
     multiple <- result_species[!only_one]
     
     # Now, take this list of multiple entries and reduce it to unique 
     # entries on a species-by-species basis.
     reduced <- 
          multiple %>%
          lapply(unique) %>% # Reduce to unique seqs
          DNAStringSetList() %>% # Reformat in order to...
          unlist(use.names = FALSE) # Unlist; use.names prevents concatenation of all the names in the sub-list
     
     add_ncbi <- append(add_ncbi, reduced)
     add_ncbi
}