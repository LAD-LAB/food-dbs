make_reference <- function(marker, filt=NULL){
     # Given a region of interest ('marker') formulated as query term, searches
     # NCBI's Nucleotide database for all sequences corresponding to that query 
     # combined with all organisms listed in template_panel$query.  Provided 
     # with an additional filter ('filt'), will only search for the region in 
     # organisms belonging to a particular category (should be one of 'plant', 
     # 'animal', 'fungus'). 
     
     # Returns a DNAStringSet of results with dada2-compatible names. 
     # Requires the following packages: BioStrings, rentrez, taxonomizr
     
     # Filter by template category if provided
     if (is.null(filt)){
          templates <- template_panel$query
     } else {
          templates <- filter(template_panel, category==filt) %>% pull(query)
     }
     
     # Now loop over templates to query Entrez for requested region of each
     # Initialize empty DNAStringSet for holding results
     results <- DNAStringSet() 
     for (template in templates){
          query <- paste0(template, ' AND biomol_genomic[PROP] AND ', marker)
          ids <- entrez_search(query, db='nucleotide', retmax = 10000, use_history=TRUE)
          
          if (ids$count==0){
               # Do nothing further if there aren't any hits
          } else {
               # Regular expression to pull FASTA header out from string
               # Don't match >, match immediately after up to first space
               ex <- '[^>]\\S*' 
               seqs <- 
                    entrez_fetch(db='nucleotide', 
                                 web_history = ids$web_history, 
                                 rettype='fasta') %>%
                    # This returns concatenated sequence strings; split apart 
                    # so we can re-name inline
                    strsplit('\n\n') %>%
                    unlist()
               # Save this to ultimately combine with taxonomy data, as want to
               # be able to identify these sequences after the fact
               accs <- str_extract(seqs, ex) 
               
               seqs <- 
                    seqs %>%
                    # Now update seqs to sequence only, stripping header
                    sub('^[^\n]*\n', '', .) %>%
                    # And removing separating \n characters
                    gsub('\n', '', .)
               
               # Can't figure out how to use web_history with entrez_link!!
               if (ids$count > 100){
                    chunk_size <- 100
                    ids_chunked <- split(ids$ids, ceiling(seq_along(ids$ids)/chunk_size))
                    links <- lapply(ids_chunked, 
                                    function(x) entrez_link(dbfrom='nucleotide', 
                                                            id=x, db='taxonomy', 
                                                            by_id=TRUE)) %>%
                         unlist(recursive=FALSE)
               } else {
                    links <- entrez_link(dbfrom='nucleotide', id=ids$ids, db='taxonomy', by_id=TRUE)
               }
               
               if (length(links$links)==1){
                    taxids <- links$links$nuccore_taxonomy %>% unlist()
               } else {
                    # Extract taxids
                    taxids <- lapply(links, function(x) x$links$nuccore_taxonomy) %>% 
                         unlist()
                    
               }
               
               revised_names <- data.frame(taxid = taxids, 
                                           accession = accs,
                                           seq = seqs,
                                           stringsAsFactors = F)
               
               taxonomy <- getTaxonomy(unique(revised_names$taxid), 
                                       sqlFile = '/Volumes/brianna/david/taxonomy/accessionTaxa.sql')
               
               # Pull taxid information before converting to dataframe, which 
               # changes formatting
               # Put taxonomy into dada2-compatible format (single string with 
               # taxonomic entries separated by semicolons as described in
               # https://benjjneb.github.io/dada2/training.html)
               rows <- trimws(row.names(taxonomy)) 
               
               taxonomy <- 
                    taxonomy %>% 
                    data.frame() %>% 
                    mutate(taxid=rows) %>% 
                    mutate(dadaname = apply(select(., superkingdom:species), 1, paste, collapse=';'))  
               
               # Should put statement here to catch error if dadaname somehow ends up as NA
               revised_names <- 
                    select(taxonomy, c(taxid, dadaname)) %>%
                    right_join(revised_names, by = 'taxid') %>%
                    # Add accession info as final level
                    mutate(name = paste(dadaname, accession, sep=';'))
               
               # Now add to DNAStringSet
               seqs <- DNAStringSet(revised_names$seq)
               names(seqs) <- revised_names$name
               
               results <- append(results, seqs)
          }
          cat(ids$count, 'sequences processed for', template, '\n')
     }
     results
}