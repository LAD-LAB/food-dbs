query_ncbi_plastome <- function(organisms, max_per_species = 10,
                                min_len = 50001, max_len = 400000,
                                max_tries = 4, retry_wait = 10){
     # Finds complete (or near-complete) chloroplast genomes for each organism
     # and returns them as a DNAStringSet, one record per genome, with the full
     # FASTA header preserved -- the same shape query_ncbi() returns, so
     # find_primer_pair() and the rest of the extension can take it unchanged.
     #
     # Why this exists: query_ncbi() adds a free-text marker term ("trnL") and
     # a 50 kb length cap. The cap keeps chromosomal assemblies out, but it also
     # removes every complete plastome (~150 kb), and the marker term misses
     # about 96% of them because plastome titles ("<species> chloroplast,
     # complete genome") never name the gene. A plant whose only trnL sequence
     # is inside a GenBank plastome is therefore invisible to both. This
     # function searches by organism and plastome title instead, and fetches one
     # record per request so a single 150 kb genome cannot push a combined
     # payload past R's string limit.
     #
     # max_per_species: cap on plastomes fetched per organism. Records are not
     #                  ranked, so a high cap only costs download time; the
     #                  caller de-duplicates identical amplicons afterwards.
     # min_len/max_len: SLEN window. The floor is the marker path's cap, so the
     #                  two paths do not overlap.

     library(Biostrings)
     library(rentrez)

     retry <- function(f, label) {
          for (attempt in seq_len(max_tries)) {
               result <- tryCatch(f(), error = function(e) {
                    cat("  ", label, "- error (attempt", attempt, "of", max_tries, "):",
                        conditionMessage(e), "\n")
                    Sys.sleep(retry_wait)
                    NULL
               })
               if (!is.null(result)) return(result)
          }
          warning(label, ": failed after ", max_tries, " attempts - skipping.")
          NULL
     }

     out <- list()
     n_found <- integer(length(organisms)); names(n_found) <- organisms

     for (k in seq_along(organisms)) {
          org   <- organisms[[k]]
          # sprintf, not paste0: paste0() renders 400000 as "4e+05", which NCBI
          # reads as a malformed length range and answers with zero hits.
          term  <- paste0('"', org, '"[ORGN] AND (chloroplast[TITL] OR plastid[TITL]) AND ',
                          sprintf('%d:%d[SLEN]', as.integer(min_len), as.integer(max_len)))
          srch  <- retry(function() entrez_search(db = 'nucleotide', term = term,
                                                  retmax = max_per_species),
                         paste('search', org))
          Sys.sleep(0.4)
          if (is.null(srch) || length(srch$ids) == 0) next
          n_found[[org]] <- srch$count

          for (id in srch$ids) {
               fa <- retry(function() entrez_fetch(db = 'nucleotide', id = id,
                                                   rettype = 'fasta'),
                           paste('fetch', id))
               Sys.sleep(0.4)
               if (is.null(fa)) next
               header <- sub('^>', '', sub('\n.*$', '', fa))
               body   <- gsub('[\n\r ]', '', sub('^[^\n]*\n', '', fa))
               s      <- DNAStringSet(body)
               names(s) <- header
               out[[length(out) + 1]] <- s
          }
          if (k %% 10 == 0)
               cat('Organism', k, 'of', length(organisms), '|', length(out), 'plastomes so far\n')
     }

     cat(sum(n_found), 'plastome records found for', sum(n_found > 0), 'of',
         length(organisms), 'organisms;', length(out), 'fetched\n')
     if (length(out) == 0) return(DNAStringSet())
     do.call(c, out)
}
