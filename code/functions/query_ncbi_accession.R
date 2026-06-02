query_ncbi_accession <- function(acc, max_tries = 4, retry_wait = 5){
     # Given an accession ID ('acc'), searches NCBI's Nucleotide database and
     # returns the corresponding integer taxon ID, or NA if not found.
     #
     # Both entrez_search and entrez_fetch are wrapped in the retry loop so
     # that transient network errors at either step are handled gracefully.
     #
     # max_tries:  number of attempts before giving up and returning NA
     # retry_wait: seconds to wait between attempts

     library(rentrez)
     library(stringr)

     for (attempt in seq_len(max_tries)) {

          result <- tryCatch({

               ids <- entrez_search(acc,
                                    db          = 'nucleotide',
                                    retmax      = 1,
                                    use_history = TRUE)

               if (ids$count == 0) return(NA_integer_)

               fetch_text <- entrez_fetch(db          = 'nucleotide',
                                          web_history = ids$web_history,
                                          rettype     = 'xml',
                                          retmax      = 1)

               taxon <- str_extract(fetch_text, 'taxon:(.*)<')
               taxon <- str_replace(taxon, 'taxon:', '')
               taxon <- str_replace(taxon, '<',     '')
               as.integer(taxon)

          }, error = function(e) {
               cat("  query_ncbi_accession(", acc, ") — error (attempt",
                   attempt, "of", max_tries, "):", conditionMessage(e), "\n")
               Sys.sleep(retry_wait)
               NULL
          })

          if (!is.null(result)) return(result)
     }

     warning("query_ncbi_accession: could not resolve taxon ID for '", acc,
             "' after ", max_tries, " attempts — returning NA")
     NA_integer_
}
