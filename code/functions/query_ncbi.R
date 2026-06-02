query_ncbi <- function(marker, organisms, chunk_size = 5,
                       max_seq_len = 50000, retmax_fetch = 500,
                       max_tries = 4, retry_wait = 10){
     # Given a region of interest ('marker') formulated as query term, searches
     # NCBI's Nucleotide database for all sequences corresponding to that query
     # combined with all organisms listed in 'organisms'.

     # Returns a DNAStringSet of results with full header preserved for
     # any downstream manipulation/renaming.

     # chunk_size:    number of species per NCBI query
     # max_seq_len:   upper sequence-length filter (SLEN); excludes chromosomal /
     #                whole-genome assemblies that can exceed R's 2^31-1 byte
     #                string limit even at small chunk sizes
     # retmax_fetch:  max sequences downloaded per entrez_fetch call; a second
     #                guard against very large payloads
     # max_tries:     number of times to retry a failed entrez_search or entrez_fetch
     # retry_wait:    seconds to wait between retries

     library(Biostrings)
     library(rentrez)

     # Helper: search NCBI with automatic retry on 502 / network errors
     search_with_retry <- function(query, chunk_label) {
          for (attempt in seq_len(max_tries)) {
               result <- tryCatch(
                    entrez_search(query,
                                  db          = 'nucleotide',
                                  retmax      = 10000,
                                  use_history = TRUE),
                    error = function(e) {
                         cat("  Chunk", chunk_label, "— search error (attempt", attempt,
                             "of", max_tries, "):", conditionMessage(e), "\n")
                         Sys.sleep(retry_wait)
                         NULL
                    }
               )
               if (!is.null(result)) return(result)
          }
          warning("Chunk ", chunk_label, ": entrez_search failed after ", max_tries,
                  " attempts — skipping this chunk.")
          NULL
     }

     # Helper: fetch FASTA from NCBI history with automatic retry on SSL / network errors
     fetch_with_retry <- function(web_history, chunk_label) {
          for (attempt in seq_len(max_tries)) {
               result <- tryCatch(
                    entrez_fetch(db          = 'nucleotide',
                                 web_history = web_history,
                                 rettype     = 'fasta',
                                 retmax      = retmax_fetch),
                    error = function(e) {
                         cat("  Chunk", chunk_label, "— fetch error (attempt", attempt,
                             "of", max_tries, "):", conditionMessage(e), "\n")
                         Sys.sleep(retry_wait)
                         NULL
                    }
               )
               if (!is.null(result)) return(result)
          }
          warning("Chunk ", chunk_label, ": entrez_fetch failed after ", max_tries,
                  " attempts — skipping this chunk.")
          NULL
     }

     # Chunk organisms to keep each fetch payload manageable
     organisms <- split(organisms, ceiling(seq_along(organisms) / chunk_size))

     organism_query <-
          organisms %>%
          sapply(dQuote) %>%
          sapply(paste0, '[ORGN]') %>%
          sapply(paste, collapse = ' OR ') %>%
          sapply(function(x){ paste0('(', x, ')') })

     # Build query: marker AND organisms AND sequence-length cap.
     # The SLEN filter excludes chromosomal / scaffold assemblies that contain
     # the marker gene but are far too large to download and process efficiently.
     slen_filter <- paste0('0:', max_seq_len, '[SLEN]')
     query <- paste(marker, organism_query, slen_filter, sep = ' AND ')

     seqs.return = NULL
     seqs.count  = 0

     for (i in seq_along(query)) {

          ids <- search_with_retry(query[[i]], chunk_label = i)
          if (is.null(ids) || ids$count == 0) next

          ex  <- '[^>]\\S*'   # regex: accession from FASTA header

          raw <- fetch_with_retry(ids$web_history, chunk_label = i)
          if (is.null(raw)) next   # skip chunk if all retries failed

          seqs <- raw %>%
               strsplit('\n{2,}') %>%
               unlist()

          accs    <- str_extract(seqs, ex)
          headers <- str_extract(seqs, '^[^\n]*')

          seqs <- seqs %>%
               sub('^[^\n]*\n', '', .) %>%
               gsub('\n', '', .)

          seqs        <- DNAStringSet(seqs)
          names(seqs) <- headers

          seqs.count  <- seqs.count + ids$count
          seqs.return <- append(seqs.return, seqs)

          cat('Chunk', i, 'of', length(query), 'processed |',
              length(seqs.return), 'sequences so far\n')

          Sys.sleep(0.5)   # respect NCBI rate limits between chunks
     }

     cat(seqs.count, 'sequences processed for', marker, '\n')
     seqs.return
}
