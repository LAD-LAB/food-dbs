fetch_lineage_ncbi <- function(taxids, batch_size = 200, max_tries = 4,
                               retry_wait = 5, verbose = TRUE){
     ### Given a vector of NCBI taxon IDs, returns a data frame with one row per
     ### input taxid and one column per taxonomic rank, in the 10-rank order the
     ### DADA2-compatible taxonomy FASTAs use:
     ###   superkingdom, phylum, class, order, family, genus, species,
     ###   subspecies, varietas, forma
     ###
     ### This is the SQL-free counterpart to extract_taxonomy(): it reads
     ### lineages from NCBI's taxonomy database over e-utils rather than from a
     ### local accessionTaxa.sql. That makes it usable on a laptop, which is the
     ### whole point of the extension workflow, at the cost of needing network
     ### access and being far slower per record. Use extract_taxonomy() for
     ### whole-database builds; use this for adding a handful of species.
     ###
     ### Rows are returned in the order of `taxids`, with NA for any taxid NCBI
     ### does not resolve. Duplicate taxids are fetched once and mapped back.

     library(rentrez)
     library(xml2)

     RANKS <- c("superkingdom", "phylum", "class", "order", "family",
                "genus", "species", "subspecies", "varietas", "forma")

     empty_row <- function() setNames(rep(NA_character_, length(RANKS)), RANKS)

     out <- as.data.frame(matrix(NA_character_, nrow = length(taxids),
                                 ncol = length(RANKS),
                                 dimnames = list(NULL, RANKS)),
                          stringsAsFactors = FALSE)
     out$taxid <- taxids
     out <- out[, c("taxid", RANKS)]

     keep <- !is.na(taxids)
     uniq <- unique(taxids[keep])
     if (length(uniq) == 0) return(out)

     fetched <- list()
     chunks  <- split(uniq, ceiling(seq_along(uniq) / batch_size))

     for (i in seq_along(chunks)) {
          xml <- NULL
          for (attempt in seq_len(max_tries)) {
               xml <- tryCatch(
                    entrez_fetch(db = "taxonomy", id = chunks[[i]], rettype = "xml"),
                    error = function(e) {
                         cat("  fetch_lineage_ncbi: batch", i, "error (attempt",
                             attempt, "of", max_tries, "):", conditionMessage(e), "\n")
                         Sys.sleep(retry_wait)
                         NULL
                    })
               if (!is.null(xml)) break
          }
          if (is.null(xml)) {
               warning("fetch_lineage_ncbi: batch ", i, " failed after ", max_tries,
                       " attempts; those taxids stay NA.")
               next
          }

          doc <- read_xml(xml)
          for (tx in xml_find_all(doc, "/TaxaSet/Taxon")) {
               id  <- xml_text(xml_find_first(tx, "./TaxId"))
               row <- empty_row()

               # Ancestors, each carrying its own rank
               anc <- xml_find_all(tx, "./LineageEx/Taxon")
               a_rank <- xml_text(xml_find_all(anc, "./Rank"))
               a_name <- xml_text(xml_find_all(anc, "./ScientificName"))

               # NCBI renamed the top-level rank "superkingdom" -> "domain"
               # (see the same fix in extract_taxonomy()). Accept either and
               # store it as superkingdom so the output schema is stable.
               a_rank[a_rank == "domain"] <- "superkingdom"

               hit <- a_rank %in% RANKS
               row[a_rank[hit]] <- a_name[hit]

               # The queried taxon itself is not part of LineageEx
               own_rank <- xml_text(xml_find_first(tx, "./Rank"))
               own_name <- xml_text(xml_find_first(tx, "./ScientificName"))
               if (own_rank == "domain") own_rank <- "superkingdom"
               if (own_rank %in% RANKS) row[[own_rank]] <- own_name

               fetched[[id]] <- row
          }

          if (verbose)
               cat(sprintf("  lineages: %d of %d\n",
                           min(i * batch_size, length(uniq)), length(uniq)))
          Sys.sleep(0.15)   # respect NCBI rate limits
     }

     for (r in RANKS) {
          out[[r]] <- vapply(as.character(taxids), function(t) {
               if (is.na(t) || is.null(fetched[[t]])) NA_character_ else fetched[[t]][[r]]
          }, character(1), USE.NAMES = FALSE)
     }

     out
}
