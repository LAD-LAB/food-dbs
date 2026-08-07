suppressMessages(library(tidyverse))
work <- "C:/Users/sab236/Box/Home Folder sab236/Private/FoodSeq Reference Database/2026 test/coverage-recheck"

classify <- function(path, marker) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      ncbi_hits = as.integer(ncbi_hits),
      n_words   = str_count(str_trim(scientific_name), "\\s+") + 1,
      rank_type = if_else(n_words >= 2, "species-level (binomial)", "higher rank (genus/family/order)")
    )
}

pl <- classify(file.path(work, "plants_trnL_recheck.csv"),  "trnL")
an <- classify(file.path(work, "animals_12S_recheck.csv"), "12S")

report <- function(df, label, marker) {
  cat(sprintf("\n=== %s (%s) ===\n", label, marker))
  tab <- df %>%
    group_by(rank_type) %>%
    summarise(
      total      = n(),
      with_data  = sum(ncbi_hits > 0, na.rm = TRUE),
      no_data    = sum(ncbi_hits == 0, na.rm = TRUE),
      pct_avail  = round(100 * with_data / n(), 1),
      .groups    = "drop"
    )
  print(tab, n = Inf)

  sp <- df %>% filter(rank_type == "species-level (binomial)")
  cat(sprintf("\n  -> SPECIES-LEVEL ONLY: %d of %d (%.1f%%) have %s data\n",
              sum(sp$ncbi_hits > 0, na.rm = TRUE), nrow(sp),
              100*sum(sp$ncbi_hits > 0, na.rm = TRUE)/nrow(sp), marker))
  invisible(NULL)
}

report(pl, "PLANTS", "trnL")
report(an, "ANIMALS", "12S")

cat("\n=== Higher-rank entries with data (these need different handling) ===\n")
bind_rows(
  pl %>% filter(rank_type != "species-level (binomial)", ncbi_hits > 0) %>%
    transmute(marker = "trnL", scientific_name, ncbi_hits),
  an %>% filter(rank_type != "species-level (binomial)", ncbi_hits > 0) %>%
    transmute(marker = "12S", scientific_name, ncbi_hits)
) %>% arrange(desc(ncbi_hits)) %>% print(n = 40)

cat("\n=== Top 20 TRUE SPECIES-LEVEL plant candidates ===\n")
pl %>% filter(rank_type == "species-level (binomial)", ncbi_hits > 0) %>%
  arrange(desc(ncbi_hits)) %>% select(scientific_name, ncbi_hits) %>% head(20) %>% print(n = Inf)

cat("\n=== Top 20 TRUE SPECIES-LEVEL animal candidates ===\n")
an %>% filter(rank_type == "species-level (binomial)", ncbi_hits > 0) %>%
  arrange(desc(ncbi_hits)) %>% select(scientific_name, ncbi_hits) %>% head(20) %>% print(n = Inf)

# rewrite candidate files split by rank type
pl %>% filter(rank_type == "species-level (binomial)", ncbi_hits > 0) %>% arrange(desc(ncbi_hits)) %>%
  write_csv(file.path(work, "CANDIDATES_plants_trnL_available.csv"))
an %>% filter(rank_type == "species-level (binomial)", ncbi_hits > 0) %>% arrange(desc(ncbi_hits)) %>%
  write_csv(file.path(work, "CANDIDATES_animals_12S_available.csv"))
bind_rows(pl, an) %>% filter(rank_type != "species-level (binomial)") %>%
  write_csv(file.path(work, "REVIEW_higher-rank-entries.csv"))
cat("\nRewrote CANDIDATES_* (species-level only) + REVIEW_higher-rank-entries.csv\n")
