suppressMessages(library(tidyverse))

work <- "C:/Users/sab236/Box/Home Folder sab236/Private/FoodSeq Reference Database/2026 test/coverage-recheck"
fd   <- "C:/Users/sab236/Box/Home Folder sab236/Private/FoodSeq Reference Database/food-dbs"

plants_missing  <- read_csv(file.path(fd, "data/outputs/plants_missing_trnL.csv"), show_col_types = FALSE)
animals_missing <- read_csv(file.path(fd, "data/outputs/animals_missing_12SV5.csv"), show_col_types = FALSE)

pl <- read_csv(file.path(work, "plants_trnL_recheck.csv"), show_col_types = FALSE) %>%
  mutate(ncbi_hits = as.integer(ncbi_hits)) %>%
  left_join(plants_missing, by = "scientific_name")

an <- read_csv(file.path(work, "animals_12S_recheck.csv"), show_col_types = FALSE) %>%
  mutate(ncbi_hits = as.integer(ncbi_hits)) %>%
  left_join(animals_missing, by = "scientific_name")

summarise_marker <- function(df, label, marker) {
  n     <- nrow(df)
  avail <- sum(df$ncbi_hits > 0, na.rm = TRUE)
  cat(sprintf("\n=== %s (%s) ===\n", label, marker))
  cat(sprintf("Species currently missing from reference : %d\n", n))
  cat(sprintf("Now have %s data at NCBI (hits > 0)      : %d (%.1f%%)\n", marker, avail, 100*avail/n))
  cat(sprintf("Still no %s data (hits = 0)              : %d (%.1f%%)\n", marker, n-avail, 100*(n-avail)/n))
  cat("\nHit-count distribution among those with data:\n")
  b <- df %>% filter(ncbi_hits > 0) %>%
    mutate(bucket = cut(ncbi_hits, c(0,1,2,5,10,50,Inf),
                        labels = c("1","2","3-5","6-10","11-50",">50"))) %>%
    count(bucket)
  print(b, n = Inf)
  invisible(NULL)
}

summarise_marker(pl, "PLANTS", "trnL")
summarise_marker(an, "ANIMALS", "12S")

# Which sources stand to gain most?
cat("\n=== Plants now-available, by source (top 12) ===\n")
pl %>% filter(ncbi_hits > 0) %>% count(source, sort = TRUE) %>% head(12) %>% print(n = Inf)

cat("\n=== Animals now-available, by source (top 12) ===\n")
an %>% filter(ncbi_hits > 0) %>% count(source, sort = TRUE) %>% head(12) %>% print(n = Inf)

cat("\n=== Top 25 plant candidates (most trnL records) ===\n")
pl %>% filter(ncbi_hits > 0) %>% arrange(desc(ncbi_hits)) %>%
  select(scientific_name, common_name, ncbi_hits, source) %>% head(25) %>% print(n = Inf)

cat("\n=== Top 25 animal candidates (most 12S records) ===\n")
an %>% filter(ncbi_hits > 0) %>% arrange(desc(ncbi_hits)) %>%
  select(scientific_name, common_name, ncbi_hits, source) %>% head(25) %>% print(n = Inf)

# Write actionable outputs
pl %>% filter(ncbi_hits > 0) %>% arrange(desc(ncbi_hits)) %>%
  write_csv(file.path(work, "CANDIDATES_plants_trnL_available.csv"))
an %>% filter(ncbi_hits > 0) %>% arrange(desc(ncbi_hits)) %>%
  write_csv(file.path(work, "CANDIDATES_animals_12S_available.csv"))
pl %>% filter(ncbi_hits == 0) %>% write_csv(file.path(work, "STILL_missing_plants_trnL.csv"))
an %>% filter(ncbi_hits == 0) %>% write_csv(file.path(work, "STILL_missing_animals_12S.csv"))

cat("\nWrote CANDIDATES_* and STILL_missing_* CSVs to coverage-recheck/\n")
