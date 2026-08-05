suppressMessages(library(tidyverse))
work <- "C:/Users/sab236/Box/Home Folder sab236/Private/FoodSeq Reference Database/2026 test/coverage-recheck"
fd   <- "C:/Users/sab236/Box/Home Folder sab236/Private/FoodSeq Reference Database/food-dbs"

want <- read_lines(file.path(work, "vert_input.txt")) %>% str_trim() %>% discard(~ .x == "")

# raw vert output has no header and may contain malformed rows; parse defensively
raw <- read_lines(file.path(work, "animals_vertebrate_check.csv.bak")) %>% discard(~ .x == "")
parsed <- tibble(line = raw) %>%
  extract(line, into = c("scientific_name", "vert_12S_hits"),
          regex = '^"(.*)",([0-9]+)$', remove = FALSE, convert = FALSE) %>%
  filter(!is.na(scientific_name)) %>%
  mutate(vert_12S_hits = as.integer(vert_12S_hits)) %>%
  # keep only rows whose species exactly matches an input species (drops malformed + bogus)
  filter(scientific_name %in% want) %>%
  distinct(scientific_name, .keep_all = TRUE) %>%
  select(scientific_name, vert_12S_hits)

stopifnot(nrow(parsed) == length(want))
cat(sprintf("Parsed %d/%d species cleanly.\n", nrow(parsed), length(want)))

phase1  <- read_csv(file.path(work, "animals_12S_recheck.csv"), show_col_types = FALSE,
                    col_types = cols(.default = col_character())) %>%
  mutate(total_12S_hits = as.integer(ncbi_hits)) %>% select(scientific_name, total_12S_hits)
missing <- read_csv(file.path(fd, "data/outputs/animals_missing_12SV5.csv"), show_col_types = FALSE)

final <- parsed %>%
  left_join(phase1,  by = "scientific_name") %>%
  left_join(missing, by = "scientific_name") %>%
  mutate(is_vertebrate = vert_12S_hits > 0) %>%
  select(scientific_name, common_name, source, total_12S_hits, vert_12S_hits, is_vertebrate)

n_v  <- sum(final$is_vertebrate)
n_iv <- sum(!final$is_vertebrate)
cat(sprintf("\n=== 12SV5 candidate list, vertebrate-filtered ===\n"))
cat(sprintf("Animal candidates from phase 1        : %d\n", nrow(final)))
cat(sprintf("TRUE VERTEBRATES (actionable)         : %d (%.1f%%)\n", n_v, 100*n_v/nrow(final)))
cat(sprintf("Invertebrates (12SV5 won't amplify)   : %d (%.1f%%)\n", n_iv, 100*n_iv/nrow(final)))

cat("\n--- Top 20 actionable vertebrate candidates ---\n")
final %>% filter(is_vertebrate) %>% arrange(desc(vert_12S_hits)) %>%
  select(scientific_name, common_name, vert_12S_hits) %>% head(20) %>% print(n = Inf)

cat("\n--- Invertebrates being excluded, by source ---\n")
final %>% filter(!is_vertebrate) %>% count(source, sort = TRUE) %>% head(10) %>% print(n = Inf)

cat("\n--- Sample of excluded invertebrates ---\n")
final %>% filter(!is_vertebrate) %>% arrange(desc(total_12S_hits)) %>%
  select(scientific_name, common_name, total_12S_hits) %>% head(15) %>% print(n = Inf)

final %>% filter(is_vertebrate)  %>% arrange(desc(vert_12S_hits)) %>%
  write_csv(file.path(work, "CANDIDATES_animals_12SV5_VERTEBRATES.csv"))
final %>% filter(!is_vertebrate) %>% arrange(desc(total_12S_hits)) %>%
  write_csv(file.path(work, "EXCLUDED_animals_invertebrates.csv"))
cat("\nWrote CANDIDATES_animals_12SV5_VERTEBRATES.csv and EXCLUDED_animals_invertebrates.csv\n")
