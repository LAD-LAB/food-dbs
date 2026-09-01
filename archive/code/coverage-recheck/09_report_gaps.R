suppressMessages(library(tidyverse))
scr  <- "C:/Users/sab236/AppData/Local/Temp/claude/C--Users-sab236-Box-Home-Folder-sab236-Private-FoodSeq-Reference-Database/5f594644-027b-4313-89aa-f4f890b18de0/scratchpad"
work <- "C:/Users/sab236/Box/Home Folder sab236/Private/FoodSeq Reference Database/2026 test/coverage-recheck"

probe <- read_csv(file.path(scr, "probe_result.csv"), show_col_types = FALSE) %>%
  filter(str_starts(status, "GAP")) %>%
  mutate(marker = if_else(group == "plant", "trnL", "12S"))
ncbi <- read_csv(file.path(scr, "gap_ncbi.csv"), show_col_types = FALSE,
                 col_types = cols(.default = col_character())) %>%
  mutate(hits = as.integer(hits))

final <- probe %>%
  left_join(ncbi, by = c("scientific_name", "marker")) %>%
  mutate(sequence_available = hits > 0) %>%
  select(scientific_name, common_name, group, marker, status, hits, sequence_available) %>%
  arrange(desc(sequence_available), group, desc(hits))

cat("=== PHASE 2: foods absent from human-foods.csv, and whether marker data exists ===\n")
cat(sprintf("Gap species identified          : %d\n", nrow(final)))
cat(sprintf("Marker sequence data AVAILABLE  : %d\n", sum(final$sequence_available, na.rm = TRUE)))
cat(sprintf("No marker data                  : %d\n", sum(!final$sequence_available, na.rm = TRUE)))

cat("\n=== ADDABLE NOW: absent from list AND sequence data exists ===\n")
final %>% filter(sequence_available) %>%
  select(scientific_name, common_name, group, hits) %>% print(n = Inf)

cat("\n=== gaps with NO marker data ===\n")
nod <- final %>% filter(!sequence_available)
if (nrow(nod) == 0) cat("  (none)\n") else print(select(nod, scientific_name, common_name, group), n = Inf)

cat("\n=== summary by group ===\n")
final %>% group_by(group) %>%
  summarise(gaps = n(), with_data = sum(sequence_available, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(gaps)) %>% print(n = Inf)

write_csv(final, file.path(work, "PHASE2_uncatalogued_foods_with_sequence_data.csv"))
cat("\nwrote PHASE2_uncatalogued_foods_with_sequence_data.csv\n")
