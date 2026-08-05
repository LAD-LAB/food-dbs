suppressMessages(library(tidyverse))
scr <- "C:/Users/sab236/AppData/Local/Temp/claude/C--Users-sab236-Box-Home-Folder-sab236-Private-FoodSeq-Reference-Database/5f594644-027b-4313-89aa-f4f890b18de0/scratchpad"

gaps <- read_csv(file.path(scr, "probe_result.csv"), show_col_types = FALSE) %>%
  filter(str_starts(status, "GAP")) %>%
  mutate(marker = if_else(group == "plant", "trnL", "12S"))

gaps %>%
  transmute(out = paste(scientific_name, marker, sep = "|")) %>%
  pull(out) %>%
  write_lines(file.path(scr, "gap_query.txt"))

cat("wrote", nrow(gaps), "gap species to gap_query.txt\n")
