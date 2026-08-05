suppressMessages(library(tidyverse))
fd  <- "C:/Users/sab236/Box/Home Folder sab236/Private/FoodSeq Reference Database/food-dbs"
scr <- "C:/Users/sab236/AppData/Local/Temp/claude/C--Users-sab236-Box-Home-Folder-sab236-Private-FoodSeq-Reference-Database/5f594644-027b-4313-89aa-f4f890b18de0/scratchpad"

foods <- read_csv(file.path(fd, "data/inputs/human-foods.csv"), show_col_types = FALSE)
norm  <- function(x) x %>% replace_na("") %>% str_to_lower() %>% str_squish()

sci_set   <- norm(foods$scientific_name)
genus_set <- word(sci_set, 1)
hay <- paste(norm(foods$scientific_name), norm(foods$common_name),
             norm(foods$alternative_names), sep = " | ")

probe <- read_delim(file.path(scr, "probe_list.txt"), delim = "|",
                    col_names = c("scientific_name", "common_name", "group"),
                    show_col_types = FALSE, trim_ws = TRUE) %>%
  filter(!is.na(scientific_name), scientific_name != "")

res <- probe %>%
  mutate(
    key           = norm(scientific_name),
    genus         = word(key, 1),
    exact_sci     = key %in% sci_set,
    genus_present = genus %in% genus_set,
    anywhere      = map_lgl(key, ~ any(str_detect(hay, fixed(.x)))),
    status = case_when(
      exact_sci     ~ "present (species)",
      anywhere      ~ "present (synonym/alt name)",
      genus_present ~ "GAP (species absent; genus present)",
      TRUE          ~ "GAP (fully absent)"
    )
  )

cat("=== PROBE: 88 additional food species ===\n")
res %>% count(status, sort = TRUE) %>% print(n = Inf)

cat("\n=== by group ===\n")
res %>% mutate(gap = str_starts(status, "GAP")) %>%
  count(group, gap) %>% pivot_wider(names_from = gap, values_from = n, values_fill = 0) %>%
  rename(present = `FALSE`, gaps = `TRUE`) %>% arrange(desc(gaps)) %>% print(n = Inf)

cat("\n=== ALL GAPS ===\n")
res %>% filter(str_starts(status, "GAP")) %>%
  select(scientific_name, common_name, group, status) %>%
  arrange(group, scientific_name) %>% print(n = Inf)

write_csv(res, file.path(scr, "probe_result.csv"))
cat("\nwrote probe_result.csv\n")
