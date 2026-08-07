suppressMessages(library(tidyverse))

fd  <- "C:/Users/sab236/Box/Home Folder sab236/Private/FoodSeq Reference Database/food-dbs"
scr <- "C:/Users/sab236/AppData/Local/Temp/claude/C--Users-sab236-Box-Home-Folder-sab236-Private-FoodSeq-Reference-Database/5f594644-027b-4313-89aa-f4f890b18de0/scratchpad"

foods <- read_csv(file.path(fd, "data/inputs/human-foods.csv"), show_col_types = FALSE)
cat("rows in human-foods.csv:", nrow(foods), "\n")
cat("columns:", paste(names(foods), collapse = ", "), "\n\n")

# normalise: lowercase, squash whitespace
norm <- function(x) x %>% replace_na("") %>% str_to_lower() %>% str_squish()

sci_set <- norm(foods$scientific_name)
# haystack of every name-bearing field, for synonym/common-name fallback
hay <- paste(norm(foods$scientific_name), norm(foods$common_name),
             norm(foods$alternative_names), sep = " | ")

recon <- read_delim(file.path(scr, "recon_list.txt"), delim = "|",
                    col_names = c("scientific_name", "common_name", "category"),
                    show_col_types = FALSE, trim_ws = TRUE) %>%
  filter(!is.na(scientific_name), scientific_name != "")

recon <- recon %>%
  mutate(
    key          = norm(scientific_name),
    exact_sci    = key %in% sci_set,
    # genus-level presence: is the genus represented at all?
    genus        = word(key, 1),
    genus_present= genus %in% word(sci_set, 1),
    anywhere     = map_lgl(key, ~ any(str_detect(hay, fixed(.x)))),
    status       = case_when(
      exact_sci ~ "present (exact scientific name)",
      anywhere  ~ "present (matched elsewhere: synonym/alt name)",
      TRUE      ~ "ABSENT"
    )
  )

cat("=== RECON RESULT ===\n")
recon %>% count(status, sort = TRUE) %>% print(n = Inf)

cat("\n=== ABSENT (by category) ===\n")
recon %>% filter(status == "ABSENT") %>% count(category) %>% print(n = Inf)

cat("\n=== ABSENT species, with genus-present flag ===\n")
recon %>% filter(status == "ABSENT") %>%
  select(scientific_name, common_name, category, genus_present) %>%
  arrange(category, scientific_name) %>% print(n = Inf)

write_csv(recon, file.path(scr, "recon_result.csv"))
cat("\nwrote recon_result.csv\n")
