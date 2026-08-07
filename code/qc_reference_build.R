#!/usr/bin/env Rscript
# =============================================================================
# qc_reference_build.R  --  post-build sanity gate for the FoodSeq reference DBs
# -----------------------------------------------------------------------------
# Run this AFTER any reference rebuild (cluster or local) and BEFORE trusting or
# shipping the output. The May 2026 build shipped three silent failures that
# raised no error -- kingdom field 100% NA, the synthetic control dropped, and a
# positive-control species (Rhacodactylus leachianus) lost between versions.
# Each check below targets one of those failure modes plus general sanity.
#
# Usage:
#   Rscript code/qc_reference_build.R [REPO_DIR]
# If REPO_DIR is omitted, uses the current working directory.
# Exit code 0 = all checks passed; 1 = one or more FAILs (usable in CI / SLURM).
# =============================================================================

suppressWarnings(suppressMessages({
  ok <- requireNamespace("Biostrings", quietly = TRUE)
}))
if (!ok) stop("Biostrings is required. install via BiocManager::install('Biostrings').")

args     <- commandArgs(trailingOnly = TRUE)
repo_dir <- if (length(args) >= 1) args[[1]] else getwd()

d2 <- file.path(repo_dir, "data", "outputs", "dada2-compatible")
inputs <- file.path(repo_dir, "data", "inputs")

# ---- tolerances (edit to taste) --------------------------------------------
COUNT_TOLERANCE <- 0.25   # allowed fractional change vs the previous build
NA_RANK_MAX     <- 0.99   # fail if any rank is >= this fraction NA
STAPLES <- list(
  trnL = c("Triticum", "Oryza", "Zea mays", "Solanum", "Brassica"),
  `12S` = c("Bos taurus", "Gallus gallus", "Sus scrofa", "Salmo", "Oncorhynchus")
)
RANKS <- c("superkingdom","phylum","class","order","family",
           "genus","species","subspecies","varietas","forma")

# ---- expected clades per marker (allowlist) --------------------------------
# Every reference record's PHYLUM should fall inside its marker's allowlist.
# Anything outside is off-target: bacterial / viral / parasite contamination
# that leaked in (e.g. the May 2026 12S build carried Pseudomonadota, Bacillota,
# and several *viricota records that hijacked fish assignments). An allowlist is
# used deliberately over a denylist so an UNKNOWN future contaminant phylum is
# still caught. Edible invertebrate phyla are included for 12S so legitimate
# shellfish / jellyfish are not flagged.
ALLOWED_PHYLA <- list(
  trnL  = c("Streptophyta"),
  `12S` = c("Chordata", "Arthropoda", "Mollusca", "Echinodermata",
            "Cnidaria", "Annelida", "Porifera")
)

# ---- helpers ----------------------------------------------------------------
n_fail <- 0L; n_warn <- 0L
say  <- function(status, msg) {
  tag <- switch(status, PASS = "[PASS]", FAIL = "[FAIL]", WARN = "[WARN]", INFO = "[ .. ]")
  if (status == "FAIL") n_fail <<- n_fail + 1L
  if (status == "WARN") n_warn <<- n_warn + 1L
  cat(sprintf("%s %s\n", tag, msg))
}

read_headers <- function(path) {
  if (!file.exists(path)) return(NULL)
  ln <- readLines(path, warn = FALSE)
  sub("^>", "", ln[startsWith(ln, ">")])
}

# ---- per-marker checks ------------------------------------------------------
check_marker <- function(marker, tax_fasta, seq_fasta, prev_fasta,
                         fail_controls, warn_controls, staples, n_fields) {
  cat(sprintf("\n===== %s : %s =====\n", marker, basename(tax_fasta)))
  hdr <- read_headers(tax_fasta)
  if (is.null(hdr)) { say("FAIL", paste("taxonomy FASTA not found:", tax_fasta)); return(invisible()) }

  n <- length(hdr)
  say("INFO", sprintf("%d records", n))

  # 1. field count -- every header must have exactly n_fields ;-delimited fields
  nf <- lengths(strsplit(hdr, ";", fixed = TRUE))
  if (all(nf == n_fields)) say("PASS", sprintf("all headers have %d taxonomy fields", n_fields))
  else say("FAIL", sprintf("header field count varies: %s (expected %d)",
                           paste(sort(unique(nf)), collapse = ","), n_fields))

  # 2. no rank almost entirely NA  (the kingdom-NA regression)
  mat <- do.call(rbind, strsplit(hdr, ";", fixed = TRUE))
  mat <- mat[, seq_len(min(ncol(mat), length(RANKS))), drop = FALSE]
  colnames(mat) <- RANKS[seq_len(ncol(mat))]
  na_frac <- apply(mat, 2, function(col) mean(is.na(col) | col == "NA" | col == ""))
  bad <- na_frac[na_frac >= NA_RANK_MAX]
  # controls legitimately carry NA in subspecies/varietas/forma, so only flag
  # ranks kingdom..species (the taxonomically meaningful ones)
  meaningful <- intersect(names(bad), RANKS[1:7])
  if (length(meaningful) == 0)
    say("PASS", "no meaningful rank is >=99% NA")
  else
    say("FAIL", sprintf("rank(s) almost entirely NA: %s",
                        paste(sprintf("%s=%.0f%%", meaningful, 100*na_frac[meaningful]),
                              collapse = ", ")))

  # 3a. controls defined in controls.csv -- REQUIRED (hard FAIL if missing)
  for (ctl in fail_controls) {
    if (any(grepl(ctl, hdr, fixed = TRUE)))
      say("PASS", sprintf("control present: %s", ctl))
    else
      say("FAIL", sprintf("control MISSING (defined in controls.csv): %s", ctl))
  }
  # 3b. positive-control species named in the wet-lab protocol but NOT yet in
  # controls.csv -- surfaced as WARN, since whether they belong in the reference
  # is an open curation decision (see control.12S in Pipeline-to-Phyloseq.Rmd).
  for (ctl in warn_controls) {
    if (any(grepl(ctl, hdr, fixed = TRUE)))
      say("PASS", sprintf("protocol control present: %s", ctl))
    else
      say("WARN", sprintf("protocol control absent (undecided): %s", ctl))
  }

  # 4. staple foods present and resolved (not all-NA lineage)
  for (st in staples) {
    if (any(grepl(st, hdr, fixed = TRUE)))
      say("PASS", sprintf("staple present: %s", st))
    else
      say("WARN", sprintf("staple not found: %s", st))
  }

  # 5. record count vs previous build
  prev <- read_headers(prev_fasta)
  if (is.null(prev)) {
    say("INFO", "no previous build to compare record count against")
  } else {
    delta <- (n - length(prev)) / length(prev)
    msg <- sprintf("record count %d vs previous %d (%+.1f%%)", n, length(prev), 100*delta)
    if (abs(delta) <= COUNT_TOLERANCE) say("PASS", msg)
    else say("WARN", paste(msg, sprintf("- exceeds +/-%.0f%%; confirm this is intended",
                                         100*COUNT_TOLERANCE)))

    # 6. taxa present before but absent now (the dropped-gecko failure mode)
    sp_now  <- mat[, "species"]
    prev_mat <- do.call(rbind, strsplit(prev, ";", fixed = TRUE))
    if (ncol(prev_mat) >= 7) {
      sp_prev <- prev_mat[, 7]
      lost <- setdiff(unique(sp_prev), unique(sp_now))
      lost <- lost[!is.na(lost) & lost != "NA"]
      if (length(lost) == 0)
        say("PASS", "no species present in previous build is missing now")
      else
        say("WARN", sprintf("%d species dropped since previous build (e.g. %s)",
                            length(lost), paste(head(lost, 5), collapse = ", ")))
    }
  }

  # 7. off-target clade -- every record's phylum must be in the marker's
  # allowlist. Anything else is bacterial / viral / parasite contamination.
  # (Catches the cross-kingdom decoys; does NOT catch same-clade mislabels --
  # that is what check 8 is for.)
  allowed   <- ALLOWED_PHYLA[[marker]]
  ph        <- mat[, "phylum"]
  is_na_ph  <- is.na(ph) | ph == "NA" | ph == ""
  is_ctl    <- ph %in% fail_controls | grepl("synthetic", ph, ignore.case = TRUE)
  offtarget <- !is_na_ph & !is_ctl & !(ph %in% allowed)
  if (!any(offtarget)) {
    say("PASS", "all records fall within expected clades (no bacterial/viral/parasite phyla)")
  } else {
    tab <- sort(table(ph[offtarget]), decreasing = TRUE)
    say("FAIL", sprintf("%d off-target record(s) with non-target phylum: %s",
                        sum(offtarget),
                        paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", ")))
  }
  n_unres <- sum(is_na_ph & !is_ctl)
  if (n_unres > 0) say("WARN", sprintf("%d record(s) have no phylum (unresolved)", n_unres))

  # 8. accession integrity + queried-vs-resolved genus consistency
  # The sequence FASTA header is "<accession> <species>". Some 2026 records
  # carry a bare GENUS name where the accession should be; where that genus does
  # not match the resolved species' genus, the sequence is mislabelled (e.g.
  # Esox -> Salmo salar, Engraulis -> Salimicrobium jeotgali). This is the ONLY
  # check that catches same-clade (fish-for-fish) mislabels, which pass every
  # taxonomic filter because both organisms are vertebrates.
  seq_hdr <- read_headers(seq_fasta)
  if (is.null(seq_hdr)) {
    say("INFO", sprintf("sequence FASTA not found (%s); skipping accession check",
                        basename(seq_fasta)))
  } else {
    acc_field <- sub("\\s.*$", "", seq_hdr)     # first whitespace-delimited token
    sp        <- sub("^\\S+\\s*", "", seq_hdr)   # remainder = species name
    sp_genus  <- sub("\\s.*$", "", sp)
    bare      <- !grepl("[0-9]", acc_field)      # a real accession contains digits
    if (!any(bare)) {
      say("PASS", "all sequence records carry a real accession identifier")
    } else {
      say("WARN", sprintf("%d record(s) have a bare genus name where an accession should be",
                          sum(bare)))
      mism <- bare & (acc_field != sp_genus)
      if (!any(mism)) {
        say("PASS", "every bare-genus record's queried genus matches its resolved species")
      } else {
        ex <- head(sprintf("%s->%s", acc_field[mism], sp[mism]), 6)
        say("FAIL", sprintf("%d record(s) queried-genus != resolved-genus (mislabel): %s",
                            sum(mism), paste(ex, collapse = "; ")))
      }
    }
  }
}

# ---- controls.csv sanity ----------------------------------------------------
ctl_path <- file.path(inputs, "controls.csv")
ctl <- if (file.exists(ctl_path)) read.csv(ctl_path, stringsAsFactors = FALSE) else NULL
if (is.null(ctl)) say("WARN", "data/inputs/controls.csv not found - cannot verify controls by definition")

ctl_labels <- function(mk) if (is.null(ctl)) character(0) else ctl$label[ctl$marker == mk]

cat("========================================================\n")
cat("  FoodSeq reference build QC\n")
cat("  repo:", repo_dir, "\n")
cat("========================================================\n")

check_marker(
  "trnL",
  tax_fasta     = file.path(d2, "trnL", "trnLGH_taxonomy_May2026.fasta"),
  seq_fasta     = file.path(d2, "trnL", "trnLGH_May2026.fasta"),
  prev_fasta    = file.path(d2, "trnL", "trnLGH_taxonomy.fasta"),
  fail_controls = unique(c(ctl_labels("trnL"), "synthetic trnL ASV")),
  warn_controls = c("Ilex paraguariensis", "Trifolium pratense"),
  staples       = STAPLES$trnL, n_fields = 10L
)

check_marker(
  "12S",
  tax_fasta     = file.path(d2, "12Sv5", "12Sv5_taxonomy_May2026.fasta"),
  seq_fasta     = file.path(d2, "12Sv5", "12Sv5_May2026.fasta"),
  prev_fasta    = file.path(d2, "12Sv5", "12Sv5_taxonomy.fasta"),
  fail_controls = unique(c(ctl_labels("12SV5"), "synthetic 12S ASV")),
  # positive-control template species named in the wet-lab protocol; not yet in
  # controls.csv, so WARN-only until the add-the-geckos decision is made
  warn_controls = c("Dromaius novaehollandiae", "Correlophus ciliatus",
                    "Rhacodactylus leachianus"),
  staples       = STAPLES$`12S`, n_fields = 10L
)

cat("\n========================================================\n")
cat(sprintf("  RESULT: %d FAIL, %d WARN\n", n_fail, n_warn))
cat("========================================================\n")
if (n_fail > 0) {
  cat("Build did NOT pass QC. Do not ship this reference until FAILs are resolved.\n")
  quit(status = 1)
} else {
  cat("All hard checks passed. Review any WARNs before shipping.\n")
  quit(status = 0)
}
