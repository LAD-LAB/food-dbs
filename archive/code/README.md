# Archived code

Superseded notebooks and scripts, kept for history. Nothing here is run by the
current pipeline (`foodseq_reference_pipeline*.Rmd`, `code/`). Several files
have hard-coded paths from the machines they were written on and will not run
as-is.

| File / folder | What it is |
|---|---|
| `trnL-reference.Rmd` | Early trnL (g/h) reference build. |
| `12SV5-reference.Rmd` | Early 12SV5 reference build. |
| `trnLcd-reference.Rmd` | Early trnL c/d build by Ashish Subramanian — see below. |
| `coverage-recheck/` | Scripts from the July 2026 coverage re-check, since consolidated into `code/Coverage recheck.Rmd`. |
| `functions/` | Helper functions used by the notebooks above. |

## `trnLcd-reference.Rmd`

Added 2026-10-28. Author Ashish Subramanian, dated 2024-11-26, run from his
own machine (hard-coded `/Users/ams292/...` paths for the food list, the manual
renaming table, the RefSeq plastid FASTA and the `accessionTaxa.sql` database).

**What it does.** Builds a trnL reference trimmed to an amplicon, by primer
matching. It takes RefSeq plastid sequences and a GenBank pull from
`query_ncbi()`, keeps those that carry both primer sites (`find_primer_pair()`),
drops UNVERIFIED records, merges them with manual additions, and looks up
taxonomy with `taxonomizr`. It follows the structure and wording of the g/h
notebook (`trnL-reference.Rmd`) with different primers.

**What it is not.** It is not the code behind the 2025 `trnLCD.fasta`
(`data/outputs/dada2-compatible/archive/trnLCD/`, committed 2025-03-03 with no
build code). Two reasons:

1. **The primers are not trnL-c / trnL-d.** The notebook's `trnLC` is
   `GATAGGTGCAGAGACTCAATGG`, which is the reverse complement of trnL-**h**
   (`CCATTGAGTCTCTGCACCTATC`). The real trnL-c is `CGAAATCGGTAGACGCTACG`, as used
   in the current pipeline's Part A-CD. `trnLD` (`GGGGATAGAGGGACTTGAAC`) is
   correct.
2. **The method differs.** An analysis in August 2026 showed the 2025
   `trnLCD.fasta` was not produced by primer matching: 644 of its 2,239
   records come from source records with no c/d primer site, and the amplicon
   boundaries sit at homologous but divergent positions inside the source
   records. That points to extraction by annotation coordinates or alignment.
   This notebook uses `find_primer_pair()`, so it cannot reproduce that file.

So it is best read as an earlier or different attempt. It does not close the
gap that the build code for the 2025 trnL c/d database is missing. The current
pipeline builds c/d in Part A-CD with alignment-augmented extraction instead.
