# food-dbs

Reference database pipeline for dietary metabarcoding using the trnL (plant) and 12SV5 (vertebrate) markers. Builds DADA2- and QIIME2-compatible reference databases by combining RefSeq sequences downloaded locally with remote GenBank queries via the NCBI API.

> **Just adding a few species to an existing reference?** Skip "Getting started" below — that's for a full cluster rebuild. See "Using Extend reference.Rmd and Coverage recheck.Rmd" instead: it runs on a laptop in minutes and needs no cluster, no SLURM, and none of the large downloads (see "Extend, or rebuild?" for the full distinction).

------------------------------------------------------------------------

## Repository structure

```         
food-dbs/
├── foodseq_reference_pipeline.Rmd  # Main pipeline (start here)
├── food-dbs.Rproj                  # RStudio project file
│
├── code/
│   ├── functions/                  # Functions sourced by the pipeline
│   │   ├── find_primer_pair.R      # In silico PCR trimming
│   │   ├── query_ncbi.R            # Batch NCBI nucleotide queries
│   │   ├── query_ncbi_accession.R  # Resolve accessions to taxon IDs
│   │   └── fetch_lineage_ncbi.R    # SQL-free taxid -> lineage lookup (used by Extend reference.Rmd)
│   ├── Descriptive-statistics.Rmd  # Summary statistics for built databases
│   ├── Extend reference.Rmd        # Add species to an existing reference without a rebuild - see "Extend, or rebuild?" below
│   ├── Parse ecoPCR.Rmd            # Parse ecoPCR output
│   ├── SQL to reference.Rmd        # Alternative taxonomy approach
│   ├── Taxa names.Rmd              # Taxon name handling
│   └── tree-building.Rmd           # Phylogenetic tree construction
│
├── data/
│   ├── inputs/
│   │   ├── human-foods.csv         # Curated list of food plant and animal species (`build_name` column maps NCBI synonyms)
│   │   ├── Manual renaming.csv     # Manual additions (merged in at step 1) and omissions (step 6) for trnL; the Brassica variety renamings are hard-coded in the pipeline, not read from here
│   │   ├── reference-additions.csv # Species list consumed by Extend reference.Rmd (blank template)
│   │   ├── reference-additions-*.csv # The lists used for past Extend runs, kept as a record of what was added and why
│   │   └── known-absent-species.csv # Food-list species deliberately absent from a build, with the reason (read by the QC coverage ratchet)
│   └── outputs/
│       ├── dada2-compatible/       # Reference FASTAs for use with DADA2
│       │   ├── trnL/               # trnL database, current release (_v2026.10)
│       │   ├── 12Sv5/              # 12SV5 database, current release (_v2026.10)
│       │   ├── trnLCD/             # trnL c/d database, current release (_v2026.10)
│       │   ├── miscellaneous/      # Other marker databases
│       │   └── archive/            # Superseded releases (aug2026 = Zenodo v2026.08, may2026, 2022) and recovery snapshots
│       ├── qiime2-compatible/      # Reference FASTAs and TSVs for use with QIIME2
│       │   ├── trnL/               # _v2026.10
│       │   ├── 12Sv5/
│       │   ├── trnLCD/
│       │   └── archive/            # superseded releases, as above
│       ├── common-names/           # trnL and 12S ASV common names CSVs, built against a pre-Aug 2026 reference
│       ├── plants_missing_trnL.csv    # Food plants without trnL coverage (current run)
│       └── animals_missing_12SV5.csv  # Food animals without 12SV5 coverage (current run)
│
└── archive/                        # Superseded pipeline files (kept for reference)
    ├── code/
    │   ├── trnL-reference.Rmd
    │   ├── 12SV5-reference.Rmd
    │   └── functions/
└── ...
```

------------------------------------------------------------------------

## Acquiring FoodSeq reference databases

This repository contains the latest FoodSeq reference databases for use with the Edible Atlas FoodSeq Handbook. FASTA files containing reference sequences can be downloaded from this repository by going to data/outputs/ and selecting either dada2-compatible or qiime2-compatible FASTA files. File names carry the release version (`_v2026.10`), matching the Zenodo version numbering; earlier releases are kept in each folder's `archive/`. Alternatively, the latest version of the reference database is available at <https://doi.org/10.5281/zenodo.22756737>

## How the pipeline works

The pipeline builds two reference databases in parallel, following the same steps for each:

### Data sources

Sequences are drawn from two sources and merged, which apply the food species list at different points:

-   **RefSeq** (local): Full plastid (trnL) or mitochondrial (12SV5) genomes downloaded from the NCBI FTP server, filtered to one sequence per species, then in silico amplified using the target primer pair for *all* RefSeq species — the result is subset down to the food species list only after amplification.
-   **GenBank** (remote): Sequences queried directly from NCBI's nucleotide database using `query_ncbi()`, which searches only for the target marker across species already in the food species list — the query itself is scoped to `human-foods.csv`, so there's no separate post-hoc filter step on this side.

See "Pipeline flowchart" below for how this looks as a diagram.

### Species lists

-   **trnL**: All species in `human-foods.csv` where `category == "plant"`
-   **12SV5**: All species in `human-foods.csv` where `category == "animal"`

### Taxonomy

Accession numbers are mapped to full taxonomic lineages using a local `taxonomizr` SQL database (a mirror of NCBI taxonomy). Any accessions added to NCBI after the SQL database was built are resolved automatically using `query_ncbi_accession()`.

### Processing steps

For each marker, the pipeline:

1.  Downloads RefSeq genomes and queries GenBank (via SLURM jobs) — the GenBank query is already scoped to the food species list; RefSeq is not yet
2.  Applies in silico PCR with `find_primer_pair()` to extract the target amplicon
3.  Subsets the RefSeq amplicons down to the food species list (GenBank needs no equivalent step — see "Data sources" above)
4.  Removes unverified sequences and sequences with degenerate nucleotides
5.  Looks up full taxonomy for all accessions
6.  Applies manual curation (trnL only): omits the records listed as `omit` in `Manual renaming.csv`, and applies a few hard-coded *Brassica oleracea* variety renamings. The file's `add` rows, mostly whole-genome records that the 50 kb query cap cannot reach, are merged in earlier, alongside the GenBank sequences
7.  Standardises sequence orientation
8.  Deduplicates (identical sequences from the same taxon collapsed to one, with RefSeq accessions preferred)
9.  Saves DADA2- and QIIME2-compatible output files

### Primers

| Marker | Forward | Reverse |
|----|----|----|
| trnL (plant) | `GGGCAATCCTGAGCCAA` (trnLg) | `CCATTGAGTCTCTGCACCTATC` (trnLh) |
| 12SV5 (vertebrate) | `TAGAACAGGCTCCTCTAG` (V5F) | `TTAGATACCCCACTATGC` (V5R) |

### Pipeline flowchart

The same diagram as `pipeline-flowchart.html` (open that file directly for the interactive version, plus a table of everything that differs between the trnL and 12SV5 runs — species list, primers, phylum allowlist, output filenames):

``` mermaid
flowchart TD

    CTL["controls.csv"]
    SETUP["SLURM setup jobs"]
    REFSEQ_DL["RefSeq plastid + RefSeq mitochondrial"]
    SQL[("accessionTaxa.sql")]
    NCBI_FB["query_ncbi_accession<br/>accessions newer than the SQL build"]

    SETUP -->|download| REFSEQ_DL
    SETUP -->|prepareDatabase| SQL

    REFSEQ_DL -->|filter one per species · trnL, 12SV5| GENOME_SP["② RefSeq genomes<br/>one per species, accession kept"]

    GENOME_SP -->|find_primer_pair, all RefSeq species| REFSEQ["③ RefSeq amplicons"]
    NCBI_SPECIES["① plants or animals from human-foods.csv"] -->|query_ncbi + find_primer_pair| NCBI["④ GenBank sequences"]

    REFSEQ -->|"① Use plants or animals from human-foods.csv to filter amplicons"| REFSEQ_FOOD["RefSeq food amplicons"]

    REFSEQ_FOOD --> COMB["⑤ combined sequences"]
    NCBI -->|combine · trnL also merges manual additions| COMB
    COMB -->|accessionToTaxa| TAX["⑥ + taxonomy"]
    SQL --> TAX
    NCBI_FB -.->|fallback for missing IDs| TAX
    TAX -->|A6b/B6b · phylum allowlist| CLEAN["⑦ on-target only"]

    CLEAN -->|⑧ QC · ⑨ orientation · ⑩ dedup · add controls| OUT["⑩ reference database"]
    CLEAN -.->|trnL only| CUR["⑦b manual curation<br/>omissions from Manual renaming.csv<br/>+ hard-coded variety renamings"]
    CUR -.->|⑧ QC · ⑨ orientation · ⑩ dedup · add controls| OUT
    CTL -.->|add controls| OUT

    QC{"qc_reference_build.R<br/>controls · rank completeness · clade consistency<br/>accession integrity · orientation · count tolerance<br/>coverage ratchet vs the 2025 reference<br/>(known-absent-species.csv · build_name)"}
    OUT --> QC
    QC -->|pass| SHIP["release<br/>files named _vYYYY.MM · previous release moved to archive/"]
    QC -->|fail, exit 1| BLOCK["do not ship"]

    NCBI_SPECIES -.->|trnL only · same species list and NCBI pull| CD["Part A-CD · trnL c/d<br/>alignment-augmented extraction:<br/>plastome introns · GenBank primers · GenBank alignment rescue"]
    GENOME_SP -.->|trnL only · plastome introns are the alignment panel| CD
    CD -.->|⑧ QC · ⑩ dedup · add controls| OUTCD["trnLCD reference database"]
    OUTCD -.-> QC

    SHIP -.->|after a release: species still missing| RECHECK["Coverage recheck.Rmd<br/>which missing species now have data at NCBI?"]
    RECHECK -.-> EXT["Extend reference.Rmd · laptop, no cluster<br/>SOURCE marker: NCBI records under 50 kb<br/>SOURCE plastome (trnL): complete chloroplast genomes<br/>forward-orients amplicons · rejects out-of-range lengths"]
    EXT -.->|new suffixed pair · gate again| QC

    style CUR fill:#F5EEDC,stroke:#C9A96E,color:#5F4A1E
    style QC fill:#C8E0C9,stroke:#2C5F2D,color:#1E3A1F
    style BLOCK fill:#F2DCD6,stroke:#BC5138,color:#5F1E1E
    style SHIP fill:#CDE3CE,stroke:#2C5F2D,color:#1E5F1E
    style CD fill:#E6EDF7,stroke:#7A93BC,color:#1E2F4F
    style OUTCD fill:#E6EDF7,stroke:#7A93BC,color:#1E2F4F
    style RECHECK fill:#EEE8F5,stroke:#9A82B8,color:#3A2A55
    style EXT fill:#EEE8F5,stroke:#9A82B8,color:#3A2A55
```

------------------------------------------------------------------------

## How to run this pipeline

There are two ways to run the FoodSeq reference sequence database creation yourself: you can either re-create the full database from scratch, or you can add reference sequences to an existing version of the database.

### Extend, or rebuild?

#### Extend (`code/Extend reference.Rmd`)

Adds a handful of newly-available species to an **existing** reference without touching any record already in it. It runs on a laptop in minutes: it queries NCBI directly for both sequences (`query_ncbi()`) and taxonomy (`fetch_lineage_ncbi()`, an e-utils-based lookup), so it doesn't need the cluster or the \~70 GB `accessionTaxa.sql` database the full pipeline depends on.

**Use it when:**

-   New species have shown up at NCBI since the last build for names already in `human-foods.csv`, or you're closing gaps identified by a coverage re-check (`data/outputs/coverage-recheck/CANDIDATES_*.csv` are already in the expected input format — a CSV with a `scientific_name` column)
-   You're adding species that were never on the target list at all (point `ADDITIONS` at any CSV with a `scientific_name` column, e.g. `data/inputs/reference-additions.csv`)
-   The only thing needed is new sequence — nothing about the existing records needs to change

**It can't help with anything that has to change an existing record**, because it only appends:

-   A taxonomy or curation fix, a mislabeled or renamed record
-   A change to the off-target filter, the dedup logic, or any other pipeline step
-   trnLCD — its extraction is alignment-based against a panel of plastome introns, not a single per-species NCBI query, so extending it needs a rebuild
-   A species whose only NCBI record is a complete genome — by default Extend applies the same 50 kb length cap as the full pipeline (`query_ncbi()`), so a complete plastid genome comes back empty even though sequence technically exists. For plants this is solved by setting `SOURCE <- "plastome"` (trnL only): it searches complete chloroplast genomes by organism and title, fetches them one at a time, and cuts the amplicon out by primer matching, with the same QC as the default path. Complete vertebrate mitogenomes (~16 kb) are already under the cap, so 12SV5 needs no equivalent

Extend never overwrites the reference it's extending — it always writes a new, separately-suffixed pair of files (e.g. `_v2026.10` → `_v2026.10_ext`), so a bad run is easy to discard and re-run. Like the full pipeline, it stores every amplicon forward-primer-first and rejects amplicons far outside the reference's length range; the QC gate fails on any record that starts with the reverse primer. Before shipping the result, run the QC gate against it with the marker as a 4th argument, e.g.:

``` bash
Rscript code/qc_reference_build.R . _v2026.10_ext _v2026.10 trnL
```

The 4th argument restricts the gate to the marker you actually extended — without it, the gate checks every marker under the new suffix and hard-fails on the ones you didn't touch, which reads as a failure even on a clean run.

#### Rebuild (`foodseq_reference_pipeline.Rmd`)

Needed for anything Extend can't do (above), and for periodic full refreshes that pick up whatever new sequence data has accumulated at NCBI more broadly since the last build. Requires the cluster setup described under **Getting started** below.

------------------------------------------------------------------------

### Getting started

#### Prerequisites

**R packages:**

``` r
# Bioconductor
BiocManager::install(c("Biostrings", "ShortRead"))

# CRAN
install.packages(c("taxonomizr", "tidyverse", "rentrez", "remotes"))

# GitHub
remotes::install_github("ammararuby/MButils")
```

**NCBI API key** (recommended): Create an NCBI account, go to Account Settings → API Key Management, and copy your key. In R, run:

``` r
rentrez::set_entrez_key("your_key_here")
```

This increases the NCBI query rate limit from 3 to 10 requests per second. NCBI also recommends running large queries on weekends or between 9 PM and 5 AM EST on weekdays.

#### Large file setup (cluster recommended)

Two large files must be obtained before running the pipeline. These are best downloaded on an HPC cluster due to their size and download time:

| File | Size | How to obtain |
|----|----|----|
| RefSeq plastid FASTA | \~15–20 GB uncompressed | NCBI FTP: `ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/` |
| RefSeq mitochondrial FASTA | \~5–10 GB uncompressed | NCBI FTP: `ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/` |
| `accessionTaxa.sql` | \~70 GB | Built with `taxonomizr::prepareDatabase()` |

SLURM job scripts for all three are generated and submitted automatically by the pipeline (sections 2a–2c of the Rmd).

#### Running the pipeline

1.  Clone this repository:

    ``` bash
    git clone https://github.com/LAD-LAB/food-dbs.git
    ```

2.  Open `foodseq_reference_pipeline.Rmd` in RStudio (or via RStudio Server on Open OnDemand)

3.  Update the paths in the **Configuration** chunk (section 0) to match your environment:

    ``` r
    SCRATCH  <- "/scratch/your_username"
    REPO_DIR <- "/path/to/food-dbs"
    SQL_PATH <- "/path/to/accessionTaxa.sql"
    QC_PREVIOUS_SUFFIX <- "_v2026.10"  # last archived build, for the QC gate
    ```

4.  Run the **Install packages** chunk (section 1) once on first use

5.  Submit the three SLURM jobs (sections 2a–2c) and monitor their progress (section 3)

6.  Once jobs are complete, run **Part A** (trnL) and/or **Part B** (12SV5) sequentially

7.  The notebook's final section, **Verify before shipping**, runs the QC gate automatically and prints the result — 0 FAIL means the build is safe to archive. If `QC_PREVIOUS_SUFFIX` was left at its default (`"-"`), this step fails loudly rather than comparing against the wrong build; set it and re-run just that chunk.

Output files are written to `data/outputs/dada2-compatible/` and `data/outputs/qiime2-compatible/` as date-less files (`trnLGH.fasta`, not `trnLGH_v2026.10.fasta`) — once the QC gate passes, archive them under a new suffix before the next rebuild overwrites them in place.

------------------------------------------------------------------------

### Using Extend reference.Rmd and Coverage recheck.Rmd

Together these grow an existing reference from a laptop, without the cluster setup above — see "Extend, or rebuild?" for when this is (and isn't) the right tool. The usual order is: find candidates, pull sequence, verify, promote.

#### 1. Find candidates — `code/Coverage recheck.Rmd`

Open in RStudio and run top to bottom. Configuration options:

-   `RUN_PHASE1` (default `TRUE`): re-checks species already on `human-foods.csv` but missing from the current build against NCBI, resolving each name to its current NCBI-accepted synonym first so a renamed species isn't misreported as missing. Writes `data/outputs/coverage-recheck/CANDIDATES_plants_trnL<suffix>.csv` and `CANDIDATES_animals_12SV5<suffix>.csv`. This is the routine thing to run before an Extend session.
-   `RUN_PHASE2` (default `TRUE`, set `FALSE` to skip): audits a curated list of \~199 globally significant foods against `human-foods.csv` for species missing from the target list entirely, not just missing sequence. An occasional audit rather than something to re-run every time.
-   `REF_SUFFIX` / `OUT_SUFFIX`: which built reference to check "missing" against, and what to suffix the output CSVs with.

Queries are checkpointed to `data/outputs/coverage-recheck/raw/`, so a killed or resumed run doesn't re-query species already recorded — Part 1 alone is a multi-hour job against \~1,500 species without an NCBI API key. A hit means NCBI has *some* record for that species, not a confirmed amplicon — `Extend reference.Rmd` is the real test.

#### 2. Pull sequence — `code/Extend reference.Rmd`

Point `ADDITIONS` at a `CANDIDATES_*.csv` from step 1 (or any CSV with a `scientific_name` column — `data/inputs/reference-additions.csv` is a blank template for ad hoc additions):

``` r
MARKER     <- "trnL"        # or "12SV5" - run once per marker
ADDITIONS  <- here('data', 'outputs', 'coverage-recheck',
                    'CANDIDATES_plants_trnL_Aug2026.csv')
SOURCE     <- "marker"      # or "plastome" (trnL only): complete chloroplast genomes
IN_SUFFIX  <- "_v2026.10"   # the reference you're extending
OUT_SUFFIX <- "_v2026.10_ext" # what to write - never overwrites IN_SUFFIX
```

Run top to bottom. Species already in the reference are skipped automatically; species with no primer-spanning record are reported at the end rather than silently dropped.

#### 3. Verify before shipping

``` bash
Rscript code/qc_reference_build.R . _v2026.10_ext _v2026.10 trnL
```

The same gate the full pipeline uses. The 4th argument restricts it to the marker you extended — without it, the gate checks every marker under the new suffix and hard-fails on the ones you didn't touch. A clean run should show the record count rise by exactly what was added and every other check still pass.

#### 4. Promote the result

If the QC gate passes, rename the `OUT_SUFFIX` files to whatever the new canonical suffix is (e.g. `_v2026.10_ext` → `_v2026.11`) so they become what the rest of the repo — and the next `Extend` or `Coverage recheck` run — treats as current.

------------------------------------------------------------------------

## Database coverage

### Food species list

The pipeline targets all species in `human-foods.csv`, a curated list of 3,806 food species assembled from 32+ literature and database sources (29 species were added in a July 2026 coverage re-check; see "Coverage over time" below for how that affects comparisons across builds).

| Category                             | Species   |
|--------------------------------------|-----------|
| Plants (species-level)               | 1,573     |
| Plants (genus/family-level entries)  | 7         |
| Animals (species-level)              | 2,121     |
| Animals (genus/family-level entries) | 70        |
| Fungi                                | 33        |
| Bacteria                             | 2         |
| **Total**                            | **3,806** |

**Synonyms (`build_name`).** NCBI sometimes files a species' sequences under a different name from the one on the food list: a rename (*Pisum sativum* is now *Lathyrus oleraceus*), a spelling variant (*Ctenopharyngodon idellus* / *idella*), or an infraspecific name whose species is in the reference (*Brassica oleracea* var. *italica*). The `build_name` column of `human-foods.csv` records the name the reference actually uses for such entries — 185 of them (131 plants, 54 animals). `scientific_name` is never changed, so the food list and its denominators are unaffected. A species counts as covered if either its `scientific_name` or its `build_name` is a species in the reference, and the QC gate applies the same rule. The mappings come from two audits in October 2026: the taxon NCBI assigns to the records it returns for each name, and a lookup of every uncovered name in NCBI Taxonomy (per-name results in `data/outputs/coverage-recheck/SYNONYM_pass_Oct2026.csv`). Hybrids are not mapped to a single parent. To add a mapping, put the reference's name in `build_name` (separate several with a comma); the QC gate warns if a `build_name` is no longer in the build.

### Coverage over time

**trnL (plants)**

| Version  | Sequences | Unique taxa | Food plants covered (exact name) | Coverage | Incl. synonyms  | Coverage |
|----------|-----------|-------------|----------------------------------|----------|-----------------|----------|
| Oct 2022 | 1,402     | 807         | 716 / 1,570                      | 46%      | —               | —        |
| 2025     | 1,991     | 1,169       | 1,060 / 1,570                    | 68%      | —               | —        |
| May 2026 | 1,991     | 1,169       | 1,060 / 1,570                    | 68%      | —               | —        |
| Aug 2026 (Zenodo v2026.08) | 2,027 | 1,196 | 1,078 / 1,573             | 69%      | 1,225 / 1,573   | 78%      |
| v2026.10 | 2,116     | 1,272       | 1,137 / 1,573                    | 72%      | 1,302 / 1,573   | 83%      |

> **Exact name vs. including synonyms.** The "exact name" columns count a food-list species only when the reference holds a species of that exact name — the method used for every earlier row, so those rows compare directly. The "incl. synonyms" columns also count species the reference holds under their `build_name` (see "Synonyms" above). Earlier builds were not recomputed that way.

> Plant denominators exclude 7 genus-only entries in `human-foods.csv` (*Eucheuma, Fragaria, Gelidium, Gracilaria, Gurania, Mentha, Xanthosoma*) that aren't species-level binomials and can't be matched against a species-level reference — applied consistently across all four rows above, the same way the 12SV5 table below excludes non-vertebrates from every row.

**12SV5 (vertebrates)**

| Version  | Sequences | Unique taxa | Food animals covered (exact name) | Coverage | Incl. synonyms  | Coverage |
|----------|-----------|-------------|-----------------------------------|----------|-----------------|----------|
| Oct 2022 | 57        | —           | —                                 | —        | —               | —        |
| 2025     | 2,991     | 2,112       | 1,099 / 2,095                     | 52%      | —               | —        |
| May 2026 | 3,390     | 2,337       | 1,168 / 2,095                     | 56%      | —               | —        |
| Aug 2026 (Zenodo v2026.08) | 2,384 | 1,720 | 982 / 2,121             | 46%      | 1,042 / 2,121   | 49%      |
| v2026.10 | 2,684     | 1,984       | 1,238 / 2,121                     | 58%      | 1,307 / 2,121   | 62%      |

> The Oct 2022 12SV5 database was a small food-filtered subset of the Schneider et al. database. The 2025 and May 2026 databases were built from scratch using the full pipeline.

> **Reading the last two rows:** `human-foods.csv` grew between the May and Aug 2026 builds (2,095 → 2,121 food-animal species, 1,570 → 1,573 food-plant species — see "Food species list" above), so the Aug 2026 and v2026.10 coverage fractions are not directly comparable to earlier rows via raw denominator. Recomputed against the *current* list, May 2026 covers 1,060 / 1,573 (67%) plants and 891 / 2,121 (42%) animals — so v2026.10 is a genuine improvement on both markers (72% and 58% by exact name; 83% and 62% counting synonyms). The released Aug 2026 12SV5 database covered only 46% of food animals; the recoveries described below account for most of the difference. The Aug 2026 rebuild added a human host-taxon control, fixed 21 mislabelled off-target records, fixed a bare-genus accession bug, and closed 43 animal / 83 plant gaps found by the July 2026 coverage re-check — see `data/outputs/coverage-recheck/README.md` for the full account.

> **September 2026 recovery.** The figures above already include a follow-up fix: 51 animal species with a confirmed valid 12SV5 amplicon — present in the May 2026 build, or found by the July coverage re-check — had been silently dropped from the original Aug 2026 rebuild. Root cause: `query_ncbi()` queries species in batches of 5 joined by `OR`, then fetches only the first 500 combined results (`retmax_fetch`), so a rare species batched with a heavily-sequenced one (chicken, pig) can be crowded out of the fetch entirely with no error. `code/Extend reference.Rmd` was used to re-query and merge these back in (+53 sequences, +42 unique taxa). One species (*Acanthurus gahhm*) remains unrecovered despite a confirmed valid amplicon — same unexplained character as the "39 unexplained" cases below.

> **September 2026 trnL control.** The trnL row above gained one sequence (2,026 → 2,027; 1,195 → 1,196 taxa) for the same reason the human and gecko 12SV5 controls were added: *Nicotiana tabacum* (tobacco) reaches dietary samples as exposure or contamination rather than as food, so it is correctly absent from `human-foods.csv` — and was therefore dropped silently when the pipeline became food-list-driven. Comparing the megaphyloseq across both reference vintages found 600,342 reads across 334 of 21,488 samples that went from a confident tobacco call to no assignment at all. It is now a real-species control in `data/inputs/controls.csv`, so every future rebuild includes it and the QC gate hard-fails if it goes missing again. Coverage percentages are unaffected — tobacco is not a food plant and does not enter the denominator.

> **October 2026 recovery and QC change.** The September check compared a build only with its immediate predecessor, so it could not see species an earlier step had already dropped. Comparing the Aug 2026 build with the 2025 reference found 202 food-list animal species that the 2025 12SV5 reference held and Aug 2026 lacked (cod, pollock and tuna among them). Every one has a valid V5 amplicon at NCBI, and all were re-added with `code/Extend reference.Rmd` (10 species, then 191, then *Sander volgensis*), taking 12SV5 from 2,384 to 2,632 sequences and from 1,720 to 1,922 unique taxa. The pre-recovery states are kept in `archive/aug2026-pre-*` under `dada2-compatible/` and `qiime2-compatible/`. The cause of the loss is established for only two of the first ten (the 500-record batch cap described above); for the rest it is not known. The same comparison for trnL found no undocumented losses. `qc_reference_build.R` now has a coverage ratchet (check 6b): every food-list species in the 2025 reference must be present in the build, directly or through `build_name`, unless `data/inputs/known-absent-species.csv` records why not — an undocumented loss is a FAIL. Seven plants (*Amaranthus cruentus, A. dubius, A. thunbergii, Cordia africana, Costus erythrophyllus, Vaccinium tenellum, Victoria amazonica*) have a g/h amplicon at NCBI that the build excludes because it contains ambiguity codes (N, K, S or R); the build drops every sequence with a non-ACGT base. A further pass on 2026-10-28 added the 58 animals that the gap check (below) found to have a valid amplicon at NCBI but no place in the build, including the yak (*Bos grunniens*), llama and guanaco, 53 of them under the name on the food list and 5 under a different NCBI name (mapped with `build_name`). That added 64 sequences and 58 taxa, taking 12SV5 to 2,696 sequences and 1,980 unique taxa; the pre-run state is in `archive/aug2026-pre-gapcheck58/`.

> **v2026.10 release (28 October 2026).** Files are now named `_v2026.10` to match the Zenodo version; the Aug 2026 release is kept in each `archive/aug2026/` folder (identical in content to the deposited files, but stored with LF line endings: the deposited 12SV5 and trnL g/h files are CRLF), and the QC gate finds it there as the previous build. Since that release, 12SV5 has grown from 2,384 to 2,684 sequences and trnL from 2,027 to 2,116. Besides the animal recoveries described above, two things were added. **(1) A name review.** The 149 food-list names that NCBI Taxonomy did not recognise were checked against GBIF and NCBI (`data/outputs/coverage-recheck/NAMES_resolution_Oct2026.csv`, one row per name). 15 misspellings were corrected in `human-foods.csv` (the original spelling is kept in the `notes` column), and 31 names were mapped to the name the reference uses through `build_name`, 5 of them rows that duplicate a correctly spelled row already on the list. 4 vertebrates and 5 plants gained sequence from the corrected names directly, and more through the plastome path below. Synonymies that GBIF asserted but that looked doubtful (e.g. *Melilotus macrorhizus*, *Omphalea megacarpa*, *Polyprion yanezi*) were left alone. 57 of the 101 unrecognised animals are molluscs, crustaceans or insects, which a Chordata-only 12SV5 reference cannot hold. **(2) A plastome path for trnL.** `Extend reference.Rmd` gained `SOURCE <- "plastome"`. A plant whose only trnL sequence sits inside a complete chloroplast genome was invisible to every earlier path: the 50 kb length cap removes the genome, and plastome titles never name the gene, so the marker search misses 96% of them. The new path searches by organism and title, and all 344 plants then still missing were checked. 73 had a plastome with a clean trnL amplicon (70–117 bp, the same range as the existing records) and were added as 83 sequences covering 76 taxa; 15 of them are filed by NCBI under another name (e.g. *Polianthes tuberosa* is *Agave amica*) and are mapped with `build_name`; two of these were among the 149 unrecognised names. Four `known-absent-species.csv` rows (*Alpinia officinarum*, *Amaranthus cruentus*, *A. dubius*, *Triticum carthlicum*) were superseded, as the QC gate flagged, and removed. **(3) An orientation fix.** The pipeline stores every amplicon forward-primer-first (steps A9 / B8), but `Extend reference.Rmd` had no such step, so any record taken from a minus-strand GenBank entry went in reverse-complemented. DADA2's `assignTaxonomy()` matches the forward strand only by default, so reads of those species could not be assigned to them. The Aug 2026 release (Zenodo v2026.08) contains 52 such 12SV5 records, all from the September recovery; the October recoveries raised the count to 355 12SV5 and 12 trnL records before this was found. All were reverse-complemented in v2026.10, 20 exact duplicates that this exposed were removed, and two 437 bp *Pseudophycis breviuscula* amplicons (the reference range is 71–180 bp) were dropped. `Extend reference.Rmd` now orients every amplicon and rejects out-of-range lengths, and the QC gate has a new check (9) that fails on any record starting with the reverse primer. `data/outputs/coverage-recheck/ORIENTATION_fix_Oct2026.csv` lists every record affected and whether it was in the deposited release. **(4) Two mislabelled trnL c/d records and line endings.** `JX390727.1` (labelled *Murraya koenigii*, in fact a peanut sequence) and `GU595140.1` (labelled *Secale cereale*, in fact rice) were in the trnL c/d database, including the 2026.08 release; both were removed (4,522 → 4,520 records, no species lost) and are now omitted by accession in step CD7. The same accessions were already omitted from g/h through `Manual renaming.csv`, but those rows are keyed to g/h amplicon sequences and cannot match c/d. Separately, the 2026.08 trnL g/h and 12SV5 files had Windows (CRLF) line endings, written by the laptop Extend workflow; all files are now LF, `Extend reference.Rmd` strips the carriage returns it would write, and the QC gate has a line-ending check (10).

### Remaining gaps (v2026.10)

Most of the gap reflects species with little or no public sequence for the marker, but not all of it. In October 2026 each of the 1,248 species not in the reference was checked against NCBI (`data/outputs/coverage-recheck/GAP_check_Oct2026.csv`, one row per species), applying the same rules the build uses: sequences up to 50 kb, both primer sites (allowing the build's mismatch tolerance), no UNVERIFIED records, no ambiguity codes, and an on-target phylum. The 58 animals that passed were then added to the reference, the unrecognised names were reviewed, and the plastome path was run for plants (see the v2026.10 note above), so the table shows what remains. Counts for species whose name was resolved are reported against the corrected name.

| Why the species is not in the reference | Plants (trnL) | Animals (12SV5) |
|---|---|---|
| Name not recognised by NCBI Taxonomy (obsolete name, unrecognised variety, or a valid name NCBI has no sequence for; 52 of the animals are invertebrates) | 17 | 78 |
| No record for the marker at NCBI under 50 kb (plants: and no complete chloroplast genome either) | 178 | 343 |
| Records exist, but none carries both primer sites (plants: includes 6 algae whose plastome lacks the trnL g/h sites) | 55 | 366 |
| Primer pair present only on UNVERIFIED records, or only on sequences with ambiguity codes | 15 | 2 |
| Non-vertebrate whose amplicon the Chordata filter excludes (crustaceans, a jellyfish) | — | 11 |
| Hybrid, not mapped to a single parent | — | 3 |
| Name resolved to an NCBI taxon, but no usable amplicon there | 6 | 11 |
| Valid amplicon at NCBI, not yet added to the reference | 0 | 0 |
| **Total missing** | **271** | **814** |

The 58 animals were 53 that NCBI files under the name on the food list and 5 that it files under a different name: *Arripis georgianus* (NCBI: *A. georgiana*), *Crenimugil crenilabis* (*Moolgarda crenilabis*), *Cypselurus agoo* (*Cheilopogon agoo*), *Meuschenia scaber* (*M. scabra*) and *Torpedo nobiliana* (*Tetronarce nobiliana*). Each of the five is mapped to its food-list name by `build_name`. The yak had hundreds of valid records, which collapse to a single distinct amplicon. For plants, the 233 species with no record under 50 kb were re-searched as complete chloroplast genomes (see the v2026.10 note): 73 were added, and the 178 now listed as having no record have neither a marker record nor a plastome.

**Plants without trnL sequences: 271 / 1,573 (17%)** (species not in the reference under either name)

| Source                                            | Missing species |
|---------------------------------------------------|-----------------|
| Lim, *Edible Medicinal and Non-Medicinal Plants*  | 95              |
| Milla 2020 Crop Origins                           | 95              |
| BP (contributor)                                  | 22              |
| Newton, *The Oldest Foods on Earth*               | 16              |
| Peters, *Edible Wild Plants of Subsaharan Africa* | 10              |
| JL (contributor)                                  | 9               |
| van Wyk, *Food Plants of the World*               | 5               |
| Other sources                                     | 19              |

**Animals without 12SV5 sequences: 814 / 2,121 (38%)** (species not in the reference under either name)

| Source                                   | Missing species |
|------------------------------------------|-----------------|
| FDA The Seafood List                     | 655             |
| SJ (contributor)                         | 50              |
| Halloran et al., *Edible Insects* †      | 42              |
| FAO Cultured Aquatic Species Fact Sheets | 20              |
| BP (contributor)                         | 20              |
| Newton, *The Oldest Foods on Earth*      | 9               |
| JL (contributor)                         | 8               |
| Other sources                            | 10              |

† Insects are not vertebrates and are therefore not targeted by the 12SV5 marker; these gaps are expected.

Full lists of species without sequence coverage are in `data/outputs/plants_missing_trnL.csv` and `data/outputs/animals_missing_12SV5.csv`, regenerated in October 2026 against the current v2026.10 build and counting `build_name` matches as covered. The counts above come from the same calculation.

------------------------------------------------------------------------

## Citation

The latest version of this database (2026.10) can be cited using the following: Brown, S., Petrone, B., Subramanian, A., Aqeel, A., Jiang, S., Superdock, D., & David, L. (2026). FoodSeq Reference Database (Version 2026.10) [Dataset]. Zenodo. <https://doi.org/10.5281/zenodo.23772049>

The files in this repository are release 2026.10, deposited on Zenodo at <https://zenodo.org/records/23772049>. To cite the database without specifying a version, use the concept DOI <https://doi.org/10.5281/zenodo.22756736>, which always resolves to the latest version. The previous release, 2026.08 (<https://doi.org/10.5281/zenodo.22756737>), is kept here in `archive/aug2026/`; it contains reverse-oriented 12SV5 records and two mislabelled trnL c/d records, corrected in 2026.10, so please cite and use 2026.10 instead.
