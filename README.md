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
│   │   ├── human-foods.csv         # Curated list of food plant and animal species
│   │   ├── Manual renaming.csv     # Manual curation edits (omissions and renamings)
│   │   └── reference-additions.csv # Species list consumed by Extend reference.Rmd
│   └── outputs/
│       ├── dada2-compatible/       # Reference FASTAs for use with DADA2
│       │   ├── trnL/               # trnL databases (current + dated versions)
│       │   ├── 12Sv5/              # 12SV5 databases (current + dated versions)
│       │   ├── miscellaneous/      # Other marker databases
│       │   └── archive/            # Oct 2022 database versions
│       ├── qiime2-compatible/      # Reference FASTAs and TSVs for use with QIIME2
│       │   ├── trnL/
│       │   └── 12Sv5/
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

This repository contains the latest FoodSeq reference databases for use with the Edible Atlas FoodSeq Handbook. FASTA files containing reference sequences can be downloaded from this repository by going to data/outputs/ and selecting either dada2-compatible or qiime2-compatible FASTA files.

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
6.  Applies manual curation edits from `Manual renaming.csv`
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
    CLEAN -.->|trnL only| CUR["⑦b manual curation<br/>Manual renaming.csv"]
    CUR -.->|⑧ QC · ⑨ orientation · ⑩ dedup · add controls| OUT
    CTL -.->|add controls| OUT

    QC{"qc_reference_build.R<br/>controls · rank completeness · clade consistency<br/>accession integrity · count tolerance"}
    OUT --> QC
    QC -->|pass| SHIP["release"]
    QC -->|fail, exit 1| BLOCK["do not ship"]

    style CUR fill:#F5EEDC,stroke:#C9A96E,color:#5F4A1E
    style QC fill:#C8E0C9,stroke:#2C5F2D,color:#1E3A1F
    style BLOCK fill:#F2DCD6,stroke:#BC5138,color:#5F1E1E
    style SHIP fill:#CDE3CE,stroke:#2C5F2D,color:#1E5F1E
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
-   A species whose only NCBI record is a complete genome — Extend applies the same 50 kb length cap as the full pipeline (`query_ncbi()`), so a complete plastid/mitochondrial genome comes back empty even though sequence technically exists (see "Remaining gaps" below)

Extend never overwrites the reference it's extending — it always writes a new, separately-suffixed pair of files (e.g. `_Aug2026` → `_Aug2026_ext`), so a bad run is easy to discard and re-run. Before shipping the result, run the QC gate against it with the marker as a 4th argument, e.g.:

``` bash
Rscript code/qc_reference_build.R . _Aug2026_ext _Aug2026 trnL
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
    QC_PREVIOUS_SUFFIX <- "_Aug2026"   # last archived build, for the QC gate
    ```

4.  Run the **Install packages** chunk (section 1) once on first use

5.  Submit the three SLURM jobs (sections 2a–2c) and monitor their progress (section 3)

6.  Once jobs are complete, run **Part A** (trnL) and/or **Part B** (12SV5) sequentially

7.  The notebook's final section, **Verify before shipping**, runs the QC gate automatically and prints the result — 0 FAIL means the build is safe to archive. If `QC_PREVIOUS_SUFFIX` was left at its default (`"-"`), this step fails loudly rather than comparing against the wrong build; set it and re-run just that chunk.

Output files are written to `data/outputs/dada2-compatible/` and `data/outputs/qiime2-compatible/` as date-less files (`trnLGH.fasta`, not `trnLGH_Aug2026.fasta`) — once the QC gate passes, archive them under a new suffix before the next rebuild overwrites them in place.

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
IN_SUFFIX  <- "_Aug2026"    # the reference you're extending
OUT_SUFFIX <- "_Aug2026_ext" # what to write - never overwrites IN_SUFFIX
```

Run top to bottom. Species already in the reference are skipped automatically; species with no primer-spanning record are reported at the end rather than silently dropped.

#### 3. Verify before shipping

``` bash
Rscript code/qc_reference_build.R . _Aug2026_ext _Aug2026 trnL
```

The same gate the full pipeline uses. The 4th argument restricts it to the marker you extended — without it, the gate checks every marker under the new suffix and hard-fails on the ones you didn't touch. A clean run should show the record count rise by exactly what was added and every other check still pass.

#### 4. Promote the result

If the QC gate passes, rename the `OUT_SUFFIX` files to whatever the new canonical suffix is (e.g. `_Aug2026_ext` → `_Sep2026`) so they become what the rest of the repo — and the next `Extend` or `Coverage recheck` run — treats as current.

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

### Coverage over time

**trnL (plants)**

| Version  | Sequences | Unique taxa | Food plants covered | Coverage |
|----------|-----------|-------------|---------------------|----------|
| Oct 2022 | 1,402     | 807         | 716 / 1,570         | 46%      |
| 2025     | 1,991     | 1,169       | 1,060 / 1,570       | 68%      |
| May 2026 | 1,991     | 1,169       | 1,060 / 1,570       | 68%      |
| Aug 2026 | 2,026     | 1,195       | 1,078 / 1,573       | 69%      |

> Plant denominators exclude 7 genus-only entries in `human-foods.csv` (*Eucheuma, Fragaria, Gelidium, Gracilaria, Gurania, Mentha, Xanthosoma*) that aren't species-level binomials and can't be matched against a species-level reference — applied consistently across all four rows above, the same way the 12SV5 table below excludes non-vertebrates from every row.

**12SV5 (vertebrates)**

| Version  | Sequences | Unique taxa | Food animals covered | Coverage |
|----------|-----------|-------------|----------------------|----------|
| Oct 2022 | 57        | —           | —                    | —        |
| 2025     | 2,991     | 2,112       | 1,099 / 2,095        | 52%      |
| May 2026 | 3,390     | 2,337       | 1,168 / 2,095        | 56%      |
| Aug 2026 | 2,384     | 1,720       | 981 / 2,121          | 46%      |

> The Oct 2022 12SV5 database was a small food-filtered subset of the Schneider et al. database. The 2025 and May 2026 databases were built from scratch using the full pipeline.

> **Reading the Aug 2026 row:** `human-foods.csv` grew between the May and Aug 2026 builds (2,095 → 2,121 food-animal species, 1,570 → 1,573 food-plant species — see "Food species list" above), so the Aug 2026 coverage fraction is not directly comparable to earlier rows via raw denominator. Recomputed against the *current* list, May 2026 covers 1,060 / 1,573 (67%) plants and 891 / 2,121 (42%) animals — so Aug 2026 is a genuine improvement on both markers (69% and 46%), not the regression the raw 56% → 46% 12SV5 comparison would otherwise suggest. The Aug 2026 rebuild added a human host-taxon control, fixed 21 mislabelled off-target records, fixed a bare-genus accession bug, and closed 43 animal / 83 plant gaps found by the July 2026 coverage re-check — see `data/outputs/coverage-recheck/README.md` for the full account.

> **September 2026 recovery.** The figures above already include a follow-up fix: 51 animal species with a confirmed valid 12SV5 amplicon — present in the May 2026 build, or found by the July coverage re-check — had been silently dropped from the original Aug 2026 rebuild. Root cause: `query_ncbi()` queries species in batches of 5 joined by `OR`, then fetches only the first 500 combined results (`retmax_fetch`), so a rare species batched with a heavily-sequenced one (chicken, pig) can be crowded out of the fetch entirely with no error. `code/Extend reference.Rmd` was used to re-query and merge these back in (+53 sequences, +42 unique taxa; verified zero regressions against the pre-recovery build via `qc_reference_build.R`). One species (*Acanthurus gahhm*) remains unrecovered despite a confirmed valid amplicon — same unexplained character as the "39 unexplained" cases below.

### Remaining gaps (Aug 2026)

Coverage gaps primarily reflect species with limited public sequence data rather than pipeline limitations, with one exception: 41 plant species recoverable only from a complete plastome (blocked by the pipeline's 50 kb length cap — see "Extend, or rebuild?" above). The largest sources of missing species (based on those listed in `human-foods.csv`) are listed below.

**Plants without trnL sequences: 495 / 1,573 (31%)**

| Source                                            | Missing species |
|---------------------------------------------------|-----------------|
| Lim, *Edible Medicinal and Non-Medicinal Plants*  | 163             |
| Milla 2020 Crop Origins                           | 144             |
| BP (contributor)                                  | 85              |
| Newton, *The Oldest Foods on Earth*               | 21              |
| JL (contributor)                                  | 16              |
| van Wyk, *Food Plants of the World*               | 12              |
| GRIN Taxonomy for Plants                          | 11              |
| Peters, *Edible Wild Plants of Subsaharan Africa* | 11              |
| Other sources                                     | 32              |

**Animals without 12SV5 sequences: 1,140 / 2,121 (54%)**

| Source                                   | Missing species |
|------------------------------------------|-----------------|
| FDA The Seafood List                     | 900             |
| SJ (contributor)                         | 64              |
| Halloran et al., *Edible Insects* †      | 42              |
| NOAA Fisheries                           | 28              |
| FAO Cultured Aquatic Species Fact Sheets | 26              |
| BP (contributor)                         | 25              |
| Other sources                            | 55              |

† Insects are not vertebrates and are therefore not targeted by the 12SV5 marker; these gaps are expected.

Full lists of species without sequence coverage are provided in `data/outputs/plants_missing_trnL.csv` and `data/outputs/animals_missing_12SV5.csv` — **not yet regenerated against the Aug 2026 build** as of this update; the counts above were computed directly from the current `human-foods.csv` and the Aug 2026 reference FASTAs, not from those (still May-2026-era) files.

------------------------------------------------------------------------

## Reference

WEPSR: World Economic Plants: A Standard Reference, by John H. Wiersema and Blanca Leon.
