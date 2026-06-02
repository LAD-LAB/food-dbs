# food-dbs

Reference database pipeline for dietary metabarcoding using the trnL (plant) and 12SV5 (vertebrate) markers. Builds DADA2- and QIIME2-compatible reference databases by combining RefSeq sequences downloaded locally with remote GenBank queries via the NCBI API.

---

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
│   │   └── query_ncbi_accession.R  # Resolve accessions to taxon IDs
│   ├── Descriptive-statistics.Rmd  # Summary statistics for built databases
│   ├── Extend reference.Rmd        # Add sequences to an existing reference
│   ├── Parse ecoPCR.Rmd            # Parse ecoPCR output
│   ├── SQL to reference.Rmd        # Alternative taxonomy approach
│   ├── Taxa names.Rmd              # Taxon name handling
│   └── tree-building.Rmd           # Phylogenetic tree construction
│
├── data/
│   ├── inputs/
│   │   ├── human-foods.csv         # Curated list of food plant and animal species
│   │   └── Manual renaming.csv     # Manual curation edits (omissions and renamings)
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

---

## How the pipeline works

The pipeline builds two reference databases in parallel, following the same steps for each:

### Data sources
Sequences are drawn from two sources and merged:
- **RefSeq** (local): Full plastid (trnL) or mitochondrial (12SV5) genomes downloaded from the NCBI FTP server, filtered to one sequence per species, then in silico amplified using the target primer pair
- **GenBank** (remote): Sequences queried directly from NCBI's nucleotide database using `query_ncbi()`, which searches for the target marker across all species in the food species list

### Species lists
- **trnL**: All species in `human-foods.csv` where `category == "plant"`
- **12SV5**: All species in `human-foods.csv` where `category == "animal"`

### Taxonomy
Accession numbers are mapped to full taxonomic lineages using a local `taxonomizr` SQL database (a mirror of NCBI taxonomy). Any accessions added to NCBI after the SQL database was built are resolved automatically using `query_ncbi_accession()`.

### Processing steps
For each marker, the pipeline:
1. Downloads RefSeq genomes and queries GenBank (via SLURM jobs)
2. Applies in silico PCR with `find_primer_pair()` to extract the target amplicon
3. Removes unverified sequences and sequences with degenerate nucleotides
4. Looks up full taxonomy for all accessions
5. Applies manual curation edits from `Manual renaming.csv`
6. Standardises sequence orientation
7. Deduplicates (identical sequences from the same taxon collapsed to one, with RefSeq accessions preferred)
8. Saves DADA2- and QIIME2-compatible output files

### Primers
| Marker | Forward | Reverse |
|--------|---------|---------|
| trnL (plant) | `GGGCAATCCTGAGCCAA` (trnLg) | `CCATTGAGTCTCTGCACCTATC` (trnLh) |
| 12SV5 (vertebrate) | `TAGAACAGGCTCCTCTAG` (V5F) | `TTAGATACCCCACTATGC` (V5R) |

---

## Database coverage

### Food species list

The pipeline targets all species in `human-foods.csv`, a curated list of 3,777 food species assembled from 32 literature and database sources.

| Category | Species |
|---|---|
| Plants | 1,577 |
| Animals (species-level) | 2,095 |
| Animals (genus/family-level entries) | 70 |
| Fungi | 33 |
| Bacteria | 2 |
| **Total** | **3,777** |

### Coverage over time

**trnL (plants)**

| Version | Sequences | Unique taxa | Food plants covered | Coverage |
|---|---|---|---|---|
| Oct 2022 | 1,402 | 807 | 716 / 1,577 | 45% |
| 2025 | 1,991 | 1,169 | 1,060 / 1,577 | 67% |
| May 2026 | 1,991 | 1,169 | 1,060 / 1,577 | 67% |

**12SV5 (vertebrates)**

| Version | Sequences | Unique taxa | Food animals covered | Coverage |
|---|---|---|---|---|
| Oct 2022 | 57 | — | — | — |
| 2025 | 2,991 | 2,112 | 1,099 / 2,095 | 52% |
| May 2026 | 3,390 | 2,337 | 1,168 / 2,095 | 56% |

> The Oct 2022 12SV5 database was a small food-filtered subset of the Schneider et al. database. The 2025 and May 2026 databases were built from scratch using the full pipeline.

### Remaining gaps (May 2026)

Coverage gaps primarily reflect species with limited public sequence data rather than pipeline limitations. The largest sources of missing species are listed below.

**Plants without trnL sequences: 517 / 1,577 (33%)**

| Source | Missing species |
|---|---|
| Lim, *Edible Medicinal and Non-Medicinal Plants* | 177 |
| Milla 2020 Crop Origins | 146 |
| BP (contributor) | 86 |
| Newton, *The Oldest Foods on Earth* | 20 |
| JL (contributor) | 18 |
| van Wyk, *Food Plants of the World* | 12 |
| Peters, *Edible Wild Plants of Subsaharan Africa* | 11 |
| GRIN Taxonomy for Plants | 11 |
| Other sources | 36 |

**Animals without 12SV5 sequences: 927 / 2,095 (44%)**

| Source | Missing species |
|---|---|
| FDA The Seafood List | 755 |
| SJ (contributor) | 53 |
| Halloran et al., *Edible Insects* † | 42 |
| FAO Cultured Aquatic Species Fact Sheets | 24 |
| BP (contributor) | 22 |
| Other sources | 31 |

† Insects are not vertebrates and are therefore not targeted by the 12SV5 marker; these gaps are expected.

Full lists of species without sequence coverage are provided in `data/outputs/plants_missing_trnL.csv` and `data/outputs/animals_missing_12SV5.csv`.

---

## Getting started

### Prerequisites

**R packages:**
```r
# Bioconductor
BiocManager::install(c("Biostrings", "ShortRead"))

# CRAN
install.packages(c("taxonomizr", "tidyverse", "rentrez", "remotes"))

# GitHub
remotes::install_github("ammararuby/MButils")
```

**NCBI API key** (recommended): Create an NCBI account, go to Account Settings → API Key Management, and copy your key. In R, run:
```r
rentrez::set_entrez_key("your_key_here")
```
This increases the NCBI query rate limit from 3 to 10 requests per second. NCBI also recommends running large queries on weekends or between 9 PM and 5 AM EST on weekdays.

### Large file setup (cluster recommended)

Two large files must be obtained before running the pipeline. These are best downloaded on an HPC cluster due to their size and download time:

| File | Size | How to obtain |
|------|------|---------------|
| RefSeq plastid FASTA | ~15–20 GB uncompressed | NCBI FTP: `ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/` |
| RefSeq mitochondrial FASTA | ~5–10 GB uncompressed | NCBI FTP: `ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/` |
| `accessionTaxa.sql` | ~70 GB | Built with `taxonomizr::prepareDatabase()` |

SLURM job scripts for all three are generated and submitted automatically by the pipeline (sections 2a–2c of the Rmd).

### Running the pipeline

1. Clone this repository:
   ```bash
   git clone https://github.com/LAD-LAB/food-dbs.git
   ```

2. Open `foodseq_reference_pipeline.Rmd` in RStudio (or via RStudio Server on Open OnDemand)

3. Update the paths in the **Configuration** chunk (section 0) to match your environment:
   ```r
   SCRATCH  <- "/scratch/your_username"
   REPO_DIR <- "/path/to/food-dbs"
   SQL_PATH <- "/path/to/accessionTaxa.sql"
   ```

4. Run the **Install packages** chunk (section 1) once on first use

5. Submit the three SLURM jobs (sections 2a–2c) and monitor their progress (section 3)

6. Once jobs are complete, run **Part A** (trnL) and/or **Part B** (12SV5) sequentially

Output files are written to `data/outputs/dada2-compatible/` and `data/outputs/qiime2-compatible/`.

---

## Reference

WEPSR: World Economic Plants: A Standard Reference, by John H. Wiersema and Blanca Leon.
