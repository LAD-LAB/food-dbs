---
editor_options: 
  markdown: 
    wrap: 72
---

# Coverage re-check — July 2026

Which food species could be **added** to the trnL / 12SV5 reference
databases, and which foods are missing from the target species list
entirely.

Generated 2026-07-20. Originally produced by nine standalone scripts;
those are now archived at
[`archive/code/coverage-recheck/`](../../../archive/code/coverage-recheck/)
for provenance. Re-run this with
[`code/Coverage recheck.Rmd`](../../../code/Coverage%20recheck.Rmd),
which consolidates them into one notebook with portable paths.

------------------------------------------------------------------------

## Why

The May 2026 build left real coverage gaps: 517 food plants with no trnL
sequence and 997 food animals with no 12SV5 sequence
(`../plants_missing_trnL.csv`, `../animals_missing_12SV5.csv`). Those
gaps were recorded when the reference was built — but GenBank grows, so
some are now closable. This re-check asks two questions:

1.  **Phase 1** — of the species already on the target list but missing
    a sequence, which now have marker data at NCBI?
2.  **Phase 2** — which foods are absent from
    `data/inputs/human-foods.csv` altogether?

------------------------------------------------------------------------

## Act on these

| File | Rows | What it is |
|------------------------|------------------------|------------------------|
| `CANDIDATES_plants_trnL_Aug2026.csv` | 203 | **Use this one.** The July plant candidates re-derived against the Aug 2026 build (see below) |
| `CANDIDATES_animals_12SV5_Aug2026.csv` | 209 | **Use this one.** The July animal candidates re-derived against the Aug 2026 build (see below) |
| `PHASE2_uncatalogued_foods_SOURCED.csv` | 40 | Foods **absent from `human-foods.csv` entirely**, each with a literature/database source. **Use this one, not the file below** (see caveat 5) |

### Re-derived against the Aug 2026 build (2026-08-25)

The July plant list was checked against the Aug 2026 reference,
resolving every name to its **current** NCBI name first. That mattered:
**82 of the 203 names (40%) were outdated** — `Allium porrum` →
`Allium ampeloprasum`, `Acacia senegal` → `Senegalia senegal`,
`Poncirus trifoliata` → `Citrus trifoliata` — so a plain name match
reports species as missing when they are already present.

| Outcome                                                               | n   |
|------------------------------------|------------------------------------|
| Added by the Aug 2026 build                                           | 83  |
| Queried; no record spans both trnL-g/h primer sites                   | 71  |
| Recoverable from a complete plastome, blocked by the 50 kb `SLEN` cap | 41  |
| Complete plastome exists, but not a land plant (kombu, wakame, nori)  | 4   |
| Complete plastome exists; phylum unresolved, check manually           | 3   |
| No usable trnL record at NCBI                                         | 1   |

Two things follow:

-   **`code/Extend reference.Rmd` cannot add any of the remaining 120.**
    All are on `human-foods.csv`, so the Aug rebuild already queried
    them through the same code path and got no amplicon. The extension
    notebook is for species *not* on the food list, or sequences
    published since the last build.
-   **41 named species are recoverable only by fixing the plastome gap**
    — the `0:50000[SLEN]` cap in `query_ncbi()` excludes complete
    chloroplast genomes, and a `trnL` text query would not find them
    anyway (96% of complete plastomes do not match the term). That is
    the single highest-yield outstanding fix.

The 209 July **animal** candidates were re-derived the same way. 60 of
209 (29%) had outdated names (`Trionyx sinensis` →
`Pelodiscus sinensis`, `Catla catla` → `Labeo catla`, `Manta birostris`
→ `Mobula birostris`).

| Outcome | n |
|------------------------------------|------------------------------------|
| Added by the Aug 2026 build | 43 |
| 12S records exist but none spans both V5 primer sites | 127 |
| **UNEXPLAINED — mitogenome exists and yields a valid amplicon, yet absent** | **39** |

**12S has neither structural gap that blocks trnL.** Vertebrate
mitogenomes are \~16 kb, so 100% fall inside the 50 kb cap, and 97%
match the free-text term `12S` (versus 4% of plastomes matching `trnL`).
Complete mitogenomes are reachable for animals.

Which makes the third row a genuine anomaly, not a known limitation.
Spot-checked *Dicentrarchus punctatus*, *Dentex dentex*, *Alosa alosa*
and *Aetobatus narinari*: each is on `human-foods.csv`, each has only
2–4 matching records (so the 500-record fetch cap is not implicated),
and each mitogenome yields a clean 135–138 bp V5 amplicon under
`find_primer_pair()` on test. *D. punctatus*'s amplicon is 1 edit from
the *D. labrax* record already in the reference. They should be in the
build and are not. **Root cause unknown — needs investigation before
these 39 are written off.**

## Supporting / reference

| File | Rows | What it is |
|------------------------|------------------------|------------------------|
| `SUPERSEDED_plants_trnL_available_Jul2026.csv` | 203 | **Superseded** by `CANDIDATES_plants_trnL_Aug2026.csv`. The original July list. Kept for provenance; do not act on it |
| `PHASE2_uncatalogued_foods_with_sequence_data.csv` | 39 | **Superseded** by `PHASE2_uncatalogued_foods_SOURCED.csv`. Kept only to document the unsourced starting point — it has no `source` column, so it does not meet the standard every row of `human-foods.csv` is held to. Do not act on it |
| `EXCLUDED_animals_invertebrates.csv` | 261 | Invertebrates removed from the 12SV5 candidates — documents what was filtered and why |
| `REVIEW_higher-rank-entries.csv` | 77 | Gap-list entries that are **not species** but genus/family/order names — need a separate decision (see caveat 3) |
| `STILL_missing_plants_trnL.csv` | 307 | Plants with still no trnL data anywhere |
| `STILL_missing_animals_12S.csv` | 460 | Animals with still no 12S data |
| `raw/` | — | Per-species NCBI hit counts, as queried. Provenance for everything above |

------------------------------------------------------------------------

## Method

One NCBI E-utilities `esearch` per species against `nucleotide`:

-   plants — `"<species>"[Organism] AND trnL`
-   animals — `"<species>"[Organism] AND 12S`
-   vertebrate filter —
    `"<species>"[Organism] AND 12S AND Vertebrata[Organism]`

`[Organism]` is hierarchical, so ANDing `Vertebrata[Organism]` returns 0
for non-vertebrates. That resolves vertebrate status in one query per
species rather than a two-step taxid lookup.

Rate-limited to \~2.5 req/sec (under NCBI's 3/sec no-API-key limit),
with `tool=`/`email=` set, and checkpointed so a run resumes rather than
re-querying. All 1,514 phase-1 species completed with zero errors.

------------------------------------------------------------------------

## Caveats — please read before acting on these lists

**1. A hit count \> 0 does NOT mean the amplicon is recoverable.** It
only means NCBI holds records *mentioning* trnL or 12S for that
organism. The record may not span the trnL-g/h or 12S-V5 primer binding
sites. The real test is the pipeline's in-silico PCR
(`find_primer_pair()`). **Treat every row as a candidate requiring
pipeline validation, not a confirmed addition** — expect attrition.

**2. 12SV5 is a vertebrate marker, so the raw animal counts were
misleading.** Every animal has a 12S rRNA gene, so invertebrates return
hits even though the V5 primers will not amplify them. Filtering cut the
animal candidates from 470 to **209** — 261 were invertebrates (shrimp,
crab, lobster, abalone, mussel, clam, barnacle, silkworm, honeybee). Use
`CANDIDATES_animals_12SV5_VERTEBRATES.csv`, not an unfiltered list.

**3. 77 gap-list entries are not species.** They are genus/family/order
names (`Siluriformes` catfish, `Cetacea` whales, `Clupeidae` herring,
`Mentha`, `Fragaria`). NCBI `[Organism]` searches expand to all
descendants, so these returned huge counts (Siluriformes: 3,540) that
say nothing about species-level coverage. They are excluded from the
candidate files and collected in `REVIEW_higher-rank-entries.csv` — each
needs a decision: pick a representative species, or keep as a
genus-level reference entry.

**4. Phase 2 is a targeted probe, not an exhaustive diff.** \~199
globally significant food species were checked against
`human-foods.csv`; 40 gaps were found. This is "gaps among species
probed", not a complete enumeration — no comprehensive external
food-species catalogue was used.

**5. Phase 2 candidates were sourced retrospectively, and not all of
them could be.** The probe list was assembled from general knowledge, so
the initial output (`PHASE2_uncatalogued_foods_with_sequence_data.csv`)
had **no source column** — unlike every one of the 3,777 rows in
`human-foods.csv`, which carries an attribution. Each candidate was
afterwards checked against the literature; results are in
`PHASE2_uncatalogued_foods_SOURCED.csv` with `source` and
`source_status` columns:

-   **28 confirmed** against a citable source. Several map onto sources
    already in the `human-foods.csv` vocabulary — FAO Cultured Aquatic
    Species Fact Sheets (American bullfrog), PFAF (hausa potato), PROTA
    (fluted pumpkin), FishBase.
-   **11 unsourced** — `source` deliberately left **blank**. These are
    plausibly eaten (`Bos indicus`, `Ovis canadensis`, `Bison bonasus`,
    `Capra ibex`, `Ovibos moschatus`, `Cervus nippon`,
    `Muntiacus reevesi`, `Odocoileus hemionus`, `Lepus americanus`,
    `Perdix perdix`), but no source naming the specific species was
    found. A blank field is deliberate: a plausible-looking but
    unverified citation would be worse than none.
-   **1 rejected** — `Carassius auratus`. FishBase states "edible but
    rarely eaten", and the original entry also conflated goldfish (*C.
    auratus*) with crucian carp (*C. carassius*). It is retained in the
    file with `source_status = REJECTED` so the decision is auditable.
    **Do not add it to `human-foods.csv`.**

Two confirmed rows carry qualifications, noted inline in
`source_status`: `Lepus europaeus` rests on a source naming "hare" only
at genus level, and the wild-harvested species (African buffalo,
warthog, Nile crocodile, the deer) are sourced to bushmeat/game-trade
literature — which documents consumption, but several are hunted under
licence or CITES restriction. That matters if the food list is meant to
reflect legal or commercial food supply rather than consumption in
general.

------------------------------------------------------------------------

## Phase 2 finding: the gap is systematic

Coverage is strong where the source list is strong and thin in one
specific area:

| Group                | Present | Gaps |
|----------------------|---------|------|
| Fish                 | 24/25   | 1    |
| Plants               | 93/96   | 3    |
| Bovids               | 1/10    | 9    |
| Birds (poultry/game) | 2/9     | 7    |
| Cervids              | 1/7     | 6    |
| Camelids             | 0/5     | 5    |
| Lagomorphs           | 0/3     | 3    |
| Game rodents         | 0/3     | 3    |

Every one of 76 staple crops checked was already present. But camel,
llama, alpaca, yak, zebu, domestic rabbit, hare, greater cane rat,
capybara, ostrich, guineafowl, species-level turkey, squab and most deer
were not. This tracks the target list's provenance — built heavily on
the FDA Seafood List and Western academic plant references, so deep on
seafood and crops, light on non-Western and non-industrial terrestrial
livestock.

**39 of the 40 gaps already have marker sequence data available** (only
*Plectranthus rotundifolius*, hausa potato, has none). After source
verification (caveat 5), **28 are confirmed with a citation, 11 remain
unsourced, and 1 was rejected** — so the number of defensible additions
is 28, not 40.

Some phase-2 entries are curation judgement calls rather than clear
omissions: `Bos indicus` is often treated as a *B. taurus* subspecies,
and turkey/red deer exist at genus level (`Meleagris`) or via a congener
(*Cervus canadensis*). The `status` column distinguishes
`GAP (fully absent)` from `GAP (species absent; genus present)`.
