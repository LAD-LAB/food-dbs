# Coverage re-check — July 2026

Which food species could be **added** to the trnL / 12SV5 reference databases, and which
foods are missing from the target species list entirely.

Generated 2026-07-20. Scripts that produce everything here live in
[`code/coverage-recheck/`](../../../code/coverage-recheck/).

---

## Why

The May 2026 build left real coverage gaps: 517 food plants with no trnL sequence and
997 food animals with no 12SV5 sequence (`../plants_missing_trnL.csv`,
`../animals_missing_12SV5.csv`). Those gaps were recorded when the reference was built —
but GenBank grows, so some are now closable. This re-check asks two questions:

1. **Phase 1** — of the species already on the target list but missing a sequence, which
   now have marker data at NCBI?
2. **Phase 2** — which foods are absent from `data/inputs/human-foods.csv` altogether?

---

## Act on these

| File | Rows | What it is |
|---|---|---|
| `CANDIDATES_plants_trnL_available.csv` | 203 | Food plants on the target list, missing from the reference, that **now have trnL records** at NCBI |
| `CANDIDATES_animals_12SV5_VERTEBRATES.csv` | 209 | Same for 12SV5, **filtered to true vertebrates** (see caveat 2) |
| `PHASE2_uncatalogued_foods_SOURCED.csv` | 40 | Foods **absent from `human-foods.csv` entirely**, each with a literature/database source. **Use this one, not the file below** (see caveat 5) |

## Supporting / reference

| File | Rows | What it is |
|---|---|---|
| `PHASE2_uncatalogued_foods_with_sequence_data.csv` | 39 | **Superseded** by `PHASE2_uncatalogued_foods_SOURCED.csv`. Kept only to document the unsourced starting point — it has no `source` column, so it does not meet the standard every row of `human-foods.csv` is held to. Do not act on it |
| `EXCLUDED_animals_invertebrates.csv` | 261 | Invertebrates removed from the 12SV5 candidates — documents what was filtered and why |
| `REVIEW_higher-rank-entries.csv` | 77 | Gap-list entries that are **not species** but genus/family/order names — need a separate decision (see caveat 3) |
| `STILL_missing_plants_trnL.csv` | 307 | Plants with still no trnL data anywhere |
| `STILL_missing_animals_12S.csv` | 460 | Animals with still no 12S data |
| `raw/` | — | Per-species NCBI hit counts, as queried. Provenance for everything above |

---

## Method

One NCBI E-utilities `esearch` per species against `nucleotide`:

- plants — `"<species>"[Organism] AND trnL`
- animals — `"<species>"[Organism] AND 12S`
- vertebrate filter — `"<species>"[Organism] AND 12S AND Vertebrata[Organism]`

`[Organism]` is hierarchical, so ANDing `Vertebrata[Organism]` returns 0 for non-vertebrates.
That resolves vertebrate status in one query per species rather than a two-step taxid lookup.

Rate-limited to ~2.5 req/sec (under NCBI's 3/sec no-API-key limit), with `tool=`/`email=`
set, and checkpointed so a run resumes rather than re-querying. All 1,514 phase-1 species
completed with zero errors.

---

## Caveats — please read before acting on these lists

**1. A hit count > 0 does NOT mean the amplicon is recoverable.**
It only means NCBI holds records *mentioning* trnL or 12S for that organism. The record may
not span the trnL-g/h or 12S-V5 primer binding sites. The real test is the pipeline's
in-silico PCR (`find_primer_pair()`). **Treat every row as a candidate requiring pipeline
validation, not a confirmed addition** — expect attrition.

**2. 12SV5 is a vertebrate marker, so the raw animal counts were misleading.**
Every animal has a 12S rRNA gene, so invertebrates return hits even though the V5 primers
will not amplify them. Filtering cut the animal candidates from 470 to **209** — 261 were
invertebrates (shrimp, crab, lobster, abalone, mussel, clam, barnacle, silkworm, honeybee).
Use `CANDIDATES_animals_12SV5_VERTEBRATES.csv`, not an unfiltered list.

**3. 77 gap-list entries are not species.**
They are genus/family/order names (`Siluriformes` catfish, `Cetacea` whales, `Clupeidae`
herring, `Mentha`, `Fragaria`). NCBI `[Organism]` searches expand to all descendants, so
these returned huge counts (Siluriformes: 3,540) that say nothing about species-level
coverage. They are excluded from the candidate files and collected in
`REVIEW_higher-rank-entries.csv` — each needs a decision: pick a representative species, or
keep as a genus-level reference entry.

**4. Phase 2 is a targeted probe, not an exhaustive diff.**
~199 globally significant food species were checked against `human-foods.csv`; 40 gaps were
found. This is "gaps among species probed", not a complete enumeration — no comprehensive
external food-species catalogue was used.

**5. Phase 2 candidates were sourced retrospectively, and not all of them could be.**
The probe list was assembled from general knowledge, so the initial output
(`PHASE2_uncatalogued_foods_with_sequence_data.csv`) had **no source column** — unlike every
one of the 3,777 rows in `human-foods.csv`, which carries an attribution. Each candidate was
afterwards checked against the literature; results are in
`PHASE2_uncatalogued_foods_SOURCED.csv` with `source` and `source_status` columns:

- **28 confirmed** against a citable source. Several map onto sources already in the
  `human-foods.csv` vocabulary — FAO Cultured Aquatic Species Fact Sheets (American bullfrog),
  PFAF (hausa potato), PROTA (fluted pumpkin), FishBase.
- **11 unsourced** — `source` deliberately left **blank**. These are plausibly eaten
  (`Bos indicus`, `Ovis canadensis`, `Bison bonasus`, `Capra ibex`, `Ovibos moschatus`,
  `Cervus nippon`, `Muntiacus reevesi`, `Odocoileus hemionus`, `Lepus americanus`,
  `Perdix perdix`), but no source naming the specific species was found. A blank field is
  deliberate: a plausible-looking but unverified citation would be worse than none.
- **1 rejected** — `Carassius auratus`. FishBase states "edible but rarely eaten", and the
  original entry also conflated goldfish (*C. auratus*) with crucian carp (*C. carassius*).
  It is retained in the file with `source_status = REJECTED` so the decision is auditable.
  **Do not add it to `human-foods.csv`.**

Two confirmed rows carry qualifications, noted inline in `source_status`: `Lepus europaeus`
rests on a source naming "hare" only at genus level, and the wild-harvested species (African
buffalo, warthog, Nile crocodile, the deer) are sourced to bushmeat/game-trade literature —
which documents consumption, but several are hunted under licence or CITES restriction. That
matters if the food list is meant to reflect legal or commercial food supply rather than
consumption in general.

---

## Phase 2 finding: the gap is systematic

Coverage is strong where the source list is strong and thin in one specific area:

| Group | Present | Gaps |
|---|---|---|
| Fish | 24/25 | 1 |
| Plants | 93/96 | 3 |
| Bovids | 1/10 | 9 |
| Birds (poultry/game) | 2/9 | 7 |
| Cervids | 1/7 | 6 |
| Camelids | 0/5 | 5 |
| Lagomorphs | 0/3 | 3 |
| Game rodents | 0/3 | 3 |

Every one of 76 staple crops checked was already present. But camel, llama, alpaca, yak,
zebu, domestic rabbit, hare, greater cane rat, capybara, ostrich, guineafowl, species-level
turkey, squab and most deer were not. This tracks the target list's provenance — built
heavily on the FDA Seafood List and Western academic plant references, so deep on seafood
and crops, light on non-Western and non-industrial terrestrial livestock.

**39 of the 40 gaps already have marker sequence data available** (only *Plectranthus
rotundifolius*, hausa potato, has none). After source verification (caveat 5), **28 are
confirmed with a citation, 11 remain unsourced, and 1 was rejected** — so the number of
defensible additions is 28, not 40.

Some phase-2 entries are curation judgement calls rather than clear omissions: `Bos indicus`
is often treated as a *B. taurus* subspecies, and turkey/red deer exist at genus level
(`Meleagris`) or via a congener (*Cervus canadensis*). The `status` column distinguishes
`GAP (fully absent)` from `GAP (species absent; genus present)`.
