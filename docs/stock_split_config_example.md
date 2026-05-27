---
title: "Stock-splitting config — plain-English walkthrough"
subtitle: "For review of `data/config/stock_split_config_example.json`"
---

# Stock-splitting config — plain-English walkthrough

**File under review:** `data/config/stock_split_config_example.json`

**Purpose of this document:** the config that drives the pipeline's
`split-stocks` stage is written in JSON, which can be unfriendly to read.
This is a faithful plain-English translation of that same file so you can
review the choices (keywords, filters, output sheets) and flag anything
you would change. Every key in the JSON is covered below — nothing is
hidden or summarised away.

A short list of **questions for you** appears at the end. Those are the
places I'd most like your judgment.

---

## 1. TL;DR

The pipeline takes the FlyBase stock report for our input genes and slices
it into a workbook of prioritised Excel sheets. This config controls three
things:

1. **Keywords.** What counts as a "relevant" publication for a stock.
   Here: `sleep` and `circadian`.
2. **Filters.** Reusable yes/no rules about a stock (e.g. "is from
   Bloomington", "has no balancers", "has at least one keyword-hit
   paper"). The config defines 23 named filters.
3. **Combinations.** Each entry in the `combinations` list is a set of
   filters that get AND-ed together to produce one output sheet in the
   final workbook. This config defines 13 sheets.

The Ref++ / Ref+ / Ref- tiering used throughout is driven purely by whether
the keywords appear in publication titles or abstracts — there is no GPT
call or functional-validity check involved at this stage.

---

## 2. Search keywords and global settings

From the `settings` block of the JSON:

| Setting | Value in this config | Plain-English meaning |
|---|---|---|
| `relevantSearchTerms` | `sleep`, `circadian` | Keywords used to tier alleles into Ref++ / Ref+ / Ref-. Matching against publication titles and abstracts is **case-insensitive**. |
| `phenotypeSimilarityTargets` | `sleep`, `circadian` (each with `embedding_text` equal to the keyword) | Targets used **only** when OpenAI embedding-based phenotype similarity is enabled (off by default). Mirrors the keywords above. |
| `maxStocksPerGene` | `999999999999999` | Per-gene cap on how many stocks survive limiting. The value is so large in practice it imposes **no** cap. |
| `maxStocksPerAllele` | `9999999999999999` | Same idea, per allele. Also effectively **no** cap. |

See the questions at the end about whether the keyword set and the
effectively-unlimited caps are intentional.

---

## 3. Columns the pipeline reads or computes

These are the data columns that the filters reference. "Raw" means the
column comes straight from the FlyBase stock report; "Computed" means the
pipeline derives it from other columns.

| Column | Source | What it means |
|---|---|---|
| `genotype` | Raw | Canonical FlyBase stock genotype string. |
| `all_stock_constructs` | Raw / aggregated | Semicolon-separated list of unique construct/insertion symbols derived from FlyBase stock components (mostly FBti / FBtp symbols). |
| `collection` | Raw | FlyBase stock collection name (e.g. Bloomington, Vienna, Kyoto). |
| `PMID` | Raw / aggregated | Semicolon-separated list of PMIDs for papers associated with the allele, or `-` if none. |
| `RNAi` | **Computed** | Unified RNAi proxy. `True` if the genotype contains `UAS` or `RNAi`, **or** if a Vienna/VDRC stock matches a knockdown family (GD, KK, or VSH). |
| `UAS` | **Computed** | `True` if the genotype string contains `UAS` (case-insensitive). |
| `sgRNA` | **Computed** | `True` if the genotype string contains `sgRNA` (case-insensitive). |
| `Balancers` | **Computed** | Semicolon-separated balancer symbols derived from FBba stock components, or `-` if none. |
| `num_Balancers` | **Computed** | How many FBba (balancer) components are associated with the stock. |
| `multiple_insertions` | **Computed** | `True` if the stock has more than one unique transgenic construct (counted from `all_stock_constructs`). |
| `FBti_count` | **Computed** | How many FBti (transgenic insertion) elements are associated with the stock. |
| `attP_count` | **Computed** | How many times `attP` appears in the genotype. |
| `ALLELE_PAPER_RELEVANCE_SCORE` | **Computed** | 0, 1, or 2. Drives the Ref-/Ref+/Ref++ tier. See section 4. |
| `custom_stock` | **Computed** | `True` if this is a phenotype-backed reagent with no matching FBst stock (a "custom phenotype reagent"). Used by the `Custom Phenotype Reagent` filter. |

**Known gap (called out in the JSON itself):** filtering by chromosome arm
is not currently available. It would require an extra gene-mapping data
source that isn't wired in yet.

---

## 4. Filters in plain English

The config defines 23 named filters. Each filter is just a yes/no check
on one of the columns above. They are grouped here by intent so you can
scan by what they do, not by their JSON key order.

### 4a. Source / stock collection

| Filter | Plain-English rule | Underlying check |
|---|---|---|
| `Bloomington` | Stock comes from the Bloomington collection. | `collection` contains `Bloomington` |
| `Non-Bloomington` | Stock comes from any collection other than Bloomington. | `collection` does **not** contain `Bloomington` |

### 4b. Publication-evidence tier (the heart of Ref++/Ref+/Ref-)

`ALLELE_PAPER_RELEVANCE_SCORE` is set once per allele as follows:

- `0` → the allele has no publication references at all (= **Ref-**)
- `1` → the allele has at least one publication, but none of those papers
  have any of the keywords (`sleep`, `circadian`) in the title or
  abstract (= **Ref+**)
- `2` → the allele has at least one publication whose title or abstract
  contains a keyword (case-insensitive). This is the strongest tier
  (= **Ref++**)

| Filter | Plain-English rule | Underlying check |
|---|---|---|
| `Ref++` | Allele has at least one publication whose title or abstract mentions a keyword. | `ALLELE_PAPER_RELEVANCE_SCORE` equals `2` |
| `Ref+` | Allele has publications, but none mention a keyword in title/abstract. | `ALLELE_PAPER_RELEVANCE_SCORE` equals `1` |
| `Ref-` | Allele has no publication references at all. | `ALLELE_PAPER_RELEVANCE_SCORE` equals `0` |
| `Allele Has Paper Refs` | Allele has at least one publication. | `PMID` does **not** equal `-` |
| `Allele Has No Paper Refs` | Allele has no publications. | `PMID` equals `-` |

### 4c. Reagent class

| Filter | Plain-English rule | Underlying check |
|---|---|---|
| `UAS` | Stock contains a UAS-responsive component (GAL4-driven). | `UAS` equals `True` |
| `RNAi` | Stock matches the unified RNAi proxy (genotype has `UAS` or `RNAi`, or it's a VDRC GD/KK/VSH knockdown). | `RNAi` equals `True` |
| `AlleleOrInsertion` | Stock does **not** match the unified RNAi proxy (i.e. it's an allele or insertion, not an RNAi reagent). | `RNAi` equals `False` |
| `Has sgRNA` | Genotype includes `sgRNA` (CRISPR-based reagent). | `sgRNA` equals `True` |
| `No sgRNA` | Genotype does **not** include `sgRNA`. | `sgRNA` equals `False` |

### 4d. Balancer status

| Filter | Plain-English rule | Underlying check |
|---|---|---|
| `No Balancers` | Stock has no balancer components (no FBba). | `num_Balancers` equals `0` |
| `Has Balancers` | Stock has at least one balancer component. | `num_Balancers` does **not** equal `0` |

### 4e. Insertion details

The `attP40` and `attP2` filters are stricter than they look: they require
the stock to be a **single-construct** stock landed at that site (the
pipeline's `insertion_site` filter type checks the genotype text **and**
that `multiple_insertions` is `False`).

| Filter | Plain-English rule | Underlying check |
|---|---|---|
| `attP40` | Single-construct stock with its insertion at the attP40 landing site. | genotype contains `attP40` AND `multiple_insertions` is `False` |
| `attP2` | Single-construct stock with its insertion at the attP2 landing site. | genotype contains `attP2` AND `multiple_insertions` is `False` |
| `other_insertion` | Stock has more than one unique transgenic construct. | `multiple_insertions` equals `True` |
| `no_other_insertion` | Stock has only one unique transgenic construct. | `multiple_insertions` equals `False` |

### 4f. RNAi vector

These two filters do a **case-sensitive** substring check on the
genotype string.

| Filter | Plain-English rule | Underlying check |
|---|---|---|
| `VALIUM10` | Genotype mentions the VALIUM10 vector. | `genotype` contains `VALIUM10` |
| `VALIUM20` | Genotype mentions the VALIUM20 vector. | `genotype` contains `VALIUM20` |

### 4g. Split-Gal4

| Filter | Plain-English rule | Underlying check |
|---|---|---|
| `Split_AD` | Split-Gal4 activation-domain half (`.AD` in the genotype). | `genotype` contains `.AD` |
| `Split_DBD` | Split-Gal4 DNA-binding-domain half (`.DBD` in the genotype). | `genotype` contains `.DBD` |

### 4h. Phenotype-only reagent

| Filter | Plain-English rule | Underlying check |
|---|---|---|
| `Custom Phenotype Reagent` | Reagent has a phenotype association in FlyBase but no FBst stock record (i.e. it's not stockable; we surface it separately). | `custom_stock` equals `True` |

---

## 5. Output sheets (`combinations`)

Each entry in `combinations` produces one Excel sheet. All filters listed
in an entry are **AND-ed together** — a stock has to pass every one of
them to land in that sheet.

This config produces **13 sheets**. Twelve of them are a Bloomington-only,
non-CRISPR RNAi matrix; the thirteenth is an all-sources catch-all.

### 5a. The Bloomington RNAi matrix (12 sheets)

All 12 of these sheets share the same base filters:

> **Bloomington collection** AND **RNAi** AND **No sgRNA**

They differ along three axes — balancer status × insertion count ×
publication tier — yielding 2 × 2 × 3 = 12 sheets:

| # | Balancer status | Insertion count | Pub tier | Sheet contents |
|---|---|---|---|---|
| 1 | No balancers | Single insertion | Ref++ | Clean Bloomington RNAi stocks with keyword-hit papers. |
| 2 | No balancers | Single insertion | Ref+ | Clean Bloomington RNAi stocks with papers but no keyword hits. |
| 3 | No balancers | Single insertion | Ref- | Clean Bloomington RNAi stocks with no papers at all. |
| 4 | Has balancers | Single insertion | Ref++ | Balanced Bloomington RNAi stocks with keyword-hit papers. |
| 5 | Has balancers | Single insertion | Ref+ | Balanced Bloomington RNAi stocks with papers but no keyword hits. |
| 6 | Has balancers | Single insertion | Ref- | Balanced Bloomington RNAi stocks with no papers. |
| 7 | No balancers | Multiple insertions | Ref++ | Multi-insertion, unbalanced Bloomington RNAi stocks with keyword-hit papers. |
| 8 | No balancers | Multiple insertions | Ref+ | Multi-insertion, unbalanced Bloomington RNAi stocks with papers but no keyword hits. |
| 9 | No balancers | Multiple insertions | Ref- | Multi-insertion, unbalanced Bloomington RNAi stocks with no papers. |
| 10 | Has balancers | Multiple insertions | Ref++ | Multi-insertion, balanced Bloomington RNAi stocks with keyword-hit papers. |
| 11 | Has balancers | Multiple insertions | Ref+ | Multi-insertion, balanced Bloomington RNAi stocks with papers but no keyword hits. |
| 12 | Has balancers | Multiple insertions | Ref- | Multi-insertion, balanced Bloomington RNAi stocks with no papers. |

Reading the matrix:

- **Sheets 1, 4, 7, 10** are the four Ref++ buckets — these are the
  highest-priority Bloomington RNAi sheets and are the natural first stop
  for ordering.
- The "No balancers / Single insertion / Ref++" sheet (#1) is the
  cleanest possible reagent: one construct, no balancer headaches, and a
  paper that already mentions sleep or circadian.
- The 12-sheet matrix is **Bloomington-only**. There is no parallel
  matrix for Vienna, NIG-Fly, or Kyoto stocks — see questions.

### 5b. The catch-all sheet (13th sheet)

| # | Combination | Sheet contents |
|---|---|---|
| 13 | `RNAi` AND `No sgRNA` | All non-CRISPR RNAi stocks across every collection (Bloomington, Vienna, Kyoto, etc.), with no balancer or publication-tier restriction. Functions as a complete non-CRISPR RNAi inventory. |

### 5c. Filters that are defined but never used in a sheet

The following filters exist in the config but no `combinations` entry
references them: `Non-Bloomington`, `UAS`, `AlleleOrInsertion`,
`Has sgRNA`, `Allele Has Paper Refs`, `Allele Has No Paper Refs`,
`attP40`, `attP2`, `VALIUM10`, `VALIUM20`, `Split_AD`, `Split_DBD`,
`Custom Phenotype Reagent`.

They're functional — they could be slotted into a new combination — but
they don't currently generate any output sheet. See questions.

---

## 6. Domain glossary

Concise definitions for the terms used above. Most of these are standard
fly-community vocabulary; a few (`Ref++`, `custom_stock`) are
in-house to this pipeline.

- **Ref++ / Ref+ / Ref-** — the pipeline's three publication-evidence
  tiers for an allele:
  - **Ref++** = at least one paper whose title/abstract contains a
    configured keyword (case-insensitive).
  - **Ref+** = papers exist, but none of them contain a keyword in
    title/abstract.
  - **Ref-** = no papers at all.
- **UAS** — Upstream Activating Sequence. A GAL4-responsive construct;
  needs to be crossed to a GAL4 driver to be expressed.
- **RNAi proxy** (as used here) — `True` if the genotype contains
  `UAS` or `RNAi`, **or** the stock is from a Vienna/VDRC knockdown
  family (GD, KK, VSH).
- **VDRC GD / KK / VSH** — three Vienna Drosophila Resource Center
  RNAi-knockdown library families. The pipeline treats Vienna stocks
  from any of these families as RNAi reagents even when the genotype
  text alone doesn't make that obvious.
- **VALIUM10 / VALIUM20** — TRiP RNAi vectors. VALIUM20 is the
  shmiR-based vector commonly used for tissue-specific knockdown;
  VALIUM10 is the earlier-generation vector.
- **attP / attP40 / attP2** — phiC31 landing sites. `attP40` (chr 2L) and
  `attP2` (chr 3L) are the two workhorse sites for transgene insertion.
  In this config the `attP40` and `attP2` filters require a
  **single-construct** stock at that site.
- **Split-Gal4 `.AD` / `.DBD`** — the two halves of a Split-Gal4 driver.
  The activation domain (`.AD`) and DNA-binding domain (`.DBD`) are
  carried on separate transgenes and reconstitute GAL4 only where both
  expression patterns overlap.
- **sgRNA** — single-guide RNA. The presence of `sgRNA` in the genotype
  text marks the stock as a CRISPR reagent. The current config
  explicitly excludes these (every sheet uses `No sgRNA`).
- **Balancer** — a chromosome that suppresses recombination and carries
  dominant markers. Counted via FBba components in the stock.
- **FlyBase ID prefixes used here:**
  - `FBst` — a FlyBase stock record.
  - `FBal` — an allele.
  - `FBtp` — a transgenic construct (the plasmid-level entity).
  - `FBti` — a transgenic insertion (a construct landed at a site).
  - `FBba` — a balancer.
  - `FBgn` — a gene.
- **Custom phenotype reagent** — a reagent that has a phenotype
  annotation in FlyBase's `genotype_phenotype_data` but no associated
  FBst stock record, so it cannot be ordered as-is. The pipeline keeps
  these visible via the `custom_stock` flag rather than dropping them.

---

## 7. Questions for you

These are the spots where I'd like your judgment most. The pipeline will
do whatever this config tells it to — your call on what we want it to do.

1. **Keyword set.** This config uses `sleep` and `circadian` only. The
   GUI-design notes suggest a default of `sleep`, `rhythm`, `circadian`,
   `locomotor`. Should we add `rhythm` and `locomotor` here too, or are
   they too noisy / too broad for this particular screen?

2. **Per-gene / per-allele caps.** Both caps are set to absurdly large
   numbers (`9.99×10¹⁴` and `9.99×10¹⁵`), which is the same as having
   **no cap at all**. Is the intent really "give me everything" for this
   screen? If not, what would a reasonable cap per gene look like
   (e.g. 5? 10?).

3. **Bloomington-only matrix.** The 12 detailed sheets are all
   Bloomington. Vienna / Kyoto / NIG-Fly stocks only appear via the 13th
   catch-all sheet. Is that intentional (Bloomington is the primary
   order source for the lab) or is there a Vienna matrix you'd want as
   well, especially for genes with no Bloomington RNAi line?

4. **Unused filters.** These filters are defined but unused in
   `combinations`: `Non-Bloomington`, `UAS`, `AlleleOrInsertion`,
   `Has sgRNA`, `Allele Has Paper Refs`, `Allele Has No Paper Refs`,
   `attP40`, `attP2`, `VALIUM10`, `VALIUM20`, `Split_AD`, `Split_DBD`,
   `Custom Phenotype Reagent`. Are any of these worth turning into a
   sheet? For example, do you want dedicated `attP40` vs `attP2`
   sub-sheets, or a UAS-driver sheet, or a Split-Gal4 sheet?

5. **No-sgRNA filter on everything.** Every combination filters out
   CRISPR (`No sgRNA`). Is that the right default for your screen, or
   should there be a parallel CRISPR-only set of sheets?

6. **Multi-insertion stocks.** The matrix gives multi-insertion stocks
   their own row of sheets, but no extra prioritisation between them.
   Are multi-insertion lines roughly equivalent to single-insertion
   ones for ordering, or should they be deprioritised (and therefore
   maybe just folded into one combined sheet)?

7. **Chromosome-arm filtering** is currently unavailable in the
   pipeline. Would having that be valuable enough to invest in (e.g.
   for picking lines for a chr-3 screen), or is it nice-to-have?

8. **Anything else** about how stocks are bucketed, prioritised, or
   labelled — please flag it. The whole point of this document is to
   make those choices visible.
