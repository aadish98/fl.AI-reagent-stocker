# fly_stocker_v2: Practical Pipeline Guide

This guide explains how to go from a gene list to a curated set of candidate fly stocks with supporting literature.

## What You Get

After a full run, you will have:

- a stock-level table linked to references
- organized output sheets based on your filtering rules
- a reference table narrowed to papers relevant to selected stocks
- optional functional-validation signals for Ref++ sheets

## End-to-End Flow (Purpose and Method)

### Flowchart 1: High-Level Overview

Use this when you want the one-minute story of what goes in, what the pipeline does, and what comes out.

![High-level end-to-end pipeline flowchart](pipeline-flowchart-high-level-standalone.png)

### Flowchart 2: Rules and Evidence Decisions

Use this when you want to understand how prioritization rules are applied and how `Ref++` decisions are handled.

![Rules and evidence decision flowchart](pipeline-flowchart-rules-and-tools-standalone.png)

## Before You Run

- Run commands from the project directory (the folder with `cli.py` and `data/JSON/`).
- Install dependencies:

```bash
pip install -r requirements.txt
```

- Put your input CSV files in one folder (example: `./gene_lists`).
- Optional environment variables:
  - `OPENAI_API_KEY` (needed if you want GPT-based validation)
  - `NCBI_API_KEY` (helps with PubMed rate limits)
  - `UNPAYWALL_TOKEN` (optional, can improve full-text retrieval)

## 5-Minute Quick Start

### Option A (recommended): one command

```bash
python -m fly_stocker_v2.cli run-full-pipeline ./gene_lists --config ./my_split_config.json --test-log
```

### Option B: two-step workflow

```bash
python -m fly_stocker_v2.cli get-allele-refs ./gene_lists
python -m fly_stocker_v2.cli split-stocks ./gene_lists/Stocks --test-log
```

Default config path:

- `data/JSON/stock_split_config_example.json`

---

## Step 1: Build Stocks + References (`get-allele-refs`)

This step reads your gene lists and creates an Excel workbook that links genes, stocks, and publications.

### Command

```bash
python -m fly_stocker_v2.cli get-allele-refs ./gene_lists
```

### Main output

- `./gene_lists/Stocks/aggregated_bdsc_stock_refs.xlsx`

### What is inside

- `Stocks`: stock-centric rows with associated metadata
- `References`: publication-centric rows connected to stock entries

### Common options

- `--config`, `-c`: path to your split config JSON (loads keyword terms)
- `--skip-fbgnid-conversion`: use this if your input already has FBgn IDs
- `--run-validation`: run GPT validation during Step 1 (off by default)
- `--soft-run`: estimate GPT calls without making them (only meaningful with validation)
- `--test-log`: write GPT query logs (only when validation runs)

---

## Step 2: Organize and Prioritize (`split-stocks`)

This step takes the Step 1 workbook, applies your config logic, and writes organized output sheets.

### Command

```bash
python -m fly_stocker_v2.cli split-stocks ./gene_lists/Stocks --config ./my_split_config.json
```

### Main output folder

- `./gene_lists/Stocks/Organized Stocks/`

For each input workbook, output is `<input_name>_aggregated.xlsx`.

### What is inside

- `Contents`: overview of generated sheets
- `Sheet1`, `Sheet2`, ...: config-driven stock groupings
- `References`: only PMIDs cited by stocks that remain in output sheets
- `Stock Sheet by Gene` (when applicable): gene-grouped stock/reference summary

### Common options

- `--config`, `-c`: split configuration file
- `--quiet`, `-q`: less console output
- `--soft-run`: estimate validation calls and stop before GPT calls
- `--test-log`: write GPT query logs to `data/logs/GPT-Queries/`
- `--max-gpt-calls-per-stock`: cap GPT calls per stock during validation

---

## How Validation Is Applied

Validation is intentionally selective.

In Step 2, a stock-reference pair is considered for GPT validation only when it passes all of the following:

1. the stock survives filter and limit rules into a sheet whose combination includes `Ref++`
2. the reference is keyword-relevant (from `settings.relevantSearchTerms`)
3. full-text and pattern checks indicate the pair is eligible

If a stock is only present in non-`Ref++` sheets, it is not sent for GPT validation.

---

## Choosing Between Commands

- Use `run-full-pipeline` when you want a single end-to-end run.
- Use separate `get-allele-refs` then `split-stocks` when you want to inspect or reuse intermediate outputs.

### Full pipeline command

```bash
python -m fly_stocker_v2.cli run-full-pipeline ./gene_lists --config ./my_split_config.json --test-log
```

Equivalent package entry-point usage:

```bash
python -m fly_stocker_v2 get-allele-refs ./gene_lists --config ./my_split_config.json
python -m fly_stocker_v2 split-stocks ./gene_lists/Stocks --config ./my_split_config.json
```

---

## Config File in Plain Terms

The JSON config controls:

- which terms define “reference relevance” (`settings.relevantSearchTerms`)
- max stocks per gene (`settings.maxStocksPerGene`)
- max stocks per allele (`settings.maxStocksPerAllele`)
- reusable filter blocks (`filters`)
- output-sheet definitions (`combinations`)
- human-readable sheet descriptions (`filterDescriptions`)

Example config:

- `data/JSON/stock_split_config_example.json`

---

## Example Run-Through (inc, tim, vri)

Illustrative counts for genes `inc`, `tim`, `vri` with default config (keywords: `["sleep", "circadian"]`, maxStocksPerGene: 5, maxStocksPerAllele: 3).

### `get-allele-refs`

```
┌─────────────────────────────────────────────────────────────────────────┐
│ INPUT: gene_list.csv  (inc, tim, vri)                                   │
└────────────────────────────────┬────────────────────────────────────────┘
                                 │
                                 ▼
              ┌──────────────────────────────────┐
              │ 1. Gene → FBgn ID conversion     │
              │    inc → FBgn0036816              │
              │    tim → FBgn0014396              │
              │    vri → FBgn0016076              │
              └───────────────┬──────────────────┘
                              │
                              ▼
              ┌──────────────────────────────────┐
              │ 2. FBgn → BDSC Stocks            │
              │    (via stockgenes.csv)           │
              │    inc: 12 stocks, 8 alleles      │
              │    tim: 35 stocks, 15 alleles     │
              │    vri: 10 stocks, 6 alleles      │
              │    ────────────────────────        │
              │    57 stocks, 29 alleles total     │
              └───────────────┬──────────────────┘
                              │
               ┌──────────────┴──────────────┐
               ▼                             ▼
┌─────────────────────────────┐  ┌──────────────────────────┐
│ 3. Alleles → References     │  │ 4. Flag & Metadata       │
│    FlyBase entity_publication│  │    UAS: 18 stocks        │
│    → FBrf → PMID            │  │    sgRNA: 3 stocks       │
│    (papers only, ≤50 genes) │  │    + chromosome info     │
│    104 unique PMIDs         │  │    + stock comments      │
│    17 FBrfs without PMID    │  └──────────────────────────┘
└──────────────┬──────────────┘
               │
               ▼
┌─────────────────────────────┐
│ 5. PubMed Metadata Fetch    │
│    104 PMIDs → Entrez API   │
│    → title, abstract,       │
│      authors, journal, date │
└──────────────┬──────────────┘
               │
               ▼
┌─────────────────────────────┐
│ 6. Keyword Scoring          │
│    "sleep" OR "circadian"   │
│    in title/abstract?       │
│    38 match, 66 don't       │
│    → per-stock counts       │
└──────────────┬──────────────┘
               │
               ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ OUTPUT: aggregated_bdsc_stock_refs.xlsx                                  │
│   Stocks sheet (57 rows) ── stock, genotype, alleles, UAS, sgRNA,       │
│                              PMIDs, keyword counts, keyword boolean      │
│   References sheet (104 rows) ── PMID, title, genes, keyword boolean    │
└─────────────────────────────────────────────────────────────────────────┘
```

### `split-stocks`

```
┌─────────────────────────────────────────────────────────────────────────┐
│ INPUT: aggregated_bdsc_stock_refs.xlsx (57 stocks, 104 refs)            │
│      + stock_split_config_example.json (18 filters, 15 combinations)    │
└────────────────────────────────┬────────────────────────────────────────┘
                                 │
                                 ▼
              ┌──────────────────────────────────────────┐
              │ 1. Compute Derived Columns               │
              │    Balancers:      14 have / 43 none      │
              │    multi-insert:    9 yes  / 48 no        │
              │    Relevance score:                       │
              │      Ref++ (score 2): 30    ── keyword hit│
              │      Ref+  (score 1): 19    ── refs only  │
              │      Ref−  (score 0):  8    ── no refs    │
              └───────────────────┬──────────────────────┘
                                  │
                                  ▼
              ┌──────────────────────────────────────────┐
              │ 2. Sort by Priority                      │
              │    score → keyword hits → total refs     │
              └───────────────────┬──────────────────────┘
                                  │
                                  ▼
              ┌──────────────────────────────────────────┐
              │ 3. Apply 15 Filter Combinations          │
              │                                          │
              │  Example (combo 1):                      │
              │    Bloomington ∩ UAS ∩ No sgRNA ∩        │
              │    No Balancers ∩ single insert ∩ Ref++  │
              │    57→18→15→11→9→6 match                 │
              │    after limits (5/gene, 3/allele): 5    │
              │    → Sheet1                              │
              │                                          │
              │  Each combo: AND filters → limit → sheet │
              │  Stocks used in earlier sheets excluded  │
              │  Empty combos produce no sheet           │
              └───────────────────┬──────────────────────┘
                                  │
                                  ▼
              ┌──────────────────────────────────────────┐
              │ 4. GPT Validation (Ref++ sheets only)    │
              │    13 Ref++ stocks, 28 (stock, PMID)     │
              │    tasks — fetch full text, check for    │
              │    allele mention, send to GPT           │
              │    → functional validity, phenotypes     │
              │    (capped at 5 calls/stock)             │
              └───────────────────┬──────────────────────┘
                                  │
                                  ▼
              ┌──────────────────────────────────────────┐
              │ 5. Filter References                     │
              │    Keep PMIDs cited by surviving stocks   │
              │    104 → 78 rows                         │
              └───────────────────┬──────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ OUTPUT: aggregated_bdsc_stock_refs_aggregated.xlsx                       │
│   Contents ──────── summary table, filter definitions, gene lists       │
│   Sheet1–Sheet15 ── filtered stock groups (35 unique stocks)            │
│   References ────── 78 rows, keyword column color-coded                 │
│   Stock Sheet by Gene ── Ref++ stocks grouped by gene, one row per      │
│                          (stock, keyword-hit PMID) pair                  │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Notes for Reliability and Reproducibility

- Full-text retrieval uses multiple sources and retries transient failures.
- Metadata and method caches are used to reduce repeated lookup work.
- Duplicate reference rows with the same PMID are merged with preference for richer identifiers (PMCID/DOI).
- Validation can short-circuit per stock once functional evidence is found.

---

## Useful Flags (Reference)

### `get-allele-refs`

- `input_dir` (required): folder containing your CSV files
- `--gene-col` (default: `flybase_gene_id`)
- `--input-gene-col` (default: `ext_gene`)
- `--batch-size`, `-b` (default: `50`)
- `--config`, `-c`
- `--skip-fbgnid-conversion`
- `--run-validation`
- `--soft-run`
- `--test-log`
- `--max-gpt-calls-per-stock` (default: `5`)

### `split-stocks`

- `input_dir` (required): folder with Step 1 Excel files (example: `./gene_lists/Stocks`)
- `--config`, `-c`
- `--quiet`, `-q`
- `--soft-run`
- `--test-log`
- `--max-gpt-calls-per-stock` (default: `5`)

### `run-full-pipeline`

- `input_dir` (required): folder containing gene-list CSV files
- shared options include `--config`, `--soft-run`, `--test-log`, `--max-gpt-calls-per-stock`
- pipeline-specific options include `--gene-col`, `--input-gene-col`, `--batch-size`, `--skip-fbgnid-conversion`, `--quiet`

---

## Quick Troubleshooting

- If validation does not run, check whether your output sheet is `Ref++` and whether `OPENAI_API_KEY` is set.
- If PMID coverage looks low, verify keyword terms in `settings.relevantSearchTerms`.
- If full-text retrieval is sparse, add `UNPAYWALL_TOKEN` and retry.
- If API usage is a concern, start with `--soft-run`.

