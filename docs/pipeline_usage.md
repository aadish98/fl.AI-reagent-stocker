# fl_ai_reagent_stocker: Practical Pipeline Guide

This guide explains how to go from a gene list to a curated set of candidate fly stocks with supporting literature, config-driven sheet organization, and optional AI validation.

## What You Get

After a full run, you will have:

- a stock-level workbook linked to references
- organized output sheets based on your JSON rules
- a reference table narrowed to papers relevant to selected stocks
- optional validation columns for `Ref++` output-sheet stocks

## Visual Overview

The maintained flowcharts live in `docs/pipeline_flowcharts.md` as editable
Mermaid diagrams:

- high-level command flow from input gene lists to final workbooks
- low-level evidence flow showing stock-level references, phenotype-sheet
  references, keyword buckets, embeddings, and validation decisions

## Before You Run

- Run commands from the project directory.
- Install dependencies:

```bash
pip install -r requirements.txt
```

- Put your input CSV files in one folder such as `./gene_lists`.
- Optional environment variables:
  - `OPENAI_API_KEY`
  - `NCBI_API_KEY`
  - `UNPAYWALL_TOKEN`
  - `OPENAI_MODEL`
  - `OPENAI_EMBEDDING_MODEL`

## 5-Minute Quick Start

### Option A: one command, normal mode (split + validate)

```bash
python -m fl_ai_reagent_stocker run-full-pipeline ./gene_lists --config ./my_split_config.json --mode normal --test-log
```

### Option B: one command, phenotype mode (Stage 1 + phenotype sheet)

```bash
python -m fl_ai_reagent_stocker run-full-pipeline ./gene_lists --config ./my_split_config.json --mode phenotype --similarity tiers
```

### Option C: explicit stages

```bash
python -m fl_ai_reagent_stocker find-stocks ./gene_lists
python -m fl_ai_reagent_stocker split-stocks ./gene_lists/Stocks --config ./my_split_config.json
python -m fl_ai_reagent_stocker phenotype-sheet ./gene_lists/Stocks --config ./my_split_config.json --similarity tiers
python -m fl_ai_reagent_stocker validate-stocks ./gene_lists/Stocks --config ./my_split_config.json --test-log
```

Default config path:

- `data/config/stock_split_config_example.json`

## Stage 1: Build Stocks + References (`find-stocks`)

This stage reads your gene lists and creates an Excel workbook that links genes, stocks, and publications.

### Command

```bash
python -m fl_ai_reagent_stocker find-stocks ./gene_lists
```

### Main output

- `./gene_lists/Stocks/aggregated_stock_refs.xlsx`

### What Stage 1 does

1. Convert gene symbols to FBgn IDs unless skipped.
2. Build gene -> allele mappings from `fbal_to_fbgn`.
3. Build allele -> construct mappings from `transgenic_construct_descriptions`.
4. Expand construct IDs to insertion IDs through `data/flybase/transgenic_insertions/fbtp_to_fbti.csv`.
5. Match allele, construct, and insertion IDs to stocks via `fbst_to_derived_stock_component.csv`.
6. Link matched components to references and PMIDs.
7. Score title/abstract keyword relevance using `settings.relevantSearchTerms`.

### Common options

- `--config`, `-c`
- `--skip-fbgnid-conversion`
- `--gene-col`
- `--input-gene-col`
- `--batch-size`

## Stage 2: Organize Stocks (`split-stocks`)

This stage takes the Stage 1 workbook, applies your config logic, and writes organized output sheets. It does not produce the phenotype sheet and does not run GPT validation.

### Command

```bash
python -m fl_ai_reagent_stocker split-stocks ./gene_lists/Stocks --config ./my_split_config.json
```

### Main output folder

- `./gene_lists/Stocks/Organized Stocks/`

For each input workbook, output is `<input_name>_aggregated.xlsx`.

### What is inside

- `Contents`
- `Sheet1`, `Sheet2`, ...
- `References`
- `Stock Sheet by Gene`

`split-stocks` no longer accepts `--soft-run`, `--OAI-embedding`,
`--simple-buckets`, or `--keyword-bucketing`. Use `phenotype-sheet` for the
phenotype reagent flow.

### Common options

- `--config`, `-c`
- `--quiet`, `-q`

## Phenotype Sheet (`phenotype-sheet`)

This command builds the gene-first `Stock Phenotype Sheet` from a Stage 1
workbook directory. It can also write an optional similarity sidecar
workbook.

### Command

```bash
python -m fl_ai_reagent_stocker phenotype-sheet ./gene_lists/Stocks --config ./my_split_config.json --similarity tiers
```

### Main output folder

- `./gene_lists/Stocks/Organized Stocks/`

For each input workbook, output is `<input_name>_aggregated.xlsx` containing the
`Stock Phenotype Sheet`. When `--similarity tiers|simple-buckets|keyword-buckets`
is set, the pipeline also writes
`<input_name>_aggregated_similarity_tiers.xlsx`.

### Reagent universe

The Stock Phenotype Sheet is driven by the complete Stage 1 **Gene Reagent
Index** (every FBal/FBtp/FBti reagent associated with each input gene), not by
stocks that survived Stage 1 ranking or JSON split limits. Stock metadata is
attached after phenotype-associated reagents are identified, using the full,
unfiltered `fbst_to_derived_stock_component.csv` so every FBst FlyBase
associates with a matched reagent appears in the sheet.

`phenotype-sheet` consumes Stage 1 workbook directories such as
`./gene_lists/Stocks` and expects a `Gene Reagent Index` sheet for full-scope
output. Workbooks without the index still run via a legacy fallback that is
best-effort and not true full scope; the command emits a stderr warning in
that case (even under `--quiet`) telling you to re-run `find-stocks` to
regenerate the index.

### `--similarity` choices

- `--similarity none` writes only the Stock Phenotype Sheet inside the
  aggregated workbook.
- `--similarity tiers` adds a sidecar workbook with:
  - `Contents` as the first tab
  - `Gene Set` as the second tab, copied from the input gene-list CSV data
    when available
  - `Stock Phenotype Sheet` as the third tab
  - fixed cosine-similarity threshold tabs after it (0.05 steps,
    `0.95-1.0`, `0.9-0.95`, ..., final `<0.05` bucket; empty buckets are
    skipped)
  - tier assignment based on the maximum cosine score across the configured
    `settings.phenotypeSimilarityTargets`
- `--similarity simple-buckets` writes the same sidecar with rule-based
  `collection / UAS / sleep-circ / balancer` tabs instead of cosine tiers:
  - `Contents` lists one row per combination with `Sheet name`, `# Stocks`,
    `# Alleles`, and `# Genes`
  - one sheet is written for every listed combination, including zero-count
    rows
  - genes, alleles, and reagents are assigned once and are never
    double-counted within or across combinations
  - column headers always show all three dimensions:
    - **UAS / non-UAS** — whether the stock resolves to the pipeline's UAS
      reagent bucket
    - **slp/ circ / Non slp/ circ** — whether any linked phenotype mentions
      "sleep" or "circadian"
    - **No bal / Has bal** — whether the stock carries at least one balancer
      chromosome
  - rows are ordered by collection first, then `UAS=true`,
    `sleep/circ=true`, and `No bal` before `Has bal`
- `--similarity keyword-buckets` writes the same sidecar with
  `Keyword Hits / No Keyword Hits` tabs instead of cosine tiers.

The default phenotype-similarity embedding model is `text-embedding-3-large`;
set `OPENAI_EMBEDDING_MODEL` to override it.

### Common options

- `--config`, `-c`
- `--quiet`, `-q`
- `--similarity {none,tiers,simple-buckets,keyword-buckets}` (default: `none`)

## Stage 3: Validate Ref++ Stocks (`validate-stocks`)

This stage reruns the same split logic, identifies `Ref++` output-sheet stocks, runs selective validation, and merges validation columns back into the organized workbook.

### Command

```bash
python -m fl_ai_reagent_stocker validate-stocks ./gene_lists/Stocks --config ./my_split_config.json --test-log
```

### Validation rules

Validation is intentionally selective:

1. the stock must survive filter and limit rules into a sheet whose combination includes `Ref++`
2. the reference must be keyword-relevant according to `settings.relevantSearchTerms`
3. full-text and pattern checks must make the pair eligible

Validation short-circuits on the first functional hit per stock.

`validate-stocks` no longer accepts `--soft-run`, `--OAI-embedding`,
`--simple-buckets`, or `--keyword-bucketing`. Use `phenotype-sheet` for the
phenotype reagent flow.

### Common options

- `--config`, `-c`
- `--quiet`, `-q`
- `--test-log`
- `--max-gpt-calls-per-stock`

## Full Pipeline

`run-full-pipeline` accepts a `--mode` choice that selects between two
end-to-end chains:

```bash
python -m fl_ai_reagent_stocker run-full-pipeline ./gene_lists --config ./my_split_config.json --mode normal --test-log
python -m fl_ai_reagent_stocker run-full-pipeline ./gene_lists --config ./my_split_config.json --mode phenotype --similarity tiers
```

- `--mode normal` runs Stage 1 + `split-stocks` + `validate-stocks`. It honors
  `--test-log` and `--max-gpt-calls-per-stock`.
- `--mode phenotype` runs Stage 1 + `phenotype-sheet` only and never invokes
  GPT validation. It honors `--similarity`.

Stage 1 args (`--gene-col`, `--input-gene-col`, `--batch-size/-b`,
`--skip-fbgnid-conversion`) apply to both modes. Validation flags are ignored
in `--mode phenotype`; `--similarity` is ignored in `--mode normal`.

Output paths stay the same as before:

- Stage 1 workbook: `./gene_lists/Stocks/aggregated_stock_refs.xlsx`
- Organized output directory: `./gene_lists/Stocks/Organized Stocks/`

## Config File Contract

The JSON config files stay in `data/config/` and keep the same effect on stock splitting:

- `settings.relevantSearchTerms` still defines keyword relevance and `Ref++`
- `settings.phenotypeSimilarityTargets` is required for phenotype-sheet cosine similarity targets
- `filters` still defines reusable filter predicates
- `combinations` still defines the sheet partitions
- `filterDescriptions` still defines user-facing explanations
- `maxStocksPerGene` and `maxStocksPerAllele` still define stock limits

Example configs:

- `data/config/stock_split_config_example.json`
- `data/config/stock_split_config_no_bloomington.json`

## Helper Scripts

Canonical helper entry points:

- `scripts/fetch_fbgn_ids.py`
- `scripts/build_fbst_derived_stock_components.py`
- `scripts/build_fbtp_to_fbti_mapping.py`
- `scripts/refresh_flybase_data.py`

### Refreshing FlyBase data

`refresh_flybase_data.py` discovers the latest TSV report families from
`https://s3ftp.flybase.org/releases/current/precomputed_files/`, downloads
them into `data/flybase/`, and prunes older versions in place.

```bash
python scripts/refresh_flybase_data.py            # TSVs only (default)
python scripts/refresh_flybase_data.py --with-xml # also refresh chado XML
python scripts/refresh_flybase_data.py --dry-run --with-xml
```

The optional chado XML families (`chado_FBst*.xml(.gz)`,
`chado_FBti*.xml(.gz)`) are downloaded from
`https://s3ftp.flybase.org/releases/current/chado-xml/` and are required by
`scripts/build_fbst_derived_stock_components.py` and
`scripts/build_fbtp_to_fbti_mapping.py`. Both build scripts discover the local
XML by prefix, so they keep working if FlyBase ever embeds a release suffix in
chado XML filenames.

## Reliability Notes

- Full-text retrieval uses multiple sources and retries transient failures.
- Stage 1 now supports FBtp -> FBti expansion through the generated mapping CSV.
- Stage 2 and Stage 3 share the same split/filter logic so config behavior stays stable across the refactor.
