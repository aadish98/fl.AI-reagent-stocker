# fl_ai_reagent_stocker

Modular Drosophila stock-processing pipeline for turning gene lists into stock workbooks, config-driven organized sheets, and optional AI-assisted validation.

## Documentation

- Practical guide (Markdown): `docs/pipeline_usage.md`
- Practical guide (HTML): `docs/pipeline_usage.html`
- Practical guide (PDF): `docs/pipeline_usage.pdf`
- Data flowchart: `docs/images/pipeline-data-flowchart.png`
- High-level flowchart: `docs/images/pipeline-flowchart-high-level-standalone.png`
- Rules/evidence flowchart: `docs/images/pipeline-flowchart-rules-and-tools-standalone.png`
- End-to-end flow diagram: `docs/images/pipeline-end-to-end-flowchart.png`
- Example split config: `data/config/stock_split_config_example.json`

## What The Pipeline Produces

After a full run, you get:

- a stock-level workbook linked to references
- organized output sheets based on your JSON rules
- a references sheet narrowed to papers cited by selected stocks
- optional validation columns for `Ref++` output-sheet stocks

## Pipeline Overview

The canonical package is `fl_ai_reagent_stocker` and it exposes three explicit stages plus an end-to-end wrapper:

1. `find-stocks`
   - Reads gene-list CSV files.
   - Converts gene symbols to FBgn IDs unless skipped.
   - Maps genes to FlyBase stocks through allele, construct, and insertion components.
   - Links stocks and components to publications and PMIDs.
   - Writes `./gene_lists/Stocks/aggregated_stock_refs.xlsx`.

2. `split-stocks`
   - Reads the Stage 1 workbook.
   - Applies `data/config/*.json` filters and combinations.
   - Writes organized output sheets to `./gene_lists/Stocks/Organized Stocks/`.
   - Does not run GPT validation.

3. `validate-stocks`
   - Rebuilds the same split outputs using the same JSON config.
   - Selectively validates only `Ref++` output-sheet stocks.
   - Merges validation columns back into the organized workbook.

4. `run-full-pipeline`
   - Runs Stages 1, 2, and 3 in sequence.

The JSON configs in `data/config/` stay in place and keep the same effect on stock splitting as before.

## Installation

From the repository root:

```bash
pip install -r requirements.txt
```

Optional environment variables:

- `OPENAI_API_KEY` for validation
- `NCBI_API_KEY` for improved PubMed rate limits
- `UNPAYWALL_TOKEN` for improved full-text retrieval coverage
- `OPENAI_MODEL` to override the default model

## Quick Start

### One command

```bash
python -m fl_ai_reagent_stocker run-full-pipeline ./gene_lists --config ./my_split_config.json --test-log
```

### Three-step workflow

```bash
python -m fl_ai_reagent_stocker find-stocks ./gene_lists
python -m fl_ai_reagent_stocker split-stocks ./gene_lists/Stocks --config ./my_split_config.json
python -m fl_ai_reagent_stocker validate-stocks ./gene_lists/Stocks --config ./my_split_config.json --test-log
```

Default config path:

- `data/config/stock_split_config_example.json`

## Command Summary

### Stage 1: `find-stocks`

```bash
python -m fl_ai_reagent_stocker find-stocks ./gene_lists
```

Main output:

- `./gene_lists/Stocks/aggregated_stock_refs.xlsx`

Common options:

- `--config`, `-c`
- `--skip-fbgnid-conversion`
- `--gene-col`
- `--input-gene-col`
- `--batch-size`

### Stage 2: `split-stocks`

```bash
python -m fl_ai_reagent_stocker split-stocks ./gene_lists/Stocks --config ./my_split_config.json
```

Main output folder:

- `./gene_lists/Stocks/Organized Stocks/`

Common options:

- `--config`, `-c`
- `--quiet`, `-q`
- `--soft-run`
- `--OAI-embedding`
- `--simple-buckets`

### Stage 3: `validate-stocks`

```bash
python -m fl_ai_reagent_stocker validate-stocks ./gene_lists/Stocks --config ./my_split_config.json --test-log
```

Validation is selective by design:

- only `Ref++` output-sheet stocks are considered
- references must be keyword-relevant
- validation short-circuits after the first functional hit per stock

Common options:

- `--config`, `-c`
- `--quiet`, `-q`
- `--soft-run`
- `--OAI-embedding`
- `--simple-buckets`
- `--test-log`
- `--max-gpt-calls-per-stock`

## Config File Contract

The stock-splitting JSON files remain under `data/config/` and retain their existing semantics:

- `settings.relevantSearchTerms` still defines `Ref++`
- `settings.phenotypeSimilarityTargets` is required for phenotype-sheet cosine similarity targets
- `filters` still define reusable predicates
- `combinations` still define the sheet partitions
- `filterDescriptions` still control user-facing sheet descriptions
- `maxStocksPerGene` and `maxStocksPerAllele` still enforce the same stock limits

## Helper Scripts

Canonical helper entry points:

- `scripts/fetch_fbgn_ids.py`
- `scripts/build_fbst_derived_stock_components.py`
- `scripts/build_fbtp_to_fbti_mapping.py`
- `scripts/refresh_flybase_data.py`

## Notes

- Stage 1 now supports FBtp -> FBti expansion via `data/flybase/transgenic_insertions/fbtp_to_fbti.csv`.
- Soft-run mode still stops before GPT execution.
- Output paths remain unchanged even though the package and command names changed.
