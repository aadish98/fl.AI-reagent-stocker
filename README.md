# fly-stocker-v2

Modular Drosophila stock-processing pipeline for going from gene lists to prioritized stock sheets with linked references and optional GPT-based functional validation.

This README consolidates the existing project documentation and points to the full detailed guides.

## Documentation Map

- Full practical guide (Markdown): `PIPELINE_USAGE.md`
- Full practical guide (HTML): `PIPELINE_USAGE.html`
- Full practical guide (PDF): `PIPELINE_USAGE.pdf`
- High-level flowchart: `pipeline-flowchart-high-level-standalone.png`
- Rules/evidence flowchart: `pipeline-flowchart-rules-and-tools-standalone.png`
- End-to-end flow diagram: `pipeline-end-to-end-flowchart.png`
- Example split config: `data/JSON/stock_split_config_example.json`

## What the Pipeline Produces

After a full run, you get:

- a stock-level table linked to references
- organized output sheets based on your filtering rules
- a references table narrowed to papers relevant to selected stocks
- optional functional-validation signals for `Ref++` stocks

## Pipeline Overview

The project has two CLI stages (plus an end-to-end wrapper command):

1. `get-allele-refs`
   - Reads gene-list CSV files.
   - Converts gene symbols to FBgn IDs (unless skipped).
   - Maps genes to BDSC stocks.
   - Links stocks/alleles to references and PMIDs.
   - Produces an aggregated workbook.

2. `split-stocks`
   - Reads Step 1 workbook(s).
   - Applies JSON-driven filters and combinations.
   - Builds organized output sheets.
   - Runs selective GPT validation for eligible `Ref++` stock-reference pairs.

3. `run-full-pipeline`
   - Runs Step 1 then Step 2 in one command.

## Installation

From the repository root:

```bash
pip install -r requirements.txt
```

## Before You Run

- Run commands from this project directory (the folder containing `cli.py` and `data/JSON/`).
- Put your input CSV files in one folder (for example, `./gene_lists`).
- Optional environment variables:
  - `OPENAI_API_KEY` for GPT validation
  - `NCBI_API_KEY` to improve PubMed API rate limits
  - `UNPAYWALL_TOKEN` to improve full-text retrieval coverage

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

## Commands

### Step 1: Build Stocks + References

```bash
python -m fly_stocker_v2.cli get-allele-refs ./gene_lists
```

Main output:

- `./gene_lists/Stocks/aggregated_bdsc_stock_refs.xlsx`

Workbook sheets:

- `Stocks`
- `References`

Common options:

- `--config`, `-c`
- `--skip-fbgnid-conversion`
- `--run-validation`
- `--soft-run`
- `--test-log`
- `--max-gpt-calls-per-stock` (default: `5`)

### Step 2: Organize and Prioritize

```bash
python -m fly_stocker_v2.cli split-stocks ./gene_lists/Stocks --config ./my_split_config.json
```

Main output folder:

- `./gene_lists/Stocks/Organized Stocks/`

For each input workbook, output is `<input_name>_aggregated.xlsx`.

Common options:

- `--config`, `-c`
- `--quiet`, `-q`
- `--soft-run`
- `--test-log`
- `--max-gpt-calls-per-stock` (default: `5`)

### Full Pipeline

```bash
python -m fly_stocker_v2.cli run-full-pipeline ./gene_lists --config ./my_split_config.json --test-log
```

Use full-pipeline mode when you want a single end-to-end run. Use separate commands when you want to inspect intermediate outputs.

## Validation Behavior (Important)

Validation is selective by design. In Step 2, a stock-reference pair is considered for GPT validation only when all are true:

1. The stock survives filtering/limits into a sheet whose combination includes `Ref++`.
2. The reference is keyword-relevant (from `settings.relevantSearchTerms` in config).
3. Full-text and pattern checks mark the pair as eligible.

If a stock only appears in non-`Ref++` sheets, it is not sent for GPT validation.

## Config File (Plain Terms)

The JSON split config controls:

- relevance keywords (`settings.relevantSearchTerms`)
- maximum stocks per gene (`settings.maxStocksPerGene`)
- maximum stocks per allele (`settings.maxStocksPerAllele`)
- reusable named filters (`filters`)
- output sheet definitions (`combinations`)
- user-facing filter descriptions (`filterDescriptions`)

See:

- `data/JSON/stock_split_config_example.json`

## Reliability and Reproducibility Notes

- Full-text retrieval uses multiple sources and retries transient failures.
- Caches are used to reduce repeated metadata/full-text method lookups.
- Duplicate rows with the same PMID are merged with preference for richer IDs (PMCID/DOI).
- Validation can stop early per stock when functional evidence is already found.

## Quick Troubleshooting

- Validation not running: confirm output sheet is `Ref++` and `OPENAI_API_KEY` is set.
- Low PMID relevance: check `settings.relevantSearchTerms` in config.
- Sparse full-text coverage: set `UNPAYWALL_TOKEN` and rerun.
- Concerned about API usage: start with `--soft-run`.

## Detailed Guides and Visuals

For the complete narrative, full example run-through, and visual flows, use:

- `PIPELINE_USAGE.md`
- `PIPELINE_USAGE.html`
- `PIPELINE_USAGE.pdf`
- `pipeline-flowchart-high-level-standalone.png`
- `pipeline-flowchart-rules-and-tools-standalone.png`
- `pipeline-end-to-end-flowchart.png`
