# CLAUDE.md

This file provides guidance to coding agents working in this repository.

## Project Overview

`fl_ai_reagent_stocker` is a modular Python pipeline that converts Drosophila gene lists into prioritized stock sheets with linked references and optional AI-based functional validation.

High-level flow:

`gene lists -> find-stocks -> split-stocks -> validate-stocks`

Outputs stay in the same locations as before:

- Stage 1: `./gene_lists/Stocks/aggregated_stock_refs.xlsx`
- Stage 2/3: `./gene_lists/Stocks/Organized Stocks/<name>_aggregated.xlsx`

## Setup

```bash
pip install -r requirements.txt
```

Environment variables:

```bash
OPENAI_API_KEY=...   # Required for validation
NCBI_API_KEY=...     # Optional, improves PubMed rate limits
UNPAYWALL_TOKEN=...  # Optional, improves full-text coverage
OPENAI_MODEL=...     # Optional, overrides default
```

## Canonical Commands

### Stage 1: find-stocks

```bash
python -m fl_ai_reagent_stocker find-stocks ./gene_lists [options]
```

Purpose:

1. Convert gene symbols -> FBgn IDs
2. Map genes -> FlyBase stocks via `data/flybase/alleles_and_stocks/fbst_to_derived_stock_component.csv`
3. Expand construct-linked insertions via `data/flybase/transgenic_insertions/fbtp_to_fbti.csv`
4. Link stocks and components -> references
5. Score keyword relevance from title/abstract metadata

### Stage 2: split-stocks

```bash
python -m fl_ai_reagent_stocker split-stocks ./gene_lists/Stocks --config ./my_config.json [options]
```

Purpose:

1. Load the Stage 1 workbook
2. Compute derived columns such as `Balancers`, `multiple_insertions`, and `ALLELE_PAPER_RELEVANCE_SCORE`
3. Apply JSON filter rules and combinations
4. Write organized workbooks without GPT validation side effects

### Stage 3: validate-stocks

```bash
python -m fl_ai_reagent_stocker validate-stocks ./gene_lists/Stocks --config ./my_config.json [options]
```

Purpose:

1. Rebuild the same split outputs using the same JSON config
2. Identify `Ref++` output-sheet stocks
3. Run selective AI validation
4. Merge validation columns back into the organized workbook

### End-to-end wrapper

```bash
python -m fl_ai_reagent_stocker run-full-pipeline ./gene_lists --config ./my_config.json [options]
```

## JSON Config Contract

The JSON files in `data/config/` remain canonical and must keep the same behavior through refactors:

- `settings.relevantSearchTerms` defines keyword relevance and `Ref++`
- `settings.phenotypeSimilarityTargets` is required for phenotype-sheet cosine similarity targets
- `filters` defines reusable column predicates
- `combinations` defines sheet partitions
- `filterDescriptions` defines user-facing sheet descriptions
- `maxStocksPerGene` and `maxStocksPerAllele` define stock limits

Do not move or rename the existing config files unless explicitly requested.

## Architecture

All canonical code lives inside the `fl_ai_reagent_stocker/` package.
Root-level `.py` files (`config.py`, `utils.py`, `pipeline_references.py`,
`pipeline_split.py`, etc.) are thin backward-compat stubs that re-export
from the package. **Do not add new code to root-level files.**

### Package layout

```
fl_ai_reagent_stocker/
├── __init__.py                  # Public API surface
├── __main__.py                  # `python -m fl_ai_reagent_stocker`
├── cli.py                       # argparse CLI
├── config.py                    # Settings, paths, ValidationStatus
├── utils.py                     # ID cleaning, keyword helpers, TSV loading
├── validation_runner.py         # Shared GPT validation logic
├── integrations/
│   ├── pubmed.py                # PubMed/Entrez client and cache
│   └── fulltext.py              # Full-text retrieval + OpenAI validator
└── pipelines/
    ├── stock_finding.py         # Stage 1: genes → stocks → references
    └── stock_splitting.py       # Stage 2/3: filters, limits, Excel output
```

### Key classes

- **`StockFindingPipeline`** (`pipelines/stock_finding.py`): Stage 1.
  Loads FlyBase data, maps `FBgn → FBal → FBtp/FBti → FBst`, fetches
  references, scores keywords, writes the Stage 1 workbook.

- **`StockSplittingPipeline`** (`pipelines/stock_splitting.py`): Stage 2
  *and* Stage 3. Loads Stage 1 output, computes derived columns, applies
  JSON filters and stock limits, writes organized workbooks. Pass
  `run_validation=True` to also run GPT validation (Stage 3).

- **`Settings`** (`config.py`): Dataclass holding API keys, paths, and
  feature flags. Loads `.env` at init time.

## Helper Scripts

Canonical helper scripts (standalone, no package imports):

- `scripts/fetch_fbgn_ids.py`
- `scripts/build_fbst_derived_stock_components.py`
- `scripts/build_fbtp_to_fbti_mapping.py`
- `scripts/refresh_flybase_data.py`

## Data Layout

Expected local data layout:

- `data/flybase/alleles_and_stocks/*.tsv(.gz)`
- `data/flybase/references/*.tsv(.gz)`
- `data/flybase/transgenic_constructs/*.tsv(.gz)`
- `data/flybase/transgenic_insertions/chado_FBti.xml.gz`
- `data/flybase/transgenic_insertions/fbtp_to_fbti.csv`
- `data/config/*.json`

The code prefers local `data/` and falls back to `../Data` for legacy layouts.
