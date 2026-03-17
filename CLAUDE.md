# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**fly-stocker-v2** is a modular Python pipeline that converts Drosophila (fruit fly) gene lists into prioritized, annotated stock sheets with linked scientific references and optional GPT-based functional validation.

Flow: gene lists → BDSC stocks → linked references → organized Excel output (with optional AI validation)

## Setup & Commands

```bash
pip install -r requirements.txt
```

Set environment variables (or `.env`):
```bash
OPENAI_API_KEY=...        # Required for GPT validation
NCBI_API_KEY=...          # Optional, improves PubMed rate limits
UNPAYWALL_TOKEN=...       # Optional, improves full-text coverage
OPENAI_MODEL=...          # Optional, overrides default (gpt-5-mini)
```

### Three CLI commands

**Pipeline 1** — Gene → Stocks → References:
```bash
python -m fly_stocker_v2.cli get-allele-refs ./gene_lists [options]
```
Output: `./gene_lists/Stocks/aggregated_bdsc_stock_refs.xlsx`

**Pipeline 2** — Organize & filter stocks by JSON rules:
```bash
python -m fly_stocker_v2.cli split-stocks ./gene_lists/Stocks --config ./my_config.json [options]
```
Output: `./gene_lists/Stocks/Organized Stocks/<name>_aggregated.xlsx`

**Full pipeline** — Run both stages end-to-end:
```bash
python -m fly_stocker_v2.cli run-full-pipeline ./gene_lists --config ./my_config.json [options]
```

### Common options
| Option | Description |
|--------|-------------|
| `--config`/`-c PATH` | JSON configuration file (keywords + filter rules) |
| `--run-validation` | Enable GPT functional validation (default: off) |
| `--soft-run` | Estimate GPT calls without executing |
| `--test-log` | Log all GPT queries to `data/logs/GPT-Queries/` |
| `--max-gpt-calls-per-stock N` | Limit GPT calls per stock (default: 5) |
| `--skip-fbgnid-conversion` | Skip FlyBase ID conversion (if input already has FBgn IDs) |
| `--quiet`/`-q` | Suppress verbose output (Pipeline 2 only) |

## Architecture

### Two-Stage Pipeline

**Stage 1 (`pipeline_references.py`)**: Data collection
1. Convert gene symbols → FlyBase IDs (FBgn)
2. Map genes → BDSC stocks via `data/BDSC/stockgenes.csv`
3. Link stocks/alleles → publications via FlyBase TSV data
4. Fetch PubMed metadata; score titles/abstracts against keywords
5. Optional GPT validation for functional evidence extraction

**Stage 2 (`pipeline_split.py`)**: Organization
1. Load Stage 1 Excel output
2. Compute derived columns: `Balancers`, `multiple_insertions`, `ALLELE_PAPER_RELEVANCE_SCORE`
3. Apply JSON filter rules → create one Excel sheet per filter combination
4. Selective GPT validation for Ref++ stocks only (keyword-relevant stocks)

### JSON-Driven Configuration

All filter logic lives in a JSON config file (see `data/JSON/stock_split_config_example.json`). No code changes are needed to define new analysis rules. The config specifies:
- **`keywords`**: Used for Ref± / Ref++ scoring of titles/abstracts
- **`filters`**: Named column-level predicates (field/type/value)
- **`filter_combinations`**: AND-ed filter sets → one output sheet each
- **`settings`**: Stock limits, sheet descriptions, GPT behavior

### Relevance Scoring Hierarchy
```
Ref++  (score 2): keyword hits in title/abstract — eligible for GPT validation
Ref+   (score 1): reference exists but no keyword match
Ref-   (score 0): no associated publications
```

### GPT Validation (selective)
Only triggered for Ref++ stocks. Validation short-circuits on first "Functionally validated" result. Output columns appended:
- `Functional Validity (by reference)?`
- `Functionally Valid Stock?`
- `PMID of [keywords] references that showed stocks functional validity`

### External Services (`external/`)
- **`pubmed.py`**: PubMed API client with CSV-based caching
- **`fulltext.py`**: Multi-source full-text retrieval (PMC → Unpaywall → CrossRef) + OpenAI validation

Both degrade gracefully when unavailable.

### Key Modules

| Module | Purpose |
|--------|---------|
| `cli.py` | Argument parsing, command routing |
| `config.py` | Settings dataclass, path resolution, validation status constants |
| `pipeline_references.py` | Stage 1 logic |
| `pipeline_split.py` | Stage 2 logic |
| `validation_runner.py` | Shared GPT validation logic used by both pipelines |
| `utils.py` | ID cleaning, FlyBase TSV loading, keyword extraction |

### Data Layout

Required static files (bundled or on network drive `/Volumes/umms-rallada/...`):
- `data/FlyBase/Alleles_And_Stocks/*.tsv(.gz)` — FlyBase allele/stock data
- `data/FlyBase/FlyBase_References/*.tsv(.gz)` — FlyBase publication mappings
- `data/BDSC/{stockgenes,bloomington,stockcomps_map_comments,balancers}.csv`

The code checks local `data/` first and falls back to the network drive (portable mode).

PubMed metadata is cached in CSV files at the network path to avoid redundant API calls across runs.
