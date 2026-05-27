# AGENTS.md

This file provides guidance to coding agents working in this repository.

## Project Overview

`fl_ai_reagent_stocker` is a modular Python pipeline that converts Drosophila
gene lists into prioritized stock sheets with linked references, optional
phenotype-similarity scoring, and optional AI-based functional validation.

High-level flows:

```
gene lists -> find-stocks -> split-stocks -> validate-stocks         (normal)
gene lists -> find-stocks -> phenotype-sheet                         (phenotype)
```

`run-full-pipeline --mode normal` runs the first chain end-to-end.
`run-full-pipeline --mode phenotype` runs the second chain and never invokes
GPT validation. `phenotype-sheet` is its own command and is no longer a flag
on `split-stocks` or `validate-stocks`.

Output locations:

- Stage 1: `./gene_lists/Stocks/aggregated_stock_refs.xlsx`
- Stage 2/3: `./gene_lists/Stocks/Organized Stocks/<name>_aggregated.xlsx`
- phenotype-sheet sidecar (when `--similarity tiers|simple-buckets|keyword-buckets`):
  `./gene_lists/Stocks/Organized Stocks/<name>_aggregated_similarity_tiers.xlsx`

## Related Documentation

These user-facing docs must move in lockstep with CLI / pipeline changes.
When you add, rename, or change the behavior of a command, flag, output
artifact, JSON config key, helper script, or environment variable, update the
relevant doc(s) in the same change:

- `README.md` — top-level overview, install, quick start, command summary,
  config contract, helper-script entry points.
- `docs/pipeline_usage.md` — practical end-to-end guide; mirrors the CLI
  surface and `--similarity` semantics in user terms.
- `docs/pipeline_flowcharts.md` — Mermaid flowcharts for the data flow,
  high-level command flow, and low-level evidence flow.
- `docs/citations.md` — how downstream users cite the data sources, APIs, and
  AI models the pipeline depends on. Update when adding/removing an external
  data source or model dependency, or when the default model changes.
- `docs/AGENTS.md` — this file; agent-facing reference for commands, package
  layout, key classes, helper scripts, tests, and data layout.

## Setup

```bash
pip install -r requirements.txt
```

Environment variables (loaded from repo-root `.env`):

```bash
OPENAI_API_KEY=...           # Required for validation and OpenAI embeddings
NCBI_API_KEY=...             # Optional, improves PubMed rate limits
UNPAYWALL_TOKEN=...          # Optional, improves full-text coverage
OPENAI_MODEL=...             # Optional, overrides default (gpt-5-mini)
OPENAI_EMBEDDING_MODEL=...   # Optional, overrides default (text-embedding-3-large)
```

## Canonical Commands

### Stage 1: find-stocks

```bash
python -m fl_ai_reagent_stocker find-stocks ./gene_lists [options]
```

Purpose:

1. Convert gene symbols -> FBgn IDs (unless `--skip-fbgnid-conversion`)
2. Map genes -> FlyBase stocks via
   `data/flybase/alleles_and_stocks/fbst_to_derived_stock_component.csv`
3. Expand construct-linked insertions via
   `data/flybase/transgenic_insertions/fbtp_to_fbti.csv`
4. Link stocks and components -> references
5. Score keyword relevance from title/abstract metadata

Common options: `--config/-c`, `--gene-col`, `--input-gene-col`,
`--batch-size/-b`, `--skip-fbgnid-conversion`.

### Stage 2: split-stocks (organization only)

```bash
python -m fl_ai_reagent_stocker split-stocks ./gene_lists/Stocks --config ./my_config.json [options]
```

Purpose:

1. Load the Stage 1 workbook
2. Compute derived columns such as `Balancers`, `multiple_insertions`, and
   `ALLELE_PAPER_RELEVANCE_SCORE`
3. Apply JSON filter rules and combinations
4. Write organized workbooks without GPT validation or phenotype-sheet side
   effects

`split-stocks` is intentionally narrow: it no longer accepts `--soft-run`,
`--OAI-embedding`, `--simple-buckets`, or `--keyword-bucketing`. Use
`phenotype-sheet` for the phenotype reagent flow.

Common options: `--config/-c`, `--quiet/-q`.

### phenotype-sheet (gene-first phenotype reagent flow)

```bash
python -m fl_ai_reagent_stocker phenotype-sheet ./gene_lists/Stocks --config ./my_config.json --similarity tiers
```

Purpose:

1. Load the Stage 1 workbook
2. Drive the `Stock Phenotype Sheet` from the **complete** Stage 1
   `Gene Reagent Index` (every FBal / FBtp / FBti reagent associated with each
   input gene) filtered against FlyBase `genotype_phenotype_data`
3. Attach stock metadata afterward from the full
   `fbst_to_derived_stock_component.csv` so the sheet is not restricted to
   stocks that survived Stage 1 ranking or JSON split limits
4. Optionally write a similarity sidecar workbook based on `--similarity`

`phenotype-sheet` consumes Stage 1 workbook directories such as
`./gene_lists/Stocks` and expects a `Gene Reagent Index` sheet for full-scope
output. Workbooks without the index still run via a legacy fallback that is
best-effort and not true full scope; the command emits a stderr warning in
that case (even under `--quiet`) telling you to re-run `find-stocks` to
regenerate the index.

`--similarity` choices:

- `none` writes only the Stock Phenotype Sheet inside the aggregated workbook.
- `tiers` adds the cosine-threshold tier sidecar workbook computed against
  `settings.phenotypeSimilarityTargets`.
- `simple-buckets` adds the rule-based
  `collection / UAS / sleep-circ / balancer` sidecar workbook.
- `keyword-buckets` adds the `Keyword Hits / No Keyword Hits` sidecar
  workbook.

Common options: `--config/-c`, `--quiet/-q`, `--similarity {none,tiers,simple-buckets,keyword-buckets}`.

### Stage 3: validate-stocks (GPT validation only)

```bash
python -m fl_ai_reagent_stocker validate-stocks ./gene_lists/Stocks --config ./my_config.json [options]
```

Purpose:

1. Rebuild the same split outputs using the same JSON config
2. Identify `Ref++` output-sheet stocks
3. Run selective AI validation (short-circuits on first functional hit per stock)
4. Merge validation columns back into the organized workbook

`validate-stocks` no longer accepts `--soft-run`, `--OAI-embedding`,
`--simple-buckets`, or `--keyword-bucketing`. Use `phenotype-sheet` for the
phenotype reagent flow.

Validation-specific options: `--test-log` (write GPT call logs under
`data/logs/gpt_queries/`), `--max-gpt-calls-per-stock` (default `5`).

### End-to-end wrapper

```bash
python -m fl_ai_reagent_stocker run-full-pipeline ./gene_lists --config ./my_config.json --mode normal [options]
python -m fl_ai_reagent_stocker run-full-pipeline ./gene_lists --config ./my_config.json --mode phenotype --similarity tiers
```

`--mode normal` runs Stage 1 + `split-stocks` + `validate-stocks` and accepts
the Stage 1 args plus `--test-log` and `--max-gpt-calls-per-stock`.
`--mode phenotype` runs Stage 1 + `phenotype-sheet` only and accepts the
Stage 1 args plus `--similarity`. Validation flags are ignored in
`--mode phenotype`; `--similarity` is ignored in `--mode normal`.

Stage 1 always writes a `Gene Reagent Index` sheet alongside `Stocks` /
`References` in `aggregated_stock_refs.xlsx`. Older workbooks without the
sheet still run via a legacy fallback that derives the reagent universe from
the `Stocks` sheet (best-effort, not true full scope).

## JSON Config Contract

The JSON files in `data/config/` are canonical and must keep the same behavior
through refactors. The default is `stock_split_config_example.json`.

- `settings.relevantSearchTerms` defines keyword relevance and `Ref++`
- `settings.phenotypeSimilarityTargets` is required for phenotype-sheet
  cosine similarity targets (each entry needs `keyword` + `embedding_text`,
  validated by `normalize_phenotype_similarity_targets` in `config.py`)
- `settings.maxStocksPerGene` and `settings.maxStocksPerAllele` define stock
  limits
- `filters` defines reusable column predicates
- `combinations` defines sheet partitions
- `filterDescriptions` defines user-facing sheet descriptions

Do not move or rename the existing config files unless explicitly requested.

## Architecture

All canonical code lives inside the `fl_ai_reagent_stocker/` package. Root-level
`.py` files (`config.py`, `utils.py`, `pipeline_references.py`,
`pipeline_split.py`, `validation_runner.py`) are thin backward-compat stubs
that re-export from the package. **Do not add new code to root-level files.**

### Package layout

```
fl_ai_reagent_stocker/
├── __init__.py                  # Public API surface (Settings, pipelines)
├── __main__.py                  # `python -m fl_ai_reagent_stocker`
├── cli.py                       # argparse CLI
├── config.py                    # Settings, paths, ValidationStatus,
│                                # normalize_phenotype_similarity_targets
├── utils.py                     # ID cleaning, keyword helpers, TSV loading,
│                                # reagent-bucket / RNAi-type helpers
├── validation_runner.py         # Shared GPT validation logic
├── integrations/
│   ├── __init__.py
│   ├── pubmed.py                # PubMed/Entrez client and PubMedCache
│   ├── fulltext.py              # Full-text retrieval + FunctionalValidator
│   └── phenotype_similarity.py  # OpenAI embedding cosine scoring + t-SNE plots
└── pipelines/
    ├── __init__.py
    ├── stock_finding.py         # Stage 1: genes → stocks → references → Gene Reagent Index
    └── stock_splitting.py       # split-stocks, phenotype-sheet, validate-stocks:
                                 # filters, limits, Excel output, gene-first
                                 # phenotype reagent flow + similarity sidecars
```

### Key classes

- **`StockFindingPipeline`** (`pipelines/stock_finding.py`): Stage 1.
  Loads FlyBase data, maps `FBgn → FBal → FBtp/FBti → FBst`, fetches
  references, scores keywords, writes the Stage 1 workbook.

- **`StockSplittingPipeline`** (`pipelines/stock_splitting.py`): Stage 2,
  Stage 3, *and* the phenotype-sheet flow. Loads Stage 1 output, computes
  derived columns, applies JSON filters and stock limits, writes organized
  workbooks. Pass `run_validation=True` to also run GPT validation
  (`validate-stocks`). Set `Settings.soft_run=True` (selected by the
  `phenotype-sheet` CLI command via `--similarity`) to switch to the
  phenotype reagent flow and emit the optional similarity sidecars.

- **`Settings`** (`config.py`): Dataclass holding API keys, model names,
  paths, and feature flags (`soft_run`, `enable_oai_embedding`,
  `simple_buckets`, `keyword_bucketing`, `enable_gpt_logging`,
  `max_gpt_calls_per_stock`, `phenotype_similarity_targets`,
  `skip_fbgnid_conversion`). Loads `.env` at init time. Default model is
  `gpt-5-mini`; default embedding model is `text-embedding-3-large`.

- **`ValidationStatus`** (`config.py`): Validation status constants,
  category mapping, and priority ordering used by Stage 3 outputs.

- **`EmbeddingSimilarityScorer`** (`integrations/phenotype_similarity.py`):
  Computes phenotype embeddings against configured targets and powers the
  similarity-tier sidecar workbook.

## Helper Scripts

Canonical helper scripts live in `scripts/` (standalone, no package imports).
Run them from the repo root.

Data builders / refreshers:

- `refresh_flybase_data.py` - pull the latest FlyBase TSV/TSV.GZ reports.
- `fetch_fbgn_ids.py` - convert gene symbols to FBgn IDs.
- `build_fbst_derived_stock_components.py` - build
  `data/flybase/alleles_and_stocks/fbst_to_derived_stock_component.csv`.
- `build_fbtp_to_fbti_mapping.py` - build
  `data/flybase/transgenic_insertions/fbtp_to_fbti.csv`.
- `build_fbsf_to_fbgn_mapping.py` - build
  `data/flybase/sequence_features/fbsf_to_fbgn.csv` from the FlyBase SQL dump.

Audits and aggregations:

- `audit_stock_phenotype_pipeline.py` - write Chado/phenotype audit JSON +
  Markdown to `audit_outputs/`.
- `aggregate_similarity_workbook_counts.py` - summarize
  `_aggregated_similarity_tiers.xlsx` Contents tables across a directory tree.

Visualizations / decks:

- `generate_flowchart.py`, `generate_pipeline_data_flowchart.py`
- `generate_openai_embedding_tsne_demo.py`, `visualize_embeddings_tsne.py`
- `create_stock_phenotype_mapping_deck.js`

## Tests

Run from the repo root:

```bash
python -m pytest tests
```

Notable suites:

- `tests/test_phenotype_reagent_capture.py`
- `tests/test_phenotype_similarity_pipeline.py`
- `tests/test_aggregate_similarity_workbooks.py`
- `tests/test_stock_phenotype_audit.py`

Tests rely on the package being importable from the repo root (they prepend
`REPO_ROOT` to `sys.path`).

## Data Layout

Expected local data layout under `data/`:

- `flybase/alleles_and_stocks/`
  - `fbst_to_derived_stock_component.csv` (built helper)
  - `fbal_to_fbgn_fb_*.tsv(.gz)`
  - `stocks_FB*.tsv.gz`
  - `dmel_classical_and_insertion_allele_descriptions_fb_*.tsv.gz`
  - `genotype_phenotype_data_fb_*.tsv(.gz)`
  - `chado_FBst.xml.gz` (gitignored)
- `flybase/references/`
  - `entity_publication_fb_*.tsv.gz` (gitignored)
  - `fbrf_pmid_pmcid_doi_fb_*.tsv.gz`
- `flybase/transgenic_constructs/`
  - `transgenic_construct_descriptions_fb_*.tsv(.gz)`
- `flybase/transgenic_insertions/`
  - `fbtp_to_fbti.csv` (built helper)
- `flybase/genes/`
  - `fb_synonym_fb_*.tsv.gz`
- `flybase/sequence_features/`
  - `fbsf_to_fbgn.csv` (built helper, gitignored dir)
- `config/*.json`
- `cache/` (gitignored) - PubMed and embedding caches:
  - `pmid_to_title_abstract.csv`
  - `pmid_to_fulltext_method.csv`
  - `phenotype_text_embeddings.csv`
  - `phenotype_target_embeddings.csv`
- `logs/gpt_queries/` (gitignored) - GPT call logs when `--test-log` is set

The code prefers the local `data/` directory and falls back to a sibling
`../Data` directory for legacy layouts (`is_portable_mode()` returns whether
the local layout is in use).
