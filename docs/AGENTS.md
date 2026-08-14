# AGENTS.md

This file provides guidance to coding agents working in this repository.

## Project Overview

`fl_ai_reagent_stocker` is a modular Python pipeline that converts Drosophila
gene lists into prioritized stock sheets with linked references,
curated-phenotype and publication-evidence tiers, and optional AI-based
functional validation and phenotype-embedding analysis.

Primary flow:

```
gene lists -> run  (find-stocks -> organize -> validate [+ embeddings])
```

`run` is the single public command. It builds the Stage 1 workbook, organizes
stocks with the JSON config, and validates `Ref++` stocks or runs embeddings
when enabled by the JSON config. The individual stages
(`find-stocks`, `split-stocks`, `phenotype-sheet`, `validate-stocks`) remain as
advanced entry points for debugging and tests.

Output locations (primary `run`, recursive over input gene-set CSVs):

- Per gene set:
  `./gene_lists/Per Gene Set Runs/<relative path>/<gene set>/Stocks/<gene set>_aggregated.xlsx`
- Combined summary at the input root: `combination_counts_summary.xlsx` / `.csv`

The primary `run` produces a clean tree: the final aggregated workbook lives
directly in `Stocks` (no `Organized Stocks` wrapper), and it does not emit
`.txt` sidecars, `*_similarity_tiers.xlsx` workbooks, `*_similarity/` plot
folders, or the intermediate unsplit `<gene set>.xlsx` unless
`settings.output.preserveUnsplitWorkbook` is `true`. `settings.embeddings.enabled`
adds cosine columns to the aggregated workbook only.

The advanced per-stage commands (`split-stocks`, `phenotype-sheet`) keep the
legacy `Stocks/Organized Stocks/<name>_aggregated.xlsx` layout, and
`phenotype-sheet` still writes the
`<name>_aggregated_similarity_tiers.xlsx` sidecar plus a
`<name>_aggregated_similarity/` plot directory when embeddings are enabled in
config.

## Related Documentation

These user-facing docs must move in lockstep with CLI / pipeline changes.
When you add, rename, or change the behavior of a command, flag, output
artifact, JSON config key, helper script, or environment variable, update the
relevant doc(s) in the same change:

- `README.md` — top-level overview, install, quick start, command summary,
  config contract, helper-script entry points.
- `docs/pipeline_usage.md` — practical end-to-end guide; mirrors the CLI
  surface and config semantics in user terms.
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
OPENAI_API_KEY=...           # Required for validation and phenotype embeddings when enabled in config
NCBI_API_KEY=...             # Optional, improves PubMed rate limits
UNPAYWALL_TOKEN=...          # Optional contact email; defaults to aadishms@umich.edu
OPENAI_MODEL=...             # Optional, overrides default (gpt-5-mini)
OPENAI_EMBEDDING_MODEL=...   # Optional, overrides default (text-embedding-3-large)
```

## Commands

### run (primary, config-driven)

```bash
python -m fl_ai_reagent_stocker run ./gene_lists --config ./my_config.json [--quiet]
```

`run` executes the full pipeline:

1. `find-stocks`: convert gene symbols -> FBgn IDs unless
   `settings.input.skipFbgnidConversion` is `true`, map genes -> FlyBase stocks, expand
   construct-linked insertions, link references, score keyword relevance, and
   write the Stage 1 workbook with `Stocks`, `References`, and the complete
   `Gene Reagent Index`.
2. organize: compute derived columns (`Balancers`, `multiple_insertions`,
   `ALLELE_PAPER_RELEVANCE_SCORE`, and `PHENOTYPE_RELEVANCE_SCORE` when a filter
   needs it), then apply the JSON filters and combinations.
3. validate: selectively GPT-validate `Ref++` output stocks and merge the
   validation columns back into the organized workbook.

`settings.embeddings.enabled` additionally runs the phenotype-sheet flow with
OpenAI embeddings to add cosine columns to the aggregated workbook. For the primary
`run`, the similarity tier sidecar and t-SNE/plot folders are suppressed; the
input is scanned recursively for gene-set CSVs; the final aggregated workbook
is written directly into each gene set's `Stocks` folder; and a combined
combination-counts summary is generated at the input root when the run finishes.
During a multi-gene-set `run`, Stage 1 reuses one `StockFindingPipeline`
instance, the gene-set-independent FlyBase auxiliary loaders in
`stock_splitting.py` are memoized per process, and `PubMedCache.load()` reuses a
parsed CSV cache by `(path, mtime)` while copying data into each instance. The
Ref++ validation pass is skipped when no `OPENAI_API_KEY` is configured or when
the required JSON config block `settings.validation.enabled` is `false`.

Options: `--config/-c`, `--quiet/-q`. Input columns, FBgn conversion, PubMed
batch size, embeddings, output preservation, and validation policy are required
JSON config fields.

### Advanced per-stage entry points

These subcommands run individual stages with the same config and outputs as
`run`:

```bash
python -m fl_ai_reagent_stocker find-stocks ./gene_lists --config ./my_config.json
python -m fl_ai_reagent_stocker split-stocks ./gene_lists/Stocks --config ./my_config.json [--quiet]
python -m fl_ai_reagent_stocker phenotype-sheet ./gene_lists/Stocks --config ./my_config.json [--quiet]
python -m fl_ai_reagent_stocker validate-stocks ./gene_lists/Stocks --config ./my_config.json [--quiet]
```

- `split-stocks` writes organized workbooks only (no validation or
  phenotype-sheet side effects).
- `phenotype-sheet` drives the gene-first `All Phenotypic Stocks Sheet` from the
  Stage 1 `Gene Reagent Index` filtered against FlyBase
  `genotype_phenotype_data`, attaching stock metadata from the full
  `fbst_to_derived_stock_component.csv`. `settings.embeddings.enabled` adds the
  cosine analysis and similarity-tier sidecar; without it, only the All
  Phenotypic Stocks Sheet is written. It expects a `Gene Reagent Index` sheet for full-scope output and
  warns to stderr (even under `--quiet`) when one is missing.
- `validate-stocks` rebuilds the split outputs and GPT-validates `Ref++`
  stocks, short-circuiting on the first functional hit per stock.

Stage 1 always writes a `Gene Reagent Index` sheet alongside `Stocks` /
`References` in `aggregated_stock_refs.xlsx`.

## JSON Config Contract

The JSON files in `data/config/` are canonical. The default is
`stock_split_config_example.json`; `stock_split_config_phenotype_example.json`
is the phenotype-aware variant. `stock_split_config_priority_example.json` and
`stock_split_config_priority_all_phenotypes.json` are the phenotype-bucket
priority configs (see `docs/stock_split_config_priority_example.md`).

- `settings.relevantSearchTerms` defines keyword relevance, backing
  `ALLELE_PAPER_RELEVANCE_SCORE` and the `Ref++` / `Ref+` / `Ref-` tiers
- `settings.phenotypeSimilarityTargets` defines the curated-phenotype targets
  (each entry needs `keyword` + `embedding_text`, validated by
  `normalize_phenotype_similarity_targets` in `config.py`). Used both for the
  `PHENOTYPE_RELEVANCE_SCORE` tiers and for the OpenAI embedding analysis.
- `settings.maxStocksPerGene` and `settings.maxStocksPerAllele` define stock
  limits
- `settings.input` is required and defines gene columns plus FBgn conversion
- `settings.pubmed` is required and defines PubMed/full-text batch size
- `settings.embeddings` is required and controls embedding analysis
- `settings.output` is required and controls output preservation
- `settings.validation` is required in every JSON config; shipped configs set
  `enabled: false` explicitly
- `settings.validation.maxGptCallsPerStock`, `minFullTextChars`,
  `gptCallDelaySeconds`, `shortCircuitOnFunctionalValidation`, and
  `enableGptLogging` are required validation behavior controls
- `settings.contentsSeparatorEvery` defines how often the `Contents` sheet
  breakdown table draws a faint separator row (default `3`)
- `filters` defines reusable column predicates
- `combinations` defines ordered sheet partitions; each `(stock_id, collection)`
  reagent is placed in only the first combination it matches
- optional `combination_names` defines Excel tab names aligned one-to-one with
  `combinations`; configured labels are stored without numbers and populated
  sheets receive sequential `1_`, `2_`, ... prefixes. Absent names preserve
  legacy non-empty `Sheet1..N` numbering
- `filterDescriptions` defines user-facing sheet descriptions

Computed score columns available to filters:

- `ALLELE_PAPER_RELEVANCE_SCORE`: `2` (`Ref++`), `1` (`Ref+`), `0` (`Ref-`).
- `PHENOTYPE_RELEVANCE_SCORE`: `2` (`Phenotype++`, a component has a curated
  phenotype matching a `phenotypeSimilarityTargets` term), `1` (`Phenotype+`,
  any curated phenotype), `0` (none). Computed by
  `compute_phenotype_relevance_score` in `pipelines/stock_splitting.py` only
  when a filter references the column.

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
    └── stock_splitting.py       # organize, phenotype-sheet, validate:
                                 # filters, limits, Excel output, gene-first
                                 # phenotype reagent flow + embedding analysis
```

### Key classes

- **`StockFindingPipeline`** (`pipelines/stock_finding.py`): Stage 1.
  Loads FlyBase data, maps `FBgn → FBal → FBtp/FBti → FBst`, fetches
  references, scores keywords, writes the Stage 1 workbook.

- **`StockSplittingPipeline`** (`pipelines/stock_splitting.py`): the organize,
  validate, *and* phenotype-sheet flows. Loads Stage 1 output, computes derived
  columns (including `PHENOTYPE_RELEVANCE_SCORE` when a filter needs it),
  applies JSON filters and stock limits, and writes organized workbooks. Pass
  `run_validation=True` to also run GPT validation. Set `Settings.soft_run=True`
  to switch to the phenotype reagent flow; combine with
  `enable_oai_embedding=True` to emit the cosine analysis and similarity-tier
  sidecar.

- **`Settings`** (`config.py`): Dataclass holding API keys, model names,
  paths, and feature flags (`soft_run`, `enable_oai_embedding`,
  `phenotype_similarity_sidecar`, `phenotype_similarity_plots` (gate the
  t-SNE/density plot folder independently of cosine scoring),
  `emit_no_pmid_report` (gate the no-PMID `.txt` sidecar),
  `organized_output_dir` (write the aggregated workbook directly into this
  directory instead of a nested `Organized Stocks`), `enable_gpt_logging`,
  `max_gpt_calls_per_stock`, `phenotype_similarity_targets`,
  `skip_fbgnid_conversion`). Loads `.env` at init time. Default model is
  `gpt-5-mini`; default embedding model is `text-embedding-3-large`.

- **`ValidationStatus`** (`config.py`): Validation status constants,
  category mapping, and priority ordering used by Stage 3 outputs.

- **`EmbeddingSimilarityScorer`** (`integrations/phenotype_similarity.py`):
  Computes phenotype embeddings against configured targets. It powers cosine
  columns, plots, and the similarity-tier sidecar when
  `settings.embeddings.enabled` is `true`.

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
- `logs/gpt_queries/` (gitignored) - GPT call logs when
  `settings.validation.enableGptLogging` is `true`

The code prefers the local `data/` directory and falls back to a sibling
`../Data` directory for legacy layouts (`is_portable_mode()` returns whether
the local layout is in use).
