# AGENTS.md

This file provides guidance to coding agents working in this repository.

Repo-root `AGENTS.md` carries the rules for handling downstream requests on
stocker output (filtering, selection, counting, and exports on reagents a run
has already produced). Follow it whenever a request is about an existing run's
reagents rather than about the pipeline itself.

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
- `docs/features.md` — the running list of features that work on reagents a
  stocker run already produced: what each one reads, the columns and rule for
  every filtering or selection step in order, and what it writes. Add an entry
  with every new `scripts/` feature, and update it in the same change that
  changes the script's steps.
- `AGENTS.md` (repo root) — rules for downstream requests on stocker output,
  plus the wording conventions for step comments and `docs/features.md`.
- `docs/SOP_Stocker_Cursor_CLI.pdf` — standalone SOP: GitHub install, how
  Stocker works, solid gene-list review, default config, and stock limits.
  Keep the FBgn conversion hard rule in lockstep with this file.

## FBgn ID conversion (hard rule)

Wrong FBgn IDs are the largest failure mode in this system.

Gene symbol → FlyBase gene ID (`FBgn` / `flybase_gene_id`) conversion
**must always** go through `scripts/fetch_fbgn_ids.py`. That script
**never overwrites** the user CSV. It writes:

- `validated_<original>.csv` — mapped rows only
- `needs-review.csv` — unmatched rows from every input, plus `source_file`
- `validated_<original>.xlsx` — review workbook; each FBgn is a hyperlink
  to `https://flybase.org/reports/<FBgn>`

Agents must **never** edit the user-generated gene list. Agents must
**never** invent, recall, autocomplete, or hand-fill FBgn IDs. Agents
may edit `validated_*.csv` or FBgn columns **only if the user explicitly
asks**. Tell the user the xlsx exists and wait for them to review it for
precision. If `needs-review.csv` has data rows, send them to
https://flybase.org/convert/id first. Run stocker **only** on
`validated_*.csv`, and only after that review. `run` skips
`needs-review.csv` and prefers a `validated_` sibling when one exists.

Shipped configs may set `settings.input.skipFbgnidConversion` to `true`
(including `stock_split_config_priority_example.json`). In that case
stocker trusts the validated CSV IDs and will not convert symbols during
`run`.

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
python -m fl_ai_reagent_stocker run ./gene_lists --config data/config/stock_split_config_priority_example.json [--quiet]
```

Unless the user specifies another config, use
`data/config/stock_split_config_priority_example.json`.

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

When the JSON config uses a `Phenotype` filter, `run` writes the gene-first
`All Phenotypic Stocks Sheet` and adds phenotype/qualifier/reference/partner-GAL4
enrichment to the matching combination sheets. This does not require embeddings.
`settings.embeddings.enabled` additionally runs OpenAI embeddings to add cosine
columns to the aggregated workbook. For the primary `run`, the similarity tier
sidecar and t-SNE/plot folders are suppressed; the input is scanned recursively
for gene-set CSVs; the final aggregated workbook is written directly into each
gene set's `Stocks` folder; and a combined combination-counts summary is
generated at the input root when the run finishes.
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
`data/config/stock_split_config_priority_example.json`; use it unless another
config is explicitly defined. `stock_split_config_example.json` is the
publication-evidence baseline; `stock_split_config_phenotype_example.json` is
the phenotype-aware variant. `stock_split_config_priority_all_phenotypes.json`
is the exhaustive phenotype-bucket variant (see
`docs/stock_split_config_priority_example.md`).

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
  reagent is placed in only the first combination it matches. When any filter
  or combination name starts with `Phenotype`, combination sheets also receive
  phenotype/qualifier/reference/GAL4 enrichment columns from the internal
  phenotype sheet (`_add_phenotype_enrichment_to_stock_sheet`), and the
  full-scope `All Phenotypic Stocks Sheet` is emitted even when embeddings are
  disabled
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

### Reagent-list ingredients

The public, tracked config is
`data/config/stock_split_config_priority_example.json`; it is the only file in
`data/config/` that users receive from a fresh clone. Treat its `filters`
object as the catalog for curating a new reagent list, then create a
priority-ordered `combinations` entry from those names.

Filterable columns and ingredients:

- **Source and reagent class:** `collection`, `custom_stock`, `RNAi`, and
  `transgenic_product_class_terms` support Bloomington/Vienna, stock-center,
  RNAi, and allele/insertion rules.
- **Stock structure:** `sgRNA`, `num_Balancers`, and
  `multiple_insertions` are derived before filters run and support CRISPR,
  balancer, and insertion-count rules.
- **Evidence:** `ALLELE_PAPER_RELEVANCE_SCORE` is always derived and supports
  `Ref++` (2), `Ref+` (1), and `Ref-` (0). `PHENOTYPE_RELEVANCE_SCORE` is
  derived only when a filter uses it and supports `Phenotype` or the finer
  `Phenotype++`/`Phenotype+` tiers.
- **Supported predicate types:** `contains`, `doesnt_contain`, `equals`,
  `doesnt_equal`, and `insertion_site`. `insertion_site` additionally requires
  a single-insertion stock.

The shipped catalog includes collection (`Bloomington`, `Vienna`), stock-source
(`Stock Center`, `Non Stock Center`), reagent (`RNAi_reagent`, `UAS`,
`AlleleOrInsertion`), evidence (`Phenotype`, `Ref++`, `Ref+`, `Ref-`), and
structure (`No sgRNA`, `Has Balancers`, `No Balancers`, `other_insertion`,
`no_other_insertion`) ingredients. A defined filter need not be used by the
default combinations; it remains available for user-authored combinations.

Combinations use first-match-wins priority. `maxStocksPerGene` and
`maxStocksPerAllele` are also enforced in that order across the full workbook.
With the default 5/3 caps, phenotype sheets can consume a gene or allele's
selection budget before later Ref sheets are considered.

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
│                                # reagent-bucket / RNAi-type helpers, shared
│                                # GAL4 symbol primitives
├── validation_runner.py         # Shared GPT validation logic
├── integrations/
│   ├── __init__.py
│   ├── pubmed.py                # PubMed/Entrez client and PubMedCache
│   ├── fulltext.py              # Full-text retrieval + FunctionalValidator
│   └── phenotype_similarity.py  # OpenAI embedding cosine scoring + t-SNE plots
├── pipelines/
│   ├── __init__.py
│   ├── stock_finding.py         # Stage 1: genes → stocks → references → Gene Reagent Index
│   └── stock_splitting.py       # organize, phenotype-sheet, validate:
│                                # filters, limits, Excel output, gene-first
│                                # phenotype reagent flow + embedding analysis
└── tx_omics_revision/           # Separate Tx-Omics Revision analysis (see below)
    ├── __init__.py              # Exposes run_pipeline
    ├── pipeline.py              # Stage orchestration + run manifest
    ├── constants.py             # Pinned paths, hashes, DAVID settings
    ├── categories.py            # The four publication-defined category CSVs
    ├── flybase.py               # Synonym-aware FlyBase identity resolution
    ├── go_analysis.py           # DAVID GO enrichment via hash-pinned scripts
    ├── keywords.py              # Keyword matching + DAVID term parsing
    ├── pathways.py              # FDR-filtered pathway-hit gene sets
    ├── overlap.py               # FBgn overlap tables + narrative report
    ├── plots.py                 # Deterministic figures
    └── io_utils.py              # Filesystem, hashing, git helpers
```

The GAL4 symbol primitives in `utils.py` (`looks_like_gal4_symbol`,
`normalized_text_looks_like_gal4_driver`, `gal4_search_keys`,
`gal4_driver_fingerprint`, `is_split_gal4_component`,
`primary_gal4_component_symbol_for_fbst`, `normalize_symbol_token`,
`compact_symbol_key`, `split_genotype_symbol_tokens`) are shared by
`pipelines/stock_splitting.py` and `scripts/export_gal4_driver_workbook.py`.
Change them in `utils.py`, not by re-defining a local copy in either consumer.

### Tx-Omics Revision subpackage

`fl_ai_reagent_stocker/tx_omics_revision/` is a self-contained analysis for a
manuscript revision, independent of the stocker pipeline and not reachable from
the `fl_ai_reagent_stocker` CLI. It builds four publication-defined gene-set
categories, resolves identities against the pinned FlyBase synonym release,
runs DAVID GO enrichment through hash-pinned external scripts, derives
FDR-filtered pathway gene sets, computes FBgn overlap tables, and emits
figures. `pipeline.run_pipeline` executes these as named stages (`preflight`,
`categories`, `category_audit`, `go`, `review_terms`, `pathways`,
`pathway_audit`, `overlap`, `plots`, `docs`) and records input hashes, package
versions, and git state in a run manifest for reproducibility. Drive it through
the `scripts/*_tx_omics_revision*.py` entry points below.

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

Canonical helper scripts live in `scripts/`. Run them from the repo root;
data-building scripts may import shared package utilities when they need the
same FlyBase, PubMed, or configuration logic as the pipeline. A script that
works on reagents a stocker run already produced also needs a test built from
real rows and an entry in `docs/features.md` (see repo-root `AGENTS.md`).

Data builders / refreshers:

- `refresh_flybase_data.py` - pull the latest FlyBase TSV/TSV.GZ reports.
- `refresh_bdsc_data.py` - download the native BDSC datasets into `data/bdsc/`,
  used by the GAL4 workbook exporter.
- `fetch_fbgn_ids.py` - convert gene symbols to FBgn IDs. Never overwrites
  the user CSV. Writes `validated_*.csv`, `needs-review.csv`, and a
  hyperlinked `validated_*.xlsx`. The only allowed way to create
  `flybase_gene_id` values. Agents must never fill those IDs by hand.
- `build_sop_stocker_pdf.py` - regenerate `docs/SOP_Stocker_Cursor_CLI.pdf`.
- `render_config_walkthrough_pdf.py` - render `docs/stock_split_config_example.md`
  to `docs/stock_split_config_example.pdf` (gitignored output).
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
- `summarize_combination_counts.py` - aggregate the "Prioritized sheet
  breakdown" tables from organized workbook `Contents` sheets across a
  directory tree, one row per gene set and config combination. Backs the
  `combination_counts_summary.xlsx` / `.csv` that `run` writes at the input root.
- `build_temperature_phenotype_gene_list.py` - match all FlyBase curated
  phenotype/qualifier text against Temperature/heat/cold/thermo keywords
  (direct text match, no embeddings), attach orderable and custom reagents,
  enrich citations through the PubMed cache, and write a gene-grouped Excel
  workbook with `Reagents` and `References` sheets. The inclusion rule is a
  cited phenotype with a boundary-aware keyword hit AND an allele/insertion
  reagent (UAS/RNAi constructs, including common VDRC knockdown allele
  families, are excluded). The `Reagents` sheet carries a `Stock Center`
  column (`BDSC`, `Vienna`, other FlyBase collection names, or `Custom`).
  The ubiquitous `w` gene (`FBgn0003996`) is manually blacklisted by both
  FBgn ID and symbol (`exclude_genes`), applied once before stock/reference
  expansion and again on the final deduplicated reagent table, so its
  `w[1118]` background allele can never generate downstream reagent or
  reference rows.

Exports and workbook post-processing:

- `export_gal4_driver_workbook.py` - take a CSV with a `GAL4 symbol` column,
  resolve orderable GAL4 driver stocks (including mixed-genotype Bloomington
  lines and cross-format GawB naming) through FlyBase derived stock components,
  and write one sheet per stock collection plus `References` and `Coverage`
  sheets. Shares the GAL4 symbol primitives in `utils.py` with
  `stock_splitting.py`; keep those two in sync.
- `refactor_aggregated_workbook_layout.py` - improve readability of an
  aggregated stock workbook (column widths, freezing, formatting) without
  changing any cell values.

Tx-Omics Revision (see the subpackage section above):

- `build_tx_omics_revision_sets.py` - run the full staged pipeline: gene sets,
  audits, overlap tables, and figures.
- `run_tx_omics_revision_go.py` - run only the hash-pinned DAVID GO enrichment
  stage.
- `plot_tx_omics_revision.py` - regenerate figures from saved evidence tables
  without re-running the analysis.

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
- `tests/test_config_defaults.py` - guards that the default split-config path
  points at the git-tracked config
- `tests/test_export_gal4_driver_workbook.py`
- `tests/test_tx_omics_revision.py`

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
- `config/*.json` — default is `stock_split_config_priority_example.json`
  (tracked); other JSON configs in this directory remain gitignored unless
  another config is explicitly defined for a run
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
