# fl_ai_reagent_stocker

Modular Drosophila stock-processing pipeline for turning gene lists into stock workbooks, config-driven organized sheets, an optional gene-first phenotype reagent sheet, and optional AI-assisted validation.

## Documentation

- Practical guide: `docs/pipeline_usage.md`
- Editable Mermaid flowcharts: `docs/pipeline_flowcharts.md`
- Agent guidance: `docs/AGENTS.md`
- How to cite the data sources, APIs, and models this pipeline uses: `docs/citations.md`
- Example split config: `data/config/stock_split_config_example.json`

## What The Pipeline Produces

Depending on which command(s) you run, you get:

- a Stage 1 stock-level workbook linked to references and a complete `Gene Reagent Index`
- organized output sheets based on your JSON rules
- a references sheet narrowed to papers cited by selected stocks
- a gene-first `Stock Phenotype Sheet` with optional similarity sidecar workbook
- optional validation columns for `Ref++` output-sheet stocks

## Pipeline Overview

The canonical package is `fl_ai_reagent_stocker` and it exposes five explicit commands:

1. `find-stocks`
   - Reads gene-list CSV files.
   - Converts gene symbols to FBgn IDs unless skipped.
   - Maps genes to FlyBase stocks through allele, construct, and insertion components.
   - Links stocks and components to publications and PMIDs.
   - Persists a complete `Gene Reagent Index` of every FBal/FBtp/FBti reagent associated with each input gene.
   - Writes `./gene_lists/Stocks/aggregated_stock_refs.xlsx`.

2. `split-stocks`
   - Reads the Stage 1 workbook.
   - Applies `data/config/*.json` filters and combinations.
   - Writes organized output sheets to `./gene_lists/Stocks/Organized Stocks/`.
   - Does not produce the phenotype sheet and does not run GPT validation.

3. `phenotype-sheet`
   - Reads the Stage 1 workbook.
   - Drives the gene-first `Stock Phenotype Sheet` from the `Gene Reagent Index`, attaching stock metadata afterward from the full `fbst_to_derived_stock_component.csv`.
   - Optionally adds a similarity sidecar workbook through `--similarity tiers|simple-buckets|keyword-buckets`.

4. `validate-stocks`
   - Rebuilds the same split outputs using the same JSON config.
   - Selectively validates only `Ref++` output-sheet stocks.
   - Merges validation columns back into the organized workbook.

5. `run-full-pipeline`
   - End-to-end wrapper from raw gene-list CSVs.
   - `--mode normal` runs Stage 1 + `split-stocks` + `validate-stocks`.
   - `--mode phenotype` runs Stage 1 + `phenotype-sheet` only and never invokes GPT validation.

The JSON configs in `data/config/` stay in place and keep the same effect on stock splitting as before.

## Installation

From the repository root:

```bash
pip install -r requirements.txt
```

Optional environment variables:

- `OPENAI_API_KEY` for validation and OpenAI embeddings
- `NCBI_API_KEY` for improved PubMed rate limits
- `UNPAYWALL_TOKEN` for improved full-text retrieval coverage
- `OPENAI_MODEL` to override the default model
- `OPENAI_EMBEDDING_MODEL` to override the default embedding model

## Quick Start

### One command, normal mode (split + validate)

```bash
python -m fl_ai_reagent_stocker run-full-pipeline ./gene_lists --config ./my_split_config.json --mode normal --test-log
```

### One command, phenotype mode (Stage 1 + phenotype sheet)

```bash
python -m fl_ai_reagent_stocker run-full-pipeline ./gene_lists --config ./my_split_config.json --mode phenotype --similarity tiers
```

### Stage-by-stage workflow

```bash
python -m fl_ai_reagent_stocker find-stocks ./gene_lists
python -m fl_ai_reagent_stocker split-stocks ./gene_lists/Stocks --config ./my_split_config.json
python -m fl_ai_reagent_stocker phenotype-sheet ./gene_lists/Stocks --config ./my_split_config.json --similarity tiers
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

The workbook contains three sheets: `Stocks`, `References`, and `Gene Reagent
Index`. The `Gene Reagent Index` sheet is a flat, complete catalog of every
FBal / FBtp / FBti reagent associated with each input gene, independent of
stock availability or downstream JSON split limits. The `phenotype-sheet`
command uses it as the reagent universe for the Stock Phenotype Sheet so
phenotype rows are no longer narrowed by Stage 1 stock matching or ranking.

Common options:

- `--config`, `-c`
- `--skip-fbgnid-conversion`
- `--gene-col`
- `--input-gene-col`
- `--batch-size`

### Stage 2: `split-stocks` (organization only)

```bash
python -m fl_ai_reagent_stocker split-stocks ./gene_lists/Stocks --config ./my_split_config.json
```

Main output folder:

- `./gene_lists/Stocks/Organized Stocks/`

`split-stocks` only performs JSON-driven stock organization. It no longer
accepts `--soft-run`, `--OAI-embedding`, `--simple-buckets`, or
`--keyword-bucketing`; use `phenotype-sheet` for that workflow.

Common options:

- `--config`, `-c`
- `--quiet`, `-q`

### `phenotype-sheet` (gene-first phenotype reagent flow)

```bash
python -m fl_ai_reagent_stocker phenotype-sheet ./gene_lists/Stocks --config ./my_split_config.json --similarity tiers
```

Main output folder:

- `./gene_lists/Stocks/Organized Stocks/`

`phenotype-sheet` consumes Stage 1 workbook directories such as
`./gene_lists/Stocks` and expects a `Gene Reagent Index` sheet for full-scope
output. Workbooks without the index still run via a legacy fallback that is
best-effort and not true full scope; the command emits a warning to stderr in
that case (even under `--quiet`) telling you to re-run `find-stocks` to
regenerate the index.

Common options:

- `--config`, `-c`
- `--quiet`, `-q`
- `--similarity {none,tiers,simple-buckets,keyword-buckets}` (default: `none`)

`--similarity` choices:

- `none` writes only the Stock Phenotype Sheet inside the aggregated workbook.
- `tiers` adds `<input>_aggregated_similarity_tiers.xlsx` with cosine-threshold
  tier tabs (0.05-wide buckets) computed against
  `settings.phenotypeSimilarityTargets`.
- `simple-buckets` adds the same sidecar with rule-based
  `collection / UAS / sleep-circ / balancer` tabs instead of cosine tiers.
- `keyword-buckets` adds the same sidecar with `Keyword Hits / No Keyword Hits`
  tabs instead of cosine tiers.

### Stage 3: `validate-stocks` (GPT validation only)

```bash
python -m fl_ai_reagent_stocker validate-stocks ./gene_lists/Stocks --config ./my_split_config.json --test-log
```

Validation is selective by design:

- only `Ref++` output-sheet stocks are considered
- references must be keyword-relevant
- validation short-circuits after the first functional hit per stock

`validate-stocks` no longer accepts `--soft-run`, `--OAI-embedding`,
`--simple-buckets`, or `--keyword-bucketing`; use `phenotype-sheet` for that
workflow.

Common options:

- `--config`, `-c`
- `--quiet`, `-q`
- `--test-log`
- `--max-gpt-calls-per-stock`

### `run-full-pipeline`

```bash
python -m fl_ai_reagent_stocker run-full-pipeline ./gene_lists --config ./my_split_config.json --mode normal --test-log
python -m fl_ai_reagent_stocker run-full-pipeline ./gene_lists --config ./my_split_config.json --mode phenotype --similarity tiers
```

`--mode normal` runs Stage 1 + `split-stocks` + `validate-stocks`.
`--mode phenotype` runs Stage 1 + `phenotype-sheet` only and never invokes
GPT validation. `--similarity` is honored only in `--mode phenotype`.
Validation flags (`--test-log`, `--max-gpt-calls-per-stock`) are honored
only in `--mode normal`.

Common options:

- `--config`, `-c`
- `--quiet`, `-q`
- Stage 1: `--gene-col`, `--input-gene-col`, `--batch-size`, `-b`, `--skip-fbgnid-conversion`
- `--mode {normal,phenotype}` (default: `normal`)
- `--similarity {none,tiers,simple-buckets,keyword-buckets}` (default: `none`)
- `--test-log`, `--max-gpt-calls-per-stock`

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

`refresh_flybase_data.py` downloads the precomputed TSV report families by
default. Pass `--with-xml` to additionally refresh the optional chado XML
families (`chado_FBst*.xml(.gz)`, `chado_FBti*.xml(.gz)`) that the build
helpers consume. The build and audit scripts discover the local XML files by
prefix, so they remain compatible if FlyBase ever embeds a release suffix in
chado XML filenames.

## Notes

- Stage 1 supports FBtp -> FBti expansion via `data/flybase/transgenic_insertions/fbtp_to_fbti.csv`.
- `phenotype-sheet` stops before GPT execution. Use `validate-stocks` to run GPT validation.
- Output paths remain unchanged across the CLI cleanup.
