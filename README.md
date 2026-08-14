# fl_ai_reagent_stocker

Modular Drosophila stock-processing pipeline. It turns gene lists into a
reference-linked stock workbook, organizes stocks into config-driven sheets
tiered by publication and curated-phenotype evidence, and optionally runs
OpenAI validation and phenotype-embedding analysis.

## Documentation

- Pipeline guide: `docs/pipeline_usage.md`
- Mermaid flowcharts: `docs/pipeline_flowcharts.md`
- Agent guidance: `docs/AGENTS.md`
- Citing data sources, APIs, and models: `docs/citations.md`
- Config reference: `docs/stock_split_config_example.md`
- Phenotypic-stock priority configs: `docs/stock_split_config_priority_example.md`

## What The Pipeline Produces

- a Stage 1 stock workbook linked to references plus a complete `Gene Reagent Index`
- organized output sheets defined by your JSON config
- a references sheet narrowed to papers cited by selected stocks
- `Ref++` validation columns when enabled by JSON config
- cosine-similarity columns when embeddings are enabled by JSON config
- a combined combination-counts summary at the input-folder root

## Quick Start

```bash
pip install -r requirements.txt

python -m fl_ai_reagent_stocker run ./gene_lists --config ./my_config.json
```

`run` scans the input folder recursively for gene-list CSVs and processes each
independently: it builds the Stage 1 workbook, organizes stocks with the JSON
config, validates `Ref++` stocks when enabled, and writes a combined summary at
the input root. JSON config controls input columns, FBgn conversion, PubMed
batch size, embeddings, output preservation, and validation policy. Default
config: `data/config/stock_split_config_example.json`.

Environment variables (optional, from a repo-root `.env`):

- `OPENAI_API_KEY` — required for validation and embeddings when enabled in config
- `NCBI_API_KEY` — higher PubMed rate limits
- `UNPAYWALL_TOKEN` — full-text contact email (defaults to `aadishms@umich.edu`)
- `OPENAI_MODEL`, `OPENAI_EMBEDDING_MODEL` — override default models

## `run` Options

- `--config`, `-c` — JSON split config
- `--quiet`, `-q` — suppress per-file progress logging

## Outputs

`run` produces a clean tree under the input folder:

- Per gene set: `./gene_lists/Per Gene Set Runs/<relative path>/<gene set>/Stocks/<gene set>_aggregated.xlsx`
- Combined summary at the input root: `combination_counts_summary.xlsx` / `.csv`

The primary `run` does not emit a nested `Organized Stocks` directory, `.txt`
sidecars, `*_similarity_tiers.xlsx` workbooks, `*_similarity/` plot folders, or
the intermediate unsplit `<gene set>.xlsx` unless
`settings.output.preserveUnsplitWorkbook` is `true`. The advanced
per-stage commands below retain the legacy `Stocks/Organized Stocks/` layout
and `phenotype-sheet` sidecar artifacts when embeddings are enabled in config.

## Configuration

Output sheets come from ordered `combinations` of named `filters`; each reagent
(`stock_id`, `collection`) lands in only the first combination it matches. An
optional `combination_names` list assigns corresponding Excel tab names
one-to-one. Populated categories are automatically prefixed `1_`, `2_`, and so
on; without names, they use legacy `Sheet1..N` names.
Two computed score columns back the standard tiers:

- `settings.input` — required input column and FBgn-conversion policy.
- `settings.pubmed` — required PubMed/full-text batch policy.
- `settings.embeddings` — required embedding-analysis policy.
- `settings.output` — required output-preservation policy.
- `settings.validation` — required explicit block in every config. Shipped
  configs set `enabled: false` so Ref++ validation is skipped even when
  `OPENAI_API_KEY` is available.
- `settings.validation.maxGptCallsPerStock` — required explicit cap for GPT
  validation calls per stock when validation is enabled (`null` means no cap).
- `ALLELE_PAPER_RELEVANCE_SCORE` → `Ref++` / `Ref+` / `Ref-`, from title/abstract
  keyword matches against `settings.relevantSearchTerms`.
- `PHENOTYPE_RELEVANCE_SCORE` → `Phenotype++` / `Phenotype+` / `0`, from FlyBase
  `genotype_phenotype_data`. `Phenotype++` requires a curated phenotype matching
  a `settings.phenotypeSimilarityTargets` term; `Phenotype+` is any curated
  phenotype. Computed only when a filter references it.

Config fields: `settings.relevantSearchTerms`,
`settings.phenotypeSimilarityTargets`, `settings.maxStocksPerGene`,
`settings.maxStocksPerAllele`, `settings.input`, `settings.pubmed`,
`settings.embeddings`, `settings.output`, `settings.validation`,
`settings.contentsSeparatorEvery`, `filters`, `combinations`,
`combination_names`, `filterDescriptions`.

Example configs:

- `data/config/stock_split_config_example.json` — publication-evidence tiers
- `data/config/stock_split_config_phenotype_example.json` — adds `Phenotype++`
  and `Phenotype+` sheets above each `Ref` group
- `data/config/stock_split_config_priority_example.json` — six named phenotype
  buckets; see `docs/stock_split_config_priority_example.md`
- `data/config/stock_split_config_priority_all_phenotypes.json` — exhaustive
  RNAi → UAS → allele coverage plus a phenotype audit bucket

See `docs/stock_split_config_example.md` for the full config reference.

## Advanced: per-stage entry points

The individual stages remain available for debugging or partial reruns:

```bash
python -m fl_ai_reagent_stocker find-stocks ./gene_lists
python -m fl_ai_reagent_stocker split-stocks ./gene_lists/Stocks --config ./my_config.json
python -m fl_ai_reagent_stocker phenotype-sheet ./gene_lists/Stocks --config ./my_config.json
python -m fl_ai_reagent_stocker validate-stocks ./gene_lists/Stocks --config ./my_config.json
```

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
- `phenotype-sheet` builds the gene-first All Phenotypic Stocks Sheet; `validate-stocks` runs the GPT validation step.
