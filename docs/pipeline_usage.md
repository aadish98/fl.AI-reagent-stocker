# fl_ai_reagent_stocker: Pipeline Guide

`fl_ai_reagent_stocker` converts Drosophila gene lists into organized,
reference-linked stock sheets. It maps genes to FlyBase stocks, tiers them by
publication and curated-phenotype evidence using a JSON config, and optionally
runs OpenAI validation and phenotype-embedding analysis when enabled by config.

## Inputs

- One folder of gene-list CSV files (for example `./gene_lists`). Each CSV
  provides gene symbols and/or FBgn IDs.
- A local FlyBase data tree under `data/flybase/` (see
  [stock_split_config_example.md](stock_split_config_example.md) and the helper
  scripts below).
- A JSON split config (defaults to
  `data/config/stock_split_config_example.json`).

Environment variables (optional, read from a repo-root `.env`):

- `OPENAI_API_KEY` — required for embeddings and Ref++ validation when enabled.
- `NCBI_API_KEY` — higher PubMed rate limits.
- `UNPAYWALL_TOKEN` — contact email for full-text retrieval.
- `OPENAI_MODEL`, `OPENAI_EMBEDDING_MODEL` — override default models.

Install dependencies once:

```bash
pip install -r requirements.txt
```

## Command

A single command runs the pipeline end to end:

```bash
python -m fl_ai_reagent_stocker run ./gene_lists --config ./my_config.json
```

The input folder is scanned recursively: every gene-list CSV found under it
(at any depth) is processed independently. A single config-provided
`settings.input.geneCol` / `settings.input.inputGeneCol` pair applies to the
whole run; if any discovered CSV is missing those columns, the run aborts before
doing heavy work.

`run` performs these steps per gene set:

1. Build the Stage 1 stock/reference workbook from the gene list.
2. Organize stocks into output sheets defined by the JSON config.
3. Validate `Ref++` output stocks when `settings.validation.enabled` is `true`
   (requires `OPENAI_API_KEY`).

Then, once all gene sets finish, it writes a combined combination-counts
summary at the input-folder root.

Set `settings.embeddings.enabled` to `true` to add OpenAI cosine-similarity
columns to each aggregated workbook. The primary `run` deliberately does **not**
emit the similarity tier sidecar workbook or the t-SNE/density plot folders
(those remain available via the advanced `phenotype-sheet` command when
embeddings are enabled in config).

Options:

- `--config`, `-c` — JSON split config (default:
  `data/config/stock_split_config_example.json`).
- `--quiet`, `-q` — suppress per-file progress logging.

## Outputs

`run` produces a clean, summary-ready tree under the input folder:

- Per gene set: `./gene_lists/Per Gene Set Runs/<relative path>/<gene set>/Stocks/<gene set>_aggregated.xlsx`
  (`Contents`, configured combination tabs or legacy `Sheet1..N`, `References`, `Stock Sheet by Gene`, plus Ref++
  validation columns when enabled, plus embedding cosine columns when enabled).
- Combined summary at the input root: `combination_counts_summary.xlsx` and
  `combination_counts_summary.csv` (one row per gene set and config
  combination, grouped by input folder).

The primary `run` intentionally does not emit:

- a nested `Organized Stocks` wrapper directory (the final workbook lives
  directly in `Stocks`),
- `.txt` sidecar reports (e.g. `references_without_pmid_fbrf.txt`),
- `*_similarity_tiers.xlsx` workbooks or `*_similarity/` plot folders,
- the intermediate unsplit `<gene set>.xlsx` workbook unless
  `settings.output.preserveUnsplitWorkbook` is `true`.

The advanced per-stage commands below retain the legacy
`Stocks/Organized Stocks/` layout and sidecar artifacts.

## Configuration

The JSON config drives stock organization. Output sheets come from
`combinations` of named `filters`, applied in order; each reagent
(`stock_id`, `collection`) is placed in only the first combination it matches.

Key fields:

- `settings.relevantSearchTerms` — keywords for the publication-evidence tiers
  (`Ref++` / `Ref+` / `Ref-`), matched case-insensitively against titles and
  abstracts.
- `settings.phenotypeSimilarityTargets` — `{keyword, embedding_text}` targets
  used for the curated-phenotype tiers and for embedding analysis.
- `settings.maxStocksPerGene`, `settings.maxStocksPerAllele` — selection caps.
- `settings.input` — required gene-column and FBgn-conversion policy.
- `settings.pubmed` — required PubMed/full-text batch-size policy.
- `settings.embeddings` — required embedding-analysis policy.
- `settings.output` — required output preservation policy.
- `settings.validation` — required explicit block in every config. Shipped
  configs set `enabled: false` so the Ref++ validation pass is skipped even when
  `OPENAI_API_KEY` is available.
- `settings.validation.maxGptCallsPerStock` — required explicit cap for GPT
  validation calls per stock when validation is enabled (`null` means no cap).
- `settings.validation.minFullTextChars`,
  `settings.validation.gptCallDelaySeconds`,
  `settings.validation.shortCircuitOnFunctionalValidation`, and
  `settings.validation.enableGptLogging` — required GPT-validation behavior
  controls.
- `settings.contentsSeparatorEvery` — draw a faint separator after this many
  rows in the `Contents` sheet breakdown table (default: `3`).
- `filters` — reusable column predicates.
- `combinations` — ordered sheet definitions.
- `combination_names` — optional Excel tab names aligned one-to-one with
  `combinations`. Names must be unique, valid Excel worksheet names of at most
  31 characters after numbering. Store unnumbered labels: populated categories
  are automatically prefixed `1_`, `2_`, and so on in output order. Empty
  combinations show `-` in Contents and create no tab; without this field,
  populated combinations use `Sheet1..N`.
- `filterDescriptions` — text shown in the `Contents` sheet.

### Evidence tiers

Two computed score columns back the standard tiers:

- `ALLELE_PAPER_RELEVANCE_SCORE` → `Ref++` (2), `Ref+` (1), `Ref-` (0),
  from title/abstract keyword matches.
- `PHENOTYPE_RELEVANCE_SCORE` → `Phenotype++` (2), `Phenotype+` (1), `0`,
  from FlyBase `genotype_phenotype_data`. A stock scores `2` when one of its
  components has a curated phenotype whose name matches a
  `phenotypeSimilarityTargets` term, `1` when it has any curated phenotype, and
  `0` otherwise. `PHENOTYPE_RELEVANCE_SCORE` is computed only when a filter
  references it.

Phenotype evidence is row-level over `(stock, collection, phenotype,
reference)`, but each score collapses to one value per stock so organized
sheets stay reagent-level over `(stock_id, collection)`.

### Example configs

- `data/config/stock_split_config_example.json` — publication-evidence tiers.
- `data/config/stock_split_config_phenotype_example.json` — adds `Phenotype++`
  and `Phenotype+` sheets above each `Ref` group, so phenotype-backed reagents
  are claimed before the publication-only tiers.
- `data/config/stock_split_config_priority_example.json` — six named,
  priority-ordered phenotype buckets. See
  [stock_split_config_priority_example.md](stock_split_config_priority_example.md).
- `data/config/stock_split_config_priority_all_phenotypes.json` — exhaustive
  RNAi → residual UAS → allele/insertion coverage, plus a final phenotype
  audit bucket.

## Finding phenotypic stocks

Use the priority configs when the goal is phenotype-backed reagents rather
than the publication-evidence Ref matrix. Combinations are evaluated top to
bottom; each `(stock_id, collection)` reagent lands in only the first matching
sheet.

```bash
python -m fl_ai_reagent_stocker run ./gene_lists \
  --config data/config/stock_split_config_priority_example.json
```

For exhaustive RNAi → UAS → allele coverage, use
`data/config/stock_split_config_priority_all_phenotypes.json`.

## Advanced: per-stage entry points

The stages are also available individually for debugging or partial reruns.
They share the same config and outputs as `run`:

```bash
python -m fl_ai_reagent_stocker find-stocks ./gene_lists
python -m fl_ai_reagent_stocker split-stocks ./gene_lists/Stocks --config ./my_config.json
python -m fl_ai_reagent_stocker phenotype-sheet ./gene_lists/Stocks --config ./my_config.json
python -m fl_ai_reagent_stocker validate-stocks ./gene_lists/Stocks --config ./my_config.json
```

`phenotype-sheet` builds the gene-first `All Phenotypic Stocks Sheet` from the Stage 1
`Gene Reagent Index`; `settings.embeddings.enabled` adds the cosine analysis and
tier sidecar.

## Helper scripts

Standalone scripts under `scripts/` (run from the repo root):

- `fetch_fbgn_ids.py` — gene symbols to FBgn IDs.
- `build_fbst_derived_stock_components.py` — build the derived stock-component CSV.
- `build_fbtp_to_fbti_mapping.py` — build the FBtp→FBti mapping CSV.
- `refresh_flybase_data.py` — download the latest FlyBase reports into
  `data/flybase/` and prune older versions:

```bash
python scripts/refresh_flybase_data.py            # TSVs only (default)
python scripts/refresh_flybase_data.py --with-xml # also refresh chado XML
```

The chado XML families are required by the two build scripts above.
