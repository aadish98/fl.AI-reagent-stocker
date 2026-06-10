# 11. Dependencies and Sequencing

The GUI depends on two upstream plans. One is implemented; one is pending.
This file captures the order of operations and what each upstream plan
guarantees for the GUI.

## Upstream Plans

| Plan | Status | Purpose for the GUI |
|---|---|---|
| [`/Users/aadishms/.cursor/plans/phenotype_sheet_rework_844d6a5c.plan.md`](/Users/aadishms/.cursor/plans/phenotype_sheet_rework_844d6a5c.plan.md) | Implemented | Provides the gene-first `All Phenotypic Stocks Sheet` and the persisted Stage 1 `Gene Reagent Index`. |
| [`/Users/aadishms/.cursor/plans/simplify-cli_b16a2d0c.plan.md`](/Users/aadishms/.cursor/plans/simplify-cli_b16a2d0c.plan.md) | Pending | Cleans up the public CLI surface so the GUI can rely on a single, clearer set of commands and flag semantics. |

## Phenotype-Sheet Rework: What The GUI Now Inherits

Implementation is verified in
[`fl_ai_reagent_stocker/pipelines/stock_finding.py`](../fl_ai_reagent_stocker/pipelines/stock_finding.py)
and
[`fl_ai_reagent_stocker/pipelines/stock_splitting.py`](../fl_ai_reagent_stocker/pipelines/stock_splitting.py).

What this gives the GUI:

- `All Phenotypic Stocks Sheet` is **gene-first**: every input-gene FBal/FBtp/FBti
  reagent is considered, regardless of whether a stock survived earlier stages.
- A `Gene Reagent Index` sheet is persisted in `aggregated_stock_refs.xlsx`
  for every input gene.
- No-stock phenotype reagents are emitted from the index when no FBst lookup
  match exists, so the keyword-hit/no-hit views are not narrowed by stock
  survival.
- Stage 2 already loads the index and threads it into the phenotype-sheet
  builder; legacy fallback exists for older Stage 1 workbooks.

GUI consequences:

- The keyword-hit/no-hit browser is more complete than it would have been
  before the rework.
- The "Input Genes Needing Fallback" panel is meaningful: a missing keyword
  hit really does mean the user has follow-up work to do, not that the
  pipeline silently dropped data.
- The optional **Non-Phenotypic Reagent Explorer**
  ([`04-phenotype-results-browser.md`](04-phenotype-results-browser.md))
  reads the same `Gene Reagent Index` to surface input-gene reagents with no
  curated phenotype row.

## Simplify-CLI: What The GUI Will Inherit

After the simplify-cli plan executes, the public CLI surface becomes:

```text
fl-ai-reagent-stocker find-stocks         GENE_LIST_DIR --config ...
fl-ai-reagent-stocker split-stocks        STOCKS_DIR    --config ...
fl-ai-reagent-stocker phenotype-sheet     STOCKS_DIR    --config ... [--similarity tiers|simple-buckets|keyword-buckets] [--embeddings|--no-embeddings]
fl-ai-reagent-stocker validate-stocks     STOCKS_DIR    --config ...
fl-ai-reagent-stocker run-full-pipeline   GENE_LIST_DIR --config ... --mode normal|phenotype
```

GUI consequences:

- The GUI calls the existing Python APIs directly (not via subprocess), so
  command-name changes don't break the GUI directly. But the `--similarity`
  enum makes documentation, smoke checks, and any optional
  shell-out commands much simpler.
- The GUI's default phenotype-first run cleanly maps to:
  - Equivalent of `phenotype-sheet ... --no-embeddings` when OpenAI is off.
  - Equivalent of `phenotype-sheet ... --similarity keyword-buckets` if the
    user wants the keyword-hit/no-hit sidecar workbook without OpenAI.
  - Equivalent of `phenotype-sheet ... --similarity tiers --embeddings` when
    the user explicitly wants cosine-tier output.
- Internal `Settings` flags (`soft_run`, `enable_oai_embedding`,
  `simple_buckets`, `keyword_bucketing`, `phenotype_similarity_sidecar`) are kept as the implementation
  substrate; the GUI does not have to care about them once the simplify-cli
  mapping is in place.

## Recommended Build Order

1. **Phenotype-sheet rework** — done.
2. **Simplify-CLI** — pending. Land first if possible so docs and tests align.
3. **GUI MVP** —
   1. App shell, setup checks, gene input
   2. Default config + keyword/priority screen
   3. Run pipeline (default phenotype-first)
   4. Phenotype-first results browser
   5. Stock planning stage with phenotype folding and duplicate audit
   6. RNAi planning Python port (the canonical Maddy outputs)
   7. Allada `.xlsx` upload + parser
   8. `gene_schema-flowchart.png` automation
   9. Downloads bundle
   10. Optional non-phenotypic reagent explorer
   11. Lab data privacy gitignore + docs

## Risks Worth Watching

- **Phenotype sheet row count growth**: gene-first behavior surfaces strictly
  more rows than the legacy narrowing. The GUI should default to keyword-hit
  + filtered views, not full-table dumps.
- **Schema regression for Maddy's PowerShell scripts**: those scripts read
  `All Phenotypic Stocks Sheet` by name and column shape; column-name changes
  must keep aliases until the Python RNAi planning module replaces the
  PowerShell flow end-to-end.
- **Streamlit fragility for senior users**: trade-off accepted for MVP; track
  user feedback to decide when to rebuild on FastAPI/React.
- **Allada parser coverage**: today's PowerShell parser only picks one
  stock-ID column per worksheet, missing several inventory boxes. The Python
  rewrite must be strictly more thorough — see
  [`07-allada-snapshot.md`](07-allada-snapshot.md).
