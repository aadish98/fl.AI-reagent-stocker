# GUI Plan: Chat Summary

This folder summarizes the GUI specification decisions reached during chat
discussion. The canonical, executable plan lives at
[`/Users/aadishms/.cursor/plans/streamlit_gui_07f929cc.plan.md`](/Users/aadishms/.cursor/plans/streamlit_gui_07f929cc.plan.md).
These markdowns are derivative notes; the plan is the source of truth.

## Audience

Non-coding users (e.g., a neuroscience PhD student running fly RNAi screens)
who currently rely on Maddy's PowerShell scripts under
[`Maddy_Playground/RNAi plan script TEST/`](../Maddy_Playground/RNAi%20plan%20script%20TEST).
They should never touch CLI flags, JSON configs, PowerShell, or column-position
tricks.

## One-Sentence Goal

Take an input gene list, surface phenotype-backed stock candidates that hit the
user's keywords first, let the user shape priorities through plain-language
controls, and emit the same final outputs Maddy's scripts produce
(`RNAi Order Form.xlsx`, `RNAi screening plan.xlsx`) plus a visual
`gene_schema-flowchart.png` audit.

## File Index

| # | File | Topic |
|---|---|---|
| 01 | [`01-vision-and-priorities.md`](01-vision-and-priorities.md) | Goal, audience, primary endpoint, product priority |
| 02 | [`02-user-flow.md`](02-user-flow.md) | High-level user flow + screen list |
| 03 | [`03-config-and-keywords.md`](03-config-and-keywords.md) | Default config, default keywords, plain-language priorities (no raw JSON) |
| 04 | [`04-phenotype-results-browser.md`](04-phenotype-results-browser.md) | Keyword-hit first, no-hit secondary, fallback, optional non-phenotypic explorer |
| 05 | [`05-stock-planning-and-dedupe.md`](05-stock-planning-and-dedupe.md) | Stock planning stage, phenotype folding into reagents, duplicate audit, count path |
| 06 | [`06-rnai-outputs.md`](06-rnai-outputs.md) | Maddy outputs ported to Python, final docs, column retention, `<keyword> phenotype` columns |
| 07 | [`07-allada-snapshot.md`](07-allada-snapshot.md) | Allada Lab Order Sheet `.xlsx` upload (no Google auth), parser robustness, reference fixture |
| 08 | [`08-schema-flowchart.md`](08-schema-flowchart.md) | `gene_schema-flowchart.png` automation and checkpoint emission |
| 09 | [`09-data-privacy.md`](09-data-privacy.md) | Gitignore for lab-sensitive data, project gene sets, generated outputs |
| 10 | [`10-deployment.md`](10-deployment.md) | Streamlit MVP decision, when to graduate to other frameworks |
| 11 | [`11-dependencies-and-sequencing.md`](11-dependencies-and-sequencing.md) | Phenotype-sheet rework (done) and simplify-cli (pending) sequencing |

## Highest-Stakes Decisions In One Place

- The GUI is **phenotype-first**: it runs the `phenotype-sheet` flow by default
  and centers the UI on keyword-hit reagents.
- Default keywords: `sleep`, `rhythm`, `circadian`, `locomotor`.
- Default split config: `data/config/stock_split_config_no_bloomington.json`.
- OpenAI is **starred premium**, never required. Keyword hit/no-hit grouping
  must work without an API key. `Cosine Similarity (<keyword>)` columns are
  replaced by interpretable `<keyword> phenotype` columns in the default path.
- Lab-stock lookup uses an uploaded `.xlsx` snapshot of the Allada Lab Order
  Sheet. **No Google auth in MVP.**
- Two final user-facing workbooks: `RNAi Order Form.xlsx` and
  `RNAi screening plan.xlsx`, plus `gene_schema-flowchart.png`.
- Every run emits `gene_schema-flowchart.png` with all count checkpoints.
- Duplicate phenotypes are **folded into the kept reagent**, not silently
  dropped, and the dropped source rows are persisted to a duplicate-tracking
  audit workbook.
- Lab-sensitive workbooks, project-specific gene sets, and generated
  screen-planning outputs must be gitignored.

## Related Plans

- Canonical GUI plan: [`/Users/aadishms/.cursor/plans/streamlit_gui_07f929cc.plan.md`](/Users/aadishms/.cursor/plans/streamlit_gui_07f929cc.plan.md)
- Phenotype-sheet rework (done): [`/Users/aadishms/.cursor/plans/phenotype_sheet_rework_844d6a5c.plan.md`](/Users/aadishms/.cursor/plans/phenotype_sheet_rework_844d6a5c.plan.md)
- CLI simplification (pending): [`/Users/aadishms/.cursor/plans/simplify-cli_b16a2d0c.plan.md`](/Users/aadishms/.cursor/plans/simplify-cli_b16a2d0c.plan.md)
