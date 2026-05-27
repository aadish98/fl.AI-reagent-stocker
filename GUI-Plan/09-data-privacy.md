# 09. Lab-Sensitive Data Privacy

## Principle

Anything that reveals in-progress lab screens or lab inventory stays out of
git. Schema documentation and synthetic test fixtures are fair game.

## What Must Be Gitignored

- Uploaded Allada Lab Order Sheet snapshots and any local lab inventory /
  order workbook, e.g.:
  - `Allada Lab Common Stocks*.xlsx`
  - `*Lab*Order*Sheet*.xlsx`
  - run-local uploaded reference workbooks under `runs/<timestamp>/`
- Project- or screen-specific gene-set inputs that imply an active screen,
  e.g.:
  - `data/gene_sets/vGAT-Screen-Maddy-04292026/` and similar
  - any other `data/gene_sets/<active-screen>/` folder, by policy
- Local Maddy_Playground workbooks:
  - generated `RNAi reagents.xlsx`, `other reagents.xlsx`,
    `RNAi Order Form.xlsx`, `RNAi screening plan.xlsx`
  - `StockerOutput.xlsx` snapshots from real screens
- Generated screen-planning outputs from any GUI run:
  - duplicate reagent tracking workbooks
  - follow-up no-stock / no-keyword-hit gene lists
  - `gene_schema-flowchart.png` from real runs
  - similarity workbooks if OpenAI embeddings were used
  - non-phenotypic reagent exports

## What Stays In Git

- Schema documentation (e.g., `docs/allada_lab_sheet_snapshot.md`).
- Tiny synthetic fixture workbooks for parser tests.
- Tiny synthetic gene-list fixtures for end-to-end tests.
- The non-sensitive default split config
  ([`data/config/stock_split_config_no_bloomington.json`](../data/config/stock_split_config_no_bloomington.json)).
- Code paths that point to ignored locations, with explicit instructions on
  where the user should drop their own real workbook locally.

## Operational Pattern For Real Workbooks

When a developer needs a real workbook for manual development:

- Document the local path the GUI expects (e.g.,
  `data/local/Allada Lab Common Stocks.xlsx`).
- Make sure that local path is gitignored.
- Never reference a real workbook from a committed test or example.

## Anti-Goals

- Do not commit real Allada lab sheets.
- Do not commit in-progress screen workbooks.
- Do not commit project-specific gene lists.
- Do not commit generated screen-planning outputs.
- Do not check in test fixtures derived from the real lab inventory.
- Do not cache real workbook contents in source-controllable directories.

## Verification Hooks

Manual checks the maintainer should run before pushing:

- `git status` shows no `.xlsx` from `runs/`, `Maddy_Playground/`, or
  `data/gene_sets/<screen>/` staged.
- `git diff` shows no real stock IDs from the Allada workbook.
- The committed parser tests use synthetic fixtures only.
- `gene_schema-flowchart.png` is not committed.
