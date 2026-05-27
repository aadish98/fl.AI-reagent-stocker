# 03. Config and Keywords

## Default Split Config

Default config the GUI loads on first run:
[`data/config/stock_split_config_no_bloomington.json`](../data/config/stock_split_config_no_bloomington.json).

Important naming note: despite its filename, this config does not exclude
Bloomington. Its description says it has "no Bloomington filter in
combinations," meaning all stock centers are included rather than restricting
to Bloomington-only. The user-facing label in the GUI should be something
like **"All stock centers"**, not "no Bloomington."

The config is treated as a **template**, not as the user's editable file:

- On run start, copy it into the per-run working directory
  (`runs/<timestamp>/`), then mutate the copy.
- Never mutate the repo template directly from the GUI.

## Default Phenotype Keywords (Pre-Selected)

`sleep`, `rhythm`, `circadian`, `locomotor`.

These match the screening domain Maddy currently runs and the targets already
in [`data/config/stock_split_config_no_bloomington.json`](../data/config/stock_split_config_no_bloomington.json).

Users can add, remove, or replace keywords. The selected keyword set drives:

- the keyword-hit / no-hit partition in the result browser
- per-keyword `<keyword> phenotype` columns in the planning workbooks
- the `Phenotype keyword hits` summary count
- the planning-stage default candidate set
- the labels in `gene_schema-flowchart.png`

## Plain-Language Priority Controls (No JSON)

The user-facing fields exposed on the keyword/priority screen are intentionally
narrow:

- `settings.relevantSearchTerms` (the keyword input)
- simple reagent-category priorities surfaced as labeled cards:
  - "Prefer RNAi first"
  - "Prefer no balancers"
  - "Group custom phenotype reagents separately"
  - "Limit stocks per gene"
- `settings.maxStocksPerGene`
- `settings.maxStocksPerAllele`
- `settings.phenotypeSimilarityTargets` — starred, only useful with OpenAI
  embeddings, hidden by default

Each card shows a one-sentence explanation of what it changes (sort vs filter,
which output buckets it affects).

The raw `filters` and `combinations` blocks of the JSON are not editable in
the MVP. Per-run sorting/filter behavior is generated from the simple cards.

Before each run, show a human-readable summary of effective settings and a
live count preview where possible.

## What Is Hidden Or Locked

- Raw `filters` / `combinations` JSON.
- CLI flag names (`--soft-run`, `--keyword-bucketing`, etc.).
- OpenAI embedding similarity targets, unless an OpenAI API key is configured
  and the user opts in.
- `phenotypeSimilarityTargets` editing UI.
- The `settings.relevantSearchTerms` / `settings.phenotypeSimilarityTargets`
  semantic distinction; the GUI surfaces a single "phenotype keywords" concept
  and treats similarity targets as a starred premium variant.

## Why This Matters

The user must consistently know what is going on:

- The keyword set chosen at the top of the run is the same set referenced in
  results, planning, RNAi outputs, and the schema flowchart.
- No silent transformations.
- No editing JSON.
- No per-screen concept inversion.
