# 01. Vision and Product Priorities

## Goal

Build a local, no-code GUI for fly geneticists to upload or paste a gene list,
enter phenotype keywords, apply simple reagent preferences/limits, and review
phenotype-backed stock candidates with clear gene/stock counts at every step.
The GUI then runs a defined stock-planning stage and emits the same practical
RNAi outputs currently produced by
[`Maddy_Playground/RNAi plan script TEST/`](../Maddy_Playground/RNAi%20plan%20script%20TEST),
without requiring PowerShell, OpenAI, or manual Excel handling.

## Primary User

Non-coding research scientists, especially fly geneticists running RNAi or
related screens. Concrete reference user: Maddy (neuroscience PhD student) who
currently runs Windows + Excel + PowerShell to produce her screening plan.

## Primary Endpoint

**Reagents whose phenotypes hit the user-selected keywords**, sourced from the
soft-run `All Phenotypic Stocks Sheet`.

That is the canonical thing the GUI is optimizing for. Everything else either
feeds it (gene list, keywords, priorities) or branches off it (no-hit
aggregation, fallback stocks, non-phenotypic explorer, RNAi outputs).

## Product Priority Order

1. **Primary**: keyword-hit phenotype-backed reagents from the
   `All Phenotypic Stocks Sheet`.
2. **Secondary aggregation material**: phenotype-backed reagents that do not
   hit the selected keywords.
3. **Tertiary**: input genes with no stocks, no phenotype-backed stocks, or no
   keyword-hit phenotype rows. Surfaced explicitly as a "needs follow-up"
   panel, not buried in a workbook.
4. **Interaction layer**: user-defined priorities (RNAi, UAS, mutant, balancer
   status, stock source, custom phenotype reagent). These sort, group, filter,
   and limit results; they do not reorganize the primary endpoint.
5. **Final workflow endpoint**: stock planning outputs equivalent to Maddy's
   scripts, especially `RNAi Order Form.xlsx` and `RNAi screening plan.xlsx`.
6. **Optional starred enhancement**: OpenAI embedding similarity and AI
   validation. Not required for the main GUI path.

## Defaults That Encode This Priority

- Default config: [`data/config/stock_split_config_no_bloomington.json`](../data/config/stock_split_config_no_bloomington.json).
- Default phenotype keywords selected for the user: `sleep`, `rhythm`,
  `circadian`, `locomotor`.
- Default run mode: `soft_run=True`, OpenAI embeddings off,
  keyword hit/no-hit grouping on.
- Default result tab: keyword-hit phenotypic reagents.
- Default planning input: keyword-hit reagents only.

## What This Vision Explicitly Excludes (MVP)

- No raw JSON editing for end users.
- No CLI flags surfaced to end users.
- No Google sign-in or Google Drive/Sheets API integration.
- No requirement for an OpenAI API key.
- No hosted deployment, multi-user auth, or job queues.
- Mixing non-phenotypic reagents into the phenotype-first browser, stock
  planning, RNAi planning, or any phenotype-derived export. Non-phenotypic
  reagents have a separate, opt-in explorer with a separate export.
