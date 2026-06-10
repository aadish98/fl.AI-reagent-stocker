# 06. RNAi Outputs

## Source Of Truth

Maddy's PowerShell scripts under
[`Maddy_Playground/RNAi plan script TEST/`](../Maddy_Playground/RNAi%20plan%20script%20TEST):

- [`Build-RNAi-Split-And-Order.ps1`](../Maddy_Playground/RNAi%20plan%20script%20TEST/Build-RNAi-Split-And-Order.ps1)
- [`Build-RNAi-Screening-Plan-Step1.ps1`](../Maddy_Playground/RNAi%20plan%20script%20TEST/Build-RNAi-Screening-Plan-Step1.ps1)

These are treated as a **behavioral spec**, not as code to preserve. The GUI's
RNAi planning module must produce the same final outputs but with named
columns, no fixed Excel positions, no PowerShell, no Excel COM automation, and
explicit auditability.

The Python landing point will be a new module, e.g.
`fl_ai_reagent_stocker/pipelines/rnai_planning.py`.

## Final User-Facing Documents

Two primary final workbooks plus a visual audit:

1. **`RNAi screening plan.xlsx`** — the enriched screening workbook the user
   actually runs the screen from.
2. **`RNAi Order Form.xlsx`** — the orderable form for known stock centers
   (Bloomington / Vienna / NIG-Fly), sorted by center.
3. **`gene_schema-flowchart.png`** — emitted beside both workbooks at the end
   of every run; see [`08-schema-flowchart.md`](08-schema-flowchart.md).

Audit / supporting outputs (downloadable, not the deliverable):

- `RNAi reagents.xlsx` (split + deduped intermediate)
- `other reagents.xlsx` (non-RNAi rows from the same input)
- Duplicate reagent tracking workbook (every dropped source row plus reason)
- Follow-up gene list (no stocks / no keyword hits)

## Pipeline Steps Inside The Module

Following the spec from Maddy's scripts but with named columns and explicit
counts:

1. **Select source rows** from `All Phenotypic Stocks Sheet`.
   - Default: keyword-hit RNAi rows.
   - Option: include no-keyword-hit RNAi rows.
2. **Classify reagents**.
   - RNAi vs non-RNAi (today driven by Maddy's column 42 / `RNAi_reagent`
     marker; the Python port must use a named column).
   - Known stock centers vs custom reagents.
   - One-hot reagent buckets (RNAi / UAS / mutant / balancer status).
3. **Deduplicate with phenotype folding**.
   - Dedupe by stock ID and effective RNAi shorthand.
   - Fold phenotype, references, and keyword-hit context into the kept row.
   - Save discarded rows to the duplicate audit workbook with reasons and
     kept-row pointers.
4. **Lab snapshot lookup** (see [`07-allada-snapshot.md`](07-allada-snapshot.md)).
   - Gate the screen-planning flow on an uploaded Allada `.xlsx` snapshot.
   - Mark already-in-lab rows with sheet/tab/location context.
5. **Reference enrichment** (Bloomington genotype/comments, Vienna/VDRC
   viability/chromosome/library, TRiP vector, RNAi control, validation flags).
6. **Generate outputs**, each with a summary/README sheet explaining what
   happened.

## Column Retention In `RNAi screening plan.xlsx`

Maddy's screening plan keeps 44 columns. Columns 1-41 come from the source
phenotype sheet via this hardcoded source-column list:

```text
@(1..11) + @(12, 43, 13, 14, 15) + @(12, 16, 17, 18, 19, 20, 21, 22, 23, 24)
+ @(26..32) + @(41, 42, 44, 45, 47, 48, 49, 50)
```

Source columns 47-50 are `Cosine Similarity (sleep / circadian / rhythm /
locomotor)`, summed into output column 42 (`SumCosine`). Output columns 43-44
are computed `ValidationFlag` and `ValidationNote`.

The Python port must:

- Map by **named columns**, not Excel column letters or numeric indices.
- Replace `Cosine Similarity (<keyword>)` with interpretable
  **`<keyword> phenotype` columns** in the default non-OpenAI path. Each
  column indicates whether the row's phenotype/reference text matched that
  keyword; these columns are used for sorting, filtering, counts, and
  planning exports.
- Include `Cosine Similarity (<keyword>)` and a summed score only as
  starred supplemental columns when OpenAI embeddings were explicitly run.
- Fail loudly if a required source column is missing.
- Provide compatibility aliases for any current workbook column names so
  drop-in replacement is feasible.

## Naming Decisions

- Final workbook column name pattern in the default path:
  `<keyword> phenotype` (e.g., `sleep phenotype`, `rhythm phenotype`,
  `circadian phenotype`, `locomotor phenotype`).
- Final workbook label for the cosine sum: not present unless OpenAI embeddings
  were run; if present, it stays `SumCosine` for compatibility.
- Final workbook columns `ValidationFlag` and `ValidationNote` are kept as-is
  and populated by the same RNAi-control assignment rules used in
  `Set-RnaiControlColumnAc`.

## What This Module Must Not Do

- No reliance on Excel COM automation.
- No reliance on Windows or PowerShell.
- No reliance on column positions.
- No silent duplicate dropping.
- No silent re-sorting hidden from the user.
- No writing back to the uploaded Allada snapshot.
- No mixing non-phenotypic reagents into RNAi outputs.
