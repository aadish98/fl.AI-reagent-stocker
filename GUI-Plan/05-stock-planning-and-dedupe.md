# 05. Stock Planning and Dedupe

## Where Planning Sits

Planning is an **explicit stage between browsing and final downloads**. It is
not implicit in browsing, and it is not an opaque step inside RNAi output
generation.

Default planning input: keyword-hit phenotype reagents only. Users can opt
into adding no-keyword-hit phenotype reagents as secondary material.

## Planning Controls (Plain Language)

All controls are mapped from the priority cards on the keyword/priority
screen — never from raw JSON:

- Reagent type priority (RNAi, UAS, mutant, balancer status, stock source,
  custom phenotype reagent).
- Stock source priority (e.g., Bloomington vs Vienna vs NIG-Fly).
- Balancer preference (no balancer first, allow balancer, etc.).
- Duplicate handling (fold + keep one row vs keep all).
- Per-gene and per-reagent limits.
- Already-in-lab handling (include / flag / exclude from order exports).

## Dedupe with Phenotype Folding

Maddy's PowerShell scripts deduplicate by stock ID and RNAi shorthand and
**drop entire rows**, losing whatever phenotype lived on the discarded rows.
The GUI must improve on this:

- Deduplicate to one **reagent record**.
- **Fold** all phenotype context from duplicates into the kept reagent record,
  preserving:
  - all phenotype terms
  - all PMIDs / references
  - all matched keyword terms
  - best/highest similarity score (if OpenAI embeddings were run)
  - the full list of duplicate source rows (for audit)
- Persist a separate **duplicate-tracking workbook** containing every dropped
  source row plus dedupe reason and a pointer to the kept reagent record.

## Visible Count Path Through Planning

Surface counts at every step so nothing happens silently:

- Candidate genes
- Candidate reagents
- Duplicate groups
- Kept reagents (after fold)
- Folded duplicate phenotypes
- Already-in-lab hits
- Final order-form rows
- Final screening-plan rows

These counts feed `gene_schema-flowchart.png` (see
[`08-schema-flowchart.md`](08-schema-flowchart.md)).

## Persisted Outputs From Planning

- The planned reagent set (orderable / screenable).
- The duplicate-tracking workbook (audit).
- The follow-up gene list (genes with no stocks or no keyword hits).
- The schema flowchart PNG.
- The non-phenotypic reagent export (only if the secondary explorer ran).

## Anti-Goals

- No silent dedupe that drops phenotype context.
- No hidden reordering by "advanced" criteria the user didn't set.
- No mixing already-in-lab reagents into the order form by default — they
  should be flagged separately, with the user choosing whether to include them.
- No relying on column positions (column 42, J, K, AC, etc.) the way Maddy's
  scripts do today. The Python port must use named columns end-to-end.
