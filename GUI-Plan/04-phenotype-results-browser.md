# 04. Phenotype-First Results Browser

The result browser is the GUI's primary surface. Everything else (planning,
RNAi outputs, downloads) is reached from here.

## Tab Order (Strict)

1. **Keyword Hits** — phenotype rows whose `Phenotype` text matches any of the
   user-selected keywords. This is the primary endpoint.
2. **No Keyword Hits** — phenotype rows that do not match any keyword. Optional
   aggregation/review material.
3. **Input Genes Needing Fallback** — input genes with no stocks, no
   phenotype-backed stock rows, or no keyword-hit phenotype rows. First-class
   panel, not buried in a workbook.
4. **Fallback Stock Candidates** — opt-in. Triggered explicitly per-gene or
   for the whole no-hit set; never the default.

## Curated Columns

Display order, scientist-first:

- Gene
- Phenotype
- Matched keyword terms
- `<keyword> phenotype` columns (one per selected keyword) — the non-OpenAI
  planning signal
- Reagent type (one-hot reagent bucket)
- Stock ID, Stock Source
- Genotype
- Balancers
- RNAi / UAS / mutant / sgRNA flags
- PMID / Reference
- Optional starred `Cosine Similarity (<keyword>)` and `Max Cosine Similarity`
  if OpenAI embeddings were run

## Filter / Sort Rules

User-defined reagent priorities (RNAi, UAS, mutant, no balancers, stock
source, custom phenotype reagent) act as **interaction layers** on top of the
phenotype-first views:

- They **sort** and **group**, they don't replace the phenotype-first
  organization.
- Filters that hide rows are reversible and visibly indicated.
- A reagent category filter must never silently remove keyword-hit phenotypic
  reagents unless the user explicitly toggled it.
- Clearing filters returns the browser to the keyword-hit-first default
  immediately.

## "Genes Needing Follow-Up" Panel

Surfaced as a dedicated panel, not just a workbook download. It includes:

- Input genes with no stocks at all.
- Input genes with stocks but no phenotype-backed stock rows.
- Input genes with phenotype rows but no keyword hits.

For each gene, the panel shows the reason and offers actions:

- "Show fallback stock candidates" — opens the fallback view for that gene.
- "Save as follow-up gene list" — exports a reusable CSV for a subsequent run.

## Summary Counts

A summary tile at the top of the browser shows:

- Input genes
- Genes with keyword-hit phenotype reagents
- Genes with phenotype rows but no keyword hit
- Genes with no stocks
- RNAi / UAS / mutant counts
- Stock count
- Reference count
- Number of output buckets

These are the same counts that feed `gene_schema-flowchart.png`.

## Optional Secondary View: Non-Phenotypic Reagent Explorer

This is **strictly separate** from the phenotype-first browser, hidden behind
an opt-in toggle, and excluded from any phenotype-derived export.

- Source: the complete Stage 1 `Gene Reagent Index` (persisted as part of the
  phenotype-sheet rework — see
  [`/Users/aadishms/.cursor/plans/phenotype_sheet_rework_844d6a5c.plan.md`](/Users/aadishms/.cursor/plans/phenotype_sheet_rework_844d6a5c.plan.md)).
- Surfaces input-gene FBal/FBtp/FBti reagents that have no row in
  `genotype_phenotype_data`.
- Columns: input gene symbol, FBgn, component ID, component type, reagent
  symbol, allele class, transgenic product class, RNAi type, plus Stage 1
  stock metadata when available.
- Filters: gene multi-select, component type, stock availability, source
  collection, RNAi type, balancer presence, multiple insertions, allele class,
  transgenic product class.
- Counts panel: total reagents in the index, with-phenotype, without-phenotype,
  filtered-in.
- Export target: a separate `Non-Phenotypic Reagents.xlsx` workbook or a
  labeled sheet in the run workbook. **Never inlined into**
  `Stock Phenotype Sheet`.
- Default state: collapsed/hidden so the primary phenotype-first journey stays
  uncluttered.

## Anti-Goals

- No mixing of non-phenotypic reagents into keyword-hit / no-keyword-hit views.
- No raw JSON config editor in this screen.
- No silent dropping of keyword-hit reagents because of an applied category
  filter.
- No "advanced" filters that change the meaning of the keyword sets.
