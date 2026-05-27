# 08. Schema Flowchart Automation

## Required Artifact

Every completed run must emit, alongside the two primary final workbooks:

`gene_schema-flowchart.png`

Together with:

- `RNAi screening plan.xlsx`
- `RNAi Order Form.xlsx`

These three files are the audit-ready bundle a scientist needs after a run.

## Why It Exists

Each screen, planning step, and RNAi output should be inspectable without
opening logs. The flowchart shows **what actually happened** in the run, not
what could have happened in general.

## What It Should Show

At minimum:

- Input gene count
- Selected phenotype keywords
- Phenotype keyword-hit count vs no-hit count
- Genes with no stocks / no phenotype-backed stocks / no keyword-hit
  phenotype rows
- Reagent priority settings the user applied
- Duplicate folding counts (groups, kept reagents, folded phenotypes)
- Allada lab-sheet match count
- Final order-form row count
- Final screening-plan row count
- Whether OpenAI embeddings were run (starred when used)
- Whether GPT validation was run (starred when used)
- Allada snapshot status (uploaded filename, sheets indexed)
- Whether lab-stock lookup was active

## Honest Outputs

The flowchart must reflect actual choices, not theoretical ones:

- "Fallback stock search: skipped" if the user did not run it.
- "Lab stock lookup: not run" if no Allada snapshot was uploaded.
- "★ OpenAI similarity: locked / not run" if OpenAI was unavailable.
- "Non-phenotypic explorer: not opened" if the user did not toggle it.

This keeps the artifact useful for troubleshooting.

## Architecture Sketch

The implementation should not recompute counts ad hoc. Each pipeline step
emits a **structured checkpoint**:

```python
{
  "step": "keyword_hit_partition",
  "label": "Phenotype keyword split",
  "inputs":  {"phenotype_rows": 1240},
  "outputs": {
    "keyword_hit_rows": 312,
    "no_keyword_hit_rows": 928,
    "genes_with_keyword_hits": 42
  },
  "settings": {
    "keywords": ["sleep", "rhythm", "circadian", "locomotor"]
  }
}
```

These checkpoints accumulate into a `RunSummary` (or `GeneSchemaSummary`)
object, then the run emits:

- `gene_schema-summary.json` — machine-readable audit/debug.
- `gene_schema-flowchart.png` — the scientist-friendly visual artifact.
- Optional `Run Summary` sheet inside the output workbooks.

Suggested module layout:

```text
fl_ai_reagent_stocker/reporting/run_summary.py
fl_ai_reagent_stocker/reporting/flowchart.py
```

`flowchart.py` consumes the summary object and writes the PNG. Use Graphviz
or `matplotlib`/`networkx` for reliable local PNG generation. Avoid Mermaid
PNG rendering, which typically requires a browser/Node toolchain.

## Display Order Inside The PNG

Roughly mirrors the user flow:

```text
Input genes
  -> Mapped genes
  -> Genes with phenotype-backed reagents
  -> Keyword-hit reagents | No-keyword-hit reagents
  -> RNAi planning candidates
  -> Duplicate groups -> Kept reagents (after fold)
  -> Already-in-lab matches
  -> Order-form rows | Screening-plan rows
```

## Failure Mode

If a checkpoint is missing or inconsistent, the GUI must:

- Mark it explicitly as "missing" rather than render a misleading number.
- Surface a warning before the user trusts the schema artifact.
