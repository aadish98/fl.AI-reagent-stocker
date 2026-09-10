# AGENTS.md

Repository-wide instructions for coding agents. The full reference — commands,
package layout, key classes, helper scripts, tests, data layout, and the FBgn
conversion hard rule — lives in [`docs/AGENTS.md`](docs/AGENTS.md). Read it
before changing pipeline code.

## Downstream requests on stocker output

This section governs any request to work on reagents, stocks, or sheets that a
stocker run has already produced: filtering, selecting, ranking, capping,
deduplicating, splitting, counting, joining, reformatting, or exporting. It
applies whether the request arrives as "drop the RNAi lines", "keep the top two
stocks per gene", "give me only the ones with a published phenotype", or
anything similar.

Do not answer these requests with ad-hoc terminal commands, throwaway snippets,
notebook cells, one-off pasted code, or hand edits to a workbook. Every such
request ships three things, together, in the same change.

### 1. A script in `scripts/`

- One script per feature, run from the repo root, following the existing
  `scripts/*.py` shape: a module docstring, an `argparse` parser, and named
  functions for each step.
- Read the stocker output workbook or `data/` file and write a **new** file.
  Never overwrite, re-save, or edit the run output in place. The one exception
  is a pass that changes only presentation — column widths, freezing, header
  formatting — and no cell value.
- Reuse the shared helpers in `fl_ai_reagent_stocker/utils.py` rather than
  re-implementing FlyBase, PubMed, or GAL4 logic locally.
- Keep the workbook's own column names on the output. Do not rename, drop, or
  reorder columns unless the request asks for it.
- If the rule depends on a threshold, a keyword list, or a cap, expose it as a
  command-line option with a documented default instead of burying the number
  in the code.

### 2. A test in `tests/`

- Name it `tests/test_<script name>.py` and load the script the way
  `tests/test_summarize_combination_counts.py` does.
- **Derive fixtures from the real files.** Copy actual rows out of the workbook
  or `data/` file the request is about, trim them to the few rows that exercise
  the rule, and save the extract under `tests/fixtures/`. Never invent column
  names, ID formats, or values; a fixture that does not match the real file
  proves nothing.
- Fixtures must stay small and must not carry unpublished or private data.
- Cover each keep/remove/select step: one case that the step keeps, one that it
  removes, and the boundary case if there is a threshold or a cap.
- `python -m pytest tests` must pass before the change is done.

### 3. An entry in `docs/features.md`

- One section per feature. Record what goes in, every step in the order it
  runs, what comes out, and the exact column names each step reads.
- Update the entry in the same change that changes the script. A script whose
  documented steps no longer match its code is a bug.

### Wording of comments and documentation

Comments above each filtering or selection step, and the steps written in
`docs/features.md`, use the vocabulary of a spreadsheet:

- Name columns exactly as they appear in the workbook: `Stock Center`,
  `flybase_gene_id`, `ALLELE_PAPER_RELEVANCE_SCORE`, `RNAi`, `genotype`.
- Describe the action with a plain verb: keep, remove, sort, group, count,
  match, combine, deduplicate, cap.
- State the rule and the reason for it, not the mechanics of the code.

```python
# Keep only rows where Stock Center is BDSC.
rows = rows[rows["Stock Center"] == "BDSC"]

# Sort by PHENOTYPE_RELEVANCE_SCORE, highest first, then keep at most 3 rows
# for each flybase_gene_id so no single gene fills the order.
rows = rows.sort_values("PHENOTYPE_RELEVANCE_SCORE", ascending=False)
rows = rows.groupby("flybase_gene_id").head(3)
```

Avoid implementation jargon in these comments and in `docs/features.md`:
"boolean mask", "vectorized", "predicate", "left join on the key", "regex", "O(n)",
"dataframe", "dedupe via hash", "short-circuit". Ordinary code that is not a
filtering or selection step needs no such comment; do not narrate obvious lines.

Never write down who the comments or the documentation are for, and never
describe a reader's background, role, or level of technical knowledge — not in
the code, not in `docs/features.md`, not in commit messages.
