# Features

Every feature that works on the reagents a stocker run has already produced is
listed here, in the order its steps run. Each entry gives the file it reads,
the columns it looks at, the rule for keeping or removing a row, and the file it
writes. A feature that selects or filters rows always writes a new file and
leaves the run output untouched; the only features that write back to a workbook
are the ones that change nothing but its layout.

Add an entry whenever a new script is added under `scripts/`, and update the
entry in the same change that changes the script's steps.

---

## Combination counts summary

**Script:** `scripts/summarize_combination_counts.py`

**What it reads:** every organized stocker workbook (`*_aggregated.xlsx`) found
in a folder and its subfolders, for example a `Per Gene Set Runs` folder.

**What it does:**

1. Find every `*_aggregated.xlsx` workbook in the folder and its subfolders.
2. In each workbook, open the `Contents` sheet and read the
   `Prioritized sheet breakdown` table.
3. From that table, take the `Sheet criteria`, `# Stocks`, `# Alleles`,
   `# Genes`, and `Sheet Name` columns.
4. Name the gene set after the folder the workbook came from, and record the
   folder path so two gene sets with the same name stay separate.
5. Put every row from every workbook into one table, one row per gene set and
   per sheet.
6. When a JSON config is supplied with `--config`, list the sheets in the order
   the config gives them, and fill in a count of `0` for any sheet a workbook
   did not produce.

**What it writes:** `combination_counts_summary.xlsx` and a matching `.csv` in
the folder that was scanned, with the columns `relative_path`, `gene_set`,
`combination`, `num_genes`, `num_stocks`, `num_alleles`, `sheet_name`, and
`workbook_path`. The `run` command writes this file automatically at the end of
a run.

**Test:** `tests/test_summarize_combination_counts.py`

---

## Aggregated workbook readability pass

**Script:** `scripts/refactor_aggregated_workbook_layout.py`

**What it reads:** one organized stocker workbook (`*_aggregated.xlsx`).

**What it does:**

1. Widen the columns that hold long text, such as `genotype`, `phenotype`,
   `title`, `author`, `journal`, `notes`, and `reference` columns.
2. Set a readable width for the remaining columns.
3. Freeze the header row so it stays visible while scrolling.
4. Apply consistent header formatting and text wrapping.

No cell value is added, removed, or edited; only the way the sheet is displayed
changes.

**What it writes:** the same workbook, with its layout updated.

---

## Adding a feature

Copy the shape of the entries above:

```markdown
## <Name of the feature>

**Script:** `scripts/<file>.py`

**What it reads:** <the workbook or data file, and the sheet>

**What it does:**

1. <One step per line, in the order it runs, naming the columns it uses and
   whether it keeps, removes, sorts, groups, counts, matches, combines,
   deduplicates, or caps.>

**What it writes:** <the new file, and the columns on it>

**Test:** `tests/test_<file>.py`
```

Steps are written the way the same operation would be described in a
spreadsheet: name the column exactly as it appears in the workbook, and say
whether the rows are kept, removed, sorted, grouped, counted, matched,
combined, deduplicated, or capped.
