
# Priority Stock-Splitting Configuration Reference

**Visual overview:** [stock_split_config_priority_example.html](stock_split_config_priority_example.html)

**File under review:** `data/config/stock_split_config_priority_example.json`

**Purpose:** public default config for phenotype-backed reagents plus
Bloomington RNAi literature follow-up. It defines eight named, priority-ordered
output buckets: six phenotype buckets followed by `Ref++` and `Ref+` RNAi
buckets. Its `filters` object is also the reusable ingredient catalog for
user-authored combinations.

**Related:** [stock_split_config_example.md](stock_split_config_example.md) (shared
Ref-tier filter glossary) · [pipeline_usage.md § Finding phenotypic stocks](pipeline_usage.md#finding-phenotypic-stocks)

---

## Priority rule

Combinations are evaluated top to bottom. Each `(stock_id, collection)` reagent
lands in only the **first** matching sheet.

| Tier | Sheets | What it captures |
|---|---|---|
| Phenotype buckets | 1–6 | Any curated FlyBase phenotype, by collection/reagent class |
| Bloomington RNAi literature buckets | 7–8 | Remaining Bloomington RNAi reagents with `Ref++` or `Ref+` publication evidence |

---

## Tier 1 — Phenotype buckets (sheets 1–6)

All entries require `Phenotype` (`PHENOTYPE_RELEVANCE_SCORE ≠ 0`).

| # | Combination | Intent |
|---|---|---|
| 1 | Bloomington + RNAi_reagent + Phenotype | Bloomington RNAi with phenotype |
| 2 | Vienna + RNAi_reagent + Phenotype | Vienna RNAi with phenotype |
| 3 | Bloomington + AlleleOrInsertion + Phenotype | Bloomington allele/insertion with phenotype |
| 4 | Stock Center + AlleleOrInsertion + Phenotype | Orderable non-RNAi with phenotype |
| 5 | Non Stock Center + RNAi_reagent + Phenotype | Custom RNAi reagent with phenotype |
| 6 | Non Stock Center + AlleleOrInsertion + Phenotype | Custom allele/insertion with phenotype |

---

## Tier 2 — Bloomington RNAi literature buckets (sheets 7–8)

These sheets are evaluated only after all phenotype buckets. `Ref++` requires
at least one configured search term in a cited paper's title or abstract;
`Ref+` has citations but no such title/abstract hit.

| # | Combination | Intent |
|---|---|---|
| 7 | Bloomington + RNAi_reagent + Ref++ | Literature-supported Bloomington RNAi not already claimed by a phenotype bucket |
| 8 | Bloomington + RNAi_reagent + Ref+ | Bloomington RNAi with citations but no configured keyword hit, not already claimed earlier |

The default 5-stocks-per-gene and 3-stocks-per-allele caps apply across every
sheet in priority order. A phenotype bucket can therefore exhaust a gene or
allele cap before either Ref bucket is considered.

---

## Settings highlights

| Setting | Value |
|---|---|
| `embeddings.enabled` | `true` |
| `validation.enabled` | `false` |
| `input.skipFbgnidConversion` | `true` |
| `maxStocksPerGene` / `maxStocksPerAllele` | `5` / `3` |

Visible tab numbers are assigned sequentially only to populated categories.

---

## Workflow

```bash
python -m fl_ai_reagent_stocker run ./gene_lists \
  --config data/config/stock_split_config_priority_example.json

python scripts/refactor_aggregated_workbook_layout.py \
  ./gene_lists/Per\ Gene\ Set\ Runs/<gene_set>/Stocks/<gene_set>_aggregated.xlsx
```
