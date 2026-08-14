
# Priority Stock-Splitting Configuration Reference

**Visual overview:** [stock_split_config_priority_example.html](stock_split_config_priority_example.html)

**File under review:** `data/config/stock_split_config_priority_example.json`

**Purpose:** simple config for finding phenotypic stocks. It defines six named,
priority-ordered phenotype buckets. For exhaustive RNAi → UAS → allele coverage,
use `data/config/stock_split_config_priority_all_phenotypes.json`.

**Related:** [stock_split_config_example.md](stock_split_config_example.md) (shared
Ref-tier filter glossary) · [pipeline_usage.md § Finding phenotypic stocks](pipeline_usage.md#finding-phenotypic-stocks)

---

## Priority rule

Combinations are evaluated top to bottom. Each `(stock_id, collection)` reagent
lands in only the **first** matching sheet.

| Tier | Sheets | What it captures |
|---|---|---|
| Phenotype buckets | 1–6 | Any curated FlyBase phenotype, by collection/reagent class |

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

## Settings highlights

| Setting | Value |
|---|---|
| `embeddings.enabled` | `true` |
| `validation.enabled` | `false` |
| `input.skipFbgnidConversion` | `true` |
| `maxStocksPerGene` / `maxStocksPerAllele` | `5` / `3` |

The exhaustive `stock_split_config_priority_all_phenotypes.json` variant
sets these limits effectively unlimited, disables functional validation, keeps
embeddings enabled, and adds a final `Pheno·Audit` coverage bucket. Visible tab
numbers are assigned sequentially only to populated categories.

---

## Workflow

```bash
python -m fl_ai_reagent_stocker run ./gene_lists \
  --config data/config/stock_split_config_priority_example.json

python scripts/refactor_aggregated_workbook_layout.py \
  ./gene_lists/Per\ Gene\ Set\ Runs/<gene_set>/Stocks/<gene_set>_aggregated.xlsx
```
