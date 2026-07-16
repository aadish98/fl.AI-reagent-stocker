
# Boss-Priority Stock-Splitting Configuration Reference

**Visual overview:** [stock_split_config_boss_priority_example.html](stock_split_config_boss_priority_example.html)

**File under review:** `data/config/stock_split_config_boss_priority_example.json`

**Purpose:** standard config for finding phenotypic stocks. Six phenotype-priority
buckets run first; remaining Bloomington `RNAi_reagent` stocks fall through to a
twelve-sheet Ref matrix. **18 output sheets** total.

**Related:** [stock_split_config_example.md](stock_split_config_example.md) (shared
Ref-tier filter glossary) · [pipeline_usage.md § Finding phenotypic stocks](pipeline_usage.md#finding-phenotypic-stocks)

---

## Priority rule

Combinations are evaluated top to bottom. Each `(stock_id, collection)` reagent
lands in only the **first** matching sheet.

| Tier | Sheets | What it captures |
|---|---|---|
| Phenotype buckets | 1–6 | Any curated FlyBase phenotype, by collection/reagent class |
| Bloomington Ref matrix | 7–18 | Remaining Bloomington `RNAi_reagent` stocks (non-CRISPR), tiered by balancer, insertion count, and Ref evidence |

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

## Tier 2 — Bloomington Ref matrix (sheets 7–18)

Base: `Bloomington` + `RNAi_reagent` + `No sgRNA`.

| Insertion | Balancers | Ref++ | Ref+ | Ref− |
|---|---|---|---|---|
| Single | No | 7 | 9 | 11 |
| Single | Has | 8 | 10 | 12 |
| Multiple | No | 13 | 15 | 17 |
| Multiple | Has | 14 | 16 | 18 |

Ref tiers come from `ALLELE_PAPER_RELEVANCE_SCORE` and keywords
`sleep`, `circadian`, `rhythm`, `locomotor`.

---

## Settings highlights

| Setting | Value |
|---|---|
| `embeddings.enabled` | `true` |
| `validation.enabled` | `false` |
| `input.skipFbgnidConversion` | `true` |
| `maxStocksPerGene` / `maxStocksPerAllele` | effectively unlimited |

---

## Workflow

```bash
python -m fl_ai_reagent_stocker run ./gene_lists \
  --config data/config/stock_split_config_boss_priority_example.json

python scripts/refactor_aggregated_workbook_layout.py \
  ./gene_lists/Per\ Gene\ Set\ Runs/<gene_set>/Stocks/<gene_set>_aggregated.xlsx
```
