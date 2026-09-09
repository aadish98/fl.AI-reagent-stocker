# LEARNINGS.md

Living notes for coding agents working in this repository.

## How to use this file

- Read this file at the start of every agent session in this repo.
- Append a dated entry after non-trivial pipeline, output, or documentation changes.
- Keep entries factual and short: context → decision → where in code.

## Template

```markdown
### YYYY-MM-DD — Short title
- Context: …
- Decision: …
- Code: `path/to/file.py` (`function_name`)
```

## Seed entries

### 2026-06-25 — Sheet1..N phenotype enrichment columns
- Context: Sheet1..N previously exposed aggregated `Phenotype` text and stock-level PMIDs separately, without empirical phenotype↔reference pairing or published GAL4 partner context.
- Decision: When the JSON config uses any Phenotype filter, build the internal phenotype sheet (even without embeddings) and merge five Sheet1..N columns per focal reagent: `Phenotype Reference-ID Map`, `Partner GAL4 FlyBase IDs`, `Partner GAL4 Symbols`, `Published GAL4/ Positive Control`, and `Published GAL4 stock id`.
- Code: `fl_ai_reagent_stocker/pipelines/stock_splitting.py` (`_build_reagent_phenotype_enrichment`, `_add_phenotype_enrichment_to_stock_sheet`)

### 2026-06-25 — Paired GAL4 enrichment entries
- Context: Researchers need both the published GAL4 symbol and the orderable partner stock number with the same collection, reference ID, and phenotype context.
- Decision: Emit parallel comma-separated `{…}` entries. Positive Control field 1 = published GAL4 symbol; field 2 = orderable driver genotype on the recommended gal4-only stock (e.g. `P{GawB}elav[C155]` for Bloomington 458); stock id field 1 = stock-center number; fields 3–5 (collection, reference ID, phenotype) must match between the two columns at each index. Use `-` for genotype/collection when no gal4-only stock resolves.
- Code: `fl_ai_reagent_stocker/pipelines/stock_splitting.py` (`_build_reagent_published_gal4_enrichment`, `_build_stock_candidate_label_to_driver_genotype`, `PARTNER_DRIVER_STOCK_GENOTYPES_COLUMN`)

### 2026-06-25 — Prefixed PMID / PMCID export
- Context: Bare numeric PMIDs/PMCIDs are hard to validate and copy into PubMed workflows.
- Decision: Format display-only reference IDs as `PMID{n}` / `PMCID{n}` at Excel write time; keep internal joins on `clean_id()`.
- Code: `fl_ai_reagent_stocker/utils.py` (`format_pmid_display`, `format_pmcid_display`, `format_reference_id_columns`)

### 2026-06-25 — Phenotype sheet build gate
- Context: Phenotype enrichment on Sheet1..N requires phenotype-sheet data, not only the embedding soft-run.
- Decision: Build the internal phenotype sheet when `soft_run=True` **or** `_config_uses_phenotype_filter(config)`.
- Code: `fl_ai_reagent_stocker/pipelines/stock_splitting.py` (`write_aggregated_excel`, `_config_uses_phenotype_filter`)

### 2026-06-25 — Join key for stock ↔ phenotype enrichment
- Context: Stock rows and phenotype-sheet rows must merge on a stable reagent label.
- Decision: Use `_format_source_stock_label(collection, stock_number)` on stock rows and `Source/ Stock #` on phenotype rows.
- Code: `fl_ai_reagent_stocker/pipelines/stock_splitting.py` (`_stock_sheet_reagent_key`, `_get_source_stock_series`)

### 2026-06-25 — Masterlist vs Sheet1..N GAL4 columns
- Context: The All Phenotypic Stocks Sheet masterlist export already exposes symbol-only `Published Gal4/ Positive control`.
- Decision: Do not change masterlist semantics. Sheet1..N uses separate enriched columns `Published GAL4/ Positive Control` and `Published GAL4 stock id`.
- Code: `fl_ai_reagent_stocker/pipelines/stock_splitting.py` (`MASTERLIST_TEMPLATE_SOURCE_COLUMNS`, `STOCK_SHEET_ENRICHMENT_COLUMNS`)

### 2026-06-25 — Avoid double masterlist reorder
- Context: Re-running `_reorder_to_masterlist_columns()` on already-masterlist sheets copied PMID values into `Column 30`.
- Decision: Skip masterlist reorder when `_is_masterlist_formatted()` is true before writing phenotype workbook sheets.
- Code: `fl_ai_reagent_stocker/pipelines/stock_splitting.py` (`_prepare_phenotype_workbook_output`, `_is_masterlist_formatted`)

### 2026-06-25 — Partner GAL4 columns: token split + FBal→FBst chain (do not match whole genotype)
- Context: `Partner driver symbols (best-effort)` and `Partner driver stock candidates` on the Stock Phenotype Sheet were wrong on real data (e.g. vGAT `*_similarity_tiers.xlsx`). Symptoms: (1) partner symbols were the **entire** `genotype_symbols` string (e.g. `A23[A23] Scer\GAL4[Act5C.PI]`) instead of the GAL4 token alone; (2) stock candidates were **always empty** (0 non-empty rows) even when symbols looked populated; (3) agents tried to “match on the whole genotype” or look up partner `FBal` IDs directly in `fbst_to_derived_stock_component.csv`.
- Root causes:
  - FlyBase `genotype_symbols` are **space-delimited** allele components (`UAS-line Scer\GAL4[driver]`). Tokenizers that only split on `;`, `/`, `|`, or `,` leave the full string as one token, so GAL4 fallback logic treats the whole genotype as a “GAL4 token”.
  - `genotype_phenotype_data` co-reagent IDs are usually **`FBal` alleles** (e.g. `FBal0366188` = `Scer\GAL4[Hand.Delta]`). The derived stock CSV indexes most GAL4 **stocks** by **`FBti` insertion** IDs, not by those `FBal` IDs — direct `FBal` lookup in the derived CSV almost always returns nothing.
  - Partner columns are **GAL4-only** (not LexA, balancers, or other co-reagents). Do not fill them when the focal stock *is* the driver and there are no co-reagent IDs.
- Decision:
  - Split genotype text with `_split_genotype_symbol_tokens`: `;|,\n` → `/` → **whitespace** so each allele is its own token; use `_extract_unmatched_gal4_tokens` only when `co_reagent_ids` is non-empty and the token contains `gal4` but is not part of the focal stock’s matched symbols.
  - Resolve partner stocks via **`FBal → FBtp → FBti → FBst`**: construct descriptions (`transgenic_construct_descriptions_fb_*.tsv`) map allele → construct; `data/flybase/transgenic_insertions/fbtp_to_fbti.csv` maps construct → insertion; `fbst_to_derived_stock_component.csv` maps insertion → stock. Use `_resolve_fbal_stock_candidates(..., gal4_only=True)` when direct component lookup on the `FBal` fails.
  - Emit stock candidates as **`(<stock #>, <collection>)`** via `_format_stock_candidate_label` (not `Source/ Stock # [FBst…]`).
  - Prefer exact genotype-aligned symbols from `_extract_genotype_id_symbol_pairs` over generic FlyBase lookup symbols (e.g. `Scer\GAL4` vs `Scer\GAL4[tim-GAL4]`).
- Verification (after any change to partner logic): on a real phenotype sheet, `Partner driver symbols` should be short GAL4 tokens only; `Partner driver stock candidates` should have **some** non-empty rows when drivers are orderable; spot-check e.g. `Scer\GAL4[Act5C.PI]` → `(28816, Bloomington)`-style labels. Regression tests: `tests/test_phenotype_reagent_capture.py` (`test_stock_phenotype_sheet_recovers_partner_driver_candidates`, `test_stock_phenotype_sheet_recovers_gal4_partner_when_symbol_pairing_is_imperfect`, `test_stock_phenotype_sheet_excludes_non_gal4_partner_reagents`).
- Code: `fl_ai_reagent_stocker/pipelines/stock_splitting.py` (`_split_genotype_symbol_tokens`, `_extract_unmatched_gal4_tokens`, `_looks_like_gal4_symbol`, `_format_stock_candidate_label`, `_resolve_fbal_stock_candidates`, `_build_stock_phenotype_sheet`; columns `PARTNER_DRIVER_SYMBOLS_COLUMN`, `PARTNER_DRIVER_STOCK_CANDIDATES_COLUMN`)

### 2026-06-25 — Sheet1..N Partner GAL4 summary + symbol fallback
- Context: `Published GAL4/ Positive Control` and `Published GAL4 stock id` on Sheet1..N were `-` when genotype detected a partner GAL4 but no gal4-only stock resolved (common for mixed `tim-GAL4 / UAS-GFP` stocks). Partner symbols lived in `Co-reagent symbols` but `Partner driver symbols (best-effort)` was cleared whenever stock resolution failed.
- Decision: Keep `partner_driver_symbols` populated from genotype-detected GAL4 even without gal4-only stock candidates; fall back to co-reagent symbols in `_build_reagent_published_gal4_enrichment`; add `Partner GAL4 FlyBase IDs` / `Partner GAL4 Symbols` summary columns to all Sheet1..N; reorder columns for screening priority via `_reorder_screening_stock_sheet_columns`.
- Code: `fl_ai_reagent_stocker/pipelines/stock_splitting.py` (`_build_stock_phenotype_sheet`, `_partner_gal4_symbols_for_row`, `_build_reagent_partner_gal4_summary`, `_build_reagent_published_gal4_enrichment`, `_reorder_screening_stock_sheet_columns`, `write_aggregated_excel`)

### 2026-06-25 — Orderable driver genotype in Published GAL4 brace entries
- Context: Sheet1..N `{…}` entries listed published partner symbols and stock ids but not the FlyBase component genotype on the recommended gal4-only stock (e.g. `P{GawB}elav[C155]` for Bloomington 458).
- Decision: Extend brace entries to five fields: `{published symbol, orderable driver genotype, collection, reference ID, phenotype}` and the parallel stock-id form. Populate field 2 from gal4-only stock component symbols via `Partner driver stock genotypes` (internal column, stripped from Excel phenotype sheets like `_reference_url`).
- Code: `fl_ai_reagent_stocker/pipelines/stock_splitting.py` (`_build_stock_candidate_label_to_driver_genotype`, `_join_aligned_driver_stock_genotypes`, `_build_reagent_published_gal4_enrichment`, `PARTNER_DRIVER_STOCK_GENOTYPES_COLUMN`)

### 2026-06-25 — FBal → FBti via insertion allele TSV on phenotype sheet
- Context: `_build_stock_phenotype_sheet` only chained FBal → FBtp → FBti via construct descriptions, but Stage 1 also links FBal → FBti through `dmel_classical_and_insertion_allele_descriptions_*.tsv`. Many phenotype partner FBals (e.g. driver alleles missing from construct TSV) failed stock resolution unless symbol fingerprint fallback matched.
- Decision: Load insertion allele descriptions in `_build_stock_phenotype_sheet`, build `fbal_to_fbtis`, and union direct FBal → FBti links with the existing FBtp chain in `_resolve_fbal_stock_candidates`.
- Code: `fl_ai_reagent_stocker/pipelines/stock_splitting.py` (`_build_stock_phenotype_sheet`, `_resolve_fbal_stock_candidates`, `_stock_candidates_for_fbti_id`)

### 2026-06-25 — Collection alias normalization (BDSC ↔ Bloomington)
- Context: Partner stock labels, brace field 3, enrichment join keys, and Bloomington config filters treated `BDSC` and `Bloomington` as different strings despite referring to the same stock center.
- Decision: Canonicalize via `normalize_collection_name()` / `canonical_stock_candidate_label()` in `utils.py`; apply when formatting/parsing candidate labels, source stock keys, brace output, gal4-only genotype lookup, collection filters, and simple-bucket summaries. Display canonical name `Bloomington`.
- Code: `fl_ai_reagent_stocker/utils.py`, `fl_ai_reagent_stocker/pipelines/stock_splitting.py`

### 2026-06-25 — GawB / cross-format GAL4 driver stock resolution
- Context: Classic driver stocks such as Bloomington **458** (`P{GawB}elav[C155]`) were excluded from gal4-only partner stock candidates because `_looks_like_gal4_symbol` required the literal substring `"gal4"`. Phenotype rows cite the same driver as `Scer\GAL4[elav-C155]`, so stock ids stayed `-` despite a clean one-insertion stock existing.
- Decision: Treat `GawB`, `P{Gal4-…}`, and related driver naming as GAL4; match drivers across symbol formats via `_gal4_driver_fingerprint`; when FBal→FBti chain lookup fails, resolve gal4-only stocks from a fingerprint index built from derived stock components (e.g. `(458, Bloomington)` for elav-C155).
- Code: `fl_ai_reagent_stocker/pipelines/stock_splitting.py` (`_looks_like_gal4_symbol`, `_gal4_driver_fingerprint`, `_resolve_gal4_stock_candidates_by_driver_symbol`, `_build_gal4_driver_stock_candidate_index`)

### 2026-07-14 — GAL4 construct-alias stock resolution
- Context: Published allele names can differ from orderable insertion names; `Scer\GAL4[Mef2.PR]` and stock 600193's `P{GAL4-Mef2.R}39-1` share FlyBase construct `FBtp0006434` but have no useful text-key overlap.
- Decision: In the GAL4 workbook exporter, resolve allele and construct symbols through the authoritative `FBal → FBtp → FBti → FBst` chain before fuzzy mixed-stock matching, while excluding generic GAL4 alias keys.
- Code: `scripts/export_gal4_driver_workbook.py` (`_build_gal4_construct_alias_fbti_index`, `_fbst_map_from_construct_aliases`, `_resolve_gal4_symbol_to_stocks`)

### 2026-09-01 — Combination-sheet phenotype enrichment landed
- Context: LEARNINGS described Sheet1..N phenotype/qualifier/reference/GAL4 columns, but `write_aggregated_excel` still wrote raw Stage 1 stock rows. The June helpers were never committed.
- Decision: Build the internal phenotype sheet when `soft_run=True` or the JSON config uses a Phenotype filter, then merge `Phenotype Reference-ID Map`, `Partner GAL4 FlyBase IDs`, `Partner GAL4 Symbols`, `Published GAL4/ Positive Control`, and `Published GAL4 stock id` onto every populated combination sheet. Qualifier is folded into the phenotype text when present. Prefix PMID/PMCID values at Excel write time.
- Code: `fl_ai_reagent_stocker/pipelines/stock_splitting.py` (`_config_uses_phenotype_filter`, `_build_reagent_phenotype_enrichment`, `_add_phenotype_enrichment_to_stock_sheet`, `_reorder_screening_stock_sheet_columns`, `write_aggregated_excel`)

### 2026-07-15 — Generic GAL4 placeholders and coverage audit
- Context: FlyBase symbols such as `nSyb.PU`, `sca.PU`, and `elav.PU` describe unknown particular constructs and therefore have no authoritative stock links; the exporter silently omitted them from collection sheets.
- Decision: When exact and construct-alias resolution fail for generic `.P?` symbols, recommend one simple current BDSC stock from the same promoter family and label it `Candidate substitution—not the exact published construct`. Emit a `Coverage` sheet with one row per input and never rely on collection sheets alone to detect omissions.
- Code: `scripts/export_gal4_driver_workbook.py` (`_generic_shorthand_driver_base`, `_gal4_promoter_keys`, `_resolve_generic_promoter_driver_fallback`, `_build_coverage_sheet`)

### 2026-09-09 — Doc audit: dead/moved code references in earlier entries

This is a correction entry, not a new decision. A documentation audit found
several `Code:` pointers above no longer match the codebase. History is left
as-is; use these corrections when following those entries:

- **"Avoid double masterlist reorder" (2026-06-25):** `_is_masterlist_formatted`
  no longer exists anywhere in the repo. The reorder-skip logic it described
  was superseded by `_prepare_phenotype_workbook_output` in
  `fl_ai_reagent_stocker/pipelines/stock_splitting.py`, which is the current
  entry point for that behavior.
- **"Sheet1..N Partner GAL4 summary + symbol fallback" (2026-06-25):**
  `_partner_gal4_symbols_for_row` does not exist in `stock_splitting.py` (or
  anywhere in the repo). The partner-symbol fallback behavior it described is
  now implemented inside `_build_reagent_published_gal4_enrichment`.
- **"Orderable driver genotype in Published GAL4 brace entries" (2026-06-25):**
  `_join_aligned_driver_stock_genotypes` does not exist. The writer this entry
  describes is `_build_stock_candidate_label_to_driver_genotype`, which is
  defined twice for two different consumers: as a local helper inside
  `_build_stock_phenotype_sheet` in `stock_splitting.py` (for Sheet1..N) and as
  a module-level function in `scripts/export_gal4_driver_workbook.py` (for the
  GAL4 driver workbook). Both populate field 2 of the
  `Published GAL4/ Positive Control` brace entries.
- **"GawB / cross-format GAL4 driver stock resolution" (2026-06-25):**
  `_resolve_gal4_stock_candidates_by_driver_symbol` and
  `_build_gal4_driver_stock_candidate_index` live only in
  `scripts/export_gal4_driver_workbook.py`; full cross-format driver *stock
  resolution* still happens there, not in `stock_splitting.py`. The
  symbol-level primitives, however, are shared: `looks_like_gal4_symbol`,
  `normalized_text_looks_like_gal4_driver`, `gal4_search_keys`, and
  `gal4_driver_fingerprint` all live in `fl_ai_reagent_stocker/utils.py`, so
  the GawB-aware matching is available to the stocker and the exporter alike.
- Code (for the corrections above): `fl_ai_reagent_stocker/utils.py` (shared GAL4 symbol primitives), `fl_ai_reagent_stocker/pipelines/stock_splitting.py` (`_prepare_phenotype_workbook_output`, `_build_reagent_published_gal4_enrichment`, `_build_stock_phenotype_sheet`, `PARTNER_DRIVER_STOCK_GENOTYPES_COLUMN`), `scripts/export_gal4_driver_workbook.py` (`_resolve_gal4_stock_candidates_by_driver_symbol`, `_build_gal4_driver_stock_candidate_index`, `_build_stock_candidate_label_to_driver_genotype`)

### 2026-09-09 — Phenotype evidence and orderable GAL4 genotypes decoupled from embeddings
- Context: A `Phenotype` filter built the internal phenotype frame but hid row-level phenotype data from combination sheets unless the separate embedding pass ran; the All Phenotypic Stocks Sheet likewise required embeddings. The orderable driver-genotype field in published-GAL4 brace entries had no writer.
- Decision: A Phenotype filter now independently emits the All Phenotypic Stocks Sheet and enriches combination sheets with `Phenotype`, `Qualifier`, reference, partner-GAL4, and published-driver context. Orderable GAL4 genotype values use keyed candidate entries internally so they remain aligned after aggregation, and populate brace field 2. GawB-aware driver matching is shared by the stocker and GAL4 exporter.
- Code: `fl_ai_reagent_stocker/utils.py` (shared GAL4 symbol utilities), `fl_ai_reagent_stocker/pipelines/stock_splitting.py` (`_build_stock_phenotype_sheet`, `_build_reagent_phenotype_enrichment`, `_build_reagent_published_gal4_enrichment`, `write_aggregated_excel`)

### 2026-09-09 — Public config ingredient catalog
- Context: The only tracked configuration had documented Ref, balancer, insertion, and CRISPR concepts without executable filter definitions, and exposed no simple literature-follow-up buckets after phenotype selection.
- Decision: The tracked priority config now defines every documented filter as a reusable ingredient, appends Bloomington RNAi `Ref++` / `Ref+` buckets after phenotype buckets, and documents settings, ingredients, ordering, and shared stock-cap behavior for Cursor users.
- Code: `data/config/stock_split_config_priority_example.json`, `docs/AGENTS.md`, `docs/stock_split_config_priority_example.md`
