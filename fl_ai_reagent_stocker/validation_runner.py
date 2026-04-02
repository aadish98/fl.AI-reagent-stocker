"""
Shared functional-validation runner used by both pipelines.
"""

import time
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
from tqdm import tqdm

from .config import ValidationStatus
from .utils import clean_id, parse_semicolon_list, generate_keyword_column_name


def run_functional_validation(
    stocks_df: pd.DataFrame,
    references_df: Optional[pd.DataFrame],
    keywords: Optional[List[str]],
    soft_run: bool,
    *,
    validator,
    fulltext_fetcher,
    pubmed_client,
    pubmed_cache,
    max_gpt_calls_per_stock: Optional[int],
    component_to_pmids: Optional[Dict[str, Set[str]]] = None,
    component_id_to_symbol: Optional[Dict[str, str]] = None,
    stock_tasks_override: Optional[Dict[str, List[Tuple[str, str]]]] = None,
) -> pd.DataFrame:
    """
    Run functional validation for stock/reference tasks.

    Two modes are supported:
    1) Standard mode (pipeline_references): provide ``component_to_pmids`` and
       leave ``stock_tasks_override`` as None.
    2) Override mode (pipeline_split): provide ``stock_tasks_override`` with
       pre-filtered (stock -> [(allele, pmid), ...]) tasks.
    """
    kw_col_name = generate_keyword_column_name(keywords) if keywords else None
    validated_col = (
        f"PMID of {' OR '.join(keywords)} references that showed stocks functional validity"
        if keywords
        else "PMID of keyword references that showed stocks functional validity"
    )

    stocks_df = stocks_df.copy()
    stocks_df["Functional Validity (by reference)?"] = ""
    stocks_df["Functionally Valid Stock?"] = ""
    stocks_df[validated_col] = ""
    stocks_df["validation_rationale"] = ""
    stocks_df["phenotypes_observed"] = ""

    if "total_refs" in stocks_df.columns:
        no_refs_mask = stocks_df["total_refs"].fillna(0).astype(int) <= 0
    else:
        no_refs_mask = (
            stocks_df["PMID"].isna()
            | (stocks_df["PMID"].astype(str).str.strip() == "")
            | (stocks_df["PMID"].astype(str).str.strip() == "-")
        )
    stocks_df.loc[no_refs_mask, "Functionally Valid Stock?"] = ValidationStatus.NO_REFS

    if not validator.is_available:
        print("    Warning: OpenAI not available. Skipping functional validation.")
        stocks_df.loc[~no_refs_mask, "Functionally Valid Stock?"] = ValidationStatus.AMBIGUOUS
        return stocks_df

    keyword_pmids: Set[str] = set()
    if kw_col_name and references_df is not None and kw_col_name in references_df.columns:
        keyword_refs = references_df[references_df[kw_col_name] == True]
        keyword_pmids = set(keyword_refs["PMID"].apply(clean_id))
        keyword_pmids.discard("")

    print(f"\n    Found {len(keyword_pmids)} keyword-matching PMIDs")

    stocks_with_refs = stocks_df[~no_refs_mask].copy()
    if len(stocks_with_refs) == 0 or len(keyword_pmids) == 0:
        stocks_df.loc[~no_refs_mask, "Functionally Valid Stock?"] = ValidationStatus.AMBIGUOUS
        return stocks_df

    print(f"    Fetching publication dates for {len(keyword_pmids)} PMIDs...")
    pmid_dates: Dict[str, dict] = pubmed_client.fetch_author_date_metadata(list(keyword_pmids))
    print(f"    Retrieved dates for {len(pmid_dates)} PMIDs")

    def get_date_key(pmid: str) -> Tuple[int, int, int]:
        date_info = pmid_dates.get(pmid, {})
        date_str = date_info.get("date_published", "") if isinstance(date_info, dict) else ""
        if not date_str:
            return (0, 0, 0)
        parts = date_str.split("-")
        try:
            year = int(parts[0]) if len(parts) > 0 and parts[0].isdigit() else 0
            month = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0
            day = int(parts[2]) if len(parts) > 2 and parts[2].isdigit() else 0
            return (year, month, day)
        except (ValueError, IndexError):
            return (0, 0, 0)

    stock_tasks: Dict[str, List[Tuple[str, str]]] = {}
    stock_is_custom: Dict[str, bool] = {}
    stock_gene_terms: Dict[str, List[str]] = {}

    if stock_tasks_override is None:
        if component_to_pmids is None:
            raise ValueError("component_to_pmids must be provided when stock_tasks_override is not set")

        id_to_sym = component_id_to_symbol or {}

        for _, row in stocks_with_refs.iterrows():
            stock_num = clean_id(row.get("stock_number", ""))
            if not stock_num:
                continue

            relevant_ids = parse_semicolon_list(row.get("relevant_component_ids", ""))
            _custom_flag = row.get("custom_stock", row.get("is_custom_stock", False))
            stock_is_custom[stock_num] = _custom_flag is True or str(_custom_flag).strip().lower() in ("true", "1", "yes")

            gene_terms: List[str] = []
            seen_terms: Set[str] = set()
            for sym in parse_semicolon_list(row.get("relevant_gene_symbols", "")):
                sym_clean = str(sym).strip()
                if not sym_clean:
                    continue
                k = sym_clean.lower()
                if k in seen_terms:
                    continue
                seen_terms.add(k)
                gene_terms.append(sym_clean)
            stock_gene_terms[stock_num] = gene_terms

            all_pairs: List[Tuple[str, str]] = []
            for cid in relevant_ids:
                allele_symbol = id_to_sym.get(cid, cid)
                for pmid in component_to_pmids.get(cid, set()):
                    pmid_clean = clean_id(pmid)
                    if pmid_clean and pmid_clean in keyword_pmids:
                        all_pairs.append((allele_symbol, pmid_clean))

            all_pairs.sort(key=lambda x: get_date_key(x[1]), reverse=True)
            stock_tasks[stock_num] = all_pairs
    else:
        stock_indexed = {}
        for _, row in stocks_with_refs.iterrows():
            stock_num = clean_id(row.get("stock_number", ""))
            if stock_num:
                stock_indexed[stock_num] = row

        for stock_num, tasks in stock_tasks_override.items():
            stock_clean = clean_id(stock_num)
            row = stock_indexed.get(stock_clean)
            if row is None:
                continue

            _custom_flag = row.get("custom_stock", row.get("is_custom_stock", False))
            stock_is_custom[stock_clean] = _custom_flag is True or str(_custom_flag).strip().lower() in ("true", "1", "yes")
            gene_terms: List[str] = []
            seen_terms: Set[str] = set()
            for sym in parse_semicolon_list(row.get("relevant_gene_symbols", "")):
                sym_clean = str(sym).strip()
                if not sym_clean:
                    continue
                k = sym_clean.lower()
                if k in seen_terms:
                    continue
                seen_terms.add(k)
                gene_terms.append(sym_clean)
            stock_gene_terms[stock_clean] = gene_terms

            deduped: Set[Tuple[str, str]] = set()
            for allele, pmid in tasks:
                pmid_clean = clean_id(pmid)
                allele_clean = str(allele).strip()
                if not allele_clean or not pmid_clean:
                    continue
                if pmid_clean in keyword_pmids:
                    deduped.add((allele_clean, pmid_clean))
            sorted_pairs = sorted(list(deduped), key=lambda x: get_date_key(x[1]), reverse=True)
            stock_tasks[stock_clean] = sorted_pairs

    total_tasks = sum(len(tasks) for tasks in stock_tasks.values())
    num_stocks_with_refs = len([k for k, v in stock_tasks.items() if v])
    if max_gpt_calls_per_stock is not None:
        print(f"\n    Unique stocks with keyword-matching refs: {num_stocks_with_refs}")
        print(f"    Total keyword-matching references to process: {total_tasks}")
        print(f"    Max OpenAI API calls per stock: {max_gpt_calls_per_stock}")
        print("    Note: Only actual GPT calls count against this limit.")
        print("          References without accessible text or stock/allele patterns are marked ambiguous and don't count.")
    else:
        print(f"\n    Unique stocks with keyword-matching refs: {num_stocks_with_refs}")
        print(f"    Total keyword-matching references: {total_tasks}")
        print("    Max OpenAI API calls per stock: unlimited")

    if soft_run:
        predicted_upper_bound = total_tasks
        if max_gpt_calls_per_stock is not None:
            predicted_upper_bound = min(
                total_tasks,
                num_stocks_with_refs * max_gpt_calls_per_stock
            )
        print("\n    SOFT RUN: Stopping before OpenAI API calls.")
        print(f"    Predicted GPT query candidates (raw): {total_tasks}")
        if max_gpt_calls_per_stock is not None:
            print(f"    Predicted GPT query upper bound (with per-stock cap): {predicted_upper_bound}")
        else:
            print("    Predicted GPT query upper bound (with per-stock cap): unlimited cap configured")
        print("    Note: Actual calls may be lower due to missing full text, missing stock/allele/gene patterns,")
        print("          and short-circuiting after 'Functionally validated' results.")
        return stocks_df

    refs_by_pmid: Dict[str, dict] = {}
    if references_df is not None and "PMID" in references_df.columns:
        for _, ref_row in references_df.iterrows():
            pmid = clean_id(ref_row.get("PMID", ""))
            if not pmid:
                continue
            # Merge duplicate PMID rows and keep whichever row has richer identifiers.
            if pmid not in refs_by_pmid:
                refs_by_pmid[pmid] = ref_row
            else:
                existing = refs_by_pmid[pmid]
                existing_pmcid = clean_id(existing.get("PMCID", "")) if pd.notna(existing.get("PMCID", "")) else ""
                existing_doi = clean_id(existing.get("DOI", "")) if pd.notna(existing.get("DOI", "")) else ""
                new_pmcid = clean_id(ref_row.get("PMCID", "")) if pd.notna(ref_row.get("PMCID", "")) else ""
                new_doi = clean_id(ref_row.get("DOI", "")) if pd.notna(ref_row.get("DOI", "")) else ""
                existing_score = int(bool(existing_pmcid)) + int(bool(existing_doi))
                new_score = int(bool(new_pmcid)) + int(bool(new_doi))
                if new_score > existing_score:
                    refs_by_pmid[pmid] = ref_row

    results: Dict[Tuple[str, str, str], dict] = {}
    validated_stocks: Set[str] = set()
    stock_gpt_call_count: Dict[str, int] = {}
    skipped_short_circuit = 0
    skipped_limit_reached = 0
    skipped_no_fulltext = 0
    no_fulltext_reasons: Dict[str, int] = {}
    skipped_no_patterns = 0
    actual_gpt_calls = 0

    if total_tasks > 0:
        print("\n    Running functional validation...")
        pbar = tqdm(total=total_tasks, desc="Validation", unit="ref", ncols=80)

        for stock_num, tasks in stock_tasks.items():
            is_custom = stock_is_custom.get(stock_num, False)
            stock_gpt_call_count[stock_num] = 0

            for allele, pmid in tasks:
                if stock_num in validated_stocks:
                    skipped_short_circuit += 1
                    pbar.update(1)
                    continue

                pmcid = None
                doi = None
                ref_row = refs_by_pmid.get(clean_id(pmid))
                if ref_row is not None:
                    pmcid_val = ref_row.get("PMCID", "")
                    if pd.notna(pmcid_val) and str(pmcid_val).strip():
                        pmcid = str(pmcid_val).strip()
                    doi_val = ref_row.get("DOI", "")
                    if pd.notna(doi_val) and str(doi_val).strip():
                        doi = str(doi_val).strip()

                full_text = ""
                fetch_reason = ""
                if hasattr(fulltext_fetcher, "fetch_with_diagnostics"):
                    full_text, _source, fetch_reason = fulltext_fetcher.fetch_with_diagnostics(pmid, pmcid, doi)
                else:
                    full_text, _source = fulltext_fetcher.fetch(pmid, pmcid, doi)
                if not full_text or len(full_text) <= 500:
                    reason_suffix = f" ({fetch_reason})" if fetch_reason else ""
                    results[(stock_num, allele, pmid)] = {
                        "functional_validity": ValidationStatus.AMBIGUOUS,
                        "phenotypes": None,
                        "confidence": 0,
                        "rationale": f"Full text not accessible{reason_suffix}",
                        "from_gpt": False,
                    }
                    skipped_no_fulltext += 1
                    if fetch_reason:
                        no_fulltext_reasons[fetch_reason] = no_fulltext_reasons.get(fetch_reason, 0) + 1
                    pbar.update(1)
                    continue

                gene_terms = stock_gene_terms.get(stock_num, [])
                patterns_found = validator.check_patterns_exist(
                    full_text, stock_num, allele, is_custom, gene_terms=gene_terms
                )
                if not patterns_found:
                    if is_custom:
                        rationale = f"Allele symbol '{allele}' and fallback gene terms not found in paper text"
                    else:
                        rationale = f"Stock number '{stock_num}', allele symbol '{allele}', and fallback gene terms not found in paper text"
                    results[(stock_num, allele, pmid)] = {
                        "functional_validity": ValidationStatus.AMBIGUOUS,
                        "phenotypes": None,
                        "confidence": 95,
                        "rationale": rationale,
                        "from_gpt": False,
                    }
                    skipped_no_patterns += 1
                    pbar.update(1)
                    continue

                if max_gpt_calls_per_stock is not None and stock_gpt_call_count[stock_num] >= max_gpt_calls_per_stock:
                    results[(stock_num, allele, pmid)] = {
                        "functional_validity": ValidationStatus.AMBIGUOUS,
                        "phenotypes": None,
                        "confidence": 0,
                        "rationale": f"GPT call limit ({max_gpt_calls_per_stock}) reached for this stock",
                        "from_gpt": False,
                    }
                    skipped_limit_reached += 1
                    pbar.update(1)
                    continue

                meta = pubmed_cache.get(pmid) or {}
                result = validator.validate(
                    stock_number=stock_num,
                    allele_symbol=allele,
                    full_text=full_text,
                    title=meta.get("title", ""),
                    abstract=meta.get("abstract", ""),
                    is_custom_stock=is_custom,
                    gene_terms=gene_terms,
                )
                time.sleep(0.5)

                results[(stock_num, allele, pmid)] = result
                stock_gpt_call_count[stock_num] += 1
                actual_gpt_calls += 1

                if result["functional_validity"] == ValidationStatus.FUNCTIONALLY_VALIDATED:
                    validated_stocks.add(stock_num)

                pbar.update(1)

        pbar.close()
        print("\n    Processing summary:")
        print(f"      Actual GPT API calls: {actual_gpt_calls}")
        if skipped_short_circuit > 0:
            print(f"      Skipped (short-circuit after validation): {skipped_short_circuit}")
        if skipped_limit_reached > 0:
            print(f"      Skipped (GPT call limit reached): {skipped_limit_reached}")
        if skipped_no_fulltext > 0:
            print(f"      Skipped (no full text accessible): {skipped_no_fulltext}")
            if no_fulltext_reasons:
                print("        Full-text miss reasons:")
                for reason, count in sorted(no_fulltext_reasons.items(), key=lambda x: (-x[1], x[0])):
                    print(f"          - {reason}: {count}")
        if skipped_no_patterns > 0:
            print(f"      Skipped (no stock/allele/gene patterns in text): {skipped_no_patterns}")

    def aggregate_and_format_results(row):
        stock_num = clean_id(row.get("stock_number", ""))
        per_ref_results: List[str] = []
        validated_pmids: List[str] = []
        all_rationales: List[str] = []
        all_phenotypes: List[str] = []
        validities: List[str] = []

        for allele, pmid in stock_tasks.get(stock_num, []):
            key = (stock_num, allele, clean_id(pmid))
            if key not in results:
                continue
            result = results[key]
            validity = result["functional_validity"]
            validities.append(validity)
            pmid_clean = clean_id(pmid)
            # Keep PMID-tagged payload format consistent across GPT-derived columns
            # so stock-sheet per-reference extraction can parse by PMID reliably.
            per_ref_results.append(f"PMID:{pmid_clean}: {validity}")
            if validity == ValidationStatus.FUNCTIONALLY_VALIDATED:
                validated_pmids.append(pmid_clean)
            if result.get("from_gpt") and result.get("rationale"):
                all_rationales.append(f"PMID:{pmid_clean}: {result['rationale']}")
            if result.get("from_gpt") and result.get("phenotypes"):
                all_phenotypes.append(f"PMID:{pmid_clean}: {result['phenotypes']}")

        if not validities:
            total_refs_val = row.get("total_refs", 0)
            try:
                has_any_refs = int(float(total_refs_val)) > 0
            except (ValueError, TypeError):
                has_any_refs = False
            overall_validity = ValidationStatus.AMBIGUOUS if has_any_refs else ValidationStatus.NO_REFS
        else:
            overall_validity = ValidationStatus.AMBIGUOUS
            for val in [
                ValidationStatus.FUNCTIONALLY_VALIDATED,
                ValidationStatus.TESTED_NO_PHENOTYPE,
                ValidationStatus.AMBIGUOUS,
            ]:
                if val in validities:
                    overall_validity = val
                    break

        return pd.Series(
            {
                "Functional Validity (by reference)?": "; ".join(per_ref_results),
                "Functionally Valid Stock?": overall_validity,
                validated_col: ", ".join(validated_pmids),
                "validation_rationale": " | ".join(all_rationales),
                "phenotypes_observed": " | ".join(all_phenotypes),
            }
        )

    if len(stocks_with_refs) > 0:
        result_cols = stocks_df.loc[~no_refs_mask].apply(aggregate_and_format_results, axis=1)
        for col in result_cols.columns:
            stocks_df.loc[~no_refs_mask, col] = result_cols[col]

    print("\n    Validation complete:")
    print(
        "      Functionally validated: "
        f"{(stocks_df['Functionally Valid Stock?'] == ValidationStatus.FUNCTIONALLY_VALIDATED).sum()}"
    )
    print(
        "      Tested, no phenotype: "
        f"{(stocks_df['Functionally Valid Stock?'] == ValidationStatus.TESTED_NO_PHENOTYPE).sum()}"
    )
    print(
        "      Ambiguous: "
        f"{(stocks_df['Functionally Valid Stock?'] == ValidationStatus.AMBIGUOUS).sum()}"
    )
    print(
        "      No references: "
        f"{(stocks_df['Functionally Valid Stock?'] == ValidationStatus.NO_REFS).sum()}"
    )

    return stocks_df
