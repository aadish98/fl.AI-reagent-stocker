"""Build FDR-filtered pathway-hit gene sets from approved GO terms."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .constants import (
    BUCKET_TO_SET_NAME,
    BUCKET_TO_STEM,
    PATHWAY_BUCKETS,
)
from .io_utils import join_sorted, write_csv
from .keywords import adjacent_unmatched, match_buckets, proposed_decision


def build_term_review(
    terms: pd.DataFrame, genes: pd.DataFrame
) -> pd.DataFrame:
    if terms.empty:
        return pd.DataFrame()
    gene_counts = (
        genes.groupby(["source_table", "category", "term_id"], dropna=False)
        .size()
        .rename("hit_gene_count")
        .reset_index()
    )
    review = terms.merge(
        gene_counts, on=["source_table", "category", "term_id"], how="left"
    )
    matched_names = []
    bucket_lists = []
    rule_lists = []
    for name in review["term_name"]:
        buckets, rules = match_buckets(name)
        bucket_lists.append(buckets)
        rule_lists.append(rules)
        if buckets:
            matched_names.append(name)
    review["candidate_buckets"] = [";".join(items) for items in bucket_lists]
    review["matched_rule"] = [";".join(items) for items in rule_lists]
    review["conflict_flag"] = ["TRUE" if len(items) > 1 else "FALSE" for items in bucket_lists]
    review["proposed_decision"] = [proposed_decision(items) for items in bucket_lists]
    review["adjacent_unmatched"] = [
        "FALSE"
        if buckets
        else ("TRUE" if adjacent_unmatched(name, matched_names) else "FALSE")
        for name, buckets in zip(review["term_name"], bucket_lists)
    ]
    review["user_approved_bucket"] = ""
    review["approval_status"] = "pending"
    keep_mask = (review["proposed_decision"] != "exclude") | (
        review["adjacent_unmatched"] == "TRUE"
    )
    return review.loc[keep_mask].sort_values(
        ["proposed_decision", "fdr_percent", "term_name"]
    ).reset_index(drop=True)


def approve_proposed_terms(review: pd.DataFrame, conflict_policy: str) -> pd.DataFrame:
    approved = review.copy()
    buckets = []
    statuses = []
    for _, rec in approved.iterrows():
        decision = rec["proposed_decision"]
        candidates = [part for part in str(rec["candidate_buckets"]).split(";") if part]
        if decision == "include" and len(candidates) == 1:
            buckets.append(candidates[0])
            statuses.append("approved_proposed")
        elif decision == "conflict":
            if conflict_policy == "include_all_matched_buckets":
                buckets.append(";".join(candidates))
                statuses.append("approved_conflict_all_matched_buckets")
            else:
                buckets.append("")
                statuses.append("pending_conflict")
        else:
            buckets.append("")
            statuses.append("excluded_proposed")
    approved["user_approved_bucket"] = buckets
    approved["approval_status"] = statuses
    return approved


def _observation_lookup(observations: pd.DataFrame) -> dict[tuple[str, str, str], pd.Series]:
    lookup = {}
    for _, rec in observations.iterrows():
        lookup[(rec["flybase_gene_id"], rec["threshold"], rec["direction"])] = rec
    return lookup


def _obs_value(lookup, fbgn: str, threshold: str, direction: str, column: str) -> str:
    rec = lookup.get((fbgn, threshold, direction))
    if rec is None:
        return ""
    return "" if rec.get(column) is None else str(rec.get(column))


def build_pathway_tables(
    approved: pd.DataFrame,
    genes: pd.DataFrame,
    observations: pd.DataFrame,
    category_ids: dict[str, set[str]],
    csw4_thresholds: dict[str, str],
) -> dict[str, pd.DataFrame]:
    approved_terms = approved[approved["user_approved_bucket"].astype(str) != ""].copy()
    lookup = _observation_lookup(observations)
    symbol_by_fbgn = (
        observations.sort_values(["threshold", "direction"])
        .drop_duplicates("flybase_gene_id")
        .set_index("flybase_gene_id")["ext_gene"]
        .to_dict()
    )
    outputs: dict[str, pd.DataFrame] = {}
    for bucket in PATHWAY_BUCKETS:
        term_hits = approved_terms[
            approved_terms["user_approved_bucket"].str.split(";").apply(lambda items: bucket in items)
        ]
        if term_hits.empty:
            outputs[bucket] = pd.DataFrame(columns=_pathway_columns())
            continue
        keys = term_hits[["source_table", "category", "term_id", "term_name", "fdr_percent"]]
        gene_cols = [col for col in genes.columns if col != "term_name"]
        merged = genes[gene_cols].merge(
            keys,
            on=["source_table", "category", "term_id"],
            how="inner",
        )
        rows = []
        for fbgn, group in merged.groupby("flybase_gene_id"):
            source_tables = sorted(set(group["source_table"]))
            directions = sorted(
                {
                    "wake" if "Wake" in table else "sleep"
                    for table in source_tables
                }
            )
            thresholds = sorted(
                {
                    "FC0.5" if "FC0.5" in table else "FC1"
                    for table in source_tables
                }
            )
            rows.append(
                {
                    "ext_gene": symbol_by_fbgn.get(fbgn, ""),
                    "source_flybase_gene_id": fbgn,
                    "flybase_gene_id": fbgn,
                    "pathway_bucket": bucket,
                    "source_tables": join_sorted(source_tables),
                    "directions": join_sorted(directions),
                    "thresholds": join_sorted(thresholds),
                    "go_term_ids": join_sorted(group["term_id"]),
                    "go_term_names": join_sorted(group["term_name"]),
                    "go_categories": join_sorted(group["category"]),
                    "go_evidence_count": str(len(group)),
                    "go_fdr_min": f"{pd.to_numeric(group['fdr_percent'], errors='coerce').min():.6g}",
                    "frequency_FC0.5_wake": _obs_value(lookup, fbgn, "FC0.5", "wake", "frequency"),
                    "frequency_FC1_wake": _obs_value(lookup, fbgn, "FC1", "wake", "frequency"),
                    "frequency_FC0.5_sleep": _obs_value(lookup, fbgn, "FC0.5", "sleep", "frequency"),
                    "frequency_FC1_sleep": _obs_value(lookup, fbgn, "FC1", "sleep", "frequency"),
                    "wake_corr_exps_FC0.5": _obs_value(lookup, fbgn, "FC0.5", "wake", "wake_corr_exps"),
                    "wake_corr_exps_FC1": _obs_value(lookup, fbgn, "FC1", "wake", "wake_corr_exps"),
                    "sleep_corr_exps_FC0.5": _obs_value(lookup, fbgn, "FC0.5", "sleep", "sleep_corr_exps"),
                    "sleep_corr_exps_FC1": _obs_value(lookup, fbgn, "FC1", "sleep", "sleep_corr_exps"),
                    "is_cycling_FC0.5_wake": _obs_value(lookup, fbgn, "FC0.5", "wake", "is_cycling"),
                    "is_cycling_FC1_wake": _obs_value(lookup, fbgn, "FC1", "wake", "is_cycling"),
                    "is_cycling_FC0.5_sleep": _obs_value(lookup, fbgn, "FC0.5", "sleep", "is_cycling"),
                    "is_cycling_FC1_sleep": _obs_value(lookup, fbgn, "FC1", "sleep", "is_cycling"),
                    "ZT_cyclers_FC0.5_wake": _obs_value(lookup, fbgn, "FC0.5", "wake", "ZT_cyclers"),
                    "ZT_cyclers_FC1_wake": _obs_value(lookup, fbgn, "FC1", "wake", "ZT_cyclers"),
                    "ZT_cyclers_FC0.5_sleep": _obs_value(lookup, fbgn, "FC0.5", "sleep", "ZT_cyclers"),
                    "ZT_cyclers_FC1_sleep": _obs_value(lookup, fbgn, "FC1", "sleep", "ZT_cyclers"),
                    "in_CSW_4plus": "TRUE" if fbgn in category_ids["CSW 4+"] else "FALSE",
                    "csw_4plus_thresholds": csw4_thresholds.get(fbgn, ""),
                    "in_Mechanistic_6": "TRUE" if fbgn in category_ids["Mechanistic"] else "FALSE",
                    "in_Homeostatic_20": "TRUE" if fbgn in category_ids["Homeostatic History/Rebound"] else "FALSE",
                    "in_HLH_4": "TRUE" if fbgn in category_ids["HLH Upstream Regulators"] else "FALSE",
                }
            )
        out = pd.DataFrame(rows).sort_values(["ext_gene", "flybase_gene_id"]).reset_index(drop=True)
        outputs[bucket] = out
    return outputs


def _pathway_columns() -> list[str]:
    return [
        "ext_gene",
        "source_flybase_gene_id",
        "flybase_gene_id",
        "pathway_bucket",
        "source_tables",
        "directions",
        "thresholds",
        "go_term_ids",
        "go_term_names",
        "go_categories",
        "go_evidence_count",
        "go_fdr_min",
        "frequency_FC0.5_wake",
        "frequency_FC1_wake",
        "frequency_FC0.5_sleep",
        "frequency_FC1_sleep",
        "wake_corr_exps_FC0.5",
        "wake_corr_exps_FC1",
        "sleep_corr_exps_FC0.5",
        "sleep_corr_exps_FC1",
        "is_cycling_FC0.5_wake",
        "is_cycling_FC1_wake",
        "is_cycling_FC0.5_sleep",
        "is_cycling_FC1_sleep",
        "ZT_cyclers_FC0.5_wake",
        "ZT_cyclers_FC1_wake",
        "ZT_cyclers_FC0.5_sleep",
        "ZT_cyclers_FC1_sleep",
        "in_CSW_4plus",
        "csw_4plus_thresholds",
        "in_Mechanistic_6",
        "in_Homeostatic_20",
        "in_HLH_4",
    ]


def reconstruct_summary_ok(pathway_df: pd.DataFrame, genes: pd.DataFrame, approved: pd.DataFrame) -> bool:
    if pathway_df.empty:
        return True
    bucket = pathway_df["pathway_bucket"].iloc[0]
    term_hits = approved[
        approved["user_approved_bucket"].str.split(";").apply(lambda items: bucket in items)
    ]
    expected = set(
        genes.merge(
            term_hits[["source_table", "category", "term_id"]],
            on=["source_table", "category", "term_id"],
            how="inner",
        )["flybase_gene_id"]
    )
    got = set(pathway_df["flybase_gene_id"])
    return expected == got
