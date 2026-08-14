"""FBgn-based overlap tables and narrative report."""

from __future__ import annotations

from collections import Counter
from itertools import combinations
from pathlib import Path

import pandas as pd

from .constants import SET_NAMES, TRHN_FBGN, TRHN_SYMBOL
from .io_utils import join_sorted, write_csv


def membership_table(sets: dict[str, pd.DataFrame]) -> pd.DataFrame:
    union: dict[str, dict[str, object]] = {}
    for set_name in SET_NAMES:
        df = sets[set_name]
        for _, rec in df.iterrows():
            fbgn = str(rec["flybase_gene_id"])
            symbol = str(rec.get("ext_gene") or rec.get("primary_symbol") or "")
            row = union.setdefault(
                fbgn,
                {
                    "flybase_gene_id": fbgn,
                    "primary_symbol": symbol,
                    **{name: False for name in SET_NAMES},
                },
            )
            row[set_name] = True
            if symbol and not row["primary_symbol"]:
                row["primary_symbol"] = symbol
    records = []
    for fbgn, row in union.items():
        members = [name for name in SET_NAMES if row[name]]
        row["membership_count"] = len(members)
        row["membership_signature"] = join_sorted(members)
        for name in SET_NAMES:
            row[name] = "TRUE" if row[name] else "FALSE"
        records.append(row)
    if not records:
        columns = ["flybase_gene_id", "primary_symbol", *SET_NAMES, "membership_count", "membership_signature"]
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(records).sort_values(
        ["membership_count", "primary_symbol", "flybase_gene_id"],
        ascending=[False, True, True],
    ).reset_index(drop=True)


def overlapping_genes(membership: pd.DataFrame) -> pd.DataFrame:
    out = membership[membership["membership_count"] >= 2].copy()
    out["set_names"] = out["membership_signature"]
    return out[
        ["flybase_gene_id", "primary_symbol", "set_names", "membership_count", "membership_signature"]
        + SET_NAMES
    ].reset_index(drop=True)


def exact_intersections(membership: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    grouped = (
        membership.groupby("membership_signature", dropna=False)
        .agg(
            intersection_size=("flybase_gene_id", "size"),
            membership_count=("membership_count", "first"),
        )
        .reset_index()
        .sort_values(["intersection_size", "membership_signature"], ascending=[False, True])
    )
    genes = membership[["membership_signature", "flybase_gene_id", "primary_symbol", "membership_count"]]
    return grouped, genes.sort_values(["membership_signature", "primary_symbol"]).reset_index(drop=True)


def pairwise_overlap(membership: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    pair_rows = []
    gene_rows = []
    for a, b in combinations(SET_NAMES, 2):
        both = membership[(membership[a] == "TRUE") & (membership[b] == "TRUE")]
        pair_rows.append(
            {
                "set_a": a,
                "set_b": b,
                "n_overlapping_genes": int(len(both)),
            }
        )
        pair_rows.append(
            {
                "set_a": b,
                "set_b": a,
                "n_overlapping_genes": int(len(both)),
            }
        )
        for _, rec in both.iterrows():
            gene_rows.append(
                {
                    "set_a": a,
                    "set_b": b,
                    "flybase_gene_id": rec["flybase_gene_id"],
                    "primary_symbol": rec["primary_symbol"],
                }
            )
            gene_rows.append(
                {
                    "set_a": b,
                    "set_b": a,
                    "flybase_gene_id": rec["flybase_gene_id"],
                    "primary_symbol": rec["primary_symbol"],
                }
            )
    pairs = pd.DataFrame(pair_rows)
    if pairs.empty:
        pairs = pd.DataFrame(columns=["set_a", "set_b", "n_overlapping_genes"])
    else:
        pairs = pairs.sort_values(["set_a", "set_b"]).reset_index(drop=True)
    genes = pd.DataFrame(gene_rows)
    if genes.empty:
        genes = pd.DataFrame(columns=["set_a", "set_b", "flybase_gene_id", "primary_symbol"])
    else:
        genes = genes.sort_values(["set_a", "set_b", "primary_symbol"]).reset_index(drop=True)
    return pairs, genes


def write_overlap_tables(membership: pd.DataFrame, evidence_dir: Path) -> dict[str, Path]:
    overlapping = overlapping_genes(membership)
    exact, exact_genes = exact_intersections(membership)
    pairs, pair_genes = pairwise_overlap(membership)
    paths = {
        "Set_Membership.csv": evidence_dir / "Set_Membership.csv",
        "Overlapping_Genes.csv": evidence_dir / "Overlapping_Genes.csv",
        "Exact_Intersections.csv": evidence_dir / "Exact_Intersections.csv",
        "Exact_Intersection_Genes.csv": evidence_dir / "Exact_Intersection_Genes.csv",
        "Pairwise_Overlap.csv": evidence_dir / "Pairwise_Overlap.csv",
        "Pairwise_Overlap_Genes.csv": evidence_dir / "Pairwise_Overlap_Genes.csv",
    }
    write_csv(membership, paths["Set_Membership.csv"])
    write_csv(overlapping, paths["Overlapping_Genes.csv"])
    write_csv(exact, paths["Exact_Intersections.csv"])
    write_csv(exact_genes, paths["Exact_Intersection_Genes.csv"])
    write_csv(pairs, paths["Pairwise_Overlap.csv"])
    write_csv(pair_genes, paths["Pairwise_Overlap_Genes.csv"])
    return paths


def overlap_report_markdown(membership: pd.DataFrame) -> str:
    overlapping = overlapping_genes(membership)
    exact, _ = exact_intersections(membership)
    pairs, _ = pairwise_overlap(membership)
    counts = Counter(membership["membership_count"].astype(int))
    lines = [
        "# Tx-Omics Revision overlap",
        "",
        f"Unique genes across the seven sets: **{len(membership)}**.",
        "",
        "Genes by number of sets:",
    ]
    for n in sorted(counts):
        lines.append(f"- exactly {n} set(s): {counts[n]}")
    lines.extend(["", "## Largest exact intersections", ""])
    for _, rec in exact.head(12).iterrows():
        lines.append(
            f"- {int(rec['intersection_size'])} genes in exactly `{rec['membership_signature']}`"
        )
    directed = pairs[pairs["set_a"] < pairs["set_b"]].sort_values(
        "n_overlapping_genes", ascending=False
    )
    lines.extend(["", "## Largest pairwise overlaps", ""])
    for _, rec in directed.head(12).iterrows():
        lines.append(
            f"- {rec['set_a']} ∩ {rec['set_b']}: {int(rec['n_overlapping_genes'])} genes"
        )

    categories = [
        "Mechanistic",
        "Homeostatic History/Rebound",
        "CSW 4+",
        "HLH Upstream Regulators",
    ]
    pathways = [
        "CSW Ribosomal/Translation",
        "CSW Mitochondrial/Metabolism",
        "CSW Immune",
    ]
    multi_category = membership[
        membership[categories].apply(lambda row: sum(val == "TRUE" for val in row) >= 2, axis=1)
    ]
    lines.extend(["", "## Genes in more than one publication-defined category", ""])
    if multi_category.empty:
        lines.append("None.")
    else:
        for _, rec in multi_category.iterrows():
            names = [name for name in categories if rec[name] == "TRUE"]
            lines.append(f"- {rec['primary_symbol']} (`{rec['flybase_gene_id']}`): {', '.join(names)}")

    category_in_pathways = membership[
        (membership[categories].apply(lambda row: any(val == "TRUE" for val in row), axis=1))
        & (membership[pathways].apply(lambda row: any(val == "TRUE" for val in row), axis=1))
    ]
    lines.extend(["", "## Publication-category genes also in a pathway set", ""])
    if category_in_pathways.empty:
        lines.append("None.")
    else:
        for _, rec in category_in_pathways.iterrows():
            names = [name for name in SET_NAMES if rec[name] == "TRUE"]
            lines.append(f"- {rec['primary_symbol']} (`{rec['flybase_gene_id']}`): {', '.join(names)}")

    multi_path = membership[
        membership[pathways].apply(lambda row: sum(val == "TRUE" for val in row) >= 2, axis=1)
    ]
    lines.extend(["", "## Genes spanning multiple pathway sets", ""])
    if multi_path.empty:
        lines.append("None.")
    else:
        for _, rec in multi_path.iterrows():
            names = [name for name in pathways if rec[name] == "TRUE"]
            lines.append(f"- {rec['primary_symbol']} (`{rec['flybase_gene_id']}`): {', '.join(names)}")

    lines.extend(["", "## High-priority gene membership", ""])
    priority = membership[
        (membership["flybase_gene_id"] == TRHN_FBGN)
        | (membership["primary_symbol"].isin(["Trhn", "AstA-R2", "RpL23", "unc79", "SIFa", "rumpel", "CG9377", "bigmax", "HLH3B", "E(spl)m3-HLH", "E(spl)mbeta-HLH"]))
    ]
    if priority.empty:
        lines.append("No high-priority genes were present in the seven-set union.")
    else:
        for _, rec in priority.iterrows():
            names = [name for name in SET_NAMES if rec[name] == "TRUE"]
            lines.append(f"- {rec['primary_symbol']} (`{rec['flybase_gene_id']}`): {', '.join(names)}")
    lines.append("")
    return "\n".join(lines)
