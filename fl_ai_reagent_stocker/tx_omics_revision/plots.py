"""Deterministic Tx-Omics Revision figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

from .constants import PATHWAY_BUCKETS, SET_NAMES
from .io_utils import write_csv


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_pathway_gene_counts(pathway_sets: dict[str, pd.DataFrame], dest: Path, data_dir: Path) -> None:
    rows = []
    for bucket in PATHWAY_BUCKETS:
        df = pathway_sets.get(bucket, pd.DataFrame())
        wake_only = sleep_only = both = 0
        for _, rec in df.iterrows():
            dirs = set(str(rec.get("directions", "")).split(";")) - {""}
            if dirs == {"wake"}:
                wake_only += 1
            elif dirs == {"sleep"}:
                sleep_only += 1
            elif "wake" in dirs and "sleep" in dirs:
                both += 1
        rows.append(
            {
                "pathway_bucket": bucket,
                "wake_only": wake_only,
                "sleep_only": sleep_only,
                "both": both,
                "unique_total": wake_only + sleep_only + both,
            }
        )
    data = pd.DataFrame(rows)
    write_csv(data, data_dir / "pathway_gene_counts.csv")
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    x = np.arange(len(data))
    ax.bar(x, data["wake_only"], label="wake only")
    ax.bar(x, data["sleep_only"], bottom=data["wake_only"], label="sleep only")
    ax.bar(
        x,
        data["both"],
        bottom=data["wake_only"] + data["sleep_only"],
        label="wake and sleep",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(data["pathway_bucket"], rotation=20, ha="right")
    ax.set_ylabel("Unique genes")
    ax.set_title("Pathway gene counts")
    for i, total in enumerate(data["unique_total"]):
        ax.text(i, total + 0.3, str(int(total)), ha="center", va="bottom", fontsize=8)
    if data["unique_total"].sum() == 0:
        ax.text(0.5, 0.5, "No genes", transform=ax.transAxes, ha="center")
    ax.legend()
    _save(fig, dest)


def plot_go_terms(approved: pd.DataFrame, dest: Path, data_dir: Path) -> None:
    rows = []
    for _, rec in approved.iterrows():
        buckets = [part for part in str(rec.get("user_approved_bucket", "")).split(";") if part]
        for bucket in buckets:
            rows.append(
                {
                    "pathway_bucket": bucket,
                    "term_name": rec["term_name"],
                    "term_id": rec["term_id"],
                    "hit_gene_count": rec.get("hit_gene_count", rec.get("count", "")),
                    "fdr_percent": rec.get("fdr_percent", ""),
                    "source_table": rec.get("source_table", ""),
                }
            )
    data = pd.DataFrame(rows)
    write_csv(data, data_dir / "go_terms_per_bucket.csv")
    if data.empty:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No approved GO terms", ha="center", va="center")
        ax.axis("off")
        _save(fig, dest)
        return
    grouped = (
        data.groupby(["pathway_bucket", "term_name"], as_index=False)
        .agg(hit_gene_count=("hit_gene_count", "max"), fdr_percent=("fdr_percent", "min"))
    )
    grouped["hit_gene_count"] = pd.to_numeric(grouped["hit_gene_count"], errors="coerce").fillna(0)
    grouped["fdr_percent"] = pd.to_numeric(grouped["fdr_percent"], errors="coerce")
    n = max(len(grouped), 1)
    fig, ax = plt.subplots(figsize=(10, max(4, 0.28 * n)))
    grouped = grouped.sort_values(["pathway_bucket", "hit_gene_count"])
    colors = {
        "Ribosomal/translation": "#4C78A8",
        "Mitochondrial/metabolism": "#F58518",
        "Immune": "#54A24B",
    }
    y = np.arange(len(grouped))
    ax.barh(y, grouped["hit_gene_count"], color=[colors.get(b, "#888") for b in grouped["pathway_bucket"]])
    labels = [
        f"{row.pathway_bucket}: {row.term_name} (FDR {row.fdr_percent:g}%)"
        if pd.notna(row.fdr_percent)
        else f"{row.pathway_bucket}: {row.term_name}"
        for row in grouped.itertuples()
    ]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Hit gene count")
    ax.set_title("Approved GO/KEGG terms by pathway bucket")
    _save(fig, dest)


def plot_pathway_venn(pathway_sets: dict[str, pd.DataFrame], dest: Path, data_dir: Path) -> None:
    from matplotlib_venn import venn3

    sets = []
    for bucket in PATHWAY_BUCKETS:
        df = pathway_sets.get(bucket, pd.DataFrame())
        sets.append(set(df["flybase_gene_id"]) if not df.empty else set())
    write_csv(
        pd.DataFrame(
            {
                "pathway_bucket": PATHWAY_BUCKETS,
                "n_genes": [len(s) for s in sets],
            }
        ),
        data_dir / "pathway_venn_sizes.csv",
    )
    fig, ax = plt.subplots(figsize=(8, 6))
    if all(len(s) == 0 for s in sets):
        ax.text(0.5, 0.5, "No pathway genes", ha="center", va="center")
        ax.axis("off")
    else:
        try:
            venn3(sets, set_labels=PATHWAY_BUCKETS, ax=ax)
            ax.set_title("Pathway gene overlap")
        except ValueError:
            ax.text(
                0.5,
                0.5,
                "Pathway sets are empty, identical, or cannot form a 3-set Venn",
                ha="center",
                va="center",
            )
            ax.axis("off")
    _save(fig, dest)


def plot_pairwise_counts(pairs: pd.DataFrame, dest: Path, data_dir: Path) -> None:
    write_csv(pairs, data_dir / "seven_set_pairwise_counts.csv")
    matrix = pd.DataFrame(0, index=SET_NAMES, columns=SET_NAMES, dtype=int)
    for _, rec in pairs.iterrows():
        matrix.loc[rec["set_a"], rec["set_b"]] = int(rec["n_overlapping_genes"])
    fig, ax = plt.subplots(figsize=(9, 7))
    im = ax.imshow(matrix.values, cmap="Blues")
    ax.set_xticks(range(len(SET_NAMES)))
    ax.set_yticks(range(len(SET_NAMES)))
    ax.set_xticklabels(SET_NAMES, rotation=40, ha="right", fontsize=8)
    ax.set_yticklabels(SET_NAMES, fontsize=8)
    for i in range(len(SET_NAMES)):
        for j in range(len(SET_NAMES)):
            ax.text(j, i, str(int(matrix.iat[i, j])), ha="center", va="center", fontsize=7)
    ax.set_title("Pairwise overlapping gene counts")
    fig.colorbar(im, ax=ax, fraction=0.046)
    _save(fig, dest)


def plot_upset(membership: pd.DataFrame, dest: Path, data_dir: Path) -> None:
    from upsetplot import UpSet, from_memberships

    write_csv(membership, data_dir / "seven_set_upset.csv")
    memberships = []
    for _, rec in membership.iterrows():
        names = [name for name in SET_NAMES if rec[name] == "TRUE"]
        memberships.append(names)
    fig = plt.figure(figsize=(10, 6))
    if not memberships:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No genes", ha="center")
        ax.axis("off")
        _save(fig, dest)
        return
    data = from_memberships(memberships)
    if not isinstance(data.index, pd.MultiIndex) or data.index.nlevels < 2:
        ax = fig.add_subplot(111)
        ax.text(
            0.5,
            0.5,
            "UpSet plot needs genes in at least two named sets",
            ha="center",
            va="center",
        )
        ax.axis("off")
        fig.suptitle("Exact seven-set intersections")
        _save(fig, dest)
        return
    UpSet(data, subset_size="count", show_counts=True).plot(fig=fig)
    fig.suptitle("Exact seven-set intersections")
    _save(fig, dest)


def plot_membership_matrix(overlapping: pd.DataFrame, dest: Path, data_dir: Path) -> None:
    write_csv(overlapping, data_dir / "overlapping_gene_membership_matrix.csv")
    if overlapping.empty:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "No overlapping genes", ha="center")
        ax.axis("off")
        _save(fig, dest)
        return
    genes = overlapping.sort_values(["membership_count", "primary_symbol"], ascending=[False, True])
    chunk = 40
    n_pages = max(1, int(np.ceil(len(genes) / chunk)))
    fig, axes = plt.subplots(n_pages, 1, figsize=(10, min(16, 0.28 * chunk * n_pages + 1.5)))
    if n_pages == 1:
        axes = [axes]
    cmap = ListedColormap(["#F2F2F2", "#1F4E79"])
    for page, ax in enumerate(axes):
        part = genes.iloc[page * chunk : (page + 1) * chunk]
        matrix = part[SET_NAMES].replace({"TRUE": 1, "FALSE": 0}).astype(int).values
        ax.imshow(matrix, cmap=cmap, aspect="auto", vmin=0, vmax=1)
        ax.set_xticks(range(len(SET_NAMES)))
        ax.set_xticklabels(SET_NAMES, rotation=40, ha="right", fontsize=7)
        ax.set_yticks(range(len(part)))
        ax.set_yticklabels(part["primary_symbol"].tolist(), fontsize=7)
        ax.set_title(f"Overlapping gene membership ({page + 1}/{n_pages})")
    fig.tight_layout()
    _save(fig, dest)


def plot_all(
    pathway_sets: dict[str, pd.DataFrame],
    approved: pd.DataFrame,
    membership: pd.DataFrame,
    overlapping: pd.DataFrame,
    pairs: pd.DataFrame,
    figures: Path,
    figure_data: Path,
) -> list[Path]:
    figure_data.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    outputs = [
        figures / "pathway_gene_counts.png",
        figures / "go_terms_per_bucket.png",
        figures / "pathway_overlap_venn.png",
        figures / "seven_set_pairwise_counts.png",
        figures / "seven_set_upset.png",
        figures / "overlapping_gene_membership_matrix.png",
    ]
    plot_pathway_gene_counts(pathway_sets, outputs[0], figure_data)
    plot_go_terms(approved, outputs[1], figure_data)
    plot_pathway_venn(pathway_sets, outputs[2], figure_data)
    plot_pairwise_counts(pairs, outputs[3], figure_data)
    plot_upset(membership, outputs[4], figure_data)
    plot_membership_matrix(overlapping, outputs[5], figure_data)
    return outputs
