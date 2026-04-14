#!/usr/bin/env python3
"""
Generate 2D t-SNE visualizations of OpenAI embeddings for biologically
relevant phrases spanning Drosophila behavior, genetics, and sleep/circadian
biology.

Usage:
    python scripts/visualize_embeddings_tsne.py

Requires:
    OPENAI_API_KEY environment variable set.

Outputs:
    docs/images/embedding_tsne_*.png
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
from dotenv import load_dotenv
from matplotlib.lines import Line2D
from openai import OpenAI
from sklearn.manifold import TSNE

load_dotenv()

# ── Phrase corpus ─────────────────────────────────────────────────────────────
# Organized into semantic categories so we can colour-code the plot.

CATEGORIES: dict[str, list[str]] = {
    "Sleep & Circadian": [
        "circadian rhythm",
        "sleep homeostasis",
        "clock gene period",
        "timeless mutant",
        "PDF neuron signaling",
        "rest-activity rhythm",
        "sleep deprivation rebound",
        "daytime sleep",
        "nighttime activity",
        "doubletime kinase",
        "cryptochrome photoreceptor",
        "morning anticipation",
        "evening anticipation",
        "free-running period",
        "light entrainment",
    ],
    "Behavior": [
        "courtship behavior",
        "aggression assay",
        "locomotor activity",
        "grooming behavior",
        "phototaxis response",
        "geotaxis climbing",
        "olfactory avoidance",
        "food preference",
        "mating latency",
        "male wing song",
        "startle response",
        "ethanol sensitivity",
        "learning and memory",
        "social clustering",
        "optogenetic activation",
    ],
    "Genetics & Molecular": [
        "UAS-GAL4 expression system",
        "RNAi knockdown",
        "CRISPR gene editing",
        "balancer chromosome",
        "P-element insertion",
        "transgenic construct",
        "enhancer trap line",
        "loss of function allele",
        "gain of function mutation",
        "homozygous lethal",
        "X-linked inheritance",
        "dosage compensation",
        "chromatin remodeling",
        "transcription factor binding",
        "epigenetic silencing",
    ],
    "Neuroscience & Anatomy": [
        "mushroom body neurons",
        "dopaminergic signaling",
        "GABAergic inhibition",
        "central complex neuropil",
        "ventral nerve cord",
        "antennal lobe glomeruli",
        "projection neuron",
        "kenyon cell calcium imaging",
        "synaptic vesicle release",
        "neuropeptide signaling",
    ],
}

CATEGORY_COLORS = {
    "Sleep & Circadian": "#5B8DB8",
    "Behavior": "#E07B54",
    "Genetics & Molecular": "#7CB872",
    "Neuroscience & Anatomy": "#C279C0",
}

EMBEDDING_MODEL = "text-embedding-3-large"


def get_embeddings(phrases: list[str], client: OpenAI) -> np.ndarray:
    """Fetch embeddings from OpenAI in a single batched call."""
    response = client.embeddings.create(input=phrases, model=EMBEDDING_MODEL)
    vectors = [item.embedding for item in response.data]
    return np.array(vectors, dtype=np.float64)


def build_label_and_color_arrays(
    categories: dict[str, list[str]],
) -> tuple[list[str], list[str], list[str]]:
    """Flatten category dict into parallel lists of phrases, labels, colours."""
    phrases, labels, colors = [], [], []
    for cat, items in categories.items():
        for item in items:
            phrases.append(item)
            labels.append(item)
            colors.append(CATEGORY_COLORS[cat])
    return phrases, labels, colors


def run_tsne(embeddings: np.ndarray, perplexity: float = 12.0) -> np.ndarray:
    tsne = TSNE(
        n_components=2,
        perplexity=perplexity,
        random_state=42,
        init="pca",
        learning_rate="auto",
        max_iter=2000,
    )
    return tsne.fit_transform(embeddings)


def plot_main(
    coords: np.ndarray,
    labels: list[str],
    colors: list[str],
    out_path: Path,
) -> None:
    """Full scatter with text labels for every point."""
    fig, ax = plt.subplots(figsize=(16, 12))
    ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=70, alpha=0.85,
               edgecolors="white", linewidths=0.5, zorder=3)

    texts = []
    for i, label in enumerate(labels):
        texts.append(
            ax.text(
                coords[i, 0], coords[i, 1], label,
                fontsize=7.5, alpha=0.92, ha="center", va="bottom",
            )
        )
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color="grey", lw=0.4))

    legend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=c,
               markersize=9, label=cat)
        for cat, c in CATEGORY_COLORS.items()
    ]
    ax.legend(handles=legend_handles, loc="upper left", fontsize=9,
              framealpha=0.9, title="Category", title_fontsize=10)

    ax.set_title(
        "t-SNE of OpenAI Embeddings — Drosophila Biology Phrases",
        fontsize=14, fontweight="bold", pad=14,
    )
    ax.set_xlabel("t-SNE 1", fontsize=10)
    ax.set_ylabel("t-SNE 2", fontsize=10)
    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.15)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved → {out_path}")


def plot_clusters_only(
    coords: np.ndarray,
    labels: list[str],
    colors: list[str],
    out_path: Path,
) -> None:
    """Cleaner view: colour-coded dots with only a few representative labels."""
    representative_indices = set()
    cat_offset = 0
    for cat, items in CATEGORIES.items():
        n = len(items)
        sub = coords[cat_offset : cat_offset + n]
        centroid = sub.mean(axis=0)
        dists = np.linalg.norm(sub - centroid, axis=1)
        representative_indices.add(cat_offset + int(np.argmin(dists)))
        representative_indices.add(cat_offset + int(np.argmax(dists)))
        cat_offset += n

    fig, ax = plt.subplots(figsize=(14, 10))
    ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=60, alpha=0.8,
               edgecolors="white", linewidths=0.5, zorder=3)

    texts = []
    for i in representative_indices:
        texts.append(
            ax.text(
                coords[i, 0], coords[i, 1], labels[i],
                fontsize=8, fontweight="bold", alpha=0.95,
                ha="center", va="bottom",
            )
        )
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))

    legend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=c,
               markersize=9, label=cat)
        for cat, c in CATEGORY_COLORS.items()
    ]
    ax.legend(handles=legend_handles, loc="upper left", fontsize=9,
              framealpha=0.9, title="Category", title_fontsize=10)

    ax.set_title(
        "t-SNE Clusters — Drosophila Biology Embeddings",
        fontsize=14, fontweight="bold", pad=14,
    )
    ax.set_xlabel("t-SNE 1", fontsize=10)
    ax.set_ylabel("t-SNE 2", fontsize=10)
    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.15)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved → {out_path}")


def plot_distance_heatmap(
    embeddings: np.ndarray,
    labels: list[str],
    out_path: Path,
) -> None:
    """Cosine-similarity heatmap across all phrases."""
    norms = embeddings / np.linalg.norm(embeddings, axis=1, keepdims=True)
    sim = norms @ norms.T

    fig, ax = plt.subplots(figsize=(18, 15))
    im = ax.imshow(sim, cmap="RdYlBu_r", vmin=0.2, vmax=1.0, aspect="auto")
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=5.5)
    ax.set_yticklabels(labels, fontsize=5.5)

    cat_boundaries, offset = [], 0
    for cat, items in CATEGORIES.items():
        mid = offset + len(items) / 2
        cat_boundaries.append((offset, mid, cat))
        offset += len(items)

    for start, _, _ in cat_boundaries[1:]:
        ax.axhline(start - 0.5, color="black", linewidth=0.8)
        ax.axvline(start - 0.5, color="black", linewidth=0.8)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Cosine Similarity", fontsize=10)
    ax.set_title(
        "Pairwise Cosine Similarity — OpenAI Embeddings",
        fontsize=14, fontweight="bold", pad=14,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved → {out_path}")


def main() -> None:
    if not os.environ.get("OPENAI_API_KEY"):
        print("ERROR: Set OPENAI_API_KEY before running.", file=sys.stderr)
        sys.exit(1)

    out_dir = Path(__file__).resolve().parent.parent / "docs" / "images"
    out_dir.mkdir(parents=True, exist_ok=True)

    phrases, labels, colors = build_label_and_color_arrays(CATEGORIES)
    print(f"Embedding {len(phrases)} phrases with {EMBEDDING_MODEL} …")

    client = OpenAI()
    embeddings = get_embeddings(phrases, client)
    print(f"  shape: {embeddings.shape}")

    print("Running t-SNE …")
    coords = run_tsne(embeddings)

    print("Generating visualizations …")
    plot_main(coords, labels, colors, out_dir / "embedding_tsne_labeled.png")
    plot_clusters_only(coords, labels, colors, out_dir / "embedding_tsne_clusters.png")
    plot_distance_heatmap(embeddings, labels, out_dir / "embedding_cosine_heatmap.png")

    print("Done.")


if __name__ == "__main__":
    main()
