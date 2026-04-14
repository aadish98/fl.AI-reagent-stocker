#!/usr/bin/env python3
"""Generate a small OpenAI embedding demo for Drosophila concepts."""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import matplotlib
import numpy as np
import pandas as pd
from adjustText import adjust_text
from dotenv import load_dotenv
from openai import OpenAI
from sklearn.manifold import TSNE
from sklearn.metrics.pairwise import cosine_similarity

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


@dataclass(frozen=True)
class EmbeddingExample:
    category: str
    text: str


DEFAULT_MODEL = "text-embedding-3-large"
DEFAULT_OUTPUT_DIR = Path("docs/images")
DEFAULT_PREFIX = "openai_embeddings_drosophila"
COLOR_BY_CATEGORY = {
    "Behavior": "#2563EB",
    "Neural": "#7C3AED",
    "Genetics": "#DB2777",
    "Development": "#EA580C",
}
EXAMPLES: Sequence[EmbeddingExample] = (
    EmbeddingExample("Behavior", "courtship song"),
    EmbeddingExample("Behavior", "sleep bout"),
    EmbeddingExample("Behavior", "negative geotaxis"),
    EmbeddingExample("Behavior", "olfactory learning"),
    EmbeddingExample("Neural", "mushroom body neuron"),
    EmbeddingExample("Neural", "dopamine signaling"),
    EmbeddingExample("Neural", "circadian pacemaker"),
    EmbeddingExample("Neural", "synaptic vesicle release"),
    EmbeddingExample("Genetics", "fruitless mutant"),
    EmbeddingExample("Genetics", "doublesex regulation"),
    EmbeddingExample("Genetics", "GAL4/UAS driver"),
    EmbeddingExample("Genetics", "RNAi knockdown"),
    EmbeddingExample("Development", "wing imaginal disc"),
    EmbeddingExample("Development", "Notch pathway"),
    EmbeddingExample("Development", "axon guidance"),
    EmbeddingExample("Development", "photoreceptor differentiation"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a t-SNE embedding demo, similarity heatmap, and summary table "
            "for short Drosophila behavior and genetics phrases."
        )
    )
    parser.add_argument(
        "--model",
        default=DEFAULT_MODEL,
        help=f"OpenAI embedding model to use (default: {DEFAULT_MODEL}).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for generated images (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--prefix",
        default=DEFAULT_PREFIX,
        help=f"Output file prefix (default: {DEFAULT_PREFIX}).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=17,
        help="Random seed for the t-SNE projection.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=320,
        help="Output resolution for PNGs.",
    )
    return parser.parse_args()


def get_client() -> OpenAI:
    load_dotenv()
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY is not set. Load it in the environment or .env.")
    return OpenAI(api_key=api_key)


def fetch_embeddings(
    client: OpenAI,
    examples: Sequence[EmbeddingExample],
    model: str,
) -> np.ndarray:
    texts = [example.text for example in examples]
    response = client.embeddings.create(model=model, input=texts)
    rows = [[float(value) for value in item.embedding] for item in response.data]
    return np.asarray(rows, dtype=float)


def project_embeddings(embeddings: np.ndarray, seed: int) -> np.ndarray:
    perplexity = min(6, len(embeddings) - 1)
    tsne = TSNE(
        n_components=2,
        perplexity=perplexity,
        init="pca",
        learning_rate="auto",
        random_state=seed,
        max_iter=2000,
    )
    return tsne.fit_transform(embeddings)


def build_summary_frame(
    examples: Sequence[EmbeddingExample],
    embeddings: np.ndarray,
    coordinates: np.ndarray,
) -> pd.DataFrame:
    similarity = cosine_similarity(embeddings)
    labels = [example.text for example in examples]
    nearest_neighbor = []
    nearest_score = []

    for row_idx in range(len(labels)):
        ranked = similarity[row_idx].copy()
        ranked[row_idx] = -1.0
        neighbor_idx = int(np.argmax(ranked))
        nearest_neighbor.append(labels[neighbor_idx])
        nearest_score.append(float(similarity[row_idx, neighbor_idx]))

    frame = pd.DataFrame(
        {
            "source_index": list(range(len(examples))),
            "category": [example.category for example in examples],
            "text": labels,
            "nearest_neighbor": nearest_neighbor,
            "nearest_neighbor_cosine": nearest_score,
            "tsne_x": coordinates[:, 0],
            "tsne_y": coordinates[:, 1],
        }
    )
    return frame


def save_tsne_plot(frame: pd.DataFrame, output_path: Path, model: str, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(12.5, 9.0))
    texts = []

    for category in COLOR_BY_CATEGORY:
        subset = frame[frame["category"] == category]
        ax.scatter(
            subset["tsne_x"],
            subset["tsne_y"],
            s=160,
            color=COLOR_BY_CATEGORY[category],
            label=category,
            alpha=0.88,
            edgecolors="white",
            linewidths=1.3,
        )
        for _, row in subset.iterrows():
            texts.append(
                ax.text(
                    row["tsne_x"],
                    row["tsne_y"],
                    row["text"],
                    fontsize=10,
                    color="#0F172A",
                )
            )

    adjust_text(
        texts,
        ax=ax,
        arrowprops={"arrowstyle": "-", "color": "#94A3B8", "lw": 0.8},
    )
    ax.set_title("OpenAI embeddings projected with t-SNE", fontsize=18, pad=16)
    ax.set_xlabel("t-SNE 1")
    ax.set_ylabel("t-SNE 2")
    ax.grid(alpha=0.18, linestyle="--", linewidth=0.6)
    ax.legend(title="Concept family", frameon=True)
    fig.text(
        0.5,
        0.02,
        f"Examples: Drosophila behavior and genetics phrases | Model: {model}",
        ha="center",
        fontsize=10,
        color="#475569",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def save_similarity_heatmap(
    frame: pd.DataFrame,
    embeddings: np.ndarray,
    output_path: Path,
    model: str,
    dpi: int,
) -> None:
    similarity = cosine_similarity(embeddings)
    ordered_frame = frame.sort_values(["category", "text"], kind="stable")
    order = ordered_frame["source_index"].tolist()
    ordered_similarity = similarity[np.ix_(order, order)]
    ordered_labels = ordered_frame["text"].tolist()

    fig, ax = plt.subplots(figsize=(12.0, 10.0))
    image = ax.imshow(ordered_similarity, cmap="magma", vmin=0.15, vmax=1.0)
    ax.set_xticks(range(len(ordered_labels)))
    ax.set_yticks(range(len(ordered_labels)))
    ax.set_xticklabels(ordered_labels, rotation=60, ha="right", fontsize=8)
    ax.set_yticklabels(ordered_labels, fontsize=8)
    ax.set_title("Cosine similarity between embedding examples", fontsize=18, pad=16)
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="Cosine similarity")
    fig.text(
        0.5,
        0.02,
        f"Rows and columns are sorted by concept family | Model: {model}",
        ha="center",
        fontsize=10,
        color="#475569",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def save_summary_table(frame: pd.DataFrame, output_path: Path, dpi: int) -> None:
    table_frame = frame.drop(columns=["source_index"]).copy()
    table_frame["nearest_neighbor_cosine"] = table_frame["nearest_neighbor_cosine"].map(
        lambda value: f"{value:.3f}"
    )
    table_frame["tsne_x"] = table_frame["tsne_x"].map(lambda value: f"{value:.2f}")
    table_frame["tsne_y"] = table_frame["tsne_y"].map(lambda value: f"{value:.2f}")
    table_frame = table_frame.rename(
        columns={
            "category": "Category",
            "text": "Phrase",
            "nearest_neighbor": "Nearest phrase",
            "nearest_neighbor_cosine": "Cosine",
            "tsne_x": "t-SNE x",
            "tsne_y": "t-SNE y",
        }
    )

    fig_height = max(6.5, 0.45 * len(table_frame) + 1.8)
    fig, ax = plt.subplots(figsize=(14.0, fig_height))
    ax.axis("off")
    ax.set_title("Embedding example summary table", fontsize=18, pad=16)

    table = ax.table(
        cellText=table_frame.values,
        colLabels=table_frame.columns,
        loc="center",
        cellLoc="left",
        colLoc="left",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1, 1.35)

    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#CBD5E1")
        cell.set_linewidth(0.6)
        if row == 0:
            cell.set_facecolor("#E2E8F0")
            cell.set_text_props(weight="bold", color="#0F172A")
        else:
            cell.set_facecolor("#F8FAFC" if row % 2 else "#FFFFFF")

    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    client = get_client()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    embeddings = fetch_embeddings(client, EXAMPLES, args.model)
    coordinates = project_embeddings(embeddings, args.seed)
    frame = build_summary_frame(EXAMPLES, embeddings, coordinates)
    frame = frame.sort_values(["category", "text"], kind="stable").reset_index(drop=True)

    csv_path = args.output_dir / f"{args.prefix}_table.csv"
    tsne_path = args.output_dir / f"{args.prefix}_tsne.png"
    heatmap_path = args.output_dir / f"{args.prefix}_similarity_heatmap.png"
    table_path = args.output_dir / f"{args.prefix}_table.png"

    frame.drop(columns=["source_index"]).to_csv(csv_path, index=False)
    save_tsne_plot(frame, tsne_path, args.model, args.dpi)
    save_similarity_heatmap(frame, embeddings, heatmap_path, args.model, args.dpi)
    save_summary_table(frame, table_path, args.dpi)

    print(f"Saved table CSV: {csv_path}")
    print(f"Saved t-SNE plot: {tsne_path}")
    print(f"Saved heatmap: {heatmap_path}")
    print(f"Saved table image: {table_path}")


if __name__ == "__main__":
    main()
