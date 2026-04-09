"""
Phenotype similarity scoring and visualization helpers.

This module provides OpenAI embedding cosine similarity scoring for FlyBase
phenotype terms used in the soft-run Stock Phenotype Sheet.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patheffects as patheffects  # noqa: E402
from adjustText import adjust_text  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize  # noqa: E402

_RING_PALETTE = ["#D6CEC8", "#B8AFA8", "#8C8278", "#645C54"]

_TSNE_CMAP_EMBEDDING = LinearSegmentedColormap.from_list(
    "tsne_embedding",
    ["#5B7C99", "#7C3AED", "#DB2777", "#B91C1C"],
)

_TSNE_FIGSIZE = (16.0, 11.5)
_TSNE_EXPORT_DPI = 360
_TSNE_MIN_FONT = 9.5
_TSNE_MAX_FONT = 17.5

_PUB_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "axes.linewidth": 0.8,
    "axes.edgecolor": "#334155",
    "axes.labelcolor": "#1E293B",
    "axes.titlepad": 14,
    "xtick.color": "#475569",
    "ytick.color": "#475569",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "figure.facecolor": "white",
    "axes.facecolor": "#FAFBFC",
    "axes.grid": False,
    "legend.frameon": True,
    "legend.framealpha": 0.85,
    "legend.edgecolor": "#CBD5E1",
    "legend.fontsize": 9,
}

try:
    from openai import OpenAI

    OPENAI_AVAILABLE = True
except ImportError:
    OPENAI_AVAILABLE = False


def sanitize_filename_fragment(value: str) -> str:
    """Return a filesystem-safe slug for target-specific filenames."""
    slug = re.sub(r"[^a-zA-Z0-9]+", "_", str(value or "").strip().lower())
    return slug.strip("_") or "target"


def round_score(value: Optional[float]) -> Optional[float]:
    """Round similarity values for stable workbook output."""
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    return round(float(value), 6)


def cosine_similarity(vector_a: Sequence[float], vector_b: Sequence[float]) -> float:
    """Compute cosine similarity for two dense vectors."""
    a = np.asarray(vector_a, dtype=float)
    b = np.asarray(vector_b, dtype=float)
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom == 0:
        return 0.0
    return float(np.dot(a, b) / denom)


def normalize_phenotype_text(
    phenotype_name: str,
) -> str:
    """Normalize the phenotype string for embedding input."""
    return str(phenotype_name or "").strip()


def normalize_qualifier_text(qualifier_names: str) -> str:
    """Normalize FlyBase qualifier names for separate sheet display."""
    qualifiers = [q.strip() for q in str(qualifier_names or "").split("|") if q.strip()]
    return ", ".join(qualifiers)


@dataclass(frozen=True)
class PhenotypeSimilarityTarget:
    keyword: str
    embedding_text: str

    @property
    def slug(self) -> str:
        return sanitize_filename_fragment(self.keyword)

    @property
    def cosine_similarity_column(self) -> str:
        return f"Cosine Similarity ({self.keyword})"


def build_similarity_targets(raw_targets: Sequence[Dict[str, str]]) -> List[PhenotypeSimilarityTarget]:
    """Convert normalized config dictionaries into typed targets."""
    return [
        PhenotypeSimilarityTarget(
            keyword=str(target["keyword"]),
            embedding_text=str(target["embedding_text"]),
        )
        for target in raw_targets
    ]


class TextEmbeddingCache:
    """CSV-backed cache for text embeddings."""

    CACHE_COLUMNS = ["model", "text", "embedding_json", "updated_at"]

    def __init__(self, cache_path: Optional[Path]):
        self.cache_path = Path(cache_path) if cache_path else None
        self._cache: Dict[Tuple[str, str], List[float]] = {}
        self._loaded = False

    def load(self) -> Dict[Tuple[str, str], List[float]]:
        if self._loaded:
            return self._cache
        if self.cache_path and self.cache_path.exists():
            try:
                cache_df = pd.read_csv(self.cache_path, dtype=str, keep_default_na=False)
                for _, row in cache_df.iterrows():
                    model = str(row.get("model", "") or "").strip()
                    text = str(row.get("text", "") or "").strip()
                    payload = str(row.get("embedding_json", "") or "").strip()
                    if not model or not text or not payload:
                        continue
                    try:
                        embedding = json.loads(payload)
                    except json.JSONDecodeError:
                        continue
                    if isinstance(embedding, list) and embedding:
                        self._cache[(model, text)] = [float(v) for v in embedding]
            except Exception as exc:
                print(f"    Warning: Could not load embedding cache {self.cache_path}: {exc}")
        self._loaded = True
        return self._cache

    def get(self, model: str, text: str) -> Optional[List[float]]:
        self.load()
        return self._cache.get((str(model).strip(), str(text).strip()))

    def set(self, model: str, text: str, embedding: Sequence[float]) -> None:
        self.load()
        model_key = str(model).strip()
        text_key = str(text).strip()
        if not model_key or not text_key:
            return
        self._cache[(model_key, text_key)] = [float(v) for v in embedding]

    def save(self) -> None:
        if not self.cache_path:
            return
        rows = []
        for (model, text), embedding in sorted(self._cache.items()):
            rows.append(
                {
                    "model": model,
                    "text": text,
                    "embedding_json": json.dumps(embedding),
                    "updated_at": datetime.now(timezone.utc).isoformat(),
                }
            )
        out_df = pd.DataFrame(rows, columns=self.CACHE_COLUMNS)
        try:
            self.cache_path.parent.mkdir(parents=True, exist_ok=True)
            out_df.to_csv(self.cache_path, index=False)
        except OSError as exc:
            print(f"    Warning: Could not persist embedding cache {self.cache_path}: {exc}")


class EmbeddingSimilarityScorer:
    """OpenAI embedding scorer with local CSV caches."""

    def __init__(
        self,
        openai_api_key: Optional[str],
        model: str,
        phenotype_cache_path: Optional[Path],
        target_cache_path: Optional[Path],
    ):
        self.model = model
        self._phenotype_cache = TextEmbeddingCache(phenotype_cache_path)
        self._target_cache = TextEmbeddingCache(target_cache_path)
        self._client = None
        self.is_available = bool(OPENAI_AVAILABLE and openai_api_key)
        if self.is_available:
            self._client = OpenAI(api_key=openai_api_key)

    def _ensure_embeddings(
        self,
        texts: Sequence[str],
        cache: TextEmbeddingCache,
    ) -> Dict[str, List[float]]:
        cache.load()
        resolved: Dict[str, List[float]] = {}
        missing: List[str] = []
        for text in texts:
            normalized = str(text or "").strip()
            if not normalized:
                continue
            cached = cache.get(self.model, normalized)
            if cached is not None:
                resolved[normalized] = cached
            else:
                missing.append(normalized)

        if missing and self.is_available and self._client is not None:
            try:
                for batch_start in range(0, len(missing), 100):
                    batch = missing[batch_start: batch_start + 100]
                    response = self._client.embeddings.create(model=self.model, input=batch)
                    for item, text in zip(response.data, batch):
                        embedding = [float(v) for v in item.embedding]
                        cache.set(self.model, text, embedding)
                        resolved[text] = embedding
                cache.save()
            except Exception:
                self.is_available = False
                print(
                    "    Warning: OpenAI embedding similarity disabled after embedding "
                    "request failed. Check OPENAI_API_KEY and embedding model settings."
                )

        for text in missing:
            cached = cache.get(self.model, text)
            if cached is not None:
                resolved[text] = cached
        return resolved

    def score_texts(
        self,
        phenotype_texts: Iterable[str],
        targets: Sequence[PhenotypeSimilarityTarget],
    ) -> Dict[str, Dict[str, Optional[float]]]:
        unique_texts = sorted({str(text or "").strip() for text in phenotype_texts if str(text or "").strip()})
        if not unique_texts:
            return {}

        empty_result = {
            text: {target.cosine_similarity_column: None for target in targets}
            for text in unique_texts
        }
        if not self.is_available:
            return empty_result

        phenotype_embeddings = self._ensure_embeddings(unique_texts, self._phenotype_cache)
        target_embeddings = self._ensure_embeddings(
            [target.embedding_text for target in targets],
            self._target_cache,
        )

        results: Dict[str, Dict[str, Optional[float]]] = {}
        for text in unique_texts:
            row_scores: Dict[str, Optional[float]] = {}
            text_embedding = phenotype_embeddings.get(text)
            for target in targets:
                target_embedding = target_embeddings.get(target.embedding_text)
                if text_embedding is None or target_embedding is None:
                    row_scores[target.cosine_similarity_column] = None
                else:
                    row_scores[target.cosine_similarity_column] = round_score(
                        cosine_similarity(text_embedding, target_embedding)
                    )
            results[text] = row_scores
        return results

    def tsne_frame_for_target(
        self,
        target: PhenotypeSimilarityTarget,
        phenotype_df: pd.DataFrame,
    ) -> pd.DataFrame:
        if not self.is_available:
            return pd.DataFrame()

        grouped = (
            phenotype_df.groupby(["Phenotype"], as_index=False)
            .size()
            .rename(columns={"size": "frequency"})
        )
        if grouped.empty:
            return pd.DataFrame()

        target_row = pd.DataFrame(
            [{"Phenotype": target.embedding_text, "frequency": max(1, int(grouped["frequency"].max()))}]
        )
        grouped = pd.concat([target_row, grouped], ignore_index=True)
        grouped = grouped.drop_duplicates(subset=["Phenotype"], keep="first")
        texts = grouped["Phenotype"].tolist()
        embeddings = self._ensure_embeddings(texts, self._phenotype_cache)
        if len(embeddings) < 1:
            return pd.DataFrame()

        valid_rows = grouped[grouped["Phenotype"].isin(embeddings.keys())].copy()
        if valid_rows.empty:
            return pd.DataFrame()
        valid_rows["is_target"] = valid_rows["Phenotype"].eq(target.embedding_text)
        target_embedding = embeddings.get(target.embedding_text)
        valid_rows[target.cosine_similarity_column] = valid_rows["Phenotype"].map(
            lambda text: round_score(
                None
                if target_embedding is None
                else cosine_similarity(embeddings[text], target_embedding)
            )
        )
        valid_rows = _select_tsne_rows(
            valid_rows,
            target_mask=valid_rows["is_target"].fillna(False),
            closest_sort_cols=[target.cosine_similarity_column, "frequency", "Phenotype"],
            closest_ascending=[False, False, True],
            dedupe_cols=["Phenotype"],
            n_closest=20,
            n_frequent=10,
        )
        matrix = np.asarray([embeddings[text] for text in valid_rows["Phenotype"]], dtype=float)
        coords = compute_tsne_or_fallback(matrix, metric="cosine")
        valid_rows["tsne_x"] = coords[:, 0]
        valid_rows["tsne_y"] = coords[:, 1]
        valid_rows = _center_tsne_on_target(valid_rows)
        return valid_rows


def compute_tsne_or_fallback(data: np.ndarray, metric: str) -> np.ndarray:
    """Compute a stable 2D layout, falling back for very small inputs."""
    n_samples = len(data)
    if n_samples == 0:
        return np.zeros((0, 2), dtype=float)
    if n_samples == 1:
        return np.asarray([[0.0, 0.0]])
    if n_samples == 2:
        return np.asarray([[-1.0, 0.0], [1.0, 0.0]])

    perplexity = min(30.0, max(2.0, float(n_samples - 1) / 3.0))
    perplexity = min(perplexity, float(n_samples - 1) - 1e-6)

    if metric == "precomputed":
        tsne = TSNE(
            n_components=2,
            metric="precomputed",
            init="random",
            learning_rate="auto",
            perplexity=perplexity,
            random_state=42,
            max_iter=1000,
        )
    else:
        tsne = TSNE(
            n_components=2,
            metric=metric,
            init="random",
            learning_rate="auto",
            perplexity=perplexity,
            random_state=42,
            max_iter=1000,
        )
    return tsne.fit_transform(data)


def _select_tsne_rows(
    df: pd.DataFrame,
    *,
    target_mask: pd.Series,
    closest_sort_cols: List[str],
    closest_ascending: List[bool],
    dedupe_cols: List[str],
    n_closest: int = 20,
    n_frequent: int = 10,
) -> pd.DataFrame:
    """Keep target plus the union of nearest and most frequent rows."""
    if df.empty:
        return df

    target_df = df[target_mask].copy()
    candidates = df[~target_mask].copy()

    closest_df = candidates.sort_values(
        by=closest_sort_cols,
        ascending=closest_ascending,
        na_position="last",
    ).head(n_closest)
    frequent_df = candidates.sort_values(
        by=["frequency", *dedupe_cols],
        ascending=[False, *([True] * len(dedupe_cols))],
        na_position="last",
    ).head(n_frequent)

    selected = pd.concat([target_df, closest_df, frequent_df], ignore_index=True)
    return selected.drop_duplicates(subset=dedupe_cols, keep="first").reset_index(drop=True)


def _center_tsne_on_target(
    df: pd.DataFrame,
    *,
    target_mask_col: str = "is_target",
) -> pd.DataFrame:
    """Shift coordinates so the target point lands at the origin."""
    if df.empty or target_mask_col not in df.columns:
        return df

    target_rows = df[df[target_mask_col].fillna(False)]
    if target_rows.empty:
        return df

    centered = df.copy()
    target_x = float(target_rows["tsne_x"].iloc[0])
    target_y = float(target_rows["tsne_y"].iloc[0])
    centered["tsne_x"] = centered["tsne_x"] - target_x
    centered["tsne_y"] = centered["tsne_y"] - target_y
    return centered


def _style_axes(ax, *, grid: bool = False, light_spines: bool = False) -> None:
    """Apply clean publication styling to an axes instance."""
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    if light_spines:
        for spine in ("bottom", "left"):
            ax.spines[spine].set_color("#CBD5E1")
            ax.spines[spine].set_linewidth(0.6)
    ax.tick_params(labelsize=8, colors="#475569")
    if grid:
        ax.grid(True, linewidth=0.4, color="#E2E8F0", alpha=0.7, zorder=0)
        ax.set_axisbelow(True)


def _set_centered_axes(ax, tsne_df: pd.DataFrame, *, padding: float = 1.35) -> None:
    """Use symmetric limits around zero so the target stays centered visually."""
    if tsne_df.empty:
        return
    max_abs_x = max(1.0, float(np.abs(tsne_df["tsne_x"]).max()))
    max_abs_y = max(1.0, float(np.abs(tsne_df["tsne_y"]).max()))
    ax.set_xlim(-max_abs_x * padding, max_abs_x * padding)
    ax.set_ylim(-max_abs_y * padding, max_abs_y * padding)
    ax.margins(x=0.03, y=0.05)


def _draw_similarity_rings(
    ax,
    tsne_df: pd.DataFrame,
    value_col: str,
    thresholds: Sequence[float] = (0.2, 0.4, 0.6, 0.8),
) -> None:
    """Draw concentric rings around the target (origin) marking similarity zones.

    Each ring's radius is determined by the farthest point whose similarity
    score >= that threshold, so the rings reflect actual t-SNE layout rather
    than arbitrary spacing.
    """
    scores = pd.to_numeric(tsne_df[value_col], errors="coerce").fillna(0.0).to_numpy()
    is_target = tsne_df["is_target"].fillna(False).to_numpy()
    xs = tsne_df["tsne_x"].to_numpy(dtype=float)
    ys = tsne_df["tsne_y"].to_numpy(dtype=float)
    dists = np.sqrt(xs ** 2 + ys ** 2)
    non_target = ~is_target

    ring_radii: Dict[float, float] = {}
    for thresh in thresholds:
        mask = (scores >= thresh) & non_target
        if mask.any():
            ring_radii[thresh] = float(dists[mask].max()) * 1.15

    if not ring_radii:
        max_dist = float(dists[non_target].max()) if non_target.any() else 1.0
        for thresh in thresholds:
            ring_radii[thresh] = max_dist * (1.0 - thresh) * 1.2 + 0.5

    sorted_thresh = sorted(ring_radii.keys())
    for i in range(len(sorted_thresh) - 1):
        lo, hi = sorted_thresh[i], sorted_thresh[i + 1]
        if ring_radii[lo] < ring_radii[hi]:
            ring_radii[lo] = ring_radii[hi] * 1.15

    n_rings = len(ring_radii)
    for idx, thresh in enumerate(sorted(ring_radii.keys())):
        radius = ring_radii[thresh]
        pal_idx = min(idx, len(_RING_PALETTE) - 1)
        edge_color = _RING_PALETTE[pal_idx]
        fill_alpha = 0.03 + 0.02 * idx
        circle = plt.Circle(
            (0, 0), radius,
            fill=True,
            facecolor=(*matplotlib.colors.to_rgb(edge_color), fill_alpha),
            edgecolor=(*matplotlib.colors.to_rgb(edge_color), 0.48),
            linewidth=1.0,
            linestyle=(0, (6, 4)),
            zorder=0,
        )
        ax.add_patch(circle)
        label_angle = math.pi * (0.28 + 0.12 * idx)
        lx = radius * math.cos(label_angle)
        ly = radius * math.sin(label_angle)
        ax.annotate(
            f"\u2265 {thresh:.1f}",
            xy=(lx, ly),
            fontsize=8.5,
            fontweight="semibold",
            color="#5B514B",
            ha="left",
            va="bottom",
            zorder=1,
            bbox=dict(
                boxstyle="round,pad=0.18",
                facecolor="white",
                edgecolor="#CBD5E1",
                linewidth=0.6,
                alpha=0.92,
            ),
        )


def _plot_tsne_text_map(
    ax,
    tsne_df: pd.DataFrame,
    value_col: str,
    cmap: LinearSegmentedColormap,
    target_keyword: str,
) -> Optional[plt.cm.ScalarMappable]:
    """Plot t-SNE where colored/sized phenotype text replaces scatter dots.

    Returns the ScalarMappable for a colorbar, or None if no data.
    """
    if tsne_df.empty or value_col not in tsne_df.columns:
        return None

    scores = pd.to_numeric(tsne_df[value_col], errors="coerce").fillna(0.0).to_numpy()
    freqs = tsne_df["frequency"].astype(float).to_numpy()
    is_target = tsne_df["is_target"].fillna(False).to_numpy()
    labels = tsne_df["Phenotype"].astype(str).to_numpy()
    xs = tsne_df["tsne_x"].to_numpy(dtype=float)
    ys = tsne_df["tsne_y"].to_numpy(dtype=float)

    norm = Normalize(vmin=0.0, vmax=max(float(scores.max()), 0.01))

    _draw_similarity_rings(ax, tsne_df, value_col)

    min_fontsize, max_fontsize = _TSNE_MIN_FONT, _TSNE_MAX_FONT
    freq_min, freq_max = float(freqs.min()), max(float(freqs.max()), 1.0)
    freq_range = freq_max - freq_min if freq_max > freq_min else 1.0

    ax.scatter(xs, ys, s=18, c="#334155", zorder=2, alpha=0.24, linewidths=0)

    texts = []
    for i in range(len(tsne_df)):
        frac = (freqs[i] - freq_min) / freq_range
        fontsize = min_fontsize + frac * (max_fontsize - min_fontsize)
        color = cmap(norm(scores[i]))

        if is_target[i]:
            txt = ax.text(
                xs[i], ys[i], labels[i],
                fontsize=max_fontsize + 2.0,
                fontweight="bold",
                color="white",
                ha="center", va="center",
                zorder=5,
                bbox=dict(
                    boxstyle="round,pad=0.36",
                    facecolor="#1F2937",
                    edgecolor="#111827",
                    linewidth=1.2,
                    alpha=0.96,
                ),
            )
            txt.set_path_effects(
                [patheffects.withStroke(linewidth=1.2, foreground="#111827", alpha=0.85)]
            )
        else:
            txt = ax.text(
                xs[i], ys[i], labels[i],
                fontsize=fontsize,
                fontweight="semibold" if scores[i] >= 0.35 else "normal",
                color=color,
                ha="center", va="center",
                zorder=4,
            )
            txt.set_path_effects(
                [patheffects.withStroke(linewidth=2.4, foreground="white", alpha=0.95)]
            )
        texts.append(txt)

    adjust_text(
        texts,
        x=xs, y=ys,
        expand=(1.7, 2.0),
        force_text=(0.9, 1.2),
        force_points=(0.6, 0.8),
        arrowprops=dict(
            arrowstyle="-",
            color="#94A3B8",
            alpha=0.45,
            lw=0.8,
            shrinkA=8,
            shrinkB=4,
        ),
        ax=ax,
        only_move={"text": "xy"},
        ensure_inside_axes=True,
    )

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    return sm


def _plot_no_data(ax, title: str, message: str) -> None:
    ax.set_title(title, fontsize=11, fontweight="semibold", color="#1E293B")
    ax.text(
        0.5, 0.5, message,
        ha="center", va="center", fontsize=10, color="#64748B",
        transform=ax.transAxes,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_similarity_outputs(
    phenotype_sheet_df: pd.DataFrame,
    targets: Sequence[PhenotypeSimilarityTarget],
    output_dir: Path,
    embedding_scorer: Optional[EmbeddingSimilarityScorer],
) -> List[Path]:
    """Generate distribution and t-SNE PNGs from the final phenotype sheet."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    written_paths: List[Path] = []

    plot_df = phenotype_sheet_df.copy()
    if plot_df.empty:
        return written_paths

    _prev_rc = {k: matplotlib.rcParams.get(k) for k in _PUB_RC}
    matplotlib.rcParams.update(_PUB_RC)

    try:
        written_paths = _generate_all_plots(
            plot_df, targets, output_dir, embedding_scorer,
        )
    finally:
        for k, v in _prev_rc.items():
            if v is not None:
                matplotlib.rcParams[k] = v

    return written_paths


def _generate_all_plots(
    plot_df: pd.DataFrame,
    targets: Sequence[PhenotypeSimilarityTarget],
    output_dir: Path,
    embedding_scorer: Optional["EmbeddingSimilarityScorer"],
) -> List[Path]:
    written_paths: List[Path] = []

    plot_df["Phenotype"] = plot_df.get("Phenotype", "").fillna("").astype(str)

    for target in targets:
        cosine_col = target.cosine_similarity_column

        fig, ax = plt.subplots(figsize=(8, 5))
        cosine_values = pd.to_numeric(plot_df.get(cosine_col, pd.Series(dtype=float)), errors="coerce").dropna().to_numpy()
        bins = np.linspace(0.0, 1.0, 21)
        if len(cosine_values) > 0:
            ax.hist(
                cosine_values, bins=bins, density=True, alpha=0.55,
                label=cosine_col, color="#0D9488", edgecolor="white", linewidth=0.5,
            )
            ax.set_title(
                f"Similarity score density\nTarget: {target.keyword}",
                fontsize=11, fontweight="semibold", color="#1E293B",
            )
            ax.set_xlabel("Similarity score", fontsize=9, color="#475569")
            ax.set_ylabel("Density", fontsize=9, color="#475569")
            ax.legend(fontsize=8, loc="upper right")
            _style_axes(ax, grid=True)
        else:
            _plot_no_data(
                ax,
                f"Similarity score density: {target.keyword}",
                "No similarity values available for this target.",
            )
        path = output_dir / f"{target.slug}_similarity_density.png"
        fig.tight_layout(pad=1.2)
        fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        written_paths.append(path)

        if embedding_scorer is not None and embedding_scorer.is_available:
            tsne_df = embedding_scorer.tsne_frame_for_target(target, plot_df)
            fig, ax = plt.subplots(figsize=_TSNE_FIGSIZE)
            if tsne_df.empty:
                _plot_no_data(
                    ax,
                    f"Cosine t-SNE: {target.keyword}",
                    "No embedding vectors available for this target.",
                )
            else:
                _set_centered_axes(ax, tsne_df, padding=1.45)
                sm = _plot_tsne_text_map(ax, tsne_df, cosine_col, _TSNE_CMAP_EMBEDDING, target.keyword)
                if sm is not None:
                    cbar = fig.colorbar(sm, ax=ax, shrink=0.82, pad=0.018, fraction=0.05)
                    cbar.set_label(cosine_col, fontsize=10.5, color="#475569")
                    cbar.ax.tick_params(labelsize=9, colors="#475569")
                    cbar.outline.set_edgecolor("#CBD5E1")
                    cbar.outline.set_linewidth(0.8)
                ax.set_title(
                    f"Cosine t-SNE (embedding space)\nTarget: {target.keyword}",
                    fontsize=14, fontweight="semibold", color="#1E293B",
                )
                ax.set_xlabel("t-SNE 1", fontsize=10.5, color="#64748B")
                ax.set_ylabel("t-SNE 2", fontsize=10.5, color="#64748B")
                _style_axes(ax, light_spines=True)
            path = output_dir / f"tsne_cosine_{target.slug}.png"
            fig.tight_layout(pad=2.2)
            fig.savefig(path, dpi=_TSNE_EXPORT_DPI, bbox_inches="tight", facecolor="white")
            plt.close(fig)
            written_paths.append(path)

    fig, axes = plt.subplots(
        len(targets), 1,
        figsize=(8, 4.5 * max(1, len(targets))),
        squeeze=False,
    )
    for row_idx, target in enumerate(targets):
        target_df = (
            plot_df.groupby(["Phenotype"], as_index=False)
            .agg(
                frequency=("Phenotype", "size"),
                cosine=(target.cosine_similarity_column, "max"),
            )
        )
        target_df["cosine"] = pd.to_numeric(target_df["cosine"], errors="coerce")

        ax = axes[row_idx][0]
        if target_df["cosine"].notna().any():
            ax.hexbin(
                target_df["cosine"], target_df["frequency"],
                gridsize=16, cmap="BuGn", mincnt=1,
                edgecolors="#CBD5E1", linewidths=0.3,
            )
            ax.set_title(
                f"Cosine similarity vs frequency: {target.keyword}",
                fontsize=10, fontweight="semibold", color="#1E293B",
            )
            ax.set_xlabel(target.cosine_similarity_column, fontsize=8.5, color="#475569")
            ax.set_ylabel("Phenotype frequency", fontsize=8.5, color="#475569")
            _style_axes(ax, grid=True)
        else:
            _plot_no_data(
                ax,
                f"Cosine similarity vs frequency: {target.keyword}",
                "OpenAI embedding similarity is unavailable for this run.",
            )
    fig.tight_layout(pad=1.5)
    path = output_dir / "similarity_vs_frequency.png"
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    written_paths.append(path)

    return written_paths
