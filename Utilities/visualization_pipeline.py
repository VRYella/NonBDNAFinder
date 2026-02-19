"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Visualization Pipeline - Streamlit-Safe On-Demand Plot Rendering             │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Visualization pipeline that operates exclusively on pre-binned summary data
    produced by ``VisualizationAccumulator``.

    The full motif table is *never* loaded into memory for plotting.  Each plot
    is rendered on demand from fixed-size arrays (``density_bins``,
    ``length_bins``, ``cooccurrence_matrix``) whose size is constant with
    respect to genome size.

    Replaces:
        - Raw scatter plotting (replaced by pre-binned density histogram)
        - Large in-memory histograms (replaced by fixed-length numpy arrays)
        - Multi-megabyte matplotlib objects held between renders

    Memory is bounded because:
        - All inputs are summary arrays of fixed length (``bin_count`` ≤ 200)
        - Figures are rendered and returned; the caller controls lifetime
        - No global state is accumulated between calls

USAGE::

    from Utilities.visualization_accumulator import VisualizationAccumulator
    from Utilities.visualization_pipeline import VisualizationPipeline

    acc = VisualizationAccumulator(seq_length=5_000_000)
    # … populate acc …

    pipeline = VisualizationPipeline(acc.get_summary())
    fig_density  = pipeline.plot_density_histogram()
    fig_lengths  = pipeline.plot_length_histogram()
    fig_classes  = pipeline.plot_class_distribution()
    fig_cooccur  = pipeline.plot_cooccurrence_heatmap()
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# Lazy matplotlib imports – avoids loading the heavy library at module import
# time (important for Streamlit cold-start performance).
_plt = None
_sns = None


def _ensure_matplotlib():
    """Lazy-import matplotlib and seaborn; return (plt, sns)."""
    global _plt, _sns
    if _plt is None:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns
        _plt = plt
        _sns = sns
    return _plt, _sns


class VisualizationPipeline:
    """
    On-demand plot renderer backed by pre-binned summary statistics.

    All plot methods accept optional ``figsize`` and ``title`` overrides and
    return a ``matplotlib.figure.Figure`` that the caller can display or save.

    Parameters
    ----------
    summary : dict
        Output of ``VisualizationAccumulator.get_summary()``.
    """

    def __init__(self, summary: Dict[str, Any]):
        self._summary = summary
        self._seq_length: int = summary.get("seq_length", 1)
        self._bin_count: int = summary.get("bin_count", 100)
        self._max_length: int = summary.get("max_length", 10_000)

    # ------------------------------------------------------------------
    # DENSITY HISTOGRAM
    # ------------------------------------------------------------------

    def plot_density_histogram(
        self,
        figsize: Tuple[int, int] = (10, 4),
        title: str = "Motif Positional Density",
        color: str = "#1f77b4",
    ):
        """
        Bar chart of motif counts per genomic position bin.

        Args:
            figsize: Figure size (width, height) in inches.
            title:   Plot title.
            color:   Bar colour.

        Returns:
            matplotlib Figure.
        """
        plt, _ = _ensure_matplotlib()
        density_bins: np.ndarray = self._summary.get(
            "density_bins", np.zeros(self._bin_count, dtype=np.int64)
        )

        bin_edges = np.linspace(0, self._seq_length, self._bin_count + 1)
        bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_width = bin_edges[1] - bin_edges[0]

        fig, ax = plt.subplots(figsize=figsize)
        ax.bar(
            bin_centres,
            density_bins,
            width=bin_width * 0.9,
            color=color,
            alpha=0.8,
            edgecolor="none",
        )
        ax.set_xlabel("Genomic Position (bp)", fontsize=11)
        ax.set_ylabel("Motif Count", fontsize=11)
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.grid(axis="y", alpha=0.3)
        plt.tight_layout()
        return fig

    # ------------------------------------------------------------------
    # LENGTH HISTOGRAM
    # ------------------------------------------------------------------

    def plot_length_histogram(
        self,
        figsize: Tuple[int, int] = (8, 4),
        title: str = "Motif Length Distribution",
        color: str = "#ff7f0e",
    ):
        """
        Bar chart of motif length distribution.

        Args:
            figsize: Figure size.
            title:   Plot title.
            color:   Bar colour.

        Returns:
            matplotlib Figure.
        """
        plt, _ = _ensure_matplotlib()
        length_bins: np.ndarray = self._summary.get(
            "length_bins", np.zeros(self._bin_count, dtype=np.int64)
        )

        bin_edges = np.linspace(0, self._max_length, self._bin_count + 1)
        bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_width = bin_edges[1] - bin_edges[0]

        fig, ax = plt.subplots(figsize=figsize)
        ax.bar(
            bin_centres,
            length_bins,
            width=bin_width * 0.9,
            color=color,
            alpha=0.8,
            edgecolor="none",
        )
        ax.set_xlabel("Motif Length (bp)", fontsize=11)
        ax.set_ylabel("Count", fontsize=11)
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.grid(axis="y", alpha=0.3)
        plt.tight_layout()
        return fig

    # ------------------------------------------------------------------
    # CLASS DISTRIBUTION BAR CHART
    # ------------------------------------------------------------------

    def plot_class_distribution(
        self,
        figsize: Tuple[int, int] = (8, 5),
        title: str = "Motif Class Distribution",
        color: str = "#2ca02c",
        top_n: int = 20,
    ):
        """
        Horizontal bar chart of motif class counts.

        Args:
            figsize: Figure size.
            title:   Plot title.
            color:   Bar colour.
            top_n:   Maximum number of classes to display.

        Returns:
            matplotlib Figure.
        """
        plt, _ = _ensure_matplotlib()
        class_counts: Dict[str, int] = self._summary.get("class_counts", {})

        if not class_counts:
            return self._empty_figure(title, figsize)

        sorted_items = sorted(class_counts.items(), key=lambda x: x[1], reverse=True)
        sorted_items = sorted_items[:top_n]
        labels = [item[0] for item in sorted_items]
        values = [item[1] for item in sorted_items]

        fig, ax = plt.subplots(figsize=figsize)
        ax.barh(labels[::-1], values[::-1], color=color, alpha=0.85)
        ax.set_xlabel("Count", fontsize=11)
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.grid(axis="x", alpha=0.3)
        plt.tight_layout()
        return fig

    # ------------------------------------------------------------------
    # SUBCLASS DISTRIBUTION BAR CHART
    # ------------------------------------------------------------------

    def plot_subclass_distribution(
        self,
        figsize: Tuple[int, int] = (9, 6),
        title: str = "Motif Subclass Distribution",
        top_n: int = 25,
    ):
        """
        Horizontal bar chart of motif subclass counts.

        Args:
            figsize: Figure size.
            title:   Plot title.
            top_n:   Maximum number of subclasses to display.

        Returns:
            matplotlib Figure.
        """
        plt, _ = _ensure_matplotlib()
        subclass_counts: Dict[str, int] = self._summary.get("subclass_counts", {})

        if not subclass_counts:
            return self._empty_figure(title, figsize)

        sorted_items = sorted(
            subclass_counts.items(), key=lambda x: x[1], reverse=True
        )
        sorted_items = sorted_items[:top_n]
        labels = [item[0] for item in sorted_items]
        values = [item[1] for item in sorted_items]

        import matplotlib
        try:
            cmap = matplotlib.colormaps.get_cmap("tab20").resampled(len(labels))
        except AttributeError:
            # Fallback for matplotlib < 3.7
            cmap = matplotlib.cm.get_cmap("tab20", len(labels))
        colors = [cmap(i) for i in range(len(labels))]

        fig, ax = plt.subplots(figsize=figsize)
        ax.barh(labels[::-1], values[::-1], color=colors[::-1], alpha=0.85)
        ax.set_xlabel("Count", fontsize=11)
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.grid(axis="x", alpha=0.3)
        plt.tight_layout()
        return fig

    # ------------------------------------------------------------------
    # CO-OCCURRENCE HEATMAP
    # ------------------------------------------------------------------

    def plot_cooccurrence_heatmap(
        self,
        figsize: Tuple[int, int] = (8, 6),
        title: str = "Motif Class Co-occurrence",
    ):
        """
        Heatmap of class-level co-occurrence counts.

        Args:
            figsize: Figure size.
            title:   Plot title.

        Returns:
            matplotlib Figure.
        """
        plt, sns = _ensure_matplotlib()
        cooccurrence: Dict[str, Dict[str, int]] = self._summary.get(
            "cooccurrence_matrix", {}
        )

        if not cooccurrence:
            return self._empty_figure(title, figsize)

        classes = sorted(cooccurrence.keys())
        n = len(classes)
        matrix = np.zeros((n, n), dtype=np.int64)
        for i, ca in enumerate(classes):
            for j, cb in enumerate(classes):
                matrix[i, j] = cooccurrence.get(ca, {}).get(cb, 0)

        fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(
            matrix,
            xticklabels=classes,
            yticklabels=classes,
            annot=n <= 15,
            fmt="d",
            cmap="YlOrRd",
            ax=ax,
            cbar_kws={"label": "Co-occurrence count"},
        )
        ax.set_title(title, fontsize=13, fontweight="bold")
        plt.xticks(rotation=45, ha="right", fontsize=9)
        plt.yticks(rotation=0, fontsize=9)
        plt.tight_layout()
        return fig

    # ------------------------------------------------------------------
    # GENERATE ALL PLOTS
    # ------------------------------------------------------------------

    def generate_all(self) -> Dict[str, Any]:
        """
        Render all standard plots and return them in a dict.

        Returns:
            Dict mapping plot name → matplotlib Figure::

                {
                    "density_histogram":     Figure,
                    "length_histogram":      Figure,
                    "class_distribution":    Figure,
                    "subclass_distribution": Figure,
                    "cooccurrence_heatmap":  Figure,
                }
        """
        return {
            "density_histogram": self.plot_density_histogram(),
            "length_histogram": self.plot_length_histogram(),
            "class_distribution": self.plot_class_distribution(),
            "subclass_distribution": self.plot_subclass_distribution(),
            "cooccurrence_heatmap": self.plot_cooccurrence_heatmap(),
        }

    # ------------------------------------------------------------------
    # HELPER
    # ------------------------------------------------------------------

    @staticmethod
    def _empty_figure(title: str, figsize: Tuple[int, int]):
        """Return a blank figure with a 'No data' message."""
        plt, _ = _ensure_matplotlib()
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(
            0.5, 0.5, "No data available",
            ha="center", va="center", fontsize=12, color="grey",
            transform=ax.transAxes,
        )
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.axis("off")
        plt.tight_layout()
        return fig
