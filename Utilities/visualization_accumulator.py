"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Visualization Accumulator - Streaming Aggregation for Plots                  │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Maintains fixed-size summary arrays that are updated incrementally as motif
    results stream in from disk, so the full motif table never needs to be loaded
    into memory for visualization.

    Tracks:
        - class_counts          : {class_name: int}
        - subclass_counts       : {subclass_name: int}
        - density_bins          : np.ndarray – motifs-per-kb histogram
        - length_bins           : np.ndarray – motif length histogram
        - cooccurrence_matrix   : {class_a: {class_b: int}}

MEMORY GUARANTEES:
    - All maintained structures have O(1) size relative to genome length
    - Arrays are fixed length (bin_count buckets)
    - Suitable for Streamlit Community Cloud (≤1 GB RAM)

USAGE::

    acc = VisualizationAccumulator(seq_length=5_000_000)
    for motifs, _ in manager.iter_chunk_results():
        acc.update(motifs)

    summary = acc.get_summary()
    # summary["class_counts"], summary["density_bins"], etc.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Any, Dict, List, Optional

import numpy as np

logger = logging.getLogger(__name__)

# Default number of histogram bins – kept constant regardless of genome size
_DEFAULT_BINS = 100


class VisualizationAccumulator:
    """
    Streaming accumulator for visualization-ready summary statistics.

    Call :meth:`update` with a list of motif dicts each time a chunk of
    results is available.  When processing is complete, call
    :meth:`get_summary` to retrieve all aggregated data structures needed
    to render plots.

    All internal arrays have fixed length (``bin_count``) so memory usage
    is constant with respect to genome size.
    """

    def __init__(
        self,
        seq_length: int = 0,
        bin_count: int = _DEFAULT_BINS,
        max_length: int = 10_000,
    ):
        """
        Initialise the accumulator.

        Args:
            seq_length:  Total length of the analysed sequence in bp.
                         Used to compute positional density bins.
            bin_count:   Number of fixed bins for histograms (default 100).
            max_length:  Maximum motif length tracked by the length histogram
                         (default 10 000 bp; longer motifs are clamped to the
                         last bin).
        """
        self.seq_length = max(seq_length, 1)
        self.bin_count = bin_count
        self.max_length = max(max_length, 1)

        # --- per-class / per-subclass counts ---
        self.class_counts: Dict[str, int] = defaultdict(int)
        self.subclass_counts: Dict[str, int] = defaultdict(int)

        # --- positional density histogram ---
        # Bins span [0, seq_length); each bin covers seq_length/bin_count bp
        self.density_bins: np.ndarray = np.zeros(bin_count, dtype=np.int64)

        # --- motif length histogram ---
        # Bins span [1, max_length]; motifs longer than max_length go in the
        # last bin.
        self.length_bins: np.ndarray = np.zeros(bin_count, dtype=np.int64)

        # --- class co-occurrence matrix ---
        # cooccurrence_matrix[class_a][class_b] = number of motifs found
        # in the same chunk where both class_a and class_b appear.
        self.cooccurrence_matrix: Dict[str, Dict[str, int]] = defaultdict(
            lambda: defaultdict(int)
        )

        self._total_motifs: int = 0

        logger.debug(
            f"VisualizationAccumulator ready "
            f"(seq_length={seq_length:,}, bins={bin_count}, max_len={max_length})"
        )

    # ------------------------------------------------------------------
    # INCREMENTAL UPDATE
    # ------------------------------------------------------------------

    def update(self, motifs: List[Dict[str, Any]]) -> None:
        """
        Incorporate a batch of motifs into the accumulator.

        This method is designed to be called repeatedly as chunks of results
        stream in from disk.  Each call is O(len(motifs)).

        Args:
            motifs: List of motif dicts.  Each dict should contain at minimum:
                    ``Class``, ``Subclass``, ``Start``, ``End`` / ``Length``,
                    and optionally ``Score``.
        """
        if not motifs:
            return

        classes_in_batch: List[str] = []

        for motif in motifs:
            cls = motif.get("Class") or "Unknown"
            subcls = motif.get("Subclass") or "Unknown"
            start = motif.get("Start", 0)

            # --- Motif length ---
            length = motif.get("Length")
            if length is None:
                end = motif.get("End", start)
                length = max(int(end) - int(start), 0)
            else:
                try:
                    length = int(length)
                except (ValueError, TypeError):
                    length = 0

            # 1. Class / subclass counts
            self.class_counts[cls] += 1
            self.subclass_counts[subcls] += 1

            # 2. Positional density bin
            try:
                pos = int(start)
            except (ValueError, TypeError):
                pos = 0
            bin_idx = min(
                int(pos * self.bin_count / self.seq_length),
                self.bin_count - 1,
            )
            self.density_bins[bin_idx] += 1

            # 3. Length histogram bin
            clamped_len = min(max(length, 0), self.max_length)
            len_bin = min(
                int(clamped_len * self.bin_count / self.max_length),
                self.bin_count - 1,
            )
            self.length_bins[len_bin] += 1

            classes_in_batch.append(cls)

        # 4. Co-occurrence: O(k²) where k = distinct classes in this chunk
        self._total_motifs += len(motifs)
        unique_classes = list(set(classes_in_batch))
        for i, ca in enumerate(unique_classes):
            for cb in unique_classes:
                self.cooccurrence_matrix[ca][cb] += 1

    # ------------------------------------------------------------------
    # SUMMARY EXPORT
    # ------------------------------------------------------------------

    def get_summary(self) -> Dict[str, Any]:
        """
        Return a snapshot of all accumulated statistics.

        Returns:
            Dict with keys:

            * ``total_motifs``        – int
            * ``class_counts``        – ``{class: count}``
            * ``subclass_counts``     – ``{subclass: count}``
            * ``density_bins``        – ``np.ndarray`` shape ``(bin_count,)``
            * ``length_bins``         – ``np.ndarray`` shape ``(bin_count,)``
            * ``cooccurrence_matrix`` – ``{class_a: {class_b: count}}``
            * ``bin_count``           – int
            * ``seq_length``          – int
            * ``max_length``          – int
        """
        return {
            "total_motifs": self._total_motifs,
            "class_counts": dict(self.class_counts),
            "subclass_counts": dict(self.subclass_counts),
            "density_bins": self.density_bins.copy(),
            "length_bins": self.length_bins.copy(),
            "cooccurrence_matrix": {
                ca: dict(cb_counts)
                for ca, cb_counts in self.cooccurrence_matrix.items()
            },
            "bin_count": self.bin_count,
            "seq_length": self.seq_length,
            "max_length": self.max_length,
        }

    def reset(self) -> None:
        """Reset all accumulators to zero (reuse the same instance)."""
        self.class_counts = defaultdict(int)
        self.subclass_counts = defaultdict(int)
        self.density_bins = np.zeros(self.bin_count, dtype=np.int64)
        self.length_bins = np.zeros(self.bin_count, dtype=np.int64)
        self.cooccurrence_matrix = defaultdict(lambda: defaultdict(int))
        self._total_motifs = 0
        logger.debug("VisualizationAccumulator reset")
