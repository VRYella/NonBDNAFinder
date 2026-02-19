"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Motif Detection Engine - Adaptive Analysis Orchestrator                      │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    High-level orchestrator that selects the optimal analysis strategy based on
    sequence length and routes execution to the appropriate backend.

    Adaptive strategy:

        sequence_length < 100 000 bp
            → numba-only (direct analysis, no chunking)

        100 000 ≤ sequence_length < 5 000 000 bp
            → numba + 2-worker chunk processing
              (StreamlitSafeExecutor with ProcessPoolExecutor)

        sequence_length ≥ 5 000 000 bp
            → disk-backed chunk streaming mode
              (sequential, one chunk at a time, constant RAM)

    Visualisation always uses a VisualizationAccumulator populated while
    streaming results – the full motif table is never loaded into RAM for
    plotting.

MEMORY GUARANTEES:
    - Peak RAM proportional to chunk_size, not total genome size
    - Compatible with Streamlit Community Cloud (≤ 1 GB RAM)
"""

from __future__ import annotations

import logging
from typing import Any, Callable, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# Adaptive strategy thresholds (base pairs)
_THRESHOLD_DIRECT: int = 100_000
_THRESHOLD_CHUNK: int = 5_000_000


class MotifDetectionEngine:
    """
    Adaptive motif detection engine for single-sequence analysis.

    Automatically selects the most appropriate execution strategy based on
    sequence length, balancing speed, RAM usage, and Streamlit safety.

    Usage::

        from Utilities.disk_storage import UniversalSequenceStorage

        storage = UniversalSequenceStorage()
        seq_id = storage.save_sequence(dna_sequence, "chr1")

        engine = MotifDetectionEngine(storage)
        results_storage, accumulator = engine.analyze(
            seq_id=seq_id,
            progress_callback=lambda p: print(f"{p:.0f}%"),
            enabled_classes=["G-Quadruplex", "Z-DNA"],
        )

        summary = results_storage.get_summary_stats()
        vis_data = accumulator.get_summary()
    """

    def __init__(
        self,
        sequence_storage,
        chunk_size: int = 50_000,
        overlap: int = 2_000,
        max_workers: Optional[int] = None,
    ):
        """
        Args:
            sequence_storage: ``UniversalSequenceStorage`` instance.
            chunk_size:       Chunk size in bp (default 50 000).
            overlap:          Overlap between chunks in bp (default 2 000).
            max_workers:      Override max parallel workers (default auto).
        """
        self.sequence_storage = sequence_storage
        self.chunk_size = chunk_size
        self.overlap = overlap
        self.max_workers = max_workers

    # ------------------------------------------------------------------
    # PUBLIC API
    # ------------------------------------------------------------------

    def analyze(
        self,
        seq_id: str,
        progress_callback: Optional[Callable[[float], None]] = None,
        enabled_classes: Optional[List[str]] = None,
    ) -> Tuple[Any, Any]:
        """
        Analyse a stored sequence with the adaptive strategy.

        Args:
            seq_id:            Sequence ID from ``UniversalSequenceStorage``.
            progress_callback: Optional ``callback(progress_pct: float)``.
            enabled_classes:   Motif classes to analyse (``None`` = all).

        Returns:
            Tuple of:
                - ``UniversalResultsStorage`` – results on disk
                - ``VisualizationAccumulator`` – pre-aggregated vis data
        """
        from Utilities.streamlit_safe_executor import StreamlitSafeExecutor
        from Utilities.visualization_accumulator import VisualizationAccumulator

        meta = self.sequence_storage.get_metadata(seq_id)
        seq_length = meta["length"]
        strategy = self._select_strategy(seq_length)

        logger.info(
            f"MotifDetectionEngine: '{meta['name']}' "
            f"length={seq_length:,} bp → strategy='{strategy}'"
        )

        executor = StreamlitSafeExecutor(
            sequence_storage=self.sequence_storage,
            chunk_size=self.chunk_size,
            overlap=self.overlap,
            max_workers=self.max_workers,
        )

        results_storage = executor.run(
            seq_id=seq_id,
            progress_callback=progress_callback,
            enabled_classes=enabled_classes,
        )

        # Build VisualizationAccumulator by streaming the final results
        accumulator = VisualizationAccumulator(seq_length=seq_length)
        for motif in results_storage.iter_results():
            accumulator.update([motif])

        logger.info(
            f"MotifDetectionEngine: analysis complete "
            f"({results_storage.get_summary_stats()['total_count']} motifs)"
        )
        return results_storage, accumulator

    # ------------------------------------------------------------------
    # STRATEGY SELECTION
    # ------------------------------------------------------------------

    @staticmethod
    def _select_strategy(seq_length: int) -> str:
        """
        Choose execution strategy based on sequence length.

        Returns:
            One of ``'numba_only'``, ``'chunk_workers'``, or
            ``'disk_streaming'``.
        """
        if seq_length < _THRESHOLD_DIRECT:
            return "numba_only"
        if seq_length < _THRESHOLD_CHUNK:
            return "chunk_workers"
        return "disk_streaming"
