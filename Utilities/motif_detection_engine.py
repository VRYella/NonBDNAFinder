"""
Motif Detection Engine: adaptive orchestrator that selects analysis strategy by sequence length.

    sequence_length < 100_000 bp  → numba_only (direct analysis)
    100_000 ≤ length < 5_000_000 → chunk_workers (2-worker ProcessPool)
    length ≥ 5_000_000 bp         → disk_streaming (sequential, constant RAM)
"""

from __future__ import annotations

import logging
from typing import Any, Callable, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

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
        self.sequence_storage = sequence_storage
        self.chunk_size = chunk_size
        self.overlap = overlap
        self.max_workers = max_workers

    # PUBLIC API

    def analyze(
        self,
        seq_id: str,
        progress_callback: Optional[Callable[[float], None]] = None,
        enabled_classes: Optional[List[str]] = None,
    ) -> Tuple[Any, Any]:
        """Analyse a stored sequence. Returns (UniversalResultsStorage, VisualizationAccumulator)."""
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

        accumulator = VisualizationAccumulator(seq_length=seq_length)
        for motif in results_storage.iter_results():
            accumulator.update([motif])

        logger.info(
            f"MotifDetectionEngine: analysis complete "
            f"({results_storage.get_summary_stats()['total_count']} motifs)"
        )
        return results_storage, accumulator

    # STRATEGY SELECTION

    @staticmethod
    def _select_strategy(seq_length: int) -> str:
        """Select 'numba_only', 'chunk_workers', or 'disk_streaming' by length."""
        if seq_length < _THRESHOLD_DIRECT:
            return "numba_only"
        if seq_length < _THRESHOLD_CHUNK:
            return "chunk_workers"
        return "disk_streaming"
