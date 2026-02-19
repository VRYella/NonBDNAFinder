"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Overlap Deduplicator - Remove Duplicate Motifs at Chunk Boundaries           │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Removes duplicate motifs that arise when the same genomic region is
    analysed by two consecutive overlapping chunks.

    Strategy: only keep motifs whose ``Start`` position is strictly less than
    ``core_end``.  Motifs at or beyond ``core_end`` are within the overlap
    region and will be re-detected (and kept) by the next chunk.

    This is equivalent to the "core-region" filtering described in the
    architecture blueprint.

USAGE::

    dedup = OverlapDeduplicator()

    # With a list of motif dicts:
    unique = dedup.filter_core(motifs, core_end=48_000)

    # With a pandas DataFrame:
    df_unique = dedup.filter_core_df(df, core_end=48_000)
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List

logger = logging.getLogger(__name__)


class OverlapDeduplicator:
    """
    Filter motifs to the core (authoritative) region of a chunk.

    The *core region* of a chunk spans ``[chunk_start, core_end)``.  Only
    motifs that **start** within this region are kept.  Motifs that start
    in the overlap zone (``[core_end, chunk_end)``) will be re-detected – and
    kept – when the next chunk is processed.

    Usage::

        dedup = OverlapDeduplicator()
        unique = dedup.filter_core(motifs, core_end=48_000)
    """

    def filter_core(
        self,
        motifs: List[Dict[str, Any]],
        core_end: int,
    ) -> List[Dict[str, Any]]:
        """
        Return only motifs whose ``Start`` is strictly less than *core_end*.

        Args:
            motifs:   List of motif dicts.  Each must contain a ``"Start"``
                      key with an integer-compatible value.
            core_end: Exclusive upper boundary of the authoritative region
                      (from ``ChunkGenerator`` chunk dict ``"core_end"``).

        Returns:
            Filtered list containing only core-region motifs.
        """
        result = []
        for motif in motifs:
            try:
                start = int(motif.get("Start", 0))
            except (ValueError, TypeError):
                start = 0
            if start < core_end:
                result.append(motif)

        logger.debug(
            f"OverlapDeduplicator: {len(motifs)} motifs → "
            f"{len(result)} kept (core_end={core_end:,})"
        )
        return result

    def filter_core_df(self, motif_df, core_end: int):
        """
        DataFrame variant: return rows where ``Start < core_end``.

        Args:
            motif_df: ``pandas.DataFrame`` with a ``"Start"`` column.
            core_end: Exclusive boundary of the authoritative region.

        Returns:
            Filtered ``pandas.DataFrame``.
        """
        if motif_df is None or len(motif_df) == 0:
            return motif_df
        return motif_df[motif_df["Start"] < core_end]
