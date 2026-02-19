"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Final Exporter - Stream-Merge Chunk Results for Full-Table Export            │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Assembles the final, deduplicated motif table *only* when the user requests
    an explicit export.  The full table is never kept in memory during analysis;
    it is only materialised here on demand by streaming and merging per-chunk
    files from a DiskChunkManager (or compatible UniversalResultsStorage).

    Two export targets are supported:

        1. Generator-based streaming  (``assemble``)
           Yields deduplicated chunks as lists of motif dicts.  Use this when
           you want to write results to disk incrementally.

        2. Full DataFrame assembly    (``to_dataframe``)
           Concatenates everything into a single pandas DataFrame.  Only call
           when the caller has confirmed sufficient RAM.

USAGE::

    from Utilities.disk_chunk_manager import DiskChunkManager
    from Utilities.overlap_deduplicator import OverlapDeduplicator

    exporter  = FinalExporter()
    dedup     = OverlapDeduplicator()

    # Stream deduplicated chunks
    for chunk_motifs in exporter.assemble(disk_manager, dedup):
        write_to_csv(chunk_motifs)

    # Or build a DataFrame (RAM-permitting)
    df = exporter.to_dataframe(disk_manager, dedup)
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Generator, List, Optional

logger = logging.getLogger(__name__)


class FinalExporter:
    """
    Assemble the complete, deduplicated motif table from disk-stored chunks.

    Results are only materialised when :meth:`assemble` or
    :meth:`to_dataframe` is called, keeping memory usage constant during
    the analysis phase.

    Usage::

        exporter = FinalExporter()
        for motifs in exporter.assemble(disk_manager, deduplicator):
            process(motifs)
    """

    def assemble(
        self,
        disk_manager,
        deduplicator,
        delete_after_read: bool = False,
    ) -> Generator[List[Dict[str, Any]], None, None]:
        """
        Stream deduplicated motif chunks from *disk_manager*.

        Yields one list of motif dicts per stored chunk, after applying
        *deduplicator*.  Chunk files may optionally be deleted after reading
        to bound disk usage.

        Compatible with:
            - :class:`~Utilities.disk_chunk_manager.DiskChunkManager`
              (``iter_chunk_results`` interface)
            - :class:`~Utilities.disk_storage.UniversalResultsStorage`
              (``iter_results`` interface – treated as a single chunk)

        Args:
            disk_manager:      Storage object providing chunk iteration.
            deduplicator:      :class:`~Utilities.overlap_deduplicator.OverlapDeduplicator`
                               instance used to filter overlap regions.
            delete_after_read: Whether to delete chunk files after reading.

        Yields:
            Filtered list of motif dicts for each chunk.
        """
        # DiskChunkManager interface: iter_chunk_results yields (motifs, meta)
        if hasattr(disk_manager, "iter_chunk_results"):
            for motifs, meta in disk_manager.iter_chunk_results(
                delete_after_read=delete_after_read
            ):
                core_end = meta.get("core_end", meta.get("chunk_end", float("inf")))
                filtered = deduplicator.filter_core(motifs, int(core_end))
                logger.debug(
                    f"FinalExporter chunk {meta.get('chunk_index', '?')}: "
                    f"{len(motifs)} raw → {len(filtered)} after dedup"
                )
                yield filtered

        # UniversalResultsStorage interface: iter_results yields individual motifs
        elif hasattr(disk_manager, "iter_results"):
            motifs = list(disk_manager.iter_results())
            # Results already deduplicated at write time; yield as-is
            yield motifs

        else:
            raise TypeError(
                f"disk_manager of type {type(disk_manager).__name__} does not "
                "provide 'iter_chunk_results' or 'iter_results'."
            )

    def to_dataframe(
        self,
        disk_manager,
        deduplicator,
        delete_after_read: bool = False,
    ):
        """
        Assemble all chunks into a single pandas DataFrame.

        .. warning::
            This loads the entire motif table into RAM.  Only call when the
            caller has verified that sufficient memory is available.

        Args:
            disk_manager:      Storage object (see :meth:`assemble`).
            deduplicator:      Deduplicator instance.
            delete_after_read: Whether to delete chunk files after reading.

        Returns:
            ``pandas.DataFrame`` containing all deduplicated motifs, or an
            empty DataFrame if no motifs were found.
        """
        import pandas as pd

        all_motifs: List[Dict[str, Any]] = []
        for chunk_motifs in self.assemble(
            disk_manager, deduplicator, delete_after_read=delete_after_read
        ):
            all_motifs.extend(chunk_motifs)

        if not all_motifs:
            logger.info("FinalExporter.to_dataframe: no motifs to assemble")
            return pd.DataFrame()

        df = pd.DataFrame(all_motifs)
        logger.info(f"FinalExporter.to_dataframe: assembled {len(df):,} motifs")
        return df
