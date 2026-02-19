"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Disk Chunk Manager - Disk-First Parallel Architecture                        │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Disk-backed chunk manager for memory-bounded parallel genome analysis.

    Instead of returning full motif tables from parallel workers, workers:
      1. Process a sequence chunk
      2. Write results to a temporary CSV file
      3. Return only lightweight metadata (file path, motif count)

    The main process then:
      - Streams chunk files via an iterator
      - Aggregates only summary statistics
      - Deletes each chunk file after it has been merged

    This keeps RAM bounded regardless of genome size.

TEMPORARY FILE STRUCTURE:
    /tmp/nonbdna_<session>/
        chunk_000.csv
        chunk_001.csv
        ...

MEMORY GUARANTEES:
    - Workers never accumulate large DataFrames in memory
    - Each chunk CSV is deleted immediately after aggregation
    - Peak RAM usage scales with chunk_size, not total genome size
"""

import csv
import gc
import json
import logging
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple

logger = logging.getLogger(__name__)

# Fieldnames written to each chunk CSV
_CSV_FIELDS = ["Class", "Subclass", "Start", "End", "Length", "Score", "ID", "Strand"]


class DiskChunkManager:
    """
    Disk-backed chunk manager for memory-bounded genome analysis.

    Workers write results to per-chunk CSV files; the manager streams
    and aggregates those files without loading the full dataset into RAM.

    Usage::

        manager = DiskChunkManager()
        try:
            chunk_id = manager.write_chunk_results(chunk_index=0, motifs=[...])
            for motifs, meta in manager.iter_chunk_results():
                # process motifs
                pass
        finally:
            manager.cleanup()
    """

    def __init__(self, base_dir: Optional[str] = None):
        """
        Initialise the manager and create a temporary working directory.

        Args:
            base_dir: Optional directory in which to create the temp folder.
                      Defaults to the system temp directory.
        """
        self.base_dir = Path(
            tempfile.mkdtemp(prefix="nonbdna_chunks_", dir=base_dir)
        )
        # Map chunk_index -> metadata dict
        self._chunk_registry: Dict[int, Dict[str, Any]] = {}
        logger.info(f"DiskChunkManager initialised at {self.base_dir}")

    # ------------------------------------------------------------------
    # WRITING
    # ------------------------------------------------------------------

    def write_chunk_results(
        self, chunk_index: int, motifs: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """
        Write motif results for one chunk to a CSV file.

        Args:
            chunk_index: Zero-based index of the chunk.
            motifs:      List of motif dicts produced by ``analyze_sequence``.

        Returns:
            Lightweight metadata dict::

                {
                    "chunk_index": int,
                    "file_path":   str,
                    "motif_count": int,
                }
        """
        file_path = self.base_dir / f"chunk_{chunk_index:04d}.csv"

        with open(file_path, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=_CSV_FIELDS, extrasaction="ignore"
            )
            writer.writeheader()
            for motif in motifs:
                writer.writerow(motif)

        meta = {
            "chunk_index": chunk_index,
            "file_path": str(file_path),
            "motif_count": len(motifs),
        }
        self._chunk_registry[chunk_index] = meta
        logger.debug(
            f"Chunk {chunk_index}: wrote {len(motifs)} motifs to {file_path}"
        )
        return meta

    # ------------------------------------------------------------------
    # READING / STREAMING
    # ------------------------------------------------------------------

    def iter_chunk_results(
        self, delete_after_read: bool = True
    ) -> Iterator[Tuple[List[Dict[str, Any]], Dict[str, Any]]]:
        """
        Iterate over chunk results in index order, streaming from disk.

        Each chunk file is read, yielded, and (optionally) deleted to keep
        disk usage bounded.

        Args:
            delete_after_read: Remove the CSV file after yielding its rows.

        Yields:
            Tuple of (motifs, metadata) for each recorded chunk.
        """
        for chunk_index in sorted(self._chunk_registry.keys()):
            meta = self._chunk_registry[chunk_index]
            file_path = Path(meta["file_path"])

            if not file_path.exists():
                logger.warning(
                    f"Chunk {chunk_index} file missing: {file_path}"
                )
                continue

            motifs: List[Dict[str, Any]] = []
            with open(file_path, newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    # Convert numeric fields back from strings
                    for int_field in ("Start", "End", "Length"):
                        if row.get(int_field):
                            try:
                                row[int_field] = int(row[int_field])
                            except ValueError:
                                pass
                    if row.get("Score"):
                        try:
                            row["Score"] = float(row["Score"])
                        except ValueError:
                            pass
                    motifs.append(dict(row))

            yield motifs, meta

            if delete_after_read:
                file_path.unlink(missing_ok=True)
                logger.debug(f"Chunk {chunk_index}: deleted {file_path}")

            gc.collect()

    def get_chunk_metadata(self) -> List[Dict[str, Any]]:
        """
        Return metadata for all registered chunks (no disk I/O).

        Returns:
            List of metadata dicts sorted by chunk_index.
        """
        return [
            self._chunk_registry[i]
            for i in sorted(self._chunk_registry.keys())
        ]

    def total_motif_count(self) -> int:
        """Return the sum of motif counts across all chunk metadata records."""
        return sum(m["motif_count"] for m in self._chunk_registry.values())

    # ------------------------------------------------------------------
    # CLEANUP
    # ------------------------------------------------------------------

    def cleanup(self):
        """Remove the entire temporary directory and all chunk files."""
        if self.base_dir.exists():
            shutil.rmtree(self.base_dir, ignore_errors=True)
            logger.info(f"DiskChunkManager: cleaned up {self.base_dir}")
        self._chunk_registry.clear()
