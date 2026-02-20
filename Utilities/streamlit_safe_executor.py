"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Streamlit-Safe Executor - Controlled Parallel Analysis                       │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Streamlit Community Cloud compatible executor with hard limits on process
    count and memory.

    Key design constraints:
        - max_workers = min(2, os.cpu_count())
        - No nested process pools
        - Chunk *indices* are passed to workers, not DataFrames
        - Workers write results to disk; only metadata is returned
        - No background threads/services that outlive the Streamlit request

    Adaptive execution strategy based on sequence length:
        < 100 000 bp   → direct (no workers, no chunking)
        100 000 – 5 M  → 2-worker ProcessPoolExecutor with chunk indices
        > 5 000 000 bp → disk-backed streaming (chunks written by workers)

MEMORY GUARANTEES:
    - Peak RAM scales with chunk_size, not genome size
    - Suitable for 1 GB RAM containers
    - No large objects passed across process boundaries
"""

from __future__ import annotations

import gc
import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Callable, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────────────────
# CONSTANTS
# ──────────────────────────────────────────────────────────────────────────────
# Hard ceiling on worker count to protect Streamlit Community Cloud
MAX_WORKERS: int = min(2, os.cpu_count() or 1)

# Thresholds for adaptive strategy selection (base pairs)
# THRESHOLD_DIRECT = 0 forces ALL sequences through chunking (50K/2K)
# This ensures consistent analysis with overlap deduplication regardless of size
THRESHOLD_DIRECT: int = 0              # always chunk (50K/2K for every sequence)
THRESHOLD_CHUNK: int = 5_000_000       # below → chunk workers; above → disk streaming

# Default chunk parameters
DEFAULT_CHUNK_SIZE: int = 50_000
DEFAULT_OVERLAP: int = 2_000


# ──────────────────────────────────────────────────────────────────────────────
# WORKER (module-level, picklable)
# ──────────────────────────────────────────────────────────────────────────────

def _chunk_worker(
    args: Tuple[str, str, int, int, Optional[List[str]], str],
) -> Dict[str, Any]:
    """
    Process a single chunk and write results to a CSV file on disk.

    This is a *module-level* function so it is picklable by
    ``ProcessPoolExecutor`` without serialising any large objects.

    Args:
        args: Tuple of
            (chunk_seq, seq_name, chunk_start, chunk_end, enabled_classes,
             tmp_dir)

    Returns:
        Lightweight metadata dict::

            {
                "chunk_index": int,
                "chunk_start": int,
                "file_path":   str,
                "motif_count": int,
            }
    """
    from Utilities.disk_chunk_manager import DiskChunkManager
    from Utilities.nonbscanner import analyze_sequence

    (
        chunk_seq,
        seq_name,
        chunk_start,
        chunk_end,
        enabled_classes,
        tmp_dir,
        chunk_index,
    ) = args

    # Analyse the chunk (no large objects returned across process boundary)
    raw_motifs = analyze_sequence(
        sequence=chunk_seq,
        sequence_name=seq_name,
        use_fast_mode=True,
        enabled_classes=enabled_classes,
    )

    # Adjust positions to genome-global coordinates
    adjusted: List[Dict[str, Any]] = []
    for motif in raw_motifs:
        m = motif.copy()
        m["Start"] = motif.get("Start", 0) + chunk_start
        m["End"] = motif.get("End", 0) + chunk_start
        if "ID" in m:
            parts = m["ID"].split("_")
            if len(parts) >= 3:
                parts[-1] = str(m["Start"])
                m["ID"] = "_".join(parts)
        adjusted.append(m)

    # Write to disk; return only metadata
    manager = DiskChunkManager(base_dir=tmp_dir)
    # Override the auto-created subdir to use the shared tmp_dir directly
    import tempfile
    from pathlib import Path
    import csv

    file_path = Path(tmp_dir) / f"chunk_{chunk_index:04d}.csv"
    _fields = ["Class", "Subclass", "Start", "End", "Length", "Score", "ID", "Strand"]
    with open(file_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_fields, extrasaction="ignore")
        writer.writeheader()
        for m in adjusted:
            writer.writerow(m)

    return {
        "chunk_index": chunk_index,
        "chunk_start": chunk_start,
        "chunk_end": chunk_end,
        "file_path": str(file_path),
        "motif_count": len(adjusted),
    }


# ──────────────────────────────────────────────────────────────────────────────
# EXECUTOR
# ──────────────────────────────────────────────────────────────────────────────

class StreamlitSafeExecutor:
    """
    Streamlit Community Cloud compatible analysis executor.

    Selects an execution strategy based on sequence length and enforces a
    hard cap of ``max_workers = min(2, os.cpu_count())`` to avoid RAM and
    CPU exhaustion inside the Streamlit container.

    Usage::

        executor = StreamlitSafeExecutor(sequence_storage)
        results_storage = executor.run(
            seq_id=seq_id,
            progress_callback=lambda p: print(f"{p:.0f}%"),
            enabled_classes=["G-Quadruplex"],
        )
    """

    def __init__(
        self,
        sequence_storage,
        chunk_size: int = DEFAULT_CHUNK_SIZE,
        overlap: int = DEFAULT_OVERLAP,
        max_workers: Optional[int] = None,
    ):
        """
        Args:
            sequence_storage: ``UniversalSequenceStorage`` instance.
            chunk_size:       Chunk size in bp (default 50 000).
            overlap:          Overlap between consecutive chunks in bp (default 2 000).
            max_workers:      Override the default worker cap.
        """
        self.sequence_storage = sequence_storage
        self.chunk_size = chunk_size
        self.overlap = overlap
        self.max_workers = max_workers if max_workers is not None else MAX_WORKERS

        logger.info(
            f"StreamlitSafeExecutor ready "
            f"(max_workers={self.max_workers}, chunk={chunk_size:,}, overlap={overlap:,})"
        )

    # ------------------------------------------------------------------
    # PUBLIC API
    # ------------------------------------------------------------------

    def run(
        self,
        seq_id: str,
        progress_callback: Optional[Callable[[float], None]] = None,
        enabled_classes: Optional[List[str]] = None,
    ):
        """
        Analyse a stored sequence using the adaptive strategy.

        Args:
            seq_id:            Sequence identifier from ``UniversalSequenceStorage``.
            progress_callback: Optional ``callback(progress_pct: float)``.
            enabled_classes:   Motif classes to analyse (None = all).

        Returns:
            ``UniversalResultsStorage`` instance with results on disk.
        """
        meta = self.sequence_storage.get_metadata(seq_id)
        seq_length = meta["length"]
        strategy = self._select_strategy(seq_length)

        logger.info(
            f"StreamlitSafeExecutor: seq_length={seq_length:,} → strategy='{strategy}'"
        )

        if strategy == "direct":
            return self._run_direct(seq_id, meta, progress_callback, enabled_classes)
        elif strategy == "chunk_workers":
            return self._run_chunk_workers(
                seq_id, meta, progress_callback, enabled_classes
            )
        else:  # disk_streaming
            return self._run_disk_streaming(
                seq_id, meta, progress_callback, enabled_classes
            )

    # ------------------------------------------------------------------
    # STRATEGY SELECTION
    # ------------------------------------------------------------------

    @staticmethod
    def _select_strategy(seq_length: int) -> str:
        """Return 'direct', 'chunk_workers', or 'disk_streaming'."""
        if seq_length < THRESHOLD_DIRECT:
            return "direct"
        if seq_length < THRESHOLD_CHUNK:
            return "chunk_workers"
        return "disk_streaming"

    # ------------------------------------------------------------------
    # DIRECT (small sequences – no chunking)
    # ------------------------------------------------------------------

    def _run_direct(
        self,
        seq_id: str,
        meta: Dict[str, Any],
        progress_callback: Optional[Callable[[float], None]],
        enabled_classes: Optional[List[str]],
    ):
        from Utilities.disk_storage import UniversalResultsStorage
        from Utilities.nonbscanner import analyze_sequence

        seq_name = meta["name"]
        results_storage = UniversalResultsStorage(
            base_dir=str(self.sequence_storage.base_dir / "results"),
            seq_id=seq_id,
        )

        chunk = self.sequence_storage.get_sequence_chunk(
            seq_id, 0, meta["length"]
        )
        motifs = analyze_sequence(
            sequence=chunk,
            sequence_name=seq_name,
            use_fast_mode=True,
            enabled_classes=enabled_classes,
        )
        results_storage.append_batch(motifs)

        if progress_callback:
            progress_callback(100.0)

        logger.info(
            f"Direct analysis: {len(motifs)} motifs for '{seq_name}'"
        )
        return results_storage

    # ------------------------------------------------------------------
    # CHUNK WORKERS (medium sequences – 2-worker pool)
    # ------------------------------------------------------------------

    def _run_chunk_workers(
        self,
        seq_id: str,
        meta: Dict[str, Any],
        progress_callback: Optional[Callable[[float], None]],
        enabled_classes: Optional[List[str]],
    ):
        """2-worker ProcessPoolExecutor, workers write to disk."""
        import tempfile
        from pathlib import Path
        from Utilities.disk_storage import UniversalResultsStorage

        results_storage = UniversalResultsStorage(
            base_dir=str(self.sequence_storage.base_dir / "results"),
            seq_id=seq_id,
        )

        # Collect chunk args (only small scalars, not DataFrames)
        chunk_args: List[Tuple] = []
        for chunk_seq, chunk_start, chunk_end in self.sequence_storage.iter_chunks(
            seq_id, self.chunk_size, self.overlap
        ):
            idx = len(chunk_args)
            chunk_args.append((
                chunk_seq,
                f"{meta['name']}_chunk{idx}",
                chunk_start,
                chunk_end,
                enabled_classes,
                str(self.sequence_storage.base_dir / "chunks"),
                idx,
            ))

        # Ensure the shared chunk directory exists
        chunk_dir = self.sequence_storage.base_dir / "chunks"
        chunk_dir.mkdir(parents=True, exist_ok=True)

        total = len(chunk_args)
        logger.info(
            f"ChunkWorkers: processing {total} chunks "
            f"(chunk_size={self.chunk_size:,}, overlap={self.overlap:,}) "
            f"for '{meta['name']}' ({meta['length']:,} bp)"
        )
        completed = 0
        chunk_metadata: List[Optional[Dict[str, Any]]] = [None] * total

        try:
            with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                future_to_idx = {
                    executor.submit(_chunk_worker, args): args[6]
                    for args in chunk_args
                }
                for future in as_completed(future_to_idx):
                    idx = future_to_idx[future]
                    result_meta = future.result()  # raises on worker error
                    chunk_metadata[idx] = result_meta
                    completed += 1
                    if progress_callback:
                        progress_callback(completed / total * 100.0)
                    logger.info(
                        f"Chunk {idx + 1}/{total} complete: "
                        f"{result_meta['motif_count']} motifs at "
                        f"[{result_meta['chunk_start']:,}–{result_meta['chunk_end']:,}] "
                        f"({completed / total * 100:.1f}% overall)"
                    )

        except Exception as exc:
            logger.warning(
                f"Parallel chunk workers failed ({exc}); "
                "falling back to sequential"
            )
            return self._run_disk_streaming(
                seq_id, meta, progress_callback, enabled_classes
            )

        # Annotate each chunk_metadata entry with core_end for deduplication.
        # core_end is the exclusive boundary of the authoritative region:
        #   non-last chunks: core_end = chunk_end - overlap
        #   last chunk:      core_end = chunk_end  (full chunk is authoritative)
        for i, cm in enumerate(chunk_metadata):
            if cm is not None:
                is_last = (i == total - 1)
                cm["core_end"] = (
                    cm["chunk_end"] if is_last
                    else cm["chunk_end"] - self.overlap
                )

        # Stream chunk files into results_storage, deduplicate at boundaries
        self._merge_chunk_files_to_storage(chunk_metadata, results_storage)

        # Clean up chunk CSV files
        for cm in chunk_metadata:
            if cm and Path(cm["file_path"]).exists():
                Path(cm["file_path"]).unlink(missing_ok=True)

        gc.collect()
        return results_storage

    # ------------------------------------------------------------------
    # DISK STREAMING (large sequences – fully disk-backed)
    # ------------------------------------------------------------------

    def _run_disk_streaming(
        self,
        seq_id: str,
        meta: Dict[str, Any],
        progress_callback: Optional[Callable[[float], None]],
        enabled_classes: Optional[List[str]],
    ):
        """
        Sequential disk-backed processing for very large sequences.

        Each chunk is processed, written to disk, and immediately merged
        into the results storage.  Only one chunk is held in RAM at a time.

        Overlap deduplication uses the core_end boundary from ChunkGenerator:
        - non-last chunks: core_end = chunk_end - overlap
        - last chunk:      core_end = chunk_end
        Only motifs with Start < core_end are kept per chunk, ensuring each
        motif is counted exactly once in its authoritative (core) region.
        """
        from Utilities.disk_storage import UniversalResultsStorage
        from Utilities.nonbscanner import analyze_sequence
        from Utilities.overlap_deduplicator import OverlapDeduplicator

        seq_name = meta["name"]
        seq_length = meta["length"]

        results_storage = UniversalResultsStorage(
            base_dir=str(self.sequence_storage.base_dir / "results"),
            seq_id=seq_id,
        )

        # First pass: count total chunks for accurate progress reporting.
        # This iterates metadata only (reads chunk boundaries, not chunk content)
        # so memory usage is constant even for very large sequences.
        total_chunks = sum(
            1
            for _ in self.sequence_storage.iter_chunks(
                seq_id, self.chunk_size, self.overlap
            )
        )

        logger.info(
            f"DiskStreaming: processing {total_chunks} chunks "
            f"(chunk_size={self.chunk_size:,}, overlap={self.overlap:,}) "
            f"for '{seq_name}' ({seq_length:,} bp)"
        )

        dedup = OverlapDeduplicator()
        total_raw = 0
        total_kept = 0
        chunk_num = 0

        # Second pass: process one chunk at a time (constant RAM, no list materialisation)
        for chunk_seq, chunk_start, chunk_end in self.sequence_storage.iter_chunks(
            seq_id, self.chunk_size, self.overlap
        ):
            chunk_num += 1
            is_last = (chunk_num == total_chunks)
            # core_end: authoritative region is [chunk_start, core_end)
            # motifs with Start >= core_end are in the overlap zone and will be
            # re-detected (and kept) by the next chunk's core region
            core_end = chunk_end if is_last else chunk_end - self.overlap

            raw_motifs = analyze_sequence(
                sequence=chunk_seq,
                sequence_name=f"{seq_name}_chunk{chunk_num}",
                use_fast_mode=True,
                enabled_classes=enabled_classes,
            )

            # Adjust positions to genome-global coordinates
            globally_positioned_motifs: List[Dict[str, Any]] = []
            for motif in raw_motifs:
                m = motif.copy()
                m["Start"] = motif.get("Start", 0) + chunk_start
                m["End"] = motif.get("End", 0) + chunk_start
                if "ID" in m:
                    parts = m["ID"].split("_")
                    if len(parts) >= 3:
                        parts[-1] = str(m["Start"])
                        m["ID"] = "_".join(parts)
                globally_positioned_motifs.append(m)

            # Rigorous core-region filtering: keep only motifs in authoritative zone
            unique_motifs = dedup.filter_core(globally_positioned_motifs, core_end)
            total_raw += len(raw_motifs)
            total_kept += len(unique_motifs)

            results_storage.append_batch(unique_motifs)

            logger.info(
                f"Chunk {chunk_num}/{total_chunks} "
                f"[{chunk_start:,}–{chunk_end:,}] core_end={core_end:,}: "
                f"{len(raw_motifs)} raw → {len(unique_motifs)} kept "
                f"({chunk_num / total_chunks * 100:.1f}% done)"
            )

            if progress_callback:
                progress_callback(chunk_num / total_chunks * 100.0)

            # Free memory immediately
            del chunk_seq, raw_motifs, globally_positioned_motifs, unique_motifs
            gc.collect()

        logger.info(
            f"DiskStreaming complete: {total_kept} motifs kept "
            f"from {total_raw} raw detections across {total_chunks} chunks"
        )
        return results_storage

    # ------------------------------------------------------------------
    # HELPERS
    # ------------------------------------------------------------------

    @staticmethod
    def _merge_chunk_files_to_storage(
        chunk_metadata: List[Optional[Dict[str, Any]]],
        results_storage,
    ) -> None:
        """
        Read chunk CSV files in order and append unique motifs to results_storage.

        Performs rigorous core-region boundary deduplication using the pre-computed
        ``core_end`` value in each chunk's metadata entry (set by _run_chunk_workers):
        - non-last chunks: core_end = chunk_end - overlap
        - last chunk:      core_end = chunk_end

        Only motifs with Start < core_end are kept, ensuring each motif is counted
        exactly once in its authoritative (core) region without any cross-chunk
        seen-key bookkeeping.
        """
        import csv
        from pathlib import Path
        from Utilities.overlap_deduplicator import OverlapDeduplicator

        dedup = OverlapDeduplicator()
        total_raw = 0
        total_kept = 0
        total_chunks = len([cm for cm in chunk_metadata if cm is not None])

        for chunk_idx, cm in enumerate(chunk_metadata):
            if cm is None:
                continue

            file_path = Path(cm["file_path"])
            if not file_path.exists():
                logger.warning(f"Missing chunk file: {file_path}")
                continue

            # core_end was computed in _run_chunk_workers
            core_end = cm.get("core_end", cm.get("chunk_end", float("inf")))

            chunk_motifs: list = []
            with open(file_path, newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
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
                    chunk_motifs.append(dict(row))

            # Rigorous core-region filtering via OverlapDeduplicator
            unique_motifs = dedup.filter_core(chunk_motifs, int(core_end))
            total_raw += len(chunk_motifs)
            total_kept += len(unique_motifs)

            logger.info(
                f"Merge chunk {chunk_idx + 1}/{total_chunks} "
                f"[core_end={core_end:,}]: "
                f"{len(chunk_motifs)} raw → {len(unique_motifs)} kept"
            )

            results_storage.append_batch(unique_motifs)
            gc.collect()
