"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ ChunkExecutor - Chunk-Based Parallel Detection Engine                        │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.2            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Replaces detector-level parallelism with chunk-level parallel single-pass
    execution.  Each worker process receives one genome chunk and runs all
    enabled detectors sequentially on a shared ``SequenceContext``, eliminating
    the nested ``ThreadPoolExecutor`` inside worker processes.

    Execution model::

        Main Process
            ↓
        Generate Chunks (with overlap)
            ↓
        ProcessPoolExecutor (max_workers = min(4, cpu_count))
            ↓
        _chunk_worker()  [per chunk]
            ↓
        SequenceContext (uppercase once, cache length)
            ↓
        All detectors run *sequentially* on shared context
            ↓
        Return motif results + per-detector timings
            ↓
        Main: merge, deduplicate, collect telemetry

    Worker outputs are written to CSV files (tiny cross-process footprint).

    Adaptive worker cap:
        workers = min(4, os.cpu_count())
    Memory-aware cap is applied via psutil when available.

KEY DIFFERENCES FROM StreamlitSafeExecutor
    - max_workers raised to min(4, cpu_count)  (was hard-capped at 2)
    - No ThreadPoolExecutor inside worker (was used for parallel detectors)
    - SequenceContext preprocessing happens exactly once per chunk
    - Per-detector timings included in every worker return payload
    - Rich progress_callback receives a telemetry dict, not just a float

BACKWARD COMPATIBILITY
    ``ChunkExecutor.run()`` returns the same ``UniversalResultsStorage`` type
    as ``StreamlitSafeExecutor.run()`` so it is a drop-in replacement.
"""

from __future__ import annotations

import csv
import gc
import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from time import perf_counter
from typing import Any, Callable, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# Adaptive worker cap: balance speed vs RAM for Community Cloud
_DEFAULT_MAX_WORKERS: int = min(4, os.cpu_count() or 1)

# Default chunk parameters (kept in sync with StreamlitSafeExecutor defaults)
_DEFAULT_CHUNK_SIZE: int = 50_000
_DEFAULT_OVERLAP: int = 2_000

# CSV fields written by the worker
_MOTIF_FIELDS = ["Class", "Subclass", "Start", "End", "Length", "Score", "ID", "Strand"]


# ──────────────────────────────────────────────────────────────────────────────
# MODULE-LEVEL WORKER (must be picklable for ProcessPoolExecutor)
# ──────────────────────────────────────────────────────────────────────────────

def _chunk_worker(
    args: Tuple[str, str, int, int, Optional[List[str]], int, str],
) -> Dict[str, Any]:
    """
    Process one genome chunk with all detectors running sequentially.

    This is a *module-level* function so it can be serialised by
    ``ProcessPoolExecutor`` without pickling any large objects.

    Each detector is timed individually via the NonBScanner progress_callback.
    Results are written to a CSV file; only lightweight metadata is returned
    across the process boundary.

    Args:
        args: Tuple of
            (chunk_seq, seq_name, chunk_start, chunk_end, enabled_classes,
             chunk_index, tmp_dir)

    Returns:
        Dict with keys::

            {
                "chunk_index":      int,
                "chunk_start":      int,
                "chunk_end":        int,
                "file_path":        str,
                "motif_count":      int,
                "elapsed":          float,   # total wall-clock for this chunk
                "bp_processed":     int,
                "detector_timings": {detector_name: elapsed, ...},
            }
    """
    from Utilities.core.sequence_context import SequenceContext
    from Utilities.nonbscanner import _get_cached_scanner

    chunk_seq, seq_name, chunk_start, chunk_end, enabled_classes, chunk_index, tmp_dir = args

    t0 = perf_counter()

    # Shared preprocessing: uppercase once, cache length
    context = SequenceContext(chunk_seq, seq_name)

    # Collect per-detector timings via the existing progress_callback hook
    detector_timings: Dict[str, float] = {}

    def _timing_cb(
        detector_name: str,
        _completed: int,
        _total: int,
        elapsed: float,
        _motif_count: int,
    ) -> None:
        detector_timings[detector_name] = elapsed

    # Run ALL detectors SEQUENTIALLY (no nested ThreadPoolExecutor)
    raw_motifs = _get_cached_scanner().analyze_sequence(
        context.sequence,      # already uppercase – avoids per-detector .upper()
        context.name,
        progress_callback=_timing_cb,
        enabled_classes=enabled_classes,
        use_parallel_detectors=False,   # single-pass sequential within chunk
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

    t1 = perf_counter()

    # Write to disk; return only lightweight metadata across process boundary
    file_path = Path(tmp_dir) / f"cexec_{chunk_index:04d}.csv"
    with open(file_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_MOTIF_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for m in adjusted:
            writer.writerow(m)

    return {
        "chunk_index": chunk_index,
        "chunk_start": chunk_start,
        "chunk_end": chunk_end,
        "file_path": str(file_path),
        "motif_count": len(adjusted),
        "elapsed": t1 - t0,
        "bp_processed": context.length,
        "detector_timings": detector_timings,
    }


# ──────────────────────────────────────────────────────────────────────────────
# CHUNK EXECUTOR
# ──────────────────────────────────────────────────────────────────────────────

class ChunkExecutor:
    """
    Chunk-parallel motif detection engine with real-time performance telemetry.

    Drop-in replacement for ``StreamlitSafeExecutor`` with the following
    improvements:

    * Workers raised to ``min(4, cpu_count)``  (was 2)
    * Detectors run **sequentially** inside each worker → no nested parallelism
    * ``SequenceContext`` centralises preprocessing (uppercase once per chunk)
    * ``progress_callback`` receives a rich telemetry dict (see below)
    * ``get_performance_summary()`` available after ``run()`` completes

    Progress callback payload::

        {
            "stage":        "detection",
            "chunk_id":     int,
            "detector":     str,        # from detector_timings merge
            "elapsed":      float,      # chunk wall-clock
            "bp_processed": int,
            "throughput":   float,      # bp/sec for this chunk
            "memory_mb":    float,
            "progress_pct": float,      # 0–100
        }

    If the callable accepts a single float argument (legacy API), only
    ``progress_pct`` is forwarded.

    Usage::

        from Utilities.disk_storage import UniversalSequenceStorage
        from Utilities.core.chunk_executor import ChunkExecutor

        storage = UniversalSequenceStorage()
        seq_id  = storage.save_sequence(dna, "chr1")

        executor = ChunkExecutor(storage)
        results_storage = executor.run(seq_id, progress_callback=cb,
                                       enabled_classes=["G-Quadruplex"])
        print(executor.get_performance_summary())
    """

    def __init__(
        self,
        sequence_storage,
        chunk_size: int = _DEFAULT_CHUNK_SIZE,
        overlap: int = _DEFAULT_OVERLAP,
        max_workers: Optional[int] = None,
    ) -> None:
        """
        Args:
            sequence_storage: ``UniversalSequenceStorage`` instance.
            chunk_size:       Chunk size in bp (default 50 000).
            overlap:          Overlap between consecutive chunks in bp (default 2 000).
            max_workers:      Override the worker cap (default ``min(4, cpu_count)``).
        """
        self.sequence_storage = sequence_storage
        self.chunk_size = chunk_size
        self.overlap = overlap
        self.max_workers = max_workers if max_workers is not None else _DEFAULT_MAX_WORKERS

        from Utilities.core.performance_monitor import PerformanceMonitor
        self._monitor = PerformanceMonitor()

        logger.info(
            f"ChunkExecutor ready (max_workers={self.max_workers}, "
            f"chunk={chunk_size:,}, overlap={overlap:,})"
        )

    # ------------------------------------------------------------------
    # PUBLIC API
    # ------------------------------------------------------------------

    def run(
        self,
        seq_id: str,
        progress_callback: Optional[Callable] = None,
        enabled_classes: Optional[List[str]] = None,
    ):
        """
        Analyse a stored sequence using chunk-parallel single-pass execution.

        Args:
            seq_id:            Sequence ID from ``UniversalSequenceStorage``.
            progress_callback: Optional callback.  Receives a telemetry dict
                               *or* a plain ``float`` (legacy) – auto-detected.
            enabled_classes:   Motif classes to analyse (``None`` = all).

        Returns:
            ``UniversalResultsStorage`` instance with results on disk.
        """
        from Utilities.disk_storage import UniversalResultsStorage

        self._monitor.start()
        meta = self.sequence_storage.get_metadata(seq_id)
        seq_length = meta["length"]

        logger.info(
            f"ChunkExecutor: '{meta['name']}' length={seq_length:,} bp"
        )

        results_storage = UniversalResultsStorage(
            base_dir=str(self.sequence_storage.base_dir / "results"),
            seq_id=seq_id,
        )

        # Ensure the shared chunk temp directory exists
        chunk_dir = self.sequence_storage.base_dir / "cexec_chunks"
        chunk_dir.mkdir(parents=True, exist_ok=True)

        # Build list of chunk args (scalars only – no large objects sent to workers)
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
                idx,
                str(chunk_dir),
            ))

        total = len(chunk_args)
        if total == 0:
            logger.warning("ChunkExecutor: no chunks generated – empty sequence?")
            return results_storage

        logger.info(f"ChunkExecutor: {total} chunks → {self.max_workers} workers")

        chunk_metadata: List[Optional[Dict[str, Any]]] = [None] * total
        completed = 0

        t_detection_start = perf_counter()

        try:
            with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                future_to_idx = {
                    executor.submit(_chunk_worker, args): args[5]   # args[5] = chunk_index
                    for args in chunk_args
                }
                for future in as_completed(future_to_idx):
                    idx = future_to_idx[future]
                    result_meta = future.result()   # raises on worker error
                    chunk_metadata[idx] = result_meta
                    completed += 1

                    # Record chunk telemetry
                    self._monitor.record_chunk(
                        chunk_id=idx,
                        elapsed=result_meta["elapsed"],
                        motif_count=result_meta["motif_count"],
                        bp_count=result_meta["bp_processed"],
                    )
                    # Record per-detector timings for this chunk
                    for det_name, det_elapsed in result_meta.get(
                        "detector_timings", {}
                    ).items():
                        self._monitor.record_detector(
                            chunk_id=idx,
                            detector_name=det_name,
                            elapsed=det_elapsed,
                            bp_count=result_meta["bp_processed"],
                        )

                    # Snapshot memory
                    mem_mb = self._monitor.snapshot_memory()
                    progress_pct = completed / total * 100.0

                    if progress_callback is not None:
                        self._emit_progress(
                            cb=progress_callback,
                            chunk_id=idx,
                            elapsed=result_meta["elapsed"],
                            bp_processed=result_meta["bp_processed"],
                            progress_pct=progress_pct,
                            memory_mb=mem_mb,
                        )

                    logger.debug(
                        f"ChunkExecutor: chunk {idx} done "
                        f"({result_meta['motif_count']} motifs, "
                        f"{result_meta['elapsed']:.3f}s)"
                    )

        except Exception as exc:
            logger.warning(
                f"ChunkExecutor: parallel execution failed ({exc}); "
                "falling back to sequential"
            )
            return self._run_sequential(
                seq_id, meta, results_storage, chunk_dir,
                progress_callback, enabled_classes, total,
            )

        t_detection_end = perf_counter()
        self._monitor.record_stage("detection", t_detection_end - t_detection_start)

        # Merge chunk files → results_storage with boundary deduplication
        t_filter_start = perf_counter()
        self._merge_chunks(chunk_metadata, results_storage)
        self._monitor.record_stage("filtering", perf_counter() - t_filter_start)

        # Clean up temp CSV files
        for cm in chunk_metadata:
            if cm and Path(cm["file_path"]).exists():
                Path(cm["file_path"]).unlink(missing_ok=True)

        gc.collect()

        stats = results_storage.get_summary_stats()
        logger.info(
            f"ChunkExecutor: complete – {stats['total_count']} motifs\n"
            + self._monitor.format_summary()
        )
        return results_storage

    def get_performance_summary(self) -> Dict[str, Any]:
        """Return the performance summary dict from the last ``run()`` call."""
        return self._monitor.get_summary()

    # ------------------------------------------------------------------
    # SEQUENTIAL FALLBACK
    # ------------------------------------------------------------------

    def _run_sequential(
        self,
        seq_id: str,
        meta: Dict[str, Any],
        results_storage,
        chunk_dir: Path,
        progress_callback: Optional[Callable],
        enabled_classes: Optional[List[str]],
        total_chunks: int,
    ):
        """
        Sequential fallback when ProcessPoolExecutor is unavailable.

        Processes each chunk in the main process, still using SequenceContext
        and per-detector timing.
        """
        from Utilities.core.sequence_context import SequenceContext
        from Utilities.nonbscanner import _get_cached_scanner

        seq_name = meta["name"]
        seen_keys: set = set()
        chunk_num = 0

        for chunk_seq, chunk_start, chunk_end in self.sequence_storage.iter_chunks(
            seq_id, self.chunk_size, self.overlap
        ):
            chunk_num += 1
            t0 = perf_counter()

            context = SequenceContext(chunk_seq, f"{seq_name}_chunk{chunk_num}")
            detector_timings: Dict[str, float] = {}

            def _timing_cb(det, _c, _t, elapsed, _m):  # noqa: E731
                detector_timings[det] = elapsed

            raw_motifs = _get_cached_scanner().analyze_sequence(
                context.sequence,
                context.name,
                progress_callback=_timing_cb,
                enabled_classes=enabled_classes,
                use_parallel_detectors=False,
            )

            unique: List[Dict[str, Any]] = []
            for motif in raw_motifs:
                m = motif.copy()
                m["Start"] = motif.get("Start", 0) + chunk_start
                m["End"] = motif.get("End", 0) + chunk_start
                if "ID" in m:
                    parts = m["ID"].split("_")
                    if len(parts) >= 3:
                        parts[-1] = str(m["Start"])
                        m["ID"] = "_".join(parts)
                key = (m.get("Class", ""), m.get("Subclass", ""), m["Start"], m["End"])
                if key not in seen_keys:
                    unique.append(m)
                    if m["Start"] >= chunk_end - self.overlap:
                        seen_keys.add(key)

            results_storage.append_batch(unique)

            elapsed = perf_counter() - t0
            self._monitor.record_chunk(
                chunk_id=chunk_num - 1,
                elapsed=elapsed,
                motif_count=len(unique),
                bp_count=context.length,
            )
            for det_name, det_elapsed in detector_timings.items():
                self._monitor.record_detector(
                    chunk_id=chunk_num - 1,
                    detector_name=det_name,
                    elapsed=det_elapsed,
                    bp_count=context.length,
                )

            mem_mb = self._monitor.snapshot_memory()
            progress_pct = chunk_num / total_chunks * 100.0

            if progress_callback is not None:
                self._emit_progress(
                    cb=progress_callback,
                    chunk_id=chunk_num - 1,
                    elapsed=elapsed,
                    bp_processed=context.length,
                    progress_pct=progress_pct,
                    memory_mb=mem_mb,
                )

            del chunk_seq, raw_motifs, unique
            gc.collect()

        return results_storage

    # ------------------------------------------------------------------
    # MERGE + DEDUPLICATION
    # ------------------------------------------------------------------

    @staticmethod
    def _merge_chunks(
        chunk_metadata: List[Optional[Dict[str, Any]]],
        results_storage,
    ) -> None:
        """
        Read per-chunk CSV files in order, deduplicate at boundaries, and
        append unique motifs to ``results_storage``.
        """
        seen_keys: set = set()

        for cm in chunk_metadata:
            if cm is None:
                continue
            file_path = Path(cm["file_path"])
            if not file_path.exists():
                logger.warning(f"ChunkExecutor: missing chunk file {file_path}")
                continue

            unique_motifs: List[Dict[str, Any]] = []
            chunk_end = cm.get("chunk_end", 0)

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
                    key = (
                        row.get("Class", ""),
                        row.get("Subclass", ""),
                        row.get("Start", 0),
                        row.get("End", 0),
                    )
                    if key not in seen_keys:
                        unique_motifs.append(dict(row))
                        if row.get("Start", 0) >= chunk_end - _DEFAULT_OVERLAP:
                            seen_keys.add(key)

            results_storage.append_batch(unique_motifs)
            gc.collect()

    # ------------------------------------------------------------------
    # PROGRESS HELPER
    # ------------------------------------------------------------------

    @staticmethod
    def _emit_progress(
        cb: Callable,
        chunk_id: int,
        elapsed: float,
        bp_processed: int,
        progress_pct: float,
        memory_mb: float,
    ) -> None:
        """
        Fire ``cb`` with either a rich telemetry dict or a plain float,
        depending on what the caller provided.
        """
        throughput = bp_processed / elapsed if elapsed > 0 else 0.0
        telemetry = {
            "stage": "detection",
            "chunk_id": chunk_id,
            "elapsed": elapsed,
            "bp_processed": bp_processed,
            "throughput": throughput,
            "memory_mb": memory_mb,
            "progress_pct": progress_pct,
        }
        try:
            cb(telemetry)
        except TypeError:
            # Legacy callback expects a single float
            try:
                cb(progress_pct)
            except Exception as exc:
                logger.debug(f"ChunkExecutor: progress_callback error: {exc}")
        except Exception as exc:
            logger.debug(f"ChunkExecutor: progress_callback error: {exc}")
