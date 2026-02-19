"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Detector Runner & Parallel Chunk Executor - Parallel Motif Detection         │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Two complementary classes that together drive parallel motif detection:

    DetectorRunner
        Runs all Non-B DNA motif detectors on a single sequence chunk.
        Writes results to a CSV file on disk and returns only lightweight
        metadata (file path, chunk bounds, motif count).
        Workers never return large objects across process boundaries.

    ParallelChunkExecutor
        Submits chunks to a ``ProcessPoolExecutor``, collects per-chunk
        metadata as futures complete, and yields results one at a time.
        Worker count is capped at ``MAX_WORKERS`` (≤ 2 by default) to stay
        within Streamlit Community Cloud limits.
        Falls back to sequential processing on any executor failure.

MEMORY GUARANTEES:
    - Only lightweight metadata (file path, ints) is returned by workers.
    - No full motif lists cross the process boundary.
    - One chunk at a time is in RAM on each worker.

USAGE::

    from Utilities.detector_runner import DetectorRunner, ParallelChunkExecutor
    from Utilities.chunk_generator import ChunkGenerator

    runner   = DetectorRunner(tmp_dir="/tmp/my_run")
    executor = ParallelChunkExecutor(detector_runner=runner, workers=2)

    gen = ChunkGenerator(sequence, chunk_size=50_000, overlap=2_000)
    for result_meta in executor.execute(gen.generate()):
        # result_meta: {"file", "chunk_start", "core_end", "chunk_index", ...}
        motifs = read_chunk_csv(result_meta["file"])
        ...
"""

from __future__ import annotations

import csv
import gc
import logging
import os
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Generator, Iterator, List, Optional

logger = logging.getLogger(__name__)

# Hard ceiling on worker count to protect Streamlit Community Cloud
MAX_WORKERS: int = min(2, os.cpu_count() or 1)

# Fieldnames for per-chunk CSV files
_CSV_FIELDS = ["Class", "Subclass", "Start", "End", "Length", "Score", "ID", "Strand"]


# ──────────────────────────────────────────────────────────────────────────────
# MODULE-LEVEL WORKER (must be picklable for ProcessPoolExecutor)
# ──────────────────────────────────────────────────────────────────────────────

def _worker_run_chunk(args: tuple) -> Dict[str, Any]:
    """
    Process a single chunk and write results to disk.

    This is a module-level function so it is picklable by
    ``ProcessPoolExecutor`` without serialising large objects.

    Args:
        args: Tuple of
            (chunk_seq, seq_name, chunk_start, chunk_end, core_end,
             enabled_classes, tmp_dir, chunk_index)

    Returns:
        Lightweight metadata dict::

            {
                "chunk_index" : int,
                "chunk_start" : int,
                "chunk_end"   : int,
                "core_end"    : int,
                "file"        : str,
                "motif_count" : int,
            }
    """
    from Utilities.nonbscanner import analyze_sequence

    (
        chunk_seq,
        seq_name,
        chunk_start,
        chunk_end,
        core_end,
        enabled_classes,
        tmp_dir,
        chunk_index,
    ) = args

    # Detect motifs (chunk-local coordinates)
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
        m["Start"] = int(motif.get("Start", 0)) + chunk_start
        m["End"] = int(motif.get("End", 0)) + chunk_start
        if "ID" in m:
            parts = m["ID"].split("_")
            if len(parts) >= 3:
                parts[-1] = str(m["Start"])
                m["ID"] = "_".join(parts)
        adjusted.append(m)

    # Write to disk; return only metadata (no large objects across process boundary)
    file_path = Path(tmp_dir) / f"chunk_{chunk_index:04d}.csv"
    with open(file_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for m in adjusted:
            writer.writerow(m)

    return {
        "chunk_index": chunk_index,
        "chunk_start": chunk_start,
        "chunk_end": chunk_end,
        "core_end": core_end,
        "file": str(file_path),
        "motif_count": len(adjusted),
    }


# ──────────────────────────────────────────────────────────────────────────────
# DETECTOR RUNNER
# ──────────────────────────────────────────────────────────────────────────────

class DetectorRunner:
    """
    Run all Non-B DNA motif detectors on a single sequence chunk.

    Writes chunk results to a CSV file on disk and returns lightweight
    metadata.  Designed to be called from within a worker process so that
    no large objects need to cross the process boundary.

    Usage::

        runner = DetectorRunner(tmp_dir="/tmp/nbdna_run")
        meta   = runner.run(chunk_data, chunk_index=0)
        # meta["file"] → path to CSV; meta["core_end"] → dedup boundary
    """

    def __init__(
        self,
        tmp_dir: Optional[str] = None,
        enabled_classes: Optional[List[str]] = None,
        seq_name: str = "seq",
    ):
        """
        Args:
            tmp_dir:         Directory for temporary CSV files.
                             Created automatically if ``None``.
            enabled_classes: Motif classes to analyse (``None`` = all).
            seq_name:        Base name used in motif IDs.
        """
        if tmp_dir is None:
            self._tmp_dir = tempfile.mkdtemp(prefix="nonbdna_detector_")
            self._owns_tmp = True
        else:
            self._tmp_dir = tmp_dir
            Path(tmp_dir).mkdir(parents=True, exist_ok=True)
            self._owns_tmp = False

        self.enabled_classes = enabled_classes
        self.seq_name = seq_name

    def run(
        self,
        chunk_data: Dict[str, Any],
        chunk_index: int = 0,
    ) -> Dict[str, Any]:
        """
        Detect motifs in *chunk_data* and write results to disk.

        Args:
            chunk_data:  Dict produced by :class:`~Utilities.chunk_generator.ChunkGenerator`.
                         Must contain ``"sequence"``, ``"start"``, ``"end"``,
                         ``"core_end"``.
            chunk_index: Zero-based chunk index used to name the output file.

        Returns:
            Lightweight metadata dict::

                {
                    "chunk_index" : int,
                    "chunk_start" : int,
                    "chunk_end"   : int,
                    "core_end"    : int,
                    "file"        : str,   # path to CSV file
                    "motif_count" : int,
                }
        """
        args = (
            chunk_data["sequence"],
            f"{self.seq_name}_chunk{chunk_index}",
            chunk_data["start"],
            chunk_data["end"],
            chunk_data["core_end"],
            self.enabled_classes,
            self._tmp_dir,
            chunk_index,
        )
        return _worker_run_chunk(args)

    def cleanup(self) -> None:
        """Remove the temporary directory (only if created by this instance)."""
        import shutil
        if self._owns_tmp and Path(self._tmp_dir).exists():
            shutil.rmtree(self._tmp_dir, ignore_errors=True)
            logger.info(f"DetectorRunner: cleaned up {self._tmp_dir}")


# ──────────────────────────────────────────────────────────────────────────────
# PARALLEL CHUNK EXECUTOR
# ──────────────────────────────────────────────────────────────────────────────

class ParallelChunkExecutor:
    """
    Execute chunk detection in parallel and yield per-chunk metadata.

    Uses a ``ProcessPoolExecutor`` with a hard cap of :attr:`MAX_WORKERS`
    workers.  If the executor fails for any reason the executor falls back
    to sequential (in-process) execution automatically.

    Usage::

        runner   = DetectorRunner(tmp_dir="/tmp/run")
        executor = ParallelChunkExecutor(detector_runner=runner, workers=2)

        for result_meta in executor.execute(chunk_generator.generate()):
            # result_meta["file"] – path to per-chunk CSV
            # result_meta["core_end"] – dedup boundary
            ...
    """

    def __init__(self, detector_runner: DetectorRunner, workers: Optional[int] = None):
        """
        Args:
            detector_runner: :class:`DetectorRunner` instance that defines
                             ``tmp_dir``, ``enabled_classes``, and ``seq_name``.
            workers:         Number of parallel workers.  Capped at
                             :attr:`MAX_WORKERS` (≤ 2 on Streamlit Cloud).
        """
        self.detector_runner = detector_runner
        requested = workers if workers is not None else MAX_WORKERS
        self.workers = min(requested, MAX_WORKERS)

        logger.info(
            f"ParallelChunkExecutor ready (workers={self.workers})"
        )

    def execute(
        self,
        chunks: Iterator[Dict[str, Any]],
    ) -> Generator[Dict[str, Any], None, None]:
        """
        Process all *chunks* and yield result metadata as chunks complete.

        Results are yielded in **completion order** (not necessarily genome
        order) to maximise throughput.  Callers that require ordered results
        should sort by ``"chunk_index"`` before merging.

        Falls back to sequential processing if the
        ``ProcessPoolExecutor`` raises any exception.

        Args:
            chunks: Iterator of chunk dicts from
                    :class:`~Utilities.chunk_generator.ChunkGenerator`.

        Yields:
            Lightweight metadata dicts (see :meth:`DetectorRunner.run`).
        """
        # Materialise chunks to build worker arg tuples
        all_args: List[tuple] = []
        for idx, chunk in enumerate(chunks):
            all_args.append((
                chunk["sequence"],
                f"{self.detector_runner.seq_name}_chunk{idx}",
                chunk["start"],
                chunk["end"],
                chunk["core_end"],
                self.detector_runner.enabled_classes,
                self.detector_runner._tmp_dir,
                idx,
            ))

        if not all_args:
            return

        try:
            with ProcessPoolExecutor(max_workers=self.workers) as executor:
                futures = {
                    executor.submit(_worker_run_chunk, args): args[7]
                    for args in all_args
                }
                for future in as_completed(futures):
                    result = future.result()  # propagates worker exceptions
                    logger.debug(
                        f"Chunk {result['chunk_index']} done: "
                        f"{result['motif_count']} motifs"
                    )
                    yield result
                    gc.collect()

        except Exception as exc:
            logger.warning(
                f"ParallelChunkExecutor: parallel execution failed "
                f"({type(exc).__name__}: {exc}); falling back to sequential"
            )
            # Sequential fallback: run each chunk in-process
            for args in all_args:
                result = _worker_run_chunk(args)
                logger.debug(
                    f"Sequential fallback chunk {result['chunk_index']}: "
                    f"{result['motif_count']} motifs"
                )
                yield result
                gc.collect()
