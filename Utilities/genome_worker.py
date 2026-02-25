"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Genome Worker – chromosome-level ProcessPoolExecutor worker                  │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.2            │
└──────────────────────────────────────────────────────────────────────────────┘

Each worker:
  • Initialises detector classes locally (no shared state across processes).
  • Processes one chromosome/sequence end-to-end.
  • For large chromosomes (> LARGE_CHR_THRESHOLD) splits into overlapping chunks,
    processes them, then deduplicates boundary motifs.
  • Writes per-chromosome results to a Parquet file (pyarrow engine).
  • Returns a lightweight summary tuple – no large objects are passed back
    between processes.

This module must remain importable from worker processes (it is a plain module,
not defined inside a Jupyter cell) so that ``spawn``-based multiprocessing works
on all platforms.
"""
from __future__ import annotations

import os
import time
import warnings
from typing import Any, Dict, List, Optional, Tuple

warnings.filterwarnings("ignore")

# ── Optional Parquet writer ────────────────────────────────────────────────────
try:
    import pyarrow as _pa
    import pyarrow.parquet as _pq
    _HAS_PYARROW = True
except ImportError:
    _HAS_PYARROW = False

import pandas as pd

# ── Expected output columns (must match rest of pipeline) ─────────────────────
_OUTPUT_COLS = [
    "Class", "Subclass", "Start", "End", "Length",
    "Score", "Strand", "Sequence_Name", "Source_File", "File_Type",
]

_COL_DEFAULTS: Dict[str, Any] = {
    "Class": "Unknown", "Subclass": "Other",
    "Start": 0, "End": 0, "Length": 0,
    "Score": 0.0, "Strand": "+",
    "Sequence_Name": "", "Source_File": "", "File_Type": "",
}

import numpy as np


def _chunk_sequence(
    seq: str,
    chunk_size: int,
    overlap: int,
) -> List[Tuple[int, int]]:
    """Return (start, end) pairs that tile *seq* with the given overlap."""
    positions: List[Tuple[int, int]] = []
    start = 0
    n = len(seq)
    while start < n:
        end = min(start + chunk_size, n)
        positions.append((start, end))
        if end == n:
            break
        start += chunk_size - overlap
    return positions


def _dedup_boundary(motifs: List[Dict], chunk_starts: List[int], overlap: int) -> List[Dict]:
    """Remove duplicate motifs produced in the overlap regions between chunks.

    A motif is a duplicate if another motif with the same (Class, Subclass,
    Start, End) already exists in the merged list.  The one with the higher
    Score is kept.
    """
    seen: Dict[Tuple, Dict] = {}
    for m in motifs:
        key = (m.get("Class", ""), m.get("Subclass", ""), m.get("Start", 0), m.get("End", 0))
        if key not in seen or m.get("Score", 0) > seen[key].get("Score", 0):
            seen[key] = m
    return sorted(seen.values(), key=lambda x: x.get("Start", 0))


def _run_detectors_locally(
    seq: str,
    seq_name: str,
    enabled_classes: Optional[List[str]],
) -> List[Dict]:
    """Initialise a fresh NonBScanner and scan *seq*.

    Detectors are created inside this call so no state is shared between
    worker processes.
    """
    from Utilities.nonbscanner import NonBScanner, CLASS_TO_DETECTOR

    scanner = NonBScanner(enable_all_detectors=True)
    return scanner.analyze_sequence(
        seq,
        seq_name,
        enabled_classes=enabled_classes,
        use_parallel_detectors=True,
    )


def process_chromosome(
    args: Tuple,
) -> Tuple[str, Optional[str], int, float]:
    """Top-level worker executed by ``ProcessPoolExecutor``.

    Parameters
    ----------
    args:
        Tuple of::

            (seq_name, seq, source_file, file_type,
             large_chr_threshold, genome_chunk_size, genome_chunk_overlap,
             enabled_classes, parquet_dir, chunk_size, chunk_overlap)

    Returns
    -------
    (seq_name, parquet_path_or_None, n_motifs, elapsed_seconds)
        *parquet_path_or_None* is ``None`` when the sequence has no motifs.
    """
    (
        seq_name,
        seq,
        source_file,
        file_type,
        large_chr_threshold,
        genome_chunk_size,
        genome_chunk_overlap,
        enabled_classes,
        parquet_dir,
        chunk_size,
        chunk_overlap,
    ) = args

    t0 = time.perf_counter()
    seq_len = len(seq)

    if seq_len < 10:
        return seq_name, None, 0, 0.0

    # ── Split large chromosomes into sub-chunks processed sequentially ─────────
    if seq_len > large_chr_threshold:
        tile_positions = _chunk_sequence(seq, genome_chunk_size, genome_chunk_overlap)
        all_motifs: List[Dict] = []
        chunk_starts: List[int] = [s for s, _ in tile_positions]
        for c_start, c_end in tile_positions:
            chunk_seq = seq[c_start:c_end]
            chunk_motifs = _run_detectors_locally(chunk_seq, seq_name, enabled_classes)
            for m in chunk_motifs:
                m["Start"] += c_start
                m["End"] += c_start
            all_motifs.extend(chunk_motifs)
        motifs = _dedup_boundary(all_motifs, chunk_starts, genome_chunk_overlap)
    else:
        motifs = _run_detectors_locally(seq, seq_name, enabled_classes)

    elapsed = time.perf_counter() - t0

    if not motifs:
        return seq_name, None, 0, elapsed

    # ── Build lightweight DataFrame ────────────────────────────────────────────
    # Use only scalar column access to avoid heavy dict-in-loop overhead.
    records: List[Tuple] = []
    for m in motifs:
        length = m.get("Length", 0)
        start  = m.get("Start",  0)
        end    = m.get("End",    0)
        if length == 0:
            length = max(0, end - start)
        records.append((
            m.get("Class",         "Unknown"),
            m.get("Subclass",      "Other"),
            start,
            end,
            length,
            m.get("Score",         0.0),
            m.get("Strand",        "+"),
            seq_name,
            source_file,
            file_type,
        ))

    df = pd.DataFrame.from_records(records, columns=_OUTPUT_COLS)

    # ── Write to Parquet ───────────────────────────────────────────────────────
    safe_name = seq_name.replace("/", "_").replace("\\", "_")[:80]
    parquet_path = os.path.join(parquet_dir, f"{safe_name}.parquet")

    if _HAS_PYARROW:
        table = _pa.Table.from_pandas(df, preserve_index=False)
        _pq.write_table(table, parquet_path, compression="snappy")
    else:
        df.to_parquet(parquet_path, index=False)

    return seq_name, parquet_path, len(df), elapsed
