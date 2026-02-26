"""SMID: Seed-based Motif Identification and Dispatch layer.

SMID is a non-invasive preprocessing orchestration layer that reduces the
genome search space before invoking the full-resolution Non-B DNA detectors.

Workflow
--------
1. SMIDSeedEngine scans the chromosome for lightweight seed patterns.
2. Seed hits are clustered with linear_1d_clustering (O(n), no sklearn).
3. Each cluster is extended and overlapping windows are merged.
4. SMIDDispatcher runs the appropriate detector on each window only.
5. Results are coordinate-translated back to chromosome space.

Nothing inside the detector classes is modified.

Quick-start
-----------
    from SMID import smid_pipeline
    from Detectors import ZDNADetector, APhilicDetector   # unchanged

    results = smid_pipeline(chromosome_sequence, sequence_name="chr1")
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Dict, List, Tuple

from SMID.seed_engine import SMIDSeedEngine
from SMID.clustering import linear_1d_clustering, extend_and_merge
from SMID.dispatcher import SMIDDispatcher
from SMID.seed_database import SMID_SEEDS


# ---------------------------------------------------------------------------
# Helper: build per-class seed parameter lookup
# ---------------------------------------------------------------------------

def _build_class_params(seeds: List[dict]) -> Dict[str, Dict[str, int]]:
    """Return {class_name: {epsilon, min_samples, max_extension}} from seeds."""
    params: Dict[str, Dict[str, int]] = {}
    for seed in seeds:
        for cls in seed["classes"]:
            # Use the strictest (smallest) extension / most permissive epsilon
            # when a class appears in multiple seeds.  The conservative choice
            # is to take the *largest* epsilon (most merging) and the *largest*
            # extension so we never miss hits.
            if cls not in params:
                params[cls] = {
                    "epsilon": seed["epsilon"],
                    "min_samples": seed["min_samples"],
                    "max_extension": seed["max_extension"],
                }
            else:
                params[cls]["epsilon"] = max(params[cls]["epsilon"], seed["epsilon"])
                params[cls]["max_extension"] = max(
                    params[cls]["max_extension"], seed["max_extension"]
                )
    return params


# Pre-compute once at import time
_CLASS_PARAMS: Dict[str, Dict[str, int]] = _build_class_params(SMID_SEEDS)


# ---------------------------------------------------------------------------
# Public pipeline entry point
# ---------------------------------------------------------------------------

def smid_pipeline(
    chromosome_seq: str,
    sequence_name: str = "seq",
    class_to_detector_map: Dict[str, Any] | None = None,
    max_workers: int | None = None,
) -> List[Dict[str, Any]]:
    """End-to-end SMID pipeline for one chromosome/contig.

    Parameters
    ----------
    chromosome_seq : str
        Full chromosome sequence (any case).
    sequence_name : str, optional
        Label for result records.
    class_to_detector_map : dict, optional
        Override the default detector class mapping (for testing).
    max_workers : int, optional
        Number of parallel worker processes for window-level parallelism.
        ``None`` uses the system default (CPU count).

    Returns
    -------
    list of dict
        Motif records from all detectors, coordinates in chromosome space.
    """
    seq_len = len(chromosome_seq)

    # Step 1 – seed scan
    engine = SMIDSeedEngine()
    seed_hits = engine.scan(chromosome_seq)

    # Step 2 & 3 – cluster and extend per class
    windows_by_class: Dict[str, List[Tuple[int, int]]] = {}
    for cls, positions in seed_hits.items():
        if len(positions) == 0:
            continue
        params = _CLASS_PARAMS.get(cls, {"epsilon": 150, "min_samples": 2, "max_extension": 300})
        clusters = linear_1d_clustering(positions, params["epsilon"], params["min_samples"])
        windows = extend_and_merge(clusters, params["max_extension"], seq_len)
        if windows:
            windows_by_class[cls] = windows

    if not windows_by_class:
        return []

    # Step 4 – dispatch (window-level parallelism)
    dispatcher = SMIDDispatcher(class_to_detector_map)

    if max_workers == 1:
        # Serial execution (useful for debugging / small sequences)
        return dispatcher.run(chromosome_seq, windows_by_class, sequence_name)

    # Flatten windows into tasks: (class, (start, end))
    tasks: List[Tuple[str, Tuple[int, int]]] = [
        (cls, win) for cls, wins in windows_by_class.items() for win in wins
    ]

    # Parallel execution — each window is an independent unit of work
    results: List[Dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_map = {
            executor.submit(
                _run_window_task,
                cls,
                win,
                chromosome_seq[win[0]: win[1] + 1],
                win[0],
                sequence_name,
                class_to_detector_map,
            ): (cls, win)
            for cls, win in tasks
        }
        for future in as_completed(future_map):
            results.extend(future.result())

    return results


# ---------------------------------------------------------------------------
# Worker function (must be top-level for multiprocessing pickling)
# ---------------------------------------------------------------------------

def _run_window_task(
    class_name: str,
    window: Tuple[int, int],
    window_seq: str,
    win_start: int,
    sequence_name: str,
    class_to_detector_map: Dict[str, Any] | None,
) -> List[Dict[str, Any]]:
    """Run one detector on one window — called by worker processes."""
    dispatcher = SMIDDispatcher(class_to_detector_map)
    # Pass the window as a 0-based local window to avoid allocating a padded
    # chromosome string.  The dispatcher translates coordinates using win_start.
    return dispatcher.run_window(
        window_seq,
        class_name,
        win_start,
        sequence_name,
    )


__all__ = [
    "smid_pipeline",
    "SMIDSeedEngine",
    "SMIDDispatcher",
    "linear_1d_clustering",
    "extend_and_merge",
    "SMID_SEEDS",
]
