"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ PerformanceMonitor - Real-Time Performance Telemetry                         │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.2            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Thread-safe performance monitor that tracks per-detector, per-chunk, and
    per-stage timings across an entire analysis run.

    Tracks:
        - Per-detector runtime (across all chunks)
        - Per-chunk runtime + motif count
        - Named stage durations (detection, filtering, visualization)
        - Peak RSS memory via psutil (optional)
        - CPU utilisation snapshot

    Usage::

        from Utilities.core.performance_monitor import PerformanceMonitor

        monitor = PerformanceMonitor()
        monitor.start()

        # Inside a chunk worker (called from the main process after worker returns)
        monitor.record_chunk(chunk_id=0, elapsed=1.2, motif_count=45, bp_count=50_000)
        monitor.record_detector(chunk_id=0, detector_name="G-Quadruplex",
                                elapsed=0.3, bp_count=50_000)

        # After detection stage
        monitor.record_stage("detection", elapsed=5.6)

        summary = monitor.get_summary()
        print(summary["throughput_bps"])
        print(summary["detector_breakdown"])
"""

from __future__ import annotations

import logging
import os
import threading
from time import perf_counter
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


class PerformanceMonitor:
    """
    Thread-safe real-time performance monitor for motif detection runs.

    All ``record_*`` methods are safe to call from multiple threads
    (e.g., when collecting results from a ``ProcessPoolExecutor``).

    Attributes are intentionally private to enforce the ``record_*`` /
    ``get_summary()`` API boundary.
    """

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._start_time: Optional[float] = None
        self._stage_records: Dict[str, float] = {}
        self._chunk_records: List[Dict[str, Any]] = []
        self._detector_records: List[Dict[str, Any]] = []
        self._peak_memory_mb: float = 0.0

    # ------------------------------------------------------------------
    # LIFECYCLE
    # ------------------------------------------------------------------

    def start(self) -> None:
        """Start the global timer.  Must be called before ``record_*`` methods."""
        self._start_time = perf_counter()
        self.snapshot_memory()

    # ------------------------------------------------------------------
    # RECORDING
    # ------------------------------------------------------------------

    def record_stage(self, stage_name: str, elapsed: float) -> None:
        """
        Record the wall-clock duration of a named analysis stage.

        Args:
            stage_name: Human-readable stage name (e.g. ``"detection"``,
                        ``"visualization"``).
            elapsed:    Duration in seconds.
        """
        with self._lock:
            self._stage_records[stage_name] = elapsed

    def record_chunk(
        self,
        chunk_id: int,
        elapsed: float,
        motif_count: int,
        bp_count: int,
    ) -> None:
        """
        Record the processing of a single genome chunk.

        Args:
            chunk_id:    Zero-based chunk index.
            elapsed:     Wall-clock time to process the chunk (seconds).
            motif_count: Number of motifs detected in this chunk.
            bp_count:    Number of base-pairs in this chunk.
        """
        with self._lock:
            self._chunk_records.append(
                {
                    "chunk_id": chunk_id,
                    "elapsed": elapsed,
                    "motif_count": motif_count,
                    "bp_count": bp_count,
                }
            )

    def record_detector(
        self,
        chunk_id: int,
        detector_name: str,
        elapsed: float,
        bp_count: int = 0,
    ) -> None:
        """
        Record one detector's runtime within a chunk.

        Args:
            chunk_id:      Zero-based chunk index.
            detector_name: Canonical detector name (e.g. ``"g_quadruplex"``).
            elapsed:       Detector wall-clock time (seconds).
            bp_count:      Base-pairs processed (same as chunk size; informational).
        """
        with self._lock:
            self._detector_records.append(
                {
                    "chunk_id": chunk_id,
                    "detector": detector_name,
                    "elapsed": elapsed,
                    "bp_count": bp_count,
                }
            )

    # ------------------------------------------------------------------
    # MEMORY / CPU
    # ------------------------------------------------------------------

    def snapshot_memory(self) -> float:
        """
        Capture current RSS memory and update peak.

        Returns:
            Current RSS memory in MB, or 0.0 if psutil is unavailable.
        """
        try:
            import psutil  # noqa: PLC0415

            proc = psutil.Process(os.getpid())
            mb = proc.memory_info().rss / 1_048_576  # bytes → MB
            with self._lock:
                if mb > self._peak_memory_mb:
                    self._peak_memory_mb = mb
            return mb
        except Exception:  # psutil not installed or permissions error
            return 0.0

    def get_cpu_percent(self) -> float:
        """
        Return current process CPU utilisation (0–100).

        Returns 0.0 if psutil is unavailable.
        """
        try:
            import psutil  # noqa: PLC0415

            return psutil.Process(os.getpid()).cpu_percent(interval=0.1)
        except Exception:
            return 0.0

    # ------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------

    def get_summary(self) -> Dict[str, Any]:
        """
        Return a snapshot of all collected performance metrics.

        Returns:
            Dict with keys::

                {
                    "total_elapsed":      float,   # seconds since start()
                    "total_bp_processed": int,
                    "throughput_bps":     float,   # bp / second
                    "peak_memory_mb":     float,
                    "chunk_count":        int,
                    "chunk_records":      list,    # per-chunk dicts
                    "stage_times":        dict,    # stage_name → elapsed
                    "detector_breakdown": dict,    # detector → {total, avg, calls}
                    "slowest_detector":   str | None,
                }
        """
        with self._lock:
            elapsed_since_start = (
                perf_counter() - self._start_time if self._start_time else 0.0
            )
            total_bp = sum(c["bp_count"] for c in self._chunk_records)

            # --- Per-detector aggregation ---------------------------------
            detector_agg: Dict[str, Dict[str, Any]] = {}
            for rec in self._detector_records:
                name = rec["detector"]
                if name not in detector_agg:
                    detector_agg[name] = {"total_elapsed": 0.0, "call_count": 0}
                detector_agg[name]["total_elapsed"] += rec["elapsed"]
                detector_agg[name]["call_count"] += 1

            breakdown = {
                name: {
                    "total_elapsed": data["total_elapsed"],
                    "avg_elapsed": (
                        data["total_elapsed"] / data["call_count"]
                        if data["call_count"] > 0
                        else 0.0
                    ),
                    "call_count": data["call_count"],
                    "pct_total": (
                        data["total_elapsed"] / elapsed_since_start * 100.0
                        if elapsed_since_start > 0
                        else 0.0
                    ),
                }
                for name, data in detector_agg.items()
            }

            slowest = (
                max(breakdown, key=lambda n: breakdown[n]["total_elapsed"])
                if breakdown
                else None
            )

            # --- Chunk distribution stats ---------------------------------
            chunk_times = [c["elapsed"] for c in self._chunk_records]
            avg_chunk_time = (
                sum(chunk_times) / len(chunk_times) if chunk_times else 0.0
            )

            return {
                "total_elapsed": elapsed_since_start,
                "total_bp_processed": total_bp,
                "throughput_bps": (
                    total_bp / elapsed_since_start if elapsed_since_start > 0 else 0.0
                ),
                "peak_memory_mb": self._peak_memory_mb,
                "chunk_count": len(self._chunk_records),
                "avg_chunk_time": avg_chunk_time,
                "chunk_records": list(self._chunk_records),
                "stage_times": dict(self._stage_records),
                "detector_breakdown": breakdown,
                "slowest_detector": slowest,
            }

    def format_summary(self) -> str:
        """
        Return a human-readable performance summary table.

        Returns:
            Multi-line string with detection breakdown and resource metrics.
        """
        s = self.get_summary()

        lines: List[str] = [
            "══════════════════════════════════════════════════",
            "  Performance Summary",
            "══════════════════════════════════════════════════",
            f"  Total runtime      : {s['total_elapsed']:.3f} s",
            f"  Base pairs         : {s['total_bp_processed']:,}",
            f"  Throughput         : {s['throughput_bps']:,.0f} bp/s",
            f"  Peak memory        : {s['peak_memory_mb']:.1f} MB",
            f"  Chunks processed   : {s['chunk_count']}",
            f"  Avg chunk time     : {s['avg_chunk_time']:.3f} s",
        ]

        if s.get("stage_times"):
            lines.append("")
            lines.append("  Stage Times:")
            for stage, t in sorted(s["stage_times"].items()):
                lines.append(f"    {stage:<20} {t:.3f} s")

        if s.get("detector_breakdown"):
            lines.append("")
            lines.append(
                "  Detection Breakdown (cumulative across all chunks):"
            )
            lines.append(
                f"    {'Detector':<22} {'Time (s)':>9}  {'% Total':>8}  {'Calls':>6}"
            )
            lines.append("    " + "-" * 50)
            for name, d in sorted(
                s["detector_breakdown"].items(),
                key=lambda kv: -kv[1]["total_elapsed"],
            ):
                lines.append(
                    f"    {name:<22} {d['total_elapsed']:>9.3f}  "
                    f"{d['pct_total']:>7.1f}%  {d['call_count']:>6}"
                )
            if s.get("slowest_detector"):
                lines.append(f"\n  Slowest detector: {s['slowest_detector']}")

        lines.append("══════════════════════════════════════════════════")
        return "\n".join(lines)
