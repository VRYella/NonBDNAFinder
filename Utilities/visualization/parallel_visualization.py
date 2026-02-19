"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ ParallelVisualization - Concurrent Plot Generation                           │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.2            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Parallelises the five standard ``VisualizationPipeline`` plots using a
    ``ThreadPoolExecutor`` (not ``ProcessPoolExecutor``) because:

        * Matplotlib is not fork-safe
        * Data is already in memory (no IPC overhead)
        * Plot rendering is CPU-light relative to motif detection

    Each plot is generated concurrently and its wall-clock time is recorded.
    The aggregated timing is included in the return value so it can be
    surfaced in the GUI performance dashboard.

USAGE::

    from Utilities.visualization_pipeline import VisualizationPipeline
    from Utilities.visualization.parallel_visualization import ParallelVisualization

    pipeline = VisualizationPipeline(accumulator.get_summary())
    runner  = ParallelVisualization(pipeline, max_workers=4)
    result  = runner.generate_all()

    # result["figures"]  → {plot_name: Figure}
    # result["timings"]  → {plot_name: elapsed_seconds}
    # result["total_elapsed"] → float
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from time import perf_counter
from typing import Any, Callable, Dict, Optional, Tuple

logger = logging.getLogger(__name__)

# Default number of concurrent plot threads.
# Four plots are generated in the standard pipeline; using 4 workers
# allows all of them to render concurrently.
_DEFAULT_VIZ_WORKERS = 4


class ParallelVisualization:
    """
    Concurrent wrapper around ``VisualizationPipeline.generate_all()``.

    All five standard plots are submitted to a ``ThreadPoolExecutor``
    simultaneously so the total visualization time approaches that of the
    single slowest plot rather than the sum of all plots.

    Per-plot timing is recorded and returned alongside the figures.

    Parameters
    ----------
    pipeline : VisualizationPipeline
        Fully-initialised pipeline instance backed by pre-binned summary data.
    max_workers : int
        Number of concurrent threads.  Defaults to ``min(4, plot_count)``.
    progress_callback : callable, optional
        Called after each plot finishes with a dict::

            {
                "stage":     "visualization",
                "component": plot_name,
                "elapsed":   float,          # seconds for this plot
            }
    """

    def __init__(
        self,
        pipeline,
        max_workers: int = _DEFAULT_VIZ_WORKERS,
        progress_callback: Optional[Callable[[Dict[str, Any]], None]] = None,
    ) -> None:
        self._pipeline = pipeline
        self._max_workers = max_workers
        self._progress_callback = progress_callback

    # ------------------------------------------------------------------
    # PUBLIC API
    # ------------------------------------------------------------------

    def generate_all(self) -> Dict[str, Any]:
        """
        Render all standard plots in parallel and return results + timings.

        Returns:
            Dict with keys::

                {
                    "figures": {
                        "density_histogram":     Figure,
                        "length_histogram":      Figure,
                        "class_distribution":    Figure,
                        "subclass_distribution": Figure,
                        "cooccurrence_heatmap":  Figure,
                    },
                    "timings": {plot_name: elapsed_seconds, ...},
                    "total_elapsed": float,
                }
        """
        plot_tasks: Dict[str, Callable] = {
            "density_histogram": self._pipeline.plot_density_histogram,
            "length_histogram": self._pipeline.plot_length_histogram,
            "class_distribution": self._pipeline.plot_class_distribution,
            "subclass_distribution": self._pipeline.plot_subclass_distribution,
            "cooccurrence_heatmap": self._pipeline.plot_cooccurrence_heatmap,
        }

        figures: Dict[str, Any] = {}
        timings: Dict[str, float] = {}
        t_start = perf_counter()

        n_workers = min(self._max_workers, len(plot_tasks))

        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            future_to_name = {
                executor.submit(self._timed_plot, name, fn): name
                for name, fn in plot_tasks.items()
            }
            for future in as_completed(future_to_name):
                plot_name = future_to_name[future]
                try:
                    fig, elapsed = future.result()
                    figures[plot_name] = fig
                    timings[plot_name] = elapsed

                    if self._progress_callback is not None:
                        try:
                            self._progress_callback(
                                {
                                    "stage": "visualization",
                                    "component": plot_name,
                                    "elapsed": elapsed,
                                }
                            )
                        except Exception as cb_exc:
                            logger.debug(
                                f"ParallelVisualization: callback error: {cb_exc}"
                            )

                    logger.debug(
                        f"ParallelVisualization: '{plot_name}' "
                        f"rendered in {elapsed:.3f}s"
                    )
                except Exception as exc:
                    logger.warning(
                        f"ParallelVisualization: '{plot_name}' failed: {exc}"
                    )
                    figures[plot_name] = None
                    timings[plot_name] = 0.0

        total_elapsed = perf_counter() - t_start
        logger.info(
            f"ParallelVisualization: all plots done in {total_elapsed:.3f}s "
            f"(sequential would be ~{sum(timings.values()):.3f}s)"
        )

        return {
            "figures": figures,
            "timings": timings,
            "total_elapsed": total_elapsed,
        }

    # ------------------------------------------------------------------
    # HELPERS
    # ------------------------------------------------------------------

    @staticmethod
    def _timed_plot(name: str, fn: Callable) -> Tuple[Any, float]:
        """Call ``fn()`` and return ``(result, elapsed_seconds)``."""
        t0 = perf_counter()
        result = fn()
        return result, perf_counter() - t0
