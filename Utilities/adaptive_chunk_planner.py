"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Adaptive Chunk Planner - Dynamic Chunk Size / Worker Count Selection         │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Decides chunk size, overlap, worker count, and execution mode dynamically
    based on genome length, available RAM, and CPU core count.

    The overlap is fixed at 2 000 bp to ensure motifs at chunk boundaries
    (longest known Non-B motif is ~2 000 bp) are always detected fully.

    Execution modes:

        "disk_stream"  – Sequential chunk processing, one chunk at a time.
                         Used when RAM is scarce (< 4 GB) or small genomes
                         do not justify parallelism overhead.

        "hybrid"       – Parallel chunk processing with up to 4 workers.
                         RAM budget must be ≥ 4 GB.

USAGE::

    from Utilities.system_resource_inspector import SystemResourceInspector

    resources = SystemResourceInspector()
    planner   = AdaptiveChunkPlanner()

    plan = planner.plan(
        genome_length = 50_000_000,
        ram_budget    = resources.get_memory_budget(),
        cpu_count     = resources.get_cpu_count(),
    )
    # plan["chunk_size"], plan["overlap"], plan["workers"], plan["mode"]
"""

from __future__ import annotations

import logging
from typing import Any, Dict

logger = logging.getLogger(__name__)


class AdaptiveChunkPlanner:
    """
    Determine optimal chunk size and worker count for a given environment.

    The overlap is always fixed at :attr:`OVERLAP` (2 000 bp) regardless of
    chunk size so that the :class:`ChunkGenerator` covers motifs that span
    chunk boundaries.

    Usage::

        planner = AdaptiveChunkPlanner()
        plan = planner.plan(50_000_000, ram_budget=2_000_000_000, cpu_count=4)
        # {'chunk_size': 25000, 'overlap': 2000, 'workers': 2, 'mode': 'disk_stream'}
    """

    #: Fixed overlap between consecutive chunks (bp).
    #: Must be ≥ the longest Non-B DNA motif (~2 000 bp).
    OVERLAP: int = 2_000

    def plan(
        self,
        genome_length: int,
        ram_budget: int,
        cpu_count: int,
    ) -> Dict[str, Any]:
        """
        Compute chunk parameters appropriate for the given resources.

        Args:
            genome_length: Total sequence length in bp.
            ram_budget:    Safe available RAM in bytes
                           (e.g. from ``SystemResourceInspector.get_memory_budget()``).
            cpu_count:     Number of logical CPU cores.

        Returns:
            Dict with keys:

            * ``chunk_size`` – int, chunk length in bp
            * ``overlap``    – int, always :attr:`OVERLAP` (2 000 bp)
            * ``workers``    – int, suggested parallel worker count
            * ``mode``       – str, ``"disk_stream"`` or ``"hybrid"``
        """
        if ram_budget < 1_000_000_000:          # < 1 GB
            chunk_size = 10_000
            workers = 1
            mode = "disk_stream"

        elif ram_budget < 4_000_000_000:        # 1 – 4 GB
            chunk_size = 25_000
            workers = min(2, cpu_count)
            mode = "disk_stream"

        else:                                   # ≥ 4 GB
            chunk_size = 50_000
            workers = min(4, cpu_count)
            mode = "hybrid"

        plan: Dict[str, Any] = {
            "chunk_size": chunk_size,
            "overlap": self.OVERLAP,
            "workers": workers,
            "mode": mode,
        }

        logger.info(
            f"AdaptiveChunkPlanner: genome={genome_length:,} bp, "
            f"ram_budget={ram_budget / 1e9:.2f} GB, cpus={cpu_count} → "
            f"chunk_size={chunk_size:,}, workers={workers}, mode={mode}"
        )
        return plan
