"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ System Resource Inspector - RAM / CPU / Disk Availability Detection          │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Detects RAM, disk, and CPU availability and determines safe operating limits.
    Used by AdaptiveChunkPlanner to size chunks and worker pools appropriately.

USAGE::

    inspector = SystemResourceInspector()
    budget = inspector.get_memory_budget()   # 60% of available RAM in bytes
    cpus   = inspector.get_cpu_count()
    disk   = inspector.get_available_disk("/tmp")
"""

from __future__ import annotations

import logging
import os
from typing import Optional

logger = logging.getLogger(__name__)

# Fraction of available RAM exposed as the safe memory budget
_MEMORY_BUDGET_FRACTION: float = 0.60


class SystemResourceInspector:
    """
    Inspect system resources and compute safe operating limits.

    All values are in bytes unless otherwise stated.

    Usage::

        inspector = SystemResourceInspector()
        print(f"RAM budget : {inspector.get_memory_budget() / 1e9:.2f} GB")
        print(f"CPU cores  : {inspector.get_cpu_count()}")
        print(f"Disk free  : {inspector.get_available_disk('/tmp') / 1e9:.2f} GB")
    """

    def get_available_ram(self) -> int:
        """
        Return currently available (free + reclaimable) RAM in bytes.

        Uses ``psutil`` if available; falls back to a conservative 512 MB
        estimate so the rest of the pipeline can still run in restricted
        environments.
        """
        try:
            import psutil
            return psutil.virtual_memory().available
        except Exception as exc:
            logger.warning(f"psutil unavailable – defaulting RAM to 512 MB: {exc}")
            return 512 * 1024 * 1024  # 512 MB fallback

    def get_total_ram(self) -> int:
        """
        Return total installed RAM in bytes.

        Falls back to 1 GB when ``psutil`` is unavailable.
        """
        try:
            import psutil
            return psutil.virtual_memory().total
        except Exception as exc:
            logger.warning(f"psutil unavailable – defaulting total RAM to 1 GB: {exc}")
            return 1 * 1024 * 1024 * 1024

    def get_cpu_count(self) -> int:
        """
        Return the number of logical CPU cores available.

        Uses ``os.cpu_count()`` as the primary source; falls back to 1.
        """
        count = os.cpu_count()
        if count is None or count < 1:
            logger.warning("os.cpu_count() returned None – defaulting to 1")
            return 1
        return count

    def get_available_disk(self, path: str = "/tmp") -> int:
        """
        Return free disk space at *path* in bytes.

        Args:
            path: Filesystem path to inspect (default ``/tmp``).

        Falls back to 1 GB when ``psutil`` / ``shutil`` is unavailable.
        """
        try:
            import shutil
            return shutil.disk_usage(path).free
        except Exception as exc:
            logger.warning(
                f"Could not determine disk space at '{path}' – defaulting to 1 GB: {exc}"
            )
            return 1 * 1024 * 1024 * 1024

    def get_memory_budget(self) -> int:
        """
        Return the safe usable RAM budget in bytes.

        Defined as ``_MEMORY_BUDGET_FRACTION`` (60 %) of currently available
        RAM so that the analysis pipeline does not starve other processes.

        Returns:
            Safe memory budget in bytes (≥ 1 byte).
        """
        budget = int(self.get_available_ram() * _MEMORY_BUDGET_FRACTION)
        return max(budget, 1)
