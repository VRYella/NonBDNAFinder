"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ SequenceContext - Shared Preprocessing Context for Detectors                 │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.2            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Encapsulates a DNA sequence with all one-time preprocessing already applied.

    Instead of each detector independently calling ``sequence.upper()``,
    a single ``SequenceContext`` is created per chunk and shared across all
    detectors that run on that chunk.

    Benefits:
        - Uppercase conversion performed exactly once per chunk
        - Sequence length cached (avoids repeated ``len()`` calls)
        - Optional prefix-GC array precomputed once for GC-heavy detectors
        - Clean interface: ``detector.scan(context)`` instead of
          ``detector.detect_motifs(raw_sequence, name)``

USAGE::

    from Utilities.core.sequence_context import SequenceContext

    ctx = SequenceContext("atcgGGGTTTaggg", name="chr1_chunk0")
    print(ctx.sequence)  # "ATCGGGGTTTAGGG"
    print(ctx.length)    # 14
"""

from __future__ import annotations

from typing import Optional
import numpy as np


class SequenceContext:
    """
    Preprocessed DNA sequence context shared across all detectors for one chunk.

    Parameters
    ----------
    sequence : str
        Raw DNA sequence string (any case).
    name : str
        Identifier for the sequence / chunk.
    precompute_gc : bool
        When ``True``, build a cumulative GC prefix array so detectors can
        query GC content for any sub-interval in O(1).  Costs O(n) time and
        memory; disabled by default.

    Attributes
    ----------
    sequence : str
        Uppercase, stripped DNA sequence ready for regex / scoring.
    name : str
        Sequence identifier passed through from the caller.
    length : int
        Length of ``sequence`` (cached – avoids repeated ``len()`` calls).
    gc_prefix : np.ndarray or None
        Cumulative GC count array of shape ``(length + 1,)`` where
        ``gc_prefix[i+1] - gc_prefix[j]`` gives the GC count in
        ``sequence[j:i+1]``.  ``None`` when ``precompute_gc=False``.
    """

    __slots__ = ("sequence", "name", "length", "gc_prefix")

    def __init__(
        self,
        sequence: str,
        name: str = "sequence",
        precompute_gc: bool = False,
    ) -> None:
        self.sequence: str = sequence.upper().strip()
        self.name: str = name
        self.length: int = len(self.sequence)
        self.gc_prefix: Optional[np.ndarray] = None

        if precompute_gc and self.length > 0:
            gc_flag = np.frombuffer(
                self.sequence.encode("ascii"), dtype=np.uint8
            )
            # G = 71, C = 67
            is_gc = (gc_flag == 71) | (gc_flag == 67)
            prefix = np.empty(self.length + 1, dtype=np.int64)
            prefix[0] = 0
            np.cumsum(is_gc, out=prefix[1:])
            self.gc_prefix = prefix

    # ------------------------------------------------------------------
    # GC helpers
    # ------------------------------------------------------------------

    def gc_content(self, start: int = 0, end: Optional[int] = None) -> float:
        """
        Fraction of G/C bases in ``sequence[start:end]``.

        Uses the precomputed prefix array when available (O(1));
        falls back to a direct count otherwise (O(n)).

        Args:
            start: Inclusive start index (default 0).
            end:   Exclusive end index (default ``length``).

        Returns:
            GC fraction in ``[0.0, 1.0]``.
        """
        if end is None:
            end = self.length
        region_len = end - start
        if region_len <= 0:
            return 0.0

        if self.gc_prefix is not None:
            gc_count = int(self.gc_prefix[end] - self.gc_prefix[start])
        else:
            sub = self.sequence[start:end]
            gc_count = sub.count("G") + sub.count("C")

        return gc_count / region_len

    # ------------------------------------------------------------------
    # Dunder helpers
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        return (
            f"SequenceContext(name={self.name!r}, length={self.length:,}, "
            f"gc_precomputed={self.gc_prefix is not None})"
        )

    def __len__(self) -> int:
        return self.length
