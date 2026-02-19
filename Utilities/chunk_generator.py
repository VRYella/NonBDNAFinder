"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Chunk Generator - Overlapping Genome Segment Iterator                        │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Yields overlapping genome segments from a raw sequence string.

    The overlap ensures that motifs spanning chunk boundaries are detected
    fully in at least one chunk.  The ``core_end`` field marks the position
    up to which a chunk's motifs are authoritative; motifs starting at or
    beyond ``core_end`` are duplicates of the next chunk and should be
    discarded by the OverlapDeduplicator.

    Example with chunk_size=50 000 and overlap=2 000:

        Chunk 1:  [    0 –  50 000]  core_end = 48 000
        Chunk 2:  [48 000 –  98 000]  core_end = 96 000
        Chunk 3:  [96 000 – 146 000]  core_end = 144 000
        ...
        Last chunk: core_end = end (full chunk is authoritative)

USAGE::

    gen = ChunkGenerator(sequence="ACGT...", chunk_size=50_000, overlap=2_000)
    for chunk in gen.generate():
        # chunk["sequence"], chunk["start"], chunk["end"], chunk["core_end"]
        run_detector(chunk)
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Generator

logger = logging.getLogger(__name__)


class ChunkGenerator:
    """
    Yield overlapping genome chunks as dicts for downstream processing.

    Each yielded dict contains:

    * ``sequence``  – str, the raw DNA subsequence for this chunk
    * ``start``     – int, genome-global start position (0-based, inclusive)
    * ``end``       – int, genome-global end position (exclusive)
    * ``core_end``  – int, the exclusive boundary of the authoritative region;
                     motifs starting at or after this position are overlap
                     duplicates and should be discarded (except for the last
                     chunk, where ``core_end == end``).

    Usage::

        gen = ChunkGenerator(sequence, chunk_size=50_000, overlap=2_000)
        for chunk in gen.generate():
            process(chunk["sequence"], chunk["start"], chunk["core_end"])
    """

    def __init__(self, genome_sequence: str, chunk_size: int, overlap: int):
        """
        Args:
            genome_sequence: Full DNA sequence string (uppercase, no whitespace).
            chunk_size:      Target chunk length in bp.
            overlap:         Overlap between consecutive chunks in bp.
                             Must be < chunk_size.

        Raises:
            ValueError: If overlap ≥ chunk_size.
        """
        if overlap >= chunk_size:
            raise ValueError(
                f"overlap ({overlap}) must be less than chunk_size ({chunk_size})"
            )
        self.seq = genome_sequence
        self.chunk_size = chunk_size
        self.overlap = overlap

    def generate(self) -> Generator[Dict[str, Any], None, None]:
        """
        Yield chunk dicts in genome order.

        Yields:
            Dict with keys ``sequence``, ``start``, ``end``, ``core_end``.
        """
        genome_length = len(self.seq)
        start = 0
        chunk_num = 0

        while start < genome_length:
            end = min(start + self.chunk_size, genome_length)
            is_last = end >= genome_length

            # core_end: motifs from [start, core_end) are authoritative.
            # For the last chunk the entire chunk is authoritative.
            core_end = end if is_last else end - self.overlap

            chunk_num += 1
            logger.debug(
                f"ChunkGenerator: chunk {chunk_num} "
                f"[{start:,}–{end:,}] core_end={core_end:,}"
            )

            yield {
                "sequence": self.seq[start:end],
                "start": start,
                "end": end,
                "core_end": core_end,
            }

            if is_last:
                break

            # Advance by (chunk_size – overlap) so next chunk overlaps by `overlap` bp
            start = end - self.overlap

        logger.info(
            f"ChunkGenerator: yielded {chunk_num} chunk(s) "
            f"for sequence length {genome_length:,}"
        )
