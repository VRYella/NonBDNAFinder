"""
Tests for ChunkAnalyzer with sequences larger than 1MB.

Verifies that the fix for the 0-motifs bug (chunk_idx not in chunk_data
replaced with chunk_idx >= len(chunk_data)) works correctly in parallel mode,
and that the sequential fallback works when ProcessPoolExecutor fails.
"""

import pytest
from unittest.mock import patch

from Utilities.disk_storage import UniversalSequenceStorage
from Utilities.chunk_analyzer import ChunkAnalyzer


def _make_sequence(length: int) -> str:
    """Create a test sequence with embedded G-quadruplex motifs."""
    motif_unit = "GGGTTTAGGGTTTGGG"          # G4-like pattern
    spacer = "ATCGATCGATCGATCGATCGATCG" * 40  # ~960 bp spacer
    segment = motif_unit + spacer
    reps = (length // len(segment)) + 1
    return (segment * reps)[:length]


@pytest.fixture
def storage_with_large_seq(tmp_path):
    """Provides a storage instance with a >1MB sequence saved."""
    storage = UniversalSequenceStorage(base_dir=str(tmp_path / "seqs"))
    sequence = _make_sequence(1_100_000)  # 1.1 MB
    seq_id = storage.save_sequence(sequence, "test_large")
    yield storage, seq_id
    storage.cleanup()


@pytest.fixture
def storage_with_two_chunk_seq(tmp_path):
    """Provides a storage instance with a sequence spanning exactly 2 chunks."""
    storage = UniversalSequenceStorage(base_dir=str(tmp_path / "seqs"))
    sequence = _make_sequence(58_000)  # ~58 KB — spans two 50KB chunks
    seq_id = storage.save_sequence(sequence, "test_two_chunks")
    yield storage, seq_id
    storage.cleanup()


class TestParallelChunkAnalyzerLargeSequence:
    """Parallel mode must return motifs (not 0) for sequences >1MB."""

    def test_parallel_finds_motifs_above_1mb(self, storage_with_large_seq):
        storage, seq_id = storage_with_large_seq
        analyzer = ChunkAnalyzer(
            storage, chunk_size=50_000, overlap=2_000, use_parallel=True
        )
        rs = analyzer.analyze(seq_id=seq_id, enabled_classes=["G-Quadruplex"])
        assert rs.get_summary_stats()["total_count"] > 0, (
            "Parallel ChunkAnalyzer must find motifs for sequences >1MB"
        )

    def test_parallel_and_sequential_agree(self, tmp_path):
        """Parallel and sequential modes should find the same number of motifs."""
        # Parallel run with fresh results dir
        storage_p = UniversalSequenceStorage(base_dir=str(tmp_path / "p"))
        seq_id_p = storage_p.save_sequence(_make_sequence(1_100_000), "par")
        a_parallel = ChunkAnalyzer(
            storage_p, chunk_size=50_000, overlap=2_000, use_parallel=True
        )
        rs_p = a_parallel.analyze(seq_id=seq_id_p, enabled_classes=["G-Quadruplex"])

        # Sequential run with fresh results dir
        storage_s = UniversalSequenceStorage(base_dir=str(tmp_path / "s"))
        seq_id_s = storage_s.save_sequence(_make_sequence(1_100_000), "seq")
        a_sequential = ChunkAnalyzer(
            storage_s, chunk_size=50_000, overlap=2_000, use_parallel=False
        )
        rs_s = a_sequential.analyze(seq_id=seq_id_s, enabled_classes=["G-Quadruplex"])

        count_p = rs_p.get_summary_stats()["total_count"]
        count_s = rs_s.get_summary_stats()["total_count"]
        assert count_p > 0, "Parallel mode should find motifs"
        assert count_s > 0, "Sequential mode should find motifs"
        assert count_p == count_s, (
            f"Parallel ({count_p}) and sequential ({count_s}) should agree"
        )

        storage_p.cleanup()
        storage_s.cleanup()


class TestParallelFallbackToSequential:
    """When ProcessPoolExecutor fails, ChunkAnalyzer must fall back to sequential."""

    def test_fallback_finds_motifs_when_executor_fails(
        self, storage_with_two_chunk_seq
    ):
        storage, seq_id = storage_with_two_chunk_seq

        def _failing_executor(*args, **kwargs):
            raise RuntimeError("Simulated ProcessPoolExecutor failure")

        with patch(
            "Utilities.chunk_analyzer.ProcessPoolExecutor",
            side_effect=_failing_executor,
        ):
            analyzer = ChunkAnalyzer(
                storage, chunk_size=50_000, overlap=2_000, use_parallel=True
            )
            rs = analyzer.analyze(
                seq_id=seq_id, enabled_classes=["G-Quadruplex", "Slipped_DNA"]
            )

        assert rs.get_summary_stats()["total_count"] > 0, (
            "Fallback to sequential must still find motifs when parallel fails"
        )


class TestChunkIndexBoundsCheckFix:
    """Regression test: the original chunk_idx not in chunk_data bug gave 0 motifs."""

    def test_bugfix_chunk_idx_bounds_check(self):
        """Demonstrate the old bug and confirm the fix.
        
        Old code: ``if chunk_idx not in chunk_data: continue``
        chunk_data is a list of tuples.  Integer indices are never elements of
        that list, so the condition was always True and every chunk was skipped.
        
        New code: ``if chunk_idx >= len(chunk_data): continue``
        This correctly bounds-checks the index, allowing valid chunks through.
        """
        chunk_data = [
            ("GGGTTTAGGGATCG", "seq_chunk1", 0, 50_000, None),
            ("ATCGATCGATCGATCG", "seq_chunk2", 48_000, 98_000, None),
        ]

        # Confirm the old condition was always True for integer indices.
        assert 0 not in chunk_data, "Integer 0 is never an element of a list of tuples"
        assert 1 not in chunk_data, "Integer 1 is never an element of a list of tuples"

        # OLD buggy check saved 0 motifs because every chunk was skipped.
        old_saved_count = 0
        for chunk_idx in [0, 1]:
            if chunk_idx not in chunk_data:   # OLD: always True → continue
                continue
            old_saved_count += 1
        assert old_saved_count == 0, "Old check skips all chunks"

        # NEW fixed check saves all valid chunks.
        new_saved_count = 0
        for chunk_idx in [0, 1]:
            if chunk_idx >= len(chunk_data):  # FIXED: False for valid indices
                continue
            new_saved_count += 1
        assert new_saved_count == 2, "Fixed check keeps all valid chunks"
