"""
Tests for the new adaptive chunking architecture components:
  - SystemResourceInspector
  - AdaptiveChunkPlanner
  - ChunkGenerator
  - OverlapDeduplicator
  - FinalExporter
  - DetectorRunner (metadata / interface only; no nonbscanner call)
  - ParallelChunkExecutor (strategy only; no nonbscanner call)
"""

import pytest
import os
from pathlib import Path


# ──────────────────────────────────────────────────────────────────────────────
# SystemResourceInspector
# ──────────────────────────────────────────────────────────────────────────────

class TestSystemResourceInspector:
    """Unit tests for SystemResourceInspector."""

    def test_available_ram_positive(self):
        from Utilities.system_resource_inspector import SystemResourceInspector
        inspector = SystemResourceInspector()
        assert inspector.get_available_ram() > 0

    def test_total_ram_positive(self):
        from Utilities.system_resource_inspector import SystemResourceInspector
        inspector = SystemResourceInspector()
        assert inspector.get_total_ram() > 0

    def test_cpu_count_at_least_one(self):
        from Utilities.system_resource_inspector import SystemResourceInspector
        inspector = SystemResourceInspector()
        assert inspector.get_cpu_count() >= 1

    def test_available_disk_positive(self, tmp_path):
        from Utilities.system_resource_inspector import SystemResourceInspector
        inspector = SystemResourceInspector()
        assert inspector.get_available_disk(str(tmp_path)) > 0

    def test_memory_budget_less_than_available_ram(self):
        from Utilities.system_resource_inspector import SystemResourceInspector
        inspector = SystemResourceInspector()
        assert inspector.get_memory_budget() <= inspector.get_available_ram()

    def test_memory_budget_positive(self):
        from Utilities.system_resource_inspector import SystemResourceInspector
        inspector = SystemResourceInspector()
        assert inspector.get_memory_budget() > 0

    def test_memory_budget_is_60pct_of_available(self):
        from Utilities.system_resource_inspector import SystemResourceInspector
        inspector = SystemResourceInspector()
        avail = inspector.get_available_ram()
        budget = inspector.get_memory_budget()
        # Budget should be approximately 60% of available RAM
        expected = int(avail * 0.60)
        assert budget == expected


# ──────────────────────────────────────────────────────────────────────────────
# AdaptiveChunkPlanner
# ──────────────────────────────────────────────────────────────────────────────

class TestAdaptiveChunkPlanner:
    """Unit tests for AdaptiveChunkPlanner."""

    def test_overlap_is_always_2000(self):
        from Utilities.adaptive_chunk_planner import AdaptiveChunkPlanner
        planner = AdaptiveChunkPlanner()
        for ram in [500_000_000, 2_000_000_000, 8_000_000_000]:
            plan = planner.plan(10_000_000, ram, 4)
            assert plan["overlap"] == 2_000

    def test_low_ram_disk_stream_single_worker(self):
        from Utilities.adaptive_chunk_planner import AdaptiveChunkPlanner
        planner = AdaptiveChunkPlanner()
        plan = planner.plan(10_000_000, ram_budget=500_000_000, cpu_count=8)
        assert plan["chunk_size"] == 10_000
        assert plan["workers"] == 1
        assert plan["mode"] == "disk_stream"

    def test_mid_ram_disk_stream_two_workers(self):
        from Utilities.adaptive_chunk_planner import AdaptiveChunkPlanner
        planner = AdaptiveChunkPlanner()
        plan = planner.plan(10_000_000, ram_budget=2_000_000_000, cpu_count=4)
        assert plan["chunk_size"] == 25_000
        assert plan["workers"] == min(2, 4)
        assert plan["mode"] == "disk_stream"

    def test_high_ram_hybrid_mode(self):
        from Utilities.adaptive_chunk_planner import AdaptiveChunkPlanner
        planner = AdaptiveChunkPlanner()
        plan = planner.plan(100_000_000, ram_budget=8_000_000_000, cpu_count=16)
        assert plan["chunk_size"] == 50_000
        assert plan["workers"] == min(4, 16)
        assert plan["mode"] == "hybrid"

    def test_workers_capped_by_cpu_count(self):
        from Utilities.adaptive_chunk_planner import AdaptiveChunkPlanner
        planner = AdaptiveChunkPlanner()
        # Only 1 CPU available; high RAM should cap workers at 1
        plan = planner.plan(100_000_000, ram_budget=8_000_000_000, cpu_count=1)
        assert plan["workers"] == 1

    def test_plan_returns_required_keys(self):
        from Utilities.adaptive_chunk_planner import AdaptiveChunkPlanner
        planner = AdaptiveChunkPlanner()
        plan = planner.plan(1_000_000, 2_000_000_000, 4)
        assert set(plan.keys()) == {"chunk_size", "overlap", "workers", "mode"}

    def test_overlap_constant_matches_class_attribute(self):
        from Utilities.adaptive_chunk_planner import AdaptiveChunkPlanner
        assert AdaptiveChunkPlanner.OVERLAP == 2_000


# ──────────────────────────────────────────────────────────────────────────────
# ChunkGenerator
# ──────────────────────────────────────────────────────────────────────────────

class TestChunkGenerator:
    """Unit tests for ChunkGenerator."""

    def _seq(self, length: int) -> str:
        return ("ACGT" * (length // 4 + 1))[:length]

    def test_single_chunk_for_short_sequence(self):
        from Utilities.chunk_generator import ChunkGenerator
        seq = self._seq(1_000)
        gen = ChunkGenerator(seq, chunk_size=50_000, overlap=2_000)
        chunks = list(gen.generate())
        assert len(chunks) == 1
        assert chunks[0]["start"] == 0
        assert chunks[0]["end"] == 1_000
        assert chunks[0]["core_end"] == 1_000   # last chunk: core_end == end

    def test_two_chunks_for_58k_sequence(self):
        from Utilities.chunk_generator import ChunkGenerator
        seq = self._seq(58_000)
        gen = ChunkGenerator(seq, chunk_size=50_000, overlap=2_000)
        chunks = list(gen.generate())
        assert len(chunks) == 2

    def test_overlap_is_applied_between_chunks(self):
        from Utilities.chunk_generator import ChunkGenerator
        seq = self._seq(100_000)
        gen = ChunkGenerator(seq, chunk_size=50_000, overlap=2_000)
        chunks = list(gen.generate())
        # Second chunk should start at (50_000 - 2_000) = 48_000
        assert chunks[1]["start"] == 48_000

    def test_core_end_excludes_overlap(self):
        from Utilities.chunk_generator import ChunkGenerator
        seq = self._seq(100_000)
        gen = ChunkGenerator(seq, chunk_size=50_000, overlap=2_000)
        chunks = list(gen.generate())
        # First chunk core_end = 50_000 - 2_000 = 48_000
        assert chunks[0]["core_end"] == 48_000

    def test_last_chunk_core_end_equals_end(self):
        from Utilities.chunk_generator import ChunkGenerator
        seq = self._seq(100_000)
        gen = ChunkGenerator(seq, chunk_size=50_000, overlap=2_000)
        chunks = list(gen.generate())
        last = chunks[-1]
        assert last["core_end"] == last["end"]

    def test_all_bases_covered(self):
        from Utilities.chunk_generator import ChunkGenerator
        seq = self._seq(120_000)
        gen = ChunkGenerator(seq, chunk_size=50_000, overlap=2_000)
        chunks = list(gen.generate())
        assert chunks[0]["start"] == 0
        assert chunks[-1]["end"] == 120_000

    def test_sequence_content_matches(self):
        from Utilities.chunk_generator import ChunkGenerator
        seq = self._seq(60_000)
        gen = ChunkGenerator(seq, chunk_size=50_000, overlap=2_000)
        chunks = list(gen.generate())
        for chunk in chunks:
            expected = seq[chunk["start"]:chunk["end"]]
            assert chunk["sequence"] == expected

    def test_overlap_must_be_less_than_chunk_size(self):
        from Utilities.chunk_generator import ChunkGenerator
        with pytest.raises(ValueError):
            ChunkGenerator(self._seq(10_000), chunk_size=2_000, overlap=2_000)

    def test_empty_sequence_yields_no_chunks(self):
        from Utilities.chunk_generator import ChunkGenerator
        gen = ChunkGenerator("", chunk_size=50_000, overlap=2_000)
        chunks = list(gen.generate())
        assert chunks == []

    def test_exact_multiple_of_chunk_size(self):
        from Utilities.chunk_generator import ChunkGenerator
        seq = self._seq(100_000)
        gen = ChunkGenerator(seq, chunk_size=50_000, overlap=2_000)
        chunks = list(gen.generate())
        assert chunks[-1]["end"] == 100_000


# ──────────────────────────────────────────────────────────────────────────────
# OverlapDeduplicator
# ──────────────────────────────────────────────────────────────────────────────

class TestOverlapDeduplicator:
    """Unit tests for OverlapDeduplicator."""

    def _motifs(self, starts):
        return [
            {"Class": "G-Quadruplex", "Subclass": "G4",
             "Start": s, "End": s + 20, "Length": 20, "Score": 0.9}
            for s in starts
        ]

    def test_keeps_motifs_before_core_end(self):
        from Utilities.overlap_deduplicator import OverlapDeduplicator
        dedup = OverlapDeduplicator()
        motifs = self._motifs([100, 200, 47_999])
        result = dedup.filter_core(motifs, core_end=48_000)
        assert len(result) == 3

    def test_removes_motifs_at_core_end(self):
        from Utilities.overlap_deduplicator import OverlapDeduplicator
        dedup = OverlapDeduplicator()
        motifs = self._motifs([48_000, 49_000])
        result = dedup.filter_core(motifs, core_end=48_000)
        assert len(result) == 0

    def test_mixed_core_and_overlap(self):
        from Utilities.overlap_deduplicator import OverlapDeduplicator
        dedup = OverlapDeduplicator()
        motifs = self._motifs([100, 47_999, 48_000, 50_000])
        result = dedup.filter_core(motifs, core_end=48_000)
        assert len(result) == 2
        assert all(m["Start"] < 48_000 for m in result)

    def test_empty_motif_list(self):
        from Utilities.overlap_deduplicator import OverlapDeduplicator
        dedup = OverlapDeduplicator()
        assert dedup.filter_core([], core_end=48_000) == []

    def test_last_chunk_with_large_core_end(self):
        from Utilities.overlap_deduplicator import OverlapDeduplicator
        dedup = OverlapDeduplicator()
        motifs = self._motifs([100, 98_000, 99_999])
        result = dedup.filter_core(motifs, core_end=100_000)
        assert len(result) == 3  # all kept when core_end == end

    def test_filter_core_df_basic(self):
        """DataFrame variant filters correctly."""
        pytest.importorskip("pandas")
        import pandas as pd
        from Utilities.overlap_deduplicator import OverlapDeduplicator

        dedup = OverlapDeduplicator()
        df = pd.DataFrame([
            {"Start": 100, "End": 120},
            {"Start": 48_000, "End": 48_020},
        ])
        result = dedup.filter_core_df(df, core_end=48_000)
        assert len(result) == 1
        assert result.iloc[0]["Start"] == 100

    def test_filter_core_df_empty(self):
        pytest.importorskip("pandas")
        import pandas as pd
        from Utilities.overlap_deduplicator import OverlapDeduplicator

        dedup = OverlapDeduplicator()
        df = pd.DataFrame()
        result = dedup.filter_core_df(df, core_end=48_000)
        assert result is not None
        assert len(result) == 0


# ──────────────────────────────────────────────────────────────────────────────
# FinalExporter
# ──────────────────────────────────────────────────────────────────────────────

class TestFinalExporter:
    """Unit tests for FinalExporter using a mock DiskChunkManager."""

    def _make_disk_manager(self, tmp_path, chunks_data):
        """
        Create a real DiskChunkManager populated with test data.

        chunks_data: list of (chunk_index, motifs, core_end) tuples.
        """
        from Utilities.disk_chunk_manager import DiskChunkManager

        mgr = DiskChunkManager(base_dir=str(tmp_path))
        for chunk_index, motifs, core_end in chunks_data:
            meta = mgr.write_chunk_results(chunk_index=chunk_index, motifs=motifs)
            # Inject core_end into stored metadata for FinalExporter to read
            mgr._chunk_registry[chunk_index]["core_end"] = core_end
        return mgr

    def _motifs(self, starts, cls="G-Quadruplex"):
        return [
            {"Class": cls, "Subclass": "G4",
             "Start": s, "End": s + 20, "Length": 20,
             "Score": 0.9, "ID": f"seq_G4_{s}", "Strand": "+"}
            for s in starts
        ]

    def test_assemble_yields_filtered_chunks(self, tmp_path):
        from Utilities.final_exporter import FinalExporter
        from Utilities.overlap_deduplicator import OverlapDeduplicator

        # Chunk 0: motifs at 100, 47_999, 48_000 (last two in overlap)
        mgr = self._make_disk_manager(tmp_path, [
            (0, self._motifs([100, 47_999, 48_000]), 48_000),
        ])

        exporter = FinalExporter()
        dedup = OverlapDeduplicator()
        results = list(exporter.assemble(mgr, dedup))

        assert len(results) == 1
        # Only motifs with Start < 48_000 are kept
        assert all(m["Start"] < 48_000 for m in results[0])
        assert len(results[0]) == 2
        mgr.cleanup()

    def test_assemble_multiple_chunks(self, tmp_path):
        from Utilities.final_exporter import FinalExporter
        from Utilities.overlap_deduplicator import OverlapDeduplicator

        mgr = self._make_disk_manager(tmp_path, [
            (0, self._motifs([100, 200, 300]), 48_000),
            (1, self._motifs([48_500, 50_000]), 96_000),
        ])

        exporter = FinalExporter()
        dedup = OverlapDeduplicator()
        results = list(exporter.assemble(mgr, dedup))
        assert len(results) == 2
        assert len(results[0]) == 3
        assert len(results[1]) == 2
        mgr.cleanup()

    def test_to_dataframe_returns_dataframe(self, tmp_path):
        pytest.importorskip("pandas")
        import pandas as pd
        from Utilities.final_exporter import FinalExporter
        from Utilities.overlap_deduplicator import OverlapDeduplicator

        mgr = self._make_disk_manager(tmp_path, [
            (0, self._motifs([100, 200, 300]), 1_000_000),
        ])

        exporter = FinalExporter()
        dedup = OverlapDeduplicator()
        df = exporter.to_dataframe(mgr, dedup)

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3
        mgr.cleanup()

    def test_to_dataframe_empty_when_no_motifs(self, tmp_path):
        pytest.importorskip("pandas")
        import pandas as pd
        from Utilities.final_exporter import FinalExporter
        from Utilities.overlap_deduplicator import OverlapDeduplicator

        mgr = self._make_disk_manager(tmp_path, [
            (0, [], 48_000),
        ])

        exporter = FinalExporter()
        dedup = OverlapDeduplicator()
        df = exporter.to_dataframe(mgr, dedup)

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0
        mgr.cleanup()

    def test_assemble_with_results_storage(self, tmp_path):
        """FinalExporter works with UniversalResultsStorage as well."""
        from Utilities.disk_storage import UniversalResultsStorage
        from Utilities.final_exporter import FinalExporter
        from Utilities.overlap_deduplicator import OverlapDeduplicator

        rs = UniversalResultsStorage(
            base_dir=str(tmp_path / "results"), seq_id="test"
        )
        for s in [100, 200, 300]:
            rs.append({"Class": "Z-DNA", "Subclass": "ZA",
                       "Start": s, "End": s + 15, "Length": 15,
                       "Score": 0.7, "ID": f"x_{s}", "Strand": "+"})

        exporter = FinalExporter()
        dedup = OverlapDeduplicator()
        results = list(exporter.assemble(rs, dedup))
        # All motifs yielded as single batch (no chunk_meta → no dedup filtering)
        assert len(results) == 1
        assert len(results[0]) == 3


# ──────────────────────────────────────────────────────────────────────────────
# DetectorRunner – interface / metadata only
# ──────────────────────────────────────────────────────────────────────────────

class TestDetectorRunnerInterface:
    """Tests for DetectorRunner initialisation and cleanup (no nonbscanner)."""

    def test_creates_tmp_dir_automatically(self):
        from Utilities.detector_runner import DetectorRunner
        runner = DetectorRunner()
        assert Path(runner._tmp_dir).is_dir()
        runner.cleanup()
        # After cleanup, the auto-created dir should be removed
        assert not Path(runner._tmp_dir).exists()

    def test_uses_provided_tmp_dir(self, tmp_path):
        from Utilities.detector_runner import DetectorRunner
        runner = DetectorRunner(tmp_dir=str(tmp_path / "chunks"))
        assert Path(runner._tmp_dir).is_dir()
        # Should NOT delete caller-provided directory on cleanup
        runner.cleanup()
        assert Path(runner._tmp_dir).exists()

    def test_seq_name_stored(self):
        from Utilities.detector_runner import DetectorRunner
        runner = DetectorRunner(seq_name="chr1")
        assert runner.seq_name == "chr1"
        runner.cleanup()


# ──────────────────────────────────────────────────────────────────────────────
# ParallelChunkExecutor – worker cap
# ──────────────────────────────────────────────────────────────────────────────

class TestParallelChunkExecutorWorkerCap:
    """Verify worker count is capped at MAX_WORKERS."""

    def test_workers_never_exceed_max(self):
        from Utilities.detector_runner import DetectorRunner, ParallelChunkExecutor, MAX_WORKERS
        runner = DetectorRunner()
        # Request many workers; should be capped
        executor = ParallelChunkExecutor(detector_runner=runner, workers=100)
        assert executor.workers <= MAX_WORKERS
        runner.cleanup()

    def test_max_workers_at_most_2(self):
        from Utilities.detector_runner import MAX_WORKERS
        assert MAX_WORKERS <= 2

    def test_default_workers_matches_max_workers(self):
        from Utilities.detector_runner import DetectorRunner, ParallelChunkExecutor, MAX_WORKERS
        runner = DetectorRunner()
        executor = ParallelChunkExecutor(detector_runner=runner)
        assert executor.workers == MAX_WORKERS
        runner.cleanup()
