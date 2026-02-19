"""
Tests for the new chunk-based parallel execution modules:

  - SequenceContext
  - PerformanceMonitor
  - ChunkExecutor
  - ParallelVisualization
  - BaseMotifDetector.scan(context) method
"""

import time
import pytest
from unittest.mock import MagicMock, patch


# ──────────────────────────────────────────────────────────────────────────────
# SequenceContext
# ──────────────────────────────────────────────────────────────────────────────

class TestSequenceContext:
    """Unit tests for SequenceContext preprocessing."""

    def test_uppercases_sequence(self):
        from Utilities.core.sequence_context import SequenceContext
        ctx = SequenceContext("atcg", "s1")
        assert ctx.sequence == "ATCG"

    def test_strips_whitespace(self):
        from Utilities.core.sequence_context import SequenceContext
        ctx = SequenceContext("  ATCG  ", "s1")
        assert ctx.sequence == "ATCG"

    def test_length_cached(self):
        from Utilities.core.sequence_context import SequenceContext
        seq = "GGGTTTAGGG" * 10
        ctx = SequenceContext(seq)
        assert ctx.length == len(seq)
        # Length should equal len(sequence) exactly
        assert len(ctx) == ctx.length

    def test_name_stored(self):
        from Utilities.core.sequence_context import SequenceContext
        ctx = SequenceContext("ATCG", name="chr1")
        assert ctx.name == "chr1"

    def test_default_name(self):
        from Utilities.core.sequence_context import SequenceContext
        ctx = SequenceContext("ATCG")
        assert ctx.name == "sequence"

    def test_gc_content_without_precompute(self):
        from Utilities.core.sequence_context import SequenceContext
        ctx = SequenceContext("GCGCAT")   # 4 GC out of 6
        assert abs(ctx.gc_content() - 4 / 6) < 1e-9

    def test_gc_content_with_precompute(self):
        from Utilities.core.sequence_context import SequenceContext
        ctx = SequenceContext("GCGCAT", precompute_gc=True)
        assert ctx.gc_prefix is not None
        assert abs(ctx.gc_content() - 4 / 6) < 1e-9

    def test_gc_content_subregion(self):
        from Utilities.core.sequence_context import SequenceContext
        # "ATCG" → 0 GC in [0:2] ("AT"), 2 GC in [2:4] ("CG")
        ctx = SequenceContext("ATCG", precompute_gc=True)
        assert ctx.gc_content(0, 2) == 0.0   # "AT"  – 0 GC in 2 bases
        assert ctx.gc_content(2, 4) == 1.0   # "CG"  – 2 GC in 2 bases

    def test_repr_contains_name_and_length(self):
        from Utilities.core.sequence_context import SequenceContext
        ctx = SequenceContext("ATCG", name="myseq")
        r = repr(ctx)
        assert "myseq" in r
        assert "4" in r

    def test_empty_sequence(self):
        from Utilities.core.sequence_context import SequenceContext
        ctx = SequenceContext("")
        assert ctx.length == 0
        assert ctx.gc_content() == 0.0


# ──────────────────────────────────────────────────────────────────────────────
# PerformanceMonitor
# ──────────────────────────────────────────────────────────────────────────────

class TestPerformanceMonitor:
    """Unit tests for PerformanceMonitor."""

    def test_start_and_elapsed(self):
        from Utilities.core.performance_monitor import PerformanceMonitor
        monitor = PerformanceMonitor()
        monitor.start()
        time.sleep(0.05)
        summary = monitor.get_summary()
        assert summary["total_elapsed"] >= 0.04

    def test_record_chunk(self):
        from Utilities.core.performance_monitor import PerformanceMonitor
        monitor = PerformanceMonitor()
        monitor.start()
        monitor.record_chunk(chunk_id=0, elapsed=1.0, motif_count=10, bp_count=50_000)
        summary = monitor.get_summary()
        assert summary["chunk_count"] == 1
        assert summary["chunk_records"][0]["motif_count"] == 10
        assert summary["chunk_records"][0]["bp_count"] == 50_000
        assert summary["total_bp_processed"] == 50_000

    def test_record_detector(self):
        from Utilities.core.performance_monitor import PerformanceMonitor
        monitor = PerformanceMonitor()
        monitor.start()
        monitor.record_detector(chunk_id=0, detector_name="g_quadruplex",
                                elapsed=0.3, bp_count=50_000)
        monitor.record_detector(chunk_id=0, detector_name="g_quadruplex",
                                elapsed=0.2, bp_count=50_000)
        summary = monitor.get_summary()
        bd = summary["detector_breakdown"]
        assert "g_quadruplex" in bd
        assert bd["g_quadruplex"]["call_count"] == 2
        assert abs(bd["g_quadruplex"]["total_elapsed"] - 0.5) < 1e-9

    def test_record_stage(self):
        from Utilities.core.performance_monitor import PerformanceMonitor
        monitor = PerformanceMonitor()
        monitor.start()
        monitor.record_stage("detection", 3.5)
        summary = monitor.get_summary()
        assert summary["stage_times"]["detection"] == 3.5

    def test_slowest_detector(self):
        from Utilities.core.performance_monitor import PerformanceMonitor
        monitor = PerformanceMonitor()
        monitor.start()
        monitor.record_detector(0, "fast_det", 0.1)
        monitor.record_detector(0, "slow_det", 5.0)
        summary = monitor.get_summary()
        assert summary["slowest_detector"] == "slow_det"

    def test_throughput_zero_before_bp_recorded(self):
        from Utilities.core.performance_monitor import PerformanceMonitor
        monitor = PerformanceMonitor()
        monitor.start()
        summary = monitor.get_summary()
        assert summary["throughput_bps"] == 0.0

    def test_throughput_positive_after_chunk(self):
        from Utilities.core.performance_monitor import PerformanceMonitor
        monitor = PerformanceMonitor()
        monitor.start()
        time.sleep(0.01)
        monitor.record_chunk(0, elapsed=1.0, motif_count=5, bp_count=100_000)
        summary = monitor.get_summary()
        assert summary["throughput_bps"] > 0

    def test_snapshot_memory_returns_float(self):
        from Utilities.core.performance_monitor import PerformanceMonitor
        monitor = PerformanceMonitor()
        mem = monitor.snapshot_memory()
        assert isinstance(mem, float)
        # May be 0.0 if psutil unavailable, but should be non-negative
        assert mem >= 0.0

    def test_format_summary_contains_key_fields(self):
        from Utilities.core.performance_monitor import PerformanceMonitor
        monitor = PerformanceMonitor()
        monitor.start()
        monitor.record_chunk(0, 1.0, 20, 50_000)
        monitor.record_detector(0, "z_dna", 0.5)
        monitor.record_stage("detection", 1.0)
        text = monitor.format_summary()
        assert "Total runtime" in text
        assert "z_dna" in text
        assert "detection" in text

    def test_thread_safe_concurrent_writes(self):
        """Multiple threads can safely call record_* concurrently."""
        import threading
        from Utilities.core.performance_monitor import PerformanceMonitor

        monitor = PerformanceMonitor()
        monitor.start()

        def _write(i):
            monitor.record_chunk(i, 0.1, 1, 1_000)
            monitor.record_detector(i, "det_a", 0.05)

        threads = [threading.Thread(target=_write, args=(i,)) for i in range(20)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        summary = monitor.get_summary()
        assert summary["chunk_count"] == 20
        assert summary["detector_breakdown"]["det_a"]["call_count"] == 20


# ──────────────────────────────────────────────────────────────────────────────
# BaseMotifDetector.scan(context)
# ──────────────────────────────────────────────────────────────────────────────

class TestBaseDetectorScan:
    """scan(context) should delegate to detect_motifs with the uppercase sequence."""

    def test_scan_calls_detect_motifs(self):
        from Utilities.core.sequence_context import SequenceContext
        from Detectors.base.base_detector import BaseMotifDetector

        # Minimal concrete subclass for testing
        class _DummyDetector(BaseMotifDetector):
            def get_patterns(self):
                return {}

            def get_motif_class_name(self):
                return "Dummy"

            def calculate_score(self, sequence, pattern_info):
                return 1.0

        det = _DummyDetector()
        ctx = SequenceContext("atcg", name="test_seq")

        # Patch detect_motifs to track calls
        called_with = {}

        def _fake_detect(sequence, name):
            called_with["sequence"] = sequence
            called_with["name"] = name
            return []

        det.detect_motifs = _fake_detect
        det.scan(ctx)

        assert called_with["sequence"] == "ATCG"
        assert called_with["name"] == "test_seq"


# ──────────────────────────────────────────────────────────────────────────────
# ChunkExecutor
# ──────────────────────────────────────────────────────────────────────────────

class TestChunkExecutorBasic:
    """Basic integration tests for ChunkExecutor (no ProcessPool)."""

    def _make_seq(self, length: int) -> str:
        unit = "GGGTTTAGGGTTT"
        spacer = "ATCGATCG" * 10
        seg = unit + spacer
        return (seg * (length // len(seg) + 1))[:length]

    def test_run_returns_results_storage(self, tmp_path):
        from Utilities.disk_storage import UniversalSequenceStorage
        from Utilities.core.chunk_executor import ChunkExecutor

        storage = UniversalSequenceStorage(base_dir=str(tmp_path / "seqs"))
        seq = self._make_seq(60_000)
        seq_id = storage.save_sequence(seq, "test_seq")

        executor = ChunkExecutor(storage, chunk_size=50_000, overlap=2_000, max_workers=1)
        rs = executor.run(seq_id=seq_id, enabled_classes=["G-Quadruplex"])

        stats = rs.get_summary_stats()
        assert stats["total_count"] >= 0   # may be 0 for very short sequences

        storage.cleanup()

    def test_performance_summary_available_after_run(self, tmp_path):
        from Utilities.disk_storage import UniversalSequenceStorage
        from Utilities.core.chunk_executor import ChunkExecutor

        storage = UniversalSequenceStorage(base_dir=str(tmp_path / "seqs"))
        seq = self._make_seq(55_000)
        seq_id = storage.save_sequence(seq, "ts")

        executor = ChunkExecutor(storage, chunk_size=50_000, overlap=2_000, max_workers=1)
        executor.run(seq_id=seq_id, enabled_classes=["G-Quadruplex"])

        perf = executor.get_performance_summary()
        assert perf is not None
        assert "total_elapsed" in perf
        assert "chunk_count" in perf

        storage.cleanup()

    def test_legacy_float_callback_works(self, tmp_path):
        """A legacy callback(float) must not raise."""
        from Utilities.disk_storage import UniversalSequenceStorage
        from Utilities.core.chunk_executor import ChunkExecutor

        storage = UniversalSequenceStorage(base_dir=str(tmp_path / "seqs"))
        seq = self._make_seq(55_000)
        seq_id = storage.save_sequence(seq, "ts")

        received = []

        def _legacy_cb(pct: float):
            received.append(pct)

        executor = ChunkExecutor(storage, chunk_size=50_000, overlap=2_000, max_workers=1)
        executor.run(seq_id=seq_id, enabled_classes=["G-Quadruplex"],
                     progress_callback=_legacy_cb)

        # Callback should have been fired at least once
        assert len(received) > 0

        storage.cleanup()

    def test_rich_dict_callback_works(self, tmp_path):
        """A rich callback(dict) must receive telemetry keys."""
        from Utilities.disk_storage import UniversalSequenceStorage
        from Utilities.core.chunk_executor import ChunkExecutor

        storage = UniversalSequenceStorage(base_dir=str(tmp_path / "seqs"))
        seq = self._make_seq(55_000)
        seq_id = storage.save_sequence(seq, "ts")

        events = []

        def _rich_cb(event):
            events.append(event)

        executor = ChunkExecutor(storage, chunk_size=50_000, overlap=2_000, max_workers=1)
        executor.run(seq_id=seq_id, enabled_classes=["G-Quadruplex"],
                     progress_callback=_rich_cb)

        assert len(events) > 0
        first = events[0]
        assert "stage" in first
        assert "progress_pct" in first
        assert first["stage"] == "detection"

        storage.cleanup()


class TestChunkExecutorFallback:
    """When ProcessPoolExecutor raises, ChunkExecutor falls back to sequential."""

    def _make_seq(self, length: int) -> str:
        unit = "GGGTTTAGGGTTT"
        spacer = "ATCGATCG" * 10
        seg = unit + spacer
        return (seg * (length // len(seg) + 1))[:length]

    def test_fallback_sequential_still_returns_results(self, tmp_path):
        from Utilities.disk_storage import UniversalSequenceStorage
        from Utilities.core.chunk_executor import ChunkExecutor

        storage = UniversalSequenceStorage(base_dir=str(tmp_path / "seqs"))
        seq = self._make_seq(60_000)
        seq_id = storage.save_sequence(seq, "ts")

        def _raise_on_create(*a, **kw):
            raise RuntimeError("Simulated pool failure")

        with patch("Utilities.core.chunk_executor.ProcessPoolExecutor",
                   side_effect=_raise_on_create):
            executor = ChunkExecutor(storage, chunk_size=50_000, overlap=2_000, max_workers=1)
            rs = executor.run(seq_id=seq_id, enabled_classes=["G-Quadruplex"])

        # Should have results (may be 0 but storage object must exist)
        assert rs is not None
        assert rs.get_summary_stats()["total_count"] >= 0

        storage.cleanup()


# ──────────────────────────────────────────────────────────────────────────────
# ParallelVisualization
# ──────────────────────────────────────────────────────────────────────────────

class TestParallelVisualization:
    """Unit tests for ParallelVisualization."""

    def _make_pipeline(self):
        """Return a VisualizationPipeline backed by minimal test data."""
        import numpy as np
        from Utilities.visualization_pipeline import VisualizationPipeline

        summary = {
            "total_motifs": 5,
            "class_counts": {"G-Quadruplex": 3, "Z-DNA": 2},
            "subclass_counts": {"G4_classic": 3, "Z-DNA_CG": 2},
            "density_bins": np.zeros(100, dtype=np.int64),
            "length_bins": np.zeros(100, dtype=np.int64),
            "cooccurrence_matrix": {"G-Quadruplex": {"Z-DNA": 1}, "Z-DNA": {"G-Quadruplex": 1}},
            "bin_count": 100,
            "seq_length": 1_000,
            "max_length": 10_000,
        }
        return VisualizationPipeline(summary)

    def test_generate_all_parallel_returns_five_figures(self):
        from Utilities.visualization.parallel_visualization import ParallelVisualization

        pipeline = self._make_pipeline()
        runner = ParallelVisualization(pipeline, max_workers=2)
        result = runner.generate_all()

        assert "figures" in result
        assert "timings" in result
        assert "total_elapsed" in result
        assert len(result["figures"]) == 5
        assert all(k in result["figures"] for k in [
            "density_histogram", "length_histogram",
            "class_distribution", "subclass_distribution",
            "cooccurrence_heatmap",
        ])

    def test_all_figures_are_not_none(self):
        from Utilities.visualization.parallel_visualization import ParallelVisualization

        pipeline = self._make_pipeline()
        runner = ParallelVisualization(pipeline, max_workers=2)
        result = runner.generate_all()

        for name, fig in result["figures"].items():
            assert fig is not None, f"Figure '{name}' should not be None"

    def test_progress_callback_fired(self):
        from Utilities.visualization.parallel_visualization import ParallelVisualization

        pipeline = self._make_pipeline()
        events = []
        runner = ParallelVisualization(
            pipeline, max_workers=2, progress_callback=lambda e: events.append(e)
        )
        runner.generate_all()

        assert len(events) == 5
        for ev in events:
            assert ev["stage"] == "visualization"
            assert "component" in ev
            assert "elapsed" in ev

    def test_generate_all_parallel_method_on_pipeline(self):
        """VisualizationPipeline.generate_all_parallel() is a convenience wrapper."""
        pipeline = self._make_pipeline()
        result = pipeline.generate_all_parallel(max_workers=2)

        assert "figures" in result
        assert len(result["figures"]) == 5

    def test_timings_keys_match_figures(self):
        from Utilities.visualization.parallel_visualization import ParallelVisualization

        pipeline = self._make_pipeline()
        runner = ParallelVisualization(pipeline, max_workers=2)
        result = runner.generate_all()

        assert set(result["timings"].keys()) == set(result["figures"].keys())

    def test_total_elapsed_is_positive(self):
        from Utilities.visualization.parallel_visualization import ParallelVisualization

        pipeline = self._make_pipeline()
        runner = ParallelVisualization(pipeline, max_workers=2)
        result = runner.generate_all()

        assert result["total_elapsed"] > 0.0


# ──────────────────────────────────────────────────────────────────────────────
# MotifDetectionEngine integration
# ──────────────────────────────────────────────────────────────────────────────

class TestMotifDetectionEngineWithChunkExecutor:
    """Engine should use ChunkExecutor for medium sequences."""

    def _make_seq(self, length: int) -> str:
        unit = "GGGTTTAGGGTTT"
        spacer = "ATCGATCG" * 10
        seg = unit + spacer
        return (seg * (length // len(seg) + 1))[:length]

    def test_engine_returns_results_and_accumulator(self, tmp_path):
        from Utilities.disk_storage import UniversalSequenceStorage
        from Utilities.motif_detection_engine import MotifDetectionEngine

        storage = UniversalSequenceStorage(base_dir=str(tmp_path / "seqs"))
        seq = self._make_seq(200_000)   # medium: triggers chunk_workers strategy
        seq_id = storage.save_sequence(seq, "eng_test")

        engine = MotifDetectionEngine(storage, chunk_size=50_000, overlap=2_000, max_workers=1)
        rs, acc = engine.analyze(seq_id=seq_id, enabled_classes=["G-Quadruplex"])

        assert rs is not None
        assert acc is not None
        assert rs.get_summary_stats()["total_count"] >= 0

        storage.cleanup()

    def test_engine_get_performance_summary_after_chunk_strategy(self, tmp_path):
        from Utilities.disk_storage import UniversalSequenceStorage
        from Utilities.motif_detection_engine import MotifDetectionEngine

        storage = UniversalSequenceStorage(base_dir=str(tmp_path / "seqs"))
        seq = self._make_seq(200_000)
        seq_id = storage.save_sequence(seq, "eng_test2")

        engine = MotifDetectionEngine(storage, chunk_size=50_000, overlap=2_000, max_workers=1)
        engine.analyze(seq_id=seq_id, enabled_classes=["G-Quadruplex"])

        perf = engine.get_performance_summary()
        # ChunkExecutor provides performance summary
        assert perf is not None
        assert "total_elapsed" in perf

        storage.cleanup()
