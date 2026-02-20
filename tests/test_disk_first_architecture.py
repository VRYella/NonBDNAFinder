"""
Tests for the disk-first parallel architecture components:
  - DiskChunkManager
  - VisualizationAccumulator
  - StreamlitSafeExecutor (strategy selection only; full integration skipped
    without pandas/nonbscanner)
  - MotifDetectionEngine (strategy selection only)
  - VisualizationPipeline
"""

import pytest
import numpy as np
from pathlib import Path


# ──────────────────────────────────────────────────────────────────────────────
# DiskChunkManager
# ──────────────────────────────────────────────────────────────────────────────

class TestDiskChunkManager:
    """Unit tests for DiskChunkManager."""

    def test_write_and_iter(self, tmp_path):
        from Utilities.disk_chunk_manager import DiskChunkManager

        mgr = DiskChunkManager(base_dir=str(tmp_path))
        motifs = [
            {"Class": "G-Quadruplex", "Subclass": "G4", "Start": 100, "End": 120,
             "Length": 20, "Score": 0.9, "ID": "seq_G4_100", "Strand": "+"},
            {"Class": "Z-DNA", "Subclass": "ZA", "Start": 200, "End": 215,
             "Length": 15, "Score": 0.7, "ID": "seq_ZA_200", "Strand": "-"},
        ]
        meta = mgr.write_chunk_results(chunk_index=0, motifs=motifs)

        assert meta["chunk_index"] == 0
        assert meta["motif_count"] == 2
        assert Path(meta["file_path"]).exists()

        results = list(mgr.iter_chunk_results(delete_after_read=False))
        assert len(results) == 1
        read_motifs, read_meta = results[0]
        assert len(read_motifs) == 2
        assert read_motifs[0]["Class"] == "G-Quadruplex"
        assert read_motifs[1]["Start"] == 200

        mgr.cleanup()

    def test_numeric_fields_restored(self, tmp_path):
        from Utilities.disk_chunk_manager import DiskChunkManager

        mgr = DiskChunkManager(base_dir=str(tmp_path))
        motifs = [
            {"Class": "R-Loop", "Subclass": "RL", "Start": 50, "End": 75,
             "Length": 25, "Score": 1.2, "ID": "seq_RL_50", "Strand": "+"},
        ]
        mgr.write_chunk_results(chunk_index=0, motifs=motifs)

        read_motifs, _ = next(iter(mgr.iter_chunk_results(delete_after_read=False)))
        assert isinstance(read_motifs[0]["Start"], int)
        assert isinstance(read_motifs[0]["End"], int)
        assert isinstance(read_motifs[0]["Score"], float)
        mgr.cleanup()

    def test_delete_after_read(self, tmp_path):
        from Utilities.disk_chunk_manager import DiskChunkManager

        mgr = DiskChunkManager(base_dir=str(tmp_path))
        motifs = [{"Class": "G-Quadruplex", "Subclass": "G4", "Start": 0,
                   "End": 10, "Length": 10, "Score": 0.5, "ID": "x", "Strand": "+"}]
        meta = mgr.write_chunk_results(chunk_index=0, motifs=motifs)
        file_path = Path(meta["file_path"])

        assert file_path.exists()
        list(mgr.iter_chunk_results(delete_after_read=True))
        assert not file_path.exists()
        mgr.cleanup()

    def test_multiple_chunks_ordered(self, tmp_path):
        from Utilities.disk_chunk_manager import DiskChunkManager

        mgr = DiskChunkManager(base_dir=str(tmp_path))
        for i in range(3):
            mgr.write_chunk_results(
                chunk_index=i,
                motifs=[{"Class": "G-Quadruplex", "Subclass": "G4",
                         "Start": i * 100, "End": i * 100 + 10, "Length": 10,
                         "Score": 0.5, "ID": f"x_{i}", "Strand": "+"}],
            )

        indices = [meta["chunk_index"] for _, meta in mgr.iter_chunk_results(delete_after_read=False)]
        assert indices == [0, 1, 2]
        assert mgr.total_motif_count() == 3
        mgr.cleanup()

    def test_cleanup_removes_directory(self, tmp_path):
        from Utilities.disk_chunk_manager import DiskChunkManager

        mgr = DiskChunkManager(base_dir=str(tmp_path))
        base = mgr.base_dir
        assert base.exists()
        mgr.cleanup()
        assert not base.exists()


# ──────────────────────────────────────────────────────────────────────────────
# VisualizationAccumulator
# ──────────────────────────────────────────────────────────────────────────────

class TestVisualizationAccumulator:
    """Unit tests for VisualizationAccumulator."""

    def _make_motifs(self, n=5, motif_class="G-Quadruplex", subcls="G4", start_offset=0):
        return [
            {"Class": motif_class, "Subclass": subcls,
             "Start": start_offset + i * 100, "End": start_offset + i * 100 + 20,
             "Length": 20, "Score": 0.8}
            for i in range(n)
        ]

    def test_class_counts_accumulated(self):
        from Utilities.visualization_accumulator import VisualizationAccumulator

        acc = VisualizationAccumulator(seq_length=10_000)
        acc.update(self._make_motifs(3, motif_class="G-Quadruplex"))
        acc.update(self._make_motifs(2, motif_class="Z-DNA"))

        summary = acc.get_summary()
        assert summary["class_counts"]["G-Quadruplex"] == 3
        assert summary["class_counts"]["Z-DNA"] == 2
        assert summary["total_motifs"] == 5

    def test_subclass_counts(self):
        from Utilities.visualization_accumulator import VisualizationAccumulator

        acc = VisualizationAccumulator(seq_length=10_000)
        acc.update(self._make_motifs(4, motif_class="G-Quadruplex", subcls="G4"))
        acc.update(self._make_motifs(2, motif_class="G-Quadruplex", subcls="G2"))

        summary = acc.get_summary()
        assert summary["subclass_counts"]["G4"] == 4
        assert summary["subclass_counts"]["G2"] == 2

    def test_density_bins_fixed_size(self):
        from Utilities.visualization_accumulator import VisualizationAccumulator

        acc = VisualizationAccumulator(seq_length=100_000, bin_count=50)
        acc.update(self._make_motifs(10))

        summary = acc.get_summary()
        assert summary["density_bins"].shape == (50,)
        assert summary["density_bins"].sum() == 10

    def test_length_bins_fixed_size(self):
        from Utilities.visualization_accumulator import VisualizationAccumulator

        acc = VisualizationAccumulator(seq_length=10_000, bin_count=20)
        acc.update(self._make_motifs(5))

        summary = acc.get_summary()
        assert summary["length_bins"].shape == (20,)
        assert summary["length_bins"].sum() == 5

    def test_cooccurrence_matrix(self):
        from Utilities.visualization_accumulator import VisualizationAccumulator

        acc = VisualizationAccumulator(seq_length=10_000)
        batch = (
            self._make_motifs(2, motif_class="G-Quadruplex") +
            self._make_motifs(2, motif_class="Z-DNA")
        )
        acc.update(batch)

        summary = acc.get_summary()
        cooc = summary["cooccurrence_matrix"]
        # Both classes appear in the same batch → co-occurrence should be > 0
        assert cooc.get("G-Quadruplex", {}).get("Z-DNA", 0) > 0
        assert cooc.get("Z-DNA", {}).get("G-Quadruplex", 0) > 0

    def test_empty_update_is_safe(self):
        from Utilities.visualization_accumulator import VisualizationAccumulator

        acc = VisualizationAccumulator(seq_length=1_000)
        acc.update([])
        summary = acc.get_summary()
        assert summary["total_motifs"] == 0
        assert summary["density_bins"].sum() == 0

    def test_reset(self):
        from Utilities.visualization_accumulator import VisualizationAccumulator

        acc = VisualizationAccumulator(seq_length=10_000)
        acc.update(self._make_motifs(5))
        acc.reset()

        summary = acc.get_summary()
        assert summary["total_motifs"] == 0
        assert summary["class_counts"] == {}
        assert summary["density_bins"].sum() == 0

    def test_missing_length_field_uses_end_minus_start(self):
        from Utilities.visualization_accumulator import VisualizationAccumulator

        acc = VisualizationAccumulator(seq_length=10_000)
        motifs = [{"Class": "G-Quadruplex", "Subclass": "G4",
                   "Start": 100, "End": 125, "Score": 0.5}]
        acc.update(motifs)

        summary = acc.get_summary()
        assert summary["total_motifs"] == 1


# ──────────────────────────────────────────────────────────────────────────────
# StreamlitSafeExecutor – strategy selection only
# ──────────────────────────────────────────────────────────────────────────────

class TestStreamlitSafeExecutorStrategySelection:
    """Verify adaptive strategy thresholds without running analysis."""

    def test_direct_below_threshold(self):
        from Utilities.streamlit_safe_executor import StreamlitSafeExecutor

        # THRESHOLD_DIRECT = 0 means all sequences use chunking (no direct path).
        # Sequences that were previously "direct" are now "chunk_workers" so that
        # 50K/2K chunking is applied universally for every sequence.
        assert StreamlitSafeExecutor._select_strategy(0) == "chunk_workers"
        assert StreamlitSafeExecutor._select_strategy(99_999) == "chunk_workers"

    def test_chunk_workers_middle_band(self):
        from Utilities.streamlit_safe_executor import StreamlitSafeExecutor

        assert StreamlitSafeExecutor._select_strategy(100_000) == "chunk_workers"
        assert StreamlitSafeExecutor._select_strategy(4_999_999) == "chunk_workers"

    def test_disk_streaming_large(self):
        from Utilities.streamlit_safe_executor import StreamlitSafeExecutor

        assert StreamlitSafeExecutor._select_strategy(5_000_000) == "disk_streaming"
        assert StreamlitSafeExecutor._select_strategy(100_000_000) == "disk_streaming"

    def test_max_workers_capped_at_2(self):
        from Utilities.streamlit_safe_executor import MAX_WORKERS

        assert MAX_WORKERS <= 2


# ──────────────────────────────────────────────────────────────────────────────
# MotifDetectionEngine – strategy selection only
# ──────────────────────────────────────────────────────────────────────────────

class TestMotifDetectionEngineStrategySelection:
    """Verify adaptive strategy thresholds without running analysis."""

    def test_numba_only_below_100k(self):
        from Utilities.motif_detection_engine import MotifDetectionEngine

        assert MotifDetectionEngine._select_strategy(50_000) == "numba_only"
        assert MotifDetectionEngine._select_strategy(99_999) == "numba_only"

    def test_chunk_workers_100k_to_5m(self):
        from Utilities.motif_detection_engine import MotifDetectionEngine

        assert MotifDetectionEngine._select_strategy(100_000) == "chunk_workers"
        assert MotifDetectionEngine._select_strategy(4_999_999) == "chunk_workers"

    def test_disk_streaming_above_5m(self):
        from Utilities.motif_detection_engine import MotifDetectionEngine

        assert MotifDetectionEngine._select_strategy(5_000_000) == "disk_streaming"
        assert MotifDetectionEngine._select_strategy(50_000_000) == "disk_streaming"


# ──────────────────────────────────────────────────────────────────────────────
# VisualizationPipeline
# ──────────────────────────────────────────────────────────────────────────────

class TestVisualizationPipeline:
    """Tests for VisualizationPipeline using synthetic summary data."""

    def _make_summary(self, n_classes=3, bin_count=10):
        from Utilities.visualization_accumulator import VisualizationAccumulator

        acc = VisualizationAccumulator(seq_length=100_000, bin_count=bin_count)
        classes = ["G-Quadruplex", "Z-DNA", "R-Loop"][:n_classes]
        for i, cls in enumerate(classes):
            motifs = [
                {"Class": cls, "Subclass": cls[:2], "Start": i * 1000 + j * 50,
                 "End": i * 1000 + j * 50 + 20, "Length": 20, "Score": 0.7}
                for j in range(5)
            ]
            acc.update(motifs)
        return acc.get_summary()

    def test_density_histogram_returns_figure(self):
        from Utilities.visualization_pipeline import VisualizationPipeline
        import matplotlib.figure

        pipeline = VisualizationPipeline(self._make_summary())
        fig = pipeline.plot_density_histogram()
        assert isinstance(fig, matplotlib.figure.Figure)

    def test_length_histogram_returns_figure(self):
        from Utilities.visualization_pipeline import VisualizationPipeline
        import matplotlib.figure

        pipeline = VisualizationPipeline(self._make_summary())
        fig = pipeline.plot_length_histogram()
        assert isinstance(fig, matplotlib.figure.Figure)

    def test_class_distribution_returns_figure(self):
        from Utilities.visualization_pipeline import VisualizationPipeline
        import matplotlib.figure

        pipeline = VisualizationPipeline(self._make_summary())
        fig = pipeline.plot_class_distribution()
        assert isinstance(fig, matplotlib.figure.Figure)

    def test_subclass_distribution_returns_figure(self):
        from Utilities.visualization_pipeline import VisualizationPipeline
        import matplotlib.figure

        pipeline = VisualizationPipeline(self._make_summary())
        fig = pipeline.plot_subclass_distribution()
        assert isinstance(fig, matplotlib.figure.Figure)

    def test_cooccurrence_heatmap_returns_figure(self):
        from Utilities.visualization_pipeline import VisualizationPipeline
        import matplotlib.figure

        pipeline = VisualizationPipeline(self._make_summary())
        fig = pipeline.plot_cooccurrence_heatmap()
        assert isinstance(fig, matplotlib.figure.Figure)

    def test_generate_all_returns_five_figures(self):
        from Utilities.visualization_pipeline import VisualizationPipeline
        import matplotlib.figure

        pipeline = VisualizationPipeline(self._make_summary())
        figs = pipeline.generate_all()
        assert set(figs.keys()) == {
            "density_histogram",
            "length_histogram",
            "class_distribution",
            "subclass_distribution",
            "cooccurrence_heatmap",
        }
        for name, fig in figs.items():
            assert isinstance(fig, matplotlib.figure.Figure), f"{name} is not a Figure"

    def test_empty_summary_returns_blank_figures(self):
        from Utilities.visualization_pipeline import VisualizationPipeline
        import matplotlib.figure

        empty_summary = {
            "total_motifs": 0,
            "class_counts": {},
            "subclass_counts": {},
            "density_bins": np.zeros(10, dtype=np.int64),
            "length_bins": np.zeros(10, dtype=np.int64),
            "cooccurrence_matrix": {},
            "bin_count": 10,
            "seq_length": 1000,
            "max_length": 10_000,
        }
        pipeline = VisualizationPipeline(empty_summary)
        fig = pipeline.plot_class_distribution()
        assert isinstance(fig, matplotlib.figure.Figure)
        fig2 = pipeline.plot_cooccurrence_heatmap()
        assert isinstance(fig2, matplotlib.figure.Figure)
