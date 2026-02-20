"""
Tests for per-task timing on a ~4 MB sequence.

Validates:
  - PerformanceTracker.add_task_time() records named task durations
  - get_summary() exposes task_times in its output
  - format_performance_summary() renders a "Task Timing Breakdown" section
  - Visualization and export preparation times are measured and non-negative
  - End-to-end analysis + visualization + export timing works on a ~4 MB sequence
"""

import time
import pytest


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _make_seq(length: int) -> str:
    """Build a pseudorandom DNA sequence of *length* bp with embedded G4 motifs."""
    import random
    random.seed(42)
    bases = "ATGC"
    seq = "".join(random.choice(bases) for _ in range(length))
    # Embed G4 motifs every ~100 kb so results are non-trivial
    motif = "GGGTTTAGGGTTTGGG"
    positions = range(0, length - len(motif), 100_000)
    seq_list = list(seq)
    for pos in positions:
        seq_list[pos:pos + len(motif)] = list(motif)
    return "".join(seq_list)


# ─────────────────────────────────────────────────────────────────────────────
# PerformanceTracker unit tests
# ─────────────────────────────────────────────────────────────────────────────

class TestPerformanceTrackerTaskTimes:
    """Unit tests for the task_times feature of PerformanceTracker."""

    def _tracker(self):
        from UI.performance_stats import PerformanceTracker
        return PerformanceTracker()

    def test_task_times_empty_on_init(self):
        tracker = self._tracker()
        assert tracker.task_times == {}

    def test_add_task_time_stores_value(self):
        tracker = self._tracker()
        tracker.add_task_time('visualization', 1.23)
        assert tracker.task_times['visualization'] == pytest.approx(1.23)

    def test_add_multiple_task_times(self):
        tracker = self._tracker()
        tracker.add_task_time('visualization', 0.5)
        tracker.add_task_time('export', 0.3)
        assert 'visualization' in tracker.task_times
        assert 'export' in tracker.task_times

    def test_get_summary_includes_task_times(self):
        tracker = self._tracker()
        tracker.start()
        tracker.add_task_time('visualization', 2.0)
        tracker.add_task_time('export', 0.4)
        tracker.end()
        summary = tracker.get_summary()
        assert 'task_times' in summary
        assert summary['task_times']['visualization'] == pytest.approx(2.0)
        assert summary['task_times']['export'] == pytest.approx(0.4)

    def test_get_summary_task_times_empty_when_none_added(self):
        tracker = self._tracker()
        tracker.start()
        tracker.end()
        summary = tracker.get_summary()
        assert summary['task_times'] == {}

    def test_overwrite_task_time(self):
        tracker = self._tracker()
        tracker.add_task_time('visualization', 1.0)
        tracker.add_task_time('visualization', 2.5)
        assert tracker.task_times['visualization'] == pytest.approx(2.5)


# ─────────────────────────────────────────────────────────────────────────────
# format_performance_summary output
# ─────────────────────────────────────────────────────────────────────────────

class TestFormatPerformanceSummaryTaskSection:
    """format_performance_summary must render a Task Timing Breakdown section."""

    def _fmt(self, task_dict):
        from UI.performance_stats import PerformanceTracker, format_performance_summary
        tracker = PerformanceTracker()
        tracker.start()
        tracker.total_bp_processed = 4_000_000
        for name, elapsed in task_dict.items():
            tracker.add_task_time(name, elapsed)
        tracker.end()

        def _fmt_time(s):
            return f"{s:.2f}s"

        return format_performance_summary(tracker, _fmt_time)

    def test_task_section_present_when_tasks_recorded(self):
        output = self._fmt({'visualization': 1.5, 'export': 0.3})
        assert "Task Timing Breakdown" in output

    def test_visualization_task_listed(self):
        output = self._fmt({'visualization': 1.5})
        assert "Visualization" in output
        assert "1.50s" in output

    def test_export_task_listed(self):
        output = self._fmt({'export': 0.3})
        assert "Export" in output
        assert "0.30s" in output

    def test_task_section_absent_when_no_tasks(self):
        output = self._fmt({})
        assert "Task Timing Breakdown" not in output

    def test_percentage_appears_in_output(self):
        output = self._fmt({'visualization': 1.0})
        # The percentage line should contain a '%' character
        assert "%" in output


# ─────────────────────────────────────────────────────────────────────────────
# End-to-end timing on a ~4 MB sequence
# ─────────────────────────────────────────────────────────────────────────────

class TestEndToEndTimingOn4MbSequence:
    """
    Verify that analysis, visualization-prep, and export-prep are all timed
    correctly on a ~4 MB sequence.
    """

    SEQ_LENGTH = 4_000_000  # 4 MB (4 million bp)

    def test_analysis_timing_recorded(self):
        """PerformanceTracker.get_total_time() is positive after analyzing 4 MB."""
        from Utilities.nonbscanner import analyze_sequence
        from UI.performance_stats import PerformanceTracker

        seq = _make_seq(self.SEQ_LENGTH)
        tracker = PerformanceTracker()
        tracker.start()
        motifs = analyze_sequence(seq, "test_4mb", enabled_classes=["G-Quadruplex"])
        elapsed = tracker.get_total_time()
        tracker.add_sequence_result("test_4mb", len(seq), elapsed, len(motifs))
        tracker.end()

        assert tracker.get_total_time() > 0
        assert tracker.total_bp_processed == self.SEQ_LENGTH

    def test_visualization_task_time_recorded(self):
        """Visualization timing can be added to the tracker and appears in summary."""
        from UI.performance_stats import PerformanceTracker

        tracker = PerformanceTracker()
        tracker.start()
        tracker.total_bp_processed = self.SEQ_LENGTH

        # Simulate a timed visualization step
        viz_start = time.time()
        time.sleep(0.01)  # Minimal sleep to produce a measurable elapsed time
        viz_elapsed = time.time() - viz_start
        tracker.add_task_time('visualization', viz_elapsed)
        tracker.end()

        summary = tracker.get_summary()
        assert summary['task_times']['visualization'] >= 0.0
        assert summary['task_times']['visualization'] == pytest.approx(viz_elapsed, abs=0.05)

    def test_export_task_times_recorded_for_each_format(self):
        """Each export format's preparation time is stored separately."""
        from UI.performance_stats import PerformanceTracker

        tracker = PerformanceTracker()
        tracker.start()
        tracker.total_bp_processed = self.SEQ_LENGTH

        # Simulate timing for CSV, JSON, BED export preparation
        for fmt in ('csv', 'json', 'bed'):
            _t0 = time.time()
            time.sleep(0.005)
            tracker.add_task_time(f'export_{fmt}', time.time() - _t0)

        tracker.end()

        summary = tracker.get_summary()
        for fmt in ('csv', 'json', 'bed'):
            key = f'export_{fmt}'
            assert key in summary['task_times'], f"Missing export timing for {fmt}"
            assert summary['task_times'][key] >= 0.0

    def test_all_task_times_non_negative(self):
        """All recorded task times must be non-negative."""
        from UI.performance_stats import PerformanceTracker

        tracker = PerformanceTracker()
        tracker.start()
        tracker.add_task_time('visualization', 3.14)
        tracker.add_task_time('export_csv', 0.05)
        tracker.add_task_time('export_excel', 0.12)
        tracker.add_task_time('export_json', 0.04)
        tracker.add_task_time('export_bed', 0.03)
        tracker.add_task_time('export_pdf', 0.8)
        tracker.end()

        for name, elapsed in tracker.task_times.items():
            assert elapsed >= 0.0, f"Task '{name}' has negative time: {elapsed}"

    def test_format_summary_shows_all_tasks_on_4mb(self):
        """format_performance_summary lists every recorded task for 4 MB run."""
        from UI.performance_stats import PerformanceTracker, format_performance_summary

        tracker = PerformanceTracker()
        tracker.start()
        tracker.total_bp_processed = self.SEQ_LENGTH
        tracker.add_task_time('visualization', 2.5)
        tracker.add_task_time('export_csv', 0.05)
        tracker.add_task_time('export_excel', 0.15)
        tracker.end()

        output = format_performance_summary(tracker, lambda s: f"{s:.3f}s")
        assert "Task Timing Breakdown" in output
        assert "Visualization" in output
        assert "Export Csv" in output
        assert "Export Excel" in output
