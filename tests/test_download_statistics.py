"""Unit tests for generate_statistics in UI/download.py.

Validates that dist_df and sub_df are always pandas DataFrames (never undefined)
regardless of whether all_motifs is empty or populated, ensuring a stable
Streamlit widget tree on every rerun.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import pandas as pd

# Patch streamlit with a pass-through cache_data so the decorator returns the
# original function unchanged, allowing direct testing without a Streamlit server.
import unittest.mock as mock

_st_mock = mock.MagicMock()
# Make cache_data a pass-through decorator
_st_mock.cache_data = lambda *a, **kw: (lambda fn: fn)
sys.modules['streamlit'] = _st_mock

# Import (or reload) after patching so the decorator is the pass-through
import importlib
import UI.download as dl_mod
importlib.reload(dl_mod)

_generate_statistics = dl_mod.generate_statistics


class TestGenerateStatisticsAlwaysDataFrame:
    """dist_df and sub_df must always be DataFrames, never undefined."""

    def test_empty_motifs_returns_dataframes(self):
        dist_df, sub_df, stats_excel = _generate_statistics((), (), (), 0)
        assert isinstance(dist_df, pd.DataFrame)
        assert isinstance(sub_df, pd.DataFrame)
        assert dist_df.empty
        assert sub_df.empty
        assert stats_excel == b""

    def test_populated_motifs_returns_dataframes(self):
        motifs = (
            {'Sequence_Name': 'seq1', 'Class': 'G-Quadruplex', 'Subclass': 'Canonical', 'Start': 1, 'End': 10, 'Length': 9},
            {'Sequence_Name': 'seq1', 'Class': 'Z-DNA', 'Subclass': 'Classic', 'Start': 20, 'End': 30, 'Length': 10},
        )
        names = ('seq1',)
        lengths = (1000,)
        dist_df, sub_df, stats_excel = _generate_statistics(motifs, names, lengths, 1)
        assert isinstance(dist_df, pd.DataFrame)
        assert isinstance(sub_df, pd.DataFrame)
        assert not dist_df.empty
        assert not sub_df.empty
        assert isinstance(stats_excel, bytes)
        assert len(stats_excel) > 0

    def test_dist_df_columns_present(self):
        motifs = (
            {'Sequence_Name': 'seq1', 'Class': 'G-Quadruplex', 'Subclass': 'Canonical', 'Start': 1, 'End': 5, 'Length': 4},
        )
        dist_df, _, _ = _generate_statistics(motifs, ('seq1',), (500,), 1)
        expected = {'Sequence Name', 'Motif Class', 'Count', 'Coverage (%)', 'Motifs per kbp', 'Average Length (bp)', 'Covered Bases (bp)'}
        assert expected.issubset(set(dist_df.columns))

    def test_sub_df_columns_present(self):
        motifs = (
            {'Sequence_Name': 'seq1', 'Class': 'G-Quadruplex', 'Subclass': 'Canonical', 'Start': 1, 'End': 5, 'Length': 4},
        )
        _, sub_df, _ = _generate_statistics(motifs, ('seq1',), (500,), 1)
        expected = {'Sequence Name', 'Motif Subclass', 'Count', 'Coverage (%)'}
        assert expected.issubset(set(sub_df.columns))

    def test_no_early_return_pattern(self):
        """Ensure dist_df and sub_df are always DataFrames (not undefined) even with empty input."""
        for motifs in [(), tuple()]:
            dist_df, sub_df, stats_excel = _generate_statistics(motifs, (), (), 0)
            assert isinstance(dist_df, pd.DataFrame), "dist_df must always be a DataFrame"
            assert isinstance(sub_df, pd.DataFrame), "sub_df must always be a DataFrame"

    def test_stats_excel_empty_when_no_data(self):
        dist_df, sub_df, stats_excel = _generate_statistics((), (), (), 0)
        assert stats_excel == b""

    def test_stats_excel_populated_when_data(self):
        motifs = (
            {'Sequence_Name': 'chr1', 'Class': 'Z-DNA', 'Subclass': 'Classic', 'Start': 10, 'End': 20, 'Length': 10},
        )
        _, _, stats_excel = _generate_statistics(motifs, ('chr1',), (2000,), 1)
        assert isinstance(stats_excel, bytes)
        assert len(stats_excel) > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
