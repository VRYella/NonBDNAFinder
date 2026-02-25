"""Unit tests for the stabilized has_results() in UI/storage_helpers.py.

Validates that has_results() returns True only when analysis_done is True
AND "results" exists in session_state, preventing false negatives during
Streamlit reruns triggered by download button clicks.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import unittest.mock as mock

# Patch streamlit before importing the module under test
_st_mock = mock.MagicMock()
_session_state = {}

def _get(key, default=None):
    return _session_state.get(key, default)

def _contains(key):
    return key in _session_state

_st_mock.session_state = mock.MagicMock()
_st_mock.session_state.get = _get
_st_mock.session_state.__contains__ = lambda self, key: _contains(key)

sys.modules['streamlit'] = _st_mock

import importlib
import UI.storage_helpers as sh_mod
importlib.reload(sh_mod)

_has_results = sh_mod.has_results


def setup_session(**kwargs):
    """Reset session state to given kwargs."""
    _session_state.clear()
    _session_state.update(kwargs)


class TestHasResultsStable:
    """has_results() must be stable across reruns."""

    def test_false_when_nothing_set(self):
        setup_session()
        assert _has_results() is False

    def test_false_when_only_analysis_done(self):
        setup_session(analysis_done=True)
        assert _has_results() is False

    def test_false_when_only_results_present(self):
        setup_session(results=[[{"Class": "Z-DNA"}]])
        assert _has_results() is False

    def test_true_when_both_set(self):
        setup_session(analysis_done=True, results=[[{"Class": "Z-DNA"}]])
        assert _has_results() is True

    def test_false_when_analysis_done_is_false(self):
        setup_session(analysis_done=False, results=[[{"Class": "Z-DNA"}]])
        assert _has_results() is False

    def test_true_with_empty_results_list_but_analysis_done(self):
        """analysis_done=True + results=[] is still True (key exists)."""
        setup_session(analysis_done=True, results=[])
        assert _has_results() is True

    def test_no_disk_storage_dependency(self):
        """Result availability must not depend on use_disk_storage."""
        setup_session(analysis_done=True, results=[[]], use_disk_storage=True, seq_ids=["s1"])
        assert _has_results() is True

    def test_no_results_storage_dependency(self):
        """Result availability must not depend on results_storage."""
        setup_session(analysis_done=True, results=[[]], results_storage={})
        assert _has_results() is True


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
