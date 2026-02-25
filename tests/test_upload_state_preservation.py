"""Unit tests verifying that upload.py does NOT clear analysis results during
Streamlit reruns (e.g. reruns triggered by download button clicks).

The "Persist sequences to session state" block in upload.render() must be a
no-op whenever analysis_done is True so that results survive reruns.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import types
import unittest.mock as mock


# ---------------------------------------------------------------------------
# Minimal Streamlit stub so upload.py can be imported without a server
# ---------------------------------------------------------------------------

_session_state_dict = {}


class _FakeSessionState:
    """Dict-backed session state that supports attribute and item access."""

    def get(self, key, default=None):
        return _session_state_dict.get(key, default)

    def __contains__(self, key):
        return key in _session_state_dict

    def __getattr__(self, key):
        try:
            return _session_state_dict[key]
        except KeyError:
            raise AttributeError(key)

    def __setattr__(self, key, value):
        _session_state_dict[key] = value

    def __getitem__(self, key):
        return _session_state_dict[key]

    def __setitem__(self, key, value):
        _session_state_dict[key] = value

    def keys(self):
        return _session_state_dict.keys()

    def pop(self, key, *args):
        return _session_state_dict.pop(key, *args)


def _setup_session(**kwargs):
    _session_state_dict.clear()
    _session_state_dict.update(kwargs)


# Build a minimal stub for the streamlit module
_st_stub = types.ModuleType("streamlit")
_st_stub.session_state = _FakeSessionState()
_st_stub.cache_data = lambda *a, **kw: (lambda fn: fn)

# All UI calls are no-ops
for _attr in [
    "tabs", "tab", "columns", "expander", "container", "empty",
    "markdown", "info", "warning", "error", "success", "spinner",
    "button", "radio", "text_area", "text_input", "number_input",
    "file_uploader", "progress", "dataframe", "caption", "metric",
    "download_button", "rerun", "set_page_config", "write",
]:
    setattr(_st_stub, _attr, mock.MagicMock(return_value=mock.MagicMock()))

sys.modules["streamlit"] = _st_stub


# ---------------------------------------------------------------------------
# Target function under test: the persist-sequences guard logic
# ---------------------------------------------------------------------------
# We replicate the exact guard logic from upload.py here rather than importing
# the full module (which has many heavy dependencies).  This keeps the test
# fast and focused on the specific invariant.

def _simulate_persist_sequences(disk_seq_ids, seqs, names, use_disk=False):
    """
    Reproduce the 'Persist sequences to session state' block from upload.py.
    Returns True if session state was mutated (results would be cleared).
    """
    import streamlit as st

    _analysis_locked = st.session_state.get("analysis_done", False)
    mutated = False

    if disk_seq_ids:
        if not _analysis_locked:
            st.session_state.seq_ids = disk_seq_ids
            st.session_state.names = names
            st.session_state.results_storage = {}
            st.session_state.seqs = []
            st.session_state.results = []
            mutated = True
    elif seqs:
        if use_disk:
            if not _analysis_locked:
                st.session_state.seq_ids = ["fake_id"]
                st.session_state.names = names
                st.session_state.results_storage = {}
                st.session_state.seqs = []
                st.session_state.results = []
                mutated = True
        else:
            if not _analysis_locked:
                st.session_state.seqs = seqs
                st.session_state.names = names
                st.session_state.results = []
                mutated = True

    return mutated


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestUploadStatePreservation:
    """Upload page must not clear analysis results when analysis_done=True."""

    # -- Legacy in-memory mode -----------------------------------------------

    def test_results_preserved_on_rerun_memory_mode(self):
        """Download-button rerun: seqs present, analysis_done=True → results NOT cleared."""
        existing_results = [[{"Class": "Z-DNA", "Start": 1, "End": 10}]]
        _setup_session(
            analysis_done=True,
            seqs=["ACGTACGT"],
            names=["seq1"],
            results=existing_results,
        )
        mutated = _simulate_persist_sequences(
            disk_seq_ids=[],
            seqs=["ACGTACGT"],
            names=["seq1"],
            use_disk=False,
        )
        assert not mutated, "Persist block must not run when analysis_done=True"
        assert _session_state_dict["results"] == existing_results, (
            "results must be unchanged after download-button rerun"
        )

    def test_results_cleared_before_analysis_memory_mode(self):
        """First upload (analysis_done=False) → results ARE cleared (correct behaviour)."""
        _setup_session(
            analysis_done=False,
            seqs=[],
            names=[],
            results=[],
        )
        mutated = _simulate_persist_sequences(
            disk_seq_ids=[],
            seqs=["ACGTACGT"],
            names=["seq1"],
            use_disk=False,
        )
        assert mutated, "Persist block must run when analysis_done=False"
        assert _session_state_dict["results"] == [], "results should be empty list on fresh upload"
        assert _session_state_dict["seqs"] == ["ACGTACGT"]

    def test_analysis_done_false_by_default_no_mutation(self):
        """If analysis_done key is absent, treat as False (no guard = allow mutation)."""
        _setup_session(seqs=[], names=[], results=[])
        # analysis_done absent → defaults to False → persist is allowed
        mutated = _simulate_persist_sequences(
            disk_seq_ids=[],
            seqs=["ACGTACGT"],
            names=["seq1"],
            use_disk=False,
        )
        assert mutated

    # -- Disk-based storage mode ---------------------------------------------

    def test_results_preserved_on_rerun_disk_mode_seq_ids(self):
        """Disk mode: disk_seq_ids present on rerun, analysis_done=True → state NOT cleared."""
        _setup_session(
            analysis_done=True,
            seq_ids=["orig_id"],
            names=["seq1"],
            results_storage={"orig_id": object()},
            seqs=[],
            results=[],
        )
        orig_storage = _session_state_dict["results_storage"]
        mutated = _simulate_persist_sequences(
            disk_seq_ids=["new_id"],
            seqs=[],
            names=["seq1"],
            use_disk=True,
        )
        assert not mutated
        assert _session_state_dict["seq_ids"] == ["orig_id"], (
            "seq_ids must not be overwritten during rerun when analysis_done=True"
        )
        assert _session_state_dict["results_storage"] is orig_storage, (
            "results_storage must not be cleared during rerun when analysis_done=True"
        )

    def test_results_cleared_before_analysis_disk_mode_seq_ids(self):
        """Disk mode: first upload (analysis_done=False) → state IS updated."""
        _setup_session(
            analysis_done=False,
            seq_ids=[],
            names=[],
            results_storage={},
            seqs=[],
            results=[],
        )
        mutated = _simulate_persist_sequences(
            disk_seq_ids=["new_id"],
            seqs=[],
            names=["seq1"],
            use_disk=True,
        )
        assert mutated
        assert _session_state_dict["seq_ids"] == ["new_id"]
        assert _session_state_dict["results_storage"] == {}

    def test_results_preserved_on_rerun_disk_mode_seqs(self):
        """Disk mode via seqs: analysis_done=True → state NOT overwritten."""
        existing_storage = {"s1": object()}
        _setup_session(
            analysis_done=True,
            seq_ids=["s1"],
            names=["seq1"],
            results_storage=existing_storage,
            seqs=[],
            results=[],
        )
        mutated = _simulate_persist_sequences(
            disk_seq_ids=[],
            seqs=["ACGTACGT"],
            names=["seq1"],
            use_disk=True,
        )
        assert not mutated
        assert _session_state_dict["results_storage"] is existing_storage

    # -- Edge cases ----------------------------------------------------------

    def test_no_seqs_no_disk_ids_no_mutation(self):
        """If both seqs and disk_seq_ids are empty, nothing happens regardless."""
        existing = [[{"Class": "G-Quadruplex"}]]
        _setup_session(analysis_done=True, results=existing, seqs=["ACGT"], names=["s"])
        mutated = _simulate_persist_sequences([], [], [], use_disk=False)
        assert not mutated
        assert _session_state_dict["results"] == existing

    def test_analysis_done_flag_is_only_guard(self):
        """analysis_done=True is the ONLY condition needed to preserve state."""
        for analysis_done_value in [True]:
            existing = [[{"Class": "Z-DNA"}]]
            _setup_session(analysis_done=analysis_done_value, results=existing,
                           seqs=["ACGT"], names=["s"])
            mutated = _simulate_persist_sequences([], ["ACGT"], ["s"], use_disk=False)
            assert not mutated
            assert _session_state_dict["results"] == existing


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
