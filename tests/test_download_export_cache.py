"""
Tests for the export-data caching logic in UI/download.py.

The cache prevents all download buttons from disappearing when a user clicks
one button (which triggers a Streamlit rerun) by storing the already-computed
bytes in st.session_state and reusing them on the next run.
"""

import pytest


# ─────────────────────────────────────────────────────────────────────────────
# Helpers – simulate the cache logic without a live Streamlit session
# ─────────────────────────────────────────────────────────────────────────────

_EXPORT_CACHE_KEY = '_export_cache'
_EXPORT_CACHE_VER_KEY = '_export_cache_ver'


def _run_cache_logic(session, seq_count, all_motifs,
                     compute_fn, safe_fn="seq"):
    """
    Reproduce the caching block from UI/download.render() using a plain dict
    as a stand-in for st.session_state.

    compute_fn(all_motifs, seq_count, safe_fn) -> dict with keys:
        csv, excel, excel_label, excel_fname, json, bed, pdf,
        excel_error, pdf_error, export_times
    """
    cache_ver = (seq_count, len(all_motifs))

    if session.get(_EXPORT_CACHE_VER_KEY) != cache_ver or _EXPORT_CACHE_KEY not in session:
        data = compute_fn(all_motifs, seq_count, safe_fn)
        session[_EXPORT_CACHE_KEY] = data
        session[_EXPORT_CACHE_VER_KEY] = cache_ver
    else:
        data = session[_EXPORT_CACHE_KEY]

    return data


# ─────────────────────────────────────────────────────────────────────────────
# Tests
# ─────────────────────────────────────────────────────────────────────────────

class TestExportCache:
    """Unit tests for the export-data caching strategy."""

    def _make_compute_fn(self):
        """Return a compute function that records how many times it is called."""
        calls = {'count': 0}

        def compute(all_motifs, seq_count, safe_fn):
            calls['count'] += 1
            return {
                'csv': b'csv_data' if all_motifs else None,
                'excel': b'excel_data' if all_motifs else None,
                'excel_label': '📊 Excel',
                'excel_fname': f'{safe_fn}_results.xlsx',
                'json': b'json_data' if all_motifs else None,
                'bed': b'bed_data' if all_motifs else None,
                'pdf': b'pdf_data' if all_motifs else None,
                'excel_error': None,
                'pdf_error': None,
                'export_times': {'csv': 0.01},
            }

        return compute, calls

    def test_cache_populated_on_first_run(self):
        """compute_fn is called once on the very first render."""
        session = {}
        motifs = [{'Class': 'G-Quadruplex'}]
        compute, calls = self._make_compute_fn()

        _run_cache_logic(session, seq_count=1, all_motifs=motifs, compute_fn=compute)

        assert calls['count'] == 1
        assert _EXPORT_CACHE_KEY in session
        assert session[_EXPORT_CACHE_VER_KEY] == (1, 1)

    def test_cache_reused_on_rerun(self):
        """compute_fn is NOT called again when seq_count and motif count are unchanged."""
        session = {}
        motifs = [{'Class': 'G-Quadruplex'}]
        compute, calls = self._make_compute_fn()

        # First run (cache miss)
        _run_cache_logic(session, seq_count=1, all_motifs=motifs, compute_fn=compute)
        # Second run (simulates rerun after clicking a download button)
        result = _run_cache_logic(session, seq_count=1, all_motifs=motifs, compute_fn=compute)

        assert calls['count'] == 1, "compute_fn should only be called once"
        assert result['csv'] == b'csv_data'

    def test_all_format_data_available_after_rerun(self):
        """All five export formats are available from the cache on a rerun."""
        session = {}
        motifs = [{'Class': 'G-Quadruplex'}, {'Class': 'Z-DNA'}]
        compute, _ = self._make_compute_fn()

        _run_cache_logic(session, seq_count=1, all_motifs=motifs, compute_fn=compute)
        result = _run_cache_logic(session, seq_count=1, all_motifs=motifs, compute_fn=compute)

        assert result['csv'] is not None
        assert result['excel'] is not None
        assert result['json'] is not None
        assert result['bed'] is not None
        assert result['pdf'] is not None

    def test_cache_invalidated_when_motif_count_changes(self):
        """Cache is recomputed when a new analysis adds or removes motifs."""
        session = {}
        motifs_v1 = [{'Class': 'G-Quadruplex'}]
        motifs_v2 = [{'Class': 'G-Quadruplex'}, {'Class': 'Z-DNA'}]
        compute, calls = self._make_compute_fn()

        _run_cache_logic(session, seq_count=1, all_motifs=motifs_v1, compute_fn=compute)
        _run_cache_logic(session, seq_count=1, all_motifs=motifs_v2, compute_fn=compute)

        assert calls['count'] == 2, "Cache should be invalidated when motif count changes"
        assert session[_EXPORT_CACHE_VER_KEY] == (1, 2)

    def test_cache_invalidated_when_seq_count_changes(self):
        """Cache is recomputed when the number of sequences changes."""
        session = {}
        motifs = [{'Class': 'G-Quadruplex'}]
        compute, calls = self._make_compute_fn()

        _run_cache_logic(session, seq_count=1, all_motifs=motifs, compute_fn=compute)
        _run_cache_logic(session, seq_count=2, all_motifs=motifs, compute_fn=compute)

        assert calls['count'] == 2

    def test_empty_motifs_yields_none_data(self):
        """With no motifs, all export data is None and buttons should not appear."""
        session = {}
        compute, calls = self._make_compute_fn()

        result = _run_cache_logic(session, seq_count=0, all_motifs=[], compute_fn=compute)

        assert calls['count'] == 1
        assert result['csv'] is None
        assert result['excel'] is None
        assert result['json'] is None
        assert result['bed'] is None
        assert result['pdf'] is None

    def test_partial_failure_does_not_lose_successful_formats(self):
        """If PDF generation fails, CSV/Excel/JSON/BED remain available in the cache."""
        session = {}

        def compute_with_pdf_failure(all_motifs, seq_count, safe_fn):
            return {
                'csv': b'csv_data',
                'excel': b'excel_data',
                'excel_label': '📊 Excel',
                'excel_fname': f'{safe_fn}_results.xlsx',
                'json': b'json_data',
                'bed': b'bed_data',
                'pdf': None,           # PDF generation failed
                'excel_error': None,
                'pdf_error': 'PDF rendering error',
                'export_times': {},
            }

        motifs = [{'Class': 'G-Quadruplex'}]
        _run_cache_logic(session, seq_count=1, all_motifs=motifs,
                         compute_fn=compute_with_pdf_failure)

        # Simulate rerun (e.g. user clicked CSV button)
        result = _run_cache_logic(session, seq_count=1, all_motifs=motifs,
                                  compute_fn=compute_with_pdf_failure)

        assert result['csv'] is not None
        assert result['excel'] is not None
        assert result['json'] is not None
        assert result['bed'] is not None
        assert result['pdf'] is None
        assert result['pdf_error'] == 'PDF rendering error'
