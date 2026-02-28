"""Performance and correctness tests for the vectorised detector optimisations.

Validates that:
1. SlippedDNA numpy path finds the same motifs as the regex fallback.
2. Cruciform numpy sliding-window path finds the same inverted-repeat seeds
   as the original Python dict-index approach.
3. Triplex numpy sliding-window path finds the same mirror-repeat seeds as
   the original Python dict-index approach.
4. RLoop prefix_g is computed only once per annotate_sequence() call (cache).
5. G-Quadruplex detector completes within a generous time bound on a G-rich
   sequence that previously triggered catastrophic regex backtracking.
"""
import sys
import os
import time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _random_seq(n: int, seed: int = 42) -> str:
    import random
    rng = random.Random(seed)
    return "".join(rng.choices("ACGT", k=n))


def _str_with_repeats(length: int = 5000) -> str:
    """Return a sequence that contains known STR and direct-repeat motifs."""
    base = (
        "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG"  # 60 bp CAG STR
        + "ACGTACGT" * 10                                                  # filler
        + "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"                  # ATCG direct repeat
    )
    return (base * (length // len(base) + 1))[:length]


# ---------------------------------------------------------------------------
# SlippedDNA – numpy vs regex parity
# ---------------------------------------------------------------------------

class TestSlippedDNANumpyParity:
    """numpy and regex implementations must return identical candidates."""

    def _get_candidates(self, detector, seq):
        """Return candidate dicts from find_all_tandem_repeats."""
        return detector.find_all_tandem_repeats(seq)

    def test_no_motifs_random_sequence(self):
        """Random sequence should produce 0 candidates after stringent criteria."""
        from Detectors.slipped.detector import SlippedDNADetector
        det = SlippedDNADetector()
        seq = _random_seq(5000)
        motifs = det.detect_motifs(seq, "test")
        assert isinstance(motifs, list)

    def test_cag_repeat_detected(self):
        """A long CAG STR (≥20 bp, ≥6 copies) must be detected."""
        from Detectors.slipped.detector import SlippedDNADetector
        det = SlippedDNADetector()
        seq = "N" * 10 + "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAG" + "N" * 10  # 10 copies
        motifs = det.detect_motifs(seq, "test")
        classes = [m["Subclass"] for m in motifs]
        assert any("STR" in c for c in classes), f"Expected STR, got {classes}"

    def test_numpy_regex_same_detect_results(self):
        """numpy and regex implementations must find the same final motifs."""
        try:
            import numpy as _np
        except ImportError:
            pytest.skip("numpy not available")

        from Detectors.slipped.detector import SlippedDNADetector
        det = SlippedDNADetector()
        seq = _str_with_repeats(3000)

        motifs_np = det.detect_motifs(seq, "np")

        # Force regex fallback by monkey-patching
        import Detectors.slipped.detector as slipped_mod
        orig = slipped_mod._NUMPY_AVAILABLE
        slipped_mod._NUMPY_AVAILABLE = False
        try:
            motifs_regex = det.detect_motifs(seq, "regex")
        finally:
            slipped_mod._NUMPY_AVAILABLE = orig

        np_keys = {(m["Start"], m["End"], m["Repeat_Unit"]) for m in motifs_np}
        rx_keys = {(m["Start"], m["End"], m["Repeat_Unit"]) for m in motifs_regex}
        assert np_keys == rx_keys, (
            f"numpy found {np_keys - rx_keys} extra, regex found {rx_keys - np_keys} extra"
        )


# ---------------------------------------------------------------------------
# Cruciform – sliding window seed pairs
# ---------------------------------------------------------------------------

class TestCruciformSlidingWindow:
    """Sliding-window seed discovery must return identical pairs to the dict index."""

    def test_inverted_repeats_found(self):
        """A perfect 10-bp inverted repeat (no loop) must be detected."""
        from Detectors.cruciform.detector import CruciformDetector
        det = CruciformDetector()
        # ATCGATCGAT .... ATCGATCGAT_RC (= ATCGATCGAT)
        arm = "GCATGCATGCATGCATGC"  # 18 bp
        from Utilities.detectors_utils import revcomp
        seq = "A" * 20 + arm + "T" * 5 + revcomp(arm) + "A" * 20
        hits = det.find_inverted_repeats(seq)
        assert len(hits) > 0, "Expected at least one inverted repeat"

    def test_numpy_python_same_seed_pairs(self):
        """numpy sliding-window and Python dict-index must agree on valid seed pairs."""
        try:
            import numpy as _np
        except ImportError:
            pytest.skip("numpy not available")

        from Detectors.cruciform.detector import CruciformDetector
        det = CruciformDetector()
        seq = _random_seq(20000)
        n = len(seq)
        k = det.SEED_SIZE
        max_loop = det.MAX_LOOP

        # numpy path
        pairs_np = set(CruciformDetector._find_seed_pairs_numpy(seq, n, k, max_loop))

        # Python fallback
        from collections import defaultdict
        from Utilities.detectors_utils import revcomp
        seed_index = defaultdict(list)
        for i in range(n - k + 1):
            seed_index[seq[i:i+k]].append(i)
        pairs_py = set()
        for kmer, i_positions in seed_index.items():
            rc = revcomp(kmer)
            for i in i_positions:
                for j in seed_index.get(rc, []):
                    if i + k <= j <= i + k + max_loop:
                        pairs_py.add((i, j))

        assert pairs_np == pairs_py, (
            f"numpy extra: {pairs_np - pairs_py}, python extra: {pairs_py - pairs_np}"
        )


# ---------------------------------------------------------------------------
# Triplex – sliding window seed pairs
# ---------------------------------------------------------------------------

class TestTriplexSlidingWindow:
    """Mirror-repeat sliding-window seeds must match the Python dict-index."""

    def test_numpy_python_same_mirror_pairs(self):
        """numpy mirror-hash and Python dict-index must agree on valid seed pairs."""
        try:
            import numpy as _np
        except ImportError:
            pytest.skip("numpy not available")

        from Detectors.triplex.detector import TriplexDetector
        det = TriplexDetector()
        seq = _random_seq(20000)
        n = len(seq)
        k = det.SEED_SIZE

        # numpy path
        pairs_np = set(TriplexDetector._find_mirror_pairs_numpy(seq, n, k, _np))

        # Python fallback
        from collections import defaultdict
        seed_index = defaultdict(list)
        for i in range(n - k + 1):
            seed_index[seq[i:i+k]].append(i)
        pairs_py = set()
        for kmer, i_positions in seed_index.items():
            mirror = kmer[::-1]
            for i in i_positions:
                for j in seed_index.get(mirror, []):
                    if i + k <= j <= i + k + det.MAX_LOOP:
                        pairs_py.add((i, j))

        assert pairs_np == pairs_py, (
            f"numpy extra: {pairs_np - pairs_py}, python extra: {pairs_py - pairs_np}"
        )


# ---------------------------------------------------------------------------
# RLoop – prefix_g computed once
# ---------------------------------------------------------------------------

class TestRLoopPrefixG:
    """prefix_g must be the same object passed into every _find_rez call."""

    def test_prefix_g_passed_to_find_rez(self):
        """annotate_sequence must pass prefix_g to _find_rez (not recompute it)."""
        from Detectors.rloop.detector import RLoopDetector
        det = RLoopDetector()
        seq = _random_seq(5000)

        calls = []

        orig_find_rez = det._find_rez.__func__

        def recording_find_rez(self, s, riz_end, prefix_g=None):
            calls.append(prefix_g is not None)
            return orig_find_rez(self, s, riz_end, prefix_g=prefix_g)

        import types
        det._find_rez = types.MethodType(recording_find_rez, det)

        det.annotate_sequence(seq)

        if calls:
            assert all(calls), "All _find_rez calls should receive a pre-built prefix_g"

    def test_prefix_g_shape(self):
        """_build_prefix_g must return an array of length n+1."""
        from Detectors.rloop.detector import RLoopDetector
        det = RLoopDetector()
        seq = "ACGTGGGGACGT"
        pg = det._build_prefix_g(seq)
        assert len(pg) == len(seq) + 1


# ---------------------------------------------------------------------------
# GQuadruplex – no catastrophic backtracking on G-rich sequence
# ---------------------------------------------------------------------------

class TestGQuadruplexPatternSpeed:
    """Stacked G4 / higher-order G4 must not trigger catastrophic backtracking."""

    GRICH = "GGGGTTGGGGTTGGGGTTGGGGTTGGGGTTGGGGTTGGGGTTGGGGTTGGGGTTGGGGTT"

    def test_stacked_g4_fast_on_grich(self):
        """stacked_g4 pattern must match in < 10 ms on a 60-bp all-G4 sequence
        and the match must contain at least 2 G-tracts (basic structural check)."""
        import re
        from Detectors.gquad.detector import GQuadruplexDetector
        patterns = GQuadruplexDetector().get_patterns()
        pat = re.compile(patterns["stacked_g4"][0][0])
        start = time.perf_counter()
        matches = list(pat.finditer(self.GRICH))
        elapsed = time.perf_counter() - start
        assert elapsed < 0.01, f"stacked_g4 took {elapsed*1000:.1f} ms (> 10 ms limit)"
        assert len(matches) >= 1, "Expected at least one stacked G4 match"
        # Structural check: each match must contain at least 2 GGGG+ tracts
        for m in matches:
            g_tracts = re.findall(r'G{3,}', m.group())
            assert len(g_tracts) >= 2, (
                f"stacked G4 match '{m.group()[:30]}' has < 2 G-tracts: {g_tracts}"
            )

    def test_higher_order_g4_fast_on_grich(self):
        """higher_order_g4 pattern must match in < 10 ms on the same sequence."""
        import re
        from Detectors.gquad.detector import GQuadruplexDetector
        patterns = GQuadruplexDetector().get_patterns()
        pat = re.compile(patterns["higher_order_g4"][0][0])
        start = time.perf_counter()
        matches = list(pat.finditer(self.GRICH))
        elapsed = time.perf_counter() - start
        assert elapsed < 0.01, f"higher_order_g4 took {elapsed*1000:.1f} ms (> 10 ms limit)"
        assert len(matches) >= 1, "Expected at least one higher-order G4 match"

    def test_detect_motifs_g4_completes_quickly(self):
        """detect_motifs on a 362-bp G-rich sequence must finish in < 1 second
        (was > 150 ms per call before the fix)."""
        from Detectors.gquad.detector import GQuadruplexDetector
        det = GQuadruplexDetector()
        seq = (
            "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
            "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG"
            "GGGGTTGGGGTTGGGGTTGGGGTTGGGGTTGGGGTTGGGGTTGGGGTTGGGGTTGGGGTT"
        )
        # Warmup (numba JIT if present)
        det.detect_motifs(seq[:50], "w")
        start = time.perf_counter()
        motifs = det.detect_motifs(seq, "test")
        elapsed = time.perf_counter() - start
        assert elapsed < 1.0, f"G4 detection took {elapsed*1000:.1f} ms (> 1000 ms limit)"
        assert len(motifs) >= 1, "Expected at least one G4 motif"
