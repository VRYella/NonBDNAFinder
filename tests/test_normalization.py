"""Unit tests for length-aware, structure-bounded normalization (1–3 scale)."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
from Detectors.base.base_detector import STRUCTURAL_LENGTH_CAPS, DISEASE_CAP_OVERRIDES
from Detectors.gquad.detector import GQuadruplexDetector
from Detectors.aphilic.detector import APhilicDetector
from Detectors.cruciform.detector import CruciformDetector
from Detectors.curved.detector import CurvedDNADetector
from Detectors.zdna.detector import ZDNADetector
from Detectors.slipped.detector import SlippedDNADetector
from Detectors.triplex.detector import TriplexDetector
from Detectors.rloop.detector import RLoopDetector
from Detectors.imotif.detector import IMotifDetector


# ---------------------------------------------------------------------------
# Structural Length Caps
# ---------------------------------------------------------------------------

class TestStructuralLengthCaps:
    def test_all_required_keys_present(self):
        expected_keys = {
            "A-philic_DNA", "Cruciform", "Curved_DNA_local", "Curved_DNA_global",
            "G-Quadruplex", "Z-DNA", "Slipped_DNA_STR", "Slipped_DNA_Direct",
            "Triplex", "R-Loop", "i-Motif",
        }
        assert expected_keys.issubset(set(STRUCTURAL_LENGTH_CAPS.keys()))

    def test_cap_values(self):
        assert STRUCTURAL_LENGTH_CAPS["A-philic_DNA"] == 300
        assert STRUCTURAL_LENGTH_CAPS["Cruciform"] == 200
        assert STRUCTURAL_LENGTH_CAPS["Curved_DNA_local"] == 50
        assert STRUCTURAL_LENGTH_CAPS["Curved_DNA_global"] == 120
        assert STRUCTURAL_LENGTH_CAPS["G-Quadruplex"] == 120
        assert STRUCTURAL_LENGTH_CAPS["Z-DNA"] == 300
        assert STRUCTURAL_LENGTH_CAPS["Slipped_DNA_STR"] == 1000
        assert STRUCTURAL_LENGTH_CAPS["Slipped_DNA_Direct"] == 500
        assert STRUCTURAL_LENGTH_CAPS["Triplex"] == 150
        assert STRUCTURAL_LENGTH_CAPS["R-Loop"] == 2000
        assert STRUCTURAL_LENGTH_CAPS["i-Motif"] == 60

    def test_disease_cap_overrides_present(self):
        assert "Slipped_DNA" in DISEASE_CAP_OVERRIDES
        assert DISEASE_CAP_OVERRIDES["Slipped_DNA"]["Huntington"] == 1000
        assert DISEASE_CAP_OVERRIDES["Slipped_DNA"]["Myotonic_Dystrophy"] == 1500


# ---------------------------------------------------------------------------
# Detector Length Cap Methods
# ---------------------------------------------------------------------------

class TestDetectorLengthCaps:
    def test_aphilic_cap(self):
        d = APhilicDetector()
        assert d.get_length_cap() == 300

    def test_cruciform_cap(self):
        d = CruciformDetector()
        assert d.get_length_cap() == 200

    def test_curved_local_cap(self):
        d = CurvedDNADetector()
        assert d.get_length_cap("Local Curvature") == 50

    def test_curved_global_cap(self):
        d = CurvedDNADetector()
        assert d.get_length_cap("Global Curvature") == 120

    def test_curved_default_cap_is_global(self):
        d = CurvedDNADetector()
        assert d.get_length_cap() == 120

    def test_gquad_cap(self):
        d = GQuadruplexDetector()
        assert d.get_length_cap() == 120

    def test_zdna_cap(self):
        d = ZDNADetector()
        assert d.get_length_cap() == 300

    def test_slipped_str_cap(self):
        d = SlippedDNADetector()
        assert d.get_length_cap("STR") == 1000

    def test_slipped_direct_cap(self):
        d = SlippedDNADetector()
        assert d.get_length_cap("Direct Repeat") == 500

    def test_slipped_default_cap_is_direct(self):
        d = SlippedDNADetector()
        assert d.get_length_cap() == 500

    def test_triplex_cap(self):
        d = TriplexDetector()
        assert d.get_length_cap() == 150

    def test_rloop_cap(self):
        d = RLoopDetector()
        assert d.get_length_cap() == 2000

    def test_imotif_cap(self):
        d = IMotifDetector()
        assert d.get_length_cap() == 60


# ---------------------------------------------------------------------------
# normalize_score formula: Score ∈ [1, 3]
# ---------------------------------------------------------------------------

class TestNormalizeScoreFormula:
    """Validate the new length-aware normalization formula on GQuadruplexDetector."""

    def setup_method(self):
        self.d = GQuadruplexDetector()
        # G4: raw_min=0.5, raw_max=4.0, length_cap=120

    def test_at_cap_max_raw_gives_3(self):
        score = self.d.normalize_score(4.0, 120)
        assert score == 3.0

    def test_above_cap_clamped_to_3(self):
        score = self.d.normalize_score(4.0, 240)
        assert score == 3.0

    def test_min_raw_gives_1(self):
        score = self.d.normalize_score(0.5, 120)
        assert score == 1.0

    def test_half_length_reduces_score(self):
        score_full = self.d.normalize_score(4.0, 120)
        score_half = self.d.normalize_score(4.0, 60)
        assert score_half < score_full
        assert score_half == pytest.approx(2.0, abs=0.01)

    def test_score_always_at_least_1(self):
        for raw in [0.0, 0.5, 1.0, 2.0, 4.0]:
            for length in [1, 10, 50, 120, 200, 1000]:
                score = self.d.normalize_score(raw, length)
                assert score >= 1.0, f"score={score} < 1.0 for raw={raw}, length={length}"

    def test_score_never_exceeds_3(self):
        for raw in [0.0, 0.5, 1.0, 2.0, 4.0, 10.0, 100.0]:
            for length in [1, 10, 50, 120, 200, 1000]:
                score = self.d.normalize_score(raw, length)
                assert score <= 3.0, f"score={score} > 3.0 for raw={raw}, length={length}"

    def test_subclass_does_not_break_formula(self):
        score = self.d.normalize_score(2.0, 60, "Canonical intramolecular G4")
        assert 1.0 <= score <= 3.0

    def test_disease_subclass_slipped(self):
        d = SlippedDNADetector()
        # STR subclass with disease expansion length
        score = d.normalize_score(2.5, 500, "STR")
        assert 1.0 <= score <= 3.0

    def test_curved_local_vs_global_differ(self):
        d = CurvedDNADetector()
        local = d.normalize_score(0.8, 50, "Local Curvature")
        glob = d.normalize_score(0.8, 50, "Global Curvature")
        # local cap=50 → length_component=1.0 → higher score
        # global cap=120 → length_component≈0.417 → lower score
        assert local > glob


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
