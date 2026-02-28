"""Tests for i-Motif detector: RC scanning and extended-loop (Relaxed) detection.

Verifies that the i-Motif detector no longer under-predicts by:
1. Detecting i-motifs on the minus strand when the G-rich strand is provided.
2. Detecting Relaxed i-motifs with loops 8-12 nt (previously missed entirely).
3. Still classifying short-loop sequences as Canonical i-motif.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
from Detectors.imotif.detector import IMotifDetector
from Utilities.detectors_utils import revcomp


@pytest.fixture
def detector():
    return IMotifDetector()


CANONICAL_SEQ = "CCCCTCCCCTCCCCTCCCC"   # Validated Gehring 1993 (forward C-rich)
EXTENDED_SEQ  = "CCCCAAAAAAAACCCCAAAAAAAACCCCAAAAAAAACCCC"  # 8-nt loops


class TestReverseComplementScanning:
    """i-Motifs should be detected even when the G-rich strand is provided."""

    def test_canonical_forward_strand_detected(self, detector):
        motifs = detector.detect_motifs(CANONICAL_SEQ, "fwd")
        assert len(motifs) == 1
        assert motifs[0]["Strand"] == "+"
        assert motifs[0]["Subclass"] == "Canonical i-motif"

    def test_g_rich_strand_detected_on_minus(self, detector):
        g_rich = revcomp(CANONICAL_SEQ)         # GGGGAGGGGAGGGGAGGGG
        motifs = detector.detect_motifs(g_rich, "rc")
        assert len(motifs) == 1, "G-rich strand should yield one i-motif on the minus strand"
        assert motifs[0]["Strand"] == "-"
        assert motifs[0]["Subclass"] == "Canonical i-motif"

    def test_minus_strand_positions_within_sequence(self, detector):
        g_rich = revcomp(CANONICAL_SEQ)
        motifs = detector.detect_motifs(g_rich, "rc")
        assert len(motifs) == 1
        m = motifs[0]
        assert 1 <= m["Start"] <= len(g_rich)
        assert m["Start"] <= m["End"] <= len(g_rich)

    def test_minus_strand_structural_features_c_rich(self, detector):
        """C-tract analysis must use the C-rich complement, not the G-rich input."""
        g_rich = revcomp(CANONICAL_SEQ)
        motifs = detector.detect_motifs(g_rich, "rc")
        assert len(motifs) == 1
        m = motifs[0]
        assert m["Num_Stems"] >= 3, "C-tract analysis should reflect the C-rich complement"

    def test_rc_of_extended_loop_detected(self, detector):
        g_rich_ext = revcomp(EXTENDED_SEQ)
        motifs = detector.detect_motifs(g_rich_ext, "rc_ext")
        assert len(motifs) == 1
        assert motifs[0]["Strand"] == "-"
        assert motifs[0]["Subclass"] == "Relaxed i-motif"


class TestExtendedLoopDetection:
    """Relaxed i-motifs (loops 8-12 nt) must be detected after the fix."""

    def test_8nt_loop_detected_as_relaxed(self, detector):
        motifs = detector.detect_motifs(EXTENDED_SEQ, "ext")
        assert len(motifs) == 1, "8-nt loop i-motif should be detected"
        assert motifs[0]["Subclass"] == "Relaxed i-motif"

    def test_12nt_loop_detected_as_relaxed(self, detector):
        seq_12 = "CCCCTTTTTTTTTTTTCCCCTTTTTTTTTTTTCCCCTTTTTTTTTTTTCCCC"
        motifs = detector.detect_motifs(seq_12, "ext12")
        assert len(motifs) == 1
        assert motifs[0]["Subclass"] == "Relaxed i-motif"

    def test_relaxed_criterion_text(self, detector):
        motifs = detector.detect_motifs(EXTENDED_SEQ, "ext")
        assert len(motifs) == 1
        assert "Relaxed" in motifs[0]["Criterion"] or "12" in motifs[0]["Criterion"]

    def test_relaxed_type_of_repeat(self, detector):
        motifs = detector.detect_motifs(EXTENDED_SEQ, "ext")
        assert len(motifs) == 1
        assert "Relaxed" in motifs[0]["Type_Of_Repeat"]


class TestCanonicalStillPreferred:
    """Short-loop sequences should still be classified as Canonical i-motif."""

    def test_1nt_loops_are_canonical(self, detector):
        motifs = detector.detect_motifs(CANONICAL_SEQ, "can")
        assert motifs[0]["Subclass"] == "Canonical i-motif"

    def test_3nt_loops_are_canonical(self, detector):
        seq = "CCCATCCCATCCCATCCC"
        motifs = detector.detect_motifs(seq, "can3")
        assert len(motifs) == 1
        assert motifs[0]["Subclass"] == "Canonical i-motif"

    def test_7nt_loops_are_canonical(self, detector):
        seq = "CCCCAAAAAAACCCCAAAAAAACCCCAAAAAAACCCC"
        motifs = detector.detect_motifs(seq, "can7")
        assert len(motifs) >= 1
        assert motifs[0]["Subclass"] == "Canonical i-motif"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
