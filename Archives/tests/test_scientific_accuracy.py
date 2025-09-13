"""
Test suite for scientific accuracy of scoring algorithms

This module validates that the implemented scoring algorithms produce
scientifically accurate results consistent with published literature.
"""

import pytest
import sys
import os

# Add the parent directory to the path to import modules
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from motifs.base import (
    g4hunter_score,
    triplex_stability_score, 
    z_dna_score,
    curvature_score
)


class TestG4HunterAccuracy:
    """Test G4Hunter algorithm accuracy against known sequences"""
    
    def test_canonical_g4_sequence(self):
        """Test on a canonical G4-forming sequence"""
        # Well-known G4 sequence from literature
        g4_seq = "GGGCCCGGGCCCGGGCCCGGG"
        score = g4hunter_score(g4_seq)
        
        # Should be positive and significant
        assert score > 0.5, f"G4 sequence should have positive score, got {score}"
        assert score <= 2.0, f"Score should be within range, got {score}"
    
    def test_c_rich_sequence(self):
        """Test on C-rich sequence (should be negative)"""
        c_rich_seq = "CCCGGGCCCGGGCCCGGGCCC"
        score = g4hunter_score(c_rich_seq)
        
        # Should be negative due to C-runs
        assert score < 0, f"C-rich sequence should have negative score, got {score}"
        assert score >= -2.0, f"Score should be within range, got {score}"
    
    def test_no_runs_sequence(self):
        """Test on sequence with no G/C runs"""
        neutral_seq = "ATGATGATGATGATGATGATG"
        score = g4hunter_score(neutral_seq)
        
        # Should be near zero (no G/C runs of 2+)
        assert abs(score) < 0.1, f"Neutral sequence should have low score, got {score}"
    
    def test_empty_sequence(self):
        """Test edge case of empty sequence"""
        assert g4hunter_score("") == 0.0
        assert g4hunter_score(None) == 0.0
    
    def test_score_range(self):
        """Test that scores stay within expected range"""
        test_sequences = [
            "GGGGGGGGGGGGGGGGGGGG",  # All G
            "CCCCCCCCCCCCCCCCCCCC",  # All C
            "ATATATATATATATATATAT",  # No G/C
            "GGTGGCGGAGGCGGAGGCGG",  # Mixed
        ]
        
        for seq in test_sequences:
            score = g4hunter_score(seq)
            assert -2.0 <= score <= 2.0, f"Score {score} out of range for {seq}"


class TestTriplexStabilityAccuracy:
    """Test triplex stability scoring accuracy"""
    
    def test_homopurine_sequence(self):
        """Test high stability for homopurine tract"""
        homopurine = "AAAAAGAAAAAGAAAAAGAA"
        score = triplex_stability_score(homopurine)
        
        assert score > 0.7, f"Homopurine should have high stability, got {score}"
        assert score <= 1.0, f"Score should be within range, got {score}"
    
    def test_homopyrimidine_sequence(self):
        """Test high stability for homopyrimidine tract"""
        homopyrimidine = "CCCCTCCCCTCCCCTCCCCT"
        score = triplex_stability_score(homopyrimidine)
        
        assert score > 0.7, f"Homopyrimidine should have high stability, got {score}"
        assert score <= 1.0, f"Score should be within range, got {score}"
    
    def test_mixed_sequence(self):
        """Test lower stability for mixed sequence"""
        mixed = "ATGCATGCATGCATGCATGC"
        score = triplex_stability_score(mixed)
        
        assert score < 0.6, f"Mixed sequence should have lower stability, got {score}"
    
    def test_length_effect(self):
        """Test that longer sequences get length bonus"""
        short_pur = "AAAAAGGGG"  # 9 bp
        long_pur = "AAAAAGAAAAGAAAAGAAAAG"  # 21 bp
        
        short_score = triplex_stability_score(short_pur)
        long_score = triplex_stability_score(long_pur)
        
        assert long_score > short_score, "Longer tracts should score higher"
    
    def test_score_range(self):
        """Test that scores stay within [0.0, 1.0]"""
        test_sequences = [
            "AAAAAAAAAAAAAAAAAAA",  # All A
            "TTTTTTTTTTTTTTTTTTT",  # All T  
            "GGGGGGGGGGGGGGGGGGG",  # All G
            "CCCCCCCCCCCCCCCCCCC",  # All C
            "ATGCATGCATGCATGCATG",  # Mixed
        ]
        
        for seq in test_sequences:
            score = triplex_stability_score(seq)
            assert 0.0 <= score <= 1.0, f"Score {score} out of range for {seq}"


class TestZDNAAccuracy:
    """Test Z-DNA propensity scoring accuracy"""
    
    def test_cg_alternating_sequence(self):
        """Test high propensity for CG alternating sequence"""
        cg_alternating = "CGCGCGCGCGCGCGCGCGCG"
        score = z_dna_score(cg_alternating)
        
        assert score > 0.8, f"CG alternating should have high Z-DNA propensity, got {score}"
        assert score <= 1.0, f"Score should be within range, got {score}"
    
    def test_gc_alternating_sequence(self):
        """Test high propensity for GC alternating sequence"""  
        gc_alternating = "GCGCGCGCGCGCGCGCGCGC"
        score = z_dna_score(gc_alternating)
        
        assert score > 0.8, f"GC alternating should have high Z-DNA propensity, got {score}"
    
    def test_ca_tg_alternating(self):
        """Test moderate propensity for CA/TG alternating"""
        ca_alternating = "CACACACACACACACACACA"
        score = z_dna_score(ca_alternating)
        
        assert 0.25 < score < 0.8, f"CA alternating should have moderate propensity, got {score}"
    
    def test_non_alternating_sequence(self):
        """Test low propensity for non-alternating sequence"""
        non_alt = "AAAAATTTTGGGGCCCCAAA"
        score = z_dna_score(non_alt)
        
        assert score < 0.3, f"Non-alternating should have low propensity, got {score}"
    
    def test_length_effect(self):
        """Test length dependency"""
        short_cg = "CGCGCG"  # 6 bp
        long_cg = "CGCGCGCGCGCGCGCGCGCG"  # 20 bp
        
        short_score = z_dna_score(short_cg)
        long_score = z_dna_score(long_cg)
        
        # Both should be high, but longer should benefit from length factor
        assert short_score > 0.4, "Short CG sequence should still score well"
        assert long_score >= short_score, "Longer sequence should score at least as well"
    
    def test_score_range(self):
        """Test that scores stay within [0.0, 1.0]"""
        test_sequences = [
            "CGCGCGCGCGCGCGCGCGCG",  # Perfect alternating
            "AAAAATTTTGGGGCCCCAAA",  # No alternation
            "ATGCATGCATGCATGCATGC",  # Some alternation
            "CG",  # Very short
        ]
        
        for seq in test_sequences:
            score = z_dna_score(seq)
            assert 0.0 <= score <= 1.0, f"Score {score} out of range for {seq}"


class TestCurvatureAccuracy:
    """Test DNA curvature prediction accuracy"""
    
    def test_a_tract_sequence(self):
        """Test high curvature for A-tract containing sequence"""
        a_tract_seq = "ATCAAAAAATCAAAAAATCA"
        score = curvature_score(a_tract_seq)
        
        assert score > 0.5, f"A-tract sequence should have high curvature, got {score}"
        assert score <= 1.0, f"Score should be within range, got {score}"
    
    def test_phased_a_tracts(self):
        """Test enhanced curvature for phased A-tracts"""
        phased = "AAAAAATCGATCAAAAAATCGATC"  # A-tracts spaced ~10 bp apart
        unphased = "AAAAAAAAAAATCGATCGATCGAT"  # A-tracts not phased
        
        phased_score = curvature_score(phased)
        unphased_score = curvature_score(unphased)
        
        assert phased_score > unphased_score, "Phased A-tracts should score higher"
    
    def test_gc_rich_sequence(self):
        """Test low curvature for GC-rich sequence"""
        gc_rich = "GCGCGCGCGCGCGCGCGCGC"
        score = curvature_score(gc_rich)
        
        assert score < 0.3, f"GC-rich sequence should have low curvature, got {score}"
    
    def test_no_a_tracts(self):
        """Test low curvature for sequence without A-tracts"""
        no_a_tracts = "TCGTCGTCGTCGTCGTCGTC"
        score = curvature_score(no_a_tracts)
        
        assert score < 0.4, f"Sequence without A-tracts should have low curvature, got {score}"
    
    def test_window_size_effect(self):
        """Test that window size affects scoring appropriately"""
        test_seq = "AAAAAATCGATCAAAAAATCGATC"
        
        score_10 = curvature_score(test_seq, window_size=10)
        score_5 = curvature_score(test_seq, window_size=5)
        
        # Both should detect curvature, window size affects resolution
        assert score_10 > 0.3, "Should detect curvature with window 10"
        assert score_5 > 0.3, "Should detect curvature with window 5"
    
    def test_score_range(self):
        """Test that scores stay within [0.0, 1.0]"""
        test_sequences = [
            "AAAAAATCGATCAAAAAATCGATC",  # Strong A-tracts
            "GCGCGCGCGCGCGCGCGCGCGCGC",  # GC-rich
            "ATCGATCGATCGATCGATCGATCG",  # Mixed
            "AAAAAAAAAAAAAAAAAAAAAAAA",  # Long A-tract
        ]
        
        for seq in test_sequences:
            score = curvature_score(seq)
            assert 0.0 <= score <= 1.0, f"Score {score} out of range for {seq}"


class TestScientificValidation:
    """Integration tests for scientific validation"""
    
    def test_known_g4_sequence_from_literature(self):
        """Test on a well-studied G4 sequence from literature"""
        # Human telomeric G4 sequence
        telomeric_g4 = "GGGCTAGGGCTAGGGCTAGGG"
        
        g4_score = g4hunter_score(telomeric_g4)
        z_score = z_dna_score(telomeric_g4)
        triplex_score = triplex_stability_score(telomeric_g4)
        
        # G4 should dominate
        assert g4_score > z_score, "G4 score should be higher than Z-DNA for G4 sequence"
        assert g4_score > 0.5, "Should have significant G4 propensity"
    
    def test_known_zdna_sequence(self):
        """Test on a known Z-DNA forming sequence"""
        # (CG)10 repeat - strong Z-DNA former
        zdna_seq = "CGCGCGCGCGCGCGCGCGCG"
        
        z_score = z_dna_score(zdna_seq)
        g4_score = g4hunter_score(zdna_seq)
        
        # Z-DNA should dominate
        assert z_score > 0.7, "Should have high Z-DNA propensity"
        assert z_score > abs(g4_score), "Z-DNA score should exceed G4 score magnitude"
    
    def test_curvature_vs_other_scores(self):
        """Test that curvature sequences don't score high for other motifs"""
        curve_seq = "AAAAAATCGATCAAAAAATCGATC"
        
        curve_score = curvature_score(curve_seq)
        g4_score = g4hunter_score(curve_seq)
        z_score = z_dna_score(curve_seq)
        
        # Curvature should dominate
        assert curve_score > 0.4, "Should have significant curvature"
        assert curve_score > abs(g4_score), "Curvature should exceed G4 score"
        assert curve_score > z_score, "Curvature should exceed Z-DNA score"
    
    def test_algorithmic_consistency(self):
        """Test that algorithms are consistent and deterministic"""
        test_seq = "GGGCCCGGGCCCAAAAAATCGATC"
        
        # Run multiple times to ensure consistency
        scores1 = [
            g4hunter_score(test_seq),
            triplex_stability_score(test_seq),
            z_dna_score(test_seq),
            curvature_score(test_seq)
        ]
        
        scores2 = [
            g4hunter_score(test_seq),
            triplex_stability_score(test_seq), 
            z_dna_score(test_seq),
            curvature_score(test_seq)
        ]
        
        assert scores1 == scores2, "Scoring should be deterministic"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])