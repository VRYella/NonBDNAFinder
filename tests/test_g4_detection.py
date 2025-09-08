"""
Unit tests for G4 detection using canonical test sequences.
"""

import pytest
from motif_detectors import G4Detector
from motifs.base import Candidate


class TestG4Detection:
    """Test G-quadruplex detection"""
    
    def setup_method(self):
        """Setup G4 detector for each test"""
        self.detector = G4Detector()
    
    def test_detector_initialization(self):
        """Test G4 detector initialization"""
        assert self.detector.class_name == 'g_quadruplex'
        assert self.detector.class_id == 6
        assert len(self.detector.patterns) > 0
    
    def test_canonical_g4_detection(self):
        """Test detection of canonical G4 sequence"""
        # Known G4-forming sequence
        test_sequence = "GGGTTTGGGTTTGGGTTTGGG"
        
        candidates = self.detector.detect(
            seq=test_sequence,
            seq_name="test_seq",
            contig="test_contig", 
            offset=0
        )
        
        # Should detect at least one G4
        assert len(candidates) > 0
        
        # Check candidate properties
        for candidate in candidates:
            assert isinstance(candidate, Candidate)
            assert candidate.class_name == 'g_quadruplex'
            assert candidate.sequence_name == 'test_seq'
            assert candidate.contig == 'test_contig'
            assert candidate.start > 0  # 1-based coordinates
            assert candidate.end >= candidate.start
            assert len(candidate.matched_seq) > 0
    
    def test_g4_scoring(self):
        """Test G4 scoring functionality"""
        test_sequence = "GGGTTTGGGTTTGGGTTTGGG"
        
        # Detect candidates
        candidates = self.detector.detect(
            seq=test_sequence,
            seq_name="test_seq",
            contig="test_contig",
            offset=0
        )
        
        # Score candidates
        scored_candidates = self.detector.score(candidates)
        
        assert len(scored_candidates) == len(candidates)
        
        for candidate in scored_candidates:
            assert candidate.raw_score is not None
            assert candidate.scoring_method == "G4Hunter_v2"
            assert isinstance(candidate.raw_score, float)
    
    def test_negative_sequence(self):
        """Test with sequence that should not form G4"""
        # AT-rich sequence with no G4 potential
        negative_sequence = "ATATATATATATATATATATATATATAT"
        
        candidates = self.detector.detect(
            seq=negative_sequence,
            seq_name="negative_test",
            contig="negative_contig",
            offset=0
        )
        
        # Should detect few or no G4s
        assert len(candidates) <= 1  # Allow for very permissive patterns
    
    def test_multiple_g4s(self):
        """Test sequence with multiple G4 structures"""
        # Sequence with two potential G4 regions
        multi_g4_sequence = "GGGTTTGGGTTTGGGTTTGGGAAAAAGGGTTTGGGTTTGGGTTTGGG"
        
        candidates = self.detector.detect(
            seq=multi_g4_sequence,
            seq_name="multi_test",
            contig="multi_contig",
            offset=0
        )
        
        # Should detect multiple candidates
        assert len(candidates) >= 2
        
        # Check that positions are different
        positions = [(c.start, c.end) for c in candidates]
        assert len(set(positions)) >= 2  # At least 2 unique positions
    
    def test_imperfect_g4(self):
        """Test detection of imperfect G4 structures"""
        # G4 with 2-G runs (imperfect)
        imperfect_sequence = "GGTTTGGGTTTGGGTTTGGG"
        
        candidates = self.detector.detect(
            seq=imperfect_sequence,
            seq_name="imperfect_test",
            contig="imperfect_contig",
            offset=0
        )
        
        # May or may not detect depending on pattern stringency
        # Just ensure no errors occur
        assert isinstance(candidates, list)
        for candidate in candidates:
            assert isinstance(candidate, Candidate)
    
    def test_long_sequence(self):
        """Test with longer sequence containing G4s"""
        # Longer sequence with embedded G4
        long_sequence = ("ATCGATCGATCGATCG" * 10 + 
                        "GGGTTTGGGTTTGGGTTTGGG" +
                        "ATCGATCGATCGATCG" * 10)
        
        candidates = self.detector.detect(
            seq=long_sequence,
            seq_name="long_test",
            contig="long_contig",
            offset=0
        )
        
        # Should find the embedded G4
        assert len(candidates) > 0
        
        # Check that coordinates make sense
        for candidate in candidates:
            assert candidate.start > 0
            assert candidate.end <= len(long_sequence)
    
    def test_offset_handling(self):
        """Test proper handling of sequence offset"""
        test_sequence = "GGGTTTGGGTTTGGGTTTGGG"
        offset = 1000
        
        candidates = self.detector.detect(
            seq=test_sequence,
            seq_name="offset_test",
            contig="offset_contig",
            offset=offset
        )
        
        # Check that offset is properly added
        for candidate in candidates:
            assert candidate.start > offset
            assert candidate.end > offset
    
    def test_edge_cases(self):
        """Test edge cases"""
        # Empty sequence
        empty_candidates = self.detector.detect("", "empty", "empty", 0)
        assert len(empty_candidates) == 0
        
        # Very short sequence
        short_candidates = self.detector.detect("GGG", "short", "short", 0)
        assert isinstance(short_candidates, list)
        
        # Single nucleotide
        single_candidates = self.detector.detect("G", "single", "single", 0)
        assert isinstance(single_candidates, list)