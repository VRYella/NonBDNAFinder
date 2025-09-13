"""
Additional unit tests for triplex detection and validation flow.
"""

import pytest
from motif_detectors import TriplexDetector, get_detector
from motifs.base import Candidate


class TestTriplexDetection:
    """Test triplex DNA detection"""
    
    def setup_method(self):
        """Setup triplex detector for each test"""
        self.detector = TriplexDetector()
    
    def test_detector_initialization(self):
        """Test triplex detector initialization"""
        assert self.detector.class_name == 'triplex'
        assert self.detector.class_id == 5
    
    def test_homopurine_detection(self):
        """Test detection of homopurine tracts"""
        # Long homopurine sequence
        purine_sequence = "A" * 20 + "G" * 10 + "A" * 5
        
        candidates = self.detector.detect(
            seq=purine_sequence,
            seq_name="purine_test",
            contig="purine_contig",
            offset=0
        )
        
        # Should detect homopurine tract
        assert len(candidates) > 0
        
        purine_candidates = [c for c in candidates if c.subclass == 'Homopurine']
        assert len(purine_candidates) > 0
    
    def test_homopyrimidine_detection(self):
        """Test detection of homopyrimidine tracts"""
        # Long homopyrimidine sequence
        pyrimidine_sequence = "C" * 15 + "T" * 10 + "C" * 5
        
        candidates = self.detector.detect(
            seq=pyrimidine_sequence,
            seq_name="pyrimidine_test", 
            contig="pyrimidine_contig",
            offset=0
        )
        
        # Should detect homopyrimidine tract
        pyrimidine_candidates = [c for c in candidates if c.subclass == 'Homopyrimidine']
        assert len(pyrimidine_candidates) > 0
    
    def test_triplex_scoring(self):
        """Test triplex stability scoring"""
        test_sequence = "A" * 20  # Pure homopurine
        
        candidates = self.detector.detect(
            seq=test_sequence,
            seq_name="score_test",
            contig="score_contig",
            offset=0
        )
        
        if candidates:
            scored_candidates = self.detector.score(candidates)
            
            for candidate in scored_candidates:
                assert candidate.raw_score is not None
                assert candidate.scoring_method == "Triplex_stability"
                assert 0.0 <= candidate.raw_score <= 1.0  # Should be a proportion
    
    def test_mixed_sequence(self):
        """Test sequence with both purine and pyrimidine regions"""
        mixed_sequence = "A" * 15 + "ATCGATCG" * 5 + "C" * 15
        
        candidates = self.detector.detect(
            seq=mixed_sequence,
            seq_name="mixed_test",
            contig="mixed_contig", 
            offset=0
        )
        
        # May detect both types
        subtypes = {c.subclass for c in candidates}
        # Just ensure no errors and valid candidates
        for candidate in candidates:
            assert candidate.subclass in ['Homopurine', 'Homopyrimidine']


class TestAnchorValidationFlow:
    """Test anchor-based detection with validation"""
    
    def test_detector_factory(self):
        """Test detector factory function"""
        # Test all detector classes can be instantiated
        detector_classes = [
            'g_quadruplex', 'curved_dna', 'slipped_dna', 'cruciform',
            'r_loop', 'triplex', 'i_motif', 'z_dna', 'hybrid', 'cluster'
        ]
        
        for class_name in detector_classes:
            detector = get_detector(class_name)
            assert detector.class_name == class_name
            assert hasattr(detector, 'detect')
            assert hasattr(detector, 'score')
    
    def test_invalid_detector_class(self):
        """Test error handling for invalid detector class"""
        with pytest.raises(ValueError):
            get_detector('nonexistent_class')
    
    def test_candidate_creation(self):
        """Test candidate object creation"""
        detector = get_detector('g_quadruplex')
        
        candidate = detector.make_candidate(
            seq="GGGTTTGGGTTTGGGTTTGGG",
            seq_name="test",
            contig="test_contig",
            offset=100,
            start=0,
            end=20,
            pattern_name="test_pattern",
            subclass="test_subclass",
            motif_id=1
        )
        
        assert isinstance(candidate, Candidate)
        assert candidate.sequence_name == "test"
        assert candidate.contig == "test_contig"
        assert candidate.start == 101  # offset + start + 1 (1-based)
        assert candidate.end == 121    # offset + end + 1 (1-based)
        assert candidate.length == 21
        assert candidate.subclass == "test_subclass"
    
    def test_gc_content_calculation(self):
        """Test GC content calculation in candidates"""
        from motifs.base import calculate_gc_content
        
        # Test sequences with known GC content
        test_cases = [
            ("GGCC", 1.0),  # 100% GC
            ("AATT", 0.0),  # 0% GC  
            ("ATCG", 0.5),  # 50% GC
            ("GGGATTT", 3/7),  # 3 GCs out of 7
        ]
        
        for seq, expected_gc in test_cases:
            calculated_gc = calculate_gc_content(seq.encode())
            assert abs(calculated_gc - expected_gc) < 0.001
    
    def test_multiple_detector_coordination(self):
        """Test that multiple detectors can work on the same sequence"""
        test_sequence = "GGGTTTGGGTTTGGGTTTGGGAAAAAAAAATGCGTAAAAAAAAAA"
        
        # Run multiple detectors
        g4_detector = get_detector('g_quadruplex')
        curved_detector = get_detector('curved_dna')
        
        g4_candidates = g4_detector.detect(test_sequence, "multi_test", "multi_contig", 0)
        curved_candidates = curved_detector.detect(test_sequence, "multi_test", "multi_contig", 0)
        
        # Both should work without interference
        assert isinstance(g4_candidates, list)
        assert isinstance(curved_candidates, list)
        
        # G4 detector should find something in the G-rich region
        assert len(g4_candidates) > 0
        
        # Curved detector should find something in the A-rich region
        # (may or may not depending on exact algorithm)
        assert isinstance(curved_candidates, list)  # At least no errors