#!/usr/bin/env python3
"""
Comprehensive tests for the enhanced overlap resolution system.

This test suite validates the overlap detection and resolution functionality
for all supported strategies and edge cases.
"""

import pytest
import sys
import os
from typing import List

# Add the project root to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from overlap_resolution import (
    EnhancedOverlapResolver, OverlapConfig, OverlapStrategy,
    resolve_motif_overlaps
)
from motifs.base import Candidate


class TestOverlapResolution:
    """Test suite for overlap resolution functionality"""
    
    def create_test_candidate(self, seq_name: str, class_name: str, subclass: str,
                             start: int, end: int, score: float = 1.0,
                             motif_id: int = 1) -> Candidate:
        """
        Create a test candidate with specified parameters.
        
        Args:
            seq_name: Sequence name
            class_name: Motif class name
            subclass: Motif subclass
            start: Start position
            end: End position
            score: Raw score
            motif_id: Motif identifier
            
        Returns:
            Candidate object for testing
        """
        return Candidate(
            sequence_name=seq_name,
            contig=seq_name,
            class_id=1,
            class_name=class_name,
            subclass=subclass,
            motif_id=motif_id,
            start=start,
            end=end,
            length=end - start + 1,
            matched_seq=b"A" * (end - start + 1),
            pattern_name=f"test_{class_name}",
            raw_score=score
        )
    
    def test_no_overlaps(self):
        """Test that non-overlapping candidates are preserved"""
        candidates = [
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 10, 20, 0.8),
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 30, 40, 0.7),
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 50, 60, 0.9)
        ]
        
        config = OverlapConfig(strategy=OverlapStrategy.HIGHEST_SCORE)
        resolver = EnhancedOverlapResolver(config)
        
        result = resolver.resolve_overlaps(candidates)
        
        assert len(result) == 3
        assert all(c in result for c in candidates)
    
    def test_highest_score_strategy(self):
        """Test overlap resolution using highest score strategy"""
        candidates = [
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 10, 25, 0.8),
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 20, 35, 0.9),  # Higher score
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 40, 50, 0.7)
        ]
        
        config = OverlapConfig(strategy=OverlapStrategy.HIGHEST_SCORE)
        resolver = EnhancedOverlapResolver(config)
        
        result = resolver.resolve_overlaps(candidates)
        
        # Should keep the candidate with score 0.9 and the non-overlapping one
        assert len(result) == 2
        scores = [c.raw_score for c in result]
        assert 0.9 in scores
        assert 0.7 in scores
        assert 0.8 not in scores
    
    def test_longest_motif_strategy(self):
        """Test overlap resolution using longest motif strategy"""
        candidates = [
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 10, 20, 0.8),  # Length 11
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 15, 30, 0.7),  # Length 16, longer
        ]
        
        config = OverlapConfig(strategy=OverlapStrategy.LONGEST_MOTIF)
        resolver = EnhancedOverlapResolver(config)
        
        result = resolver.resolve_overlaps(candidates)
        
        # Should keep the longer candidate
        assert len(result) == 1
        assert result[0].length == 16
        assert result[0].end == 30
    
    def test_scientific_priority_strategy(self):
        """Test overlap resolution using scientific priority"""
        candidates = [
            self.create_test_candidate("seq1", "curved_dna", "A_phased", 10, 25, 0.9),
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 20, 35, 0.7),  # Higher priority
        ]
        
        config = OverlapConfig(
            strategy=OverlapStrategy.SCIENTIFIC_PRIORITY,
            same_class_only=False  # Allow cross-class resolution
        )
        resolver = EnhancedOverlapResolver(config)
        
        result = resolver.resolve_overlaps(candidates)
        
        # Should keep G4 due to higher scientific priority
        assert len(result) == 1
        assert result[0].class_name == "g_quadruplex"
    
    def test_same_class_only_grouping(self):
        """Test that same_class_only setting works correctly"""
        candidates = [
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 10, 25, 0.8),
            self.create_test_candidate("seq1", "triplex", "homopurine", 20, 35, 0.9),  # Different class
        ]
        
        config = OverlapConfig(
            strategy=OverlapStrategy.HIGHEST_SCORE,
            same_class_only=True  # Only resolve within same class
        )
        resolver = EnhancedOverlapResolver(config)
        
        result = resolver.resolve_overlaps(candidates)
        
        # Should keep both since they're different classes
        assert len(result) == 2
    
    def test_cross_class_resolution(self):
        """Test cross-class overlap resolution"""
        candidates = [
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 10, 25, 0.8),
            self.create_test_candidate("seq1", "triplex", "homopurine", 20, 35, 0.9),  # Different class, higher score
        ]
        
        config = OverlapConfig(
            strategy=OverlapStrategy.HIGHEST_SCORE,
            same_class_only=False  # Allow cross-class resolution
        )
        resolver = EnhancedOverlapResolver(config)
        
        result = resolver.resolve_overlaps(candidates)
        
        # Should keep triplex due to higher score
        assert len(result) == 1
        assert result[0].class_name == "triplex"
    
    def test_merge_compatible_strategy(self):
        """Test merging of compatible overlapping motifs"""
        candidates = [
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 10, 25, 0.8),
            self.create_test_candidate("seq1", "g_quadruplex", "relaxed", 20, 35, 0.7),
        ]
        
        config = OverlapConfig(
            strategy=OverlapStrategy.MERGE_COMPATIBLE,
            merge_threshold=0.3  # Low threshold to allow merging
        )
        resolver = EnhancedOverlapResolver(config)
        
        result = resolver.resolve_overlaps(candidates)
        
        # Should create a merged candidate
        assert len(result) == 1
        merged = result[0]
        assert merged.start == 10  # Start of first candidate
        assert merged.end == 35    # End of second candidate
        assert "merged" in merged.subclass
    
    def test_keep_all_strategy(self):
        """Test keeping all overlapping motifs with annotation"""
        candidates = [
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 10, 25, 0.8),
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 20, 35, 0.9),
        ]
        
        config = OverlapConfig(strategy=OverlapStrategy.KEEP_ALL)
        resolver = EnhancedOverlapResolver(config)
        
        result = resolver.resolve_overlaps(candidates)
        
        # Should keep both candidates but mark overlaps
        assert len(result) == 2
        for candidate in result:
            assert candidate.overlap_classes is not None
            assert len(candidate.overlap_classes) > 0
    
    def test_minimal_overlap_threshold(self):
        """Test minimum overlap percentage threshold"""
        candidates = [
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 10, 20, 0.8),  # Length 11
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 20, 30, 0.9),  # Length 11, overlap=1
        ]
        
        # Overlap is 1bp out of 11bp = ~9%, below 10% threshold
        config = OverlapConfig(
            strategy=OverlapStrategy.HIGHEST_SCORE,
            min_overlap_percent=0.1  # 10% minimum
        )
        resolver = EnhancedOverlapResolver(config)
        
        result = resolver.resolve_overlaps(candidates)
        
        # Should keep both due to minimal overlap
        assert len(result) == 2
    
    def test_significant_overlap_threshold(self):
        """Test overlap resolution with significant overlap"""
        candidates = [
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 10, 25, 0.8),  # Length 16
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 20, 35, 0.9),  # Length 16, overlap=6
        ]
        
        # Overlap is 6bp out of 16bp = 37.5%, above 10% threshold
        config = OverlapConfig(
            strategy=OverlapStrategy.HIGHEST_SCORE,
            min_overlap_percent=0.1  # 10% minimum
        )
        resolver = EnhancedOverlapResolver(config)
        
        result = resolver.resolve_overlaps(candidates)
        
        # Should resolve overlap, keeping higher score
        assert len(result) == 1
        assert result[0].raw_score == 0.9
    
    def test_multiple_sequence_handling(self):
        """Test handling of candidates from multiple sequences"""
        candidates = [
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 10, 25, 0.8),
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 20, 35, 0.9),
            self.create_test_candidate("seq2", "g_quadruplex", "canonical", 10, 25, 0.7),  # Different sequence
            self.create_test_candidate("seq2", "g_quadruplex", "canonical", 20, 35, 0.6),  # Different sequence
        ]
        
        config = OverlapConfig(strategy=OverlapStrategy.HIGHEST_SCORE)
        resolver = EnhancedOverlapResolver(config)
        
        result = resolver.resolve_overlaps(candidates)
        
        # Should resolve overlaps within each sequence separately
        assert len(result) == 2
        seq_names = [c.sequence_name for c in result]
        assert "seq1" in seq_names
        assert "seq2" in seq_names
    
    def test_equal_scores_tiebreaking(self):
        """Test tiebreaking when candidates have equal scores"""
        candidates = [
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 10, 20, 0.8),  # Length 11
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 15, 30, 0.8),  # Length 16, longer
        ]
        
        config = OverlapConfig(strategy=OverlapStrategy.HIGHEST_SCORE)
        resolver = EnhancedOverlapResolver(config)
        
        result = resolver.resolve_overlaps(candidates)
        
        # Should keep longer motif as tiebreaker
        assert len(result) == 1
        assert result[0].length == 16
    
    def test_convenience_function(self):
        """Test the convenience function for overlap resolution"""
        candidates = [
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 10, 25, 0.8),
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 20, 35, 0.9),
        ]
        
        result = resolve_motif_overlaps(
            candidates,
            strategy=OverlapStrategy.HIGHEST_SCORE
        )
        
        assert len(result) == 1
        assert result[0].raw_score == 0.9
    
    def test_empty_input_handling(self):
        """Test handling of empty candidate list"""
        config = OverlapConfig(strategy=OverlapStrategy.HIGHEST_SCORE)
        resolver = EnhancedOverlapResolver(config)
        
        result = resolver.resolve_overlaps([])
        
        assert result == []
    
    def test_single_candidate_handling(self):
        """Test handling of single candidate"""
        candidates = [
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 10, 25, 0.8)
        ]
        
        config = OverlapConfig(strategy=OverlapStrategy.HIGHEST_SCORE)
        resolver = EnhancedOverlapResolver(config)
        
        result = resolver.resolve_overlaps(candidates)
        
        assert len(result) == 1
        assert result[0] == candidates[0]
    
    def test_complex_overlap_scenario(self):
        """Test complex scenario with multiple overlapping motifs"""
        candidates = [
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 10, 30, 0.9),   # Best score
            self.create_test_candidate("seq1", "g_quadruplex", "relaxed", 15, 25, 0.7),    # Fully contained
            self.create_test_candidate("seq1", "g_quadruplex", "bulged", 25, 40, 0.8),     # Partial overlap
            self.create_test_candidate("seq1", "g_quadruplex", "canonical", 50, 60, 0.6),  # No overlap
        ]
        
        config = OverlapConfig(strategy=OverlapStrategy.HIGHEST_SCORE)
        resolver = EnhancedOverlapResolver(config)
        
        result = resolver.resolve_overlaps(candidates)
        
        # Should keep the highest scoring overlapping motif (0.9) and the non-overlapping one (0.6)
        assert len(result) == 2
        scores = sorted([c.raw_score for c in result])
        assert scores == [0.6, 0.9]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])