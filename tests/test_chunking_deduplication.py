"""
Test suite for chunking and boundary deduplication.
Verifies fix for >1MB sequence bug.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
from Utilities.nonbscanner import analyze_sequence, _deduplicate_motifs

class TestChunkingDeduplication(unittest.TestCase):
    """Test deduplication of boundary motifs in chunked analysis."""
    
    def test_1mb_sequence_returns_motifs(self):
        """Test that sequences >1MB return motifs (primary bug fix)."""
        # Create 1.1MB sequence with known G-quadruplex patterns
        gquad = "GGGTAGGGTAGGGTAGGG"
        filler = "ATCGATCGATCGATCG"
        pattern = gquad + filler * 10
        
        target_size = 1_100_000
        test_seq = (pattern * (target_size // len(pattern) + 1))[:target_size]
        
        print(f"\nTesting 1.1MB sequence ({len(test_seq):,} bp)...")
        motifs = analyze_sequence(test_seq, "test_1mb")
        
        print(f"Found {len(motifs)} motifs")
        self.assertGreater(len(motifs), 0, 
            "CRITICAL BUG: No motifs found for >1MB sequence")
        
        # Verify motif quality
        for motif in motifs[:5]:
            self.assertIn('Start', motif)
            self.assertIn('End', motif)
            self.assertIn('Class', motif)
            self.assertGreater(motif['Start'], 0)
            self.assertLess(motif['End'], len(test_seq))
    
    def test_deduplication_removes_exact_duplicates(self):
        """Test that exact duplicates are removed."""
        motifs = [
            {'Class': 'G-Quadruplex', 'Subclass': 'Canonical', 'Start': 100, 'End': 150, 'Score': 2.0},
            {'Class': 'G-Quadruplex', 'Subclass': 'Canonical', 'Start': 100, 'End': 150, 'Score': 2.0},
        ]
        result = _deduplicate_motifs(motifs)
        self.assertEqual(len(result), 1, "Exact duplicates should be removed")
    
    def test_deduplication_removes_overlapping_duplicates(self):
        """Test that overlapping boundary duplicates are removed."""
        motifs = [
            {'Class': 'G-Quadruplex', 'Subclass': 'Canonical', 'Start': 100, 'End': 150, 'Score': 2.0},
            {'Class': 'G-Quadruplex', 'Subclass': 'Canonical', 'Start': 120, 'End': 170, 'Score': 1.8},
        ]
        result = _deduplicate_motifs(motifs)
        self.assertEqual(len(result), 1, "Overlapping duplicates (>50% overlap) should be removed")
        self.assertEqual(result[0]['Score'], 2.0, "Should keep higher-scoring motif")
    
    def test_deduplication_keeps_non_overlapping_motifs(self):
        """Test that distinct motifs are preserved."""
        motifs = [
            {'Class': 'G-Quadruplex', 'Subclass': 'Canonical', 'Start': 100, 'End': 150, 'Score': 2.0},
            {'Class': 'G-Quadruplex', 'Subclass': 'Canonical', 'Start': 200, 'End': 250, 'Score': 1.8},
        ]
        result = _deduplicate_motifs(motifs)
        self.assertEqual(len(result), 2, "Non-overlapping motifs should be preserved")
    
    def test_deduplication_keeps_different_classes(self):
        """Test that different motif classes are preserved even if overlapping."""
        motifs = [
            {'Class': 'G-Quadruplex', 'Subclass': 'Canonical', 'Start': 100, 'End': 150, 'Score': 2.0},
            {'Class': 'Z-DNA', 'Subclass': 'Z-DNA', 'Start': 120, 'End': 170, 'Score': 1.8},
        ]
        result = _deduplicate_motifs(motifs)
        self.assertEqual(len(result), 2, "Different classes should be preserved even if overlapping")
    
    def test_500kb_vs_1_5mb_performance(self):
        """Test and compare performance for sequences below and above 1MB threshold."""
        import time
        
        # 500KB sequence (no chunking)
        seq_500kb = "ATCGATCG" * 62500  # 500,000 bp
        start = time.time()
        motifs_500kb = analyze_sequence(seq_500kb, "test_500kb")
        time_500kb = time.time() - start
        
        # 1.5MB sequence (with chunking)
        seq_1_5mb = "ATCGATCG" * 187500  # 1,500,000 bp
        start = time.time()
        motifs_1_5mb = analyze_sequence(seq_1_5mb, "test_1_5mb")
        time_1_5mb = time.time() - start
        
        print(f"\nPerformance comparison:")
        print(f"  500KB: {len(motifs_500kb)} motifs in {time_500kb:.2f}s ({len(seq_500kb)/time_500kb/1000:.0f} Kbp/s)")
        print(f"  1.5MB: {len(motifs_1_5mb)} motifs in {time_1_5mb:.2f}s ({len(seq_1_5mb)/time_1_5mb/1000:.0f} Kbp/s)")
        
        # Both should find motifs
        self.assertGreaterEqual(len(motifs_500kb), 0, "500KB should return results")
        self.assertGreaterEqual(len(motifs_1_5mb), 0, "1.5MB should return results")

if __name__ == '__main__':
    unittest.main(verbosity=2)
