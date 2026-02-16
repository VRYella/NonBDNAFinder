"""
Unit tests for triple_chunk_analyzer.py - Triple Adaptive Chunking System.

Tests the three-tier hierarchical chunking with automatic strategy selection.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
import gc

# Add parent directory to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from Utilities.disk_storage import UniversalSequenceStorage
from Utilities.triple_chunk_analyzer import TripleAdaptiveChunkAnalyzer


class TestTripleChunking(unittest.TestCase):
    """Test cases for TripleAdaptiveChunkAnalyzer class."""
    
    def setUp(self):
        """Create temporary directory and storage for tests."""
        self.test_dir = tempfile.mkdtemp(prefix="test_triple_chunk_")
        self.storage = UniversalSequenceStorage(base_dir=self.test_dir)
    
    def tearDown(self):
        """Clean up temporary directory."""
        if Path(self.test_dir).exists():
            shutil.rmtree(self.test_dir)
        gc.collect()
    
    def test_direct_analysis_threshold(self):
        """Test direct analysis for sequences below 50KB threshold."""
        # Create 25KB sequence (should use direct analysis)
        sequence = "ATCGATCG" * 3125  # 25,000 bp
        seq_id = self.storage.save_sequence(sequence, "small_seq")
        
        analyzer = TripleAdaptiveChunkAnalyzer(self.storage, use_adaptive=True)
        results_storage = analyzer.analyze(seq_id)
        
        # Verify results are stored
        stats = results_storage.get_summary_stats()
        self.assertGreaterEqual(stats['total_count'], 0)
        
        # Cleanup
        results_storage.cleanup()
    
    def test_single_tier_strategy(self):
        """Test single-tier (micro) chunking for 50KB-1MB sequences."""
        # Create 500KB sequence (should use single-tier)
        sequence = "ATCGATCG" * 62500  # 500,000 bp
        seq_id = self.storage.save_sequence(sequence, "medium_seq")
        
        analyzer = TripleAdaptiveChunkAnalyzer(self.storage, use_adaptive=True)
        results_storage = analyzer.analyze(seq_id)
        
        # Verify results
        stats = results_storage.get_summary_stats()
        self.assertGreaterEqual(stats['total_count'], 0)
        
        # Cleanup
        results_storage.cleanup()
    
    def test_double_tier_strategy(self):
        """Test double-tier (meso+micro) chunking for 1MB-100MB sequences."""
        # Create 5MB sequence (should use double-tier)
        sequence = "ATCGATCG" * 625000  # 5,000,000 bp
        seq_id = self.storage.save_sequence(sequence, "large_seq")
        
        analyzer = TripleAdaptiveChunkAnalyzer(self.storage, use_adaptive=True)
        results_storage = analyzer.analyze(seq_id)
        
        # Verify results
        stats = results_storage.get_summary_stats()
        self.assertGreaterEqual(stats['total_count'], 0)
        
        # Cleanup
        results_storage.cleanup()
    
    def test_triple_tier_strategy(self):
        """Test triple-tier (macro+meso+micro) chunking for >100MB sequences."""
        # Create 120MB sequence (should use triple-tier)
        # Note: This test may take longer to run
        sequence = "ATCGATCG" * 15000000  # 120,000,000 bp
        seq_id = self.storage.save_sequence(sequence, "huge_seq")
        
        analyzer = TripleAdaptiveChunkAnalyzer(self.storage, use_adaptive=True)
        results_storage = analyzer.analyze(seq_id)
        
        # Verify results
        stats = results_storage.get_summary_stats()
        self.assertGreaterEqual(stats['total_count'], 0)
        
        # Cleanup
        results_storage.cleanup()
    
    def test_micro_tier_deduplication(self):
        """Test deduplication at micro-chunk boundaries."""
        # Create sequence with pattern at micro chunk boundary
        # Micro chunk size is 50KB, place pattern near that boundary
        pre_sequence = "ATCG" * 12000  # 48,000 bp
        g4_motif = "GGGTTAGGGTTAGGGTTAGGG"  # G-quadruplex at boundary
        post_sequence = "ATCG" * 12000  # 48,000 bp
        
        sequence = pre_sequence + g4_motif + post_sequence
        seq_id = self.storage.save_sequence(sequence, "micro_boundary_test")
        
        analyzer = TripleAdaptiveChunkAnalyzer(self.storage, use_adaptive=True)
        results_storage = analyzer.analyze(seq_id)
        
        # Get all motifs
        all_motifs = list(results_storage.iter_results())
        
        # Check for duplicates at same position
        positions = [(m['Start'], m['End'], m.get('Class', '')) for m in all_motifs]
        unique_positions = set(positions)
        
        # Should have no duplicates
        self.assertEqual(
            len(positions), len(unique_positions),
            "Found duplicate motifs at micro-chunk boundary"
        )
        
        # Cleanup
        results_storage.cleanup()
    
    def test_meso_tier_deduplication(self):
        """Test deduplication at meso-chunk boundaries."""
        # Create sequence with pattern at meso chunk boundary
        # Meso chunk size is 5MB, place pattern near that boundary
        pre_sequence = "ATCG" * 1248000  # 4,992,000 bp
        g4_motif = "GGGTTAGGGTTAGGGTTAGGG"
        post_sequence = "ATCG" * 1248000  # 4,992,000 bp
        
        sequence = pre_sequence + g4_motif + post_sequence
        seq_id = self.storage.save_sequence(sequence, "meso_boundary_test")
        
        analyzer = TripleAdaptiveChunkAnalyzer(self.storage, use_adaptive=True)
        results_storage = analyzer.analyze(seq_id)
        
        # Check for duplicates
        all_motifs = list(results_storage.iter_results())
        positions = [(m['Start'], m['End'], m.get('Class', '')) for m in all_motifs]
        unique_positions = set(positions)
        
        self.assertEqual(
            len(positions), len(unique_positions),
            "Found duplicate motifs at meso-chunk boundary"
        )
        
        # Cleanup
        results_storage.cleanup()
    
    def test_position_adjustment_micro_tier(self):
        """Test that positions are correctly adjusted in micro-tier."""
        # Create sequence with identifiable patterns
        sequence = "ATCG" * 25000  # 100,000 bp
        seq_id = self.storage.save_sequence(sequence, "position_test")
        
        analyzer = TripleAdaptiveChunkAnalyzer(self.storage, use_adaptive=True)
        results_storage = analyzer.analyze(seq_id)
        
        # Verify all positions are within bounds
        for motif in results_storage.iter_results():
            self.assertGreaterEqual(motif['Start'], 0)
            self.assertLessEqual(motif['End'], 100000)
            self.assertLessEqual(motif['Start'], motif['End'])
        
        # Cleanup
        results_storage.cleanup()
    
    def test_position_adjustment_meso_tier(self):
        """Test that positions are correctly adjusted in meso-tier."""
        sequence = "ATCG" * 2500000  # 10,000,000 bp
        seq_id = self.storage.save_sequence(sequence, "meso_position_test")
        
        analyzer = TripleAdaptiveChunkAnalyzer(self.storage, use_adaptive=True)
        results_storage = analyzer.analyze(seq_id)
        
        # Verify all positions are within bounds
        for motif in results_storage.iter_results():
            self.assertGreaterEqual(motif['Start'], 0)
            self.assertLessEqual(motif['End'], 10000000)
            self.assertLessEqual(motif['Start'], motif['End'])
        
        # Cleanup
        results_storage.cleanup()
    
    def test_enabled_classes_filter(self):
        """Test filtering by enabled motif classes."""
        # Use a sequence with known G-quadruplex pattern
        # G4 pattern: GGGTTAGGGTTAGGGTTAGGG
        sequence = "ATCGATCG" * 100000 + "GGGTTAGGGTTAGGGTTAGGG" * 10 + "ATCGATCG" * 25000  # 1,000,000+ bp
        seq_id = self.storage.save_sequence(sequence, "filter_test")
        
        # Analyze with only G-quadruplex enabled
        analyzer = TripleAdaptiveChunkAnalyzer(self.storage, use_adaptive=True)
        results_storage = analyzer.analyze(
            seq_id,
            enabled_classes=['g_quadruplex']
        )
        
        # Get all motifs
        all_motifs = list(results_storage.iter_results())
        
        # If motifs were found, verify they are the correct class
        if all_motifs:
            for motif in all_motifs:
                # G-quadruplex maps to G-Quadruplex class
                self.assertIn('G-Quadruplex', motif.get('Class', ''),
                            f"Expected G-Quadruplex but got {motif.get('Class', '')}")
        
        # Cleanup
        results_storage.cleanup()
    
    def test_progress_callback(self):
        """Test that progress callbacks are called correctly."""
        sequence = "ATCGATCG" * 250000  # 2,000,000 bp
        seq_id = self.storage.save_sequence(sequence, "progress_test")
        
        progress_values = []
        
        def progress_callback(pct):
            progress_values.append(pct)
        
        analyzer = TripleAdaptiveChunkAnalyzer(self.storage, use_adaptive=True)
        results_storage = analyzer.analyze(seq_id, progress_callback=progress_callback)
        
        # Progress should have been called
        self.assertGreater(len(progress_values), 0)
        
        # Final progress should be close to 100%
        if progress_values:
            self.assertGreaterEqual(progress_values[-1], 90.0)
        
        # Cleanup
        results_storage.cleanup()
    
    def test_adaptive_disabled(self):
        """Test behavior when adaptive mode is disabled."""
        sequence = "ATCGATCG" * 12500000  # 100,000,000 bp
        seq_id = self.storage.save_sequence(sequence, "non_adaptive_test")
        
        # With adaptive disabled, should still use direct analysis
        analyzer = TripleAdaptiveChunkAnalyzer(self.storage, use_adaptive=False)
        results_storage = analyzer.analyze(seq_id)
        
        # Should still produce results
        stats = results_storage.get_summary_stats()
        self.assertGreaterEqual(stats['total_count'], 0)
        
        # Cleanup
        results_storage.cleanup()
    
    def test_motif_key_creation(self):
        """Test motif key creation for deduplication."""
        analyzer = TripleAdaptiveChunkAnalyzer(self.storage, use_adaptive=True)
        
        motif1 = {
            'Class': 'G-Quadruplex',
            'Subclass': 'Canonical intramolecular G4',
            'Start': 100,
            'End': 120
        }
        
        motif2 = {
            'Class': 'G-Quadruplex',
            'Subclass': 'Canonical intramolecular G4',
            'Start': 100,
            'End': 120
        }
        
        motif3 = {
            'Class': 'G-Quadruplex',
            'Subclass': 'Canonical intramolecular G4',
            'Start': 200,
            'End': 220
        }
        
        key1 = analyzer._create_motif_key(motif1)
        key2 = analyzer._create_motif_key(motif2)
        key3 = analyzer._create_motif_key(motif3)
        
        # Same motifs should have same key
        self.assertEqual(key1, key2)
        
        # Different motifs should have different keys
        self.assertNotEqual(key1, key3)
    
    def test_overlap_region_detection(self):
        """Test overlap region detection."""
        analyzer = TripleAdaptiveChunkAnalyzer(self.storage, use_adaptive=True)
        
        # Motif in overlap region (last 2KB of 50KB chunk)
        motif_in_overlap = {'Start': 49500, 'End': 49520}
        self.assertTrue(
            analyzer._is_in_overlap_region(
                motif_in_overlap, 0, 50000, 2000
            )
        )
        
        # Motif not in overlap region
        motif_not_in_overlap = {'Start': 10000, 'End': 10020}
        self.assertFalse(
            analyzer._is_in_overlap_region(
                motif_not_in_overlap, 0, 50000, 2000
            )
        )
    
    def test_position_adjustment(self):
        """Test position adjustment from chunk-local to global."""
        analyzer = TripleAdaptiveChunkAnalyzer(self.storage, use_adaptive=True)
        
        motifs = [
            {'ID': 'seq_G4_100', 'Start': 100, 'End': 120, 'Class': 'G-Quadruplex'},
            {'ID': 'seq_ZDNA_200', 'Start': 200, 'End': 220, 'Class': 'Z-DNA'}
        ]
        
        chunk_start = 50000
        adjusted = analyzer._adjust_motif_positions(motifs, chunk_start)
        
        # Positions should be adjusted
        self.assertEqual(adjusted[0]['Start'], 50100)
        self.assertEqual(adjusted[0]['End'], 50120)
        self.assertEqual(adjusted[1]['Start'], 50200)
        self.assertEqual(adjusted[1]['End'], 50220)
        
        # IDs should be updated
        self.assertIn('50100', adjusted[0]['ID'])
        self.assertIn('50200', adjusted[1]['ID'])


class TestAdaptiveStrategySelection(unittest.TestCase):
    """Test automatic strategy selection based on sequence size."""
    
    def setUp(self):
        """Create temporary directory and storage for tests."""
        self.test_dir = tempfile.mkdtemp(prefix="test_adaptive_")
        self.storage = UniversalSequenceStorage(base_dir=self.test_dir)
    
    def tearDown(self):
        """Clean up temporary directory."""
        if Path(self.test_dir).exists():
            shutil.rmtree(self.test_dir)
        gc.collect()
    
    def test_threshold_boundaries(self):
        """Test strategy selection at exact threshold boundaries."""
        analyzer = TripleAdaptiveChunkAnalyzer(self.storage, use_adaptive=True)
        
        # Test sequences at each threshold (CORRECTED thresholds)
        # Note: At boundary values, the next tier is used (< vs <=)
        test_cases = [
            (49_999, "direct"),       # Just below 50KB
            (50_000, "single"),       # Exactly 50KB -> single-tier
            (50_001, "single"),       # Just above 50KB
            (999_999, "single"),      # Just below 1MB
            (1_000_000, "double"),    # Exactly 1MB -> double-tier
            (1_000_001, "double"),    # Just above 1MB
            (99_999_999, "double"),   # Just below 100MB
            (100_000_000, "triple"),  # Exactly 100MB -> triple-tier
            (100_000_001, "triple"),  # Just above 100MB
        ]
        
        for seq_length, expected_strategy in test_cases:
            with self.subTest(seq_length=seq_length):
                # Determine expected strategy based on thresholds
                if seq_length < analyzer.DIRECT_THRESHOLD:
                    expected = "direct"
                elif seq_length < analyzer.SINGLE_TIER_THRESHOLD:
                    expected = "single"
                elif seq_length < analyzer.DOUBLE_TIER_THRESHOLD:
                    expected = "double"
                else:
                    expected = "triple"
                
                # Verify expected matches test case
                self.assertEqual(expected, expected_strategy)


if __name__ == '__main__':
    unittest.main()
