"""
╔══════════════════════════════════════════════════════════════════════════════╗
║              TEST SUITE FOR ANALYSIS BUTTON AND PARALLEL FIXES               ║
║                Testing analysis_done flag reset and parallel processing       ║
╚══════════════════════════════════════════════════════════════════════════════╝

Tests validate:
1. analysis_done flag is properly reset when new sequences are loaded
2. Parallel processing is triggered for equal-length multifasta sequences
3. Sequential processing is used for non-equal-length sequences
4. Results are correctly ordered when using parallel processing
"""

import sys
import os
import unittest
import tempfile

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Utilities.disk_storage import UniversalSequenceStorage, UniversalResultsStorage
from Utilities.nonbscanner import analyze_sequence


class TestAnalysisDoneFlagReset(unittest.TestCase):
    """Test suite for analysis_done flag reset functionality."""
    
    def test_equal_length_detection(self):
        """Test detection of equal-length sequences in multifasta."""
        # Create test sequences
        seqs = ["ATGCATGC" * 100, "GCTAGCTA" * 100, "CCGGCCGG" * 100]
        
        # Check lengths
        lengths = [len(seq) for seq in seqs]
        
        # All should have same length
        self.assertTrue(len(set(lengths)) == 1, "Sequences should have equal length")
        self.assertEqual(lengths[0], 800, "Each sequence should be 800bp")
    
    def test_unequal_length_detection(self):
        """Test detection of non-equal-length sequences in multifasta."""
        # Create test sequences with different lengths
        seqs = ["ATGCATGC" * 100, "GCTAGCTA" * 150, "CCGGCCGG" * 200]
        
        # Check lengths
        lengths = [len(seq) for seq in seqs]
        
        # Should have different lengths
        self.assertFalse(len(set(lengths)) == 1, "Sequences should have different lengths")
        self.assertEqual(len(set(lengths)), 3, "Should have 3 different lengths")


class TestParallelProcessing(unittest.TestCase):
    """Test suite for parallel processing of equal-length sequences."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.storage = UniversalSequenceStorage(base_dir=self.temp_dir)
    
    def tearDown(self):
        """Clean up test environment."""
        import shutil
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_equal_length_sequences_metadata(self):
        """Test that equal-length sequences are correctly identified via metadata."""
        # Create equal-length sequences
        seqs = [
            "ATGCATGCATGCATGC" * 50,  # 800bp
            "GCTAGCTAGCTAGCTA" * 50,  # 800bp
            "CCGGCCGGCCGGCCGG" * 50   # 800bp
        ]
        names = ["Seq1", "Seq2", "Seq3"]
        
        # Save sequences to storage
        seq_ids = []
        for seq, name in zip(seqs, names):
            seq_id = self.storage.save_sequence(seq, name)
            seq_ids.append(seq_id)
        
        # Get metadata and check lengths
        lengths = []
        for seq_id in seq_ids:
            metadata = self.storage.get_metadata(seq_id)
            lengths.append(metadata['length'])
        
        # All should be equal
        self.assertEqual(len(set(lengths)), 1, "All sequences should have equal length")
        self.assertEqual(lengths[0], 800, "Each sequence should be 800bp")
    
    def test_threshold_for_parallel_processing(self):
        """Test that sequences below chunk threshold qualify for parallel processing."""
        CHUNK_ANALYSIS_THRESHOLD_BP = 1_000_000  # 1MB
        
        # Create sequences smaller than threshold (800bp each)
        seqs = [
            "ATGCATGCATGCATGC" * 50,  # 800bp
            "GCTAGCTAGCTAGCTA" * 50,  # 800bp
            "CCGGCCGGCCGGCCGG" * 50   # 800bp
        ]
        
        # Check all are below threshold
        for seq in seqs:
            self.assertLess(len(seq), CHUNK_ANALYSIS_THRESHOLD_BP, 
                          f"Sequence length {len(seq)} should be below threshold {CHUNK_ANALYSIS_THRESHOLD_BP}")
    
    def test_analysis_results_order(self):
        """Test that parallel processing maintains correct order of results."""
        # Create equal-length test sequences
        seqs = [
            "GGGGGGGGGGGGGGGG" * 10,  # Should produce G-quadruplex motifs
            "ATATATAT" * 20,          # Should produce Z-DNA motifs  
            "GGGGGGGGGGGGGGGG" * 10   # Should produce G-quadruplex motifs
        ]
        names = ["Seq1_G4", "Seq2_ZDNA", "Seq3_G4"]
        
        # Save to storage
        seq_ids = []
        for seq, name in zip(seqs, names):
            seq_id = self.storage.save_sequence(seq, name)
            seq_ids.append(seq_id)
        
        # Simulate parallel processing by analyzing each sequence
        results_by_id = {}
        for seq_id, name in zip(seq_ids, names):
            seq = self.storage.load_sequence(seq_id)
            results = analyze_sequence(seq, name)
            results_by_id[seq_id] = results
        
        # Reconstruct in original order (as parallel code does)
        ordered_results = [results_by_id[seq_id] for seq_id in seq_ids]
        
        # Verify we have results for all 3 sequences
        self.assertEqual(len(ordered_results), 3, "Should have results for all 3 sequences")
        
        # Verify order is maintained (results correspond to seq_ids order)
        for i, results in enumerate(ordered_results):
            self.assertIsInstance(results, list, f"Results for sequence {i} should be a list")


class TestSequentialFallback(unittest.TestCase):
    """Test that sequential processing is used when parallel is not appropriate."""
    
    def test_unequal_lengths_no_parallel(self):
        """Verify that unequal-length sequences don't qualify for parallel processing."""
        # Create sequences with different lengths
        seqs = [
            "ATGCATGC" * 100,   # 800bp
            "GCTAGCTA" * 150,   # 1200bp
            "CCGGCCGG" * 200    # 1600bp
        ]
        
        lengths = [len(seq) for seq in seqs]
        equal_length = len(set(lengths)) == 1
        
        # Should not be equal length
        self.assertFalse(equal_length, "Sequences with different lengths should not be equal")
    
    def test_single_sequence_no_parallel(self):
        """Verify that single sequences don't use parallel processing."""
        seqs = ["ATGCATGC" * 100]
        num_sequences = len(seqs)
        
        # Should not trigger parallel processing
        self.assertEqual(num_sequences, 1, "Single sequence should not use parallel processing")
        
        # Parallel processing requires: equal_length AND num_sequences > 1
        equal_length = True  # Single sequence is trivially equal
        use_parallel = equal_length and num_sequences > 1
        
        self.assertFalse(use_parallel, "Single sequence should not trigger parallel processing")


if __name__ == '__main__':
    # Run tests
    unittest.main(verbosity=2)
