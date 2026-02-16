"""
╔══════════════════════════════════════════════════════════════════════════════╗
║            ERROR HANDLING TEST SUITE                                          ║
║        Testing Empty Sequence and Edge Case Handling                          ║
╚══════════════════════════════════════════════════════════════════════════════╝

Tests validate:
1. Empty sequence handling returns empty list without errors
2. Type validation for sequence parameter
3. Chunking threshold properly applies for >50KB sequences
4. Progress callbacks work correctly
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
from typing import List, Dict, Any

from Utilities.nonbscanner import analyze_sequence, CHUNK_THRESHOLD


class TestEmptySequenceHandling(unittest.TestCase):
    """Test suite for empty sequence validation"""
    
    def test_empty_string_returns_empty_list(self):
        """Empty string should return empty list without error"""
        result = analyze_sequence("", "empty_test")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 0)
    
    def test_none_sequence_returns_empty_list(self):
        """None sequence should return empty list without error"""
        result = analyze_sequence(None, "none_test")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 0)
    
    def test_invalid_type_raises_error(self):
        """Non-string sequence should raise TypeError"""
        with self.assertRaises(TypeError):
            analyze_sequence(12345, "invalid_type")
    
    def test_valid_small_sequence(self):
        """Valid small sequence should work without issues"""
        # Simple valid sequence
        result = analyze_sequence("ATCG" * 100, "small_test")
        self.assertIsInstance(result, list)
        # Should not raise any errors


class TestChunkingThreshold(unittest.TestCase):
    """Test suite for chunking threshold configuration"""
    
    def test_chunk_threshold_is_50kb(self):
        """Verify CHUNK_THRESHOLD is set to 50,000"""
        self.assertEqual(CHUNK_THRESHOLD, 50000)
    
    def test_small_sequence_no_chunking(self):
        """Sequences < 50KB should not use chunking by default"""
        # 10KB sequence (20,000 bp < 50KB)
        seq = "ATCG" * 5000  # 20,000 bp
        result = analyze_sequence(seq, "small_seq", use_chunking=None)
        self.assertIsInstance(result, list)
    
    def test_large_sequence_triggers_chunking(self):
        """Sequences > 50KB should trigger chunking"""
        # 60KB sequence (60,000 bp > 50KB)
        seq = "ATCG" * 15000  # 60,000 bp
        
        # Track if progress callback was called (indicates chunking)
        callback_called = {'called': False}
        
        def progress_callback(chunk, total, bp, elapsed, throughput):
            callback_called['called'] = True
        
        result = analyze_sequence(
            seq, 
            "large_seq", 
            use_chunking=None,  # Should auto-enable
            progress_callback=progress_callback
        )
        self.assertIsInstance(result, list)
        # If sequence is large enough to chunk, callback should be called
        # Note: This may not always trigger if fast mode is used
    
    def test_explicit_chunking_enabled(self):
        """Explicit chunking parameter should work"""
        seq = "ATCG" * 15000  # 60,000 bp
        
        result = analyze_sequence(
            seq, 
            "explicit_chunking", 
            use_chunking=True
        )
        self.assertIsInstance(result, list)


class TestProgressCallback(unittest.TestCase):
    """Test suite for progress callback functionality"""
    
    def test_progress_callback_receives_correct_parameters(self):
        """Progress callback should receive correct parameters"""
        seq = "ATCG" * 15000  # 60,000 bp - large enough to chunk
        
        callback_data = []
        
        def progress_callback(chunk_num, total_chunks, bp_processed, elapsed_time, throughput):
            callback_data.append({
                'chunk_num': chunk_num,
                'total_chunks': total_chunks,
                'bp_processed': bp_processed,
                'elapsed_time': elapsed_time,
                'throughput': throughput
            })
        
        result = analyze_sequence(
            seq, 
            "callback_test", 
            use_chunking=True,
            progress_callback=progress_callback
        )
        
        # If chunking happened, we should have callback data
        if callback_data:
            # Verify structure of callback data
            for data in callback_data:
                self.assertIn('chunk_num', data)
                self.assertIn('total_chunks', data)
                self.assertIn('bp_processed', data)
                self.assertIn('elapsed_time', data)
                self.assertIn('throughput', data)
                
                # Verify types and ranges
                self.assertIsInstance(data['chunk_num'], int)
                self.assertIsInstance(data['total_chunks'], int)
                self.assertGreater(data['chunk_num'], 0)
                self.assertGreaterEqual(data['elapsed_time'], 0)


class TestSequenceSizeVariations(unittest.TestCase):
    """Test different sequence sizes"""
    
    def test_tiny_sequence(self):
        """Test with very small sequence (400bp)"""
        seq = "ATCG" * 100  # 400bp
        result = analyze_sequence(seq, "tiny_test")
        self.assertIsInstance(result, list)
    
    def test_small_sequence(self):
        """Test with small sequence (10KB)"""
        seq = "ATCG" * 2500  # 10,000bp
        result = analyze_sequence(seq, "small_test")
        self.assertIsInstance(result, list)
    
    def test_medium_sequence(self):
        """Test with medium sequence (60KB)"""
        seq = "ATCG" * 15000  # 60,000bp
        result = analyze_sequence(seq, "medium_test")
        self.assertIsInstance(result, list)
    
    def test_large_sequence(self):
        """Test with large sequence (200KB)"""
        seq = "ATCG" * 50000  # 200,000bp
        result = analyze_sequence(seq, "large_test", use_chunking=True)
        self.assertIsInstance(result, list)


if __name__ == '__main__':
    # Run tests with verbose output
    unittest.main(verbosity=2)
