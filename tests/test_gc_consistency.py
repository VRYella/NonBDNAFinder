"""
╔══════════════════════════════════════════════════════════════════════════════╗
║            GC CONTENT CALCULATION CONSISTENCY TEST                            ║
║        Verify all GC% calculations use the same method                       ║
╚══════════════════════════════════════════════════════════════════════════════╝

Tests validate:
1. All GC content functions produce identical results
2. Consistency across utilities.py, disk_storage.py, and detectors_utils.py
3. Proper handling of edge cases (empty, mixed case, all GC, no GC)
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
import tempfile
import shutil

from Utilities.detectors_utils import calc_gc_content
from Utilities.utilities import gc_content
from Utilities.disk_storage import UniversalSequenceStorage


class TestGCConsistency(unittest.TestCase):
    """Test that all GC content calculation methods produce identical results."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create temporary directory for disk storage tests
        self.temp_dir = tempfile.mkdtemp()
        self.storage = UniversalSequenceStorage(self.temp_dir)
    
    def tearDown(self):
        """Clean up test fixtures."""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_empty_sequence(self):
        """Test that empty sequence returns 0.0 GC%."""
        seq = ""
        result1 = calc_gc_content(seq)
        result2 = gc_content(seq)
        
        self.assertEqual(result1, 0.0, "calc_gc_content should return 0.0 for empty sequence")
        self.assertEqual(result2, 0.0, "gc_content should return 0.0 for empty sequence")
        self.assertEqual(result1, result2, "Both methods should return identical results")
    
    def test_50_percent_gc(self):
        """Test sequence with 50% GC content."""
        seq = "ATCG"
        result1 = calc_gc_content(seq)
        result2 = gc_content(seq)
        
        self.assertEqual(result1, 50.0, "calc_gc_content should return 50.0")
        self.assertEqual(result2, 50.0, "gc_content should return 50.0")
        self.assertEqual(result1, result2, "Both methods should return identical results")
    
    def test_100_percent_gc(self):
        """Test sequence with 100% GC content."""
        seq = "GGCC"
        result1 = calc_gc_content(seq)
        result2 = gc_content(seq)
        
        self.assertEqual(result1, 100.0, "calc_gc_content should return 100.0")
        self.assertEqual(result2, 100.0, "gc_content should return 100.0")
        self.assertEqual(result1, result2, "Both methods should return identical results")
    
    def test_0_percent_gc(self):
        """Test sequence with 0% GC content."""
        seq = "AAATTT"
        result1 = calc_gc_content(seq)
        result2 = gc_content(seq)
        
        self.assertEqual(result1, 0.0, "calc_gc_content should return 0.0")
        self.assertEqual(result2, 0.0, "gc_content should return 0.0")
        self.assertEqual(result1, result2, "Both methods should return identical results")
    
    def test_mixed_case(self):
        """Test that mixed case sequences are handled correctly."""
        seq = "ATATatgc"
        result1 = calc_gc_content(seq)
        result2 = gc_content(seq)
        
        self.assertEqual(result1, 25.0, "calc_gc_content should return 25.0 for mixed case")
        self.assertEqual(result2, 25.0, "gc_content should return 25.0 for mixed case")
        self.assertEqual(result1, result2, "Both methods should return identical results")
    
    def test_real_sequence(self):
        """Test with a realistic DNA sequence."""
        seq = "ATCGTAGCTGCAGTCGATCGTAGCTAGCTAGCTAGCTAGCTAG"
        result1 = calc_gc_content(seq)
        result2 = gc_content(seq)
        
        self.assertIsInstance(result1, float, "calc_gc_content should return a float")
        self.assertIsInstance(result2, float, "gc_content should return a float")
        self.assertEqual(result1, result2, "Both methods should return identical results")
        self.assertGreaterEqual(result1, 0.0, "GC% should be >= 0")
        self.assertLessEqual(result1, 100.0, "GC% should be <= 100")
    
    def test_disk_storage_consistency(self):
        """Test that disk storage uses the same GC calculation."""
        test_seq = "ATCGTAGCTGCAGTCGATCG"
        expected_gc = calc_gc_content(test_seq)
        
        # Save sequence to disk storage
        seq_id = self.storage.save_sequence(test_seq, "test_seq")
        
        # Retrieve metadata
        metadata = self.storage.get_metadata(seq_id)
        
        # Check that GC content in metadata matches
        self.assertAlmostEqual(
            metadata['gc_content'], 
            expected_gc, 
            places=5,
            msg="Disk storage GC% should match calc_gc_content()"
        )
    
    def test_all_three_methods_identical(self):
        """Test that all three GC calculation methods produce identical results."""
        test_cases = [
            ("ATCG", 50.0),
            ("GGCC", 100.0),
            ("AAATTT", 0.0),
            ("ATATatgc", 25.0),
            ("GCGCGCGC", 100.0),
            ("", 0.0),
        ]
        
        for seq, expected in test_cases:
            with self.subTest(seq=seq):
                result1 = calc_gc_content(seq)
                result2 = gc_content(seq)
                
                # For disk storage, we need a non-empty sequence
                if seq:
                    seq_id = self.storage.save_sequence(seq, f"test_{seq[:5]}")
                    metadata = self.storage.get_metadata(seq_id)
                    result3 = metadata['gc_content']
                    
                    self.assertEqual(result1, result2, 
                                   f"calc_gc_content and gc_content differ for '{seq}'")
                    self.assertAlmostEqual(result1, result3, places=5,
                                         msg=f"disk_storage GC differs for '{seq}'")
                else:
                    # Empty sequence - only test the two functions
                    self.assertEqual(result1, result2, 
                                   f"calc_gc_content and gc_content differ for empty sequence")
                
                self.assertAlmostEqual(result1, expected, places=5,
                                     msg=f"Expected {expected}% for '{seq}', got {result1}%")


if __name__ == '__main__':
    unittest.main()
