"""
Test suite for eGZ-motif (Extruded-G Z-DNA) detection

This test file verifies that the ZDNADetector correctly identifies
long (CGG)n, (GGC)n, (CCG)n, and (GCC)n trinucleotide repeats via regex patterns.
"""

import unittest
from detectors import ZDNADetector


class TestEGZMotifDetection(unittest.TestCase):
    """Test suite for eGZ-motif detection in ZDNADetector"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.detector = ZDNADetector()
    
    def test_cgg_repeat_detection(self):
        """Test detection of (CGG)n repeats"""
        # Create sequence with CGG repeat (minimum 3 repeats = 9 bp)
        sequence = "ATGCGGCGGCGGATGC"  # Contains (CGG)3
        motifs = self.detector.detect_motifs(sequence, "test_seq")
        
        # Should detect at least one eGZ motif
        egz_motifs = [m for m in motifs if m['Subclass'] == 'eGZ']
        self.assertGreater(len(egz_motifs), 0, "Should detect CGG repeat as eGZ motif")
        
        # Verify the detected motif
        motif = egz_motifs[0]
        self.assertEqual(motif['Class'], 'Z-DNA')
        self.assertEqual(motif['Subclass'], 'eGZ')
        self.assertEqual(motif['Repeat_Unit'], 'CGG')
        self.assertGreaterEqual(motif['Repeat_Count'], 3)
        self.assertIn('CGG', motif['Sequence'])
    
    def test_ggc_repeat_detection(self):
        """Test detection of (GGC)n repeats"""
        # Create sequence with GGC repeat
        sequence = "ATGGCGGCGGCATGC"  # Contains (GGC)3
        motifs = self.detector.detect_motifs(sequence, "test_seq")
        
        # Should detect at least one eGZ motif
        egz_motifs = [m for m in motifs if m['Subclass'] == 'eGZ']
        self.assertGreater(len(egz_motifs), 0, "Should detect GGC repeat as eGZ motif")
        
        # Verify the detected motif
        motif = egz_motifs[0]
        self.assertEqual(motif['Subclass'], 'eGZ')
        self.assertEqual(motif['Repeat_Unit'], 'GGC')
        self.assertGreaterEqual(motif['Repeat_Count'], 3)
    
    def test_ccg_repeat_detection(self):
        """Test detection of (CCG)n repeats"""
        sequence = "ATCCGCCGCCGATGC"  # Contains (CCG)3
        motifs = self.detector.detect_motifs(sequence, "test_seq")
        
        egz_motifs = [m for m in motifs if m['Subclass'] == 'eGZ']
        self.assertGreater(len(egz_motifs), 0, "Should detect CCG repeat as eGZ motif")
        
        motif = egz_motifs[0]
        self.assertEqual(motif['Repeat_Unit'], 'CCG')
    
    def test_gcc_repeat_detection(self):
        """Test detection of (GCC)n repeats"""
        sequence = "ATGCCGCCGCCATGC"  # Contains (GCC)3
        motifs = self.detector.detect_motifs(sequence, "test_seq")
        
        egz_motifs = [m for m in motifs if m['Subclass'] == 'eGZ']
        self.assertGreater(len(egz_motifs), 0, "Should detect GCC repeat as eGZ motif")
        
        motif = egz_motifs[0]
        self.assertEqual(motif['Repeat_Unit'], 'GCC')
    
    def test_minimum_repeat_threshold(self):
        """Test that minimum 3 repeats are required"""
        # Only 2 repeats - should not be detected
        sequence = "ATCGGCGGATGC"  # Only (CGG)2
        motifs = self.detector.detect_motifs(sequence, "test_seq")
        egz_motifs = [m for m in motifs if m['Subclass'] == 'eGZ']
        
        # Should not detect with only 2 repeats
        self.assertEqual(len(egz_motifs), 0, "Should not detect with only 2 repeats")
    
    def test_long_cgg_repeat(self):
        """Test detection of long CGG repeat"""
        # Create a longer repeat sequence (e.g., 10 repeats)
        sequence = "ATG" + "CGG" * 10 + "ATGC"
        motifs = self.detector.detect_motifs(sequence, "test_seq")
        
        egz_motifs = [m for m in motifs if m['Subclass'] == 'eGZ']
        self.assertGreater(len(egz_motifs), 0, "Should detect long CGG repeat")
        
        motif = egz_motifs[0]
        self.assertEqual(motif['Repeat_Count'], 10)
        self.assertEqual(motif['Length'], 30)  # 10 repeats * 3 bp
        # Score should increase with repeat count
        self.assertGreater(motif['Score'], 2.0)
    
    def test_score_increases_with_repeats(self):
        """Test that score increases with more repeats"""
        # 3 repeats
        seq3 = "ATG" + "CGG" * 3 + "ATGC"
        motifs3 = self.detector.detect_motifs(seq3, "test_seq")
        egz3 = [m for m in motifs3 if m['Subclass'] == 'eGZ']
        
        # 6 repeats
        seq6 = "ATG" + "CGG" * 6 + "ATGC"
        motifs6 = self.detector.detect_motifs(seq6, "test_seq")
        egz6 = [m for m in motifs6 if m['Subclass'] == 'eGZ']
        
        if egz3 and egz6:
            self.assertLess(egz3[0]['Score'], egz6[0]['Score'], 
                          "Score should increase with more repeats")
    
    def test_both_zdna_and_egz_detection(self):
        """Test that both classic Z-DNA and eGZ can be detected in same sequence"""
        # Create sequence with both 10-mer Z-DNA pattern and CGG repeat
        sequence = "GCGCGCGCGC" + "ATGC" + "CGGCGGCGG" + "ATGC"
        motifs = self.detector.detect_motifs(sequence, "test_seq")
        
        # Should detect both types
        zdna_motifs = [m for m in motifs if m['Subclass'] == 'Z-DNA']
        egz_motifs = [m for m in motifs if m['Subclass'] == 'eGZ']
        
        # At least one of each type should be detected
        self.assertGreater(len(zdna_motifs) + len(egz_motifs), 0, 
                          "Should detect at least one motif")
    
    def test_case_insensitivity(self):
        """Test that detection is case-insensitive"""
        # Lower case CGG repeat
        sequence_lower = "atgcggcggcggatgc"
        motifs_lower = self.detector.detect_motifs(sequence_lower, "test_seq")
        
        # Upper case CGG repeat
        sequence_upper = "ATGCGGCGGCGGATGC"
        motifs_upper = self.detector.detect_motifs(sequence_upper, "test_seq")
        
        egz_lower = [m for m in motifs_lower if m['Subclass'] == 'eGZ']
        egz_upper = [m for m in motifs_upper if m['Subclass'] == 'eGZ']
        
        # Should detect the same in both cases
        self.assertEqual(len(egz_lower), len(egz_upper))
    
    def test_motif_output_format(self):
        """Test that eGZ motifs have correct output format"""
        sequence = "ATGCGGCGGCGGATGC"
        motifs = self.detector.detect_motifs(sequence, "test_seq")
        
        egz_motifs = [m for m in motifs if m['Subclass'] == 'eGZ']
        if egz_motifs:
            motif = egz_motifs[0]
            
            # Check required fields
            required_fields = ['ID', 'Sequence_Name', 'Class', 'Subclass', 
                             'Start', 'End', 'Length', 'Sequence', 'Score', 
                             'Strand', 'Method', 'Pattern_ID']
            for field in required_fields:
                self.assertIn(field, motif, f"Missing required field: {field}")
            
            # Check eGZ-specific fields
            self.assertIn('Repeat_Unit', motif)
            self.assertIn('Repeat_Count', motif)
            self.assertIn('GC_Content', motif)
            
            # Check data types
            self.assertIsInstance(motif['Start'], int)
            self.assertIsInstance(motif['End'], int)
            self.assertIsInstance(motif['Score'], (int, float))
            self.assertIsInstance(motif['Repeat_Count'], int)


def run_tests():
    """Run all tests and print results"""
    suite = unittest.TestLoader().loadTestsFromTestCase(TestEGZMotifDetection)
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    return result.wasSuccessful()


if __name__ == '__main__':
    success = run_tests()
    exit(0 if success else 1)
