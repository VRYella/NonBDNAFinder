"""
╔══════════════════════════════════════════════════════════════════════════════╗
║            SEQUENCE PREPROCESSOR TEST SUITE                                   ║
║        Testing Gold Standard Genomic Sequence Preprocessing                  ║
╚══════════════════════════════════════════════════════════════════════════════╝

Tests validate:
1. GC calculation excludes ambiguous bases (gold standard)
2. FASTA header extraction and removal
3. Invalid character detection with position reporting
4. Multi-line sequence concatenation
5. GC balance classification
6. Comprehensive character counting
7. Validation status determination
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
from Utilities.sequence_preprocessor import (
    preprocess_sequence,
    PreprocessingResult,
    format_preprocessing_report
)


class TestPreprocessingBasics(unittest.TestCase):
    """Test basic preprocessing functionality."""
    
    def test_simple_sequence(self):
        """Test preprocessing of simple DNA sequence."""
        seq = "ATCG"
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.sequence, "ATCG")
        self.assertEqual(result.length, 4)
        self.assertEqual(result.valid_bases, 4)
        self.assertEqual(result.validation_status, "valid")
        self.assertEqual(len(result.errors), 0)
        self.assertEqual(len(result.warnings), 0)
    
    def test_lowercase_normalization(self):
        """Test that lowercase sequences are normalized to uppercase."""
        seq = "atcg"
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.sequence, "ATCG")
        self.assertEqual(result.validation_status, "valid")
    
    def test_mixed_case(self):
        """Test mixed case normalization."""
        seq = "AtCgAtCg"
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.sequence, "ATCGATCG")
        self.assertEqual(result.length, 8)
    
    def test_empty_sequence(self):
        """Test that empty sequence is handled properly."""
        seq = ""
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.validation_status, "error")
        self.assertTrue(any("Empty sequence" in err for err in result.errors))


class TestGCCalculation(unittest.TestCase):
    """Test GC% calculation with correct formula."""
    
    def test_gc_calculation_excludes_ambiguous_bases(self):
        """
        Test GC% uses (A+T+G+C) denominator, not total length.
        This is the gold standard: (G+C) / (A+T+G+C) × 100
        """
        seq = "ATCGNNNN"
        result = preprocess_sequence(seq)
        
        # Should be 50% (2 GC / 4 ATGC), NOT 25% (2 GC / 8 total)
        self.assertAlmostEqual(result.gc_percentage, 50.0, places=2)
        self.assertEqual(result.valid_bases, 4)
        self.assertEqual(result.length, 8)
        self.assertEqual(result.validation_status, "warning")  # N's present
    
    def test_gc_at_sum_to_100(self):
        """Test that GC% + AT% = 100%."""
        test_cases = [
            "ATCG",
            "GGGGCCCC",
            "AAAATTTT",
            "ATCGNNNN",
        ]
        
        for seq in test_cases:
            with self.subTest(seq=seq):
                result = preprocess_sequence(seq)
                total = result.gc_percentage + result.at_percentage
                self.assertAlmostEqual(total, 100.0, places=2,
                                     msg=f"GC% + AT% should equal 100% for {seq}")
    
    def test_100_percent_gc(self):
        """Test sequence with 100% GC content."""
        seq = "GGGGCCCC"
        result = preprocess_sequence(seq)
        
        self.assertAlmostEqual(result.gc_percentage, 100.0, places=2)
        self.assertAlmostEqual(result.at_percentage, 0.0, places=2)
    
    def test_0_percent_gc(self):
        """Test sequence with 0% GC content."""
        seq = "AAAATTTT"
        result = preprocess_sequence(seq)
        
        self.assertAlmostEqual(result.gc_percentage, 0.0, places=2)
        self.assertAlmostEqual(result.at_percentage, 100.0, places=2)
    
    def test_only_n_bases(self):
        """Test sequence with only N bases."""
        seq = "NNNN"
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.validation_status, "error")
        self.assertTrue(any("No valid ATGC bases" in err for err in result.errors))


class TestFASTAHandling(unittest.TestCase):
    """Test FASTA format handling."""
    
    def test_fasta_header_removal(self):
        """Test that FASTA header is extracted and removed from sequence."""
        fasta = ">header\nATCG\nATCG"
        result = preprocess_sequence(fasta)
        
        self.assertEqual(result.header, "HEADER")  # Headers are uppercased
        self.assertEqual(result.sequence, "ATCGATCG")
        self.assertNotIn('>', result.sequence)
    
    def test_multiline_concatenation(self):
        """Test that multi-line FASTA sequences are concatenated correctly."""
        fasta = ">seq\nATCG\nATCG\nATCG"
        result = preprocess_sequence(fasta)
        
        self.assertEqual(result.sequence, "ATCGATCGATCG")
        self.assertEqual(result.length, 12)
        self.assertNotIn('\n', result.sequence)
    
    def test_fasta_with_description(self):
        """Test FASTA header with description."""
        fasta = ">seq1 chromosome 1 complete sequence\nATCGATCG"
        result = preprocess_sequence(fasta)
        
        self.assertEqual(result.header, "SEQ1 CHROMOSOME 1 COMPLETE SEQUENCE")  # Headers are uppercased
        self.assertEqual(result.sequence, "ATCGATCG")
    
    def test_no_header(self):
        """Test sequence without FASTA header."""
        seq = "ATCGATCG"
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.header, "")
        self.assertEqual(result.sequence, "ATCGATCG")
    
    def test_multiple_headers(self):
        """Test that only first header is kept."""
        fasta = ">header1\nATCG\n>header2\nATCG"
        result = preprocess_sequence(fasta)
        
        self.assertEqual(result.header, "HEADER1")  # Headers are uppercased
        # Note: Second header line will be ignored as a sequence line starting with '>'
        # This is intentional - multi-FASTA should be handled by other tools


class TestInvalidCharacterDetection(unittest.TestCase):
    """Test detection of invalid characters."""
    
    def test_invalid_character_detection(self):
        """Test that invalid characters are detected."""
        seq = "ATCG123XZ"  # Using X and Z which are not valid IUPAC codes
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.validation_status, "error")
        self.assertIn('1', result.invalid_characters)
        self.assertIn('2', result.invalid_characters)
        self.assertIn('3', result.invalid_characters)
        self.assertIn('X', result.invalid_characters)
        self.assertIn('Z', result.invalid_characters)
    
    def test_invalid_char_position_reporting(self):
        """Test that positions of invalid characters are reported."""
        seq = "AT1CG2"
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.validation_status, "error")
        self.assertEqual(result.invalid_characters['1'], [2])
        self.assertEqual(result.invalid_characters['2'], [5])
    
    def test_max_10_positions_per_char(self):
        """Test that at most 10 positions per invalid char are stored."""
        seq = "A" + "X" * 20 + "T"  # 20 X's
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.validation_status, "error")
        self.assertEqual(len(result.invalid_characters['X']), 10)


class TestCharacterCounting(unittest.TestCase):
    """Test comprehensive character counting."""
    
    def test_character_counts(self):
        """Test that character counts are accurate."""
        seq = "AAATTTGGGCCC"
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.character_counts['A'], 3)
        self.assertEqual(result.character_counts['T'], 3)
        self.assertEqual(result.character_counts['G'], 3)
        self.assertEqual(result.character_counts['C'], 3)
        self.assertEqual(result.character_counts['N'], 0)
    
    def test_n_counting(self):
        """Test that N bases are counted separately."""
        seq = "ATCGNNN"
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.character_counts['N'], 3)
        self.assertEqual(result.valid_bases, 4)
        self.assertEqual(result.length, 7)
    
    def test_valid_bases_excludes_n(self):
        """Test that valid_bases count excludes N."""
        seq = "AAATTTGGGCCCNNN"
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.valid_bases, 12)
        self.assertEqual(result.character_counts['N'], 3)


class TestGCBalanceClassification(unittest.TestCase):
    """Test GC balance classification."""
    
    def test_gc_rich_classification(self):
        """Test that GC-rich sequences are classified correctly."""
        seq = "GGGGCCCC"  # 100% GC
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.gc_balance, "GC-rich")
    
    def test_at_rich_classification(self):
        """Test that AT-rich sequences are classified correctly."""
        seq = "AAAATTTT"  # 0% GC
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.gc_balance, "AT-rich")
    
    def test_balanced_classification(self):
        """Test that balanced sequences are classified correctly."""
        seq = "ATCGATCG"  # 50% GC
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.gc_balance, "Balanced")
    
    def test_gc_balance_thresholds(self):
        """Test GC balance classification thresholds."""
        # GC% > 60% = GC-rich
        seq_gc_rich = "GGGGGGCCCCCCAATT"  # 75% GC
        result = preprocess_sequence(seq_gc_rich)
        self.assertEqual(result.gc_balance, "GC-rich")
        
        # GC% < 40% = AT-rich
        seq_at_rich = "AAAAAAATTTTTGGCC"  # 25% GC
        result = preprocess_sequence(seq_at_rich)
        self.assertEqual(result.gc_balance, "AT-rich")
        
        # 40% <= GC% <= 60% = Balanced
        seq_balanced = "AATTGGCC"  # 50% GC
        result = preprocess_sequence(seq_balanced)
        self.assertEqual(result.gc_balance, "Balanced")


class TestValidationStatus(unittest.TestCase):
    """Test validation status determination."""
    
    def test_valid_status(self):
        """Test that clean sequences get 'valid' status."""
        seq = "ATCGATCG"
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.validation_status, "valid")
        self.assertEqual(len(result.warnings), 0)
        self.assertEqual(len(result.errors), 0)
    
    def test_warning_status_with_n(self):
        """Test that sequences with N get 'warning' status."""
        seq = "ATCGNNNN"
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.validation_status, "warning")
        self.assertTrue(len(result.warnings) > 0)
        self.assertTrue(any("ambiguous" in w.lower() for w in result.warnings))
    
    def test_error_status_with_invalid_chars(self):
        """Test that sequences with invalid chars get 'error' status."""
        seq = "ATCG123"
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.validation_status, "error")
        self.assertTrue(len(result.errors) > 0)


class TestFormattingReport(unittest.TestCase):
    """Test report formatting."""
    
    def test_format_report(self):
        """Test that formatting function produces output."""
        seq = "ATCGATCG"
        result = preprocess_sequence(seq)
        report = format_preprocessing_report(result)
        
        self.assertIsInstance(report, str)
        self.assertIn("PREPROCESSING REPORT", report)
        self.assertIn("Length:", report)
        self.assertIn("GC Content:", report)
    
    def test_report_includes_warnings(self):
        """Test that report includes warnings when present."""
        seq = "ATCGNNNN"
        result = preprocess_sequence(seq)
        report = format_preprocessing_report(result)
        
        self.assertIn("WARNINGS", report)
        self.assertIn("ambiguous", report.lower())
    
    def test_report_includes_errors(self):
        """Test that report includes errors when present."""
        seq = "ATCG123"
        result = preprocess_sequence(seq)
        report = format_preprocessing_report(result)
        
        self.assertIn("ERRORS", report)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and boundary conditions."""
    
    def test_very_long_sequence(self):
        """Test that long sequences are handled correctly."""
        seq = "ATCG" * 10000  # 40,000 bp
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.length, 40000)
        self.assertEqual(result.valid_bases, 40000)
        self.assertAlmostEqual(result.gc_percentage, 50.0, places=2)
    
    def test_whitespace_handling(self):
        """Test that whitespace is handled correctly."""
        seq = "  ATCG  \n  ATCG  "
        result = preprocess_sequence(seq)
        
        self.assertEqual(result.sequence, "ATCGATCG")
    
    def test_various_n_positions(self):
        """Test N's at various positions."""
        test_cases = [
            ("NNNATCG", 50.0),   # N's at start
            ("ATCGNNN", 50.0),   # N's at end
            ("ATNNCG", 50.0),    # N's in middle
            ("NANTNCGN", 50.0),  # N's scattered
        ]
        
        for seq, expected_gc in test_cases:
            with self.subTest(seq=seq):
                result = preprocess_sequence(seq)
                self.assertAlmostEqual(result.gc_percentage, expected_gc, places=2)
                self.assertEqual(result.validation_status, "warning")


if __name__ == '__main__':
    unittest.main()
