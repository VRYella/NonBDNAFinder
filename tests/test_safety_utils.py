"""
╔══════════════════════════════════════════════════════════════════════════════╗
║            SAFETY UTILITIES TEST SUITE                                        ║
║        Testing Safe List/Array Access Utilities                              ║
╚══════════════════════════════════════════════════════════════════════════════╝

Tests validate:
1. Safe list access with bounds checking
2. Safe min/max operations on empty lists
3. Safe DataFrame access with bounds checking
4. Index filtering for valid ranges
5. Sequence validation and filtering
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
import pandas as pd
from Utilities.safety import (
    safe_list_access,
    safe_min_max,
    safe_dataframe_access,
    filter_valid_indices,
    validate_sequence_input
)


class TestSafeListAccess(unittest.TestCase):
    """Test suite for safe_list_access function"""
    
    def test_valid_index(self):
        """Test accessing valid index returns correct value"""
        lst = [10, 20, 30, 40, 50]
        self.assertEqual(safe_list_access(lst, 2), 30)
        self.assertEqual(safe_list_access(lst, 0), 10)
        self.assertEqual(safe_list_access(lst, 4), 50)
    
    def test_out_of_bounds_positive(self):
        """Test out of bounds positive index returns default"""
        lst = [1, 2, 3]
        self.assertEqual(safe_list_access(lst, 10, default=-1), -1)
        self.assertEqual(safe_list_access(lst, 3, default=0), 0)
    
    def test_out_of_bounds_negative(self):
        """Test negative index returns default (we only allow positive indices)"""
        lst = [1, 2, 3]
        self.assertEqual(safe_list_access(lst, -1, default=0), 0)
        self.assertEqual(safe_list_access(lst, -5, default=0), 0)
    
    def test_empty_list(self):
        """Test accessing empty list returns default"""
        lst = []
        self.assertEqual(safe_list_access(lst, 0, default=-1), -1)
    
    def test_non_list_input(self):
        """Test non-list input returns default"""
        self.assertEqual(safe_list_access("not a list", 0, default=None), None)
        self.assertEqual(safe_list_access(42, 0, default=None), None)
    
    def test_tuple_input(self):
        """Test tuple input works correctly"""
        tpl = (10, 20, 30)
        self.assertEqual(safe_list_access(tpl, 1), 20)
        self.assertEqual(safe_list_access(tpl, 5, default=-1), -1)


class TestSafeMinMax(unittest.TestCase):
    """Test suite for safe_min_max function"""
    
    def test_min_on_valid_list(self):
        """Test min operation on valid list"""
        lst = [5, 2, 8, 1, 9]
        self.assertEqual(safe_min_max(lst, 'min'), 1)
    
    def test_max_on_valid_list(self):
        """Test max operation on valid list"""
        lst = [5, 2, 8, 1, 9]
        self.assertEqual(safe_min_max(lst, 'max'), 9)
    
    def test_min_on_empty_list(self):
        """Test min on empty list returns default"""
        self.assertEqual(safe_min_max([], 'min', default=0), 0)
        self.assertEqual(safe_min_max([], 'min', default=100), 100)
    
    def test_max_on_empty_list(self):
        """Test max on empty list returns default"""
        self.assertEqual(safe_min_max([], 'max', default=0), 0)
        self.assertEqual(safe_min_max([], 'max', default=100), 100)
    
    def test_single_element(self):
        """Test min/max on single element list"""
        lst = [42]
        self.assertEqual(safe_min_max(lst, 'min'), 42)
        self.assertEqual(safe_min_max(lst, 'max'), 42)
    
    def test_negative_numbers(self):
        """Test min/max with negative numbers"""
        lst = [-5, -2, -8, -1]
        self.assertEqual(safe_min_max(lst, 'min'), -8)
        self.assertEqual(safe_min_max(lst, 'max'), -1)


class TestSafeDataFrameAccess(unittest.TestCase):
    """Test suite for safe_dataframe_access function"""
    
    def test_valid_index(self):
        """Test accessing valid DataFrame row"""
        df = pd.DataFrame({
            'A': [1, 2, 3],
            'B': [4, 5, 6]
        })
        row = safe_dataframe_access(df, 1)
        self.assertIsNotNone(row)
        self.assertEqual(row['A'], 2)
        self.assertEqual(row['B'], 5)
    
    def test_out_of_bounds(self):
        """Test out of bounds returns default"""
        df = pd.DataFrame({'A': [1, 2, 3]})
        self.assertIsNone(safe_dataframe_access(df, 10))
        self.assertIsNone(safe_dataframe_access(df, 3))
    
    def test_negative_index(self):
        """Test negative index returns default"""
        df = pd.DataFrame({'A': [1, 2, 3]})
        self.assertIsNone(safe_dataframe_access(df, -1))
    
    def test_empty_dataframe(self):
        """Test accessing empty DataFrame returns default"""
        df = pd.DataFrame()
        self.assertIsNone(safe_dataframe_access(df, 0))
    
    def test_non_dataframe_input(self):
        """Test non-DataFrame input returns default"""
        self.assertIsNone(safe_dataframe_access("not a dataframe", 0))
        self.assertIsNone(safe_dataframe_access([1, 2, 3], 0))


class TestFilterValidIndices(unittest.TestCase):
    """Test suite for filter_valid_indices function"""
    
    def test_all_valid_indices(self):
        """Test when all indices are valid"""
        indices = [0, 1, 2, 3, 4]
        result = filter_valid_indices(indices, max_length=10)
        self.assertEqual(result, [0, 1, 2, 3, 4])
    
    def test_some_invalid_indices(self):
        """Test filtering out invalid indices"""
        indices = [0, 5, 10, 15, 20]
        result = filter_valid_indices(indices, max_length=12)
        self.assertEqual(result, [0, 5, 10])
    
    def test_all_invalid_indices(self):
        """Test when all indices are invalid"""
        indices = [10, 20, 30]
        result = filter_valid_indices(indices, max_length=5)
        self.assertEqual(result, [])
    
    def test_negative_indices(self):
        """Test that negative indices are filtered out"""
        indices = [-5, -1, 0, 5, 10]
        result = filter_valid_indices(indices, max_length=8)
        self.assertEqual(result, [0, 5])
    
    def test_empty_list(self):
        """Test filtering empty list"""
        result = filter_valid_indices([], max_length=10)
        self.assertEqual(result, [])
    
    def test_boundary_values(self):
        """Test boundary values"""
        indices = [0, 4, 5]
        result = filter_valid_indices(indices, max_length=5)
        # 5 is at the boundary (exclusive), should be filtered
        self.assertEqual(result, [0, 4])


class TestValidateSequenceInput(unittest.TestCase):
    """Test suite for validate_sequence_input function"""
    
    def test_all_valid_sequences(self):
        """Test when all sequences are valid"""
        sequences = [
            "ATCGATCGATCG",
            "GGGGCCCCTTTTAAAA",
            "NNNNATCGNNN"
        ]
        valid, issues = validate_sequence_input(sequences, min_length=5)
        self.assertEqual(len(valid), 3)
        self.assertEqual(len(issues), 0)
    
    def test_empty_sequence_filtered(self):
        """Test that empty sequences are filtered"""
        sequences = ["ATCG", "", "GCTA"]
        valid, issues = validate_sequence_input(sequences, min_length=2)
        self.assertEqual(len(valid), 2)
        self.assertGreater(len(issues), 0)
        self.assertIn("empty", issues[0].lower())
    
    def test_short_sequence_warning(self):
        """Test that short sequences generate warnings but are kept"""
        sequences = ["AT", "ATCGATCGATCG"]
        valid, issues = validate_sequence_input(sequences, min_length=10)
        # Short sequence should still be in valid list (just warning)
        self.assertEqual(len(valid), 2)
        self.assertGreater(len(issues), 0)
        self.assertIn("short", issues[0].lower())
    
    def test_invalid_characters_filtered(self):
        """Test that sequences with invalid characters are filtered"""
        sequences = [
            "ATCGATCG",
            "ATCG123XYZ",
            "GCTAGCTA"
        ]
        valid, issues = validate_sequence_input(sequences, min_length=5)
        self.assertEqual(len(valid), 2)
        self.assertGreater(len(issues), 0)
        self.assertIn("invalid", issues[0].lower())
    
    def test_mixed_case_sequences(self):
        """Test that mixed case sequences are valid"""
        sequences = ["ATCGatcg", "GcTaGcTa"]
        valid, issues = validate_sequence_input(sequences, min_length=5)
        self.assertEqual(len(valid), 2)
    
    def test_sequences_with_n(self):
        """Test that sequences with N are valid"""
        sequences = ["NNNNATCGNNNN", "NNNNNNNN"]
        valid, issues = validate_sequence_input(sequences, min_length=5)
        self.assertEqual(len(valid), 2)
    
    def test_all_invalid(self):
        """Test when all sequences are invalid"""
        sequences = ["", "ABC", "123"]
        valid, issues = validate_sequence_input(sequences, min_length=5)
        self.assertEqual(len(valid), 0)
        self.assertGreater(len(issues), 0)


# ═══════════════════════════════════════════════════════════════════════════════
# TEST RUNNER
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)
