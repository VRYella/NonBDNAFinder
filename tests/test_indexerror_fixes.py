"""
╔══════════════════════════════════════════════════════════════════════════════╗
║            INDEX ERROR FIXES INTEGRATION TEST SUITE                           ║
║        Testing IndexError Prevention in Real-World Scenarios                  ║
╚══════════════════════════════════════════════════════════════════════════════╝

Tests validate:
1. False positive analysis handles invalid indices gracefully
2. False negative analysis handles invalid indices gracefully
3. System doesn't crash when indices are out of bounds
4. Appropriate warnings are logged for invalid data
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
import pandas as pd
from NBSTVALIDATION.validation_analysis_subclass import (
    analyze_false_positives,
    analyze_false_negatives
)


class TestFalsePositivesIndexSafety(unittest.TestCase):
    """Test suite for false positives analysis with invalid indices"""
    
    def test_valid_indices_only(self):
        """Test FP analysis with all valid indices"""
        nbf_motifs = [
            {'Start': 1, 'End': 10, 'Length': 9, 'Score': 0.8},
            {'Start': 20, 'End': 30, 'Length': 10, 'Score': 0.9},
            {'Start': 50, 'End': 65, 'Length': 15, 'Score': 0.7}
        ]
        
        result = {
            'nbf_only_indices': [0, 1, 2]
        }
        
        # Should process all motifs successfully
        fp_analysis = analyze_false_positives(nbf_motifs, result)
        
        self.assertEqual(fp_analysis['count'], 3)
        self.assertGreater(fp_analysis['mean_length'], 0)
        self.assertGreater(fp_analysis['mean_score'], 0)
    
    def test_some_invalid_indices(self):
        """Test FP analysis handles out-of-bounds indices gracefully"""
        nbf_motifs = [
            {'Start': 1, 'End': 10, 'Length': 9, 'Score': 0.8},
            {'Start': 20, 'End': 30, 'Length': 10, 'Score': 0.9}
        ]
        
        result = {
            'nbf_only_indices': [0, 1, 5, 10, 15]  # 5, 10, 15 are invalid
        }
        
        # Should not crash, only count valid indices
        fp_analysis = analyze_false_positives(nbf_motifs, result)
        
        # Should only count the 2 valid motifs
        self.assertEqual(fp_analysis['count'], 2)
        self.assertGreater(fp_analysis['mean_length'], 0)
    
    def test_all_invalid_indices(self):
        """Test FP analysis when all indices are invalid"""
        nbf_motifs = [
            {'Start': 1, 'End': 10, 'Length': 9}
        ]
        
        result = {
            'nbf_only_indices': [10, 20, 30]  # All invalid
        }
        
        # Should return empty result gracefully without crashing
        fp_analysis = analyze_false_positives(nbf_motifs, result)
        
        self.assertEqual(fp_analysis['count'], 0)
        self.assertEqual(fp_analysis['mean_length'], 0)
        self.assertEqual(fp_analysis['mean_score'], 0)
    
    def test_empty_indices_list(self):
        """Test FP analysis with empty indices list"""
        nbf_motifs = [
            {'Start': 1, 'End': 10, 'Length': 9}
        ]
        
        result = {
            'nbf_only_indices': []
        }
        
        # Should handle empty list gracefully
        fp_analysis = analyze_false_positives(nbf_motifs, result)
        
        self.assertEqual(fp_analysis['count'], 0)
        self.assertEqual(fp_analysis['common_characteristics'], 'N/A')
    
    def test_negative_indices(self):
        """Test FP analysis filters negative indices"""
        nbf_motifs = [
            {'Start': 1, 'End': 10, 'Length': 9, 'Score': 0.8},
            {'Start': 20, 'End': 30, 'Length': 10, 'Score': 0.9}
        ]
        
        result = {
            'nbf_only_indices': [-1, 0, 1, -5]  # Negative indices should be filtered
        }
        
        # Should only use valid positive indices
        fp_analysis = analyze_false_positives(nbf_motifs, result)
        
        self.assertEqual(fp_analysis['count'], 2)


class TestFalseNegativesIndexSafety(unittest.TestCase):
    """Test suite for false negatives analysis with invalid indices"""
    
    def test_valid_indices_only(self):
        """Test FN analysis with all valid indices"""
        nbst_data = {
            'Start': [10, 50, 100],
            'Stop': [20, 65, 120],
            'Class': ['GQ', 'GQ', 'curved']
        }
        nbst_df = pd.DataFrame(nbst_data)
        
        result = {
            'nbst_only_indices': [0, 1, 2]
        }
        
        # Should process all rows successfully
        fn_analysis = analyze_false_negatives(nbst_df, result)
        
        self.assertEqual(fn_analysis['count'], 3)
        self.assertGreater(fn_analysis['mean_length'], 0)
    
    def test_some_invalid_indices(self):
        """Test FN analysis handles out-of-bounds indices gracefully"""
        nbst_data = {
            'Start': [10, 50],
            'Stop': [20, 65],
            'Class': ['GQ', 'curved']
        }
        nbst_df = pd.DataFrame(nbst_data)
        
        result = {
            'nbst_only_indices': [0, 1, 5, 10]  # 5 and 10 are invalid
        }
        
        # Should not crash, only count valid indices
        fn_analysis = analyze_false_negatives(nbst_df, result)
        
        # Should only count the 2 valid rows
        self.assertEqual(fn_analysis['count'], 2)
        self.assertGreater(fn_analysis['mean_length'], 0)
    
    def test_all_invalid_indices(self):
        """Test FN analysis when all indices are invalid"""
        nbst_data = {
            'Start': [10],
            'Stop': [20],
            'Class': ['GQ']
        }
        nbst_df = pd.DataFrame(nbst_data)
        
        result = {
            'nbst_only_indices': [10, 20, 30]  # All invalid
        }
        
        # Should return empty result gracefully without crashing
        fn_analysis = analyze_false_negatives(nbst_df, result)
        
        self.assertEqual(fn_analysis['count'], 0)
        self.assertEqual(fn_analysis['mean_length'], 0)
        self.assertIn('out of bounds', fn_analysis['potential_causes'])
    
    def test_empty_indices_list(self):
        """Test FN analysis with empty indices list"""
        nbst_data = {
            'Start': [10],
            'Stop': [20],
            'Class': ['GQ']
        }
        nbst_df = pd.DataFrame(nbst_data)
        
        result = {
            'nbst_only_indices': []
        }
        
        # Should handle empty list gracefully
        fn_analysis = analyze_false_negatives(nbst_df, result)
        
        self.assertEqual(fn_analysis['count'], 0)
        self.assertEqual(fn_analysis['potential_causes'], 'N/A')
    
    def test_empty_dataframe(self):
        """Test FN analysis with empty DataFrame"""
        nbst_df = pd.DataFrame()
        
        result = {
            'nbst_only_indices': [0, 1, 2]
        }
        
        # Should handle empty DataFrame gracefully
        fn_analysis = analyze_false_negatives(nbst_df, result)
        
        self.assertEqual(fn_analysis['count'], 0)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and boundary conditions"""
    
    def test_fp_with_missing_score_field(self):
        """Test FP analysis when motifs lack Score field"""
        nbf_motifs = [
            {'Start': 1, 'End': 10, 'Length': 9},  # No Score field
            {'Start': 20, 'End': 30, 'Length': 10}
        ]
        
        result = {
            'nbf_only_indices': [0, 1]
        }
        
        # Should handle missing Score field gracefully
        fp_analysis = analyze_false_positives(nbf_motifs, result)
        
        self.assertEqual(fp_analysis['count'], 2)
        self.assertEqual(fp_analysis['mean_score'], 0)
    
    def test_fn_with_end_column_instead_of_stop(self):
        """Test FN analysis when DataFrame uses 'End' instead of 'Stop'"""
        nbst_data = {
            'Start': [10, 50],
            'End': [20, 65],  # Using 'End' instead of 'Stop'
            'Class': ['GQ', 'curved']
        }
        nbst_df = pd.DataFrame(nbst_data)
        
        result = {
            'nbst_only_indices': [0, 1]
        }
        
        # Should handle both column names
        fn_analysis = analyze_false_negatives(nbst_df, result)
        
        self.assertEqual(fn_analysis['count'], 2)
        self.assertGreater(fn_analysis['mean_length'], 0)
    
    def test_very_large_index(self):
        """Test with extremely large index values"""
        nbf_motifs = [
            {'Start': 1, 'End': 10, 'Length': 9, 'Score': 0.8}
        ]
        
        result = {
            'nbf_only_indices': [0, 1000000, 9999999]
        }
        
        # Should filter out large indices without error
        fp_analysis = analyze_false_positives(nbf_motifs, result)
        
        self.assertEqual(fp_analysis['count'], 1)


# ═══════════════════════════════════════════════════════════════════════════════
# TEST RUNNER
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)
