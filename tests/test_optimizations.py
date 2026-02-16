"""
Unit tests for the Aho-Corasick matcher and optimization components.
"""

import unittest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Utilities.ac_matcher import (
    AhoCorasickMatcher, 
    PatternGroup, 
    merge_pattern_groups,
    create_simple_matcher
)


class TestAhoCorasickMatcher(unittest.TestCase):
    """Test the Aho-Corasick multi-pattern matcher."""
    
    def test_basic_matching(self):
        """Test basic pattern matching functionality."""
        matcher = AhoCorasickMatcher()
        matcher.add_pattern("GGGG", detector="g_quadruplex", subclass="G4")
        matcher.add_pattern("CCCC", detector="i_motif", subclass="iM")
        matcher.build()
        
        sequence = "ATGGGGCCCCTA"
        matches = list(matcher.search(sequence))
        
        # Should find both patterns
        self.assertEqual(len(matches), 2)
        
        # Check GGGG match
        gggg_match = [m for m in matches if m[2] == "GGGG"][0]
        self.assertEqual(gggg_match[0], 2)  # Start position
        self.assertEqual(gggg_match[1], 6)  # End position (exclusive)
        self.assertEqual(gggg_match[3]['detector'], 'g_quadruplex')
        
        # Check CCCC match
        cccc_match = [m for m in matches if m[2] == "CCCC"][0]
        self.assertEqual(cccc_match[0], 6)
        self.assertEqual(cccc_match[1], 10)
        self.assertEqual(cccc_match[3]['detector'], 'i_motif')
    
    def test_case_insensitive(self):
        """Test case-insensitive matching."""
        matcher = AhoCorasickMatcher(case_sensitive=False)
        matcher.add_pattern("ATCG", test=True)
        matcher.build()
        
        # Should match both upper and lower case
        matches_upper = list(matcher.search("ATCGATCG"))
        matches_lower = list(matcher.search("atcgatcg"))
        matches_mixed = list(matcher.search("AtCgAtCg"))
        
        self.assertEqual(len(matches_upper), 2)
        self.assertEqual(len(matches_lower), 2)
        self.assertEqual(len(matches_mixed), 2)
    
    def test_overlapping_patterns(self):
        """Test detection of overlapping patterns."""
        matcher = AhoCorasickMatcher()
        matcher.add_pattern("GGG", short_g4=True)
        matcher.add_pattern("GGGG", long_g4=True)
        matcher.build()
        
        sequence = "ATGGGGTA"
        matches = list(matcher.search(sequence))
        
        # Should find both the GGG and GGGG patterns
        self.assertGreaterEqual(len(matches), 2)
    
    def test_empty_sequence(self):
        """Test handling of empty sequence."""
        matcher = AhoCorasickMatcher()
        matcher.add_pattern("ATCG")
        matcher.build()
        
        matches = list(matcher.search(""))
        self.assertEqual(len(matches), 0)
    
    def test_no_matches(self):
        """Test sequence with no pattern matches."""
        matcher = AhoCorasickMatcher()
        matcher.add_pattern("GGGG")
        matcher.build()
        
        sequence = "ATATATATATAT"  # No GGGG pattern
        matches = list(matcher.search(sequence))
        self.assertEqual(len(matches), 0)
    
    def test_pattern_group(self):
        """Test PatternGroup functionality."""
        group = PatternGroup("g_quadruplex")
        group.add_pattern("GGGG", subclass="canonical_G4")
        group.add_pattern("GGG", subclass="weak_G4")
        
        self.assertEqual(len(group), 2)
        
        matcher = group.to_matcher()
        matcher.build()
        
        sequence = "ATGGGGTA"
        matches = list(matcher.search(sequence))
        
        # Should find patterns from the group
        self.assertGreater(len(matches), 0)
        
        # All matches should have detector name
        for match in matches:
            self.assertEqual(match[3]['detector'], 'g_quadruplex')
    
    def test_merge_pattern_groups(self):
        """Test merging multiple pattern groups."""
        g4_group = PatternGroup("g_quadruplex")
        g4_group.add_pattern("GGGG")
        
        im_group = PatternGroup("i_motif")
        im_group.add_pattern("CCCC")
        
        matcher = merge_pattern_groups([g4_group, im_group])
        matcher.build()
        
        sequence = "ATGGGGCCCCTA"
        matches = list(matcher.search(sequence))
        
        # Should find patterns from both groups
        self.assertEqual(len(matches), 2)
        
        detectors = set(m[3]['detector'] for m in matches)
        self.assertEqual(detectors, {'g_quadruplex', 'i_motif'})
    
    def test_create_simple_matcher(self):
        """Test convenience function for creating simple matchers."""
        matcher = create_simple_matcher(
            ["GGGG", "CCCC", "ATCG"],
            detector="test_detector"
        )
        
        sequence = "ATCGGGGGCCCCATCG"
        matches = list(matcher.search(sequence))
        
        # Should find all three patterns
        self.assertGreaterEqual(len(matches), 3)
        
        # All should have the same detector
        for match in matches:
            self.assertEqual(match[3]['detector'], 'test_detector')
    
    def test_stats(self):
        """Test matcher statistics."""
        matcher = AhoCorasickMatcher()
        matcher.add_pattern("GGGG")
        matcher.add_pattern("CCCC")
        
        stats = matcher.get_stats()
        
        self.assertEqual(stats['num_patterns'], 2)
        self.assertFalse(stats['built'])
        
        matcher.build()
        stats = matcher.get_stats()
        
        self.assertTrue(stats['built'])
    
    def test_error_handling(self):
        """Test error handling for invalid operations."""
        matcher = AhoCorasickMatcher()
        
        # Cannot search before building
        with self.assertRaises(RuntimeError):
            list(matcher.search("ATCG"))
        
        # Cannot build with no patterns
        with self.assertRaises(RuntimeError):
            matcher.build()
        
        # Cannot add empty pattern
        with self.assertRaises(ValueError):
            matcher.add_pattern("")


class TestOptimizationInfo(unittest.TestCase):
    """Test optimization information functions."""
    
    def test_get_optimization_info(self):
        """Test getting optimization information."""
        from Utilities.nonbscanner_optimized import get_optimization_info
        
        info = get_optimization_info()
        
        # Should have required keys
        self.assertIn('ac_available', info)
        self.assertIn('nonbscanner_available', info)
        self.assertIn('optimization_enabled', info)
        self.assertIn('expected_speedup', info)
        
        # Values should be boolean or string
        self.assertIsInstance(info['ac_available'], bool)
        self.assertIsInstance(info['nonbscanner_available'], bool)
        self.assertIsInstance(info['optimization_enabled'], bool)
        self.assertIsInstance(info['expected_speedup'], str)


def run_tests():
    """Run all tests and return results."""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestAhoCorasickMatcher))
    suite.addTests(loader.loadTestsFromTestCase(TestOptimizationInfo))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result


if __name__ == '__main__':
    result = run_tests()
    
    # Exit with appropriate code
    sys.exit(0 if result.wasSuccessful() else 1)
