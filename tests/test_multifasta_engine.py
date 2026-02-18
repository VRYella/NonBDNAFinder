"""
╔══════════════════════════════════════════════════════════════════════════════╗
║              MULTIFASTA ENGINE TEST SUITE                                     ║
║          Testing aggregation-first engine for >20 sequences                   ║
╚══════════════════════════════════════════════════════════════════════════════╝

Tests validate:
1. Global class distribution across sequences
2. Class sequence coverage percentages
3. Motif density per sequence calculation
4. Positional conservation for equal-length sequences
5. Memory safety with large datasets
6. Edge cases (empty sequences, single sequence, etc.)
"""

import sys
import os
import unittest

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Utilities.multifasta_engine import MultiFastaEngine


class TestMultiFastaEngine(unittest.TestCase):
    """Test suite for MultiFASTA aggregation engine."""
    
    def setUp(self):
        """Set up test data."""
        # Create sample annotations for multiple sequences
        self.annotations_by_sequence = {
            'Seq1': [
                {'Start': 10, 'End': 20, 'Class': 'G-Quadruplex', 'Subclass': 'G4_Canonical', 'Length': 11},
                {'Start': 30, 'End': 40, 'Class': 'Z-DNA', 'Subclass': 'ZDNA_CG_Alternating', 'Length': 11},
                {'Start': 50, 'End': 60, 'Class': 'G-Quadruplex', 'Subclass': 'G4_Telomeric', 'Length': 11},
            ],
            'Seq2': [
                {'Start': 15, 'End': 25, 'Class': 'G-Quadruplex', 'Subclass': 'G4_Canonical', 'Length': 11},
                {'Start': 35, 'End': 45, 'Class': 'i-Motif', 'Subclass': 'iMotif_Canonical', 'Length': 11},
            ],
            'Seq3': [
                {'Start': 12, 'End': 22, 'Class': 'Z-DNA', 'Subclass': 'ZDNA_CG_Alternating', 'Length': 11},
                {'Start': 32, 'End': 42, 'Class': 'Z-DNA', 'Subclass': 'ZDNA_AT_Alternating', 'Length': 11},
                {'Start': 52, 'End': 62, 'Class': 'G-Quadruplex', 'Subclass': 'G4_Canonical', 'Length': 11},
            ]
        }
        
        self.sequence_lengths = {
            'Seq1': 1000,
            'Seq2': 1000,
            'Seq3': 1000
        }
    
    def test_global_class_distribution(self):
        """Test global class distribution across all sequences."""
        engine = MultiFastaEngine(self.annotations_by_sequence, self.sequence_lengths)
        
        distribution = engine.global_class_distribution()
        
        # Check structure
        self.assertIsInstance(distribution, dict)
        
        # Verify expected classes are present
        self.assertIn('G-Quadruplex', distribution)
        self.assertIn('Z-DNA', distribution)
        self.assertIn('i-Motif', distribution)
        
        # Verify counts
        # Seq1: 2 G-Quadruplex, 1 Z-DNA
        # Seq2: 1 G-Quadruplex, 1 i-Motif
        # Seq3: 1 G-Quadruplex, 2 Z-DNA
        # Total: 4 G-Quadruplex, 3 Z-DNA, 1 i-Motif
        self.assertEqual(distribution['G-Quadruplex'], 4)
        self.assertEqual(distribution['Z-DNA'], 3)
        self.assertEqual(distribution['i-Motif'], 1)
        
        # Verify total
        total = sum(distribution.values())
        expected_total = sum(len(motifs) for motifs in self.annotations_by_sequence.values())
        self.assertEqual(total, expected_total)
    
    def test_class_sequence_coverage(self):
        """Test class sequence coverage percentages."""
        engine = MultiFastaEngine(self.annotations_by_sequence, self.sequence_lengths)
        
        coverage = engine.class_sequence_coverage()
        
        # Check structure
        self.assertIsInstance(coverage, dict)
        
        # Verify percentages
        # G-Quadruplex: in all 3 sequences = 100%
        # Z-DNA: in Seq1 and Seq3 = 66.67%
        # i-Motif: in Seq2 only = 33.33%
        self.assertAlmostEqual(coverage['G-Quadruplex'], 100.0, places=1)
        self.assertAlmostEqual(coverage['Z-DNA'], 66.67, places=1)
        self.assertAlmostEqual(coverage['i-Motif'], 33.33, places=1)
        
        # Verify all values are percentages
        for cls, pct in coverage.items():
            self.assertGreaterEqual(pct, 0.0)
            self.assertLessEqual(pct, 100.0)
    
    def test_motif_density_per_sequence(self):
        """Test motif density calculation per sequence."""
        engine = MultiFastaEngine(self.annotations_by_sequence, self.sequence_lengths)
        
        densities = engine.motif_density_per_sequence()
        
        # Check structure
        self.assertIsInstance(densities, list)
        self.assertEqual(len(densities), 3)
        
        # Verify all densities are positive numbers
        for density in densities:
            self.assertIsInstance(density, (int, float))
            self.assertGreater(density, 0)
        
        # Calculate expected densities
        # Seq1: 3 motifs / 1kb = 3.0
        # Seq2: 2 motifs / 1kb = 2.0
        # Seq3: 3 motifs / 1kb = 3.0
        expected_densities = [3.0, 2.0, 3.0]
        
        for actual, expected in zip(densities, expected_densities):
            self.assertAlmostEqual(actual, expected, places=1)
    
    def test_positional_conservation(self):
        """Test positional conservation for equal-length sequences."""
        engine = MultiFastaEngine(self.annotations_by_sequence, self.sequence_lengths)
        
        seq_length = 1000
        conservation = engine.positional_conservation(seq_length)
        
        # Check structure
        self.assertIsInstance(conservation, list)
        
        # With 1000bp sequences and 1000bp windows, should have 1 bin
        self.assertEqual(len(conservation), 1)
        
        # All 3 sequences have motifs in the first 1000bp, so conservation should be 3
        self.assertEqual(conservation[0], 3)
    
    def test_positional_conservation_multiple_windows(self):
        """Test positional conservation with multiple windows."""
        # Create annotations spanning multiple windows
        annotations = {
            'Seq1': [
                {'Start': 100, 'End': 200, 'Class': 'A', 'Subclass': 'A1'},
                {'Start': 1500, 'End': 1600, 'Class': 'B', 'Subclass': 'B1'},
            ],
            'Seq2': [
                {'Start': 150, 'End': 250, 'Class': 'A', 'Subclass': 'A1'},
                {'Start': 2500, 'End': 2600, 'Class': 'C', 'Subclass': 'C1'},
            ],
            'Seq3': [
                {'Start': 1400, 'End': 1500, 'Class': 'B', 'Subclass': 'B1'},
            ]
        }
        
        sequence_lengths = {'Seq1': 3000, 'Seq2': 3000, 'Seq3': 3000}
        engine = MultiFastaEngine(annotations, sequence_lengths)
        
        conservation = engine.positional_conservation(3000)
        
        # Should have 3 bins (0-999, 1000-1999, 2000-2999)
        self.assertEqual(len(conservation), 3)
        
        # Bin 0: Seq1 and Seq2 have motifs
        self.assertEqual(conservation[0], 2)
        
        # Bin 1: Seq1 and Seq3 have motifs
        self.assertEqual(conservation[1], 2)
        
        # Bin 2: Only Seq2 has motifs
        self.assertEqual(conservation[2], 1)
    
    def test_empty_annotations(self):
        """Test handling of empty annotations."""
        empty_annotations = {
            'Seq1': [],
            'Seq2': [],
            'Seq3': []
        }
        
        engine = MultiFastaEngine(empty_annotations, self.sequence_lengths)
        
        # Global distribution should be empty
        distribution = engine.global_class_distribution()
        self.assertEqual(len(distribution), 0)
        
        # Coverage should be empty
        coverage = engine.class_sequence_coverage()
        self.assertEqual(len(coverage), 0)
        
        # Densities should all be 0
        densities = engine.motif_density_per_sequence()
        self.assertEqual(len(densities), 3)
        for density in densities:
            self.assertEqual(density, 0)
        
        # Conservation should be all zeros
        conservation = engine.positional_conservation(1000)
        self.assertTrue(all(c == 0 for c in conservation))
    
    def test_single_sequence(self):
        """Test handling of single sequence."""
        single_seq = {
            'OnlySeq': [
                {'Start': 10, 'End': 20, 'Class': 'G-Quadruplex', 'Subclass': 'G4_Canonical'}
            ]
        }
        
        single_length = {'OnlySeq': 1000}
        engine = MultiFastaEngine(single_seq, single_length)
        
        # Should work with single sequence
        distribution = engine.global_class_distribution()
        self.assertEqual(distribution['G-Quadruplex'], 1)
        
        coverage = engine.class_sequence_coverage()
        self.assertEqual(coverage['G-Quadruplex'], 100.0)
        
        densities = engine.motif_density_per_sequence()
        self.assertEqual(len(densities), 1)
        self.assertAlmostEqual(densities[0], 1.0, places=1)
    
    def test_large_dataset_memory_safety(self):
        """Test memory safety with large dataset (>20 sequences)."""
        # Create 25 sequences with multiple motifs each
        large_annotations = {}
        large_lengths = {}
        
        for i in range(25):
            seq_id = f'Seq{i}'
            large_annotations[seq_id] = [
                {'Start': j * 100, 'End': j * 100 + 50, 'Class': f'Class{j % 5}', 
                 'Subclass': f'Sub{j % 10}', 'Length': 51}
                for j in range(10)
            ]
            large_lengths[seq_id] = 10000
        
        # Should not crash or consume excessive memory
        engine = MultiFastaEngine(large_annotations, large_lengths)
        
        distribution = engine.global_class_distribution()
        self.assertEqual(sum(distribution.values()), 25 * 10)  # 25 sequences * 10 motifs
        
        coverage = engine.class_sequence_coverage()
        self.assertEqual(len(coverage), 5)  # 5 different classes
        
        densities = engine.motif_density_per_sequence()
        self.assertEqual(len(densities), 25)  # 25 sequences
        
        conservation = engine.positional_conservation(10000)
        self.assertEqual(len(conservation), 10)  # 10 windows (10000bp / 1000bp)
    
    def test_different_sequence_lengths(self):
        """Test handling of sequences with different lengths."""
        different_lengths = {
            'Seq1': 500,
            'Seq2': 2000,
            'Seq3': 1500
        }
        
        engine = MultiFastaEngine(self.annotations_by_sequence, different_lengths)
        
        # Densities should be calculated correctly for each length
        densities = engine.motif_density_per_sequence()
        
        # Seq1: 3 motifs / 0.5kb = 6.0
        # Seq2: 2 motifs / 2kb = 1.0
        # Seq3: 3 motifs / 1.5kb = 2.0
        self.assertAlmostEqual(densities[0], 6.0, places=1)
        self.assertAlmostEqual(densities[1], 1.0, places=1)
        self.assertAlmostEqual(densities[2], 2.0, places=1)


class TestMultiFastaEngineEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""
    
    def test_missing_class_field(self):
        """Test handling of motifs with missing Class field."""
        annotations = {
            'Seq1': [
                {'Start': 10, 'End': 20},  # Missing Class
                {'Start': 30, 'End': 40, 'Class': 'G-Quadruplex'}
            ]
        }
        
        lengths = {'Seq1': 1000}
        engine = MultiFastaEngine(annotations, lengths)
        
        # Should handle None class gracefully
        distribution = engine.global_class_distribution()
        self.assertIn(None, distribution)
        self.assertIn('G-Quadruplex', distribution)
    
    def test_zero_length_sequence(self):
        """Test handling of zero-length sequence."""
        annotations = {
            'Seq1': [
                {'Start': 10, 'End': 20, 'Class': 'Test', 'Subclass': 'Test_Sub', 'Length': 11}
            ]
        }
        
        lengths = {'Seq1': 0}
        engine = MultiFastaEngine(annotations, lengths)
        
        # Should not crash with division by zero
        densities = engine.motif_density_per_sequence()
        self.assertEqual(len(densities), 1)
        # Density should be 0 or handled gracefully
        self.assertIsInstance(densities[0], (int, float))
    
    def test_motif_spanning_multiple_bins(self):
        """Test motif that spans multiple conservation bins."""
        annotations = {
            'Seq1': [
                {'Start': 900, 'End': 1100, 'Class': 'A', 'Subclass': 'A1', 'Length': 201}  # Spans bins 0 and 1
            ]
        }
        
        lengths = {'Seq1': 2000}
        engine = MultiFastaEngine(annotations, lengths)
        
        conservation = engine.positional_conservation(2000)
        
        # Should mark both bins as covered
        self.assertEqual(conservation[0], 1)
        self.assertEqual(conservation[1], 1)


if __name__ == '__main__':
    # Run tests
    unittest.main(verbosity=2)
