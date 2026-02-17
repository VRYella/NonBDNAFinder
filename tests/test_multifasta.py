"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    MULTIFASTA FUNCTIONALITY TEST SUITE                        ║
║          Testing MultiFASTA visualization and export capabilities             ║
╚══════════════════════════════════════════════════════════════════════════════╝

Tests validate:
1. Equal-length detection logic
2. Positional occurrence computation
3. Excel export structure with multiple sheets
4. Unified visualization generation
5. Edge cases (empty sequences, single motif, large datasets)
"""

import sys
import os
import tempfile
import unittest

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Utilities.multifasta_visualizer import (
    MultiFastaVisualizer,
    prepare_multifasta_excel_data
)
from Utilities.utilities import export_multifasta_to_excel
import pandas as pd


class TestMultiFastaVisualizer(unittest.TestCase):
    """Test suite for MultiFASTA visualizer functionality."""
    
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
            'Seq1': 100,
            'Seq2': 100,
            'Seq3': 100
        }
        
        # Different lengths for testing
        self.different_lengths = {
            'Seq1': 100,
            'Seq2': 150,
            'Seq3': 200
        }
    
    def test_equal_length_detection(self):
        """Test detection of equal-length sequences."""
        visualizer = MultiFastaVisualizer(self.annotations_by_sequence)
        
        # Equal length sequences
        sequences_equal = {'Seq1': 'A' * 100, 'Seq2': 'C' * 100, 'Seq3': 'G' * 100}
        self.assertTrue(visualizer.all_sequences_equal_length(sequences_equal))
        
        # Different length sequences
        sequences_diff = {'Seq1': 'A' * 100, 'Seq2': 'C' * 150, 'Seq3': 'G' * 200}
        self.assertFalse(visualizer.all_sequences_equal_length(sequences_diff))
        
        # Empty sequences
        sequences_empty = {}
        self.assertFalse(visualizer.all_sequences_equal_length(sequences_empty))
        
        # Single sequence
        sequences_single = {'Seq1': 'A' * 100}
        self.assertTrue(visualizer.all_sequences_equal_length(sequences_single))
    
    def test_positional_occurrence_computation(self):
        """Test positional occurrence computation."""
        visualizer = MultiFastaVisualizer(self.annotations_by_sequence)
        seq_length = 100
        
        # Compute positional occurrence
        result = visualizer.compute_positional_occurrence(seq_length, by_class=True, by_subclass=True)
        
        # Check result structure
        self.assertIn('Overall', result)
        self.assertIn('Class', result)
        self.assertIn('Subclass', result)
        
        # Check overall counts
        self.assertIsInstance(result['Overall'], dict)
        self.assertGreater(len(result['Overall']), 0)
        
        # Check class-level counts
        self.assertIn('G-Quadruplex', result['Class'])
        self.assertIn('Z-DNA', result['Class'])
        
        # Verify position 10-20 (Seq1 G-Quadruplex)
        g4_counts = result['Class']['G-Quadruplex']
        self.assertGreater(g4_counts.get(10, 0), 0)
        self.assertGreater(g4_counts.get(20, 0), 0)
        
        # Verify counts are within valid range
        for pos, count in result['Overall'].items():
            self.assertGreaterEqual(count, 1)
            self.assertLessEqual(count, len(self.annotations_by_sequence))
    
    def test_unified_summary_generation(self):
        """Test unified summary statistics generation."""
        visualizer = MultiFastaVisualizer(self.annotations_by_sequence)
        
        summary = visualizer.generate_unified_summary()
        
        # Check summary structure
        self.assertIn('sequence_stats', summary)
        self.assertIn('class_distribution', summary)
        self.assertIn('subclass_distribution', summary)
        
        # Verify sequence statistics
        self.assertEqual(len(summary['sequence_stats']), 3)
        
        # Check each sequence has proper stats
        for stat in summary['sequence_stats']:
            self.assertIn('FASTA_ID', stat)
            self.assertIn('Total_Motifs', stat)
            self.assertIn('Class_Counts', stat)
        
        # Verify class distribution
        self.assertIn('G-Quadruplex', summary['class_distribution'])
        self.assertIn('Z-DNA', summary['class_distribution'])
        
        # Check total motif count
        total_motifs = sum(summary['class_distribution'].values())
        expected_total = sum(len(motifs) for motifs in self.annotations_by_sequence.values())
        self.assertEqual(total_motifs, expected_total)
    
    def test_excel_data_preparation(self):
        """Test Excel data preparation for export."""
        # Test with equal-length sequences
        sheet_data = prepare_multifasta_excel_data(
            self.annotations_by_sequence,
            self.sequence_lengths,
            equal_length=True,
            seq_length=100
        )
        
        # Check all required sheets are present
        self.assertIn('All_Motifs', sheet_data)
        self.assertIn('Sequence_Summary', sheet_data)
        self.assertIn('Class_Summary', sheet_data)
        self.assertIn('Positional_Occurrence', sheet_data)  # Should be present for equal length
        
        # Verify All_Motifs structure
        all_motifs = sheet_data['All_Motifs']
        self.assertGreater(len(all_motifs), 0)
        for motif in all_motifs:
            self.assertIn('FASTA_ID', motif)
            self.assertIn('Start', motif)
            self.assertIn('End', motif)
            self.assertIn('Class', motif)
        
        # Verify Sequence_Summary structure
        seq_summary = sheet_data['Sequence_Summary']
        self.assertEqual(len(seq_summary), 3)
        for summary in seq_summary:
            self.assertIn('FASTA_ID', summary)
            self.assertIn('Length', summary)
            self.assertIn('Total_Motifs', summary)
            self.assertIn('Motifs_per_kb', summary)
        
        # Verify Class_Summary structure
        class_summary = sheet_data['Class_Summary']
        self.assertGreater(len(class_summary), 0)
        for summary in class_summary:
            self.assertIn('Class', summary)
            self.assertIn('Subclass', summary)
            self.assertIn('Total_Count', summary)
        
        # Verify Positional_Occurrence structure
        pos_occurrence = sheet_data['Positional_Occurrence']
        self.assertGreater(len(pos_occurrence), 0)
        for entry in pos_occurrence:
            self.assertIn('Position', entry)
            self.assertIn('Class', entry)
            self.assertIn('Subclass', entry)
            self.assertIn('Count', entry)
    
    def test_excel_data_preparation_different_lengths(self):
        """Test Excel data preparation with different sequence lengths."""
        # Test with different-length sequences
        sheet_data = prepare_multifasta_excel_data(
            self.annotations_by_sequence,
            self.different_lengths,
            equal_length=False,
            seq_length=None
        )
        
        # Check required sheets
        self.assertIn('All_Motifs', sheet_data)
        self.assertIn('Sequence_Summary', sheet_data)
        self.assertIn('Class_Summary', sheet_data)
        
        # Positional_Occurrence should NOT be present for different lengths
        self.assertNotIn('Positional_Occurrence', sheet_data)
    
    def test_excel_export_integration(self):
        """Test full Excel export with file creation."""
        with tempfile.NamedTemporaryFile(suffix='.xlsx', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            # Export to Excel
            result = export_multifasta_to_excel(
                self.annotations_by_sequence,
                self.sequence_lengths,
                tmp_path,
                equal_length=True,
                seq_length=100
            )
            
            # Check success message
            self.assertIn("success", result.lower())
            
            # Verify file was created
            self.assertTrue(os.path.exists(tmp_path))
            self.assertGreater(os.path.getsize(tmp_path), 0)
            
            # Read and verify Excel structure
            excel_file = pd.ExcelFile(tmp_path, engine='openpyxl')
            sheet_names = excel_file.sheet_names
            
            # Check required sheets
            self.assertIn('All_Motifs', sheet_names)
            self.assertIn('Sequence_Summary', sheet_names)
            self.assertIn('Class_Summary', sheet_names)
            self.assertIn('Positional_Occurrence', sheet_names)
            
            # Verify All_Motifs sheet structure
            df_motifs = pd.read_excel(tmp_path, sheet_name='All_Motifs')
            self.assertIn('FASTA_ID', df_motifs.columns)
            self.assertIn('Class', df_motifs.columns)
            self.assertIn('Start', df_motifs.columns)
            self.assertIn('End', df_motifs.columns)
            
            # Verify Sequence_Summary sheet
            df_seq = pd.read_excel(tmp_path, sheet_name='Sequence_Summary')
            self.assertIn('FASTA_ID', df_seq.columns)
            self.assertIn('Length', df_seq.columns)
            self.assertIn('Total_Motifs', df_seq.columns)
            self.assertIn('Motifs_per_kb', df_seq.columns)
            self.assertEqual(len(df_seq), 3)
            
        finally:
            # Clean up
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
    
    def test_visualization_generation(self):
        """Test visualization generation (non-rendering check)."""
        visualizer = MultiFastaVisualizer(self.annotations_by_sequence)
        
        # Test class distribution plot generation
        try:
            fig = visualizer.generate_class_distribution_plot()
            self.assertIsNotNone(fig)
        except Exception as e:
            self.fail(f"Class distribution plot failed: {e}")
        
        # Test density heatmap generation
        try:
            fig = visualizer.generate_density_heatmap(self.sequence_lengths)
            self.assertIsNotNone(fig)
        except Exception as e:
            self.fail(f"Density heatmap failed: {e}")
        
        # Test positional panels generation
        try:
            figures = visualizer.generate_positional_panels(100, smooth=False)
            self.assertIsInstance(figures, dict)
            self.assertIn('overall', figures)
        except Exception as e:
            self.fail(f"Positional panels failed: {e}")
    
    def test_empty_annotations(self):
        """Test handling of empty annotations."""
        empty_annotations = {
            'Seq1': [],
            'Seq2': [],
            'Seq3': []
        }
        
        visualizer = MultiFastaVisualizer(empty_annotations)
        
        # Should not crash with empty data
        summary = visualizer.generate_unified_summary()
        self.assertEqual(len(summary['sequence_stats']), 3)
        self.assertEqual(len(summary['class_distribution']), 0)
        
        # Positional occurrence with empty data
        result = visualizer.compute_positional_occurrence(100)
        self.assertEqual(len(result['Overall']), 0)
    
    def test_single_sequence(self):
        """Test handling of single sequence (edge case)."""
        single_seq = {
            'OnlySeq': [
                {'Start': 10, 'End': 20, 'Class': 'G-Quadruplex', 'Subclass': 'G4_Canonical', 'Length': 11}
            ]
        }
        
        visualizer = MultiFastaVisualizer(single_seq)
        summary = visualizer.generate_unified_summary()
        
        self.assertEqual(len(summary['sequence_stats']), 1)
        self.assertEqual(summary['sequence_stats'][0]['Total_Motifs'], 1)


class TestMultiFastaEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""
    
    def test_large_position_range(self):
        """Test handling of large position ranges."""
        annotations = {
            'LargeSeq': [
                {'Start': 1, 'End': 1000000, 'Class': 'Test', 'Subclass': 'TestSub', 'Length': 1000000}
            ]
        }
        
        visualizer = MultiFastaVisualizer(annotations)
        # Should handle large ranges without memory issues
        result = visualizer.compute_positional_occurrence(1000000)
        self.assertIsInstance(result, dict)
    
    def test_overlapping_motifs(self):
        """Test handling of overlapping motifs in positional occurrence."""
        annotations = {
            'Seq1': [
                {'Start': 10, 'End': 30, 'Class': 'A', 'Subclass': 'A1', 'Length': 21},
                {'Start': 20, 'End': 40, 'Class': 'B', 'Subclass': 'B1', 'Length': 21}
            ]
        }
        
        visualizer = MultiFastaVisualizer(annotations)
        result = visualizer.compute_positional_occurrence(100)
        
        # Position 20-30 should have count >= 2 (overlapping)
        self.assertGreaterEqual(result['Overall'].get(25, 0), 2)


if __name__ == '__main__':
    # Run tests
    unittest.main(verbosity=2)
