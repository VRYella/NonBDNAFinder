"""
╔══════════════════════════════════════════════════════════════════════════════╗
║     TEST: STACKED BAR CLASS-SUBCLASS VISUALIZATION                          ║
║     Validates new stacked bar plot implementation                             ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
import matplotlib.pyplot as plt
from typing import List, Dict, Any

from Utilities.visualization.stacked_bar_class_subclass import (
    plot_stacked_bar_class_subclass,
    plot_nested_pie_chart
)


class TestStackedBarVisualization(unittest.TestCase):
    """Test suite for stacked bar class-subclass visualization."""
    
    def setUp(self):
        """Create test motif data."""
        self.motifs = [
            {'Class': 'G-Quadruplex', 'Subclass': 'Telomeric', 'Start': 0, 'End': 24},
            {'Class': 'G-Quadruplex', 'Subclass': 'Telomeric', 'Start': 50, 'End': 74},
            {'Class': 'G-Quadruplex', 'Subclass': 'Non-Telomeric', 'Start': 100, 'End': 124},
            {'Class': 'Z-DNA', 'Subclass': 'Extended CG', 'Start': 200, 'End': 230},
            {'Class': 'Z-DNA', 'Subclass': 'Extended CG', 'Start': 250, 'End': 280},
            {'Class': 'Z-DNA', 'Subclass': 'Extended CG', 'Start': 300, 'End': 330},
            {'Class': 'Z-DNA', 'Subclass': 'Extended CG', 'Start': 350, 'End': 380},
            {'Class': 'Cruciform', 'Subclass': 'Inverted Repeat', 'Start': 400, 'End': 450},
            {'Class': 'Cruciform', 'Subclass': 'Inverted Repeat', 'Start': 500, 'End': 550},
            {'Class': 'i-Motif', 'Subclass': 'C-Rich', 'Start': 600, 'End': 620},
            {'Class': 'i-Motif', 'Subclass': 'C-Rich', 'Start': 650, 'End': 670},
            {'Class': 'i-Motif', 'Subclass': 'C-Rich', 'Start': 700, 'End': 720},
            {'Class': 'Triplex', 'Subclass': 'Mirror Repeat', 'Start': 800, 'End': 830},
        ]
    
    def test_basic_functionality(self):
        """Test basic stacked bar plot generation."""
        fig = plot_stacked_bar_class_subclass(self.motifs)
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_backward_compatible_wrapper(self):
        """Test that plot_nested_pie_chart wrapper works."""
        fig = plot_nested_pie_chart(self.motifs)
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_empty_motifs(self):
        """Test handling of empty motif list."""
        fig = plot_stacked_bar_class_subclass([])
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_exclude_classes(self):
        """Test exclusion of specific classes."""
        fig = plot_stacked_bar_class_subclass(
            self.motifs,
            exclude_classes=['Z-DNA']
        )
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_all_excluded(self):
        """Test handling when all motifs are excluded."""
        all_classes = ['G-Quadruplex', 'Z-DNA', 'Cruciform', 'i-Motif', 'Triplex']
        fig = plot_stacked_bar_class_subclass(
            self.motifs,
            exclude_classes=all_classes
        )
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_custom_title_and_figsize(self):
        """Test custom title and figure size."""
        fig = plot_stacked_bar_class_subclass(
            self.motifs,
            title="Custom Title",
            figsize=(10, 5)
        )
        self.assertIsInstance(fig, plt.Figure)
        self.assertEqual(fig.get_figwidth(), 10)
        self.assertEqual(fig.get_figheight(), 5)
        plt.close(fig)
    
    def test_single_class(self):
        """Test plot with motifs from a single class."""
        single_class_motifs = [m for m in self.motifs if m['Class'] == 'G-Quadruplex']
        fig = plot_stacked_bar_class_subclass(single_class_motifs)
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_single_subclass_per_class(self):
        """Test plot where each class has only one subclass."""
        single_subclass_motifs = [
            {'Class': 'G-Quadruplex', 'Subclass': 'Telomeric', 'Start': 0, 'End': 24},
            {'Class': 'Z-DNA', 'Subclass': 'Extended CG', 'Start': 200, 'End': 230},
            {'Class': 'Cruciform', 'Subclass': 'Inverted Repeat', 'Start': 400, 'End': 450},
        ]
        fig = plot_stacked_bar_class_subclass(single_subclass_motifs)
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)
    
    def test_many_subclasses(self):
        """Test plot with many subclasses to verify legend limiting."""
        many_subclasses = []
        for i in range(15):
            many_subclasses.append({
                'Class': 'G-Quadruplex',
                'Subclass': f'Subclass_{i}',
                'Start': i * 100,
                'End': i * 100 + 24
            })
        fig = plot_stacked_bar_class_subclass(many_subclasses)
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)


if __name__ == '__main__':
    # Run with verbose output
    unittest.main(verbosity=2)
