"""
╔══════════════════════════════════════════════════════════════════════════════╗
║            VISUALIZATION INTEGRATION TEST SUITE                              ║
║        Testing Detector Output Compatibility with Visualizations              ║
╚══════════════════════════════════════════════════════════════════════════════╝

Tests validate:
1. All detector outputs have visualization-compatible structure
2. Motif class colors are properly mapped
3. Advanced parameters are reported for specialized visualizations
4. Score normalization is consistent across detectors
5. Position data can be plotted on linear/Manhattan tracks
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
from typing import Dict, List, Any, Set
import pandas as pd

# Import all detectors
from Detectors.gquad.detector import GQuadruplexDetector
from Detectors.zdna.detector import ZDNADetector
from Detectors.cruciform.detector import CruciformDetector
from Detectors.triplex.detector import TriplexDetector
from Detectors.imotif.detector import IMotifDetector
from Detectors.rloop.detector import RLoopDetector
from Detectors.slipped.detector import SlippedDNADetector
from Detectors.curved.detector import CurvedDNADetector
from Detectors.aphilic.detector import APhilicDetector

# Import visualization standards
from Utilities.config.colors import UNIFIED_MOTIF_COLORS, VISUALIZATION_MOTIF_COLORS
from Utilities.visualization.standards import (
    NATURE_MOTIF_COLORS,
    MetricFilter,
    FigurePanel
)
from Utilities.config.motif_taxonomy import (
    VALID_CLASSES,
    VALID_SUBCLASSES,
    CLASS_TO_SUBCLASSES
)


# ═══════════════════════════════════════════════════════════════════════════════
# COMPREHENSIVE TEST SEQUENCES
# ═══════════════════════════════════════════════════════════════════════════════

TEST_SEQUENCES = {
    'G4_TELOMERIC': 'TTAGGGTTAGGGTTAGGGTTAGGG',
    'G4_CANONICAL': 'GGGATGGGCTGGGAAGGG',
    'ZDNA_EXTENDED': 'CGCGCGCGCGCGCGCGCGCGCGCGCGCGCG',
    'EGZ_CGG': 'CGGCGGCGGCGGCGG',
    'CRUCIFORM': 'ATCGATCGATCGGGGCGATCGATCGAT',
    'TRIPLEX_MIRROR': 'GAAGAAGAAGAAAAGAAGAAGAAG',
    'STICKY_GAA': 'GAAGAAGAAGAAGAAGAA',
    'IMOTIF': 'CCCCTCCCCTCCCCTCCCC',
    # R-Loop: G-rich RIZ + linker region + G-rich REZ (QmRLFS model)
    'RLOOP': 'GGGGCGGGGGCGGGGCGGGGGCGGGG' + ('A' * 100) + 'GGGGCGGGG',
    'STR_CAG': 'CAGCAGCAGCAGCAGCAGCAGCAG',
    'CURVED': 'AAAAAATAAAAATAAAAATAAAAAAAAAA',
    'APHILIC': 'AACGAACGAACGAACGAACGAACGAACG',
    'MIXED': 'TTAGGGTTAGGGTTAGGGTTAGGG' + 'CGCGCGCGCGCGCG' + 'CCCCTCCCCTCCCCTCCCC',
}


# ═══════════════════════════════════════════════════════════════════════════════
# DETECTOR PARAMETER SPECIFICATIONS
# These define the expected advanced parameters each detector should report
# ═══════════════════════════════════════════════════════════════════════════════

DETECTOR_EXPECTED_PARAMS = {
    'G-Quadruplex': {
        'core': {'ID', 'Class', 'Subclass', 'Start', 'End', 'Length', 'Sequence', 'Score', 'Strand', 'Method', 'Pattern_ID'},
        'detector_specific': {'Sequence_Name'}  # Additional base field
    },
    'Z-DNA': {
        'core': {'ID', 'Class', 'Subclass', 'Start', 'End', 'Length', 'Sequence', 'Score', 'Strand', 'Method', 'Pattern_ID'},
        'detector_specific': {'GC_Content', 'Contributing_10mers', 'Mean_10mer_Score', 'CG_Dinucleotides', 
                    'AT_Dinucleotides', 'Alternating_CG_Regions', 'Alternating_AT_Regions',
                    'Repeat_Unit', 'Repeat_Count'}
    },
    'Cruciform': {
        'core': {'ID', 'Class', 'Subclass', 'Start', 'End', 'Length', 'Sequence', 'Score', 'Strand', 'Method', 'Pattern_ID'},
        'detector_specific': {'Left_Arm', 'Right_Arm', 'Loop_Seq', 'Arm_Length', 'Loop_Length', 'Stem_Length',
                    'GC_Total', 'GC_Left_Arm', 'GC_Right_Arm', 'GC_Loop', 'Mismatches', 'Match_Fraction', 'DeltaG'}
    },
    'Triplex': {
        'core': {'ID', 'Class', 'Subclass', 'Start', 'End', 'Length', 'Sequence', 'Score', 'Strand', 'Method', 'Pattern_ID'},
        'detector_specific': {'Repeat_Unit', 'Copy_Number'}  # For Sticky DNA subclass
    },
    'i-Motif': {
        'core': {'ID', 'Class', 'Subclass', 'Start', 'End', 'Length', 'Sequence', 'Score', 'Strand', 'Method', 'Pattern_ID'},
        'detector_specific': {'Stems', 'Loops', 'Num_Stems', 'Num_Loops', 'Stem_Lengths', 'Loop_Lengths',
                    'GC_Total', 'GC_Stems', 'Stem_Length', 'Loop_Length'}
    },
    'R-Loop': {
        'core': {'ID', 'Class', 'Subclass', 'Start', 'End', 'Length', 'Sequence', 'Score', 'Strand', 'Method', 'Pattern_ID'},
        # R-Loop specific fields - always included in output (REZ values may be 0 if REZ zone not detected)
        'detector_specific': {'Model', 'RIZ_Length', 'RIZ_Perc_G', 'REZ_Length', 'REZ_Perc_G'}
    },
    'Slipped_DNA': {
        'core': {'ID', 'Class', 'Subclass', 'Start', 'End', 'Length', 'Sequence', 'Score', 'Strand', 'Method', 'Pattern_ID'},
        'detector_specific': {'Repeat_Unit', 'Unit_Size', 'Copy_Number', 'Purity', 'Slippage_Energy_Score'}
    },
    'Curved_DNA': {
        'core': {'ID', 'Class', 'Subclass', 'Start', 'End', 'Length', 'Sequence', 'Score', 'Strand', 'Method', 'Pattern_ID'},
        'detector_specific': {'A_Tracts', 'T_Tracts', 'Num_A_Tracts', 'Num_T_Tracts', 'A_Tract_Lengths', 
                    'T_Tract_Lengths', 'GC_Content', 'AT_Content', 'Center_Positions', 'Tract_Type', 'Tract_Length'}
    },
    'A-philic_DNA': {
        'core': {'ID', 'Class', 'Subclass', 'Start', 'End', 'Length', 'Sequence', 'Score', 'Strand', 'Method', 'Pattern_ID'},
        'detector_specific': set()
    }
}


# ═══════════════════════════════════════════════════════════════════════════════
# COLOR MAPPING TESTS
# ═══════════════════════════════════════════════════════════════════════════════

class TestColorMapping(unittest.TestCase):
    """Tests for color mapping consistency."""
    
    def test_unified_colors_include_all_classes(self):
        """Test that UNIFIED_MOTIF_COLORS includes all valid classes."""
        for class_name in VALID_CLASSES:
            self.assertIn(class_name, UNIFIED_MOTIF_COLORS,
                f"Class '{class_name}' missing from UNIFIED_MOTIF_COLORS")
    
    def test_visualization_colors_match_unified(self):
        """Test that VISUALIZATION_MOTIF_COLORS matches UNIFIED_MOTIF_COLORS."""
        for class_name in VALID_CLASSES:
            unified = UNIFIED_MOTIF_COLORS.get(class_name)
            vis = VISUALIZATION_MOTIF_COLORS.get(class_name)
            self.assertEqual(unified, vis,
                f"Color mismatch for '{class_name}': unified={unified}, vis={vis}")
    
    def test_nature_colors_match_unified(self):
        """Test that NATURE_MOTIF_COLORS matches UNIFIED_MOTIF_COLORS."""
        for class_name in VALID_CLASSES:
            unified = UNIFIED_MOTIF_COLORS.get(class_name)
            nature = NATURE_MOTIF_COLORS.get(class_name)
            self.assertEqual(unified, nature,
                f"Color mismatch for '{class_name}': unified={unified}, nature={nature}")
    
    def test_colors_are_valid_hex(self):
        """Test that all colors are valid hex codes."""
        import re
        hex_pattern = re.compile(r'^#[0-9A-Fa-f]{6}$')
        
        for class_name, color in UNIFIED_MOTIF_COLORS.items():
            self.assertTrue(hex_pattern.match(color),
                f"Invalid hex color for '{class_name}': {color}")


# ═══════════════════════════════════════════════════════════════════════════════
# DETECTOR OUTPUT FORMAT TESTS
# ═══════════════════════════════════════════════════════════════════════════════

class TestDetectorOutputFormat(unittest.TestCase):
    """Tests for detector output format compatibility with visualizations."""
    
    @classmethod
    def setUpClass(cls):
        """Initialize all detectors once."""
        cls.detectors = {
            'G-Quadruplex': GQuadruplexDetector(),
            'Z-DNA': ZDNADetector(),
            'Cruciform': CruciformDetector(),
            'Triplex': TriplexDetector(),
            'i-Motif': IMotifDetector(),
            'R-Loop': RLoopDetector(),
            'Slipped_DNA': SlippedDNADetector(),
            'Curved_DNA': CurvedDNADetector(),
            'A-philic_DNA': APhilicDetector(),
        }
        
        # Map detectors to appropriate test sequences
        cls.test_mappings = {
            'G-Quadruplex': ['G4_TELOMERIC', 'G4_CANONICAL'],
            'Z-DNA': ['ZDNA_EXTENDED', 'EGZ_CGG'],
            'Cruciform': ['CRUCIFORM'],
            'Triplex': ['TRIPLEX_MIRROR', 'STICKY_GAA'],
            'i-Motif': ['IMOTIF'],
            'R-Loop': ['RLOOP'],
            'Slipped_DNA': ['STR_CAG'],
            'Curved_DNA': ['CURVED'],
            'A-philic_DNA': ['APHILIC'],
        }
    
    def test_all_detectors_have_core_fields(self):
        """Test all detectors output required core fields."""
        for class_name, detector in self.detectors.items():
            test_seq_keys = self.test_mappings.get(class_name, ['MIXED'])
            
            for seq_key in test_seq_keys:
                sequence = TEST_SEQUENCES.get(seq_key, TEST_SEQUENCES['MIXED'])
                motifs = detector.detect_motifs(sequence, f"test_{class_name}")
                
                expected_params = DETECTOR_EXPECTED_PARAMS.get(class_name, {})
                core_fields = expected_params.get('core', set())
                
                for motif in motifs:
                    missing_core = core_fields - set(motif.keys())
                    self.assertEqual(len(missing_core), 0,
                        f"{class_name}: Missing core fields {missing_core} in motif")
    
    def test_motif_classes_match_detector(self):
        """Test that motif Class field matches detector class."""
        for class_name, detector in self.detectors.items():
            test_seq_keys = self.test_mappings.get(class_name, ['MIXED'])
            
            for seq_key in test_seq_keys:
                sequence = TEST_SEQUENCES.get(seq_key, TEST_SEQUENCES['MIXED'])
                motifs = detector.detect_motifs(sequence, f"test_{class_name}")
                
                expected_class = detector.get_motif_class_name()
                
                for motif in motifs:
                    self.assertEqual(motif.get('Class'), expected_class,
                        f"{class_name}: Motif class mismatch")
    
    def test_subclasses_belong_to_class(self):
        """Test that all subclasses are valid for their parent class."""
        for class_name, detector in self.detectors.items():
            test_seq_keys = self.test_mappings.get(class_name, ['MIXED'])
            
            for seq_key in test_seq_keys:
                sequence = TEST_SEQUENCES.get(seq_key, TEST_SEQUENCES['MIXED'])
                motifs = detector.detect_motifs(sequence, f"test_{class_name}")
                
                parent_class = detector.get_motif_class_name()
                valid_subclasses = CLASS_TO_SUBCLASSES.get(parent_class, [])
                
                for motif in motifs:
                    subclass = motif.get('Subclass')
                    self.assertIn(subclass, valid_subclasses,
                        f"{class_name}: Invalid subclass '{subclass}' for class '{parent_class}'")


# ═══════════════════════════════════════════════════════════════════════════════
# ADVANCED PARAMETER TESTS
# ═══════════════════════════════════════════════════════════════════════════════

class TestAdvancedParameters(unittest.TestCase):
    """Tests for detector-specific advanced parameters."""
    
    def setUp(self):
        self.detectors = {
            'G-Quadruplex': GQuadruplexDetector(),
            'Z-DNA': ZDNADetector(),
            'Cruciform': CruciformDetector(),
            'Triplex': TriplexDetector(),
            'i-Motif': IMotifDetector(),
            'R-Loop': RLoopDetector(),
            'Slipped_DNA': SlippedDNADetector(),
            'Curved_DNA': CurvedDNADetector(),
            'A-philic_DNA': APhilicDetector(),
        }
    
    def test_gquad_reports_pattern_priority(self):
        """Test G-Quadruplex reports pattern-based information."""
        detector = self.detectors['G-Quadruplex']
        sequence = TEST_SEQUENCES['G4_TELOMERIC']
        motifs = detector.detect_motifs(sequence, "test")
        
        self.assertGreater(len(motifs), 0)
        # G4 subclasses should be descriptive
        for motif in motifs:
            self.assertIn(motif['Subclass'], CLASS_TO_SUBCLASSES['G-Quadruplex'])
    
    def test_zdna_reports_gc_content(self):
        """Test Z-DNA reports GC content and 10-mer information."""
        detector = self.detectors['Z-DNA']
        sequence = TEST_SEQUENCES['ZDNA_EXTENDED']
        motifs = detector.detect_motifs(sequence, "test")
        
        # Z-DNA detection is score-based
        for motif in motifs:
            if motif.get('Subclass') == 'Z-DNA':
                self.assertIn('GC_Content', motif)
                self.assertIn('Contributing_10mers', motif)
    
    def test_cruciform_reports_thermodynamic_data(self):
        """Test Cruciform reports thermodynamic (deltaG) information."""
        detector = self.detectors['Cruciform']
        sequence = TEST_SEQUENCES['CRUCIFORM']
        motifs = detector.detect_motifs(sequence, "test")
        
        for motif in motifs:
            # Cruciform should report structural data
            self.assertIn('Arm_Length', motif)
            self.assertIn('Loop_Length', motif)
            self.assertIn('DeltaG', motif)
    
    def test_triplex_sticky_dna_reports_repeats(self):
        """Test Triplex Sticky DNA reports repeat information."""
        detector = self.detectors['Triplex']
        sequence = TEST_SEQUENCES['STICKY_GAA']
        motifs = detector.detect_motifs(sequence, "test")
        
        self.assertGreater(len(motifs), 0)
        for motif in motifs:
            if motif.get('Subclass') == 'Sticky DNA':
                self.assertIn('Repeat_Unit', motif)
                self.assertIn('Copy_Number', motif)
    
    def test_imotif_reports_stem_loop_structure(self):
        """Test i-Motif reports C-tract stem/loop structure."""
        detector = self.detectors['i-Motif']
        sequence = TEST_SEQUENCES['IMOTIF']
        motifs = detector.detect_motifs(sequence, "test")
        
        self.assertGreater(len(motifs), 0)
        for motif in motifs:
            self.assertIn('Stems', motif)
            self.assertIn('Loops', motif)
            self.assertIn('Num_Stems', motif)
    
    def test_rloop_reports_riz_rez_data(self):
        """Test R-Loop reports RIZ/REZ zone information."""
        detector = self.detectors['R-Loop']
        sequence = TEST_SEQUENCES['RLOOP']
        motifs = detector.detect_motifs(sequence, "test")
        
        for motif in motifs:
            self.assertIn('RIZ_Length', motif)
            self.assertIn('RIZ_Perc_G', motif)
            # REZ may not always be found
            self.assertIn('REZ_Length', motif)
    
    def test_slipped_reports_repeat_metrics(self):
        """Test Slipped DNA reports repeat unit and slippage score."""
        detector = self.detectors['Slipped_DNA']
        sequence = TEST_SEQUENCES['STR_CAG']
        motifs = detector.detect_motifs(sequence, "test")
        
        for motif in motifs:
            self.assertIn('Repeat_Unit', motif)
            self.assertIn('Unit_Size', motif)
            self.assertIn('Copy_Number', motif)
            self.assertIn('Purity', motif)
            self.assertIn('Slippage_Score', motif)  # Changed from Slippage_Energy_Score
    
    def test_curved_reports_tract_data(self):
        """Test Curved DNA reports A/T-tract information."""
        detector = self.detectors['Curved_DNA']
        sequence = TEST_SEQUENCES['CURVED']
        motifs = detector.detect_motifs(sequence, "test")
        
        for motif in motifs:
            self.assertIn('AT_Content', motif)
            self.assertIn('GC_Content', motif)


# ═══════════════════════════════════════════════════════════════════════════════
# DATAFRAME COMPATIBILITY TESTS
# ═══════════════════════════════════════════════════════════════════════════════

class TestDataFrameCompatibility(unittest.TestCase):
    """Tests for pandas DataFrame compatibility of detector outputs."""
    
    @classmethod
    def setUpClass(cls):
        """Run all detectors on mixed sequence and collect motifs."""
        cls.all_motifs = []
        
        detectors = [
            GQuadruplexDetector(),
            ZDNADetector(),
            CruciformDetector(),
            TriplexDetector(),
            IMotifDetector(),
            RLoopDetector(),
            SlippedDNADetector(),
            CurvedDNADetector(),
            APhilicDetector(),
        ]
        
        for detector in detectors:
            motifs = detector.detect_motifs(TEST_SEQUENCES['MIXED'], "test_mixed")
            cls.all_motifs.extend(motifs)
    
    def test_can_create_dataframe(self):
        """Test that motifs can be converted to DataFrame."""
        if not self.all_motifs:
            self.skipTest("No motifs detected in mixed sequence")
        
        df = pd.DataFrame(self.all_motifs)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertGreater(len(df), 0)
    
    def test_core_columns_present(self):
        """Test that core columns are present in DataFrame."""
        if not self.all_motifs:
            self.skipTest("No motifs detected in mixed sequence")
        
        df = pd.DataFrame(self.all_motifs)
        core_columns = {'ID', 'Class', 'Subclass', 'Start', 'End', 'Length', 'Score'}
        
        for col in core_columns:
            self.assertIn(col, df.columns, f"Missing column: {col}")
    
    def test_numeric_columns_are_numeric(self):
        """Test that numeric columns have correct types."""
        if not self.all_motifs:
            self.skipTest("No motifs detected in mixed sequence")
        
        df = pd.DataFrame(self.all_motifs)
        
        numeric_columns = {'Start', 'End', 'Length', 'Score'}
        
        for col in numeric_columns:
            if col in df.columns:
                # Should be convertible to numeric
                try:
                    pd.to_numeric(df[col], errors='raise')
                except ValueError:
                    self.fail(f"Column '{col}' contains non-numeric values")
    
    def test_groupby_class(self):
        """Test that motifs can be grouped by Class."""
        if not self.all_motifs:
            self.skipTest("No motifs detected in mixed sequence")
        
        df = pd.DataFrame(self.all_motifs)
        grouped = df.groupby('Class')
        
        # Should be able to iterate over groups
        for class_name, group_df in grouped:
            self.assertIn(class_name, VALID_CLASSES)
            self.assertGreater(len(group_df), 0)
    
    def test_groupby_subclass(self):
        """Test that motifs can be grouped by Subclass."""
        if not self.all_motifs:
            self.skipTest("No motifs detected in mixed sequence")
        
        df = pd.DataFrame(self.all_motifs)
        grouped = df.groupby('Subclass')
        
        for subclass_name, group_df in grouped:
            self.assertIn(subclass_name, VALID_SUBCLASSES)


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE PANEL COMPATIBILITY TESTS
# ═══════════════════════════════════════════════════════════════════════════════

class TestFigurePanelCompatibility(unittest.TestCase):
    """Tests for detector output compatibility with Nature-ready figure panels."""
    
    @classmethod
    def setUpClass(cls):
        """Collect motifs from all detectors."""
        cls.all_motifs = []
        
        test_seq = TEST_SEQUENCES['G4_TELOMERIC'] + 'N' * 50 + TEST_SEQUENCES['IMOTIF']
        
        detectors = [
            GQuadruplexDetector(),
            ZDNADetector(),
            CruciformDetector(),
            TriplexDetector(),
            IMotifDetector(),
            RLoopDetector(),
            SlippedDNADetector(),
            CurvedDNADetector(),
            APhilicDetector(),
        ]
        
        for detector in detectors:
            motifs = detector.detect_motifs(test_seq, "test_panel")
            cls.all_motifs.extend(motifs)
    
    def test_nested_pie_chart_data(self):
        """Test data structure for nested pie chart (composition)."""
        if not self.all_motifs:
            self.skipTest("No motifs detected")
        
        df = pd.DataFrame(self.all_motifs)
        
        # Should be able to count by Class and Subclass
        class_counts = df.groupby('Class').size()
        subclass_counts = df.groupby(['Class', 'Subclass']).size()
        
        self.assertGreater(len(class_counts), 0)
        self.assertGreater(len(subclass_counts), 0)
    
    def test_manhattan_plot_data(self):
        """Test data structure for Manhattan plot (position)."""
        if not self.all_motifs:
            self.skipTest("No motifs detected")
        
        df = pd.DataFrame(self.all_motifs)
        
        # Need Start, End, Class for Manhattan plot
        required_cols = {'Start', 'End', 'Class'}
        self.assertTrue(required_cols.issubset(df.columns))
        
        # Positions should be plottable
        self.assertTrue(df['Start'].dtype in ['int64', 'float64'])
        self.assertTrue(df['End'].dtype in ['int64', 'float64'])
    
    def test_score_distribution_data(self):
        """Test data structure for score distribution plot."""
        if not self.all_motifs:
            self.skipTest("No motifs detected")
        
        df = pd.DataFrame(self.all_motifs)
        
        # Need Score and Class for score distribution
        self.assertIn('Score', df.columns)
        self.assertIn('Class', df.columns)
        
        # Scores should be numeric and finite
        scores = df['Score'].dropna()
        self.assertTrue(all(isinstance(s, (int, float)) for s in scores))
        self.assertFalse(pd.isna(scores).any())  # No NaN
    
    def test_length_kde_data(self):
        """Test data structure for length KDE plot."""
        if not self.all_motifs:
            self.skipTest("No motifs detected")
        
        df = pd.DataFrame(self.all_motifs)
        
        # Need Length and Class for KDE
        self.assertIn('Length', df.columns)
        self.assertIn('Class', df.columns)
        
        # Lengths should be positive integers
        lengths = df['Length'].dropna()
        self.assertTrue(all(l > 0 for l in lengths))


# ═══════════════════════════════════════════════════════════════════════════════
# METRIC FILTER COMPATIBILITY TESTS
# ═══════════════════════════════════════════════════════════════════════════════

class TestMetricFilterCompatibility(unittest.TestCase):
    """Tests for MetricFilter compatibility with detector outputs."""
    
    def test_core_metrics_present_in_outputs(self):
        """Test that MetricFilter.CORE_METRICS are present in detector outputs."""
        detectors = [
            ('G-Quadruplex', GQuadruplexDetector(), TEST_SEQUENCES['G4_TELOMERIC']),
            ('Z-DNA', ZDNADetector(), TEST_SEQUENCES['ZDNA_EXTENDED']),
            ('Triplex', TriplexDetector(), TEST_SEQUENCES['STICKY_GAA']),
            ('i-Motif', IMotifDetector(), TEST_SEQUENCES['IMOTIF']),
            ('Slipped_DNA', SlippedDNADetector(), TEST_SEQUENCES['STR_CAG']),
            ('Curved_DNA', CurvedDNADetector(), TEST_SEQUENCES['CURVED']),
        ]
        
        for name, detector, sequence in detectors:
            motifs = detector.detect_motifs(sequence, "test")
            
            if motifs:
                motif = motifs[0]
                # Check core metrics are present
                for metric in MetricFilter.CORE_METRICS:
                    if metric == 'Sequence_Name':
                        continue  # Optional
                    self.assertIn(metric, motif,
                        f"{name}: Missing core metric '{metric}'")
    
    def test_conditional_metrics_match_class(self):
        """Test that conditional metrics are class-appropriate."""
        conditional = MetricFilter.CONDITIONAL_METRICS
        
        # G-Quadruplex should not have Repeat_Unit from Slipped_DNA
        g4_metrics = conditional.get('G-Quadruplex', [])
        slipped_metrics = conditional.get('Slipped_DNA', [])
        
        # These should be different
        self.assertNotIn('Repeat_Unit', g4_metrics)
        self.assertIn('Repeat_Unit', slipped_metrics)


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    unittest.main(verbosity=2)
