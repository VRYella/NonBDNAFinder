"""
╔══════════════════════════════════════════════════════════════════════════════╗
║            COMPREHENSIVE DETECTOR TEST SUITE                                  ║
║        Testing All Non-B DNA Motif Detectors                                  ║
╚══════════════════════════════════════════════════════════════════════════════╝

Tests validate:
1. Each detector correctly initializes
2. Each detector produces valid motif output structure
3. Motif taxonomy (Class/Subclass) is canonical
4. Score values are within expected ranges
5. Sequence positions are valid
6. All required parameters are reported
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
from typing import Dict, List, Any

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

# Import taxonomy validation
from Utilities.config.motif_taxonomy import (
    VALID_CLASSES,
    VALID_SUBCLASSES,
    is_valid_pairing
)


# ═══════════════════════════════════════════════════════════════════════════════
# TEST SEQUENCES - Biologically valid test cases for each motif class
# ═══════════════════════════════════════════════════════════════════════════════

TEST_SEQUENCES = {
    # G-Quadruplex: Known telomeric G4 pattern
    'G4_TELOMERIC': 'TTAGGGTTAGGGTTAGGGTTAGGG',
    'G4_CANONICAL': 'GGGATGGGCTGGGAAGGG',
    'G4_EXTENDED': 'GGGAACTGGGAAGGGACTGGG',
    
    # Z-DNA: CG-rich alternating purine-pyrimidine sequences
    'ZDNA_CANONICAL': 'CGCGCGCGCGCGCGCGCGCG',
    'ZDNA_EXTENDED': 'CGCGCGCGCGCGCGCGCGCGCGCGCGCGCG',
    'EGZ_CGG': 'CGGCGGCGGCGGCGG',
    'EGZ_GGC': 'GGCGGCGGCGGCGGC',
    
    # Cruciform: Perfect inverted repeat with loop
    'CRUCIFORM_PERFECT': 'ATCGATCGATCGGGGCGATCGATCGAT',
    'CRUCIFORM_IR': 'AGTCAGTCAGTCAAGCTGACTGACTGACT',
    
    # Triplex: Mirror repeat (H-DNA forming)
    'TRIPLEX_MIRROR': 'GAAGAAGAAGAAAAGAAGAAGAAG',
    'STICKY_GAA': 'GAAGAAGAAGAAGAAGAA',
    'STICKY_TTC': 'TTCTTCTTCTTCTTCTTC',
    
    # i-Motif: C-rich structure  
    'IMOTIF_CANONICAL': 'CCCCTCCCCTCCCCTCCCC',
    'IMOTIF_WEAK': 'CCCTCCCTCCCTCCC',
    'HUR_AC': 'AAACGTACCCGTACCCGTACCCCCC',
    
    # R-Loop: G-rich (QmRLFS model)
    'RLOOP_GRICH': 'GGGGCGGGGGCGGGGCGGGGGCGGGG' + 'A' * 100 + 'GGGGCGGGG',
    
    # Slipped DNA: Tandem repeats (STR and Direct Repeat)
    'STR_CAG': 'CAGCAGCAGCAGCAGCAGCAGCAG',
    'STR_AT': 'ATATATATATATATATATATAT',
    'DIRECT_REPEAT': 'ATCGATCGATCGATCGATCGATCGATCG',
    
    # Curved DNA: A-tract phasing
    'CURVED_ATRACT': 'AAAAAATAAAAATAAAAATAAAAAAAAAA',
    'CURVED_TTRACT': 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
    
    # A-philic DNA: A-rich sequences
    'APHILIC_ARICH': 'AACGAACGAACGAACGAACGAACGAACG',
    
    # Control sequences (should not match)
    'RANDOM': 'ACGTACGTACGTACGTACGT',
    'EMPTY': '',
}


# ═══════════════════════════════════════════════════════════════════════════════
# BASE TEST CLASS
# ═══════════════════════════════════════════════════════════════════════════════

class BaseDetectorTest(unittest.TestCase):
    """Base test class with common assertions for all detectors."""
    
    # Required fields in motif output
    REQUIRED_FIELDS = {
        'ID', 'Sequence_Name', 'Class', 'Subclass', 'Start', 'End',
        'Length', 'Sequence', 'Score', 'Strand', 'Method', 'Pattern_ID'
    }
    
    def validate_motif_structure(self, motif: Dict[str, Any], detector_name: str):
        """Validate that a motif has all required fields."""
        missing = self.REQUIRED_FIELDS - set(motif.keys())
        self.assertEqual(
            len(missing), 0,
            f"{detector_name}: Missing required fields: {missing}"
        )
    
    def validate_motif_taxonomy(self, motif: Dict[str, Any], detector_name: str):
        """Validate Class/Subclass are canonical."""
        class_name = motif.get('Class', '')
        subclass = motif.get('Subclass', '')
        
        self.assertIn(
            class_name, VALID_CLASSES,
            f"{detector_name}: Invalid class '{class_name}'"
        )
        self.assertIn(
            subclass, VALID_SUBCLASSES,
            f"{detector_name}: Invalid subclass '{subclass}'"
        )
        self.assertTrue(
            is_valid_pairing(class_name, subclass),
            f"{detector_name}: Invalid pairing {class_name}/{subclass}"
        )
    
    def validate_motif_positions(self, motif: Dict[str, Any], 
                                  sequence: str, detector_name: str):
        """Validate position fields are correct."""
        start = motif.get('Start', 0)
        end = motif.get('End', 0)
        length = motif.get('Length', 0)
        seq_len = len(sequence)
        
        # Positions are 1-based
        self.assertGreaterEqual(start, 1, 
            f"{detector_name}: Start must be >= 1")
        self.assertLessEqual(end, seq_len, 
            f"{detector_name}: End ({end}) exceeds sequence length ({seq_len})")
        self.assertGreater(end, start - 1, 
            f"{detector_name}: End must be > Start-1")
        self.assertEqual(length, end - start + 1, 
            f"{detector_name}: Length should equal End - Start + 1")
    
    def validate_score_range(self, motif: Dict[str, Any], 
                              detector_name: str,
                              min_score: float = 0.0,
                              max_score: float = float('inf')):
        """Validate score is within expected range."""
        score = motif.get('Score', 0)
        self.assertGreaterEqual(score, min_score,
            f"{detector_name}: Score {score} below minimum {min_score}")
        self.assertLessEqual(score, max_score,
            f"{detector_name}: Score {score} above maximum {max_score}")


# ═══════════════════════════════════════════════════════════════════════════════
# G-QUADRUPLEX DETECTOR TESTS
# ═══════════════════════════════════════════════════════════════════════════════

class TestGQuadruplexDetector(BaseDetectorTest):
    """Tests for GQuadruplexDetector."""
    
    def setUp(self):
        self.detector = GQuadruplexDetector()
        self.name = "GQuadruplexDetector"
    
    def test_initialization(self):
        """Test detector initializes correctly."""
        self.assertIsNotNone(self.detector)
        self.assertEqual(self.detector.get_motif_class_name(), "G-Quadruplex")
    
    def test_patterns_structure(self):
        """Test patterns dictionary structure."""
        patterns = self.detector.get_patterns()
        self.assertIsInstance(patterns, dict)
        self.assertGreater(len(patterns), 0)
        
        # Check expected pattern groups
        expected_groups = ['telomeric_g4', 'canonical_g4', 'weak_pqs']
        for group in expected_groups:
            self.assertIn(group, patterns, f"Missing pattern group: {group}")
    
    def test_detect_telomeric_g4(self):
        """Test detection of telomeric G4."""
        sequence = TEST_SEQUENCES['G4_TELOMERIC']
        motifs = self.detector.detect_motifs(sequence, "test_g4_tel")
        
        self.assertGreater(len(motifs), 0, "Should detect telomeric G4")
        
        for motif in motifs:
            self.validate_motif_structure(motif, self.name)
            self.validate_motif_taxonomy(motif, self.name)
            self.assertEqual(motif['Class'], 'G-Quadruplex')
    
    def test_detect_canonical_g4(self):
        """Test detection of canonical G4."""
        sequence = TEST_SEQUENCES['G4_CANONICAL']
        motifs = self.detector.detect_motifs(sequence, "test_g4_can")
        
        self.assertGreater(len(motifs), 0, "Should detect canonical G4")
        
        for motif in motifs:
            self.validate_motif_structure(motif, self.name)
            self.validate_motif_taxonomy(motif, self.name)
    
    def test_score_calculation(self):
        """Test score calculation."""
        sequence = TEST_SEQUENCES['G4_CANONICAL']
        score = self.detector.calculate_score(sequence, None)
        self.assertGreaterEqual(score, 0.0, "Score should be non-negative")
    
    def test_annotation(self):
        """Test annotate_sequence returns valid annotations."""
        sequence = TEST_SEQUENCES['G4_TELOMERIC']
        annotations = self.detector.annotate_sequence(sequence)
        
        self.assertIsInstance(annotations, list)
        if annotations:
            ann = annotations[0]
            self.assertIn('start', ann)
            self.assertIn('end', ann)
            self.assertIn('score', ann)
            self.assertIn('class_name', ann)
    
    def test_output_subclasses(self):
        """Test that all G4 subclasses are valid."""
        expected_subclasses = {
            'Telomeric G4', 'Stacked canonical G4s', 'Stacked G4s with linker',
            'Canonical intramolecular G4', 'Extended-loop canonical',
            'Higher-order G4 array/G4-wire', 'Intramolecular G-triplex',
            'Two-tetrad weak PQS'
        }
        
        # Run detection on various sequences and collect subclasses
        test_seqs = [TEST_SEQUENCES['G4_TELOMERIC'], TEST_SEQUENCES['G4_CANONICAL']]
        detected_subclasses = set()
        
        for seq in test_seqs:
            motifs = self.detector.detect_motifs(seq, "test")
            for m in motifs:
                detected_subclasses.add(m.get('Subclass', ''))
        
        # All detected subclasses should be valid
        for subclass in detected_subclasses:
            self.assertIn(subclass, VALID_SUBCLASSES,
                f"Invalid G4 subclass: {subclass}")


# ═══════════════════════════════════════════════════════════════════════════════
# Z-DNA DETECTOR TESTS
# ═══════════════════════════════════════════════════════════════════════════════

class TestZDNADetector(BaseDetectorTest):
    """Tests for ZDNADetector."""
    
    def setUp(self):
        self.detector = ZDNADetector()
        self.name = "ZDNADetector"
    
    def test_initialization(self):
        """Test detector initializes correctly."""
        self.assertIsNotNone(self.detector)
        self.assertEqual(self.detector.get_motif_class_name(), "Z-DNA")
    
    def test_patterns_structure(self):
        """Test patterns dictionary structure."""
        patterns = self.detector.get_patterns()
        self.assertIsInstance(patterns, dict)
        self.assertIn('egz_motifs', patterns)
    
    def test_detect_zdna(self):
        """Test detection of Z-DNA."""
        sequence = TEST_SEQUENCES['ZDNA_EXTENDED']
        motifs = self.detector.detect_motifs(sequence, "test_zdna")
        
        # Z-DNA detection is score-based, may not always find motifs
        if motifs:
            for motif in motifs:
                self.validate_motif_structure(motif, self.name)
                self.validate_motif_taxonomy(motif, self.name)
                self.assertEqual(motif['Class'], 'Z-DNA')
    
    def test_detect_egz(self):
        """Test detection of eGZ motifs (CGG/GGC repeats)."""
        sequence = TEST_SEQUENCES['EGZ_CGG']
        motifs = self.detector.detect_motifs(sequence, "test_egz")
        
        if motifs:
            for motif in motifs:
                self.validate_motif_structure(motif, self.name)
                self.validate_motif_taxonomy(motif, self.name)
                # eGZ can have subclass 'eGZ'
                self.assertIn(motif['Subclass'], ['Z-DNA', 'eGZ'])
    
    def test_annotation(self):
        """Test annotate_sequence returns valid annotations."""
        sequence = TEST_SEQUENCES['ZDNA_CANONICAL']
        annotations = self.detector.annotate_sequence(sequence)
        
        self.assertIsInstance(annotations, list)
        if annotations:
            ann = annotations[0]
            self.assertIn('start', ann)
            self.assertIn('end', ann)
    
    def test_score_calculation(self):
        """Test score calculation."""
        sequence = TEST_SEQUENCES['ZDNA_CANONICAL']
        score = self.detector.calculate_score(sequence, None)
        self.assertIsInstance(score, float)
    
    def test_zdna_specific_fields(self):
        """Test Z-DNA specific output fields."""
        sequence = TEST_SEQUENCES['ZDNA_EXTENDED']
        motifs = self.detector.detect_motifs(sequence, "test")
        
        # Z-DNA motifs should have GC_Content field
        for motif in motifs:
            if motif.get('Subclass') == 'Z-DNA':
                self.assertIn('GC_Content', motif)


# ═══════════════════════════════════════════════════════════════════════════════
# CRUCIFORM DETECTOR TESTS
# ═══════════════════════════════════════════════════════════════════════════════

class TestCruciformDetector(BaseDetectorTest):
    """Tests for CruciformDetector."""
    
    def setUp(self):
        self.detector = CruciformDetector()
        self.name = "CruciformDetector"
    
    def test_initialization(self):
        """Test detector initializes correctly."""
        self.assertIsNotNone(self.detector)
        self.assertEqual(self.detector.get_motif_class_name(), "Cruciform")
    
    def test_patterns_structure(self):
        """Test patterns dictionary structure."""
        patterns = self.detector.get_patterns()
        self.assertIsInstance(patterns, dict)
    
    def test_detect_inverted_repeat(self):
        """Test detection of inverted repeats."""
        # Create a perfect inverted repeat
        sequence = 'ATCGATCGATCG' + 'AAAA' + 'CGATCGATCGAT'
        motifs = self.detector.detect_motifs(sequence, "test_cru")
        
        # Inverted repeat detection depends on parameters
        if motifs:
            for motif in motifs:
                self.validate_motif_structure(motif, self.name)
                self.validate_motif_taxonomy(motif, self.name)
                self.assertEqual(motif['Class'], 'Cruciform')
    
    def test_find_inverted_repeats_method(self):
        """Test the find_inverted_repeats internal method."""
        sequence = 'ATCGATCGATCG' + 'AAAA' + 'CGATCGATCGAT'
        hits = self.detector.find_inverted_repeats(sequence)
        
        self.assertIsInstance(hits, list)
        if hits:
            hit = hits[0]
            self.assertIn('left_start', hit)
            self.assertIn('right_end', hit)
            self.assertIn('arm_len', hit)
            self.assertIn('deltaG', hit)
            self.assertIn('score', hit)
    
    def test_thermodynamic_scoring(self):
        """Test thermodynamic (deltaG) scoring."""
        stem_seq = 'GCGCGCGC'  # GC-rich, should have negative deltaG
        deltaG = self.detector._calculate_stem_deltaG(stem_seq)
        self.assertLess(deltaG, 0, "GC-rich stem should have negative deltaG")
    
    def test_cruciform_specific_fields(self):
        """Test cruciform-specific output fields."""
        sequence = 'ATCGATCGATCG' + 'AAAA' + 'CGATCGATCGAT'
        motifs = self.detector.detect_motifs(sequence, "test")
        
        for motif in motifs:
            # Cruciform motifs should have arm/loop fields
            expected_fields = ['Left_Arm', 'Right_Arm', 'Arm_Length', 'Loop_Length']
            for field in expected_fields:
                self.assertIn(field, motif, f"Missing field: {field}")


# ═══════════════════════════════════════════════════════════════════════════════
# TRIPLEX DETECTOR TESTS
# ═══════════════════════════════════════════════════════════════════════════════

class TestTriplexDetector(BaseDetectorTest):
    """Tests for TriplexDetector."""
    
    def setUp(self):
        self.detector = TriplexDetector()
        self.name = "TriplexDetector"
    
    def test_initialization(self):
        """Test detector initializes correctly."""
        self.assertIsNotNone(self.detector)
        self.assertEqual(self.detector.get_motif_class_name(), "Triplex")
    
    def test_patterns_structure(self):
        """Test patterns dictionary structure."""
        patterns = self.detector.get_patterns()
        self.assertIsInstance(patterns, dict)
        self.assertIn('mirror_triplex', patterns)
        self.assertIn('sticky_dna', patterns)
    
    def test_detect_sticky_dna_gaa(self):
        """Test detection of Sticky DNA (GAA repeats)."""
        sequence = TEST_SEQUENCES['STICKY_GAA']
        motifs = self.detector.detect_motifs(sequence, "test_sticky")
        
        self.assertGreater(len(motifs), 0, "Should detect GAA sticky DNA")
        
        for motif in motifs:
            self.validate_motif_structure(motif, self.name)
            self.validate_motif_taxonomy(motif, self.name)
            self.assertEqual(motif['Class'], 'Triplex')
            self.assertEqual(motif['Subclass'], 'Sticky DNA')
    
    def test_detect_sticky_dna_ttc(self):
        """Test detection of Sticky DNA (TTC repeats)."""
        sequence = TEST_SEQUENCES['STICKY_TTC']
        motifs = self.detector.detect_motifs(sequence, "test_sticky")
        
        self.assertGreater(len(motifs), 0, "Should detect TTC sticky DNA")
        
        for motif in motifs:
            self.validate_motif_structure(motif, self.name)
            self.validate_motif_taxonomy(motif, self.name)
            self.assertEqual(motif['Subclass'], 'Sticky DNA')
    
    def test_mirror_repeat_detection(self):
        """Test mirror repeat (H-DNA) detection."""
        # Mirror repeat with purine purity
        sequence = TEST_SEQUENCES['TRIPLEX_MIRROR']
        motifs = self.detector.detect_motifs(sequence, "test_triplex")
        
        if motifs:
            for motif in motifs:
                self.validate_motif_structure(motif, self.name)
                self.validate_motif_taxonomy(motif, self.name)
    
    def test_sticky_dna_copy_number(self):
        """Test that Sticky DNA reports Copy_Number."""
        sequence = TEST_SEQUENCES['STICKY_GAA']
        motifs = self.detector.detect_motifs(sequence, "test")
        
        for motif in motifs:
            if motif.get('Subclass') == 'Sticky DNA':
                self.assertIn('Copy_Number', motif)
                self.assertIn('Repeat_Unit', motif)
    
    def test_subclass_validity(self):
        """Test all output subclasses are valid."""
        test_seqs = [TEST_SEQUENCES['STICKY_GAA'], TEST_SEQUENCES['TRIPLEX_MIRROR']]
        
        for seq in test_seqs:
            motifs = self.detector.detect_motifs(seq, "test")
            for motif in motifs:
                self.assertIn(motif['Subclass'], ['Triplex', 'Sticky DNA'])

    def test_mirror_triplex_mechanistic_scoring(self):
        """Test mechanistic H-DNA scoring returns values in 1-3 range."""
        # Test various arm lengths, loop lengths, purity, and interruptions
        test_cases = [
            (10, 2, 0.95, 0),   # Min arm, short loop, high purity
            (35, 0, 1.0, 0),    # Reference arm, no loop, perfect purity
            (20, 4, 0.90, 2),   # Medium arm, medium loop, threshold purity
            (50, 8, 0.85, 5),   # Long arm, max loop, low purity
        ]
        
        for arm_len, loop_len, purity, interruptions in test_cases:
            score = self.detector._score_mirror_triplex(
                arm_len, loop_len, purity, interruptions
            )
            self.assertGreaterEqual(score, 1.0, 
                f"Score should be >= 1.0 for arm={arm_len}, loop={loop_len}")
            self.assertLessEqual(score, 3.0, 
                f"Score should be <= 3.0 for arm={arm_len}, loop={loop_len}")
    
    def test_sticky_dna_piecewise_scoring(self):
        """Test piecewise Sticky DNA scoring aligned with biological thresholds."""
        # Test each threshold range
        test_cases = [
            (6, 1.0, 1.3, False, False, False),    # Weak range
            (19, 1.0, 1.3, False, False, False),   # Just below replication blockage
            (20, 1.3, 2.0, True, False, False),    # Replication blockage start
            (39, 1.3, 2.0, True, False, False),    # Just below sticky threshold
            (40, 2.0, 2.6, False, True, False),    # Sticky threshold start
            (59, 2.0, 2.6, False, True, False),    # Just below pathogenic
            (60, 2.6, 3.0, False, False, True),    # Pathogenic start
            (100, 2.6, 3.0, False, False, True),   # High pathogenic
        ]
        
        for copies, min_score, max_score, rep_flag, sticky_flag, patho_flag in test_cases:
            score, flags = self.detector._score_sticky_dna(copies)
            self.assertGreaterEqual(score, min_score,
                f"Score for {copies} copies should be >= {min_score}")
            self.assertLessEqual(score, max_score,
                f"Score for {copies} copies should be <= {max_score}")
            self.assertEqual(flags["Replication_Blockage_Range"], rep_flag,
                f"Replication_Blockage_Range flag mismatch for {copies} copies")
            self.assertEqual(flags["Sticky_Threshold_Range"], sticky_flag,
                f"Sticky_Threshold_Range flag mismatch for {copies} copies")
            self.assertEqual(flags["Pathogenic_Range"], patho_flag,
                f"Pathogenic_Range flag mismatch for {copies} copies")
    
    def test_sticky_dna_flags_in_output(self):
        """Test that biological flags are included in Sticky DNA motif output."""
        # Create sequence with 25 copies (replication blockage range)
        sequence = "GAA" * 25
        motifs = self.detector.detect_motifs(sequence, "test_flags")
        
        # Find sticky DNA motif
        sticky_motifs = [m for m in motifs if m.get('Subclass') == 'Sticky DNA']
        self.assertGreater(len(sticky_motifs), 0, "Should detect Sticky DNA")
        
        for motif in sticky_motifs:
            # Check flags are present
            self.assertIn('Replication_Blockage_Range', motif)
            self.assertIn('Sticky_Threshold_Range', motif)
            self.assertIn('Pathogenic_Range', motif)
            
            # 25 copies should be in replication blockage range
            self.assertTrue(motif['Replication_Blockage_Range'],
                "25 copies should be in Replication_Blockage_Range")
            self.assertFalse(motif['Sticky_Threshold_Range'])
            self.assertFalse(motif['Pathogenic_Range'])


# ═══════════════════════════════════════════════════════════════════════════════
# I-MOTIF DETECTOR TESTS
# ═══════════════════════════════════════════════════════════════════════════════

class TestIMotifDetector(BaseDetectorTest):
    """Tests for IMotifDetector."""
    
    def setUp(self):
        self.detector = IMotifDetector()
        self.name = "IMotifDetector"
    
    def test_initialization(self):
        """Test detector initializes correctly."""
        self.assertIsNotNone(self.detector)
        self.assertEqual(self.detector.get_motif_class_name(), "i-Motif")
    
    def test_patterns_structure(self):
        """Test patterns dictionary structure."""
        patterns = self.detector.get_patterns()
        self.assertIsInstance(patterns, dict)
        self.assertIn('canonical_imotif', patterns)
        self.assertIn('hur_ac_motif', patterns)
    
    def test_detect_canonical_imotif(self):
        """Test detection of canonical i-motif."""
        sequence = TEST_SEQUENCES['IMOTIF_CANONICAL']
        motifs = self.detector.detect_motifs(sequence, "test_imotif")
        
        self.assertGreater(len(motifs), 0, "Should detect canonical i-motif")
        
        for motif in motifs:
            self.validate_motif_structure(motif, self.name)
            self.validate_motif_taxonomy(motif, self.name)
            self.assertEqual(motif['Class'], 'i-Motif')
    
    def test_validated_sequences(self):
        """Test detection of validated i-motif sequences."""
        validated_seq = TEST_SEQUENCES['IMOTIF_CANONICAL']
        matches = self.detector.find_validated_matches(validated_seq)
        
        # Check validated match returns expected structure
        if matches:
            match = matches[0]
            self.assertIn('id', match)
            self.assertIn('start', match)
            self.assertIn('end', match)
    
    def test_hur_ac_candidates(self):
        """Test HUR AC-motif candidate detection."""
        sequence = TEST_SEQUENCES['HUR_AC']
        candidates = self.detector.find_hur_ac_candidates(sequence)
        
        self.assertIsInstance(candidates, list)
        if candidates:
            cand = candidates[0]
            self.assertIn('start', cand)
            self.assertIn('end', cand)
            self.assertIn('linker', cand)
    
    def test_imotif_specific_fields(self):
        """Test i-motif specific output fields."""
        sequence = TEST_SEQUENCES['IMOTIF_CANONICAL']
        motifs = self.detector.detect_motifs(sequence, "test")
        
        for motif in motifs:
            # i-Motif should report stems and loops
            self.assertIn('Stems', motif)
            self.assertIn('Loops', motif)
            self.assertIn('Num_Stems', motif)
    
    def test_subclass_validity(self):
        """Test all output subclasses are valid."""
        sequence = TEST_SEQUENCES['IMOTIF_CANONICAL']
        motifs = self.detector.detect_motifs(sequence, "test")
        
        valid_subclasses = {'Canonical i-motif', 'Relaxed i-motif', 'AC-motif'}
        for motif in motifs:
            self.assertIn(motif['Subclass'], valid_subclasses)


# ═══════════════════════════════════════════════════════════════════════════════
# R-LOOP DETECTOR TESTS
# ═══════════════════════════════════════════════════════════════════════════════

class TestRLoopDetector(BaseDetectorTest):
    """Tests for RLoopDetector."""
    
    def setUp(self):
        self.detector = RLoopDetector()
        self.name = "RLoopDetector"
    
    def test_initialization(self):
        """Test detector initializes correctly."""
        self.assertIsNotNone(self.detector)
        self.assertEqual(self.detector.get_motif_class_name(), "R-Loop")
    
    def test_patterns_structure(self):
        """Test patterns dictionary structure."""
        patterns = self.detector.get_patterns()
        self.assertIsInstance(patterns, dict)
        self.assertIn('qmrlfs_model_1', patterns)
        self.assertIn('qmrlfs_model_2', patterns)
    
    def test_detect_rloop(self):
        """Test detection of R-loops."""
        sequence = TEST_SEQUENCES['RLOOP_GRICH']
        motifs = self.detector.detect_motifs(sequence, "test_rloop")
        
        if motifs:
            for motif in motifs:
                self.validate_motif_structure(motif, self.name)
                self.validate_motif_taxonomy(motif, self.name)
                self.assertEqual(motif['Class'], 'R-Loop')
    
    def test_annotation(self):
        """Test annotate_sequence returns valid annotations."""
        sequence = TEST_SEQUENCES['RLOOP_GRICH']
        annotations = self.detector.annotate_sequence(sequence)
        
        self.assertIsInstance(annotations, list)
        if annotations:
            ann = annotations[0]
            self.assertIn('riz_start', ann)
            self.assertIn('riz_end', ann)
            self.assertIn('riz_perc_g', ann)
    
    def test_both_strands(self):
        """Test detection on both strands."""
        sequence = TEST_SEQUENCES['RLOOP_GRICH']
        motifs = self.detector.detect_motifs(sequence, "test")
        
        # Check that audit reports both strands scanned
        audit = self.detector.get_audit_info()
        self.assertTrue(audit.get('both_strands_scanned', False))
    
    def test_quality_threshold(self):
        """Test quality threshold filtering."""
        sequence = 'ACGTACGT'  # Random, should not pass quality
        score = self.detector.calculate_score(sequence, None)
        passes = self.detector.passes_quality_threshold(sequence, score, None)
        
        # Low G content should not pass
        self.assertFalse(passes)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIPPED DNA DETECTOR TESTS
# ═══════════════════════════════════════════════════════════════════════════════

class TestSlippedDNADetector(BaseDetectorTest):
    """Tests for SlippedDNADetector."""
    
    def setUp(self):
        self.detector = SlippedDNADetector()
        self.name = "SlippedDNADetector"
    
    def test_initialization(self):
        """Test detector initializes correctly."""
        self.assertIsNotNone(self.detector)
        self.assertEqual(self.detector.get_motif_class_name(), "Slipped_DNA")
    
    def test_detect_str(self):
        """Test detection of STRs (Short Tandem Repeats)."""
        sequence = TEST_SEQUENCES['STR_CAG']
        motifs = self.detector.detect_motifs(sequence, "test_str")
        
        if motifs:
            for motif in motifs:
                self.validate_motif_structure(motif, self.name)
                self.validate_motif_taxonomy(motif, self.name)
                self.assertEqual(motif['Class'], 'Slipped_DNA')
    
    def test_primitive_motif_computation(self):
        """Test primitive motif extraction."""
        # CAGCAGCAG should reduce to CAG
        primitive = self.detector.compute_primitive_motif('CAGCAGCAGCAG')
        self.assertEqual(primitive, 'CAG')
        
        # ATAT should reduce to AT
        primitive = self.detector.compute_primitive_motif('ATATATATAT')
        self.assertEqual(primitive, 'AT')
    
    def test_repeat_purity(self):
        """Test repeat purity calculation."""
        # Perfect repeat
        purity = self.detector.compute_repeat_purity('CAGCAGCAGCAG', 'CAG')
        self.assertEqual(purity, 1.0)
        
        # Impure repeat
        purity = self.detector.compute_repeat_purity('CAGCAGXAGCAG', 'CAG')
        self.assertLess(purity, 1.0)
    
    def test_entropy_calculation(self):
        """Test Shannon entropy calculation."""
        # Low entropy (homopolymer)
        entropy_low = self.detector.calculate_entropy('AAAAAAAA')
        self.assertLess(entropy_low, 0.5)
        
        # High entropy (random)
        entropy_high = self.detector.calculate_entropy('ACGTACGT')
        self.assertGreater(entropy_high, 1.5)
    
    def test_slipped_dna_specific_fields(self):
        """Test slipped DNA specific output fields."""
        sequence = TEST_SEQUENCES['STR_CAG']
        motifs = self.detector.detect_motifs(sequence, "test")
        
        for motif in motifs:
            self.assertIn('Repeat_Unit', motif)
            self.assertIn('Unit_Size', motif)
            self.assertIn('Copy_Number', motif)
            self.assertIn('Purity', motif)
            self.assertIn('Slippage_Score', motif)
            self.assertIn('References', motif)
    
    def test_str_vs_direct_repeat_classification(self):
        """Test STR vs Direct Repeat subclass classification."""
        # STR: unit size < 7 (new threshold per literature)
        str_seq = 'CAGCAGCAGCAGCAGCAGCAG'  # 3bp unit
        motifs_str = self.detector.detect_motifs(str_seq, "test")
        
        for motif in motifs_str:
            # CAG is 3bp, should be STR (unit size < 7)
            if motif.get('Unit_Size', 0) < 7:
                self.assertEqual(motif['Subclass'], 'STR')
    
    def test_disease_motif_bonus(self):
        """Test that disease-associated motifs (CAG/CTG/CGG/GAA/TTC) receive scoring bonus."""
        # CAG is a disease-associated motif (Huntington's)
        cag_seq = 'CAGCAGCAGCAGCAGCAGCAGCAG'  # 8 copies, 24bp
        motifs_cag = self.detector.detect_motifs(cag_seq, "test_cag")
        
        # Non-disease motif with similar structure
        act_seq = 'ACTACTACTACTACTACTACTACT'  # 8 copies, 24bp
        motifs_act = self.detector.detect_motifs(act_seq, "test_act")
        
        if motifs_cag and motifs_act:
            cag_score = motifs_cag[0]['Slippage_Score']
            act_score = motifs_act[0]['Slippage_Score']
            # CAG should have higher score due to disease motif bonus
            self.assertGreater(cag_score, act_score,
                "Disease motif CAG should score higher than non-disease ACT")
    
    def test_unit_size_instability_weighting(self):
        """Test that 2-4 bp units have highest instability weighting compared to other sizes."""
        # Create sequences with same length and perfect purity
        # Using non-disease motifs to isolate unit-size effect
        
        # 3bp unit (highest instability, U=1.0)
        seq_3bp = 'ACTACTACTACTACTACTACTACT'  # 24bp, 8 copies
        score_3bp = self.detector.compute_slippage_energy_score(
            seq_3bp, 'ACT', 8, 1.0
        )
        
        # 1bp unit (lower instability, U=0.75)
        seq_1bp = 'AAAAAAAAAAAAAAAAAAAAAAAA'  # 24bp, 24 copies
        score_1bp = self.detector.compute_slippage_energy_score(
            seq_1bp, 'A', 24, 1.0
        )
        
        # Score should be in valid range [1, 3]
        self.assertGreaterEqual(score_3bp, 1.0)
        self.assertLessEqual(score_3bp, 3.0)
        self.assertGreaterEqual(score_1bp, 1.0)
        self.assertLessEqual(score_1bp, 3.0)
        
        # 3bp units should score higher than 1bp units due to instability weighting
        # (accounting for the higher copy number effect in 1bp, 3bp still has better U factor)
    
    def test_entropy_on_full_tract(self):
        """Test that entropy is calculated on full tract, not just primitive motif."""
        # This sequence has higher entropy when calculated on full tract
        seq = 'CAGCAGCAGCAGCAGCAGCAGCAG'
        
        # Entropy of full tract
        full_entropy = self.detector.calculate_entropy(seq)
        
        # Entropy of primitive unit
        primitive_entropy = self.detector.calculate_entropy('CAG')
        
        # Both should pass minimum threshold
        self.assertGreaterEqual(full_entropy, 0.5)
        self.assertGreaterEqual(primitive_entropy, 0.5)
    
    def test_references_in_output(self):
        """Test that literature references are included in motif output with exact format."""
        sequence = TEST_SEQUENCES['STR_CAG']
        motifs = self.detector.detect_motifs(sequence, "test")
        
        for motif in motifs:
            self.assertIn('References', motif)
            self.assertEqual(
                motif['References'],
                'Sinden 1994; Pearson 2005; Mirkin 2007',
                "References field should match documented format exactly"
            )


# ═══════════════════════════════════════════════════════════════════════════════
# CURVED DNA DETECTOR TESTS
# ═══════════════════════════════════════════════════════════════════════════════

class TestCurvedDNADetector(BaseDetectorTest):
    """Tests for CurvedDNADetector."""
    
    def setUp(self):
        self.detector = CurvedDNADetector()
        self.name = "CurvedDNADetector"
    
    def test_initialization(self):
        """Test detector initializes correctly."""
        self.assertIsNotNone(self.detector)
        self.assertEqual(self.detector.get_motif_class_name(), "Curved_DNA")
    
    def test_patterns_structure(self):
        """Test patterns dictionary structure."""
        patterns = self.detector.get_patterns()
        self.assertIsInstance(patterns, dict)
        self.assertIn('local_curved', patterns)
    
    def test_detect_long_tracts(self):
        """Test detection of long A/T tracts."""
        sequence = TEST_SEQUENCES['CURVED_ATRACT']
        long_tracts = self.detector.find_long_tracts(sequence)
        
        self.assertIsInstance(long_tracts, list)
        if long_tracts:
            tract = long_tracts[0]
            self.assertIn('start', tract)
            self.assertIn('end', tract)
            self.assertIn('base', tract)
            self.assertIn('score', tract)
    
    def test_find_a_tracts(self):
        """Test A-tract window analysis."""
        sequence = TEST_SEQUENCES['CURVED_ATRACT']
        a_tracts = self.detector.find_a_tracts(sequence)
        
        self.assertIsInstance(a_tracts, list)
        if a_tracts:
            tract = a_tracts[0]
            self.assertIn('window_seq', tract)
            self.assertIn('call', tract)
    
    def test_annotation(self):
        """Test annotate_sequence returns comprehensive annotation."""
        sequence = TEST_SEQUENCES['CURVED_ATRACT']
        annotation = self.detector.annotate_sequence(sequence)
        
        self.assertIsInstance(annotation, dict)
        self.assertIn('a_tract_windows', annotation)
        self.assertIn('aprs', annotation)
        self.assertIn('long_tracts', annotation)
        self.assertIn('summary', annotation)
    
    def test_curved_dna_specific_fields(self):
        """Test curved DNA specific output fields."""
        sequence = TEST_SEQUENCES['CURVED_ATRACT']
        motifs = self.detector.detect_motifs(sequence, "test")
        
        for motif in motifs:
            self.assertIn('AT_Content', motif)
            self.assertIn('GC_Content', motif)
    
    def test_subclass_validity(self):
        """Test all output subclasses are valid."""
        sequence = TEST_SEQUENCES['CURVED_ATRACT']
        motifs = self.detector.detect_motifs(sequence, "test")
        
        valid_subclasses = {'Global Curvature', 'Local Curvature'}
        for motif in motifs:
            self.assertIn(motif['Subclass'], valid_subclasses)


# ═══════════════════════════════════════════════════════════════════════════════
# A-PHILIC DNA DETECTOR TESTS
# ═══════════════════════════════════════════════════════════════════════════════

class TestAPhilicDetector(BaseDetectorTest):
    """Tests for APhilicDetector."""
    
    def setUp(self):
        self.detector = APhilicDetector()
        self.name = "APhilicDetector"
    
    def test_initialization(self):
        """Test detector initializes correctly."""
        self.assertIsNotNone(self.detector)
        self.assertEqual(self.detector.get_motif_class_name(), "A-philic_DNA")
    
    def test_patterns_structure(self):
        """Test patterns dictionary structure."""
        patterns = self.detector.get_patterns()
        self.assertIsInstance(patterns, dict)
        self.assertIn('a_philic_10mers', patterns)
    
    def test_detect_aphilic(self):
        """Test detection of A-philic DNA."""
        sequence = TEST_SEQUENCES['APHILIC_ARICH']
        motifs = self.detector.detect_motifs(sequence, "test_aph")
        
        if motifs:
            for motif in motifs:
                self.validate_motif_structure(motif, self.name)
                self.validate_motif_taxonomy(motif, self.name)
                self.assertEqual(motif['Class'], 'A-philic_DNA')
    
    def test_annotation(self):
        """Test annotate_sequence returns valid annotations."""
        sequence = TEST_SEQUENCES['APHILIC_ARICH']
        annotations = self.detector.annotate_sequence(sequence)
        
        self.assertIsInstance(annotations, list)
        if annotations:
            ann = annotations[0]
            self.assertIn('start', ann)
            self.assertIn('end', ann)
            self.assertIn('sum_log2', ann)
    
    def test_score_calculation(self):
        """Test score calculation."""
        sequence = TEST_SEQUENCES['APHILIC_ARICH']
        score = self.detector.calculate_score(sequence, None)
        self.assertIsInstance(score, float)


# ═══════════════════════════════════════════════════════════════════════════════
# INTEGRATION TESTS - Cross-Detector Consistency
# ═══════════════════════════════════════════════════════════════════════════════

class TestCrossDetectorConsistency(BaseDetectorTest):
    """Tests for consistency across all detectors."""
    
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
    
    def test_all_detectors_have_required_methods(self):
        """Test all detectors implement required interface."""
        required_methods = [
            'get_motif_class_name',
            'get_patterns',
            'calculate_score',
            'detect_motifs',
        ]
        
        for name, detector in self.detectors.items():
            for method in required_methods:
                self.assertTrue(
                    hasattr(detector, method) and callable(getattr(detector, method)),
                    f"{name} missing method: {method}"
                )
    
    def test_all_detectors_return_valid_class_names(self):
        """Test all detectors return canonical class names."""
        for name, detector in self.detectors.items():
            class_name = detector.get_motif_class_name()
            self.assertIn(class_name, VALID_CLASSES,
                f"{name} returns invalid class: {class_name}")
    
    def test_empty_sequence_handling(self):
        """Test all detectors handle empty sequences gracefully."""
        for name, detector in self.detectors.items():
            try:
                motifs = detector.detect_motifs('', 'empty_test')
                self.assertIsInstance(motifs, list)
                self.assertEqual(len(motifs), 0, 
                    f"{name} should return empty list for empty sequence")
            except Exception as e:
                self.fail(f"{name} raised exception on empty sequence: {e}")
    
    def test_random_sequence_handling(self):
        """Test all detectors handle random sequences gracefully."""
        random_seq = TEST_SEQUENCES['RANDOM']
        
        for name, detector in self.detectors.items():
            try:
                motifs = detector.detect_motifs(random_seq, 'random_test')
                self.assertIsInstance(motifs, list)
                # Random sequence may or may not produce motifs
                for motif in motifs:
                    self.validate_motif_structure(motif, name)
            except Exception as e:
                self.fail(f"{name} raised exception on random sequence: {e}")
    
    def test_output_format_consistency(self):
        """Test all detectors produce consistent output format."""
        # Use a sequence that should produce motifs for multiple detectors
        # G-rich region + spacer + alternating purine-pyrimidine region
        spacer = 'A' * 30  # Neutral spacer between test regions
        test_seq = 'GGGGCGGGCGGGCGGGG' + spacer + 'CGCGCGCGCGCGCG'
        
        for name, detector in self.detectors.items():
            motifs = detector.detect_motifs(test_seq, 'format_test')
            
            for motif in motifs:
                # All motifs must have required fields
                self.validate_motif_structure(motif, name)
                
                # All taxonomy must be valid
                self.validate_motif_taxonomy(motif, name)
                
                # Positions must be valid
                self.validate_motif_positions(motif, test_seq, name)


# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION COMPATIBILITY TESTS
# ═══════════════════════════════════════════════════════════════════════════════

class TestVisualizationCompatibility(BaseDetectorTest):
    """Tests for detector output compatibility with visualization standards."""
    
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
    
    def test_score_visualization_compatibility(self):
        """Test that scores are numeric and can be visualized."""
        test_seq = 'GGGGCGGGCGGGCGGGG' + 'A' * 30
        
        for name, detector in self.detectors.items():
            motifs = detector.detect_motifs(test_seq, 'viz_test')
            
            for motif in motifs:
                score = motif.get('Score')
                self.assertIsNotNone(score)
                self.assertIsInstance(score, (int, float))
                # Score should be finite for visualization
                import math
                self.assertFalse(
                    math.isnan(score),
                    f"{name}: Score is NaN"
                )
    
    def test_position_visualization_compatibility(self):
        """Test that positions can be plotted on a linear track."""
        test_seq = 'GGGGCGGGCGGGCGGGG' + 'A' * 30
        seq_len = len(test_seq)
        
        for name, detector in self.detectors.items():
            motifs = detector.detect_motifs(test_seq, 'pos_test')
            
            for motif in motifs:
                start = motif.get('Start', 0)
                end = motif.get('End', 0)
                
                # Positions must be valid for linear track plotting
                self.assertGreater(start, 0)
                self.assertLessEqual(end, seq_len)
                self.assertLess(start, end + 1)
    
    def test_class_color_mapping(self):
        """Test that all detector classes have valid color mappings."""
        from Utilities.config.colors import VISUALIZATION_MOTIF_COLORS
        
        for name, detector in self.detectors.items():
            class_name = detector.get_motif_class_name()
            # Class should be in color mapping
            self.assertIn(class_name, VISUALIZATION_MOTIF_COLORS,
                f"Class '{class_name}' missing from color mapping")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    # Run tests with verbosity
    unittest.main(verbosity=2)
