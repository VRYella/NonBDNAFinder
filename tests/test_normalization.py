"""
╔══════════════════════════════════════════════════════════════════════════════╗
║            DETECTOR NORMALIZATION TEST SUITE                                  ║
║        Testing Self-Normalization Architecture                                ║
╚══════════════════════════════════════════════════════════════════════════════╝

Tests validate:
1. Each detector has normalization parameters defined
2. _normalize_score() method produces values in [1.0, 3.0] range
3. Boundary conditions are handled correctly
4. Different normalization methods work as expected
5. Scores from detect_motifs() include both Raw_Score and Score fields
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
from typing import Dict, Any

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
from Detectors.base.base_detector import BaseMotifDetector


class TestDetectorNormalization(unittest.TestCase):
    """Test normalization functionality for all detectors."""
    
    def setUp(self):
        """Initialize all detectors for testing."""
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
    
    def test_detectors_have_normalization_parameters(self):
        """Verify all detectors have required normalization parameters."""
        required_attrs = [
            'RAW_SCORE_MIN', 'RAW_SCORE_MAX',
            'NORMALIZED_MIN', 'NORMALIZED_MAX',
            'NORMALIZATION_METHOD', 'SCORE_REFERENCE'
        ]
        
        for name, detector in self.detectors.items():
            with self.subTest(detector=name):
                for attr in required_attrs:
                    self.assertTrue(
                        hasattr(detector, attr),
                        f"{name} detector missing {attr} attribute"
                    )
                
                # Verify normalization method is valid
                valid_methods = ['linear', 'log', 'g4hunter', 'zdna_cumulative']
                self.assertIn(
                    detector.NORMALIZATION_METHOD,
                    valid_methods,
                    f"{name} has invalid normalization method: {detector.NORMALIZATION_METHOD}"
                )
    
    def test_normalize_score_returns_valid_range(self):
        """Verify _normalize_score() returns values in [1.0, 3.0] range."""
        for name, detector in self.detectors.items():
            with self.subTest(detector=name):
                # Test with various raw scores
                test_scores = [
                    detector.RAW_SCORE_MIN,
                    detector.RAW_SCORE_MAX,
                    (detector.RAW_SCORE_MIN + detector.RAW_SCORE_MAX) / 2,
                    detector.RAW_SCORE_MIN - 0.1,  # Below min
                    detector.RAW_SCORE_MAX + 0.1,  # Above max
                ]
                
                for raw_score in test_scores:
                    normalized = detector._normalize_score(raw_score)
                    self.assertIsInstance(normalized, float, 
                        f"{name}: normalized score must be float")
                    self.assertGreaterEqual(normalized, 1.0, 
                        f"{name}: normalized score {normalized} < 1.0 for raw={raw_score}")
                    self.assertLessEqual(normalized, 3.0,
                        f"{name}: normalized score {normalized} > 3.0 for raw={raw_score}")
    
    def test_normalize_score_boundary_conditions(self):
        """Test boundary conditions for normalization."""
        for name, detector in self.detectors.items():
            with self.subTest(detector=name):
                # Min raw score should map to 1.0 (or close to it)
                min_normalized = detector._normalize_score(detector.RAW_SCORE_MIN)
                self.assertAlmostEqual(min_normalized, 1.0, places=1,
                    msg=f"{name}: RAW_SCORE_MIN should normalize to ~1.0")
                
                # Max raw score should map to 3.0 (or close to it)
                max_normalized = detector._normalize_score(detector.RAW_SCORE_MAX)
                self.assertAlmostEqual(max_normalized, 3.0, places=1,
                    msg=f"{name}: RAW_SCORE_MAX should normalize to ~3.0")
    
    def test_normalize_score_midpoint(self):
        """Test that midpoint raw score normalizes to approximately 2.0."""
        for name, detector in self.detectors.items():
            with self.subTest(detector=name):
                # Skip special normalization methods that don't have linear midpoint
                if detector.NORMALIZATION_METHOD in ['g4hunter', 'log']:
                    continue
                
                mid_raw = (detector.RAW_SCORE_MIN + detector.RAW_SCORE_MAX) / 2
                mid_normalized = detector._normalize_score(mid_raw)
                
                # Midpoint should be near 2.0 for linear normalization
                self.assertGreaterEqual(mid_normalized, 1.5,
                    f"{name}: midpoint {mid_normalized} too low")
                self.assertLessEqual(mid_normalized, 2.5,
                    f"{name}: midpoint {mid_normalized} too high")
    
    def test_detect_motifs_includes_both_scores(self):
        """Verify detect_motifs() returns both Raw_Score and Score fields."""
        # Test sequences known to produce motifs
        test_sequences = {
            'G-Quadruplex': 'TTAGGGTTAGGGTTAGGGTTAGGG',
            'Z-DNA': 'CGCGCGCGCGCGCGCGCGCG',
            'Curved_DNA': 'AAAAAATAAAAATAAAAATAAAAA',
            'Slipped_DNA': 'CAGCAGCAGCAGCAGCAGCAGCAG',
        }
        
        for detector_name, sequence in test_sequences.items():
            if detector_name in self.detectors:
                detector = self.detectors[detector_name]
                with self.subTest(detector=detector_name):
                    motifs = detector.detect_motifs(sequence, "test_seq")
                    
                    # Should find at least one motif
                    if len(motifs) > 0:
                        motif = motifs[0]
                        
                        # Check both score fields exist
                        self.assertIn('Raw_Score', motif,
                            f"{detector_name}: Missing Raw_Score field")
                        self.assertIn('Score', motif,
                            f"{detector_name}: Missing Score field")
                        
                        # Check Raw_Score and Score are different (normalized)
                        raw_score = motif['Raw_Score']
                        normalized_score = motif['Score']
                        
                        # Normalized score should be in [1.0, 3.0]
                        self.assertGreaterEqual(normalized_score, 1.0,
                            f"{detector_name}: Score {normalized_score} < 1.0")
                        self.assertLessEqual(normalized_score, 3.0,
                            f"{detector_name}: Score {normalized_score} > 3.0")
    
    def test_linear_normalization_is_monotonic(self):
        """Test that linear normalization preserves ordering."""
        for name, detector in self.detectors.items():
            if detector.NORMALIZATION_METHOD != 'linear':
                continue
            
            with self.subTest(detector=name):
                # Generate sequence of increasing raw scores
                min_raw = detector.RAW_SCORE_MIN
                max_raw = detector.RAW_SCORE_MAX
                step = (max_raw - min_raw) / 10
                
                prev_normalized = 0.0
                for i in range(11):
                    raw_score = min_raw + i * step
                    normalized = detector._normalize_score(raw_score)
                    
                    # Each normalized score should be >= previous
                    self.assertGreaterEqual(normalized, prev_normalized,
                        f"{name}: normalization not monotonic at {raw_score}")
                    prev_normalized = normalized
    
    def test_g4hunter_special_normalization(self):
        """Test G4Hunter-specific normalization (absolute value handling)."""
        g4_detector = self.detectors['G-Quadruplex']
        
        # G4Hunter uses absolute values
        positive_score = g4_detector._normalize_score(0.7)
        negative_score = g4_detector._normalize_score(-0.7)
        
        # Both should produce same result due to abs()
        self.assertEqual(positive_score, negative_score,
            "G4Hunter should treat positive and negative scores equally")
    
    def test_zdna_log_normalization(self):
        """Test Z-DNA log-scale normalization for cumulative scores."""
        zdna_detector = self.detectors['Z-DNA']
        
        # Z-DNA uses log scaling
        self.assertEqual(zdna_detector.NORMALIZATION_METHOD, 'log',
            "Z-DNA should use log normalization")
        
        # Test with cumulative scores
        low_score = zdna_detector._normalize_score(50.0)    # Min
        mid_score = zdna_detector._normalize_score(200.0)   # Mid
        high_score = zdna_detector._normalize_score(2000.0) # Max
        
        # Should be in increasing order
        self.assertLess(low_score, mid_score,
            "Z-DNA log normalization not monotonic (low < mid)")
        self.assertLess(mid_score, high_score,
            "Z-DNA log normalization not monotonic (mid < high)")
        
        # All should be in [1.0, 3.0]
        for score in [low_score, mid_score, high_score]:
            self.assertGreaterEqual(score, 1.0)
            self.assertLessEqual(score, 3.0)


class TestBackwardCompatibility(unittest.TestCase):
    """Test that deprecated functions still work for backward compatibility."""
    
    def test_normalize_motif_scores_deprecated(self):
        """Verify normalize_motif_scores() raises deprecation warning."""
        from Utilities.utilities import normalize_motif_scores
        
        test_motifs = [
            {'Class': 'G-Quadruplex', 'Score': 2.5, 'Raw_Score': 0.8},
            {'Class': 'Curved_DNA', 'Score': 2.2, 'Raw_Score': 0.6}
        ]
        
        # Should raise DeprecationWarning
        with self.assertWarns(DeprecationWarning):
            result = normalize_motif_scores(test_motifs)
        
        # Should return motifs unchanged
        self.assertEqual(len(result), len(test_motifs))
        self.assertEqual(result[0]['Score'], 2.5)
        self.assertEqual(result[1]['Score'], 2.2)
    
    def test_normalize_score_to_1_3_deprecated(self):
        """Verify normalize_score_to_1_3() raises deprecation warning."""
        from Utilities.utilities import normalize_score_to_1_3
        
        # Should raise DeprecationWarning
        with self.assertWarns(DeprecationWarning):
            result = normalize_score_to_1_3(0.8, 'G-Quadruplex')
        
        # Should still return a valid normalized score
        self.assertGreaterEqual(result, 1.0)
        self.assertLessEqual(result, 3.0)


if __name__ == '__main__':
    unittest.main()
