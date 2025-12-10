"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    COMPREHENSIVE DETECTOR TEST SUITE                          ║
║          Unit Tests for All Non-B DNA Motif Detectors                        ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: test_detectors.py
VERSION: 2024.1
LICENSE: MIT

DESCRIPTION:
    Comprehensive test suite for all motif detector classes ensuring:
    - Accuracy of detection algorithms
    - Performance benchmarks
    - Edge case handling
    - Consistency across detectors
    
TEST COVERAGE:
    ✓ All 9 detector classes
    ✓ Known motif sequences (positive tests)
    ✓ Negative controls (should not detect)
    ✓ Edge cases (empty, very long sequences)
    ✓ Performance benchmarks
    ✓ Output format consistency
"""

import unittest
import time
import random
from typing import List, Dict, Any
import sys
import os

# Import detector classes
from detectors import (
    CurvedDNADetector,
    SlippedDNADetector,
    CruciformDetector,
    RLoopDetector,
    TriplexDetector,
    GQuadruplexDetector,
    IMotifDetector,
    ZDNADetector,
    APhilicDetector
)


class TestDetectorBase(unittest.TestCase):
    """Base class for detector tests with common utilities"""
    
    def assertMotifFormat(self, motifs: List[Dict[str, Any]]):
        """Verify motifs follow standard format"""
        required_fields = ['Start', 'End', 'Class', 'Subclass', 'Score']
        for motif in motifs:
            for field in required_fields:
                self.assertIn(field, motif, f"Motif missing required field: {field}")
            self.assertIsInstance(motif['Start'], int)
            self.assertIsInstance(motif['End'], int)
            self.assertGreater(motif['End'], motif['Start'])
            self.assertGreater(motif['Score'], 0)
    
    def benchmark_detector(self, detector, sequence: str, iterations: int = 10) -> float:
        """Benchmark detector performance"""
        start_time = time.time()
        for _ in range(iterations):
            detector.detect_motifs(sequence)
        elapsed = time.time() - start_time
        avg_time = elapsed / iterations
        bp_per_sec = len(sequence) / avg_time
        return bp_per_sec


class TestCurvedDNADetector(TestDetectorBase):
    """Test suite for Curved DNA (A-tract) detection"""
    
    def setUp(self):
        self.detector = CurvedDNADetector()
    
    def test_known_atract_sequence(self):
        """Test detection of known A-tract sequence"""
        # A-tract: phased A/T runs (5+ bp) spaced ~10 bp apart cause curvature
        sequence = "AAAAATGCATAAAAATGCATAAAAATGCA"
        motifs = self.detector.detect_motifs(sequence)
        self.assertGreater(len(motifs), 0, "Should detect A-tract curvature")
        self.assertMotifFormat(motifs)
        for motif in motifs:
            self.assertEqual(motif['Class'], 'Curved_DNA')
    
    def test_no_curvature_sequence(self):
        """Test that random sequence doesn't produce false positives"""
        sequence = "GCGCTAGCTAGCTAGCTAGCT"
        motifs = self.detector.detect_motifs(sequence)
        # Random sequences might have short runs, but should have low scores
        if motifs:
            for motif in motifs:
                self.assertLess(motif['Score'], 2.0, "Random sequence should have low curvature score")
    
    def test_empty_sequence(self):
        """Test handling of empty sequence"""
        motifs = self.detector.detect_motifs("")
        self.assertEqual(len(motifs), 0)
    
    def test_performance_benchmark(self):
        """Benchmark performance on 10kb sequence"""
        sequence = "ATGC" * 2500  # 10kb
        bp_per_sec = self.benchmark_detector(self.detector, sequence, iterations=5)
        self.assertGreater(bp_per_sec, 1000, f"Performance too low: {bp_per_sec:.0f} bp/s")
        print(f"  Curved DNA: {bp_per_sec:,.0f} bp/s")


class TestGQuadruplexDetector(TestDetectorBase):
    """Test suite for G-Quadruplex detection"""
    
    def setUp(self):
        self.detector = GQuadruplexDetector()
    
    def test_canonical_g4(self):
        """Test canonical G4 pattern: GGG(N1-7)GGG(N1-7)GGG(N1-7)GGG"""
        sequence = "GGGTAGGGTGGGTAGGG"  # Classic G4 motif
        motifs = self.detector.detect_motifs(sequence)
        self.assertGreater(len(motifs), 0, "Should detect canonical G4")
        self.assertMotifFormat(motifs)
        # Verify at least one is classified as G-Quadruplex
        classes = [m['Class'] for m in motifs]
        self.assertIn('G-Quadruplex', classes)
    
    def test_no_g4_in_atrich(self):
        """Test that AT-rich sequences don't produce G4 false positives"""
        sequence = "ATATATATATATATAT"
        motifs = self.detector.detect_motifs(sequence)
        # Should not detect G4 in AT-rich sequence
        g4_motifs = [m for m in motifs if m['Class'] == 'G-Quadruplex']
        self.assertEqual(len(g4_motifs), 0)
    
    def test_performance_benchmark(self):
        """Benchmark G4 detection performance"""
        sequence = "ATGC" * 2500  # 10kb
        bp_per_sec = self.benchmark_detector(self.detector, sequence, iterations=5)
        self.assertGreater(bp_per_sec, 1000, f"Performance too low: {bp_per_sec:.0f} bp/s")
        print(f"  G-Quadruplex: {bp_per_sec:,.0f} bp/s")


class TestCruciformDetector(TestDetectorBase):
    """Test suite for Cruciform (inverted repeat) detection"""
    
    def setUp(self):
        self.detector = CruciformDetector()
    
    def test_perfect_inverted_repeat(self):
        """Test detection of perfect palindromic sequence"""
        # Create perfect inverted repeat: ATGC...GCAT
        arm = "ATGCTAGCTAG"
        sequence = arm + "AA" + self._reverse_complement(arm)
        motifs = self.detector.detect_motifs(sequence)
        self.assertGreater(len(motifs), 0, "Should detect inverted repeat")
        self.assertMotifFormat(motifs)
    
    def test_no_cruciform_in_random(self):
        """Test that random sequence doesn't trigger false positives"""
        sequence = "ATCGATCGATCGATCG"  # No inverted repeat
        motifs = self.detector.detect_motifs(sequence)
        # May detect short imperfect matches, but should have low scores
        if motifs:
            for motif in motifs:
                self.assertLess(motif['Score'], 2.0)
    
    def _reverse_complement(self, seq: str) -> str:
        """Helper: generate reverse complement"""
        comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(comp.get(b, b) for b in reversed(seq))
    
    def test_performance_benchmark(self):
        """Benchmark cruciform detection"""
        sequence = "ATGC" * 2500  # 10kb
        bp_per_sec = self.benchmark_detector(self.detector, sequence, iterations=3)
        self.assertGreater(bp_per_sec, 500, f"Performance too low: {bp_per_sec:.0f} bp/s")
        print(f"  Cruciform: {bp_per_sec:,.0f} bp/s")


class TestSlippedDNADetector(TestDetectorBase):
    """Test suite for Slipped DNA (direct repeat) detection"""
    
    def setUp(self):
        self.detector = SlippedDNADetector()
    
    def test_simple_tandem_repeat(self):
        """Test detection of simple tandem repeat"""
        sequence = "ATGATGATGATGATG"  # ATG repeated 5 times
        motifs = self.detector.detect_motifs(sequence)
        self.assertGreater(len(motifs), 0, "Should detect tandem repeat")
        self.assertMotifFormat(motifs)
    
    def test_str_detection(self):
        """Test Short Tandem Repeat (STR) detection"""
        sequence = "CAG" * 20  # Long CAG repeat (Huntington's disease-associated)
        motifs = self.detector.detect_motifs(sequence)
        self.assertGreater(len(motifs), 0, "Should detect STR")
    
    def test_performance_benchmark(self):
        """Benchmark slipped DNA detection - should be very fast"""
        sequence = "ATGC" * 2500  # 10kb
        bp_per_sec = self.benchmark_detector(self.detector, sequence, iterations=5)
        self.assertGreater(bp_per_sec, 10000, f"Performance too low: {bp_per_sec:.0f} bp/s")
        print(f"  Slipped DNA: {bp_per_sec:,.0f} bp/s")


class TestZDNADetector(TestDetectorBase):
    """Test suite for Z-DNA detection"""
    
    def setUp(self):
        self.detector = ZDNADetector()
    
    def test_alternating_pur_pyr(self):
        """Test Z-DNA detection in alternating purine-pyrimidine"""
        sequence = "CGCGCGCGCGCG"  # Classic Z-DNA forming sequence
        motifs = self.detector.detect_motifs(sequence)
        self.assertGreater(len(motifs), 0, "Should detect Z-DNA potential")
        self.assertMotifFormat(motifs)
    
    def test_performance_benchmark(self):
        """Benchmark Z-DNA detection"""
        sequence = "ATGC" * 2500  # 10kb
        bp_per_sec = self.benchmark_detector(self.detector, sequence, iterations=5)
        self.assertGreater(bp_per_sec, 1000, f"Performance too low: {bp_per_sec:.0f} bp/s")
        print(f"  Z-DNA: {bp_per_sec:,.0f} bp/s")


class TestPerformanceRegression(TestDetectorBase):
    """Test suite to ensure no performance regressions"""
    
    def test_all_detectors_on_large_sequence(self):
        """Test all detectors on a large 50kb sequence"""
        # Generate realistic sequence (40% GC content)
        random.seed(42)
        bases = 'ATGC'
        weights = [0.3, 0.3, 0.2, 0.2]  # A, T, G, C
        sequence = ''.join(random.choices(bases, weights=weights, k=50000))
        
        detectors = {
            'Curved DNA': CurvedDNADetector(),
            'G-Quadruplex': GQuadruplexDetector(),
            'i-Motif': IMotifDetector(),
            'Z-DNA': ZDNADetector(),
            'Slipped DNA': SlippedDNADetector(),
            'Cruciform': CruciformDetector(),
            'Triplex': TriplexDetector(),
            'R-Loop': RLoopDetector(),
            'A-philic': APhilicDetector()
        }
        
        print("\n" + "="*60)
        print("PERFORMANCE BENCHMARK - 50kb sequence")
        print("="*60)
        
        total_time = 0
        for name, detector in detectors.items():
            start = time.time()
            motifs = detector.detect_motifs(sequence)
            elapsed = time.time() - start
            total_time += elapsed
            bp_per_sec = len(sequence) / elapsed if elapsed > 0 else 0
            print(f"{name:15s}: {bp_per_sec:>8,.0f} bp/s | {len(motifs):>4d} motifs | {elapsed:.3f}s")
        
        print("="*60)
        overall_rate = len(sequence) * len(detectors) / total_time if total_time > 0 else 0
        print(f"Overall: {overall_rate:,.0f} bp/s (all detectors)")
        print("="*60)
        
        # Overall performance should be reasonable
        self.assertGreater(overall_rate, 1000, "Overall performance too low")


class TestEdgeCases(TestDetectorBase):
    """Test edge cases and error handling"""
    
    def test_very_short_sequence(self):
        """Test all detectors on very short sequence"""
        sequence = "ATGC"
        detectors = [
            CurvedDNADetector(),
            GQuadruplexDetector(),
            IMotifDetector(),
            ZDNADetector(),
            SlippedDNADetector(),
        ]
        
        for detector in detectors:
            motifs = detector.detect_motifs(sequence)
            # Should not crash, may return 0 or few motifs
            self.assertIsInstance(motifs, list)
            if motifs:
                self.assertMotifFormat(motifs)
    
    def test_very_long_sequence(self):
        """Test detectors on very long sequence (100kb)"""
        sequence = "ATGC" * 25000  # 100kb
        detector = GQuadruplexDetector()  # Test one fast detector
        
        start = time.time()
        motifs = detector.detect_motifs(sequence)
        elapsed = time.time() - start
        
        # Should complete in reasonable time
        self.assertLess(elapsed, 60, "Detection took too long on 100kb sequence")
        self.assertIsInstance(motifs, list)
    
    def test_non_standard_bases(self):
        """Test handling of non-standard bases (N, ambiguous codes)"""
        sequence = "ATGCNNATGCNNATGC"
        detector = GQuadruplexDetector()
        motifs = detector.detect_motifs(sequence)
        # Should not crash
        self.assertIsInstance(motifs, list)


def run_test_suite():
    """Run the complete test suite with detailed output"""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestCurvedDNADetector))
    suite.addTests(loader.loadTestsFromTestCase(TestGQuadruplexDetector))
    suite.addTests(loader.loadTestsFromTestCase(TestCruciformDetector))
    suite.addTests(loader.loadTestsFromTestCase(TestSlippedDNADetector))
    suite.addTests(loader.loadTestsFromTestCase(TestZDNADetector))
    suite.addTests(loader.loadTestsFromTestCase(TestPerformanceRegression))
    suite.addTests(loader.loadTestsFromTestCase(TestEdgeCases))
    
    # Run tests with verbose output
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Print summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Success rate: {(result.testsRun - len(result.failures) - len(result.errors)) / result.testsRun * 100:.1f}%")
    print("="*60)
    
    return result.wasSuccessful()


if __name__ == '__main__':
    success = run_test_suite()
    sys.exit(0 if success else 1)
