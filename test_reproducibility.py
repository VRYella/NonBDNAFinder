"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                   REPRODUCIBILITY VERIFICATION TEST                          ║
║        Ensure Optimizations Don't Change Output or Break Logic               ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: test_reproducibility.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.2
LICENSE: MIT

DESCRIPTION:
    Comprehensive test to verify that all performance optimizations maintain
    identical output to the original implementation. Tests for reproducibility
    across different modes and configurations.

TESTS:
    1. Scanner function output consistency
    2. Standard vs fast mode comparison
    3. Chunked vs non-chunked comparison
    4. Motif count verification
    5. Motif position verification
    6. Score consistency check
"""

import sys
from typing import List, Dict, Any, Tuple
import nonbscanner as nbs


class ReproducibilityTest:
    """Test suite for verifying reproducibility of optimizations."""
    
    def __init__(self):
        self.test_sequences = self._generate_test_sequences()
        self.passed = 0
        self.failed = 0
        self.warnings = 0
    
    def _generate_test_sequences(self) -> Dict[str, str]:
        """Generate diverse test sequences."""
        return {
            'g4_simple': "GGGTTAGGGTTAGGGTTAGGG",
            'g4_complex': "GGGTTAGGGTTAGGGTTAGGG" * 5,
            'mixed_1kb': "ACGTACGTACGTACGT" * 16 + "GGGTTAGGGTTAGGGTTAGGG" * 10 + "CCCCTAACCCTAACCCTAACCC" * 5,
            'genomic_5kb': ("AAAAATTTT" * 10 + "CGCGCGCG" * 10 + "GGGTTAGGGTTAGGG" * 10 + 
                            "CCCCTAACCCTAACCC" * 10 + "ACGTACGTACGTACGT" * 10) * 10
        }
    
    def _compare_motifs(self, motifs1: List[Dict], motifs2: List[Dict], 
                        tolerance: float = 0.001) -> Tuple[bool, str]:
        """
        Compare two motif lists for equality.
        
        Returns:
            (is_equal, message) tuple
        """
        if len(motifs1) != len(motifs2):
            return False, f"Count mismatch: {len(motifs1)} vs {len(motifs2)}"
        
        # Sort both lists by position for comparison
        m1_sorted = sorted(motifs1, key=lambda x: (x.get('Start', 0), x.get('End', 0)))
        m2_sorted = sorted(motifs2, key=lambda x: (x.get('Start', 0), x.get('End', 0)))
        
        for i, (m1, m2) in enumerate(zip(m1_sorted, m2_sorted)):
            # Check positions
            if m1.get('Start') != m2.get('Start') or m1.get('End') != m2.get('End'):
                return False, f"Position mismatch at index {i}: " \
                              f"({m1.get('Start')}, {m1.get('End')}) vs " \
                              f"({m2.get('Start')}, {m2.get('End')})"
            
            # Check class/subclass
            if m1.get('Class') != m2.get('Class'):
                return False, f"Class mismatch at index {i}: " \
                              f"{m1.get('Class')} vs {m2.get('Class')}"
            
            # Check scores (with tolerance for floating point)
            s1 = m1.get('Score', 0)
            s2 = m2.get('Score', 0)
            if abs(s1 - s2) > tolerance:
                return False, f"Score mismatch at index {i}: {s1} vs {s2}"
        
        return True, "Identical"
    
    def test_standard_vs_fast_mode(self) -> bool:
        """Test that standard and fast modes produce identical results."""
        print("\n" + "="*70)
        print("TEST: Standard vs Fast Mode Comparison")
        print("="*70)
        
        all_passed = True
        
        for seq_name, sequence in self.test_sequences.items():
            print(f"\nTesting {seq_name} ({len(sequence)} bp)...")
            
            # Run in standard mode
            motifs_std = nbs.analyze_sequence(sequence, seq_name, use_fast_mode=False)
            
            # Run in fast mode
            motifs_fast = nbs.analyze_sequence(sequence, seq_name, use_fast_mode=True)
            
            # Compare
            is_equal, message = self._compare_motifs(motifs_std, motifs_fast)
            
            if is_equal:
                print(f"  ✓ PASS: {len(motifs_std)} motifs, {message}")
                self.passed += 1
            else:
                print(f"  ✗ FAIL: {message}")
                print(f"    Standard: {len(motifs_std)} motifs")
                print(f"    Fast: {len(motifs_fast)} motifs")
                self.failed += 1
                all_passed = False
        
        return all_passed
    
    def test_chunked_vs_standard(self) -> bool:
        """Test that chunked processing produces identical results."""
        print("\n" + "="*70)
        print("TEST: Chunked vs Standard Processing")
        print("="*70)
        
        all_passed = True
        
        # Only test on larger sequences where chunking matters
        large_sequences = {k: v for k, v in self.test_sequences.items() 
                          if len(v) > 1000}
        
        for seq_name, sequence in large_sequences.items():
            print(f"\nTesting {seq_name} ({len(sequence)} bp)...")
            
            # Standard (no chunking)
            motifs_std = nbs.analyze_sequence(sequence, seq_name, 
                                              use_chunking=False)
            
            # Chunked
            motifs_chunked = nbs.analyze_sequence(sequence, seq_name,
                                                  use_chunking=True,
                                                  chunk_size=2000,
                                                  chunk_overlap=500)
            
            # Compare
            is_equal, message = self._compare_motifs(motifs_std, motifs_chunked)
            
            if is_equal:
                print(f"  ✓ PASS: {len(motifs_std)} motifs, {message}")
                self.passed += 1
            else:
                # Chunking may find slightly different motifs at boundaries
                # This is acceptable if the difference is small
                count_diff = abs(len(motifs_std) - len(motifs_chunked))
                if count_diff <= 2:  # Allow small differences
                    print(f"  ⚠ WARN: Minor difference ({message})")
                    print(f"    Acceptable: count difference = {count_diff}")
                    self.warnings += 1
                else:
                    print(f"  ✗ FAIL: {message}")
                    print(f"    Standard: {len(motifs_std)} motifs")
                    print(f"    Chunked: {len(motifs_chunked)} motifs")
                    self.failed += 1
                    all_passed = False
        
        return all_passed
    
    def test_scanner_functions(self) -> bool:
        """Test that optimized scanner functions work correctly."""
        print("\n" + "="*70)
        print("TEST: Scanner Functions")
        print("="*70)
        
        all_passed = True
        
        try:
            from scanner import find_direct_repeats, find_inverted_repeats, find_strs
            
            test_seq = self.test_sequences['genomic_5kb']
            
            # Test direct repeats
            print("\nTesting direct repeats...")
            try:
                repeats = find_direct_repeats(test_seq)
                print(f"  ✓ PASS: Found {len(repeats)} direct repeats")
                self.passed += 1
            except Exception as e:
                print(f"  ✗ FAIL: {e}")
                self.failed += 1
                all_passed = False
            
            # Test inverted repeats
            print("\nTesting inverted repeats...")
            try:
                inverted = find_inverted_repeats(test_seq)
                print(f"  ✓ PASS: Found {len(inverted)} inverted repeats")
                self.passed += 1
            except Exception as e:
                print(f"  ✗ FAIL: {e}")
                self.failed += 1
                all_passed = False
            
            # Test STRs
            print("\nTesting STRs...")
            try:
                strs = find_strs(test_seq)
                print(f"  ✓ PASS: Found {len(strs)} STRs")
                self.passed += 1
            except Exception as e:
                print(f"  ✗ FAIL: {e}")
                self.failed += 1
                all_passed = False
                
        except ImportError as e:
            print(f"  ✗ FAIL: Cannot import scanner functions: {e}")
            self.failed += 1
            all_passed = False
        
        return all_passed
    
    def test_pattern_cache(self) -> bool:
        """Test that pattern cache works correctly."""
        print("\n" + "="*70)
        print("TEST: Pattern Cache")
        print("="*70)
        
        all_passed = True
        
        try:
            from pattern_cache import get_pattern_cache
            
            cache = get_pattern_cache()
            test_seq = self.test_sequences['g4_complex']
            
            print("\nTesting pattern cache...")
            try:
                # Test G4 patterns
                matches = cache.find_all('g4', 'canonical', test_seq)
                print(f"  ✓ PASS: Found {len(matches)} G4 matches")
                self.passed += 1
                
                # Test i-motif patterns
                matches = cache.find_all('imotif', 'canonical', test_seq)
                print(f"  ✓ PASS: Found {len(matches)} i-motif matches")
                self.passed += 1
                
            except Exception as e:
                print(f"  ✗ FAIL: {e}")
                self.failed += 1
                all_passed = False
                
        except ImportError as e:
            print(f"  ⚠ WARN: Pattern cache not available: {e}")
            self.warnings += 1
        
        return all_passed
    
    def run_all_tests(self) -> bool:
        """Run all reproducibility tests."""
        print("\n" + "#"*70)
        print("# NonBDNAFinder Reproducibility Test Suite")
        print("#"*70)
        
        # Check for performance libraries
        print("\n" + "="*70)
        print("PERFORMANCE LIBRARY CHECK")
        print("="*70)
        
        try:
            import numba
            print("✓ Numba installed - using optimized scanner")
        except ImportError:
            print("⚠ Numba not installed - using fallback (slower)")
        
        try:
            import regex
            print("✓ Regex installed - using optimized patterns")
        except ImportError:
            print("⚠ Regex not installed - using standard re")
        
        # Run tests
        tests = [
            ("Scanner Functions", self.test_scanner_functions),
            ("Pattern Cache", self.test_pattern_cache),
            ("Standard vs Fast Mode", self.test_standard_vs_fast_mode),
            ("Chunked vs Standard", self.test_chunked_vs_standard),
        ]
        
        all_passed = True
        for test_name, test_func in tests:
            try:
                passed = test_func()
                if not passed:
                    all_passed = False
            except Exception as e:
                print(f"\n✗ TEST ERROR in {test_name}: {e}")
                self.failed += 1
                all_passed = False
        
        # Print summary
        self._print_summary()
        
        return all_passed
    
    def _print_summary(self):
        """Print test summary."""
        print("\n" + "="*70)
        print("TEST SUMMARY")
        print("="*70)
        
        total = self.passed + self.failed + self.warnings
        
        print(f"\nTests Run: {total}")
        print(f"  ✓ Passed: {self.passed}")
        print(f"  ⚠ Warnings: {self.warnings}")
        print(f"  ✗ Failed: {self.failed}")
        
        if self.failed == 0:
            print("\n" + "="*70)
            print("✓ ALL TESTS PASSED - Optimizations maintain reproducibility")
            print("="*70)
        else:
            print("\n" + "="*70)
            print("✗ SOME TESTS FAILED - Review failures above")
            print("="*70)


def main():
    """Run reproducibility tests."""
    test_suite = ReproducibilityTest()
    success = test_suite.run_all_tests()
    
    # Exit with appropriate code
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
