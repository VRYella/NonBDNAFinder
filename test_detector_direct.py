#!/usr/bin/env python3
"""
Direct Detector Test for Slipped DNA (No Pipeline Dependencies)
===============================================================

Tests the SlippedDNADetector directly without requiring the full
NonBScanner pipeline or external dependencies.
"""

import sys
sys.path.insert(0, '/home/runner/work/NonBDNAFinder/NonBDNAFinder')

from detectors import SlippedDNADetector

def test_real_world_sequences():
    """Test with real-world disease-relevant sequences."""
    print("\n" + "="*70)
    print("REAL-WORLD SEQUENCE TESTING")
    print("="*70)
    
    detector = SlippedDNADetector()
    
    test_cases = [
        {
            'name': "Huntington's Disease (CAG)20 (Normal)",
            'sequence': "ATGAAGGCCTTC" + "CAG" * 20 + "CCGCCGCCGCCG",
            'expected_unit': 'CAG',
            'expected_copies': 20,
        },
        {
            'name': "Huntington's Disease (CAG)50 (Expanded)",
            'sequence': "ATGAAGGCCTTC" + "CAG" * 50 + "CCGCCGCCGCCG",
            'expected_unit': 'CAG',
            'expected_copies': 50,
        },
        {
            'name': "Fragile X (CGG)30",
            'sequence': "GCGCGCGC" + "CGG" * 30 + "GCATGCATGCAT",
            'expected_unit': 'CGG',
            'expected_copies': 30,
        },
        {
            'name': "Myotonic Dystrophy (CTG)25",
            'sequence': "ATGCATGC" + "CTG" * 25 + "GCGCGCGC",
            'expected_unit': 'CTG',
            'expected_copies': 25,
        },
    ]
    
    passed = 0
    for test in test_cases:
        print(f"\nTest: {test['name']}")
        print(f"  Sequence length: {len(test['sequence'])} bp")
        
        motifs = detector.detect_motifs(test['sequence'], test['name'])
        
        if motifs:
            motif = motifs[0]  # Should only be one
            unit = motif.get('Repeat_Unit', 'N/A')
            copies = motif.get('Copy_Number', 0)
            score = motif.get('Slippage_Energy_Score', 0)
            purity = motif.get('Purity', 0)
            
            unit_match = unit == test['expected_unit']
            copies_match = copies == test['expected_copies']
            
            print(f"  Detected: ({unit})×{copies}, Purity={purity:.3f}, Score={score:.2f}")
            
            if unit_match and copies_match:
                print(f"  ✓ PASS: Correct unit and copy number")
                passed += 1
            else:
                print(f"  ✗ FAIL: Expected ({test['expected_unit']})×{test['expected_copies']}")
        else:
            print(f"  ✗ FAIL: No motifs detected")
    
    print(f"\nResult: {passed}/{len(test_cases)} tests passed")
    return passed == len(test_cases)

def test_edge_cases():
    """Test edge cases and boundary conditions."""
    print("\n" + "="*70)
    print("EDGE CASE TESTING")
    print("="*70)
    
    detector = SlippedDNADetector()
    
    test_cases = [
        {
            'name': "Exactly 20bp (minimum)",
            'sequence': "CAG" * 7 + "C",  # 21 bp (just above minimum)
            'should_detect': True,
        },
        {
            'name': "Just below 20bp",
            'sequence': "CAG" * 6,  # 18 bp (below minimum)
            'should_detect': False,
        },
        {
            'name': "Interrupted repeat (low purity)",
            'sequence': "CAGCAGCATCAGCAGCAGCAGCAGCAG",  # One interruption
            'should_detect': False,  # Should fail purity check
        },
        {
            'name': "Perfect long repeat",
            'sequence': "GATA" * 20,  # 80 bp, perfect repeat
            'should_detect': True,
        },
        {
            'name': "Homopolymer (low entropy)",
            'sequence': "A" * 50,  # Should fail entropy check
            'should_detect': False,
        },
    ]
    
    passed = 0
    for test in test_cases:
        print(f"\nTest: {test['name']}")
        print(f"  Sequence: {test['sequence'][:30]}..." if len(test['sequence']) > 30 else f"  Sequence: {test['sequence']}")
        
        motifs = detector.detect_motifs(test['sequence'], test['name'])
        detected = len(motifs) > 0
        
        if detected == test['should_detect']:
            status = "detected" if detected else "rejected"
            print(f"  ✓ PASS: Correctly {status}")
            passed += 1
        else:
            expected = "detection" if test['should_detect'] else "rejection"
            actual = "detected" if detected else "rejected"
            print(f"  ✗ FAIL: Expected {expected}, but {actual}")
    
    print(f"\nResult: {passed}/{len(test_cases)} tests passed")
    return passed == len(test_cases)

def test_performance():
    """Test performance on longer sequences."""
    print("\n" + "="*70)
    print("PERFORMANCE TESTING")
    print("="*70)
    
    import time
    
    detector = SlippedDNADetector()
    
    # Generate test sequences of varying sizes
    test_sizes = [100, 1000, 10000]
    
    for size in test_sizes:
        # Create sequence with periodic repeats
        repeat_unit = "CAGCAGCAG"
        num_repeats = size // len(repeat_unit)
        test_seq = repeat_unit * num_repeats
        
        print(f"\nTesting {len(test_seq)} bp sequence...")
        
        start_time = time.time()
        motifs = detector.detect_motifs(test_seq, f"perf_test_{size}")
        elapsed = time.time() - start_time
        
        throughput = len(test_seq) / elapsed if elapsed > 0 else 0
        
        print(f"  Time: {elapsed:.3f}s")
        print(f"  Throughput: {throughput:,.0f} bp/s")
        print(f"  Motifs detected: {len(motifs)}")
        
        if elapsed < 5.0:  # Should be fast even for 10kb
            print(f"  ✓ PASS: Acceptable performance")
        else:
            print(f"  ⚠ WARNING: Slower than expected")
    
    return True  # Performance test always passes (just informative)

def test_output_format():
    """Verify output format is publication-ready."""
    print("\n" + "="*70)
    print("OUTPUT FORMAT VALIDATION")
    print("="*70)
    
    detector = SlippedDNADetector()
    
    test_seq = "CAG" * 20  # Simple test case
    motifs = detector.detect_motifs(test_seq, "format_test")
    
    if not motifs:
        print("✗ FAIL: No motifs detected")
        return False
    
    motif = motifs[0]
    
    # Check for required publication fields
    required_fields = {
        'Class': str,
        'Subclass': str,
        'Start': int,
        'End': int,
        'Length': int,
        'Sequence': str,
        'Repeat_Unit': str,
        'Unit_Size': int,
        'Copy_Number': (int, float),
        'Purity': float,
        'Slippage_Energy_Score': float,
    }
    
    print("\nChecking required fields:")
    all_present = True
    for field, expected_type in required_fields.items():
        if field in motif:
            value = motif[field]
            if isinstance(expected_type, tuple):
                type_ok = isinstance(value, expected_type)
            else:
                type_ok = isinstance(value, expected_type)
            
            status = "✓" if type_ok else "✗"
            print(f"  {status} {field}: {value} ({type(value).__name__})")
            
            if not type_ok:
                all_present = False
        else:
            print(f"  ✗ {field}: MISSING")
            all_present = False
    
    # Check unified subclass
    if motif.get('Subclass') == 'Slipped_DNA':
        print("\n✓ PASS: Uses unified 'Slipped_DNA' subclass")
    else:
        print(f"\n✗ FAIL: Subclass is '{motif.get('Subclass')}', expected 'Slipped_DNA'")
        all_present = False
    
    # Check score range
    score = motif.get('Slippage_Energy_Score', 0)
    if 1.0 <= score <= 3.0:
        print(f"✓ PASS: Score {score:.2f} in valid range [1.0-3.0]")
    else:
        print(f"✗ FAIL: Score {score:.2f} outside valid range [1.0-3.0]")
        all_present = False
    
    return all_present

def main():
    """Run all direct detector tests."""
    print("\n" + "#"*70)
    print("# SLIPPED DNA DETECTOR DIRECT TESTING")
    print("# (No pipeline dependencies)")
    print("#"*70)
    
    results = {
        "Real-World Sequences": test_real_world_sequences(),
        "Edge Cases": test_edge_cases(),
        "Output Format": test_output_format(),
        "Performance": test_performance(),
    }
    
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    
    for test_name, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {test_name}")
    
    total_passed = sum(results.values())
    total_tests = len(results)
    print(f"\nOverall: {total_passed}/{total_tests} test groups passed")
    
    if total_passed == total_tests:
        print("\n✓ All tests PASSED!")
        return 0
    else:
        print(f"\n✗ {total_tests - total_passed} test group(s) FAILED")
        return 1

if __name__ == "__main__":
    sys.exit(main())
