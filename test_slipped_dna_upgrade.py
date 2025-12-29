#!/usr/bin/env python3
"""
Test Script for Mechanism-Driven Slipped DNA Detection Upgrade
================================================================

Tests the new unified, non-redundant slipped DNA detection:
1. Primitive motif computation
2. Repeat purity calculation
3. Stringent entry criteria (≥20bp, ≥90% purity, min copies)
4. Redundancy elimination (one call per locus)
5. Mechanistic slippage scoring (1-3 scale)
6. Publication-ready output format
"""

import sys
sys.path.insert(0, '/home/runner/work/NonBDNAFinder/NonBDNAFinder')

from detectors import SlippedDNADetector

def test_primitive_motif():
    """Test primitive motif computation."""
    print("\n" + "="*70)
    print("TEST 1: Primitive Motif Computation")
    print("="*70)
    
    detector = SlippedDNADetector()
    
    test_cases = [
        ("CAGCAGCAGCAG", "CAG"),  # Should find CAG, not CAGCAG
        ("ATATATATAT", "AT"),      # Should find AT, not ATAT
        ("CACACACACACA", "CA"),    # Should find CA
        ("TATATATATATA", "TA"),    # Should find TA
        ("AGAGAGAGAGAG", "AG"),    # Should find AG
    ]
    
    passed = 0
    for sequence, expected in test_cases:
        result = detector.compute_primitive_motif(sequence)
        status = "PASS" if result == expected else "FAIL"
        print(f"{status}: {sequence[:20]}... → primitive={result} (expected={expected})")
        if status == "PASS":
            passed += 1
    
    print(f"\nResult: {passed}/{len(test_cases)} tests passed")
    return passed == len(test_cases)

def test_repeat_purity():
    """Test repeat purity calculation."""
    print("\n" + "="*70)
    print("TEST 2: Repeat Purity Calculation")
    print("="*70)
    
    detector = SlippedDNADetector()
    
    test_cases = [
        ("CAGCAGCAGCAG", "CAG", 1.0),     # Perfect purity
        ("CAGCAGCAACAG", "CAG", 0.917),   # One mismatch (11/12)
        ("ATATATATAT", "AT", 1.0),        # Perfect purity
        ("ATACATATAC", "AT", 0.8),        # 8/10 matches (corrected)
    ]
    
    passed = 0
    for sequence, unit, expected_purity in test_cases:
        result = detector.compute_repeat_purity(sequence, unit)
        tolerance = 0.01
        status = "PASS" if abs(result - expected_purity) < tolerance else "FAIL"
        print(f"{status}: seq={sequence[:20]}... unit={unit} → purity={result:.3f} (expected={expected_purity:.3f})")
        if status == "PASS":
            passed += 1
    
    print(f"\nResult: {passed}/{len(test_cases)} tests passed")
    return passed == len(test_cases)

def test_redundancy_elimination():
    """Test that redundant STR calls are eliminated (one per locus)."""
    print("\n" + "="*70)
    print("TEST 3: Redundancy Elimination (One Call Per Locus)")
    print("="*70)
    
    detector = SlippedDNADetector()
    
    # Test sequence: (CAG)10 - previously would generate STR_1, STR_2, STR_3, STR_6 etc.
    # Now should generate only ONE call with primitive motif CAG
    test_seq = "CAG" * 10  # 30 bp, perfect (CAG)10 repeat
    
    print(f"Test sequence: {test_seq} (length={len(test_seq)})")
    
    motifs = detector.detect_motifs(test_seq, "test_redundancy")
    
    print(f"\nNumber of motifs detected: {len(motifs)}")
    
    if len(motifs) == 1:
        print("PASS: Exactly one motif detected (redundancy eliminated)")
        motif = motifs[0]
        print(f"  Repeat_Unit: {motif.get('Repeat_Unit', 'N/A')}")
        print(f"  Unit_Size: {motif.get('Unit_Size', 'N/A')}")
        print(f"  Copy_Number: {motif.get('Copy_Number', 'N/A')}")
        print(f"  Purity: {motif.get('Purity', 'N/A')}")
        print(f"  Slippage_Energy_Score: {motif.get('Slippage_Energy_Score', 'N/A')}")
        return True
    else:
        print(f"FAIL: Expected 1 motif, got {len(motifs)}")
        for i, motif in enumerate(motifs):
            print(f"  Motif {i+1}: Unit={motif.get('Repeat_Unit', 'N/A')}, " +
                  f"Size={motif.get('Unit_Size', 'N/A')}, " +
                  f"Copies={motif.get('Copy_Number', 'N/A')}")
        return False

def test_stringent_criteria():
    """Test stringent entry criteria enforcement."""
    print("\n" + "="*70)
    print("TEST 4: Stringent Entry Criteria (≥20bp, ≥90% purity, min copies)")
    print("="*70)
    
    detector = SlippedDNADetector()
    
    # Test 1: Too short (should be rejected)
    short_seq = "CAGCAG"  # 6 bp, only 2 copies
    motifs_short = detector.detect_motifs(short_seq, "test_short")
    print(f"Test 1 (too short, 6bp): {len(motifs_short)} motifs (expected 0)")
    pass1 = len(motifs_short) == 0
    
    # Test 2: Too few copies (should be rejected)
    few_copies = "CAGCAG"  # 6 bp, only 2 copies
    motifs_few = detector.detect_motifs(few_copies, "test_few")
    print(f"Test 2 (too few copies, 2): {len(motifs_few)} motifs (expected 0)")
    pass2 = len(motifs_few) == 0
    
    # Test 3: Low purity (should be rejected if <90%)
    low_purity = "CAGCAGCATCAGCAGCAG"  # Interrupted repeat
    motifs_impure = detector.detect_motifs(low_purity, "test_impure")
    print(f"Test 3 (low purity): {len(motifs_impure)} motifs (expected 0)")
    pass3 = len(motifs_impure) == 0
    
    # Test 4: Valid slipped DNA (should pass)
    valid = "CAG" * 10  # 30 bp, 10 copies, 100% purity
    motifs_valid = detector.detect_motifs(valid, "test_valid")
    print(f"Test 4 (valid, 30bp, 10 copies, 100% purity): {len(motifs_valid)} motifs (expected 1)")
    pass4 = len(motifs_valid) == 1
    
    print(f"\nResult: {sum([pass1, pass2, pass3, pass4])}/4 tests passed")
    return all([pass1, pass2, pass3, pass4])

def test_unified_subclass():
    """Test that all slipped DNA uses unified 'Slipped_DNA' subclass."""
    print("\n" + "="*70)
    print("TEST 5: Unified Subclass (No STR vs Direct_Repeat distinction)")
    print("="*70)
    
    detector = SlippedDNADetector()
    
    # Test with different k values
    test_sequences = [
        ("CAG" * 10, "k=3 (STR range)"),
        ("CAGGAT" * 5, "k=6 (STR range)"),
        ("ACGTACGTACGT" * 3, "k=12 (Direct repeat range)"),
    ]
    
    all_unified = True
    for seq, desc in test_sequences:
        motifs = detector.detect_motifs(seq, "test_unified")
        if motifs:
            subclass = motifs[0].get('Subclass', 'N/A')
            is_unified = subclass == 'Slipped_DNA'
            status = "PASS" if is_unified else "FAIL"
            print(f"{status}: {desc} → Subclass='{subclass}' (expected='Slipped_DNA')")
            if not is_unified:
                all_unified = False
        else:
            print(f"SKIP: {desc} → No motifs detected (may not meet stringent criteria)")
    
    return all_unified

def test_slippage_score_range():
    """Test that slippage energy scores are in [1-3] range."""
    print("\n" + "="*70)
    print("TEST 6: Slippage Energy Score Range (1-3)")
    print("="*70)
    
    detector = SlippedDNADetector()
    
    # Test sequences with varying characteristics
    test_sequences = [
        ("CAG" * 7, "Short, pure repeat (21bp)"),
        ("CAG" * 15, "Long, pure repeat (45bp)"),
        ("CAGGAT" * 8, "Longer unit (k=6, 48bp)"),
    ]
    
    all_in_range = True
    for seq, desc in test_sequences:
        motifs = detector.detect_motifs(seq, "test_score")
        if motifs:
            score = motifs[0].get('Slippage_Energy_Score', 0)
            in_range = 1.0 <= score <= 3.0
            status = "PASS" if in_range else "FAIL"
            print(f"{status}: {desc} → Score={score:.2f} (expected 1.0-3.0)")
            if not in_range:
                all_in_range = False
        else:
            print(f"SKIP: {desc} → No motifs detected")
    
    return all_in_range

def main():
    """Run all tests."""
    print("\n" + "#"*70)
    print("# MECHANISM-DRIVEN SLIPPED DNA DETECTION VALIDATION")
    print("#"*70)
    
    results = {
        "Primitive Motif": test_primitive_motif(),
        "Repeat Purity": test_repeat_purity(),
        "Redundancy Elimination": test_redundancy_elimination(),
        "Stringent Criteria": test_stringent_criteria(),
        "Unified Subclass": test_unified_subclass(),
        "Score Range (1-3)": test_slippage_score_range(),
    }
    
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    for test_name, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {test_name}")
    
    total_passed = sum(results.values())
    total_tests = len(results)
    print(f"\nOverall: {total_passed}/{total_tests} test groups passed")
    
    if total_passed == total_tests:
        print("\n✓ All tests PASSED - Mechanism-driven slipped DNA detection working correctly!")
        return 0
    else:
        print(f"\n✗ {total_tests - total_passed} test group(s) FAILED")
        return 1

if __name__ == "__main__":
    sys.exit(main())
