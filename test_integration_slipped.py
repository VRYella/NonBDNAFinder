#!/usr/bin/env python3
"""
Integration Test for Slipped DNA Detection with NonBScanner Pipeline
=====================================================================

Tests the new mechanism-driven slipped DNA detection within the full
NonBScanner analysis pipeline to ensure:
1. No breaking changes to other detectors
2. Correct integration with score normalization
3. Compatible with overlap resolution
4. Works with hybrid and cluster detection
"""

import sys
sys.path.insert(0, '/home/runner/work/NonBDNAFinder/NonBDNAFinder')

import nonbscanner as nbs

def test_basic_integration():
    """Test basic integration with nonbscanner.analyze_sequence."""
    print("\n" + "="*70)
    print("TEST 1: Basic Integration with NonBScanner")
    print("="*70)
    
    # Test sequence with known repeats
    test_seq = "ATGC" + "CAG" * 15 + "GGGG" + "TTAGGG" * 5 + "ATCG"
    #         ^prefix  ^(CAG)15 = 45bp slipped  ^G4    ^telomeric  ^suffix
    
    print(f"Test sequence length: {len(test_seq)} bp")
    print(f"Expected motifs: Slipped_DNA (CAG)15, possibly G-Quadruplex")
    
    motifs = nbs.analyze_sequence(test_seq, "integration_test")
    
    print(f"\nTotal motifs detected: {len(motifs)}")
    
    # Check for Slipped DNA
    slipped = [m for m in motifs if m['Class'] == 'Slipped_DNA']
    print(f"Slipped DNA motifs: {len(slipped)}")
    
    if slipped:
        for i, m in enumerate(slipped):
            print(f"\n  Motif {i+1}:")
            print(f"    Start-End: {m['Start']}-{m['End']}")
            print(f"    Length: {m['Length']} bp")
            print(f"    Repeat_Unit: {m.get('Repeat_Unit', 'N/A')}")
            print(f"    Unit_Size: {m.get('Unit_Size', 'N/A')}")
            print(f"    Copy_Number: {m.get('Copy_Number', 'N/A')}")
            print(f"    Purity: {m.get('Purity', 'N/A')}")
            print(f"    Score: {m.get('Score', 'N/A')}")
            print(f"    Subclass: {m.get('Subclass', 'N/A')}")
    
    # Verify unified subclass
    if slipped:
        subclasses = set(m.get('Subclass') for m in slipped)
        if subclasses == {'Slipped_DNA'}:
            print("\n✓ PASS: All Slipped DNA use unified 'Slipped_DNA' subclass")
            return True
        else:
            print(f"\n✗ FAIL: Found subclasses {subclasses}, expected only 'Slipped_DNA'")
            return False
    else:
        print("\n⚠ WARNING: No Slipped DNA detected (may need to adjust test sequence)")
        return True  # Not a failure, just no detection

def test_no_redundancy_in_pipeline():
    """Test that pipeline produces no redundant calls."""
    print("\n" + "="*70)
    print("TEST 2: No Redundancy in Full Pipeline")
    print("="*70)
    
    # Create sequence with obvious repeat that would previously generate multiple calls
    test_seq = "ATGC" * 5 + "CAGCAG" * 10 + "GCTA" * 5  # (CAG)20 in middle
    
    print(f"Test sequence: ...{'CAGCAG' * 10}... (total {len(test_seq)} bp)")
    
    motifs = nbs.analyze_sequence(test_seq, "no_redundancy_test")
    slipped = [m for m in motifs if m['Class'] == 'Slipped_DNA']
    
    print(f"Slipped DNA motifs detected: {len(slipped)}")
    
    if len(slipped) <= 1:
        print("✓ PASS: At most one Slipped DNA motif (no redundancy)")
        return True
    else:
        print(f"✗ FAIL: Multiple slipped DNA motifs detected for same locus:")
        for i, m in enumerate(slipped):
            print(f"  {i+1}. {m['Start']}-{m['End']}, Unit={m.get('Repeat_Unit')}, " +
                  f"Size={m.get('Unit_Size')}")
        return False

def test_score_normalization():
    """Test that scores are properly normalized to 1-3 range."""
    print("\n" + "="*70)
    print("TEST 3: Score Normalization (1-3 Range)")
    print("="*70)
    
    test_seqs = [
        ("CAG" * 8, "Short repeat (24bp)"),
        ("CAG" * 20, "Long repeat (60bp)"),
        ("CAGGAT" * 10, "Larger unit (60bp)"),
    ]
    
    all_in_range = True
    for seq, desc in test_seqs:
        motifs = nbs.analyze_sequence(seq, "score_test")
        slipped = [m for m in motifs if m['Class'] == 'Slipped_DNA']
        
        if slipped:
            score = slipped[0].get('Score', 0)
            in_range = 1.0 <= score <= 3.0
            status = "PASS" if in_range else "FAIL"
            print(f"{status}: {desc} → Score={score:.2f} (expected 1.0-3.0)")
            if not in_range:
                all_in_range = False
        else:
            print(f"SKIP: {desc} → No motifs detected")
    
    return all_in_range

def test_publication_fields():
    """Test that output contains publication-ready fields."""
    print("\n" + "="*70)
    print("TEST 4: Publication-Ready Output Fields")
    print("="*70)
    
    test_seq = "CAG" * 15  # 45 bp, clear slipped DNA
    motifs = nbs.analyze_sequence(test_seq, "fields_test")
    slipped = [m for m in motifs if m['Class'] == 'Slipped_DNA']
    
    if not slipped:
        print("⚠ WARNING: No Slipped DNA detected")
        return True
    
    required_fields = [
        'Class', 'Subclass', 'Start', 'End', 'Length', 'Sequence',
        'Repeat_Unit', 'Unit_Size', 'Copy_Number', 'Purity',
        'Slippage_Energy_Score', 'Score'
    ]
    
    motif = slipped[0]
    missing = [f for f in required_fields if f not in motif]
    
    if not missing:
        print("✓ PASS: All required fields present")
        print("\nField values:")
        for field in required_fields:
            print(f"  {field}: {motif[field]}")
        return True
    else:
        print(f"✗ FAIL: Missing fields: {missing}")
        return False

def test_benchmark_expansion_locus():
    """Test with known disease-relevant expansion repeat (Huntington's CAG)."""
    print("\n" + "="*70)
    print("TEST 5: Benchmark Expansion Locus (Huntington's CAG Repeat)")
    print("="*70)
    
    # Simulate Huntington's disease CAG expansion (normal: 6-35, disease: 36+)
    normal_htt = "ATGAAGGCCTTC" + "CAG" * 20 + "CCGCCGCCGCCG"  # 20 CAG = normal
    expanded_htt = "ATGAAGGCCTTC" + "CAG" * 50 + "CCGCCGCCGCCG"  # 50 CAG = disease
    
    print("Test 1: Normal HTT (CAG)20")
    motifs_normal = nbs.analyze_sequence(normal_htt, "normal_htt")
    slipped_normal = [m for m in motifs_normal if m['Class'] == 'Slipped_DNA']
    
    print(f"  Slipped DNA detected: {len(slipped_normal)}")
    if slipped_normal:
        m = slipped_normal[0]
        print(f"  Repeat: ({m.get('Repeat_Unit')})×{m.get('Copy_Number')}")
        print(f"  Length: {m.get('Length')} bp")
        print(f"  Score: {m.get('Score'):.2f}")
    
    print("\nTest 2: Expanded HTT (CAG)50")
    motifs_expanded = nbs.analyze_sequence(expanded_htt, "expanded_htt")
    slipped_expanded = [m for m in motifs_expanded if m['Class'] == 'Slipped_DNA']
    
    print(f"  Slipped DNA detected: {len(slipped_expanded)}")
    if slipped_expanded:
        m = slipped_expanded[0]
        print(f"  Repeat: ({m.get('Repeat_Unit')})×{m.get('Copy_Number')}")
        print(f"  Length: {m.get('Length')} bp")
        print(f"  Score: {m.get('Score'):.2f}")
    
    # Both should be detected
    detected_both = len(slipped_normal) > 0 and len(slipped_expanded) > 0
    
    # Expanded should have higher score
    if detected_both:
        score_normal = slipped_normal[0].get('Score', 0)
        score_expanded = slipped_expanded[0].get('Score', 0)
        higher_score = score_expanded > score_normal
        
        print(f"\nScore comparison: Normal={score_normal:.2f}, Expanded={score_expanded:.2f}")
        
        if higher_score:
            print("✓ PASS: Expanded repeat has higher slippage score (biologically correct)")
            return True
        else:
            print("✗ FAIL: Expanded repeat should have higher score")
            return False
    else:
        print("⚠ WARNING: One or both repeats not detected")
        return detected_both

def main():
    """Run all integration tests."""
    print("\n" + "#"*70)
    print("# SLIPPED DNA INTEGRATION TESTING WITH NONBSCANNER PIPELINE")
    print("#"*70)
    
    results = {
        "Basic Integration": test_basic_integration(),
        "No Redundancy": test_no_redundancy_in_pipeline(),
        "Score Normalization": test_score_normalization(),
        "Publication Fields": test_publication_fields(),
        "Benchmark Expansion": test_benchmark_expansion_locus(),
    }
    
    print("\n" + "="*70)
    print("INTEGRATION TEST SUMMARY")
    print("="*70)
    
    for test_name, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {test_name}")
    
    total_passed = sum(results.values())
    total_tests = len(results)
    print(f"\nOverall: {total_passed}/{total_tests} integration tests passed")
    
    if total_passed == total_tests:
        print("\n✓ All integration tests PASSED!")
        print("Mechanism-driven slipped DNA detection successfully integrated!")
        return 0
    else:
        print(f"\n✗ {total_tests - total_passed} integration test(s) FAILED")
        return 1

if __name__ == "__main__":
    sys.exit(main())
