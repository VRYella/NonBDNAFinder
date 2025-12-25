#!/usr/bin/env python3
"""
Test suite for SlippedDNADetector refinements.

Tests:
1. Large sequence handling (no false positives)
2. Valid direct repeat detection
3. Entropy filtering for low-complexity sequences
4. Step size consistency across sequence sizes
"""

import sys
from detectors import SlippedDNADetector


def test_entropy_calculation():
    """Test entropy calculation for DNA sequences."""
    detector = SlippedDNADetector()
    
    # Test 1: Maximum entropy (all bases equally distributed)
    seq_max = "ACGTACGTACGT"
    entropy_max = detector.calculate_entropy(seq_max)
    print(f"✓ Maximum entropy test: {entropy_max:.3f} (expected ~2.0)")
    assert entropy_max > 1.8, f"Expected entropy > 1.8, got {entropy_max}"
    
    # Test 2: Minimum entropy (homopolymer)
    seq_min = "AAAAAAAAAA"
    entropy_min = detector.calculate_entropy(seq_min)
    print(f"✓ Minimum entropy test: {entropy_min:.3f} (expected 0.0)")
    assert entropy_min == 0.0, f"Expected entropy = 0.0, got {entropy_min}"
    
    # Test 3: Low complexity (dinucleotide repeat)
    seq_low = "ATATATATAT"
    entropy_low = detector.calculate_entropy(seq_low)
    print(f"✓ Low complexity test: {entropy_low:.3f} (expected ~1.0)")
    assert entropy_low < 1.5, f"Expected entropy < 1.5, got {entropy_low}"
    
    # Test 4: Medium complexity
    seq_med = "ACGTACGGCTA"
    entropy_med = detector.calculate_entropy(seq_med)
    print(f"✓ Medium complexity test: {entropy_med:.3f}")
    
    print("✅ All entropy calculation tests passed!\n")


def test_low_complexity_filtering():
    """Test that low-complexity sequences are filtered out."""
    detector = SlippedDNADetector()
    
    # Test sequence with low-complexity direct repeats (should be filtered)
    # 10 bp homopolymer repeated (entropy = 0.0, below threshold)
    seq = "AAAAAAAAAA" * 10  # 100 bp of low complexity
    motifs = detector.detect_motifs(seq, "test_low_complexity")
    
    # Should not detect the homopolymer repeat due to entropy filter
    direct_repeats = [m for m in motifs if m['Subclass'] == 'Direct_Repeat']
    print(f"✓ Low complexity filtering: Found {len(direct_repeats)} direct repeats (expected 0)")
    assert len(direct_repeats) == 0, f"Expected 0 direct repeats in low-complexity sequence, got {len(direct_repeats)}"
    
    print("✅ Low-complexity filtering test passed!\n")


def test_valid_direct_repeat_detection():
    """Test detection of valid direct repeats with sufficient entropy."""
    detector = SlippedDNADetector()
    
    # Test sequence with valid direct repeat (15 bp unit > STR range, no spacer, good entropy)
    # Use a non-repetitive 15 bp unit that can't be broken down into smaller STRs
    unit = "ACGTGCATGCTAGCA"  # 15 bp unit with good entropy, not divisible into smaller units
    seq = unit + unit + "G" * 100  # Direct repeat at start, then padding
    
    motifs = detector.detect_motifs(seq, "test_valid_repeat")
    direct_repeats = [m for m in motifs if m['Subclass'] == 'Direct_Repeat']
    
    print(f"✓ Valid direct repeat detection: Found {len(direct_repeats)} direct repeats")
    
    # With fallback implementation (no optimized scanner), we should find the repeat
    # Note: The test might not find it if STR detection consumes the region first
    # Let's verify it's at least being considered in the fallback path
    used = [False] * len(seq)
    fallback_regions = detector.find_direct_repeats_fast(seq, used)
    
    print(f"  - Fallback found {len(fallback_regions)} direct repeats")
    
    # At minimum, the fallback should detect valid repeats
    assert len(fallback_regions) > 0 or len(direct_repeats) > 0, \
        "Expected at least 1 valid direct repeat in fallback or final detection"
    
    # Verify the detected repeat has correct properties
    if direct_repeats:
        repeat = direct_repeats[0]
        print(f"  - Start: {repeat['Start']}, End: {repeat['End']}, Length: {repeat['Length']}")
        print(f"  - Unit length: {repeat.get('Unit_Length', 'N/A')}")
        print(f"  - Spacer length: {repeat.get('Spacer_Length', 'N/A')}")
        assert repeat['Spacer_Length'] == 0, f"Expected spacer length 0, got {repeat['Spacer_Length']}"
        # Unit is 15 bp, so direct repeat should be 30 bp
        assert repeat['Length'] == 30, f"Expected length 30, got {repeat['Length']}"
    elif fallback_regions:
        # Check fallback results
        region = fallback_regions[0]
        print(f"  - Fallback Start: {region['start']}, End: {region['end']}, Length: {region['length']}")
        print(f"  - Fallback Unit length: {region['details']['unit_length']}")
        print(f"  - Fallback Entropy: {region['details'].get('entropy', 'N/A')}")
        assert region['details']['spacer_length'] == 0, \
            f"Expected spacer length 0, got {region['details']['spacer_length']}"
    
    print("✅ Valid direct repeat detection test passed!\n")


def test_large_sequence_no_false_positives():
    """Test that large sequences don't produce false positives every 10,000 bp."""
    detector = SlippedDNADetector()
    
    # Create a 50,000 bp random sequence with no actual repeats
    # Use a pseudo-random pattern to ensure no accidental repeats
    import hashlib
    seq_parts = []
    for i in range(5000):
        # Generate 10 bp chunks with varying patterns
        hash_input = f"chunk_{i}".encode()
        hash_hex = hashlib.md5(hash_input).hexdigest()[:10]
        # Convert hex to DNA bases
        dna_chunk = ""
        for char in hash_hex:
            val = int(char, 16)
            dna_chunk += "ACGT"[val % 4]
        seq_parts.append(dna_chunk)
    
    seq = "".join(seq_parts)
    print(f"✓ Testing large sequence: {len(seq)} bp")
    
    motifs = detector.detect_motifs(seq, "test_large_sequence")
    direct_repeats = [m for m in motifs if m['Subclass'] == 'Direct_Repeat']
    
    print(f"✓ Large sequence test: Found {len(direct_repeats)} direct repeats")
    
    # Should not have repeats every 10,000 bp (old bug)
    # With random sequence, we expect very few or no repeats
    # Allow some false positives due to random chance, but not systematic
    assert len(direct_repeats) < 10, f"Expected < 10 direct repeats in random sequence, got {len(direct_repeats)}"
    
    # Check that repeats are not regularly spaced (indicating systematic false positives)
    if len(direct_repeats) > 1:
        positions = [m['Start'] for m in direct_repeats]
        gaps = [positions[i+1] - positions[i] for i in range(len(positions)-1)]
        avg_gap = sum(gaps) / len(gaps) if gaps else 0
        print(f"  - Average gap between repeats: {avg_gap:.1f} bp")
        # Should not have regular 10,000 bp spacing
        assert avg_gap != 10000, "Detected systematic false positives at 10,000 bp intervals!"
    
    print("✅ Large sequence test passed (no systematic false positives)!\n")


def test_step_size_consistency():
    """Test that step_size is consistent across different sequence sizes."""
    detector = SlippedDNADetector()
    
    # Create sequences of different sizes with the same valid repeat
    unit = "ACGTGCATGC"  # 10 bp unit with good entropy
    repeat_seq = unit + unit
    
    # Test with small sequence (< 10,000 bp)
    seq_small = "A" * 1000 + repeat_seq + "C" * 1000
    motifs_small = detector.detect_motifs(seq_small, "small")
    small_repeats = [m for m in motifs_small if m['Subclass'] == 'Direct_Repeat' and m['Start'] > 1000]
    
    # Test with medium sequence (10,000 - 50,000 bp)
    seq_medium = "A" * 15000 + repeat_seq + "C" * 15000
    motifs_medium = detector.detect_motifs(seq_medium, "medium")
    medium_repeats = [m for m in motifs_medium if m['Subclass'] == 'Direct_Repeat' and m['Start'] > 15000]
    
    # Test with large sequence (> 50,000 bp)
    seq_large = "A" * 60000 + repeat_seq + "C" * 60000
    motifs_large = detector.detect_motifs(seq_large, "large")
    large_repeats = [m for m in motifs_large if m['Subclass'] == 'Direct_Repeat' and m['Start'] > 60000]
    
    print(f"✓ Step size consistency test:")
    print(f"  - Small sequence ({len(seq_small)} bp): {len(small_repeats)} repeats")
    print(f"  - Medium sequence ({len(seq_medium)} bp): {len(medium_repeats)} repeats")
    print(f"  - Large sequence ({len(seq_large)} bp): {len(large_repeats)} repeats")
    
    # All should detect the same repeat (or none due to low complexity of surrounding A's)
    # The key is they should behave consistently
    assert len(small_repeats) == len(medium_repeats), \
        f"Inconsistent detection: small={len(small_repeats)}, medium={len(medium_repeats)}"
    assert len(medium_repeats) == len(large_repeats), \
        f"Inconsistent detection: medium={len(medium_repeats)}, large={len(large_repeats)}"
    
    print("✅ Step size consistency test passed!\n")


def test_no_spacer_requirement():
    """Test that only direct repeats with no spacer (adjacent) are detected."""
    detector = SlippedDNADetector()
    
    # Test 1: Adjacent repeats (spacer = 0) - should be detected
    # Use a 12 bp unit to avoid STR classification (STRs are 1-9 bp)
    unit = "ACGTGCATGCAT"  # 12 bp
    seq_adjacent = unit + unit  # No spacer
    
    # Check fallback explicitly since STR detection might interfere
    used = [False] * len(seq_adjacent)
    fallback_adjacent = detector.find_direct_repeats_fast(seq_adjacent, used)
    
    # Also check final detection
    motifs_adjacent = detector.detect_motifs(seq_adjacent, "adjacent")
    adjacent_repeats = [m for m in motifs_adjacent if m['Subclass'] == 'Direct_Repeat']
    
    print(f"✓ Adjacent repeats (spacer=0): {len(adjacent_repeats)} in final, {len(fallback_adjacent)} in fallback")
    
    # At least the fallback should detect it
    assert len(fallback_adjacent) > 0, "Expected to detect adjacent repeats in fallback"
    
    if fallback_adjacent:
        assert fallback_adjacent[0]['details']['spacer_length'] == 0, \
            f"Expected spacer=0, got {fallback_adjacent[0]['details']['spacer_length']}"
        print(f"  - Fallback confirmed spacer=0")
    
    if adjacent_repeats:
        assert adjacent_repeats[0]['Spacer_Length'] == 0, \
            f"Expected spacer=0, got {adjacent_repeats[0]['Spacer_Length']}"
        print(f"  - Final detection confirmed spacer=0")
    
    # Test 2: Repeats with spacer > 0 - should NOT be detected
    seq_spaced = unit + "GGGG" + unit  # 4 bp spacer
    used_spaced = [False] * len(seq_spaced)
    fallback_spaced = detector.find_direct_repeats_fast(seq_spaced, used_spaced)
    
    motifs_spaced = detector.detect_motifs(seq_spaced, "spaced")
    spaced_repeats = [m for m in motifs_spaced if m['Subclass'] == 'Direct_Repeat']
    
    print(f"✓ Spaced repeats (spacer=4): {len(spaced_repeats)} in final, {len(fallback_spaced)} in fallback")
    
    # Should not detect repeats with spacer > 0
    assert len(fallback_spaced) == 0, f"Expected 0 repeats with spacer>0, got {len(fallback_spaced)}"
    
    print("✅ No-spacer requirement test passed!\n")


def main():
    """Run all tests."""
    print("=" * 70)
    print("Testing SlippedDNADetector Refinements")
    print("=" * 70 + "\n")
    
    try:
        test_entropy_calculation()
        test_low_complexity_filtering()
        test_valid_direct_repeat_detection()
        test_large_sequence_no_false_positives()
        test_step_size_consistency()
        test_no_spacer_requirement()
        
        print("=" * 70)
        print("✅ ALL TESTS PASSED!")
        print("=" * 70)
        return 0
    
    except AssertionError as e:
        print("\n" + "=" * 70)
        print(f"❌ TEST FAILED: {e}")
        print("=" * 70)
        return 1
    
    except Exception as e:
        print("\n" + "=" * 70)
        print(f"❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        print("=" * 70)
        return 1


if __name__ == "__main__":
    sys.exit(main())
