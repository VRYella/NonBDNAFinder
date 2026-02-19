"""
End-to-end test to verify motif detection works for sequences > 1MB.
This test specifically addresses the reported issue where motif detection
was not working for sequences longer than 1MB.
"""
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from Utilities.nonbscanner import analyze_sequence


def test_1mb_sequence_detection():
    """Test that motif detection works for a sequence > 1MB."""
    print("Creating a 1.1MB test sequence...")
    
    # Create a 1.1MB sequence with some recognizable patterns
    # Include G-quadruplex motifs: GGG sequences
    gquad_pattern = "GGGTAGGGTAGGGTAGGG"  # G-quadruplex pattern
    filler = "ATCGATCGATCGATCG"
    
    # Calculate how many times to repeat to get > 1MB
    target_size = 1_100_000  # 1.1 MB
    pattern = gquad_pattern + filler * 10  # Mix in some patterns
    repetitions = target_size // len(pattern) + 1
    
    test_sequence = pattern * repetitions
    test_sequence = test_sequence[:target_size]  # Trim to exact size
    
    seq_length = len(test_sequence)
    print(f"Test sequence length: {seq_length:,} bp ({seq_length / 1_000_000:.2f} MB)")
    
    if seq_length <= 1_000_000:
        print("ERROR: Test sequence is not > 1MB!")
        return False
    
    print("Running motif detection...")
    try:
        motifs = analyze_sequence(test_sequence, "test_1mb_sequence")
        print(f"✓ Motif detection completed successfully!")
        print(f"  Found {len(motifs)} motifs")
        
        # Check for critical bug
        if len(motifs) == 0:
            print("❌ CRITICAL BUG: No motifs detected!")
            print("   This indicates the deduplication bug is still present")
            return False
        
        # Show some statistics
        if motifs:
            motif_classes = {}
            for motif in motifs:
                cls = motif.get('Class', 'Unknown')
                motif_classes[cls] = motif_classes.get(cls, 0) + 1
            
            print(f"  Motif classes detected:")
            for cls, count in sorted(motif_classes.items()):
                print(f"    - {cls}: {count}")
        
        return True
        
    except Exception as e:
        print(f"✗ ERROR: Motif detection failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_500kb_sequence_detection():
    """Test that motif detection works for a 500KB sequence (< 1MB, should not chunk)."""
    print("\nCreating a 500KB test sequence...")
    
    # Create a 500KB sequence
    target_size = 500_000
    pattern = "ATCGATCGATCGATCG"
    repetitions = target_size // len(pattern) + 1
    test_sequence = (pattern * repetitions)[:target_size]
    
    seq_length = len(test_sequence)
    print(f"Test sequence length: {seq_length:,} bp ({seq_length / 1_000:.0f} KB)")
    
    print("Running motif detection...")
    try:
        motifs = analyze_sequence(test_sequence, "test_500kb_sequence")
        print(f"✓ Motif detection completed successfully!")
        print(f"  Found {len(motifs)} motifs")
        return True
    except Exception as e:
        print(f"✗ ERROR: Motif detection failed: {e}")
        return False


if __name__ == "__main__":
    print("="*70)
    print("Testing motif detection for sequences > 1MB")
    print("="*70)
    
    # Test 500KB (should not chunk)
    result1 = test_500kb_sequence_detection()
    
    # Test 1.1MB (should chunk)
    result2 = test_1mb_sequence_detection()
    
    print("\n" + "="*70)
    if result1 and result2:
        print("✓ All tests PASSED!")
        print("  Motif detection now works correctly for sequences of all sizes.")
        sys.exit(0)
    else:
        print("✗ Some tests FAILED!")
        sys.exit(1)
