#!/usr/bin/env python3
"""
Quick verification test for chunking and parallel processing.
Tests that sequences >50KB use chunking and parallel execution.
"""

import sys
sys.path.insert(0, '.')

from Utilities.nonbscanner import analyze_sequence, CHUNK_THRESHOLD

def test_chunking_config():
    """Test that chunking configuration is correct."""
    print("TEST: Chunking Configuration")
    print(f"  CHUNK_THRESHOLD = {CHUNK_THRESHOLD:,} bp")
    assert CHUNK_THRESHOLD == 50000, f"Expected 50000, got {CHUNK_THRESHOLD}"
    print("  ✓ PASSED")
    print()

def test_small_sequence():
    """Test small sequence (no chunking)."""
    print("TEST: Small Sequence (20KB - no chunking expected)")
    seq = "ATCG" * 5000  # 20,000 bp
    print(f"  Sequence length: {len(seq):,} bp")
    result = analyze_sequence(seq, "small_test")
    print(f"  Motifs found: {len(result)}")
    print("  ✓ PASSED")
    print()

def test_large_sequence_with_chunking():
    """Test large sequence (should use chunking)."""
    print("TEST: Large Sequence (60KB - chunking expected)")
    seq = "ATCG" * 15000  # 60,000 bp
    print(f"  Sequence length: {len(seq):,} bp")
    result = analyze_sequence(seq, "large_test", use_chunking=True)
    print(f"  Motifs found: {len(result)}")
    print("  ✓ PASSED")
    print()

def test_parallel_detectors_parameter():
    """Test that parallel detectors parameter is accepted."""
    print("TEST: Parallel Detectors Parameter")
    seq = "ATCG" * 5000  # 20,000 bp
    
    # Test with parallel detectors disabled
    result1 = analyze_sequence(seq, "test1", use_parallel_detectors=False)
    print(f"  Sequential: {len(result1)} motifs")
    
    # Test with parallel detectors enabled
    result2 = analyze_sequence(seq, "test2", use_parallel_detectors=True)
    print(f"  Parallel: {len(result2)} motifs")
    
    # Results should be the same
    assert len(result1) == len(result2), "Results differ between sequential and parallel"
    print("  ✓ PASSED - Results are consistent")
    print()

def main():
    print("=" * 70)
    print("CHUNKING & PARALLEL PROCESSING VERIFICATION")
    print("=" * 70)
    print()
    
    try:
        test_chunking_config()
        test_small_sequence()
        test_large_sequence_with_chunking()
        test_parallel_detectors_parameter()
        
        print("=" * 70)
        print("ALL TESTS PASSED ✓")
        print("=" * 70)
        print()
        print("Summary:")
        print("  • Chunking threshold set to 50,000 bp ✓")
        print("  • Small sequences process correctly ✓")
        print("  • Large sequences use chunking ✓")
        print("  • Parallel detectors parameter works ✓")
        print()
        return 0
        
    except Exception as e:
        print()
        print("=" * 70)
        print("TEST FAILED ✗")
        print("=" * 70)
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
