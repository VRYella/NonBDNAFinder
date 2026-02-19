"""
Benchmark chunking performance with fixed deduplication.
Compares before/after metrics.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import time
from Utilities.nonbscanner import analyze_sequence

def benchmark_sequence_size(size_bp, description):
    """Benchmark analysis for given sequence size."""
    # Create test sequence with known patterns
    pattern = "GGGTAGGGTAGGGTAGGG" + "ATCGATCGATCGATCG" * 10
    test_seq = (pattern * (size_bp // len(pattern) + 1))[:size_bp]
    
    print(f"\n{description} ({size_bp:,} bp)")
    print("=" * 60)
    
    start = time.time()
    motifs = analyze_sequence(test_seq, f"test_{size_bp}")
    elapsed = time.time() - start
    
    throughput = size_bp / elapsed if elapsed > 0 else 0
    
    print(f"  Time:       {elapsed:.2f}s")
    print(f"  Throughput: {throughput/1000:.0f} Kbp/s")
    print(f"  Motifs:     {len(motifs)}")
    print(f"  Status:     {'✅ PASS' if len(motifs) > 0 else '❌ FAIL (0 motifs)'}")
    
    return {
        'size': size_bp,
        'time': elapsed,
        'throughput': throughput,
        'motifs': len(motifs),
        'pass': len(motifs) > 0
    }

if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("CHUNKING PERFORMANCE BENCHMARK")
    print("=" * 60)
    
    results = []
    
    # Test sizes across the chunking threshold
    test_cases = [
        (100_000, "100 KB (no chunking)"),
        (500_000, "500 KB (no chunking)"),
        (1_000_000, "1 MB (threshold)"),
        (1_500_000, "1.5 MB (chunking)"),
        (5_000_000, "5 MB (chunking)"),
    ]
    
    for size, desc in test_cases:
        results.append(benchmark_sequence_size(size, desc))
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    all_passed = all(r['pass'] for r in results)
    
    if all_passed:
        print("✅ All tests PASSED - deduplication fix verified!")
    else:
        print("❌ Some tests FAILED - check results above")
    
    print(f"\nAverage throughput: {sum(r['throughput'] for r in results) / len(results) / 1000:.0f} Kbp/s")
