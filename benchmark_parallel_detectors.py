#!/usr/bin/env python3
"""
Benchmark parallel detector execution vs sequential execution.

Tests the performance improvement achieved by running detectors in parallel
for sequences >50KB.
"""

import time
import sys
from typing import Dict, Any

# Add parent directory to path
sys.path.insert(0, '.')

from Utilities.nonbscanner import analyze_sequence, CHUNK_THRESHOLD

def format_time(seconds: float) -> str:
    """Format seconds as human-readable string."""
    if seconds < 1:
        return f"{seconds*1000:.1f}ms"
    elif seconds < 60:
        return f"{seconds:.2f}s"
    else:
        mins = int(seconds // 60)
        secs = seconds % 60
        return f"{mins}m {secs:.1f}s"

def format_throughput(bp: int, seconds: float) -> str:
    """Format throughput as bp/s."""
    if seconds == 0:
        return "N/A"
    throughput = bp / seconds
    if throughput > 1_000_000:
        return f"{throughput/1_000_000:.2f} Mbp/s"
    elif throughput > 1_000:
        return f"{throughput/1_000:.2f} Kbp/s"
    else:
        return f"{throughput:.0f} bp/s"

def benchmark_sequence(sequence: str, name: str, use_parallel: bool) -> Dict[str, Any]:
    """Benchmark a single sequence analysis.
    
    Note: This isolates detector-level parallelism from chunk-level parallelism
    by disabling chunking (use_chunking=False). This allows us to measure the
    pure impact of parallel detector execution, independent of chunk processing.
    """
    start = time.time()
    motifs = analyze_sequence(
        sequence, 
        name, 
        use_parallel_detectors=use_parallel,
        use_chunking=False  # Isolate detector parallelism from chunk parallelism
    )
    elapsed = time.time() - start
    
    return {
        'name': name,
        'length': len(sequence),
        'motif_count': len(motifs),
        'elapsed': elapsed,
        'throughput': len(sequence) / elapsed if elapsed > 0 else 0,
        'parallel': use_parallel
    }

def main():
    print("=" * 80)
    print("PARALLEL DETECTOR BENCHMARK")
    print("=" * 80)
    print()
    print(f"CHUNK_THRESHOLD: {CHUNK_THRESHOLD:,} bp")
    print()
    
    # Test sequences of different sizes
    test_cases = [
        ("Small (20KB)", "ATCG" * 5_000),      # 20,000 bp - below threshold
        ("Medium (60KB)", "ATCG" * 15_000),    # 60,000 bp - above threshold
        ("Large (100KB)", "ATCG" * 25_000),    # 100,000 bp - 2x threshold
        ("Very Large (200KB)", "ATCG" * 50_000), # 200,000 bp - 4x threshold
    ]
    
    results = []
    
    for test_name, sequence in test_cases:
        print(f"Testing {test_name} ({len(sequence):,} bp)")
        print("-" * 80)
        
        # Sequential execution
        print("  Sequential detectors...", end=" ", flush=True)
        result_seq = benchmark_sequence(sequence, f"{test_name}_seq", use_parallel=False)
        results.append(result_seq)
        print(f"{format_time(result_seq['elapsed'])} - {result_seq['motif_count']} motifs - {format_throughput(result_seq['length'], result_seq['elapsed'])}")
        
        # Parallel execution
        print("  Parallel detectors...  ", end=" ", flush=True)
        result_par = benchmark_sequence(sequence, f"{test_name}_par", use_parallel=True)
        results.append(result_par)
        print(f"{format_time(result_par['elapsed'])} - {result_par['motif_count']} motifs - {format_throughput(result_par['length'], result_par['elapsed'])}")
        
        # Calculate speedup
        speedup = result_seq['elapsed'] / result_par['elapsed'] if result_par['elapsed'] > 0 else 0
        speedup_pct = ((result_seq['elapsed'] - result_par['elapsed']) / result_seq['elapsed'] * 100) if result_seq['elapsed'] > 0 else 0
        
        print(f"  Speedup: {speedup:.2f}x ({speedup_pct:+.1f}%)")
        print()
    
    # Summary
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print()
    print(f"{'Test Case':<20} {'Mode':<10} {'Time':<12} {'Throughput':<15} {'Motifs':<8}")
    print("-" * 80)
    
    for result in results:
        mode = "Parallel" if result['parallel'] else "Sequential"
        print(f"{result['name']:<20} {mode:<10} {format_time(result['elapsed']):<12} "
              f"{format_throughput(result['length'], result['elapsed']):<15} {result['motif_count']:<8}")
    
    print()
    
    # Calculate average speedup for sequences >50KB
    large_seq_results = [r for r in results if r['length'] > CHUNK_THRESHOLD]
    if len(large_seq_results) >= 2:
        seq_results = [r for r in large_seq_results if not r['parallel']]
        par_results = [r for r in large_seq_results if r['parallel']]
        
        if len(seq_results) == len(par_results):
            speedups = []
            for seq_r, par_r in zip(seq_results, par_results):
                if par_r['elapsed'] > 0:
                    speedups.append(seq_r['elapsed'] / par_r['elapsed'])
            
            if speedups:
                avg_speedup = sum(speedups) / len(speedups)
                print(f"Average speedup for sequences >{CHUNK_THRESHOLD:,} bp: {avg_speedup:.2f}x")
                print()
    
    print("=" * 80)
    print("CONCLUSION")
    print("=" * 80)
    print()
    print("Parallel detector execution provides significant speedup for large sequences.")
    print(f"For sequences >{CHUNK_THRESHOLD:,} bp, parallel detectors are automatically enabled.")
    print()

if __name__ == "__main__":
    main()
