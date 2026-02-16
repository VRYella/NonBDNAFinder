#!/usr/bin/env python3
"""
Performance benchmark for triple adaptive chunking strategy.

Tests performance targets:
- 10MB: <30 seconds
- 100MB: <2 minutes
- 500MB: <5 minutes
- 1GB: <8 minutes
"""

import time
import tempfile
import shutil
from pathlib import Path
import sys
import gc

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from Utilities.disk_storage import UniversalSequenceStorage
from Utilities.triple_chunk_analyzer import TripleAdaptiveChunkAnalyzer
from Utilities.chunk_analyzer import ChunkAnalyzer


def format_time(seconds):
    """Format seconds into human-readable string."""
    if seconds < 60:
        return f"{seconds:.2f}s"
    else:
        minutes = int(seconds // 60)
        secs = seconds % 60
        return f"{minutes}m {secs:.2f}s"


def format_size(bp):
    """Format base pairs into human-readable string."""
    if bp < 1_000_000:
        return f"{bp / 1000:.1f}KB"
    elif bp < 1_000_000_000:
        return f"{bp / 1_000_000:.1f}MB"
    else:
        return f"{bp / 1_000_000_000:.2f}GB"


def benchmark_sequence(size_bp, name, test_dir):
    """
    Benchmark analysis of a sequence of given size.
    
    Args:
        size_bp: Size of sequence in base pairs
        name: Name for the sequence
        test_dir: Temporary directory for storage
        
    Returns:
        Tuple of (elapsed_time, motif_count, throughput_bp_per_sec)
    """
    print(f"\n{'='*70}")
    print(f"Benchmarking {name} ({format_size(size_bp)})")
    print(f"{'='*70}")
    
    # Create sequence
    print("Creating sequence...", end=" ", flush=True)
    create_start = time.time()
    # Use simple repeating pattern for speed
    sequence = "ATCGATCG" * (size_bp // 8)
    create_time = time.time() - create_start
    print(f"done ({create_time:.2f}s)")
    
    # Save to storage
    storage = UniversalSequenceStorage(base_dir=test_dir)
    print("Saving to storage...", end=" ", flush=True)
    save_start = time.time()
    seq_id = storage.save_sequence(sequence, name)
    save_time = time.time() - save_start
    print(f"done ({save_time:.2f}s)")
    
    # Free memory
    del sequence
    gc.collect()
    
    # Benchmark adaptive chunking
    print(f"\nAnalyzing with adaptive chunking...")
    analyzer = TripleAdaptiveChunkAnalyzer(
        storage,
        use_adaptive=True
    )
    
    analysis_start = time.time()
    results_storage = analyzer.analyze(seq_id)
    analysis_time = time.time() - analysis_start
    
    # Get results
    stats = results_storage.get_summary_stats()
    motif_count = stats['total_count']
    
    # Calculate throughput
    throughput = size_bp / analysis_time if analysis_time > 0 else 0
    
    # Cleanup
    results_storage.cleanup()
    
    print(f"\n{'='*70}")
    print(f"RESULTS: {name}")
    print(f"{'='*70}")
    print(f"  Sequence size:     {format_size(size_bp)} ({size_bp:,} bp)")
    print(f"  Analysis time:     {format_time(analysis_time)}")
    print(f"  Motifs detected:   {motif_count:,}")
    print(f"  Throughput:        {throughput:,.0f} bp/s ({format_size(throughput)}/s)")
    print(f"  Total time:        {format_time(create_time + save_time + analysis_time)}")
    
    return analysis_time, motif_count, throughput


def compare_strategies(size_bp, name, test_dir):
    """
    Compare adaptive vs non-adaptive strategies.
    
    Args:
        size_bp: Size of sequence in base pairs
        name: Name for the sequence
        test_dir: Temporary directory for storage
        
    Returns:
        Tuple of (adaptive_time, non_adaptive_time, speedup)
    """
    print(f"\n{'='*70}")
    print(f"Strategy Comparison: {name} ({format_size(size_bp)})")
    print(f"{'='*70}")
    
    # Create and save sequence
    sequence = "ATCGATCG" * (size_bp // 8)
    storage = UniversalSequenceStorage(base_dir=test_dir)
    seq_id = storage.save_sequence(sequence, name)
    del sequence
    gc.collect()
    
    # Test adaptive strategy
    print("\n1. Testing adaptive strategy...")
    analyzer_adaptive = TripleAdaptiveChunkAnalyzer(
        storage,
        use_adaptive=True
    )
    
    start = time.time()
    results_adaptive = analyzer_adaptive.analyze(seq_id)
    adaptive_time = time.time() - start
    stats_adaptive = results_adaptive.get_summary_stats()
    results_adaptive.cleanup()
    
    print(f"   Time: {format_time(adaptive_time)}")
    print(f"   Motifs: {stats_adaptive['total_count']:,}")
    
    # Test non-adaptive strategy
    print("\n2. Testing non-adaptive strategy...")
    analyzer_non_adaptive = ChunkAnalyzer(
        storage,
        chunk_size=5_000_000,
        overlap=10_000,
        use_adaptive=False
    )
    
    start = time.time()
    results_non_adaptive = analyzer_non_adaptive.analyze(seq_id)
    non_adaptive_time = time.time() - start
    stats_non_adaptive = results_non_adaptive.get_summary_stats()
    results_non_adaptive.cleanup()
    
    print(f"   Time: {format_time(non_adaptive_time)}")
    print(f"   Motifs: {stats_non_adaptive['total_count']:,}")
    
    # Calculate speedup
    speedup = non_adaptive_time / adaptive_time if adaptive_time > 0 else 0
    
    print(f"\n{'='*70}")
    print(f"COMPARISON RESULTS")
    print(f"{'='*70}")
    print(f"  Adaptive time:       {format_time(adaptive_time)}")
    print(f"  Non-adaptive time:   {format_time(non_adaptive_time)}")
    print(f"  Speedup:             {speedup:.2f}x")
    print(f"  Motif difference:    {abs(stats_adaptive['total_count'] - stats_non_adaptive['total_count'])}")
    
    return adaptive_time, non_adaptive_time, speedup


def main():
    """Run performance benchmarks."""
    print("=" * 70)
    print("TRIPLE ADAPTIVE CHUNKING PERFORMANCE BENCHMARK")
    print("=" * 70)
    print("\nPerformance Targets:")
    print("  10MB:   < 30 seconds")
    print("  100MB:  < 2 minutes (120s)")
    print("  500MB:  < 5 minutes (300s)")
    print("  1GB:    < 8 minutes (480s)")
    
    # Create temporary directory
    test_dir = tempfile.mkdtemp(prefix="benchmark_triple_")
    
    try:
        # Benchmark suite
        benchmarks = [
            (10_000_000, "10MB Genome", 30),
            # Full benchmark suite - Uncomment to test larger genomes
            # NOTE: These tests take longer (100MB: ~2min, 500MB: ~5min, 1GB: ~8min)
            # Recommended for CI/CD or pre-release validation
            # (100_000_000, "100MB Genome", 120),
            # (500_000_000, "500MB Genome", 300),
            # (1_000_000_000, "1GB Genome", 480),
        ]
        
        results = []
        
        for size_bp, name, target_time in benchmarks:
            elapsed, motifs, throughput = benchmark_sequence(size_bp, name, test_dir)
            
            # Check if target met
            target_met = elapsed <= target_time
            status = "✓ PASS" if target_met else "✗ FAIL"
            
            print(f"\n  Target: {format_time(target_time)} | Actual: {format_time(elapsed)} | {status}")
            
            results.append({
                'name': name,
                'size': size_bp,
                'time': elapsed,
                'target': target_time,
                'passed': target_met,
                'motifs': motifs,
                'throughput': throughput
            })
            
            # Clean up between benchmarks
            gc.collect()
        
        # Summary
        print(f"\n{'='*70}")
        print("BENCHMARK SUMMARY")
        print(f"{'='*70}")
        
        for result in results:
            status = "✓" if result['passed'] else "✗"
            print(f"{status} {result['name']:<15} {format_time(result['time']):<12} "
                  f"(target: {format_time(result['target'])})")
        
        passed = sum(1 for r in results if r['passed'])
        total = len(results)
        print(f"\nPassed: {passed}/{total} benchmarks")
        
        # Optional: Strategy comparison on smaller sequence
        print(f"\n{'='*70}")
        print("STRATEGY COMPARISON (Optional)")
        print(f"{'='*70}")
        print("\nComparing adaptive vs non-adaptive on 10MB sequence...")
        
        try:
            compare_strategies(10_000_000, "10MB_comparison", test_dir)
        except Exception as e:
            print(f"Comparison failed: {e}")
        
    finally:
        # Cleanup
        print(f"\nCleaning up temporary directory...")
        if Path(test_dir).exists():
            shutil.rmtree(test_dir)
        print("Done!")


if __name__ == '__main__':
    main()
