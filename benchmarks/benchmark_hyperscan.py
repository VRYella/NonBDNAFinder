#!/usr/bin/env python3
"""
Hyperscan Backend Benchmark
============================

This benchmark compares the performance of Hyperscan-accelerated pattern matching
against pure Python regex matching for Non-B DNA motif detection.

Requirements:
    - libhs (Hyperscan library) must be installed on the system
    - pip install hyperscan

Usage:
    python benchmarks/benchmark_hyperscan.py [--size SIZE] [--iterations ITERATIONS]

Example:
    python benchmarks/benchmark_hyperscan.py --size 100000 --iterations 5
"""

import argparse
import random
import sys
import time
from typing import Dict, List, Tuple

# Add parent directory to path for imports
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def generate_random_sequence(length: int) -> str:
    """Generate a random DNA sequence of specified length."""
    return ''.join(random.choice('ACGT') for _ in range(length))


def generate_sequence_with_motifs(length: int) -> str:
    """
    Generate a DNA sequence with embedded motifs for testing.
    
    Includes:
    - G-quadruplex patterns (GGGG...GGGG)
    - A-tracts (AAAAA)
    - Alternating patterns for Z-DNA
    """
    seq = list(generate_random_sequence(length))
    
    # Insert some G-quadruplex-like patterns
    g4_pattern = "GGGTTAGGGTTAGGGTTAGGG"
    for i in range(0, length - len(g4_pattern) - 100, length // 10):
        pos = i + random.randint(0, 50)
        if pos + len(g4_pattern) < length:
            for j, base in enumerate(g4_pattern):
                seq[pos + j] = base
    
    # Insert some A-tracts
    a_tract = "AAAAAAAAAA"
    for i in range(length // 20, length - len(a_tract), length // 8):
        pos = i + random.randint(0, 30)
        if pos + len(a_tract) < length:
            for j, base in enumerate(a_tract):
                seq[pos + j] = base
    
    # Insert alternating CG for Z-DNA potential
    cg_pattern = "CGCGCGCGCGCG"
    for i in range(length // 15, length - len(cg_pattern), length // 6):
        pos = i + random.randint(0, 40)
        if pos + len(cg_pattern) < length:
            for j, base in enumerate(cg_pattern):
                seq[pos + j] = base
    
    return ''.join(seq)


def get_test_patterns() -> List[Tuple[str, str, Dict]]:
    """
    Get a set of test patterns for benchmarking.
    
    Note: Uses possessive quantifiers ({n,}+) to prevent catastrophic
    backtracking on long sequences. Requires Python 3.11+.
    
    Returns:
        List of (regex, pattern_id, metadata) tuples
    """
    return [
        # G-quadruplex patterns (possessive quantifiers to prevent backtracking)
        (r'G{3,}+[ACGT]{1,7}+G{3,}+[ACGT]{1,7}+G{3,}+[ACGT]{1,7}+G{3,}+', 
         'G4_1', {'class': 'G-Quadruplex', 'subclass': 'Canonical', 'score': 0.9}),
        (r'G{4,}', 'G_RUN', {'class': 'G-Quadruplex', 'subclass': 'G-run', 'score': 0.5}),
        
        # i-Motif patterns (possessive quantifiers to prevent backtracking)
        (r'C{3,}+[ACGT]{1,7}+C{3,}+[ACGT]{1,7}+C{3,}+[ACGT]{1,7}+C{3,}+', 
         'IM_1', {'class': 'i-Motif', 'subclass': 'Canonical', 'score': 0.9}),
        (r'C{4,}', 'C_RUN', {'class': 'i-Motif', 'subclass': 'C-run', 'score': 0.5}),
        
        # A-tract patterns (curved DNA)
        (r'A{4,}', 'A_TRACT', {'class': 'Curved_DNA', 'subclass': 'A-tract', 'score': 0.7}),
        (r'T{4,}', 'T_TRACT', {'class': 'Curved_DNA', 'subclass': 'T-tract', 'score': 0.7}),
        
        # Z-DNA patterns
        (r'(?:CG){5,}', 'ZDNA_CG', {'class': 'Z-DNA', 'subclass': 'CG-repeat', 'score': 0.8}),
        (r'(?:GC){5,}', 'ZDNA_GC', {'class': 'Z-DNA', 'subclass': 'GC-repeat', 'score': 0.8}),
        (r'(?:CA){5,}', 'ZDNA_CA', {'class': 'Z-DNA', 'subclass': 'CA-repeat', 'score': 0.6}),
        
        # Simple repeats
        (r'(?:AT){6,}', 'AT_REPEAT', {'class': 'STR', 'subclass': 'Dinucleotide', 'score': 0.6}),
        (r'(?:TA){6,}', 'TA_REPEAT', {'class': 'STR', 'subclass': 'Dinucleotide', 'score': 0.6}),
    ]


def benchmark_regex(sequence: str, patterns: List[Tuple[str, str, Dict]], 
                    iterations: int = 3) -> Tuple[float, int]:
    """
    Benchmark pure Python regex matching.
    
    Note: Pattern compilation is done outside the timing loop to measure
    only the matching performance, consistent with the Hyperscan benchmark.
    
    Returns:
        Tuple of (average_time_seconds, total_matches)
    """
    import re
    
    # Pre-compile patterns (not included in timing)
    compiled_patterns = [
        (re.compile(regex, re.IGNORECASE), pattern_id, metadata)
        for regex, pattern_id, metadata in patterns
    ]
    
    total_matches = 0
    times = []
    
    for _ in range(iterations):
        start = time.perf_counter()
        matches = 0
        
        for compiled, pattern_id, metadata in compiled_patterns:
            for match in compiled.finditer(sequence):
                matches += 1
        
        elapsed = time.perf_counter() - start
        times.append(elapsed)
        total_matches = matches  # Same each iteration
    
    avg_time = sum(times) / len(times)
    return avg_time, total_matches


def benchmark_hyperscan(sequence: str, patterns: List[Tuple[str, str, Dict]], 
                        iterations: int = 3) -> Tuple[float, int]:
    """
    Benchmark Hyperscan pattern matching.
    
    Returns:
        Tuple of (average_time_seconds, total_matches)
    """
    from scanner_backends.hyperscan_backend import HyperscanScanner, is_hyperscan_available
    
    if not is_hyperscan_available():
        raise RuntimeError("Hyperscan is not available")
    
    scanner = HyperscanScanner()
    
    # Compilation time (not included in benchmark)
    scanner.compile_patterns(patterns)
    
    times = []
    total_matches = 0
    
    for _ in range(iterations):
        start = time.perf_counter()
        matches = scanner.scan(sequence)
        elapsed = time.perf_counter() - start
        
        times.append(elapsed)
        total_matches = len(matches)
    
    avg_time = sum(times) / len(times)
    return avg_time, total_matches


def benchmark_fallback_regex(sequence: str, patterns: List[Tuple[str, str, Dict]], 
                             iterations: int = 3) -> Tuple[float, int]:
    """
    Benchmark the fallback regex implementation from hyperscan_backend.
    
    Note: This benchmarks the fallback_regex_scan function which includes
    pattern compilation time internally. This represents the actual
    performance when using the fallback path.
    
    Returns:
        Tuple of (average_time_seconds, total_matches)
    """
    from scanner_backends.hyperscan_backend import fallback_regex_scan
    
    times = []
    total_matches = 0
    
    for _ in range(iterations):
        start = time.perf_counter()
        matches = fallback_regex_scan(sequence, patterns)
        elapsed = time.perf_counter() - start
        
        times.append(elapsed)
        total_matches = len(matches)
    
    avg_time = sum(times) / len(times)
    return avg_time, total_matches


def format_time(seconds: float) -> str:
    """Format time in human-readable format."""
    if seconds < 0.001:
        return f"{seconds * 1000000:.2f} µs"
    elif seconds < 1:
        return f"{seconds * 1000:.2f} ms"
    else:
        return f"{seconds:.3f} s"


def format_throughput(sequence_length: int, seconds: float) -> str:
    """Format throughput in bases per second."""
    if seconds == 0:
        return "N/A"
    bps = sequence_length / seconds
    if bps >= 1e9:
        return f"{bps / 1e9:.2f} Gbp/s"
    elif bps >= 1e6:
        return f"{bps / 1e6:.2f} Mbp/s"
    elif bps >= 1e3:
        return f"{bps / 1e3:.2f} Kbp/s"
    else:
        return f"{bps:.2f} bp/s"


def run_benchmark(size: int = 100000, iterations: int = 5):
    """
    Run the complete benchmark suite.
    
    Args:
        size: Sequence length in base pairs
        iterations: Number of iterations per benchmark
    """
    from scanner_backends import is_hyperscan_available
    
    print("=" * 70)
    print("Hyperscan Backend Benchmark")
    print("=" * 70)
    print(f"\nSequence size: {size:,} bp")
    print(f"Iterations: {iterations}")
    print()
    
    # Generate test sequence
    print("Generating test sequence with embedded motifs...")
    random.seed(42)  # For reproducibility
    sequence = generate_sequence_with_motifs(size)
    patterns = get_test_patterns()
    print(f"Number of patterns: {len(patterns)}")
    print()
    
    results = {}
    
    # Benchmark pure Python regex
    print("Running pure Python regex benchmark...")
    try:
        regex_time, regex_matches = benchmark_regex(sequence, patterns, iterations)
        results['regex'] = {'time': regex_time, 'matches': regex_matches}
        print(f"  Time: {format_time(regex_time)}")
        print(f"  Matches: {regex_matches}")
        print(f"  Throughput: {format_throughput(size, regex_time)}")
    except Exception as e:
        print(f"  Error: {e}")
        results['regex'] = None
    print()
    
    # Benchmark fallback regex from hyperscan_backend
    print("Running fallback_regex_scan benchmark...")
    try:
        fallback_time, fallback_matches = benchmark_fallback_regex(sequence, patterns, iterations)
        results['fallback'] = {'time': fallback_time, 'matches': fallback_matches}
        print(f"  Time: {format_time(fallback_time)}")
        print(f"  Matches: {fallback_matches}")
        print(f"  Throughput: {format_throughput(size, fallback_time)}")
    except Exception as e:
        print(f"  Error: {e}")
        results['fallback'] = None
    print()
    
    # Benchmark Hyperscan (if available)
    print("Running Hyperscan benchmark...")
    if is_hyperscan_available():
        try:
            hs_time, hs_matches = benchmark_hyperscan(sequence, patterns, iterations)
            results['hyperscan'] = {'time': hs_time, 'matches': hs_matches}
            print(f"  Time: {format_time(hs_time)}")
            print(f"  Matches: {hs_matches}")
            print(f"  Throughput: {format_throughput(size, hs_time)}")
        except Exception as e:
            print(f"  Error: {e}")
            results['hyperscan'] = None
    else:
        print("  Hyperscan not available - skipping")
        print("  Install with: pip install hyperscan")
        print("  (requires libhs system library)")
        results['hyperscan'] = None
    print()
    
    # Summary comparison
    print("=" * 70)
    print("Summary")
    print("=" * 70)
    
    if results.get('regex') and results.get('hyperscan'):
        speedup = results['regex']['time'] / results['hyperscan']['time']
        print(f"\nHyperscan vs Pure Regex: {speedup:.1f}x faster")
    
    if results.get('fallback') and results.get('hyperscan'):
        speedup = results['fallback']['time'] / results['hyperscan']['time']
        print(f"Hyperscan vs Fallback Regex: {speedup:.1f}x faster")
    
    if not results.get('hyperscan'):
        print("\nNote: Hyperscan benchmark was not run.")
        print("To enable Hyperscan:")
        print("  1. Install libhs (Intel Hyperscan library)")
        print("     - Ubuntu/Debian: apt-get install libhyperscan-dev")
        print("     - macOS: brew install hyperscan")
        print("  2. Install Python bindings: pip install hyperscan")
    
    print()
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark Hyperscan vs Python regex for Non-B DNA detection"
    )
    parser.add_argument(
        "--size", "-s",
        type=int,
        default=100000,
        help="Sequence size in base pairs (default: 100000)"
    )
    parser.add_argument(
        "--iterations", "-i",
        type=int,
        default=5,
        help="Number of iterations per benchmark (default: 5)"
    )
    
    args = parser.parse_args()
    
    run_benchmark(
        size=args.size,
        iterations=args.iterations
    )


if __name__ == "__main__":
    main()
