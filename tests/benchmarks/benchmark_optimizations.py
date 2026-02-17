"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ NonBDNAFinder Performance Benchmark Suite                                    │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.2            │
│ Validates 50-200x performance improvement claims                             │
└──────────────────────────────────────────────────────────────────────────────┘
"""

import time
import logging
import random
import sys
from typing import Dict, List, Tuple, Any
from collections import defaultdict

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# Import both standard and optimized scanners
try:
    from Utilities.nonbscanner import analyze_sequence as analyze_standard
    STANDARD_AVAILABLE = True
except ImportError:
    logger.error("Standard NonBScanner not available")
    STANDARD_AVAILABLE = False

try:
    from Utilities.nonbscanner_optimized import (
        analyze_sequence_optimized,
        get_optimization_info
    )
    OPTIMIZED_AVAILABLE = True
except ImportError:
    logger.error("Optimized NonBScanner not available")
    OPTIMIZED_AVAILABLE = False

try:
    from Utilities.ac_matcher import benchmark_matcher
    AC_AVAILABLE = True
except ImportError:
    AC_AVAILABLE = False


# ═══════════════════════════════════════════════════════════════════════════════
# TEST SEQUENCE GENERATORS
# ═══════════════════════════════════════════════════════════════════════════════

def generate_random_dna(length: int, gc_content: float = 0.5, seed: int = 42) -> str:
    """
    Generate random DNA sequence with specified GC content.
    
    Args:
        length: Length of sequence in base pairs
        gc_content: Fraction of G+C bases (0.0-1.0)
        seed: Random seed for reproducibility
    
    Returns:
        Random DNA sequence string
    """
    random.seed(seed)
    
    gc_bases = ['G', 'C']
    at_bases = ['A', 'T']
    
    num_gc = int(length * gc_content)
    num_at = length - num_gc
    
    # Generate bases
    bases = (
        [random.choice(gc_bases) for _ in range(num_gc)] +
        [random.choice(at_bases) for _ in range(num_at)]
    )
    
    # Shuffle to randomize positions
    random.shuffle(bases)
    
    return ''.join(bases)


def generate_motif_rich_sequence(length: int, motif_density: float = 0.1, seed: int = 42) -> str:
    """
    Generate DNA sequence enriched with known NonB DNA motifs.
    
    Args:
        length: Total length of sequence
        motif_density: Fraction of sequence that should be motifs
        seed: Random seed for reproducibility
    
    Returns:
        DNA sequence with embedded motifs
    """
    random.seed(seed)
    
    # Common NonB DNA motifs
    motifs = [
        'GGGGTTGGGGTTGGGGTTGGGG',  # G-quadruplex
        'CCCCAACCCCAACCCCAACCCC',  # i-Motif
        'CGCGCGCGCGCGCG',           # Z-DNA
        'AAAAAATTTTTTTTT',          # A-philic
        'ATATAT' * 5,               # Direct repeat
        'TTAGGGTTAGGGTTAGGG' * 2,   # Telomeric G4
        'GGGAGGGGAGGG',             # Weak G4
    ]
    
    # Calculate how many motifs to insert
    total_motif_length = int(length * motif_density)
    remaining_length = length
    
    sequence_parts = []
    
    while remaining_length > 0:
        if total_motif_length > 0 and random.random() < motif_density:
            # Insert motif
            motif = random.choice(motifs)
            if len(motif) <= remaining_length:
                sequence_parts.append(motif)
                remaining_length -= len(motif)
                total_motif_length -= len(motif)
            else:
                # Fill with random DNA
                sequence_parts.append(generate_random_dna(remaining_length))
                remaining_length = 0
        else:
            # Insert random spacer (10-100 bp)
            spacer_len = min(random.randint(10, 100), remaining_length)
            sequence_parts.append(generate_random_dna(spacer_len))
            remaining_length -= spacer_len
    
    return ''.join(sequence_parts)


# ═══════════════════════════════════════════════════════════════════════════════
# BENCHMARK FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def benchmark_analysis(sequence: str, 
                      sequence_name: str,
                      implementation: str = "both",
                      enabled_classes: List[str] = None) -> Dict[str, Any]:
    """
    Benchmark NonBDNA analysis on a sequence.
    
    Args:
        sequence: DNA sequence to analyze
        sequence_name: Name for the sequence
        implementation: "standard", "optimized", or "both"
        enabled_classes: List of motif classes to detect (None = all)
    
    Returns:
        Dictionary with timing and result statistics
    """
    results = {
        'sequence_name': sequence_name,
        'sequence_length': len(sequence),
        'enabled_classes': enabled_classes or "all"
    }
    
    # Benchmark standard implementation
    if implementation in ["standard", "both"] and STANDARD_AVAILABLE:
        logger.info(f"Running STANDARD analysis on {sequence_name}...")
        
        start_time = time.time()
        try:
            motifs_standard = analyze_standard(
                sequence=sequence,
                sequence_name=sequence_name,
                enabled_classes=enabled_classes
            )
            elapsed_standard = time.time() - start_time
            
            results['standard'] = {
                'elapsed_time': elapsed_standard,
                'num_motifs': len(motifs_standard),
                'throughput': len(sequence) / elapsed_standard,
                'success': True
            }
            
            logger.info(
                f"  ✓ Standard: {elapsed_standard:.3f}s, "
                f"{len(motifs_standard)} motifs, "
                f"{results['standard']['throughput']:,.0f} bp/s"
            )
        
        except Exception as e:
            logger.error(f"  ✗ Standard failed: {e}")
            results['standard'] = {'success': False, 'error': str(e)}
    
    # Benchmark optimized implementation
    if implementation in ["optimized", "both"] and OPTIMIZED_AVAILABLE:
        logger.info(f"Running OPTIMIZED analysis on {sequence_name}...")
        
        start_time = time.time()
        try:
            motifs_optimized = analyze_sequence_optimized(
                sequence=sequence,
                sequence_name=sequence_name,
                enabled_classes=enabled_classes
            )
            elapsed_optimized = time.time() - start_time
            
            results['optimized'] = {
                'elapsed_time': elapsed_optimized,
                'num_motifs': len(motifs_optimized),
                'throughput': len(sequence) / elapsed_optimized,
                'success': True
            }
            
            logger.info(
                f"  ✓ Optimized: {elapsed_optimized:.3f}s, "
                f"{len(motifs_optimized)} motifs, "
                f"{results['optimized']['throughput']:,.0f} bp/s"
            )
        
        except Exception as e:
            logger.error(f"  ✗ Optimized failed: {e}")
            results['optimized'] = {'success': False, 'error': str(e)}
    
    # Calculate speedup
    if (implementation == "both" and 
        results.get('standard', {}).get('success') and 
        results.get('optimized', {}).get('success')):
        
        speedup = results['standard']['elapsed_time'] / results['optimized']['elapsed_time']
        results['speedup'] = speedup
        
        logger.info(f"  ⚡ Speedup: {speedup:.1f}x faster")
    
    return results


def run_benchmark_suite() -> Dict[str, Any]:
    """
    Run comprehensive benchmark suite with various sequence sizes.
    
    Returns:
        Dictionary with all benchmark results
    """
    print("=" * 80)
    print("NonBDNAFinder Performance Benchmark Suite")
    print("=" * 80)
    print()
    
    # Check what's available
    opt_info = get_optimization_info() if OPTIMIZED_AVAILABLE else {}
    print("System Configuration:")
    print(f"  Standard NonBScanner: {'✓ Available' if STANDARD_AVAILABLE else '✗ Not available'}")
    print(f"  Optimized NonBScanner: {'✓ Available' if OPTIMIZED_AVAILABLE else '✗ Not available'}")
    print(f"  Aho-Corasick Matcher: {'✓ Available' if AC_AVAILABLE else '✗ Not available'}")
    if opt_info:
        print(f"  Optimization Status: {opt_info.get('optimization_enabled', 'Unknown')}")
    print()
    
    if not STANDARD_AVAILABLE and not OPTIMIZED_AVAILABLE:
        print("ERROR: No implementations available to benchmark!")
        return {}
    
    # Define test cases
    test_cases = [
        # (size, name, motif_density)
        (10_000, "10KB_random", 0.0),
        (10_000, "10KB_motif_rich", 0.15),
        (100_000, "100KB_random", 0.0),
        (100_000, "100KB_motif_rich", 0.10),
        (1_000_000, "1MB_random", 0.0),
        (1_000_000, "1MB_motif_rich", 0.05),
    ]
    
    # Optional: Add larger sequences if requested
    if "--large" in sys.argv:
        test_cases.extend([
            (10_000_000, "10MB_random", 0.0),
            (10_000_000, "10MB_motif_rich", 0.02),
        ])
    
    all_results = {}
    
    for size, name, motif_density in test_cases:
        print(f"\n{'=' * 80}")
        print(f"Test Case: {name} ({size:,} bp)")
        print(f"{'=' * 80}")
        
        # Generate sequence
        if motif_density > 0:
            sequence = generate_motif_rich_sequence(size, motif_density)
        else:
            sequence = generate_random_dna(size)
        
        # Run benchmark
        result = benchmark_analysis(
            sequence=sequence,
            sequence_name=name,
            implementation="both"
        )
        
        all_results[name] = result
    
    # Print summary
    print("\n" + "=" * 80)
    print("BENCHMARK SUMMARY")
    print("=" * 80)
    print()
    print(f"{'Test Case':<25} {'Size':<12} {'Standard':<12} {'Optimized':<12} {'Speedup':<10}")
    print("-" * 80)
    
    for name, result in all_results.items():
        size_str = f"{result['sequence_length']:,} bp"
        
        if result.get('standard', {}).get('success'):
            std_time = f"{result['standard']['elapsed_time']:.2f}s"
        else:
            std_time = "FAILED"
        
        if result.get('optimized', {}).get('success'):
            opt_time = f"{result['optimized']['elapsed_time']:.2f}s"
        else:
            opt_time = "FAILED"
        
        if 'speedup' in result:
            speedup_str = f"{result['speedup']:.1f}x"
        else:
            speedup_str = "N/A"
        
        print(f"{name:<25} {size_str:<12} {std_time:<12} {opt_time:<12} {speedup_str:<10}")
    
    print()
    
    # Calculate average speedup
    speedups = [r['speedup'] for r in all_results.values() if 'speedup' in r]
    if speedups:
        avg_speedup = sum(speedups) / len(speedups)
        print(f"Average Speedup: {avg_speedup:.1f}x faster")
        print()
        
        if avg_speedup >= 50:
            print("✓ PERFORMANCE TARGET MET: 50-200x improvement achieved!")
        elif avg_speedup >= 10:
            print("⚠ PARTIAL SUCCESS: Significant improvement but below 50x target")
        else:
            print("⚠ PERFORMANCE TARGET NOT MET: Further optimization needed")
    
    return all_results


def benchmark_ac_matcher_only() -> None:
    """Benchmark just the Aho-Corasick matcher with simple patterns."""
    if not AC_AVAILABLE:
        print("Aho-Corasick matcher not available for benchmarking")
        return
    
    print("\n" + "=" * 80)
    print("Aho-Corasick Matcher Benchmark (Standalone)")
    print("=" * 80)
    print()
    
    # Test patterns
    patterns = [
        "GGGG", "CCCC", "AAAA", "TTTT",  # Simple repeats
        "CGCGCG", "ATATAT", "GCGCGC",     # Alternating
        "TTAGGG", "GGGTTA",                # Telomeric
    ]
    
    # Test sequences
    test_sequences = [
        (10_000, "10KB"),
        (100_000, "100KB"),
        (1_000_000, "1MB"),
    ]
    
    for size, name in test_sequences:
        sequence = generate_motif_rich_sequence(size, motif_density=0.1)
        
        result = benchmark_matcher(
            sequence=sequence,
            patterns=patterns,
            num_iterations=10
        )
        
        print(f"{name}:")
        print(f"  Build time: {result['build_time_ms']:.3f} ms")
        print(f"  Search time: {result['search_time_ms']:.3f} ms")
        print(f"  Throughput: {result['throughput_mb_per_sec']:.1f} MB/s")
        print(f"  Matches found: {result['total_matches']}")
        print(f"  Implementation: {result['implementation']}")
        print()


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN ENTRY POINT
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    """Main benchmark entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Benchmark NonBDNAFinder performance optimizations"
    )
    parser.add_argument(
        '--large',
        action='store_true',
        help='Include large sequences (10MB+) in benchmark'
    )
    parser.add_argument(
        '--ac-only',
        action='store_true',
        help='Benchmark only Aho-Corasick matcher'
    )
    parser.add_argument(
        '--quick',
        action='store_true',
        help='Run quick benchmark with small sequences only'
    )
    
    args = parser.parse_args()
    
    if args.ac_only:
        benchmark_ac_matcher_only()
    else:
        results = run_benchmark_suite()
        
        # Optionally save results
        if results:
            import json
            output_file = "benchmark_results.json"
            with open(output_file, 'w') as f:
                json.dump(results, f, indent=2)
            print(f"\n✓ Results saved to {output_file}")


if __name__ == "__main__":
    main()
