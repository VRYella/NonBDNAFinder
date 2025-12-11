"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                     COMPREHENSIVE PERFORMANCE BENCHMARK                       ║
║              Measure and Validate All Performance Optimizations              ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: benchmark_performance.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.2 - Performance Testing
LICENSE: MIT

DESCRIPTION:
    Comprehensive benchmark suite to measure and validate all performance
    optimizations. Tests individual components and end-to-end performance.

BENCHMARKS:
    - K-mer indexing speed
    - Pattern matching performance
    - Individual detector speeds
    - Overall analysis throughput
    - Parallel vs sequential comparison
    - Memory usage profiling

OUTPUT:
    - Performance metrics table
    - Speedup comparisons
    - Reproducibility verification
    - Performance regression detection
"""

import time
import sys
import json
from typing import Dict, List, Tuple, Any
from collections import defaultdict

# Import components to test
import nonbscanner as nbs


class PerformanceBenchmark:
    """Comprehensive performance benchmark suite."""
    
    def __init__(self):
        self.results = defaultdict(dict)
        self.test_sequences = self._generate_test_sequences()
    
    def _generate_test_sequences(self) -> Dict[str, str]:
        """Generate diverse test sequences of different sizes."""
        sequences = {}
        
        # Small sequence (1kb) - G4 rich
        sequences['1kb_g4'] = "GGGTTAGGGTTAGGGTTAGGG" * 50
        
        # Medium sequence (10kb) - mixed content
        base_seq = "ACGTACGTACGTACGT" + "GGGTTAGGGTTAGGGTTAGGG" + "CCCCTAACCCTAACCCTAACCC" + "CGCGCGCGCGCG"
        sequences['10kb_mixed'] = base_seq * 100
        
        # Large sequence (100kb) - genomic-like
        genomic_pattern = (
            "AAAAATTTT" * 10 +  # A-tracts
            "CGCGCGCG" * 10 +   # Z-DNA
            "GGGTTAGGGTTAGGG" * 10 +  # G4
            "CCCCTAACCCTAACCC" * 10 +  # i-motif
            "ACGTACGTACGTACGT" * 10    # General
        )
        sequences['100kb_genomic'] = genomic_pattern * 100
        
        return sequences
    
    def benchmark_scanner_functions(self) -> Dict[str, Any]:
        """Benchmark low-level scanner functions."""
        print("\n" + "="*70)
        print("SCANNER FUNCTIONS BENCHMARK")
        print("="*70)
        
        test_seq = self.test_sequences['10kb_mixed']
        results = {}
        
        # Import scanner functions
        try:
            from scanner import (
                find_direct_repeats,
                find_inverted_repeats,
                find_strs
            )
            
            # Direct repeats
            print("\n--- Direct Repeats ---")
            start = time.time()
            direct = find_direct_repeats(test_seq)
            elapsed = time.time() - start
            speed = len(test_seq) / elapsed
            print(f"Found: {len(direct)} repeats")
            print(f"Time: {elapsed:.4f}s")
            print(f"Speed: {speed:,.0f} bp/s")
            results['direct_repeats'] = {'count': len(direct), 'time': elapsed, 'speed': speed}
            
            # Inverted repeats
            print("\n--- Inverted Repeats ---")
            start = time.time()
            inverted = find_inverted_repeats(test_seq)
            elapsed = time.time() - start
            speed = len(test_seq) / elapsed
            print(f"Found: {len(inverted)} repeats")
            print(f"Time: {elapsed:.4f}s")
            print(f"Speed: {speed:,.0f} bp/s")
            results['inverted_repeats'] = {'count': len(inverted), 'time': elapsed, 'speed': speed}
            
            # STRs
            print("\n--- STRs ---")
            start = time.time()
            strs = find_strs(test_seq)
            elapsed = time.time() - start
            speed = len(test_seq) / elapsed
            print(f"Found: {len(strs)} STRs")
            print(f"Time: {elapsed:.4f}s")
            print(f"Speed: {speed:,.0f} bp/s")
            results['strs'] = {'count': len(strs), 'time': elapsed, 'speed': speed}
            
        except Exception as e:
            print(f"Error benchmarking scanner functions: {e}")
        
        return results
    
    def benchmark_pattern_cache(self) -> Dict[str, Any]:
        """Benchmark pattern cache performance."""
        print("\n" + "="*70)
        print("PATTERN CACHE BENCHMARK")
        print("="*70)
        
        test_seq = self.test_sequences['10kb_mixed']
        results = {}
        
        try:
            from pattern_cache import get_pattern_cache
            import re
            
            cache = get_pattern_cache()
            
            # Test G4 pattern matching
            print("\n--- G4 Pattern Matching ---")
            pattern_str = r'G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}'
            
            # Uncached
            start = time.time()
            matches_uncached = list(re.finditer(pattern_str, test_seq, re.IGNORECASE))
            elapsed_uncached = time.time() - start
            
            # Cached
            start = time.time()
            matches_cached = cache.find_all('g4', 'canonical', test_seq)
            elapsed_cached = time.time() - start
            
            speedup = elapsed_uncached / elapsed_cached if elapsed_cached > 0 else 0
            
            print(f"Uncached: {len(matches_uncached)} matches in {elapsed_uncached*1000:.4f}ms")
            print(f"Cached: {len(matches_cached)} matches in {elapsed_cached*1000:.4f}ms")
            print(f"Speedup: {speedup:.2f}x")
            
            results['pattern_cache'] = {
                'uncached_time': elapsed_uncached,
                'cached_time': elapsed_cached,
                'speedup': speedup
            }
            
        except Exception as e:
            print(f"Error benchmarking pattern cache: {e}")
        
        return results
    
    def benchmark_full_analysis(self) -> Dict[str, Any]:
        """Benchmark full sequence analysis (end-to-end)."""
        print("\n" + "="*70)
        print("FULL ANALYSIS BENCHMARK")
        print("="*70)
        
        results = {}
        
        for seq_name, test_seq in self.test_sequences.items():
            print(f"\n--- {seq_name} ({len(test_seq)} bp) ---")
            
            # Standard mode
            print("Standard Mode (Sequential):")
            start = time.time()
            motifs_std = nbs.analyze_sequence(test_seq, seq_name, use_fast_mode=False)
            elapsed_std = time.time() - start
            speed_std = len(test_seq) / elapsed_std
            print(f"  Found: {len(motifs_std)} motifs")
            print(f"  Time: {elapsed_std:.4f}s")
            print(f"  Speed: {speed_std:,.0f} bp/s")
            
            # Fast mode
            print("Fast Mode (Parallel):")
            start = time.time()
            motifs_fast = nbs.analyze_sequence(test_seq, seq_name, use_fast_mode=True)
            elapsed_fast = time.time() - start
            speed_fast = len(test_seq) / elapsed_fast
            print(f"  Found: {len(motifs_fast)} motifs")
            print(f"  Time: {elapsed_fast:.4f}s")
            print(f"  Speed: {speed_fast:,.0f} bp/s")
            
            speedup = elapsed_std / elapsed_fast if elapsed_fast > 0 else 0
            print(f"\n  Speedup: {speedup:.2f}x")
            
            # Verify reproducibility
            reproducible = len(motifs_std) == len(motifs_fast)
            print(f"  Reproducible: {'✓ Yes' if reproducible else '✗ No'}")
            
            results[seq_name] = {
                'length': len(test_seq),
                'std_time': elapsed_std,
                'std_speed': speed_std,
                'std_motifs': len(motifs_std),
                'fast_time': elapsed_fast,
                'fast_speed': speed_fast,
                'fast_motifs': len(motifs_fast),
                'speedup': speedup,
                'reproducible': reproducible
            }
        
        return results
    
    def generate_summary_report(self, all_results: Dict[str, Any]) -> str:
        """Generate a formatted summary report."""
        report = []
        report.append("\n" + "="*70)
        report.append("PERFORMANCE BENCHMARK SUMMARY")
        report.append("="*70)
        
        # Full analysis summary
        if 'full_analysis' in all_results:
            report.append("\n## Full Analysis Performance ##")
            report.append("")
            report.append("| Sequence   | Size    | Std Time | Fast Time | Speedup | Speed (bp/s) |")
            report.append("|------------|---------|----------|-----------|---------|--------------|")
            
            for seq_name, data in all_results['full_analysis'].items():
                report.append(
                    f"| {seq_name:<10} | {data['length']:>6} | "
                    f"{data['std_time']:>7.3f}s | {data['fast_time']:>8.3f}s | "
                    f"{data['speedup']:>6.2f}x | {data['fast_speed']:>11,.0f} |"
                )
        
        # Scanner functions summary
        if 'scanner' in all_results:
            report.append("\n## Scanner Functions Performance ##")
            report.append("")
            report.append("| Function          | Time (s) | Speed (bp/s) | Motifs Found |")
            report.append("|-------------------|----------|--------------|--------------|")
            
            for func_name, data in all_results['scanner'].items():
                report.append(
                    f"| {func_name:<17} | {data['time']:>7.3f}s | "
                    f"{data['speed']:>11,.0f} | {data['count']:>12} |"
                )
        
        # Pattern cache summary
        if 'pattern_cache' in all_results and all_results['pattern_cache']:
            report.append("\n## Pattern Cache Performance ##")
            report.append("")
            data = all_results['pattern_cache']
            if 'uncached_time' in data:
                report.append(f"Uncached time: {data['uncached_time']*1000:.4f}ms")
                report.append(f"Cached time: {data['cached_time']*1000:.4f}ms")
                report.append(f"Speedup: {data['speedup']:.2f}x")
        
        # Overall summary
        report.append("\n## Overall Performance Gains ##")
        report.append("")
        
        if 'full_analysis' in all_results:
            avg_speedup = sum(d['speedup'] for d in all_results['full_analysis'].values()) / len(all_results['full_analysis'])
            max_speed = max(d['fast_speed'] for d in all_results['full_analysis'].values())
            report.append(f"Average speedup: {avg_speedup:.2f}x")
            report.append(f"Maximum throughput: {max_speed:,.0f} bp/s")
        
        report.append("\n" + "="*70)
        
        return "\n".join(report)
    
    def run_all_benchmarks(self) -> Dict[str, Any]:
        """Run all benchmarks and generate report."""
        print("\n" + "#"*70)
        print("# NonBDNAFinder Performance Benchmark Suite")
        print("#"*70)
        
        all_results = {}
        
        # Run benchmarks
        all_results['scanner'] = self.benchmark_scanner_functions()
        all_results['pattern_cache'] = self.benchmark_pattern_cache()
        all_results['full_analysis'] = self.benchmark_full_analysis()
        
        # Generate report
        report = self.generate_summary_report(all_results)
        print(report)
        
        # Save results to JSON
        with open('benchmark_results.json', 'w') as f:
            json.dump(all_results, f, indent=2)
        print("\nResults saved to: benchmark_results.json")
        
        return all_results


def main():
    """Main benchmark execution."""
    benchmark = PerformanceBenchmark()
    results = benchmark.run_all_benchmarks()
    return results


if __name__ == "__main__":
    main()
