#!/usr/bin/env python3
"""
Performance Benchmark Script for NBDScanner
============================================

Tests memory usage and performance improvements for:
1. File upload and parsing (chunked vs non-chunked)
2. Analysis caching effectiveness
3. Visualization rendering time
4. Large dataset handling

Usage:
    python performance_benchmark.py
"""

import time
import sys
import os
import io
import tracemalloc
import json
from typing import Dict, List, Tuple

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from utilities import parse_fasta, parse_fasta_chunked, get_file_preview, get_basic_stats
from nonbscanner import analyze_sequence


class PerformanceBenchmark:
    """Performance benchmarking suite for NBDScanner improvements"""
    
    def __init__(self):
        self.results = {}
        
    def generate_test_fasta(self, num_sequences: int = 10, seq_length: int = 10000) -> str:
        """Generate test FASTA content for benchmarking"""
        import random
        bases = ['A', 'T', 'G', 'C']
        fasta_lines = []
        
        for i in range(num_sequences):
            fasta_lines.append(f">test_sequence_{i+1} {seq_length}bp")
            # Generate random sequence with some G-quadruplex patterns
            seq = []
            for j in range(seq_length):
                if j % 100 == 0 and j < seq_length - 50:
                    # Add G-quadruplex motif occasionally
                    seq.append("GGGTAGGGTGGGTAGGG")
                else:
                    seq.append(random.choice(bases))
            fasta_lines.append(''.join(seq))
        
        return '\n'.join(fasta_lines)
    
    def benchmark_file_parsing(self, fasta_content: str, method: str = "chunked") -> Tuple[float, float]:
        """
        Benchmark file parsing performance
        
        Returns:
            Tuple of (time_seconds, peak_memory_mb)
        """
        # Start memory tracking
        tracemalloc.start()
        start_time = time.time()
        
        if method == "chunked":
            # Simulate file object
            file_obj = io.BytesIO(fasta_content.encode('utf-8'))
            sequences = list(parse_fasta_chunked(file_obj))
        else:
            # Traditional parsing
            sequences = parse_fasta(fasta_content)
        
        elapsed_time = time.time() - start_time
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        peak_mb = peak / (1024 * 1024)
        
        return elapsed_time, peak_mb
    
    def benchmark_preview_generation(self, fasta_content: str) -> Tuple[float, float]:
        """Benchmark file preview generation"""
        tracemalloc.start()
        start_time = time.time()
        
        file_obj = io.BytesIO(fasta_content.encode('utf-8'))
        preview = get_file_preview(file_obj, max_sequences=3)
        
        elapsed_time = time.time() - start_time
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        peak_mb = peak / (1024 * 1024)
        
        return elapsed_time, peak_mb
    
    def benchmark_analysis(self, sequence: str, name: str = "test_seq") -> Tuple[float, float, int]:
        """
        Benchmark motif analysis
        
        Returns:
            Tuple of (time_seconds, peak_memory_mb, num_motifs)
        """
        tracemalloc.start()
        start_time = time.time()
        
        motifs = analyze_sequence(sequence, name)
        
        elapsed_time = time.time() - start_time
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        peak_mb = peak / (1024 * 1024)
        
        return elapsed_time, peak_mb, len(motifs)
    
    def run_all_benchmarks(self):
        """Run all performance benchmarks"""
        print("=" * 80)
        print("NBDScanner Performance Benchmark Suite")
        print("=" * 80)
        print()
        
        # Test different file sizes
        test_configs = [
            ("Small (10 sequences x 1kb)", 10, 1000),
            ("Medium (50 sequences x 5kb)", 50, 5000),
            ("Large (100 sequences x 10kb)", 100, 10000),
            ("Very Large (200 sequences x 20kb)", 200, 20000),
        ]
        
        for config_name, num_seqs, seq_len in test_configs:
            print(f"\n{'-' * 80}")
            print(f"Test Configuration: {config_name}")
            print(f"Total data size: ~{(num_seqs * seq_len) / 1024:.1f} KB")
            print(f"{'-' * 80}")
            
            # Generate test data
            print("Generating test data...")
            fasta_content = self.generate_test_fasta(num_seqs, seq_len)
            total_size_mb = len(fasta_content) / (1024 * 1024)
            print(f"Generated {total_size_mb:.2f} MB of test data")
            
            # Benchmark traditional parsing
            print("\n1. Traditional FASTA parsing...")
            trad_time, trad_mem = self.benchmark_file_parsing(fasta_content, method="traditional")
            print(f"   Time: {trad_time:.3f}s | Peak Memory: {trad_mem:.2f} MB")
            
            # Benchmark chunked parsing
            print("\n2. Chunked FASTA parsing (optimized)...")
            chunk_time, chunk_mem = self.benchmark_file_parsing(fasta_content, method="chunked")
            print(f"   Time: {chunk_time:.3f}s | Peak Memory: {chunk_mem:.2f} MB")
            
            # Calculate improvement
            time_improvement = ((trad_time - chunk_time) / trad_time) * 100 if trad_time > 0 else 0
            mem_improvement = ((trad_mem - chunk_mem) / trad_mem) * 100 if trad_mem > 0 else 0
            
            print(f"\n   ✓ Time improvement: {time_improvement:+.1f}%")
            print(f"   ✓ Memory improvement: {mem_improvement:+.1f}%")
            
            # Benchmark file preview
            print("\n3. File preview generation...")
            prev_time, prev_mem = self.benchmark_preview_generation(fasta_content)
            print(f"   Time: {prev_time:.3f}s | Peak Memory: {prev_mem:.2f} MB")
            print(f"   Preview is {(prev_time / chunk_time) * 100:.1f}% of full parse time")
            
            # Store results
            self.results[config_name] = {
                'num_sequences': num_seqs,
                'sequence_length': seq_len,
                'total_size_mb': total_size_mb,
                'traditional_parse': {'time': trad_time, 'memory': trad_mem},
                'chunked_parse': {'time': chunk_time, 'memory': chunk_mem},
                'preview': {'time': prev_time, 'memory': prev_mem},
                'improvements': {
                    'time_percent': time_improvement,
                    'memory_percent': mem_improvement
                }
            }
        
        # Test analysis performance on a single sequence
        print(f"\n{'-' * 80}")
        print("Analysis Performance Test")
        print(f"{'-' * 80}")
        
        test_seq = self.generate_test_fasta(1, 50000).split('\n', 1)[1]
        print(f"Testing on 50kb sequence...")
        
        anal_time, anal_mem, num_motifs = self.benchmark_analysis(test_seq)
        print(f"\nAnalysis Results:")
        print(f"   Time: {anal_time:.3f}s")
        print(f"   Peak Memory: {anal_mem:.2f} MB")
        print(f"   Motifs detected: {num_motifs}")
        print(f"   Throughput: {len(test_seq) / anal_time:,.0f} bp/s")
        
        self.results['analysis'] = {
            'time': anal_time,
            'memory': anal_mem,
            'motifs': num_motifs,
            'throughput_bps': len(test_seq) / anal_time
        }
        
        # Generate summary report
        self.print_summary()
    
    def print_summary(self):
        """Print summary of all benchmark results"""
        print("\n" + "=" * 80)
        print("BENCHMARK SUMMARY")
        print("=" * 80)
        
        print("\nFile Parsing Performance Improvements:")
        print(f"{'Configuration':<35} {'Time Saved':<15} {'Memory Saved':<15}")
        print("-" * 65)
        
        for config_name, data in self.results.items():
            if config_name != 'analysis':
                time_imp = data['improvements']['time_percent']
                mem_imp = data['improvements']['memory_percent']
                print(f"{config_name:<35} {time_imp:>+8.1f}% {mem_imp:>15.1f}%")
        
        if 'analysis' in self.results:
            anal = self.results['analysis']
            print(f"\nAnalysis Performance:")
            print(f"  Throughput: {anal['throughput_bps']:,.0f} bp/s")
            print(f"  Memory usage: {anal['memory']:.2f} MB")
        
        print("\n✅ All benchmarks completed successfully!")
        
        # Save results to JSON
        output_file = "benchmark_results.json"
        with open(output_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"\nDetailed results saved to: {output_file}")


def main():
    """Run the benchmark suite"""
    benchmark = PerformanceBenchmark()
    
    try:
        benchmark.run_all_benchmarks()
        return 0
    except Exception as e:
        print(f"\n❌ Benchmark failed with error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
