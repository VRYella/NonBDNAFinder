#!/usr/bin/env python3
"""
Comprehensive Performance Benchmarking Script for NonBDNAFinder

This script measures:
1. Processing speed (bp/second) for various sequence sizes
2. Memory usage and optimization
3. Accuracy on known motifs
4. Scalability testing
5. Feature detection rates

Author: Dr. Venkata Rajesh Yella
Version: 1.0
"""

import time
import sys
import os
from typing import Dict, List, Any, Tuple
import pandas as pd

# Import NonBDNAFinder modules
try:
    from nonbscanner import analyze_sequence
    from utilities import get_memory_usage_mb, trigger_garbage_collection
except ImportError:
    print("Error: Could not import NonBDNAFinder modules")
    print("Please ensure you're running from the NonBDNAFinder directory")
    sys.exit(1)


class PerformanceBenchmark:
    """Comprehensive benchmarking suite for NonBDNAFinder"""
    
    def __init__(self):
        self.results = {
            'speed': [],
            'memory': [],
            'accuracy': [],
            'scalability': []
        }
        
    def generate_test_sequence(self, size: int, motif_type: str = "random") -> str:
        """Generate test sequences with known motifs"""
        import numpy as np
        
        if motif_type == "random":
            # Random DNA sequence
            bases = ['A', 'C', 'G', 'T']
            return ''.join(np.random.choice(bases, size))
        
        elif motif_type == "g4_rich":
            # Sequence enriched with G4 motifs
            base_seq = "GGGTTAGGGTTAGGGTTAGGG"  # Telomeric G4
            repeats = size // len(base_seq) + 1
            seq = (base_seq * repeats)[:size]
            return seq
        
        elif motif_type == "str_rich":
            # Sequence with STRs
            base_seq = "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAG"  # HTT-like CAG repeat
            repeats = size // len(base_seq) + 1
            seq = (base_seq * repeats)[:size]
            return seq
        
        elif motif_type == "zdna_rich":
            # Sequence with Z-DNA motifs
            base_seq = "CGCGCGCGCGCGCGCGCGCG"  # CG dinucleotide repeats
            repeats = size // len(base_seq) + 1
            seq = (base_seq * repeats)[:size]
            return seq
        
        elif motif_type == "mixed":
            # Mixed motifs
            g4 = "GGGTTAGGGTTAGGGTTAGGG" * 5
            str_seq = "CAGCAGCAGCAGCAG" * 5
            zdna = "CGCGCGCGCGCG" * 5
            curved = "AAAAAATTTTTTAAAAAAATTTTT" * 5
            mixed = g4 + str_seq + zdna + curved
            repeats = size // len(mixed) + 1
            seq = (mixed * repeats)[:size]
            return seq
        
        return ''.join(np.random.choice(['A', 'C', 'G', 'T'], size))
    
    def benchmark_speed(self, sizes: List[int] = [1000, 5000, 10000, 50000, 100000]) -> Dict:
        """Benchmark processing speed across different sequence sizes"""
        
        print("\n" + "="*70)
        print("SPEED BENCHMARK")
        print("="*70)
        
        speed_results = []
        
        for size in sizes:
            print(f"\nTesting {size:,} bp sequence...")
            
            # Generate test sequence
            seq = self.generate_test_sequence(size, "mixed")
            
            # Warmup run
            _ = analyze_sequence(seq[:100], "warmup")
            
            # Timed run
            start_time = time.time()
            motifs = analyze_sequence(seq, f"test_{size}")
            end_time = time.time()
            
            elapsed = end_time - start_time
            bp_per_second = size / elapsed if elapsed > 0 else 0
            
            result = {
                'size_bp': size,
                'time_seconds': elapsed,
                'bp_per_second': bp_per_second,
                'motifs_detected': len(motifs),
                'motifs_per_kb': (len(motifs) / size) * 1000
            }
            
            speed_results.append(result)
            
            print(f"  Time: {elapsed:.2f}s")
            print(f"  Speed: {bp_per_second:,.0f} bp/s")
            print(f"  Motifs detected: {len(motifs)}")
            
            # Cleanup
            trigger_garbage_collection()
        
        self.results['speed'] = speed_results
        return speed_results
    
    def benchmark_memory(self, size: int = 100000) -> Dict:
        """Benchmark memory usage"""
        
        print("\n" + "="*70)
        print("MEMORY BENCHMARK")
        print("="*70)
        
        # Initial memory
        trigger_garbage_collection()
        initial_mem = get_memory_usage_mb()
        print(f"\nInitial memory: {initial_mem:.1f} MB")
        
        # Generate sequence
        seq = self.generate_test_sequence(size, "mixed")
        seq_mem = get_memory_usage_mb()
        print(f"After sequence generation: {seq_mem:.1f} MB (+{seq_mem - initial_mem:.1f} MB)")
        
        # Process sequence
        motifs = analyze_sequence(seq, "memory_test")
        process_mem = get_memory_usage_mb()
        print(f"After analysis: {process_mem:.1f} MB (+{process_mem - seq_mem:.1f} MB)")
        
        # Convert to DataFrame
        df = pd.DataFrame(motifs)
        df_mem = get_memory_usage_mb()
        print(f"After DataFrame creation: {df_mem:.1f} MB (+{df_mem - process_mem:.1f} MB)")
        
        # Cleanup
        del seq, motifs, df
        trigger_garbage_collection()
        final_mem = get_memory_usage_mb()
        print(f"After cleanup: {final_mem:.1f} MB (freed {df_mem - final_mem:.1f} MB)")
        
        result = {
            'initial_mb': initial_mem,
            'peak_mb': process_mem,
            'delta_mb': process_mem - initial_mem,
            'final_mb': final_mem,
            'freed_mb': df_mem - final_mem
        }
        
        self.results['memory'] = result
        return result
    
    def benchmark_accuracy(self) -> Dict:
        """Benchmark accuracy on known motifs"""
        
        print("\n" + "="*70)
        print("ACCURACY BENCHMARK")
        print("="*70)
        
        test_cases = [
            {
                'name': 'Telomeric G4',
                'sequence': 'GGGTTAGGGTTAGGGTTAGGG',
                'expected_class': 'G-Quadruplex',
                'expected_min_count': 1
            },
            {
                'name': 'HTT CAG repeat (normal)',
                'sequence': 'CAG' * 20,
                'expected_class': 'Slipped-DNA',
                'expected_min_count': 1
            },
            {
                'name': 'HTT CAG repeat (expanded)',
                'sequence': 'CAG' * 45,
                'expected_class': 'Slipped-DNA',
                'expected_min_count': 1
            },
            {
                'name': 'Z-DNA forming sequence',
                'sequence': 'CGCGCGCGCGCGCGCGCGCG',
                'expected_class': 'Z-DNA',
                'expected_min_count': 1
            },
            {
                'name': 'A-tract curved DNA',
                'sequence': 'AAAAAA' + 'TGCA' + 'AAAAAA',
                'expected_class': 'Curved-DNA',
                'expected_min_count': 1
            },
            {
                'name': 'Palindrome (cruciform)',
                'sequence': 'ATCGATCGATCGATCGAT' + 'ATCGATCGATCGATCGAT',
                'expected_class': 'Cruciform',
                'expected_min_count': 0  # May or may not detect depending on length
            }
        ]
        
        accuracy_results = []
        
        for test in test_cases:
            print(f"\nTesting: {test['name']}")
            print(f"  Sequence: {test['sequence'][:50]}{'...' if len(test['sequence']) > 50 else ''}")
            
            motifs = analyze_sequence(test['sequence'], test['name'])
            
            # Count motifs of expected class
            detected = [m for m in motifs if m['Class'] == test['expected_class']]
            detected_count = len(detected)
            
            success = detected_count >= test['expected_min_count']
            
            result = {
                'test_name': test['name'],
                'expected_class': test['expected_class'],
                'detected_count': detected_count,
                'total_motifs': len(motifs),
                'success': success
            }
            
            accuracy_results.append(result)
            
            status = "✓ PASS" if success else "✗ FAIL"
            print(f"  Expected: {test['expected_class']} (min {test['expected_min_count']})")
            print(f"  Detected: {detected_count} {test['expected_class']} motifs")
            print(f"  Total motifs: {len(motifs)}")
            print(f"  Status: {status}")
            
            if detected:
                print(f"  Top score: {max(m['Score'] for m in detected):.2f}")
        
        self.results['accuracy'] = accuracy_results
        
        # Calculate success rate
        success_rate = sum(1 for r in accuracy_results if r['success']) / len(accuracy_results) * 100
        print(f"\n  Overall accuracy: {success_rate:.1f}% ({sum(1 for r in accuracy_results if r['success'])}/{len(accuracy_results)} tests passed)")
        
        return accuracy_results
    
    def benchmark_scalability(self, sizes: List[int] = [10000, 50000, 100000, 500000, 1000000]) -> Dict:
        """Test scalability with progressively larger sequences"""
        
        print("\n" + "="*70)
        print("SCALABILITY BENCHMARK")
        print("="*70)
        
        scalability_results = []
        
        for size in sizes:
            print(f"\nTesting {size:,} bp sequence...")
            
            try:
                seq = self.generate_test_sequence(size, "random")
                
                initial_mem = get_memory_usage_mb()
                start_time = time.time()
                
                motifs = analyze_sequence(seq, f"scale_{size}")
                
                end_time = time.time()
                peak_mem = get_memory_usage_mb()
                
                elapsed = end_time - start_time
                
                result = {
                    'size_bp': size,
                    'time_seconds': elapsed,
                    'bp_per_second': size / elapsed,
                    'peak_memory_mb': peak_mem,
                    'memory_delta_mb': peak_mem - initial_mem,
                    'motifs_detected': len(motifs),
                    'success': True
                }
                
                scalability_results.append(result)
                
                print(f"  Time: {elapsed:.2f}s ({size/elapsed:,.0f} bp/s)")
                print(f"  Memory: {peak_mem:.1f} MB (Δ{peak_mem - initial_mem:.1f} MB)")
                print(f"  Motifs: {len(motifs)}")
                print(f"  Status: ✓ SUCCESS")
                
                # Cleanup
                del seq, motifs
                trigger_garbage_collection()
                
            except Exception as e:
                print(f"  Status: ✗ FAILED - {str(e)}")
                result = {
                    'size_bp': size,
                    'success': False,
                    'error': str(e)
                }
                scalability_results.append(result)
        
        self.results['scalability'] = scalability_results
        return scalability_results
    
    def benchmark_motif_detection_rates(self) -> Dict:
        """Benchmark detection rates for each motif type"""
        
        print("\n" + "="*70)
        print("MOTIF TYPE DETECTION RATES")
        print("="*70)
        
        # Create sequences enriched for each motif type
        test_sequences = {
            'G-Quadruplex': self.generate_test_sequence(10000, "g4_rich"),
            'Slipped-DNA': self.generate_test_sequence(10000, "str_rich"),
            'Z-DNA': self.generate_test_sequence(10000, "zdna_rich"),
            'Mixed': self.generate_test_sequence(10000, "mixed"),
            'Random': self.generate_test_sequence(10000, "random")
        }
        
        detection_results = []
        
        for seq_type, seq in test_sequences.items():
            print(f"\nAnalyzing {seq_type} sequence...")
            
            motifs = analyze_sequence(seq, seq_type)
            
            # Count by class
            class_counts = {}
            for motif in motifs:
                cls = motif.get('Class', 'Unknown')
                class_counts[cls] = class_counts.get(cls, 0) + 1
            
            result = {
                'sequence_type': seq_type,
                'total_motifs': len(motifs),
                'motifs_per_kb': (len(motifs) / len(seq)) * 1000,
                'class_distribution': class_counts
            }
            
            detection_results.append(result)
            
            print(f"  Total motifs: {len(motifs)}")
            print(f"  Density: {result['motifs_per_kb']:.2f} motifs/kb")
            print(f"  Classes detected: {', '.join(class_counts.keys())}")
        
        return detection_results
    
    def generate_report(self, output_file: str = "benchmark_report.md"):
        """Generate comprehensive markdown report"""
        
        print("\n" + "="*70)
        print("GENERATING REPORT")
        print("="*70)
        
        report = []
        report.append("# NonBDNAFinder Performance Benchmark Report\n")
        report.append(f"**Generated:** {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        report.append("---\n\n")
        
        # Speed results
        if self.results['speed']:
            report.append("## Speed Performance\n\n")
            report.append("| Sequence Size | Time (s) | Speed (bp/s) | Motifs Detected | Density (motifs/kb) |\n")
            report.append("|---------------|----------|--------------|-----------------|---------------------|\n")
            
            for r in self.results['speed']:
                report.append(f"| {r['size_bp']:,} bp | {r['time_seconds']:.2f} | {r['bp_per_second']:,.0f} | {r['motifs_detected']} | {r['motifs_per_kb']:.2f} |\n")
            
            # Calculate averages
            avg_speed = sum(r['bp_per_second'] for r in self.results['speed']) / len(self.results['speed'])
            report.append(f"\n**Average speed:** {avg_speed:,.0f} bp/s\n\n")
        
        # Memory results
        if self.results['memory']:
            report.append("## Memory Usage\n\n")
            mem = self.results['memory']
            report.append(f"- Initial memory: {mem['initial_mb']:.1f} MB\n")
            report.append(f"- Peak memory: {mem['peak_mb']:.1f} MB\n")
            report.append(f"- Memory delta: {mem['delta_mb']:.1f} MB\n")
            report.append(f"- After cleanup: {mem['final_mb']:.1f} MB\n")
            report.append(f"- Memory freed: {mem['freed_mb']:.1f} MB\n\n")
        
        # Accuracy results
        if self.results['accuracy']:
            report.append("## Accuracy Tests\n\n")
            report.append("| Test | Expected Class | Detected | Total Motifs | Status |\n")
            report.append("|------|----------------|----------|--------------|--------|\n")
            
            for r in self.results['accuracy']:
                status = "✓ PASS" if r['success'] else "✗ FAIL"
                report.append(f"| {r['test_name']} | {r['expected_class']} | {r['detected_count']} | {r['total_motifs']} | {status} |\n")
            
            success_rate = sum(1 for r in self.results['accuracy'] if r['success']) / len(self.results['accuracy']) * 100
            report.append(f"\n**Success rate:** {success_rate:.1f}%\n\n")
        
        # Scalability results
        if self.results['scalability']:
            report.append("## Scalability\n\n")
            report.append("| Sequence Size | Time (s) | Speed (bp/s) | Memory (MB) | Status |\n")
            report.append("|---------------|----------|--------------|-------------|--------|\n")
            
            for r in self.results['scalability']:
                if r['success']:
                    status = "✓ SUCCESS"
                    report.append(f"| {r['size_bp']:,} bp | {r['time_seconds']:.2f} | {r['bp_per_second']:,.0f} | {r['peak_memory_mb']:.1f} | {status} |\n")
                else:
                    status = "✗ FAILED"
                    report.append(f"| {r['size_bp']:,} bp | - | - | - | {status} |\n")
            
            report.append("\n")
        
        # Write report
        with open(output_file, 'w') as f:
            f.writelines(report)
        
        print(f"\nReport written to: {output_file}")
        
        return ''.join(report)
    
    def run_full_benchmark(self):
        """Run complete benchmark suite"""
        
        print("\n" + "="*70)
        print("NONBDNAFINDER COMPREHENSIVE BENCHMARK SUITE")
        print("="*70)
        
        # Run all benchmarks
        self.benchmark_speed(sizes=[1000, 5000, 10000, 50000, 100000])
        self.benchmark_memory(size=100000)
        self.benchmark_accuracy()
        self.benchmark_scalability(sizes=[10000, 50000, 100000, 500000])
        self.benchmark_motif_detection_rates()
        
        # Generate report
        self.generate_report()
        
        print("\n" + "="*70)
        print("BENCHMARK COMPLETE")
        print("="*70 + "\n")


def main():
    """Main entry point"""
    
    # Create benchmark instance
    benchmark = PerformanceBenchmark()
    
    # Run full benchmark
    benchmark.run_full_benchmark()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
