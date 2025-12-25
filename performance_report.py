#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    PERFORMANCE BENCHMARKING REPORT GENERATOR                  ║
║              Comprehensive Timing Analysis for NonBDNAFinder                 ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: performance_report.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2025.1
LICENSE: MIT

DESCRIPTION:
    Generates comprehensive performance reports for NonBDNAFinder by testing
    all detector functions with a 100,000 nucleotide sequence. Measures timing,
    memory usage, and throughput for each component.

FEATURES:
    - Individual detector timing
    - Memory usage monitoring
    - Throughput calculation (bp/s)
    - Detailed motif detection statistics
    - HTML/PDF report generation
    - Visualization of timing data

USAGE:
    python performance_report.py [--sequence-file FILE] [--output OUTPUT]
    
EXAMPLE:
    python performance_report.py --output performance_report.html
"""

import time
import sys
import os
import argparse
import tracemalloc
import json
from datetime import datetime
from typing import Dict, List, Any, Tuple
import random

# Ensure local modules can be imported
_current_dir = os.path.dirname(os.path.abspath(__file__))
if _current_dir not in sys.path:
    sys.path.insert(0, _current_dir)

# Import NonBDNAFinder modules
from nonbscanner import NonBScanner, DETECTOR_DISPLAY_NAMES, _get_cached_scanner
from detectors import (
    CurvedDNADetector,
    SlippedDNADetector,
    CruciformDetector,
    RLoopDetector,
    TriplexDetector,
    GQuadruplexDetector,
    IMotifDetector,
    ZDNADetector,
    APhilicDetector
)


# =============================================================================
# TEST SEQUENCE GENERATION
# =============================================================================

def generate_test_sequence(length: int = 100000, seed: int = 42) -> str:
    """
    Generate a realistic test DNA sequence with varied composition.
    
    Creates a sequence with:
    - Random AT/GC-rich regions
    - G-quadruplex motifs
    - Direct/inverted repeats
    - Curved DNA (A-tracts)
    
    Args:
        length: Length of sequence in nucleotides
        seed: Random seed for reproducibility
        
    Returns:
        DNA sequence string
    """
    random.seed(seed)
    nucleotides = ['A', 'T', 'G', 'C']
    sequence = []
    
    # Generate varied regions
    pos = 0
    while pos < length:
        region_type = random.choice(['random', 'at_rich', 'gc_rich', 'g4', 'repeat'])
        region_len = min(random.randint(100, 1000), length - pos)
        
        if region_type == 'random':
            # Random sequence
            region = ''.join(random.choices(nucleotides, k=region_len))
        elif region_type == 'at_rich':
            # AT-rich region (60-80% AT)
            region = ''.join(random.choices(['A', 'T', 'G', 'C'], 
                                           weights=[0.35, 0.35, 0.15, 0.15], 
                                           k=region_len))
        elif region_type == 'gc_rich':
            # GC-rich region (60-80% GC)
            region = ''.join(random.choices(['A', 'T', 'G', 'C'], 
                                           weights=[0.15, 0.15, 0.35, 0.35], 
                                           k=region_len))
        elif region_type == 'g4':
            # G-quadruplex motif
            g4_pattern = 'GGG' + ''.join(random.choices(['A', 'T', 'G', 'C'], k=5)) + \
                        'GGG' + ''.join(random.choices(['A', 'T', 'G', 'C'], k=5)) + \
                        'GGG' + ''.join(random.choices(['A', 'T', 'G', 'C'], k=5)) + 'GGG'
            remaining = region_len - len(g4_pattern)
            region = g4_pattern + ''.join(random.choices(nucleotides, k=max(0, remaining)))
        elif region_type == 'repeat':
            # Direct repeat
            unit_len = random.randint(5, 15)
            unit = ''.join(random.choices(nucleotides, k=unit_len))
            repeats = region_len // unit_len
            region = (unit * repeats)[:region_len]
        
        sequence.append(region)
        pos += region_len
    
    return ''.join(sequence)[:length]


# =============================================================================
# DETECTOR BENCHMARKING
# =============================================================================

class DetectorBenchmark:
    """Benchmark individual detector performance."""
    
    def __init__(self, sequence: str):
        self.sequence = sequence
        self.sequence_length = len(sequence)
        self.results = {}
        
    def benchmark_detector(self, detector_name: str, detector_class) -> Dict[str, Any]:
        """
        Benchmark a single detector.
        
        Args:
            detector_name: Name of detector
            detector_class: Detector class to instantiate
            
        Returns:
            Dictionary with timing and statistics
        """
        print(f"  Benchmarking {detector_name}...")
        
        # Start timing and memory tracking
        tracemalloc.start()
        start_time = time.time()
        start_memory = tracemalloc.get_traced_memory()[0]
        
        # Run detector
        try:
            detector = detector_class()
            motifs = detector.detect_motifs(self.sequence, "test_sequence")
            
            # End timing
            end_time = time.time()
            end_memory = tracemalloc.get_traced_memory()[0]
            
            elapsed = end_time - start_time
            memory_used = (end_memory - start_memory) / (1024 * 1024)  # Convert to MB
            throughput = self.sequence_length / elapsed if elapsed > 0 else 0
            
            # Calculate motif statistics
            motif_count = len(motifs)
            total_bp_covered = sum(m.get('Length', 0) for m in motifs)
            coverage_pct = (total_bp_covered / self.sequence_length) * 100 if self.sequence_length > 0 else 0
            
            # Subclass distribution
            subclasses = {}
            for motif in motifs:
                subclass = motif.get('Subclass', 'Unknown')
                subclasses[subclass] = subclasses.get(subclass, 0) + 1
            
            result = {
                'detector_name': detector_name,
                'elapsed_time': elapsed,
                'memory_mb': memory_used,
                'throughput_bps': throughput,
                'motif_count': motif_count,
                'total_bp_covered': total_bp_covered,
                'coverage_percent': coverage_pct,
                'subclasses': subclasses,
                'success': True,
                'error': None
            }
            
        except Exception as e:
            result = {
                'detector_name': detector_name,
                'elapsed_time': 0,
                'memory_mb': 0,
                'throughput_bps': 0,
                'motif_count': 0,
                'total_bp_covered': 0,
                'coverage_percent': 0,
                'subclasses': {},
                'success': False,
                'error': str(e)
            }
            print(f"    ERROR: {e}")
        
        finally:
            tracemalloc.stop()
        
        return result
    
    def benchmark_all_detectors(self) -> Dict[str, Dict[str, Any]]:
        """
        Benchmark all detectors.
        
        Returns:
            Dictionary mapping detector names to results
        """
        print(f"\nBenchmarking all detectors on {self.sequence_length:,} bp sequence...")
        print("=" * 70)
        
        detectors = {
            'Curved DNA': CurvedDNADetector,
            'Slipped DNA': SlippedDNADetector,
            'Cruciform': CruciformDetector,
            'R-Loop': RLoopDetector,
            'Triplex': TriplexDetector,
            'G-Quadruplex': GQuadruplexDetector,
            'i-Motif': IMotifDetector,
            'Z-DNA': ZDNADetector,
            'A-philic DNA': APhilicDetector
        }
        
        results = {}
        for name, detector_class in detectors.items():
            result = self.benchmark_detector(name, detector_class)
            results[name] = result
            
            if result['success']:
                print(f"    ✓ {name}: {result['elapsed_time']:.3f}s, "
                      f"{result['throughput_bps']:,.0f} bp/s, "
                      f"{result['motif_count']} motifs")
            else:
                print(f"    ✗ {name}: FAILED")
        
        print("=" * 70)
        return results
    
    def benchmark_full_pipeline(self) -> Dict[str, Any]:
        """
        Benchmark the full analysis pipeline.
        
        Returns:
            Dictionary with pipeline timing and statistics
        """
        print("\nBenchmarking full analysis pipeline...")
        print("=" * 70)
        
        # Start timing and memory tracking
        tracemalloc.start()
        start_time = time.time()
        start_memory = tracemalloc.get_traced_memory()[0]
        
        try:
            # Run full pipeline
            scanner = _get_cached_scanner()
            motifs = scanner.analyze_sequence(self.sequence, "test_sequence")
            
            # End timing
            end_time = time.time()
            end_memory = tracemalloc.get_traced_memory()[0]
            
            elapsed = end_time - start_time
            memory_used = (end_memory - start_memory) / (1024 * 1024)  # Convert to MB
            throughput = self.sequence_length / elapsed if elapsed > 0 else 0
            
            # Calculate statistics
            motif_count = len(motifs)
            total_bp_covered = sum(m.get('Length', 0) for m in motifs)
            coverage_pct = (total_bp_covered / self.sequence_length) * 100 if self.sequence_length > 0 else 0
            
            # Class and subclass distribution
            classes = {}
            subclasses = {}
            for motif in motifs:
                cls = motif.get('Class', 'Unknown')
                subclass = motif.get('Subclass', 'Unknown')
                classes[cls] = classes.get(cls, 0) + 1
                subclasses[subclass] = subclasses.get(subclass, 0) + 1
            
            result = {
                'elapsed_time': elapsed,
                'memory_mb': memory_used,
                'throughput_bps': throughput,
                'motif_count': motif_count,
                'total_bp_covered': total_bp_covered,
                'coverage_percent': coverage_pct,
                'classes': classes,
                'subclasses': subclasses,
                'success': True,
                'error': None
            }
            
            print(f"  ✓ Pipeline: {elapsed:.3f}s, {throughput:,.0f} bp/s, {motif_count} motifs")
            
        except Exception as e:
            result = {
                'elapsed_time': 0,
                'memory_mb': 0,
                'throughput_bps': 0,
                'motif_count': 0,
                'total_bp_covered': 0,
                'coverage_percent': 0,
                'classes': {},
                'subclasses': {},
                'success': False,
                'error': str(e)
            }
            print(f"  ✗ Pipeline: FAILED - {e}")
        
        finally:
            tracemalloc.stop()
        
        print("=" * 70)
        return result


# =============================================================================
# REPORT GENERATION
# =============================================================================

def generate_html_report(detector_results: Dict[str, Dict[str, Any]], 
                         pipeline_result: Dict[str, Any],
                         sequence_length: int,
                         output_file: str = "performance_report.html") -> None:
    """
    Generate an HTML performance report.
    
    Args:
        detector_results: Results from detector benchmarking
        pipeline_result: Results from pipeline benchmarking
        sequence_length: Length of test sequence
        output_file: Output HTML file path
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Calculate summary statistics
    total_detector_time = sum(r['elapsed_time'] for r in detector_results.values() if r['success'])
    total_motifs = sum(r['motif_count'] for r in detector_results.values() if r['success'])
    avg_throughput = sum(r['throughput_bps'] for r in detector_results.values() if r['success']) / len([r for r in detector_results.values() if r['success']])
    
    # Sort detectors by time
    sorted_detectors = sorted(detector_results.items(), key=lambda x: x[1]['elapsed_time'], reverse=True)
    
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>NonBDNAFinder Performance Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        h1 {{
            margin: 0;
            font-size: 2.5em;
        }}
        .timestamp {{
            margin-top: 10px;
            font-size: 0.9em;
            opacity: 0.9;
        }}
        .summary {{
            background: white;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .summary h2 {{
            margin-top: 0;
            color: #667eea;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }}
        .stat-card {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}
        .stat-label {{
            font-size: 0.9em;
            color: #666;
            margin-bottom: 5px;
        }}
        .stat-value {{
            font-size: 1.8em;
            font-weight: bold;
            color: #333;
        }}
        table {{
            width: 100%;
            background: white;
            border-collapse: collapse;
            margin: 20px 0;
            border-radius: 10px;
            overflow: hidden;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        th {{
            background: #667eea;
            color: white;
            padding: 15px;
            text-align: left;
            font-weight: 600;
        }}
        td {{
            padding: 12px 15px;
            border-bottom: 1px solid #eee;
        }}
        tr:hover {{
            background-color: #f8f9fa;
        }}
        .success {{
            color: #28a745;
            font-weight: bold;
        }}
        .failure {{
            color: #dc3545;
            font-weight: bold;
        }}
        .section {{
            background: white;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .section h2 {{
            margin-top: 0;
            color: #667eea;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
        }}
        .progress-bar {{
            width: 100%;
            height: 20px;
            background-color: #e9ecef;
            border-radius: 10px;
            overflow: hidden;
            margin: 10px 0;
        }}
        .progress-fill {{
            height: 100%;
            background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
            transition: width 0.3s ease;
        }}
        .chart {{
            margin: 20px 0;
        }}
        .bar {{
            background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
            height: 30px;
            margin: 5px 0;
            border-radius: 5px;
            position: relative;
            transition: width 0.3s ease;
        }}
        .bar-label {{
            position: absolute;
            left: 10px;
            top: 50%;
            transform: translateY(-50%);
            color: white;
            font-weight: 600;
            font-size: 0.9em;
        }}
        .footer {{
            text-align: center;
            padding: 20px;
            color: #666;
            font-size: 0.9em;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🔬 NonBDNAFinder Performance Report</h1>
        <div class="timestamp">Generated: {timestamp}</div>
        <div class="timestamp">Test Sequence: {sequence_length:,} nucleotides</div>
    </div>
    
    <div class="summary">
        <h2>📊 Summary Statistics</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-label">Total Detector Time</div>
                <div class="stat-value">{total_detector_time:.2f}s</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Pipeline Time</div>
                <div class="stat-value">{pipeline_result['elapsed_time']:.2f}s</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Total Motifs</div>
                <div class="stat-value">{total_motifs:,}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Avg Throughput</div>
                <div class="stat-value">{avg_throughput:,.0f} bp/s</div>
            </div>
        </div>
    </div>
    
    <div class="section">
        <h2>⏱️ Individual Detector Performance</h2>
        <table>
            <thead>
                <tr>
                    <th>Detector</th>
                    <th>Time (s)</th>
                    <th>Throughput (bp/s)</th>
                    <th>Memory (MB)</th>
                    <th>Motifs Found</th>
                    <th>Coverage (%)</th>
                    <th>Status</th>
                </tr>
            </thead>
            <tbody>
"""
    
    for detector_name, result in sorted_detectors:
        status_class = "success" if result['success'] else "failure"
        status_text = "✓ Success" if result['success'] else "✗ Failed"
        
        html_content += f"""
                <tr>
                    <td><strong>{detector_name}</strong></td>
                    <td>{result['elapsed_time']:.3f}</td>
                    <td>{result['throughput_bps']:,.0f}</td>
                    <td>{result['memory_mb']:.2f}</td>
                    <td>{result['motif_count']:,}</td>
                    <td>{result['coverage_percent']:.2f}%</td>
                    <td class="{status_class}">{status_text}</td>
                </tr>
"""
    
    html_content += """
            </tbody>
        </table>
    </div>
    
    <div class="section">
        <h2>📈 Timing Comparison</h2>
        <div class="chart">
"""
    
    # Add bar chart for timing
    max_time = max(r['elapsed_time'] for r in detector_results.values() if r['success']) if detector_results else 1
    for detector_name, result in sorted_detectors:
        if result['success']:
            width_pct = (result['elapsed_time'] / max_time) * 100
            html_content += f"""
            <div>
                <div style="font-size: 0.9em; color: #666; margin-bottom: 3px;">{detector_name}</div>
                <div class="bar" style="width: {width_pct}%;">
                    <span class="bar-label">{result['elapsed_time']:.3f}s</span>
                </div>
            </div>
"""
    
    html_content += """
        </div>
    </div>
    
    <div class="section">
        <h2>🔬 Full Pipeline Performance</h2>
        <table>
            <thead>
                <tr>
                    <th>Metric</th>
                    <th>Value</th>
                </tr>
            </thead>
            <tbody>
                <tr>
                    <td>Total Time</td>
                    <td><strong>{pipeline_time:.3f} seconds</strong></td>
                </tr>
                <tr>
                    <td>Throughput</td>
                    <td><strong>{pipeline_throughput:,.0f} bp/s</strong></td>
                </tr>
                <tr>
                    <td>Memory Used</td>
                    <td><strong>{pipeline_memory:.2f} MB</strong></td>
                </tr>
                <tr>
                    <td>Total Motifs</td>
                    <td><strong>{pipeline_motifs:,}</strong></td>
                </tr>
                <tr>
                    <td>Coverage</td>
                    <td><strong>{pipeline_coverage:.2f}%</strong></td>
                </tr>
                <tr>
                    <td>Status</td>
                    <td class="{pipeline_status_class}"><strong>{pipeline_status_text}</strong></td>
                </tr>
            </tbody>
        </table>
    </div>
    
    <div class="footer">
        NonBDNAFinder v2025.1 | Performance Report Generator<br>
        Dr. Venkata Rajesh Yella | KL University
    </div>
</body>
</html>
""".format(
        pipeline_time=pipeline_result['elapsed_time'],
        pipeline_throughput=pipeline_result['throughput_bps'],
        pipeline_memory=pipeline_result['memory_mb'],
        pipeline_motifs=pipeline_result['motif_count'],
        pipeline_coverage=pipeline_result['coverage_percent'],
        pipeline_status_class="success" if pipeline_result['success'] else "failure",
        pipeline_status_text="✓ Success" if pipeline_result['success'] else "✗ Failed"
    )
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"\n✓ HTML report saved to: {output_file}")


def generate_json_report(detector_results: Dict[str, Dict[str, Any]], 
                        pipeline_result: Dict[str, Any],
                        sequence_length: int,
                        output_file: str = "performance_report.json") -> None:
    """
    Generate a JSON performance report.
    
    Args:
        detector_results: Results from detector benchmarking
        pipeline_result: Results from pipeline benchmarking
        sequence_length: Length of test sequence
        output_file: Output JSON file path
    """
    report = {
        'timestamp': datetime.now().isoformat(),
        'sequence_length': sequence_length,
        'detector_results': detector_results,
        'pipeline_result': pipeline_result,
        'summary': {
            'total_detector_time': sum(r['elapsed_time'] for r in detector_results.values() if r['success']),
            'total_motifs_found': sum(r['motif_count'] for r in detector_results.values() if r['success']),
            'average_throughput': sum(r['throughput_bps'] for r in detector_results.values() if r['success']) / len([r for r in detector_results.values() if r['success']]) if detector_results else 0
        }
    }
    
    with open(output_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"✓ JSON report saved to: {output_file}")


# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    """Main function to run performance benchmarking."""
    parser = argparse.ArgumentParser(
        description='Generate performance report for NonBDNAFinder',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s
  %(prog)s --output my_report.html
  %(prog)s --sequence-file test.fasta --output report.html
  %(prog)s --length 50000 --output report.html
        """
    )
    parser.add_argument('--sequence-file', type=str, help='Input FASTA file (optional)')
    parser.add_argument('--length', type=int, default=100000, 
                       help='Length of test sequence if generating (default: 100000)')
    parser.add_argument('--output', type=str, default='performance_report.html',
                       help='Output HTML report file (default: performance_report.html)')
    parser.add_argument('--json', type=str, help='Output JSON report file (optional)')
    parser.add_argument('--seed', type=int, default=42, 
                       help='Random seed for sequence generation (default: 42)')
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("NonBDNAFinder Performance Benchmarking")
    print("=" * 70)
    
    # Get test sequence
    if args.sequence_file:
        print(f"\nLoading sequence from: {args.sequence_file}")
        # Simple FASTA parser
        with open(args.sequence_file, 'r') as f:
            lines = f.readlines()
            sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
        print(f"Loaded sequence: {len(sequence):,} bp")
    else:
        print(f"\nGenerating test sequence: {args.length:,} bp")
        sequence = generate_test_sequence(args.length, args.seed)
        print(f"Generated sequence: {len(sequence):,} bp")
    
    # Run benchmarks
    benchmark = DetectorBenchmark(sequence)
    
    # Benchmark individual detectors
    detector_results = benchmark.benchmark_all_detectors()
    
    # Benchmark full pipeline
    pipeline_result = benchmark.benchmark_full_pipeline()
    
    # Generate reports
    print("\nGenerating reports...")
    generate_html_report(detector_results, pipeline_result, len(sequence), args.output)
    
    if args.json:
        generate_json_report(detector_results, pipeline_result, len(sequence), args.json)
    
    print("\n✅ Performance benchmarking complete!")
    print("=" * 70)


if __name__ == '__main__':
    main()
