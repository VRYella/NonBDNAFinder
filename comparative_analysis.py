#!/usr/bin/env python3
"""
Comparative Genome Analysis Module
==================================

Performs comparative analysis across multiple genomes to identify:
- Motif distribution differences
- Statistical significance of patterns
- Cross-genome comparative metrics
- Enrichment and depletion patterns

Usage:
    python comparative_analysis.py [--input-dir INPUT_DIR] [--output-dir OUTPUT_DIR]
"""

import os
import sys
import json
import argparse
from pathlib import Path
from typing import List, Dict, Any
from collections import Counter, defaultdict
import numpy as np
from scipy import stats

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def load_batch_results(input_dir: str) -> Dict[str, Any]:
    """
    Load batch analysis results from directory.
    
    Args:
        input_dir: Directory containing batch results
        
    Returns:
        Dictionary with genome results
    """
    summary_file = os.path.join(input_dir, 'batch_summary.json')
    
    if not os.path.exists(summary_file):
        raise FileNotFoundError(f"Batch summary not found: {summary_file}")
    
    with open(summary_file, 'r') as f:
        summary = json.load(f)
    
    # Load individual genome motif files
    results = {}
    for genome_info in summary['genomes']:
        genome_file = genome_info['genome_file']
        safe_name = genome_file.replace('.fna', '').replace('.fasta', '').replace(' ', '_')
        motifs_file = os.path.join(input_dir, f"{safe_name}_motifs.json")
        
        if os.path.exists(motifs_file):
            with open(motifs_file, 'r') as f:
                motifs = json.load(f)
            
            results[genome_file] = {
                'info': genome_info,
                'motifs': motifs
            }
    
    return results


def calculate_genome_statistics(genome_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate comprehensive statistics for a single genome.
    
    Args:
        genome_data: Genome data with motifs
        
    Returns:
        Dictionary of statistics
    """
    motifs = genome_data['motifs']
    sequences = genome_data['info']['sequences']
    
    # Calculate total genome length
    total_length = sum(seq['length'] for seq in sequences)
    
    # Count motifs by class and subclass
    class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    subclass_counts = Counter(m.get('Subclass', 'Unknown') for m in motifs)
    
    # Calculate coverage
    total_motif_length = sum(m.get('Length', 0) for m in motifs)
    coverage_percent = (total_motif_length / total_length * 100) if total_length > 0 else 0
    
    # Motif density (per kb)
    density_per_kb = (len(motifs) / total_length * 1000) if total_length > 0 else 0
    
    # Average motif length
    avg_length = np.mean([m.get('Length', 0) for m in motifs]) if motifs else 0
    
    # Score distribution
    scores = [m.get('Score', 0) for m in motifs if m.get('Score', 0) > 0]
    avg_score = np.mean(scores) if scores else 0
    
    # Class diversity (Shannon entropy)
    class_diversity = calculate_shannon_entropy(class_counts)
    
    return {
        'total_length': total_length,
        'total_motifs': len(motifs),
        'unique_classes': len(class_counts),
        'unique_subclasses': len(subclass_counts),
        'coverage_percent': coverage_percent,
        'density_per_kb': density_per_kb,
        'avg_motif_length': avg_length,
        'avg_score': avg_score,
        'class_diversity': class_diversity,
        'class_counts': dict(class_counts),
        'subclass_counts': dict(subclass_counts)
    }


def calculate_shannon_entropy(counts: Counter) -> float:
    """
    Calculate Shannon entropy for diversity measurement.
    
    Args:
        counts: Counter object with class/subclass counts
        
    Returns:
        Shannon entropy value
    """
    total = sum(counts.values())
    if total == 0:
        return 0.0
    
    entropy = 0.0
    for count in counts.values():
        if count > 0:
            proportion = count / total
            entropy -= proportion * np.log2(proportion)
    
    return entropy


def compare_genomes(results: Dict[str, Any]) -> Dict[str, Any]:
    """
    Perform cross-genome comparative analysis.
    
    Args:
        results: Dictionary of genome results
        
    Returns:
        Comparative analysis results
    """
    # Calculate statistics for each genome
    genome_stats = {}
    for genome_name, genome_data in results.items():
        genome_stats[genome_name] = calculate_genome_statistics(genome_data)
    
    # Collect all motif classes across genomes
    all_classes = set()
    for stats in genome_stats.values():
        all_classes.update(stats['class_counts'].keys())
    
    # Calculate normalized class frequencies
    class_frequencies = {}
    for class_name in all_classes:
        frequencies = []
        for genome_name, stats in genome_stats.items():
            count = stats['class_counts'].get(class_name, 0)
            total = stats['total_motifs']
            freq = (count / total) if total > 0 else 0
            frequencies.append({
                'genome': genome_name,
                'frequency': freq,
                'count': count
            })
        
        class_frequencies[class_name] = frequencies
    
    # Identify enriched/depleted classes
    enrichment_patterns = analyze_enrichment(class_frequencies, genome_stats)
    
    # Cross-genome correlations
    correlations = calculate_cross_genome_correlations(genome_stats, all_classes)
    
    return {
        'genome_statistics': genome_stats,
        'class_frequencies': class_frequencies,
        'enrichment_patterns': enrichment_patterns,
        'correlations': correlations,
        'summary': generate_comparative_summary(genome_stats)
    }


def analyze_enrichment(class_frequencies: Dict[str, List], 
                       genome_stats: Dict[str, Any]) -> Dict[str, Any]:
    """
    Analyze enrichment/depletion patterns across genomes.
    
    Args:
        class_frequencies: Class frequencies across genomes
        genome_stats: Statistics for each genome
        
    Returns:
        Enrichment analysis results
    """
    enrichment = {}
    
    for class_name, frequencies in class_frequencies.items():
        freq_values = [f['frequency'] for f in frequencies]
        
        if len(freq_values) < 2:
            continue
        
        # Calculate statistics
        mean_freq = np.mean(freq_values)
        std_freq = np.std(freq_values)
        cv = (std_freq / mean_freq) if mean_freq > 0 else 0
        
        # Identify outliers (> 2 SD from mean)
        enriched_genomes = []
        depleted_genomes = []
        
        for freq_data in frequencies:
            z_score = ((freq_data['frequency'] - mean_freq) / std_freq) if std_freq > 0 else 0
            
            if z_score > 2.0:
                enriched_genomes.append({
                    'genome': freq_data['genome'],
                    'frequency': freq_data['frequency'],
                    'z_score': z_score
                })
            elif z_score < -2.0:
                depleted_genomes.append({
                    'genome': freq_data['genome'],
                    'frequency': freq_data['frequency'],
                    'z_score': z_score
                })
        
        enrichment[class_name] = {
            'mean_frequency': mean_freq,
            'std_frequency': std_freq,
            'coefficient_variation': cv,
            'enriched_genomes': enriched_genomes,
            'depleted_genomes': depleted_genomes
        }
    
    return enrichment


def calculate_cross_genome_correlations(genome_stats: Dict[str, Any],
                                        all_classes: set) -> Dict[str, Any]:
    """
    Calculate correlations between genome features.
    
    Args:
        genome_stats: Statistics for each genome
        all_classes: Set of all motif classes
        
    Returns:
        Correlation analysis results
    """
    # Create feature matrix
    genome_names = list(genome_stats.keys())
    
    # Features to correlate
    features = ['total_length', 'total_motifs', 'density_per_kb', 
                'coverage_percent', 'class_diversity', 'avg_score']
    
    feature_matrix = []
    for genome_name in genome_names:
        genome_stat = genome_stats[genome_name]
        feature_vector = [genome_stat[f] for f in features]
        feature_matrix.append(feature_vector)
    
    feature_matrix = np.array(feature_matrix)
    
    # Calculate correlation matrix
    if len(genome_names) > 1 and feature_matrix.shape[0] > 1:
        correlations = {}
        for i, feat1 in enumerate(features):
            for j, feat2 in enumerate(features):
                if i < j:  # Only upper triangle
                    corr, pval = stats.pearsonr(feature_matrix[:, i], feature_matrix[:, j])
                    correlations[f"{feat1}_vs_{feat2}"] = {
                        'correlation': corr,
                        'p_value': pval,
                        'significant': pval < 0.05
                    }
    else:
        correlations = {}
    
    return correlations


def generate_comparative_summary(genome_stats: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate summary statistics across all genomes.
    
    Args:
        genome_stats: Statistics for each genome
        
    Returns:
        Summary statistics
    """
    # Aggregate statistics
    total_length = sum(s['total_length'] for s in genome_stats.values())
    total_motifs = sum(s['total_motifs'] for s in genome_stats.values())
    
    avg_density = np.mean([s['density_per_kb'] for s in genome_stats.values()])
    avg_coverage = np.mean([s['coverage_percent'] for s in genome_stats.values()])
    avg_diversity = np.mean([s['class_diversity'] for s in genome_stats.values()])
    
    # Find extremes
    max_density_genome = max(genome_stats.items(), key=lambda x: x[1]['density_per_kb'])
    min_density_genome = min(genome_stats.items(), key=lambda x: x[1]['density_per_kb'])
    
    max_diversity_genome = max(genome_stats.items(), key=lambda x: x[1]['class_diversity'])
    min_diversity_genome = min(genome_stats.items(), key=lambda x: x[1]['class_diversity'])
    
    return {
        'total_genomes': len(genome_stats),
        'total_sequence_length': total_length,
        'total_motifs_detected': total_motifs,
        'average_density_per_kb': avg_density,
        'average_coverage_percent': avg_coverage,
        'average_class_diversity': avg_diversity,
        'highest_density_genome': {
            'name': max_density_genome[0],
            'density': max_density_genome[1]['density_per_kb']
        },
        'lowest_density_genome': {
            'name': min_density_genome[0],
            'density': min_density_genome[1]['density_per_kb']
        },
        'highest_diversity_genome': {
            'name': max_diversity_genome[0],
            'diversity': max_diversity_genome[1]['class_diversity']
        },
        'lowest_diversity_genome': {
            'name': min_diversity_genome[0],
            'diversity': min_diversity_genome[1]['class_diversity']
        }
    }


def save_comparative_results(comparative_results: Dict[str, Any], output_dir: str):
    """
    Save comparative analysis results.
    
    Args:
        comparative_results: Comparative analysis results
        output_dir: Output directory
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Custom JSON encoder to handle numpy types
    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.integer, np.floating)):
                return float(obj)
            elif isinstance(obj, np.bool_):
                return bool(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            return super().default(obj)
    
    # Save complete results as JSON
    results_file = os.path.join(output_dir, 'comparative_analysis.json')
    with open(results_file, 'w') as f:
        json.dump(comparative_results, f, indent=2, cls=NumpyEncoder)
    
    print(f"Saved comparative analysis: {results_file}")
    
    # Generate text report
    report_file = os.path.join(output_dir, 'comparative_report.txt')
    with open(report_file, 'w') as f:
        generate_text_report(comparative_results, f)
    
    print(f"Saved text report: {report_file}")


def generate_text_report(results: Dict[str, Any], output_file):
    """Generate human-readable text report"""
    summary = results['summary']
    
    output_file.write("="*80 + "\n")
    output_file.write("COMPARATIVE GENOME ANALYSIS REPORT\n")
    output_file.write("="*80 + "\n\n")
    
    # Summary section
    output_file.write("SUMMARY STATISTICS\n")
    output_file.write("-"*80 + "\n")
    output_file.write(f"Total genomes analyzed: {summary['total_genomes']}\n")
    output_file.write(f"Total sequence length: {summary['total_sequence_length']:,} bp\n")
    output_file.write(f"Total motifs detected: {summary['total_motifs_detected']:,}\n")
    output_file.write(f"Average motif density: {summary['average_density_per_kb']:.2f} motifs/kb\n")
    output_file.write(f"Average coverage: {summary['average_coverage_percent']:.2f}%\n")
    output_file.write(f"Average class diversity: {summary['average_class_diversity']:.3f}\n\n")
    
    # Extremes
    output_file.write("NOTABLE PATTERNS\n")
    output_file.write("-"*80 + "\n")
    output_file.write(f"Highest density: {summary['highest_density_genome']['name']} "
                     f"({summary['highest_density_genome']['density']:.2f} motifs/kb)\n")
    output_file.write(f"Lowest density: {summary['lowest_density_genome']['name']} "
                     f"({summary['lowest_density_genome']['density']:.2f} motifs/kb)\n")
    output_file.write(f"Highest diversity: {summary['highest_diversity_genome']['name']} "
                     f"(H={summary['highest_diversity_genome']['diversity']:.3f})\n")
    output_file.write(f"Lowest diversity: {summary['lowest_diversity_genome']['name']} "
                     f"(H={summary['lowest_diversity_genome']['diversity']:.3f})\n\n")
    
    # Per-genome statistics
    output_file.write("PER-GENOME STATISTICS\n")
    output_file.write("-"*80 + "\n")
    for genome_name, stats in results['genome_statistics'].items():
        output_file.write(f"\n{genome_name}:\n")
        output_file.write(f"  Length: {stats['total_length']:,} bp\n")
        output_file.write(f"  Motifs: {stats['total_motifs']:,}\n")
        output_file.write(f"  Density: {stats['density_per_kb']:.2f} motifs/kb\n")
        output_file.write(f"  Coverage: {stats['coverage_percent']:.2f}%\n")
        output_file.write(f"  Classes: {stats['unique_classes']}\n")
        output_file.write(f"  Diversity: {stats['class_diversity']:.3f}\n")


def main():
    """Main function for comparative analysis"""
    parser = argparse.ArgumentParser(
        description='Comparative analysis of batch genome results'
    )
    parser.add_argument(
        '--input-dir',
        default='batch_results',
        help='Input directory with batch results (default: batch_results)'
    )
    parser.add_argument(
        '--output-dir',
        default='comparative_results',
        help='Output directory for comparative results (default: comparative_results)'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("COMPARATIVE GENOME ANALYSIS")
    print("="*80)
    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {args.output_dir}\n")
    
    # Load batch results
    print("Loading batch results...")
    try:
        results = load_batch_results(args.input_dir)
        print(f"Loaded {len(results)} genome results\n")
    except Exception as e:
        print(f"ERROR loading batch results: {str(e)}")
        return 1
    
    # Perform comparative analysis
    print("Performing comparative analysis...")
    comparative_results = compare_genomes(results)
    print("✓ Comparative analysis complete\n")
    
    # Save results
    print("Saving results...")
    save_comparative_results(comparative_results, args.output_dir)
    print(f"\n✓ All results saved to: {args.output_dir}/\n")
    
    # Print summary
    summary = comparative_results['summary']
    print("="*80)
    print("ANALYSIS SUMMARY")
    print("="*80)
    print(f"Genomes: {summary['total_genomes']}")
    print(f"Total motifs: {summary['total_motifs_detected']:,}")
    print(f"Average density: {summary['average_density_per_kb']:.2f} motifs/kb")
    print(f"Average diversity: {summary['average_class_diversity']:.3f}")
    print("="*80)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
