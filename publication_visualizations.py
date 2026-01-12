#!/usr/bin/env python3
"""
Publication Visualizations Module
=================================

Generates publication-quality visualizations for comparative genome analysis:
- Multi-genome comparison plots
- Statistical distribution figures
- Heatmaps and clustering visualizations

Usage:
    python publication_visualizations.py [--comparative-dir DIR] [--output-dir OUTPUT_DIR]
"""

import os
import sys
import json
import argparse
from pathlib import Path
from typing import Dict, Any, List
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# Publication-quality plot styling
def setup_publication_style():
    """Setup matplotlib for publication-quality figures"""
    plt.rcParams.update({
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'font.size': 10,
        'font.family': 'sans-serif',
        'axes.labelsize': 11,
        'axes.titlesize': 12,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'figure.titlesize': 14,
        'axes.linewidth': 1.2,
        'grid.linewidth': 0.8,
        'lines.linewidth': 2,
    })


def plot_density_comparison(comparative_results: Dict, output_dir: str):
    """
    Create bar plot comparing motif density across genomes.
    
    Args:
        comparative_results: Comparative analysis results
        output_dir: Output directory
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Extract data
    genome_names = []
    densities = []
    
    for genome_name, stats in comparative_results['genome_statistics'].items():
        # Shorten genome names for display
        short_name = genome_name.replace('.fna', '').replace('.fasta', '').replace('_', ' ')
        genome_names.append(short_name[:25])  # Limit to 25 chars
        densities.append(stats['density_per_kb'])
    
    # Sort by density
    sorted_data = sorted(zip(genome_names, densities), key=lambda x: x[1], reverse=True)
    genome_names, densities = zip(*sorted_data)
    
    # Create bar plot
    colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(genome_names)))
    bars = ax.barh(genome_names, densities, color=colors, edgecolor='black', linewidth=1)
    
    # Styling
    ax.set_xlabel('Motif Density (motifs/kb)', fontweight='bold')
    ax.set_ylabel('Genome', fontweight='bold')
    ax.set_title('Non-B DNA Motif Density Across Genomes', fontweight='bold', pad=15)
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'density_comparison.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved: density_comparison.png")


def plot_class_distribution_heatmap(comparative_results: Dict, output_dir: str):
    """
    Create heatmap showing motif class distribution across genomes.
    
    Args:
        comparative_results: Comparative analysis results
        output_dir: Output directory
    """
    # Collect data
    genome_names = list(comparative_results['genome_statistics'].keys())
    
    # Get all classes
    all_classes = set()
    for stats in comparative_results['genome_statistics'].values():
        all_classes.update(stats['class_counts'].keys())
    
    # Remove Hybrid and Cluster classes for cleaner visualization
    all_classes = sorted([c for c in all_classes if c not in ['Hybrid', 'Non-B_DNA_Clusters']])
    
    # Create matrix of normalized frequencies
    matrix = []
    for genome_name in genome_names:
        stats = comparative_results['genome_statistics'][genome_name]
        total_motifs = stats['total_motifs']
        
        row = []
        for class_name in all_classes:
            count = stats['class_counts'].get(class_name, 0)
            frequency = (count / total_motifs * 100) if total_motifs > 0 else 0
            row.append(frequency)
        
        matrix.append(row)
    
    matrix = np.array(matrix)
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Shorten names for display
    short_genome_names = [gn.replace('.fna', '').replace('.fasta', '')[:20] 
                          for gn in genome_names]
    short_class_names = [cn.replace('_', ' ')[:20] for cn in all_classes]
    
    # Plot heatmap
    im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto')
    
    # Set ticks
    ax.set_xticks(np.arange(len(all_classes)))
    ax.set_yticks(np.arange(len(genome_names)))
    ax.set_xticklabels(short_class_names, rotation=45, ha='right')
    ax.set_yticklabels(short_genome_names)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Frequency (%)', rotation=270, labelpad=20, fontweight='bold')
    
    # Title
    ax.set_title('Motif Class Distribution Across Genomes', 
                 fontweight='bold', pad=15)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'class_distribution_heatmap.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved: class_distribution_heatmap.png")


def plot_diversity_vs_density(comparative_results: Dict, output_dir: str):
    """
    Create scatter plot of class diversity vs motif density.
    
    Args:
        comparative_results: Comparative analysis results
        output_dir: Output directory
    """
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Extract data
    densities = []
    diversities = []
    genome_names = []
    
    for genome_name, stats in comparative_results['genome_statistics'].items():
        densities.append(stats['density_per_kb'])
        diversities.append(stats['class_diversity'])
        short_name = genome_name.replace('.fna', '').replace('.fasta', '')[:15]
        genome_names.append(short_name)
    
    # Create scatter plot
    scatter = ax.scatter(densities, diversities, s=150, alpha=0.7, 
                        c=range(len(densities)), cmap='viridis',
                        edgecolors='black', linewidth=1.5)
    
    # Add genome labels
    for i, name in enumerate(genome_names):
        ax.annotate(name, (densities[i], diversities[i]), 
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=8, alpha=0.8)
    
    # Styling
    ax.set_xlabel('Motif Density (motifs/kb)', fontweight='bold')
    ax.set_ylabel('Class Diversity (Shannon Entropy)', fontweight='bold')
    ax.set_title('Relationship Between Motif Density and Class Diversity',
                 fontweight='bold', pad=15)
    ax.grid(alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'diversity_vs_density.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved: diversity_vs_density.png")


def plot_enrichment_patterns(comparative_results: Dict, output_dir: str):
    """
    Create visualization of enrichment/depletion patterns.
    
    Args:
        comparative_results: Comparative analysis results
        output_dir: Output directory
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Collect enrichment data
    class_names = []
    mean_frequencies = []
    std_frequencies = []
    has_enrichment = []
    
    for class_name, enrichment in comparative_results['enrichment_patterns'].items():
        if class_name not in ['Hybrid', 'Non-B_DNA_Clusters']:
            class_names.append(class_name.replace('_', ' ')[:25])
            mean_frequencies.append(enrichment['mean_frequency'] * 100)
            std_frequencies.append(enrichment['std_frequency'] * 100)
            has_enrichment.append(len(enrichment['enriched_genomes']) > 0 or 
                                 len(enrichment['depleted_genomes']) > 0)
    
    # Sort by mean frequency
    sorted_data = sorted(zip(class_names, mean_frequencies, std_frequencies, has_enrichment),
                        key=lambda x: x[1], reverse=True)
    class_names, mean_frequencies, std_frequencies, has_enrichment = zip(*sorted_data)
    
    # Create bar plot with error bars
    colors = ['#d62728' if enrich else '#1f77b4' for enrich in has_enrichment]
    
    y_pos = np.arange(len(class_names))
    bars = ax.barh(y_pos, mean_frequencies, xerr=std_frequencies,
                   color=colors, alpha=0.7, edgecolor='black', linewidth=1,
                   error_kw={'linewidth': 1.5, 'elinewidth': 1.5})
    
    # Styling
    ax.set_yticks(y_pos)
    ax.set_yticklabels(class_names)
    ax.set_xlabel('Mean Frequency (%) ± SD', fontweight='bold')
    ax.set_ylabel('Motif Class', fontweight='bold')
    ax.set_title('Motif Class Frequencies with Enrichment Patterns',
                 fontweight='bold', pad=15)
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#d62728', alpha=0.7, edgecolor='black', label='Has enrichment/depletion'),
        Patch(facecolor='#1f77b4', alpha=0.7, edgecolor='black', label='No enrichment/depletion')
    ]
    ax.legend(handles=legend_elements, loc='lower right')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'enrichment_patterns.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved: enrichment_patterns.png")


def plot_genome_statistics_radar(comparative_results: Dict, output_dir: str):
    """
    Create radar plot comparing genome statistics.
    
    Args:
        comparative_results: Comparative analysis results
        output_dir: Output directory
    """
    # Select up to 6 genomes for clarity
    genome_names = list(comparative_results['genome_statistics'].keys())[:6]
    
    if len(genome_names) < 2:
        print("  ⚠ Skipping radar plot (need at least 2 genomes)")
        return
    
    # Features for radar plot (normalized 0-1)
    features = ['density_per_kb', 'coverage_percent', 'class_diversity', 'avg_score']
    feature_labels = ['Density', 'Coverage', 'Diversity', 'Avg Score']
    
    # Normalize data
    normalized_data = {}
    for feature in features:
        values = [comparative_results['genome_statistics'][gn][feature] 
                 for gn in genome_names]
        max_val = max(values)
        min_val = min(values)
        range_val = max_val - min_val if max_val != min_val else 1
        
        normalized_data[feature] = [(v - min_val) / range_val for v in values]
    
    # Create radar plot
    angles = np.linspace(0, 2 * np.pi, len(features), endpoint=False).tolist()
    angles += angles[:1]  # Complete the circle
    
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(genome_names)))
    
    for i, genome_name in enumerate(genome_names):
        values = [normalized_data[feature][i] for feature in features]
        values += values[:1]  # Complete the circle
        
        short_name = genome_name.replace('.fna', '').replace('.fasta', '')[:20]
        ax.plot(angles, values, 'o-', linewidth=2, label=short_name, color=colors[i])
        ax.fill(angles, values, alpha=0.15, color=colors[i])
    
    # Styling
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(feature_labels, fontsize=11, fontweight='bold')
    ax.set_ylim(0, 1)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(['0.25', '0.5', '0.75', '1.0'], fontsize=9)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.set_title('Comparative Genome Statistics (Normalized)',
                 fontweight='bold', pad=20, fontsize=14)
    
    # Legend
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'genome_statistics_radar.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved: genome_statistics_radar.png")


def plot_correlation_matrix(comparative_results: Dict, output_dir: str):
    """
    Create correlation matrix heatmap for genome features.
    
    Args:
        comparative_results: Comparative analysis results
        output_dir: Output directory
    """
    if not comparative_results['correlations']:
        print("  ⚠ Skipping correlation matrix (no correlations computed)")
        return
    
    # Extract correlation data
    features = ['total_length', 'total_motifs', 'density_per_kb', 
                'coverage_percent', 'class_diversity', 'avg_score']
    
    # Build correlation matrix
    n_features = len(features)
    corr_matrix = np.eye(n_features)
    
    for i, feat1 in enumerate(features):
        for j, feat2 in enumerate(features):
            if i < j:
                key = f"{feat1}_vs_{feat2}"
                if key in comparative_results['correlations']:
                    corr = comparative_results['correlations'][key]['correlation']
                    corr_matrix[i, j] = corr
                    corr_matrix[j, i] = corr
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    
    im = ax.imshow(corr_matrix, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    
    # Set ticks
    feature_labels = [f.replace('_', ' ').title() for f in features]
    ax.set_xticks(np.arange(len(features)))
    ax.set_yticks(np.arange(len(features)))
    ax.set_xticklabels(feature_labels, rotation=45, ha='right')
    ax.set_yticklabels(feature_labels)
    
    # Add correlation values
    for i in range(len(features)):
        for j in range(len(features)):
            text = ax.text(j, i, f'{corr_matrix[i, j]:.2f}',
                          ha='center', va='center', 
                          color='white' if abs(corr_matrix[i, j]) > 0.5 else 'black',
                          fontsize=9, fontweight='bold')
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Correlation Coefficient', rotation=270, labelpad=20, fontweight='bold')
    
    ax.set_title('Feature Correlation Matrix', fontweight='bold', pad=15)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'correlation_matrix.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved: correlation_matrix.png")


def generate_all_visualizations(comparative_dir: str, output_dir: str):
    """
    Generate all publication visualizations.
    
    Args:
        comparative_dir: Directory with comparative results
        output_dir: Output directory for visualizations
    """
    # Load comparative results
    comparative_file = os.path.join(comparative_dir, 'comparative_analysis.json')
    
    with open(comparative_file, 'r') as f:
        comparative_results = json.load(f)
    
    # Setup publication style
    setup_publication_style()
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate visualizations
    print("\nGenerating visualizations:")
    
    plot_density_comparison(comparative_results, output_dir)
    plot_class_distribution_heatmap(comparative_results, output_dir)
    plot_diversity_vs_density(comparative_results, output_dir)
    plot_enrichment_patterns(comparative_results, output_dir)
    plot_genome_statistics_radar(comparative_results, output_dir)
    plot_correlation_matrix(comparative_results, output_dir)


def main():
    """Main function for visualization generation"""
    parser = argparse.ArgumentParser(
        description='Generate publication-quality visualizations'
    )
    parser.add_argument(
        '--comparative-dir',
        default='comparative_results',
        help='Comparative results directory (default: comparative_results)'
    )
    parser.add_argument(
        '--output-dir',
        default='publication_figures',
        help='Output directory for figures (default: publication_figures)'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("PUBLICATION VISUALIZATIONS GENERATOR")
    print("="*80)
    print(f"Input directory: {args.comparative_dir}")
    print(f"Output directory: {args.output_dir}")
    
    try:
        generate_all_visualizations(args.comparative_dir, args.output_dir)
        
        print("\n" + "="*80)
        print("VISUALIZATION GENERATION COMPLETE")
        print("="*80)
        print(f"All figures saved to: {args.output_dir}/")
        print("\nGenerated figures:")
        print("  • density_comparison.png - Motif density across genomes")
        print("  • class_distribution_heatmap.png - Class distribution heatmap")
        print("  • diversity_vs_density.png - Diversity vs density scatter plot")
        print("  • enrichment_patterns.png - Enrichment/depletion patterns")
        print("  • genome_statistics_radar.png - Multi-feature radar plot")
        print("  • correlation_matrix.png - Feature correlation heatmap")
        print("="*80)
        
        return 0
        
    except Exception as e:
        print(f"\nERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
