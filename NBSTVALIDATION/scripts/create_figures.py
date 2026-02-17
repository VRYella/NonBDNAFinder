#!/usr/bin/env python3
"""
NBST Validation - Comprehensive Analysis & Visualizations
===========================================================
Generates publication-quality figures and tables for NonBDNAFinder vs NBST comparison.

Author: Automated Analysis
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Matplotlib configuration
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import seaborn as sns

# Style settings
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Directory structure: scripts are in NBSTVALIDATION/scripts/
BASE_DIR = Path(__file__).parent.parent  # NBSTVALIDATION directory
DATA_DIR = BASE_DIR / "data"
FIGURES_DIR = BASE_DIR / "figures"
TABLES_DIR = BASE_DIR / "tables"

# ============================================================================
# LOAD DATA
# ============================================================================

def load_all_data():
    """Load all data files."""
    nbf_df = pd.read_csv(DATA_DIR / "nonbdnafinder_results.tsv", sep='\t')
    comparison_df = pd.read_csv(DATA_DIR / "comparison_summary.tsv", sep='\t')
    
    # Load NBST files
    nbst_data = {
        "curved": pd.read_csv(DATA_DIR / "693fc40d26a53_curved.tsv", sep='\t'),
        "GQ": pd.read_csv(DATA_DIR / "693fc40d26a53_GQ.tsv", sep='\t'),
        "Z": pd.read_csv(DATA_DIR / "693fc40d26a53_Z.tsv", sep='\t'),
        "STR": pd.read_csv(DATA_DIR / "693fc40d26a53_Slipped_STR.tsv", sep='\t'),
        "DR": pd.read_csv(DATA_DIR / "693fc40d26a53_slipped_DR.tsv", sep='\t'),
        "MR": pd.read_csv(DATA_DIR / "693fc40d26a53_MR.tsv", sep='\t'),
    }
    
    return nbf_df, comparison_df, nbst_data


# ============================================================================
# FIGURE 1: Class Distribution Comparison
# ============================================================================

def create_class_distribution_figure(nbf_df, nbst_data):
    """Create side-by-side class distribution comparison."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # NonBDNAFinder distribution (exclude Hybrid and Cluster for fair comparison)
    nbf_core = nbf_df[~nbf_df['Class'].isin(['Hybrid', 'Non-B_DNA_Clusters'])]
    nbf_counts = nbf_core['Class'].value_counts()
    
    # Colors for classes
    colors = plt.cm.Set3(np.linspace(0, 1, len(nbf_counts)))
    
    # NBF pie chart
    axes[0].pie(nbf_counts.values, labels=nbf_counts.index, autopct='%1.1f%%',
                colors=colors, startangle=90, explode=[0.02]*len(nbf_counts))
    axes[0].set_title('NonBDNAFinder\n(n={} core motifs)'.format(len(nbf_core)), fontsize=12, fontweight='bold')
    
    # NBST counts
    nbst_counts = {
        'G-Quadruplex': len(nbst_data['GQ']),
        'Slipped_DNA': len(nbst_data['STR']) + len(nbst_data['DR']),
        'Z-DNA': len(nbst_data['Z']),
        'Curved_DNA': len(nbst_data['curved']),
        'Mirror_Repeat': len(nbst_data['MR']),
    }
    nbst_series = pd.Series(nbst_counts)
    nbst_colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854']
    
    # NBST pie chart
    axes[1].pie(nbst_series.values, labels=nbst_series.index, autopct='%1.1f%%',
                colors=nbst_colors[:len(nbst_series)], startangle=90, explode=[0.02]*len(nbst_series))
    axes[1].set_title('NBST Benchmark\n(n={} motifs)'.format(sum(nbst_counts.values())), fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'figure1_class_distribution.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: figure1_class_distribution.png")


# ============================================================================
# FIGURE 2: Performance Metrics Bar Chart
# ============================================================================

def create_performance_metrics_figure(comparison_df):
    """Create bar chart of precision, recall, F1 scores."""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    metrics = comparison_df[['nbst_name', 'precision', 'recall', 'f1']].copy()
    
    x = np.arange(len(metrics))
    width = 0.25
    
    bars1 = ax.bar(x - width, metrics['precision'], width, label='Precision', color='#2ecc71', alpha=0.8)
    bars2 = ax.bar(x, metrics['recall'], width, label='Recall', color='#3498db', alpha=0.8)
    bars3 = ax.bar(x + width, metrics['f1'], width, label='F1 Score', color='#e74c3c', alpha=0.8)
    
    ax.set_xlabel('Non-B DNA Class', fontsize=12, fontweight='bold')
    ax.set_ylabel('Score', fontsize=12, fontweight='bold')
    ax.set_title('NonBDNAFinder Performance vs NBST Benchmark', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(metrics['nbst_name'], rotation=45, ha='right')
    ax.legend(loc='upper right')
    ax.set_ylim(0, 1.1)
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, label='Threshold')
    
    # Add value labels on bars
    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.annotate(f'{height:.2f}',
                           xy=(bar.get_x() + bar.get_width() / 2, height),
                           xytext=(0, 3), textcoords="offset points",
                           ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'figure2_performance_metrics.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: figure2_performance_metrics.png")


# ============================================================================
# FIGURE 3: Genomic Coverage Map
# ============================================================================

def create_genomic_coverage_figure(nbf_df, nbst_data):
    """Create genomic coverage visualization."""
    fig, ax = plt.subplots(figsize=(16, 8))
    
    seq_len = 40523  # From our analysis
    
    # Class colors
    class_colors = {
        'G-Quadruplex': '#e74c3c',
        'Curved_DNA': '#2ecc71',
        'Slipped_DNA': '#3498db',
        'Cruciform': '#9b59b6',
        'Z-DNA': '#f39c12',
        'R-Loop': '#1abc9c',
        'Triplex': '#e67e22',
        'i-Motif': '#34495e',
        'A-philic_DNA': '#95a5a6',
        'Hybrid': '#7f8c8d',
        'Non-B_DNA_Clusters': '#bdc3c7',
    }
    
    # Y positions for each track
    y_positions = {}
    y_pos = 0
    for cls in class_colors:
        y_positions[cls] = y_pos
        y_pos += 1
    
    # Plot NBF motifs
    for _, row in nbf_df.iterrows():
        cls = row['Class']
        if cls in y_positions:
            start = row['Start']
            end = row['End']
            rect = Rectangle((start, y_positions[cls] - 0.4), end - start, 0.8,
                            facecolor=class_colors[cls], edgecolor='none', alpha=0.7)
            ax.add_patch(rect)
    
    # Add class labels
    for cls, y in y_positions.items():
        ax.text(-500, y, cls, ha='right', va='center', fontsize=9)
    
    ax.set_xlim(-3000, seq_len + 500)
    ax.set_ylim(-1, len(y_positions))
    ax.set_xlabel('Genomic Position (bp)', fontsize=12, fontweight='bold')
    ax.set_title('NonBDNAFinder Motif Coverage Map (40,523 bp sequence)', fontsize=14, fontweight='bold')
    ax.set_yticks([])
    
    # Add legend
    legend_patches = [mpatches.Patch(color=c, label=cls, alpha=0.7) for cls, c in class_colors.items()]
    ax.legend(handles=legend_patches, loc='upper right', ncol=3, fontsize=8)
    
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'figure3_genomic_coverage.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: figure3_genomic_coverage.png")


# ============================================================================
# FIGURE 4: Detection Count Comparison
# ============================================================================

def create_detection_count_figure(comparison_df):
    """Create grouped bar chart of detection counts."""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    x = np.arange(len(comparison_df))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, comparison_df['nbst_count'], width, label='NBST', color='#3498db', alpha=0.8)
    bars2 = ax.bar(x + width/2, comparison_df['nbf_count'], width, label='NonBDNAFinder', color='#e74c3c', alpha=0.8)
    
    ax.set_xlabel('Non-B DNA Class', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Motifs', fontsize=12, fontweight='bold')
    ax.set_title('Detection Count Comparison: NBST vs NonBDNAFinder', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(comparison_df['nbst_name'], rotation=45, ha='right')
    ax.legend()
    
    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.annotate(f'{int(height)}',
                       xy=(bar.get_x() + bar.get_width() / 2, height),
                       xytext=(0, 3), textcoords="offset points",
                       ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'figure4_detection_counts.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: figure4_detection_counts.png")


# ============================================================================
# FIGURE 5: Overlap Analysis Heatmap
# ============================================================================

def create_overlap_heatmap(nbf_df):
    """Create heatmap of class overlap counts."""
    # Calculate overlaps
    classes = nbf_df['Class'].unique()
    overlap_matrix = pd.DataFrame(0, index=classes, columns=classes)
    
    motifs = nbf_df.to_dict('records')
    
    for i, m1 in enumerate(motifs):
        for j, m2 in enumerate(motifs):
            if i >= j:
                continue
            
            start1, end1 = m1['Start'], m1['End']
            start2, end2 = m2['Start'], m2['End']
            
            # Check overlap
            if not (end1 <= start2 or end2 <= start1):
                overlap_matrix.loc[m1['Class'], m2['Class']] += 1
                if m1['Class'] != m2['Class']:
                    overlap_matrix.loc[m2['Class'], m1['Class']] += 1
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Filter to non-zero rows/cols
    non_zero = overlap_matrix.loc[(overlap_matrix.sum(axis=1) > 0), (overlap_matrix.sum(axis=0) > 0)]
    
    sns.heatmap(non_zero, annot=True, fmt='d', cmap='YlOrRd', 
                linewidths=0.5, ax=ax, cbar_kws={'label': 'Number of Overlaps'})
    
    ax.set_title('Inter-Class Overlap Matrix\n(Overlaps between different Non-B DNA structures)', 
                fontsize=14, fontweight='bold')
    ax.set_xlabel('Class', fontsize=12)
    ax.set_ylabel('Class', fontsize=12)
    
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'figure5_overlap_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: figure5_overlap_heatmap.png")


# ============================================================================
# FIGURE 6: Subclass Distribution
# ============================================================================

def create_subclass_distribution(nbf_df):
    """Create subclass distribution for major classes."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    major_classes = ['G-Quadruplex', 'Curved_DNA', 'R-Loop', 'Cruciform']
    
    for ax, cls in zip(axes.flat, major_classes):
        subset = nbf_df[nbf_df['Class'] == cls]
        if len(subset) > 0:
            subclass_counts = subset['Subclass'].value_counts()
            colors = plt.cm.Set3(np.linspace(0, 1, len(subclass_counts)))
            
            bars = ax.barh(subclass_counts.index, subclass_counts.values, color=colors)
            ax.set_xlabel('Count', fontsize=10)
            ax.set_title(f'{cls}\n(n={len(subset)} motifs)', fontsize=12, fontweight='bold')
            
            # Add value labels
            for bar in bars:
                width = bar.get_width()
                ax.annotate(f'{int(width)}',
                           xy=(width, bar.get_y() + bar.get_height()/2),
                           xytext=(3, 0), textcoords="offset points",
                           ha='left', va='center', fontsize=9)
    
    plt.suptitle('Subclass Distribution by Major Class', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'figure6_subclass_distribution.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: figure6_subclass_distribution.png")


# ============================================================================
# TABLE GENERATION
# ============================================================================

def create_summary_tables(nbf_df, comparison_df, nbst_data):
    """Create summary tables for the report."""
    
    # Table 1: Overall comparison
    table1 = comparison_df[['nbst_name', 'class', 'nbst_count', 'nbf_count', 'tp', 'fp', 'fn', 
                            'precision', 'recall', 'f1']].copy()
    table1.columns = ['NBST Type', 'NBF Class', 'NBST Count', 'NBF Count', 'True Pos', 
                     'False Pos', 'False Neg', 'Precision', 'Recall', 'F1']
    table1.to_csv(TABLES_DIR / 'table1_comparison_metrics.csv', index=False)
    print("Saved: table1_comparison_metrics.csv")
    
    # Table 2: Class distribution
    class_counts = nbf_df['Class'].value_counts().reset_index()
    class_counts.columns = ['Class', 'Count']
    class_counts['Percentage'] = (class_counts['Count'] / class_counts['Count'].sum() * 100).round(2)
    class_counts.to_csv(TABLES_DIR / 'table2_class_distribution.csv', index=False)
    print("Saved: table2_class_distribution.csv")
    
    # Table 3: What NBF detects that NBST doesn't
    extra_detections = {
        'R-Loop': len(nbf_df[nbf_df['Class'] == 'R-Loop']),
        'Triplex': len(nbf_df[nbf_df['Class'] == 'Triplex']),
        'i-Motif': len(nbf_df[nbf_df['Class'] == 'i-Motif']),
        'A-philic_DNA': len(nbf_df[nbf_df['Class'] == 'A-philic_DNA']),
        'Hybrid': len(nbf_df[nbf_df['Class'] == 'Hybrid']),
        'Non-B_DNA_Clusters': len(nbf_df[nbf_df['Class'] == 'Non-B_DNA_Clusters']),
    }
    extra_df = pd.DataFrame(list(extra_detections.items()), columns=['Class', 'Count'])
    extra_df['Description'] = [
        'G-rich regions prone to R-loop formation',
        'Purine/pyrimidine tracts that can form triple helices',
        'C-rich sequences complementary to G-quadruplexes',
        'AT-rich regions with DNA bending properties',
        'Overlapping regions with multiple structure potential',
        'Dense clusters of multiple Non-B DNA types',
    ]
    extra_df.to_csv(TABLES_DIR / 'table3_nbf_unique_detections.csv', index=False)
    print("Saved: table3_nbf_unique_detections.csv")
    
    return table1, class_counts, extra_df


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 80)
    print("GENERATING VALIDATION FIGURES AND TABLES")
    print("=" * 80)
    
    # Load data
    print("\n[1] Loading data...")
    nbf_df, comparison_df, nbst_data = load_all_data()
    print(f"    NBF motifs: {len(nbf_df)}")
    print(f"    Comparison records: {len(comparison_df)}")
    
    # Create figures
    print("\n[2] Creating figures...")
    
    create_class_distribution_figure(nbf_df, nbst_data)
    create_performance_metrics_figure(comparison_df)
    create_genomic_coverage_figure(nbf_df, nbst_data)
    create_detection_count_figure(comparison_df)
    create_overlap_heatmap(nbf_df)
    create_subclass_distribution(nbf_df)
    
    # Create tables
    print("\n[3] Creating summary tables...")
    create_summary_tables(nbf_df, comparison_df, nbst_data)
    
    print("\n" + "=" * 80)
    print("ALL FIGURES AND TABLES GENERATED SUCCESSFULLY")
    print("=" * 80)


if __name__ == "__main__":
    main()
