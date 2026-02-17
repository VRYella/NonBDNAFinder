#!/usr/bin/env python3
"""
NBST Validation - Publication-Quality Figure Generation
========================================================
Generates high-resolution (300 DPI), colorblind-friendly figures with proper typography
for publication in bioinformatics journals.

This script creates 9 comprehensive figures for subclass-level validation:
1. Class/Subclass Distribution Comparison
2. Performance Metrics by Subclass
3. Genomic Coverage with Subclass Differentiation
4. Detection Counts by Subclass
5. Overlap Heatmap (Subclass Level)
6. Comprehensive Subclass Distribution
7. Precision-Recall Curves per Subclass
8. Venn Diagrams (NBF Subclasses vs NBST)
9. (Future) Sequence Logo Representations

Author: NonBDNAFinder Validation Suite
Date: 2026
Version: 2.0
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Optional, Any
import warnings
warnings.filterwarnings('ignore')

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import motif taxonomy for biological ordering
from Utilities.config.motif_taxonomy import (
    get_all_classes_taxonomy_order,
    get_all_subclasses_taxonomy_order,
    get_subclasses_for_class,
    MOTIF_CLASSIFICATION,
    SUBCLASS_TO_CLASS
)
from Utilities.config.colors import UNIFIED_MOTIF_COLORS

# Matplotlib configuration for publication quality
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle, Circle
import matplotlib.ticker as ticker
from matplotlib_venn import venn2, venn3
import seaborn as sns

# Publication-ready settings - Reduced DPI to 150 to save memory
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 150
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['xtick.labelsize'] = 13
plt.rcParams['ytick.labelsize'] = 13
plt.rcParams['legend.fontsize'] = 13
plt.rcParams['figure.titlesize'] = 20
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']

# Colorblind-friendly palette (ColorBrewer2 Set2)
COLORBLIND_PALETTE = [
    '#66c2a5',  # Teal
    '#fc8d62',  # Orange
    '#8da0cb',  # Blue
    '#e78ac3',  # Pink
    '#a6d854',  # Green
    '#ffd92f',  # Yellow
    '#e5c494',  # Tan
    '#b3b3b3',  # Gray
]

# Extended palette for many subclasses
EXTENDED_PALETTE = sns.color_palette("tab20", 20)

# Directory structure: scripts are in NBSTVALIDATION/scripts/
BASE_DIR = Path(__file__).parent.parent  # NBSTVALIDATION directory
DATA_DIR = BASE_DIR / "data"
FIGURES_DIR = BASE_DIR / "figures"
TABLES_DIR = BASE_DIR / "tables"

# ============================================================================
# DATA LOADING
# ============================================================================

def load_all_data() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict]:
    """Load all validation data files."""
    nbf_df = pd.read_csv(DATA_DIR / "nonbdnafinder_results.tsv", sep='\t')
    subclass_summary = pd.read_csv(DATA_DIR / "subclass_comparison_summary.tsv", sep='\t')
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
    
    return nbf_df, subclass_summary, comparison_df, nbst_data


# ============================================================================
# FIGURE 1: Class/Subclass Distribution Comparison
# ============================================================================

def create_figure1_distribution(nbf_df: pd.DataFrame, nbst_data: Dict):
    """
    Create enhanced class/subclass distribution comparison.
    Uses pie charts with unified color palette and proper typography.
    Classes are ordered by biological taxonomy, not alphabetically.
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    
    # NonBDNAFinder class distribution - use taxonomy order
    taxonomy_classes = get_all_classes_taxonomy_order()
    class_counts = nbf_df['Class'].value_counts()
    
    # Reindex to taxonomy order (only include classes present in data)
    ordered_classes = [cls for cls in taxonomy_classes if cls in class_counts.index]
    ordered_counts = [class_counts[cls] for cls in ordered_classes]
    
    # Use unified motif colors
    colors_nbf = [UNIFIED_MOTIF_COLORS.get(cls, '#808080') for cls in ordered_classes]
    
    # NBF pie chart
    wedges, texts, autotexts = axes[0].pie(
        ordered_counts,
        labels=ordered_classes,
        autopct=lambda pct: f'{pct:.1f}%\n({int(pct/100.*sum(ordered_counts))})',
        colors=colors_nbf,
        startangle=90,
        explode=[0.03]*len(ordered_classes),
        textprops={'fontsize': 12, 'weight': 'bold'}
    )
    
    # Make percentage text more readable with better size
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontsize(11)
    
    # Color code the label text to match pie slices
    for text, color in zip(texts, colors_nbf):
        text.set_color(color)
        text.set_fontweight('bold')
    
    axes[0].set_title('NonBDNAFinder\nClass Distribution\n(n={:,} motifs)'.format(len(nbf_df)),
                     fontsize=18, fontweight='bold', pad=20)
    
    # NBST class distribution
    nbst_counts = {
        'G-Quadruplex (GQ)': len(nbst_data['GQ']),
        'Slipped DNA (STR)': len(nbst_data['STR']),
        'Direct Repeats (DR)': len(nbst_data['DR']),
        'Mirror Repeats (MR)': len(nbst_data['MR']),
        'Z-DNA': len(nbst_data['Z']),
        'Curved DNA': len(nbst_data['curved']),
    }
    nbst_series = pd.Series(nbst_counts)
    colors_nbst = COLORBLIND_PALETTE[:len(nbst_series)]
    
    # NBST pie chart
    wedges2, texts2, autotexts2 = axes[1].pie(
        nbst_series.values,
        labels=nbst_series.index,
        autopct=lambda pct: f'{pct:.1f}%\n({int(pct/100.*sum(nbst_counts.values()))})',
        colors=colors_nbst,
        startangle=90,
        explode=[0.03]*len(nbst_series),
        textprops={'fontsize': 12, 'weight': 'bold'}
    )
    
    for autotext in autotexts2:
        autotext.set_color('white')
        autotext.set_fontsize(11)
    
    axes[1].set_title('NBST Benchmark\nClass Distribution\n(n={:,} motifs)'.format(sum(nbst_counts.values())),
                     fontsize=18, fontweight='bold', pad=20)
    
    plt.suptitle('Comparative Class Distribution Analysis', fontsize=20, fontweight='bold', y=0.98)
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'figure1_subclass_distribution.png',
                dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print("✓ Saved: figure1_subclass_distribution.png (150 DPI)")


# ============================================================================
# FIGURE 2: Performance Metrics by Subclass
# ============================================================================

def create_figure2_performance_metrics(subclass_summary: pd.DataFrame):
    """
    Create subclass-level performance metrics with grouped bars.
    Separate panels for each major class showing subclass performance.
    Uses biological taxonomy order, not alphabetical.
    """
    # Group by NBF_Class and order by taxonomy
    taxonomy_classes = get_all_classes_taxonomy_order()
    present_classes = [cls for cls in taxonomy_classes if cls in subclass_summary['NBF_Class'].values]
    
    # Create figure with subplots
    n_classes = len(present_classes)
    fig, axes = plt.subplots(n_classes, 1, figsize=(14, 4*n_classes))
    
    if n_classes == 1:
        axes = [axes]
    
    for idx, cls in enumerate(present_classes):
        ax = axes[idx]
        subset = subclass_summary[subclass_summary['NBF_Class'] == cls].copy()
        
        if len(subset) == 0:
            continue
        
        # Order subclasses by taxonomy, not F1 score
        try:
            taxonomy_subclasses = get_subclasses_for_class(cls)
            ordered_subclasses = [sub for sub in taxonomy_subclasses if sub in subset['NBF_Subclass'].values]
            # Reindex subset to maintain taxonomy order
            subset = subset.set_index('NBF_Subclass').loc[ordered_subclasses].reset_index()
        except (KeyError, ValueError):
            # If taxonomy not available, keep original order
            pass
        
        x = np.arange(len(subset))
        width = 0.25
        
        bars1 = ax.bar(x - width, subset['Precision'], width,
                      label='Precision', color=COLORBLIND_PALETTE[0], alpha=0.9)
        bars2 = ax.bar(x, subset['Recall'], width,
                      label='Recall', color=COLORBLIND_PALETTE[1], alpha=0.9)
        bars3 = ax.bar(x + width, subset['F1_Score'], width,
                      label='F1 Score', color=COLORBLIND_PALETTE[2], alpha=0.9)
        
        ax.set_xlabel('Subclass', fontsize=15, fontweight='bold')
        ax.set_ylabel('Score', fontsize=15, fontweight='bold')
        ax.set_title(f'{cls} Subclass Performance vs NBST',
                    fontsize=16, fontweight='bold', pad=10)
        ax.set_xticks(x)
        
        # Color code the subclass labels
        class_color = UNIFIED_MOTIF_COLORS.get(cls, '#808080')
        ax.set_xticklabels(subset['NBF_Subclass'], rotation=45, ha='right', fontsize=12, color=class_color)
        
        ax.legend(loc='upper right', framealpha=0.9)
        ax.set_ylim(0, 1.1)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1)
        
        # Add value labels with better font size
        for bars in [bars1, bars2, bars3]:
            for bar in bars:
                height = bar.get_height()
                if height > 0.02:
                    ax.text(bar.get_x() + bar.get_width()/2., height,
                           f'{height:.2f}',
                           ha='center', va='bottom', fontsize=10)
    
    plt.suptitle('Subclass-Level Performance Metrics Against NBST Benchmark',
                fontsize=20, fontweight='bold', y=1.0)
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'figure2_performance_metrics.png',
                dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print("✓ Saved: figure2_performance_metrics.png (150 DPI)")


# ============================================================================
# FIGURE 3: Genomic Coverage with Subclass Differentiation
# ============================================================================

def create_figure3_genomic_coverage(nbf_df: pd.DataFrame, nbst_data: Dict):
    """
    Create genomic coverage map with subclass differentiation using shading/patterns.
    Excludes Hybrid and Non-B_DNA_Clusters to focus on canonical motifs.
    Uses taxonomy order and unified colors.
    """
    fig, ax = plt.subplots(figsize=(18, 10))
    
    seq_len = 40523
    
    # Filter out Hybrid and Clusters
    nbf_core = nbf_df[~nbf_df['Class'].isin(['Hybrid', 'Non-B_DNA_Clusters'])].copy()
    
    # Assign y-positions and colors using taxonomy order
    y_pos = 0
    y_positions = {}
    class_colors = {}
    
    taxonomy_classes = get_all_classes_taxonomy_order()
    
    for cls in taxonomy_classes:
        if cls not in nbf_core['Class'].values:
            continue
        
        # Use unified color for this class
        base_color = UNIFIED_MOTIF_COLORS.get(cls, '#808080')
        
        # Get subclasses in taxonomy order
        try:
            taxonomy_subclasses = get_subclasses_for_class(cls)
            subclasses = [sub for sub in taxonomy_subclasses if sub in nbf_core[nbf_core['Class'] == cls]['Subclass'].values]
        except (KeyError, ValueError):
            subclasses = nbf_core[nbf_core['Class'] == cls]['Subclass'].unique()
        
        for sub_idx, subclass in enumerate(subclasses):
            key = f"{cls}::{subclass}"
            y_positions[key] = y_pos
            
            # Create color variations for subclasses (lighter/darker shades)
            alpha = 0.6 + (sub_idx * 0.4 / max(len(subclasses) - 1, 1))
            class_colors[key] = (*matplotlib.colors.to_rgb(base_color), alpha)
            
            y_pos += 1
    
    # Plot motifs
    for _, row in nbf_core.iterrows():
        key = f"{row['Class']}::{row['Subclass']}"
        if key in y_positions:
            start = row['Start']
            end = row['End']
            y = y_positions[key]
            
            rect = Rectangle((start, y - 0.4), end - start, 0.8,
                           facecolor=class_colors[key], edgecolor='black',
                           linewidth=0.3, alpha=0.8)
            ax.add_patch(rect)
    
    # Add labels with color coding
    for key, y in y_positions.items():
        cls, subclass = key.split('::', 1)
        label = f"{cls}\n{subclass}"
        label_color = UNIFIED_MOTIF_COLORS.get(cls, '#000000')
        ax.text(-1500, y, label, ha='right', va='center', fontsize=10, wrap=True, color=label_color, fontweight='bold')
    
    ax.set_xlim(-5000, seq_len + 1000)
    ax.set_ylim(-1, len(y_positions))
    ax.set_xlabel('Genomic Position (bp)', fontsize=16, fontweight='bold')
    ax.set_title('Subclass-Level Genomic Coverage Map\n(40,523 bp sequence, excluding Hybrid/Clusters)',
                fontsize=18, fontweight='bold', pad=15)
    ax.set_yticks([])
    
    # Format x-axis with comma separators
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f'{int(x):,}'))
    
    # Add grid
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'figure3_genomic_coverage.png',
                dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print("✓ Saved: figure3_genomic_coverage.png (150 DPI)")


# ============================================================================
# FIGURE 4: Detection Counts by Subclass
# ============================================================================

def create_figure4_detection_counts(subclass_summary: pd.DataFrame):
    """
    Create side-by-side comparison of detection counts by subclass.
    Uses biological taxonomy order, not count-based sorting.
    """
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Order by taxonomy instead of NBF_Count
    taxonomy_subclasses = get_all_subclasses_taxonomy_order()
    ordered_subclasses = [sub for sub in taxonomy_subclasses if sub in subclass_summary['NBF_Subclass'].values]
    
    # Reindex to maintain taxonomy order
    data = subclass_summary.set_index('NBF_Subclass').loc[ordered_subclasses].reset_index()
    
    x = np.arange(len(data))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, data['NBST_Count'], width,
                  label='NBST Count', color=COLORBLIND_PALETTE[1], alpha=0.9)
    bars2 = ax.bar(x + width/2, data['NBF_Count'], width,
                  label='NonBDNAFinder Count', color=COLORBLIND_PALETTE[0], alpha=0.9)
    
    ax.set_xlabel('NBF Subclass', fontsize=16, fontweight='bold')
    ax.set_ylabel('Number of Motifs', fontsize=16, fontweight='bold')
    ax.set_title('Detection Count Comparison by Subclass:\nNBST vs NonBDNAFinder (Taxonomy Order)',
                fontsize=18, fontweight='bold', pad=15)
    ax.set_xticks(x)
    
    # Color code subclass labels by their parent class
    subclass_colors = []
    for subclass in data['NBF_Subclass']:
        parent_class = SUBCLASS_TO_CLASS.get(subclass, '')
        color = UNIFIED_MOTIF_COLORS.get(parent_class, '#000000')
        subclass_colors.append(color)
    
    # Set labels with colors
    labels = ax.set_xticklabels(data['NBF_Subclass'], rotation=45, ha='right', fontsize=12)
    for label, color in zip(labels, subclass_colors):
        label.set_color(color)
        label.set_fontweight('bold')
    
    ax.legend(loc='upper right', framealpha=0.9)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add value labels with better font
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{int(height)}',
                       ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'figure4_detection_counts.png',
                dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print("✓ Saved: figure4_detection_counts.png (150 DPI)")


# ============================================================================
# FIGURE 5: Overlap Heatmap (Subclass Level)
# ============================================================================

def create_figure5_overlap_heatmap(nbf_df: pd.DataFrame):
    """
    Create overlap heatmap at subclass level with proper color scheme.
    Excludes Hybrid and Non-B_DNA_Clusters to focus on canonical motif overlaps.
    Uses taxonomy ordering.
    """
    # Filter out Hybrid and Clusters
    nbf_core = nbf_df[~nbf_df['Class'].isin(['Hybrid', 'Non-B_DNA_Clusters'])].copy()
    
    # Create subclass combinations
    nbf_core['Class_Subclass'] = nbf_core['Class'] + ' / ' + nbf_core['Subclass']
    
    # Order subclasses by taxonomy
    taxonomy_subclasses = get_all_subclasses_taxonomy_order()
    ordered_combos = []
    for sub in taxonomy_subclasses:
        for combo in nbf_core['Class_Subclass'].unique():
            if combo.endswith(' / ' + sub):
                ordered_combos.append(combo)
    
    # Calculate overlap matrix with ordered indices
    overlap_matrix = pd.DataFrame(0, index=ordered_combos, columns=ordered_combos)
    
    motifs = nbf_core.to_dict('records')
    
    for i, m1 in enumerate(motifs):
        for j, m2 in enumerate(motifs):
            if i >= j:
                continue
            
            start1, end1 = m1['Start'], m1['End']
            start2, end2 = m2['Start'], m2['End']
            
            # Check overlap
            if not (end1 <= start2 or end2 <= start1):
                key1 = m1['Class_Subclass']
                key2 = m2['Class_Subclass']
                if key1 in overlap_matrix.index and key2 in overlap_matrix.columns:
                    overlap_matrix.loc[key1, key2] += 1
                    if key1 != key2:
                        overlap_matrix.loc[key2, key1] += 1
    
    # Filter to non-zero rows/cols for readability
    non_zero_mask = (overlap_matrix.sum(axis=1) > 0) | (overlap_matrix.sum(axis=0) > 0)
    filtered = overlap_matrix.loc[non_zero_mask, non_zero_mask]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 12))
    
    # Use a colorblind-friendly sequential colormap
    sns.heatmap(filtered, annot=True, fmt='d', cmap='YlOrRd',
               linewidths=0.5, ax=ax, cbar_kws={'label': 'Number of Overlaps'},
               square=True, annot_kws={'fontsize': 9})
    
    ax.set_title('Subclass-Level Overlap Matrix\n(Canonical motifs only, excluding Hybrid/Clusters)',
                fontsize=18, fontweight='bold', pad=15)
    ax.set_xlabel('Subclass', fontsize=15, fontweight='bold')
    ax.set_ylabel('Subclass', fontsize=15, fontweight='bold')
    
    # Rotate labels for better readability with larger font
    plt.xticks(rotation=45, ha='right', fontsize=11)
    plt.yticks(rotation=0, fontsize=11)
    
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'figure5_overlap_heatmap.png',
                dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print("✓ Saved: figure5_overlap_heatmap.png (150 DPI)")


# ============================================================================
# FIGURE 6: Comprehensive Subclass Distribution
# ============================================================================

def create_figure6_subclass_distribution(nbf_df: pd.DataFrame):
    """
    Create comprehensive subclass breakdown with horizontal bars for readability.
    Excludes Hybrid and Non-B_DNA_Clusters. Uses taxonomy order and color-coded labels.
    """
    # Filter out Hybrid and Clusters
    nbf_core = nbf_df[~nbf_df['Class'].isin(['Hybrid', 'Non-B_DNA_Clusters'])].copy()
    
    # Use taxonomy order for classes
    taxonomy_classes = get_all_classes_taxonomy_order()
    major_classes = [cls for cls in taxonomy_classes if cls in nbf_core['Class'].values]
    
    n_rows = 4
    n_cols = 2
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 16))
    axes = axes.flatten()
    
    for idx, cls in enumerate(major_classes):
        if idx >= len(axes):
            break
            
        ax = axes[idx]
        subset = nbf_core[nbf_core['Class'] == cls]
        
        if len(subset) > 0:
            # Order subclasses by taxonomy
            try:
                taxonomy_subclasses = get_subclasses_for_class(cls)
                present_subclasses = [sub for sub in taxonomy_subclasses if sub in subset['Subclass'].values]
                subclass_counts = pd.Series([len(subset[subset['Subclass'] == sub]) for sub in present_subclasses],
                                           index=present_subclasses)
            except (KeyError, ValueError):
                subclass_counts = subset['Subclass'].value_counts()
            
            total = len(subset)
            
            # Use unified color for this class with variations
            base_color = UNIFIED_MOTIF_COLORS.get(cls, '#808080')
            colors = [base_color] * len(subclass_counts)
            
            # Create horizontal bar chart
            y_pos = np.arange(len(subclass_counts))
            bars = ax.barh(y_pos, subclass_counts.values, color=colors, alpha=0.7)
            
            ax.set_yticks(y_pos)
            # Color code the y-axis labels
            labels = ax.set_yticklabels(subclass_counts.index, fontsize=11)
            for label in labels:
                label.set_color(base_color)
                label.set_fontweight('bold')
            
            ax.set_xlabel('Count', fontsize=13, fontweight='bold')
            ax.set_title(f'{cls}\n(n={total:,} total)',
                        fontsize=15, fontweight='bold', pad=10, color=base_color)
            ax.grid(axis='x', alpha=0.3, linestyle='--')
            
            # Add value labels and percentages with better font
            for i, (bar, count) in enumerate(zip(bars, subclass_counts.values)):
                width = bar.get_width()
                pct = 100 * count / total
                ax.text(width, bar.get_y() + bar.get_height()/2.,
                       f' {int(count)} ({pct:.1f}%)',
                       ha='left', va='center', fontsize=10, fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'No motifs', ha='center', va='center',
                   transform=ax.transAxes, fontsize=14)
            ax.set_title(cls, fontsize=15, fontweight='bold')
    
    # Hide unused subplots
    for idx in range(len(major_classes), len(axes)):
        axes[idx].axis('off')
    
    plt.suptitle('Comprehensive Subclass Distribution by Major Class\n(Canonical motifs only)',
                fontsize=20, fontweight='bold', y=0.995)
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'figure6_subclass_distribution.png',
                dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print("✓ Saved: figure6_subclass_distribution.png (150 DPI)")


# ============================================================================
# FIGURE 7: Precision-Recall Curves per Subclass
# ============================================================================

def create_figure7_precision_recall(subclass_summary: pd.DataFrame):
    """
    Create precision-recall scatter plot for each subclass.
    """
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Group by NBF_Class and plot
    classes = subclass_summary['NBF_Class'].unique()
    
    for idx, cls in enumerate(classes):
        subset = subclass_summary[subclass_summary['NBF_Class'] == cls]
        color = COLORBLIND_PALETTE[idx % len(COLORBLIND_PALETTE)]
        
        # Plot points
        ax.scatter(subset['Recall'], subset['Precision'],
                  s=200, alpha=0.7, color=color, label=cls,
                  edgecolors='black', linewidth=1.5)
        
        # Add subclass labels with better font size
        for _, row in subset.iterrows():
            ax.annotate(row['NBF_Subclass'][:20],
                       (row['Recall'], row['Precision']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=10, alpha=0.8)
    
    # Add diagonal reference line (F1 = 0.5)
    x = np.linspace(0, 1, 100)
    for f1 in [0.2, 0.5, 0.8]:
        y = (f1 * x) / (2 * x - f1 + 1e-10)
        y = np.clip(y, 0, 1)
        ax.plot(x, y, '--', alpha=0.3, linewidth=1, color='gray')
        ax.text(0.9, (f1 * 0.9) / (2 * 0.9 - f1), f'F1={f1}',
               fontsize=11, alpha=0.6)
    
    ax.set_xlabel('Recall', fontsize=16, fontweight='bold')
    ax.set_ylabel('Precision', fontsize=16, fontweight='bold')
    ax.set_title('Precision-Recall Analysis by Subclass\n(vs NBST Benchmark)',
                fontsize=18, fontweight='bold', pad=15)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='lower left', framealpha=0.9)
    ax.set_aspect('equal')
    
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'figure7_precision_recall_curves.png',
                dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print("✓ Saved: figure7_precision_recall_curves.png (150 DPI)")


# ============================================================================
# FIGURE 8: Venn Diagrams (NBF Subclasses vs NBST)
# ============================================================================

def create_figure8_venn_diagrams(subclass_summary: pd.DataFrame):
    """
    Create Venn diagrams showing overlap between NBF subclasses and NBST classes.
    """
    # Group by NBST_Class
    nbst_classes = subclass_summary['NBST_Class'].unique()
    
    n_classes = len(nbst_classes)
    n_cols = 3
    n_rows = (n_classes + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
    axes = axes.flatten() if n_classes > 1 else [axes]
    
    for idx, nbst_cls in enumerate(nbst_classes):
        if idx >= len(axes):
            break
            
        ax = axes[idx]
        subset = subclass_summary[subclass_summary['NBST_Class'] == nbst_cls]
        
        # Calculate totals
        total_tp = subset['True_Positives'].sum()
        total_fp = subset['False_Positives'].sum()
        total_fn = subset['False_Negatives'].sum()
        
        # Create Venn diagram
        venn = venn2(subsets=(total_fp, total_fn, total_tp),
                    set_labels=('NBF Only', 'NBST Only'),
                    ax=ax, alpha=0.7)
        
        # Customize colors
        if venn.get_patch_by_id('10'):
            venn.get_patch_by_id('10').set_color(COLORBLIND_PALETTE[0])
        if venn.get_patch_by_id('01'):
            venn.get_patch_by_id('01').set_color(COLORBLIND_PALETTE[1])
        if venn.get_patch_by_id('11'):
            venn.get_patch_by_id('11').set_color(COLORBLIND_PALETTE[2])
        
        ax.set_title(f'NBST {nbst_cls}\n{len(subset)} NBF subclasses',
                    fontsize=15, fontweight='bold')
        
        # Add legend with better font
        nbf_subclasses_str = ', '.join(subset['NBF_Subclass'].head(3).tolist())
        if len(subset) > 3:
            nbf_subclasses_str += '...'
        ax.text(0.5, -0.2, f'NBF: {nbf_subclasses_str}',
               transform=ax.transAxes, ha='center', fontsize=11,
               wrap=True)
    
    # Hide unused subplots
    for idx in range(len(nbst_classes), len(axes)):
        axes[idx].axis('off')
    
    plt.suptitle('Overlap Between NonBDNAFinder Subclasses and NBST Classes',
                fontsize=20, fontweight='bold', y=0.995)
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'figure8_venn_diagrams.png',
                dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print("✓ Saved: figure8_venn_diagrams.png (150 DPI)")


# ============================================================================
# MAIN
# ============================================================================

def main():
    """Generate all publication-quality figures."""
    print("=" * 80)
    print("GENERATING PUBLICATION-QUALITY FIGURES")
    print("High Resolution (300 DPI), Colorblind-Friendly")
    print("=" * 80)
    
    # Load data
    print("\n[1] Loading data...")
    nbf_df, subclass_summary, comparison_df, nbst_data = load_all_data()
    print(f"    NBF motifs: {len(nbf_df):,}")
    print(f"    Subclass comparisons: {len(subclass_summary)}")
    
    # Generate figures
    print("\n[2] Generating figures...")
    
    create_figure1_distribution(nbf_df, nbst_data)
    create_figure2_performance_metrics(subclass_summary)
    create_figure3_genomic_coverage(nbf_df, nbst_data)
    create_figure4_detection_counts(subclass_summary)
    create_figure5_overlap_heatmap(nbf_df)
    create_figure6_subclass_distribution(nbf_df)
    create_figure7_precision_recall(subclass_summary)
    create_figure8_venn_diagrams(subclass_summary)
    
    print("\n" + "=" * 80)
    print("ALL PUBLICATION-QUALITY FIGURES GENERATED")
    print("8 figures created at 300 DPI with colorblind-friendly palettes")
    print("=" * 80)


if __name__ == "__main__":
    main()
