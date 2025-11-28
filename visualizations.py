"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    NBDSCANNER VISUALIZATION SUITE                             ║
║        Comprehensive Plotting and Visualization Functions                    ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: visualizations.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.1
LICENSE: MIT

DESCRIPTION:
    Comprehensive plotting and visualization functions for Non-B DNA motif analysis.
    Includes both static and interactive visualizations with scientific styling.

VISUALIZATION TABLE:
┌──────────────────────────────────────────────────────────────────────────────┐
│ Category     │ Functions                 │ Description                       │
├──────────────┼───────────────────────────┼───────────────────────────────────┤
│ Distribution │ plot_motif_distribution   │ Class/subclass distribution plots │
│ Coverage     │ plot_coverage_map         │ Sequence coverage visualization   │
│ Statistics   │ plot_length_distribution  │ Length distributions              │
│ Hierarchy    │ plot_nested_pie_chart     │ Class-subclass hierarchy          │
│ Export       │ save_all_plots            │ Batch plot export                 │
└──────────────────────────────────────────────────────────────────────────────┘

ADVANCED VISUALIZATIONS (Publication-Quality):
┌──────────────────────────────────────────────────────────────────────────────┐
│ - plot_genome_landscape_track      │ Horizontal genomic ruler with glyphs │
│ - plot_sliding_window_heat_ribbon  │ 1D density heatmap with score        │
│ - plot_ridge_plots_length_by_class │ Stacked density ridges (joyplots)    │
│ - plot_sunburst_treemap            │ Hierarchical composition             │
│ - plot_hexbin_start_vs_score       │ 2D hexbin with marginals             │
│ - plot_upset_intersection          │ UpSet plot for motif overlaps        │
│ - plot_score_violin_beeswarm       │ Score distributions with points      │
│ - plot_cluster_hotspot_map         │ Cluster regions with annotations     │
└──────────────────────────────────────────────────────────────────────────────┘

FEATURES:
    - Publication-quality static plots
    - Interactive visualizations (Plotly)
    - Colorblind-friendly palettes
    - SVG/PNG export at 300 DPI
    - Customizable styling
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Tuple, Union
from collections import Counter, defaultdict
import warnings
warnings.filterwarnings("ignore")

# Advanced visualizations feature disabled (removed from codebase for simplification)
ADVANCED_VIZ_AVAILABLE = False

# Try to import plotly for interactive plots
try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

# =============================================================================
# STYLING & CONFIGURATION
# =============================================================================

# Scientific color palette for motif classes
MOTIF_CLASS_COLORS = {
    'Curved_DNA': '#E91E63',          # Material Pink - vibrant and scientific
    'Slipped_DNA': '#FF9800',         # Material Orange - warm and visible
    'Cruciform': '#2196F3',           # Material Blue - professional
    'R-Loop': '#4CAF50',              # Material Green - natural
    'Triplex': '#9C27B0',             # Material Purple - elegant
    'G-Quadruplex': '#FFC107',        # Material Amber - distinctive
    'i-Motif': '#FF5722',             # Material Deep Orange - bold
    'Z-DNA': '#673AB7',               # Material Deep Purple - sophisticated
    'A-philic_DNA': '#00BCD4',        # Material Cyan - fresh and modern
    'Hybrid': '#9E9E9E',              # Material Gray - neutral
    'Non-B_DNA_Clusters': '#607D8B'   # Material Blue Gray - professional
}

# Enhanced scientific styling configuration for publication-quality plots
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'font.size': 11,
    'axes.titlesize': 14,
    'axes.titleweight': 'bold',
    'axes.labelsize': 12,
    'axes.labelweight': 'semibold',
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'legend.frameon': True,
    'legend.framealpha': 0.95,
    'legend.edgecolor': '#E3F2FD',
    'figure.titlesize': 16,
    'figure.titleweight': 'bold',
    'axes.grid': True,
    'grid.alpha': 0.25,
    'grid.linestyle': '--',
    'grid.linewidth': 0.6,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.spines.left': True,
    'axes.spines.bottom': True,
    'axes.linewidth': 1.2,
    'axes.edgecolor': '#0D47A1',
    'figure.facecolor': 'white',
    'axes.facecolor': '#FAFBFC',
})

# Constants for enrichment analysis visualization
INFINITE_FOLD_ENRICHMENT_CAP = 100  # Cap for infinite fold enrichment values in plots

def set_scientific_style():
    """Apply scientific publication-ready styling"""
    sns.set_style("whitegrid")
    sns.set_palette("husl")

# =============================================================================
# DISTRIBUTION PLOTS
# =============================================================================

def plot_motif_distribution(motifs: List[Dict[str, Any]], 
                           by: str = 'Class',
                           title: Optional[str] = None,
                           figsize: Tuple[int, int] = (10, 6)) -> plt.Figure:
    """
    Plot distribution of motifs by class or subclass
    
    Args:
        motifs: List of motif dictionaries
        by: Group by 'Class' or 'Subclass'
        title: Custom plot title
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center', 
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title or 'Motif Distribution')
        return fig
    
    # Count motifs by specified grouping
    counts = Counter(m.get(by, 'Unknown') for m in motifs)
    
    # Prepare data
    categories = list(counts.keys())
    values = list(counts.values())
    
    # Get colors
    if by == 'Class':
        colors = [MOTIF_CLASS_COLORS.get(cat, '#808080') for cat in categories]
    else:
        colors = sns.color_palette("husl", len(categories))
    
    # Create plot with dynamic sizing based on number of categories
    if len(categories) > 15:
        figsize = (max(12, len(categories) * 0.6), 6)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    bars = ax.bar(range(len(categories)), values, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    
    # Customize plot
    ax.set_xlabel(f'Motif {by}', fontsize=11, fontweight='bold')
    ax.set_ylabel('Count', fontsize=11, fontweight='bold')
    ax.set_title(title or f'Distribution of Motifs by {by}', fontsize=13, fontweight='bold')
    ax.set_xticks(range(len(categories)))
    
    # Adjust label rotation and alignment based on number of categories
    if len(categories) > 10:
        ax.set_xticklabels(categories, rotation=60, ha='right', fontsize=8)
    else:
        ax.set_xticklabels(categories, rotation=45, ha='right', fontsize=9)
    
    # Add count labels on bars (only if not too crowded)
    if len(categories) <= 20:
        for bar, count in zip(bars, values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + max(values) * 0.01,
                    str(count), ha='center', va='bottom', fontweight='bold', fontsize=9)
    
    plt.tight_layout()
    return fig

def plot_class_subclass_sunburst(motifs: List[Dict[str, Any]], 
                                 title: str = "Motif Class-Subclass Distribution") -> Union[plt.Figure, Any]:
    """
    Create sunburst plot showing class-subclass hierarchy
    
    Args:
        motifs: List of motif dictionaries
        title: Plot title
        
    Returns:
        Plotly figure if available, otherwise matplotlib figure
    """
    if not motifs:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center', 
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    if not PLOTLY_AVAILABLE:
        # Fallback to matplotlib nested pie chart
        return plot_nested_pie_chart(motifs, title)
    
    # Build hierarchical data for sunburst
    class_subclass_counts = defaultdict(lambda: defaultdict(int))
    
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        subclass_name = motif.get('Subclass', 'Unknown')
        class_subclass_counts[class_name][subclass_name] += 1
    
    # Prepare data for plotly
    ids = []
    labels = []
    parents = []
    values = []
    colors = []
    
    # Add classes (inner ring)
    for class_name, subclasses in class_subclass_counts.items():
        total_class_count = sum(subclasses.values())
        ids.append(class_name)
        labels.append(f"{class_name}<br>({total_class_count})")
        parents.append("")
        values.append(total_class_count)
        colors.append(MOTIF_CLASS_COLORS.get(class_name, '#808080'))
    
    # Add subclasses (outer ring)
    for class_name, subclasses in class_subclass_counts.items():
        for subclass_name, count in subclasses.items():
            ids.append(f"{class_name}_{subclass_name}")
            labels.append(f"{subclass_name}<br>({count})")
            parents.append(class_name)
            values.append(count)
            # Lighter shade of class color for subclasses
            base_color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
            colors.append(base_color + '80')  # Add transparency
    
    fig = go.Figure(go.Sunburst(
        ids=ids,
        labels=labels,
        parents=parents,
        values=values,
        branchvalues="total",
        marker=dict(colors=colors, line=dict(color="#FFFFFF", width=2)),
        hovertemplate='<b>%{label}</b><br>Count: %{value}<extra></extra>',
    ))
    
    fig.update_layout(
        title=title,
        font_size=10,
        width=600,
        height=600
    )
    
    return fig

def plot_nested_pie_chart(motifs: List[Dict[str, Any]], 
                         title: str = "Motif Distribution") -> plt.Figure:
    """
    Create nested donut chart with improved text placement to avoid overlaps
    
    Args:
        motifs: List of motif dictionaries
        title: Plot title
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    # Count by class and subclass
    class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    class_subclass_counts = defaultdict(lambda: defaultdict(int))
    
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        subclass_name = motif.get('Subclass', 'Unknown')
        class_subclass_counts[class_name][subclass_name] += 1
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Inner donut (classes)
    class_names = list(class_counts.keys())
    class_values = list(class_counts.values())
    class_colors = [MOTIF_CLASS_COLORS.get(name, '#808080') for name in class_names]
    
    # Create inner donut with better spacing
    wedges1, texts1, autotexts1 = ax.pie(
        class_values, 
        labels=class_names,
        colors=class_colors,
        radius=0.65,
        autopct=lambda pct: f'{pct:.1f}%' if pct > 5 else '',  # Only show % if > 5%
        pctdistance=0.80,
        startangle=90,
        wedgeprops=dict(width=0.35, edgecolor='white', linewidth=2)  # Donut style
    )
    
    # Outer donut (subclasses) with improved text placement
    all_subclass_counts = []
    all_subclass_colors = []
    all_subclass_labels = []
    
    for class_name in class_names:
        subclass_dict = class_subclass_counts[class_name]
        base_color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        
        for subclass_name, count in subclass_dict.items():
            all_subclass_counts.append(count)
            # Truncate long subclass names to avoid overlap
            label = subclass_name if len(subclass_name) <= 15 else subclass_name[:12] + '...'
            all_subclass_labels.append(label)
            # Create lighter shades for subclasses
            all_subclass_colors.append(base_color)
    
    # Only show subclass labels if we have room (< 25 subclasses)
    if len(all_subclass_labels) > 25:
        all_subclass_labels = ['' for _ in all_subclass_labels]
    
    wedges2, texts2 = ax.pie(
        all_subclass_counts,
        labels=all_subclass_labels,
        colors=all_subclass_colors,
        radius=1.0,
        labeldistance=1.15,
        startangle=90,
        wedgeprops=dict(width=0.35, edgecolor='white', linewidth=1.5),  # Donut style
        textprops={'fontsize': 7}
    )
    
    ax.set_title(title, fontsize=16, pad=20, fontweight='bold')
    
    # Improve text positioning to avoid overlap
    for text in texts1:
        text.set_fontsize(9)
        text.set_fontweight('bold')
    
    for text in texts2:
        text.set_fontsize(7)
    
    # Style percentage labels
    for autotext in autotexts1:
        autotext.set_fontsize(9)
        autotext.set_fontweight('bold')
        autotext.set_color('white')
    
    return fig

# =============================================================================
# COVERAGE & POSITIONAL PLOTS
# =============================================================================

def plot_coverage_map(motifs: List[Dict[str, Any]], 
                     sequence_length: int,
                     title: Optional[str] = None,
                     figsize: Tuple[int, int] = (12, 8)) -> plt.Figure:
    """
    Plot motif coverage map showing positions along sequence
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of the analyzed sequence
        title: Custom plot title
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center', 
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title or 'Motif Coverage Map')
        return fig
    
    # Group motifs by class
    class_motifs = defaultdict(list)
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        class_motifs[class_name].append(motif)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    y_pos = 0
    class_positions = {}
    
    for class_name, class_motif_list in class_motifs.items():
        class_positions[class_name] = y_pos
        color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        
        for motif in class_motif_list:
            start = motif.get('Start', 0) - 1  # Convert to 0-based
            end = motif.get('End', start + 1)
            length = end - start
            
            # Draw motif as rectangle
            rect = patches.Rectangle(
                (start, y_pos - 0.4), length, 0.8,
                facecolor=color, edgecolor='black', alpha=0.7, linewidth=0.5
            )
            ax.add_patch(rect)
            
            # Add subclass label if space allows
            if length > sequence_length * 0.02:  # Only label if wide enough
                subclass = motif.get('Subclass', '')
                ax.text(start + length/2, y_pos, subclass, 
                       ha='center', va='center', fontsize=8, 
                       rotation=0 if length > sequence_length * 0.1 else 90)
        
        y_pos += 1
    
    # Customize plot
    ax.set_xlim(0, sequence_length)
    ax.set_ylim(-0.5, len(class_motifs) - 0.5)
    ax.set_xlabel('Sequence Position (bp)')
    ax.set_ylabel('Motif Class')
    ax.set_title(title or f'Motif Coverage Map ({sequence_length} bp)')
    
    # Set y-axis labels
    ax.set_yticks(list(class_positions.values()))
    ax.set_yticklabels(list(class_positions.keys()))
    
    # Add sequence ruler
    ruler_ticks = np.arange(0, sequence_length + 1, max(1, sequence_length // 10))
    ax.set_xticks(ruler_ticks)
    ax.set_xticklabels([f'{int(x):,}' for x in ruler_ticks])
    
    plt.tight_layout()
    return fig

def plot_density_heatmap(motifs: List[Dict[str, Any]], 
                        sequence_length: int,
                        window_size: int = 1000,
                        title: Optional[str] = None,
                        figsize: Tuple[int, int] = (12, 6)) -> plt.Figure:
    """
    Plot motif density heatmap along sequence
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of the analyzed sequence
        window_size: Window size for density calculation
        title: Custom plot title
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center', 
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title or 'Motif Density Heatmap')
        return fig
    
    # Calculate windows
    num_windows = max(1, sequence_length // window_size)
    windows = np.linspace(0, sequence_length, num_windows + 1)
    
    # Get unique classes
    classes = sorted(set(m.get('Class', 'Unknown') for m in motifs))
    
    # Calculate density matrix
    density_matrix = np.zeros((len(classes), num_windows))
    
    for i, class_name in enumerate(classes):
        class_motifs = [m for m in motifs if m.get('Class') == class_name]
        
        for j in range(num_windows):
            window_start = windows[j]
            window_end = windows[j + 1]
            
            # Count motifs in window
            count = 0
            for motif in class_motifs:
                motif_start = motif.get('Start', 0)
                motif_end = motif.get('End', 0)
                
                # Check if motif overlaps with window
                if not (motif_end <= window_start or motif_start >= window_end):
                    count += 1
            
            density_matrix[i, j] = count
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=figsize)
    
    im = ax.imshow(density_matrix, cmap='YlOrRd', aspect='auto', interpolation='nearest')
    
    # Customize plot
    ax.set_xlabel(f'Sequence Position (windows of {window_size:,} bp)')
    ax.set_ylabel('Motif Class')
    ax.set_title(title or f'Motif Density Heatmap (Window size: {window_size:,} bp)')
    
    # Set ticks and labels
    ax.set_yticks(range(len(classes)))
    ax.set_yticklabels(classes)
    
    x_ticks = np.arange(0, num_windows, max(1, num_windows // 10))
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f'{int(windows[i]):,}' for i in x_ticks])
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Motif Count')
    
    plt.tight_layout()
    return fig

# =============================================================================
# STATISTICAL PLOTS
# =============================================================================

def plot_score_distribution(motifs: List[Dict[str, Any]], 
                           by_class: bool = True,
                           title: Optional[str] = None,
                           figsize: Tuple[int, int] = (10, 6)) -> plt.Figure:
    """
    Plot distribution of motif scores
    
    Args:
        motifs: List of motif dictionaries
        by_class: Whether to separate by motif class
        title: Custom plot title
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center', 
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title or 'Score Distribution')
        return fig
    
    # Extract scores (use Score field, fallback to Normalized_Score for backward compatibility)
    scores_data = []
    for motif in motifs:
        score = motif.get('Score', motif.get('Normalized_Score'))
        if isinstance(score, (int, float)):
            if by_class:
                scores_data.append({
                    'Score': score,
                    'Class': motif.get('Class', 'Unknown')
                })
            else:
                scores_data.append(score)
    
    if not scores_data:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No score data available', ha='center', va='center', 
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title or 'Score Distribution')
        return fig
    
    fig, ax = plt.subplots(figsize=figsize)
    
    if by_class and isinstance(scores_data[0], dict):
        # Create DataFrame for seaborn
        df = pd.DataFrame(scores_data)
        
        # Box plot by class
        sns.boxplot(data=df, x='Class', y='Score', ax=ax, palette='Set2')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        ax.set_ylabel('Normalized Score')
        ax.set_xlabel('Motif Class')
    else:
        # Simple histogram
        ax.hist(scores_data, bins=20, alpha=0.7, edgecolor='black')
        ax.set_xlabel('Normalized Score')
        ax.set_ylabel('Frequency')
    
    ax.set_title(title or 'Motif Score Distribution (Normalized)')
    plt.tight_layout()
    return fig

def plot_length_distribution(motifs: List[Dict[str, Any]], 
                           by_class: bool = True,
                           title: Optional[str] = None,
                           figsize: Tuple[int, int] = (10, 6)) -> plt.Figure:
    """
    Plot distribution of motif lengths
    
    Args:
        motifs: List of motif dictionaries
        by_class: Whether to separate by motif class
        title: Custom plot title
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center', 
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title or 'Length Distribution')
        return fig
    
    # Extract lengths
    length_data = []
    for motif in motifs:
        length = motif.get('Length')
        if isinstance(length, int) and length > 0:
            if by_class:
                length_data.append({
                    'Length': length,
                    'Class': motif.get('Class', 'Unknown')
                })
            else:
                length_data.append(length)
    
    if not length_data:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No length data available', ha='center', va='center', 
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title or 'Length Distribution')
        return fig
    
    fig, ax = plt.subplots(figsize=figsize)
    
    if by_class and isinstance(length_data[0], dict):
        # Create DataFrame for seaborn
        df = pd.DataFrame(length_data)
        
        # Violin plot by class
        sns.violinplot(data=df, x='Class', y='Length', ax=ax, palette='Set3')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        ax.set_ylabel('Motif Length (bp)')
        ax.set_xlabel('Motif Class')
    else:
        # Simple histogram
        ax.hist(length_data, bins=20, alpha=0.7, edgecolor='black')
        ax.set_xlabel('Motif Length (bp)')
        ax.set_ylabel('Frequency')
    
    ax.set_title(title or 'Motif Length Distribution')
    plt.tight_layout()
    return fig

# =============================================================================
# COMPARISON PLOTS
# =============================================================================

def plot_class_comparison(results: Dict[str, List[Dict[str, Any]]], 
                         metric: str = 'count',
                         title: Optional[str] = None,
                         figsize: Tuple[int, int] = (12, 8)) -> plt.Figure:
    """
    Compare motif classes across multiple sequences/samples
    
    Args:
        results: Dictionary of {sample_name: motifs_list}
        metric: Comparison metric ('count', 'coverage', 'density')
        title: Custom plot title
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    if not results:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No data to display', ha='center', va='center', 
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title or 'Class Comparison')
        return fig
    
    # Collect all classes across samples
    all_classes = set()
    for motifs in results.values():
        all_classes.update(m.get('Class', 'Unknown') for m in motifs)
    all_classes = sorted(all_classes)
    
    # Calculate metrics for each sample
    comparison_data = []
    
    for sample_name, motifs in results.items():
        class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
        
        for class_name in all_classes:
            count = class_counts.get(class_name, 0)
            
            if metric == 'count':
                value = count
            elif metric == 'coverage':
                # Calculate coverage percentage (simplified)
                class_motifs = [m for m in motifs if m.get('Class') == class_name]
                covered_length = sum(m.get('Length', 0) for m in class_motifs)
                # Assume 10kb sequence for percentage calculation
                value = (covered_length / 10000) * 100
            elif metric == 'density':
                # Motifs per kb (assume 10kb sequence)
                value = count / 10
            else:
                value = count
            
            comparison_data.append({
                'Sample': sample_name,
                'Class': class_name,
                'Value': value
            })
    
    # Create DataFrame and pivot for heatmap
    df = pd.DataFrame(comparison_data)
    pivot_df = df.pivot(index='Sample', columns='Class', values='Value').fillna(0)
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=figsize)
    
    sns.heatmap(pivot_df, annot=True, fmt='.1f', cmap='YlOrRd', 
                ax=ax, cbar_kws={'label': f'Motif {metric.title()}'})
    
    ax.set_title(title or f'Motif Class Comparison by {metric.title()}')
    ax.set_xlabel('Motif Class')
    ax.set_ylabel('Sample')
    
    plt.tight_layout()
    return fig

# =============================================================================
# INTERACTIVE PLOTS (PLOTLY)
# =============================================================================

def create_interactive_coverage_plot(motifs: List[Dict[str, Any]], 
                                   sequence_length: int,
                                   title: str = "Interactive Motif Coverage") -> Union[Any, plt.Figure]:
    """
    Create interactive coverage plot using Plotly
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of the analyzed sequence
        title: Plot title
        
    Returns:
        Plotly figure if available, otherwise matplotlib figure
    """
    if not PLOTLY_AVAILABLE:
        return plot_coverage_map(motifs, sequence_length, title)
    
    if not motifs:
        fig = go.Figure()
        fig.add_annotation(text="No motifs to display", 
                          xref="paper", yref="paper",
                          x=0.5, y=0.5, showarrow=False)
        fig.update_layout(title=title)
        return fig
    
    fig = go.Figure()
    
    # Group motifs by class
    class_motifs = defaultdict(list)
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        class_motifs[class_name].append(motif)
    
    y_pos = 0
    for class_name, class_motif_list in class_motifs.items():
        color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        
        x_starts = []
        x_ends = []
        y_positions = []
        hover_texts = []
        
        for motif in class_motif_list:
            start = motif.get('Start', 0)
            end = motif.get('End', start + 1)
            
            x_starts.append(start)
            x_ends.append(end)
            y_positions.extend([y_pos - 0.4, y_pos - 0.4, y_pos + 0.4, y_pos + 0.4, None])
            
            hover_text = f"{class_name}<br>{motif.get('Subclass', '')}<br>" + \
                        f"Position: {start}-{end}<br>Length: {motif.get('Length', 0)} bp<br>" + \
                        f"Score: {motif.get('Score', 'N/A')}"
            hover_texts.append(hover_text)
        
        # Add rectangles for motifs
        for i, (start, end) in enumerate(zip(x_starts, x_ends)):
            fig.add_shape(
                type="rect",
                x0=start, y0=y_pos - 0.4,
                x1=end, y1=y_pos + 0.4,
                fillcolor=color,
                opacity=0.7,
                line=dict(color="black", width=1)
            )
            
            # Add invisible scatter points for hover
            fig.add_trace(go.Scatter(
                x=[(start + end) / 2],
                y=[y_pos],
                mode='markers',
                marker=dict(size=0.1, opacity=0),
                hoverinfo='text',
                hovertext=hover_texts[i],
                showlegend=False
            ))
        
        # Add class label
        fig.add_trace(go.Scatter(
            x=[sequence_length * 1.02],
            y=[y_pos],
            mode='text',
            text=[class_name],
            textposition="middle right",
            showlegend=False,
            hoverinfo='skip'
        ))
        
        y_pos += 1
    
    fig.update_layout(
        title=title,
        xaxis_title="Sequence Position (bp)",
        yaxis_title="Motif Class",
        xaxis=dict(range=[0, sequence_length * 1.15]),
        yaxis=dict(range=[-0.5, len(class_motifs) - 0.5], showticklabels=False),
        showlegend=False,
        height=max(400, len(class_motifs) * 60)
    )
    
    return fig

# =============================================================================
# BATCH EXPORT FUNCTIONS
# =============================================================================

def save_all_plots(motifs: List[Dict[str, Any]], 
                   sequence_length: int,
                   output_dir: str = "plots",
                   file_format: str = "png",
                   dpi: int = 300) -> Dict[str, str]:
    """
    Generate and save all standard plots
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of analyzed sequence
        output_dir: Output directory for plots
        file_format: File format ('png', 'pdf', 'svg')
        dpi: Resolution for raster formats
        
    Returns:
        Dictionary of {plot_name: file_path}
    """
    import os
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    saved_files = {}
    
    # List of plots to generate (diverse, non-repetitive visualization types)
    plots_to_generate = [
        ("motif_distribution_class", lambda: plot_motif_distribution(motifs, by='Class')),
        ("coverage_map", lambda: plot_coverage_map(motifs, sequence_length)),
        ("density_heatmap", lambda: plot_density_heatmap(motifs, sequence_length)),
        ("score_distribution", lambda: plot_score_distribution(motifs, by_class=True)),
        ("length_distribution", lambda: plot_length_distribution(motifs, by_class=True)),
        ("nested_donut_chart", lambda: plot_nested_pie_chart(motifs))
    ]
    
    # Add advanced visualizations if available
    if ADVANCED_VIZ_AVAILABLE:
        plots_to_generate.extend([
            ("sunburst_hierarchy", lambda: plot_sunburst_treemap(motifs, plot_type='sunburst')),
            ("score_violin_beeswarm", lambda: plot_score_violin_beeswarm(motifs))
        ])
    
    for plot_name, plot_func in plots_to_generate:
        try:
            fig = plot_func()
            filename = f"{plot_name}.{file_format}"
            filepath = os.path.join(output_dir, filename)
            
            if hasattr(fig, 'savefig'):  # Matplotlib figure
                fig.savefig(filepath, format=file_format, dpi=dpi, bbox_inches='tight')
                plt.close(fig)
            else:  # Plotly figure
                if file_format.lower() == 'png':
                    fig.write_image(filepath)
                elif file_format.lower() == 'html':
                    fig.write_html(filepath)
            
            saved_files[plot_name] = filepath
            print(f"✓ Saved {plot_name} to {filepath}")
            
        except Exception as e:
            print(f"✗ Failed to generate {plot_name}: {e}")
    
    return saved_files


# =============================================================================
# ENHANCED SCIENTIFIC VISUALIZATION FUNCTIONS
# =============================================================================

def plot_class_analysis_comprehensive(motifs: List[Dict[str, Any]], 
                                     figsize: Tuple[int, int] = (16, 12)) -> plt.Figure:
    """
    Comprehensive class-level analysis with multiple subplots.
    Shows distribution, statistics, and comparison of all 11 Non-B DNA classes.
    Highlights which classes were detected and which were not.
    
    Args:
        motifs: List of motif dictionaries
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object with multiple subplots
    """
    set_scientific_style()
    
    # Define all 11 Non-B DNA classes
    all_classes = [
        'Curved_DNA', 'Slipped_DNA', 'Cruciform', 'R-Loop', 'Triplex',
        'G-Quadruplex', 'i-Motif', 'Z-DNA', 'A-philic_DNA', 
        'Hybrid', 'Non-B_DNA_Clusters'
    ]
    
    # Count motifs by class
    detected_classes = Counter(m.get('Class', 'Unknown') for m in motifs)
    
    # Identify detected vs not detected
    detected = [cls for cls in all_classes if detected_classes.get(cls, 0) > 0]
    not_detected = [cls for cls in all_classes if detected_classes.get(cls, 0) == 0]
    
    # Create figure with subplots
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)
    
    # 1. Main distribution bar chart
    ax1 = fig.add_subplot(gs[0, :])
    counts = [detected_classes.get(cls, 0) for cls in all_classes]
    colors = [MOTIF_CLASS_COLORS.get(cls, '#808080') for cls in all_classes]
    bars = ax1.bar(range(len(all_classes)), counts, color=colors, alpha=0.8, 
                   edgecolor='black', linewidth=1.5)
    
    ax1.set_xlabel('Non-B DNA Class', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax1.set_title('Distribution of All 11 Non-B DNA Classes', fontsize=14, fontweight='bold')
    ax1.set_xticks(range(len(all_classes)))
    ax1.set_xticklabels(all_classes, rotation=45, ha='right', fontsize=10)
    ax1.grid(axis='y', alpha=0.3)
    
    # Add count labels on bars
    for bar, count in zip(bars, counts):
        if count > 0:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{count}', ha='center', va='bottom', fontweight='bold', fontsize=9)
    
    # 2. Detected vs Not Detected pie chart
    ax2 = fig.add_subplot(gs[1, 0])
    detection_counts = [len(detected), len(not_detected)]
    detection_labels = [f'Detected\n({len(detected)} classes)', 
                       f'Not Detected\n({len(not_detected)} classes)']
    colors_pie = ['#4CAF50', '#FF5722']
    ax2.pie(detection_counts, labels=detection_labels, autopct='%1.1f%%',
            colors=colors_pie, startangle=90, textprops={'fontsize': 11, 'fontweight': 'bold'})
    ax2.set_title('Class Detection Status', fontsize=12, fontweight='bold')
    
    # 3. Statistics table for detected classes
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axis('off')
    
    if detected:
        # Calculate statistics for detected classes
        class_stats = []
        for cls in detected[:5]:  # Show top 5
            cls_motifs = [m for m in motifs if m.get('Class') == cls]
            count = len(cls_motifs)
            avg_length = np.mean([m.get('Length', 0) for m in cls_motifs]) if cls_motifs else 0
            avg_score = np.mean([m.get('Score', 0) for m in cls_motifs]) if cls_motifs else 0
            class_stats.append([cls[:15], count, f'{avg_length:.1f}', f'{avg_score:.3f}'])
        
        table = ax3.table(cellText=class_stats,
                         colLabels=['Class', 'Count', 'Avg Len', 'Avg Score'],
                         cellLoc='left', loc='center',
                         colWidths=[0.4, 0.2, 0.2, 0.2])
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 2)
        
        # Style header
        for i in range(4):
            table[(0, i)].set_facecolor('#E0E0E0')
            table[(0, i)].set_text_props(weight='bold')
        
        ax3.set_title('Top Detected Classes (Statistics)', fontsize=12, fontweight='bold', pad=20)
    else:
        ax3.text(0.5, 0.5, 'No classes detected', ha='center', va='center',
                transform=ax3.transAxes, fontsize=12)
    
    # 4. List of not detected classes
    ax4 = fig.add_subplot(gs[2, :])
    ax4.axis('off')
    
    if not_detected:
        not_detected_text = 'Classes NOT Detected:\n' + ', '.join(not_detected)
        ax4.text(0.5, 0.5, not_detected_text, ha='center', va='center',
                transform=ax4.transAxes, fontsize=11, 
                bbox=dict(boxstyle='round', facecolor='#FFEBEE', alpha=0.8),
                wrap=True)
    else:
        ax4.text(0.5, 0.5, 'All 11 Non-B DNA classes detected! ✓', 
                ha='center', va='center', transform=ax4.transAxes, fontsize=12,
                fontweight='bold', color='green')
    
    plt.suptitle(f'Comprehensive Class Analysis ({len(motifs)} total motifs)', 
                fontsize=16, fontweight='bold', y=0.98)
    
    return fig


def plot_subclass_analysis_comprehensive(motifs: List[Dict[str, Any]], 
                                        figsize: Tuple[int, int] = (18, 14)) -> plt.Figure:
    """
    Comprehensive subclass-level analysis showing all detected subclasses
    organized by their parent class.
    
    Args:
        motifs: List of motif dictionaries
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object with subclass analysis
    """
    set_scientific_style()
    
    # Group motifs by class and subclass
    class_subclass_counts = defaultdict(lambda: defaultdict(int))
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        subclass_name = motif.get('Subclass', 'Unknown')
        class_subclass_counts[class_name][subclass_name] += 1
    
    # Prepare data for visualization
    all_subclasses = []
    all_counts = []
    all_classes = []
    
    for class_name in sorted(class_subclass_counts.keys()):
        for subclass_name, count in sorted(class_subclass_counts[class_name].items()):
            all_subclasses.append(f"{class_name}:{subclass_name}")
            all_counts.append(count)
            all_classes.append(class_name)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, 
                                   gridspec_kw={'height_ratios': [2, 1]})
    
    # 1. Subclass distribution bar chart
    colors = [MOTIF_CLASS_COLORS.get(cls, '#808080') for cls in all_classes]
    x_pos = range(len(all_subclasses))
    bars = ax1.barh(x_pos, all_counts, color=colors, alpha=0.8, 
                    edgecolor='black', linewidth=0.5)
    
    ax1.set_yticks(x_pos)
    ax1.set_yticklabels(all_subclasses, fontsize=9)
    ax1.set_xlabel('Count', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Class:Subclass', fontsize=12, fontweight='bold')
    ax1.set_title('Distribution of All Detected Subclasses', fontsize=14, fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)
    
    # Add count labels on bars
    for bar, count in zip(bars, all_counts):
        width = bar.get_width()
        ax1.text(width, bar.get_y() + bar.get_height()/2.,
                f' {count}', ha='left', va='center', fontweight='bold', fontsize=8)
    
    # 2. Subclass summary by class
    ax2.axis('off')
    
    # Create summary text
    summary_lines = ['Subclass Summary by Class:\n']
    for class_name in sorted(class_subclass_counts.keys()):
        subclasses = class_subclass_counts[class_name]
        n_subclasses = len(subclasses)
        total_count = sum(subclasses.values())
        summary_lines.append(f'{class_name}: {n_subclasses} subclass(es), {total_count} motifs')
    
    summary_text = '\n'.join(summary_lines)
    ax2.text(0.1, 0.5, summary_text, ha='left', va='center',
            transform=ax2.transAxes, fontsize=10, family='monospace',
            bbox=dict(boxstyle='round', facecolor='#F5F5F5', alpha=0.8))
    
    plt.suptitle(f'Comprehensive Subclass Analysis ({len(motifs)} total motifs)', 
                fontsize=16, fontweight='bold', y=0.99)
    plt.tight_layout()
    
    return fig


def plot_score_statistics_by_class(motifs: List[Dict[str, Any]], 
                                   figsize: Tuple[int, int] = (14, 8)) -> plt.Figure:
    """
    Advanced statistical visualization of scores by class.
    Shows box plots, violin plots, and statistical annotations.
    
    Args:
        motifs: List of motif dictionaries
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object with score statistics
    """
    set_scientific_style()
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No motifs to analyze', ha='center', va='center',
               transform=ax.transAxes, fontsize=14)
        return fig
    
    # Prepare data
    df_data = []
    for motif in motifs:
        df_data.append({
            'Class': motif.get('Class', 'Unknown'),
            'Score': motif.get('Score', 0),
            'Length': motif.get('Length', 0)
        })
    df = pd.DataFrame(df_data)
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, 
                                   gridspec_kw={'height_ratios': [2, 1]})
    
    # 1. Violin plot with box plot overlay
    classes = sorted(df['Class'].unique())
    positions = range(len(classes))
    
    # Create violin plot
    parts = ax1.violinplot([df[df['Class'] == cls]['Score'].values for cls in classes],
                          positions=positions, widths=0.7, showmeans=True, showmedians=True)
    
    # Color violins by class
    for i, pc in enumerate(parts['bodies']):
        cls = classes[i]
        color = MOTIF_CLASS_COLORS.get(cls, '#808080')
        pc.set_facecolor(color)
        pc.set_alpha(0.6)
    
    # Overlay box plots
    bp = ax1.boxplot([df[df['Class'] == cls]['Score'].values for cls in classes],
                     positions=positions, widths=0.3, patch_artist=True,
                     boxprops=dict(facecolor='white', alpha=0.7),
                     medianprops=dict(color='red', linewidth=2))
    
    ax1.set_xticks(positions)
    ax1.set_xticklabels(classes, rotation=45, ha='right', fontsize=10)
    ax1.set_ylabel('Score', fontsize=12, fontweight='bold')
    ax1.set_title('Score Distribution by Class (Violin + Box Plot)', fontsize=14, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3)
    
    # Add statistical annotations
    for i, cls in enumerate(classes):
        scores = df[df['Class'] == cls]['Score'].values
        if len(scores) > 0:
            mean_score = np.mean(scores)
            std_score = np.std(scores)
            ax1.text(i, ax1.get_ylim()[1] * 0.95, 
                    f'μ={mean_score:.2f}\nσ={std_score:.2f}',
                    ha='center', va='top', fontsize=8,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    # 2. Statistical summary table
    ax2.axis('off')
    
    # Calculate statistics
    stats_data = []
    for cls in classes:
        cls_scores = df[df['Class'] == cls]['Score'].values
        if len(cls_scores) > 0:
            stats_data.append([
                cls[:15],
                len(cls_scores),
                f'{np.mean(cls_scores):.3f}',
                f'{np.median(cls_scores):.3f}',
                f'{np.std(cls_scores):.3f}',
                f'{np.min(cls_scores):.3f}',
                f'{np.max(cls_scores):.3f}'
            ])
    
    if stats_data:
        table = ax2.table(cellText=stats_data,
                         colLabels=['Class', 'N', 'Mean', 'Median', 'Std', 'Min', 'Max'],
                         cellLoc='center', loc='center',
                         colWidths=[0.25, 0.1, 0.13, 0.13, 0.13, 0.13, 0.13])
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 1.8)
        
        # Style header
        for i in range(7):
            table[(0, i)].set_facecolor('#E0E0E0')
            table[(0, i)].set_text_props(weight='bold')
    
    plt.suptitle('Score Statistics Analysis', fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    return fig


def plot_length_statistics_by_class(motifs: List[Dict[str, Any]], 
                                   figsize: Tuple[int, int] = (14, 10)) -> plt.Figure:
    """
    Advanced visualization of motif length distributions by class.
    
    Args:
        motifs: List of motif dictionaries
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object with length statistics
    """
    set_scientific_style()
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No motifs to analyze', ha='center', va='center',
               transform=ax.transAxes, fontsize=14)
        return fig
    
    # Prepare data
    df_data = []
    for motif in motifs:
        df_data.append({
            'Class': motif.get('Class', 'Unknown'),
            'Length': motif.get('Length', 0)
        })
    df = pd.DataFrame(df_data)
    
    # Create figure with subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=figsize, 
                                       gridspec_kw={'height_ratios': [2, 1.5, 1]})
    
    # 1. Histogram with KDE overlay for each class
    classes = sorted(df['Class'].unique())
    
    for cls in classes:
        cls_lengths = df[df['Class'] == cls]['Length'].values
        if len(cls_lengths) > 1:
            color = MOTIF_CLASS_COLORS.get(cls, '#808080')
            ax1.hist(cls_lengths, bins=20, alpha=0.4, label=cls, color=color, edgecolor='black')
    
    ax1.set_xlabel('Length (bp)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax1.set_title('Length Distribution by Class', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=9, ncol=2)
    ax1.grid(axis='y', alpha=0.3)
    
    # 2. Box plot comparison
    bp = ax2.boxplot([df[df['Class'] == cls]['Length'].values for cls in classes],
                     labels=classes, patch_artist=True, vert=True)
    
    # Color boxes by class
    for i, (patch, cls) in enumerate(zip(bp['boxes'], classes)):
        color = MOTIF_CLASS_COLORS.get(cls, '#808080')
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    
    ax2.set_xticklabels(classes, rotation=45, ha='right', fontsize=10)
    ax2.set_ylabel('Length (bp)', fontsize=12, fontweight='bold')
    ax2.set_title('Length Comparison (Box Plot)', fontsize=12, fontweight='bold')
    ax2.grid(axis='y', alpha=0.3)
    
    # 3. Statistical summary table
    ax3.axis('off')
    
    stats_data = []
    for cls in classes:
        cls_lengths = df[df['Class'] == cls]['Length'].values
        if len(cls_lengths) > 0:
            stats_data.append([
                cls[:15],
                len(cls_lengths),
                f'{np.mean(cls_lengths):.1f}',
                f'{np.median(cls_lengths):.1f}',
                f'{np.std(cls_lengths):.1f}',
                f'{np.min(cls_lengths):.0f}',
                f'{np.max(cls_lengths):.0f}'
            ])
    
    if stats_data:
        table = ax3.table(cellText=stats_data,
                         colLabels=['Class', 'N', 'Mean', 'Median', 'Std', 'Min', 'Max'],
                         cellLoc='center', loc='center',
                         colWidths=[0.25, 0.1, 0.13, 0.13, 0.13, 0.13, 0.13])
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 1.8)
        
        # Style header
        for i in range(7):
            table[(0, i)].set_facecolor('#E0E0E0')
            table[(0, i)].set_text_props(weight='bold')
    
    plt.suptitle('Length Statistics Analysis', fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    return fig


# =============================================================================
# TESTING & EXAMPLES
# =============================================================================

def test_visualizations():
    """Test visualization functions with example data"""
    print("Testing NBDScanner visualizations...")
    
    # Create example motif data
    example_motifs = [
        {'Class': 'G-Quadruplex', 'Subclass': 'Canonical G4', 'Start': 1, 'End': 21, 'Length': 21, 'Score': 0.85},
        {'Class': 'G-Quadruplex', 'Subclass': 'Relaxed G4', 'Start': 45, 'End': 60, 'Length': 16, 'Score': 0.72},
        {'Class': 'Curved_DNA', 'Subclass': 'A-tract', 'Start': 80, 'End': 95, 'Length': 16, 'Score': 0.65},
        {'Class': 'Z-DNA', 'Subclass': 'CG alternating', 'Start': 120, 'End': 135, 'Length': 16, 'Score': 0.90},
        {'Class': 'i-Motif', 'Subclass': 'Canonical i-motif', 'Start': 160, 'End': 180, 'Length': 21, 'Score': 0.78}
    ]
    
    sequence_length = 200
    
    print(f"\nTesting with {len(example_motifs)} example motifs:")
    for motif in example_motifs:
        print(f"  {motif['Class']} at {motif['Start']}-{motif['End']}")
    
    # Test basic plots
    try:
        fig1 = plot_motif_distribution(example_motifs, by='Class')
        plt.close(fig1)
        print("✓ Motif distribution plot: PASS")
    except Exception as e:
        print(f"✗ Motif distribution plot: FAIL - {e}")
    
    try:
        fig2 = plot_coverage_map(example_motifs, sequence_length)
        plt.close(fig2)
        print("✓ Coverage map plot: PASS")
    except Exception as e:
        print(f"✗ Coverage map plot: FAIL - {e}")
    
    try:
        fig3 = plot_score_distribution(example_motifs)
        plt.close(fig3)
        print("✓ Score distribution plot: PASS")
    except Exception as e:
        print(f"✗ Score distribution plot: FAIL - {e}")
    
    try:
        fig4 = plot_nested_pie_chart(example_motifs)
        plt.close(fig4)
        print("✓ Nested pie chart: PASS")
    except Exception as e:
        print(f"✗ Nested pie chart: FAIL - {e}")
    
    print(f"\n✓ Visualization testing completed")
    print(f"Plotly available: {'Yes' if PLOTLY_AVAILABLE else 'No'}")


# =============================================================================
# ENHANCED STATISTICS VISUALIZATIONS: DENSITY AND ENRICHMENT
# =============================================================================

def plot_density_comparison(genomic_density: Dict[str, float],
                            positional_density: Dict[str, float],
                            title: str = "Motif Density Analysis",
                            figsize: Tuple[int, int] = (14, 6)) -> plt.Figure:
    """
    Plot comparison of genomic density (coverage %) and positional density (motifs/kbp).
    
    Args:
        genomic_density: Dictionary of class -> genomic density (%)
        positional_density: Dictionary of class -> positional density (motifs/unit)
        title: Plot title
        figsize: Figure size
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    # Remove 'Overall' for class-specific comparison
    classes = [k for k in genomic_density.keys() if k != 'Overall']
    if not classes:
        classes = list(genomic_density.keys())
    
    genomic_vals = [genomic_density.get(c, 0) for c in classes]
    positional_vals = [positional_density.get(c, 0) for c in classes]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Genomic Density (Coverage %)
    colors1 = [MOTIF_CLASS_COLORS.get(c, '#808080') for c in classes]
    bars1 = ax1.barh(classes, genomic_vals, color=colors1, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax1.set_xlabel('Genomic Density (Coverage %)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Motif Class', fontsize=11, fontweight='bold')
    ax1.set_title('A. Genomic Density (σ_G)', fontsize=12, fontweight='bold', pad=10)
    
    # Add value labels
    for i, (bar, val) in enumerate(zip(bars1, genomic_vals)):
        if val > 0:
            ax1.text(val + max(genomic_vals) * 0.01, i, f'{val:.3f}%', 
                    va='center', fontsize=9, fontweight='bold')
    
    # Positional Density (Frequency)
    colors2 = [MOTIF_CLASS_COLORS.get(c, '#808080') for c in classes]
    bars2 = ax2.barh(classes, positional_vals, color=colors2, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('Positional Density (motifs/kbp)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Motif Class', fontsize=11, fontweight='bold')
    ax2.set_title('B. Positional Density (λ)', fontsize=12, fontweight='bold', pad=10)
    
    # Add value labels
    for i, (bar, val) in enumerate(zip(bars2, positional_vals)):
        if val > 0:
            ax2.text(val + max(positional_vals) * 0.01, i, f'{val:.2f}', 
                    va='center', fontsize=9, fontweight='bold')
    
    plt.suptitle(title, fontsize=14, fontweight='bold', y=1.00)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    
    return fig


def plot_enrichment_analysis(enrichment_results: Dict[str, Dict[str, Any]],
                             title: str = "Motif Enrichment Analysis",
                             figsize: Tuple[int, int] = (14, 8)) -> plt.Figure:
    """
    Plot enrichment analysis results with fold enrichment and p-values.
    
    Args:
        enrichment_results: Dictionary from calculate_enrichment_with_shuffling
        title: Plot title
        figsize: Figure size
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    # Extract data (exclude 'Overall' for class-specific view)
    classes = [k for k in enrichment_results.keys() if k != 'Overall']
    if not classes:
        classes = list(enrichment_results.keys())
    
    fold_enrichments = []
    p_values = []
    observed_densities = []
    background_means = []
    
    for cls in classes:
        result = enrichment_results[cls]
        fe = result.get('fold_enrichment', 0)
        # Handle both string 'Inf' and float infinity values robustly
        if fe == 'Inf' or (isinstance(fe, float) and np.isinf(fe)):
            fe = INFINITE_FOLD_ENRICHMENT_CAP  # Cap infinite values for visualization
        fold_enrichments.append(fe)
        p_values.append(result.get('p_value', 1.0))
        observed_densities.append(result.get('observed_density', 0))
        background_means.append(result.get('background_mean', 0))
    
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)
    
    # 1. Fold Enrichment
    ax1 = fig.add_subplot(gs[0, 0])
    colors = [MOTIF_CLASS_COLORS.get(c, '#808080') for c in classes]
    bars1 = ax1.barh(classes, fold_enrichments, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax1.axvline(x=1.0, color='red', linestyle='--', linewidth=2, label='No enrichment (FE=1)')
    ax1.set_xlabel('Fold Enrichment', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Motif Class', fontsize=11, fontweight='bold')
    ax1.set_title('A. Fold Enrichment', fontsize=12, fontweight='bold')
    ax1.legend(loc='best', fontsize=9)
    
    # Add value labels
    for i, (bar, val) in enumerate(zip(bars1, fold_enrichments)):
        label_text = f'{val:.2f}' if val < INFINITE_FOLD_ENRICHMENT_CAP else 'Inf'
        ax1.text(val + max(fold_enrichments) * 0.01, i, label_text, 
                va='center', fontsize=9, fontweight='bold')
    
    # 2. P-values
    ax2 = fig.add_subplot(gs[0, 1])
    # Color code by significance
    p_colors = ['green' if p < 0.05 else 'orange' if p < 0.1 else 'red' for p in p_values]
    bars2 = ax2.barh(classes, p_values, color=p_colors, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax2.axvline(x=0.05, color='green', linestyle='--', linewidth=2, label='p=0.05')
    ax2.axvline(x=0.1, color='orange', linestyle='--', linewidth=1.5, label='p=0.1')
    ax2.set_xlabel('P-value', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Motif Class', fontsize=11, fontweight='bold')
    ax2.set_title('B. Statistical Significance', fontsize=12, fontweight='bold')
    ax2.set_xlim(0, max(1.0, max(p_values) * 1.1))
    ax2.legend(loc='best', fontsize=9)
    
    # Add value labels
    for i, (bar, val) in enumerate(zip(bars2, p_values)):
        ax2.text(val + 0.02, i, f'{val:.3f}', va='center', fontsize=9, fontweight='bold')
    
    # 3. Observed vs Background Density
    ax3 = fig.add_subplot(gs[1, :])
    x = np.arange(len(classes))
    width = 0.35
    
    bars3a = ax3.bar(x - width/2, observed_densities, width, label='Observed', 
                    color='steelblue', alpha=0.8, edgecolor='black', linewidth=0.5)
    bars3b = ax3.bar(x + width/2, background_means, width, label='Background (Mean)', 
                    color='coral', alpha=0.8, edgecolor='black', linewidth=0.5)
    
    ax3.set_xlabel('Motif Class', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Density (%)', fontsize=11, fontweight='bold')
    ax3.set_title('C. Observed vs. Background Density Comparison', fontsize=12, fontweight='bold')
    ax3.set_xticks(x)
    ax3.set_xticklabels(classes, rotation=45, ha='right', fontsize=9)
    ax3.legend(loc='best', fontsize=10)
    ax3.grid(axis='y', alpha=0.3)
    
    plt.suptitle(title, fontsize=14, fontweight='bold', y=0.99)
    
    return fig


def plot_enrichment_summary_table(enrichment_results: Dict[str, Dict[str, Any]],
                                  title: str = "Enrichment Summary Statistics") -> plt.Figure:
    """
    Create a summary table visualization for enrichment results.
    
    Args:
        enrichment_results: Dictionary from calculate_enrichment_with_shuffling
        title: Plot title
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    # Prepare data for table
    classes = [k for k in enrichment_results.keys() if k != 'Overall']
    if not classes:
        return None
    
    table_data = []
    for cls in classes:
        result = enrichment_results[cls]
        fe = result.get('fold_enrichment', 0)
        fe_str = f"{fe:.2f}" if fe != 'Inf' else 'Inf'
        
        row = [
            cls,
            result.get('observed_count', 0),
            f"{result.get('observed_density', 0):.4f}%",
            f"{result.get('background_mean', 0):.4f}%",
            fe_str,
            f"{result.get('p_value', 1.0):.4f}",
            '***' if result.get('p_value', 1.0) < 0.001 else 
            '**' if result.get('p_value', 1.0) < 0.01 else 
            '*' if result.get('p_value', 1.0) < 0.05 else 'ns'
        ]
        table_data.append(row)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, max(6, len(classes) * 0.4)))
    ax.axis('tight')
    ax.axis('off')
    
    # Create table
    headers = ['Class', 'Count', 'Observed\nDensity', 'Background\nMean', 
               'Fold\nEnrichment', 'P-value', 'Sig.']
    
    table = ax.table(cellText=table_data, colLabels=headers, 
                    cellLoc='center', loc='center',
                    colWidths=[0.20, 0.10, 0.15, 0.15, 0.15, 0.12, 0.08])
    
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    # Style header
    for i, header in enumerate(headers):
        cell = table[(0, i)]
        cell.set_facecolor('#2196F3')
        cell.set_text_props(weight='bold', color='white', fontsize=11)
    
    # Style rows with alternating colors
    for i in range(1, len(table_data) + 1):
        for j in range(len(headers)):
            cell = table[(i, j)]
            if i % 2 == 0:
                cell.set_facecolor('#F5F5F5')
            else:
                cell.set_facecolor('white')
            
            # Highlight significant results
            if j == 6:  # Significance column
                if table_data[i-1][j] in ['***', '**', '*']:
                    cell.set_facecolor('#C8E6C9')
                    cell.set_text_props(weight='bold', color='green')
    
    plt.title(title, fontsize=14, fontweight='bold', pad=20)
    plt.tight_layout()
    
    return fig


if __name__ == "__main__":
    test_visualizations()