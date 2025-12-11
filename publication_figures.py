#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║              NATURE-LEVEL MULTI-PANEL PUBLICATION FIGURES                    ║
║        Comprehensive Multi-Panel Figures for Scientific Publications        ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: publication_figures.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.1.5
LICENSE: MIT

DESCRIPTION:
    Publication-ready multi-panel figures combining multiple visualizations
    following Nature Methods/Genetics submission guidelines.
    
FEATURES:
    - Multi-panel layouts (2x2, 2x3, 3x3, custom)
    - Automatic panel labeling (A, B, C, D, ...)
    - Consistent styling across all panels
    - Nature column widths (89mm, 120mm, 183mm)
    - 300 DPI vector output (PDF/SVG)
    - Comprehensive figure captions
    
USAGE:
    >>> from publication_figures import create_figure_1_overview
    >>> fig = create_figure_1_overview(motifs, sequence_length)
    >>> fig.savefig('figure1.pdf', dpi=300, bbox_inches='tight')
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from typing import List, Dict, Any, Tuple, Optional
from collections import Counter, defaultdict
import warnings
warnings.filterwarnings("ignore")

# Import scipy for KDE plots
try:
    from scipy import stats
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

# Import seaborn for violin/box plots
try:
    import seaborn as sns
    SEABORN_AVAILABLE = True
except ImportError:
    SEABORN_AVAILABLE = False

# Import visualization functions and utilities
from visualizations import (
    MOTIF_CLASS_COLORS,
    PUBLICATION_DPI,
    FIGURE_SIZES,
    set_scientific_style,
    _apply_nature_style,
    _format_display_name
)

# =============================================================================
# NATURE PUBLICATION CONSTANTS
# =============================================================================

# Nature journal column widths in inches (converted from mm)
NATURE_SINGLE_COLUMN = 3.5  # 89mm
NATURE_ONE_AND_HALF = 4.7   # 120mm
NATURE_DOUBLE_COLUMN = 7.2  # 183mm

# Constants for visualization configuration
CIRCOS_EXCLUDED_CLASSES = ['Non-B_DNA_Clusters', 'Hybrid']  # Classes to exclude from certain visualizations
MAX_CLASSES_FOR_SCATTER = 6  # Maximum classes to show in scatter plots for clarity
MAX_CLASSES_FOR_COMPARISON = 5  # Maximum classes for comparative plots
MAX_CLASSES_FOR_COOCCURRENCE = 7  # Maximum classes for co-occurrence matrix

# Panel label styling (Nature convention)
PANEL_LABEL_STYLE = {
    'fontsize': 12,
    'fontweight': 'bold',
    'family': 'Arial',
    'va': 'top',
    'ha': 'left'
}

# Figure caption template
CAPTION_TEMPLATE = """
Figure {fig_num}. {title}

{description}

{panel_descriptions}

All scale bars represent {scale}. Error bars indicate {error_type}.
n = {n_samples} independent measurements. Statistical significance: 
*P < 0.05, **P < 0.01, ***P < 0.001 (two-tailed Student's t-test).
"""


# =============================================================================
# PANEL LABELING AND LAYOUT UTILITIES
# =============================================================================

def add_panel_label(ax, label: str, x: float = -0.1, y: float = 1.05):
    """
    Add panel label (A, B, C, etc.) to subplot following Nature style.
    
    Args:
        ax: Matplotlib axis object
        label: Panel label (e.g., 'A', 'B', 'C')
        x: X position in axis coordinates
        y: Y position in axis coordinates
    """
    ax.text(x, y, label, transform=ax.transAxes, **PANEL_LABEL_STYLE)


def create_multipanel_figure(n_rows: int, n_cols: int, 
                             figsize: Optional[Tuple[float, float]] = None,
                             width_ratios: Optional[List[float]] = None,
                             height_ratios: Optional[List[float]] = None,
                             **kwargs) -> Tuple[plt.Figure, np.ndarray]:
    """
    Create a multi-panel figure with consistent styling.
    
    Args:
        n_rows: Number of rows
        n_cols: Number of columns
        figsize: Figure size (width, height) in inches. Auto if None.
        width_ratios: Relative widths of columns
        height_ratios: Relative heights of rows
        **kwargs: Additional arguments for gridspec
        
    Returns:
        Tuple of (figure, axes_array)
    """
    set_scientific_style('nature')
    
    # Auto-calculate figure size if not provided
    if figsize is None:
        if n_cols == 1:
            width = NATURE_SINGLE_COLUMN
        elif n_cols == 2:
            width = NATURE_DOUBLE_COLUMN
        else:
            width = NATURE_DOUBLE_COLUMN
        
        # Height proportional to rows (2.5 inches per row)
        height = n_rows * 2.5
        figsize = (width, height)
    
    # Create figure with gridspec
    fig = plt.figure(figsize=figsize, dpi=PUBLICATION_DPI)
    
    gs = gridspec.GridSpec(n_rows, n_cols, figure=fig,
                          width_ratios=width_ratios,
                          height_ratios=height_ratios,
                          hspace=kwargs.get('hspace', 0.4),
                          wspace=kwargs.get('wspace', 0.3))
    
    # Create axes array
    axes = np.empty((n_rows, n_cols), dtype=object)
    for i in range(n_rows):
        for j in range(n_cols):
            axes[i, j] = fig.add_subplot(gs[i, j])
    
    return fig, axes


# =============================================================================
# FIGURE 1: COMPREHENSIVE MOTIF OVERVIEW
# =============================================================================

def create_figure_1_overview(motifs: List[Dict[str, Any]], 
                            sequence_length: int,
                            title: str = "Comprehensive Non-B DNA Motif Analysis",
                            save_path: Optional[str] = None) -> plt.Figure:
    """
    Create Figure 1: Comprehensive overview combining distribution, 
    hierarchical, coverage, and statistics.
    
    Panel Layout (2x2):
    - A: Motif class distribution (bar chart)
    - B: Class-subclass hierarchy (nested donut)
    - C: Genomic coverage track
    - D: Score distributions by class (box plot)
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of analyzed sequence
        title: Figure title
        save_path: Optional path to save figure
        
    Returns:
        Matplotlib figure object (publication-ready at 300 DPI)
    """
    set_scientific_style('nature')
    
    # Create 2x2 multi-panel figure
    fig = plt.figure(figsize=(NATURE_DOUBLE_COLUMN, 7), dpi=PUBLICATION_DPI)
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.35)
    
    # Panel A: Motif class distribution
    ax_a = fig.add_subplot(gs[0, 0])
    add_panel_label(ax_a, 'A')
    
    # Get counts for all classes
    ALL_CLASSES = [
        'Curved_DNA', 'Slipped_DNA', 'Cruciform', 'R-Loop', 'Triplex',
        'G-Quadruplex', 'i-Motif', 'Z-DNA', 'A-philic_DNA', 'Hybrid', 
        'Non-B_DNA_Clusters'
    ]
    counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    values = [counts.get(cls, 0) for cls in ALL_CLASSES]
    colors = [MOTIF_CLASS_COLORS.get(cls, '#808080') for cls in ALL_CLASSES]
    
    bars = ax_a.bar(range(len(ALL_CLASSES)), values, color=colors, 
                    edgecolor='black', linewidth=0.5, width=0.8, alpha=0.85)
    ax_a.set_xlabel('Motif Class', fontweight='normal')
    ax_a.set_ylabel('Count', fontweight='normal')
    ax_a.set_title('Motif Class Distribution', fontweight='bold', fontsize=9)
    ax_a.set_xticks(range(len(ALL_CLASSES)))
    display_labels = [_format_display_name(cls) for cls in ALL_CLASSES]
    ax_a.set_xticklabels(display_labels, rotation=45, ha='right', fontsize=7)
    
    # Add count labels on bars
    max_val = max(values) if max(values) > 0 else 1
    for bar, count in zip(bars, values):
        if count > 0:
            height = bar.get_height()
            ax_a.text(bar.get_x() + bar.get_width()/2., height + max_val * 0.02,
                     str(count), ha='center', va='bottom', fontsize=6, fontweight='bold')
    
    _apply_nature_style(ax_a)
    
    # Panel B: Nested donut chart (class-subclass hierarchy)
    ax_b = fig.add_subplot(gs[0, 1])
    add_panel_label(ax_b, 'B')
    
    # Count by class and subclass
    class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    class_subclass_counts = defaultdict(lambda: defaultdict(int))
    
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        subclass_name = motif.get('Subclass', 'Unknown')
        class_subclass_counts[class_name][subclass_name] += 1
    
    # Inner donut (classes)
    class_names = list(class_counts.keys())
    class_values = list(class_counts.values())
    class_colors = [MOTIF_CLASS_COLORS.get(name, '#808080') for name in class_names]
    
    wedges1, texts1, autotexts1 = ax_b.pie(
        class_values, 
        labels=None,
        colors=class_colors,
        radius=0.65,
        autopct=lambda pct: f'{pct:.0f}%' if pct > 5 else '',
        pctdistance=0.80,
        startangle=90,
        wedgeprops=dict(width=0.35, edgecolor='white', linewidth=1)
    )
    
    # Add class labels
    for i, (wedge, class_name) in enumerate(zip(wedges1, class_names)):
        angle = (wedge.theta2 + wedge.theta1) / 2
        x = 0.5 * np.cos(np.radians(angle))
        y = 0.5 * np.sin(np.radians(angle))
        display_name = _format_display_name(class_name)
        ax_b.text(x, y, display_name, ha='center', va='center', 
                 fontsize=6, fontweight='bold')
    
    # Outer donut (subclasses) 
    all_subclass_counts = []
    all_subclass_colors = []
    
    for class_name in class_names:
        subclass_dict = class_subclass_counts[class_name]
        base_color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        
        for subclass_name, count in subclass_dict.items():
            all_subclass_counts.append(count)
            all_subclass_colors.append(base_color)
    
    wedges2, texts2 = ax_b.pie(
        all_subclass_counts,
        labels=None,
        colors=all_subclass_colors,
        radius=1.0,
        startangle=90,
        wedgeprops=dict(width=0.35, edgecolor='white', linewidth=0.8)
    )
    
    ax_b.set_title('Class-Subclass Hierarchy', fontweight='bold', fontsize=9)
    
    # Style percentage labels
    for autotext in autotexts1:
        autotext.set_fontsize(6)
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    
    # Panel C: Genomic coverage track  
    ax_c = fig.add_subplot(gs[1, 0])
    add_panel_label(ax_c, 'C')
    
    # Group motifs by class
    class_motifs = defaultdict(list)
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        class_motifs[class_name].append(motif)
    
    y_pos = 0
    class_positions = {}
    
    for class_name, class_motif_list in sorted(class_motifs.items()):
        class_positions[class_name] = y_pos
        color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        
        for motif in class_motif_list:
            start = motif.get('Start', 0) - 1
            end = motif.get('End', 0)
            length = end - start
            
            rect = mpatches.Rectangle(
                (start, y_pos - 0.35), length, 0.7,
                facecolor=color, edgecolor='black', linewidth=0.3, alpha=0.85
            )
            ax_c.add_patch(rect)
        
        y_pos += 1
    
    ax_c.set_xlim(0, sequence_length)
    ax_c.set_ylim(-0.5, len(class_motifs) - 0.5)
    ax_c.set_xlabel('Genomic Position (bp)', fontweight='normal', fontsize=8)
    ax_c.set_ylabel('Motif Class', fontweight='normal', fontsize=8)
    ax_c.set_title('Genomic Coverage Map', fontweight='bold', fontsize=9)
    
    ax_c.set_yticks(list(class_positions.values()))
    display_labels = [_format_display_name(label) for label in class_positions.keys()]
    ax_c.set_yticklabels(display_labels, fontsize=7)
    ax_c.ticklabel_format(style='sci', axis='x', scilimits=(3, 3))
    
    _apply_nature_style(ax_c)
    
    # Panel D: Length distribution by class
    ax_d = fig.add_subplot(gs[1, 1])
    add_panel_label(ax_d, 'D')
    
    # Extract length data
    length_data = []
    for motif in motifs:
        length = motif.get('Length')
        class_name = motif.get('Class', 'Unknown')
        if isinstance(length, int) and length > 0:
            length_data.append({
                'Length': length,
                'Class': class_name
            })
    
    if length_data and SEABORN_AVAILABLE:
        df_len = pd.DataFrame(length_data)
        classes_with_data = sorted(df_len['Class'].unique())
        colors_len = [MOTIF_CLASS_COLORS.get(cls, '#808080') for cls in classes_with_data]
        
        # Use violin plot for distribution
        sns.violinplot(data=df_len, x='Class', y='Length', ax=ax_d, palette=colors_len,
                      linewidth=0.8, inner='box', cut=0)
        display_labels_len = [_format_display_name(c) for c in classes_with_data]
        ax_d.set_xticklabels(display_labels_len, rotation=45, ha='right', fontsize=7)
        ax_d.set_ylabel('Length (bp)', fontweight='normal', fontsize=8)
        ax_d.set_xlabel('Motif Class', fontweight='normal', fontsize=8)
        ax_d.set_title('Length Distribution by Class', fontweight='bold', fontsize=9)
        _apply_nature_style(ax_d)
    else:
        # Fallback if no data
        ax_d.text(0.5, 0.5, 'No length data available', ha='center', va='center',
                 transform=ax_d.transAxes, fontsize=9)
        ax_d.axis('off')
    
    # Overall title
    fig.suptitle(title, fontsize=11, fontweight='bold', y=0.98)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    if save_path:
        fig.savefig(save_path, dpi=PUBLICATION_DPI, bbox_inches='tight',
                   facecolor='white', format='pdf')
        print(f"✓ Saved Figure 1 to {save_path}")
    
    return fig


# =============================================================================
# FIGURE 2: STATISTICAL CHARACTERIZATION
# =============================================================================

def create_figure_2_statistics(motifs: List[Dict[str, Any]],
                               title: str = "Statistical Characterization of Non-B DNA Motifs",
                               save_path: Optional[str] = None) -> plt.Figure:
    """
    Create Figure 2: Statistical analysis combining score and length distributions,
    density analysis, and correlations.
    
    Panel Layout (2x2):
    - A: Score distribution by class (violin plot)
    - B: Length distribution by class (box plot)
    - C: Score vs Length correlation (scatter)
    - D: Summary statistics table
    
    Args:
        motifs: List of motif dictionaries
        title: Figure title
        save_path: Optional path to save figure
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style('nature')
    
    fig = plt.figure(figsize=(NATURE_DOUBLE_COLUMN, 7), dpi=PUBLICATION_DPI)
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.35)
    
    # Extract data
    score_data = []
    length_data = []
    for motif in motifs:
        score = motif.get('Score', motif.get('Normalized_Score'))
        length = motif.get('Length')
        class_name = motif.get('Class', 'Unknown')
        
        if isinstance(score, (int, float)) and isinstance(length, int) and length > 0:
            score_data.append({
                'Score': score,
                'Length': length,
                'Class': class_name
            })
    
    df = pd.DataFrame(score_data)
    
    # Panel A: Score distribution by class (violin plot)
    ax_a = fig.add_subplot(gs[0, 0])
    add_panel_label(ax_a, 'A')
    
    if not df.empty and SEABORN_AVAILABLE:
        classes = sorted(df['Class'].unique())
        colors = [MOTIF_CLASS_COLORS.get(cls, '#808080') for cls in classes]
        
        sns.violinplot(data=df, x='Class', y='Score', ax=ax_a, palette=colors,
                      linewidth=0.8, inner='box')
        display_labels = [_format_display_name(c) for c in classes]
        ax_a.set_xticklabels(display_labels, rotation=45, ha='right', fontsize=7)
        ax_a.set_ylabel('Motif Score', fontweight='normal')
        ax_a.set_xlabel('Motif Class', fontweight='normal')
        ax_a.set_title('Score Distribution by Class', fontweight='bold', fontsize=9)
        _apply_nature_style(ax_a)
    
    # Panel B: Length distribution by class (box plot)
    ax_b = fig.add_subplot(gs[0, 1])
    add_panel_label(ax_b, 'B')
    
    if not df.empty and SEABORN_AVAILABLE:
        sns.boxplot(data=df, x='Class', y='Length', ax=ax_b, palette=colors,
                   linewidth=0.8, fliersize=2)
        display_labels = [_format_display_name(c) for c in classes]
        ax_b.set_xticklabels(display_labels, rotation=45, ha='right', fontsize=7)
        ax_b.set_ylabel('Motif Length (bp)', fontweight='normal')
        ax_b.set_xlabel('Motif Class', fontweight='normal')
        ax_b.set_title('Length Distribution by Class', fontweight='bold', fontsize=9)
        _apply_nature_style(ax_b)
    
    # Panel C: Score vs Length scatter
    ax_c = fig.add_subplot(gs[1, 0])
    add_panel_label(ax_c, 'C')
    
    if not df.empty:
        for class_name in classes:
            class_df = df[df['Class'] == class_name]
            color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
            display_name = _format_display_name(class_name)
            ax_c.scatter(class_df['Length'], class_df['Score'], 
                        c=color, s=20, alpha=0.6, edgecolors='black', 
                        linewidth=0.3, label=display_name)
        
        ax_c.set_xlabel('Motif Length (bp)', fontweight='normal')
        ax_c.set_ylabel('Motif Score', fontweight='normal')
        ax_c.set_title('Score-Length Correlation', fontweight='bold', fontsize=9)
        ax_c.legend(loc='best', fontsize=6, framealpha=0.9, ncol=2)
        ax_c.grid(alpha=0.3, linestyle='--', linewidth=0.5)
        _apply_nature_style(ax_c)
    
    # Panel D: Summary statistics table
    ax_d = fig.add_subplot(gs[1, 1])
    add_panel_label(ax_d, 'D')
    ax_d.axis('off')
    
    if not df.empty:
        # Calculate statistics
        stats_data = []
        for class_name in classes[:7]:  # Top 7 classes
            class_df = df[df['Class'] == class_name]
            if len(class_df) > 0:
                display_name = _format_display_name(class_name)
                if len(display_name) > 15:
                    display_name = display_name[:12] + '...'
                
                stats_data.append([
                    display_name,
                    len(class_df),
                    f"{class_df['Score'].mean():.3f}",
                    f"{class_df['Length'].mean():.1f}"
                ])
        
        table = ax_d.table(
            cellText=stats_data,
            colLabels=['Class', 'n', 'Score\n(mean)', 'Length\n(mean bp)'],
            cellLoc='left',
            loc='center',
            colWidths=[0.35, 0.15, 0.25, 0.25]
        )
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 2.5)
        
        # Style header
        for i in range(4):
            table[(0, i)].set_facecolor('#E0E0E0')
            table[(0, i)].set_text_props(weight='bold', fontsize=8)
        
        # Alternate row colors
        for i in range(1, len(stats_data) + 1):
            for j in range(4):
                if i % 2 == 0:
                    table[(i, j)].set_facecolor('#F5F5F5')
        
        ax_d.set_title('Summary Statistics', fontweight='bold', fontsize=9, pad=20)
    
    fig.suptitle(title, fontsize=11, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    if save_path:
        fig.savefig(save_path, dpi=PUBLICATION_DPI, bbox_inches='tight',
                   facecolor='white', format='pdf')
        print(f"✓ Saved Figure 2 to {save_path}")
    
    return fig


# =============================================================================
# FIGURE 3: GENOMIC LANDSCAPE AND DENSITY ANALYSIS
# =============================================================================

def create_figure_3_landscape(motifs: List[Dict[str, Any]], 
                             sequence_length: int,
                             window_size: int = None,
                             title: str = "Genomic Landscape and Density Analysis",
                             save_path: Optional[str] = None) -> plt.Figure:
    """
    Create Figure 3: Genomic landscape combining density profiles,
    Manhattan plot, and cumulative distribution.
    
    Panel Layout (3x1):
    - A: Density heatmap by class
    - B: Manhattan plot showing hotspots
    - C: Cumulative distribution curves
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of analyzed sequence
        window_size: Window size for density calculation
        title: Figure title
        save_path: Optional path to save figure
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style('nature')
    
    if window_size is None:
        window_size = max(1000, sequence_length // 50)
    
    fig = plt.figure(figsize=(NATURE_DOUBLE_COLUMN, 9), dpi=PUBLICATION_DPI)
    gs = gridspec.GridSpec(3, 1, figure=fig, hspace=0.4, height_ratios=[1.2, 1.2, 1])
    
    # Panel A: Density heatmap
    ax_a = fig.add_subplot(gs[0])
    add_panel_label(ax_a, 'A', x=-0.08)
    
    # Calculate density matrix
    num_windows = max(1, sequence_length // window_size)
    windows = np.linspace(0, sequence_length, num_windows + 1)
    
    # Exclude synthetic classes for cleaner visualization
    CIRCOS_EXCLUDED_CLASSES = ['Non-B_DNA_Clusters', 'Hybrid']
    classes = sorted(set(m.get('Class', 'Unknown') for m in motifs 
                        if m.get('Class') not in CIRCOS_EXCLUDED_CLASSES))
    
    if not classes:
        classes = sorted(set(m.get('Class', 'Unknown') for m in motifs))
    
    density_matrix = np.zeros((len(classes), num_windows))
    
    for i, class_name in enumerate(classes):
        class_motifs = [m for m in motifs if m.get('Class') == class_name]
        
        for j in range(num_windows):
            window_start = windows[j]
            window_end = windows[j + 1]
            
            count = 0
            for motif in class_motifs:
                motif_start = motif.get('Start', 0) - 1
                motif_end = motif.get('End', 0)
                
                if not (motif_end <= window_start or motif_start >= window_end):
                    count += 1
            
            density_matrix[i, j] = count
    
    im = ax_a.imshow(density_matrix, cmap='viridis', aspect='auto', interpolation='nearest')
    ax_a.set_xlabel('Position (kb)', fontweight='normal')
    ax_a.set_ylabel('Motif Class', fontweight='normal')
    ax_a.set_title('Density Heatmap Across Genome', fontweight='bold', fontsize=9)
    
    ax_a.set_yticks(range(len(classes)))
    display_classes = [_format_display_name(cls) for cls in classes]
    ax_a.set_yticklabels(display_classes, fontsize=7)
    
    x_ticks = np.arange(0, num_windows, max(1, num_windows // 5))
    ax_a.set_xticks(x_ticks)
    ax_a.set_xticklabels([f'{int(windows[i]/1000)}' for i in x_ticks], fontsize=7)
    
    cbar = plt.colorbar(im, ax=ax_a, shrink=0.8)
    cbar.set_label('Count', fontsize=8, fontweight='bold')
    cbar.ax.tick_params(labelsize=7)
    
    # Panel B: Manhattan plot
    ax_b = fig.add_subplot(gs[1])
    add_panel_label(ax_b, 'B', x=-0.08)
    
    positions_kb = []
    densities = []
    colors_list = []
    
    for class_name in classes:
        class_motifs = [m for m in motifs if m.get('Class') == class_name]
        color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        
        for i in range(num_windows):
            window_start = i * window_size
            window_end = (i + 1) * window_size
            window_center_kb = (window_start + window_end) / 2 / 1000
            
            window_motifs = []
            for motif in class_motifs:
                motif_start = motif.get('Start', 0) - 1
                motif_end = motif.get('End', 0)
                if not (motif_end <= window_start or motif_start >= window_end):
                    window_motifs.append(motif)
            
            if window_motifs:
                density = len(window_motifs) / (window_size / 1000)
                positions_kb.append(window_center_kb)
                densities.append(density)
                colors_list.append(color)
    
    ax_b.scatter(positions_kb, densities, c=colors_list, s=20, alpha=0.6, 
                edgecolors='black', linewidth=0.3)
    
    ax_b.set_xlabel('Genomic Position (kb)', fontweight='normal')
    ax_b.set_ylabel('Density (motifs/kb)', fontweight='normal')
    ax_b.set_title('Manhattan Plot - Density Hotspots', fontweight='bold', fontsize=9)
    ax_b.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    _apply_nature_style(ax_b)
    
    # Panel C: Cumulative distribution
    ax_c = fig.add_subplot(gs[2])
    add_panel_label(ax_c, 'C', x=-0.08)
    
    for class_name in classes[:5]:  # Top 5 classes for clarity
        class_motifs = [m for m in motifs if m.get('Class') == class_name]
        class_motifs_sorted = sorted(class_motifs, key=lambda m: m.get('Start', 0))
        
        positions = [0]
        cumulative = [0]
        
        for i, motif in enumerate(class_motifs_sorted, 1):
            positions.append(motif.get('Start', 0) / 1000)
            cumulative.append(i)
        
        positions.append(sequence_length / 1000)
        cumulative.append(len(class_motifs))
        
        color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        display_name = class_name.replace('_', ' ')
        ax_c.plot(positions, cumulative, color=color, linewidth=1.5, 
                 label=display_name, alpha=0.8)
    
    ax_c.set_xlabel('Genomic Position (kb)', fontweight='normal')
    ax_c.set_ylabel('Cumulative Count', fontweight='normal')
    ax_c.set_title('Cumulative Distribution', fontweight='bold', fontsize=9)
    ax_c.legend(loc='upper left', fontsize=7, framealpha=0.9)
    ax_c.grid(alpha=0.3, linestyle='--', linewidth=0.5)
    _apply_nature_style(ax_c)
    
    fig.suptitle(title, fontsize=11, fontweight='bold', y=0.99)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    
    if save_path:
        fig.savefig(save_path, dpi=PUBLICATION_DPI, bbox_inches='tight',
                   facecolor='white', format='pdf')
        print(f"✓ Saved Figure 3 to {save_path}")
    
    return fig


# =============================================================================
# FIGURE 4: DETAILED SUBCLASS ANALYSIS
# =============================================================================

def create_figure_4_subclass_analysis(motifs: List[Dict[str, Any]],
                                     sequence_length: int,
                                     title: str = "Detailed Subclass-Level Analysis",
                                     save_path: Optional[str] = None) -> plt.Figure:
    """
    Create Figure 4: Detailed subclass-level analysis with enhanced visualizations.
    
    Panel Layout (2x2):
    - A: Subclass distribution (bar chart)
    - B: Subclass density heatmap
    - C: Length distribution by subclass (KDE)
    - D: Motif co-occurrence matrix
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of analyzed sequence
        title: Figure title
        save_path: Optional path to save figure
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style('nature')
    
    fig = plt.figure(figsize=(NATURE_DOUBLE_COLUMN, 7.5), dpi=PUBLICATION_DPI)
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.35)
    
    # Panel A: Subclass distribution
    ax_a = fig.add_subplot(gs[0, 0])
    add_panel_label(ax_a, 'A')
    
    # Get top 15 subclasses
    subclass_counts = Counter()
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        subclass_name = motif.get('Subclass', 'Unknown')
        key = f"{class_name}:{subclass_name}"
        subclass_counts[key] += 1
    
    top_subclasses = subclass_counts.most_common(15)
    labels = [sc[0] for sc in top_subclasses]
    values = [sc[1] for sc in top_subclasses]
    
    # Get colors based on parent class
    colors = []
    for label in labels:
        if ':' in label:
            parent_class = label.split(':')[0]
            colors.append(MOTIF_CLASS_COLORS.get(parent_class, '#808080'))
        else:
            colors.append('#808080')
    
    bars = ax_a.barh(range(len(labels)), values, color=colors, 
                     edgecolor='black', linewidth=0.5, alpha=0.85)
    
    ax_a.set_yticks(range(len(labels)))
    display_labels = []
    for label in labels:
        parts = label.split(':')
        if len(parts) == 2:
            display_labels.append(f"{_format_display_name(parts[0])}: {parts[1]}")
        else:
            display_labels.append(_format_display_name(label))
    ax_a.set_yticklabels(display_labels, fontsize=7)
    ax_a.set_xlabel('Count', fontweight='normal')
    ax_a.set_ylabel('Subclass', fontweight='normal')
    ax_a.set_title('Top 15 Subclasses', fontweight='bold', fontsize=9)
    ax_a.invert_yaxis()
    
    # Add count labels
    max_val = max(values) if values else 1
    for i, (bar, count) in enumerate(zip(bars, values)):
        ax_a.text(count + max_val * 0.02, i, str(count), 
                 va='center', fontsize=7, fontweight='bold')
    
    _apply_nature_style(ax_a)
    
    # Panel B: Subclass density heatmap (simplified)
    ax_b = fig.add_subplot(gs[0, 1])
    add_panel_label(ax_b, 'B')
    
    # Calculate density for top classes
    window_size = max(1000, sequence_length // 20)
    num_windows = max(1, sequence_length // window_size)
    
    # Use top 8 subclasses for clarity
    top_8_subclasses = [sc[0] for sc in subclass_counts.most_common(8)]
    density_matrix = np.zeros((len(top_8_subclasses), num_windows))
    
    for i, subclass_key in enumerate(top_8_subclasses):
        class_name, subclass_name = subclass_key.split(':') if ':' in subclass_key else (subclass_key, '')
        subclass_motifs = [m for m in motifs 
                          if m.get('Class') == class_name and 
                          (not subclass_name or m.get('Subclass') == subclass_name)]
        
        for j in range(num_windows):
            window_start = j * window_size
            window_end = (j + 1) * window_size
            
            count = sum(1 for m in subclass_motifs 
                       if not (m.get('End', 0) <= window_start or 
                              m.get('Start', 0) - 1 >= window_end))
            density_matrix[i, j] = count
    
    im = ax_b.imshow(density_matrix, cmap='YlOrRd', aspect='auto', interpolation='nearest')
    
    ax_b.set_yticks(range(len(top_8_subclasses)))
    display_sc = []
    for sc in top_8_subclasses:
        parts = sc.split(':')
        if len(parts) == 2:
            display_sc.append(f"{parts[1][:12]}")  # Short subclass name
        else:
            display_sc.append(sc.replace('_', ' ')[:12])
    ax_b.set_yticklabels(display_sc, fontsize=7)
    
    x_ticks = np.arange(0, num_windows, max(1, num_windows // 4))
    ax_b.set_xticks(x_ticks)
    ax_b.set_xticklabels([f'{int(i*window_size/1000)}' for i in x_ticks], fontsize=7)
    
    ax_b.set_xlabel('Position (kb)', fontweight='normal')
    ax_b.set_ylabel('Subclass', fontweight='normal')
    ax_b.set_title('Density Heatmap', fontweight='bold', fontsize=9)
    
    cbar = plt.colorbar(im, ax=ax_b, shrink=0.8)
    cbar.set_label('Count', fontsize=7, fontweight='bold')
    cbar.ax.tick_params(labelsize=6)
    
    # Panel C: Length KDE by class
    ax_c = fig.add_subplot(gs[1, 0])
    add_panel_label(ax_c, 'C')
    
    # Get unique classes (not subclasses for clarity)
    CIRCOS_EXCLUDED = ['Non-B_DNA_Clusters', 'Hybrid']
    classes = sorted(set(m.get('Class', 'Unknown') for m in motifs 
                        if m.get('Class') not in CIRCOS_EXCLUDED))
    
    if not classes:
        classes = sorted(set(m.get('Class', 'Unknown') for m in motifs))
    
    # Plot KDE for top 5 classes
    for class_name in classes[:5]:
        class_motifs = [m for m in motifs if m.get('Class') == class_name]
        lengths = [m.get('Length', 0) for m in class_motifs if m.get('Length', 0) > 0]
        
        if len(lengths) > 2:
            color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
            display_name = class_name.replace('_', ' ')
            
            try:
                from scipy import stats
                kde = stats.gaussian_kde(lengths)
                x_range = np.linspace(min(lengths), max(lengths), 200)
                density = kde(x_range)
                ax_c.plot(x_range, density, color=color, linewidth=1.5, 
                         label=display_name, alpha=0.8)
                ax_c.fill_between(x_range, density, alpha=0.2, color=color)
            except:
                # Fallback to histogram
                ax_c.hist(lengths, bins=20, alpha=0.3, color=color, 
                         label=display_name, density=True, edgecolor='black', linewidth=0.5)
    
    ax_c.set_xlabel('Motif Length (bp)', fontweight='normal')
    ax_c.set_ylabel('Density', fontweight='normal')
    ax_c.set_title('Length Distribution (KDE)', fontweight='bold', fontsize=9)
    ax_c.legend(loc='upper right', fontsize=7, framealpha=0.9)
    ax_c.grid(alpha=0.3, linestyle='--', linewidth=0.5)
    _apply_nature_style(ax_c)
    
    # Panel D: Co-occurrence matrix (simplified)
    ax_d = fig.add_subplot(gs[1, 1])
    add_panel_label(ax_d, 'D')
    
    # Use top 6 classes for matrix
    top_classes = classes[:6]
    n_classes = len(top_classes)
    cooccur_matrix = np.zeros((n_classes, n_classes))
    
    # Calculate co-occurrences (within 100bp)
    for i, class_i in enumerate(top_classes):
        motifs_i = [m for m in motifs if m.get('Class') == class_i]
        
        for j, class_j in enumerate(top_classes):
            motifs_j = [m for m in motifs if m.get('Class') == class_j]
            
            count = 0
            for mi in motifs_i:
                start_i = mi.get('Start', 0)
                end_i = mi.get('End', 0)
                
                for mj in motifs_j:
                    start_j = mj.get('Start', 0)
                    end_j = mj.get('End', 0)
                    
                    # Check overlap or proximity (within 100bp)
                    distance = max(0, max(start_i, start_j) - min(end_i, end_j))
                    if distance <= 100:
                        count += 1
            
            cooccur_matrix[i, j] = count
    
    im = ax_d.imshow(cooccur_matrix, cmap='Blues', aspect='auto', interpolation='nearest')
    
    ax_d.set_xticks(range(n_classes))
    ax_d.set_yticks(range(n_classes))
    
    display_classes = [_format_display_name(c) for c in top_classes]
    ax_d.set_xticklabels(display_classes, rotation=45, ha='right', fontsize=7)
    ax_d.set_yticklabels(display_classes, fontsize=7)
    
    # Add values to cells
    for i in range(n_classes):
        for j in range(n_classes):
            value = int(cooccur_matrix[i, j])
            if value > 0:
                text_color = 'white' if value > cooccur_matrix.max() / 2 else 'black'
                ax_d.text(j, i, str(value), ha='center', va='center', 
                         color=text_color, fontsize=7, fontweight='bold')
    
    ax_d.set_title('Co-occurrence Matrix', fontweight='bold', fontsize=9)
    ax_d.set_xlabel('Motif Class', fontweight='normal', fontsize=8)
    ax_d.set_ylabel('Motif Class', fontweight='normal', fontsize=8)
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax_d, shrink=0.8)
    cbar.set_label('Count', fontsize=7, fontweight='bold')
    cbar.ax.tick_params(labelsize=6)
    
    fig.suptitle(title, fontsize=11, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    if save_path:
        fig.savefig(save_path, dpi=PUBLICATION_DPI, bbox_inches='tight',
                   facecolor='white', format='pdf')
        print(f"✓ Saved Figure 4 to {save_path}")
    
    return fig


# =============================================================================
# FIGURE 5: MOTIF QUALITY AND CHARACTERISTICS
# =============================================================================

def create_figure_5_quality_metrics(motifs: List[Dict[str, Any]],
                                    sequence_length: int,
                                    title: str = "Motif Quality and Characteristic Features",
                                    save_path: Optional[str] = None) -> plt.Figure:
    """
    Create Figure 5: Comprehensive quality metrics showing score distributions,
    length characteristics, density patterns, and positional biases.
    
    Panel Layout (2x2):
    - A: Score distributions across classes (box plot with overlaid points)
    - B: Length vs Score correlation (scatter with regression)
    - C: Positional density heatmap across genome windows
    - D: Class abundance and diversity metrics (bar + line combo)
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of analyzed sequence
        title: Figure title
        save_path: Optional path to save figure
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style('nature')
    
    fig = plt.figure(figsize=(NATURE_DOUBLE_COLUMN, 7.5), dpi=PUBLICATION_DPI)
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.35)
    
    # Extract data
    classes = sorted(set(m.get('Class', 'Unknown') for m in motifs 
                        if m.get('Class') not in CIRCOS_EXCLUDED_CLASSES))
    
    if not classes:
        classes = sorted(set(m.get('Class', 'Unknown') for m in motifs))
    
    # Panel A: Score distributions with data points
    ax_a = fig.add_subplot(gs[0, 0])
    add_panel_label(ax_a, 'A')
    
    score_data = []
    for motif in motifs:
        score = motif.get('Score', motif.get('Normalized_Score'))
        class_name = motif.get('Class', 'Unknown')
        if isinstance(score, (int, float)) and class_name in classes:
            score_data.append({
                'Score': score,
                'Class': class_name
            })
    
    if score_data and SEABORN_AVAILABLE:
        df_score = pd.DataFrame(score_data)
        colors = [MOTIF_CLASS_COLORS.get(cls, '#808080') for cls in classes]
        
        # Box plot with overlaid strip plot
        sns.boxplot(data=df_score, x='Class', y='Score', ax=ax_a, palette=colors,
                   linewidth=0.8, fliersize=0, width=0.6)
        sns.stripplot(data=df_score, x='Class', y='Score', ax=ax_a, 
                     color='black', alpha=0.3, size=2, jitter=True)
        
        display_labels = [_format_display_name(c) for c in classes]
        ax_a.set_xticklabels(display_labels, rotation=45, ha='right', fontsize=7)
        ax_a.set_ylabel('Motif Score', fontweight='normal', fontsize=8)
        ax_a.set_xlabel('Motif Class', fontweight='normal', fontsize=8)
        ax_a.set_title('Score Distribution with Data Points', fontweight='bold', fontsize=9)
        ax_a.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
        _apply_nature_style(ax_a)
    
    # Panel B: Length vs Score correlation
    ax_b = fig.add_subplot(gs[0, 1])
    add_panel_label(ax_b, 'B')
    
    corr_data = []
    for motif in motifs:
        score = motif.get('Score', motif.get('Normalized_Score'))
        length = motif.get('Length')
        class_name = motif.get('Class', 'Unknown')
        if (isinstance(score, (int, float)) and isinstance(length, int) and 
            length > 0 and class_name in classes):
            corr_data.append({
                'Score': score,
                'Length': length,
                'Class': class_name
            })
    
    if corr_data:
        df_corr = pd.DataFrame(corr_data)
        
        # Scatter plot by class (show top MAX_CLASSES_FOR_SCATTER for clarity)
        for class_name in classes[:MAX_CLASSES_FOR_SCATTER]:
            class_df = df_corr[df_corr['Class'] == class_name]
            if len(class_df) > 0:
                color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
                display_name = _format_display_name(class_name)
                ax_b.scatter(class_df['Length'], class_df['Score'], 
                           c=color, s=25, alpha=0.6, edgecolors='black', 
                           linewidth=0.3, label=display_name)
        
        # Add regression line
        if len(df_corr) > 1:
            z = np.polyfit(df_corr['Length'], df_corr['Score'], 1)
            p = np.poly1d(z)
            x_line = np.linspace(df_corr['Length'].min(), df_corr['Length'].max(), 100)
            ax_b.plot(x_line, p(x_line), 'r--', linewidth=1.5, alpha=0.8,
                     label=f'Fit: y={z[0]:.4f}x+{z[1]:.2f}')
            
            # Calculate R²
            correlation = np.corrcoef(df_corr['Length'], df_corr['Score'])[0, 1]
            r_squared = correlation ** 2
            ax_b.text(0.05, 0.95, f'R² = {r_squared:.3f}\nρ = {correlation:.3f}', 
                     transform=ax_b.transAxes,
                     fontsize=8, verticalalignment='top',
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        ax_b.set_xlabel('Motif Length (bp)', fontweight='normal', fontsize=8)
        ax_b.set_ylabel('Motif Score', fontweight='normal', fontsize=8)
        ax_b.set_title('Length-Score Relationship', fontweight='bold', fontsize=9)
        ax_b.legend(loc='lower right', fontsize=6, framealpha=0.9, ncol=2)
        ax_b.grid(alpha=0.3, linestyle='--', linewidth=0.5)
        _apply_nature_style(ax_b)
    
    # Panel C: Positional density heatmap
    ax_c = fig.add_subplot(gs[1, 0])
    add_panel_label(ax_c, 'C')
    
    # Calculate density matrix
    window_size = max(1000, sequence_length // 30)
    num_windows = max(1, sequence_length // window_size)
    
    density_matrix = np.zeros((len(classes), num_windows))
    
    for i, class_name in enumerate(classes):
        class_motifs = [m for m in motifs if m.get('Class') == class_name]
        
        for j in range(num_windows):
            window_start = j * window_size
            window_end = (j + 1) * window_size
            
            count = sum(1 for m in class_motifs 
                       if not (m.get('End', 0) <= window_start or 
                              m.get('Start', 0) - 1 >= window_end))
            density_matrix[i, j] = count
    
    im = ax_c.imshow(density_matrix, cmap='RdYlBu_r', aspect='auto', 
                     interpolation='bilinear')
    
    ax_c.set_yticks(range(len(classes)))
    display_classes = [_format_display_name(c) for c in classes]
    ax_c.set_yticklabels(display_classes, fontsize=7)
    
    x_ticks = np.arange(0, num_windows, max(1, num_windows // 5))
    ax_c.set_xticks(x_ticks)
    ax_c.set_xticklabels([f'{int(i*window_size/1000)}' for i in x_ticks], fontsize=7)
    
    ax_c.set_xlabel('Position (kb)', fontweight='normal', fontsize=8)
    ax_c.set_ylabel('Motif Class', fontweight='normal', fontsize=8)
    ax_c.set_title('Positional Density Heatmap', fontweight='bold', fontsize=9)
    
    cbar = plt.colorbar(im, ax=ax_c, shrink=0.8)
    cbar.set_label('Count', fontsize=7, fontweight='bold')
    cbar.ax.tick_params(labelsize=6)
    
    # Panel D: Class abundance and diversity
    ax_d = fig.add_subplot(gs[1, 1])
    add_panel_label(ax_d, 'D')
    
    # Calculate metrics
    class_counts = Counter(m.get('Class', 'Unknown') for m in motifs 
                          if m.get('Class') in classes)
    
    sorted_classes = sorted(classes, key=lambda x: class_counts.get(x, 0), reverse=True)
    counts = [class_counts.get(c, 0) for c in sorted_classes]
    
    # Calculate diversity (Shannon entropy normalized)
    total = sum(counts)
    if total > 0:
        diversity_values = []
        for c in sorted_classes:
            count = class_counts.get(c, 0)
            if count > 0:
                p = count / total
                diversity_values.append(-p * np.log(p))
            else:
                diversity_values.append(0)
    else:
        diversity_values = [0] * len(sorted_classes)
    
    # Create twin axes for dual y-axis
    ax_d_twin = ax_d.twinx()
    
    x_pos = np.arange(len(sorted_classes))
    colors_d = [MOTIF_CLASS_COLORS.get(c, '#808080') for c in sorted_classes]
    
    # Bar plot for counts
    bars = ax_d.bar(x_pos, counts, color=colors_d, alpha=0.7, 
                    edgecolor='black', linewidth=0.5, label='Count')
    
    # Line plot for diversity
    line = ax_d_twin.plot(x_pos, diversity_values, 'o-', color='darkred', 
                          linewidth=2, markersize=6, label='Diversity (H)')
    
    ax_d.set_xlabel('Motif Class', fontweight='normal', fontsize=8)
    ax_d.set_ylabel('Motif Count', fontweight='normal', fontsize=8, color='black')
    ax_d_twin.set_ylabel('Shannon Diversity', fontweight='normal', fontsize=8, color='darkred')
    ax_d.set_title('Abundance and Diversity Metrics', fontweight='bold', fontsize=9)
    
    ax_d.set_xticks(x_pos)
    # Format labels efficiently (avoid double computation)
    display_labels_d = []
    for c in sorted_classes:
        formatted_name = _format_display_name(c)
        if len(formatted_name) > 10:
            display_labels_d.append(formatted_name[:10] + '...')
        else:
            display_labels_d.append(formatted_name)
    ax_d.set_xticklabels(display_labels_d, rotation=45, ha='right', fontsize=7)
    
    ax_d.tick_params(axis='y', labelcolor='black')
    ax_d_twin.tick_params(axis='y', labelcolor='darkred')
    
    # Combined legend
    lines_labels = [bars] + line
    labels = ['Count', 'Diversity']
    ax_d.legend(lines_labels, labels, loc='upper right', fontsize=7, framealpha=0.9)
    
    ax_d.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    _apply_nature_style(ax_d)
    
    fig.suptitle(title, fontsize=11, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    if save_path:
        fig.savefig(save_path, dpi=PUBLICATION_DPI, bbox_inches='tight',
                   facecolor='white', format='pdf')
        print(f"✓ Saved Figure 5 to {save_path}")
    
    return fig


# =============================================================================
# FIGURE 6: COMPARATIVE CROSS-CLASS ANALYSIS
# =============================================================================

def create_figure_6_comparative_analysis(motifs: List[Dict[str, Any]],
                                        sequence_length: int,
                                        title: str = "Comparative Cross-Class Analysis",
                                        save_path: Optional[str] = None) -> plt.Figure:
    """
    Create Figure 6: Comparative analysis showing relationships and patterns
    across different motif classes.
    
    Panel Layout (2x2):
    - A: Multi-class cumulative distribution curves
    - B: Class co-occurrence network/matrix
    - C: Relative enrichment comparison (normalized densities)
    - D: Length distribution comparison (overlapping KDE)
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of analyzed sequence
        title: Figure title
        save_path: Optional path to save figure
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style('nature')
    
    fig = plt.figure(figsize=(NATURE_DOUBLE_COLUMN, 7.5), dpi=PUBLICATION_DPI)
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.35)
    
    # Get classes
    classes = sorted(set(m.get('Class', 'Unknown') for m in motifs 
                        if m.get('Class') not in CIRCOS_EXCLUDED_CLASSES))
    
    if not classes:
        classes = sorted(set(m.get('Class', 'Unknown') for m in motifs))
    
    # Panel A: Cumulative distributions (show top MAX_CLASSES_FOR_COMPARISON)
    ax_a = fig.add_subplot(gs[0, 0])
    add_panel_label(ax_a, 'A')
    
    for class_name in classes[:MAX_CLASSES_FOR_COMPARISON]:
        class_motifs = [m for m in motifs if m.get('Class') == class_name]
        class_motifs_sorted = sorted(class_motifs, key=lambda m: m.get('Start', 0))
        
        if len(class_motifs_sorted) > 0:
            positions = [0]
            cumulative = [0]
            
            for i, motif in enumerate(class_motifs_sorted, 1):
                positions.append(motif.get('Start', 0) / 1000)
                cumulative.append(i)
            
            positions.append(sequence_length / 1000)
            cumulative.append(len(class_motifs))
            
            color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
            display_name = _format_display_name(class_name)
            ax_a.plot(positions, cumulative, color=color, linewidth=2, 
                     label=display_name, alpha=0.8)
    
    ax_a.set_xlabel('Genomic Position (kb)', fontweight='normal', fontsize=8)
    ax_a.set_ylabel('Cumulative Motif Count', fontweight='normal', fontsize=8)
    ax_a.set_title('Cumulative Distribution Comparison', fontweight='bold', fontsize=9)
    ax_a.legend(loc='upper left', fontsize=7, framealpha=0.9, ncol=2)
    ax_a.grid(alpha=0.3, linestyle='--', linewidth=0.5)
    _apply_nature_style(ax_a)
    
    # Panel B: Co-occurrence heatmap (use MAX_CLASSES_FOR_COOCCURRENCE)
    ax_b = fig.add_subplot(gs[0, 1])
    add_panel_label(ax_b, 'B')
    
    # Use top MAX_CLASSES_FOR_COOCCURRENCE classes for matrix
    top_classes = classes[:MAX_CLASSES_FOR_COOCCURRENCE]
    n_classes = len(top_classes)
    cooccur_matrix = np.zeros((n_classes, n_classes))
    
    # Calculate co-occurrences (within 50bp)
    for i, class_i in enumerate(top_classes):
        motifs_i = [m for m in motifs if m.get('Class') == class_i]
        
        for j, class_j in enumerate(top_classes):
            motifs_j = [m for m in motifs if m.get('Class') == class_j]
            
            count = 0
            for mi in motifs_i:
                start_i = mi.get('Start', 0)
                end_i = mi.get('End', 0)
                
                for mj in motifs_j:
                    start_j = mj.get('Start', 0)
                    end_j = mj.get('End', 0)
                    
                    # Calculate actual inter-motif distance
                    # For overlapping motifs: distance = 0
                    # For non-overlapping motifs: distance = gap between them
                    if end_i <= start_j:
                        distance = start_j - end_i
                    elif end_j <= start_i:
                        distance = start_i - end_j
                    else:
                        distance = 0  # Overlapping
                    
                    if distance <= 50:
                        count += 1
            
            cooccur_matrix[i, j] = count
    
    # Normalize by diagonal (self-occurrences)
    norm_matrix = np.zeros_like(cooccur_matrix)
    for i in range(n_classes):
        for j in range(n_classes):
            if cooccur_matrix[i, i] > 0 and cooccur_matrix[j, j] > 0:
                norm_matrix[i, j] = cooccur_matrix[i, j] / np.sqrt(cooccur_matrix[i, i] * cooccur_matrix[j, j])
    
    im = ax_b.imshow(norm_matrix, cmap='YlOrRd', aspect='auto', 
                     interpolation='nearest', vmin=0, vmax=1)
    
    ax_b.set_xticks(range(n_classes))
    ax_b.set_yticks(range(n_classes))
    
    display_classes = [_format_display_name(c)[:8] for c in top_classes]
    ax_b.set_xticklabels(display_classes, rotation=45, ha='right', fontsize=7)
    ax_b.set_yticklabels(display_classes, fontsize=7)
    
    ax_b.set_title('Normalized Co-occurrence', fontweight='bold', fontsize=9)
    ax_b.set_xlabel('Motif Class', fontweight='normal', fontsize=8)
    ax_b.set_ylabel('Motif Class', fontweight='normal', fontsize=8)
    
    cbar = plt.colorbar(im, ax=ax_b, shrink=0.8)
    cbar.set_label('Similarity', fontsize=7, fontweight='bold')
    cbar.ax.tick_params(labelsize=6)
    
    # Panel C: Relative enrichment (density comparison)
    ax_c = fig.add_subplot(gs[1, 0])
    add_panel_label(ax_c, 'C')
    
    # Calculate densities
    class_densities = {}
    for class_name in classes:
        class_motifs = [m for m in motifs if m.get('Class') == class_name]
        total_coverage = sum(m.get('Length', 0) for m in class_motifs)
        density = (total_coverage / sequence_length) * 100  # Percentage
        class_densities[class_name] = density
    
    # Sort by density
    sorted_classes = sorted(classes, key=lambda x: class_densities.get(x, 0), reverse=True)
    densities = [class_densities.get(c, 0) for c in sorted_classes]
    
    # Calculate relative enrichment (fold over mean)
    mean_density = np.mean(densities) if densities else 0
    enrichments = [d / mean_density if mean_density > 0 else 0 for d in densities]
    
    x_pos = np.arange(len(sorted_classes))
    colors_c = [MOTIF_CLASS_COLORS.get(c, '#808080') for c in sorted_classes]
    
    bars = ax_c.barh(x_pos, enrichments, color=colors_c, alpha=0.8, 
                     edgecolor='black', linewidth=0.5)
    
    # Add reference line at 1.0 (mean)
    ax_c.axvline(x=1.0, color='red', linestyle='--', linewidth=1.5, 
                label='Mean enrichment', alpha=0.8)
    
    ax_c.set_yticks(x_pos)
    display_labels_c = [_format_display_name(c) for c in sorted_classes]
    ax_c.set_yticklabels(display_labels_c, fontsize=7)
    ax_c.set_xlabel('Fold Enrichment (vs. mean)', fontweight='normal', fontsize=8)
    ax_c.set_ylabel('Motif Class', fontweight='normal', fontsize=8)
    ax_c.set_title('Relative Coverage Enrichment', fontweight='bold', fontsize=9)
    ax_c.legend(loc='lower right', fontsize=7, framealpha=0.9)
    ax_c.invert_yaxis()
    
    # Add value labels
    max_enrich = max(enrichments) if enrichments else 1
    for i, (bar, val) in enumerate(zip(bars, enrichments)):
        if val > 0:
            ax_c.text(val + max_enrich * 0.02, i, f'{val:.2f}x', 
                     va='center', fontsize=7, fontweight='bold')
    
    ax_c.grid(axis='x', alpha=0.3, linestyle='--', linewidth=0.5)
    _apply_nature_style(ax_c)
    
    # Panel D: Length distribution comparison (KDE) - show top MAX_CLASSES_FOR_COMPARISON
    ax_d = fig.add_subplot(gs[1, 1])
    add_panel_label(ax_d, 'D')
    
    if SCIPY_AVAILABLE:
        for class_name in classes[:MAX_CLASSES_FOR_COMPARISON]:
            class_motifs = [m for m in motifs if m.get('Class') == class_name]
            lengths = [m.get('Length', 0) for m in class_motifs if m.get('Length', 0) > 0]
            
            if len(lengths) > 2:
                color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
                display_name = _format_display_name(class_name)
                
                try:
                    kde = stats.gaussian_kde(lengths)
                    x_range = np.linspace(min(lengths), max(lengths), 200)
                    density = kde(x_range)
                    ax_d.plot(x_range, density, color=color, linewidth=2, 
                             label=display_name, alpha=0.8)
                    ax_d.fill_between(x_range, density, alpha=0.2, color=color)
                except (ValueError, np.linalg.LinAlgError):
                    # KDE can fail for very small or uniform datasets
                    pass
        
        ax_d.set_xlabel('Motif Length (bp)', fontweight='normal', fontsize=8)
        ax_d.set_ylabel('Probability Density', fontweight='normal', fontsize=8)
        ax_d.set_title('Length Distribution Overlay', fontweight='bold', fontsize=9)
        ax_d.legend(loc='upper right', fontsize=7, framealpha=0.9)
        ax_d.grid(alpha=0.3, linestyle='--', linewidth=0.5)
        _apply_nature_style(ax_d)
    else:
        ax_d.text(0.5, 0.5, 'SciPy required for KDE', ha='center', va='center',
                 transform=ax_d.transAxes, fontsize=9)
        ax_d.axis('off')
    
    fig.suptitle(title, fontsize=11, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    if save_path:
        fig.savefig(save_path, dpi=PUBLICATION_DPI, bbox_inches='tight',
                   facecolor='white', format='pdf')
        print(f"✓ Saved Figure 6 to {save_path}")
    
    return fig


# =============================================================================
# SUPPLEMENTARY FIGURE S1: COMPREHENSIVE 6-PANEL OVERVIEW
# =============================================================================

def create_supplementary_figure_s1(motifs: List[Dict[str, Any]],
                                   sequence_length: int,
                                   title: str = "Supplementary Figure S1: Comprehensive Analysis Overview",
                                   save_path: Optional[str] = None) -> plt.Figure:
    """
    Create Supplementary Figure S1: Comprehensive 6-panel overview showing
    all major aspects of motif analysis in a single figure.
    
    Panel Layout (3x2):
    - A: Class distribution bar chart
    - B: Hierarchical donut chart
    - C: Genomic coverage track
    - D: Score vs length scatter
    - E: Density heatmap
    - F: Summary statistics table
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of analyzed sequence
        title: Figure title
        save_path: Optional path to save figure
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style('nature')
    
    fig = plt.figure(figsize=(NATURE_DOUBLE_COLUMN, 10), dpi=PUBLICATION_DPI)
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.45, wspace=0.35)
    
    # Get classes
    ALL_CLASSES = [
        'Curved_DNA', 'Slipped_DNA', 'Cruciform', 'R-Loop', 'Triplex',
        'G-Quadruplex', 'i-Motif', 'Z-DNA', 'A-philic_DNA', 'Hybrid', 
        'Non-B_DNA_Clusters'
    ]
    
    # Panel A: Class distribution
    ax_a = fig.add_subplot(gs[0, 0])
    add_panel_label(ax_a, 'A')
    
    counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    values = [counts.get(cls, 0) for cls in ALL_CLASSES]
    colors = [MOTIF_CLASS_COLORS.get(cls, '#808080') for cls in ALL_CLASSES]
    
    bars = ax_a.bar(range(len(ALL_CLASSES)), values, color=colors, 
                    edgecolor='black', linewidth=0.5, width=0.75, alpha=0.85)
    
    ax_a.set_xticks(range(len(ALL_CLASSES)))
    display_labels = [_format_display_name(cls)[:10] for cls in ALL_CLASSES]
    ax_a.set_xticklabels(display_labels, rotation=45, ha='right', fontsize=6)
    ax_a.set_xlabel('Motif Class', fontweight='normal', fontsize=8)
    ax_a.set_ylabel('Count', fontweight='normal', fontsize=8)
    ax_a.set_title('A. Class Distribution', fontweight='bold', fontsize=9, loc='left')
    ax_a.grid(axis='y', alpha=0.3)
    _apply_nature_style(ax_a)
    
    # Panel B: Hierarchical donut
    ax_b = fig.add_subplot(gs[0, 1])
    add_panel_label(ax_b, 'B')
    
    class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    class_names = list(class_counts.keys())
    class_values = list(class_counts.values())
    class_colors = [MOTIF_CLASS_COLORS.get(name, '#808080') for name in class_names]
    
    wedges, texts, autotexts = ax_b.pie(
        class_values, 
        colors=class_colors,
        autopct=lambda pct: f'{pct:.0f}%' if pct > 5 else '',
        startangle=90,
        wedgeprops=dict(width=0.4, edgecolor='white', linewidth=1)
    )
    
    for autotext in autotexts:
        autotext.set_fontsize(6)
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    
    ax_b.set_title('B. Class Hierarchy', fontweight='bold', fontsize=9, loc='left', pad=10)
    
    # Panel C: Genomic coverage
    ax_c = fig.add_subplot(gs[1, :])
    add_panel_label(ax_c, 'C', x=-0.05)
    
    class_motifs = defaultdict(list)
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        class_motifs[class_name].append(motif)
    
    y_pos = 0
    for class_name, class_motif_list in sorted(class_motifs.items()):
        color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        
        for motif in class_motif_list:
            start = motif.get('Start', 0) - 1
            end = motif.get('End', 0)
            length = end - start
            
            rect = mpatches.Rectangle(
                (start, y_pos - 0.35), length, 0.7,
                facecolor=color, edgecolor='black', linewidth=0.2, alpha=0.85
            )
            ax_c.add_patch(rect)
        
        y_pos += 1
    
    ax_c.set_xlim(0, sequence_length)
    ax_c.set_ylim(-0.5, len(class_motifs) - 0.5)
    ax_c.set_xlabel('Genomic Position (bp)', fontweight='normal', fontsize=8)
    ax_c.set_ylabel('Motif Class', fontweight='normal', fontsize=8)
    ax_c.set_title('C. Genomic Coverage Track', fontweight='bold', fontsize=9, loc='left')
    
    ax_c.set_yticks(range(len(class_motifs)))
    display_labels_c = [_format_display_name(c) for c in sorted(class_motifs.keys())]
    ax_c.set_yticklabels(display_labels_c, fontsize=6)
    ax_c.ticklabel_format(style='sci', axis='x', scilimits=(3, 3))
    _apply_nature_style(ax_c)
    
    # Panel D: Score vs Length
    ax_d = fig.add_subplot(gs[2, 0])
    add_panel_label(ax_d, 'D')
    
    scatter_data = []
    for motif in motifs:
        score = motif.get('Score', motif.get('Normalized_Score'))
        length = motif.get('Length')
        class_name = motif.get('Class', 'Unknown')
        if isinstance(score, (int, float)) and isinstance(length, int) and length > 0:
            scatter_data.append({
                'Score': score,
                'Length': length,
                'Class': class_name
            })
    
    if scatter_data:
        df_scatter = pd.DataFrame(scatter_data)
        
        for class_name in df_scatter['Class'].unique()[:5]:
            class_df = df_scatter[df_scatter['Class'] == class_name]
            color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
            ax_d.scatter(class_df['Length'], class_df['Score'], 
                       c=color, s=15, alpha=0.6, edgecolors='black', linewidth=0.2)
        
        ax_d.set_xlabel('Length (bp)', fontweight='normal', fontsize=8)
        ax_d.set_ylabel('Score', fontweight='normal', fontsize=8)
        ax_d.set_title('D. Score-Length Relationship', fontweight='bold', fontsize=9, loc='left')
        ax_d.grid(alpha=0.3, linestyle='--', linewidth=0.5)
        _apply_nature_style(ax_d)
    
    # Panel E: Density heatmap
    ax_e = fig.add_subplot(gs[2, 1])
    add_panel_label(ax_e, 'E')
    
    window_size = max(1000, sequence_length // 20)
    num_windows = max(1, sequence_length // window_size)
    
    classes_viz = sorted(set(m.get('Class', 'Unknown') for m in motifs 
                            if m.get('Class') not in CIRCOS_EXCLUDED_CLASSES))[:MAX_CLASSES_FOR_SCATTER]
    
    density_matrix = np.zeros((len(classes_viz), num_windows))
    
    for i, class_name in enumerate(classes_viz):
        class_motifs_e = [m for m in motifs if m.get('Class') == class_name]
        
        for j in range(num_windows):
            window_start = j * window_size
            window_end = (j + 1) * window_size
            
            count = sum(1 for m in class_motifs_e 
                       if not (m.get('End', 0) <= window_start or 
                              m.get('Start', 0) - 1 >= window_end))
            density_matrix[i, j] = count
    
    im = ax_e.imshow(density_matrix, cmap='plasma', aspect='auto', interpolation='nearest')
    
    ax_e.set_yticks(range(len(classes_viz)))
    display_classes_e = [_format_display_name(c)[:10] for c in classes_viz]
    ax_e.set_yticklabels(display_classes_e, fontsize=6)
    
    x_ticks = np.arange(0, num_windows, max(1, num_windows // 4))
    ax_e.set_xticks(x_ticks)
    ax_e.set_xticklabels([f'{int(i*window_size/1000)}' for i in x_ticks], fontsize=6)
    
    ax_e.set_xlabel('Position (kb)', fontweight='normal', fontsize=8)
    ax_e.set_ylabel('Class', fontweight='normal', fontsize=8)
    ax_e.set_title('E. Density Heatmap', fontweight='bold', fontsize=9, loc='left')
    
    cbar_e = plt.colorbar(im, ax=ax_e, shrink=0.6)
    cbar_e.set_label('Count', fontsize=7)
    cbar_e.ax.tick_params(labelsize=6)
    
    fig.suptitle(title, fontsize=12, fontweight='bold', y=0.99)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    
    if save_path:
        fig.savefig(save_path, dpi=PUBLICATION_DPI, bbox_inches='tight',
                   facecolor='white', format='pdf')
        print(f"✓ Saved Supplementary Figure S1 to {save_path}")
    
    return fig


# =============================================================================
# EXAMPLE USAGE AND TESTING
# =============================================================================

def generate_example_figures(motifs: List[Dict[str, Any]], 
                            sequence_length: int,
                            output_dir: str = "/tmp/publication_figures",
                            include_supplementary: bool = True):
    """
    Generate all publication-ready multi-panel figures.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of analyzed sequence
        output_dir: Directory to save figures
        include_supplementary: Whether to include supplementary figures
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    print("\n" + "=" * 80)
    print("GENERATING PUBLICATION-READY MULTI-PANEL FIGURES")
    print("=" * 80)
    
    figures_generated = []
    
    # Main Figures
    print("\n[MAIN FIGURES]")
    
    # Figure 1: Overview
    print("\n[1/6] Generating Figure 1: Comprehensive Overview...")
    try:
        fig1 = create_figure_1_overview(motifs, sequence_length,
                                        save_path=f"{output_dir}/figure1_overview.pdf")
        plt.close(fig1)
        figures_generated.append("Figure 1")
        print("      ✓ Saved Figure 1")
    except Exception as e:
        print(f"      ✗ Error: {e}")
    
    # Figure 2: Statistics
    print("\n[2/6] Generating Figure 2: Statistical Characterization...")
    try:
        fig2 = create_figure_2_statistics(motifs,
                                          save_path=f"{output_dir}/figure2_statistics.pdf")
        plt.close(fig2)
        figures_generated.append("Figure 2")
        print("      ✓ Saved Figure 2")
    except Exception as e:
        print(f"      ✗ Error: {e}")
    
    # Figure 3: Landscape
    print("\n[3/6] Generating Figure 3: Genomic Landscape...")
    try:
        fig3 = create_figure_3_landscape(motifs, sequence_length,
                                         save_path=f"{output_dir}/figure3_landscape.pdf")
        plt.close(fig3)
        figures_generated.append("Figure 3")
        print("      ✓ Saved Figure 3")
    except Exception as e:
        print(f"      ✗ Error: {e}")
    
    # Figure 4: Subclass Analysis
    print("\n[4/6] Generating Figure 4: Subclass Analysis...")
    try:
        fig4 = create_figure_4_subclass_analysis(motifs, sequence_length,
                                                save_path=f"{output_dir}/figure4_subclass.pdf")
        plt.close(fig4)
        figures_generated.append("Figure 4")
        print("      ✓ Saved Figure 4")
    except Exception as e:
        print(f"      ✗ Error: {e}")
    
    # Figure 5: Quality Metrics
    print("\n[5/6] Generating Figure 5: Quality Metrics...")
    try:
        fig5 = create_figure_5_quality_metrics(motifs, sequence_length,
                                              save_path=f"{output_dir}/figure5_quality.pdf")
        plt.close(fig5)
        figures_generated.append("Figure 5")
        print("      ✓ Saved Figure 5")
    except Exception as e:
        print(f"      ✗ Error: {e}")
    
    # Figure 6: Comparative Analysis
    print("\n[6/6] Generating Figure 6: Comparative Analysis...")
    try:
        fig6 = create_figure_6_comparative_analysis(motifs, sequence_length,
                                                   save_path=f"{output_dir}/figure6_comparative.pdf")
        plt.close(fig6)
        figures_generated.append("Figure 6")
        print("      ✓ Saved Figure 6")
    except Exception as e:
        print(f"      ✗ Error: {e}")
    
    # Supplementary Figures
    if include_supplementary:
        print("\n[SUPPLEMENTARY FIGURES]")
        
        print("\n[S1] Generating Supplementary Figure S1: Comprehensive Overview...")
        try:
            figS1 = create_supplementary_figure_s1(motifs, sequence_length,
                                                  save_path=f"{output_dir}/figureS1_comprehensive.pdf")
            plt.close(figS1)
            figures_generated.append("Figure S1")
            print("      ✓ Saved Supplementary Figure S1")
        except Exception as e:
            print(f"      ✗ Error: {e}")
    
    print("\n" + "=" * 80)
    print(f"✅ Generated {len(figures_generated)} figures:")
    for fig_name in figures_generated:
        print(f"   • {fig_name}")
    print(f"\n✅ All figures saved to: {output_dir}")
    print("=" * 80)
    
    return figures_generated


if __name__ == "__main__":
    print("Publication Figures Module - Nature-Level Multi-Panel Figures")
    print("Import this module to create publication-ready figures.")
