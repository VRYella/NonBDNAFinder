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
    ax_c = fig.add_subplot(gs[1, :])
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
    ax_c.set_xlabel('Genomic Position (bp)', fontweight='normal')
    ax_c.set_ylabel('Motif Class', fontweight='normal')
    ax_c.set_title('Genomic Coverage Map', fontweight='bold', fontsize=9)
    
    ax_c.set_yticks(list(class_positions.values()))
    display_labels = [_format_display_name(label) for label in class_positions.keys()]
    ax_c.set_yticklabels(display_labels, fontsize=7)
    ax_c.ticklabel_format(style='sci', axis='x', scilimits=(3, 3))
    
    _apply_nature_style(ax_c)
    
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
# EXAMPLE USAGE AND TESTING
# =============================================================================

def generate_example_figures(motifs: List[Dict[str, Any]], 
                            sequence_length: int,
                            output_dir: str = "/tmp/publication_figures"):
    """
    Generate all publication-ready multi-panel figures.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of analyzed sequence
        output_dir: Directory to save figures
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    print("\n" + "=" * 80)
    print("GENERATING PUBLICATION-READY MULTI-PANEL FIGURES")
    print("=" * 80)
    
    # Figure 1: Overview
    print("\nGenerating Figure 1: Comprehensive Overview...")
    fig1 = create_figure_1_overview(motifs, sequence_length,
                                    save_path=f"{output_dir}/figure1_overview.pdf")
    
    # Figure 2: Statistics
    print("\nGenerating Figure 2: Statistical Characterization...")
    fig2 = create_figure_2_statistics(motifs,
                                      save_path=f"{output_dir}/figure2_statistics.pdf")
    
    # Figure 3: Landscape
    print("\nGenerating Figure 3: Genomic Landscape...")
    fig3 = create_figure_3_landscape(motifs, sequence_length,
                                     save_path=f"{output_dir}/figure3_landscape.pdf")
    
    # Figure 4: Subclass Analysis
    print("\nGenerating Figure 4: Subclass Analysis...")
    fig4 = create_figure_4_subclass_analysis(motifs, sequence_length,
                                            save_path=f"{output_dir}/figure4_subclass.pdf")
    
    print("\n" + "=" * 80)
    print(f"✅ All figures saved to: {output_dir}")
    print("=" * 80)


if __name__ == "__main__":
    print("Publication Figures Module - Nature-Level Multi-Panel Figures")
    print("Import this module to create publication-ready figures.")
