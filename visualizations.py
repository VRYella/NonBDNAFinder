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

# Try to import plotly for interactive plots
try:
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

# =============================================================================
# NATURE-LEVEL PUBLICATION STYLING & CONFIGURATION
# =============================================================================
# Following Nature Methods/Nature Genetics style guidelines:
# - Sans-serif fonts (Arial/Helvetica preferred)
# - High DPI (300+ for print)
# - Clean, minimal design with proper spacing
# - Colorblind-friendly palettes

# Default DPI for publication quality (Nature requires 300 DPI minimum)
PUBLICATION_DPI = 300

# Nature-level color palette for motif classes (colorblind-friendly)
# Based on Wong, B. (2011) Nature Methods colorblind-safe palette
# Each motif class has a unique color for distinguishability
MOTIF_CLASS_COLORS = {
    'Curved_DNA': '#CC79A7',          # Reddish Purple
    'Slipped_DNA': '#E69F00',         # Orange
    'Cruciform': '#56B4E9',           # Sky Blue
    'R-Loop': '#009E73',              # Bluish Green
    'Triplex': '#F0E442',             # Yellow
    'G-Quadruplex': '#0072B2',        # Blue
    'i-Motif': '#D55E00',             # Vermillion
    'Z-DNA': '#882255',               # Wine (distinct from Curved_DNA)
    'A-philic_DNA': '#44AA99',        # Teal (distinct from Cruciform)
    'Hybrid': '#999999',              # Gray - neutral
    'Non-B_DNA_Clusters': '#666666'   # Dark Gray - professional
}

# Nature-level scientific styling configuration for publication-quality plots
# Reference: Nature author guidelines for figure preparation
_NATURE_STYLE_PARAMS = {
    # Typography - Arial/Helvetica as per Nature guidelines
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,  # Nature recommends 5-7pt, we use 8pt for readability
    
    # Title and labels
    'axes.titlesize': 9,
    'axes.titleweight': 'bold',
    'axes.labelsize': 8,
    'axes.labelweight': 'normal',
    
    # Tick labels
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    
    # Legend
    'legend.fontsize': 7,
    'legend.frameon': False,
    'legend.borderpad': 0.4,
    
    # Figure
    'figure.titlesize': 10,
    'figure.titleweight': 'bold',
    'figure.dpi': PUBLICATION_DPI,
    'figure.facecolor': 'white',
    'savefig.dpi': PUBLICATION_DPI,
    'savefig.format': 'pdf',  # Vector format for publication
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    
    # Axes - minimal, clean design
    'axes.grid': False,  # Nature typically uses minimal grids
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.spines.left': True,
    'axes.spines.bottom': True,
    'axes.linewidth': 0.8,
    'axes.edgecolor': 'black',
    'axes.facecolor': 'white',
    
    # Lines
    'lines.linewidth': 1.0,
    'lines.markersize': 4,
    
    # PDF/SVG output
    'pdf.fonttype': 42,  # TrueType for editability
    'ps.fonttype': 42,
}

# Apply Nature-level styling at module load
plt.rcParams.update(_NATURE_STYLE_PARAMS)

# Constants for enrichment analysis visualization
INFINITE_FOLD_ENRICHMENT_CAP = 100  # Cap for infinite fold enrichment values in plots

# Figure size constants (in inches) following Nature guidelines
# Nature column widths: single=89mm(3.5in), 1.5=120mm(4.7in), double=183mm(7.2in)
FIGURE_SIZES = {
    'single_column': (3.5, 2.8),      # Single column
    'one_and_half': (4.7, 3.5),       # 1.5 column
    'double_column': (7.2, 4.5),      # Double column
    'square': (3.5, 3.5),             # Square figure
    'wide': (7.2, 3.0),               # Wide panoramic
}


def set_scientific_style(style: str = 'nature'):
    """
    Apply scientific publication-ready styling.
    
    Args:
        style: Style preset ('nature', 'default', 'presentation')
               - 'nature': Nature/Science journal style (default)
               - 'default': Standard seaborn whitegrid
               - 'presentation': Larger fonts for presentations
    """
    if style == 'nature':
        plt.rcParams.update(_NATURE_STYLE_PARAMS)
        sns.set_style("white")
        sns.set_context("paper")
    elif style == 'presentation':
        plt.rcParams.update(_NATURE_STYLE_PARAMS)
        plt.rcParams.update({
            'font.size': 12,
            'axes.titlesize': 14,
            'axes.labelsize': 12,
            'xtick.labelsize': 11,
            'ytick.labelsize': 11,
            'legend.fontsize': 10,
        })
        sns.set_style("whitegrid")
        sns.set_context("talk")
    else:
        sns.set_style("whitegrid")
        sns.set_palette("husl")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def _apply_nature_style(ax):
    """Apply Nature journal style to matplotlib axis (remove top/right spines).
    
    | Parameter | Type | Description            |
    |-----------|------|------------------------|
    | ax        | Axes | Matplotlib axis object |
    """
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# =============================================================================
# DISTRIBUTION PLOTS (Nature-Level Quality)
# =============================================================================

def plot_motif_distribution(motifs: List[Dict[str, Any]], 
                           by: str = 'Class',
                           title: Optional[str] = None,
                           figsize: Tuple[float, float] = None,
                           style: str = 'nature') -> plt.Figure:
    """
    Plot distribution of motifs by class or subclass.
    
    Publication-quality bar chart following Nature Methods guidelines.
    Shows all classes/subclasses even when count is 0, ensuring comprehensive
    visualization of both detected and undetected motif types.
    
    Args:
        motifs: List of motif dictionaries.
        by: Group by 'Class' or 'Subclass'.
        title: Custom plot title.
        figsize: Figure size (width, height) in inches. Uses Nature standard if None.
        style: Style preset ('nature', 'default', 'presentation').
        
    Returns:
        Matplotlib figure object (publication-ready at 300 DPI)
    """
    set_scientific_style(style)
    
    # Use Nature-appropriate figure size
    if figsize is None:
        figsize = FIGURE_SIZES['double_column'] if by == 'Subclass' else FIGURE_SIZES['one_and_half']
    
    # Define ALL expected classes and subclasses (always show these)
    ALL_CLASSES = [
        'Curved_DNA', 'Slipped_DNA', 'Cruciform', 'R-Loop', 'Triplex',
        'G-Quadruplex', 'i-Motif', 'Z-DNA', 'A-philic_DNA', 'Hybrid', 'Non-B_DNA_Clusters'
    ]
    
    ALL_SUBCLASSES = [
        'Global Curvature', 'Local Curvature',  # Curved DNA
        'Direct Repeat', 'STR',  # Slipped DNA
        'Inverted Repeats',  # Cruciform
        'R-loop formation sites', 'QmRLFS-m1', 'QmRLFS-m2',  # R-Loop
        'Triplex', 'Sticky DNA',  # Triplex
        'Canonical G4', 'Relaxed G4', 'Bulged G4', 'Bipartite G4',  # G-Quadruplex
        'Multimeric G4', 'Imperfect G4', 'G-Triplex intermediate',
        'Canonical i-motif', 'Relaxed i-motif', 'AC-motif',  # i-Motif
        'Z-DNA', 'eGZ (Extruded-G) DNA',  # Z-DNA
        'A-philic DNA',  # A-philic
    ]
    
    # Count motifs by specified grouping
    counts = Counter(m.get(by, 'Unknown') for m in motifs) if motifs else Counter()
    
    # Prepare data with all categories
    if by == 'Class':
        categories = ALL_CLASSES
    else:
        categories = ALL_SUBCLASSES
    
    # Get counts (0 if not present)
    values = [counts.get(cat, 0) for cat in categories]
    
    # Get colors (colorblind-friendly)
    if by == 'Class':
        colors = [MOTIF_CLASS_COLORS.get(cat, '#808080') for cat in categories]
    else:
        # Map subclasses to their parent class colors
        subclass_to_class = {
            'Global Curvature': 'Curved_DNA', 'Local Curvature': 'Curved_DNA',
            'Direct Repeat': 'Slipped_DNA', 'STR': 'Slipped_DNA',
            'Inverted Repeats': 'Cruciform',
            'R-loop formation sites': 'R-Loop', 'QmRLFS-m1': 'R-Loop', 'QmRLFS-m2': 'R-Loop',
            'Triplex': 'Triplex', 'Sticky DNA': 'Triplex',
            'Canonical G4': 'G-Quadruplex', 'Relaxed G4': 'G-Quadruplex',
            'Bulged G4': 'G-Quadruplex', 'Bipartite G4': 'G-Quadruplex',
            'Multimeric G4': 'G-Quadruplex', 'Imperfect G4': 'G-Quadruplex',
            'G-Triplex intermediate': 'G-Quadruplex',
            'Canonical i-motif': 'i-Motif', 'Relaxed i-motif': 'i-Motif', 'AC-motif': 'i-Motif',
            'Z-DNA': 'Z-DNA', 'eGZ (Extruded-G) DNA': 'Z-DNA',
            'A-philic DNA': 'A-philic_DNA',
        }
        colors = [MOTIF_CLASS_COLORS.get(subclass_to_class.get(cat, ''), '#808080') for cat in categories]
    
    # Create figure with high DPI
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    # Create bars with Nature-style aesthetics
    bars = ax.bar(range(len(categories)), values, color=colors, 
                  edgecolor='black', linewidth=0.5, width=0.8)
    
    # Customize axes (Nature style - minimal, clean)
    ax.set_xlabel(f'Motif {by}', fontweight='normal')
    ax.set_ylabel('Count', fontweight='normal')
    if title:
        # Replace underscores with spaces in title
        display_title = title.replace('_', ' ')
        ax.set_title(display_title, fontweight='bold', pad=10)
    ax.set_xticks(range(len(categories)))
    
    # Replace underscores with spaces in category labels
    # Apply 45° rotation for all categories (Nature style - consistent readability)
    display_categories = [cat.replace('_', ' ') for cat in categories]
    ax.set_xticklabels(display_categories, rotation=45, ha='right')
    
    # Add count labels on bars (small font, positioned above)
    if len(categories) <= 20:
        max_val = max(values) if max(values) > 0 else 1
        for bar, count in zip(bars, values):
            if count > 0:  # Only label non-zero bars
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + max_val * 0.02,
                        str(count), ha='center', va='bottom', fontsize=6)
    
    # Apply Nature journal style
    _apply_nature_style(ax)
    
    plt.tight_layout()
    return fig


def plot_class_subclass_sunburst(motifs: List[Dict[str, Any]], 
                                 title: str = "Motif Class-Subclass Distribution",
                                 figsize: Tuple[float, float] = None) -> Union[plt.Figure, Any]:
    """
    Create sunburst plot showing class-subclass hierarchy.
    
    Publication-quality hierarchical visualization.
    
    Args:
        motifs: List of motif dictionaries
        title: Plot title
        figsize: Figure size (width, height) in inches
        
    Returns:
        Plotly figure if available, otherwise matplotlib figure
    """
    if figsize is None:
        figsize = FIGURE_SIZES['square']
        
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center', 
                transform=ax.transAxes)
        ax.axis('off')
        if title:
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
        marker=dict(colors=colors, line=dict(color="#FFFFFF", width=1)),
        hovertemplate='<b>%{label}</b><br>Count: %{value}<extra></extra>',
    ))
    
    # Publication-quality layout
    fig.update_layout(
        title=dict(text=title, font=dict(size=10, family='Arial')),
        font=dict(size=8, family='Arial'),
        width=int(figsize[0] * 100),
        height=int(figsize[1] * 100),
        margin=dict(t=30, l=10, r=10, b=10)
    )
    
    return fig


def plot_nested_pie_chart(motifs: List[Dict[str, Any]], 
                         title: str = "Motif Distribution",
                         figsize: Tuple[float, float] = None) -> plt.Figure:
    """
    Create nested donut chart with improved text placement to avoid overlapping labels.
    
    Publication-quality hierarchical pie chart following Nature guidelines.
    
    Args:
        motifs: List of motif dictionaries
        title: Plot title
        figsize: Figure size (width, height) in inches
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style('nature')
    
    if figsize is None:
        figsize = FIGURE_SIZES['square']
    
    # Count by class and subclass
    class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    class_subclass_counts = defaultdict(lambda: defaultdict(int))
    
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        subclass_name = motif.get('Subclass', 'Unknown')
        class_subclass_counts[class_name][subclass_name] += 1
    
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    # Inner donut (classes)
    class_names = list(class_counts.keys())
    class_values = list(class_counts.values())
    class_colors = [MOTIF_CLASS_COLORS.get(name, '#808080') for name in class_names]
    
    # Create inner donut with Nature-style clean design
    # Use labels=None to manually place labels later for better control
    wedges1, texts1, autotexts1 = ax.pie(
        class_values, 
        labels=None,  # We'll add labels manually
        colors=class_colors,
        radius=0.65,
        autopct=lambda pct: f'{pct:.0f}%' if pct > 5 else '',
        pctdistance=0.80,
        startangle=90,
        wedgeprops=dict(width=0.35, edgecolor='white', linewidth=1)
    )
    
    # Manually add class labels with better positioning (replace underscores with spaces)
    for i, (wedge, class_name) in enumerate(zip(wedges1, class_names)):
        angle = (wedge.theta2 + wedge.theta1) / 2
        x = 0.5 * np.cos(np.radians(angle))
        y = 0.5 * np.sin(np.radians(angle))
        # Replace underscores with spaces in labels
        display_name = class_name.replace('_', ' ')
        ax.text(x, y, display_name, ha='center', va='center', fontsize=7, fontweight='bold')
    
    # Outer donut (subclasses)
    all_subclass_counts = []
    all_subclass_colors = []
    all_subclass_labels = []
    
    for class_name in class_names:
        subclass_dict = class_subclass_counts[class_name]
        base_color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        
        for subclass_name, count in subclass_dict.items():
            all_subclass_counts.append(count)
            # Truncate long names for clean appearance and replace underscores with spaces
            # Use consistent truncation length (15 chars max, including ellipsis)
            display_name = subclass_name.replace('_', ' ')
            MAX_LABEL_LENGTH = 15
            label = display_name if len(display_name) <= MAX_LABEL_LENGTH else display_name[:MAX_LABEL_LENGTH-1] + '…'
            all_subclass_labels.append(label)
            all_subclass_colors.append(base_color)
    
    # Use smarter labeling strategy to avoid overlap
    # For many subclasses, hide labels and rely on legend instead
    if len(all_subclass_labels) > 10:
        # Hide outer ring labels when there are too many
        wedges2, texts2 = ax.pie(
            all_subclass_counts,
            labels=None,  # No labels for cleaner appearance
            colors=all_subclass_colors,
            radius=1.0,
            startangle=90,
            wedgeprops=dict(width=0.35, edgecolor='white', linewidth=0.8)
        )
        
        # Add a legend for subclasses instead
        # Note: All subclasses of the same parent class share the same color by design.
        # This ensures visual grouping in the nested donut chart.
        # Track first occurrence of each unique label for the legend.
        seen_labels = {}
        legend_handles = []
        legend_labels = []
        for i, (label, color) in enumerate(zip(all_subclass_labels, all_subclass_colors)):
            if label not in seen_labels:
                seen_labels[label] = color
                legend_handles.append(plt.Rectangle((0,0),1,1, fc=color, ec='white', lw=0.5))
                legend_labels.append(label)
                if len(legend_labels) >= 10:  # Limit to 10
                    break
        
        ax.legend(legend_handles, legend_labels, loc='center left', bbox_to_anchor=(1, 0.5),
                 fontsize=6, frameon=False, title='Top Subclasses')
    else:
        # For fewer subclasses, show labels with improved spacing
        wedges2, texts2 = ax.pie(
            all_subclass_counts,
            labels=all_subclass_labels,
            colors=all_subclass_colors,
            radius=1.0,
            labeldistance=1.15,  # Push labels further out to avoid overlap
            startangle=90,
            wedgeprops=dict(width=0.35, edgecolor='white', linewidth=0.8),
            textprops={'fontsize': 6}
        )
        
        # Adjust label positions to avoid overlap
        for text in texts2:
            text.set_fontsize(6)
            # Add slight rotation for better readability
            angle = text.get_rotation()
            if 90 < angle < 270:
                text.set_rotation(angle - 180)
    
    if title:
        # Replace underscores with spaces in title
        display_title = title.replace('_', ' ')
        ax.set_title(display_title, fontweight='bold', pad=10)
    
    # Style percentage labels
    for autotext in autotexts1:
        autotext.set_fontsize(6)
        autotext.set_color('white')
    
    return fig


# =============================================================================
# COVERAGE & POSITIONAL PLOTS (Nature-Level Quality)
# =============================================================================

def plot_coverage_map(motifs: List[Dict[str, Any]], 
                     sequence_length: int,
                     title: Optional[str] = None,
                     figsize: Tuple[float, float] = None,
                     style: str = 'nature') -> plt.Figure:
    """
    Plot motif coverage map showing positions along sequence.
    
    Publication-quality genomic track visualization.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of the analyzed sequence
        title: Custom plot title
        figsize: Figure size (width, height) in inches
        style: Style preset ('nature', 'default', 'presentation')
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style(style)
    
    if figsize is None:
        figsize = FIGURE_SIZES['wide']
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center', 
                transform=ax.transAxes)
        ax.axis('off')
        if title:
            ax.set_title(title)
        return fig
    
    # Group motifs by class
    class_motifs = defaultdict(list)
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        class_motifs[class_name].append(motif)
    
    # Create figure with high DPI
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    y_pos = 0
    class_positions = {}
    
    for class_name, class_motif_list in class_motifs.items():
        class_positions[class_name] = y_pos
        color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        
        for motif in class_motif_list:
            start = motif.get('Start', 0) - 1  # Convert to 0-based
            end = motif.get('End', start + 1)
            length = end - start
            
            # Draw motif as rectangle (Nature style - clean edges)
            rect = patches.Rectangle(
                (start, y_pos - 0.35), length, 0.7,
                facecolor=color, edgecolor='black', linewidth=0.3
            )
            ax.add_patch(rect)
        
        y_pos += 1
    
    # Customize axes (Nature style)
    ax.set_xlim(0, sequence_length)
    ax.set_ylim(-0.5, len(class_motifs) - 0.5)
    ax.set_xlabel('Position (bp)')
    ax.set_ylabel('Motif Class')
    if title:
        # Replace underscores with spaces in title
        display_title = title.replace('_', ' ')
        ax.set_title(display_title, fontweight='bold', pad=10)
    
    # Set y-axis labels with underscores replaced by spaces
    ax.set_yticks(list(class_positions.values()))
    display_labels = [label.replace('_', ' ') for label in class_positions.keys()]
    ax.set_yticklabels(display_labels)
    
    # Clean x-axis ticks
    ax.ticklabel_format(style='sci', axis='x', scilimits=(3, 3))
    
    # Apply Nature journal style
    _apply_nature_style(ax)
    
    plt.tight_layout()
    return fig


def plot_density_heatmap(motifs: List[Dict[str, Any]], 
                        sequence_length: int,
                        window_size: int = 1000,
                        title: Optional[str] = None,
                        figsize: Tuple[float, float] = None,
                        style: str = 'nature') -> plt.Figure:
    """
    Plot motif density heatmap along sequence.
    
    Publication-quality density visualization.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of the analyzed sequence
        window_size: Window size for density calculation
        title: Custom plot title
        figsize: Figure size (width, height) in inches
        style: Style preset
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style(style)
    
    if figsize is None:
        figsize = FIGURE_SIZES['wide']
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center', 
                transform=ax.transAxes)
        ax.axis('off')
        if title:
            ax.set_title(title)
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
    
    # Create heatmap with publication quality
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    # Use colorblind-friendly colormap
    im = ax.imshow(density_matrix, cmap='viridis', aspect='auto', interpolation='nearest')
    
    # Customize axes (Nature style)
    ax.set_xlabel(f'Position (kb)')
    ax.set_ylabel('Motif Class')
    if title:
        # Replace underscores with spaces in title
        display_title = title.replace('_', ' ')
        ax.set_title(display_title, fontweight='bold', pad=10)
    
    # Set ticks and labels with underscores replaced by spaces
    ax.set_yticks(range(len(classes)))
    display_classes = [cls.replace('_', ' ') for cls in classes]
    ax.set_yticklabels(display_classes)
    
    # Clean x-axis with kb units
    x_ticks = np.arange(0, num_windows, max(1, num_windows // 5))
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f'{int(windows[i]/1000)}' for i in x_ticks])
    
    # Add colorbar with proper label
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Count', fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    
    plt.tight_layout()
    return fig


# =============================================================================
# STATISTICAL PLOTS (Nature-Level Quality)
# =============================================================================

def plot_score_distribution(motifs: List[Dict[str, Any]], 
                           by_class: bool = True,
                           title: Optional[str] = None,
                           figsize: Tuple[float, float] = None,
                           style: str = 'nature') -> plt.Figure:
    """
    Plot distribution of motif scores.
    
    Publication-quality box plot visualization.
    
    Args:
        motifs: List of motif dictionaries
        by_class: Whether to separate by motif class
        title: Custom plot title
        figsize: Figure size (width, height) in inches
        style: Style preset
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style(style)
    
    if figsize is None:
        figsize = FIGURE_SIZES['one_and_half']
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center', 
                transform=ax.transAxes)
        ax.axis('off')
        if title:
            ax.set_title(title)
        return fig
    
    # Extract scores
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
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No score data available', ha='center', va='center', 
                transform=ax.transAxes)
        ax.axis('off')
        if title:
            display_title = title.replace('_', ' ')
            ax.set_title(display_title)
        return fig
    
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    if by_class and isinstance(scores_data[0], dict):
        # Create DataFrame for seaborn
        df = pd.DataFrame(scores_data)
        
        # Box plot by class (Nature style - clean, minimal)
        colors = [MOTIF_CLASS_COLORS.get(cls, '#808080') for cls in df['Class'].unique()]
        sns.boxplot(data=df, x='Class', y='Score', ax=ax, palette=colors,
                   linewidth=0.8, fliersize=2)
        # Replace underscores with spaces in x-tick labels
        ax.set_xticklabels([label.get_text().replace('_', ' ') for label in ax.get_xticklabels()], 
                          rotation=45, ha='right')
        ax.set_ylabel('Score')
        ax.set_xlabel('Motif Class')
    else:
        # Simple histogram
        ax.hist(scores_data, bins=20, edgecolor='black', linewidth=0.5, color='#56B4E9')
        ax.set_xlabel('Score')
        ax.set_ylabel('Frequency')
    
    if title:
        # Replace underscores with spaces in title
        display_title = title.replace('_', ' ')
        ax.set_title(display_title, fontweight='bold', pad=10)
    
    # Apply Nature journal style
    _apply_nature_style(ax)
    
    plt.tight_layout()
    return fig


def plot_length_distribution(motifs: List[Dict[str, Any]], 
                           by_class: bool = True,
                           title: Optional[str] = None,
                           figsize: Tuple[float, float] = None,
                           style: str = 'nature') -> plt.Figure:
    """
    Plot distribution of motif lengths.
    
    Publication-quality violin plot visualization.
    
    Args:
        motifs: List of motif dictionaries
        by_class: Whether to separate by motif class
        title: Custom plot title
        figsize: Figure size (width, height) in inches
        style: Style preset
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style(style)
    
    if figsize is None:
        figsize = FIGURE_SIZES['one_and_half']
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center', 
                transform=ax.transAxes)
        ax.axis('off')
        if title:
            ax.set_title(title)
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
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No length data available', ha='center', va='center', 
                transform=ax.transAxes)
        ax.axis('off')
        if title:
            display_title = title.replace('_', ' ')
            ax.set_title(display_title)
        return fig
    
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    if by_class and isinstance(length_data[0], dict):
        # Create DataFrame for seaborn
        df = pd.DataFrame(length_data)
        
        # Violin plot by class (Nature style - clean)
        colors = [MOTIF_CLASS_COLORS.get(cls, '#808080') for cls in df['Class'].unique()]
        sns.violinplot(data=df, x='Class', y='Length', ax=ax, palette=colors,
                      linewidth=0.8, inner='box')
        # Replace underscores with spaces in x-tick labels
        ax.set_xticklabels([label.get_text().replace('_', ' ') for label in ax.get_xticklabels()], 
                          rotation=45, ha='right')
        ax.set_ylabel('Length (bp)')
        ax.set_xlabel('Motif Class')
    else:
        # Simple histogram
        ax.hist(length_data, bins=20, edgecolor='black', linewidth=0.5, color='#009E73')
        ax.set_xlabel('Length (bp)')
        ax.set_ylabel('Frequency')
    
    if title:
        # Replace underscores with spaces in title
        display_title = title.replace('_', ' ')
        ax.set_title(display_title, fontweight='bold', pad=10)
    
    # Apply Nature journal style
    _apply_nature_style(ax)
    
    plt.tight_layout()
    return fig


# =============================================================================
# COMPARISON PLOTS (Nature-Level Quality)
# =============================================================================

def plot_class_comparison(results: Dict[str, List[Dict[str, Any]]], 
                         metric: str = 'count',
                         title: Optional[str] = None,
                         figsize: Tuple[float, float] = None,
                         style: str = 'nature') -> plt.Figure:
    """
    Compare motif classes across multiple sequences/samples.
    
    Publication-quality heatmap comparison.
    
    Args:
        results: Dictionary of {sample_name: motifs_list}
        metric: Comparison metric ('count', 'coverage', 'density')
        title: Custom plot title
        figsize: Figure size (width, height) in inches
        style: Style preset
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style(style)
    
    if figsize is None:
        figsize = FIGURE_SIZES['double_column']
    
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
# BATCH EXPORT FUNCTIONS (Publication-Ready)
# =============================================================================

def save_all_plots(motifs: List[Dict[str, Any]], 
                   sequence_length: int,
                   output_dir: str = "plots",
                   file_format: str = "pdf",
                   dpi: int = 300,
                   style: str = 'nature') -> Dict[str, str]:
    """
    Generate and save all standard plots in publication quality.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of analyzed sequence
        output_dir: Output directory for plots
        file_format: File format ('pdf', 'png', 'svg') - pdf recommended for publication
        dpi: Resolution for raster formats (default 300 DPI for publication)
        style: Style preset ('nature', 'default')
        
    Returns:
        Dictionary of {plot_name: file_path}
    """
    import os
    
    # Apply style
    set_scientific_style(style)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    saved_files = {}
    
    # Use publication DPI if not specified
    if dpi < PUBLICATION_DPI:
        dpi = PUBLICATION_DPI
    
    # List of plots to generate (diverse, non-repetitive visualization types)
    plots_to_generate = [
        ("motif_distribution_class", lambda: plot_motif_distribution(motifs, by='Class', style=style)),
        ("coverage_map", lambda: plot_coverage_map(motifs, sequence_length, style=style)),
        ("density_heatmap", lambda: plot_density_heatmap(motifs, sequence_length, style=style)),
        ("score_distribution", lambda: plot_score_distribution(motifs, by_class=True, style=style)),
        ("length_distribution", lambda: plot_length_distribution(motifs, by_class=True, style=style)),
        ("nested_donut_chart", lambda: plot_nested_pie_chart(motifs))
    ]
    
    for plot_name, plot_func in plots_to_generate:
        try:
            fig = plot_func()
            filename = f"{plot_name}.{file_format}"
            filepath = os.path.join(output_dir, filename)
            
            if hasattr(fig, 'savefig'):  # Matplotlib figure
                # Publication-quality save options
                fig.savefig(filepath, format=file_format, dpi=dpi, 
                           bbox_inches='tight', pad_inches=0.05,
                           facecolor='white', edgecolor='none')
                plt.close(fig)
            else:  # Plotly figure
                if file_format.lower() == 'png':
                    fig.write_image(filepath, scale=2)  # Higher resolution
                elif file_format.lower() in ('pdf', 'svg'):
                    fig.write_image(filepath)
                elif file_format.lower() == 'html':
                    fig.write_html(filepath)
                else:
                    # Unsupported format - try default write_image
                    print(f"⚠ Unsupported format '{file_format}' for Plotly, attempting save anyway")
                    try:
                        fig.write_image(filepath)
                    except Exception as fmt_err:
                        raise ValueError(f"Unsupported file format: {file_format}. Use 'png', 'pdf', 'svg', or 'html'.") from fmt_err
            
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
    
    # Sort classes alphabetically
    classes = sorted(classes)
    
    genomic_vals = [genomic_density.get(c, 0) for c in classes]
    positional_vals = [positional_density.get(c, 0) for c in classes]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, dpi=PUBLICATION_DPI)
    
    # Genomic Density (Coverage %)
    colors1 = [MOTIF_CLASS_COLORS.get(c, '#808080') for c in classes]
    
    # Replace underscores with spaces in labels
    display_classes = [c.replace('_', ' ') for c in classes]
    
    bars1 = ax1.barh(display_classes, genomic_vals, color=colors1, alpha=0.8, 
                     edgecolor='black', linewidth=0.5)
    ax1.set_xlabel('Genomic Density (Coverage %)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Motif Class', fontsize=11, fontweight='bold')
    ax1.set_title('A. Genomic Density (σ_G)', fontsize=12, fontweight='bold', pad=10)
    
    # Add value labels with safe max calculation
    max_genomic = max(genomic_vals) if genomic_vals and max(genomic_vals) > 0 else 1
    for i, (bar, val) in enumerate(zip(bars1, genomic_vals)):
        if val > 0:
            ax1.text(val + max_genomic * 0.01, i, f'{val:.3f}%', 
                    va='center', fontsize=9, fontweight='bold')
    
    # Positional Density (Frequency)
    colors2 = [MOTIF_CLASS_COLORS.get(c, '#808080') for c in classes]
    bars2 = ax2.barh(display_classes, positional_vals, color=colors2, alpha=0.8, 
                     edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('Positional Density (motifs/kbp)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Motif Class', fontsize=11, fontweight='bold')
    ax2.set_title('B. Positional Density (λ)', fontsize=12, fontweight='bold', pad=10)
    
    # Add value labels with safe max calculation
    max_positional = max(positional_vals) if positional_vals and max(positional_vals) > 0 else 1
    for i, (bar, val) in enumerate(zip(bars2, positional_vals)):
        if val > 0:
            ax2.text(val + max_positional * 0.01, i, f'{val:.2f}', 
                    va='center', fontsize=9, fontweight='bold')
    
    # Apply Nature journal style
    _apply_nature_style(ax1)
    _apply_nature_style(ax2)
    
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


# =============================================================================
# CIRCOS PLOT FOR NON-B DNA MOTIF DENSITY
# =============================================================================

# Motif classes to exclude from Circos visualization (dynamic classes)
CIRCOS_EXCLUDED_CLASSES = ['Hybrid', 'Non-B_DNA_Clusters']


def plot_circos_motif_density(motifs: List[Dict[str, Any]], 
                               sequence_length: int,
                               title: str = "Non-B DNA Motif Density Circos Plot",
                               figsize: Tuple[int, int] = (12, 12),
                               window_size: int = None) -> plt.Figure:
    """
    Create a circular Circos-style plot showing non-B DNA motif class density.
    
    The plot shows:
    - Outer ring: Sequence position ruler
    - Inner rings: One ring per motif class showing density
    - Center: Summary statistics
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total length of analyzed sequence in bp
        title: Plot title
        figsize: Figure size (width, height)
        window_size: Window size for density calculation (auto-calculated if None)
        
    Returns:
        Matplotlib figure object with Circos-style visualization
    """
    set_scientific_style()
    
    if not motifs or sequence_length == 0:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    # Auto-calculate window size if not provided
    if window_size is None:
        window_size = max(100, sequence_length // 50)
    
    # Calculate number of windows
    num_windows = max(1, sequence_length // window_size)
    
    # Get unique classes (excluding Hybrid and Clusters for cleaner visualization)
    classes = sorted(set(m.get('Class', 'Unknown') for m in motifs 
                        if m.get('Class') not in CIRCOS_EXCLUDED_CLASSES))
    
    if not classes:
        classes = sorted(set(m.get('Class', 'Unknown') for m in motifs))
    
    # Calculate density per window per class
    class_densities = {}
    for class_name in classes:
        densities = []
        class_motifs = [m for m in motifs if m.get('Class') == class_name]
        
        for i in range(num_windows):
            window_start = i * window_size
            window_end = (i + 1) * window_size
            
            # Count motifs in this window
            count = 0
            for motif in class_motifs:
                motif_start = motif.get('Start', 0) - 1  # 0-based
                motif_end = motif.get('End', 0)
                # Check overlap with window
                if not (motif_end <= window_start or motif_start >= window_end):
                    count += 1
            
            # Convert to density (motifs per kb)
            window_kb = window_size / 1000
            densities.append(count / window_kb if window_kb > 0 else 0)
        
        class_densities[class_name] = densities
    
    # Create figure
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='polar')
    
    # Calculate angles for each window
    theta = np.linspace(0, 2 * np.pi, num_windows, endpoint=False)
    
    # Width of each bar
    width = 2 * np.pi / num_windows * 0.8
    
    # Ring configuration
    ring_width = 0.12
    inner_radius = 0.3
    
    # Plot each class as a ring
    for i, class_name in enumerate(classes):
        densities = class_densities[class_name]
        
        # Normalize densities for this class (0-1 scale for bar height)
        max_density = max(densities) if max(densities) > 0 else 1
        normalized = [d / max_density for d in densities]
        
        # Calculate ring position
        ring_bottom = inner_radius + i * ring_width
        
        # Get color for this class
        color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        
        # Plot bars for this ring
        heights = [n * ring_width * 0.9 for n in normalized]
        # Replace underscores with spaces in legend label
        display_name = class_name.replace('_', ' ')
        bars = ax.bar(theta, heights, width=width, bottom=ring_bottom,
                     color=color, alpha=0.7, edgecolor='white', linewidth=0.5,
                     label=f'{display_name} (max: {max_density:.1f}/kb)')
    
    # Add outer position ruler
    outer_radius = inner_radius + len(classes) * ring_width + 0.05
    ruler_theta = np.linspace(0, 2 * np.pi, 12, endpoint=False)
    ruler_labels = [f'{int(i * sequence_length / 12 / 1000)}kb' for i in range(12)]
    
    ax.set_xticks(ruler_theta)
    ax.set_xticklabels(ruler_labels, fontsize=9, fontweight='bold')
    
    # Remove radial labels
    ax.set_yticklabels([])
    
    # Set limits
    ax.set_ylim(0, outer_radius + 0.1)
    
    # Add legend
    ax.legend(loc='center', bbox_to_anchor=(0.5, 0.5), fontsize=8, 
             framealpha=0.9, ncol=1)
    
    # Add title (replace underscores with spaces)
    display_title = title.replace('_', ' ')
    fig.suptitle(display_title, fontsize=14, fontweight='bold', y=0.98)
    
    # Add center statistics
    total_motifs = len([m for m in motifs if m.get('Class') not in CIRCOS_EXCLUDED_CLASSES])
    center_text = f"Total: {total_motifs}\n{len(classes)} classes\n{sequence_length/1000:.1f} kb"
    ax.text(0, 0, center_text, ha='center', va='center', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    return fig


def plot_radial_class_density(motifs: List[Dict[str, Any]], 
                               sequence_length: int,
                               title: str = "Radial Motif Class Density",
                               figsize: Tuple[int, int] = (10, 10)) -> plt.Figure:
    """
    Create a radial bar chart showing motif density per class.
    
    A simpler alternative to full Circos plot, showing aggregate density
    per motif class in a radial/polar layout.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total length of analyzed sequence in bp
        title: Plot title
        figsize: Figure size
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    if not motifs or sequence_length == 0:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    # Calculate density per class (motifs per kb)
    sequence_kb = sequence_length / 1000
    class_counts = Counter(m.get('Class', 'Unknown') for m in motifs
                          if m.get('Class') not in CIRCOS_EXCLUDED_CLASSES)
    
    if not class_counts:
        class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    
    classes = list(class_counts.keys())
    densities = [class_counts[c] / sequence_kb for c in classes]
    
    # Create polar plot
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='polar')
    
    # Calculate angles
    num_classes = len(classes)
    theta = np.linspace(0, 2 * np.pi, num_classes, endpoint=False)
    width = 2 * np.pi / num_classes * 0.7
    
    # Get colors
    colors = [MOTIF_CLASS_COLORS.get(c, '#808080') for c in classes]
    
    # Plot bars
    bars = ax.bar(theta, densities, width=width, color=colors, alpha=0.8,
                 edgecolor='white', linewidth=2)
    
    # Add class labels (replace underscores with spaces)
    ax.set_xticks(theta)
    display_classes = [cls.replace('_', ' ') for cls in classes]
    ax.set_xticklabels(display_classes, fontsize=10, fontweight='bold')
    
    # Add value labels on bars
    for angle, density, bar in zip(theta, densities, bars):
        ax.text(angle, density + max(densities) * 0.05, f'{density:.1f}',
               ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # Style (replace underscores with spaces in title)
    ax.set_ylabel('Density (motifs/kb)', labelpad=30, fontsize=11, fontweight='bold')
    display_title = title.replace('_', ' ')
    ax.set_title(display_title, fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    return fig


def plot_stacked_density_track(motifs: List[Dict[str, Any]], 
                                sequence_length: int,
                                title: str = "Stacked Motif Density Track",
                                figsize: Tuple[int, int] = (14, 6),
                                window_size: int = None) -> plt.Figure:
    """
    Create a stacked area chart showing motif density along the sequence.
    
    A linear (non-circular) alternative to Circos plot showing how different
    motif classes are distributed along the sequence.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total length of analyzed sequence in bp
        title: Plot title
        figsize: Figure size
        window_size: Window size for density calculation (auto if None)
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    if not motifs or sequence_length == 0:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    # Auto-calculate window size
    if window_size is None:
        window_size = max(100, sequence_length // 100)
    
    num_windows = max(1, sequence_length // window_size)
    
    # Get unique classes
    classes = sorted(set(m.get('Class', 'Unknown') for m in motifs
                        if m.get('Class') not in CIRCOS_EXCLUDED_CLASSES))
    
    if not classes:
        classes = sorted(set(m.get('Class', 'Unknown') for m in motifs))
    
    # Calculate density per window per class
    positions = np.arange(num_windows) * window_size / 1000  # In kb
    class_densities = {}
    
    for class_name in classes:
        densities = []
        class_motifs = [m for m in motifs if m.get('Class') == class_name]
        
        for i in range(num_windows):
            window_start = i * window_size
            window_end = (i + 1) * window_size
            
            count = 0
            for motif in class_motifs:
                motif_start = motif.get('Start', 0) - 1
                motif_end = motif.get('End', 0)
                if not (motif_end <= window_start or motif_start >= window_end):
                    count += 1
            
            densities.append(count)
        
        class_densities[class_name] = densities
    
    # Create stacked area plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Stack the densities
    colors = [MOTIF_CLASS_COLORS.get(c, '#808080') for c in classes]
    
    # Create arrays for stacking (replace underscores with spaces in labels)
    density_arrays = [np.array(class_densities[c]) for c in classes]
    display_classes = [cls.replace('_', ' ') for cls in classes]
    
    ax.stackplot(positions, *density_arrays, labels=display_classes, colors=colors, alpha=0.8)
    
    # Styling (replace underscores with spaces in title)
    ax.set_xlabel('Position (kb)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Motif Count per Window', fontsize=12, fontweight='bold')
    display_title = title.replace('_', ' ')
    ax.set_title(display_title, fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9, framealpha=0.9)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    return fig


# =============================================================================
# SUBCLASS-LEVEL DENSITY AND ENRICHMENT VISUALIZATIONS
# =============================================================================

def plot_density_comparison_by_subclass(genomic_density: Dict[str, float],
                                        positional_density: Dict[str, float],
                                        title: str = "Motif Density Analysis (by Subclass)",
                                        figsize: Tuple[int, int] = (16, 10)) -> plt.Figure:
    """
    Plot comparison of genomic density (coverage %) and positional density (motifs/kbp)
    at the subclass level.
    
    Args:
        genomic_density: Dictionary of 'Class:Subclass' -> genomic density (%)
        positional_density: Dictionary of 'Class:Subclass' -> positional density (motifs/unit)
        title: Plot title
        figsize: Figure size
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    # Remove 'Overall' for subclass-specific comparison
    subclasses = [k for k in genomic_density.keys() if k != 'Overall' and ':' in k]
    if not subclasses:
        # Fallback to regular class-level if no subclass data
        subclasses = [k for k in genomic_density.keys() if k != 'Overall']
    
    # Sort by class, then by subclass
    subclasses.sort()
    
    genomic_vals = [genomic_density.get(c, 0) for c in subclasses]
    positional_vals = [positional_density.get(c, 0) for c in subclasses]
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, dpi=PUBLICATION_DPI)
    
    # Genomic Density (Coverage %)
    y_pos = np.arange(len(subclasses))
    
    # Get colors based on parent class
    colors1 = []
    for subclass in subclasses:
        if ':' in subclass:
            parent_class = subclass.split(':')[0]
            colors1.append(MOTIF_CLASS_COLORS.get(parent_class, '#808080'))
        else:
            colors1.append(MOTIF_CLASS_COLORS.get(subclass, '#808080'))
    
    bars1 = ax1.barh(y_pos, genomic_vals, color=colors1, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax1.set_xlabel('Genomic Density (Coverage %)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Motif Subclass', fontsize=11, fontweight='bold')
    ax1.set_title('A. Genomic Density (σ_G) by Subclass', fontsize=12, fontweight='bold', pad=10)
    ax1.set_yticks(y_pos)
    
    # Replace underscores and format labels
    display_labels = [label.replace('_', ' ').replace(':', ': ') for label in subclasses]
    ax1.set_yticklabels(display_labels, fontsize=8)
    
    # Add value labels
    max_val = max(genomic_vals) if genomic_vals else 1
    for i, (bar, val) in enumerate(zip(bars1, genomic_vals)):
        if val > 0:
            ax1.text(val + max_val * 0.01, i, f'{val:.3f}%', 
                    va='center', fontsize=7, fontweight='bold')
    
    # Positional Density (Frequency)
    bars2 = ax2.barh(y_pos, positional_vals, color=colors1, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('Positional Density (motifs/kbp)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Motif Subclass', fontsize=11, fontweight='bold')
    ax2.set_title('B. Positional Density (λ) by Subclass', fontsize=12, fontweight='bold', pad=10)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(display_labels, fontsize=8)
    
    # Add value labels
    max_val2 = max(positional_vals) if positional_vals else 1
    for i, (bar, val) in enumerate(zip(bars2, positional_vals)):
        if val > 0:
            ax2.text(val + max_val2 * 0.01, i, f'{val:.2f}', 
                    va='center', fontsize=7, fontweight='bold')
    
    plt.suptitle(title, fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    
    return fig


def plot_enrichment_analysis_by_subclass(enrichment_results: Dict[str, Dict[str, Any]],
                                         title: str = "Motif Enrichment Analysis (by Subclass)",
                                         figsize: Tuple[int, int] = (16, 12)) -> plt.Figure:
    """
    Plot enrichment analysis results at subclass level with fold enrichment and p-values.
    
    Args:
        enrichment_results: Dictionary from calculate_enrichment_with_shuffling with by_subclass=True
        title: Plot title
        figsize: Figure size
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    # Extract data (exclude 'Overall' for subclass-specific view)
    subclasses = [k for k in enrichment_results.keys() if k != 'Overall' and ':' in k]
    if not subclasses:
        # Fallback to class level
        subclasses = [k for k in enrichment_results.keys() if k != 'Overall']
    
    # Sort alphabetically
    subclasses.sort()
    
    fold_enrichments = []
    p_values = []
    observed_densities = []
    background_means = []
    
    for subclass in subclasses:
        result = enrichment_results[subclass]
        fe = result.get('fold_enrichment', 0)
        # Handle both string 'Inf' and float infinity values robustly
        if fe == 'Inf' or (isinstance(fe, float) and np.isinf(fe)):
            fe = INFINITE_FOLD_ENRICHMENT_CAP  # Cap infinite values for visualization
        fold_enrichments.append(fe)
        p_values.append(result.get('p_value', 1.0))
        observed_densities.append(result.get('observed_density', 0))
        background_means.append(result.get('background_mean', 0))
    
    # Create figure with subplots
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 1, hspace=0.35, height_ratios=[2, 2, 1.5])
    
    y_pos = np.arange(len(subclasses))
    
    # Get colors based on parent class
    colors = []
    for subclass in subclasses:
        if ':' in subclass:
            parent_class = subclass.split(':')[0]
            colors.append(MOTIF_CLASS_COLORS.get(parent_class, '#808080'))
        else:
            colors.append(MOTIF_CLASS_COLORS.get(subclass, '#808080'))
    
    # 1. Fold Enrichment
    ax1 = fig.add_subplot(gs[0])
    bars1 = ax1.barh(y_pos, fold_enrichments, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax1.axvline(x=1.0, color='red', linestyle='--', linewidth=2, label='No enrichment (FE=1)')
    ax1.set_xlabel('Fold Enrichment', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Motif Subclass', fontsize=11, fontweight='bold')
    ax1.set_title('A. Fold Enrichment by Subclass', fontsize=12, fontweight='bold')
    ax1.set_yticks(y_pos)
    
    # Format labels
    display_labels = [label.replace('_', ' ').replace(':', ': ') for label in subclasses]
    ax1.set_yticklabels(display_labels, fontsize=8)
    ax1.legend(loc='best', fontsize=9)
    
    # Add value labels
    max_fe = max(fold_enrichments) if fold_enrichments else 1
    for i, (bar, val) in enumerate(zip(bars1, fold_enrichments)):
        label_text = f'{val:.2f}' if val < INFINITE_FOLD_ENRICHMENT_CAP else 'Inf'
        ax1.text(val + max_fe * 0.01, i, label_text, 
                va='center', fontsize=7, fontweight='bold')
    
    # 2. P-values
    ax2 = fig.add_subplot(gs[1])
    # Color code by significance
    p_colors = ['green' if p < 0.05 else 'orange' if p < 0.1 else 'red' for p in p_values]
    bars2 = ax2.barh(y_pos, p_values, color=p_colors, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax2.axvline(x=0.05, color='green', linestyle='--', linewidth=2, label='p=0.05')
    ax2.axvline(x=0.1, color='orange', linestyle='--', linewidth=1.5, label='p=0.1')
    ax2.set_xlabel('P-value', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Motif Subclass', fontsize=11, fontweight='bold')
    ax2.set_title('B. Statistical Significance by Subclass', fontsize=12, fontweight='bold')
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(display_labels, fontsize=8)
    ax2.set_xlim(0, max(1.0, max(p_values) * 1.1))
    ax2.legend(loc='best', fontsize=9)
    
    # Add value labels
    for i, (bar, val) in enumerate(zip(bars2, p_values)):
        ax2.text(val + 0.02, i, f'{val:.3f}', va='center', fontsize=7, fontweight='bold')
    
    # 3. Observed vs Background Density (grouped bar chart)
    ax3 = fig.add_subplot(gs[2])
    x = np.arange(len(subclasses))
    width = 0.35
    
    bars3a = ax3.bar(x - width/2, observed_densities, width, label='Observed', 
                    color='steelblue', alpha=0.8, edgecolor='black', linewidth=0.5)
    bars3b = ax3.bar(x + width/2, background_means, width, label='Background (Mean)', 
                    color='coral', alpha=0.8, edgecolor='black', linewidth=0.5)
    
    ax3.set_xlabel('Motif Subclass', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Density (%)', fontsize=11, fontweight='bold')
    ax3.set_title('C. Observed vs. Background Density Comparison', fontsize=12, fontweight='bold')
    ax3.set_xticks(x)
    ax3.set_xticklabels(display_labels, rotation=45, ha='right', fontsize=7)
    ax3.legend(loc='best', fontsize=10)
    ax3.grid(axis='y', alpha=0.3)
    
    plt.suptitle(title, fontsize=14, fontweight='bold', y=0.995)
    
    return fig


def plot_subclass_density_heatmap(motifs: List[Dict[str, Any]], 
                                  sequence_length: int,
                                  window_size: int = 1000,
                                  title: str = "Subclass Density Heatmap",
                                  figsize: Tuple[int, int] = (16, 12)) -> plt.Figure:
    """
    Create a heatmap showing density of each subclass across the sequence.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of the analyzed sequence
        window_size: Window size for density calculation
        title: Plot title
        figsize: Figure size
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style()
    
    if not motifs or sequence_length == 0:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    # Calculate windows
    num_windows = max(1, sequence_length // window_size)
    windows = np.linspace(0, sequence_length, num_windows + 1)
    
    # Get unique subclasses
    subclass_groups = defaultdict(list)
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        subclass_name = motif.get('Subclass', 'Unknown')
        key = f"{class_name}:{subclass_name}"
        subclass_groups[key].append(motif)
    
    subclasses = sorted(subclass_groups.keys())
    
    # Calculate density matrix
    density_matrix = np.zeros((len(subclasses), num_windows))
    
    for i, subclass_key in enumerate(subclasses):
        subclass_motifs = subclass_groups[subclass_key]
        
        for j in range(num_windows):
            window_start = windows[j]
            window_end = windows[j + 1]
            
            # Count motifs in window
            count = 0
            for motif in subclass_motifs:
                motif_start = motif.get('Start', 0) - 1  # 0-based
                motif_end = motif.get('End', 0)
                
                # Check if motif overlaps with window
                if not (motif_end <= window_start or motif_start >= window_end):
                    count += 1
            
            density_matrix[i, j] = count
    
    # Create heatmap with publication quality
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    # Use colorblind-friendly colormap
    im = ax.imshow(density_matrix, cmap='viridis', aspect='auto', interpolation='nearest')
    
    # Customize axes (Nature style)
    ax.set_xlabel(f'Position (kb)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Motif Subclass', fontsize=12, fontweight='bold')
    
    # Replace underscores with spaces in title
    display_title = title.replace('_', ' ')
    ax.set_title(display_title, fontsize=14, fontweight='bold', pad=10)
    
    # Set ticks and labels with underscores replaced by spaces
    ax.set_yticks(range(len(subclasses)))
    display_subclasses = [sc.replace('_', ' ').replace(':', ': ') for sc in subclasses]
    ax.set_yticklabels(display_subclasses, fontsize=8)
    
    # Clean x-axis with kb units
    x_ticks = np.arange(0, num_windows, max(1, num_windows // 10))
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f'{int(windows[i]/1000)}' for i in x_ticks])
    
    # Add colorbar with proper label
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Motif Count', fontsize=10, fontweight='bold')
    cbar.ax.tick_params(labelsize=8)
    
    plt.tight_layout()
    return fig


if __name__ == "__main__":
    test_visualizations()