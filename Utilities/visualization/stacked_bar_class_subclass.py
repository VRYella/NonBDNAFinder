# ═══════════════════════════════════════════════════════════════════════════
# FUNCTION          | plot_stacked_bar_class_subclass, plot_nested_pie_chart
# PURPOSE           | Stacked bar: Class(X) → Subclass(stacked Y) visualization
# INPUTS            | motifs(List[Dict]), title(str), figsize(tuple), exclude_classes(List[str])
# OUTPUTS           | plt.Figure object
# COLOR SCHEME      | UNIFIED_MOTIF_COLORS base + alpha gradient (0.4-1.0) for subclasses
# ORDERING          | Biological taxonomy via get_all_classes_taxonomy_order()
# ANNOTATIONS       | Count on segment (if ≥3), total above bar (n=X)
# COMPATIBILITY     | plot_nested_pie_chart() → wrapper redirect to stacked bar
# ═══════════════════════════════════════════════════════════════════════════

import matplotlib.pyplot as plt; import numpy as np; from collections import Counter; from matplotlib.colors import to_rgba
from Utilities.config.colors import UNIFIED_MOTIF_COLORS; from Utilities.config.motif_taxonomy import get_all_classes_taxonomy_order, get_subclasses_for_class

def plot_stacked_bar_class_subclass(motifs, title="Class → Subclass Composition", figsize=(14, 8), exclude_classes=None):
    if not motifs: fig, ax = plt.subplots(figsize=(8, 6)); ax.text(0.5, 0.5, 'No motifs', ha='center', va='center', fontsize=16); ax.axis('off'); return fig
    motifs = [m for m in motifs if m.get('Class') not in (exclude_classes or [])]
    if not motifs: fig, ax = plt.subplots(figsize=(8, 6)); ax.text(0.5, 0.5, 'All excluded', ha='center', va='center', fontsize=16); ax.axis('off'); return fig
    
    class_totals = Counter(); class_subclass = {}
    for m in motifs: cls, sub = m.get('Class', 'Unknown'), m.get('Subclass', 'Other'); class_totals[cls] += 1; class_subclass.setdefault(cls, Counter())[sub] += 1
    
    taxonomy = get_all_classes_taxonomy_order(); ordered = [c for c in taxonomy if c in class_totals] + [c for c in class_totals if c not in taxonomy]
    fig, ax = plt.subplots(figsize=figsize); x_pos = np.arange(len(ordered)); bottom = np.zeros(len(ordered))
    
    for cls_idx, cls in enumerate(ordered):
        base_rgba = to_rgba(UNIFIED_MOTIF_COLORS.get(cls, '#808080')); subs = class_subclass.get(cls, {})
        try: ordered_subs = [s for s in get_subclasses_for_class(cls) if s in subs]
        except: ordered_subs = sorted(subs.keys())
        
        for idx, sub in enumerate(ordered_subs):
            cnt = subs[sub]; alpha = 0.4 + 0.6 * (idx / max(1, len(ordered_subs) - 1)) if len(ordered_subs) > 1 else 1.0
            ax.bar(x_pos[cls_idx], cnt, 0.7, bottom=bottom[cls_idx], color=(*base_rgba[:3], alpha), edgecolor='white', linewidth=1.5)
            if cnt >= 3: ax.text(x_pos[cls_idx], bottom[cls_idx] + cnt/2, f"{cnt}", ha='center', va='center', fontsize=9, fontweight='bold', color='white')
            bottom[cls_idx] += cnt
        ax.text(x_pos[cls_idx], bottom[cls_idx] + 1, f"n={class_totals[cls]}", ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    ax.set_xticks(x_pos); ax.set_xticklabels(ordered, rotation=45, ha='right', fontsize=12, fontweight='bold')
    for lbl, cls in zip(ax.get_xticklabels(), ordered): lbl.set_color(UNIFIED_MOTIF_COLORS.get(cls, '#808080'))
    ax.set_ylabel('Motif Count', fontsize=14, fontweight='bold'); ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.grid(axis='y', alpha=0.3, linestyle='--')
    plt.tight_layout(); return fig

def plot_nested_pie_chart(motifs, title="Class → Subclass", **kwargs): return plot_stacked_bar_class_subclass(motifs, title=title, **kwargs)
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║         STACKED BAR CHART: CLASS → SUBCLASS VISUALIZATION                   ║
║     Publication-Quality Stacked Bar Plot for Motif Hierarchy                 ║
╚══════════════════════════════════════════════════════════════════════════════╝

Replaces nested pie chart with clearer stacked bar visualization.
Follows Nature/Science figure standards for publication-ready output.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import Counter
from typing import List, Dict, Any, Optional, Tuple
from matplotlib.colors import to_rgba

from Utilities.config.colors import UNIFIED_MOTIF_COLORS
from Utilities.config.motif_taxonomy import (
    get_all_classes_taxonomy_order,
    get_subclasses_for_class,
    SUBCLASS_TO_CLASS
)


def plot_stacked_bar_class_subclass(
    motifs: List[Dict[str, Any]],
    title: str = "Motif Distribution: Class → Subclass",
    figsize: Optional[Tuple[float, float]] = None,
    exclude_classes: Optional[List[str]] = None
) -> plt.Figure:
    """
    Create publication-quality stacked bar plot showing Class → Subclass hierarchy.
    
    Replaces nested pie chart with clearer bar visualization suitable for scientific
    publications. Shows motif counts stacked by subclasses within each class.
    
    Args:
        motifs: List of motif dictionaries with 'Class' and 'Subclass' keys
        title: Plot title
        figsize: Figure size (width, height) in inches. Default: (12, 6)
        exclude_classes: Optional list of class names to exclude from plot
        
    Returns:
        Matplotlib figure object (publication-ready)
        
    Design Features:
        - X-axis: Classes ordered by biological taxonomy
        - Y-axis: Motif counts stacked by subclasses
        - Colors: UNIFIED_MOTIF_COLORS with alpha gradients (0.4-1.0) for subclasses
        - Annotations: Count labels on segments (if count ≥ 3), total above each bar
        - Legend: Class → Subclass mapping (limited to 20 entries)
        - Styling: Color-coded X-axis labels, minimal spines, light Y-grid
    """
    # Default figure size for horizontal layout
    if figsize is None: figsize = (12, 6)
    
    # Handle empty motif list
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize); ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center', fontsize=14)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off'); return fig
    
    # Extract exclude list
    exclude_classes = exclude_classes or []
    
    # ═══════════════════════════════════════════════════════════════════════════
    # DATA AGGREGATION: Count by Class and Subclass
    # ═══════════════════════════════════════════════════════════════════════════
    class_subclass_counts = {}  # {class: {subclass: count}}
    for motif in motifs:
        cls = motif.get('Class', 'Unknown'); subcls = motif.get('Subclass', 'Unknown')
        if cls in exclude_classes: continue
        if cls not in class_subclass_counts: class_subclass_counts[cls] = {}
        class_subclass_counts[cls][subcls] = class_subclass_counts[cls].get(subcls, 0) + 1
    
    # Handle case where all motifs are excluded
    if not class_subclass_counts:
        fig, ax = plt.subplots(figsize=figsize); ax.text(0.5, 0.5, 'All motifs excluded', ha='center', va='center', fontsize=14)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off'); return fig
    
    # Calculate class totals for ordering
    class_totals = {cls: sum(subclass_counts.values()) for cls, subclass_counts in class_subclass_counts.items()}
    
    # ═══════════════════════════════════════════════════════════════════════════
    # TAXONOMY ORDERING: Order classes by biological taxonomy
    # ═══════════════════════════════════════════════════════════════════════════
    taxonomy_classes = get_all_classes_taxonomy_order(); ordered_classes = [cls for cls in taxonomy_classes if cls in class_totals]
    ordered_classes += [cls for cls in class_totals if cls not in ordered_classes]  # Fallback for any missing
    
    # ═══════════════════════════════════════════════════════════════════════════
    # COLOR SCHEME: Generate colors for subclasses with alpha gradients
    # ═══════════════════════════════════════════════════════════════════════════
    subclass_color_map = {}  # {class::subclass: rgba}
    all_subclasses = []  # For legend ordering
    
    for cls in ordered_classes:
        subclass_counts = class_subclass_counts[cls]; base_color = UNIFIED_MOTIF_COLORS.get(cls, '#808080')
        base_rgba = to_rgba(base_color); n_subs = len(subclass_counts)
        
        # Order subclasses by taxonomy
        try:
            taxonomy_subclasses = get_subclasses_for_class(cls); ordered_subclasses = [sub for sub in taxonomy_subclasses if sub in subclass_counts]
            ordered_subclasses += [sub for sub in subclass_counts if sub not in ordered_subclasses]  # Fallback
        except (KeyError, ValueError):
            ordered_subclasses = sorted(subclass_counts.keys())  # Alphabetical fallback
        
        # Assign alpha gradients: lighter (0.4) to darker (1.0)
        for idx, subcls in enumerate(ordered_subclasses):
            alpha = 0.4 + 0.6 * (idx / max(1, n_subs - 1)) if n_subs > 1 else 1.0; subclass_color_map[f"{cls}::{subcls}"] = (*base_rgba[:3], alpha)
            all_subclasses.append((cls, subcls))
    
    # ═══════════════════════════════════════════════════════════════════════════
    # PLOT CONSTRUCTION: Build stacked bar chart
    # ═══════════════════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=figsize, dpi=72)
    
    x_positions = np.arange(len(ordered_classes)); bar_width = 0.6; bottoms = np.zeros(len(ordered_classes))
    
    # Track which subclasses to include in legend (limit to 20)
    legend_entries = []; legend_added = set()
    
    # Build stacked bars
    for cls_idx, cls in enumerate(ordered_classes):
        subclass_counts = class_subclass_counts[cls]
        
        # Order subclasses same as color generation
        try:
            taxonomy_subclasses = get_subclasses_for_class(cls); ordered_subclasses = [sub for sub in taxonomy_subclasses if sub in subclass_counts]
            ordered_subclasses += [sub for sub in subclass_counts if sub not in ordered_subclasses]
        except (KeyError, ValueError):
            ordered_subclasses = sorted(subclass_counts.keys())
        
        for subcls in ordered_subclasses:
            count = subclass_counts[subcls]; color = subclass_color_map[f"{cls}::{subcls}"]
            
            # Draw bar segment
            bar = ax.bar(x_positions[cls_idx], count, bar_width, bottom=bottoms[cls_idx], color=color, edgecolor='white', linewidth=1.5)
            
            # Add count annotation if significant (count ≥ 3)
            if count >= 3:
                y_pos = bottoms[cls_idx] + count / 2; ax.text(x_positions[cls_idx], y_pos, str(count), ha='center', va='center', fontsize=10, fontweight='bold', color='white')
            
            # Add to legend (limit to 20 entries)
            if len(legend_entries) < 20:
                legend_key = f"{cls}::{subcls}"
                if legend_key not in legend_added:
                    legend_entries.append((f"{cls} → {subcls}", color)); legend_added.add(legend_key)
            
            bottoms[cls_idx] += count
    
    # ═══════════════════════════════════════════════════════════════════════════
    # ANNOTATIONS: Add total counts above each bar
    # ═══════════════════════════════════════════════════════════════════════════
    for idx, cls in enumerate(ordered_classes):
        total = class_totals[cls]; ax.text(x_positions[idx], total + 0.02 * ax.get_ylim()[1], f'n={total}', ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    # ═══════════════════════════════════════════════════════════════════════════
    # STYLING: Apply publication-quality formatting
    # ═══════════════════════════════════════════════════════════════════════════
    # Set axis labels and title
    ax.set_xlabel('Motif Class', fontsize=14, fontweight='bold'); ax.set_ylabel('Motif Count', fontsize=14, fontweight='bold')
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    
    # X-axis: Color-coded labels by class color
    ax.set_xticks(x_positions); ax.set_xticklabels([cls.replace('_', ' ') for cls in ordered_classes], fontsize=12, rotation=45, ha='right')
    for idx, (tick_label, cls) in enumerate(zip(ax.get_xticklabels(), ordered_classes)):
        color = UNIFIED_MOTIF_COLORS.get(cls, '#000000'); tick_label.set_color(color); tick_label.set_fontweight('bold')
    
    # Y-axis: Clean styling
    ax.tick_params(axis='y', labelsize=12); ax.yaxis.grid(True, alpha=0.3, linestyle='--', linewidth=0.8)
    
    # Remove top and right spines for clean appearance
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.spines['left'].set_linewidth(1.5); ax.spines['bottom'].set_linewidth(1.5)
    
    # Add legend
    if legend_entries:
        handles = [plt.Rectangle((0, 0), 1, 1, fc=color, ec='white', lw=1.5) for _, color in legend_entries]; labels = [label for label, _ in legend_entries]
        ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=10, frameon=True, title='Class → Subclass', title_fontsize=11, framealpha=0.95, edgecolor='lightgray')
    
    # Tight layout to prevent label cutoff
    plt.tight_layout()
    
    return fig


def plot_nested_pie_chart(
    motifs: List[Dict[str, Any]],
    title: str = "Motif Distribution",
    figsize: Optional[Tuple[float, float]] = None,
    exclude_classes: Optional[List[str]] = None
) -> plt.Figure:
    """
    Backward-compatible wrapper for plot_stacked_bar_class_subclass().
    
    This function maintains compatibility with existing code that calls
    plot_nested_pie_chart() but now generates a stacked bar chart instead.
    
    Args:
        motifs: List of motif dictionaries with 'Class' and 'Subclass' keys
        title: Plot title
        figsize: Figure size (width, height) in inches
        exclude_classes: Optional list of class names to exclude from plot
        
    Returns:
        Matplotlib figure object (publication-ready stacked bar chart)
    """
    # Redirect to new stacked bar implementation
    if title == "Motif Distribution": title = "Motif Distribution: Class → Subclass"
    return plot_stacked_bar_class_subclass(motifs, title=title, figsize=figsize, exclude_classes=exclude_classes)
