"""
Enhanced Visualization Functions for Nobel-Level Quality
=========================================================

This module provides enhanced publication-quality visualizations
with improved aesthetics, color schemes, and scientific accuracy.

Author: Dr. Venkata Rajesh Yella
Version: 2025.1
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from typing import List, Dict, Any, Optional, Tuple

# Enhanced color palettes for publication-quality figures
ENHANCED_SCIENTIFIC_PALETTE = {
    'vibrant_blue': '#0077BB',
    'vibrant_cyan': '#33BBEE',
    'vibrant_teal': '#009988',
    'vibrant_orange': '#EE7733',
    'vibrant_red': '#CC3311',
    'vibrant_magenta': '#EE3377',
    'vibrant_purple': '#AA3377',
    'vibrant_grey': '#BBBBBB',
}

# Nobel-level color schemes (colorblind-friendly, high contrast)
NOBEL_QUALITY_COLORS = {
    'Curved_DNA': '#E69F00',       # Orange (high visibility)
    'Z-DNA': '#56B4E9',            # Sky blue (distinct)
    'Slipped_DNA': '#009E73',      # Bluish green (clear)
    'R-Loop': '#F0E442',           # Yellow (striking)
    'Cruciform': '#0072B2',        # Blue (professional)
    'Triplex': '#D55E00',          # Vermillion (bold)
    'G-Quadruplex': '#CC79A7',     # Reddish purple (elegant)
    'G4': '#CC79A7',               # Alias for G-Quadruplex
    'i-Motif': '#999933',          # Olive (balanced)
    'IMotif': '#999933',           # Alias
    'A-philic_DNA': '#882255',    # Wine (sophisticated)
    'A-philic DNA': '#882255',     # Alias
    'Sticky_DNA': '#44AA99',      # Teal (fresh)
    'Hybrid': '#DDAA33',          # Tan (neutral)
    'Non-B_DNA_Clusters': '#666666', # Gray (functional)
}

def enhance_plot_aesthetics(ax, title=None, xlabel=None, ylabel=None,
                           spine_width=1.5, tick_length=5):
    """
    Apply Nobel-level aesthetic enhancements to matplotlib axes.
    
    Features:
    - Thicker spines for better visibility in print
    - Longer ticks for easier reading
    - Optimized font sizes
    - Professional grid styling
    """
    # Spine styling
    for spine in ax.spines.values():
        spine.set_linewidth(spine_width)
        spine.set_color('#2C3E50')  # Dark blue-grey
    
    # Tick styling
    ax.tick_params(axis='both', which='major', 
                   labelsize=11, length=tick_length,
                   width=spine_width, colors='#2C3E50')
    ax.tick_params(axis='both', which='minor', 
                   length=tick_length*0.6, width=spine_width*0.8)
    
    # Grid styling (subtle but present)
    ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.8, color='#7F8C8D')
    ax.set_axisbelow(True)  # Grid behind data
    
    # Labels with enhanced styling
    if title:
        ax.set_title(title, fontsize=14, fontweight='bold', 
                    pad=15, color='#2C3E50')
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=12, fontweight='600', 
                     labelpad=10, color='#2C3E50')
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=12, fontweight='600', 
                     labelpad=10, color='#2C3E50')
    
    return ax


def create_enhanced_colormap(name='viridis', reverse=False):
    """
    Create enhanced colormaps for density and heatmap visualizations.
    
    Returns publication-quality colormaps with better contrast.
    """
    if name == 'density':
        # Enhanced density colormap (white → blue → purple)
        colors = ['#FFFFFF', '#E3F2FD', '#90CAF9', '#42A5F5', 
                  '#1E88E5', '#1565C0', '#0D47A1', '#4A148C']
    elif name == 'heatmap':
        # Enhanced heatmap (yellow → orange → red → purple)
        colors = ['#FFF9E6', '#FFEB99', '#FFD54F', '#FFA726',
                  '#FF7043', '#E64A19', '#C62828', '#6A1B9A']
    elif name == 'diverging':
        # Enhanced diverging (blue → white → red)
        colors = ['#0D47A1', '#42A5F5', '#90CAF9', '#E3F2FD',
                  '#FFFFFF', '#FFEBEE', '#EF9A9A', '#E57373', '#C62828']
    else:
        # Default to matplotlib colormap
        return plt.get_cmap(name)
    
    if reverse:
        colors = colors[::-1]
    
    return LinearSegmentedColormap.from_list(name, colors, N=256)


def add_publication_legend(ax, handles=None, labels=None, 
                          ncol=1, loc='best', frameon=True):
    """
    Add publication-quality legend with enhanced styling.
    """
    if handles is None or labels is None:
        handles, labels = ax.get_legend_handles_labels()
    
    legend = ax.legend(handles, labels, 
                      loc=loc, 
                      ncol=ncol,
                      frameon=frameon,
                      fancybox=False,
                      shadow=False,
                      framealpha=0.95,
                      edgecolor='#2C3E50',
                      facecolor='white',
                      fontsize=10,
                      title_fontsize=11,
                      borderpad=0.8,
                      labelspacing=0.6,
                      handlelength=2.0,
                      handleheight=1.0)
    
    # Enhance legend frame
    if frameon:
        legend.get_frame().set_linewidth(1.2)
    
    return legend


def save_publication_figure(fig, filename, dpi=300, bbox_inches='tight',
                           transparent=False):
    """
    Save figure in publication-ready format.
    
    Features:
    - 300 DPI for Nature/Science journals
    - Tight bounding box
    - Multiple format support
    """
    # Save as PNG (standard)
    fig.savefig(filename, dpi=dpi, bbox_inches=bbox_inches,
                transparent=transparent, facecolor='white')
    
    # Also save as PDF for vector graphics
    pdf_filename = filename.replace('.png', '.pdf')
    fig.savefig(pdf_filename, format='pdf', bbox_inches=bbox_inches,
                transparent=transparent)
    
    print(f"✅ Saved publication-quality figures:")
    print(f"   - {filename} (raster, {dpi} DPI)")
    print(f"   - {pdf_filename} (vector)")


def enhanced_bar_plot(data, labels, colors=None, title=None,
                     xlabel=None, ylabel=None, figsize=(10, 6)):
    """
    Create enhanced bar plot with Nobel-level quality.
    """
    fig, ax = plt.subplots(figsize=figsize, dpi=100)
    
    # Use enhanced colors if not provided
    if colors is None:
        colors = [NOBEL_QUALITY_COLORS.get(label, '#808080') 
                 for label in labels]
    
    # Create bars with enhanced styling
    bars = ax.bar(range(len(data)), data, color=colors, 
                  edgecolor='#2C3E50', linewidth=1.2, alpha=0.85)
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}',
                   ha='center', va='bottom', fontsize=10,
                   fontweight='600', color='#2C3E50')
    
    # Set x-axis labels with rotation for readability
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha='right')
    
    # Apply enhanced aesthetics
    enhance_plot_aesthetics(ax, title=title, xlabel=xlabel, ylabel=ylabel)
    
    plt.tight_layout()
    return fig


# Export key functions
__all__ = [
    'NOBEL_QUALITY_COLORS',
    'ENHANCED_SCIENTIFIC_PALETTE',
    'enhance_plot_aesthetics',
    'create_enhanced_colormap',
    'add_publication_legend',
    'save_publication_figure',
    'enhanced_bar_plot'
]
