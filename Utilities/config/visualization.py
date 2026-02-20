"""
Visualization configuration for NBDScanner.

This module contains settings for plots and charts:
- DPI and figure sizes
- Plot style settings
- Motif class colors
- Grid and legend settings
"""

from Utilities.config.colors import VISUALIZATION_PALETTE

# ==================== VISUALIZATION PARAMETERS ====================
# Control all visualization and plot settings
# Modify these to customize the appearance of charts and graphs
VISUALIZATION_CONFIG = {
    # Plot resolution (DPI - dots per inch)
    # 72: Standard screen resolution (fastest rendering, lowest memory)
    # 150: Good for screen viewing
    # 300: Publication quality (Nature, Science journals)
    # 600: High-quality print (posters, large format)
    'dpi': 72,
    
    # Default figure sizes (width, height) in inches
    # 1 inch ≈ 2.54 cm; figures scale proportionally
    'figure_sizes': {
        'small': (8, 5),       # Compact plots for quick visualization
        'medium': (10, 6),     # Standard plots for reports
        'large': (12, 8),      # Large plots for presentations
        'wide': (14, 6),       # Wide plots for timeline/distribution
        'tall': (8, 10),       # Tall plots for vertical data
    },
    
    # Plot style (affects overall appearance)
    # Options: 'nature', 'default', 'presentation', 'seaborn'
    'default_style': 'seaborn',
    
    # Colors for motif classes (Wong 2011 colorblind-friendly palette)
    # These colors ensure accessibility for colorblind users
    # All colors now reference the centralized VISUALIZATION_PALETTE
    'motif_class_colors': {
        'Curved DNA': VISUALIZATION_PALETTE['chart_1'],           # Orange
        'Z-DNA': VISUALIZATION_PALETTE['chart_2'],                # Sky blue
        'Slipped DNA': VISUALIZATION_PALETTE['chart_3'],          # Bluish green
        'R-Loop': VISUALIZATION_PALETTE['chart_4'],               # Yellow
        'Cruciform': VISUALIZATION_PALETTE['chart_5'],            # Blue
        'Triplex DNA': VISUALIZATION_PALETTE['chart_6'],          # Vermillion
        'G-Quadruplex': VISUALIZATION_PALETTE['chart_7'],         # Reddish purple
        'G4': VISUALIZATION_PALETTE['chart_7'],                   # Reddish purple (alias)
        'i-Motif': VISUALIZATION_PALETTE['chart_8'],              # Olive
        'A-philic DNA': VISUALIZATION_PALETTE['chart_9'],         # Wine
        'Sticky DNA': VISUALIZATION_PALETTE['chart_10'],          # Teal
        'Hybrid': VISUALIZATION_PALETTE['chart_11'],              # Tan (for overlapping motifs)
        'Non-B DNA Clusters': VISUALIZATION_PALETTE['chart_12'],   # Gray
        'Non-B_DNA_Clusters': VISUALIZATION_PALETTE['chart_12'],   # Gray (alternate name)
    },
    
    # Grid and layout settings for plots
    'show_grid': True,                # Show background grid on plots
    'grid_alpha': 0.5,                # Grid transparency (0-1)
    'grid_style': '--',               # Grid line style ('--', '-', ':', '-.')
    
    # Legend settings for plot labels
    'legend_fontsize': 9,             # Font size for legend text
    'legend_location': 'upper right',        # Legend position ('best', 'upper right', etc.)
    'legend_frameon': False,           # Show frame around legend
    'legend_framealpha': 0.9,         # Legend background transparency
}
