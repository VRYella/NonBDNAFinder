"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    UNIFIED VISUALIZATION MODULE                               ║
║          Consolidated Non-B DNA Motif Visualization Functions                 ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: visualization/__init__.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.1
LICENSE: MIT

DESCRIPTION:
    Unified entry point for all NonBDNAFinder visualization functionality.
    This module consolidates visualization functions and standards into a single
    coherent package for easier imports and maintenance.

USAGE:
    # Import all visualization functions
    from visualization import plot_motif_distribution, plot_nested_pie_chart
    
    # Import visualization standards
    from visualization import NATURE_MOTIF_COLORS, PlotDominance, FigurePanel
    
    # Import density calculation functions
    from visualization import calculate_genomic_density, calculate_positional_density
"""

# =============================================================================
# VISUALIZATION STANDARDS - Nature-ready visualization configuration
# =============================================================================
from Utilities.visualization.standards import (
    # Color schemes
    NATURE_MOTIF_COLORS,
    SUBCLASS_TONE_ADJUSTMENTS,
    
    # Plot management classes
    PlotConcept,
    PlotDominance,
    FigurePanel,
    MetricFilter,
    LabelPolicy,
    UILayout,
    ValidationThresholds,
    
    # Transparency notes
    TRANSPARENCY_NOTE,
    SUPPLEMENTARY_NOTE,
    
    # Helper functions
    get_plot_mapping,
    should_show_plot,
    get_nature_style_params,
)

# =============================================================================
# STACKED BAR VISUALIZATION - New Class → Subclass plot
# =============================================================================
from Utilities.visualization.stacked_bar_class_subclass import (
    plot_stacked_bar_class_subclass,
    plot_nested_pie_chart,  # Backward-compatible wrapper
)

# =============================================================================
# DENSITY CALCULATION FUNCTIONS
# =============================================================================
from Utilities.utilities import (
    calculate_genomic_density,
    calculate_positional_density,
)

# =============================================================================
# PLOTTING FUNCTIONS - Comprehensive visualization suite
# =============================================================================
from Utilities.visualization.stacked_bar_class_subclass import plot_stacked_bar_class_subclass, plot_nested_pie_chart

from Utilities.utilities import (
    # Core distribution plots
    plot_motif_distribution,
    plot_class_subclass_sunburst,
    
    # Coverage and density visualization
    plot_coverage_map,
    plot_density_heatmap,
    plot_density_comparison,
    plot_density_comparison_by_subclass,
    
    # Score and length distributions
    plot_score_distribution,
    plot_length_distribution,
    plot_motif_length_kde,
    plot_motif_length_histogram,
    plot_score_violin,
    
    # Comprehensive analysis plots
    plot_class_comparison,
    plot_class_analysis_comprehensive,
    plot_subclass_analysis_comprehensive,
    plot_score_statistics_by_class,
    plot_length_statistics_by_class,
    
    # Genome-scale visualizations
    plot_manhattan_motif_density,
    plot_linear_motif_track,
    plot_genome_landscape_track,
    plot_cumulative_motif_distribution,
    
    # Circular and radial plots
    plot_circos_motif_density,
    plot_radial_class_density,
    
    # Track-based visualizations
    plot_stacked_density_track,
    plot_sliding_window_heat_ribbon,
    
    # Relationship and co-occurrence
    plot_motif_cooccurrence_matrix,
    plot_gc_content_correlation,
    
    # Cluster analysis
    plot_cluster_size_distribution,
    
    # Enrichment analysis
    plot_enrichment_analysis,
    plot_enrichment_summary_table,
    plot_enrichment_analysis_by_subclass,
    
    # Subclass-specific plots
    plot_subclass_density_heatmap,
    
    # Advanced visualization functions (Nature publication standard)
    plot_structural_heatmap,
    plot_gc_motif_correlation,
    plot_motif_network,
    plot_chromosome_density,
    plot_spacer_loop_variation,
    plot_motif_clustering_distance,
    plot_structural_competition_upset,
)

# =============================================================================
# MODULE METADATA
# =============================================================================
__version__ = "2024.1"
__author__ = "Dr. Venkata Rajesh Yella"

__all__ = [
    # Visualization Standards
    'NATURE_MOTIF_COLORS',
    'SUBCLASS_TONE_ADJUSTMENTS',
    'PlotConcept',
    'PlotDominance',
    'FigurePanel',
    'MetricFilter',
    'LabelPolicy',
    'UILayout',
    'ValidationThresholds',
    'TRANSPARENCY_NOTE',
    'SUPPLEMENTARY_NOTE',
    'get_plot_mapping',
    'should_show_plot',
    'get_nature_style_params',
    
    # Density calculations
    'calculate_genomic_density',
    'calculate_positional_density',
    
    # Core distribution plots
    'plot_motif_distribution',
    'plot_stacked_bar_class_subclass',
    'plot_nested_pie_chart',  # Now redirects to stacked bar
    'plot_class_subclass_sunburst',
    
    # Coverage and density visualization
    'plot_coverage_map',
    'plot_density_heatmap',
    'plot_density_comparison',
    'plot_density_comparison_by_subclass',
    
    # Score and length distributions
    'plot_score_distribution',
    'plot_length_distribution',
    'plot_motif_length_kde',
    'plot_motif_length_histogram',
    'plot_score_violin',
    
    # Comprehensive analysis plots
    'plot_class_comparison',
    'plot_class_analysis_comprehensive',
    'plot_subclass_analysis_comprehensive',
    'plot_score_statistics_by_class',
    'plot_length_statistics_by_class',
    
    # Genome-scale visualizations
    'plot_manhattan_motif_density',
    'plot_linear_motif_track',
    'plot_genome_landscape_track',
    'plot_cumulative_motif_distribution',
    
    # Circular and radial plots
    'plot_circos_motif_density',
    'plot_radial_class_density',
    
    # Track-based visualizations
    'plot_stacked_density_track',
    'plot_sliding_window_heat_ribbon',
    
    # Relationship and co-occurrence
    'plot_motif_cooccurrence_matrix',
    'plot_gc_content_correlation',
    
    # Cluster analysis
    'plot_cluster_size_distribution',
    
    # Enrichment analysis
    'plot_enrichment_analysis',
    'plot_enrichment_summary_table',
    'plot_enrichment_analysis_by_subclass',
    
    # Subclass-specific plots
    'plot_subclass_density_heatmap',
    
    # Advanced visualization functions (Nature publication standard)
    'plot_structural_heatmap',
    'plot_gc_motif_correlation',
    'plot_motif_network',
    'plot_chromosome_density',
    'plot_spacer_loop_variation',
    'plot_motif_clustering_distance',
    'plot_structural_competition_upset',
]
