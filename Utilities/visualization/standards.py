"""
╔══════════════════════════════════════════════════════════════════════════════╗
║              NATURE-READY VISUALIZATION STANDARDS MODULE                      ║
║          Publication-First Visualization & Redundancy Elimination             ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: visualization_standards.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2025.1
LICENSE: MIT

DESCRIPTION:
    Centralized configuration for Nature-ready visualization standards.
    Implements automatic redundancy elimination, plot dominance rules,
    and journal-ready figure panel layouts.

CORE PRINCIPLES:
    1. Eliminate redundant plots automatically
    2. Expose only biologically interpretable metrics
    3. Guarantee non-overlapping, uncluttered figures
    4. Produce journal-ready panels with minimal user tuning
    5. Preserve all advanced outputs for export/supplementary use

REFERENCE:
    Design follows standards from:
    - Nature journal figure guidelines
    - Broad Institute visualization best practices
    - EMBL-EBI data visualization standards
"""

from typing import Dict, List, Set, Tuple, Optional
from enum import Enum

# =============================================================================
# MOTIF CLASS COLOR SCHEME - UNIFIED ACROSS ALL VISUALIZATIONS
# =============================================================================
# Single source of truth: config/colors.py → UNIFIED_MOTIF_COLORS
# All 11 motif classes have unique, colorblind-friendly colors.
# This module re-exports as NATURE_MOTIF_COLORS for backward compatibility.

from Utilities.config.colors import VISUALIZATION_MOTIF_COLORS as NATURE_MOTIF_COLORS

# Subclass color variations (tonal only - lighter/darker versions)
SUBCLASS_TONE_ADJUSTMENTS = {
    'lighter': 0.3,   # Add 30% lightness
    'darker': -0.2,   # Subtract 20% lightness
    'saturate': 0.15  # Add 15% saturation
}

# =============================================================================
# PLOT REDUNDANCY ELIMINATION RULES
# =============================================================================

class PlotConcept(Enum):
    """Biological concepts that should have only ONE visualization."""
    COMPOSITION = "composition"          # What structures exist
    POSITION = "position"                # Where structures are located
    DENSITY = "density"                  # How abundant are structures
    LENGTH = "length"                    # Physical size distribution
    SCORE = "score"                      # Quality/confidence
    CLUSTERING = "clustering"            # Co-localization patterns
    COOCCURRENCE = "cooccurrence"       # Which structures appear together
    CORRELATION = "correlation"          # Relationships (e.g., GC content)

class PlotDominance:
    """
    Plot dominance rules: For each biological concept, 
    only ONE visualization is shown in main UI.
    """
    
    # Concept → Allowed plot(s)
    DOMINANT_PLOTS = {
        PlotConcept.COMPOSITION: [
            "nested_pie_chart",          # Class → Subclass hierarchy
            "motif_distribution"         # Bar plots for class and subclass counts
        ],
        PlotConcept.POSITION: [
            "manhattan_plot",            # Large genomes (>50kb)
            "linear_track"               # Small sequences (<50kb)
        ],
        PlotConcept.DENSITY: [
            "density_comparison"         # Combined genomic & positional
        ],
        PlotConcept.LENGTH: [
            "length_kde"                 # Smooth distributions by class
        ],
        PlotConcept.SCORE: [
            "score_distribution"         # Histogram by class
        ],
        PlotConcept.CLUSTERING: [
            "cluster_size_distribution"  # Cluster statistics
        ],
        PlotConcept.COOCCURRENCE: [
            "cooccurrence_matrix"        # Heatmap of motif pairs
        ],
        PlotConcept.CORRELATION: [
            "gc_correlation"             # GC content scatter
        ]
    }
    
    # Plots automatically hidden (redundant)
    HIDDEN_PLOTS = {
        "plot_class_subclass_sunburst",  # Redundant with nested pie
        "plot_coverage_map",             # Redundant with Manhattan/linear
        "plot_density_heatmap",          # Redundant with density_comparison
        "plot_circos_motif_density",     # Too complex, low information
        "plot_radial_class_density",     # Redundant with density_comparison
        "plot_stacked_density_track",    # Redundant with Manhattan
        "plot_cumulative_motif_distribution",  # Low biological value
        "plot_genome_landscape_track",   # Redundant with linear track
        "plot_sliding_window_heat_ribbon"  # Too complex
    }
    
    @classmethod
    def get_allowed_plots(cls, concept: PlotConcept) -> List[str]:
        """Get list of allowed plots for a biological concept."""
        return cls.DOMINANT_PLOTS.get(concept, [])
    
    @classmethod
    def is_plot_allowed(cls, plot_name: str) -> bool:
        """Check if a plot is allowed in main UI."""
        return plot_name not in cls.HIDDEN_PLOTS

# =============================================================================
# NATURE-READY FIGURE PANEL LAYOUT
# =============================================================================

class FigurePanel:
    """
    Canonical figure panel layout for publication.
    Enforced regardless of genome size or analysis scope.
    """
    
    # Figure 1: Global Non-B DNA Landscape
    FIGURE_1 = {
        'title': 'Global Non-B DNA Landscape',
        'panels': {
            'A': {
                'plot': 'nested_pie_chart',
                'description': 'Motif composition (Class → Subclass)',
                'required': True
            },
            'B': {
                'plot': 'manhattan_plot',  # OR linear_track for small sequences
                'description': 'Genome-scale localization',
                'required': True,
                'size_dependent': True
            },
            'C': {
                'plot': 'density_comparison',
                'description': 'Genome coverage (% per class)',
                'required': True
            }
        },
        'purpose': 'What structures exist, and where?'
    }
    
    # Figure 2: Structural Clustering & Co-Occurrence
    FIGURE_2 = {
        'title': 'Structural Clustering & Co-Occurrence',
        'panels': {
            'D': {
                'plot': 'cluster_size_distribution',
                'description': 'Cluster size distribution',
                'required': False,  # Only if clusters exist
                'conditional': 'has_clusters'
            },
            'E': {
                'plot': 'cooccurrence_matrix',
                'description': 'Motif co-occurrence (upper triangle)',
                'required': True
            }
        },
        'purpose': 'Which structures co-localize and form regulatory hotspots?'
    }
    
    # Figure 3: Structural Constraints (Optional/Toggle)
    FIGURE_3 = {
        'title': 'Structural Constraints',
        'panels': {
            'F': {
                'plot': 'length_kde',
                'description': 'Length distributions by class',
                'required': False,
                'toggle': True
            }
        },
        'purpose': 'What are the physical constraints shaping each motif class?',
        'optional': True
    }
    
    @classmethod
    def get_required_plots(cls) -> List[str]:
        """Get list of always-required plots for main UI."""
        required = []
        for fig in [cls.FIGURE_1, cls.FIGURE_2]:
            for panel_data in fig['panels'].values():
                if panel_data.get('required', True):
                    required.append(panel_data['plot'])
        return required
    
    @classmethod
    def get_optional_plots(cls) -> List[str]:
        """Get list of optional/toggle plots."""
        optional = []
        for fig in [cls.FIGURE_1, cls.FIGURE_2, cls.FIGURE_3]:
            for panel_data in fig['panels'].values():
                if not panel_data.get('required', True) or panel_data.get('toggle', False):
                    optional.append(panel_data['plot'])
        return optional

# =============================================================================
# METRIC DE-DUPLICATION RULES
# =============================================================================

class MetricFilter:
    """
    Rules for which metrics to display vs hide.
    Only normalized, cross-class comparable metrics are shown.
    """
    
    # Core metrics (always shown)
    CORE_METRICS = {
        'Score',              # Normalized 1-3 scale
        'Start',              # Position
        'End',                # Position
        'Length',             # Size
        'Strand',             # Orientation
        'Class',              # Category
        'Subclass'            # Sub-category
    }
    
    # Hidden raw/intermediate values (available in exports only)
    HIDDEN_METRICS = {
        'Raw_DeltaG',         # Raw thermodynamic
        'Unnormalized_Score', # Pre-normalization
        'Intermediate_Energy', # Calculation intermediate
        'Raw_Count',          # Pre-filtered count
        'Percentile'          # Redundant with Score
    }
    
    # Motif-specific columns (conditionally shown)
    CONDITIONAL_METRICS = {
        'G-Quadruplex': ['Num_Tracts', 'Loop_Length', 'Priority'],
        'Slipped_DNA': ['Repeat_Unit', 'Unit_Length', 'Repeat_Count'],
        'Cruciform': ['Arm_Length', 'Loop_Length', 'Num_Stems'],
        'R-Loop': ['GC_Skew', 'RIZ_Length', 'REZ_Length'],
        'Non-B_DNA_Clusters': ['Motif_Count', 'Class_Diversity']
    }
    
    @classmethod
    def should_show_metric(cls, metric_name: str, motif_class: Optional[str] = None) -> bool:
        """Determine if a metric should be displayed."""
        # Always show core metrics
        if metric_name in cls.CORE_METRICS:
            return True
        
        # Never show hidden metrics in UI
        if metric_name in cls.HIDDEN_METRICS:
            return False
        
        # Show motif-specific metrics only for relevant class
        if motif_class:
            allowed = cls.CONDITIONAL_METRICS.get(motif_class, [])
            return metric_name in allowed
        
        return False

# =============================================================================
# OVERLAP SUPPRESSION CONFIGURATION
# =============================================================================

class LabelPolicy:
    """Configuration for label rendering and overlap suppression."""
    
    # Label rendering rules
    RENDER_INDIVIDUAL_LABELS = False      # Never label individual motifs
    RENDER_CLUSTER_CENTROIDS = True       # Only cluster centers
    RENDER_TOP_QUANTILE = True            # Only top 5% density hotspots
    
    # Overlap suppression
    MIN_LABEL_DISTANCE_BP = 1000         # Minimum distance between labels
    MAX_LABELS_PER_PLOT = 20             # Maximum labels on any plot
    LABEL_PRIORITY_THRESHOLD = 0.95      # Top 5% by density/score
    
    # Annotation settings
    ANNOTATION_FONTSIZE = 7              # Small, professional
    ANNOTATION_ALPHA = 0.7               # Semi-transparent
    USE_COLLISION_DETECTION = True       # Programmatic suppression

# =============================================================================
# UI COMPACTNESS SETTINGS
# =============================================================================

class UILayout:
    """Settings for compact, information-dense layout."""
    
    # Vertical spacing
    SECTION_PADDING_PX = 8              # Reduced from default
    SUBSECTION_PADDING_PX = 4           # Minimal spacing
    
    # Layout patterns
    USE_TWO_COLUMN_LAYOUT = True        # Text + plot side-by-side
    REMOVE_EMPTY_SPACERS = True         # No unnecessary whitespace
    AUTO_WRAP_LABELS = True             # Prevent overflow
    
    # Card/metric heights
    METRIC_CARD_HEIGHT_PX = 80          # Standardized
    PLOT_CARD_HEIGHT_PX = 400           # Consistent
    TABLE_MAX_ROWS = 10                 # Compact tables

# =============================================================================
# SCIENTIFIC TRANSPARENCY BADGE
# =============================================================================

TRANSPARENCY_NOTE = """
📊 **Scientific Transparency**: Only biologically interpretable, non-redundant 
metrics are displayed. Full results including raw thermodynamic values and all 
intermediate calculations are available in CSV/Excel exports.
"""

SUPPLEMENTARY_NOTE = """
💡 **Supplementary Figures**: Additional visualizations (circos plots, radial 
tracks, cumulative distributions) are available on request for supplementary 
materials.
"""

# =============================================================================
# VALIDATION THRESHOLDS
# =============================================================================

class ValidationThresholds:
    """Quality checks for Nature-ready figures."""
    
    MAX_COLORS_PER_FIGURE = 6           # Nature guideline
    MAX_TEXT_OVERLAP_PCT = 0.0          # Zero tolerance
    MIN_DPI = 72                        # Screen display quality
    REQUIRED_FORMATS = ['png', 'pdf']   # Raster + vector
    
    @classmethod
    def validate_color_count(cls, colors: List[str]) -> bool:
        """Ensure figure doesn't exceed color limit."""
        unique_colors = set(colors)
        return len(unique_colors) <= cls.MAX_COLORS_PER_FIGURE
    
    @classmethod
    def validate_figure_quality(cls, dpi: int, formats: List[str]) -> bool:
        """Ensure figure meets quality standards."""
        dpi_ok = dpi >= cls.MIN_DPI
        formats_ok = all(fmt in formats for fmt in cls.REQUIRED_FORMATS)
        return dpi_ok and formats_ok

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_plot_mapping() -> Dict[str, str]:
    """
    Map old plot function names to new standardized names.
    Used for backward compatibility during transition.
    """
    return {
        # Keep these
        'plot_nested_pie_chart': 'nested_pie_chart',
        'plot_manhattan_motif_density': 'manhattan_plot',
        'plot_linear_motif_track': 'linear_track',
        'plot_density_comparison': 'density_comparison',
        'plot_cluster_size_distribution': 'cluster_size_distribution',
        'plot_motif_cooccurrence_matrix': 'cooccurrence_matrix',
        'plot_gc_content_correlation': 'gc_correlation',
        'plot_motif_length_kde': 'length_kde',
        'plot_score_distribution': 'score_distribution',
        'plot_motif_distribution': 'motif_distribution',
        
        # Hide these (redundant)
        'plot_coverage_map': 'HIDDEN',
        'plot_density_heatmap': 'HIDDEN',
        'plot_circos_motif_density': 'HIDDEN',
        'plot_cumulative_motif_distribution': 'HIDDEN',
        'plot_length_distribution': 'HIDDEN'  # Use KDE instead
    }

def should_show_plot(plot_name: str, sequence_length: int, has_clusters: bool = False) -> bool:
    """
    Determine if a plot should be shown based on context.
    
    Args:
        plot_name: Name of the plot function
        sequence_length: Length of sequence in bp
        has_clusters: Whether clusters were detected
    
    Returns:
        True if plot should be displayed, False otherwise
    """
    mapping = get_plot_mapping()
    canonical_name = mapping.get(plot_name, plot_name)
    
    # Hidden plots never shown
    if canonical_name == 'HIDDEN':
        return False
    
    # Size-dependent plots
    if canonical_name == 'manhattan_plot':
        return sequence_length > 50000
    elif canonical_name == 'linear_track':
        return sequence_length <= 50000
    
    # Conditional plots
    if canonical_name == 'cluster_size_distribution':
        return has_clusters
    
    # Default: show if not hidden
    return PlotDominance.is_plot_allowed(plot_name)

def get_nature_style_params() -> Dict[str, any]:
    """
    Get Nature journal style parameters for matplotlib.
    
    Returns:
        Dictionary of matplotlib rcParams
    """
    return {
        # Typography
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8,
        
        # Axes
        'axes.titlesize': 9,
        'axes.titleweight': 'bold',
        'axes.labelsize': 8,
        'axes.labelweight': 'normal',
        'axes.linewidth': 0.8,
        'axes.grid': False,
        'axes.spines.top': False,
        'axes.spines.right': False,
        
        # Ticks
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
        'figure.dpi': 72,
        'figure.facecolor': 'white',
        'savefig.dpi': 72,
        'savefig.format': 'pdf',
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05
    }
