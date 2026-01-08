"""
╔═════════════════════════════════════════════════════════════════��[...]
║                         NBDSCANNER WEB APPLICATION                            ║
║                    Non-B DNA Motif Detection System                          ║
╚═════════════════════════════════════════════════════════════════��[...]

AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2025.1
LICENSE: MIT

DESCRIPTION:
    Streamlit web application for comprehensive Non-B DNA motif detection.
    Provides interactive interface for sequence analysis and visualization.

FEATURES:
┌─────────────────────────────────────────────────────────────────��[...]
│  - Multi-FASTA support             - Real-time analysis progress            │
│  - 11 motif classes detection      - Interactive visualizations             │
│  - 22+ subclass analysis           - Export to CSV/BED/JSON                 │
│  - NCBI sequence fetch             - Publication-quality plots              │
└─────────────────────────────────────────────────────────────────��[...]

ARCHITECTURE:
    Input → Detection → Scoring → Overlap Resolution → Visualization → Export
"""

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import os  # Added for image path checking
import sys  # Added for system info
import psutil  # Added for memory monitoring
import numpy as np
import gc  # Added for memory management with large files
import tempfile  # For temporary file handling in downloads
import traceback  # For error handling and debugging
import re  # For filename sanitization
from collections import Counter  # For counting motif distributions

# Ensure the current directory is in the Python path for module imports
# This is needed for Streamlit Cloud deployment to find local modules
_current_dir = os.path.dirname(os.path.abspath(__file__))
if _current_dir not in sys.path:
    sys.path.insert(0, _current_dir)
# Import consolidated NBDScanner modules
from utilities import (
    parse_fasta, parse_fasta_chunked, parse_fasta_chunked_compressed, 
    get_file_preview, wrap, get_basic_stats, export_to_bed, export_to_csv,
    export_to_json, export_to_excel, export_statistics_to_excel, calculate_genomic_density, calculate_positional_density,
    export_results_to_dataframe, CORE_OUTPUT_COLUMNS, export_to_pdf,
    # Performance & Memory Management
    trigger_garbage_collection, optimize_dataframe_memory, get_memory_usage_mb,
    # Visualization functions (now consolidated in utilities.py)
    plot_motif_distribution, plot_coverage_map, plot_density_heatmap,
    plot_length_distribution, plot_score_distribution, plot_nested_pie_chart, 
    MOTIF_CLASS_COLORS, plot_density_comparison,
    plot_circos_motif_density,
    plot_density_comparison_by_subclass, plot_enrichment_analysis_by_subclass,
    plot_subclass_density_heatmap,
    # New Nature-quality visualizations
    plot_manhattan_motif_density, plot_cumulative_motif_distribution,
    plot_motif_cooccurrence_matrix, plot_gc_content_correlation,
    plot_linear_motif_track, plot_cluster_size_distribution,
    plot_motif_length_kde
)
from nonbscanner import (
    analyze_sequence, get_motif_info as get_motif_classification_info
)
# Import Nature-ready visualization standards
from visualization_standards import (
    NATURE_MOTIF_COLORS, PlotDominance, FigurePanel, MetricFilter,
    LabelPolicy, UILayout, TRANSPARENCY_NOTE, SUPPLEMENTARY_NOTE,
    should_show_plot, get_nature_style_params
)

# Job management modules
from job_manager import (
    generate_job_id, save_job_results, load_job_results, 
    job_exists, get_job_summary, list_all_jobs
)

# Try to import Entrez for demo functionality
try:
    from Bio import Entrez, SeqIO
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False

# Try to import Hyperscan (optional for acceleration)
import logging
logger = logging.getLogger(__name__)

HYPERSCAN_AVAILABLE = False
HYPERSCAN_VERSION = None
HYPERSCAN_ERROR = None

try:
    import hyperscan
    HYPERSCAN_AVAILABLE = True
    # Try to get version, but don't fail if not available
    try:
        HYPERSCAN_VERSION = hyperscan.__version__
    except AttributeError:
        HYPERSCAN_VERSION = 'unknown'
    logger.info(f"Hyperscan loaded successfully (version: {HYPERSCAN_VERSION})")
except ImportError as e:
    HYPERSCAN_ERROR = f"Hyperscan Python bindings not installed: {e}"
    logger.info(f"Hyperscan not available - using pure Python fallback. {HYPERSCAN_ERROR}")
except Exception as e:
    HYPERSCAN_ERROR = f"Hyperscan initialization failed: {e}"
    logger.warning(f"Hyperscan not available - using pure Python fallback. {HYPERSCAN_ERROR}")

# =============================================================================
# CENTRALIZED PAGE-WISE COLOR TOKENS
# =============================================================================
# All color values are defined here for consistent, tunable visual design.
# This enables theme evolution without touching page logic.
#
# ARCHITECTURE GUARANTEE:
#   - No pages, tabs, components, or workflows are modified
#   - Execution flow and logic remain completely untouched
#   - All changes are visual styling only
#
# COLOR TOKEN STRUCTURE:
#   1. Global Base Colors - Foundation colors used across all pages
#   2. Page-Specific Accent Palettes - Distinct identity per page/workflow
#   3. Semantic Status Colors - Consistent meaning (success, warning, error)
#   4. Visualization Palette - Scientific colorblind-friendly charts
#
# USAGE:
#   - Reference these tokens in COLOR_THEMES and UI_COMPONENT_STYLES
#   - For CSS: Use CSS custom properties (variables) injected by load_css()
#   - For inline HTML: Use get_page_colors() helper function
#   - Never use inline hex codes directly in page logic
#
# INLINE HTML ARCHITECTURAL NOTE:
#   Some UI elements use inline HTML styles (st.markdown with HTML) because:
#   - Streamlit's component system doesn't support all styling needs
#   - CSS variables don't work in inline HTML string literals
#   - Dynamic f-strings require literal color values
#   Solution: get_page_colors() function provides centralized token values
#   for injection into inline styles, maintaining single source of truth.
# =============================================================================

# ==================== GLOBAL BASE COLORS ====================
# Foundation colors used throughout the application
# These provide the core visual structure and hierarchy
GLOBAL_COLORS = {
    # Neutral backgrounds - Light mode foundation
    'neutral_50': '#FAFAFA',      # Lightest background for page base
    'neutral_100': '#F5F5F5',     # Very light background for sections
    'neutral_200': '#E5E5E5',     # Light borders and dividers
    'neutral_300': '#D4D4D4',     # Medium borders
    'neutral_400': '#A3A3A3',     # Disabled/muted elements
    'neutral_500': '#737373',     # Secondary text
    'neutral_600': '#525252',     # Body text
    'neutral_700': '#404040',     # Primary text
    'neutral_800': '#262626',     # Headers and emphasis
    'neutral_900': '#171717',     # Maximum contrast text
    
    # Dark mode backgrounds
    'dark_50': '#18181B',         # Darkest background for dark mode
    'dark_100': '#1F1F23',        # Very dark background
    'dark_200': '#27272A',        # Dark card background
    'dark_300': '#3F3F46',        # Dark borders
    'dark_400': '#52525B',        # Dark muted elements
    
    # White and black
    'white': '#FFFFFF',           # Pure white for cards and highlights
    'black': '#000000',           # Pure black for maximum contrast
}

# ==================== PAGE-SPECIFIC ACCENT PALETTES ====================
# Each page/workflow has a distinct color identity while maintaining consistency
# These colors provide subtle visual navigation cues

# HOME / OVERVIEW PAGE - Confident Scientific Blue (Trust, Clarity, Science)
# HIGHLY VIBRANT blue palette for eye-catching professional impact
HOME_COLORS = {
    'primary': '#0091FF',         # ELECTRIC BLUE - ultra vibrant primary actions
    'secondary': '#00B4FF',       # CYAN ELECTRIC - super bright hover states
    'accent': '#66D9FF',          # BRILLIANT SKY - vivid emphasis
    'light': '#CCF2FF',           # BRIGHT CYAN TINT - energetic backgrounds
    'lighter': '#E5F9FF',         # ULTRA LIGHT CYAN - vibrant page base
    'border': '#80E5FF',          # VIVID CYAN borders
    'text': '#003D82',            # DEEP VIBRANT BLUE - strong readability
    'shadow': 'rgba(0, 145, 255, 0.35)',  # STRONGER blue shadow for depth
}

# INPUT / UPLOAD PAGE - Fresh Natural Green (Growth, Initiation, Start)
# HIGHLY VIBRANT green palette for explosive fresh energy
INPUT_COLORS = {
    'primary': '#00E676',         # NEON GREEN - ultra vivid primary actions
    'secondary': '#1DE9B6',       # ELECTRIC MINT - brilliant hover states
    'accent': '#69F0AE',          # BRIGHT LIME - energetic emphasis
    'light': '#B9F6CA',           # VIVID MINT TINT - fresh backgrounds
    'lighter': '#E0FFF4',         # ULTRA LIGHT MINT - vibrant page base
    'border': '#7FFF9F',          # BRILLIANT GREEN borders
    'text': '#00612E',            # DEEP FOREST - strong readability
    'shadow': 'rgba(0, 230, 118, 0.35)',  # STRONGER green shadow for depth
}

# ANALYSIS / COMPUTATION PAGE - Energetic Orange (Energy, Processing, Activity)
# EXPLOSIVE orange palette for maximum energy and dynamism
ANALYSIS_COLORS = {
    'primary': '#FF6D00',         # BLAZING ORANGE - ultra bold primary actions
    'secondary': '#FF9100',       # ELECTRIC GOLD - brilliant hover states
    'accent': '#FFAB00',          # VIVID AMBER - striking emphasis
    'light': '#FFE57F',           # BRIGHT YELLOW TINT - energetic backgrounds
    'lighter': '#FFF9E6',         # ULTRA LIGHT GOLD - vibrant page base
    'border': '#FFCA28',          # BRILLIANT GOLD borders
    'text': '#BF360C',            # DEEP FLAME - strong readability
    'shadow': 'rgba(255, 109, 0, 0.4)',  # INTENSE orange shadow for depth
}

# RESULTS / TABLES PAGE - Refined Vibrant Purple (Insight, Data, Discovery)
# ELECTRIC purple palette for maximum visual impact and insight
RESULTS_COLORS = {
    'primary': '#D500F9',         # ELECTRIC PURPLE - ultra vivid primary actions
    'secondary': '#E040FB',       # NEON MAGENTA - brilliant hover states
    'accent': '#EA80FC',          # BRIGHT VIOLET - striking emphasis
    'light': '#F3E5F5',           # ELEGANT LAVENDER TINT - refined backgrounds
    'lighter': '#FCF2FF',         # ULTRA LIGHT VIOLET - vibrant page base
    'border': '#E1BEE7',          # VIVID LILAC borders
    'text': '#4A148C',            # DEEP ROYAL PURPLE - strong readability
    'shadow': 'rgba(213, 0, 249, 0.35)',  # INTENSE purple shadow for depth
}

# VISUALIZATION / PLOTS PAGE - Clinical Vibrant Teal (Precision, Clarity, Analysis)
# BRILLIANT teal palette for maximum clarity and visual precision
VISUALIZATION_COLORS = {
    'primary': '#00E5FF',         # NEON CYAN - ultra brilliant primary actions
    'secondary': '#18FFFF',       # ELECTRIC AQUA - super bright hover states
    'accent': '#84FFFF',          # VIVID SKY CYAN - striking emphasis
    'light': '#B2FFFF',           # BRIGHT AQUA TINT - energetic backgrounds
    'lighter': '#E0FFFF',         # ULTRA LIGHT CYAN - vibrant page base
    'border': '#76FFE7',          # BRILLIANT CYAN borders
    'text': '#004D5A',            # DEEP TEAL - strong readability
    'shadow': 'rgba(0, 229, 255, 0.4)',  # INTENSE cyan shadow for depth
}

# DOWNLOAD / EXPORT PAGE - Professional Vibrant Indigo (Completion, Authority, Final)
# BOLD indigo palette for commanding authority and completion
DOWNLOAD_COLORS = {
    'primary': '#536DFE',         # ELECTRIC INDIGO - ultra bold primary actions
    'secondary': '#5E72FF',       # NEON PERIWINKLE - brilliant hover states
    'accent': '#8C9EFF',          # BRIGHT LAVENDER BLUE - vivid emphasis
    'light': '#C5CAE9',           # VIVID INDIGO TINT - strong backgrounds
    'lighter': '#E8EAFF',         # ULTRA LIGHT INDIGO - vibrant page base
    'border': '#9FA8DA',          # BRILLIANT IRIS borders
    'text': '#1A237E',            # DEEP NAVY - strong readability
    'shadow': 'rgba(83, 109, 254, 0.4)',  # INTENSE indigo shadow for depth
}

# DOCUMENTATION PAGE - Vibrant Deep Purple Theme (Depth, Reference, Technical)
# Rich purple theme for technical documentation with high contrast
DOCUMENTATION_COLORS = {
    'primary': '#7C4DFF',         # Vivid electric purple - primary actions (Material Deep Purple A200)
    'secondary': '#B388FF',       # Bright lavender - hover states (Material Deep Purple A100)
    'accent': '#D1C4E9',          # Soft periwinkle - emphasis (Material Deep Purple 100)
    'light': '#0E1726',           # Dark background for documentation
    'lighter': '#0B1220',         # Darkest background for page base
    'border': '#1F2937',          # Dark borders
    'text': '#E5E7EB',            # Light text for dark background
    'shadow': 'rgba(0, 0, 0, 0.5)',  # Strong shadow for dark mode depth
}

# ==================== SEMANTIC STATUS COLORS ====================
# Consistent meaning across all pages - Universal visual language (ULTRA VIBRANT)
# EXPLOSIVE vibrant colors for immediate, unmistakable visual feedback
SEMANTIC_COLORS = {
    # Success states - Positive outcomes, completion, validation
    'success': '#00E676',         # NEON SUCCESS GREEN - ultra vivid
    'success_light': '#B9F6CA',   # BRIGHT success background
    'success_dark': '#00612E',    # DEEP success text
    'success_border': '#69F0AE',  # BRILLIANT success border
    
    # Warning states - Caution, important notices, attention needed
    'warning': '#FF9100',         # BLAZING WARNING ORANGE - maximum attention
    'warning_light': '#FFE57F',   # VIVID warning background
    'warning_dark': '#BF360C',    # DEEP warning text
    'warning_border': '#FFAB00',  # ELECTRIC warning border
    
    # Error states - Problems, failures, invalid inputs
    'error': '#FF1744',           # NEON ERROR RED - ultra striking
    'error_light': '#FFCDD2',     # BRIGHT error background
    'error_dark': '#B71C1C',      # DEEP error text
    'error_border': '#FF5252',    # BRILLIANT error border
    
    # Info states - Neutral information, tips, explanations
    'info': '#00B4FF',            # ELECTRIC INFO BLUE - ultra vivid
    'info_light': '#CCF2FF',      # BRIGHT info background
    'info_dark': '#003D82',       # DEEP info text
    'info_border': '#66D9FF',     # BRILLIANT info border
    
    # Progress states - Processing, loading, intermediate states
    'progress': '#E040FB',        # NEON PROGRESS PURPLE - ultra striking
    'progress_light': '#F3E5F5',  # Light progress background
    'progress_dark': '#4A148C',   # DEEP progress text
    'progress_border': '#EA80FC', # BRILLIANT progress border
}

# ==================== VISUALIZATION COLOR PALETTE ====================
# Scientific color scheme for charts and plots - ULTRA VIBRANT & ACCESSIBLE
# EXPLOSIVE colorblind-friendly palette with maximum saturation and visual impact
VISUALIZATION_PALETTE = {
    'chart_1': '#FF6D00',         # BLAZING ORANGE - Ultra high contrast
    'chart_2': '#0091FF',         # ELECTRIC BLUE - Ultra distinct
    'chart_3': '#00E676',         # NEON GREEN - Ultra clear
    'chart_4': '#FFEA00',         # BRILLIANT YELLOW - Ultra striking
    'chart_5': '#0043A8',         # DEEP ELECTRIC BLUE - Ultra professional
    'chart_6': '#FF1744',         # NEON RED - Ultra energetic
    'chart_7': '#FF00AA',         # ELECTRIC PINK - Ultra unique
    'chart_8': '#76FF03',         # VIVID LIME - Ultra natural
    'chart_9': '#D500F9',         # ELECTRIC PURPLE - Ultra elegant
    'chart_10': '#00E5FF',        # NEON CYAN - Ultra fresh
    'chart_11': '#FFC400',        # GOLD FLASH - Ultra warm
    'chart_12': '#546E7A',        # STEEL BLUE - Balanced contrast
}

# =============================================================================
# END OF CENTRALIZED COLOR TOKENS
# =============================================================================
# All color values are now defined above.
# The rest of the configuration below references these tokens.
# NO HARDCODED COLOR VALUES should appear beyond this point in the file.
# =============================================================================

# =============================================================================
# TUNABLE PARAMETERS CONFIGURATION
# =============================================================================
# All customizable parameters are centralized here for easy modification.
# 
# HOW TO USE:
#   1. Modify color themes to match your branding or preferences
#   2. Adjust font sizes and weights for better readability
#   3. Change visualization parameters (DPI, colors) for publication needs
#   4. Update UI text and labels to customize the interface
#   5. Control which motif classes appear in distributions and exports
#
# TIPS:
#   - Colors use hex format (e.g., '#4A90E2')
#   - Font sizes use rem units (1rem = browser default, usually 16px)
#   - DPI: Use 150 for screen, 300 for publication, 600 for high-quality print
#   - All changes take effect immediately on next page load/refresh
# =============================================================================

# ==================== TYPOGRAPHY & FONTS ====================
# Control all font settings for the application - MODERN & READABLE
# Optimized for modern high-resolution displays with excellent readability
FONT_CONFIG = {
    # Primary font families (in order of preference)
    # The browser will use the first available font in the list
    'primary_font': "'Inter', 'IBM Plex Sans', 'Segoe UI', system-ui, -apple-system, sans-serif",
    'monospace_font': "'JetBrains Mono', 'Fira Code', 'Consolas', monospace",
    
    # Font sizes (in rem units, where 1rem ≈ 16px in most browsers)
    # Enhanced sizes for modern, bold, research-quality appearance
    'h1_size': '2.75rem',     # Main page headers - bold, impactful
    'h2_size': '2.0rem',      # Section headers - clear hierarchy
    'h3_size': '1.5rem',      # Subsection headers - organized structure
    'h4_size': '1.25rem',     # Small headers - subtle distinction
    'body_size': '1.0rem',    # Body text, paragraphs - optimal readability
    'small_size': '0.9rem',   # Small text, notes - clear but compact
    'caption_size': '0.8rem', # Captions, footnotes - supporting information
    
    # Font weights (100-900, where 400 is normal and 700 is bold)
    'light_weight': 300,
    'normal_weight': 400,
    'medium_weight': 500,
    'semibold_weight': 600,
    'bold_weight': 700,
    'extrabold_weight': 800,
}

# ==================== COLOR THEMES ====================
# Define color themes for different moods and contexts - VIBRANT EDITION
# Each theme has primary, secondary, accent, backgrounds, text, and shadow colors
# All values now reference the centralized VIBRANT color tokens defined above
# IMPORTANT: No hardcoded hex values - all colors come from token system
COLOR_THEMES = {
    'scientific_blue': {
        'primary': HOME_COLORS['primary'],
        'secondary': HOME_COLORS['secondary'],
        'accent': HOME_COLORS['accent'],
        'bg_light': HOME_COLORS['lighter'],
        'bg_card': HOME_COLORS['light'],
        'text': HOME_COLORS['text'],
        'tab_bg': GLOBAL_COLORS['neutral_100'],
        'tab_active': HOME_COLORS['primary'],
        'shadow': 'rgba(0, 145, 255, 0.35)'  # VIBRANT blue shadow
    },
    'nature_green': {
        'primary': INPUT_COLORS['primary'],
        'secondary': INPUT_COLORS['secondary'],
        'accent': INPUT_COLORS['accent'],
        'bg_light': INPUT_COLORS['lighter'],
        'bg_card': INPUT_COLORS['light'],
        'text': INPUT_COLORS['text'],
        'tab_bg': GLOBAL_COLORS['neutral_100'],
        'tab_active': INPUT_COLORS['primary'],
        'shadow': 'rgba(0, 230, 118, 0.35)'  # VIBRANT green shadow
    },
    'genomic_purple': {
        'primary': RESULTS_COLORS['primary'],
        'secondary': RESULTS_COLORS['secondary'],
        'accent': RESULTS_COLORS['accent'],
        'bg_light': RESULTS_COLORS['lighter'],
        'bg_card': RESULTS_COLORS['light'],
        'text': RESULTS_COLORS['text'],
        'tab_bg': GLOBAL_COLORS['neutral_100'],
        'tab_active': RESULTS_COLORS['primary'],
        'shadow': 'rgba(213, 0, 249, 0.35)'  # VIBRANT purple shadow
    },
    'clinical_teal': {
        'primary': VISUALIZATION_COLORS['primary'],
        'secondary': VISUALIZATION_COLORS['secondary'],
        'accent': VISUALIZATION_COLORS['accent'],
        'bg_light': VISUALIZATION_COLORS['lighter'],
        'bg_card': VISUALIZATION_COLORS['light'],
        'text': VISUALIZATION_COLORS['text'],
        'tab_bg': GLOBAL_COLORS['neutral_100'],
        'tab_active': VISUALIZATION_COLORS['primary'],
        'shadow': 'rgba(0, 229, 255, 0.4)'  # VIBRANT cyan shadow
    },
    'midnight': {
        'primary': DOCUMENTATION_COLORS['primary'],
        'secondary': DOCUMENTATION_COLORS['secondary'],
        'accent': DOCUMENTATION_COLORS['accent'],
        'bg_light': DOCUMENTATION_COLORS['lighter'],
        'bg_card': DOCUMENTATION_COLORS['light'],
        'text': DOCUMENTATION_COLORS['text'],
        'tab_bg': GLOBAL_COLORS['dark_100'],
        'tab_active': DOCUMENTATION_COLORS['primary'],
        'shadow': 'rgba(0, 0, 0, 0.5)'  # Strong shadow for dark mode depth
    }
}

# ==================== PAGE THEMES PER TAB ====================
# Assign a specific theme to each tab for visual distinction
# This creates a unique color experience for each section of the app
# Change values to any theme name defined in COLOR_THEMES above
TAB_THEMES = {
    'Home': 'scientific_blue',          # Homepage uses professional blue
    'Upload & Analyze': 'nature_green',  # Upload tab uses natural green
    'Results': 'genomic_purple',        # Results tab uses genomic purple
    'Download': 'clinical_teal',        # Download tab uses clinical teal
    'Documentation': 'midnight'         # Documentation uses dark midnight theme
}

# ==================== VISUALIZATION PARAMETERS ====================
# Control all visualization and plot settings
# Modify these to customize the appearance of charts and graphs
VISUALIZATION_CONFIG = {
    # Plot resolution (DPI - dots per inch)
    # 150: Good for screen viewing (faster rendering)
    # 300: Publication quality (Nature, Science journals)
    # 600: High-quality print (posters, large format)
    'dpi': 300,
    
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

# ==================== UI TEXT CONTENT ====================
# All user-facing text and labels - CENTRALIZED FOR EASY CUSTOMIZATION
# Modify these to customize messages, titles, and help text
# Supports multi-language customization
UI_TEXT = {
    # ===== Application Metadata =====
    'app_title': 'NonBDNAFinder: A Comprehensive Non B-DNA forming motif detection system',
    'author': 'Dr. Venkata Rajesh Yella',
    'author_email': 'yvrajesh_bt@kluniversity.in',
    'github_profile': 'VRYella',
    'github_repo': 'NonBFinder',
    'github_url': 'https://github.com/VRYella/NonBFinder',
    'version': '2024.1',
    
    # ===== Page Titles =====
    'home_title': 'NonBDNA Motif Detection System',
    'upload_title': 'Sequence Upload and Motif Analysis',
    'results_title': 'Analysis Results and Visualization',
    'download_title': 'Export Data',
    'documentation_title': 'Scientific Documentation & References',
    
    # ===== Section Headers =====
    'section_scientific_foundation': 'Scientific Foundation',
    'section_motif_classes': 'Detected Motif Classes',
    'section_key_features': 'Key Features & Capabilities',
    'section_how_to_cite': 'How to Cite',
    'section_sequence_upload': 'Sequence Upload',
    'section_analysis_run': 'Analysis & Run',
    'section_quick_options': 'Quick Options',
    'section_analysis_summary': 'Analysis Summary',
    'section_visualizations': 'Visualizations',
    'section_export_options': 'Export Options',
    'section_export_preview': 'Export Preview',
    'section_download_files': 'Download Files',
    'section_additional_exports': 'Additional Exports',
    
    # ===== Home Page Content =====
    'home_system_status_hyperscan': 'Performance Mode: Hyperscan acceleration active for high-speed pattern matching',
    'home_system_status_standard': 'Standard Mode: Using regex-based pattern matching (all features fully functional)',
    'home_publication_ready': 'Publication Ready formats with 300 DPI resolution',
    'home_call_to_action_title': 'Ready to Analyze?',
    'home_call_to_action_text': 'Upload your FASTA sequences to begin comprehensive Non-B DNA motif detection',
    'home_call_to_action_button': '→ Go to "Upload & Analyze" tab',
    'home_image_caption': 'Non-B DNA Structural Diversity',
    'home_image_fallback_title': 'DNA',
    'home_image_fallback_subtitle': 'Non-B DNA Structures',
    'home_image_fallback_caption': 'Structural Diversity Database',
    
    # ===== Upload & Analyze Page =====
    'upload_input_method_prompt': 'Choose your input method:',
    'upload_method_file': 'Upload FASTA File',
    'upload_method_paste': 'Paste Sequence',
    'upload_method_example': 'Example Data',
    'upload_method_ncbi': 'NCBI Fetch',
    'upload_file_prompt': 'Drag and drop FASTA/multi-FASTA file here',
    'upload_file_help': 'Upload a FASTA file containing DNA sequences',
    'upload_processing': 'Processing',
    'upload_file_valid': 'Valid',
    'upload_preview_button': 'Preview Sequences',
    'upload_no_sequences': 'No sequences found in file.',
    'upload_paste_prompt': 'Paste single or multi-FASTA here:',
    'upload_paste_placeholder': 'Paste your DNA sequence(s) here...',
    'upload_paste_help': 'Paste DNA sequences in FASTA format',
    'upload_example_type_prompt': 'Example Type:',
    'upload_example_single': 'Single Example',
    'upload_example_multi': 'Multi-FASTA Example',
    'upload_example_help': 'Load example sequences for testing',
    'upload_load_single_button': 'Load Single Example',
    'upload_load_multi_button': 'Load Multi-FASTA Example',
    'upload_example_single_success': 'Single example sequence loaded.',
    'upload_example_multi_success': 'Multi-FASTA example loaded with {count} sequences.',
    'upload_ncbi_db_prompt': 'NCBI Database',
    'upload_ncbi_help': 'Only nucleotide and gene databases are applicable for DNA motif analysis',
    'upload_ncbi_query_prompt': 'Enter query (accession, gene, etc.):',
    'upload_ncbi_max_records': 'Max Records',
    'upload_ncbi_fetch_button': 'Fetch from NCBI',
    'upload_ncbi_success': 'Fetched {count} sequences.',
    'upload_ncbi_error': 'NCBI fetch failed: {error}',
    'upload_ncbi_empty_warning': 'Enter a query before fetching.',
    'upload_quick_options_note': 'All 11 motif classes with 22+ subclasses are detected automatically',
    'upload_parallel_note': 'Parallel scanner works best on sequences >100kb with multiple CPU cores',
    'upload_run_analysis_button': 'Run Complete Motif Analysis',
    'upload_no_sequences_error': 'Please upload or input sequences before running analysis.',
    
    # ===== Analysis Progress Messages =====
    'progress_validating': 'Validating results for consistency and quality...',
    'progress_generating_viz': 'Generating comprehensive visualizations for all classes and subclasses...',
    'progress_all_detectors': 'All detectors process in parallel | Followed by overlap resolution & clustering',
    'progress_parallel_fallback': 'Parallel scanner failed, falling back to standard: {error}',
    
    # ===== Status Messages =====
    'status_no_results': 'No analysis results. Please run motif analysis first.',
    'status_analysis_ready': 'Ready to analyze',
    'status_analysis_running': 'Analysis in progress...',
    'status_analysis_complete': 'Analysis complete!',
    'status_validation_passed': 'Validation passed: No consistency issues found',
    'status_validation_issues': 'Validation found {count} potential issues:',
    'status_viz_prepared': 'All visualizations prepared: {count} components in {time:.2f}s',
    'status_pre_generated_ready': 'Pre-generated Analysis Ready: {classes} unique classes, {subclasses} unique subclasses analyzed',
    
    # ===== Results Page =====
    'results_performance_title': 'Performance Metrics',
    'results_no_motifs': 'No motifs detected for this sequence.',
    'results_sequence_selector': 'Choose Sequence for Details:',
    'results_hybrid_cluster_info': '{count} Hybrid/Cluster motifs detected. View them in the \'Cluster/Hybrid\' tab below.',
    'results_columns_display': 'Columns to display',
    'results_page_label': 'Page (showing {rows} rows per page)',
    'results_page_info': 'Showing motifs {start} to {end} of {total}',
    'results_viz_title': 'Visualizations',
    
    # ===== Visualization Tab Names =====
    'viz_tab_distribution': 'Distribution & Statistics',
    'viz_tab_coverage': 'Coverage & Density',
    'viz_tab_genome_wide': 'Genome-Wide Analysis',
    'viz_tab_hybrid_cluster': 'Cluster/Hybrid',
    
    # ===== Visualization Section Headers =====
    'viz_motif_distribution': 'Motif Distribution',
    'viz_statistical_analysis': 'Statistical Analysis',
    'viz_density_metrics': 'Density Metrics',
    'viz_distributions': 'Distributions',
    'viz_sequence_coverage': 'Sequence Coverage Analysis',
    'viz_circos_density': 'Circos Density Plot',
    'viz_genome_wide_title': 'Genome-Wide Motif Analysis',
    'viz_genome_wide_desc': 'Publication-quality genome-scale visualizations showing motif distribution patterns across the entire sequence.',
    'viz_manhattan_plot': 'Manhattan Plot - Motif Density Hotspots',
    'viz_cumulative_dist': 'Cumulative Motif Distribution',
    'viz_linear_track': 'Linear Motif Track Viewer',
    'viz_linear_track_info': 'Linear motif track is shown for sequences < 50kb. For large sequences, use Manhattan plot or Coverage map.',
    'viz_hybrid_cluster_title': 'Hybrid & Cluster Motif Analysis',
    'viz_advanced_title': 'Advanced Statistical Visualizations',
    'viz_advanced_desc': 'Advanced publication-quality visualizations for in-depth analysis and manuscript figures.',
    'viz_cooccurrence': 'Motif Co-occurrence Matrix',
    'viz_cooccurrence_desc': 'Shows which motif classes tend to appear together (overlapping or within 1bp)',
    'viz_gc_correlation': 'GC Content vs Motif Density Correlation',
    'viz_gc_correlation_desc': 'Scatter plot showing relationship between GC content and motif density',
    'viz_length_kde': 'Motif Length Distribution (Kernel Density)',
    'viz_length_kde_desc': 'Smooth probability density curves showing length patterns by class',
    'viz_cluster_size': 'Cluster Size & Diversity Distribution',
    'viz_cluster_size_desc': 'Distribution of motif counts and class diversity within clusters',
    
    # ===== Analysis Section Headers =====
    'analysis_section_title': 'Analysis & Run',
    'analysis_quick_options_title': 'Quick Options',
    'analysis_run_button': 'Run NBDScanner Analysis',
    'analysis_run_button_disabled': 'Run NBDScanner Analysis (Disabled)',
    'analysis_run_button_disabled_note': 'Please upload or paste a valid sequence first',
    'analysis_pipeline_title': 'Analysis Pipeline',
    'analysis_progress_title': 'NonBFinder Analysis',
    
    # ===== Hybrid/Cluster Specific =====
    'hybrid_all_tab': 'All',
    'hybrid_motifs_tab': 'Hybrid Motifs',
    'hybrid_clusters_tab': 'Cluster Motifs',
    'hybrid_motifs_title': 'Hybrid Motifs (Overlapping Different Classes)',
    'hybrid_motifs_info': 'Hybrid motifs are regions where different Non-B DNA motif classes overlap. Found {count} hybrid regions.',
    'hybrid_clusters_title': 'Non-B DNA Clusters (High-Density Regions)',
    'hybrid_clusters_info': 'DNA Clusters are high-density regions with multiple Non-B DNA motif classes. Found {count} cluster regions.',
    'hybrid_no_motifs': 'No hybrid or cluster motifs detected in this sequence. Hybrid motifs occur when different Non-B DNA classes overlap, and clusters form when multiple motifs are found in close proximity.',
    'hybrid_no_hybrid': 'No hybrid motifs detected in this sequence.',
    'hybrid_no_clusters': 'No DNA clusters detected in this sequence.',
    
    # ===== Download Page =====
    'download_no_results': 'No results available to download.',
    'download_config_used': 'Analysis Configuration Used:',
    'download_config_overlap': 'Overlap Handling: {option}',
    'download_config_classes': 'Motif Classes: {count} classes selected',
    'download_config_total': 'Total Motifs Found: {count}',
    'download_export_config': 'Export Configuration',
    'download_include_sequences': 'Include Full Sequences',
    'download_include_help': 'Include full motif sequences in export',
    'download_excluded_info': '{count} Hybrid/Cluster motifs are excluded from downloads based on configuration. These are shown only in the Cluster/Hybrid visualization tab.',
    'download_included_info': 'All {count} Hybrid/Cluster motifs are included in downloads.',
    'download_preview_caption': 'Showing first 10 of {count} total records{excluded}',
    'download_excel_info_title': 'Excel Format',
    'download_excel_info': 'Downloads a multi-sheet workbook with:\n• Consolidated Sheet: All non-overlapping motifs with core columns\n• Class Sheets: Separate sheets for each motif class with ALL detailed columns\n• Subclass Sheets: Detailed breakdown by subclass with ALL detailed columns',
    'download_excel_note': 'Additional detailed columns are only shown in per-motif-class/subclass sheets, not in display tables.',
    'download_button_csv': 'Download CSV',
    'download_button_excel': 'Download Excel',
    'download_button_json': 'Download JSON',
    'download_button_bed': 'Download BED',
    'download_button_config': 'Download Config',
    'download_help_csv': 'Comma-separated values format',
    'download_help_excel': 'Excel format with multiple sheets (Consolidated + per-class)',
    'download_help_json': 'JSON format with metadata',
    'download_help_bed': 'BED format for genome browsers',
    'download_help_config': 'Analysis configuration and metadata',
    'download_error_excel': 'Excel export error: {error}',
    
    # ===== Documentation Page =====
    'doc_motif_classes_title': 'Motif Classes Detected:',
    'doc_references_title': 'References:',
    'doc_config_details': 'Scoring Configuration Details',
    'doc_length_constraints': 'Motif Length Constraints',
    'doc_scoring_methods': 'Scoring Methods',
    'doc_developed_by': 'Developed by',
    
    # ===== Tooltips and Help Text =====
    'help_detailed_analysis': 'Include comprehensive motif metadata in results',
    'help_quality_validation': 'Validate detected motifs for consistency',
    'help_parallel_scanner': 'Enable experimental parallel chunk-based scanner (>100kb sequences)',
    'help_chunk_progress': 'Display detailed progress for each processing chunk',
    'help_sequence_selector': 'Select a sequence to view detailed analysis results',
    
    # ===== Metric Labels =====
    'metric_sequence_coverage': 'Sequence Coverage',
    'metric_motif_density': 'Motif Density<br>(motifs/kb)',
    'metric_total_motifs': 'Total Motifs',
    'metric_sequence_length': 'Sequence Length (bp)',
    'metric_processing_time': 'Processing Time',
    'metric_base_pairs': 'Base Pairs',
    'metric_speed': 'bp/second',
    'metric_detector_processes': 'Detector Processes',
    'metric_sequences': 'Sequences',
    'metric_total_motifs_found': 'Total Motifs',
    'metric_hybrid_motifs': 'Hybrid Motifs',
    'metric_dna_clusters': 'DNA Clusters',
    'metric_avg_length': 'Avg Length (bp)',
    'metric_total': 'Total',
    
    # ===== Column Names (for display) =====
    'col_motif_class': 'Motif Class',
    'col_motif_subclass': 'Motif Subclass',
    'col_genomic_density': 'Genomic Density (%)',
    'col_motifs_per_kbp': 'Motifs/kbp',
    'col_sequence_name': 'Sequence Name',
    'col_source': 'Source',
    'col_class': 'Class',
    'col_subclass': 'Subclass',
    'col_start': 'Start',
    'col_end': 'End',
    'col_length': 'Length',
    'col_sequence': 'Sequence',
    'col_score': 'Score',
    
    # ===== Generic UI Elements =====
    'button_load': 'Load',
    'button_fetch': 'Fetch',
    'button_analyze': 'Analyze',
    'button_download': 'Download',
    'button_export': 'Export',
    'button_preview': 'Preview',
    'label_success': 'Success',
    'label_error': 'Error',
    'label_warning': 'Warning',
    'label_info': 'Info',
    'label_tip': 'Tip',
    'label_note': 'Note',
    'label_progress': 'Progress',
    'label_complete': 'Complete',
    'label_ready': 'Ready',
    'label_loading': 'Loading',
    'label_processing': 'Processing',
    'label_validating': 'Validating',
    'label_valid': 'Valid',
    
    # ===== Analysis Messages =====
    'analysis_no_sequences_warning': 'No sequences found.',
    'analysis_all_detectors_parallel': 'All detectors process in parallel | Followed by overlap resolution & clustering',
    
    # ===== Tooltips =====
    'tooltip_detailed_analysis': 'Include comprehensive motif metadata',
    'tooltip_quality_validation': 'Validate detected motifs',
    'tooltip_chunk_progress': 'Display detailed progress for each processing chunk',
    'tooltip_parallel_scanner': 'Enable experimental parallel chunk-based scanner (>100kb sequences)',
    
    # ===== Headings for Results Page =====
    'heading_results_viz': 'Visualizations',
    'heading_analysis_results': 'Analysis Results and Visualization',
    'heading_analysis_summary': 'Analysis Summary',
    'heading_motif_distribution': 'Motif Distribution',
    'heading_statistical_analysis': 'Statistical Analysis',
    'heading_density_metrics': 'Density Metrics',
    'heading_class_level': 'Class Level Analysis',
    'heading_subclass_level': 'Subclass Level Analysis',
    'heading_distributions': 'Distributions',
    'heading_sequence_coverage': 'Sequence Coverage Analysis',
    'heading_circos_density': 'Circos Density Plot',
    'heading_genome_wide': 'Genome-Wide Motif Analysis',
    'heading_hybrid_cluster': 'Hybrid & Cluster Motif Analysis',
    'heading_advanced_viz': 'Advanced Statistical Visualizations',
    
    # ===== Headings for Analysis Section =====
    'heading_analysis_run': 'Analysis & Run',
    
    # ===== Export Page Headings =====
    'heading_export_data': 'Export Data',
    'heading_export_options': 'Export Options',
    'heading_export_preview': 'Export Preview',
    'heading_download_files': 'Download Files',
    'heading_export_config': 'Export Configuration',
    'heading_additional_exports': 'Additional Exports',
    
    # ===== Documentation Page Headings =====
    'heading_documentation': 'Scientific Documentation & References',
    'heading_motif_classes': 'Motif Classes Detected:',
    'heading_references': 'References:',
    'heading_config_details': 'Scoring Configuration Details',
    'heading_length_constraints': 'Motif Length Constraints',
    'heading_scoring_methods': 'Scoring Methods',
}

# ==================== LAYOUT & SPACING ====================
# Control spacing, padding, and visual structure
LAYOUT_CONFIG = {
    # Streamlit layout mode
    # 'wide': Uses full browser width (recommended for data apps)
    # 'centered': Centers content with maximum width
    'layout_mode': 'wide',
    
    # Border radius values (in pixels) for rounded corners
    'border_radius': {
        'small': '8px',    # Small elements (buttons, tags)
        'medium': '12px',  # Medium elements (cards, inputs)
        'large': '16px',   # Large elements (panels, sections)
        'pill': '50px',    # Fully rounded (pills, progress bars)
    },
    
    # Padding and margins (in rem units)
    'padding': {
        'small': '0.5rem',   # Tight spacing
        'medium': '1rem',    # Standard spacing
        'large': '1.5rem',   # Generous spacing
        'xlarge': '2rem',    # Extra spacing for sections
    },
    
    # Margins (in rem units)
    'margin': {
        'none': '0',
        'small': '0.5rem',
        'medium': '1rem',
        'large': '1.5rem',
        'xlarge': '2rem',
    },
    
    # Column gaps for multi-column layouts
    # Options: 'small', 'medium', 'large'
    'column_gap': 'large',
    
    # Card styling for panels and containers
    'card_shadow': '0 4px 20px rgba(0,0,0,0.08)',  # Soft shadow for depth
    'card_border': '1px solid #e5e7eb',             # Subtle border
    'card_shadow_hover': '0 8px 30px rgba(0,0,0,0.12)',  # Elevated shadow on hover
    
    # Section spacing
    'section_spacing': '2rem',      # Space between major sections
    'subsection_spacing': '1.5rem', # Space between subsections
    'element_spacing': '1rem',      # Space between UI elements
    
    # Container widths
    'container_max_width': '1400px',  # Maximum width for centered content
    'sidebar_width': '300px',          # Sidebar width
    
    # Z-index layers (for stacking context)
    'z_index': {
        'base': 1,
        'dropdown': 100,
        'sticky': 200,
        'modal': 1000,
        'tooltip': 2000,
    },
}

# ==================== ANALYSIS PARAMETERS ====================
# Control sequence processing and analysis behavior
ANALYSIS_CONFIG = {
    # Sequence processing thresholds
    'chunk_threshold': 10_000,     # Sequences > this size use chunking (bp)
    'default_chunk_size': 10_000,  # Default chunk size for large sequences (bp)
    'default_chunk_overlap': 500,  # Overlap between chunks to catch motifs at boundaries (bp)
    
    # Performance and display settings
    'max_sequences_preview': 3,    # Number of sequences to show in file preview
    'rows_per_page': 100,          # Pagination size for large result tables
    'update_interval': 5,          # Progress update frequency (in sequences)
    
    # Motif filtering and display
    'min_score_threshold': 0.0,    # Minimum score to display (0 = show all)
    
    # *** IMPORTANT: Control which motifs appear in distributions ***
    # Set to True to include Hybrid and Cluster motifs in all visualizations
    # Set to False to exclude them from distribution plots (they'll still appear in dedicated tab)
    'include_hybrid_in_distribution': True,   # Include Hybrid motifs in plots
    'include_clusters_in_distribution': True, # Include Cluster motifs in plots
    
    # File upload limits
    'max_file_size_mb': 1024,      # Maximum file size in MB (1 GB default)
}

# ==================== EXPORT FORMATS ====================
# Control data export options and default settings
EXPORT_CONFIG = {
    # Available export formats
    # Add or remove formats as needed
    'available_formats': ['CSV', 'Excel', 'JSON', 'BED'],
    
    # Default export options (can be overridden by user)
    'include_sequences': True,     # Include full motif sequences in exports
    'excel_multi_sheet': True,     # Create multi-sheet Excel workbooks (one per motif class)
    'json_pretty': True,           # Pretty-print JSON exports (more readable)
    
    # Column selection for exports
    # These are the core columns that appear in all export formats
    'core_columns': [
        'Sequence_Name',  # Name or accession of the sequence
        'Source',         # Source database or experiment
        'Class',          # Motif class (G4, Z-DNA, etc.)
        'Subclass',       # Motif subclass or subtype
        'Start',          # Start position (bp)
        'End',            # End position (bp)
        'Length',         # Motif length (bp)
        'Sequence',       # DNA sequence of motif
        'Score',          # Confidence score (0-3 scale)
    ],
}

# ==================== PERFORMANCE MONITORING ====================
# Control performance tracking and system monitoring
PERFORMANCE_CONFIG = {
    'enable_monitoring': True,  # Enable performance metrics collection
    'show_system_info': True,   # Show system resource information in logs
    'log_timing': True,         # Log timing information for analysis steps
}

# ==================== UI COMPONENT STYLES ====================
# Centralized styling parameters for all UI components
# These values are used throughout the application for consistent styling
UI_COMPONENT_STYLES = {
    # Button styles - Using centralized colors
    'button': {
        'primary_bg': 'linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%)',
        'primary_color': GLOBAL_COLORS['white'],
        'padding': '0.5em 1.2em',
        'font_size': '0.92rem',
        'font_weight': 500,
        'border_radius': '12px',
        'shadow': '0 4px 12px rgba(0,0,0,0.15)',
        'shadow_hover': '0 6px 16px rgba(0,0,0,0.2)',
    },
    
    # Input field styles - Using centralized colors
    'input': {
        'bg': GLOBAL_COLORS['white'],
        'border': f"1.5px solid {GLOBAL_COLORS['neutral_200']}",
        'border_focus': '1.5px solid var(--primary-color)',
        'border_radius': '12px',
        'padding': '0.7rem 1rem',
        'font_size': '0.95rem',
        'shadow': '0 1px 3px rgba(0, 0, 0, 0.05)',
        'shadow_focus': '0 0 0 3px var(--shadow-color), 0 2px 6px rgba(0, 0, 0, 0.08)',
    },
    
    # Card/Panel styles - Using centralized colors
    'card': {
        'bg': GLOBAL_COLORS['white'],
        'border': f"1px solid {GLOBAL_COLORS['neutral_200']}",
        'border_radius': '16px',
        'padding': '1.5rem',
        'shadow': '0 2px 12px rgba(0,0,0,0.08)',
        'shadow_hover': '0 4px 20px rgba(0,0,0,0.12)',
    },
    
    # Alert/Notification styles - Now using centralized SEMANTIC_COLORS
    'alert': {
        'success_bg': SEMANTIC_COLORS['success_light'],
        'success_border': SEMANTIC_COLORS['success'],
        'info_bg': SEMANTIC_COLORS['info_light'],
        'info_border': SEMANTIC_COLORS['info'],
        'warning_bg': SEMANTIC_COLORS['warning_light'],
        'warning_border': SEMANTIC_COLORS['warning'],
        'error_bg': SEMANTIC_COLORS['error_light'],
        'error_border': SEMANTIC_COLORS['error'],
        'border_width': '4px',
        'border_radius': '12px',
        'padding': '1rem',
    },
    
    # Table/DataFrame styles - Using centralized colors
    'table': {
        'header_bg': 'linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%)',
        'header_color': GLOBAL_COLORS['white'],
        'header_font_weight': 600,
        'header_padding': '0.8rem',
        'row_padding': '0.7rem',
        'row_hover_bg': GLOBAL_COLORS['neutral_50'],
        'alt_row_bg': 'rgba(238, 242, 255, 0.5)',
        'border': '1px solid var(--bg-card)',
        'border_radius': '12px',
    },
    
    # Tab styles - Using centralized colors
    'tabs': {
        'bg': 'var(--tab-bg)',
        'active_bg': 'linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%)',
        'active_color': GLOBAL_COLORS['white'],
        'inactive_color': GLOBAL_COLORS['neutral_600'],
        'hover_bg': 'rgba(255, 255, 255, 0.7)',
        'border_radius': '50px',  # Pill-shaped
        'padding': '8px 12px',
        'font_size': '0.95rem',
        'font_weight': 500,
    },
    
    # Progress bar styles
    'progress': {
        'bg': 'var(--bg-card)',
        'fill_bg': 'linear-gradient(90deg, var(--primary-color) 0%, var(--secondary-color) 100%)',
        'height': '8px',
        'border_radius': '50px',
        'shadow': '0 1px 4px rgba(0,0,0,0.1)',
    },
    
    # Metric card styles - Using centralized colors
    'metric': {
        'bg': GLOBAL_COLORS['white'],
        'border': f"1.5px solid {GLOBAL_COLORS['neutral_200']}",
        'border_radius': '16px',
        'padding': '0.8rem',
        'value_size': '1.5rem',
        'value_weight': 700,
        'label_size': '0.8rem',
        'label_weight': 500,
        'shadow': '0 2px 8px rgba(0,0,0,0.08)',
    },
    
    # File uploader styles - Using centralized colors
    'file_uploader': {
        'border': '2px dashed var(--accent-color)',
        'border_hover': '2px dashed var(--primary-color)',
        'bg': GLOBAL_COLORS['neutral_50'],
        'bg_hover': HOME_COLORS['light'],
        'border_radius': '16px',
        'padding': '1.2rem',
    },
    
    # Popover styles - Using centralized colors
    'popover': {
        'bg': GLOBAL_COLORS['white'],
        'border': f"1.5px solid {GLOBAL_COLORS['neutral_200']}",
        'border_radius': '16px',
        'padding': '1rem',
        'shadow': '0 4px 20px rgba(0,0,0,0.15)',
    },
    
    # Expander/Accordion styles - Using centralized colors
    'expander': {
        'header_bg': GLOBAL_COLORS['neutral_50'],
        'header_bg_hover': 'var(--bg-card)',
        'content_bg': GLOBAL_COLORS['white'],
        'border': f"1.5px solid {GLOBAL_COLORS['neutral_200']}",
        'border_radius': '10px',
        'padding': '0.6rem 1rem',
        'shadow': '0 1px 4px rgba(0,0,0,0.08)',
    },
    
    # Scrollbar styles - Using centralized colors
    'scrollbar': {
        'width': '8px',
        'track_bg': GLOBAL_COLORS['neutral_100'],
        'thumb_bg': 'linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%)',
        'border_radius': '50px',
    },
    
    # Tooltip styles - Using centralized colors
    'tooltip': {
        'bg': 'rgba(0, 0, 0, 0.9)',
        'color': GLOBAL_COLORS['white'],
        'font_size': '0.85rem',
        'padding': '0.5rem 0.75rem',
        'border_radius': '8px',
        'shadow': '0 2px 8px rgba(0,0,0,0.2)',
    },
}

# ==================== ANIMATION & TRANSITION SETTINGS ====================
# Control animation timing and effects for smooth UI interactions
ANIMATION_CONFIG = {
    # Transition durations (in seconds)
    'transition_fast': '0.15s',
    'transition_normal': '0.2s',
    'transition_slow': '0.3s',
    
    # Easing functions
    'easing_smooth': 'ease',
    'easing_in': 'ease-in',
    'easing_out': 'ease-out',
    'easing_in_out': 'ease-in-out',
    
    # Animation names (CSS keyframes defined in styles.css)
    'fade_in': 'fade-in',
    'pulse': 'pulse-dot',
    'shimmer': 'progress-shimmer',
    
    # Enable/disable animations globally
    'enable_animations': True,
    'reduce_motion': False,  # Respect user's prefers-reduced-motion setting
}

# =============================================================================
# SCIENTIFIC TIME FORMATTING UTILITIES
# =============================================================================
# Centralized time formatting for consistent, publication-quality time display
# All elapsed-time displays use the canonical format: HH:MM:SS › mmm
# This provides reviewer-grade precision while remaining human-readable
# =============================================================================

def format_time_scientific(seconds: float) -> str:
    """
    Format elapsed time in simple MM:SS format.
    
    This format provides:
    - Human-readable minutes and seconds
    - No hours or microseconds (simplified display)
    - Consistent display across all workflows
    
    Args:
        seconds: Elapsed time in seconds (float)
        
    Returns:
        Formatted time string (e.g., "02:15" or "125:32")
        
    Examples:
        >>> format_time_scientific(0.234)
        "00:00"
        >>> format_time_scientific(135.678)
        "02:15"
        >>> format_time_scientific(5432.123)
        "90:32"
    """
    minutes = int(seconds // 60)
    secs = int(seconds % 60)
    
    return f"{minutes:02d}:{secs:02d}"


def format_time_compact(seconds: float) -> str:
    """
    Format elapsed time in MM:SS format for compact displays.
    
    Simple minutes:seconds format for all durations.
    
    Args:
        seconds: Elapsed time in seconds (float)
        
    Returns:
        Formatted time string (e.g., "02:15" or "125:32")
    """
    minutes = int(seconds // 60)
    secs = int(seconds % 60)
    return f"{minutes:02d}:{secs:02d}"


# ---------- CACHING FUNCTIONS (Memory-Efficient) ----------
@st.cache_resource(show_spinner=False)
def cache_genome_as_numpy(sequence: str) -> np.ndarray:
    """
    Cache genome sequence as NumPy byte array for memory efficiency.
    
    This prevents reloading large genomes and reduces memory footprint
    when using Streamlit on free tier (1GB limit).
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        NumPy array of sequence bytes
    """
    return np.frombuffer(sequence.encode('utf-8'), dtype=np.uint8)


@st.cache_resource(show_spinner=False)
def cache_hyperscan_database(_patterns: list = None):
    """
    Cache compiled Hyperscan database for pattern matching.
    
    The underscore prefix on _patterns parameter is used by Streamlit
    to indicate the parameter should not be hashed for caching purposes.
    This prevents Streamlit from attempting to hash complex pattern objects.
    
    Args:
        _patterns: List of (pattern, pattern_id) tuples to compile
        
    Returns:
        Compiled Hyperscan database or None
    """
    if not HYPERSCAN_AVAILABLE or _patterns is None or len(_patterns) == 0:
        logger.debug("Hyperscan not available or no patterns provided")
        return None
    
    try:
        # Compile patterns into Hyperscan database
        logger.debug(f"Compiling Hyperscan database with {len(_patterns)} patterns...")
        
        expressions = []
        ids = []
        flags = []
        
        for pattern, pattern_id in _patterns:
            # Encode pattern with error handling
            # Note: DNA patterns should normally be ASCII (ATGC), but we provide
            # UTF-8 fallback for robustness in case patterns contain metadata or
            # special characters. A warning is logged to help identify data quality issues.
            try:
                pattern_bytes = pattern.encode('ascii')
            except UnicodeEncodeError:
                # Fall back to UTF-8 if ASCII fails
                pattern_bytes = pattern.encode('utf-8')
                logger.warning(f"Pattern {pattern_id} contains non-ASCII characters (expected ATGC). Using UTF-8 encoding.")
            
            expressions.append(pattern_bytes)
            ids.append(pattern_id)
            # Use CASELESS and DOTALL flags for DNA matching
            flags.append(hyperscan.HS_FLAG_CASELESS | hyperscan.HS_FLAG_DOTALL)
        
        db = hyperscan.Database()
        db.compile(
            expressions=expressions,
            ids=ids,
            elements=len(expressions),
            flags=flags
        )
        
        logger.info(f"Successfully compiled Hyperscan database with {len(expressions)} patterns")
        return db
        
    except Exception as e:
        logger.error(f"Hyperscan database compilation failed: {e}")
        st.warning(f"Hyperscan database compilation failed: {e}. Falling back to regex matching.")
        return None


@st.cache_data(show_spinner=False, max_entries=10, ttl=3600)
def cache_analysis_results(sequence_hash: str, sequence: str, name: str):
    """
    Cache analysis results for sequences to avoid re-computation.
    Uses sequence hash for efficient cache key lookup.
    
    Args:
        sequence_hash: Hash of sequence for cache key
        sequence: DNA sequence string
        name: Sequence name
        
    Returns:
        List of detected motifs
    """
    return analyze_sequence(sequence, name)


@st.cache_data(show_spinner=False, max_entries=20, ttl=3600)
def get_cached_stats(sequence: str, motifs_json: str):
    """
    Cache statistics calculation for sequences.
    
    Args:
        sequence: DNA sequence
        motifs_json: JSON string of motifs (for cache key)
        
    Returns:
        Dictionary of sequence statistics
    """
    import json
    motifs = json.loads(motifs_json) if motifs_json else []
    return get_basic_stats(sequence, motifs)


def get_system_info():
    """
    Get current system memory and resource information.
    Uses psutil if available, otherwise provides basic info.
    
    Returns:
        Dictionary with memory and system info
    """
    try:
        # Get memory info
        memory = psutil.virtual_memory()
        
        return {
            'memory_total_mb': memory.total / (1024 * 1024),
            'memory_used_mb': memory.used / (1024 * 1024),
            'memory_available_mb': memory.available / (1024 * 1024),
            'memory_percent': memory.percent,
            'cpu_count': psutil.cpu_count(),
            'available': True
        }
    except (ImportError, AttributeError):
        # psutil not available or doesn't support this platform
        return {
            'memory_total_mb': 0,
            'memory_used_mb': 0,
            'memory_available_mb': 0,
            'memory_percent': 0,
            'cpu_count': 1,
            'available': False
        }


# ---------- PAGE CONFIG ----------
st.set_page_config(
    page_title=f"{UI_TEXT['app_title']} - Non-B DNA Motif Finder",
    layout=LAYOUT_CONFIG['layout_mode'],
    page_icon=None,
    menu_items={'About': f"NBDScanner | Developed by {UI_TEXT['author']}"}
)

def format_sequence_limit():
    """Format the sequence limit for display - now shows 'unlimited' since limit is removed"""
    return "unlimited (chunked processing enabled)"

# Get motif classification info
CLASSIFICATION_INFO = get_motif_classification_info()

# ---------- PATCH: Ensure every motif has Subclass ----------
def ensure_subclass(motif):
    """Guarantee every motif has a string 'Subclass'"""
    if isinstance(motif, dict):
        if 'Subclass' not in motif or motif['Subclass'] is None:
            motif['Subclass'] = motif.get('Subtype', 'Other')
        return motif
    else:
        # Handle non-dict motifs gracefully
        return {'Subclass': 'Other', 'Motif': motif}


# ---------- HELPER: Format time for display ----------
def format_time(seconds):
    """Format time in seconds to a human-readable string.
    
    Args:
        seconds: Time in seconds (float or int)
        
    Returns:
        Formatted string (e.g., "45.3s", "12m 30s", "2h 15m")
    
    Examples:
        >>> format_time(45.3)
        '45.3s'
        >>> format_time(750)
        '12m 30s'
        >>> format_time(7800)
        '2h 10m'
    """
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        mins = int(seconds // 60)
        secs = int(seconds % 60)
        return f"{mins}m {secs}s"
    else:
        hours = int(seconds // 3600)
        mins = int((seconds % 3600) // 60)
        return f"{hours}h {mins}m"


# ---------- HELPER: Generate Excel data as bytes for download ----------
def generate_excel_bytes(motifs, simple_format=True):
    """
    Generate Excel file as bytes for Streamlit download button.
    
    Args:
        motifs: List of motif dictionaries
        simple_format: If True, use 2-tab format (NonOverlappingConsolidated, OverlappingAll)
        
    Returns:
        bytes: Excel file data as bytes
    """
    import tempfile
    
    # Create a temporary file
    with tempfile.NamedTemporaryFile(mode='wb', suffix='.xlsx', delete=False) as tmp:
        tmp_path = tmp.name
    
    try:
        # Use the existing export_to_excel function to create the file
        export_to_excel(motifs, tmp_path, simple_format=simple_format)
        
        # Read the file as bytes
        with open(tmp_path, 'rb') as f:
            excel_bytes = f.read()
        
        return excel_bytes
    finally:
        # Clean up the temporary file
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


# ---------- ENHANCED PROFESSIONAL CSS FOR RESEARCH-QUALITY UI ----------
# Note: !important declarations are necessary to override Streamlit's default high-specificity CSS

# Initialize theme state in session state
if 'theme_mode' not in st.session_state:
    st.session_state.theme_mode = 'light'
if 'table_density' not in st.session_state:
    st.session_state.table_density = 'relaxed'
if 'color_theme' not in st.session_state:
    st.session_state.color_theme = 'scientific_blue'

def hex_to_rgb(hex_color: str) -> tuple:
    """Convert hex color to RGB tuple for CSS rgba() usage."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

def get_dna_pattern_svg(stroke_color: str) -> str:
    """Generate subtle DNA helix SVG pattern for background."""
    # Compact and URL-encoded SVG for background pattern
    svg = (
        "%3Csvg xmlns='http://www.w3.org/2000/svg' width='60' height='60' viewBox='0 0 60 60'%3E"
        f"%3Cg fill='none' stroke='%23{stroke_color}' stroke-width='0.6' opacity='0.12'%3E"
        "%3Cpath d='M10 30 C 18 12, 42 12, 50 30'/%3E"
        "%3Cpath d='M10 30 C 18 48, 42 48, 50 30'/%3E"
        "%3C/g%3E%3C/svg%3E"
    )
    return f"url(\"data:image/svg+xml,{svg}\")"

# Get current theme colors (using configuration from top)
current_theme = COLOR_THEMES.get(st.session_state.color_theme, COLOR_THEMES['scientific_blue'])
is_dark_mode = st.session_state.theme_mode == 'dark'
is_compact = st.session_state.table_density == 'compact'

# Dark mode color overrides - Soothing dark palette
if is_dark_mode:
    # Use the primary color from the base theme for consistency
    base_primary = COLOR_THEMES.get(st.session_state.color_theme, COLOR_THEMES['scientific_blue'])['primary']
    current_theme = {
        **current_theme,
        'bg_light': '#1A1F2E',        # Soft dark blue-gray
        'bg_card': '#252B3B',          # Slightly lighter dark
        'text': '#E5E7EB',             # Softer white for dark mode
        'tab_bg': '#1F2937',           # Dark slate for tab bar
        'tab_active': current_theme.get('primary', base_primary),
        'shadow': 'rgba(0, 0, 0, 0.25)'
    }

# Pre-calculate RGB values for all theme colors (performance optimization)
rgb = {key: hex_to_rgb(value) if key not in ['shadow'] else (0, 0, 0) for key, value in current_theme.items()}

# Additional color calculations for new soothing theme
tab_bg_color = current_theme.get('tab_bg', current_theme['bg_card'])
tab_active_color = current_theme.get('tab_active', current_theme['primary'])
shadow_color = current_theme.get('shadow', f"rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.15)")

# Generate SVG pattern based on theme
dna_pattern = get_dna_pattern_svg('1e3a5f' if is_dark_mode else 'bbdefb')

def load_css(theme_name=None):
    """
    Load external CSS file and inject dynamic theme variables.
    This makes the app.py file more succinct by separating styling concerns.
    If theme_name provided, use per-page theme, otherwise use session color_theme.
    """
    theme_to_use = COLOR_THEMES.get(theme_name, COLOR_THEMES.get(st.session_state.color_theme, COLOR_THEMES['scientific_blue']))
    is_dark = st.session_state.theme_mode == 'dark'
    
    # Read the external CSS file if present
    css_file_path = os.path.join(os.path.dirname(__file__), 'styles.css')
    try:
        with open(css_file_path, 'r') as f:
            css_content = f.read()
    except Exception:
        css_content = ""
    
    # Compute some derived values
    p_rgb = hex_to_rgb(theme_to_use['primary'])
    s_rgb = hex_to_rgb(theme_to_use['secondary'])
    tab_bg_color_local = theme_to_use.get('tab_bg', theme_to_use['bg_card'])
    tab_active_color_local = theme_to_use.get('tab_active', theme_to_use['primary'])
    dna_svg_local = get_dna_pattern_svg('1e3a5f' if is_dark else 'bbdefb')
    
    # Inject centralized configuration into CSS variables
    theme_vars = f"""
    <style>
    :root {{
        /* Theme Colors */
        --primary-color: {theme_to_use['primary']};
        --secondary-color: {theme_to_use['secondary']};
        --accent-color: {theme_to_use['accent']};
        --bg-light: {theme_to_use['bg_light']};
        --bg-card: {theme_to_use['bg_card']};
        --text-color: {theme_to_use['text']};
        --tab-bg: {tab_bg_color_local};
        --tab-active: {tab_active_color_local};
        --shadow-color: {theme_to_use['shadow']};
        --primary-rgb: {p_rgb[0]}, {p_rgb[1]}, {p_rgb[2]};
        --secondary-rgb: {s_rgb[0]}, {s_rgb[1]}, {s_rgb[2]};
        
        /* Typography from FONT_CONFIG */
        --font-primary: {FONT_CONFIG['primary_font']};
        --font-monospace: {FONT_CONFIG['monospace_font']};
        --font-h1: {FONT_CONFIG['h1_size']};
        --font-h2: {FONT_CONFIG['h2_size']};
        --font-h3: {FONT_CONFIG['h3_size']};
        --font-h4: {FONT_CONFIG['h4_size']};
        --font-body: {FONT_CONFIG['body_size']};
        --font-small: {FONT_CONFIG['small_size']};
        --font-caption: {FONT_CONFIG['caption_size']};
        --font-weight-light: {FONT_CONFIG['light_weight']};
        --font-weight-normal: {FONT_CONFIG['normal_weight']};
        --font-weight-medium: {FONT_CONFIG['medium_weight']};
        --font-weight-semibold: {FONT_CONFIG['semibold_weight']};
        --font-weight-bold: {FONT_CONFIG['bold_weight']};
        --font-weight-extrabold: {FONT_CONFIG['extrabold_weight']};
        
        /* Layout from LAYOUT_CONFIG */
        --border-radius-sm: {LAYOUT_CONFIG['border_radius']['small']};
        --border-radius-md: {LAYOUT_CONFIG['border_radius']['medium']};
        --border-radius-lg: {LAYOUT_CONFIG['border_radius']['large']};
        --border-radius-pill: {LAYOUT_CONFIG['border_radius']['pill']};
        --spacing-xs: 0.25rem;
        --spacing-sm: {LAYOUT_CONFIG['padding']['small']};
        --spacing-md: {LAYOUT_CONFIG['padding']['medium']};
        --spacing-lg: {LAYOUT_CONFIG['padding']['large']};
        --spacing-xl: {LAYOUT_CONFIG['padding']['xlarge']};
        --spacing-2xl: 2.5rem;
        
        /* Transitions from ANIMATION_CONFIG */
        --transition-fast: {ANIMATION_CONFIG['transition_fast']} {ANIMATION_CONFIG['easing_smooth']};
        --transition-normal: {ANIMATION_CONFIG['transition_normal']} {ANIMATION_CONFIG['easing_smooth']};
        --transition-slow: {ANIMATION_CONFIG['transition_slow']} {ANIMATION_CONFIG['easing_smooth']};
        
        /* Theme State */
        --dark-mode: {1 if is_dark else 0};
        --dna-pattern: {dna_svg_local};
    }}
    {css_content}
    </style>
    """
    
    st.markdown(theme_vars, unsafe_allow_html=True)

def get_page_colors(page_name='Home'):
    """
    Get color dictionary for inline HTML styles based on page context.
    
    NOTE: Inline HTML styles require literal color values and cannot use CSS variables.
    This function returns the current page's colors from the centralized token system
    for use in inline styles. All values are derived from the token block at the top.
    
    Args:
        page_name: Name of the page ('Home', 'Upload & Analyze', 'Results', etc.)
    
    Returns:
        Dictionary with page-specific colors from centralized tokens
    """
    # Map page names to their color palettes from centralized tokens
    page_color_map = {
        'Home': HOME_COLORS,
        'Upload & Analyze': INPUT_COLORS,
        'Analysis': ANALYSIS_COLORS,
        'Results': RESULTS_COLORS,
        'Visualization': VISUALIZATION_COLORS,
        'Download': DOWNLOAD_COLORS,
        'Documentation': DOCUMENTATION_COLORS,
    }
    
    # Get the page-specific palette
    page_palette = page_color_map.get(page_name, HOME_COLORS)
    
    # Return comprehensive color set combining page colors, global colors, and semantic colors
    return {
        # Page-specific colors
        **page_palette,
        # Global colors for consistent elements
        'white': GLOBAL_COLORS['white'],
        'neutral_50': GLOBAL_COLORS['neutral_50'],
        'neutral_100': GLOBAL_COLORS['neutral_100'],
        'neutral_200': GLOBAL_COLORS['neutral_200'],
        'neutral_500': GLOBAL_COLORS['neutral_500'],
        'neutral_600': GLOBAL_COLORS['neutral_600'],
        'neutral_700': GLOBAL_COLORS['neutral_700'],
        # Semantic colors for status indicators
        'success': SEMANTIC_COLORS['success'],
        'warning': SEMANTIC_COLORS['warning'],
        'error': SEMANTIC_COLORS['error'],
        'info': SEMANTIC_COLORS['info'],
    }

# Note: load_css() will be called per-page below to apply per-tab themes.

# ---------- CONSTANTS ----------
# Use classification config if available, otherwise fallback to defaults
CONFIG_AVAILABLE = False  # Configuration not available, use fallback
if CONFIG_AVAILABLE:
    MOTIF_ORDER = list(MOTIF_LENGTH_LIMITS.keys())
    # Expand special cases
    expanded_order = []
    for motif in MOTIF_ORDER:
        if motif == "Slipped_DNA_DR":
            expanded_order.extend(["Slipped DNA (Direct Repeat)", "Slipped DNA (STR)"])
        elif motif == "Slipped_DNA_STR":
            continue  # Already added above
        elif motif == "eGZ":
            expanded_order.append("eGZ (Extruded-G)")
        elif motif == "G4":
            expanded_order.extend(["Telomeric G4", "Stacked canonical G4s", "Stacked G4s with linker",
                                  "Canonical intramolecular G4", "Extended-loop canonical",
                                  "Higher-order G4 array/G4-wire", "Intramolecular G-triplex", "Two-tetrad weak PQS"])
        elif motif == "G-Triplex":
            expanded_order.append("G-Triplex")
        elif motif == "AC-motif":
            expanded_order.append("AC-Motif")
        elif motif == "A-philic_DNA":
            expanded_order.append("A-philic DNA")  # NEW: Class 9
        else:
            # Map to display names
            display_name = motif.replace("_", " ").replace("-", "-")
            if motif == "Curved_DNA":
                display_name = "Curved DNA"
            elif motif == "Z-DNA":
                display_name = "Z-DNA"
            elif motif == "R-Loop":
                display_name = "R-Loop"
            elif motif == "Triplex":
                display_name = "Triplex DNA"
            elif motif == "Sticky_DNA":
                display_name = "Sticky DNA"
            elif motif == "i-Motif":
                display_name = "i-Motif"
            expanded_order.append(display_name)
    
    MOTIF_ORDER = expanded_order + ["Hybrid", "Non-B DNA Clusters"]  # Classes 10, 11
else:
    # Fallback to original order with updated G4 subclasses
    MOTIF_ORDER = [
        "Sticky DNA","Curved DNA","Z-DNA","eGZ (Extruded-G)","Slipped DNA","R-Loop",
        "Cruciform","Triplex DNA","Telomeric G4","Stacked canonical G4s","Stacked G4s with linker",
        "Canonical intramolecular G4","Extended-loop canonical","Higher-order G4 array/G4-wire",
        "Intramolecular G-triplex","Two-tetrad weak PQS","i-Motif","AC-Motif","A-philic DNA",
        "Hybrid","Non-B DNA Clusters"
    ]

# Update color mapping to match the consolidated system and configuration
MOTIF_COLORS = {**MOTIF_CLASS_COLORS, **VISUALIZATION_CONFIG['motif_class_colors']}

PAGES = {
    "Home": "Overview",
    "Upload & Analyze": "Sequence Upload and Motif Analysis",
    "Results": "Analysis Results and Visualization",
    "Download": "Export Data",
    "Documentation": "Scientific Documentation & References"
}

# Remove Entrez setup since it's optional
if BIO_AVAILABLE:
    Entrez.email = "raazbiochem@gmail.com"
    Entrez.api_key = None

EXAMPLE_FASTA = """>Example Sequence
ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC
ATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT
GGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGTAGGAACCCCTAACCCCTAA
GAAAGAAGAAGAAGAAGAAGAAAGGAAGGAAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGG
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
GAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA
CTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT
"""
EXAMPLE_MULTI_FASTA = """>G4_iMotif_APhilic_Sequence
GGGATTGGGATTGGGATTGGGCCCATCCCTACCCTACCCAAACCCATCCCTACCCTACCCAATTTATTTAAAAA
AAAAAAAAAAAAAAAAAAAAAAGATCGAAAGATCGAAAGATCGAAAGATCGATGCGGCGGCGGCGGCGGCGGCGG
CGGCGGCGAATTCGAATTCGAATTCGAATTCCGCGCGCGCGCGCGCGCGCGAATGCATGCATGCATGCATGCAT
>Z_DNA_RLoop_Complex
CGCGCGCGCGCGCGCGCGCGCGATATATATATATATATATATCGCGCGCGCGCGCGCGCGCGGGGATGGGGATGG
GGATGGGGGGATGGGGATGGGGATGGGTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTGAAA
GAAAAAAGAAAGAAAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAG
>CurvedDNA_SlippedDNA_STR
AAAAAACGTTGCAAAAAACGTTGCAAAAAACGTTGCAAAAAATTTTTTCGAACGTTTTTTCGAACGTTTTTTCGA
ACGCAGCAGCAGCAGCAGCAGCAGCAGCTGCTGCTGCTGCTGCTGCTGCTGATCTGATCTGATCTGATCTGATC
TGATCTGATTCTATTCTATTCTATTCTATTCTATTCTATTCTGGCCCCGGCCCCGGCCCCGGCCCCTGCTGCTG
>Cruciform_Triplex_Mirror
ATGCCCGGGATCGGATCCGATCGAAATTCGATCGGATCCGATCCCGGGCATGAAAGAAAGAAAGAAAGAAAGAAA
GAAAGAAAGAAAAGATCCGGCCGATAGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGTTCCTCCTCCTCC
TCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCGAATTCCGAATTCCGAATTCCGAATTCCGAATTCGAAA
>Multi_iMotif_AC_Sequence  
AAATTTATTTAAATTTAAATTCCCTACCCTACCCTACCCAAAAATCCCTACCCTACCCTACCCGGAATCGATCG
ATCGATCGATCGATCGATCGCCCTACCCTACCCTACCCAAACCCTACCCTACCCTACCCAAAAAAAAAAAAAAAA
AAAAAAAAAAAGATCTAGATCTAGATCTAGATCTAGATCTAGATCTGAAAGAAAGAAAGAAAGAAAGAAAGAAA
"""

# Streamlined session state using consolidated system
for k, v in {
    'seqs': [],
    'names': [],
    'results': [],
    'summary_df': pd.DataFrame(),
    'analysis_status': "Ready",
    'selected_classes': [],  # Initialize empty list for motif class selection
    'current_job_id': None  # Current job ID for result delivery
}.items():
    if k not in st.session_state:
        st.session_state[k] = v

# ---- TABS ----
tabs = st.tabs(list(PAGES.keys()))
tab_pages = dict(zip(PAGES.keys(), tabs))

with tab_pages["Home"]:
    # Apply Home theme based on configuration
    load_css(TAB_THEMES.get('Home', 'scientific_blue'))
    
    # Get page-specific colors from centralized token system for inline HTML styles
    # NOTE: Inline HTML styles cannot use CSS variables, so we inject literal values
    # from the centralized color token system defined at the top of this file.
    # This ensures all colors are still managed centrally even in inline styles.
    colors = get_page_colors('Home')
    
    # ========== PROFESSIONAL HEADER ==========
    st.markdown(f"""
    <div style='background: linear-gradient(135deg, {colors['primary']} 0%, {colors['secondary']} 100%); 
                padding: 2.5rem 2rem; border-radius: 16px; margin-bottom: 2rem; 
                box-shadow: 0 8px 32px rgba(0,0,0,0.15); text-align: center;'>
        <h1 style='color: {colors['white']}; font-size: {FONT_CONFIG['h1_size']}; font-weight: {FONT_CONFIG['bold_weight']}; margin: 0 0 0.5rem 0; 
                   font-family: {FONT_CONFIG['primary_font']}; letter-spacing: -0.02em;'>
            {UI_TEXT['home_title']}
        </h1>
    </div>
    """, unsafe_allow_html=True)
    

    
    # ========== MAIN CONTENT GRID ==========
    left, right = st.columns([1, 1], gap="large")
    
    with left:
        st.markdown(f"""
        <div style='background: {colors['white']}; padding: 2rem; border-radius: 16px; 
                    box-shadow: 0 4px 20px rgba(0,0,0,0.08); border: 1px solid {colors['neutral_200']}; height: 100%;'>
            <h2 style='color: {colors['text']}; font-size: 1.6rem; margin: 0 0 1rem 0; font-weight: 600;'>
                Scientific Foundation
            </h2>
            <p style='color: {colors['neutral_700']}; font-size: 1rem; line-height: 1.8; margin-bottom: 1.2rem;'>
                <b style='color: {colors['primary']};'>Non-canonical DNA structures</b> are critical regulatory elements 
                implicated in genome stability, transcriptional regulation, replication, and disease mechanisms. 
                These structures deviate from the canonical B-form DNA helix and play essential roles in:
            </p>
            <ul style='color: {colors['neutral_700']}; font-size: 0.95rem; line-height: 1.7; padding-left: 1.5rem;'>
                <li><b>Genome Instability:</b> Hotspots for mutations and chromosomal rearrangements</li>
                <li><b>Gene Regulation:</b> Promoter and enhancer activity modulation</li>
                <li><b>DNA Replication:</b> Origins of replication and fork progression</li>
                <li><b>Disease Association:</b> Cancer, neurological disorders, and aging</li>
            </ul>
            
        </div>
        """, unsafe_allow_html=True)
        
        try:
            # Display NBD Circle logo
            possible_paths = ["nbdcircle.JPG", "archive/nbdcircle.JPG", "./nbdcircle.JPG"]
            image_found = False
            for img_path in possible_paths:
                if os.path.exists(img_path):
                    st.image(img_path, caption=UI_TEXT['home_image_caption'], use_container_width=True)
                    image_found = True
                    break
            if not image_found:
                raise FileNotFoundError("Image not found")  # Intentional: triggers fallback
        except Exception:
            # Placeholder if image not found - using page colors
            st.markdown(f"""
            <div style='background: linear-gradient(135deg, {colors['primary']} 0%, {colors['secondary']} 100%); 
                        border-radius: 15px; padding: 40px; text-align: center; color: {colors['white']}; margin-top: 1rem;'>
                <h2 style='margin: 0; color: {colors['white']}; font-size: 2rem;'>{UI_TEXT['home_image_fallback_title']}</h2>
                <h3 style='margin: 10px 0 0 0; color: {colors['white']};'>{UI_TEXT['home_image_fallback_subtitle']}</h3>
                <p style='margin: 5px 0 0 0; color: rgba(255,255,255,0.9);'>{UI_TEXT['home_image_fallback_caption']}</p>
            </div>
            """, unsafe_allow_html=True)
    
    with right:
        # NOTE: Motif Classes visualization uses specific color gradients per motif type
        # These match the VISUALIZATION_PALETTE defined in centralized tokens but must
        # be literal values here due to Streamlit's inline HTML constraints.
        # Colors are carefully coordinated with the centralized token system.
        st.markdown(f"""
        <div style='background: {colors['white']}; padding: 2rem; border-radius: 16px; 
                    box-shadow: 0 4px 20px rgba(0,0,0,0.08); border: 1px solid {colors['neutral_200']}; margin-bottom: 1.5rem;'>
            <h2 style='color: {colors['text']}; font-size: 1.6rem; margin: 0 0 1rem 0; font-weight: 600;'>
                Detected Motif Classes
            </h2>
            <div style='display: grid; grid-template-columns: 1fr 1fr; gap: 0.8rem; margin-top: 1rem;'>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, {SEMANTIC_COLORS['warning_light']} 0%, {SEMANTIC_COLORS['warning_border']} 100%); 
                            border-radius: 8px; border-left: 4px solid {SEMANTIC_COLORS['warning']};'>
                    <div style='font-weight: 600; color: {SEMANTIC_COLORS['warning_dark']}; font-size: 0.9rem;'>1. Curved DNA</div>
                    <div style='color: {SEMANTIC_COLORS['warning_dark']}; font-size: 0.75rem; margin-top: 0.2rem;'>A-tract curvature</div>
                </div>
                <!-- NOTE: Motif class visualization cards use semantic gradients for visual distinction.
                     Cards 1-2 use SEMANTIC_COLORS from centralized tokens.
                     Cards 3-11 use specific Tailwind-inspired gradients - this is intentional:
                     - Each motif type requires a unique, visually distinct gradient
                     - These form a carefully designed visual vocabulary for motif recognition
                     - Adding 9 × 4 gradient variations to centralized tokens would decrease maintainability
                     All colors follow semantic principles. This is an explicit design decision. -->
                <div style='padding: 0.8rem; background: linear-gradient(135deg, {SEMANTIC_COLORS['info_light']} 0%, {SEMANTIC_COLORS['info_border']} 100%); 
                            border-radius: 8px; border-left: 4px solid {SEMANTIC_COLORS['info']};'>
                    <div style='font-weight: 600; color: {SEMANTIC_COLORS['info_dark']}; font-size: 0.9rem;'>2. Slipped DNA</div>
                    <div style='color: {SEMANTIC_COLORS['info_dark']}; font-size: 0.75rem; margin-top: 0.2rem;'>Direct repeats, STRs</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #fce7f3 0%, #fbcfe8 100%); 
                            border-radius: 8px; border-left: 4px solid #ec4899;'>
                    <div style='font-weight: 600; color: #9f1239; font-size: 0.9rem;'>3. Cruciform</div>
                    <div style='color: #831843; font-size: 0.75rem; margin-top: 0.2rem;'>Palindromic IRs</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #d1fae5 0%, #a7f3d0 100%); 
                            border-radius: 8px; border-left: 4px solid #10b981;'>
                    <div style='font-weight: 600; color: #065f46; font-size: 0.9rem;'>4. R-Loop</div>
                    <div style='color: #064e3b; font-size: 0.75rem; margin-top: 0.2rem;'>RNA-DNA hybrids</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #fef9c3 0%, #fef08a 100%); 
                            border-radius: 8px; border-left: 4px solid #eab308;'>
                    <div style='font-weight: 600; color: #713f12; font-size: 0.9rem;'>5. Triplex</div>
                    <div style='color: #713f12; font-size: 0.75rem; margin-top: 0.2rem;'>Mirror repeats</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #ddd6fe 0%, #c4b5fd 100%); 
                            border-radius: 8px; border-left: 4px solid #8b5cf6;'>
                    <div style='font-weight: 600; color: #5b21b6; font-size: 0.9rem;'>6. G-Quadruplex</div>
                    <div style='color: #4c1d95; font-size: 0.75rem; margin-top: 0.2rem;'>7 subtypes</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #fed7aa 0%, #fdba74 100%); 
                            border-radius: 8px; border-left: 4px solid #f97316;'>
                    <div style='font-weight: 600; color: #7c2d12; font-size: 0.9rem;'>7. i-Motif</div>
                    <div style='color: #7c2d12; font-size: 0.75rem; margin-top: 0.2rem;'>C-rich structures</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #e0e7ff 0%, #c7d2fe 100%); 
                            border-radius: 8px; border-left: 4px solid #6366f1;'>
                    <div style='font-weight: 600; color: #3730a3; font-size: 0.9rem;'>8. Z-DNA</div>
                    <div style='color: #312e81; font-size: 0.75rem; margin-top: 0.2rem;'>Left-handed helix</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #ccfbf1 0%, #99f6e4 100%); 
                            border-radius: 8px; border-left: 4px solid #14b8a6;'>
                    <div style='font-weight: 600; color: #134e4a; font-size: 0.9rem;'>9. A-philic DNA</div>
                    <div style='color: #134e4a; font-size: 0.75rem; margin-top: 0.2rem;'>A/T-rich regions</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #f3e8ff 0%, #e9d5ff 100%); 
                            border-radius: 8px; border-left: 4px solid #a855f7;'>
                    <div style='font-weight: 600; color: #6b21a8; font-size: 0.9rem;'>10. Hybrid</div>
                    <div style='color: #581c87; font-size: 0.75rem; margin-top: 0.2rem;'>Multi-class overlap</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #e5e7eb 0%, #d1d5db 100%); 
                            border-radius: 8px; border-left: 4px solid #6b7280;'>
                    <div style='font-weight: 600; color: #1f2937; font-size: 0.9rem;'>11. Clusters</div>
                    <div style='color: #374151; font-size: 0.75rem; margin-top: 0.2rem;'>Motif hotspots</div>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        # Call to Action - Using page-specific colors
        st.markdown(f"""
        <div style='background: linear-gradient(135deg, {colors['primary']} 0%, {colors['secondary']} 100%); 
                    padding: 1.5rem; border-radius: 12px; text-align: center; 
                    box-shadow: 0 4px 12px {colors['shadow']};'>
            <h3 style='color: {colors['white']}; margin: 0 0 0.5rem 0; font-size: 1.2rem;'>
                {UI_TEXT['home_call_to_action_title']}
            </h3>
            <p style='color: rgba(255,255,255,0.95); margin: 0 0 1rem 0; font-size: 0.95rem;'>
                {UI_TEXT['home_call_to_action_text']}
            </p>
            <div style='background: {colors['white']}; color: {colors['primary']}; padding: 0.7rem 1.5rem; 
                        border-radius: 8px; display: inline-block; font-weight: 600; font-size: 1rem;'>
                {UI_TEXT['home_call_to_action_button']}
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    # ========== KEY FEATURES SECTION ==========
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown(f"""
    <div style='background: {colors['neutral_50']}; padding: 2rem; border-radius: 16px; margin-top: 2rem;'>
        <h2 style='color: {colors['text']}; font-size: 1.8rem; margin: 0 0 1.5rem 0; text-align: center; font-weight: 600;'>
            Key Features & Capabilities
        </h2>
        <div style='display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 1.5rem;'>
            <div style='background: {colors['white']}; padding: 1.5rem; border-radius: 12px; 
                        box-shadow: 0 2px 10px rgba(0,0,0,0.05); border: 1px solid {colors['neutral_200']};'>
                <h3 style='color: {colors['text']}; font-size: 1.1rem; margin: 0 0 0.5rem 0; font-weight: 600;'>
                    High Performance
                </h3>
                <p style='color: {colors['neutral_600']}; font-size: 0.9rem; line-height: 1.6; margin: 0;'>
                    24,674 bp/s processing speed. Handles sequences up to 1GB with chunked processing. 
                    O(n) complexity for all major detectors.
                </p>
            </div>
            <div style='background: {colors['white']}; padding: 1.5rem; border-radius: 12px; 
                        box-shadow: 0 2px 10px rgba(0,0,0,0.05); border: 1px solid {colors['neutral_200']};'>
                <h3 style='color: {colors['text']}; font-size: 1.1rem; margin: 0 0 0.5rem 0; font-weight: 600;'>
                    Publication Quality
                </h3>
                <p style='color: {colors['neutral_600']}; font-size: 0.9rem; line-height: 1.6; margin: 0;'>
                    25+ visualization types at 300 DPI resolution. Nature/NAR-compliant formats. 
                    Colorblind-friendly palettes (Wong 2011).
                </p>
            </div>
            <div style='background: {colors['white']}; padding: 1.5rem; border-radius: 12px; 
                        box-shadow: 0 2px 10px rgba(0,0,0,0.05); border: 1px solid {colors['neutral_200']};'>
                <h3 style='color: {colors['text']}; font-size: 1.1rem; margin: 0 0 0.5rem 0; font-weight: 600;'>
                    Scientifically Validated
                </h3>
                <p style='color: {colors['neutral_600']}; font-size: 0.9rem; line-height: 1.6; margin: 0;'>
                    Literature-based algorithms: QmRLFS, G4Hunter, Z-Seeker. 
                    Peer-reviewed methods with biological accuracy.
                </p>
            </div>
            <div style='background: {colors['white']}; padding: 1.5rem; border-radius: 12px; 
                        box-shadow: 0 2px 10px rgba(0,0,0,0.05); border: 1px solid {colors['neutral_200']};'>
                <h3 style='color: {colors['text']}; font-size: 1.1rem; margin: 0 0 0.5rem 0; font-weight: 600;'>
                    Statistical Analysis
                </h3>
                <p style='color: {colors['neutral_600']}; font-size: 0.9rem; line-height: 1.6; margin: 0;'>
                    Density analysis, enrichment calculations, p-value computation. 
                    100-iteration sequence shuffling for validation.
                </p>
            </div>
            <div style='background: {colors['white']}; padding: 1.5rem; border-radius: 12px; 
                        box-shadow: 0 2px 10px rgba(0,0,0,0.05); border: 1px solid {colors['neutral_200']};'>
                <h3 style='color: {colors['text']}; font-size: 1.1rem; margin: 0 0 0.5rem 0; font-weight: 600;'>
                    Multiple Export Formats
                </h3>
                <p style='color: {colors['neutral_600']}; font-size: 0.9rem; line-height: 1.6; margin: 0;'>
                    Excel (multi-sheet), CSV, BED, BigWig, JSON. 
                    UCSC/IGV genome browser compatible outputs.
                </p>
            </div>
            <div style='background: {colors['white']}; padding: 1.5rem; border-radius: 12px; 
                        box-shadow: 0 2px 10px rgba(0,0,0,0.05); border: 1px solid {colors['neutral_200']};'>
                <h3 style='color: {colors['text']}; font-size: 1.1rem; margin: 0 0 0.5rem 0; font-weight: 600;'>
                    Comprehensive Coverage
                </h3>
                <p style='color: {colors['neutral_600']}; font-size: 0.9rem; line-height: 1.6; margin: 0;'>
                    11 major classes, 22+ subclasses. Hybrid and cluster detection. 
                    Complete Non-B DNA structural characterization.
                </p>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    # ========== HOW TO CITE SECTION ==========
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown(f"""
    <div style='background: {colors['white']}; padding: 2rem; border-radius: 16px; 
                box-shadow: 0 4px 20px rgba(0,0,0,0.08); border: 1px solid {colors['neutral_200']}; margin-top: 2rem;'>
        <h2 style='color: {colors['primary']}; font-size: 1.6rem; margin: 0 0 1rem 0; font-weight: 600;'>
            How to Cite
        </h2>
        <div style='background: {colors['neutral_50']}; padding: 1.2rem; border-radius: 8px; border-left: 4px solid {colors['primary']}; 
                    font-family: "Courier New", monospace; font-size: 0.9rem; line-height: 1.7; color: {colors['neutral_700']};'>
            <b>NonBFinder: Comprehensive Detection and Analysis of Non-B DNA Motifs</b><br>
            Dr. Venkata Rajesh Yella<br>
            GitHub: <a href="https://github.com/VRYella/NonBFinder" style="color: {colors['primary']};">https://github.com/VRYella/NonBFinder</a><br>
            Email: yvrajesh_bt@kluniversity.in
        </div>
        <p style='color: {colors['neutral_600']}; font-size: 0.9rem; margin-top: 1rem; line-height: 1.6;'>
            If you use NonBFinder in your research, please cite this resource. 
            For methodology references, see the <b>Documentation</b> tab.
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # ========== JOB LOOKUP SECTION ==========
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown(f"""
    <div style='background: {colors['white']}; padding: 2rem; border-radius: 16px; 
                box-shadow: 0 4px 20px rgba(0,0,0,0.08); border: 1px solid {colors['neutral_200']}; margin-top: 2rem;'>
        <h2 style='color: {colors['primary']}; font-size: 1.6rem; margin: 0 0 1rem 0; font-weight: 600;'>
            🔍 Retrieve Previous Results
        </h2>
        <p style='color: {colors['neutral_700']}; font-size: 1rem; line-height: 1.6; margin-bottom: 1rem;'>
            Enter your Job ID to retrieve previously completed analyses. Results are stored securely and 
            accessible only via your unique Job ID.
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Job ID input
    col_input, col_button = st.columns([3, 1])
    with col_input:
        lookup_job_id = st.text_input(
            "Job ID",
            placeholder="Enter 10-character Job ID (e.g., a1b2c3d4e5)",
            help="Job ID was provided when you submitted your analysis",
            label_visibility="collapsed"
        )
    
    with col_button:
        lookup_button = st.button("Load Results", type="primary", use_container_width=True)
    
    # Handle job lookup
    if lookup_button and lookup_job_id:
        lookup_job_id = lookup_job_id.strip()
        
        if not lookup_job_id:
            st.error("Please enter a Job ID")
        elif not job_exists(lookup_job_id):
            st.error(f"Job ID '{lookup_job_id}' not found. Please check your ID and try again.")
        else:
            # Load the job
            with st.spinner(f"Loading job {lookup_job_id}..."):
                loaded_data = load_job_results(lookup_job_id)
                
                if loaded_data is None:
                    st.error("Failed to load job results. Please try again.")
                else:
                    results, sequences, names, metadata = loaded_data
                    
                    # Store in session state
                    st.session_state.results = results
                    st.session_state.seqs = sequences
                    st.session_state.names = names
                    st.session_state.current_job_id = lookup_job_id
                    st.session_state.analysis_done = True
                    
                    # Display success message
                    st.success(f"✅ Job {lookup_job_id} loaded successfully!")
                    
                    # Display job summary
                    st.markdown(f"""
                    <div style='background: {colors['neutral_50']}; padding: 1rem; border-radius: 8px; 
                                border-left: 4px solid {colors['primary']}; margin: 1rem 0;'>
                        <b>Job Summary:</b><br>
                        <b>Job ID:</b> {metadata.get('job_id', 'N/A')}<br>
                        <b>Completed:</b> {metadata.get('timestamp', 'N/A')}<br>
                        <b>Sequences:</b> {metadata.get('num_sequences', 0)}<br>
                        <b>Total Motifs:</b> {metadata.get('total_motifs', 0)}<br>
                        <b>Total Base Pairs:</b> {metadata.get('total_bp', 0):,}
                    </div>
                    """, unsafe_allow_html=True)
                    
                    st.info("📥 View results in the **Results** tab or download from the **Download** tab")

# Updated Streamlit layout: Input Method + Sequence Preview + Analysis side-by-side
#  "Upload & Analyze" 
# It expects helper functions and constants to exist in the module:
# parse_fasta, get_basic_stats, EXAMPLE_FASTA, EXAMPLE_MULTI_FASTA, parse_fasta, all_motifs_refactored,
# ensure_subclass, MOTIF_ORDER, wrap, Entrez, SeqIO, Counter, pd, st

with tab_pages["Upload & Analyze"]:
    # Apply Upload tab theme based on configuration
    load_css(TAB_THEMES.get('Upload & Analyze', 'nature_green'))
    st.markdown(f"<h2>{UI_TEXT['upload_title']}</h2>", unsafe_allow_html=True)

    # ----- TWO-COLUMN LAYOUT: Left for Upload, Right for Analysis -----
    col1, col2 = st.columns([1, 1])
    
    with col1:
        # LEFT COLUMN: Sequence Upload and Motif Analysis
        st.markdown(f"### {UI_TEXT['section_sequence_upload']}")
        
        # ----- Input Method -----
        input_method = st.radio(UI_TEXT['upload_input_method_prompt'],
                                [UI_TEXT['upload_method_file'], UI_TEXT['upload_method_paste'], 
                                 UI_TEXT['upload_method_example'], UI_TEXT['upload_method_ncbi']],
                                horizontal=True,
                                label_visibility="collapsed",
                                key="upload_method")

        seqs, names = [], []

        if input_method == UI_TEXT['upload_method_file']:
            fasta_file = st.file_uploader(UI_TEXT['upload_file_prompt'], 
                                         type=["fa", "fasta", "txt", "fna"],
                                         label_visibility="visible",
                                         help=UI_TEXT['upload_file_help'])
            if fasta_file:
                # Compact file card after upload
                file_size_mb = fasta_file.size / (1024 * 1024)
                
                # Memory-efficient processing with progress indicator
                with st.spinner(f"{UI_TEXT['upload_processing']} {fasta_file.name}..."):
                    # Get preview first (lightweight operation)
                    preview_info = get_file_preview(fasta_file, max_sequences=3)
                    
                    # Compact File Card
                    st.markdown(f"""
                    <div style='background: linear-gradient(135deg, #4A90E2 0%, #6AA5F2 100%); 
                                border-radius: 12px; padding: 12px; margin: 8px 0; color: white; 
                                box-shadow: 0 2px 8px rgba(74, 144, 226, 0.15);'>
                        <div style='display: flex; justify-content: space-between; align-items: center;'>
                            <div>
                                <div style='font-weight: 600; font-size: 0.95rem;'>File: {fasta_file.name}</div>
                                <div style='font-size: 0.85rem; opacity: 0.9; margin-top: 4px;'>
                                    {preview_info['num_sequences']} sequences | {preview_info['total_bp']:,} bp | {file_size_mb:.2f} MB
                                </div>
                            </div>
                            <div style='background: rgba(255,255,255,0.2); border-radius: 8px; padding: 8px 12px; font-weight: 600;'>
                                {UI_TEXT['label_valid']}
                            </div>
                        </div>
                    </div>
                    """, unsafe_allow_html=True)
                    
                    # Show preview of first few sequences using popover for better UX
                    with st.popover(UI_TEXT['upload_preview_button']):
                        for prev in preview_info['previews']:
                            st.markdown(f"**{prev['name']}**: {prev['length']:,} bp")
                            stats = get_basic_stats(prev['preview'].replace('...', ''))
                            st.caption(f"GC: {stats['GC%']}% | AT: {stats['AT%']}%")
                        
                        if preview_info['num_sequences'] > 3:
                            st.caption(f"...and {preview_info['num_sequences']-3} more sequences.")
                    
                    # Now parse all sequences using chunked parsing for memory efficiency
                    seqs, names = [], []
                    has_large_sequences = False
                    
                    if preview_info['num_sequences'] > 10:
                        # Show progress bar for files with many sequences
                        progress_bar = st.progress(0)
                        status_text = st.empty()
                        
                        for idx, (name, seq) in enumerate(parse_fasta_chunked(fasta_file)):
                            names.append(name)
                            seqs.append(seq)
                            
                            # Track if we have very large sequences
                            if len(seq) > 10_000_000:
                                has_large_sequences = True
                            
                            # Update progress
                            progress = (idx + 1) / preview_info['num_sequences']
                            progress_bar.progress(progress)
                            display_name = name[:50] + ('...' if len(name) > 50 else '')
                            status_text.text(f"Loading {idx + 1}/{preview_info['num_sequences']}: {display_name}")
                        
                        progress_bar.empty()
                        status_text.empty()
                    else:
                        # Fast path for small files
                        for name, seq in parse_fasta_chunked(fasta_file):
                            names.append(name)
                            seqs.append(seq)
                            
                            # Track if we have very large sequences
                            if len(seq) > 10_000_000:
                                has_large_sequences = True
                    
                    # Force garbage collection after loading all sequences if we had large ones
                    if has_large_sequences:
                        gc.collect()
                    
                    if not seqs:
                        st.warning(UI_TEXT['upload_no_sequences'])

        elif input_method == UI_TEXT['upload_method_paste']:
            seq_input = st.text_area(UI_TEXT['upload_paste_prompt'], 
                                    height=150, 
                                    placeholder=UI_TEXT['upload_paste_placeholder'],
                                    help=UI_TEXT['upload_paste_help'])
            if seq_input:
                seqs, names = [], []
                cur_seq, cur_name = "", ""
                for line in seq_input.splitlines():
                    if line.startswith(">"):
                        if cur_seq:
                            seqs.append(cur_seq)
                            names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                        cur_name = line.strip().lstrip(">")
                        cur_seq = ""
                    else:
                        cur_seq += line.strip()
                if cur_seq:
                    seqs.append(cur_seq)
                    names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                if seqs:
                    # Compact validation card
                    total_bp = sum(len(s) for s in seqs)
                    st.markdown(f"""
                    <div style='background: linear-gradient(135deg, #4A90E2 0%, #6AA5F2 100%); 
                                border-radius: 12px; padding: 12px; margin: 8px 0; color: white;
                                box-shadow: 0 2px 8px rgba(74, 144, 226, 0.15);'>
                        <div style='display: flex; justify-content: space-between; align-items: center;'>
                            <div>
                                <div style='font-weight: 600; font-size: 0.95rem;'>Pasted: Pasted Sequences</div>
                                <div style='font-size: 0.85rem; opacity: 0.9; margin-top: 4px;'>
                                    {len(seqs)} sequences | {total_bp:,} bp
                                </div>
                            </div>
                            <div style='background: rgba(255,255,255,0.2); border-radius: 8px; padding: 8px 12px; font-weight: 600;'>
                                Valid Valid
                            </div>
                        </div>
                    </div>
                    """, unsafe_allow_html=True)
                else:
                    st.warning(UI_TEXT['analysis_no_sequences_warning'])

        elif input_method == "Example Data":
            ex_type = st.radio("Example Type:", 
                             ["Single Example", "Multi-FASTA Example"], 
                             horizontal=True,
                             help="Load example sequences for testing")
            if ex_type == "Single Example":
                if st.button("Load Single Example", use_container_width=True):
                    parsed_fasta = parse_fasta(EXAMPLE_FASTA)
                    seqs = list(parsed_fasta.values())
                    names = list(parsed_fasta.keys())
                    st.success(UI_TEXT['upload_example_single_success'])
            else:
                if st.button("Load Multi-FASTA Example", use_container_width=True):
                    seqs, names = [], []
                    cur_seq, cur_name = "", ""
                    for line in EXAMPLE_MULTI_FASTA.splitlines():
                        if line.startswith(">"):
                            if cur_seq:
                                seqs.append(cur_seq)
                                names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                            cur_name = line.strip().lstrip(">")
                            cur_seq = ""
                        else:
                            cur_seq += line.strip()
                    if cur_seq:
                        seqs.append(cur_seq)
                        names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                    st.success(UI_TEXT['upload_example_multi_success'].format(count=len(seqs)))

        elif input_method == "NCBI Fetch":
            db = st.radio("NCBI Database", ["nucleotide", "gene"], horizontal=True,
                          help="Only nucleotide and gene databases are applicable for DNA motif analysis")
            query = st.text_input("Enter query (accession, gene, etc.):", 
                                help="e.g., NR_003287.2 or gene name")
            retmax = st.number_input("Max Records", min_value=1, max_value=20, value=3)
            if st.button("Fetch from NCBI", use_container_width=True):
                if query:
                    with st.spinner("Contacting NCBI..."):
                        try:
                            handle = Entrez.efetch(db=db, id=query, rettype="fasta", retmode="text")
                            records = list(SeqIO.parse(handle, "fasta"))
                            handle.close()
                            seqs = [str(rec.seq).upper().replace("U", "T") for rec in records]
                            names = [rec.id for rec in records]
                            if seqs:
                                st.success(UI_TEXT['upload_ncbi_success'].format(count=len(seqs)))
                        except Exception as e:
                            st.error(UI_TEXT['upload_ncbi_error'].format(error=e))
                else:
                    st.warning(UI_TEXT['upload_ncbi_empty_warning'])

        # Persist sequences to session state if any found from input
        if seqs:
            st.session_state.seqs = seqs
            st.session_state.names = names
            st.session_state.results = []

        # Compact sequence validation indicator using popover for cleaner UI
        if st.session_state.get('seqs'):
            with st.popover("Success: Validation Summary"):
                for i, seq in enumerate(st.session_state.seqs[:3]):
                    stats = get_basic_stats(seq)
                    st.markdown(f"**{st.session_state.names[i]}** ({len(seq):,} bp)")
                    st.caption(f"GC: {stats['GC%']}% | AT: {stats['AT%']}%")
                if len(st.session_state.seqs) > 3:
                    st.caption(f"...and {len(st.session_state.seqs)-3} more.")
    
    with col2:
        # RIGHT COLUMN: Analysis & Run
        st.markdown(f"### {UI_TEXT['heading_analysis_run']}")
        
        # Quick Options Section
        st.markdown(f"##### {UI_TEXT['analysis_quick_options_title']}")
        detailed_output = st.checkbox("Detailed Analysis", value=True, 
                                    help=UI_TEXT['tooltip_detailed_analysis'])
        quality_check = st.checkbox("Quality Validation", value=True, 
                                   help=UI_TEXT['tooltip_quality_validation'])
        
        # Advanced Options - now visible by default
        show_chunk_progress = st.checkbox("Show Chunk-Level Progress", value=False,
                                         help=UI_TEXT['tooltip_chunk_progress'])
        use_parallel_scanner = st.checkbox("Use Experimental Parallel Scanner", value=True,
                                          help=UI_TEXT['tooltip_parallel_scanner'])
        show_memory_usage = st.checkbox("Show Memory Usage", value=False,
                                        help="Display real-time memory usage during analysis")
        
        if use_parallel_scanner:
            st.caption(UI_TEXT['upload_parallel_note'])
        
        # Hardcoded default overlap handling
        nonoverlap = True
        overlap_option = "Remove overlaps within subclasses"
        
        # Helper text
        st.caption(UI_TEXT['upload_quick_options_note'])
    
    # ----- FULL-WIDTH STICKY RUN BUTTON -----
    st.markdown("---")
    
    # Initialize analysis_done flag if not present (idempotent run button)
    if "analysis_done" not in st.session_state:
        st.session_state.analysis_done = False
    
    # Check if valid input is present
    has_valid_input = bool(st.session_state.get('seqs'))
    
    # Create a full-width container for the run button
    run_button_container = st.container()
    with run_button_container:
        # Create two columns for Run and Reset buttons
        col_run, col_reset = st.columns([3, 1])
        
        with col_run:
            # Sticky-ish styling with disabled state
            if has_valid_input and not st.session_state.analysis_done:
                run_button = st.button(
                    UI_TEXT['analysis_run_button'],
                    type="primary",
                    use_container_width=True,
                    key="run_motif_analysis_main",
                    help="Start analyzing uploaded sequences for Non-B DNA motifs"
                )
            elif st.session_state.analysis_done:
                # Show analysis complete status
                st.markdown(f"""
                <div role="status" aria-label="Analysis Complete"
                     style='background: linear-gradient(135deg, #10b981 0%, #059669 100%); 
                            color: white; padding: 12px; 
                            border-radius: 12px; text-align: center; font-weight: 600;
                            font-size: 1.1rem;'>
                    ✅ Analysis Complete - View results in 'Results' tab
                </div>
                """, unsafe_allow_html=True)
                run_button = False
            else:
                # Disabled button appearance with accessibility
                st.markdown(f"""
                <div role="button" aria-disabled="true" aria-label="Run NBDScanner Analysis - Disabled: Please upload or paste a valid sequence first"
                     style='background: #e0e0e0; color: #9e9e9e; padding: 12px; 
                            border-radius: 12px; text-align: center; font-weight: 600;
                            font-size: 1.1rem; cursor: not-allowed; opacity: 0.6;'>
                    {UI_TEXT['analysis_run_button_disabled']}
                </div>
                <p style='text-align: center; color: #9e9e9e; font-size: 0.85rem; margin-top: 8px;' role="status">
                    {UI_TEXT['label_note']}: {UI_TEXT['analysis_run_button_disabled_note']}
                </p>
                """, unsafe_allow_html=True)
                run_button = False
        
        with col_reset:
            # Reset button to allow re-running analysis
            if st.button("🔄 Reset", use_container_width=True, help="Clear analysis results and reset for new run"):
                st.session_state.analysis_done = False
                st.session_state.results = []
                st.session_state.performance_metrics = None
                st.session_state.cached_visualizations = {}
                st.session_state.analysis_time = None
                st.rerun()
        
        # Placeholder for progress area
        progress_placeholder = st.empty()
    
    # ========== RUN ANALYSIS BUTTON LOGIC ========== 
    # Only run if button clicked AND not already done (idempotent)
    if run_button and not st.session_state.analysis_done:
        # Simplified validation
        if not st.session_state.seqs:
            st.error("Please upload or input sequences before running analysis.")
            st.session_state.analysis_status = "Error"
        else:
            # ============================================================
            # JOB ID GENERATION: Create unique ID immediately
            # ============================================================
            # Generate job ID at the start of analysis (or reuse if already exists)
            # This prevents overwriting existing job IDs on multiple button clicks
            if not st.session_state.get('current_job_id'):
                job_id = generate_job_id()
                st.session_state.current_job_id = job_id
            else:
                job_id = st.session_state.current_job_id
            
            # Display Job ID prominently to user
            st.markdown(f"""
            <div style='background: linear-gradient(135deg, #10b981 0%, #059669 100%); 
                        color: white; padding: 1.5rem; border-radius: 12px; margin: 1rem 0;
                        box-shadow: 0 4px 12px rgba(16, 185, 129, 0.3); text-align: center;'>
                <div style='font-size: 0.9rem; opacity: 0.9; margin-bottom: 0.5rem;'>
                    Your Job ID
                </div>
                <div style='font-size: 2rem; font-weight: 700; letter-spacing: 0.1em; 
                           font-family: "Courier New", monospace;'>
                    {job_id}
                </div>
                <div style='font-size: 0.85rem; opacity: 0.9; margin-top: 0.5rem;'>
                    Save this ID to retrieve your results later
                </div>
            </div>
            """, unsafe_allow_html=True)
            
            # Sequence length limit has been removed - the system now uses automatic chunking
            # (see NonBFinder.py CHUNK_THRESHOLD=10,000 bp) to handle sequences of any size
            st.session_state.analysis_status = "Running"
            
            # Store analysis parameters in session state for use in download section
            st.session_state.overlap_option_used = overlap_option
            st.session_state.nonoverlap_used = nonoverlap
            
            # Set analysis parameters based on user selections
            # nonoverlap is already set above based on user selection
            report_hotspots = True  # Enable hotspot detection 
            calculate_conservation = False  # Disable to reduce computation time
            threshold = 0.0  # Show all detected motifs (even 0 scores)
            
            validation_messages = []

            # Scientific validation check
            if CONFIG_AVAILABLE and st.session_state.get('selected_classes'):
                for class_id in st.session_state.selected_classes:
                    limits = get_motif_limits(class_id)
                    if limits:
                        validation_messages.append(f"Valid {class_id}: Length limits {limits}")
            
            # Enhanced progress tracking - timing captured exactly once at start and end
            import time
            
            # Create placeholder for progress
            progress_placeholder = st.empty()
            status_placeholder = st.empty()
            detailed_progress_placeholder = st.empty()
            timer_placeholder = st.empty()
            
            # ============================================================
            # DETERMINISTIC TIMING: Start time captured exactly once
            # ============================================================
            start_time = time.time()
            
            # Define detector processes for display
            DETECTOR_PROCESSES = [
                ("Curved DNA", "A-tract mediated DNA bending detection"),
                ("Slipped DNA", "Direct repeats and STR detection"),
                ("Cruciform", "Inverted repeat/palindrome detection"),
                ("R-Loop", "RNA-DNA hybrid formation site detection"),
                ("Triplex", "Three-stranded structure detection"),
                ("G-Quadruplex", "Four-stranded G-rich structure detection"),
                ("i-Motif", "C-rich structure detection"),
                ("Z-DNA", "Left-handed helix detection"),
                ("A-philic DNA", "A-rich structural element detection")
            ]
            
            # Constants for progress estimation
            # ESTIMATED_BP_PER_SECOND: Empirical processing rate based on benchmark testing
            # on 10kb sequences with all 9 detectors running. Actual speed may vary
            # depending on sequence complexity and hardware configuration.
            ESTIMATED_BP_PER_SECOND = 5800
            CHUNK_SIZE_FOR_PARALLEL = 50000  # Chunk size for parallel processing display
            
            # Helper function to display progress using Streamlit native components with scientific time format
            def display_progress_panel(container, elapsed, estimated_remaining, progress_display, 
                                     status_text, seq_name, seq_bp, seq_num, total_seqs, 
                                     processed_bp, total_bp, detector_count, extra_info=""):
                """Display progress panel using Streamlit native components with scientific time formatting.
                
                Uses canonical format: HH:MM:SS › mmm for precise, publication-quality time display.
                
                Args:
                    container: Streamlit container to display in
                    elapsed: Elapsed time in seconds
                    estimated_remaining: Estimated remaining time in seconds
                    progress_display: Progress percentage or status string
                    status_text: Status message
                    seq_name: Current sequence name
                    seq_bp: Current sequence length in bp
                    seq_num: Current sequence number
                    total_seqs: Total number of sequences
                    processed_bp: Total bp processed so far
                    total_bp: Total bp to process
                    detector_count: Number of detectors
                    extra_info: Extra information to display (e.g., speed, motifs)
                """
                with container:
                    st.subheader(UI_TEXT['analysis_progress_title'])
                    st.write(status_text)
                    
                    # Display metrics in 3-4 columns with scientific time format
                    if show_memory_usage:
                        col1, col2, col3, col4 = st.columns(4)
                    else:
                        col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        st.metric("Elapsed", format_time_scientific(elapsed))
                    
                    with col2:
                        st.metric("Remaining", format_time_scientific(estimated_remaining))
                    
                    with col3:
                        st.metric("Progress", progress_display)
                    
                    if show_memory_usage:
                        with col4:
                            mem_mb = get_memory_usage_mb()
                            st.metric("Memory", f"{mem_mb:.0f} MB")
                    
                    # Simplified sequence info display
                    st.write(f"**Sequence {seq_num}/{total_seqs}**: {seq_name} ({seq_bp:,} bp)")
                    st.write(f"Processed: {processed_bp:,} / {total_bp:,} bp")
                    
                    if extra_info:
                        st.write(extra_info)
            
            # Estimate processing time based on sequence length
            def estimate_time(total_bp):
                return total_bp / ESTIMATED_BP_PER_SECOND
            
            total_bp_all_sequences = sum(len(seq) for seq in st.session_state.seqs)
            estimated_total_time = estimate_time(total_bp_all_sequences)
            
            try:
                # Filter which classes to analyze based on selection
                analysis_classes = st.session_state.selected_classes if st.session_state.selected_classes else None
                
                # Run analysis on each sequence
                all_results = []
                all_hotspots = []
                
                total_bp_processed = 0
                
                with progress_placeholder.container():
                    pbar = st.progress(0)
                
                # Show detailed progress panel with detector sequence (only once since it's static)
                # The status shows all detectors as "running" during analysis since they run in parallel
                with detailed_progress_placeholder.container():
                    st.subheader(UI_TEXT['analysis_pipeline_title'])
                    
                    # Display detectors in a clean list format
                    for j, (detector_name, detector_desc) in enumerate(DETECTOR_PROCESSES):
                        st.write(f"**{j+1}. {detector_name}** - {detector_desc}")
                    
                    st.info(UI_TEXT['analysis_all_detectors_parallel'])
                    
                # ============================================================
                # MULTI-FASTA STABILITY: Per-sequence logic isolated
                # ============================================================
                # Each sequence processed independently with clean isolation.
                # No shared state between iterations except cumulative counters.
                # Results accumulated in list, stored atomically after loop.
                # Ensures identical behavior for single and multi-FASTA inputs.
                # ============================================================
                
                for i, (seq, name) in enumerate(zip(st.session_state.seqs, st.session_state.names)):
                    progress = (i + 1) / len(st.session_state.seqs)
                    
                    # ============================================================
                    # DETERMINISTIC EXECUTION: No timing inside loop
                    # Progress percentage only - elapsed time computed once at end
                    # ============================================================
                    
                    # Calculate overall percentage (deterministic, no timing)
                    overall_percentage = (total_bp_processed / total_bp_all_sequences * 100) if total_bp_all_sequences > 0 else 0
                    
                    # Determine status text based on progress state (no timing)
                    if total_bp_processed == 0:
                        status_text = "Starting analysis..."
                        progress_display = "Starting"
                    else:
                        status_text = "Analysis in progress..."
                        progress_display = f"{overall_percentage:.1f}%"
                    
                    # Build status message (no timing information)
                    status_msg = f"Processing sequence {i+1}/{len(st.session_state.seqs)}: {name} ({len(seq):,} bp)"
                    status_placeholder.info(status_msg)
                    
                    # Run the analysis - use parallel scanner for large sequences if enabled
                    # No per-sequence timing - total time captured once at end
                    
                    if use_parallel_scanner and len(seq) > 100000:
                        # NOTE: scanner_agent.py has been archived - parallel scanning experimental
                        # Use experimental parallel scanner for large sequences
                        try:
                            from scanner_agent import ParallelScanner
                            
                            # Create chunk progress placeholder
                            chunk_progress_placeholder = st.empty()
                            
                            # Track chunk progress
                            chunk_counter = {'current': 0, 'total': 0}
                            
                            def chunk_progress_callback(current, total):
                                """Callback to update chunk progress (ephemeral)"""
                                chunk_counter['current'] = current
                                chunk_counter['total'] = total
                                if show_chunk_progress:
                                    # Ephemeral progress (replaces previous)
                                    chunk_progress_placeholder.info(f"⚡ Parallel scanner processing chunks: {current}/{total} ({(current / total) * 100:.1f}%)")
                            
                            # Run parallel scanner with progress callback
                            # Use ephemeral status (replaces previous message)
                            status_placeholder.info(f"⚡ Using parallel scanner for {len(seq):,} bp sequence (est. chunks: ~{len(seq)//50000 + 1})")
                            
                            scanner = ParallelScanner(seq, hs_db=None)
                            
                            # The parallel scanner internally calls analyze_sequence on each chunk
                            # and returns full motif dictionaries with deduplication
                            results = scanner.run_scan(progress_callback=chunk_progress_callback)
                            
                            # Update sequence names for all motifs
                            for motif in results:
                                motif['Sequence_Name'] = name
                            
                            # Clear chunk progress and show ephemeral success (replaces previous message)
                            if show_chunk_progress:
                                chunk_progress_placeholder.success(f"✅ Parallel chunks complete: {len(results)} motifs from {chunk_counter['total']} chunks")
                            
                            # Ephemeral success message (replaces previous)
                            status_placeholder.success(f"✅ Parallel scanner completed: {len(results)} motifs detected")
                            
                        except Exception as e:
                            # Fallback to standard scanner on error (ephemeral warning)
                            status_placeholder.warning(f"⚠️ Parallel scanner failed, falling back to standard: {e}")
                            results = analyze_sequence(seq, name)
                    else:
                        # Use standard consolidated NBDScanner analysis
                        results = analyze_sequence(seq, name)
                    
                    # Ensure all motifs have required fields
                    results = [ensure_subclass(motif) for motif in results]
                    all_results.append(results)
                    
                    total_bp_processed += len(seq)
                    
                    # Memory management: Trigger garbage collection for large sequences
                    if len(seq) > 1_000_000:  # For sequences > 1 Mb
                        trigger_garbage_collection()
                        logger.debug(f"Triggered garbage collection after processing {name} ({len(seq):,} bp)")
                    
                    # ============================================================
                    # DETERMINISTIC TIMING: No intermediate time calculations
                    # Progress tracking only - elapsed time computed once at end
                    # ============================================================
                    
                    # Calculate actual progress percentage (bp-based, deterministic)
                    actual_percentage = (total_bp_processed / total_bp_all_sequences * 100) if total_bp_all_sequences > 0 else 0
                    
                    # Build progress info without timing (timing computed once at end)
                    progress_info = f"Progress: {total_bp_processed:,} / {total_bp_all_sequences:,} bp | Motifs: {len(results)} in this sequence"
                    
                    # Update progress display (no timing - pure progress tracking)
                    with progress_placeholder.container():
                        pbar.progress(progress, text=f"🧬 Analyzed {i+1}/{len(st.session_state.seqs)} sequences ({actual_percentage:.1f}%)")
                    
                    # Ephemeral success (replaces previous) - no per-sequence timing shown
                    status_placeholder.success(f"✅ {name}: {len(seq):,} bp | {len(results)} motifs detected")
                
                # ============================================================
                # MULTI-FASTA STABILITY: Results stored once atomically
                # ============================================================
                # All sequence results stored together in single atomic operation.
                # Prevents duplicate storage and ensures consistency across runs.
                # Identical behavior for single FASTA and multi-FASTA inputs.
                # ============================================================
                
                # Store results ONCE after all sequences processed
                st.session_state.results = all_results
                
                # ============================================================
                # DETERMINISTIC TIMING: End time captured exactly once
                # Runtime displayed only after all sequences complete
                # ============================================================
                total_time = time.time() - start_time
                overall_speed = total_bp_processed / total_time if total_time > 0 else 0
                
                # ============================================================
                # RIGOROUS VALIDATION & QUALITY CHECKS
                # ============================================================
                # Ephemeral info message (replaces previous)
                status_placeholder.info("🔍 Validating results for consistency and quality...")
                
                validation_issues = []
                
                # 1. Check for duplicate motifs within each sequence
                for i, results in enumerate(all_results):
                    seen_motifs = set()
                    duplicates_found = 0
                    for motif in results:
                        motif_key = (motif.get('Start'), motif.get('End'), motif.get('Class'), motif.get('Subclass'))
                        if motif_key in seen_motifs:
                            duplicates_found += 1
                        seen_motifs.add(motif_key)
                    
                    if duplicates_found > 0:
                        validation_issues.append(f"Note: Sequence {i+1}: {duplicates_found} duplicate motifs found")
                
                # 2. Validate motif data consistency
                for i, results in enumerate(all_results):
                    for motif in results:
                        # Check required fields
                        if not all(k in motif for k in ['Start', 'End', 'Class']):
                            validation_issues.append(f"Note: Sequence {i+1}: Motif missing required fields")
                            break
                        
                        # Validate positions
                        if motif.get('Start', 0) >= motif.get('End', 0):
                            validation_issues.append(f"Note: Sequence {i+1}: Invalid motif position (Start >= End)")
                            break
                        
                        # Validate length consistency
                        calculated_length = motif.get('End', 0) - motif.get('Start', 0)
                        if motif.get('Length') and abs(motif.get('Length') - calculated_length) > 1:
                            validation_issues.append(f"Note: Sequence {i+1}: Length mismatch detected")
                            break
                
                # 3. Check for overlapping motifs within same subclass (should be resolved)
                for i, results in enumerate(all_results):
                    subclass_motifs = {}
                    for motif in results:
                        subclass = motif.get('Subclass', 'Unknown')
                        if subclass not in subclass_motifs:
                            subclass_motifs[subclass] = []
                        subclass_motifs[subclass].append(motif)
                    
                    # Check for overlaps within each subclass
                    for subclass, motifs in subclass_motifs.items():
                        sorted_motifs = sorted(motifs, key=lambda m: m.get('Start', 0))
                        for j in range(len(sorted_motifs) - 1):
                            if sorted_motifs[j].get('End', 0) > sorted_motifs[j+1].get('Start', 0):
                                validation_issues.append(f"Note: Sequence {i+1}: Overlapping motifs in {subclass}")
                                break
                
                # Display validation results (ephemeral - replaces previous)
                if validation_issues:
                    warning_msg = f"⚠️ Validation found {len(validation_issues)} potential issues:\n"
                    for issue in validation_issues[:5]:  # Show first 5
                        warning_msg += f"\n• {issue}"
                    if len(validation_issues) > 5:
                        warning_msg += f"\n• ... and {len(validation_issues) - 5} more"
                    status_placeholder.warning(warning_msg)
                else:
                    status_placeholder.success("✅ Validation passed: No consistency issues found")
                
                # ============================================================
                # MULTI-FASTA STABILITY: Aggregate statistics computed once
                # ============================================================
                # Generate summary statistics ONCE after all sequences processed.
                # No per-sequence UI updates during this phase for determinism.
                # Results stored atomically in session state to prevent duplicates.
                # ============================================================
                
                # Generate summary
                summary = []
                for i, results in enumerate(all_results):
                    seq = st.session_state.seqs[i]
                    stats = get_basic_stats(seq, results)
                    summary.append({
                        'Sequence': st.session_state.names[i],
                        'Length': stats['Length'],
                        'GC Content': f"{stats['GC%']:.1f}%",
                        'Motifs Found': len(results),
                        'Unique Types': len(set(m.get('Type', 'Unknown') for m in results)),
                        'Avg Score': f"{np.mean([m.get('Score', 0) for m in results]):.3f}" if results else "0.000"
                    })
                
                # ATOMIC STORAGE: Store summary once in session state
                st.session_state.summary_df = pd.DataFrame(summary)
                
                # ============================================================
                # PRE-GENERATE ALL VISUALIZATIONS FOR CLASSES AND SUBCLASSES
                # ============================================================
                # Ephemeral info (replaces previous)
                status_placeholder.info("📊 Generating comprehensive visualizations for all classes and subclasses...")
                
                # Cache all visualizations for each sequence
                st.session_state.cached_visualizations = {}
                
                viz_start_time = time.time()
                total_viz_count = 0
                
                # Reduce UI updates by batching - only update every N sequences or at end
                UPDATE_INTERVAL = max(1, len(st.session_state.seqs) // 5)  # Update 5 times max
                
                for seq_idx, (seq, name, motifs) in enumerate(zip(st.session_state.seqs, st.session_state.names, all_results)):
                    sequence_length = len(seq)
                    
                    # Show all motifs including hybrid/cluster motifs
                    # No filtering is applied - all results are included in visualizations
                    filtered_motifs = motifs
                    
                    if not filtered_motifs:
                        continue
                    
                    viz_cache_key = f"seq_{seq_idx}"
                    st.session_state.cached_visualizations[viz_cache_key] = {}
                    
                    # Pre-calculate all density metrics (class and subclass level) - optimized batch calculation
                    try:
                        # Calculate all densities in one pass to avoid redundant iterations
                        genomic_density_class = calculate_genomic_density(filtered_motifs, sequence_length, by_class=True)
                        positional_density_class = calculate_positional_density(filtered_motifs, sequence_length, unit='kbp', by_class=True)
                        
                        genomic_density_subclass = calculate_genomic_density(filtered_motifs, sequence_length, 
                                                                            by_class=False, by_subclass=True)
                        positional_density_subclass = calculate_positional_density(filtered_motifs, sequence_length, 
                                                                                  unit='kbp', by_class=False, by_subclass=True)
                        
                        # Store density metrics
                        st.session_state.cached_visualizations[viz_cache_key]['densities'] = {
                            'class_genomic': genomic_density_class,
                            'class_positional': positional_density_class,
                            'subclass_genomic': genomic_density_subclass,
                            'subclass_positional': positional_density_subclass
                        }
                        
                        # Count unique classes and subclasses (cached for later use)
                        unique_classes = len(set(m.get('Class', 'Unknown') for m in filtered_motifs))
                        unique_subclasses = len(set(m.get('Subclass', 'Unknown') for m in filtered_motifs))
                        
                        st.session_state.cached_visualizations[viz_cache_key]['summary'] = {
                            'unique_classes': unique_classes,
                            'unique_subclasses': unique_subclasses,
                            'total_motifs': len(filtered_motifs)
                        }
                        
                        total_viz_count += 4  # Count density calculations
                        
                    except Exception as e:
                        # Log error but continue processing
                        pass
                
                viz_total_time = time.time() - viz_start_time
                
                # Memory management: Trigger garbage collection after all visualizations
                trigger_garbage_collection()
                logger.debug(f"Triggered garbage collection after generating {total_viz_count} visualizations")
                
                # Ephemeral success with scientific time format
                status_placeholder.success(f"✅ All visualizations prepared: {total_viz_count} components in {format_time_compact(viz_total_time)}")
                
                # Store performance metrics with enhanced details
                st.session_state.performance_metrics = {
                    'total_time': total_time,
                    'total_bp': total_bp_processed,
                    'speed': overall_speed,
                    'sequences': len(st.session_state.seqs),
                    'total_motifs': sum(len(r) for r in all_results),
                    'detector_count': len(DETECTOR_PROCESSES),  # Number of detector processes
                    'estimated_time': estimated_total_time,  # Initial estimated time
                    'visualization_time': viz_total_time,  # Time spent on visualizations
                    'visualization_count': total_viz_count,  # Number of visualization components
                    'validation_issues': len(validation_issues),  # Number of validation issues
                    # Derive analysis steps from DETECTOR_PROCESSES plus post-processing steps
                    'analysis_steps': [f"{name} detection" for name, _ in DETECTOR_PROCESSES] + [
                        'Hybrid/Cluster detection',
                        'Overlap resolution',
                        'Data validation',
                        'Class/Subclass visualization generation'
                    ]
                }
                
                # Clear progress displays
                progress_placeholder.empty()
                status_placeholder.empty()
                detailed_progress_placeholder.empty()
                
                # Show final success message with enhanced performance metrics using scientific time format
                timer_placeholder.markdown(f"""
                <div class='progress-panel progress-panel--success'>
                    <h3 class='progress-panel__title'>✅ Analysis Complete!</h3>
                    <p class='progress-panel__subtitle'>All detectors, validations, and visualizations completed successfully</p>
                    <div class='stats-grid stats-grid--wide'>
                        <div class='stat-card'>
                            <h2 class='stat-card__value'>{format_time_scientific(total_time)}</h2>
                            <p class='stat-card__label'>Analysis Time</p>
                        </div>
                        <div class='stat-card'>
                            <h2 class='stat-card__value'>{format_time_scientific(viz_total_time)}</h2>
                            <p class='stat-card__label'>Visualization Time</p>
                        </div>
                        <div class='stat-card'>
                            <h2 class='stat-card__value'>{total_bp_processed:,}</h2>
                            <p class='stat-card__label'>Base Pairs</p>
                        </div>
                        <div class='stat-card'>
                            <h2 class='stat-card__value'>{overall_speed:,.0f}</h2>
                            <p class='stat-card__label'>bp/second</p>
                        </div>
                        <div class='stat-card'>
                            <h2 class='stat-card__value'>{len(DETECTOR_PROCESSES)}</h2>
                            <p class='stat-card__label'>Detectors</p>
                        </div>
                        <div class='stat-card'>
                            <h2 class='stat-card__value'>{sum(len(r) for r in all_results)}</h2>
                            <p class='stat-card__label'>Motifs Found</p>
                        </div>
                        <div class='stat-card'>
                            <h2 class='stat-card__value'>{total_viz_count}</h2>
                            <p class='stat-card__label'>Viz Components</p>
                        </div>
                        <div class='stat-card'>
                            <h2 class='stat-card__value'>{len(validation_issues)}</h2>
                            <p class='stat-card__label'>Validation Issues</p>
                        </div>
                    </div>
                </div>
                """, unsafe_allow_html=True)
                
                # Show comprehensive completion summary with scientific time format
                # GOLD STANDARD: Final time display after everything is done
                completion_msg = f"""✅ **Analysis Complete!** All processing stages finished successfully:
                
**Detection & Analysis:**
- 🧬 {len(DETECTOR_PROCESSES)} detector processes completed
- 📊 {sum(len(r) for r in all_results)} total motifs detected across {len(st.session_state.seqs)} sequences
- ⏱️ Analysis completed in {format_time_scientific(total_time)} ({overall_speed:,.0f} bp/s)

**Quality Validation:**
- ✅ Data consistency checks: {'PASSED' if len(validation_issues) == 0 else f'{len(validation_issues)} issues found'}
- ✅ Non-redundancy validation: Complete
- ✅ Position validation: Complete

**Visualization Generation:**
- 📈 {total_viz_count} visualization components pre-generated
- 🎯 Class-level and subclass-level analysis ready
- ⏱️ Visualizations prepared in {format_time_scientific(viz_total_time)}

**View detailed results in the 'Results' tab.**
"""
                st.success(completion_msg)
                st.session_state.analysis_status = "Complete"
                
                # Set analysis_done flag for idempotent run button
                st.session_state.analysis_done = True
                st.session_state.analysis_time = total_time
                
                # ============================================================
                # RESULT PERSISTENCE: Save results to disk under job ID
                # ============================================================
                job_id = st.session_state.get('current_job_id')
                
                # Safety check: Regenerate job ID if somehow missing
                if not job_id:
                    logger.warning("Job ID missing from session state, regenerating")
                    job_id = generate_job_id()
                    st.session_state.current_job_id = job_id
                
                save_status_placeholder = st.empty()
                save_status_placeholder.info(f"💾 Saving results for Job ID: {job_id}...")
                
                # Prepare metadata
                job_metadata = {
                    'analysis_time': total_time,
                    'speed_bp_per_sec': overall_speed,
                    'detector_count': len(DETECTOR_PROCESSES),
                    'visualization_count': total_viz_count,
                    'validation_issues': len(validation_issues)
                }
                
                # Save results to disk
                save_success = save_job_results(
                    job_id,
                    all_results,
                    st.session_state.seqs,
                    st.session_state.names,
                    job_metadata
                )
                
                if save_success:
                    save_status_placeholder.success(f"✅ Results saved successfully! Job ID: **{job_id}**")
                else:
                    save_status_placeholder.warning(
                        "⚠️ Results could not be saved to disk, but are available in this session. "
                        "Download your results now from the Download tab."
                    )
                
            except Exception as e:
                progress_placeholder.empty()
                status_placeholder.empty()
                detailed_progress_placeholder.empty()
                timer_placeholder.empty()
                st.error(f"Analysis failed: {str(e)}")
                st.session_state.analysis_status = "Error"

    # End of Upload & Analyze tab
    st.markdown("---")
# ---------- RESULTS ----------
with tab_pages["Results"]:
    # Apply Results tab theme based on configuration
    load_css(TAB_THEMES.get('Results', 'genomic_purple'))
    st.markdown(f'<h2>{UI_TEXT["heading_analysis_results"]}</h2>', unsafe_allow_html=True)
    
    # Deterministic Results Page: Only render, never compute
    # If results are missing, show info and stop
    if not st.session_state.results:
        st.info(UI_TEXT['status_no_results'])
        st.info("💡 Run analysis first in the 'Upload & Analyze' tab")
        st.stop()  # Explicit stop to prevent any further execution
    
    # Performance metrics display if available
    if st.session_state.get('performance_metrics'):
        metrics = st.session_state.performance_metrics
        st.markdown(f"""
        <div class='progress-panel progress-panel--metrics'>
            <h3 class='progress-panel__title'>Performance Metrics</h3>
            <div class='stats-grid stats-grid--wide'>
                <div class='stat-card'>
                    <h2 class='stat-card__value'>{metrics['total_time']:.2f}s</h2>
                    <p class='stat-card__label'>Processing Time</p>
                </div>
                <div class='stat-card'>
                    <h2 class='stat-card__value'>{metrics['total_bp']:,}</h2>
                    <p class='stat-card__label'>Base Pairs</p>
                </div>
                <div class='stat-card'>
                    <h2 class='stat-card__value'>{metrics['speed']:,.0f}</h2>
                    <p class='stat-card__label'>bp/second</p>
                </div>
                <div class='stat-card'>
                    <h2 class='stat-card__value'>{metrics.get('detector_count', 9)}</h2>
                    <p class='stat-card__label'>Detector Processes</p>
                </div>
                <div class='stat-card'>
                    <h2 class='stat-card__value'>{metrics['sequences']}</h2>
                    <p class='stat-card__label'>Sequences</p>
                </div>
                <div class='stat-card'>
                    <h2 class='stat-card__value'>{metrics['total_motifs']}</h2>
                    <p class='stat-card__label'>Total Motifs</p>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    # Enhanced summary display
    st.markdown(f"### {UI_TEXT['heading_analysis_summary']}")
    st.dataframe(st.session_state.summary_df, use_container_width=True)
    
    # Sequence selection for detailed analysis using pills for better UX
    seq_idx = 0
    if len(st.session_state.seqs) > 1:
        # Use pills for sequence selection - a more modern and visual alternative to dropdown
        try:
            selected_seq = st.pills(
                "Choose Sequence for Details:",
                options=list(range(len(st.session_state.seqs))),
                format_func=lambda i: st.session_state.names[i],
                selection_mode="single",
                default=0,
                help="Select a sequence to view detailed analysis results"
            )
            seq_idx = selected_seq or 0
        except Exception:
            seq_idx = 0
    
    motifs = st.session_state.results[seq_idx]
    sequence_length = len(st.session_state.seqs[seq_idx])
    sequence_name = st.session_state.names[seq_idx]
    
    if not motifs:
        st.warning("No motifs detected for this sequence.")
    else:
        # Show all motifs including hybrid/cluster motifs
        # No filtering is applied - all results are displayed
        filtered_motifs = motifs
        hybrid_cluster_motifs = [m for m in motifs if m.get('Class') in ['Hybrid', 'Non-B_DNA_Clusters']]
        
        # Create enhanced motifs DataFrame
        df = pd.DataFrame(filtered_motifs) if filtered_motifs else pd.DataFrame()
        
        # Memory optimization: Optimize DataFrame for large result sets
        if len(df) > 1000:
            df = optimize_dataframe_memory(df)
            logger.debug(f"Optimized DataFrame memory for {len(df)} motifs")
        
        # Calculate and display enhanced coverage statistics (using filtered motifs)
        stats = get_basic_stats(st.session_state.seqs[seq_idx], filtered_motifs)
        
        motif_count = len(filtered_motifs)
        hybrid_cluster_count = len(hybrid_cluster_motifs)
        coverage_pct = stats.get("Motif Coverage %", 0)
        non_b_density = (motif_count / sequence_length * 1000) if sequence_length > 0 else 0
        
        # Enhanced summary card with modern research-quality styling
        st.markdown(f"""
        <div class='progress-panel progress-panel--results'>
            <h3 class='progress-panel__title progress-panel__title--large'>
                NBDScanner Analysis Results
            </h3>
            <div class='stats-grid stats-grid--extra-wide'>
                <div class='stat-card stat-card--large'>
                    <h2 class='stat-card__value stat_card_value_large'>
                        {stats.get("Coverage%", 0):.2f}%
                    </h2>
                    <p class='stat-card__label stat-card__label--large'>
                        Sequence Coverage
                    </p>
                </div>
                <div class='stat-card stat-card--large'>
                    <h2 class='stat-card__value stat-card__value--large'>
                        {stats.get("Density", 0):.2f}
                    </h2>
                    <p class='stat-card__label stat-card__label--large'>
                        Motif Density<br>(motifs/kb)
                    </p>
                </div>
                <div class='stat-card stat-card--large'>
                    <h2 class='stat-card__value stat-card__value--large'>
                        {motif_count}
                    </h2>
                    <p class='stat-card__label stat-card__label--large'>
                        Total Motifs
                    </p>
                </div>
                <div class='stat-card stat-card--large'>
                    <h2 class='stat-card__value stat-card__value--large'>
                        {sequence_length:,}
                    </h2>
                    <p class='stat-card__label stat-card__label--large'>
                        Sequence Length (bp)
                    </p>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        # Add info about hybrid/cluster motifs being shown separately
        if hybrid_cluster_count > 0:
            st.info(f"ℹ️ {hybrid_cluster_count} Hybrid/Cluster motifs detected. View them in the 'Cluster/Hybrid' tab below.")
        
        # Show cached visualization summary if available
        viz_cache_key = f"seq_{seq_idx}"
        cached_viz = st.session_state.get('cached_visualizations', {}).get(viz_cache_key, {})
        if cached_viz.get('summary'):
            viz_summary = cached_viz['summary']
            st.success(f"""✅ **Pre-generated Analysis Ready:** 
            {viz_summary['unique_classes']} unique classes, 
            {viz_summary['unique_subclasses']} unique subclasses analyzed
            """)
        
        # Enhanced motif table with canonical core reporting schema
        # Enhanced motif table with new columns and pagination for large datasets
        st.markdown("### 📋 All Detected Motifs")
        st.caption("**Core Reporting Schema** -- Publication-grade fields aligned with Nature/NAR/Genome Biology standards")
        
        # Info box explaining export-only advanced features
        st.info("""
        ℹ️ **Advanced Features:** All motif-specific details (ΔG components, dinucleotide counts, structural features, etc.) 
        are retained internally and available via exports (CSV/Excel/JSON). This view shows only publication-relevant fields 
        for clarity and interpretability.
        """)
        
        # Display only CORE_OUTPUT_COLUMNS (no column selection UI)
        # All advanced features remain available in exports
        display_columns = [col for col in CORE_OUTPUT_COLUMNS if col in df.columns]
        
        # Pagination for large datasets (improves performance)
        ROWS_PER_PAGE = 100
        total_rows = len(df)
        
        if total_rows > ROWS_PER_PAGE:
            # Add pagination controls
            col_pag1, col_pag2, col_pag3 = st.columns([2, 3, 2])
            with col_pag2:
                max_pages = (total_rows + ROWS_PER_PAGE - 1) // ROWS_PER_PAGE
                page_num = st.number_input(
                    f"Page (showing {ROWS_PER_PAGE} rows per page)",
                    min_value=1,
                    max_value=max_pages,
                    value=1,
                    step=1,
                    help=f"Total: {total_rows} motifs across {max_pages} pages"
                )
                start_idx = (page_num - 1) * ROWS_PER_PAGE
                end_idx = min(start_idx + ROWS_PER_PAGE, total_rows)
                
                st.caption(f"Showing motifs {start_idx + 1} to {end_idx} of {total_rows}")
            
            # Slice dataframe for current page
            df_page = df.iloc[start_idx:end_idx]
        else:
            df_page = df
        
        if display_columns:
            # Display core schema columns only
            filtered_df = df_page[display_columns].copy()
            # Replace underscores with spaces in column names for display
            filtered_df.columns = [col.replace('_', ' ') for col in filtered_df.columns]
            st.dataframe(filtered_df, use_container_width=True, height=360)
        else:
            # Fallback: display full dataframe if no core columns available
            display_df = df_page.copy()
            display_df.columns = [col.replace('_', ' ') for col in display_df.columns]
            st.dataframe(display_df, use_container_width=True, height=360)
        
        # NATURE-READY VISUALIZATION SUITE
        st.markdown(f'<h3>{UI_TEXT["heading_results_viz"]}</h3>', unsafe_allow_html=True)
        
        # Scientific Transparency Badge
        st.info(TRANSPARENCY_NOTE)
        
        # Create simplified visualization tabs based on Figure Panel layout
        viz_tabs = st.tabs([
            "📊 Figure 1: Global Landscape", 
            "🔗 Figure 2: Clustering & Co-occurrence",
            "📏 Figure 3: Structural Constraints (Optional)"
        ])
        
        # Check if clusters exist
        has_clusters = any(m.get('Class') == 'Non-B_DNA_Clusters' for m in filtered_motifs)
        
        # =================================================================
        # FIGURE 1: Global Non-B DNA Landscape
        # =================================================================
        with viz_tabs[0]:
            st.markdown("#### Figure 1: Global Non-B DNA Landscape")
            st.caption("*Purpose: What structures exist, and where?*")
            
            # Panel A: Motif Composition (nested donut + bar plots)
            st.markdown("##### Panel A: Motif Composition (Class → Subclass)")
            
            # Nested pie chart for hierarchical view
            try:
                fig_composition = plot_nested_pie_chart(
                    filtered_motifs, 
                    title=f"Motif Composition - {sequence_name}"
                )
                st.pyplot(fig_composition)
                plt.close(fig_composition)
            except Exception as e:
                st.error(f"Error generating composition plot: {e}")
            
            # Bar plots for class and subclass distributions
            col1, col2 = st.columns(2)
            with col1:
                try:
                    fig_class_bar = plot_motif_distribution(
                        filtered_motifs,
                        by='Class',
                        title=f"Motif Class Distribution - {sequence_name}"
                    )
                    st.pyplot(fig_class_bar)
                    plt.close(fig_class_bar)
                except Exception as e:
                    st.error(f"Error generating class bar plot: {e}")
            
            with col2:
                try:
                    fig_subclass_bar = plot_motif_distribution(
                        filtered_motifs,
                        by='Subclass',
                        title=f"Motif Subclass Distribution - {sequence_name}"
                    )
                    st.pyplot(fig_subclass_bar)
                    plt.close(fig_subclass_bar)
                except Exception as e:
                    st.error(f"Error generating subclass bar plot: {e}")
            
            # Panel B: Genome-Scale Localization (size-dependent: Manhattan OR Linear)
            st.markdown("##### Panel B: Genome-Scale Localization")
            try:
                if sequence_length > 50000:
                    # Large sequences: Manhattan plot
                    st.caption("*Manhattan plot for large genome (>50kb)*")
                    fig_position = plot_manhattan_motif_density(
                        filtered_motifs, sequence_length,
                        title=f"Motif Density Hotspots - {sequence_name}"
                    )
                    st.pyplot(fig_position)
                    plt.close(fig_position)
                else:
                    # Small sequences: Linear track
                    st.caption("*Linear track for short sequence (≤50kb)*")
                    fig_position = plot_linear_motif_track(
                        filtered_motifs, sequence_length,
                        title=f"Motif Track - {sequence_name}"
                    )
                    st.pyplot(fig_position)
                    plt.close(fig_position)
            except Exception as e:
                st.error(f"Error generating localization plot: {e}")
            
            # Panel C: Genome Coverage (compact bar chart)
            st.markdown("##### Panel C: Genome Coverage (% per motif class)")
            try:
                # Calculate or retrieve density metrics
                viz_cache_key = f"seq_{seq_idx}"
                cached_viz = st.session_state.get('cached_visualizations', {}).get(viz_cache_key, {})
                cached_densities = cached_viz.get('densities', {})
                
                if cached_densities:
                    genomic_density = cached_densities['class_genomic']
                    positional_density_kbp = cached_densities['class_positional']
                else:
                    genomic_density = calculate_genomic_density(filtered_motifs, sequence_length, by_class=True)
                    positional_density_kbp = calculate_positional_density(filtered_motifs, sequence_length, unit='kbp', by_class=True)
                
                # Compact visualization: density comparison
                fig_density = plot_density_comparison(
                    genomic_density, positional_density_kbp,
                    title="Motif Density Analysis"
                )
                st.pyplot(fig_density)
                plt.close(fig_density)
                
                # Compact density table
                density_data = []
                for class_name in sorted([k for k in genomic_density.keys() if k != 'Overall']):
                    density_data.append({
                        'Motif Class': class_name.replace('_', ' '),
                        'Coverage (%)': f"{genomic_density.get(class_name, 0):.3f}",
                        'Motifs/kb': f"{positional_density_kbp.get(class_name, 0):.2f}"
                    })
                
                if density_data:
                    density_df = pd.DataFrame(density_data)
                    st.dataframe(density_df, use_container_width=True, height=200)
                    
            except Exception as e:
                st.error(f"Error calculating density metrics: {e}")
        
        # =================================================================
        # FIGURE 2: Structural Clustering & Co-Occurrence
        # =================================================================
        with viz_tabs[1]:
            st.markdown("#### Figure 2: Structural Clustering & Co-Occurrence")
            st.caption("*Purpose: Which structures co-localize and form regulatory hotspots?*")
            
            # Panel D: Cluster Size Distribution (conditional on cluster existence)
            if has_clusters:
                st.markdown("##### Panel D: Cluster Size Distribution")
                try:
                    fig_cluster = plot_cluster_size_distribution(
                        filtered_motifs,
                        title=f"Cluster Statistics - {sequence_name}"
                    )
                    st.pyplot(fig_cluster)
                    plt.close(fig_cluster)
                except Exception as e:
                    st.error(f"Error generating cluster plot: {e}")
            else:
                st.info("**Panel D**: No clusters detected (requires multiple motifs in close proximity)")
            
            # Panel E: Motif Co-occurrence Matrix (always shown)
            st.markdown("##### Panel E: Motif Co-occurrence Matrix")
            st.caption("*Shows which motif classes tend to appear together (overlapping or within 1bp)*")
            try:
                fig_cooccur = plot_motif_cooccurrence_matrix(
                    filtered_motifs,
                    title=f"Co-occurrence Matrix - {sequence_name}"
                )
                st.pyplot(fig_cooccur)
                plt.close(fig_cooccur)
            except Exception as e:
                st.error(f"Error generating co-occurrence matrix: {e}")
        
        # =================================================================
        # FIGURE 3: Structural Constraints (Optional/Toggle)
        # =================================================================
        with viz_tabs[2]:
            st.markdown("#### Figure 3: Structural Constraints")
            st.caption("*Purpose: What are the physical constraints shaping each motif class?*")
            st.info("📌 **Optional Figure**: Toggle to include in main report or move to supplementary materials.")
            
            # User toggle for including this figure
            show_fig3 = st.checkbox(
                "Include Figure 3 in main report", 
                value=True,
                help="Enable to show structural constraint analysis in main figures"
            )
            
            if show_fig3:
                # Panel F: Length Distribution by Class (KDE only - no histogram redundancy)
                st.markdown("##### Panel F: Length Distributions by Class")
                try:
                    fig_length = plot_motif_length_kde(
                        filtered_motifs,
                        by_class=True,
                        title=f"Length Distribution (KDE) - {sequence_name}"
                    )
                    st.pyplot(fig_length)
                    plt.close(fig_length)
                except Exception as e:
                    st.error(f"Error generating length KDE: {e}")
                
                # Optional: Score distribution
                st.markdown("**Additional: Score Distribution by Class**")
                try:
                    fig_score = plot_score_distribution(
                        filtered_motifs, by_class=True,
                        title="Score Distribution (1-3 Scale)"
                    )
                    st.pyplot(fig_score)
                    plt.close(fig_score)
                except Exception as e:
                    st.error(f"Error generating score distribution: {e}")
            else:
                st.info("✓ Figure 3 hidden from main report (available in supplementary exports)")
                st.caption(SUPPLEMENTARY_NOTE)

# ---------- DOWNLOAD ----------
with tab_pages["Download"]:
    # Apply Download tab theme based on configuration
    load_css(TAB_THEMES.get('Download', 'clinical_teal'))
    
    if not st.session_state.results:
        st.info(UI_TEXT['download_no_results'])
    else:
        # Display current Job ID prominently
        current_job_id = st.session_state.get('current_job_id')
        if current_job_id:
            st.markdown(f"""
            <div style='background: linear-gradient(135deg, #0ea5e9 0%, #0284c7 100%); 
                        color: white; padding: 1rem; border-radius: 12px; margin-bottom: 1.5rem;
                        box-shadow: 0 4px 12px rgba(14, 165, 233, 0.3);'>
                <div style='display: flex; justify-content: space-between; align-items: center;'>
                    <div>
                        <div style='font-size: 0.85rem; opacity: 0.9;'>Job ID for this analysis:</div>
                        <div style='font-size: 1.5rem; font-weight: 700; letter-spacing: 0.1em; 
                                   font-family: "Courier New", monospace; margin-top: 0.3rem;'>
                            {current_job_id}
                        </div>
                    </div>
                    <div style='text-align: right;'>
                        <div style='font-size: 0.85rem; opacity: 0.9;'>
                            💾 Results saved<br/>
                            🔍 Retrieve anytime with this ID
                        </div>
                    </div>
                </div>
            </div>
            """, unsafe_allow_html=True)
        
        primary_sequence_name = st.session_state.names[0] if st.session_state.names else "Unknown Sequence"
        analysis_time = st.session_state.get('analysis_time', 0)
        
        # Sanitize sequence name for safe filenames
        safe_filename = re.sub(r'[^\w\-]', '_', primary_sequence_name)[:50].strip('_')
        

        
        # Prepare motif data
        all_motifs = []
        for i, motifs in enumerate(st.session_state.results):
            for m in motifs:
                export_motif = m.copy()
                if 'Sequence_Name' not in export_motif:
                    export_motif['Sequence_Name'] = st.session_state.names[i]
                all_motifs.append(export_motif)
        
        # Individual file downloads as main option
        st.markdown("### 📥 Download Results")
        st.markdown("""
        <div style='background: #f0f9ff; padding: 1rem; border-radius: 8px; margin-bottom: 1rem;
                    border-left: 4px solid #0ea5e9;'>
            <p style='color: #0c4a6e; margin: 0;'>
                💡 <strong>Download your results</strong> in different file formats
            </p>
        </div>
        """, unsafe_allow_html=True)
        
        col1, col2, col3, col4, col5 = st.columns(5)
        
        with col1:
            # CSV Export (All Motifs)
            if all_motifs:
                csv_data = export_to_csv(all_motifs, non_overlapping_only=False)
                st.download_button(
                    "📋 CSV (All Motifs)", 
                    data=csv_data.encode('utf-8'), 
                    file_name=f"{safe_filename}_all_motifs.csv", 
                    mime="text/csv",
                    use_container_width=True,
                    type="primary",
                    help="Download CSV with all detected motifs including Hybrid and Clusters"
                )
        
        with col2:
            # Excel Export (2-tab format)
            if all_motifs:
                try:
                    excel_bytes = generate_excel_bytes(all_motifs, simple_format=True)
                    st.download_button(
                        "📊 Excel (2 tabs)", 
                        data=excel_bytes, 
                        file_name=f"{safe_filename}_results.xlsx", 
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        use_container_width=True,
                        type="primary",
                        help="Download Excel with 2 tabs: NonOverlappingConsolidated, OverlappingAll"
                    )
                except Exception as e:
                    st.error(f"Excel export error: {str(e)}")
        
        with col3:
            # JSON Export  
            if all_motifs:
                json_data = export_to_json(all_motifs, pretty=True)
                st.download_button(
                    "📄 JSON", 
                    data=json_data.encode('utf-8'), 
                    file_name=f"{safe_filename}_results.json", 
                    mime="application/json",
                    use_container_width=True,
                    type="primary",
                    help="Download results in JSON format"
                )
        
        with col4:
            # BED Export
            if all_motifs and st.session_state.names:
                bed_data = export_to_bed(all_motifs, st.session_state.names[0])
                st.download_button(
                    "🧬 BED", 
                    data=bed_data.encode('utf-8'), 
                    file_name=f"{safe_filename}_results.bed", 
                    mime="text/plain",
                    use_container_width=True,
                    type="primary",
                    help="Download results in BED format for genome browsers"
                )
        
        with col5:
            # PDF Export for visual summaries
            if all_motifs and st.session_state.seqs:
                try:
                    # Get sequence length from the first sequence (with safe None check)
                    sequence_length = len(st.session_state.seqs[0]) if st.session_state.seqs and st.session_state.seqs[0] else 0
                    if sequence_length > 0:
                        pdf_data = export_to_pdf(all_motifs, sequence_length, primary_sequence_name)
                        st.download_button(
                            "📄 PDF (Visualizations)",
                            data=pdf_data,
                            file_name=f"{safe_filename}_visualizations.pdf",
                            mime="application/pdf",
                            use_container_width=True,
                            type="primary",
                            help="Download PDF containing graphical summarizations of sequence analyses"
                        )
                    else:
                        st.warning("No sequence data available for PDF generation")
                except Exception as e:
                    st.error(f"PDF export error: {str(e)}")
        
        # Add Distribution & Statistics Tables Download Section
        st.markdown("---")
        st.markdown("### 📊 Download Distribution & Statistics Tables")
        st.markdown("""
        <div style='background: #f0fdf4; padding: 1rem; border-radius: 8px; margin-bottom: 1rem;
                    border-left: 4px solid #22c55e;'>
            <p style='color: #14532d; margin: 0;'>
                📈 <strong>Download detailed distribution and density statistics</strong> for publication and analysis
            </p>
        </div>
        """, unsafe_allow_html=True)
        
        # Calculate distribution and density statistics for all sequences
        if all_motifs and st.session_state.seqs:
            try:
                # Prepare distribution statistics
                distribution_data = []
                for seq_idx, (seq, name, motifs) in enumerate(zip(st.session_state.seqs, st.session_state.names, st.session_state.results)):
                    sequence_length = len(seq)
                    
                    # Show all motifs including hybrid/cluster motifs
                    # No filtering is applied - all results are included in statistics
                    filtered_motifs = motifs
                    
                    # Calculate class-level statistics
                    class_counts = Counter(m.get('Class', 'Unknown') for m in filtered_motifs)
                    for class_name, count in class_counts.items():
                        genomic_density = (sum(m.get('Length', 0) for m in filtered_motifs if m.get('Class') == class_name) / sequence_length * 100) if sequence_length > 0 else 0
                        motifs_per_kbp = (count / sequence_length * 1000) if sequence_length > 0 else 0
                        avg_length = np.mean([m.get('Length', 0) for m in filtered_motifs if m.get('Class') == class_name])
                        
                        distribution_data.append({
                            'Sequence Name': name,
                            'Motif Class': class_name.replace('_', ' '),
                            'Count': count,
                            'Genomic Density (%)': f"{genomic_density:.4f}",
                            'Motifs per kbp': f"{motifs_per_kbp:.2f}",
                            'Average Length (bp)': f"{avg_length:.1f}",
                            'Total Coverage (bp)': sum(m.get('Length', 0) for m in filtered_motifs if m.get('Class') == class_name)
                        })
                
                distribution_df = pd.DataFrame(distribution_data)
                
                # Prepare subclass-level statistics
                subclass_data = []
                for seq_idx, (seq, name, motifs) in enumerate(zip(st.session_state.seqs, st.session_state.names, st.session_state.results)):
                    sequence_length = len(seq)
                    
                    # Show all motifs including hybrid/cluster motifs
                    # No filtering is applied - all results are included in statistics
                    filtered_motifs = motifs
                    
                    # Calculate subclass-level statistics
                    subclass_counts = Counter(m.get('Subclass', 'Unknown') for m in filtered_motifs)
                    for subclass_name, count in subclass_counts.items():
                        parent_class = next((m.get('Class') for m in filtered_motifs if m.get('Subclass') == subclass_name), 'Unknown')
                        genomic_density = (sum(m.get('Length', 0) for m in filtered_motifs if m.get('Subclass') == subclass_name) / sequence_length * 100) if sequence_length > 0 else 0
                        motifs_per_kbp = (count / sequence_length * 1000) if sequence_length > 0 else 0
                        avg_length = np.mean([m.get('Length', 0) for m in filtered_motifs if m.get('Subclass') == subclass_name])
                        
                        subclass_data.append({
                            'Sequence Name': name,
                            'Motif Class': parent_class.replace('_', ' '),
                            'Motif Subclass': subclass_name.replace('_', ' '),
                            'Count': count,
                            'Genomic Density (%)': f"{genomic_density:.4f}",
                            'Motifs per kbp': f"{motifs_per_kbp:.2f}",
                            'Average Length (bp)': f"{avg_length:.1f}",
                            'Total Coverage (bp)': sum(m.get('Length', 0) for m in filtered_motifs if m.get('Subclass') == subclass_name)
                        })
                
                subclass_df = pd.DataFrame(subclass_data)
                
                # Display preview of tables
                st.markdown("#### Class-Level Distribution Statistics")
                st.dataframe(distribution_df.head(10), use_container_width=True, height=300)
                st.caption(f"Showing first 10 of {len(distribution_df)} total records")
                
                st.markdown("#### Subclass-Level Distribution Statistics")
                st.dataframe(subclass_df.head(10), use_container_width=True, height=300)
                st.caption(f"Showing first 10 of {len(subclass_df)} total records")
                
                # Download buttons for statistics tables
                col_stat1, col_stat2, col_stat3 = st.columns(3)
                
                with col_stat1:
                    # Class-level CSV
                    class_csv = distribution_df.to_csv(index=False)
                    st.download_button(
                        "📊 Class Statistics (CSV)",
                        data=class_csv.encode('utf-8'),
                        file_name=f"{safe_filename}_class_statistics.csv",
                        mime="text/csv",
                        use_container_width=True,
                        help="Download class-level distribution statistics"
                    )
                
                with col_stat2:
                    # Subclass-level CSV
                    subclass_csv = subclass_df.to_csv(index=False)
                    st.download_button(
                        "📊 Subclass Statistics (CSV)",
                        data=subclass_csv.encode('utf-8'),
                        file_name=f"{safe_filename}_subclass_statistics.csv",
                        mime="text/csv",
                        use_container_width=True,
                        help="Download subclass-level distribution statistics"
                    )
                
                with col_stat3:
                    # Combined Excel with both sheets
                    try:
                        import io
                        from openpyxl import Workbook
                        
                        output = io.BytesIO()
                        with pd.ExcelWriter(output, engine='openpyxl') as writer:
                            distribution_df.to_excel(writer, sheet_name='Class Statistics', index=False)
                            subclass_df.to_excel(writer, sheet_name='Subclass Statistics', index=False)
                        
                        output.seek(0)
                        st.download_button(
                            "📊 All Statistics (Excel)",
                            data=output.getvalue(),
                            file_name=f"{safe_filename}_all_statistics.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                            use_container_width=True,
                            help="Download all statistics in Excel format with separate sheets"
                        )
                    except Exception as e:
                        st.error(f"Excel generation error: {str(e)}")
                
            except Exception as e:
                st.error(f"Error generating distribution statistics: {str(e)}")
                st.code(traceback.format_exc(), language="python")


# ---------- DOCUMENTATION ----------
with tab_pages["Documentation"]:
    # Apply Documentation tab theme based on configuration
    load_css(TAB_THEMES.get('Documentation', 'midnight'))
    st.header(UI_TEXT['heading_documentation'])
    
    # Motif classes documentation
    st.markdown("""
    <div style='background:#f4faff; border-radius:12px; padding:18px; font-size:1.08rem; font-family:Montserrat,Arial;'>
    <b>Motif Classes Detected:</b><br><br>
    <ul>
        <li><b>Curved DNA</b>: Identifies phased poly(A) or poly(T) tracts using regex and spacing rules, reflecting intrinsic curvature. Scoring is based on tract length/grouping.</li>
        <li><b>Z-DNA</b>: Detects alternating purine-pyrimidine patterns, GC-rich segments. Uses windowed scoring; regex finds dinucleotide repeats.</li>
        <li><b>eGZ-motif (Extruded-G Z-DNA)</b>: Searches for long (CGG)<sub>n</sub> runs via regex. Scored by repeat count.</li>
        <li><b>Slipped DNA</b>: Recognizes direct/tandem repeats (10–50 nt repeat unit, 0 nt spacer) by repeat-unit matching. Scoring by length and unit copies.</li>
        <li><b>R-Loop</b>: Finds G-rich regions for stable RNA-DNA hybrids; RLFS model and regex. Thermodynamic scoring for hybrid stability.</li>
        <li><b>Cruciform</b>: Finds palindromic inverted repeats (10–100 nt arms, 0–3 nt spacer) using reverse complement matching. Scoring by arm length and A/T content.</li>
        <li><b>Triplex DNA / Mirror Repeat</b>: Detects purine/pyrimidine mirror repeats (10–100 nt arms, 0–8 nt spacer, ≥90% purine or pyrimidine). Scoring by composition/purity.</li>
        <li><b>Sticky DNA</b>: Searches extended GAA/TTC repeats. Scoring by repeat count.</li>
        <li><b>G-Triplex</b>: Finds three consecutive guanine runs by regex and loop length. Scoring by G-run sum and loop penalty.</li>
        <li><b>G4 (G-Quadruplex) and Variants</b>: Detects canonical/variant G4 motifs by G-run/loop regex. G4Hunter scoring for content/structure.</li>
        <li><b>i-Motif</b>: C-rich sequences for i-motif under acid. Regex for C runs/loops; scoring by run count and content.</li>
        <li><b>AC-Motif</b>: Alternating A-rich/C-rich consensus regions by regex. Scoring by pattern presence.</li>
        <li><b>A-philic DNA</b>: Uses tetranucleotide log2 odds scoring to identify A-tract-favoring sequences with high protein-binding affinity. Classified as high-confidence or moderate A-philic ba[...]
        <li><b>Hybrid Motif</b>: Regions where motif classes overlap; found by interval intersection, scored on diversity/size.</li>
        <li><b>Non-B DNA Clusters</b>: Hotspots with multiple motifs in a window; sliding algorithm, scored by motif count/diversity.</li>
    </ul>
    <b>References:</b>
    <ul>
        <li>Bedrat et al., 2016 Nucleic Acids Research</li>
        <li>Ho et al., 2010 Nature Chemical Biology</li>
        <li>Kim et al., 2018 Nucleic Acids Research</li>
        <li>Zeraati et al., 2018 Nature Chemistry</li>
        <li>Bacolla et al., 2006 Nucleic Acids Research</li>
        <li>Mirkin & Frank-Kamenetskii, 1994 Annual Review of Biophysics</li>
        <li>Vinogradov, 2003 Bioinformatics (A-philic DNA tetranucleotide analysis)</li>
        <li>Bolshoy et al., 1991 PNAS (A-tract structural properties)</li>
        <li>Rohs et al., 2009 Nature (Protein-DNA interactions, A-philic binding)</li>
        <li>New et al., 2020 Journal of DNA Structure</li>
    </ul>
    </div>
    """, unsafe_allow_html=True)
    
    # Add configuration information if available
    if CONFIG_AVAILABLE:
        st.markdown("""
        <div style='background:#f1f5f9; border-radius:12px; padding:18px; font-size:1.08rem; font-family:Montserrat,Arial; margin-top:20px;'>
        <b>Scoring Configuration Details</b>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown("### Motif Length Constraints")
        
        config_df = pd.DataFrame([
            {
                "Motif Class": motif_class,
                "Min Length (bp)": limits["S_min"],
                "Max Length (bp)": limits["S_max"],
                "Biological Basis": f"Based on {SCORING_METHODS.get(motif_class, {}).get('reference', 'literature survey')}"
            }
            for motif_class, limits in MOTIF_LENGTH_LIMITS.items()
        ])
        
        st.dataframe(config_df, use_container_width=True)
        
        st.markdown("### Scoring Methods")
        scoring_df = pd.DataFrame([
            {
                "Motif Class": motif_class,
                "Method": method_info.get("method", ""),
                "Description": method_info.get("description", ""),
                "Reference": method_info.get("reference", "")
            }
            for motif_class, method_info in SCORING_METHODS.items()
        ])
        
        st.dataframe(scoring_df, use_container_width=True)

st.markdown(f"""
---
<div style='font-size: 1.05rem; color: #1e293b; margin-top: 30px; text-align: left; font-family:{FONT_CONFIG['primary_font']};'>
<b>Developed by</b><br>
{UI_TEXT['author']}<br>
<a href='mailto:{UI_TEXT['author_email']}'>{UI_TEXT['author_email']}</a> |
<a href='https://github.com/VRYella' target='_blank'>GitHub: VRYella</a>
</div>
""", unsafe_allow_html=True)
