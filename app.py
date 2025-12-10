"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                         NBDSCANNER WEB APPLICATION                            ║
║                    Non-B DNA Motif Detection System                          ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: app.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.1
LICENSE: MIT

DESCRIPTION:
    Streamlit web application for comprehensive Non-B DNA motif detection.
    Provides interactive interface for sequence analysis and visualization.

FEATURES:
┌─────────────────────────────────────────────────────────────────────────────┐
│  - Multi-FASTA support             - Real-time analysis progress            │
│  - 11 motif classes detection      - Interactive visualizations             │
│  - 22+ subclass analysis           - Export to CSV/BED/JSON                 │
│  - NCBI sequence fetch             - Publication-quality plots              │
└─────────────────────────────────────────────────────────────────────────────┘

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
# Import consolidated NBDScanner modules
from utilities import (
    parse_fasta, parse_fasta_chunked, get_file_preview, wrap, get_basic_stats, export_to_bed, export_to_csv,
    export_to_json, export_to_excel, calculate_genomic_density, calculate_positional_density
)
from nonbscanner import (
    analyze_sequence, get_motif_info as get_motif_classification_info
)
from utilities import export_results_to_dataframe
from visualizations import (
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

# Try to import Entrez for demo functionality
try:
    from Bio import Entrez, SeqIO
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False

# Try to import Hyperscan (optional for acceleration)
try:
    import hyperscan
    HYPERSCAN_AVAILABLE = True
except ImportError:
    HYPERSCAN_AVAILABLE = False

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
    Cache compiled Hyperscan database for reuse across sessions.
    
    Note: This is a placeholder for when Hyperscan patterns are defined.
    Currently returns None as Hyperscan is not yet integrated.
    
    Note: The underscore prefix on _patterns is a Streamlit convention
    to indicate the parameter should not be hashed for caching purposes.
    This prevents Streamlit from attempting to hash complex pattern objects.
    
    Args:
        _patterns: List of regex patterns to compile (underscore prefix prevents hashing)
        
    Returns:
        Compiled Hyperscan database or None
    """
    if not HYPERSCAN_AVAILABLE or _patterns is None:
        return None
    
    try:
        # Compile patterns into Hyperscan database
        # This would be implemented when Hyperscan patterns are defined
        db = None  # Placeholder
        return db
    except Exception as e:
        st.warning(f"Hyperscan database compilation failed: {e}")
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
    page_title="NBDScanner - Non-B DNA Motif Finder",
    layout="wide",
    page_icon=None,
    menu_items={'About': "NBDScanner | Developed by Dr. Venkata Rajesh Yella"}
)

# =============================================================================
# SERVER/STREAMLIT VERSION: SEQUENCE LENGTH LIMIT REMOVED
# =============================================================================
# The previous maximum sequence length limit (1 Mbp) has been REMOVED.
# The system now uses chunking with 10,000 nucleotide chunks and 500 bp overlap
# to handle sequences of any size efficiently.
#
# Chunking is automatically applied to large sequences for optimal performance.
# See nonbscanner.py for chunking configuration:
# - CHUNK_THRESHOLD: 10,000 bp (sequences larger than this are chunked)
# - DEFAULT_CHUNK_SIZE: 10,000 bp per chunk
# - DEFAULT_CHUNK_OVERLAP: 500 bp overlap between chunks


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
def generate_excel_bytes(motifs):
    """
    Generate Excel file as bytes for Streamlit download button.
    
    Args:
        motifs: List of motif dictionaries
        
    Returns:
        bytes: Excel file data as bytes
    """
    import tempfile
    
    # Create a temporary file
    with tempfile.NamedTemporaryFile(mode='wb', suffix='.xlsx', delete=False) as tmp:
        tmp_path = tmp.name
    
    try:
        # Use the existing export_to_excel function to create the file
        export_to_excel(motifs, tmp_path)
        
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
    return f"url(\"data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='60' height='60' viewBox='0 0 60 60'%3E%3Cg fill='none' stroke='%23{stroke_color}' stroke-width='0.5' opacity='0.3'%3E%3Cpath d='M10 10 Q 30 0, 50 10 T 90 10'/%3E%3Ccircle cx='10' cy='10' r='2'/%3E%3Ccircle cx='50' cy='10' r='2'/%3E%3C/g%3E%3C/svg%3E\")"

# Color theme palettes for scientific styling - Soothing UI Theme
# Using CSS variables approach with soft, modern colors optimized for readability
COLOR_THEMES = {
    'scientific_blue': {
         "primary":    "#4A90E2",   #// Clean scientific blue (balanced, not too saturated)
         "secondary":  "#6AA5F2",   #// Lighter complementary blue
         "accent":     "#AFCBFF",   #// Subtle highlight blue for emphasis
         "bg_light":   "#F7FAFC",   #// Very soft cool white (best for reading)
         "bg_card":    "#E8F0FF",   #// Gentle blue-tinted card/panel background
         "text":       "#1E252F",   #// Deep neutral blue-gray (high readability)
         "tab_bg":     "#EDF2F7",   #// Cool scientific neutral
         "tab_active": "#4A90E2",   #// Matches primary for consistency
         "shadow":     "rgba(74, 144, 226, 0.15)" #// Soft blue shadow for depth
    },
    'nature_green': {
        'primary': '#6DBB7A',        # Soft green
        'secondary': '#8FD19C',      # Lighter green
        'accent': '#B8E5C0',         # Very soft green accent
        'bg_light': '#F7FBF8',       # Very light green background
        'bg_card': '#E8F5E9',        # Soft green card background
        'text': '#1F2937',           # Darker gray for better readability
        'tab_bg': '#F0F9F1',
        'tab_active': '#6DBB7A',
        'shadow': 'rgba(109, 187, 122, 0.15)'
    },
    'genomic_purple': {
        'primary': '#9B8FD9',        # Soft purple
        'secondary': '#B3A9E5',      # Lighter purple
        'accent': '#D4CFF0',         # Very soft purple accent
        'bg_light': '#FAF9FD',       # Very light purple background
        'bg_card': '#EDE7F6',        # Soft purple card background
        'text': '#1F2937',           # Darker gray for better readability
        'tab_bg': '#F5F3FA',
        'tab_active': '#9B8FD9',
        'shadow': 'rgba(155, 143, 217, 0.15)'
    },
    'clinical_teal': {
        'primary': '#5BBFB4',        # Soft teal
        'secondary': '#7CCFC6',      # Lighter teal
        'accent': '#A8E4DE',         # Very soft teal accent
        'bg_light': '#F6FBFB',       # Very light teal background
        'bg_card': '#E0F2F1',        # Soft teal card background
        'text': '#1F2937',           # Darker gray for better readability
        'tab_bg': '#EBF7F6',
        'tab_active': '#5BBFB4',
        'shadow': 'rgba(91, 191, 180, 0.15)'
    }
}

# Get current theme colors
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

def load_css():
    """
    Load external CSS file and inject dynamic theme variables.
    This makes the app.py file more succinct by separating styling concerns.
    """
    # Read the external CSS file
    css_file_path = os.path.join(os.path.dirname(__file__), 'styles.css')
    with open(css_file_path, 'r') as f:
        css_content = f.read()
    
    # Inject dynamic theme variables
    theme_vars = f"""
    <style>
    :root {{
        --primary-color: {current_theme['primary']};
        --secondary-color: {current_theme['secondary']};
        --accent-color: {current_theme['accent']};
        --bg-light: {current_theme['bg_light']};
        --bg-card: {current_theme['bg_card']};
        --text-color: {current_theme['text']};
        --tab-bg: {tab_bg_color};
        --tab-active: {tab_active_color};
        --shadow-color: {shadow_color};
        --primary-rgb: {rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]};
        --secondary-rgb: {rgb['secondary'][0]}, {rgb['secondary'][1]}, {rgb['secondary'][2]};
        --dark-mode: {1 if is_dark_mode else 0};
        --is-compact: {1 if is_compact else 0};
        --dna-pattern: {dna_pattern};
    }}
    {css_content}
    </style>
    """
    
    st.markdown(theme_vars, unsafe_allow_html=True)

# Load the CSS with dynamic theme variables
load_css()




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
            expanded_order.extend(["G4", "Relaxed G4", "Bulged G4", "Bipartite G4", "Multimeric G4"])
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
    # Fallback to original order with A-philic DNA added
    MOTIF_ORDER = [
        "Sticky DNA","Curved DNA","Z-DNA","eGZ (Extruded-G)","Slipped DNA","R-Loop",
        "Cruciform","Triplex DNA","G-Triplex","G4","Relaxed G4","Bulged G4","Bipartite G4",
        "Multimeric G4","i-Motif","AC-Motif","A-philic DNA","Hybrid","Non-B DNA Clusters"
    ]

# Update color mapping to match the consolidated system
MOTIF_COLORS = MOTIF_CLASS_COLORS  # Use the consolidated colors

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
    'selected_classes': []  # Initialize empty list for motif class selection
}.items():
    if k not in st.session_state:
        st.session_state[k] = v

# ---- TABS ----
tabs = st.tabs(list(PAGES.keys()))
tab_pages = dict(zip(PAGES.keys(), tabs))

with tab_pages["Home"]:
    st.markdown("<h1>Non-B DNA Motif Finder</h1>", unsafe_allow_html=True)
    left, right = st.columns([1,1])
    with left:
        try:
            # Display NBD Circle logo - check multiple possible locations
            possible_paths = ["nbdcircle.JPG", "archive/nbdcircle.JPG", "./nbdcircle.JPG"]
            image_found = False
            for img_path in possible_paths:
                if os.path.exists(img_path):
                    st.image(img_path, use_container_width=True)
                    image_found = True
                    break
            if not image_found:
                raise FileNotFoundError("Image not found")
        except Exception:
            # If image doesn't exist or cannot be loaded, show placeholder
            st.markdown("""
            <div style='background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                        border-radius: 15px; padding: 40px; text-align: center; color: white;'>
                <h2 style='margin: 0; color: white;'>DNA</h2>
                <h3 style='margin: 10px 0 0 0; color: white;'>NBD Finder</h3>
                <p style='margin: 5px 0 0 0; color: #E8E8E8;'>Non-B DNA Detection</p>
            </div>
            """, unsafe_allow_html=True)
    with right:
        st.markdown("""
        <div style='font-family: Inter, system-ui, sans-serif; font-size:0.95rem; color:#263238; 
                    line-height:1.5; padding:1.2rem; background:white; border-radius:12px; 
                    box-shadow:0 2px 8px rgba(0,0,0,0.06); border: 1px solid rgba(25, 118, 210, 0.1);'>
        <p style='margin-top:0; margin-bottom:0.8rem;'><b style='color:#0d47a1;'>Non-canonical DNA structures</b> play key roles in genome stability, regulation, and evolution.</p>
        <p style='margin-bottom:0.8rem;'>Detects and analyzes <b style='color:#1976d2;'>11 major classes with 22+ subclasses</b> of Non-B DNA motifs.</p>
        <p style='margin-bottom:0.5rem;'><b style='color:#0d47a1;'>Motif Classes:</b></p>
        <div style='color:#37474f; font-size:0.88rem; line-height:1.6; padding-left:0.8rem; column-count:2; column-gap:1rem;'>
            <div style='margin-bottom:0.2rem;'><b>1.</b> Curved DNA</div>
            <div style='margin-bottom:0.2rem;'><b>2.</b> Slipped DNA</div>
            <div style='margin-bottom:0.2rem;'><b>3.</b> Cruciform DNA</div>
            <div style='margin-bottom:0.2rem;'><b>4.</b> R-loop</div>
            <div style='margin-bottom:0.2rem;'><b>5.</b> Triplex</div>
            <div style='margin-bottom:0.2rem;'><b>6.</b> G-Quadruplex</div>
            <div style='margin-bottom:0.2rem;'><b>7.</b> i-Motif</div>
            <div style='margin-bottom:0.2rem;'><b>8.</b> Z-DNA</div>
            <div style='margin-bottom:0.2rem;'><b>9.</b> A-philic DNA</div>
            <div style='margin-bottom:0.2rem;'><b>10.</b> Hybrid</div>
            <div style='margin-bottom:0.2rem;'><b>11.</b> Non-B DNA Clusters</div>
        </div>
        <p style='margin-top:0.8rem; margin-bottom:0; color:#1976d2; font-weight:500; font-size:0.9rem;'>Upload FASTA files to begin analysis</p>
        </div>
        """, unsafe_allow_html=True)

# Updated Streamlit layout: Input Method + Sequence Preview + Analysis side-by-side
#  "Upload & Analyze" 
# It expects helper functions and constants to exist in the module:
# parse_fasta, get_basic_stats, EXAMPLE_FASTA, EXAMPLE_MULTI_FASTA, parse_fasta, all_motifs_refactored,
# ensure_subclass, MOTIF_ORDER, wrap, Entrez, SeqIO, Counter, pd, st

with tab_pages["Upload & Analyze"]:
    st.markdown("<h2>Sequence Upload and Motif Analysis</h2>", unsafe_allow_html=True)
   # st.markdown('<span style="font-family:Montserrat,Arial; font-size:0.95rem;">Supports multi-FASTA and single FASTA. Paste, upload, select example, or fetch from NCBI.</span>', unsafe_allow_html=True)
    st.caption("Supported formats: .fa, .fasta, .txt, .fna | Limit: 1GB/file (optimized for large genomic sequences).")

    # ----- TWO-COLUMN LAYOUT: Left for Upload, Right for Analysis -----
    col1, col2 = st.columns([1.2, 1])
    
    with col1:
        # LEFT COLUMN: Sequence Upload and Motif Analysis
        st.markdown("### 📤 Sequence Upload")
        
        # ----- Input Method -----
        input_method = st.radio("Choose your input method:",
                                ["Upload FASTA File", "Paste Sequence", "Example Data", "NCBI Fetch"],
                                horizontal=True,
                                label_visibility="collapsed",
                                key="upload_method")

        seqs, names = [], []

        if input_method == "Upload FASTA File":
            fasta_file = st.file_uploader("Drag and drop FASTA/multi-FASTA file here", 
                                         type=["fa", "fasta", "txt", "fna"],
                                         label_visibility="visible",
                                         help="Upload a FASTA file containing DNA sequences")
            if fasta_file:
                # Compact file card after upload
                file_size_mb = fasta_file.size / (1024 * 1024)
                
                # Memory-efficient processing with progress indicator
                with st.spinner(f"Processing {fasta_file.name}..."):
                    # Get preview first (lightweight operation)
                    preview_info = get_file_preview(fasta_file, max_sequences=3)
                    
                    # Compact File Card
                    st.markdown(f"""
                    <div style='background: linear-gradient(135deg, #4A90E2 0%, #6AA5F2 100%); 
                                border-radius: 12px; padding: 12px; margin: 8px 0; color: white; 
                                box-shadow: 0 2px 8px rgba(74, 144, 226, 0.15);'>
                        <div style='display: flex; justify-content: space-between; align-items: center;'>
                            <div>
                                <div style='font-weight: 600; font-size: 0.95rem;'>📁 {fasta_file.name}</div>
                                <div style='font-size: 0.85rem; opacity: 0.9; margin-top: 4px;'>
                                    {preview_info['num_sequences']} sequences | {preview_info['total_bp']:,} bp | {file_size_mb:.2f} MB
                                </div>
                            </div>
                            <div style='background: rgba(255,255,255,0.2); border-radius: 8px; padding: 8px 12px; font-weight: 600;'>
                                ✓ Valid
                            </div>
                        </div>
                    </div>
                    """, unsafe_allow_html=True)
                    
                    # Show preview of first few sequences using popover for better UX
                    with st.popover("🔍 Preview Sequences"):
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
                        st.warning("No sequences found in file.")

        elif input_method == "Paste Sequence":
            seq_input = st.text_area("Paste single or multi-FASTA here:", 
                                    height=150, 
                                    placeholder="Paste your DNA sequence(s) here...",
                                    help="Paste DNA sequences in FASTA format")
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
                                <div style='font-weight: 600; font-size: 0.95rem;'>📝 Pasted Sequences</div>
                                <div style='font-size: 0.85rem; opacity: 0.9; margin-top: 4px;'>
                                    {len(seqs)} sequences | {total_bp:,} bp
                                </div>
                            </div>
                            <div style='background: rgba(255,255,255,0.2); border-radius: 8px; padding: 8px 12px; font-weight: 600;'>
                                ✓ Valid
                            </div>
                        </div>
                    </div>
                    """, unsafe_allow_html=True)
                else:
                    st.warning("No sequences found.")

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
                    st.success("✅ Single example sequence loaded.")
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
                    st.success(f"✅ Multi-FASTA example loaded with {len(seqs)} sequences.")

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
                                st.success(f"✅ Fetched {len(seqs)} sequences.")
                        except Exception as e:
                            st.error(f"NCBI fetch failed: {e}")
                else:
                    st.warning("Enter a query before fetching.")

        # Persist sequences to session state if any found from input
        if seqs:
            st.session_state.seqs = seqs
            st.session_state.names = names
            st.session_state.results = []

        # Compact sequence validation indicator using popover for cleaner UI
        if st.session_state.get('seqs'):
            with st.popover("✅ Validation Summary"):
                for i, seq in enumerate(st.session_state.seqs[:3]):
                    stats = get_basic_stats(seq)
                    st.markdown(f"**{st.session_state.names[i]}** ({len(seq):,} bp)")
                    st.caption(f"GC: {stats['GC%']}% | AT: {stats['AT%']}%")
                if len(st.session_state.seqs) > 3:
                    st.caption(f"...and {len(st.session_state.seqs)-3} more.")
    
    with col2:
        # RIGHT COLUMN: Analysis & Run
        st.markdown("### ⚙️ Analysis & Run")
        
        # Quick Options Section
        st.markdown("##### Quick Options")
        detailed_output = st.checkbox("Detailed Analysis", value=True, 
                                    help="Include comprehensive motif metadata")
        quality_check = st.checkbox("Quality Validation", value=True, 
                                   help="Validate detected motifs")
        
        # Advanced Options - now visible by default
        show_chunk_progress = st.checkbox("Show Chunk-Level Progress", value=False,
                                         help="Display detailed progress for each processing chunk")
        use_parallel_scanner = st.checkbox("Use Experimental Parallel Scanner", value=True,
                                          help="Enable experimental parallel chunk-based scanner (>100kb sequences)")
        
        if use_parallel_scanner:
            st.caption("⚠️ Parallel scanner works best on sequences >100kb with multiple CPU cores")
        
        # Hardcoded default overlap handling
        nonoverlap = True
        overlap_option = "Remove overlaps within subclasses"
        
        # Helper text
        st.caption("💡 All 11 motif classes with 22+ subclasses are detected automatically")
    
    # ----- FULL-WIDTH STICKY RUN BUTTON -----
    st.markdown("---")
    
    # Check if valid input is present
    has_valid_input = bool(st.session_state.get('seqs'))
    
    # Create a full-width container for the run button
    run_button_container = st.container()
    with run_button_container:
        # Sticky-ish styling with disabled state
        if has_valid_input:
            run_button = st.button(
                "🚀 Run NBDScanner Analysis",
                type="primary",
                use_container_width=True,
                key="run_motif_analysis_main",
                help="Start analyzing uploaded sequences for Non-B DNA motifs"
            )
        else:
            # Disabled button appearance with accessibility
            st.markdown(f"""
            <div role="button" aria-disabled="true" aria-label="Run NBDScanner Analysis - Disabled: Please upload or paste a valid sequence first"
                 style='background: #e0e0e0; color: #9e9e9e; padding: 12px; 
                        border-radius: 12px; text-align: center; font-weight: 600;
                        font-size: 1.1rem; cursor: not-allowed; opacity: 0.6;'>
                🚀 Run NBDScanner Analysis (Disabled)
            </div>
            <p style='text-align: center; color: #9e9e9e; font-size: 0.85rem; margin-top: 8px;' role="status">
                ⚠️ Please upload or paste a valid sequence first
            </p>
            """, unsafe_allow_html=True)
            run_button = False
        
        # Placeholder for progress area
        progress_placeholder = st.empty()
    
    # ========== RUN ANALYSIS BUTTON LOGIC ========== 
    if run_button:
        # Simplified validation
        if not st.session_state.seqs:
            st.error("Please upload or input sequences before running analysis.")
            st.session_state.analysis_status = "Error"
        else:
            # Sequence length limit has been removed - the system now uses automatic chunking
            # (see nonbscanner.py CHUNK_THRESHOLD=10,000 bp) to handle sequences of any size
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
                        validation_messages.append(f"✓ {class_id}: Length limits {limits}")
            
            # Enhanced progress tracking with timer
            import time
            
            # Create placeholder for timer and progress
            timer_placeholder = st.empty()
            progress_placeholder = st.empty()
            status_placeholder = st.empty()
            detailed_progress_placeholder = st.empty()
            
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
            
            # Helper function to display progress using Streamlit native components
            def display_progress_panel(container, elapsed, estimated_remaining, progress_display, 
                                     status_text, seq_name, seq_bp, seq_num, total_seqs, 
                                     processed_bp, total_bp, detector_count, extra_info=""):
                """Display progress panel using Streamlit native components instead of HTML.
                
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
                    st.subheader("🧬 NonBScanner Analysis")
                    st.write(status_text)
                    
                    # Display metrics in 4 columns
                    col1, col2, col3, col4 = st.columns(4)
                    
                    with col1:
                        st.metric("⏱️ Elapsed", format_time(elapsed))
                    
                    with col2:
                        st.metric("⏳ Remaining", format_time(estimated_remaining))
                    
                    with col3:
                        st.metric("📊 Progress", progress_display)
                    
                    with col4:
                        st.metric("🔬 Detectors", str(detector_count))
                    
                    # Sequence info
                    st.write(f"**📄 Sequence {seq_num}/{total_seqs}**: {seq_name} *({seq_bp:,} bp)*")
                    st.write(f"⚡ Processed: {processed_bp:,} / {total_bp:,} bp")
                    
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
                    st.subheader("🔬 Analysis Pipeline")
                    
                    # Display detectors in a clean list format
                    for j, (detector_name, detector_desc) in enumerate(DETECTOR_PROCESSES):
                        st.write(f"**{j+1}. {detector_name}** - {detector_desc}")
                    
                    st.info("✓ All detectors process in parallel | 🔄 Followed by overlap resolution & clustering")
                    
                for i, (seq, name) in enumerate(zip(st.session_state.seqs, st.session_state.names)):
                    progress = (i + 1) / len(st.session_state.seqs)
                    
                    # Calculate elapsed time and estimated remaining
                    elapsed = time.time() - start_time
                    
                    # Calculate estimated remaining time (ensure non-negative)
                    if elapsed > 0 and total_bp_processed > 0:
                        current_rate = total_bp_processed / elapsed
                        remaining_bp = total_bp_all_sequences - total_bp_processed
                        estimated_remaining = max(0, remaining_bp / current_rate) if current_rate > 0 else 0
                    else:
                        estimated_remaining = max(0, estimated_total_time - elapsed)
                    
                    # Calculate overall percentage
                    overall_percentage = (total_bp_processed / total_bp_all_sequences * 100) if total_bp_all_sequences > 0 else 0
                    
                    # Determine status text based on progress state
                    if total_bp_processed == 0:
                        status_text = "Starting analysis..."
                        progress_display = "Starting"
                    else:
                        status_text = "Analysis in progress..."
                        progress_display = f"{overall_percentage:.1f}%"
                    
                    # Build extra info for chunk progress if enabled
                    extra_info = ""
                    if show_chunk_progress and use_parallel_scanner:
                        est_chunks = max(1, (len(seq) + CHUNK_SIZE_FOR_PARALLEL - 1) // CHUNK_SIZE_FOR_PARALLEL)
                        extra_info = f"📦 Chunks: {est_chunks}"
                    
                    # Display progress using Streamlit native components
                    display_progress_panel(
                        timer_placeholder,
                        elapsed, estimated_remaining, progress_display,
                        status_text, name, len(seq), i+1, len(st.session_state.seqs),
                        total_bp_processed, total_bp_all_sequences, len(DETECTOR_PROCESSES),
                        extra_info
                    )
                    
                    # Run the analysis - use parallel scanner for large sequences if enabled
                    seq_start = time.time()
                    
                    if use_parallel_scanner and len(seq) > 100000:
                        # Use experimental parallel scanner for large sequences
                        try:
                            from scanner_agent import ParallelScanner
                            
                            # Create chunk progress placeholder
                            chunk_progress_placeholder = st.empty()
                            
                            def chunk_progress_callback(current, total):
                                """Callback to update chunk progress"""
                                if show_chunk_progress:
                                    chunk_percent = (current / total) * 100
                                    chunk_progress_placeholder.info(f"📦 Chunks: {current}/{total} ({chunk_percent:.1f}%)")
                            
                            # Run parallel scanner
                            scanner = ParallelScanner(seq, hs_db=None)
                            raw_motifs = scanner.run_scan(progress_callback=chunk_progress_callback)
                            
                            # Convert raw motifs to full motif format by running through scoring
                            # For now, just use the standard analyzer
                            results = analyze_sequence(seq, name)
                            
                            # Clear chunk progress
                            if show_chunk_progress:
                                chunk_progress_placeholder.success(f"✅ Chunks complete: {len(raw_motifs)} motifs")
                            
                        except Exception as e:
                            st.warning(f"Parallel scanner failed, falling back to standard: {e}")
                            results = analyze_sequence(seq, name)
                    else:
                        # Use standard consolidated NBDScanner analysis
                        results = analyze_sequence(seq, name)
                    
                    seq_time = time.time() - seq_start
                    
                    # Ensure all motifs have required fields
                    results = [ensure_subclass(motif) for motif in results]
                    all_results.append(results)
                    
                    total_bp_processed += len(seq)
                    
                    # Calculate processing speed and update elapsed time
                    elapsed = time.time() - start_time
                    speed = total_bp_processed / elapsed if elapsed > 0 else 0
                    
                    # Recalculate estimated remaining time with actual speed
                    if speed > 0:
                        remaining_bp = total_bp_all_sequences - total_bp_processed
                        estimated_remaining = max(0, remaining_bp / speed)
                    else:
                        estimated_remaining = 0
                    
                    # Calculate actual progress percentage
                    actual_percentage = (total_bp_processed / total_bp_all_sequences * 100) if total_bp_all_sequences > 0 else 0
                    
                    # Build extra info for completion
                    completion_info = f"⚡ Total: {total_bp_processed:,} / {total_bp_all_sequences:,} bp | 🎯 {len(results)} motifs | 🚀 {speed:,.0f} bp/s"
                    
                    # Update timer display with actual progress using native components
                    display_progress_panel(
                        timer_placeholder,
                        elapsed, estimated_remaining, f"{actual_percentage:.1f}%",
                        f"✅ Sequence {i+1}/{len(st.session_state.seqs)} completed",
                        name, len(seq), i+1, len(st.session_state.seqs),
                        total_bp_processed, total_bp_all_sequences, len(DETECTOR_PROCESSES),
                        completion_info
                    )
                    
                    with progress_placeholder.container():
                        pbar.progress(progress, text=f"🔬 Analyzed {i+1}/{len(st.session_state.seqs)} sequences")
                    
                    status_placeholder.success(f"✅ {name}: {len(seq):,} bp in {seq_time:.2f}s ({len(seq)/seq_time:.0f} bp/s) | 🎯 {len(results)} motifs")
                
                # Store results
                st.session_state.results = all_results
                
                # Final timing statistics
                total_time = time.time() - start_time
                overall_speed = total_bp_processed / total_time if total_time > 0 else 0
                
                # ============================================================
                # RIGOROUS VALIDATION & QUALITY CHECKS
                # ============================================================
                with status_placeholder.container():
                    st.info("🔍 Validating results for consistency and quality...")
                
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
                        validation_issues.append(f"⚠️ Sequence {i+1}: {duplicates_found} duplicate motifs found")
                
                # 2. Validate motif data consistency
                for i, results in enumerate(all_results):
                    for motif in results:
                        # Check required fields
                        if not all(k in motif for k in ['Start', 'End', 'Class']):
                            validation_issues.append(f"⚠️ Sequence {i+1}: Motif missing required fields")
                            break
                        
                        # Validate positions
                        if motif.get('Start', 0) >= motif.get('End', 0):
                            validation_issues.append(f"⚠️ Sequence {i+1}: Invalid motif position (Start >= End)")
                            break
                        
                        # Validate length consistency
                        calculated_length = motif.get('End', 0) - motif.get('Start', 0)
                        if motif.get('Length') and abs(motif.get('Length') - calculated_length) > 1:
                            validation_issues.append(f"⚠️ Sequence {i+1}: Length mismatch detected")
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
                                validation_issues.append(f"⚠️ Sequence {i+1}: Overlapping motifs in {subclass}")
                                break
                
                # Display validation results
                if validation_issues:
                    with status_placeholder.container():
                        st.warning(f"Validation found {len(validation_issues)} potential issues:")
                        for issue in validation_issues[:5]:  # Show first 5
                            st.write(issue)
                        if len(validation_issues) > 5:
                            st.write(f"... and {len(validation_issues) - 5} more")
                else:
                    with status_placeholder.container():
                        st.success("✅ Validation passed: No consistency issues found")
                
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
                
                st.session_state.summary_df = pd.DataFrame(summary)
                
                # ============================================================
                # PRE-GENERATE ALL VISUALIZATIONS FOR CLASSES AND SUBCLASSES
                # ============================================================
                with status_placeholder.container():
                    st.info("📊 Generating comprehensive visualizations for all classes and subclasses...")
                
                viz_progress = st.empty()
                viz_status = st.empty()
                
                # Cache all visualizations for each sequence
                st.session_state.cached_visualizations = {}
                
                viz_start_time = time.time()
                total_viz_count = 0
                
                for seq_idx, (seq, name, motifs) in enumerate(zip(st.session_state.seqs, st.session_state.names, all_results)):
                    sequence_length = len(seq)
                    
                    # Filter motifs (exclude hybrid/cluster for main visualizations)
                    filtered_motifs = [m for m in motifs if m.get('Class') not in ['Hybrid', 'Non-B_DNA_Clusters']]
                    
                    if not filtered_motifs:
                        continue
                    
                    viz_cache_key = f"seq_{seq_idx}"
                    st.session_state.cached_visualizations[viz_cache_key] = {}
                    
                    # Calculate all class and subclass statistics
                    with viz_status.container():
                        st.write(f"Analyzing sequence {seq_idx + 1}/{len(st.session_state.seqs)}: {name}")
                    
                    # Pre-calculate all density metrics (class and subclass level)
                    try:
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
                        
                        total_viz_count += 4  # Count density calculations
                        
                    except Exception as e:
                        with viz_status.container():
                            st.warning(f"Could not calculate density metrics for {name}: {e}")
                    
                    # Count unique classes and subclasses
                    unique_classes = len(set(m.get('Class', 'Unknown') for m in filtered_motifs))
                    unique_subclasses = len(set(m.get('Subclass', 'Unknown') for m in filtered_motifs))
                    
                    st.session_state.cached_visualizations[viz_cache_key]['summary'] = {
                        'unique_classes': unique_classes,
                        'unique_subclasses': unique_subclasses,
                        'total_motifs': len(filtered_motifs)
                    }
                    
                    with viz_progress.container():
                        viz_elapsed = time.time() - viz_start_time
                        st.write(f"✅ Processed {seq_idx + 1}/{len(st.session_state.seqs)} sequences")
                        st.write(f"📊 Generated {total_viz_count} visualization components in {viz_elapsed:.2f}s")
                
                viz_total_time = time.time() - viz_start_time
                
                with viz_status.container():
                    st.success(f"✅ All visualizations prepared: {total_viz_count} components in {viz_total_time:.2f}s")
                
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
                viz_progress.empty()
                viz_status.empty()
                
                # Show final success message with enhanced performance metrics
                timer_placeholder.markdown(f"""
                <div class='progress-panel progress-panel--success'>
                    <h3 class='progress-panel__title'>🎉 Analysis Complete!</h3>
                    <p class='progress-panel__subtitle'>All detectors, validations, and visualizations completed successfully</p>
                    <div class='stats-grid stats-grid--wide'>
                        <div class='stat-card'>
                            <h2 class='stat-card__value'>⏱️ {total_time:.2f}s</h2>
                            <p class='stat-card__label'>Analysis Time</p>
                        </div>
                        <div class='stat-card'>
                            <h2 class='stat-card__value'>📊 {viz_total_time:.2f}s</h2>
                            <p class='stat-card__label'>Visualization Time</p>
                        </div>
                        <div class='stat-card'>
                            <h2 class='stat-card__value'>🧬 {total_bp_processed:,}</h2>
                            <p class='stat-card__label'>Base Pairs</p>
                        </div>
                        <div class='stat-card'>
                            <h2 class='stat-card__value'>🚀 {overall_speed:,.0f}</h2>
                            <p class='stat-card__label'>bp/second</p>
                        </div>
                        <div class='stat-card'>
                            <h2 class='stat-card__value'>🔬 {len(DETECTOR_PROCESSES)}</h2>
                            <p class='stat-card__label'>Detectors</p>
                        </div>
                        <div class='stat-card'>
                            <h2 class='stat-card__value'>🎯 {sum(len(r) for r in all_results)}</h2>
                            <p class='stat-card__label'>Motifs Found</p>
                        </div>
                        <div class='stat-card'>
                            <h2 class='stat-card__value'>📈 {total_viz_count}</h2>
                            <p class='stat-card__label'>Viz Components</p>
                        </div>
                        <div class='stat-card'>
                            <h2 class='stat-card__value'>✅ {len(validation_issues)}</h2>
                            <p class='stat-card__label'>Validation Issues</p>
                        </div>
                    </div>
                </div>
                """, unsafe_allow_html=True)
                
                # Show comprehensive completion summary
                completion_msg = f"""✅ **Analysis Complete!** All processing stages finished successfully:
                
**Detection & Analysis:**
- 🔬 {len(DETECTOR_PROCESSES)} detector processes completed
- 🎯 {sum(len(r) for r in all_results)} total motifs detected across {len(st.session_state.seqs)} sequences
- ⏱️ Analysis completed in {total_time:.2f}s ({overall_speed:,.0f} bp/s)

**Quality Validation:**
- ✅ Data consistency checks: {'PASSED' if len(validation_issues) == 0 else f'{len(validation_issues)} issues found'}
- ✅ Non-redundancy validation: Complete
- ✅ Position validation: Complete

**Visualization Generation:**
- 📊 {total_viz_count} visualization components pre-generated
- 📈 Class-level and subclass-level analysis ready
- ⏱️ Visualizations prepared in {viz_total_time:.2f}s

**View detailed results in the 'Analysis Results and Visualization' tab.**
"""
                st.success(completion_msg)
                st.session_state.analysis_status = "Complete"
                
            except Exception as e:
                timer_placeholder.empty()
                progress_placeholder.empty()
                status_placeholder.empty()
                detailed_progress_placeholder.empty()
                st.error(f"Analysis failed: {str(e)}")
                st.session_state.analysis_status = "Error"

    # End of Upload & Analyze tab
    st.markdown("---")
# ---------- RESULTS ----------
with tab_pages["Results"]:
    st.markdown('<h2>Analysis Results and Visualization</h2>', unsafe_allow_html=True)
    if not st.session_state.results:
        st.info("No analysis results. Please run motif analysis first.")
    else:
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
        st.markdown("### Analysis Summary")
        st.dataframe(st.session_state.summary_df, use_container_width=True)
        
        # Sequence selection for detailed analysis using pills for better UX
        seq_idx = 0
        if len(st.session_state.seqs) > 1:
            # Use pills for sequence selection - a more modern and visual alternative to dropdown
            selected_seq = st.pills(
                "Choose Sequence for Details:",
                options=list(range(len(st.session_state.seqs))),
                format_func=lambda i: st.session_state.names[i],
                selection_mode="single",
                default=0,
                help="Select a sequence to view detailed analysis results"
            )
            seq_idx = selected_seq or 0
        
        motifs = st.session_state.results[seq_idx]
        sequence_length = len(st.session_state.seqs[seq_idx])
        sequence_name = st.session_state.names[seq_idx]
        
        if not motifs:
            st.warning("No motifs detected for this sequence.")
        else:
            # Filter motifs for main display (exclude hybrid and cluster)
            filtered_motifs = [m for m in motifs if m.get('Class') not in ['Hybrid', 'Non-B_DNA_Clusters']]
            hybrid_cluster_motifs = [m for m in motifs if m.get('Class') in ['Hybrid', 'Non-B_DNA_Clusters']]
            
            # Create enhanced motifs DataFrame (only non-hybrid/cluster motifs)
            df = pd.DataFrame(filtered_motifs) if filtered_motifs else pd.DataFrame()
            
            # Calculate and display enhanced coverage statistics (using only filtered motifs)
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
                        <h2 class='stat-card__value stat-card__value--large'>
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
                st.info(f"{hybrid_cluster_count} Hybrid/Cluster motifs detected. View them in the 'Cluster/Hybrid' tab below.")
            
            # Show cached visualization summary if available
            viz_cache_key = f"seq_{seq_idx}"
            cached_viz = st.session_state.get('cached_visualizations', {}).get(viz_cache_key, {})
            if cached_viz.get('summary'):
                viz_summary = cached_viz['summary']
                st.success(f"""📊 **Pre-generated Analysis Ready:** 
                {viz_summary['unique_classes']} unique classes, 
                {viz_summary['unique_subclasses']} unique subclasses analyzed
                """)
            
            # Enhanced motif table with new columns and pagination for large datasets
            # Display table without header as per requirements
            
            # Column selection for display using pills for a more modern, visual approach
            available_columns = df.columns.tolist()
            # Filter out specific Normalized_Score column
            available_columns = [col for col in available_columns if col != 'Normalized_Score']
            
            # Use pills for column selection - multi-selection mode for better UX
            default_cols = [col for col in ['Class', 'Subclass', 'Start', 'End', 'Length', 'Sequence', 'Score'] if col in available_columns]
            display_columns = st.pills(
                "Select columns to display:",
                options=available_columns,
                selection_mode="multi",
                default=default_cols,
                help="Click to toggle columns on/off for display"
            )
            # Ensure display_columns is always a list (st.pills with multi mode returns list or None)
            if display_columns is None:
                display_columns = []
            
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
                # Filter out sequence and sequence name columns for cleaner display
                filtered_df = df_page[display_columns].copy()
                # Replace underscores with spaces in column names for display
                filtered_df.columns = [col.replace('_', ' ') for col in filtered_df.columns]
                st.dataframe(filtered_df, use_container_width=True, height=360)
            else:
                # Replace underscores with spaces in column names for display
                display_df = df_page.copy()
                display_df.columns = [col.replace('_', ' ') for col in display_df.columns]
                st.dataframe(display_df, use_container_width=True, height=360)
            
            # CONSOLIDATED VISUALIZATION SUITE
            st.markdown('<h3>Visualizations</h3>', unsafe_allow_html=True)
            
            # Create 4 visualization tabs
            viz_tabs = st.tabs(["Distribution & Statistics", "Coverage & Density", 
                               "Genome-Wide Analysis", "Cluster/Hybrid"])
            
            with viz_tabs[0]:  # Distribution & Statistics (merged)
                st.markdown("##### Motif Distribution")
                try:
                    # Class distribution
                    fig1 = plot_motif_distribution(filtered_motifs, by='Class', title=f"Motif Classes - {sequence_name}")
                    st.pyplot(fig1)
                    plt.close(fig1)
                    
                    # Subclass distribution
                    fig2 = plot_motif_distribution(filtered_motifs, by='Subclass', title=f"Motif Subclasses - {sequence_name}")
                    st.pyplot(fig2) 
                    plt.close(fig2)
                    
                    # Nested pie chart
                    fig3 = plot_nested_pie_chart(filtered_motifs, title=f"Class-Subclass Distribution - {sequence_name}")
                    st.pyplot(fig3)
                    plt.close(fig3)
                except Exception as e:
                    st.error(f"Error generating distribution plots: {e}")
                
                # Statistical Analysis section (from old Statistics tab)
                st.markdown("##### Statistical Analysis")
                
                # Density Metrics section
                st.markdown("###### Density Metrics")
                try:
                    # Check if we have cached density metrics from analysis
                    viz_cache_key = f"seq_{seq_idx}"
                    cached_viz = st.session_state.get('cached_visualizations', {}).get(viz_cache_key, {})
                    cached_densities = cached_viz.get('densities', {})
                    
                    if cached_densities:
                        # Use cached density calculations
                        genomic_density = cached_densities['class_genomic']
                        positional_density_kbp = cached_densities['class_positional']
                        genomic_density_subclass = cached_densities['subclass_genomic']
                        positional_density_subclass = cached_densities['subclass_positional']
                        st.info("📊 Using pre-calculated density metrics from analysis phase")
                    else:
                        # Calculate fresh if not cached
                        genomic_density = calculate_genomic_density(filtered_motifs, sequence_length, by_class=True)
                        positional_density_kbp = calculate_positional_density(filtered_motifs, sequence_length, unit='kbp', by_class=True)
                        
                        genomic_density_subclass = calculate_genomic_density(filtered_motifs, sequence_length, 
                                                                            by_class=False, by_subclass=True)
                        positional_density_subclass = calculate_positional_density(filtered_motifs, sequence_length, 
                                                                                  unit='kbp', by_class=False, by_subclass=True)
                    
                    positional_density_mbp = calculate_positional_density(filtered_motifs, sequence_length, unit='Mbp', by_class=True)
                    
                    # Overall metrics in compact layout
                    col1, col2 = st.columns(2)
                    with col1:
                        st.metric("Genomic Density", f"{genomic_density.get('Overall', 0):.4f}%")
                    with col2:
                        st.metric("Motifs/kbp", f"{positional_density_kbp.get('Overall', 0):.2f}")
                    
                    # Density Analysis Level Selection
                    density_level = st.radio(
                        "Density Analysis Level:",
                        ["Class Level", "Subclass Level"],
                        horizontal=True,
                        help="Choose whether to analyze density at class or subclass level"
                    )
                    
                    if density_level == "Class Level":
                        # Per-class density table
                        density_data = []
                        for class_name in sorted([k for k in genomic_density.keys() if k != 'Overall']):
                            density_data.append({
                                'Motif Class': class_name.replace('_', ' '),
                                'Genomic Density (%)': f"{genomic_density.get(class_name, 0):.4f}",
                                'Motifs/kbp': f"{positional_density_kbp.get(class_name, 0):.2f}"
                            })
                        
                        if density_data:
                            density_df = pd.DataFrame(density_data)
                            st.dataframe(density_df, use_container_width=True, height=200)
                        
                        fig_density = plot_density_comparison(genomic_density, positional_density_kbp,
                                                              title="Motif Density Analysis (Class Level)")
                        st.pyplot(fig_density)
                        plt.close(fig_density)
                    
                    else:  # Subclass Level
                        # Use cached subclass densities if available
                        if not cached_densities:
                            # Calculate subclass-level densities only if not cached
                            genomic_density_subclass = calculate_genomic_density(filtered_motifs, sequence_length, 
                                                                                 by_class=False, by_subclass=True)
                            positional_density_subclass = calculate_positional_density(filtered_motifs, sequence_length, 
                                                                                       unit='kbp', by_class=False, by_subclass=True)
                        
                        # Per-subclass density table
                        subclass_density_data = []
                        for subclass_key in sorted([k for k in genomic_density_subclass.keys() if k != 'Overall']):
                            # Format the subclass key for display
                            display_key = subclass_key.replace('_', ' ').replace(':', ': ')
                            subclass_density_data.append({
                                'Motif Subclass': display_key,
                                'Genomic Density (%)': f"{genomic_density_subclass.get(subclass_key, 0):.4f}",
                                'Motifs/kbp': f"{positional_density_subclass.get(subclass_key, 0):.2f}"
                            })
                        
                        if subclass_density_data:
                            subclass_density_df = pd.DataFrame(subclass_density_data)
                            st.dataframe(subclass_density_df, use_container_width=True, height=300)
                        
                        # Subclass density visualization
                        fig_density_subclass = plot_density_comparison_by_subclass(
                            genomic_density_subclass, positional_density_subclass,
                            title="Motif Density Analysis (Subclass Level)"
                        )
                        st.pyplot(fig_density_subclass)
                        plt.close(fig_density_subclass)
                        
                        # Additional subclass density heatmap
                        if len(filtered_motifs) > 0:
                            st.markdown("**Subclass Density Heatmap Along Sequence**")
                            fig_heatmap = plot_subclass_density_heatmap(
                                filtered_motifs, sequence_length,
                                window_size=max(500, sequence_length // 50),
                                title="Subclass Density Distribution"
                            )
                            st.pyplot(fig_heatmap)
                            plt.close(fig_heatmap)
                    
                except Exception as e:
                    st.error(f"Error calculating density metrics: {e}")
                    with st.expander("🔍 Show detailed error trace"):
                        import traceback
                        st.code(traceback.format_exc(), language="python")
                
                # Length/Score distributions
                st.markdown("###### Distributions")
                try:
                    fig6 = plot_length_distribution(filtered_motifs, by_class=True, 
                                                   title="Length Distribution by Motif Class") 
                    st.pyplot(fig6)
                    plt.close(fig6)
                    
                    fig7 = plot_score_distribution(filtered_motifs, by_class=True,
                                                  title="Score Distribution by Motif Class (1-3 Scale)")
                    st.pyplot(fig7)
                    plt.close(fig7)
                except Exception as e:
                    st.error(f"Error generating distribution plots: {e}")
            
            with viz_tabs[1]:  # Coverage & Density (merged Coverage Map and Circos)
                st.markdown("##### Sequence Coverage Analysis")
                try:
                    # Coverage map
                    fig4 = plot_coverage_map(filtered_motifs, sequence_length, title=f"Motif Coverage - {sequence_name}")
                    st.pyplot(fig4)
                    plt.close(fig4)
                    
                    # Density heatmap
                    fig5 = plot_density_heatmap(filtered_motifs, sequence_length, 
                                               window_size=max(100, sequence_length // 20),
                                               title=f"Motif Density - {sequence_name}")
                    st.pyplot(fig5)
                    plt.close(fig5)
                    
                    # Circos plot (only one - the main one)
                    st.markdown("##### Circos Density Plot")
                    fig_circos = plot_circos_motif_density(filtered_motifs, sequence_length,
                                                          title=f"Circos Density - {sequence_name}")
                    st.pyplot(fig_circos)
                    plt.close(fig_circos)
                    
                except Exception as e:
                    st.error(f"Error generating coverage plots: {e}")
            
            with viz_tabs[2]:  # Genome-Wide Visualizations (was tab 3)
                st.markdown("##### Genome-Wide Motif Analysis")
                st.info("📊 Publication-quality genome-scale visualizations showing motif distribution patterns across the entire sequence.")
                
                try:
                    # Manhattan plot for motif density hotspots
                    st.markdown("**Manhattan Plot - Motif Density Hotspots**")
                    fig_manhattan = plot_manhattan_motif_density(
                        filtered_motifs, sequence_length,
                        title=f"Manhattan Plot - {sequence_name}"
                    )
                    st.pyplot(fig_manhattan)
                    plt.close(fig_manhattan)
                    
                    # Cumulative motif distribution
                    st.markdown("**Cumulative Motif Distribution**")
                    fig_cumulative = plot_cumulative_motif_distribution(
                        filtered_motifs, sequence_length,
                        title=f"Cumulative Distribution - {sequence_name}",
                        by_class=True
                    )
                    st.pyplot(fig_cumulative)
                    plt.close(fig_cumulative)
                    
                    # Linear motif track (for smaller regions)
                    if sequence_length <= 50000:  # Only for sequences < 50kb
                        st.markdown("**Linear Motif Track Viewer**")
                        fig_linear = plot_linear_motif_track(
                            filtered_motifs, sequence_length,
                            title=f"Motif Track - {sequence_name}"
                        )
                        st.pyplot(fig_linear)
                        plt.close(fig_linear)
                    else:
                        st.info("💡 Linear motif track is shown for sequences < 50kb. For large sequences, use Manhattan plot or Coverage map.")
                    
                except Exception as e:
                    st.error(f"Error generating genome-wide plots: {e}")
                    with st.expander("🔍 Show detailed error trace"):
                        import traceback
                        st.code(traceback.format_exc(), language="python")
            
            with viz_tabs[3]:  # Cluster/Hybrid & Advanced Visualizations (merged from tabs 4 and 5)
                st.markdown("##### Hybrid & Cluster Motif Analysis")
                
                if hybrid_cluster_motifs:
                    # Separate hybrid and cluster motifs
                    hybrid_only = [m for m in hybrid_cluster_motifs if m.get('Class') == 'Hybrid']
                    cluster_only = [m for m in hybrid_cluster_motifs if m.get('Class') == 'Non-B_DNA_Clusters']
                    
                    hc_df = pd.DataFrame(hybrid_cluster_motifs)
                    
                    # Summary metrics with detailed breakdown
                    st.markdown(f"""
                    <div class='progress-panel progress-panel--hybrid'>
                        <h3 class='progress-panel__title progress-panel__title--large'>
                            Hybrid & Cluster Motif Summary
                        </h3>
                        <div class='stats-grid stats-grid--extra-wide'>
                            <div class='stat-card stat-card--large'>
                                <h2 class='stat-card__value stat-card__value--large'>
                                    {len(hybrid_only)}
                                </h2>
                                <p class='stat-card__label stat-card__label--large'>
                                    Hybrid Motifs
                                </p>
                            </div>
                            <div class='stat-card stat-card--large'>
                                <h2 class='stat-card__value stat-card__value--large'>
                                    {len(cluster_only)}
                                </h2>
                                <p class='stat-card__label stat-card__label--large'>
                                    DNA Clusters
                                </p>
                            </div>
                            <div class='stat-card stat-card--large'>
                                <h2 class='stat-card__value stat-card__value--large'>
                                    {int(hc_df['Length'].mean()) if 'Length' in hc_df.columns else 0}
                                </h2>
                                <p class='stat-card__label stat-card__label--large'>
                                    Avg Length (bp)
                                </p>
                            </div>
                            <div class='stat-card stat-card--large'>
                                <h2 class='stat-card__value stat-card__value--large'>
                                    {len(hybrid_cluster_motifs)}
                                </h2>
                                <p class='stat-card__label stat-card__label--large'>
                                    Total
                                </p>
                            </div>
                        </div>
                    </div>
                    """, unsafe_allow_html=True)
                    
                    # Create sub-tabs for hybrid and cluster motifs
                    hc_subtabs = st.tabs(["All", "Hybrid Motifs", "Cluster Motifs"])
                    
                    with hc_subtabs[0]:  # All motifs
                        st.markdown("##### All Hybrid & Cluster Motifs")
                        display_cols = ['Class', 'Subclass', 'Start', 'End', 'Length', 'Score']
                        available_display_cols = [col for col in display_cols if col in hc_df.columns]
                        display_hc_df = hc_df[available_display_cols].copy()
                        display_hc_df.columns = [col.replace('_', ' ') for col in display_hc_df.columns]
                        st.dataframe(display_hc_df, use_container_width=True, height=300)
                        
                        # Position map for all hybrid/cluster motifs (with height cap to prevent excessive memory usage)
                        fig_height = min(20, max(4, len(hybrid_cluster_motifs) * 0.3))
                        fig, ax = plt.subplots(figsize=(12, fig_height))
                        colors_map = {'Hybrid': '#ff6b6b', 'Non-B_DNA_Clusters': '#4ecdc4'}
                        for i, motif in enumerate(hybrid_cluster_motifs):
                            color = colors_map.get(motif.get('Class'), '#95a5a6')
                            ax.barh(i, motif['End'] - motif['Start'], left=motif['Start'], 
                                   height=0.8, color=color, alpha=0.7, edgecolor='black', linewidth=0.5)
                        ax.set_xlabel("Position (bp)")
                        ax.set_ylabel("Motif Index")
                        ax.set_title("Hybrid & Cluster Position Map")
                        handles = [plt.Rectangle((0,0),1,1, color=c, alpha=0.7) for c in colors_map.values()]
                        ax.legend(handles, list(colors_map.keys()), loc='upper right')
                        ax.grid(axis='x', alpha=0.3)
                        plt.tight_layout()
                        st.pyplot(fig)
                        plt.close(fig)
                    
                    with hc_subtabs[1]:  # Hybrid motifs only
                        st.markdown("##### Hybrid Motifs (Overlapping Different Classes)")
                        if hybrid_only:
                            hybrid_df = pd.DataFrame(hybrid_only)
                            
                            # Show component classes information
                            st.info(f"**Hybrid motifs** are regions where different Non-B DNA motif classes overlap. Found {len(hybrid_only)} hybrid regions.")
                            
                            # Extended columns for hybrid motifs
                            extended_cols = ['Subclass', 'Start', 'End', 'Length', 'Score', 'Component_Classes']
                            available_cols = [col for col in extended_cols if col in hybrid_df.columns]
                            display_df = hybrid_df[available_cols].copy()
                            display_df.columns = [col.replace('_', ' ') for col in display_df.columns]
                            st.dataframe(display_df, use_container_width=True, height=300)
                            
                            # Hybrid motif position visualization (with height cap)
                            fig_height = min(15, max(3, len(hybrid_only) * 0.4))
                            fig, ax = plt.subplots(figsize=(12, fig_height))
                            for i, motif in enumerate(hybrid_only):
                                ax.barh(i, motif['End'] - motif['Start'], left=motif['Start'], 
                                       height=0.8, color='#ff6b6b', alpha=0.7, edgecolor='#c0392b', linewidth=1)
                                # Add label with component classes
                                subclass = motif.get('Subclass', '')
                                label_text = f" {subclass[:30]}{'...' if len(subclass) > 30 else ''}" if subclass else ""
                                ax.text(motif['Start'], i, label_text, va='center', fontsize=8)
                            ax.set_xlabel("Position (bp)")
                            ax.set_ylabel("Hybrid Index")
                            ax.set_title("Hybrid Motif Positions (Overlapping Classes)")
                            ax.grid(axis='x', alpha=0.3)
                            plt.tight_layout()
                            st.pyplot(fig)
                            plt.close(fig)
                        else:
                            st.info("No hybrid motifs detected in this sequence.")
                    
                    with hc_subtabs[2]:  # Cluster motifs only
                        st.markdown("##### Non-B DNA Clusters (High-Density Regions)")
                        if cluster_only:
                            cluster_df = pd.DataFrame(cluster_only)
                            
                            # Show cluster information
                            st.info(f"**DNA Clusters** are high-density regions with multiple Non-B DNA motif classes. Found {len(cluster_only)} cluster regions.")
                            
                            # Extended columns for cluster motifs
                            extended_cols = ['Subclass', 'Start', 'End', 'Length', 'Score', 'Motif_Count', 'Class_Diversity']
                            available_cols = [col for col in extended_cols if col in cluster_df.columns]
                            display_df = cluster_df[available_cols].copy()
                            display_df.columns = [col.replace('_', ' ') for col in display_df.columns]
                            st.dataframe(display_df, use_container_width=True, height=300)
                            
                            # Cluster position visualization with diversity color coding (with height cap)
                            fig_height = min(15, max(3, len(cluster_only) * 0.4))
                            fig, ax = plt.subplots(figsize=(12, fig_height))
                            for i, motif in enumerate(cluster_only):
                                diversity = motif.get('Class_Diversity', 2)
                                # Color intensity based on diversity
                                alpha = min(0.4 + diversity * 0.1, 0.9)
                                ax.barh(i, motif['End'] - motif['Start'], left=motif['Start'], 
                                       height=0.8, color='#4ecdc4', alpha=alpha, edgecolor='#16a085', linewidth=1)
                                # Add label with motif count
                                motif_count = motif.get('Motif_Count', 0)
                                ax.text(motif['End'], i, f" {motif_count} motifs", va='center', fontsize=8)
                            ax.set_xlabel("Position (bp)")
                            ax.set_ylabel("Cluster Index")
                            ax.set_title("Non-B DNA Cluster Positions (Color intensity = Class diversity)")
                            ax.grid(axis='x', alpha=0.3)
                            plt.tight_layout()
                            st.pyplot(fig)
                            plt.close(fig)
                        else:
                            st.info("No DNA clusters detected in this sequence.")
                    
                    # Advanced Statistical Visualizations Section
                    st.markdown("##### Advanced Statistical Visualizations")
                    st.info("📈 Advanced publication-quality visualizations for in-depth analysis and manuscript figures.")
                    
                    try:
                        # Motif co-occurrence matrix
                        st.markdown("**Motif Co-occurrence Matrix**")
                        st.caption("Shows which motif classes tend to appear together (overlapping or within 1bp)")
                        fig_cooccur = plot_motif_cooccurrence_matrix(
                            filtered_motifs,
                            title=f"Co-occurrence Matrix - {sequence_name}"
                        )
                        st.pyplot(fig_cooccur)
                        plt.close(fig_cooccur)
                        
                        # GC content correlation (if sequence available)
                        if st.session_state.seqs[seq_idx]:
                            st.markdown("**GC Content vs Motif Density Correlation**")
                            st.caption("Scatter plot showing relationship between GC content and motif density")
                            fig_gc = plot_gc_content_correlation(
                                filtered_motifs, st.session_state.seqs[seq_idx],
                                title=f"GC Correlation - {sequence_name}"
                            )
                            st.pyplot(fig_gc)
                            plt.close(fig_gc)
                        
                        # Motif length KDE
                        st.markdown("**Motif Length Distribution (Kernel Density)**")
                        st.caption("Smooth probability density curves showing length patterns by class")
                        fig_kde = plot_motif_length_kde(
                            filtered_motifs,
                            by_class=True,
                            title=f"Length KDE - {sequence_name}"
                        )
                        st.pyplot(fig_kde)
                        plt.close(fig_kde)
                        
                        # Cluster size distribution (if clusters exist)
                        cluster_motifs_check = [m for m in filtered_motifs if m.get('Class') == 'Non-B_DNA_Clusters']
                        if cluster_motifs_check:
                            st.markdown("**Cluster Size & Diversity Distribution**")
                            st.caption("Distribution of motif counts and class diversity within clusters")
                            fig_cluster_dist = plot_cluster_size_distribution(
                                filtered_motifs,
                                title=f"Cluster Statistics - {sequence_name}"
                            )
                            st.pyplot(fig_cluster_dist)
                            plt.close(fig_cluster_dist)
                        
                    except Exception as e:
                        st.error(f"Error generating advanced plots: {e}")
                        with st.expander("🔍 Show detailed error trace"):
                            import traceback
                            st.code(traceback.format_exc(), language="python")
                else:
                    st.info("No hybrid or cluster motifs detected in this sequence. Hybrid motifs occur when different Non-B DNA classes overlap, and clusters form when multiple motifs are found in close proximity.")

# ---------- DOWNLOAD ----------
with tab_pages["Download"]:
    st.header("Export Data")
    if not st.session_state.results:
        st.info("No results available to download.")
    else:
        st.markdown("### Export Options")
        
        # Show analysis configuration used
        overlap_option_used = st.session_state.get('overlap_option_used', 'Remove overlaps within subclasses (default)')
        st.info(f"""
        **Analysis Configuration Used:**
        • Overlap Handling: {overlap_option_used}
        • Motif Classes: {len(st.session_state.get('selected_classes', []))} classes selected
        • Total Motifs Found: {sum(len(motifs) for motifs in st.session_state.results)}
        """)
        
        # Export configuration
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**Export Configuration**")
            
        with col2:
            include_sequences = st.checkbox("Include Full Sequences", value=True,
                                           help="Include full motif sequences in export")
        
        # Prepare data for export using consolidated utilities (exclude Hybrid and Cluster)
        all_motifs = []
        hybrid_cluster_motifs_export = []
        for i, motifs in enumerate(st.session_state.results):
            for m in motifs:
                export_motif = m.copy()
                if 'Sequence_Name' not in export_motif:
                    export_motif['Sequence_Name'] = st.session_state.names[i]
                
                # Separate hybrid/cluster from regular motifs
                if export_motif.get('Class') in ['Hybrid', 'Non-B_DNA_Clusters']:
                    hybrid_cluster_motifs_export.append(export_motif)
                else:
                    all_motifs.append(export_motif)
        
        # Info about excluded motifs
        if hybrid_cluster_motifs_export:
            st.info(f"{len(hybrid_cluster_motifs_export)} Hybrid/Cluster motifs are excluded from downloads. These are shown only in the Cluster/Hybrid visualization tab.")
        
        # Display preview
        if all_motifs:
            df_preview = export_results_to_dataframe(all_motifs).head(10)
            st.markdown("### Export Preview")
            st.dataframe(df_preview, use_container_width=True)
            st.caption(f"Showing first 10 of {len(all_motifs)} total records (Hybrid/Cluster motifs excluded)")
        
        # Export buttons using consolidated functions
        st.markdown("### Download Files")
        
        # Add helpful info about Excel format
        st.info("""
        **Excel Format**: Downloads a multi-sheet workbook with:
        • **Consolidated Sheet**: All non-overlapping motifs
        • **Class Sheets**: Separate sheets for each motif class (G-Quadruplex, Z-DNA, etc.)
        • **Subclass Sheets**: Detailed breakdown by subclass (when multiple subclasses exist)
        """)
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            # CSV Export
            if all_motifs:
                csv_data = export_to_csv(all_motifs)
                st.download_button(
                    "Download CSV", 
                    data=csv_data.encode('utf-8'), 
                    file_name="nbdscanner_results.csv", 
                    mime="text/csv",
                    use_container_width=True,
                    help="Comma-separated values format"
                )
        
        with col2:
            # Excel Export
            if all_motifs:
                try:
                    excel_bytes = generate_excel_bytes(all_motifs)
                    st.download_button(
                        "Download Excel", 
                        data=excel_bytes, 
                        file_name="nbdscanner_results.xlsx", 
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        use_container_width=True,
                        help="Excel format with multiple sheets (Consolidated + per-class)"
                    )
                except Exception as e:
                    st.error(f"Excel export error: {str(e)}")
        
        with col3:
            # JSON Export  
            if all_motifs:
                json_data = export_to_json(all_motifs, pretty=True)
                st.download_button(
                    "Download JSON", 
                    data=json_data.encode('utf-8'), 
                    file_name="nbdscanner_results.json", 
                    mime="application/json",
                    use_container_width=True,
                    help="JSON format with metadata"
                )
        
        with col4:
            # BED Export
            if all_motifs and st.session_state.names:
                bed_data = export_to_bed(all_motifs, st.session_state.names[0])
                st.download_button(
                    "Download BED", 
                    data=bed_data.encode('utf-8'), 
                    file_name="nbdscanner_results.bed", 
                    mime="text/plain",
                    use_container_width=True,
                    help="BED format for genome browsers"
                )
        
        # Additional export options in a new row
        st.markdown("### Additional Exports")
        col5, col6, col7, col8 = st.columns(4)
        
        with col5:
            # Config Export
            if all_motifs:
                import json
                # Create a basic configuration summary
                config_summary = {
                    "analysis_info": {
                        "tool": "NBDScanner",
                        "version": "2024.1",
                        "motif_classes_detected": list(set(m.get('Class', 'Unknown') for m in all_motifs)),
                        "total_sequences_analyzed": len(st.session_state.names) if hasattr(st.session_state, 'names') else 0,
                        "total_motifs_found": len(all_motifs)
                    },
                    "detection_parameters": {
                        "all_classes_enabled": True,
                        "consolidated_system": True
                    }
                }
                config_json = json.dumps(config_summary, indent=2)
                st.download_button(
                    "Download Config", 
                    data=config_json, 
                    file_name="analysis_configuration.json", 
                    mime="application/json",
                    use_container_width=True,
                    help="Analysis configuration and metadata"
                )
        
        # Prepare data for genome browser exports
        all_seq_data = []
        for i, motifs in enumerate(st.session_state.results):
            all_seq_data.append({
                'sequence_name': st.session_state.names[i],
                'motifs': motifs
            })


# ---------- DOCUMENTATION ----------
with tab_pages["Documentation"]:
    st.header("Scientific Documentation & References")
    
    # Motif detection parameters
    st.markdown("""
    <div style='background:#e3f2fd; border-radius:12px; padding:18px; font-size:1.08rem; font-family:Montserrat,Arial; margin-bottom:20px;'>
    <b>Motif Detection Parameters:</b><br><br>
    <table style='width:100%; border-collapse: collapse;'>
        <tr style='background:#1976d2; color:white;'>
            <th style='padding:10px; text-align:left;'>Motif Type</th>
            <th style='padding:10px; text-align:left;'>Arm/Unit Length</th>
            <th style='padding:10px; text-align:left;'>Spacer/Loop</th>
            <th style='padding:10px; text-align:left;'>Additional Requirements</th>
        </tr>
        <tr style='background:#f5f5f5;'>
            <td style='padding:10px;'><b>Cruciform DNA</b></td>
            <td style='padding:10px;'>10–100 nt</td>
            <td style='padding:10px;'>0–3 nt</td>
            <td style='padding:10px;'>Reverse complement arms</td>
        </tr>
        <tr>
            <td style='padding:10px;'><b>Slipped DNA</b></td>
            <td style='padding:10px;'>10–50 nt repeat</td>
            <td style='padding:10px;'>0 nt (adjacent)</td>
            <td style='padding:10px;'>Direct repeat units</td>
        </tr>
        <tr style='background:#f5f5f5;'>
            <td style='padding:10px;'><b>Triplex DNA</b></td>
            <td style='padding:10px;'>10–100 nt mirrored</td>
            <td style='padding:10px;'>0–8 nt</td>
            <td style='padding:10px;'>≥90% Purine or Pyrimidine</td>
        </tr>
    </table>
    </div>
    """, unsafe_allow_html=True)
    
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
        <li><b>A-philic DNA</b>: Uses tetranucleotide log2 odds scoring to identify A-tract-favoring sequences with high protein-binding affinity. Classified as high-confidence or moderate A-philic based on score thresholds and strong tetranucleotide counts.</li>
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

st.markdown("""
---
<div style='font-size: 1.05rem; color: #1e293b; margin-top: 30px; text-align: left; font-family:Montserrat,Arial;'>
<b>Developed by</b><br>
Dr. Venkata Rajesh Yella<br>
<a href='mailto:yvrajesh_bt@kluniversity.in'>yvrajesh_bt@kluniversity.in</a> |
<a href='https://github.com/VRYella' target='_blank'>GitHub: VRYella</a>
</div>
""", unsafe_allow_html=True)
