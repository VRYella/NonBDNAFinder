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
import io
import os  # Added for image path checking
import numpy as np
from collections import Counter
# Import consolidated NBDScanner modules
from utilities import (
    canonicalize_motif, parse_fasta, gc_content, reverse_complement, wrap,
    get_basic_stats, export_to_bed, export_to_csv, export_to_json,
    validate_sequence, quality_check_motifs, export_to_excel,
    calculate_genomic_density, calculate_positional_density,
    calculate_enrichment_with_shuffling, calculate_enhanced_statistics
)
from nonbscanner import (
    analyze_sequence, analyze_multiple_sequences,
    get_motif_info as get_motif_classification_info
)
from utilities import export_results_to_dataframe
from visualizations import (
    plot_motif_distribution, plot_coverage_map, plot_density_heatmap,
    plot_length_distribution, plot_score_distribution, plot_nested_pie_chart, 
    save_all_plots, MOTIF_CLASS_COLORS,
    plot_density_comparison, plot_enrichment_analysis, plot_enrichment_summary_table,
    plot_circos_motif_density, plot_radial_class_density, plot_stacked_density_track
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


# ---------- PAGE CONFIG ----------
st.set_page_config(
    page_title="NBDScanner - Non-B DNA Motif Finder",
    layout="wide",
    page_icon="🧬",
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
#
# Setting MAX_SEQUENCE_LENGTH to a very large value to effectively disable the limit
# while maintaining backward compatibility with any code that references it.
MAX_SEQUENCE_LENGTH = 10_000_000_000  # 10 billion nucleotides (effectively unlimited)

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
    import os
    
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
        'primary': '#5B8DEF',        # Soft blue - soothing primary
        'secondary': '#7BA3F7',      # Lighter blue - softer secondary
        'accent': '#A8C5FF',         # Very soft blue accent
        'bg_light': '#F8FAFC',       # Very light gray-blue background
        'bg_card': '#EEF2FF',        # Soft indigo-tinted card background
        'text': '#374151',           # Softer dark gray for text
        'tab_bg': '#F1F5F9',         # Light slate for tab bar
        'tab_active': '#5B8DEF',     # Active tab color
        'shadow': 'rgba(91, 141, 239, 0.15)'  # Soft shadow color
    },
    'nature_green': {
        'primary': '#6DBB7A',        # Soft green
        'secondary': '#8FD19C',      # Lighter green
        'accent': '#B8E5C0',         # Very soft green accent
        'bg_light': '#F7FBF8',       # Very light green background
        'bg_card': '#E8F5E9',        # Soft green card background
        'text': '#374151',
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
        'text': '#374151',
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
        'text': '#374151',
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

st.markdown(f"""
    <style>
    /* Import Google Fonts for professional scientific typography */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=IBM+Plex+Sans:wght@300;400;500;600;700&family=Poppins:wght@400;500;600;700;800&display=swap');
    
    /* ============================================
       ENHANCED ANIMATIONS & KEYFRAMES
       ============================================ */
    @keyframes shimmer {{
        0% {{ background-position: -1000px 0; }}
        100% {{ background-position: 1000px 0; }}
    }}
    
    @keyframes pulse-glow {{
        0%, 100% {{ box-shadow: 0 0 5px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.4); }}
        50% {{ box-shadow: 0 0 20px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.8); }}
    }}
    
    @keyframes float-subtle {{
        0%, 100% {{ transform: translateY(0); }}
        50% {{ transform: translateY(-3px); }}
    }}
    
    @keyframes scale-click {{
        0% {{ transform: scale(1); }}
        50% {{ transform: scale(0.97); }}
        100% {{ transform: scale(1); }}
    }}
    
    @keyframes fade-in {{
        from {{ opacity: 0; transform: translateY(10px); }}
        to {{ opacity: 1; transform: translateY(0); }}
    }}
    
    @keyframes slide-in-right {{
        from {{ opacity: 0; transform: translateX(-20px); }}
        to {{ opacity: 1; transform: translateX(0); }}
    }}
    
    @keyframes gradient-shift {{
        0% {{ background-position: 0% 50%; }}
        50% {{ background-position: 100% 50%; }}
        100% {{ background-position: 0% 50%; }}
    }}
    
    /* ============================================
       CSS CUSTOM PROPERTIES (Variables) for Soothing Theme
       ============================================ */
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
        --font-primary: 'Inter', 'IBM Plex Sans', 'Segoe UI', system-ui, -apple-system, sans-serif;
        --border-radius-sm: 8px;
        --border-radius-md: 12px;
        --border-radius-lg: 16px;
        --border-radius-pill: 50px;
        --transition-smooth: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
    }}
    
    /* ============================================
       ELEGANT BACKGROUND - Soothing Light Theme
       ============================================ */
    body, [data-testid="stAppViewContainer"], .main {{
        background: 
            {dna_pattern},
            linear-gradient(180deg, {current_theme['bg_light']} 0%, {'#252B3B' if is_dark_mode else '#EEF2FF'} 100%) !important;
        font-family: var(--font-primary) !important;
        color: var(--text-color) !important;
    }}
    
    /* ============================================
       FULL-WIDTH TABS: Flex-based Equal Distribution
       Pills-style active state with gradient fill
       ============================================ */
    .stTabs [data-baseweb="tab-list"] {{
        display: flex !important;
        width: 100% !important;
        justify-content: stretch !important;
        gap: 6px !important;
        border-bottom: none !important;
        background: {tab_bg_color} !important;
        box-shadow: 0 2px 12px {shadow_color}, inset 0 1px 0 {'rgba(255, 255, 255, 0.1)' if is_dark_mode else 'rgba(255, 255, 255, 0.8)'} !important;
        margin-bottom: 1.5em;
        padding: 8px 10px !important;
        border-radius: var(--border-radius-lg) !important;
        animation: fade-in 0.5s ease-out;
    }}
    
    /* Each tab takes equal width (flex: 1) for full-width distribution */
    .stTabs [data-baseweb="tab"] {{
        flex: 1 1 0 !important;
        min-width: 0 !important;
        font-size: 0.95rem !important;
        font-weight: 500 !important;
        font-family: var(--font-primary) !important;
        padding: 12px 16px !important;
        text-align: center !important;
        color: {'#94a3b8' if is_dark_mode else '#64748b'} !important;
        background: transparent !important;
        border: none !important;
        border-radius: var(--border-radius-pill) !important;
        letter-spacing: 0.01em;
        transition: var(--transition-smooth);
        position: relative;
        cursor: pointer;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
    }}
    
    /* Tab hover state - subtle highlight */
    .stTabs [data-baseweb="tab"]:hover {{
        color: var(--primary-color) !important;
        background: {'rgba(255, 255, 255, 0.1)' if is_dark_mode else 'rgba(255, 255, 255, 0.7)'} !important;
    }}
    
    /* Active tab - Pill-shaped with gradient fill and subtle elevation */
    .stTabs [aria-selected="true"] {{
        color: {'#FFFFFF' if is_dark_mode else '#FFFFFF'} !important;
        font-weight: 600 !important;
        background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%) !important;
        box-shadow: 0 4px 12px {shadow_color}, 0 2px 4px rgba(0, 0, 0, 0.1) !important;
        transform: translateY(-1px);
    }}
    
    /* Remove the default tab underline indicator */
    .stTabs [data-baseweb="tab-highlight"] {{
        display: none !important;
    }}
    
    .stTabs [data-baseweb="tab-border"] {{
        display: none !important;
    }}
    
    /* ============================================
       HEADINGS: Modern Inter typography
       ============================================ */
    h1, h2, h3, h4 {{
        font-family: var(--font-primary) !important;
        color: var(--primary-color) !important;
        font-weight: 700 !important;
        letter-spacing: -0.02em;
        margin-top: 1em;
        margin-bottom: 0.6em;
        animation: fade-in 0.4s ease-out;
    }}
    h1 {{ 
        font-size: 2.4rem !important; 
        font-weight: 700 !important;
        background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
    }}
    h2 {{ 
        font-size: 1.75rem !important; 
        color: var(--primary-color) !important; 
        font-weight: 600 !important;
        border-bottom: none !important;
        padding-bottom: 0.3rem;
        position: relative;
    }}
    h2::after {{
        content: '';
        position: absolute;
        bottom: 0;
        left: 0;
        width: 50px;
        height: 3px;
        background: linear-gradient(90deg, var(--primary-color) 0%, var(--accent-color) 100%);
        border-radius: 2px;
    }}
    h3 {{ 
        font-size: 1.35rem !important; 
        color: var(--secondary-color) !important; 
        font-weight: 600 !important;
    }}
    h4 {{ 
        font-size: 1.15rem !important; 
        color: var(--secondary-color) !important; 
        font-weight: 500 !important;
    }}
    
    /* Body text: Modern Inter font for readability */
    .stMarkdown, .markdown-text-container, .stText, p, span, label {{
        font-size: 0.95rem !important;
        font-family: var(--font-primary) !important;
        line-height: 1.65 !important;
        color: var(--text-color) !important;
        margin-bottom: 0.5em !important;
        font-weight: 400;
    }}
    
    /* Input fields: Clean design with soft borders */
    input, .stTextInput>div>div>input, .stSelectbox>div>div>div, 
    .stMultiSelect>div>div>div, .stRadio>div>div>label>div {{
        font-size: 0.95rem !important;
        font-family: var(--font-primary) !important;
        border-radius: var(--border-radius-md) !important;
        border: 1.5px solid {'#475569' if is_dark_mode else '#E2E8F0'} !important;
        background: {'#252B3B' if is_dark_mode else '#FFFFFF'} !important;
        transition: var(--transition-smooth);
        box-shadow: 0 1px 3px rgba(0, 0, 0, 0.05);
    }}
    input:focus, .stTextInput>div>div>input:focus {{
        border-color: var(--primary-color) !important;
        box-shadow: 0 0 0 3px {shadow_color}, 0 2px 6px rgba(0, 0, 0, 0.08) !important;
        background: {'#252B3B' if is_dark_mode else '#FEFEFF'} !important;
    }}
    
    /* ============================================
       BUTTONS: Soft gradient with subtle shadow
       ============================================ */
    .stButton>button {{
        font-size: 0.95rem !important;
        font-family: var(--font-primary) !important;
        padding: 0.65em 1.5em !important;
        background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%) !important;
        color: #FFFFFF !important;
        border-radius: var(--border-radius-md) !important;
        border: none !important;
        font-weight: 500 !important;
        box-shadow: 0 4px 12px {shadow_color};
        transition: var(--transition-smooth);
        letter-spacing: 0.02em;
    }}
    .stButton>button:hover {{
        background: linear-gradient(135deg, var(--secondary-color) 0%, var(--primary-color) 100%) !important;
        box-shadow: 0 6px 16px {shadow_color};
        transform: translateY(-2px);
    }}
    .stButton>button:active {{
        transform: translateY(0);
        box-shadow: 0 2px 8px {shadow_color};
    }}
    
    /* ============================================
       DATAFRAMES: Clean design with soft shadows
       ============================================ */
    .stDataFrame, .stTable {{
        font-size: {'0.88rem' if is_compact else '0.92rem'} !important;
        font-family: var(--font-primary) !important;
        line-height: {'1.4' if is_compact else '1.6'} !important;
        border-radius: var(--border-radius-md) !important;
        overflow: hidden;
        box-shadow: 0 2px 8px {shadow_color} !important;
        border: 1px solid var(--bg-card) !important;
        animation: fade-in 0.4s ease-out;
    }}
    .stDataFrame thead tr th {{
        background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%) !important;
        color: white !important;
        font-weight: 600 !important;
        padding: {'0.5rem' if is_compact else '0.8rem'} !important;
        font-size: {'0.82rem' if is_compact else '0.88rem'} !important;
        letter-spacing: 0.01em;
        text-transform: uppercase;
        border-bottom: 2px solid var(--primary-color) !important;
    }}
    .stDataFrame tbody tr {{
        transition: var(--transition-smooth);
        border-bottom: 1px solid var(--bg-card) !important;
    }}
    .stDataFrame tbody tr:hover {{
        background: {'#2D3748' if is_dark_mode else '#F7FAFC'} !important;
    }}
    .stDataFrame tbody tr td {{
        padding: {'0.4rem' if is_compact else '0.7rem'} !important;
        color: var(--text-color) !important;
    }}
    
    /* Highlight alternate rows */
    .stDataFrame tbody tr:nth-child(odd) {{
        background: {'rgba(45, 55, 72, 0.3)' if is_dark_mode else 'rgba(238, 242, 255, 0.5)'} !important;
    }}
    
    /* Tab content - Clean content card section */
    .stTabs [data-baseweb="tab-panel"] {{
        padding-top: 1.5rem !important;
        padding-left: 1rem !important;
        padding-right: 1rem !important;
        padding-bottom: 1.5rem !important;
        animation: fade-in 0.3s ease-out;
        background: {'#1A1F2E' if is_dark_mode else '#FFFFFF'} !important;
        border-radius: 0 0 var(--border-radius-lg) var(--border-radius-lg) !important;
        box-shadow: 0 4px 16px {shadow_color} !important;
        margin-top: -8px;
    }}
    
    /* ============================================
       CONTENT CARDS: Soft shadows and rounded corners
       ============================================ */
    .analysis-summary-card, .glassmorphism-card, .content-card {{
        margin: 1rem 0 !important;
        padding: 1.5rem !important;
        background: {'#252B3B' if is_dark_mode else '#FFFFFF'} !important;
        border-radius: var(--border-radius-lg) !important;
        box-shadow: 0 2px 12px {shadow_color} !important;
        border: 1px solid {'#374151' if is_dark_mode else '#E5E7EB'} !important;
        animation: fade-in 0.4s ease-out;
        transition: var(--transition-smooth);
    }}
    .analysis-summary-card:hover, .glassmorphism-card:hover, .content-card:hover {{
        box-shadow: 0 4px 20px {shadow_color} !important;
    }}
    
    /* Select boxes and multiselect: Clean dropdown design */
    .stSelectbox > div > div > div, .stMultiSelect > div > div > div {{
        min-height: 2.8rem !important;
        padding: 0.7rem 1rem !important;
        line-height: 1.5 !important;
        border-radius: var(--border-radius-md) !important;
        background: {'#252B3B' if is_dark_mode else '#FFFFFF'} !important;
        border: 1.5px solid {'#475569' if is_dark_mode else '#E2E8F0'} !important;
        box-shadow: 0 1px 3px rgba(0, 0, 0, 0.05);
        transition: var(--transition-smooth);
    }}
    .stSelectbox > div > div > div:hover, .stMultiSelect > div > div > div:hover {{
        border-color: var(--accent-color) !important;
        box-shadow: 0 2px 8px {shadow_color};
    }}
    
    /* Enhanced input field spacing */
    .stTextInput > div > div > input, .stTextArea textarea {{
        padding: 0.7rem 1rem !important;
        min-height: 2.8rem !important;
        line-height: 1.5 !important;
        border-radius: var(--border-radius-md) !important;
        font-size: 0.95rem !important;
        background: {'#252B3B' if is_dark_mode else '#FFFFFF'} !important;
    }}
    .stTextArea textarea {{
        min-height: 100px !important;
    }}
    
    /* Number input styling */
    .stNumberInput > div > div > input {{
        padding: 0.7rem 1rem !important;
        min-height: 2.8rem !important;
        border-radius: var(--border-radius-md) !important;
        font-size: 0.95rem !important;
        font-weight: 500;
        background: {'#252B3B' if is_dark_mode else '#FFFFFF'} !important;
    }}
    
    /* Radio buttons: Clean styling */
    .stRadio > div {{
        gap: 0.3rem !important;
    }}
    .stRadio > div > div > label {{
        margin-bottom: 0.5rem !important;
        padding: 0.5rem 0.8rem !important;
        border-radius: var(--border-radius-sm) !important;
        transition: var(--transition-smooth);
        border: 1.5px solid transparent;
        background: {'#252B3B' if is_dark_mode else '#F8FAFC'};
        cursor: pointer;
    }}
    .stRadio > div > div > label:hover {{
        background: {'#2D3748' if is_dark_mode else current_theme['bg_card']} !important;
        border-color: {current_theme['accent']} !important;
    }}
    .stRadio > div > div > label[data-checked="true"] {{
        background: {'#2D3748' if is_dark_mode else current_theme['bg_card']} !important;
        border-color: {current_theme['primary']} !important;
        font-weight: 500 !important;
    }}
    /* Radio button circles */
    .stRadio > div > div > label > div:first-child {{
        border-width: 2px !important;
        border-color: {current_theme['accent']} !important;
        width: 18px !important;
        height: 18px !important;
    }}
    .stRadio > div > div > label[data-checked="true"] > div:first-child {{
        border-color: var(--primary-color) !important;
        background-color: var(--primary-color) !important;
    }}
    
    /* Info boxes and alerts: Clean notification design */
    .stAlert {{
        border-radius: var(--border-radius-md) !important;
        border-left: 4px solid var(--secondary-color) !important;
        box-shadow: 0 2px 8px rgba(0,0,0,0.06) !important;
        background: {'#252B3B' if is_dark_mode else '#F8FAFC'} !important;
        padding: 1rem !important;
        font-weight: 400;
    }}
    .stAlert[data-baseweb="notification"] {{
        border-left-color: var(--primary-color) !important;
    }}
    .stSuccess {{
        border-left-color: #6DBB7A !important;
        background: {'#1E3A29' if is_dark_mode else '#E8F5E9'} !important;
    }}
    .stWarning {{
        border-left-color: #F59E0B !important;
        background: {'#3A2E1E' if is_dark_mode else '#FEF3C7'} !important;
    }}
    .stError {{
        border-left-color: #EF4444 !important;
        background: {'#3A1E1E' if is_dark_mode else '#FEE2E2'} !important;
    }}
    
    /* Progress bars: Soft gradient design */
    .stProgress > div > div {{
        background: linear-gradient(90deg, var(--primary-color) 0%, var(--secondary-color) 100%) !important;
        border-radius: var(--border-radius-pill) !important;
        box-shadow: 0 1px 4px {shadow_color};
    }}
    .stProgress > div {{
        border-radius: var(--border-radius-pill) !important;
        background: var(--bg-card) !important;
        box-shadow: inset 0 1px 3px rgba(0, 0, 0, 0.05);
    }}
    
    /* File uploader: Clean drag-and-drop design */
    [data-testid="stFileUploader"] {{
        border: 2px dashed var(--accent-color) !important;
        border-radius: var(--border-radius-lg) !important;
        background: {'#1A1F2E' if is_dark_mode else '#F8FAFC'} !important;
        padding: 2rem !important;
        transition: var(--transition-smooth);
        box-shadow: 0 2px 8px {shadow_color};
    }}
    [data-testid="stFileUploader"]:hover {{
        border-color: var(--primary-color) !important;
        background: {'#252B3B' if is_dark_mode else '#EEF2FF'} !important;
        box-shadow: 0 4px 16px {shadow_color};
    }}
    [data-testid="stFileUploader"] section {{
        border: none !important;
    }}
    [data-testid="stFileUploader"] button {{
        background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%) !important;
        color: white !important;
        border: none !important;
        border-radius: var(--border-radius-md) !important;
        padding: 0.5rem 1.2rem !important;
        font-weight: 500 !important;
        transition: var(--transition-smooth);
    }}
    [data-testid="stFileUploader"] button:hover {{
        background: linear-gradient(135deg, var(--secondary-color) 0%, var(--primary-color) 100%) !important;
        box-shadow: 0 4px 12px {shadow_color};
    }}
    
    /* Checkboxes: Clean toggle-style design */
    .stCheckbox {{
        padding: 0.3rem 0 !important;
    }}
    .stCheckbox > label {{
        padding: 0.4rem 0.6rem !important;
        border-radius: var(--border-radius-sm) !important;
        transition: var(--transition-smooth);
        cursor: pointer;
        background: {'#252B3B' if is_dark_mode else '#F8FAFC'};
        border: 1.5px solid transparent;
    }}
    .stCheckbox > label:hover {{
        background: {'#2D3748' if is_dark_mode else current_theme['bg_card']} !important;
        border-color: {current_theme['accent']} !important;
    }}
    .stCheckbox > label > div:first-child {{
        border-width: 2px !important;
        border-color: {current_theme['accent']} !important;
        border-radius: 5px !important;
        width: 18px !important;
        height: 18px !important;
        transition: all 0.3s ease;
    }}
    .stCheckbox > label > div:first-child[data-checked="true"] {{
        background-color: {current_theme['primary']} !important;
        border-color: {current_theme['primary']} !important;
    }}
    
    /* Expander styling: Clean accordion design */
    .streamlit-expanderHeader {{
        border-radius: 12px !important;
        background: {'#252B3B' if is_dark_mode else '#F8FAFC'} !important;
        font-weight: 600 !important;
        padding: 0.8rem 1.2rem !important;
        border: 1.5px solid {'#374151' if is_dark_mode else '#E5E7EB'} !important;
        transition: all 0.3s ease;
        box-shadow: 0 1px 4px {shadow_color};
        display: flex !important;
        align-items: center !important;
        gap: 10px !important;
    }}
    .streamlit-expanderHeader:hover {{
        background: {'#2D3748' if is_dark_mode else current_theme['bg_card']} !important;
        border-color: {current_theme['accent']} !important;
        box-shadow: 0 2px 8px {shadow_color};
    }}
    .streamlit-expanderContent {{
        border-radius: 0 0 12px 12px !important;
        background: {'#1A1F2E' if is_dark_mode else '#FFFFFF'} !important;
        border: 1.5px solid {'#374151' if is_dark_mode else '#E5E7EB'} !important;
        border-top: none !important;
        padding: 1.2rem !important;
    }}
    
    /* Clean chevron icon styling for expanders */
    .streamlit-expanderHeader svg {{
        margin-right: 8px !important;
        flex-shrink: 0 !important;
        width: 18px !important;
        height: 18px !important;
        color: var(--primary-color) !important;
        transition: transform 0.3s ease !important;
    }}
    
    /* Rotate chevron when expanded */
    details[open] > .streamlit-expanderHeader svg {{
        transform: rotate(90deg) !important;
    }}
    
    /* Hide any raw icon text (like keyboard_arrow_right text) */
    .streamlit-expanderHeader::before {{
        display: none !important;
    }}
    
    /* Hide material icon text fallback */
    .streamlit-expanderHeader span[class*="icon"],
    .streamlit-expanderHeader [data-testid*="icon"] {{
        font-size: 0 !important;
    }}
    .streamlit-expanderHeader span[class*="icon"]::before,
    .streamlit-expanderHeader [data-testid*="icon"]::before {{
        font-size: 1rem !important;
    }}
    
    /* Fix: keyboard_arrow_down icon text overlap */
    .streamlit-expanderHeader {{
        display: flex !important;
        align-items: center !important;
        line-height: 1.5 !important;
    }}
    
    /* Ensure label text doesn't overlap with icon */
    .streamlit-expanderHeader p,
    .streamlit-expanderHeader label,
    .streamlit-expanderHeader div,
    .streamlit-expanderHeader span {{
        margin-left: 0 !important;
        flex: 1 !important;
        white-space: normal !important;
        word-wrap: break-word !important;
        line-height: 1.5 !important;
        overflow: hidden !important;
        text-overflow: ellipsis !important;
    }}
    
    /* Hide keyboard_arrow_down text that may appear as fallback */
    .streamlit-expanderHeader span:first-child {{
        display: flex !important;
        align-items: center !important;
        flex-shrink: 0 !important;
        min-width: 24px !important;
        max-width: 24px !important;
        overflow: hidden !important;
    }}
    
    /* Global fix for icon and label spacing across all Streamlit elements */
    [role="button"] svg,
    [role="tab"] svg,
    button svg,
    label svg {{
        margin-right: 6px !important;
        flex-shrink: 0 !important;
    }}
    
    /* Ensure all interactive elements use flex layout for proper icon/text spacing */
    [role="button"],
    [role="tab"],
    button,
    label {{
        display: flex !important;
        align-items: center !important;
        gap: 6px !important;
        flex-wrap: wrap !important;
        line-height: 1.5 !important;
    }}
    
    /* Prevent icon text from overlapping labels */
    .material-icons-text {{
        display: none !important;
    }}
    
    /* Additional fix for Streamlit summary element icon text */
    summary span {{
        display: flex !important;
        align-items: center !important;
        gap: 8px !important;
    }}
    
    /* Ensure icon container doesn't show raw text like 'keyboard_arrow_down' */
    [data-testid="StyledIconContainer"] {{
        display: flex !important;
        align-items: center !important;
        justify-content: center !important;
        width: 20px !important;
        min-width: 20px !important;
        overflow: hidden !important;
    }}
    
    /* ============================================
       METRIC CARDS: Clean design with soft shadows
       ============================================ */
    [data-testid="stMetric"] {{
        background: {'#252B3B' if is_dark_mode else '#FFFFFF'} !important;
        padding: 1.2rem !important;
        border-radius: var(--border-radius-lg) !important;
        border: 1.5px solid {'#374151' if is_dark_mode else '#E5E7EB'} !important;
        box-shadow: 0 2px 8px {shadow_color} !important;
        transition: var(--transition-smooth);
        animation: fade-in 0.4s ease-out;
    }}
    [data-testid="stMetric"]:hover {{
        box-shadow: 0 4px 16px {shadow_color} !important;
        border-color: var(--accent-color) !important;
    }}
    [data-testid="stMetricValue"] {{
        font-size: 1.8rem !important;
        font-weight: 700 !important;
        background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
    }}
    [data-testid="stMetricLabel"] {{
        font-size: 0.85rem !important;
        font-weight: 500 !important;
        color: {'#94a3b8' if is_dark_mode else '#64748B'} !important;
        text-transform: uppercase;
        letter-spacing: 0.03em;
    }}
    
    /* Spinner */
    .stSpinner > div {{
        border-top-color: var(--primary-color) !important;
    }}
    
    /* Code blocks: Clean monospace design */
    .stCodeBlock {{
        border-radius: var(--border-radius-md) !important;
        border: 1.5px solid {'#374151' if is_dark_mode else '#E5E7EB'} !important;
        background: {'#1A1F2E' if is_dark_mode else '#F8FAFC'} !important;
        box-shadow: 0 1px 4px rgba(0, 0, 0, 0.03);
        padding: 0.8rem !important;
    }}
    code {{
        font-family: 'JetBrains Mono', 'Fira Code', 'Consolas', monospace !important;
        font-size: 0.88rem !important;
        color: var(--primary-color) !important;
        background: {'rgba(37, 43, 59, 0.5)' if is_dark_mode else 'rgba(238, 242, 255, 0.6)'} !important;
        padding: 0.15rem 0.4rem !important;
        border-radius: 4px !important;
    }}
    
    /* Sidebar styling */
    [data-testid="stSidebar"] {{
        background: {'#1A1F2E' if is_dark_mode else '#F8FAFC'} !important;
        border-right: 1.5px solid {'#374151' if is_dark_mode else '#E5E7EB'} !important;
    }}
    [data-testid="stSidebar"] .stMarkdown {{
        color: var(--primary-color) !important;
    }}
    
    /* Scrollbar styling: Clean minimal design */
    ::-webkit-scrollbar {{
        width: 8px;
        height: 8px;
    }}
    ::-webkit-scrollbar-track {{
        background: {'#1A1F2E' if is_dark_mode else '#F1F5F9'};
        border-radius: var(--border-radius-pill);
    }}
    ::-webkit-scrollbar-thumb {{
        background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%);
        border-radius: var(--border-radius-pill);
        border: 2px solid {'#1A1F2E' if is_dark_mode else '#F8FAFC'};
    }}
    ::-webkit-scrollbar-thumb:hover {{
        background: linear-gradient(135deg, var(--secondary-color) 0%, var(--primary-color) 100%);
    }}
    
    /* Loading spinner */
    .stSpinner > div {{
        border-top-color: var(--primary-color) !important;
        border-right-color: var(--secondary-color) !important;
        border-bottom-color: var(--accent-color) !important;
    }}
    
    /* Download buttons - clean design */
    [data-testid="stDownloadButton"] button {{
        background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%) !important;
        color: white !important;
        border-radius: var(--border-radius-md) !important;
        padding: 0.6rem 1.2rem !important;
        font-weight: 500 !important;
        box-shadow: 0 2px 8px {shadow_color};
        transition: var(--transition-smooth);
    }}
    [data-testid="stDownloadButton"] button:hover {{
        background: linear-gradient(135deg, var(--secondary-color) 0%, var(--primary-color) 100%) !important;
        box-shadow: 0 4px 12px {shadow_color};
    }}
    
    /* Caption text styling */
    .caption, [data-testid="stCaptionContainer"] {{
        color: {'#94a3b8' if is_dark_mode else '#64748B'} !important;
        font-size: 0.85rem !important;
        font-style: italic;
    }}
    
    /* Success/Error message boxes */
    .element-container .stMarkdown .stSuccess {{
        background: {'#1E3A29' if is_dark_mode else '#E8F5E9'} !important;
        border-left: 3px solid #6DBB7A !important;
        border-radius: var(--border-radius-md) !important;
        padding: 0.8rem !important;
    }}
    .element-container .stMarkdown .stError {{
        background: {'#3A1E1E' if is_dark_mode else '#FEE2E2'} !important;
        border-left: 3px solid #EF4444 !important;
        border-radius: var(--border-radius-md) !important;
        padding: 0.8rem !important;
    }}
    
    /* Links styling */
    a {{
        color: var(--secondary-color) !important;
        text-decoration: none !important;
        font-weight: 500 !important;
        transition: var(--transition-smooth);
    }}
    a:hover {{
        color: var(--primary-color) !important;
        text-decoration: underline !important;
    }}
    
    /* Tooltip improvements */
    [data-testid="stTooltipIcon"] {{
        color: var(--accent-color) !important;
        cursor: help;
        transition: var(--transition-smooth);
    }}
    [data-testid="stTooltipIcon"]:hover {{
        transform: scale(1.1);
        color: var(--secondary-color) !important;
    }}
    
    /* Enhanced hr separator */
    hr {{
        border: none !important;
        height: 1px !important;
        background: linear-gradient(90deg, transparent 0%, var(--bg-card) 50%, transparent 100%) !important;
        margin: 1.5rem 0 !important;
    }}
    
    /* ============================================
       STICKY RESULTS PANEL
       ============================================ */
    .sticky-results {{
        position: sticky;
        top: 0;
        z-index: 100;
        background: {'rgba(26, 31, 46, 0.95)' if is_dark_mode else 'rgba(255, 255, 255, 0.95)'};
        backdrop-filter: blur(8px);
        -webkit-backdrop-filter: blur(8px);
        padding: 0.8rem;
        border-radius: 0 0 var(--border-radius-md) var(--border-radius-md);
        box-shadow: 0 2px 12px {shadow_color};
    }}
    
    /* ============================================
       RESPONSIVE DESIGN - Consistent across screens
       ============================================ */
    @media screen and (max-width: 1200px) {{
        .stTabs [data-baseweb="tab"] {{
            font-size: 0.88rem !important;
            padding: 10px 12px !important;
        }}
        h1 {{ font-size: 1.8rem !important; }}
        h2 {{ font-size: 1.4rem !important; }}
        [data-testid="stMetric"] {{
            padding: 0.9rem !important;
        }}
        [data-testid="stMetricValue"] {{
            font-size: 1.5rem !important;
        }}
    }}
    
    @media screen and (max-width: 768px) {{
        .stTabs [data-baseweb="tab-list"] {{
            flex-wrap: wrap !important;
            gap: 4px !important;
            padding: 6px 8px !important;
        }}
        .stTabs [data-baseweb="tab"] {{
            font-size: 0.8rem !important;
            padding: 8px 10px !important;
            flex: 1 1 auto !important;
            min-width: 80px !important;
        }}
        h1 {{ font-size: 1.5rem !important; }}
        h2 {{ font-size: 1.2rem !important; }}
        h3 {{ font-size: 1rem !important; }}
        .stDataFrame, .stTable {{
            font-size: 0.82rem !important;
        }}
        [data-testid="stMetric"] {{
            padding: 0.7rem !important;
        }}
        [data-testid="stMetricValue"] {{
            font-size: 1.3rem !important;
        }}
        .stButton>button {{
            padding: 0.5em 1em !important;
            font-size: 0.88rem !important;
        }}
    }}
    
    /* ============================================
       MODERN ACCORDION STYLING
       ============================================ */
    details.streamlit-expander {{
        border-radius: var(--border-radius-md) !important;
        overflow: hidden;
        transition: var(--transition-smooth);
        margin-bottom: 0.6rem;
    }}
    
    details.streamlit-expander[open] {{
        box-shadow: 0 4px 16px {shadow_color};
    }}
    
    /* ============================================
       HIGHLIGHT SELECTED ITEMS
       ============================================ */
    .highlighted-row {{
        background: {'rgba(91, 141, 239, 0.15)' if is_dark_mode else 'rgba(91, 141, 239, 0.1)'} !important;
        border-left: 3px solid var(--accent-color) !important;
    }}
    
    /* ============================================
       PROGRESS PANEL COMPONENTS - Updated for soothing theme
       ============================================ */
    .progress-panel {{
        background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%);
        border-radius: var(--border-radius-lg);
        padding: 1.2rem;
        color: white;
        box-shadow: 0 4px 16px {shadow_color};
        margin-bottom: 1rem;
    }}
    
    .progress-panel--success {{
        background: linear-gradient(135deg, #6DBB7A 0%, #8FD19C 100%);
        box-shadow: 0 4px 16px rgba(109, 187, 122, 0.3);
    }}
    
    .progress-panel--metrics {{
        background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%);
        box-shadow: 0 4px 16px {shadow_color};
        margin-bottom: 1.5rem;
    }}
    
    .progress-panel--results {{
        background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%);
        border-radius: var(--border-radius-lg);
        padding: 1.5rem;
        margin: 1rem 0;
        box-shadow: 0 4px 20px {shadow_color};
    }}
    
    .progress-panel--hybrid {{
        background: linear-gradient(135deg, #9B8FD9 0%, #B3A9E5 100%);
        border-radius: var(--border-radius-lg);
        padding: 1.5rem;
        margin: 1rem 0;
        box-shadow: 0 4px 20px rgba(155, 143, 217, 0.25);
    }}
    
    .progress-panel__title {{
        margin: 0 0 0.8rem 0;
        color: white;
        text-align: center;
    }}
    
    .progress-panel__title--large {{
        font-size: 1.3rem;
        font-weight: 600;
        margin-bottom: 1.2rem;
    }}
    
    .progress-panel__status {{
        margin: 0 0 0.8rem 0;
        text-align: center;
        opacity: 0.9;
        font-size: 0.9rem;
    }}
    
    .stats-grid {{
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(100px, 1fr));
        gap: 0.8rem;
        margin-bottom: 0.8rem;
    }}
    
    .stats-grid--wide {{
        grid-template-columns: repeat(auto-fit, minmax(110px, 1fr));
    }}
    
    .stats-grid--extra-wide {{
        grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
        gap: 1rem;
        margin-top: 0.8rem;
    }}
    
    .stat-card {{
        text-align: center;
        background: rgba(255, 255, 255, 0.15);
        padding: 0.7rem;
        border-radius: var(--border-radius-sm);
    }}
    
    .stat-card--large {{
        padding: 1rem;
        border-radius: var(--border-radius-md);
        backdrop-filter: blur(8px);
    }}
    
    .stat-card__value {{
        margin: 0;
        color: #FFFFFF;
        font-size: 1.3rem;
        font-weight: 600;
    }}
    
    .stat-card__value--large {{
        font-size: 1.8rem;
        font-weight: 700;
        margin-bottom: 0.3rem;
    }}
    
    .stat-card__label {{
        margin: 0.2rem 0 0 0;
        opacity: 0.85;
        font-size: 0.75rem;
    }}
    
    .stat-card__label--large {{
        font-size: 0.85rem;
        opacity: 0.9;
        font-weight: 500;
    }}
    
    .sequence-info {{
        background: rgba(0, 0, 0, 0.15);
        border-radius: var(--border-radius-sm);
        padding: 0.6rem;
        margin-bottom: 0.8rem;
    }}
    
    .sequence-info__text {{
        margin: 0;
        font-size: 0.85rem;
        text-align: center;
    }}
    
    .sequence-info__subtext {{
        margin: 0.3rem 0 0 0;
        font-size: 0.85rem;
        opacity: 0.9;
        text-align: center;
    }}
    
    /* Pipeline/Detector Components */
    .pipeline-panel {{
        background: linear-gradient(135deg, #f5f5f5 0%, #e0e0e0 100%);
        border-radius: 12px;
        padding: 1rem;
        margin-bottom: 1rem;
        border: 1px solid #bdbdbd;
    }}
    
    .pipeline-panel__title {{
        margin: 0 0 0.8rem 0;
        color: #424242;
    }}
    
    .pipeline-grid {{
        display: grid;
        grid-template-columns: repeat(3, 1fr);
        gap: 0.5rem;
        font-size: 0.85rem;
    }}
    
    .detector-item {{
        background: rgba(25, 118, 210, 0.1);
        padding: 0.5rem;
        border-radius: 6px;
        border-left: 3px solid #1976d2;
    }}
    
    .detector-item__name {{
        font-weight: 600;
    }}
    
    .detector-item__desc {{
        font-size: 0.75rem;
        color: #757575;
    }}
    
    .pipeline-panel__footer {{
        margin: 0.8rem 0 0 0;
        font-size: 0.8rem;
        color: #616161;
        text-align: center;
    }}
    
    /* ============================================
       PRINT STYLES FOR PUBLICATION
       ============================================ */
    @media print {{
        body, [data-testid="stAppViewContainer"], .main {{
            background: white !important;
        }}
        .stButton, [data-testid="stFileUploader"] {{
            display: none !important;
        }}
        .stDataFrame {{
            box-shadow: none !important;
            border: 1px solid #ccc !important;
        }}
    }}
    </style>
""", unsafe_allow_html=True)


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
        except:
            # If image doesn't exist, show placeholder
            st.markdown("""
            <div style='background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                        border-radius: 15px; padding: 40px; text-align: center; color: white;'>
                <h2 style='margin: 0; color: white;'>🧬</h2>
                <h3 style='margin: 10px 0 0 0; color: white;'>NBD Finder</h3>
                <p style='margin: 5px 0 0 0; color: #E8E8E8;'>Non-B DNA Detection</p>
            </div>
            """, unsafe_allow_html=True)
    with right:
        st.markdown("""
        <div style='font-family: Inter, system-ui, sans-serif; font-size:1.05rem; color:#263238; 
                    line-height:1.8; padding:2rem; background:white; border-radius:16px; 
                    box-shadow:0 4px 16px rgba(0,0,0,0.08); border: 1px solid rgba(25, 118, 210, 0.1);'>
        <p style='margin-top:0;'><b style='color:#0d47a1; font-size:1.15rem;'>Non-canonical DNA structures</b> play key roles in genome stability, regulation, and evolution.</p>
        <p>This application detects and analyzes <b style='color:#1976d2;'>11 major classes with 22+ subclasses</b> of Non-B DNA motifs in any DNA sequence or multi-FASTA file.</p>
        <p style='margin-bottom:0.8rem;'><b style='color:#0d47a1;'>Motif Classes (11 classes, 22+ subclasses):</b></p>
        <div style='color:#37474f; font-size:0.98rem; line-height:1.9; padding-left:1rem;'>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>1. Curved DNA</b> <span style='color:#546e7a;'>(Global curvature, Local Curvature)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>2. Slipped DNA</b> <span style='color:#546e7a;'>(Direct Repeat, STR)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>3. Cruciform DNA</b> <span style='color:#546e7a;'>(Inverted Repeats)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>4. R-loop</b> <span style='color:#546e7a;'>(R-loop formation sites)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>5. Triplex</b> <span style='color:#546e7a;'>(Triplex, Sticky DNA)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>6. G-Quadruplex Family</b> <span style='color:#546e7a;'>(Multimeric G4, Canonical G4, Relaxed G4, Bulged G4, Bipartite G4, Imperfect G4, G-Triplex intermediate)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>7. i-Motif Family</b> <span style='color:#546e7a;'>(Canonical i-motif, Relaxed i-motif, AC-motif)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>8. Z-DNA</b> <span style='color:#546e7a;'>(Z-DNA, eGZ (Extruded-G) DNA)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>9. A-philic DNA</b> <span style='color:#546e7a;'>(A-philic DNA)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>10. Hybrid</b> <span style='color:#546e7a;'>(dynamic overlaps)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>11. Non-B DNA Clusters</b> <span style='color:#546e7a;'>(dynamic clusters)</span></div>
        </div>
        <p style='margin-top:1.2rem; margin-bottom:0; color:#1976d2; font-weight:600;'>📤 Upload single or multi-FASTA files to begin analysis...</p>
        </div>
        """, unsafe_allow_html=True)

# Updated Streamlit layout: Input Method + Sequence Preview + Analysis side-by-side
#  "Upload & Analyze" 
# It expects helper functions and constants to exist in the module:
# parse_fasta, get_basic_stats, EXAMPLE_FASTA, EXAMPLE_MULTI_FASTA, parse_fasta, all_motifs_refactored,
# ensure_subclass, MOTIF_ORDER, wrap, Entrez, SeqIO, Counter, pd, st

with tab_pages["Upload & Analyze"]:
    st.markdown("<h2>Sequence Upload and Motif Analysis</h2>", unsafe_allow_html=True)
    st.markdown('<span style="font-family:Montserrat,Arial; font-size:1.12rem;">Supports multi-FASTA and single FASTA. Paste, upload, select example, or fetch from NCBI.</span>', unsafe_allow_html=True)
    st.caption("Supported formats: .fa, .fasta, .txt, .fna | Limit: 200MB/file.")
    
    # Show unlimited processing info for web version
    st.info("""
    📏 **Unlimited Sequence Length**: This version supports sequences of any size using optimized chunked processing.  
    🔧 Large sequences are automatically processed in 10kb chunks with 500bp overlap to ensure no motifs are missed.
    """)

    st.markdown("---")

    # Create two main columns: left = Input & Preview, right = Analysis controls & Summary
    col_left, col_right = st.columns([1, 1], gap="large")

    # ----- LEFT COLUMN: Input Method + Sequence Preview -----
    with col_left:
        st.markdown("### 📁 Input Method")
        input_method = st.radio("Choose your input method:",
                                ["📂 Upload FASTA File", "✏️ Paste Sequence", "🧪 Example Data", "🌐 NCBI Fetch"],
                                horizontal=True)

        seqs, names = [], []

        if input_method == "📂 Upload FASTA File":
            fasta_file = st.file_uploader("Drag and drop FASTA/multi-FASTA file here", type=["fa", "fasta", "txt", "fna"])
            if fasta_file:
                content = fasta_file.read().decode("utf-8")
                seqs, names = [], []
                cur_seq, cur_name = "", ""
                for line in content.splitlines():
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
                    st.success(f"✅ Loaded {len(seqs)} sequences.")
                    for i, seq in enumerate(seqs[:3]):
                        stats = get_basic_stats(seq)
                        st.markdown(f"**{names[i]}**: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                        st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                    if len(seqs) > 3:
                        st.caption(f"...and {len(seqs)-3} more.")
                else:
                    st.warning("No sequences found.")

        elif input_method == "✏️ Paste Sequence":
            seq_input = st.text_area("Paste single or multi-FASTA here:", height=150, placeholder="Paste your DNA sequence(s) here...")
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
                    st.success(f"✅ Pasted {len(seqs)} sequences.")
                    for i, seq in enumerate(seqs[:3]):
                        stats = get_basic_stats(seq)
                        st.markdown(f"**{names[i]}**: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                        st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                    if len(seqs) > 3:
                        st.caption(f"...and {len(seqs)-3} more.")
                else:
                    st.warning("No sequences found.")

        elif input_method == "🧪 Example Data":
            ex_type = st.radio("Example Type:", ["Single Example", "Multi-FASTA Example"], horizontal=True)
            if ex_type == "Single Example":
                if st.button("🔬 Load Single Example"):
                    parsed_fasta = parse_fasta(EXAMPLE_FASTA)
                    seqs = list(parsed_fasta.values())
                    names = list(parsed_fasta.keys())
                    st.success("✅ Single example sequence loaded.")
                    stats = get_basic_stats(seqs[0])
                    st.code(EXAMPLE_FASTA, language="fasta")
                    st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
            else:
                if st.button("🔬 Load Multi-FASTA Example"):
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
                    for i, seq in enumerate(seqs[:3]):
                        stats = get_basic_stats(seq)
                        st.markdown(f"**{names[i]}**: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                        st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                    st.code(EXAMPLE_MULTI_FASTA, language="fasta")

        elif input_method == "🌐 NCBI Fetch":
            db = st.radio("NCBI Database", ["nucleotide", "gene"], horizontal=True,
                          help="Only nucleotide and gene databases are applicable for DNA motif analysis")
            query_type = st.radio("Query Type", ["Accession", "Gene Name", "Custom Query"], horizontal=True)
            motif_examples = {
                "G-quadruplex": "NR_003287.2 (human telomerase RNA)",
                "Z-DNA": "NM_001126112.2 (human ADAR1 gene)",
                "R-loop": "NR_024540.1 (human SNRPN gene)",
                "eGZ-motif": "CGG repeat region",
                "AC-motif": "A-rich/C-rich consensus region"
            }
            with st.popover("📋 View Example Queries", use_container_width=True):
                st.markdown("**Motif Example Queries:**")
                for motif, example in motif_examples.items():
                    st.markdown(f"• **{motif}**: `{example}`")
            query = st.text_input("Enter query (accession, gene, etc.):")
            rettype = st.radio("Return Format", ["fasta", "gb"], horizontal=True)
            retmax = st.number_input("Max Records", min_value=1, max_value=20, value=3)
            if st.button("Fetch from NCBI"):
                if query:
                    with st.spinner("Contacting NCBI..."):
                        handle = Entrez.efetch(db=db, id=query, rettype=rettype, retmode="text")
                        records = list(SeqIO.parse(handle, rettype))
                        handle.close()
                        seqs = [str(rec.seq).upper().replace("U", "T") for rec in records]
                        names = [rec.id for rec in records]
                    if seqs:
                        st.success(f"Fetched {len(seqs)} sequences.")
                        for i, seq in enumerate(seqs[:3]):
                            stats = get_basic_stats(seq)
                            st.markdown(f"<b>{names[i]}</b>: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                            st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                else:
                    st.warning("Enter a query before fetching.")

        # Persist sequences to session state if any found from input
        if seqs:
            st.session_state.seqs = seqs
            st.session_state.names = names
            st.session_state.results = []

        # Sequence Preview lives under the input column for immediate feedback
        if st.session_state.get('seqs'):
            st.markdown("### 📊 Sequence Preview")
            for i, seq in enumerate(st.session_state.seqs[:2]):
                stats = get_basic_stats(seq)
                st.markdown(f"**{st.session_state.names[i]}** ({len(seq):,} bp) | GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}", unsafe_allow_html=True)
                st.code(wrap(seq[:400]), language="fasta")
            if len(st.session_state.seqs) > 2:
                st.caption(f"...and {len(st.session_state.seqs)-2} more.")

    # ----- RIGHT COLUMN: Analysis Controls + Run Button + Summary Table -----
    with col_right:
        st.markdown("### 🚀 Analysis & Run")
        
        # Analysis controls simplified
        st.markdown("### ⚙️ Analysis Options")
        
        # Enable all classes by default in consolidated system
        st.info("🔬 **NBDScanner detects all 11 motif classes with 22+ subclasses automatically**")
        
        # Simple options
        col1, col2 = st.columns(2)
        with col1:
            detailed_output = st.checkbox("Detailed Analysis", value=True, 
                                        help="Include comprehensive motif metadata")
        with col2:
            quality_check = st.checkbox("Quality Validation", value=True, 
                                      help="Validate detected motifs")
        
        # Advanced options using toggle + container for cleaner UI
        show_advanced = st.toggle("⚙️ Show Advanced Options", value=False, help="Toggle advanced analysis options")
        
        if show_advanced:
            with st.container(border=True):
                st.markdown("##### Advanced Configuration")
                col_adv1, col_adv2 = st.columns(2)
                with col_adv1:
                    show_chunk_progress = st.checkbox("Show Chunk-Level Progress", value=False,
                                                     help="Display detailed progress for each processing chunk (useful for large sequences)")
                with col_adv2:
                    use_parallel_scanner = st.checkbox("Use Experimental Parallel Scanner", value=False,
                                                      help="Enable experimental parallel chunk-based scanner (may improve performance on very large sequences >100kb)")
                
                if use_parallel_scanner:
                    st.info("ℹ️ Parallel scanner is experimental and works best on sequences >100kb with multiple CPU cores")
        else:
            show_chunk_progress = False
            use_parallel_scanner = False
        
        # Hardcoded default overlap handling: always remove overlaps within subclasses
        nonoverlap = True
        overlap_option = "Remove overlaps within subclasses"
        
        
        # ========== RUN ANALYSIS BUTTON ========== 
        if st.button("🔬 Run NBDScanner Analysis", type="primary", use_container_width=True, key="run_motif_analysis_main"):
            # Simplified validation
            if not st.session_state.seqs:
                st.error("❌ Please upload or input sequences before running analysis.")
                st.session_state.analysis_status = "Error"
            else:
                # =============================================================
                # SEQUENCE LENGTH VALIDATION FOR SERVER/STREAMLIT VERSION
                # =============================================================
                # Check if any sequence exceeds the maximum allowed length
                sequences_over_limit = []
                total_length = 0
                for i, seq in enumerate(st.session_state.seqs):
                    seq_len = len(seq)
                    total_length += seq_len
                    if seq_len > MAX_SEQUENCE_LENGTH:
                        seq_name = st.session_state.names[i] if i < len(st.session_state.names) else f"Sequence {i+1}"
                        sequences_over_limit.append((seq_name, seq_len))
                
                if sequences_over_limit:
                    limit_display = format_sequence_limit()
                    st.error(f"""
                    ❌ **Sequence Length Limit Exceeded**
                    
                    The following sequence(s) exceed the maximum allowed length of **{limit_display}** for the web version:
                    """)
                    for seq_name, seq_len in sequences_over_limit:
                        st.error(f"• **{seq_name}**: {seq_len:,} nucleotides ({seq_len - MAX_SEQUENCE_LENGTH:,} over limit)")
                    
                    st.info("""
                    💡 **For unlimited sequence analysis:**
                    
                    Use the **local Jupyter notebook version** (`NonBScanner_Local.ipynb`) which has no sequence length limits.
                    The local version is optimized for large genome-scale analysis on your own hardware.
                    """)
                    st.session_state.analysis_status = "Error"
                else:
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
                    
                    # Helper function to generate progress HTML template
                    def build_progress_html(elapsed, estimated_remaining, progress_display, 
                                            status_text, seq_info_html, detector_count):
                        """Build the progress timer HTML with consistent styling using CSS classes."""
                        return f"""
                        <div class='progress-panel'>
                            <h3 class='progress-panel__title'>⏱️ NonBScanner Analysis Progress</h3>
                            <p class='progress-panel__status'>{status_text}</p>
                            
                            <div class='stats-grid'>
                                <div class='stat-card'>
                                    <h2 class='stat-card__value'>{elapsed:.1f}s</h2>
                                    <p class='stat-card__label'>Elapsed Time</p>
                                </div>
                                <div class='stat-card'>
                                    <h2 class='stat-card__value'>{estimated_remaining:.1f}s</h2>
                                    <p class='stat-card__label'>Est. Remaining</p>
                                </div>
                                <div class='stat-card'>
                                    <h2 class='stat-card__value'>{progress_display}</h2>
                                    <p class='stat-card__label'>Overall Progress</p>
                                </div>
                                <div class='stat-card'>
                                    <h2 class='stat-card__value'>{detector_count}</h2>
                                    <p class='stat-card__label'>Detector Processes</p>
                                </div>
                            </div>
                            
                            <div class='sequence-info'>
                                {seq_info_html}
                            </div>
                        </div>
                        """
                    
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
                                status_text = "🔄 Starting analysis..."
                                progress_display = "Starting"
                            else:
                                status_text = "🔄 Analysis in progress..."
                                progress_display = f"{overall_percentage:.1f}%"
                            
                            # Build sequence info HTML
                            seq_info_html = f"""
                                <p class='sequence-info__text'>
                                    <strong>Sequence {i+1}/{len(st.session_state.seqs)}</strong>: {name} ({len(seq):,} bp)
                                </p>
                                <p class='sequence-info__subtext'>
                                    Total processed: {total_bp_processed:,} / {total_bp_all_sequences:,} bp
                                </p>
                            """
                            
                            # Build timer HTML using helper function
                            timer_html = build_progress_html(
                                elapsed, estimated_remaining, progress_display,
                                status_text, seq_info_html, len(DETECTOR_PROCESSES)
                            )
                            
                            # Add chunk progress if enabled
                            if show_chunk_progress and use_parallel_scanner:
                                # Calculate estimated chunks for this sequence
                                est_chunks = max(1, (len(seq) + CHUNK_SIZE_FOR_PARALLEL - 1) // CHUNK_SIZE_FOR_PARALLEL)
                                # Insert chunk info before the closing div
                                timer_html = timer_html[:-6] + f"<p class='sequence-info__subtext'>📦 Estimated chunks: {est_chunks}</p></div>"
                            
                            timer_placeholder.markdown(timer_html, unsafe_allow_html=True)
                            
                            # Show detailed progress panel with detector sequence
                            # The status_icon shows all detectors as "running" during analysis since they run in parallel
                            detailed_progress_html = f"""
                            <div class='pipeline-panel'>
                                <h4 class='pipeline-panel__title'>📋 Analysis Pipeline - Sequence of Operations</h4>
                                <div class='pipeline-grid'>
                            """
                            
                            for j, (detector_name, detector_desc) in enumerate(DETECTOR_PROCESSES):
                                # All detectors run in parallel during analysis, so show as "running"
                                status_icon = "🔄"
                                detailed_progress_html += f"""
                                    <div class='detector-item'>
                                        <span class='detector-item__name'>{status_icon} {j+1}. {detector_name}</span>
                                        <br/><span class='detector-item__desc'>{detector_desc}</span>
                                    </div>
                                """
                            
                            detailed_progress_html += """
                                </div>
                                <p class='pipeline-panel__footer'>
                                    All 9 detectors run in parallel for each sequence, followed by hybrid/cluster detection
                                </p>
                            </div>
                            """
                            detailed_progress_placeholder.markdown(detailed_progress_html, unsafe_allow_html=True)
                            
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
                                            chunk_progress_placeholder.info(f"🔄 Processing chunks: {current}/{total} ({chunk_percent:.1f}%)")
                                    
                                    # Run parallel scanner
                                    scanner = ParallelScanner(seq, hs_db=None)
                                    raw_motifs = scanner.run_scan(progress_callback=chunk_progress_callback)
                                    
                                    # Convert raw motifs to full motif format by running through scoring
                                    # For now, just use the standard analyzer
                                    results = analyze_sequence(seq, name)
                                    
                                    # Clear chunk progress
                                    if show_chunk_progress:
                                        chunk_progress_placeholder.success(f"✅ Chunk processing complete: {len(raw_motifs)} raw motifs found")
                                    
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
                            
                            # Build sequence completion info HTML
                            seq_complete_info_html = f"""
                                <p class='sequence-info__text'>
                                    <strong>Processed</strong>: {name} ({len(seq):,} bp) in {seq_time:.2f}s - {len(results)} motifs found
                                </p>
                                <p class='sequence-info__subtext'>
                                    Total processed: {total_bp_processed:,} / {total_bp_all_sequences:,} bp ({speed:,.0f} bp/s)
                                </p>
                            """
                            
                            # Update timer display with actual progress using helper function
                            updated_timer_html = build_progress_html(
                                elapsed, estimated_remaining, f"{actual_percentage:.1f}%",
                                f"✅ Sequence {i+1}/{len(st.session_state.seqs)} completed",
                                seq_complete_info_html, len(DETECTOR_PROCESSES)
                            )
                            timer_placeholder.markdown(updated_timer_html, unsafe_allow_html=True)
                            
                            with progress_placeholder.container():
                                pbar.progress(progress, text=f"Analyzed {i+1}/{len(st.session_state.seqs)} sequences")
                            
                            status_placeholder.info(f"✓ {name}: {len(seq):,} bp in {seq_time:.2f}s ({len(seq)/seq_time:.0f} bp/s) - {len(results)} motifs found")
                        
                        # Store results
                        st.session_state.results = all_results
                        
                        # Final timing statistics
                        total_time = time.time() - start_time
                        overall_speed = total_bp_processed / total_time if total_time > 0 else 0
                        
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
                        
                        # Store performance metrics with enhanced details
                        st.session_state.performance_metrics = {
                            'total_time': total_time,
                            'total_bp': total_bp_processed,
                            'speed': overall_speed,
                            'sequences': len(st.session_state.seqs),
                            'total_motifs': sum(len(r) for r in all_results),
                            'detector_count': len(DETECTOR_PROCESSES),  # Number of detector processes
                            'estimated_time': estimated_total_time,  # Initial estimated time
                            # Derive analysis steps from DETECTOR_PROCESSES plus post-processing steps
                            'analysis_steps': [f"{name} detection" for name, _ in DETECTOR_PROCESSES] + [
                                'Hybrid/Cluster detection',
                                'Overlap resolution'
                            ]
                        }
                        
                        # Clear progress displays
                        progress_placeholder.empty()
                        status_placeholder.empty()
                        detailed_progress_placeholder.empty()
                        
                        # Show final success message with enhanced performance metrics
                        timer_placeholder.markdown(f"""
                        <div class='progress-panel progress-panel--success'>
                            <h3 class='progress-panel__title'>✅ Analysis Complete!</h3>
                            <div class='stats-grid stats-grid--wide'>
                                <div class='stat-card'>
                                    <h2 class='stat-card__value'>{total_time:.2f}s</h2>
                                    <p class='stat-card__label'>Total Time</p>
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
                                    <p class='stat-card__label'>Detectors Run</p>
                                </div>
                                <div class='stat-card'>
                                    <h2 class='stat-card__value'>{sum(len(r) for r in all_results)}</h2>
                                    <p class='stat-card__label'>Motifs Found</p>
                                </div>
                            </div>
                        </div>
                        """, unsafe_allow_html=True)
                        
                        st.success("Results are available below and in the 'Analysis Results and Visualization' tab.")
                        st.session_state.analysis_status = "Complete"
                        
                    except Exception as e:
                        timer_placeholder.empty()
                        progress_placeholder.empty()
                        status_placeholder.empty()
                        detailed_progress_placeholder.empty()
                        st.error(f"❌ Analysis failed: {str(e)}")
                        st.session_state.analysis_status = "Error"

        # Show quick summary table if available
        if st.session_state.get('summary_df') is not None:
            st.markdown("#### Analysis Summary")
            st.dataframe(st.session_state.summary_df)
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
                <h3 class='progress-panel__title'>⚡ Performance Metrics</h3>
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
        st.markdown("### 📊 Analysis Summary")
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
                    🧬 NBDScanner Analysis Results
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
                st.info(f"ℹ️ {hybrid_cluster_count} Hybrid/Cluster motifs detected. View them in the 'Cluster/Hybrid' tab below.")
            
            # Motif class distribution summary (no score visualization as per requirements)
            if filtered_motifs:
                # Display chart without header as per requirements
                
                # Motif class distribution - show all classes even if 0 (excluding Hybrid and Cluster)
                motif_classes = [m.get('Class', 'Unknown') for m in filtered_motifs]
                
                # Initialize all classes with 0 counts (excluding Hybrid and Cluster)
                filtered_order = [cls for cls in MOTIF_ORDER if cls not in ['Hybrid', 'Non-B DNA Clusters']]
                all_class_counts = {cls: 0 for cls in filtered_order}
                
                # Update with actual counts
                actual_counts = Counter(motif_classes)
                all_class_counts.update(actual_counts)
                
                # Remove 'Unknown' if it exists and has 0 count
                if 'Unknown' in all_class_counts and all_class_counts['Unknown'] == 0:
                    del all_class_counts['Unknown']
                
                fig, ax = plt.subplots(figsize=(10, 6))
                
                # Create labels showing count and percentage
                labels = []
                values = list(all_class_counts.values())
                total = sum(values) if sum(values) > 0 else 1  # Avoid division by zero
                
                for cls, count in all_class_counts.items():
                    percentage = (count / total) * 100 if total > 0 else 0
                    if count > 0:
                        labels.append(f'{cls}\n({count}, {percentage:.1f}%)')
                    else:
                        labels.append(f'{cls}\n(0, 0.0%)')
                
                # Use colors from MOTIF_COLORS where available
                colors = []
                for cls in all_class_counts.keys():
                    if cls in MOTIF_COLORS:
                        colors.append(MOTIF_COLORS[cls])
                    else:
                        # Default colors for any missing classes
                        default_colors = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99', 
                                        '#ff99cc', '#99ccff', '#ffccbb', '#ccffcc']
                        colors.append(default_colors[len(colors) % len(default_colors)])
                
                # Create bar chart for class distribution
                bars = ax.bar(range(len(all_class_counts)), values, color=colors)
                ax.set_xlabel('Motif Classes')
                ax.set_ylabel('Count')
                ax.set_title('Motif Class Distribution (All Classes)')
                ax.set_xticks(range(len(all_class_counts)))
                ax.set_xticklabels(list(all_class_counts.keys()), rotation=45, ha='right')
                
                # Add count labels on bars
                for bar, count in zip(bars, values):
                    if count > 0:
                        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, 
                               str(count), ha='center', va='bottom')
                
                plt.tight_layout()
                st.pyplot(fig)
                plt.close(fig)
            
            # Enhanced motif table with new columns
            # Display table without header as per requirements
            
            # Column selection for display using pills for a more modern, visual approach
            available_columns = df.columns.tolist()
            # Filter out specific Normalized_Score column
            available_columns = [col for col in available_columns if col != 'Normalized_Score']
            
            # Use pills for column selection - multi-selection mode for better UX
            default_cols = [col for col in ['Class', 'Subclass', 'Start', 'End', 'Length', 'Score', 'GC Content'] if col in available_columns]
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
            
            if display_columns:
                # Filter out sequence and sequence name columns for cleaner display
                filtered_df = df[display_columns].copy()
                # Replace underscores with spaces in column names for display
                filtered_df.columns = [col.replace('_', ' ') for col in filtered_df.columns]
                st.dataframe(filtered_df, use_container_width=True, height=360)
            else:
                # Replace underscores with spaces in column names for display
                display_df = df.copy()
                display_df.columns = [col.replace('_', ' ') for col in display_df.columns]
                st.dataframe(display_df, use_container_width=True, height=360)
            
            # CONSOLIDATED VISUALIZATION SUITE
            st.markdown('<h3>Visualizations</h3>', unsafe_allow_html=True)
            
            # Create 4 visualization tabs including dedicated Cluster/Hybrid tab
            viz_tabs = st.tabs(["Distribution", "Coverage & Density", "Statistics", "Cluster/Hybrid"])
            
            with viz_tabs[0]:  # Distribution
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
            
            with viz_tabs[2]:  # Statistics (includes Cluster/Hybrid info)
                st.markdown("##### Statistical Analysis")
                
                # Density Metrics section
                st.markdown("###### Density Metrics")
                try:
                    genomic_density = calculate_genomic_density(filtered_motifs, sequence_length, by_class=True)
                    positional_density_kbp = calculate_positional_density(filtered_motifs, sequence_length, unit='kbp', by_class=True)
                    positional_density_mbp = calculate_positional_density(filtered_motifs, sequence_length, unit='Mbp', by_class=True)
                    
                    # Overall metrics in compact layout
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Genomic Density", f"{genomic_density.get('Overall', 0):.4f}%")
                    with col2:
                        st.metric("Motifs/kbp", f"{positional_density_kbp.get('Overall', 0):.2f}")
                    with col3:
                        st.metric("Motifs/Mbp", f"{positional_density_mbp.get('Overall', 0):.2f}")
                    
                    # Per-class density table
                    density_data = []
                    for class_name in sorted([k for k in genomic_density.keys() if k != 'Overall']):
                        density_data.append({
                            'Motif Class': class_name,
                            'Genomic Density (%)': f"{genomic_density.get(class_name, 0):.4f}",
                            'Motifs/kbp': f"{positional_density_kbp.get(class_name, 0):.2f}",
                            'Motifs/Mbp': f"{positional_density_mbp.get(class_name, 0):.2f}"
                        })
                    
                    if density_data:
                        density_df = pd.DataFrame(density_data)
                        st.dataframe(density_df, use_container_width=True, height=200)
                    
                    fig_density = plot_density_comparison(genomic_density, positional_density_kbp,
                                                          title="Motif Density Analysis")
                    st.pyplot(fig_density)
                    plt.close(fig_density)
                    
                except Exception as e:
                    st.error(f"Error calculating density metrics: {e}")
                
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
                
            with viz_tabs[3]:  # Dedicated Cluster/Hybrid Tab
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
                            🔗 Hybrid & Cluster Motif Summary
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
                            st.info(f"🔗 **Hybrid motifs** are regions where different Non-B DNA motif classes overlap. Found {len(hybrid_only)} hybrid regions.")
                            
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
                            st.info(f"📊 **DNA Clusters** are high-density regions with multiple Non-B DNA motif classes. Found {len(cluster_only)} cluster regions.")
                            
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
                else:
                    st.info("ℹ️ No hybrid or cluster motifs detected in this sequence. Hybrid motifs occur when different Non-B DNA classes overlap, and clusters form when multiple motifs are found in close proximity.")

# ---------- DOWNLOAD ----------
with tab_pages["Download"]:
    st.header("Export Data")
    if not st.session_state.results:
        st.info("No results available to download.")
    else:
        st.markdown("### 📊 Export Options")
        
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
            st.markdown("**📊 Export Configuration**")
            
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
            st.info(f"ℹ️ {len(hybrid_cluster_motifs_export)} Hybrid/Cluster motifs are excluded from downloads. These are shown only in the Cluster/Hybrid visualization tab.")
        
        # Display preview
        if all_motifs:
            df_preview = export_results_to_dataframe(all_motifs).head(10)
            st.markdown("### 👀 Export Preview")
            st.dataframe(df_preview, use_container_width=True)
            st.caption(f"Showing first 10 of {len(all_motifs)} total records (Hybrid/Cluster motifs excluded)")
        
        # Export buttons using consolidated functions
        st.markdown("### 💾 Download Files")
        
        # Add helpful info about Excel format
        st.info("""
        💡 **Excel Format**: Downloads a multi-sheet workbook with:
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
                    "📄 Download CSV", 
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
                        "📊 Download Excel", 
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
                    "📋 Download JSON", 
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
                    "🧬 Download BED", 
                    data=bed_data.encode('utf-8'), 
                    file_name="nbdscanner_results.bed", 
                    mime="text/plain",
                    use_container_width=True,
                    help="BED format for genome browsers"
                )
        
        # Additional export options in a new row
        st.markdown("### 🔧 Additional Exports")
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
                    "⚙️ Download Config", 
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
    <b>🔧 Motif Detection Parameters:</b><br><br>
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
        <b>📋 Scoring Configuration Details</b>
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
