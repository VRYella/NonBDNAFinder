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

# Color theme palettes for scientific styling
COLOR_THEMES = {
    'scientific_blue': {
        'primary': '#0d47a1',
        'secondary': '#1976d2',
        'accent': '#42a5f5',
        'bg_light': '#f5f9ff',
        'bg_card': '#e3f2fd',
        'text': '#1e293b'
    },
    'nature_green': {
        'primary': '#1b5e20',
        'secondary': '#388e3c',
        'accent': '#66bb6a',
        'bg_light': '#f1f8e9',
        'bg_card': '#dcedc8',
        'text': '#1e293b'
    },
    'genomic_purple': {
        'primary': '#4a148c',
        'secondary': '#7b1fa2',
        'accent': '#ba68c8',
        'bg_light': '#f3e5f5',
        'bg_card': '#e1bee7',
        'text': '#1e293b'
    },
    'clinical_teal': {
        'primary': '#006064',
        'secondary': '#00838f',
        'accent': '#26c6da',
        'bg_light': '#e0f7fa',
        'bg_card': '#b2ebf2',
        'text': '#1e293b'
    }
}

# Get current theme colors
current_theme = COLOR_THEMES.get(st.session_state.color_theme, COLOR_THEMES['scientific_blue'])
is_dark_mode = st.session_state.theme_mode == 'dark'
is_compact = st.session_state.table_density == 'compact'

# Dark mode color overrides
if is_dark_mode:
    current_theme = {
        **current_theme,
        'bg_light': '#0f172a',
        'bg_card': '#1e293b',
        'text': '#e2e8f0'
    }

# Pre-calculate RGB values for all theme colors (performance optimization)
rgb = {key: hex_to_rgb(value) for key, value in current_theme.items()}

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
       ELEGANT BACKGROUND PATTERN WITH DNA MOTIF
       ============================================ */
    body, [data-testid="stAppViewContainer"], .main {{
        background: 
            {dna_pattern},
            linear-gradient(135deg, {current_theme['bg_light']} 0%, {'#1e293b' if is_dark_mode else '#e8f4fd'} 50%, {current_theme['bg_light']} 100%) !important;
        font-family: 'Poppins', 'Inter', 'IBM Plex Sans', 'Segoe UI', system-ui, -apple-system, sans-serif !important;
        {'color: #e2e8f0 !important;' if is_dark_mode else ''}
    }}
    
    /* ============================================
       TABS: Modern elegant design with glass effect
       ============================================ */
    .stTabs [data-baseweb="tab-list"] {{
        width: 100% !important;
        justify-content: center !important;
        gap: 8px !important;
        border-bottom: none !important;
        background: {'linear-gradient(135deg, rgba(30, 41, 59, 0.9) 0%, rgba(51, 65, 85, 0.85) 100%)' if is_dark_mode else f"linear-gradient(135deg, rgba(255, 255, 255, 0.95) 0%, rgba({rgb['bg_card'][0]}, {rgb['bg_card'][1]}, {rgb['bg_card'][2]}, 0.9) 100%)"} !important;
        backdrop-filter: blur(20px) !important;
        -webkit-backdrop-filter: blur(20px) !important;
        box-shadow: 0 8px 32px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.15), 0 2px 8px rgba(0, 0, 0, 0.05) !important;
        margin-bottom: 2em;
        padding: 12px 16px !important;
        border-radius: 20px !important;
        animation: fade-in 0.6s ease-out;
    }}
    
    .stTabs [data-baseweb="tab"] {{
        font-size: 1.05rem !important;
        font-weight: 600 !important;
        font-family: 'Poppins', 'Inter', sans-serif !important;
        flex: 0 1 auto !important;
        min-width: 120px !important;
        padding: 14px 24px !important;
        text-align: center;
        color: {'#94a3b8' if is_dark_mode else '#64748b'} !important;
        background: transparent !important;
        border: none !important;
        border-radius: 14px !important;
        letter-spacing: 0.02em;
        transition: all 0.4s cubic-bezier(0.4, 0, 0.2, 1);
        position: relative;
        overflow: hidden;
    }}
    
    .stTabs [data-baseweb="tab"]::before {{
        content: '';
        position: absolute;
        inset: 0;
        background: linear-gradient(135deg, {current_theme['primary']} 0%, {current_theme['secondary']} 50%, {current_theme['accent']} 100%);
        opacity: 0;
        transition: opacity 0.4s ease;
        border-radius: 14px;
        z-index: -1;
    }}
    
    .stTabs [data-baseweb="tab"]:hover {{
        color: {current_theme['primary']} !important;
        transform: translateY(-2px);
        background: {'rgba(255, 255, 255, 0.05)' if is_dark_mode else f"rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.08)"} !important;
    }}
    
    .stTabs [aria-selected="true"] {{
        color: white !important;
        font-weight: 700 !important;
        transform: translateY(-2px);
        background: linear-gradient(135deg, {current_theme['primary']} 0%, {current_theme['secondary']} 100%) !important;
        box-shadow: 0 6px 20px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.35), 0 2px 8px rgba({rgb['secondary'][0]}, {rgb['secondary'][1]}, {rgb['secondary'][2]}, 0.2) !important;
    }}
    
    .stTabs [aria-selected="true"]::before {{
        opacity: 1;
    }}
    
    /* ============================================
       HEADINGS: Modern typography with gradients
       ============================================ */
    h1, h2, h3, h4 {{
        font-family: 'Poppins', 'IBM Plex Sans', 'Inter', 'Segoe UI', system-ui, sans-serif !important;
        color: {current_theme['primary']} !important;
        font-weight: 700 !important;
        letter-spacing: -0.01em;
        margin-top: 1.2em;
        margin-bottom: 0.7em;
        animation: fade-in 0.4s ease-out;
    }}
    h1 {{ 
        font-size: 2.8rem !important; 
        font-weight: 800 !important;
        background: linear-gradient(135deg, {current_theme['primary']} 0%, {current_theme['secondary']} 50%, {current_theme['accent']} 100%);
        background-size: 200% 200%;
        animation: gradient-shift 4s ease infinite, fade-in 0.5s ease-out;
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
        letter-spacing: -0.02em;
    }}
    h2 {{ 
        font-size: 1.9rem !important; 
        color: {current_theme['primary']} !important; 
        font-weight: 700 !important;
        border-bottom: none !important;
        padding-bottom: 0.4rem;
        position: relative;
    }}
    h2::after {{
        content: '';
        position: absolute;
        bottom: 0;
        left: 0;
        width: 60px;
        height: 4px;
        background: linear-gradient(90deg, {current_theme['primary']} 0%, {current_theme['accent']} 100%);
        border-radius: 2px;
    }}
    h3 {{ 
        font-size: 1.5rem !important; 
        color: {current_theme['secondary']} !important; 
        font-weight: 600 !important;
    }}
    h4 {{ 
        font-size: 1.25rem !important; 
        color: {current_theme['secondary']} !important; 
        font-weight: 600 !important;
    }}
    
    /* Body text: Modern readability with elegant contrast */
    .stMarkdown, .markdown-text-container, .stText, p, span, label {{
        font-size: 1.02rem !important;
        font-family: 'Poppins', 'Inter', 'Segoe UI', system-ui, sans-serif !important;
        line-height: 1.7 !important;
        color: {current_theme['text']} !important;
        margin-bottom: 0.6em !important;
        font-weight: 400;
    }}
    
    /* Input fields: Premium scientific design with elegant borders */
    input, .stTextInput>div>div>input, .stSelectbox>div>div>div, 
    .stMultiSelect>div>div>div, .stRadio>div>div>label>div {{
        font-size: 1.0rem !important;
        font-family: 'Inter', 'IBM Plex Sans', 'Segoe UI', system-ui, sans-serif !important;
        border-radius: 10px !important;
        border: 2px solid {'#475569' if is_dark_mode else '#e1e8ed'} !important;
        background: {'#1e293b' if is_dark_mode else '#ffffff'} !important;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.02);
    }}
    input:focus, .stTextInput>div>div>input:focus {{
        border-color: {current_theme['secondary']} !important;
        box-shadow: 0 0 0 4px rgba({rgb['secondary'][0]}, {rgb['secondary'][1]}, {rgb['secondary'][2]}, 0.12), 0 4px 8px rgba({rgb['secondary'][0]}, {rgb['secondary'][1]}, {rgb['secondary'][2]}, 0.08) !important;
        background: {'#1e293b' if is_dark_mode else '#fefeff'} !important;
        transform: translateY(-1px);
    }}
    
    /* ============================================
       BUTTONS: Enhanced with click animations
       ============================================ */
    .stButton>button {{
        font-size: 1.05rem !important;
        font-family: 'Inter', 'IBM Plex Sans', 'Segoe UI', system-ui, sans-serif !important;
        padding: 0.75em 1.8em !important;
        background: linear-gradient(135deg, {current_theme['primary']} 0%, {current_theme['secondary']} 50%, {current_theme['accent']} 100%) !important;
        color: #fff !important;
        border-radius: 12px !important;
        border: none !important;
        font-weight: 600 !important;
        box-shadow: 0 6px 18px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.35), 0 2px 6px rgba({rgb['secondary'][0]}, {rgb['secondary'][1]}, {rgb['secondary'][2]}, 0.2);
        transition: all 0.35s cubic-bezier(0.4, 0, 0.2, 1);
        letter-spacing: 0.03em;
        position: relative;
        overflow: hidden;
    }}
    .stButton>button::before {{
        content: '';
        position: absolute;
        top: 0;
        left: -100%;
        width: 100%;
        height: 100%;
        background: linear-gradient(90deg, transparent, rgba(255,255,255,0.3), transparent);
        transition: left 0.5s;
    }}
    .stButton>button:hover::before {{
        left: 100%;
    }}
    .stButton>button:hover {{
        background: linear-gradient(135deg, {current_theme['primary']} 0%, {current_theme['secondary']} 100%) !important;
        box-shadow: 0 8px 24px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.45), 0 4px 12px rgba({rgb['secondary'][0]}, {rgb['secondary'][1]}, {rgb['secondary'][2]}, 0.3);
        transform: translateY(-3px) scale(1.02);
    }}
    .stButton>button:active {{
        transform: translateY(-1px) scale(0.97);
        box-shadow: 0 4px 12px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.35);
        animation: scale-click 0.15s ease;
    }}
    
    /* ============================================
       DATAFRAMES: Enhanced with compact/relaxed mode
       ============================================ */
    .stDataFrame, .stTable {{
        font-size: {'0.88rem' if is_compact else '0.96rem'} !important;
        font-family: 'IBM Plex Sans', 'Inter', 'Segoe UI', system-ui, sans-serif !important;
        line-height: {'1.4' if is_compact else '1.65'} !important;
        border-radius: 12px !important;
        overflow: hidden;
        box-shadow: 0 4px 12px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.08), 0 2px 4px rgba(0, 0, 0, 0.03) !important;
        border: 1px solid {current_theme['bg_card']} !important;
        animation: fade-in 0.4s ease-out;
    }}
    .stDataFrame thead tr th {{
        background: linear-gradient(135deg, {current_theme['primary']} 0%, {current_theme['secondary']} 100%) !important;
        color: white !important;
        font-weight: 700 !important;
        padding: {'0.6rem' if is_compact else '1rem'} !important;
        font-size: {'0.85rem' if is_compact else '0.98rem'} !important;
        letter-spacing: 0.02em;
        text-transform: uppercase;
        border-bottom: 3px solid {current_theme['primary']} !important;
    }}
    .stDataFrame tbody tr {{
        transition: all 0.25s ease;
        border-bottom: 1px solid {current_theme['bg_card']} !important;
    }}
    .stDataFrame tbody tr:hover {{
        background: {'linear-gradient(90deg, #1e293b 0%, #334155 100%)' if is_dark_mode else f"linear-gradient(90deg, #f0f7ff 0%, {current_theme['bg_card']} 100%)"} !important;
        transform: scale(1.005);
        box-shadow: 0 2px 8px rgba({rgb['secondary'][0]}, {rgb['secondary'][1]}, {rgb['secondary'][2]}, 0.1);
    }}
    .stDataFrame tbody tr td {{
        padding: {'0.5rem' if is_compact else '0.85rem'} !important;
        color: {current_theme['text']} !important;
    }}
    
    /* Highlight alternate rows */
    .stDataFrame tbody tr:nth-child(odd) {{
        background: {'rgba(30, 41, 59, 0.3)' if is_dark_mode else 'rgba(227, 242, 253, 0.3)'} !important;
    }}
    
    /* Tab content spacing */
    .stTabs [data-baseweb="tab-panel"] {{
        padding-top: 2.5rem !important;
        padding-left: 1.5rem !important;
        padding-right: 1.5rem !important;
        padding-bottom: 2rem !important;
        animation: fade-in 0.4s ease-out;
    }}
    
    /* ============================================
       GLASSMORPHISM CARDS: For summary displays
       ============================================ */
    .analysis-summary-card, .glassmorphism-card {{
        margin: 1.5rem 0 !important;
        padding: 2rem !important;
        background: {'rgba(30, 41, 59, 0.8)' if is_dark_mode else 'rgba(255, 255, 255, 0.85)'} !important;
        backdrop-filter: blur(10px) !important;
        -webkit-backdrop-filter: blur(10px) !important;
        border-radius: 16px !important;
        box-shadow: 0 8px 32px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.12) !important;
        border: 1px solid {'rgba(148, 163, 184, 0.2)' if is_dark_mode else f"rgba({rgb['secondary'][0]}, {rgb['secondary'][1]}, {rgb['secondary'][2]}, 0.15)"} !important;
        animation: fade-in 0.5s ease-out;
        transition: all 0.3s ease;
    }}
    .analysis-summary-card:hover, .glassmorphism-card:hover {{
        transform: translateY(-2px);
        box-shadow: 0 12px 40px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.18) !important;
    }}
    
    /* Select boxes and multiselect: Premium dropdown design */
    .stSelectbox > div > div > div, .stMultiSelect > div > div > div {{
        min-height: 3.5rem !important;
        padding: 1rem 1.2rem !important;
        line-height: 1.65 !important;
        border-radius: 12px !important;
        background: {'linear-gradient(135deg, #1e293b 0%, #334155 100%)' if is_dark_mode else f"linear-gradient(135deg, #ffffff 0%, {current_theme['bg_light']} 100%)"} !important;
        border: 2px solid {'#475569' if is_dark_mode else '#e1e8ed'} !important;
        box-shadow: 0 2px 6px rgba(0, 0, 0, 0.03);
        transition: all 0.3s ease;
    }}
    .stSelectbox > div > div > div:hover, .stMultiSelect > div > div > div:hover {{
        border-color: {current_theme['accent']} !important;
        box-shadow: 0 4px 12px rgba({rgb['accent'][0]}, {rgb['accent'][1]}, {rgb['accent'][2]}, 0.12);
        transform: translateY(-1px);
    }}
    
    /* Enhanced input field spacing with elegant styling */
    .stTextInput > div > div > input, .stTextArea textarea {{
        padding: 1rem 1.3rem !important;
        min-height: 3.5rem !important;
        line-height: 1.6 !important;
        border-radius: 12px !important;
        font-size: 1.02rem !important;
        background: {'linear-gradient(135deg, #1e293b 0%, #334155 100%)' if is_dark_mode else f"linear-gradient(135deg, #ffffff 0%, {current_theme['bg_light']} 100%)"} !important;
    }}
    .stTextArea textarea {{
        min-height: 120px !important;
    }}
    
    /* Number input styling with scientific precision */
    .stNumberInput > div > div > input {{
        padding: 1rem 1.3rem !important;
        min-height: 3.5rem !important;
        border-radius: 12px !important;
        font-size: 1.02rem !important;
        font-weight: 600;
        background: {'linear-gradient(135deg, #1e293b 0%, #334155 100%)' if is_dark_mode else f"linear-gradient(135deg, #ffffff 0%, {current_theme['bg_light']} 100%)"} !important;
    }}
    
    /* Radio buttons: Enhanced with elegant hover states and borders */
    .stRadio > div {{
        gap: 0.4rem !important;
    }}
    .stRadio > div > div > label {{
        margin-bottom: 0.7rem !important;
        padding: 0.65rem 1rem !important;
        border-radius: 10px !important;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        border: 2px solid transparent;
        background: {'linear-gradient(135deg, #1e293b 0%, #334155 100%)' if is_dark_mode else f"linear-gradient(135deg, #ffffff 0%, {current_theme['bg_light']} 100%)"};
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.03);
        cursor: pointer;
    }}
    .stRadio > div > div > label:hover {{
        background: {'linear-gradient(135deg, #334155 0%, #475569 100%)' if is_dark_mode else f"linear-gradient(135deg, {current_theme['bg_card']} 0%, {current_theme['bg_light']} 100%)"} !important;
        border-color: {current_theme['bg_card']} !important;
        transform: translateX(4px);
        box-shadow: 0 4px 8px rgba({rgb['secondary'][0]}, {rgb['secondary'][1]}, {rgb['secondary'][2]}, 0.08);
    }}
    .stRadio > div > div > label[data-checked="true"] {{
        background: {'linear-gradient(135deg, #334155 0%, #475569 100%)' if is_dark_mode else f"linear-gradient(135deg, {current_theme['bg_card']} 0%, {current_theme['bg_card']} 100%)"} !important;
        border-color: {current_theme['secondary']} !important;
        font-weight: 600 !important;
        box-shadow: 0 4px 12px rgba({rgb['secondary'][0]}, {rgb['secondary'][1]}, {rgb['secondary'][2]}, 0.15);
    }}
    /* Radio button circles */
    .stRadio > div > div > label > div:first-child {{
        border-width: 2px !important;
        border-color: {current_theme['accent']} !important;
        width: 22px !important;
        height: 22px !important;
    }}
    .stRadio > div > div > label[data-checked="true"] > div:first-child {{
        border-color: {current_theme['primary']} !important;
        background-color: {current_theme['primary']} !important;
    }}
    
    /* Number input styling */
    .stNumberInput > div > div > input {{
        padding: 0.85rem 1.1rem !important;
        min-height: 3rem !important;
        border-radius: 8px !important;
    }}
    
    /* Info boxes and alerts: Scientific notification design */
    .stAlert {{
        border-radius: 12px !important;
        border-left: 5px solid {current_theme['secondary']} !important;
        box-shadow: 0 4px 12px rgba(0,0,0,0.08) !important;
        background: {'linear-gradient(90deg, #1e293b 0%, #334155 100%)' if is_dark_mode else f"linear-gradient(90deg, {current_theme['bg_card']} 0%, {current_theme['bg_light']} 100%)"} !important;
        padding: 1.2rem !important;
        font-weight: 500;
    }}
    .stAlert[data-baseweb="notification"] {{
        border-left-color: {current_theme['primary']} !important;
    }}
    .stSuccess {{
        border-left-color: #2e7d32 !important;
        background: {'linear-gradient(90deg, #1e3a29 0%, #1e293b 100%)' if is_dark_mode else 'linear-gradient(90deg, #e8f5e9 0%, #f1f8f4 100%)'} !important;
    }}
    .stWarning {{
        border-left-color: #f57c00 !important;
        background: {'linear-gradient(90deg, #3a2e1e 0%, #1e293b 100%)' if is_dark_mode else 'linear-gradient(90deg, #fff3e0 0%, #fff8f1 100%)'} !important;
    }}
    .stError {{
        border-left-color: #c62828 !important;
        background: {'linear-gradient(90deg, #3a1e1e 0%, #1e293b 100%)' if is_dark_mode else 'linear-gradient(90deg, #ffebee 0%, #fff5f5 100%)'} !important;
    }}
    
    /* Progress bars: Elegant animated design */
    .stProgress > div > div {{
        background: linear-gradient(90deg, {current_theme['primary']} 0%, {current_theme['secondary']} 50%, {current_theme['accent']} 100%) !important;
        border-radius: 12px !important;
        box-shadow: 0 2px 8px rgba({rgb['secondary'][0]}, {rgb['secondary'][1]}, {rgb['secondary'][2]}, 0.3);
        animation: shimmer 2s infinite;
    }}
    @keyframes shimmer {{
        0% {{ background-position: -1000px 0; }}
        100% {{ background-position: 1000px 0; }}
    }}
    .stProgress > div {{
        border-radius: 12px !important;
        background: {current_theme['bg_card']} !important;
        box-shadow: inset 0 2px 4px rgba(0, 0, 0, 0.05);
    }}
    
    /* File uploader: Premium drag-and-drop design */
    [data-testid="stFileUploader"] {{
        border: 3px dashed {current_theme['accent']} !important;
        border-radius: 16px !important;
        background: {'linear-gradient(135deg, rgba(30, 41, 59, 0.3) 0%, rgba(51, 65, 85, 0.2) 100%)' if is_dark_mode else f"linear-gradient(135deg, rgba({rgb['bg_card'][0]}, {rgb['bg_card'][1]}, {rgb['bg_card'][2]}, 0.3) 0%, rgba({rgb['bg_card'][0]}, {rgb['bg_card'][1]}, {rgb['bg_card'][2]}, 0.2) 100%)"} !important;
        padding: 2.5rem !important;
        transition: all 0.35s cubic-bezier(0.4, 0, 0.2, 1);
        box-shadow: 0 4px 12px rgba({rgb['accent'][0]}, {rgb['accent'][1]}, {rgb['accent'][2]}, 0.1);
    }}
    [data-testid="stFileUploader"]:hover {{
        border-color: {current_theme['primary']} !important;
        background: {'linear-gradient(135deg, rgba(30, 41, 59, 0.5) 0%, rgba(51, 65, 85, 0.4) 100%)' if is_dark_mode else f"linear-gradient(135deg, rgba({rgb['bg_card'][0]}, {rgb['bg_card'][1]}, {rgb['bg_card'][2]}, 0.5) 0%, rgba({rgb['bg_card'][0]}, {rgb['bg_card'][1]}, {rgb['bg_card'][2]}, 0.4) 100%)"} !important;
        transform: scale(1.01);
        box-shadow: 0 6px 20px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.15);
    }}
    [data-testid="stFileUploader"] section {{
        border: none !important;
    }}
    [data-testid="stFileUploader"] button {{
        background: linear-gradient(135deg, {current_theme['secondary']} 0%, {current_theme['accent']} 100%) !important;
        color: white !important;
        border: none !important;
        border-radius: 10px !important;
        padding: 0.6rem 1.5rem !important;
        font-weight: 600 !important;
        transition: all 0.3s ease;
    }}
    [data-testid="stFileUploader"] button:hover {{
        background: linear-gradient(135deg, {current_theme['primary']} 0%, {current_theme['secondary']} 100%) !important;
        transform: translateY(-2px);
        box-shadow: 0 6px 16px rgba({rgb['secondary'][0]}, {rgb['secondary'][1]}, {rgb['secondary'][2]}, 0.3);
    }}
    
    /* Checkboxes: Modern toggle-style design */
    .stCheckbox {{
        padding: 0.5rem 0 !important;
    }}
    .stCheckbox > label {{
        padding: 0.6rem 0.8rem !important;
        border-radius: 10px !important;
        transition: all 0.3s ease;
        cursor: pointer;
        background: {'linear-gradient(135deg, #1e293b 0%, #334155 100%)' if is_dark_mode else f"linear-gradient(135deg, #ffffff 0%, {current_theme['bg_light']} 100%)"};
        border: 2px solid transparent;
    }}
    .stCheckbox > label:hover {{
        background: {'linear-gradient(135deg, #334155 0%, #475569 100%)' if is_dark_mode else f"linear-gradient(135deg, {current_theme['bg_card']} 0%, {current_theme['bg_light']} 100%)"} !important;
        border-color: {current_theme['bg_card']} !important;
        transform: translateX(3px);
    }}
    .stCheckbox > label > div:first-child {{
        border-width: 2px !important;
        border-color: {current_theme['accent']} !important;
        border-radius: 6px !important;
        width: 22px !important;
        height: 22px !important;
        transition: all 0.3s ease;
    }}
    .stCheckbox > label > div:first-child[data-checked="true"] {{
        background-color: {current_theme['primary']} !important;
        border-color: {current_theme['primary']} !important;
    }}
    
    /* Expander styling: Scientific accordion design with clean chevron */
    .streamlit-expanderHeader {{
        border-radius: 12px !important;
        background: {'linear-gradient(135deg, #334155 0%, #475569 100%)' if is_dark_mode else f"linear-gradient(135deg, {current_theme['bg_card']} 0%, {current_theme['bg_light']} 100%)"} !important;
        font-weight: 700 !important;
        padding: 1rem 1.5rem !important;
        border: 2px solid {current_theme['bg_card']} !important;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        box-shadow: 0 2px 6px rgba({rgb['secondary'][0]}, {rgb['secondary'][1]}, {rgb['secondary'][2]}, 0.08);
        display: flex !important;
        align-items: center !important;
        gap: 12px !important;
    }}
    .streamlit-expanderHeader:hover {{
        background: {'linear-gradient(135deg, #475569 0%, #334155 100%)' if is_dark_mode else f"linear-gradient(135deg, {current_theme['bg_card']} 0%, {current_theme['bg_card']} 100%)"} !important;
        border-color: {current_theme['secondary']} !important;
        transform: translateX(4px);
        box-shadow: 0 4px 12px rgba({rgb['secondary'][0]}, {rgb['secondary'][1]}, {rgb['secondary'][2]}, 0.15);
    }}
    .streamlit-expanderContent {{
        border-radius: 0 0 12px 12px !important;
        background: {'#1e293b' if is_dark_mode else '#fefeff'} !important;
        border: 2px solid {current_theme['bg_card']} !important;
        border-top: none !important;
        padding: 1.5rem !important;
    }}
    
    /* Clean chevron icon styling for expanders - replaces keyboard_arrow icons */
    .streamlit-expanderHeader svg {{
        margin-right: 10px !important;
        flex-shrink: 0 !important;
        width: 20px !important;
        height: 20px !important;
        color: {current_theme['primary']} !important;
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
        font-size: 1.2rem !important;
    }}
    
    /* Fix: keyboard_arrow_down icon text overlap with flex container and proper line-height */
    .streamlit-expanderHeader {{
        display: flex !important;
        align-items: center !important;
        line-height: 1.5 !important;
    }}
    
    /* Ensure label text doesn't overlap with icon - enhanced flex solution */
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
        width: 24px !important;
        min-width: 24px !important;
        overflow: hidden !important;
    }}
    
    /* ============================================
       METRIC CARDS: Glassmorphism with animations
       ============================================ */
    [data-testid="stMetric"] {{
        background: {'rgba(30, 41, 59, 0.7)' if is_dark_mode else 'rgba(255, 255, 255, 0.8)'} !important;
        backdrop-filter: blur(12px) !important;
        -webkit-backdrop-filter: blur(12px) !important;
        padding: 1.5rem !important;
        border-radius: 16px !important;
        border: 2px solid {'rgba(148, 163, 184, 0.2)' if is_dark_mode else current_theme['bg_card']} !important;
        box-shadow: 0 8px 32px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.1) !important;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        animation: fade-in 0.5s ease-out;
    }}
    [data-testid="stMetric"]:hover {{
        transform: translateY(-4px) scale(1.02);
        box-shadow: 0 12px 40px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.2) !important;
        border-color: {current_theme['accent']} !important;
        animation: pulse-glow 2s infinite;
    }}
    [data-testid="stMetricValue"] {{
        font-size: 2.2rem !important;
        font-weight: 800 !important;
        background: linear-gradient(135deg, {current_theme['primary']} 0%, {current_theme['secondary']} 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
    }}
    [data-testid="stMetricLabel"] {{
        font-size: 0.95rem !important;
        font-weight: 600 !important;
        color: {'#94a3b8' if is_dark_mode else '#455a64'} !important;
        text-transform: uppercase;
        letter-spacing: 0.05em;
    }}
    
    /* Spinner */
    .stSpinner > div {{
        border-top-color: {current_theme['secondary']} !important;
    }}
    
    /* Code blocks: Scientific monospace design */
    .stCodeBlock {{
        border-radius: 12px !important;
        border: 2px solid {current_theme['bg_card']} !important;
        background: {'linear-gradient(135deg, #1e293b 0%, #0f172a 100%)' if is_dark_mode else f"linear-gradient(135deg, #f8fbff 0%, {current_theme['bg_light']} 100%)"} !important;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.03);
        padding: 1rem !important;
    }}
    code {{
        font-family: 'JetBrains Mono', 'Fira Code', 'Consolas', monospace !important;
        font-size: 0.92rem !important;
        color: {current_theme['primary']} !important;
        background: {'rgba(30, 41, 59, 0.4)' if is_dark_mode else f"rgba({rgb['bg_card'][0]}, {rgb['bg_card'][1]}, {rgb['bg_card'][2]}, 0.4)"} !important;
        padding: 0.2rem 0.5rem !important;
        border-radius: 6px !important;
    }}
    
    /* Sidebar styling for navigation */
    [data-testid="stSidebar"] {{
        background: {'linear-gradient(180deg, #1e293b 0%, #0f172a 100%)' if is_dark_mode else f"linear-gradient(180deg, {current_theme['bg_card']} 0%, {current_theme['bg_light']} 100%)"} !important;
        border-right: 2px solid {'#334155' if is_dark_mode else current_theme['bg_card']} !important;
    }}
    [data-testid="stSidebar"] .stMarkdown {{
        color: {current_theme['primary']} !important;
    }}
    
    /* Scrollbar styling: Elegant scientific design */
    ::-webkit-scrollbar {{
        width: 12px;
        height: 12px;
    }}
    ::-webkit-scrollbar-track {{
        background: {'linear-gradient(135deg, #1e293b 0%, #0f172a 100%)' if is_dark_mode else f"linear-gradient(135deg, #f1f5f9 0%, {current_theme['bg_card']} 100%)"};
        border-radius: 12px;
        border: 1px solid {current_theme['bg_card']};
    }}
    ::-webkit-scrollbar-thumb {{
        background: linear-gradient(135deg, {current_theme['primary']} 0%, {current_theme['secondary']} 50%, {current_theme['accent']} 100%);
        border-radius: 12px;
        border: 2px solid {'#1e293b' if is_dark_mode else '#f8fbff'};
        box-shadow: 0 2px 6px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.2);
    }}
    ::-webkit-scrollbar-thumb:hover {{
        background: linear-gradient(135deg, {current_theme['primary']} 0%, {current_theme['secondary']} 100%);
        box-shadow: 0 3px 10px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.3);
    }}
    
    /* Loading spinner */
    .stSpinner > div {{
        border-top-color: {current_theme['primary']} !important;
        border-right-color: {current_theme['secondary']} !important;
        border-bottom-color: {current_theme['accent']} !important;
    }}
    
    /* Download buttons enhanced */
    [data-testid="stDownloadButton"] button {{
        background: linear-gradient(135deg, {current_theme['primary']} 0%, {current_theme['secondary']} 100%) !important;
        color: white !important;
        border-radius: 12px !important;
        padding: 0.7rem 1.6rem !important;
        font-weight: 600 !important;
        box-shadow: 0 4px 12px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.25);
        transition: all 0.3s ease;
    }}
    [data-testid="stDownloadButton"] button:hover {{
        background: linear-gradient(135deg, {current_theme['primary']} 0%, {current_theme['secondary']} 100%) !important;
        transform: translateY(-2px);
        box-shadow: 0 6px 18px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.35);
    }}
    
    /* Caption text styling */
    .caption, [data-testid="stCaptionContainer"] {{
        color: {'#94a3b8' if is_dark_mode else '#607d8b'} !important;
        font-size: 0.9rem !important;
        font-style: italic;
    }}
    
    /* Success/Error message boxes */
    .element-container .stMarkdown .stSuccess {{
        background: {'linear-gradient(135deg, #1e3a29 0%, #1e293b 100%)' if is_dark_mode else 'linear-gradient(135deg, #e8f5e9 0%, #f1f8f4 100%)'} !important;
        border-left: 4px solid #2e7d32 !important;
        border-radius: 10px !important;
        padding: 1rem !important;
    }}
    .element-container .stMarkdown .stError {{
        background: {'linear-gradient(135deg, #3a1e1e 0%, #1e293b 100%)' if is_dark_mode else 'linear-gradient(135deg, #ffebee 0%, #fff5f5 100%)'} !important;
        border-left: 4px solid #c62828 !important;
        border-radius: 10px !important;
        padding: 1rem !important;
    }}
    
    /* Links styling */
    a {{
        color: {current_theme['secondary']} !important;
        text-decoration: none !important;
        font-weight: 600 !important;
        transition: all 0.2s ease;
    }}
    a:hover {{
        color: {current_theme['primary']} !important;
        text-decoration: underline !important;
    }}
    
    /* Tooltip improvements with advanced popovers */
    [data-testid="stTooltipIcon"] {{
        color: {current_theme['accent']} !important;
        cursor: help;
        transition: all 0.2s ease;
    }}
    [data-testid="stTooltipIcon"]:hover {{
        transform: scale(1.2);
        color: {current_theme['secondary']} !important;
    }}
    
    /* Enhanced hr separator */
    hr {{
        border: none !important;
        height: 2px !important;
        background: linear-gradient(90deg, transparent 0%, {current_theme['bg_card']} 50%, transparent 100%) !important;
        margin: 2rem 0 !important;
    }}
    
    /* ============================================
       STICKY RESULTS PANEL
       ============================================ */
    .sticky-results {{
        position: sticky;
        top: 0;
        z-index: 100;
        background: {'rgba(15, 23, 42, 0.95)' if is_dark_mode else 'rgba(255, 255, 255, 0.95)'};
        backdrop-filter: blur(10px);
        -webkit-backdrop-filter: blur(10px);
        padding: 1rem;
        border-radius: 0 0 12px 12px;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.1);
    }}
    
    /* ============================================
       RESPONSIVE DESIGN FOR TABLETS & LAPTOPS
       ============================================ */
    @media screen and (max-width: 1200px) {{
        .stTabs [data-baseweb="tab"] {{
            font-size: 0.95rem !important;
            padding: 12px 8px !important;
        }}
        h1 {{ font-size: 2rem !important; }}
        h2 {{ font-size: 1.5rem !important; }}
        [data-testid="stMetric"] {{
            padding: 1rem !important;
        }}
        [data-testid="stMetricValue"] {{
            font-size: 1.8rem !important;
        }}
    }}
    
    @media screen and (max-width: 768px) {{
        .stTabs [data-baseweb="tab"] {{
            font-size: 0.85rem !important;
            padding: 10px 6px !important;
        }}
        h1 {{ font-size: 1.6rem !important; }}
        h2 {{ font-size: 1.3rem !important; }}
        h3 {{ font-size: 1.1rem !important; }}
        .stDataFrame, .stTable {{
            font-size: 0.85rem !important;
        }}
        [data-testid="stMetric"] {{
            padding: 0.8rem !important;
        }}
        [data-testid="stMetricValue"] {{
            font-size: 1.5rem !important;
        }}
        .stButton>button {{
            padding: 0.6em 1.2em !important;
            font-size: 0.95rem !important;
        }}
    }}
    
    /* ============================================
       MODERN ACCORDION STYLING
       ============================================ */
    details.streamlit-expander {{
        border-radius: 12px !important;
        overflow: hidden;
        transition: all 0.3s ease;
        margin-bottom: 0.8rem;
    }}
    
    details.streamlit-expander[open] {{
        box-shadow: 0 8px 24px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.15);
    }}
    
    /* ============================================
       HIGHLIGHT SELECTED ITEMS
       ============================================ */
    .highlighted-row {{
        background: linear-gradient(90deg, rgba({rgb['accent'][0]}, {rgb['accent'][1]}, {rgb['accent'][2]}, 0.2) 0%, rgba({rgb['accent'][0]}, {rgb['accent'][1]}, {rgb['accent'][2]}, 0.1) 100%) !important;
        border-left: 4px solid {current_theme['accent']} !important;
    }}
    
    /* ============================================
       PROGRESS PANEL COMPONENTS
       ============================================ */
    .progress-panel {{
        background: linear-gradient(135deg, #1976d2 0%, #42a5f5 100%);
        border-radius: 12px;
        padding: 1.5rem;
        color: white;
        box-shadow: 0 4px 12px rgba(25, 118, 210, 0.3);
        margin-bottom: 1rem;
    }}
    
    .progress-panel--success {{
        background: linear-gradient(135deg, #2e7d32 0%, #4caf50 100%);
        box-shadow: 0 4px 12px rgba(76, 175, 80, 0.3);
    }}
    
    .progress-panel--metrics {{
        background: linear-gradient(135deg, #0d47a1 0%, #1976d2 100%);
        box-shadow: 0 4px 16px rgba(13, 71, 161, 0.25);
        margin-bottom: 2rem;
    }}
    
    .progress-panel--results {{
        background: linear-gradient(135deg, #1976d2 0%, #1565c0 100%);
        border-radius: 16px;
        padding: 2rem;
        margin: 1.5rem 0;
        box-shadow: 0 8px 24px rgba(25, 118, 210, 0.25);
    }}
    
    .progress-panel--hybrid {{
        background: linear-gradient(135deg, #9c27b0 0%, #673ab7 100%);
        border-radius: 16px;
        padding: 2rem;
        margin: 1rem 0;
        box-shadow: 0 8px 24px rgba(156, 39, 176, 0.25);
    }}
    
    .progress-panel__title {{
        margin: 0 0 1rem 0;
        color: white;
        text-align: center;
    }}
    
    .progress-panel__status {{
        margin: 0 0 1rem 0;
        text-align: center;
        opacity: 0.9;
        font-size: 0.95rem;
    }}
    
    .stats-grid {{
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(120px, 1fr));
        gap: 1rem;
        margin-bottom: 1rem;
    }}
    
    .stats-grid--wide {{
        grid-template-columns: repeat(auto-fit, minmax(130px, 1fr));
    }}
    
    .stats-grid--extra-wide {{
        grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
        gap: 1.5rem;
        margin-top: 1rem;
    }}
    
    .stat-card {{
        text-align: center;
        background: rgba(255, 255, 255, 0.15);
        padding: 0.8rem;
        border-radius: 8px;
    }}
    
    .stat-card--large {{
        padding: 1.2rem;
        border-radius: 12px;
        backdrop-filter: blur(10px);
    }}
    
    .stat-card__value {{
        margin: 0;
        color: #FFD700;
        font-size: 1.5rem;
    }}
    
    .stat-card__value--large {{
        font-size: 2.2rem;
        font-weight: 800;
        margin-bottom: 0.5rem;
    }}
    
    .stat-card__label {{
        margin: 0.3rem 0 0 0;
        opacity: 0.9;
        font-size: 0.8rem;
    }}
    
    .stat-card__label--large {{
        font-size: 0.95rem;
        opacity: 0.95;
        font-weight: 500;
    }}
    
    .sequence-info {{
        background: rgba(0, 0, 0, 0.2);
        border-radius: 8px;
        padding: 0.8rem;
        margin-bottom: 1rem;
    }}
    
    .sequence-info__text {{
        margin: 0;
        font-size: 0.9rem;
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
                <h3 class='progress-panel__title' style='font-size: 1.5rem; font-weight: 700; margin-bottom: 1.5rem;'>
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
                        <h3 class='progress-panel__title' style='font-size: 1.5rem; font-weight: 700; margin-bottom: 1.5rem;'>
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
