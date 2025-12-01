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
import streamlit.components.v1 as components
import pandas as pd
import matplotlib.pyplot as plt
import io
import os  # Added for image path checking
import numpy as np
from collections import Counter
import html as html_module  # For escaping user input
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
    get_motif_info as get_motif_classification_info,
    CHUNK_THRESHOLD, DEFAULT_CHUNK_SIZE, DEFAULT_CHUNK_OVERLAP,
    get_last_detector_timings, get_detector_display_names
)
from utilities import export_results_to_dataframe
from visualizations import (
    plot_motif_distribution, plot_coverage_map, plot_density_heatmap,
    plot_length_distribution, plot_score_distribution, plot_nested_pie_chart, 
    save_all_plots, MOTIF_CLASS_COLORS,
    plot_density_comparison, plot_enrichment_analysis, plot_enrichment_summary_table
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

# Try to import scanner backends for parallel processing
try:
    from scanner_backends import get_best_backend, NUMBA_AVAILABLE
    from scanner_backends.parallel_worker import parallel_scan, suggest_num_workers
    SCANNER_BACKENDS_AVAILABLE = True
except ImportError:
    SCANNER_BACKENDS_AVAILABLE = False
    NUMBA_AVAILABLE = False
    parallel_scan = None
    suggest_num_workers = lambda x: 1

# ---------- HTML RENDERING HELPERS ----------
# NOTE: The root cause of HTML being displayed as text (escaped) is that 
# st.markdown() escapes HTML by default. For trusted HTML content, use 
# unsafe_allow_html=True. For untrusted user input, do NOT use unsafe_allow_html 
# to avoid XSS vulnerabilities.

import re

# Regex pattern for validating CSS color values
_CSS_COLOR_PATTERN = re.compile(r'^#[0-9A-Fa-f]{6}$|^#[0-9A-Fa-f]{3}$')

def _sanitize_css_color(color: str, default: str = '#607D8B') -> str:
    """
    Validate and sanitize a CSS color value to prevent CSS injection.
    
    Only accepts valid hex color codes (#RGB or #RRGGBB format).
    Returns a safe default color if the input is invalid.
    
    Args:
        color: Color value to validate
        default: Safe default color to return if validation fails
        
    Returns:
        Sanitized color value (hex format)
    """
    if isinstance(color, str) and _CSS_COLOR_PATTERN.match(color):
        return color
    return default


def render_html(html_content: str) -> None:
    """
    Render trusted HTML content in Streamlit using st.markdown.
    
    This is a convenience wrapper that automatically sets unsafe_allow_html=True.
    Only use this for trusted, internally-generated HTML content.
    Do NOT use this for user-provided content to avoid XSS vulnerabilities.
    
    Args:
        html_content: Trusted HTML string to render
    """
    st.markdown(html_content, unsafe_allow_html=True)


def render_html_component(html_content: str, height: int = 100, scrolling: bool = False) -> None:
    """
    Render HTML using Streamlit's components.html for complex HTML content.
    
    This uses an iframe to render HTML, providing better isolation and support
    for JavaScript. Only use for trusted content.
    
    Args:
        html_content: Trusted HTML string to render
        height: Height of the component in pixels
        scrolling: Whether to allow scrolling in the component
    """
    components.html(html_content, height=height, scrolling=scrolling)


def render_motif_class_badges(motif_classes: list, colors: dict = None) -> str:
    """
    Generate HTML badges/pills for a list of motif classes.
    
    Creates styled inline badges for displaying motif classes in a visually
    appealing format. Returns HTML string to be rendered with render_html().
    
    Args:
        motif_classes: List of motif class names to display
        colors: Optional dict mapping class names to colors. 
                If None, uses MOTIF_CLASS_COLORS from visualizations.
                
    Returns:
        HTML string containing styled badges
    """
    if colors is None:
        colors = MOTIF_CLASS_COLORS
    
    badges = []
    for cls in motif_classes:
        # Escape the class name to prevent XSS from motif class names
        safe_cls = html_module.escape(str(cls))
        # Sanitize the color value to prevent CSS injection
        raw_color = colors.get(cls, '#607D8B')
        safe_color = _sanitize_css_color(raw_color)
        badge_html = f'''<span style="
            display: inline-block;
            padding: 4px 12px;
            margin: 2px 4px;
            background-color: {safe_color};
            color: white;
            border-radius: 16px;
            font-size: 0.9rem;
            font-weight: 600;
            font-family: 'Inter', 'IBM Plex Sans', sans-serif;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        ">{safe_cls}</span>'''
        badges.append(badge_html)
    
    return '<div style="display: flex; flex-wrap: wrap; gap: 4px; margin: 8px 0;">' + ''.join(badges) + '</div>'


def display_motif_classes_styled(motif_classes: list, title: str = None, colors: dict = None) -> None:
    """
    Display a list of motif classes as styled HTML badges in Streamlit.
    
    This is a convenience function that combines render_motif_class_badges
    and render_html for easy use.
    
    Args:
        motif_classes: List of motif class names to display
        title: Optional title to display above the badges
        colors: Optional dict mapping class names to colors
    """
    if title:
        render_html(f'<p style="margin-bottom: 8px; font-weight: 600;">{html_module.escape(title)}</p>')
    
    badges_html = render_motif_class_badges(motif_classes, colors)
    render_html(badges_html)


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


def analyze_sequence_optimized(sequence: str, sequence_name: str = "sequence", 
                               use_parallel: bool = False, progress_bar=None) -> list:
    """
    Analyze sequence with optional parallel processing for large sequences.
    
    For sequences >200KB, uses parallel scanning with shared memory if available.
    Falls back to standard sequential scanning otherwise.
    
    Args:
        sequence: DNA sequence string
        sequence_name: Name for the sequence
        use_parallel: Enable parallel processing (auto-enabled for large sequences)
        progress_bar: Optional Streamlit progress bar to update
        
    Returns:
        List of detected motifs
    """
    seq_len = len(sequence)
    
    # Auto-enable parallel for large sequences if available
    if SCANNER_BACKENDS_AVAILABLE and seq_len > 200_000 and (use_parallel or seq_len > 500_000):
        try:
            # Use parallel scanning with progress callback
            def update_progress(completed, total):
                if progress_bar is not None:
                    progress_bar.progress(completed / total, text=f"Scanning chunks: {completed}/{total}")
            
            num_workers = suggest_num_workers(seq_len)
            motifs = parallel_scan(
                sequence,
                sequence_name=sequence_name,
                num_workers=num_workers,
                chunk_size=100_000,
                progress_callback=update_progress if progress_bar else None
            )
            return motifs
        except Exception as e:
            # Fall back to standard scanning on error
            st.warning(f"Parallel scanning failed, using standard mode: {e}")
    
    # Standard sequential scanning
    return analyze_sequence(sequence, sequence_name)


# ---------- PAGE CONFIG ----------
st.set_page_config(
    page_title="NBDScanner - Non-B DNA Motif Finder",
    layout="wide",
    page_icon="🧬",
    menu_items={'About': "NBDScanner | Developed by Dr. Venkata Rajesh Yella"}
)

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

st.markdown(f"""
    <style>
    /* Google Fonts - optimized for performance with display=swap */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
    
    /* ============================================
       MINIMAL ANIMATIONS FOR PERFORMANCE
       ============================================ */
    @keyframes fade-in {{
        from {{ opacity: 0; }}
        to {{ opacity: 1; }}
    }}
    
    /* ============================================
       MEDIUM BACKGROUND WITH OPTIMAL READABILITY
       ============================================ */
    body, [data-testid="stAppViewContainer"], .main {{
        background: {'linear-gradient(180deg, #1a202c 0%, #2d3748 100%)' if is_dark_mode else 'linear-gradient(180deg, #f0f4f8 0%, #e2e8f0 100%)'} !important;
        color: {'#f7fafc' if is_dark_mode else '#1a202c'} !important;
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif !important;
        font-weight: 500 !important;
        -webkit-font-smoothing: antialiased !important;
        -moz-osx-font-smoothing: grayscale !important;
    }}
    
    /* ============================================
       NAV BAR: Clean, readable tabs with medium background
       ============================================ */
    
    /* Tab container with medium background for readability */
    .stTabs [data-baseweb="tab-list"] {{
        background: {'linear-gradient(90deg, rgba(45,55,72,0.95) 0%, rgba(74,85,104,0.9) 50%)' if is_dark_mode else 'linear-gradient(90deg, rgba(247,250,252,0.98) 0%, rgba(237,242,247,0.95) 50%)'} !important;
        border-bottom: 2px solid {'rgba(160,174,192,0.2)' if is_dark_mode else 'rgba(160,174,192,0.3)'} !important;
        box-shadow: 0 4px 12px rgba(0,0,0,0.08) !important;
        padding: 12px !important;
        border-radius: 12px !important;
        margin: 0 auto 1.5rem auto !important;
        max-width: 1200px !important;
    }}
    
    /* Tab pills with bold text for readability */
    .stTabs [data-baseweb="tab"] {{
        display: flex !important;
        align-items: center !important;
        justify-content: center !important;
        min-width: 140px !important;
        padding: 12px 20px !important;
        margin: 0 6px !important;
        background: {'linear-gradient(180deg, rgba(74,85,104,0.95), rgba(45,55,72,0.9))' if is_dark_mode else 'linear-gradient(180deg, rgba(255,255,255,0.98), rgba(247,250,252,0.95))'} !important;
        color: {'#f7fafc' if is_dark_mode else '#2d3748'} !important;
        border-radius: 10px !important;
        border: 1px solid {'rgba(160,174,192,0.25)' if is_dark_mode else 'rgba(160,174,192,0.35)'} !important;
        box-shadow: 0 2px 6px rgba(0,0,0,0.06) !important;
        font-weight: 600 !important;
        font-size: 0.95rem !important;
        letter-spacing: 0.01em !important;
        transition: all 0.25s ease !important;
    }}
    
    /* Ensure tab text is properly aligned */
    .stTabs [data-baseweb="tab"] > * {{
        display: inline-flex !important;
        align-items: center !important;
        justify-content: center !important;
        gap: 6px !important;
    }}
    
    /* Selected tab: clear visual distinction with medium background */
    .stTabs [aria-selected="true"] {{
        background: {'linear-gradient(180deg, #4a5568 0%, #2d3748 100%)' if is_dark_mode else 'linear-gradient(180deg, #ebf8ff 0%, #bee3f8 100%)'} !important;
        color: {'#fff' if is_dark_mode else '#2c5282'} !important;
        border: 1px solid {'rgba(99,179,237,0.5)' if is_dark_mode else 'rgba(66,153,225,0.4)'} !important;
        box-shadow: 0 4px 12px {'rgba(99,179,237,0.2)' if is_dark_mode else 'rgba(66,153,225,0.15)'} !important;
        font-weight: 700 !important;
        transform: translateY(-2px) !important;
    }}
    
    /* Hover state: subtle feedback */
    .stTabs [data-baseweb="tab"]:hover {{
        transform: translateY(-1px) !important;
        box-shadow: 0 4px 10px rgba(0,0,0,0.1) !important;
        background: {'linear-gradient(180deg, rgba(74,85,104,1), rgba(45,55,72,0.95))' if is_dark_mode else 'linear-gradient(180deg, rgba(255,255,255,1), rgba(237,242,247,0.98))'} !important;
    }}
    
    /* Hide unnecessary pseudo-elements */
    .stTabs [data-baseweb="tab"]::after,
    .stTabs [data-baseweb="tab"][aria-selected="true"]::after {{
        display: none !important;
    }}
    
    /* Responsive adjustments */
    @media screen and (max-width: 900px) {{
        .stTabs [data-baseweb="tab"] {{ min-width: 100px !important; padding: 10px 14px !important; font-size: 0.9rem !important; }}
        .stTabs [data-baseweb="tab-list"] {{ padding: 10px !important; }}
    }}
    
    /* Tab panel content */
    .stTabs [data-baseweb="tab-panel"] {{
        padding-top: 1.5rem !important;
    }}
    
    /* ============================================
       HEADINGS: Bold, readable typography
       ============================================ */
    h1, h2, h3, h4 {{
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif !important;
        color: {current_theme['primary']} !important;
        font-weight: 700 !important;
        letter-spacing: -0.01em !important;
        margin-top: 1.2em !important;
        margin-bottom: 0.6em !important;
        line-height: 1.3 !important;
    }}
    h1 {{ 
        font-size: 2.25rem !important; 
        font-weight: 800 !important;
        background: linear-gradient(135deg, {current_theme['primary']} 0%, {current_theme['secondary']} 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
    }}
    h2 {{ 
        font-size: 1.65rem !important; 
        color: {current_theme['primary']} !important; 
        font-weight: 700 !important;
        border-bottom: 2px solid {'rgba(160,174,192,0.3)' if is_dark_mode else 'rgba(160,174,192,0.4)'};
        padding-bottom: 0.4rem;
    }}
    h3 {{ 
        font-size: 1.35rem !important; 
        color: {current_theme['secondary']} !important; 
        font-weight: 700 !important;
    }}
    h4 {{ 
        font-size: 1.15rem !important; 
        color: {current_theme['secondary']} !important; 
        font-weight: 600 !important;
    }}
    
    /* Body text: Optimal readability with medium weight */
    .stMarkdown, .markdown-text-container, .stText, p, span, label {{
        font-size: 1rem !important;
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif !important;
        line-height: 1.65 !important;
        color: {current_theme['text']} !important;
        margin-bottom: 0.6em !important;
        font-weight: 500 !important;
    }}
    
    /* Input fields: Clean design with medium background */
    input, .stTextInput>div>div>input, .stSelectbox>div>div>div, 
    .stMultiSelect>div>div>div, .stRadio>div>div>label>div {{
        font-size: 0.95rem !important;
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif !important;
        border-radius: 8px !important;
        border: 1.5px solid {'#4a5568' if is_dark_mode else '#cbd5e0'} !important;
        background: {'#2d3748' if is_dark_mode else '#fff'} !important;
        transition: border-color 0.2s ease, box-shadow 0.2s ease;
        font-weight: 500 !important;
    }}
    input:focus, .stTextInput>div>div>input:focus {{
        border-color: {current_theme['secondary']} !important;
        box-shadow: 0 0 0 3px rgba({rgb['secondary'][0]}, {rgb['secondary'][1]}, {rgb['secondary'][2]}, 0.15) !important;
        outline: none !important;
    }}
    
    /* ============================================
       BUTTONS: Bold text, medium background
       ============================================ */
    .stButton>button {{
        font-size: 1rem !important;
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif !important;
        padding: 0.7em 1.6em !important;
        background: linear-gradient(135deg, {current_theme['primary']} 0%, {current_theme['secondary']} 100%) !important;
        color: #fff !important;
        border-radius: 10px !important;
        border: none !important;
        font-weight: 600 !important;
        box-shadow: 0 4px 12px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.25);
        transition: transform 0.2s ease, box-shadow 0.2s ease;
        cursor: pointer;
    }}
    .stButton>button:hover {{
        transform: translateY(-2px) !important;
        box-shadow: 0 6px 16px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.35);
    }}
    .stButton>button:active {{
        transform: translateY(0) !important;
        box-shadow: 0 2px 8px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.25);
    }}
    
    /* ============================================
       DATAFRAMES: Readable with medium backgrounds
       ============================================ */
    .stDataFrame, .stTable {{
        font-size: {'0.875rem' if is_compact else '0.925rem'} !important;
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif !important;
        line-height: {'1.45' if is_compact else '1.55'} !important;
        border-radius: 10px !important;
        overflow: hidden;
        box-shadow: 0 2px 8px rgba(0,0,0,0.08) !important;
        border: 1px solid {'#4a5568' if is_dark_mode else '#e2e8f0'} !important;
    }}
    .stDataFrame thead tr th {{
        background: {'#2d3748' if is_dark_mode else f"{current_theme['primary']}"} !important;
        color: #fff !important;
        font-weight: 600 !important;
        padding: {'0.65rem' if is_compact else '0.9rem'} !important;
        font-size: {'0.85rem' if is_compact else '0.9rem'} !important;
        letter-spacing: 0.01em;
        text-transform: none;
        border-bottom: 2px solid {'#4a5568' if is_dark_mode else current_theme['secondary']} !important;
    }}
    .stDataFrame tbody tr {{
        transition: background-color 0.15s ease;
        border-bottom: 1px solid {'#4a5568' if is_dark_mode else '#e2e8f0'} !important;
    }}
    .stDataFrame tbody tr:hover {{
        background: {'#374151' if is_dark_mode else '#f7fafc'} !important;
    }}
    .stDataFrame tbody tr td {{
        padding: {'0.5rem' if is_compact else '0.85rem'} !important;
        color: {current_theme['text']} !important;
    }}
    
    /* Highlight alternate rows */
    .stDataFrame tbody tr:nth-child(odd) {{
        background: {'rgba(55, 65, 81, 0.4)' if is_dark_mode else 'rgba(247, 250, 252, 0.6)'} !important;
    }}
    .stDataFrame tbody td {{
        padding: {'0.55rem' if is_compact else '0.75rem'} !important;
        font-weight: 500 !important;
        color: {current_theme['text']} !important;
    }}
    
    /* Tab content spacing - streamlined */
    .stTabs [data-baseweb="tab-panel"] {{
        padding: 1.5rem 1rem 1.5rem 1rem !important;
    }}
    
    /* ============================================
       CARDS: Clean medium background
       ============================================ */
    .analysis-summary-card, .glassmorphism-card {{
        margin: 1.25rem 0 !important;
        padding: 1.5rem !important;
        background: {'rgba(45, 55, 72, 0.9)' if is_dark_mode else 'rgba(255, 255, 255, 0.95)'} !important;
        border-radius: 12px !important;
        box-shadow: 0 4px 16px rgba(0,0,0,0.08) !important;
        border: 1px solid {'rgba(74, 85, 104, 0.4)' if is_dark_mode else 'rgba(226, 232, 240, 0.8)'} !important;
    }}
    .analysis-summary-card:hover, .glassmorphism-card:hover {{
        box-shadow: 0 6px 20px rgba(0,0,0,0.12) !important;
    }}
    
    /* Select boxes and multiselect: Clean design */
    .stSelectbox > div > div > div, .stMultiSelect > div > div > div {{
        min-height: 2.75rem !important;
        padding: 0.6rem 1rem !important;
        line-height: 1.5 !important;
        border-radius: 8px !important;
        background: {'#2d3748' if is_dark_mode else '#fff'} !important;
        border: 1.5px solid {'#4a5568' if is_dark_mode else '#cbd5e0'} !important;
        font-weight: 500 !important;
        transition: border-color 0.2s ease;
    }}
    .stSelectbox > div > div > div:hover, .stMultiSelect > div > div > div:hover {{
        border-color: {current_theme['secondary']} !important;
    }}
    
    /* Text input and textarea: Clean styling */
    .stTextInput > div > div > input, .stTextArea textarea {{
        padding: 0.7rem 1rem !important;
        min-height: 2.75rem !important;
        line-height: 1.5 !important;
        border-radius: 8px !important;
        font-size: 0.95rem !important;
        font-weight: 500 !important;
        background: {'#2d3748' if is_dark_mode else '#fff'} !important;
    }}
    .stTextArea textarea {{
        min-height: 100px !important;
    }}
    
    /* Number input styling */
    .stNumberInput > div > div > input {{
        padding: 0.7rem 1rem !important;
        min-height: 2.75rem !important;
        border-radius: 8px !important;
        font-size: 0.95rem !important;
        font-weight: 600 !important;
        background: {'#2d3748' if is_dark_mode else '#fff'} !important;
    }}
    
    /* Radio buttons: Clean with bold selected state */
    .stRadio > div {{
        gap: 0.35rem !important;
    }}
    .stRadio > div > div > label {{
        margin-bottom: 0.5rem !important;
        padding: 0.55rem 0.9rem !important;
        border-radius: 8px !important;
        transition: background-color 0.2s ease;
        border: 1.5px solid transparent;
        background: {'#374151' if is_dark_mode else '#f7fafc'};
        cursor: pointer;
        font-weight: 500 !important;
    }}
    .stRadio > div > div > label:hover {{
        background: {'#4a5568' if is_dark_mode else '#edf2f7'} !important;
    }}
    .stRadio > div > div > label[data-checked="true"] {{
        background: {'#4a5568' if is_dark_mode else current_theme['bg_card']} !important;
        border-color: {current_theme['secondary']} !important;
        font-weight: 600 !important;
    }}
    /* Radio button circles */
    .stRadio > div > div > label > div:first-child {{
        border-width: 2px !important;
        border-color: {'#718096' if is_dark_mode else '#a0aec0'} !important;
        width: 18px !important;
        height: 18px !important;
    }}
    .stRadio > div > div > label[data-checked="true"] > div:first-child {{
        border-color: {current_theme['primary']} !important;
        background-color: {current_theme['primary']} !important;
    }}
    
    /* Info boxes and alerts: Clean, readable design */
    .stAlert {{
        border-radius: 8px !important;
        border-left: 4px solid {current_theme['secondary']} !important;
        box-shadow: 0 2px 6px rgba(0,0,0,0.06) !important;
        background: {'#374151' if is_dark_mode else '#f7fafc'} !important;
        padding: 1rem !important;
        font-weight: 500 !important;
    }}
    .stAlert[data-baseweb="notification"] {{
        border-left-color: {current_theme['primary']} !important;
    }}
    .stSuccess {{
        border-left-color: #38a169 !important;
        background: {'#1c4532' if is_dark_mode else '#f0fff4'} !important;
    }}
    .stWarning {{
        border-left-color: #dd6b20 !important;
        background: {'#452c19' if is_dark_mode else '#fffaf0'} !important;
    }}
    .stError {{
        border-left-color: #e53e3e !important;
        background: {'#452424' if is_dark_mode else '#fff5f5'} !important;
    }}
    
    /* ============================================
       PROGRESS BARS: Clean, performant design
       ============================================ */
    .stProgress {{
        margin: 1rem 0 !important;
    }}
    .stProgress > div > div {{
        background: linear-gradient(90deg, {current_theme['primary']} 0%, {current_theme['secondary']} 100%) !important;
        border-radius: 8px !important;
        box-shadow: 0 2px 4px rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.2) !important;
        min-height: 10px !important;
        transition: width 0.3s ease;
    }}
    .stProgress > div {{
        border-radius: 8px !important;
        background: {'#4a5568' if is_dark_mode else '#e2e8f0'} !important;
        min-height: 12px !important;
        overflow: hidden;
    }}
    
    /* Progress bar text */
    .stProgress + div p {{
        font-size: 0.95rem !important;
        font-weight: 600 !important;
        color: {current_theme['secondary']} !important;
        margin-top: 0.5rem !important;
    }}
    
    /* File uploader: Clean design */
    [data-testid="stFileUploader"] {{
        border: 2px dashed {'#4a5568' if is_dark_mode else '#a0aec0'} !important;
        border-radius: 10px !important;
        background: {'#2d3748' if is_dark_mode else '#f7fafc'} !important;
        padding: 1.5rem !important;
        transition: border-color 0.2s ease;
    }}
    [data-testid="stFileUploader"]:hover {{
        border-color: {current_theme['secondary']} !important;
    }}
    [data-testid="stFileUploader"] section {{
        border: none !important;
    }}
    [data-testid="stFileUploader"] button {{
        background: {current_theme['secondary']} !important;
        color: white !important;
        border: none !important;
        border-radius: 8px !important;
        padding: 0.5rem 1.2rem !important;
        font-weight: 600 !important;
    }}
    [data-testid="stFileUploader"] button:hover {{
        background: {current_theme['primary']} !important;
    }}
    
    /* Checkboxes: Clean design */
    .stCheckbox {{
        padding: 0.35rem 0 !important;
    }}
    .stCheckbox > label {{
        padding: 0.5rem 0.7rem !important;
        border-radius: 6px !important;
        transition: background-color 0.2s ease;
        cursor: pointer;
        background: {'#374151' if is_dark_mode else '#f7fafc'};
        font-weight: 500 !important;
    }}
    .stCheckbox > label:hover {{
        background: {'#4a5568' if is_dark_mode else '#edf2f7'} !important;
    }}
    .stCheckbox > label > div:first-child {{
        border-width: 2px !important;
        border-color: {'#718096' if is_dark_mode else '#a0aec0'} !important;
        border-radius: 4px !important;
        width: 18px !important;
        height: 18px !important;
    }}
    .stCheckbox > label > div:first-child[data-checked="true"] {{
        background-color: {current_theme['primary']} !important;
        border-color: {current_theme['primary']} !important;
    }}
    
    /* Expander styling: Clean accordion design */
    .streamlit-expanderHeader {{
        border-radius: 8px !important;
        background: {'#374151' if is_dark_mode else '#f7fafc'} !important;
        font-weight: 600 !important;
        padding: 0.85rem 1.2rem !important;
        border: 1.5px solid {'#4a5568' if is_dark_mode else '#e2e8f0'} !important;
        transition: background-color 0.2s ease;
        display: flex !important;
        align-items: center !important;
        gap: 10px !important;
    }}
    .streamlit-expanderHeader:hover {{
        background: {'#4a5568' if is_dark_mode else '#edf2f7'} !important;
    }}
    .streamlit-expanderContent {{
        border-radius: 0 0 8px 8px !important;
        background: {'#2d3748' if is_dark_mode else '#fff'} !important;
        border: 1.5px solid {'#4a5568' if is_dark_mode else '#e2e8f0'} !important;
        border-top: none !important;
        padding: 1.2rem !important;
    }}
    
    /* Chevron icon for expanders */
    .streamlit-expanderHeader svg {{
        margin-right: 8px !important;
        flex-shrink: 0 !important;
        width: 18px !important;
        height: 18px !important;
        color: {current_theme['primary']} !important;
        transition: transform 0.2s ease !important;
    }}
    
    /* Rotate chevron when expanded */
    details[open] > .streamlit-expanderHeader svg {{
        transform: rotate(90deg) !important;
    }}
    
    /* Hide raw icon text fallback */
    .streamlit-expanderHeader::before {{
        display: none !important;
    }}
    
    /* Expander header layout */
    .streamlit-expanderHeader {{
        display: flex !important;
        align-items: center !important;
        line-height: 1.4 !important;
    }}
    
    /* Expander text styling */
    .streamlit-expanderHeader p,
    .streamlit-expanderHeader label,
    .streamlit-expanderHeader div,
    .streamlit-expanderHeader span {{
        margin-left: 0 !important;
        flex: 1 !important;
        line-height: 1.4 !important;
        font-weight: 600 !important;
    }}
    
    /* Hide keyboard_arrow_down text fallback */
    .streamlit-expanderHeader span:first-child {{
        display: flex !important;
        align-items: center !important;
        flex-shrink: 0 !important;
        min-width: 20px !important;
        max-width: 20px !important;
        overflow: hidden !important;
    }}
    
    /* Icon and label spacing */
    [role="button"] svg,
    [role="tab"] svg,
    button svg,
    label svg {{
        margin-right: 5px !important;
        flex-shrink: 0 !important;
    }}
    
    /* Interactive elements layout */
    [role="button"],
    [role="tab"],
    button,
    label {{
        display: flex !important;
        align-items: center !important;
        gap: 5px !important;
        line-height: 1.4 !important;
    }}
    
    /* Icon container styling */
    [data-testid="StyledIconContainer"] {{
        display: flex !important;
        align-items: center !important;
        justify-content: center !important;
        width: 20px !important;
        min-width: 20px !important;
        overflow: hidden !important;
    }}
    
    /* ============================================
       METRIC CARDS: Clean, readable design
       ============================================ */
    [data-testid="stMetric"] {{
        background: {'#2d3748' if is_dark_mode else '#fff'} !important;
        padding: 1.2rem !important;
        border-radius: 10px !important;
        border: 1.5px solid {'#4a5568' if is_dark_mode else '#e2e8f0'} !important;
        box-shadow: 0 2px 6px rgba(0,0,0,0.06) !important;
    }}
    [data-testid="stMetric"]:hover {{
        box-shadow: 0 4px 12px rgba(0,0,0,0.1) !important;
    }}
    [data-testid="stMetricValue"] {{
        font-size: 1.85rem !important;
        font-weight: 700 !important;
        color: {current_theme['primary']} !important;
    }}
    [data-testid="stMetricLabel"] {{
        font-size: 0.875rem !important;
        font-weight: 600 !important;
        color: {'#a0aec0' if is_dark_mode else '#4a5568'} !important;
        text-transform: uppercase;
        letter-spacing: 0.03em;
    }}
    
    /* Spinner */
    .stSpinner > div {{
        border-top-color: {current_theme['secondary']} !important;
    }}
    
    /* Code blocks: Clean design */
    .stCodeBlock {{
        border-radius: 8px !important;
        border: 1.5px solid {'#4a5568' if is_dark_mode else '#e2e8f0'} !important;
        background: {'#1a202c' if is_dark_mode else '#f7fafc'} !important;
        padding: 0.9rem !important;
    }}
    code {{
        font-family: 'JetBrains Mono', 'Fira Code', 'Consolas', monospace !important;
        font-size: 0.875rem !important;
        color: {current_theme['primary']} !important;
        background: {'rgba(45, 55, 72, 0.5)' if is_dark_mode else 'rgba(226, 232, 240, 0.6)'} !important;
        padding: 0.15rem 0.4rem !important;
        border-radius: 4px !important;
        font-weight: 500 !important;
    }}
    
    /* Sidebar styling */
    [data-testid="stSidebar"] {{
        background: {'#1a202c' if is_dark_mode else '#f7fafc'} !important;
        border-right: 1.5px solid {'#4a5568' if is_dark_mode else '#e2e8f0'} !important;
    }}
    [data-testid="stSidebar"] .stMarkdown {{
        color: {current_theme['text']} !important;
    }}
    
    /* Scrollbar: Simple design */
    ::-webkit-scrollbar {{
        width: 10px;
        height: 10px;
    }}
    ::-webkit-scrollbar-track {{
        background: {'#2d3748' if is_dark_mode else '#edf2f7'};
        border-radius: 8px;
    }}
    ::-webkit-scrollbar-thumb {{
        background: {'#4a5568' if is_dark_mode else '#a0aec0'};
        border-radius: 8px;
    }}
    ::-webkit-scrollbar-thumb:hover {{
        background: {'#718096' if is_dark_mode else '#718096'};
    }}
    
    /* Download buttons */
    [data-testid="stDownloadButton"] button {{
        background: {current_theme['primary']} !important;
        color: white !important;
        border-radius: 8px !important;
        padding: 0.6rem 1.4rem !important;
        font-weight: 600 !important;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }}
    [data-testid="stDownloadButton"] button:hover {{
        background: {current_theme['secondary']} !important;
    }}
    
    /* Caption text */
    .caption, [data-testid="stCaptionContainer"] {{
        color: {'#a0aec0' if is_dark_mode else '#718096'} !important;
        font-size: 0.85rem !important;
        font-style: italic;
        font-weight: 500 !important;
    }}
    
    /* Success/Error message boxes */
    .element-container .stMarkdown .stSuccess {{
        background: {'#1c4532' if is_dark_mode else '#f0fff4'} !important;
        border-left: 4px solid #38a169 !important;
        border-radius: 8px !important;
        padding: 0.9rem !important;
    }}
    .element-container .stMarkdown .stError {{
        background: {'#452424' if is_dark_mode else '#fff5f5'} !important;
        border-left: 4px solid #e53e3e !important;
        border-radius: 8px !important;
        padding: 0.9rem !important;
    }}
    
    /* Links */
    a {{
        color: {current_theme['secondary']} !important;
        text-decoration: none !important;
        font-weight: 600 !important;
    }}
    a:hover {{
        color: {current_theme['primary']} !important;
        text-decoration: underline !important;
    }}
    
    /* Tooltip */
    [data-testid="stTooltipIcon"] {{
        color: {current_theme['accent']} !important;
        cursor: help;
    }}
    [data-testid="stTooltipIcon"]:hover {{
        color: {current_theme['secondary']} !important;
    }}
    
    /* Separator */
    hr {{
        border: none !important;
        height: 1.5px !important;
        background: {'#4a5568' if is_dark_mode else '#e2e8f0'} !important;
        margin: 1.5rem 0 !important;
    }}
    
    /* Sticky results panel */
    .sticky-results {{
        position: sticky;
        top: 0;
        z-index: 100;
        background: {'rgba(26, 32, 44, 0.95)' if is_dark_mode else 'rgba(255, 255, 255, 0.95)'};
        padding: 0.9rem;
        border-radius: 0 0 10px 10px;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.08);
    }}
    
    /* ============================================
       RESPONSIVE DESIGN
       ============================================ */
    @media screen and (max-width: 1200px) {{
        .stTabs [data-baseweb="tab"] {{
            font-size: 0.875rem !important;
            padding: 10px 14px !important;
        }}
        .stTabs [data-baseweb="tab-list"] {{
            border-radius: 10px !important;
            gap: 5px !important;
        }}
        h1 {{ font-size: 1.85rem !important; }}
        h2 {{ font-size: 1.4rem !important; }}
        [data-testid="stMetric"] {{
            padding: 1rem !important;
        }}
        [data-testid="stMetricValue"] {{
            font-size: 1.6rem !important;
        }}
    }}
    
    @media screen and (max-width: 768px) {{
        .stTabs [data-baseweb="tab-list"] {{
            flex-wrap: wrap !important;
            border-radius: 8px !important;
            gap: 4px !important;
            padding: 6px !important;
        }}
        .stTabs [data-baseweb="tab"] {{
            font-size: 0.8rem !important;
            padding: 8px 10px !important;
            flex: 1 1 45% !important;
            min-width: 45% !important;
        }}
        h1 {{ font-size: 1.5rem !important; }}
        h2 {{ font-size: 1.25rem !important; }}
        h3 {{ font-size: 1.1rem !important; }}
        .stDataFrame, .stTable {{
            font-size: 0.8rem !important;
        }}
        [data-testid="stMetric"] {{
            padding: 0.8rem !important;
        }}
        [data-testid="stMetricValue"] {{
            font-size: 1.4rem !important;
        }}
        .stButton>button {{
            padding: 0.5em 1em !important;
            font-size: 0.9rem !important;
        }}
    }}
    
    /* Accordion styling */
    details.streamlit-expander {{
        border-radius: 8px !important;
        overflow: hidden;
        margin-bottom: 0.6rem;
    }}
    
    details.streamlit-expander[open] {{
        box-shadow: 0 4px 12px rgba(0,0,0,0.08);
    }}
    
    /* Highlighted rows */
    .highlighted-row {{
        background: rgba({rgb['accent'][0]}, {rgb['accent'][1]}, {rgb['accent'][2]}, 0.15) !important;
        border-left: 3px solid {current_theme['accent']} !important;
    }}
    
    /* Print styles */
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
    "Disease Analysis": "Disease-Associated Non-B DNA",
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
                        st.markdown(f"**{html_module.escape(names[i])}**: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
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
                        st.markdown(f"**{html_module.escape(names[i])}**: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
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
                        st.markdown(f"**{html_module.escape(names[i])}**: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
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
            with st.popover("Motif Example Queries"):
                for motif, example in motif_examples.items():
                    st.write(f"**{motif}**: `{example}`")
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
                            st.markdown(f"<b>{html_module.escape(names[i])}</b>: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
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
                st.markdown(f"**{html_module.escape(st.session_state.names[i])}** ({len(seq):,} bp) | GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}", unsafe_allow_html=True)
                st.code(wrap(seq[:400]), language="fasta")
            if len(st.session_state.seqs) > 2:
                st.caption(f"...and {len(st.session_state.seqs)-2} more.")

    # ----- RIGHT COLUMN: Analysis Controls + Run Button + Summary Table -----
    with col_right:
        # Elegant Analysis Control Panel
        st.markdown("""
        <div style='background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                    border-radius: 20px; padding: 2rem; margin-bottom: 1.5rem;
                    box-shadow: 0 10px 40px rgba(102, 126, 234, 0.3);
                    position: relative; overflow: hidden;'>
            <div style='position: absolute; top: -30px; right: -30px; width: 120px; height: 120px;
                        background: radial-gradient(circle, rgba(255,255,255,0.15) 0%, transparent 70%);
                        border-radius: 50%;'></div>
            <div style='position: relative; z-index: 1;'>
                <h2 style='margin: 0 0 0.5rem 0; color: white; font-size: 1.6rem; font-weight: 800;
                           display: flex; align-items: center; gap: 0.5rem;'>
                    🚀 Analysis & Run
                </h2>
                <p style='color: rgba(255,255,255,0.9); font-size: 0.95rem; margin: 0;'>
                    Configure and execute Non-B DNA motif detection
                </p>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        # Analysis Options Card
        st.markdown("""
        <div style='background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
                    border-radius: 16px; padding: 1.5rem; margin-bottom: 1.5rem;
                    border: 1px solid rgba(0,0,0,0.08);
                    box-shadow: 0 4px 15px rgba(0,0,0,0.05);'>
            <h4 style='margin: 0 0 1rem 0; color: #1a1a2e; display: flex; align-items: center; gap: 0.5rem;'>
                <span style='background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                             -webkit-background-clip: text; -webkit-text-fill-color: transparent;
                             font-size: 1.2rem;'>⚙️</span>
                Analysis Configuration
            </h4>
        </div>
        """, unsafe_allow_html=True)
        
        # Detection info with elegant badge
        st.markdown("""
        <div style='background: linear-gradient(135deg, #e0f7fa 0%, #b2ebf2 100%);
                    border-radius: 12px; padding: 1rem; margin-bottom: 1rem;
                    border-left: 4px solid #00bcd4;
                    display: flex; align-items: center; gap: 0.75rem;'>
            <span style='font-size: 1.5rem;'>🔬</span>
            <div>
                <p style='margin: 0; color: #006064; font-weight: 700; font-size: 0.95rem;'>
                    Comprehensive Detection Enabled
                </p>
                <p style='margin: 0.2rem 0 0 0; color: #00838f; font-size: 0.85rem;'>
                    11 motif classes • 22+ subclasses • Parallel processing
                </p>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        # Simple options with modern toggle-like appearance
        col1, col2 = st.columns(2)
        with col1:
            detailed_output = st.checkbox("📊 Detailed Analysis", value=True, 
                                        help="Include comprehensive motif metadata")
        with col2:
            quality_check = st.checkbox("✓ Quality Validation", value=True, 
                                      help="Validate detected motifs")
        
        # Chunking and parallel programming are now enabled by default
        # to improve performance for all users without requiring configuration.
        # Large sequences (>10kb) automatically use chunking, and progress is always shown.
        show_chunk_progress = True
        use_parallel_scanner = True
        
        # Hardcoded default overlap handling: always remove overlaps within subclasses
        nonoverlap = True
        overlap_option = "Remove overlaps within subclasses"
        
        
        # ========== MODERN NBDSCANNER ANALYSIS STYLES ==========
        st.markdown("""
        <style>
        .nbdscanner-analysis-block {
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
            border-radius: 20px;
            padding: 2rem;
            margin: 1.5rem 0;
            box-shadow: 0 15px 50px rgba(15, 52, 96, 0.4),
                        0 0 0 1px rgba(255, 255, 255, 0.05);
            position: relative;
            overflow: hidden;
        }
        .nbdscanner-analysis-block::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 3px;
            background: linear-gradient(90deg, #667eea, #764ba2, #f64f59, #c471ed, #12c2e9);
            background-size: 300% 100%;
            animation: gradient-flow 4s ease infinite;
        }
        @keyframes gradient-flow {
            0%, 100% { background-position: 0% 50%; }
            50% { background-position: 100% 50%; }
        }
        .nbdscanner-title {
            font-family: 'Inter', 'IBM Plex Sans', sans-serif;
            font-size: 1.5rem;
            font-weight: 800;
            background: linear-gradient(135deg, #ffffff 0%, #e0e0e0 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            text-align: center;
            margin-bottom: 1.5rem;
            letter-spacing: 0.02em;
        }
        .progress-container {
            background: rgba(255, 255, 255, 0.03);
            border-radius: 16px;
            padding: 1.5rem;
            margin: 1.5rem 0;
            border: 1px solid rgba(255, 255, 255, 0.08);
            backdrop-filter: blur(10px);
        }
        .progress-bar-outer {
            background: rgba(0, 0, 0, 0.4);
            border-radius: 10px;
            height: 24px;
            overflow: hidden;
            margin: 1rem 0;
            box-shadow: inset 0 2px 6px rgba(0, 0, 0, 0.4);
            position: relative;
        }
        .progress-bar-inner {
            height: 100%;
            background: linear-gradient(90deg, 
                #667eea 0%, #764ba2 25%, #f64f59 50%, #c471ed 75%, #12c2e9 100%);
            background-size: 300% 100%;
            border-radius: 10px;
            transition: width 0.5s cubic-bezier(0.4, 0, 0.2, 1);
            box-shadow: 0 0 30px rgba(102, 126, 234, 0.5);
            animation: shimmer-gradient 3s linear infinite;
            position: relative;
        }
        @keyframes shimmer-gradient {
            0% { background-position: 100% 50%; }
            100% { background-position: 0% 50%; }
        }
        .progress-bar-inner::after {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            bottom: 0;
            background: linear-gradient(90deg, transparent, rgba(255,255,255,0.3), transparent);
            animation: shine-sweep 2s ease-in-out infinite;
        }
        @keyframes shine-sweep {
            0% { transform: translateX(-100%); }
            100% { transform: translateX(100%); }
        }
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(4, 1fr);
            gap: 1rem;
            margin-top: 1.5rem;
        }
        .stat-item {
            background: linear-gradient(135deg, rgba(255, 255, 255, 0.08) 0%, rgba(255, 255, 255, 0.02) 100%);
            border-radius: 16px;
            padding: 1.2rem;
            text-align: center;
            border: 1px solid rgba(255, 255, 255, 0.08);
            transition: all 0.3s ease;
        }
        .stat-item:hover {
            transform: translateY(-3px);
            border-color: rgba(102, 126, 234, 0.3);
            box-shadow: 0 8px 25px rgba(102, 126, 234, 0.15);
        }
        .stat-value {
            font-family: 'Inter', sans-serif;
            font-size: 1.6rem;
            font-weight: 800;
            background: linear-gradient(135deg, #12c2e9 0%, #c471ed 50%, #f64f59 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            margin: 0;
        }
        .stat-label {
            font-size: 0.75rem;
            color: rgba(255, 255, 255, 0.7);
            margin: 0.4rem 0 0 0;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.1em;
        }
        .step-indicator {
            background: linear-gradient(90deg, rgba(102, 126, 234, 0.15) 0%, transparent 100%);
            border-left: 4px solid;
            border-image: linear-gradient(180deg, #667eea, #764ba2) 1;
            border-radius: 0 12px 12px 0;
            padding: 1rem 1.5rem;
            margin: 1rem 0;
        }
        .step-text {
            font-family: 'Inter', sans-serif;
            font-weight: 700;
            color: rgba(255, 255, 255, 0.95);
            font-size: 1rem;
            letter-spacing: 0.02em;
        }
        .completed-badge {
            background: linear-gradient(135deg, #00b09b 0%, #96c93d 100%);
            border-radius: 16px;
            padding: 1.2rem 2rem;
            text-align: center;
            margin-top: 1.5rem;
            box-shadow: 0 8px 30px rgba(0, 176, 155, 0.35);
        }
        .completed-text {
            font-family: 'Inter', sans-serif;
            font-size: 1.2rem;
            font-weight: 800;
            color: #ffffff;
            text-shadow: 0 1px 2px rgba(0, 0, 0, 0.2);
        }
        </style>
        """, unsafe_allow_html=True)
        
        if st.button("🔬 Run NBDScanner Analysis", type="primary", use_container_width=True, key="run_motif_analysis_main"):
            # Simplified validation
            if not st.session_state.seqs:
                st.error("❌ Please upload or input sequences before running analysis.")
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
                
                # Create placeholder for combined progress block
                combined_progress_placeholder = st.empty()
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
                
                # Estimate processing time based on sequence length
                def estimate_time(total_bp):
                    return total_bp / ESTIMATED_BP_PER_SECOND
                
                # Helper function to render the combined progress block
                def render_progress_block(progress_pct, current_step, elapsed_sec, motifs_found, speed_bps, is_complete=False, detector_timings=None):
                    """Render the combined NBDScanner Run Button + Progress + Stats Block"""
                    # Calculate progress bar blocks (using Unicode block characters)
                    filled_blocks = int(progress_pct / 12.5)  # 8 total blocks
                    empty_blocks = 8 - filled_blocks
                    progress_visual = "█" * filled_blocks + "▒" * empty_blocks
                    
                    # Build detector timings HTML if available
                    detector_timings_html = ""
                    if is_complete and detector_timings:
                        detector_display_names = get_detector_display_names()
                        # Sort by detector name for consistent display order
                        sorted_detectors = sorted(detector_timings.items(), key=lambda x: x[0])
                        
                        detector_timings_html = """
                        <div style='margin-top: 1.5rem; padding: 1rem; background: rgba(255,255,255,0.05); 
                                    border-radius: 12px; border: 1px solid rgba(255,255,255,0.1);'>
                            <div style='font-size: 0.85rem; font-weight: 700; color: rgba(255,255,255,0.9); 
                                        margin-bottom: 0.8rem; text-transform: uppercase; letter-spacing: 0.1em;'>
                                ⏱️ Detector Timings
                            </div>
                            <div style='display: grid; grid-template-columns: repeat(3, 1fr); gap: 0.5rem;'>
                        """
                        for det_name, det_time in sorted_detectors:
                            display_name = detector_display_names.get(det_name, det_name)
                            detector_timings_html += f"""
                                <div style='background: rgba(102, 126, 234, 0.15); padding: 0.5rem 0.8rem; 
                                            border-radius: 8px; border-left: 3px solid #667eea;'>
                                    <span style='font-weight: 600; color: rgba(255,255,255,0.95); font-size: 0.8rem;'>
                                        {html_module.escape(display_name)}
                                    </span>
                                    <span style='float: right; color: #4CAF50; font-weight: 700; font-size: 0.8rem;'>
                                        {det_time:.3f}s
                                    </span>
                                </div>
                            """
                        detector_timings_html += """
                            </div>
                        </div>
                        """
                    
                    if is_complete:
                        block_html = f"""
                        <div class='nbdscanner-analysis-block'>
                            <div class='nbdscanner-title'>🧬 NBDScanner Analysis</div>
                            <div class='progress-container'>
                                <div style='display: flex; justify-content: space-between; align-items: center;'>
                                    <span style='color: #ffffff; font-weight: 700; font-size: 0.9rem;'>Progress:</span>
                                    <span style='font-family: monospace; color: #4CAF50; font-size: 1.2rem; letter-spacing: 2px;'>{progress_visual}</span>
                                    <span style='color: #4CAF50; font-weight: 800;'>100%</span>
                                </div>
                                <div class='progress-bar-outer'>
                                    <div class='progress-bar-inner' style='width: 100%;'></div>
                                </div>
                            </div>
                            <div class='completed-badge'>
                                <span class='completed-text'>✓ Analysis Completed</span>
                            </div>
                            <div class='stats-grid'>
                                <div class='stat-item'>
                                    <p class='stat-value'>{elapsed_sec:.1f}s</p>
                                    <p class='stat-label'>Total Time</p>
                                </div>
                                <div class='stat-item'>
                                    <p class='stat-value'>{motifs_found:,}</p>
                                    <p class='stat-label'>Motifs Found</p>
                                </div>
                                <div class='stat-item'>
                                    <p class='stat-value'>{speed_bps:,.0f}</p>
                                    <p class='stat-label'>bp/sec</p>
                                </div>
                                <div class='stat-item'>
                                    <p class='stat-value'>9</p>
                                    <p class='stat-label'>Detectors</p>
                                </div>
                            </div>
                            {detector_timings_html}
                        </div>
                        """
                    else:
                        block_html = f"""
                        <div class='nbdscanner-analysis-block'>
                            <div class='nbdscanner-title'>🧬 NBDScanner Analysis</div>
                            <div class='progress-container'>
                                <div style='display: flex; justify-content: space-between; align-items: center;'>
                                    <span style='color: #ffffff; font-weight: 700; font-size: 0.9rem;'>Progress:</span>
                                    <span style='font-family: monospace; color: #FFD700; font-size: 1.2rem; letter-spacing: 2px;'>{progress_visual}</span>
                                    <span style='color: #FFD700; font-weight: 800;'>{progress_pct:.0f}%</span>
                                </div>
                                <div class='progress-bar-outer'>
                                    <div class='progress-bar-inner' style='width: {progress_pct}%;'></div>
                                </div>
                            </div>
                            <div class='step-indicator'>
                                <span class='step-text'>🔄 Step: {current_step}</span>
                            </div>
                            <div class='stats-grid'>
                                <div class='stat-item'>
                                    <p class='stat-value'>{elapsed_sec:.1f}s</p>
                                    <p class='stat-label'>Elapsed Time</p>
                                </div>
                                <div class='stat-item'>
                                    <p class='stat-value'>{motifs_found:,}</p>
                                    <p class='stat-label'>Motifs Found</p>
                                </div>
                                <div class='stat-item'>
                                    <p class='stat-value'>{speed_bps:,.0f}</p>
                                    <p class='stat-label'>bp/sec</p>
                                </div>
                                <div class='stat-item'>
                                    <p class='stat-value'>9</p>
                                    <p class='stat-label'>Detectors</p>
                                </div>
                            </div>
                        </div>
                        """
                    return block_html
                
                total_bp_all_sequences = sum(len(seq) for seq in st.session_state.seqs)
                estimated_total_time = estimate_time(total_bp_all_sequences)
                
                try:
                    # Filter which classes to analyze based on selection
                    analysis_classes = st.session_state.selected_classes if st.session_state.selected_classes else None
                    
                    # Run analysis on each sequence
                    all_results = []
                    all_hotspots = []
                    
                    total_bp_processed = 0
                    
                    # Track total motifs found so far
                    total_motifs_so_far = 0
                    
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
                        
                        # Calculate current speed
                        current_speed = total_bp_processed / elapsed if elapsed > 0 else 0
                        
                        # Determine current step name (escape for HTML safety)
                        safe_name = html_module.escape(name)
                        current_step_name = f"Analyzing {safe_name} ({i+1}/{len(st.session_state.seqs)})"
                        
                        # Render the combined progress block
                        progress_html = render_progress_block(
                            progress_pct=overall_percentage,
                            current_step=current_step_name,
                            elapsed_sec=elapsed,
                            motifs_found=total_motifs_so_far,
                            speed_bps=current_speed,
                            is_complete=False
                        )
                        combined_progress_placeholder.markdown(progress_html, unsafe_allow_html=True)
                        
                        # Show detailed detector status panel
                        detailed_progress_html = f"""
                        <div style='background: linear-gradient(135deg, #f5f5f5 0%, #e0e0e0 100%); 
                                    border-radius: 12px; padding: 1rem; margin: 0.8rem 0;
                                    border: 1px solid #bdbdbd;'>
                            <h4 style='margin: 0 0 0.8rem 0; color: #424242;'>📋 Detection Pipeline</h4>
                            <p style='margin: 0; font-size: 0.85rem; color: #616161;'>
                                <strong>Sequence:</strong> {safe_name} ({len(seq):,} bp)<br/>
                                <strong>Processed:</strong> {total_bp_processed:,} / {total_bp_all_sequences:,} bp
                            </p>
                            <div style='display: grid; grid-template-columns: repeat(3, 1fr); gap: 0.4rem; margin-top: 0.8rem; font-size: 0.8rem;'>
                        """
                        
                        for j, (detector_name, detector_desc) in enumerate(DETECTOR_PROCESSES):
                            # All detectors run in parallel during analysis
                            status_icon = "🔄"
                            bg_color = "rgba(25, 118, 210, 0.1)"
                            detailed_progress_html += f"""
                                <div style='background: {bg_color}; padding: 0.4rem; border-radius: 6px; 
                                            border-left: 3px solid #1976d2;'>
                                    <span style='font-weight: 600;'>{status_icon} {detector_name}</span>
                                </div>
                            """
                        
                        detailed_progress_html += """
                            </div>
                        </div>
                        """
                        detailed_progress_placeholder.markdown(detailed_progress_html, unsafe_allow_html=True)
                        
                        # Run the analysis - use chunking for large sequences (>10,000 bp)
                        seq_start = time.time()
                        
                        # Create a chunk progress placeholder for large sequences
                        chunk_progress_placeholder = st.empty()
                        
                        # Define progress callback for chunked processing
                        def chunking_progress_callback(current_chunk, total_chunks, bp_done):
                            """Callback to update chunk progress for large sequences"""
                            if show_chunk_progress and total_chunks > 1:
                                chunk_percent = (current_chunk / total_chunks) * 100
                                chunk_progress_placeholder.info(
                                    f"🔄 Processing chunks: {current_chunk}/{total_chunks} ({chunk_percent:.1f}%) - "
                                    f"{bp_done:,} bp processed"
                                )
                        
                        if len(seq) > CHUNK_THRESHOLD:
                            # Large sequence: use chunking with progress tracking
                            if show_chunk_progress:
                                num_chunks = (len(seq) // DEFAULT_CHUNK_SIZE) + 1
                                chunk_progress_placeholder.info(
                                    f"📦 Large sequence detected ({len(seq):,} bp). "
                                    f"Splitting into ~{num_chunks} chunks of {DEFAULT_CHUNK_SIZE:,} bp..."
                                )
                            
                            try:
                                results = analyze_sequence(
                                    seq, 
                                    name, 
                                    use_chunking=True,
                                    chunk_size=DEFAULT_CHUNK_SIZE,
                                    chunk_overlap=DEFAULT_CHUNK_OVERLAP,
                                    progress_callback=chunking_progress_callback if show_chunk_progress else None
                                )
                                
                                # Clear chunk progress
                                if show_chunk_progress:
                                    chunk_progress_placeholder.success(
                                        f"✅ Chunk processing complete: {len(results)} motifs found in {len(seq):,} bp"
                                    )
                            except Exception as e:
                                st.warning(f"Chunked analysis failed, trying standard mode: {e}")
                                results = analyze_sequence(seq, name, use_chunking=False)
                        elif use_parallel_scanner and len(seq) > 100000:
                            # Use experimental parallel scanner for large sequences (legacy path)
                            try:
                                from scanner_agent import ParallelScanner
                                
                                def chunk_progress_callback(current, total):
                                    """Callback to update chunk progress"""
                                    if show_chunk_progress:
                                        chunk_percent = (current / total) * 100
                                        chunk_progress_placeholder.info(f"🔄 Processing chunks: {current}/{total} ({chunk_percent:.1f}%)")
                                
                                # Run parallel scanner
                                scanner = ParallelScanner(seq, hs_db=None)
                                raw_motifs = scanner.run_scan(progress_callback=chunk_progress_callback)
                                
                                # Convert raw motifs to full motif format by running through scoring
                                results = analyze_sequence(seq, name, use_chunking=False)
                                
                                # Clear chunk progress
                                if show_chunk_progress:
                                    chunk_progress_placeholder.success(f"✅ Chunk processing complete: {len(raw_motifs)} raw motifs found")
                                
                            except Exception as e:
                                st.warning(f"Parallel scanner failed, falling back to standard: {e}")
                                results = analyze_sequence(seq, name, use_chunking=False)
                        else:
                            # Standard consolidated NBDScanner analysis (for smaller sequences)
                            results = analyze_sequence(seq, name, use_chunking=False)
                        
                        seq_time = time.time() - seq_start
                        
                        # Ensure all motifs have required fields
                        results = [ensure_subclass(motif) for motif in results]
                        all_results.append(results)
                        
                        total_bp_processed += len(seq)
                        total_motifs_so_far += len(results)
                        
                        # Calculate processing speed
                        elapsed = time.time() - start_time
                        speed = total_bp_processed / elapsed if elapsed > 0 else 0
                        
                        # Update the combined progress block after sequence completion
                        progress_pct = (total_bp_processed / total_bp_all_sequences * 100) if total_bp_all_sequences > 0 else 100
                        current_step_name = f"Completed {safe_name}"
                        progress_html = render_progress_block(
                            progress_pct=progress_pct,
                            current_step=current_step_name,
                            elapsed_sec=elapsed,
                            motifs_found=total_motifs_so_far,
                            speed_bps=speed,
                            is_complete=False
                        )
                        combined_progress_placeholder.markdown(progress_html, unsafe_allow_html=True)
                        
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
                    
                    # Get detector timings from the last analysis
                    detector_timings = get_last_detector_timings()
                    
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
                        ],
                        'detector_timings': detector_timings  # Individual detector timings
                    }
                    
                    # Clear progress displays
                    progress_placeholder.empty()
                    status_placeholder.empty()
                    detailed_progress_placeholder.empty()
                    timer_placeholder.empty()
                    
                    # Show final completed progress block
                    total_motifs_found = sum(len(r) for r in all_results)
                    completion_html = render_progress_block(
                        progress_pct=100,
                        current_step="All detectors completed",
                        elapsed_sec=total_time,
                        motifs_found=total_motifs_found,
                        speed_bps=overall_speed,
                        is_complete=True,
                        detector_timings=detector_timings
                    )
                    combined_progress_placeholder.markdown(completion_html, unsafe_allow_html=True)
                    
                    st.success("Results are available below and in the 'Analysis Results and Visualization' tab.")
                    st.session_state.analysis_status = "Complete"
                    
                except Exception as e:
                    timer_placeholder.empty()
                    progress_placeholder.empty()
                    status_placeholder.empty()
                    detailed_progress_placeholder.empty()
                    combined_progress_placeholder.empty()
                    st.error(f"❌ Analysis failed: {str(e)}")
                    st.session_state.analysis_status = "Error"

        # Show quick summary table if available
        if st.session_state.get('summary_df') is not None:
            st.markdown("""
            <div style='background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
                        border-radius: 16px; padding: 1.5rem; margin-top: 1.5rem;
                        border: 1px solid rgba(0,0,0,0.08);
                        box-shadow: 0 4px 15px rgba(0,0,0,0.05);'>
                <h4 style='margin: 0 0 1rem 0; color: #1a1a2e; display: flex; align-items: center; gap: 0.5rem;'>
                    <span style='font-size: 1.2rem;'>📊</span>
                    Analysis Summary
                </h4>
            </div>
            """, unsafe_allow_html=True)
            st.dataframe(st.session_state.summary_df, use_container_width=True)
    # End of Upload & Analyze tab
    st.markdown("---")
# ---------- RESULTS ----------
with tab_pages["Results"]:
    st.markdown('<h2>Analysis Results and Visualization</h2>', unsafe_allow_html=True)
    if not st.session_state.results:
        st.info("No analysis results. Please run motif analysis first.")
    else:
        # Performance metrics display if available - modern glassmorphism design
        if st.session_state.get('performance_metrics'):
            metrics = st.session_state.performance_metrics
            st.markdown(f"""
            <div style='background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
                        border-radius: 20px; padding: 2rem; margin-bottom: 2rem;
                        box-shadow: 0 15px 50px rgba(15, 52, 96, 0.35);
                        position: relative; overflow: hidden;'>
                <div style='position: absolute; top: 0; left: 0; right: 0; height: 3px;
                            background: linear-gradient(90deg, #667eea, #764ba2, #f64f59, #c471ed, #12c2e9);
                            background-size: 300% 100%;'></div>
                <h3 style='margin: 0 0 1.5rem 0; text-align: center;
                           background: linear-gradient(135deg, #ffffff 0%, #e0e0e0 100%);
                           -webkit-background-clip: text; -webkit-text-fill-color: transparent;
                           font-size: 1.4rem; font-weight: 800;'>
                    ⚡ NBDScanner Performance Metrics
                </h3>
                <div style='display: grid; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr)); gap: 1rem;'>
                    <div style='text-align: center; background: rgba(255,255,255,0.05); 
                                padding: 1.2rem; border-radius: 14px; border: 1px solid rgba(255,255,255,0.08);'>
                        <p style='margin: 0; font-size: 1.8rem; font-weight: 800;
                                  background: linear-gradient(135deg, #12c2e9 0%, #c471ed 100%);
                                  -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                            {metrics['total_time']:.2f}s</p>
                        <p style='margin: 0.3rem 0 0 0; color: rgba(255,255,255,0.7); font-size: 0.75rem;
                                  text-transform: uppercase; letter-spacing: 0.1em;'>Processing Time</p>
                    </div>
                    <div style='text-align: center; background: rgba(255,255,255,0.05); 
                                padding: 1.2rem; border-radius: 14px; border: 1px solid rgba(255,255,255,0.08);'>
                        <p style='margin: 0; font-size: 1.8rem; font-weight: 800;
                                  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                                  -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                            {metrics['total_bp']:,}</p>
                        <p style='margin: 0.3rem 0 0 0; color: rgba(255,255,255,0.7); font-size: 0.75rem;
                                  text-transform: uppercase; letter-spacing: 0.1em;'>Base Pairs</p>
                    </div>
                    <div style='text-align: center; background: rgba(255,255,255,0.05); 
                                padding: 1.2rem; border-radius: 14px; border: 1px solid rgba(255,255,255,0.08);'>
                        <p style='margin: 0; font-size: 1.8rem; font-weight: 800;
                                  background: linear-gradient(135deg, #f64f59 0%, #c471ed 100%);
                                  -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                            {metrics['speed']:,.0f}</p>
                        <p style='margin: 0.3rem 0 0 0; color: rgba(255,255,255,0.7); font-size: 0.75rem;
                                  text-transform: uppercase; letter-spacing: 0.1em;'>bp/second</p>
                    </div>
                    <div style='text-align: center; background: rgba(255,255,255,0.05); 
                                padding: 1.2rem; border-radius: 14px; border: 1px solid rgba(255,255,255,0.08);'>
                        <p style='margin: 0; font-size: 1.8rem; font-weight: 800;
                                  background: linear-gradient(135deg, #00b09b 0%, #96c93d 100%);
                                  -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                            {metrics.get('detector_count', 9)}</p>
                        <p style='margin: 0.3rem 0 0 0; color: rgba(255,255,255,0.7); font-size: 0.75rem;
                                  text-transform: uppercase; letter-spacing: 0.1em;'>Detectors</p>
                    </div>
                    <div style='text-align: center; background: rgba(255,255,255,0.05); 
                                padding: 1.2rem; border-radius: 14px; border: 1px solid rgba(255,255,255,0.08);'>
                        <p style='margin: 0; font-size: 1.8rem; font-weight: 800;
                                  background: linear-gradient(135deg, #ff6b6b 0%, #ffa502 100%);
                                  -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                            {metrics['sequences']}</p>
                        <p style='margin: 0.3rem 0 0 0; color: rgba(255,255,255,0.7); font-size: 0.75rem;
                                  text-transform: uppercase; letter-spacing: 0.1em;'>Sequences</p>
                    </div>
                    <div style='text-align: center; background: rgba(255,255,255,0.05); 
                                padding: 1.2rem; border-radius: 14px; border: 1px solid rgba(255,255,255,0.08);'>
                        <p style='margin: 0; font-size: 1.8rem; font-weight: 800;
                                  background: linear-gradient(135deg, #a8e6cf 0%, #dcedc1 100%);
                                  -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                            {metrics['total_motifs']}</p>
                        <p style='margin: 0.3rem 0 0 0; color: rgba(255,255,255,0.7); font-size: 0.75rem;
                                  text-transform: uppercase; letter-spacing: 0.1em;'>Total Motifs</p>
                    </div>
                </div>
            </div>
            """, unsafe_allow_html=True)
            
            # Display detector timings if available
            detector_timings = metrics.get('detector_timings', {})
            if detector_timings:
                detector_display_names = get_detector_display_names()
                sorted_detectors = sorted(detector_timings.items(), key=lambda x: x[0])
                
                detector_timings_html = """
                <div style='background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
                            border-radius: 16px; padding: 1.5rem; margin-bottom: 2rem;
                            box-shadow: 0 10px 40px rgba(15, 52, 96, 0.25);'>
                    <h4 style='margin: 0 0 1rem 0; text-align: center;
                               color: rgba(255,255,255,0.95); font-size: 1.1rem; font-weight: 700;'>
                        ⏱️ Individual Detector Timings
                    </h4>
                    <div style='display: grid; grid-template-columns: repeat(3, 1fr); gap: 0.6rem;'>
                """
                for det_name, det_time in sorted_detectors:
                    display_name = detector_display_names.get(det_name, det_name)
                    detector_timings_html += f"""
                        <div style='background: rgba(102, 126, 234, 0.15); padding: 0.7rem 1rem; 
                                    border-radius: 10px; border-left: 3px solid #667eea;
                                    display: flex; justify-content: space-between; align-items: center;'>
                            <span style='font-weight: 600; color: rgba(255,255,255,0.95); font-size: 0.85rem;'>
                                {html_module.escape(display_name)}
                            </span>
                            <span style='color: #4CAF50; font-weight: 700; font-size: 0.85rem;'>
                                {det_time:.3f}s
                            </span>
                        </div>
                    """
                detector_timings_html += """
                    </div>
                </div>
                """
                st.markdown(detector_timings_html, unsafe_allow_html=True)
        
        # Enhanced summary display with modern styling
        st.markdown("""
        <div style='background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
                    border-radius: 16px; padding: 1.5rem; margin-bottom: 1.5rem;
                    border: 1px solid rgba(0,0,0,0.08);
                    box-shadow: 0 4px 15px rgba(0,0,0,0.05);'>
            <h3 style='margin: 0; color: #1a1a2e; display: flex; align-items: center; gap: 0.5rem;'>
                <span style='font-size: 1.3rem;'>📊</span>
                Analysis Summary
            </h3>
        </div>
        """, unsafe_allow_html=True)
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
            
            # Enhanced summary card with modern glassmorphism styling
            st.markdown(f"""
            <div style='background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
                        border-radius: 24px; padding: 2.5rem; margin: 1.5rem 0;
                        box-shadow: 0 20px 60px rgba(15, 52, 96, 0.4);
                        position: relative; overflow: hidden;'>
                <div style='position: absolute; top: 0; left: 0; right: 0; height: 4px;
                            background: linear-gradient(90deg, #667eea, #764ba2, #f64f59, #c471ed, #12c2e9);
                            background-size: 300% 100%;'></div>
                <div style='position: absolute; top: -100px; right: -100px; width: 300px; height: 300px;
                            background: radial-gradient(circle, rgba(102,126,234,0.15) 0%, transparent 70%);
                            border-radius: 50%;'></div>
                <h3 style='margin: 0 0 2rem 0; text-align: center; font-size: 1.6rem; font-weight: 800;
                           background: linear-gradient(135deg, #ffffff 0%, #e0e0e0 100%);
                           -webkit-background-clip: text; -webkit-text-fill-color: transparent;
                           position: relative; z-index: 1;'>
                    🧬 NBDScanner Analysis Results
                </h3>
                <div style='display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); 
                            gap: 1.5rem; position: relative; z-index: 1;'>
                    <div style='text-align: center; background: rgba(255,255,255,0.05); 
                                padding: 1.5rem; border-radius: 16px; border: 1px solid rgba(255,255,255,0.08);
                                transition: all 0.3s ease;'>
                        <p style='margin: 0 0 0.5rem 0; font-size: 2.4rem; font-weight: 800;
                                  background: linear-gradient(135deg, #12c2e9 0%, #c471ed 100%);
                                  -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                            {stats.get("Coverage%", 0):.2f}%
                        </p>
                        <p style='margin: 0; font-size: 0.85rem; color: rgba(255,255,255,0.7);
                                  text-transform: uppercase; letter-spacing: 0.1em; font-weight: 600;'>
                            Sequence Coverage
                        </p>
                    </div>
                    <div style='text-align: center; background: rgba(255,255,255,0.05); 
                                padding: 1.5rem; border-radius: 16px; border: 1px solid rgba(255,255,255,0.08);
                                transition: all 0.3s ease;'>
                        <p style='margin: 0 0 0.5rem 0; font-size: 2.4rem; font-weight: 800;
                                  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                                  -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                            {stats.get("Density", 0):.2f}
                        </p>
                        <p style='margin: 0; font-size: 0.85rem; color: rgba(255,255,255,0.7);
                                  text-transform: uppercase; letter-spacing: 0.1em; font-weight: 600;'>
                            Motif Density<br><span style='font-size: 0.7rem;'>(motifs/kb)</span>
                        </p>
                    </div>
                    <div style='text-align: center; background: rgba(255,255,255,0.05); 
                                padding: 1.5rem; border-radius: 16px; border: 1px solid rgba(255,255,255,0.08);
                                transition: all 0.3s ease;'>
                        <p style='margin: 0 0 0.5rem 0; font-size: 2.4rem; font-weight: 800;
                                  background: linear-gradient(135deg, #f64f59 0%, #c471ed 100%);
                                  -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                            {motif_count}
                        </p>
                        <p style='margin: 0; font-size: 0.85rem; color: rgba(255,255,255,0.7);
                                  text-transform: uppercase; letter-spacing: 0.1em; font-weight: 600;'>
                            Total Motifs
                        </p>
                    </div>
                    <div style='text-align: center; background: rgba(255,255,255,0.05); 
                                padding: 1.5rem; border-radius: 16px; border: 1px solid rgba(255,255,255,0.08);
                                transition: all 0.3s ease;'>
                        <p style='margin: 0 0 0.5rem 0; font-size: 2.4rem; font-weight: 800;
                                  background: linear-gradient(135deg, #00b09b 0%, #96c93d 100%);
                                  -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                            {sequence_length:,}
                        </p>
                        <p style='margin: 0; font-size: 0.85rem; color: rgba(255,255,255,0.7);
                                  text-transform: uppercase; letter-spacing: 0.1em; font-weight: 600;'>
                            Sequence Length<br><span style='font-size: 0.7rem;'>(bp)</span>
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
            st.markdown('<h3>📊 NBDScanner Visualizations</h3>', unsafe_allow_html=True)
            
            # Create tabs for different visualization categories including Cluster/Hybrid tab
            viz_tabs = st.tabs(["📈 Distribution", "🗺️ Coverage Map", "📊 Statistics", "🔗 Cluster/Hybrid"])
            
            with viz_tabs[0]:  # Distribution
                st.subheader("Motif Distribution Analysis")
                try:
                    fig1 = plot_motif_distribution(filtered_motifs, by='Class', title=f"Motif Classes - {sequence_name}")
                    st.pyplot(fig1)
                    plt.close(fig1)
                    
                    fig2 = plot_motif_distribution(filtered_motifs, by='Subclass', title=f"Motif Subclasses - {sequence_name}")
                    st.pyplot(fig2) 
                    plt.close(fig2)
                    
                    # Pie chart
                    fig3 = plot_nested_pie_chart(filtered_motifs, title=f"Class-Subclass Distribution - {sequence_name}")
                    st.pyplot(fig3)
                    plt.close(fig3)
                except Exception as e:
                    st.error(f"Error generating distribution plots: {e}")
            
            with viz_tabs[1]:  # Coverage Map
                st.subheader("Sequence Coverage Analysis")
                try:
                    # Coverage map showing motif positions
                    st.markdown("**Motif Position Map**")
                    fig4 = plot_coverage_map(filtered_motifs, sequence_length, title=f"Motif Coverage - {sequence_name}")
                    st.pyplot(fig4)
                    plt.close(fig4)
                    
                    # Add density heatmap for comprehensive coverage analysis
                    st.markdown("**Motif Density Heatmap**")
                    fig5 = plot_density_heatmap(filtered_motifs, sequence_length, 
                                               window_size=max(100, sequence_length // 20),
                                               title=f"Motif Density - {sequence_name}")
                    st.pyplot(fig5)
                    plt.close(fig5)
                except Exception as e:
                    st.error(f"Error generating coverage map: {e}")
            
            with viz_tabs[2]:  # Statistics  
                st.subheader("📊 Statistical Analysis")
                
                # Tabs for different statistical analyses
                stat_tabs = st.tabs(["📏 Density Metrics", "✨ Enrichment Analysis", "📈 Distributions"])
                
                with stat_tabs[0]:  # Density Metrics
                    st.markdown("### Density Calculations")
                    st.markdown("""
                    <div style='background: #e3f2fd; border-left: 4px solid #2196F3; padding: 15px; border-radius: 5px; margin: 10px 0;'>
                    <p><b>Genomic Density (σ_G):</b> Percentage of sequence covered by motifs</p>
                    <p><b>Positional Density (λ):</b> Number of motifs per unit length (kbp or Mbp)</p>
                    </div>
                    """, unsafe_allow_html=True)
                    
                    try:
                        # Calculate density metrics
                        genomic_density = calculate_genomic_density(filtered_motifs, sequence_length, by_class=True)
                        positional_density_kbp = calculate_positional_density(filtered_motifs, sequence_length, unit='kbp', by_class=True)
                        positional_density_mbp = calculate_positional_density(filtered_motifs, sequence_length, unit='Mbp', by_class=True)
                        
                        # Display overall metrics
                        st.markdown("#### Overall Metrics")
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric("Genomic Density", f"{genomic_density.get('Overall', 0):.4f}%")
                        with col2:
                            st.metric("Density (motifs/kbp)", f"{positional_density_kbp.get('Overall', 0):.2f}")
                        with col3:
                            st.metric("Density (motifs/Mbp)", f"{positional_density_mbp.get('Overall', 0):.2f}")
                        
                        # Display per-class density table
                        st.markdown("#### Density by Motif Class")
                        density_data = []
                        for class_name in sorted([k for k in genomic_density.keys() if k != 'Overall']):
                            density_data.append({
                                'Motif Class': class_name,
                                'Genomic Density (%)': f"{genomic_density.get(class_name, 0):.4f}",
                                'Motifs per kbp': f"{positional_density_kbp.get(class_name, 0):.2f}",
                                'Motifs per Mbp': f"{positional_density_mbp.get(class_name, 0):.2f}"
                            })
                        
                        if density_data:
                            density_df = pd.DataFrame(density_data)
                            st.dataframe(density_df, use_container_width=True, height=300)
                        
                        # Visualize density comparison
                        st.markdown("#### Density Visualization")
                        fig_density = plot_density_comparison(genomic_density, positional_density_kbp,
                                                              title="Motif Density Analysis")
                        st.pyplot(fig_density)
                        plt.close(fig_density)
                        
                    except Exception as e:
                        st.error(f"Error calculating density metrics: {e}")
                        import traceback
                        st.error(traceback.format_exc())
                
                with stat_tabs[1]:  # Enrichment Analysis
                    st.markdown("### Enrichment Analysis")
                    st.markdown("""
                    <div style='background: #fff3e0; border-left: 4px solid #ff9800; padding: 15px; border-radius: 5px; margin: 10px 0;'>
                    <p><b>Fold Enrichment:</b> Ratio of observed motif density to background (shuffled sequences)</p>
                    <p><b>P-value:</b> Statistical significance (proportion of shuffled sequences with equal or higher density)</p>
                    <p><b>Method:</b> 100 iterations of sequence shuffling to generate background distribution</p>
                    </div>
                    """, unsafe_allow_html=True)
                    
                    # Option to run enrichment analysis
                    if st.button("🔬 Run Enrichment Analysis (100 shuffles)", key="run_enrichment"):
                        with st.spinner("Performing enrichment analysis... This may take a few minutes."):
                            progress_bar = st.progress(0)
                            status_text = st.empty()
                            
                            def progress_callback(current, total):
                                progress = current / total
                                progress_bar.progress(progress)
                                status_text.text(f"Processing shuffle {current}/{total}...")
                            
                            try:
                                enrichment_results = calculate_enrichment_with_shuffling(
                                    filtered_motifs, 
                                    st.session_state.seqs[seq_idx],
                                    n_shuffles=100,
                                    by_class=True,
                                    progress_callback=progress_callback
                                )
                                
                                # Store results in session state
                                st.session_state.enrichment_results = enrichment_results
                                progress_bar.empty()
                                status_text.empty()
                                st.success("✅ Enrichment analysis completed!")
                                
                            except Exception as e:
                                st.error(f"Error during enrichment analysis: {e}")
                                import traceback
                                st.error(traceback.format_exc())
                    
                    # Display enrichment results if available
                    if hasattr(st.session_state, 'enrichment_results') and st.session_state.enrichment_results:
                        enrichment_results = st.session_state.enrichment_results
                        
                        st.markdown("#### Enrichment Results")
                        
                        # Create summary table
                        enrich_data = []
                        for class_name in sorted([k for k in enrichment_results.keys() if k != 'Overall']):
                            result = enrichment_results[class_name]
                            fe = result.get('fold_enrichment', 0)
                            fe_str = f"{fe:.2f}" if fe != 'Inf' else 'Inf'
                            p_val = result.get('p_value', 1.0)
                            sig = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
                            
                            enrich_data.append({
                                'Motif Class': class_name,
                                'Count': result.get('observed_count', 0),
                                'Observed Density (%)': f"{result.get('observed_density', 0):.4f}",
                                'Background Mean (%)': f"{result.get('background_mean', 0):.4f}",
                                'Fold Enrichment': fe_str,
                                'P-value': f"{p_val:.4f}",
                                'Significance': sig
                            })
                        
                        if enrich_data:
                            enrich_df = pd.DataFrame(enrich_data)
                            st.dataframe(enrich_df, use_container_width=True, height=300)
                            
                            st.markdown("""
                            <div style='background: #f5f5f5; padding: 10px; border-radius: 5px; margin: 10px 0; font-size: 0.9em;'>
                            <b>Significance levels:</b> *** p < 0.001, ** p < 0.01, * p < 0.05, ns = not significant
                            </div>
                            """, unsafe_allow_html=True)
                        
                        # Visualizations
                        st.markdown("#### Enrichment Visualizations")
                        
                        try:
                            fig_enrich = plot_enrichment_analysis(enrichment_results,
                                                                 title="Motif Enrichment Analysis")
                            st.pyplot(fig_enrich)
                            plt.close(fig_enrich)
                            
                            # Summary table visualization
                            fig_table = plot_enrichment_summary_table(enrichment_results,
                                                                     title="Enrichment Summary Statistics")
                            if fig_table:
                                st.pyplot(fig_table)
                                plt.close(fig_table)
                        except Exception as e:
                            st.error(f"Error generating enrichment plots: {e}")
                    else:
                        st.info("👆 Click the button above to run enrichment analysis")
                
                with stat_tabs[2]:  # Distributions
                    st.markdown("### Distribution Analysis")
                    try:
                        # Length distribution by class - violin plot for better distribution visualization
                        st.markdown("**Motif Length Distribution by Class**")
                        fig6 = plot_length_distribution(filtered_motifs, by_class=True, 
                                                       title="Length Distribution by Motif Class") 
                        st.pyplot(fig6)
                        plt.close(fig6)
                        
                        # Score distribution by class - box plot for score comparison
                        st.markdown("**Motif Score Distribution by Class**")
                        fig7 = plot_score_distribution(filtered_motifs, by_class=True,
                                                      title="Score Distribution by Motif Class")
                        st.pyplot(fig7)
                        plt.close(fig7)
                    except Exception as e:
                        st.error(f"Error generating distribution plots: {e}")
            
            with viz_tabs[3]:  # Cluster/Hybrid
                st.subheader("🔗 Hybrid and Cluster Motifs")
                
                if not hybrid_cluster_motifs:
                    st.info("No Hybrid or Cluster motifs detected for this sequence.")
                else:
                    st.markdown(f"""
                    <div style='background: #f0f8ff; border-left: 4px solid #1e88e5; padding: 15px; border-radius: 5px; margin: 10px 0;'>
                    <h4 style='margin-top: 0; color: #1565c0;'>About Hybrid and Cluster Motifs</h4>
                    <p><b>Hybrid Motifs:</b> Regions where different non-B DNA classes overlap (30-70% overlap between classes).</p>
                    <p><b>Cluster Motifs:</b> High-density regions containing multiple non-B DNA motifs from different classes.</p>
                    <p>These motifs represent complex genomic regions with potential biological significance.</p>
                    </div>
                    """, unsafe_allow_html=True)
                    
                    # Create DataFrame for hybrid/cluster motifs
                    hc_df = pd.DataFrame(hybrid_cluster_motifs)
                    
                    # Summary statistics
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        hybrid_count = len([m for m in hybrid_cluster_motifs if m.get('Class') == 'Hybrid'])
                        st.metric("Hybrid Motifs", hybrid_count)
                    with col2:
                        cluster_count = len([m for m in hybrid_cluster_motifs if m.get('Class') == 'Non-B_DNA_Clusters'])
                        st.metric("Cluster Motifs", cluster_count)
                    with col3:
                        avg_length = int(hc_df['Length'].mean()) if 'Length' in hc_df.columns else 0
                        st.metric("Avg Length (bp)", avg_length)
                    
                    # Detailed table
                    st.markdown("### 📋 Detailed Cluster/Hybrid Table")
                    display_cols = ['Class', 'Subclass', 'Start', 'End', 'Length', 'Score']
                    if 'Component_Classes' in hc_df.columns:
                        display_cols.append('Component_Classes')
                    if 'Motif_Count' in hc_df.columns:
                        display_cols.append('Motif_Count')
                    if 'Class_Diversity' in hc_df.columns:
                        display_cols.append('Class_Diversity')
                    
                    # Filter to only show available columns
                    available_display_cols = [col for col in display_cols if col in hc_df.columns]
                    display_hc_df = hc_df[available_display_cols].copy()
                    display_hc_df.columns = [col.replace('_', ' ') for col in display_hc_df.columns]
                    st.dataframe(display_hc_df, use_container_width=True, height=300)
                    
                    # Visualizations for hybrid/cluster
                    st.markdown("### 📊 Visualizations")
                    
                    viz_col1, viz_col2 = st.columns(2)
                    
                    with viz_col1:
                        # Position map
                        st.markdown("**Position Map**")
                        fig, ax = plt.subplots(figsize=(10, 4))
                        colors_map = {'Hybrid': '#ff6b6b', 'Non-B_DNA_Clusters': '#4ecdc4'}
                        for i, motif in enumerate(hybrid_cluster_motifs):
                            color = colors_map.get(motif.get('Class'), '#95a5a6')
                            ax.barh(i, motif['End'] - motif['Start'], left=motif['Start'], 
                                   height=0.8, color=color, alpha=0.7, label=motif.get('Class'))
                        ax.set_xlabel("Sequence Position (bp)")
                        ax.set_ylabel("Motif Index")
                        ax.set_title("Hybrid/Cluster Position Map")
                        # Remove duplicate labels
                        handles, labels = ax.get_legend_handles_labels()
                        by_label = dict(zip(labels, handles))
                        ax.legend(by_label.values(), by_label.keys())
                        plt.tight_layout()
                        st.pyplot(fig)
                        plt.close(fig)
                    
                    with viz_col2:
                        # Raw score distribution (no normalized scores as per requirements)
                        st.markdown("**Raw Score Distribution**")
                        fig, ax = plt.subplots(figsize=(10, 4))
                        hybrid_scores = [m.get('Score', 0) for m in hybrid_cluster_motifs if m.get('Class') == 'Hybrid']
                        cluster_scores = [m.get('Score', 0) for m in hybrid_cluster_motifs if m.get('Class') == 'Non-B_DNA_Clusters']
                        
                        if hybrid_scores:
                            ax.hist(hybrid_scores, bins=10, alpha=0.6, color='#ff6b6b', label='Hybrid')
                        if cluster_scores:
                            ax.hist(cluster_scores, bins=10, alpha=0.6, color='#4ecdc4', label='Cluster')
                        
                        ax.set_xlabel("Raw Score")
                        ax.set_ylabel("Frequency")
                        ax.set_title("Raw Score Distribution")
                        ax.legend()
                        plt.tight_layout()
                        st.pyplot(fig)
                        plt.close(fig)

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


# ---------- DISEASE ANALYSIS ----------
with tab_pages["Disease Analysis"]:
    st.header("🧬 Sequence Analysis - Disease Risk Assessment")
    
    # Show input data analysis if available
    if st.session_state.results:
        # Collect all motifs from input data
        input_motifs = []
        for i, motifs in enumerate(st.session_state.results):
            for m in motifs:
                m_copy = m.copy()
                m_copy['Sequence_Name'] = st.session_state.names[i]
                input_motifs.append(m_copy)
        
        if input_motifs:
            # Count disease-relevant motifs
            disease_classes = ['G-Quadruplex', 'Slipped_DNA', 'Cruciform', 'Triplex', 'Z-DNA', 'R-Loop']
            disease_counts = {cls: 0 for cls in disease_classes}
            for m in input_motifs:
                cls = m.get('Class', 'Unknown')
                if cls in disease_counts:
                    disease_counts[cls] += 1
            
            total_disease_relevant = sum(disease_counts.values())
            total_motifs = len(input_motifs)
            
            # Modern header card with disease summary
            st.markdown(f"""
            <div style='background: linear-gradient(135deg, #1a237e 0%, #283593 50%, #3949ab 100%); 
                        border-radius: 20px; padding: 2.5rem; margin-bottom: 2rem;
                        box-shadow: 0 12px 40px rgba(26, 35, 126, 0.35);
                        position: relative; overflow: hidden;'>
                <div style='position: absolute; top: -50px; right: -50px; width: 200px; height: 200px; 
                            background: radial-gradient(circle, rgba(255,255,255,0.1) 0%, transparent 70%);
                            border-radius: 50%;'></div>
                <div style='position: relative; z-index: 1;'>
                    <h2 style='margin: 0 0 0.5rem 0; color: #ffffff; font-size: 1.8rem; font-weight: 800;
                               text-shadow: 0 2px 4px rgba(0,0,0,0.2);'>
                        🧬 Disease Risk Assessment Report
                    </h2>
                    <p style='color: rgba(255,255,255,0.85); font-size: 1.1rem; margin-bottom: 1.5rem;'>
                        Comprehensive analysis of disease-associated Non-B DNA structures in your sequences
                    </p>
                    <div style='display: grid; grid-template-columns: repeat(3, 1fr); gap: 1rem;'>
                        <div style='background: rgba(255,255,255,0.15); backdrop-filter: blur(10px);
                                    border-radius: 12px; padding: 1.2rem; text-align: center;
                                    border: 1px solid rgba(255,255,255,0.2);'>
                            <p style='color: #ffd54f; font-size: 2.2rem; font-weight: 800; margin: 0;
                                      text-shadow: 0 0 20px rgba(255,213,79,0.5);'>{total_disease_relevant}</p>
                            <p style='color: rgba(255,255,255,0.9); font-size: 0.85rem; margin: 0.3rem 0 0 0;
                                      font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em;'>
                                Disease-Relevant Motifs</p>
                        </div>
                        <div style='background: rgba(255,255,255,0.15); backdrop-filter: blur(10px);
                                    border-radius: 12px; padding: 1.2rem; text-align: center;
                                    border: 1px solid rgba(255,255,255,0.2);'>
                            <p style='color: #81d4fa; font-size: 2.2rem; font-weight: 800; margin: 0;'>{total_motifs}</p>
                            <p style='color: rgba(255,255,255,0.9); font-size: 0.85rem; margin: 0.3rem 0 0 0;
                                      font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em;'>
                                Total Motifs Analyzed</p>
                        </div>
                        <div style='background: rgba(255,255,255,0.15); backdrop-filter: blur(10px);
                                    border-radius: 12px; padding: 1.2rem; text-align: center;
                                    border: 1px solid rgba(255,255,255,0.2);'>
                            <p style='color: #a5d6a7; font-size: 2.2rem; font-weight: 800; margin: 0;'>
                                {len(st.session_state.seqs)}</p>
                            <p style='color: rgba(255,255,255,0.9); font-size: 0.85rem; margin: 0.3rem 0 0 0;
                                      font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em;'>
                                Sequences Scanned</p>
                        </div>
                    </div>
                </div>
            </div>
            """, unsafe_allow_html=True)
            
            # Disease-relevant motif breakdown with elegant cards
            st.markdown("""
            <h3 style='color: #1a237e; margin: 2rem 0 1.5rem 0; font-weight: 700;'>
                🔬 Disease-Associated Motif Breakdown
            </h3>
            """, unsafe_allow_html=True)
            
            # Create 3x2 grid for disease motif cards
            disease_info = {
                'G-Quadruplex': {
                    'icon': '🧬', 'color': '#e91e63', 'gradient': 'linear-gradient(135deg, #e91e63 0%, #c2185b 100%)',
                    'diseases': 'Fragile X Syndrome, ALS/FTD, Multiple Cancers',
                    'mechanism': 'Transcription blockage, replication stalling'
                },
                'Slipped_DNA': {
                    'icon': '🔄', 'color': '#9c27b0', 'gradient': 'linear-gradient(135deg, #9c27b0 0%, #7b1fa2 100%)',
                    'diseases': "Huntington's Disease, Myotonic Dystrophy",
                    'mechanism': 'Repeat expansion, trinucleotide instability'
                },
                'Cruciform': {
                    'icon': '✚', 'color': '#673ab7', 'gradient': 'linear-gradient(135deg, #673ab7 0%, #512da8 100%)',
                    'diseases': 'Chromosomal Translocations, Lymphomas',
                    'mechanism': 'Double-strand breaks, aberrant recombination'
                },
                'Triplex': {
                    'icon': '🔺', 'color': '#3f51b5', 'gradient': 'linear-gradient(135deg, #3f51b5 0%, #303f9f 100%)',
                    'diseases': "Friedreich's Ataxia, Werner Syndrome",
                    'mechanism': 'Transcription inhibition, replication blocks'
                },
                'Z-DNA': {
                    'icon': '⚡', 'color': '#2196f3', 'gradient': 'linear-gradient(135deg, #2196f3 0%, #1976d2 100%)',
                    'diseases': 'Autoimmune Disorders, Systemic Lupus',
                    'mechanism': 'Immune activation, gene regulation'
                },
                'R-Loop': {
                    'icon': '🔗', 'color': '#00bcd4', 'gradient': 'linear-gradient(135deg, #00bcd4 0%, #0097a7 100%)',
                    'diseases': 'Neurodegeneration, Genome Instability',
                    'mechanism': 'Replication-transcription conflicts'
                }
            }
            
            col1, col2, col3 = st.columns(3)
            cols = [col1, col2, col3]
            
            # Risk level thresholds for disease association assessment
            DISEASE_RISK_HIGH_THRESHOLD = 5  # More than 5 motifs = High Risk
            DISEASE_TEXT_TRUNCATE_LENGTH = 50  # Max characters for disease text display
            
            for idx, (cls, info) in enumerate(disease_info.items()):
                with cols[idx % 3]:
                    count = disease_counts.get(cls, 0)
                    risk_level = "High" if count > DISEASE_RISK_HIGH_THRESHOLD else "Moderate" if count > 0 else "Low"
                    risk_color = "#f44336" if risk_level == "High" else "#ff9800" if risk_level == "Moderate" else "#4caf50"
                    
                    # Truncate disease text for display
                    disease_text = info["diseases"][:DISEASE_TEXT_TRUNCATE_LENGTH]
                    disease_ellipsis = '...' if len(info["diseases"]) > DISEASE_TEXT_TRUNCATE_LENGTH else ''
                    
                    st.markdown(f"""
                    <div style='background: {info["gradient"]}; border-radius: 16px; padding: 1.5rem;
                                margin-bottom: 1rem; box-shadow: 0 8px 24px rgba(0,0,0,0.15);
                                transition: transform 0.3s ease, box-shadow 0.3s ease;'>
                        <div style='display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.8rem;'>
                            <span style='font-size: 1.8rem;'>{info["icon"]}</span>
                            <span style='background: {risk_color}; color: white; padding: 0.3rem 0.8rem;
                                        border-radius: 20px; font-size: 0.75rem; font-weight: 700;
                                        text-transform: uppercase;'>{risk_level} Risk</span>
                        </div>
                        <h4 style='color: white; margin: 0 0 0.3rem 0; font-size: 1.1rem; font-weight: 700;'>
                            {cls.replace("_", " ")}</h4>
                        <p style='color: rgba(255,255,255,0.95); font-size: 2rem; font-weight: 800; margin: 0.5rem 0;
                                  text-shadow: 0 2px 4px rgba(0,0,0,0.2);'>{count} <span style='font-size: 1rem;
                                  font-weight: 500;'>motifs</span></p>
                        <p style='color: rgba(255,255,255,0.85); font-size: 0.8rem; margin: 0.5rem 0 0 0;
                                  line-height: 1.4;'>
                            <strong>Associated:</strong> {disease_text}{disease_ellipsis}
                        </p>
                    </div>
                    """, unsafe_allow_html=True)
            
            # Detailed motifs table section
            st.markdown("""
            <h3 style='color: #1a237e; margin: 2rem 0 1rem 0; font-weight: 700;'>
                📋 Detailed Disease-Relevant Motifs
            </h3>
            """, unsafe_allow_html=True)
            
            # Filter only disease-relevant motifs
            disease_relevant_motifs = [m for m in input_motifs if m.get('Class') in disease_classes]
            
            if disease_relevant_motifs:
                dr_df = pd.DataFrame(disease_relevant_motifs)
                display_cols = ['Sequence_Name', 'Class', 'Subclass', 'Start', 'End', 'Length', 'Score']
                available_cols = [c for c in display_cols if c in dr_df.columns]
                
                # Add disease association column
                dr_df['Disease_Association'] = dr_df['Class'].map({
                    'G-Quadruplex': 'Fragile X, ALS/FTD, Cancer',
                    'Slipped_DNA': "Huntington's, Myotonic Dystrophy",
                    'Cruciform': 'Chromosomal translocations',
                    'Triplex': "Friedreich's Ataxia",
                    'Z-DNA': 'Autoimmune disorders',
                    'R-Loop': 'Neurodegeneration'
                })
                
                if 'Disease_Association' in dr_df.columns:
                    available_cols.append('Disease_Association')
                
                st.dataframe(dr_df[available_cols], use_container_width=True, height=350)
                st.caption(f"Showing {len(disease_relevant_motifs)} disease-relevant motifs from your analysis")
            else:
                st.info("No disease-relevant motifs found in the analyzed sequences.")
                
        else:
            st.info("Analysis complete but no motifs were detected. Try analyzing a different sequence.")
    else:
        # Show placeholder when no results available
        st.markdown("""
        <div style='background: linear-gradient(135deg, #f5f5f5 0%, #eeeeee 100%); 
                    border-radius: 20px; padding: 3rem; text-align: center;
                    border: 2px dashed #bdbdbd; margin: 2rem 0;'>
            <div style='font-size: 4rem; margin-bottom: 1rem;'>🧬</div>
            <h3 style='color: #424242; margin: 0 0 1rem 0;'>No Analysis Results Available</h3>
            <p style='color: #757575; font-size: 1.1rem; max-width: 500px; margin: 0 auto 1.5rem auto;'>
                Upload and analyze DNA sequences in the <strong>'Upload & Analyze'</strong> tab 
                to see disease-relevant Non-B DNA motifs from your data.
            </p>
            <div style='display: inline-flex; gap: 0.5rem; flex-wrap: wrap; justify-content: center;'>
                <span style='background: #e3f2fd; color: #1565c0; padding: 0.5rem 1rem; 
                            border-radius: 20px; font-size: 0.85rem; font-weight: 600;'>G-Quadruplex</span>
                <span style='background: #fce4ec; color: #c2185b; padding: 0.5rem 1rem; 
                            border-radius: 20px; font-size: 0.85rem; font-weight: 600;'>Cruciform</span>
                <span style='background: #e8f5e9; color: #2e7d32; padding: 0.5rem 1rem; 
                            border-radius: 20px; font-size: 0.85rem; font-weight: 600;'>Triplex</span>
                <span style='background: #fff3e0; color: #ef6c00; padding: 0.5rem 1rem; 
                            border-radius: 20px; font-size: 0.85rem; font-weight: 600;'>Z-DNA</span>
            </div>
        </div>
        """, unsafe_allow_html=True)


# ---------- DOCUMENTATION ----------
with tab_pages["Documentation"]:
    st.header("Scientific Documentation & References")
    
    # Motif classes documentation
    st.markdown("""
    <div style='background:#f4faff; border-radius:12px; padding:18px; font-size:1.08rem; font-family:Montserrat,Arial;'>
    <b>Motif Classes Detected:</b><br><br>
    <ul>
        <li><b>Curved DNA</b>: Identifies phased poly(A) or poly(T) tracts using regex and spacing rules, reflecting intrinsic curvature. Scoring is based on tract length/grouping.</li>
        <li><b>Z-DNA</b>: Detects alternating purine-pyrimidine patterns, GC-rich segments. Uses windowed scoring; regex finds dinucleotide repeats.</li>
        <li><b>eGZ-motif (Extruded-G Z-DNA)</b>: Searches for long (CGG)<sub>n</sub> runs via regex. Scored by repeat count.</li>
        <li><b>Slipped DNA</b>: Recognizes direct/tandem repeats by repeat-unit matching and regex. Scoring by length and unit copies.</li>
        <li><b>R-Loop</b>: Finds G-rich regions for stable RNA-DNA hybrids; RLFS model and regex. Thermodynamic scoring for hybrid stability.</li>
        <li><b>Cruciform</b>: Finds palindromic inverted repeats with spacers, regex and reverse complement. Scoring by arm length and A/T content.</li>
        <li><b>Triplex DNA / Mirror Repeat</b>: Detects purine/pyrimidine mirror repeats/triplex motifs. Regex identifies units; scoring by composition/purity.</li>
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
