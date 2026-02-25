"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Upload & Analyze Page - Sequence Upload and Motif Analysis                   │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘
"""
# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
import streamlit as st
import time
import os
import gc
import hashlib
import numpy as np
import pandas as pd
import logging
import traceback

from Utilities.config.text import UI_TEXT
from Utilities.config.typography import FONT_CONFIG
from Utilities.config.themes import TAB_THEMES
from Utilities.config.analysis import ANALYSIS_CONFIG, GC_BALANCE_MIN, GC_BALANCE_MAX
from Utilities.config.motif_taxonomy import (
    MOTIF_CLASSIFICATION,
    get_all_classes,
    build_motif_selector_data,
    get_enabled_from_selector_data
)
from Utilities.config.colors import CLASS_COLORS, CLASS_COLOR_NAMES
from UI.css import load_css, get_page_colors
from UI.formatters import format_time_scientific, format_time_human, format_time
from UI.headers import render_section_heading
from UI.performance_stats import PerformanceTracker, format_performance_summary, format_motif_distribution
import io as _io

from Utilities.utilities import (
    parse_fasta_chunked,
    get_file_preview,
    get_basic_stats,
    trigger_garbage_collection,
    parse_fasta,
    get_memory_usage_mb,
    calculate_genomic_density,
    calculate_positional_density,
    _NON_IUPAC_RE,
    plot_motif_distribution,
    plot_linear_motif_track,
    plot_linear_subclass_track,
    plot_density_comparison,
    plot_motif_length_kde,
    plot_score_violin,
    plot_nested_pie_chart,
    plot_structural_heatmap,
    plot_motif_network,
    plot_motif_cooccurrence_matrix,
    plot_chromosome_density,
    plot_motif_clustering_distance,
    plot_spacer_loop_variation,
    plot_cluster_size_distribution,
    plot_structural_competition_upset,
)
from Utilities.nonbscanner import analyze_sequence
from Utilities.job_manager import save_job_results, generate_job_id
from Utilities.disk_storage import UniversalSequenceStorage, UniversalResultsStorage
from Utilities.chunk_analyzer import ChunkAnalyzer
from Utilities.detectors_utils import calc_gc_content, _count_bases
from Utilities.multifasta_engine import analyze_sequences_parallel

# Import progress tracking for Streamlit UI (optional, graceful degradation)
try:
    from Utilities.streamlit_progress import analyze_with_streamlit_progress
    STREAMLIT_PROGRESS_AVAILABLE = True
except ImportError:
    STREAMLIT_PROGRESS_AVAILABLE = False
    logger.warning("Streamlit progress tracking not available")

try: from Bio import Entrez, SeqIO; BIO_AVAILABLE = True
except ImportError: BIO_AVAILABLE = False

logger = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════════
# TUNABLE PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════════
CONFIG_AVAILABLE = False
GRID_COLUMNS = 6; GRID_GAP = "0.10rem"; ROW_GAP = "0.10rem"; DOT_SIZE = 5; GLOW_SIZE = 5
# Chunk analysis threshold: always use ChunkAnalyzer for all sequences (0 = no minimum)
CHUNK_ANALYSIS_THRESHOLD_BP = 0  # Always chunk with 50KB/2KB approach
# Parallel processing threshold: multi-FASTA with >= this many sequences will use parallel processing
PARALLEL_PROCESSING_THRESHOLD = 2  # Use parallel processing for 2+ sequences
SUBMOTIF_ABBREVIATIONS = {
    'Global Curvature': 'Global Curv', 'Local Curvature': 'Local Curv',
    'Direct Repeat': 'DR', 'STR': 'STR', 'Cruciform forming IRs': 'Cruciform',
    'R-loop formation sites': 'R-loop', 'Triplex': 'Triplex', 'Sticky DNA': 'Sticky DNA',
    'Telomeric G4': 'Telo G4', 'Stacked G4': 'Stacked G4',
    'Canonical intramolecular G4': 'Intra G4', 'Extended-loop canonical': 'Ext. Loop G4',
    'Higher-order G4 array/G4-wire': 'G4 Array', 'Intramolecular G-triplex': 'G-triplex', 'Two-tetrad weak PQS': 'Weak PQS', 'Bulged G4': 'Bulged G4',
    'Canonical i-motif': 'i-Motif', 'Relaxed i-motif': 'Relaxed iM', 'AC-motif': 'AC Motif',
    'Z-DNA': 'Z-DNA', 'eGZ': 'eGZ', 'A-philic DNA': 'A-DNA',
    'Dynamic overlaps': 'Hybrid', 'Dynamic clusters': 'Cluster',
}
# ═══════════════════════════════════════════════════════════════════════════════


def _fig_to_bytes(fig) -> bytes:
    """Convert a matplotlib Figure to PNG bytes and close it.

    Used during visualization pre-generation to store rendered figures as
    compact bytes in session state so the Results tab can display them
    instantly without re-rendering.
    """
    import matplotlib.pyplot as _plt
    buf = _io.BytesIO()
    fig.savefig(buf, format='png', dpi=80, bbox_inches='tight')
    _plt.close(fig)
    buf.seek(0)
    return buf.getvalue()


def get_abbreviated_label(subclass: str) -> str:
    """
    Get abbreviated label for display in the submotif selector grid.
    
    Returns the abbreviated label for a given subclass name if a mapping exists
    in SUBMOTIF_ABBREVIATIONS, otherwise returns the original subclass name.
    
    Args:
        subclass: Full canonical subclass name from MOTIF_CLASSIFICATION
        
    Returns:
        Abbreviated label for compact display, or original name if not mapped
    """
    return SUBMOTIF_ABBREVIATIONS.get(subclass, subclass)


def ensure_subclass(motif):
    """Guarantee every motif has a string 'Subclass'"""
    if isinstance(motif, dict):
        if 'Subclass' not in motif or motif['Subclass'] is None:
            motif['Subclass'] = motif.get('Subtype', 'Other')
        return motif
    else:
        # Handle non-dict motifs gracefully
        return {'Subclass': 'Other', 'Motif': motif}


def render_sequence_stats_card(idx: int, name: str, length: int, gc_pct: float, at_pct: float, 
                               gradient_colors: str = "135deg, #667eea 0%, #764ba2 100%") -> str:
    """
    Generate styled HTML card for sequence statistics.
    
    Args:
        idx: Sequence index (1-based)
        name: Sequence name
        length: Sequence length in bp
        gc_pct: GC content percentage
        at_pct: AT content percentage
        gradient_colors: CSS gradient colors
        
    Returns:
        HTML string for the stats card
    """
    # Truncate long sequence names
    display_name = name[:60] + ('...' if len(name) > 60 else '')
    # Determine GC balance indicator using configurable thresholds
    gc_balance = 'OK' if GC_BALANCE_MIN <= gc_pct <= GC_BALANCE_MAX else '!'
    
    return f"""
    <div style='background: linear-gradient({gradient_colors}); 
                border-radius: 8px; padding: 16px; margin: 8px 0; color: white;
                box-shadow: 0 2px 6px rgba(102, 126, 234, 0.2);'>
        <div style='font-weight: 600; font-size: 1.05rem; margin-bottom: 8px;'>
            {idx}. {display_name}
        </div>
        <div style='display: grid; grid-template-columns: repeat(4, 1fr); gap: 10px; 
                    font-size: 0.95rem; margin-top: 8px;'>
            <div style='background: rgba(255,255,255,0.15); padding: 8px; border-radius: 4px; text-align: center;'>
                <div style='font-weight: 700; font-size: 1.1rem;'>{length:,}</div>
                <div style='opacity: 0.9; font-size: 0.82rem;'>Base Pairs</div>
            </div>
            <div style='background: rgba(255,255,255,0.15); padding: 8px; border-radius: 4px; text-align: center;'>
                <div style='font-weight: 700; font-size: 1.1rem;'>{gc_pct:.2f}%</div>
                <div style='opacity: 0.9; font-size: 0.82rem;'>GC Content</div>
            </div>
            <div style='background: rgba(255,255,255,0.15); padding: 8px; border-radius: 4px; text-align: center;'>
                <div style='font-weight: 700; font-size: 1.1rem;'>{at_pct:.2f}%</div>
                <div style='opacity: 0.9; font-size: 0.82rem;'>AT Content</div>
            </div>
            <div style='background: rgba(255,255,255,0.15); padding: 8px; border-radius: 4px; text-align: center;'>
                <div style='font-weight: 700; font-size: 1.1rem;'>{gc_balance}</div>
                <div style='opacity: 0.9; font-size: 0.82rem;'>GC Balance</div>
            </div>
        </div>
    </div>
    """


# Helper function to format time in compact format
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


def calculate_gc_percentage(sequences: list) -> float:
    """
    Calculate average GC percentage across multiple sequences.
    
    Delegates to calc_gc_content() from detectors_utils.py to ensure
    consistency with all other GC calculations in the codebase.
    
    Args:
        sequences: List of DNA sequence strings
        
    Returns:
        Average GC percentage (0.0 if no sequences/base pairs)
        
    Note:
        Uses standardized calc_gc_content() method for consistency.
        See: Utilities/detectors_utils.py
    """
    if not sequences:
        return 0.0
    
    # Calculate total base pairs and weighted GC in a single pass
    total_bp = 0
    total_gc_weighted = 0.0
    for seq in sequences:
        seq_len = len(seq)
        if seq_len > 0:
            total_bp += seq_len
            total_gc_weighted += calc_gc_content(seq) * seq_len
    
    if total_bp == 0:
        return 0.0
    
    return total_gc_weighted / total_bp


def clear_analysis_placeholders(progress_placeholder, status_placeholder,
                                detailed_progress_placeholder, timer_placeholder):
    """
    Helper function to clear all analysis UI placeholders during error handling.
    
    Args:
        progress_placeholder: Streamlit placeholder for progress bar
        status_placeholder: Streamlit placeholder for status messages
        detailed_progress_placeholder: Streamlit placeholder for detailed progress info
        timer_placeholder: Streamlit placeholder for timer/elapsed time display
    """
    progress_placeholder.empty()
    status_placeholder.empty()
    detailed_progress_placeholder.empty()
    timer_placeholder.empty()


def show_technical_details():
    """
    Show technical error details in collapsible expander.
    
    Must be called within an exception handler context where traceback.format_exc()
    can capture the active exception information.
    
    Example:
        try:
            risky_operation()
        except Exception:
            show_technical_details()  # Shows traceback of caught exception
    """
    with st.expander("🔧 Technical Details (for debugging)"):
        st.code(traceback.format_exc())


# ═══════════════════════════════════════════════════════════════════════════════
# EXAMPLE DATA
# Load from examples/ folder; fall back to inline strings if files are missing
# ═══════════════════════════════════════════════════════════════════════════════
_EXAMPLES_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")


def _load_example_file(filename: str, fallback: str) -> str:
    """Return file contents from examples/ dir, or fallback string if unavailable."""
    path = os.path.join(_EXAMPLES_DIR, filename)
    try:
        with open(path, "r") as fh:
            return fh.read()
    except OSError:
        logger.warning("Example file not found: %s — using inline fallback", path)
        return fallback


# Example FASTA data
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
ATCGATCGATCGATCGATCGATCGCCCTACCCTACCCTACCCAAACCCTACCCTACCCTACCCAAAAAAAAAAAAA
AAAAAAAAAAAGATCTAGATCTAGATCTAGATCTAGATCTAGATCTGAAAGAAAGAAAGAAAGAAAGAAAGAAA
"""


def render():
    """Render the Upload & Analyze tab content."""
    # Apply Upload tab theme based on configuration
    load_css(TAB_THEMES.get('Upload & Analyze', 'nature_green'))
    
    # ============================================================
    # INITIALIZE CRITICAL SESSION STATE VARIABLES
    # Ensure selected_classes and selected_subclasses exist before any logic checks them
    # Initialize with all classes/subclasses enabled by default (matching checkbox defaults)
    # This prevents race conditions on server deployments (Streamlit Community Cloud)
    # ============================================================
    if 'selected_classes' not in st.session_state:
        # Initialize with all classes from taxonomy
        st.session_state.selected_classes = list(get_all_classes())
    if 'selected_subclasses' not in st.session_state:
        # Initialize with all subclasses from taxonomy
        all_subclasses = []
        for class_id in sorted(MOTIF_CLASSIFICATION.keys()):
            entry = MOTIF_CLASSIFICATION[class_id]
            all_subclasses.extend(entry['subclasses'])
        st.session_state.selected_subclasses = all_subclasses
    
    # Uniform section heading with page-specific color
    render_section_heading("Upload & Analyze Sequences", page="Upload & Analyze")

    # ============================================================
    # TWO-COLUMN LAYOUT: SEQUENCE CONTEXT (LEFT) | DETECTION SCOPE (RIGHT)
    # Clean scientific layout: sequence input left, motif selection right
    # Analysis execution moves to full-width section below
    # ============================================================

    left_col, right_col = st.columns([0.7, 1.8])

    # ============================================================
    # LEFT COLUMN — SEQUENCE CONTEXT
    # Compact input + QC summary for scientific workflow
    # ============================================================
    with left_col:
        st.markdown("""
        <div style='background: linear-gradient(135deg, #d1fae5 0%, #a7f3d0 100%); 
                    padding: 0.5rem 1rem; border-radius: 6px; margin-bottom: 0.75rem;
                    border-left: 3px solid #10b981;'>
            <h4 style='margin: 0; color: #065f46; font-size: 0.95rem; font-weight: 600;'>
                Sequence Context
            </h4>
        </div>
        """, unsafe_allow_html=True)
        
        # ----- Input Method with colorful vertical buttons -----
        # Inject comprehensive CSS for colorful radio buttons
        st.markdown("""
        <style>
        /* ===== SEQUENCE CONTEXT - Colorful Radio Buttons ===== */
        /* Style each radio option label - these are inside a radiogroup */
        label[data-baseweb="radio"] {
            padding: 0.6rem 1rem !important;
            border-radius: 8px !important;
            transition: all 0.2s ease !important;
            margin-bottom: 0.35rem !important;
            width: 100% !important;
        }
        
        /* Color 1: Upload FASTA File - Green */
        label[data-baseweb="radio"]:has(input[value="0"]) {
            background: linear-gradient(135deg, #d1fae5 0%, #a7f3d0 100%) !important;
            border: 2px solid #10b981 !important;
        }
        label[data-baseweb="radio"]:has(input[value="0"]):hover {
            background: linear-gradient(135deg, #a7f3d0 0%, #6ee7b7 100%) !important;
            box-shadow: 0 2px 8px rgba(16, 185, 129, 0.3) !important;
        }
        label[data-baseweb="radio"]:has(input[value="0"]) p {
            color: #065f46 !important;
            font-weight: 600 !important;
        }
        
        /* Color 2: Paste FASTA Sequence - Blue */
        label[data-baseweb="radio"]:has(input[value="1"]) {
            background: linear-gradient(135deg, #dbeafe 0%, #bfdbfe 100%) !important;
            border: 2px solid #3b82f6 !important;
        }
        label[data-baseweb="radio"]:has(input[value="1"]):hover {
            background: linear-gradient(135deg, #bfdbfe 0%, #93c5fd 100%) !important;
            box-shadow: 0 2px 8px rgba(59, 130, 246, 0.3) !important;
        }
        label[data-baseweb="radio"]:has(input[value="1"]) p {
            color: #1e40af !important;
            font-weight: 600 !important;
        }
        
        /* Color 3: Example Data - Amber/Yellow */
        label[data-baseweb="radio"]:has(input[value="2"]) {
            background: linear-gradient(135deg, #fef3c7 0%, #fde68a 100%) !important;
            border: 2px solid #f59e0b !important;
        }
        label[data-baseweb="radio"]:has(input[value="2"]):hover {
            background: linear-gradient(135deg, #fde68a 0%, #fcd34d 100%) !important;
            box-shadow: 0 2px 8px rgba(245, 158, 11, 0.3) !important;
        }
        label[data-baseweb="radio"]:has(input[value="2"]) p {
            color: #92400e !important;
            font-weight: 600 !important;
        }
        
        /* Color 4: NCBI Fetch - Purple */
        label[data-baseweb="radio"]:has(input[value="3"]) {
            background: linear-gradient(135deg, #ede9fe 0%, #ddd6fe 100%) !important;
            border: 2px solid #8b5cf6 !important;
        }
        label[data-baseweb="radio"]:has(input[value="3"]):hover {
            background: linear-gradient(135deg, #ddd6fe 0%, #c4b5fd 100%) !important;
            box-shadow: 0 2px 8px rgba(139, 92, 246, 0.3) !important;
        }
        label[data-baseweb="radio"]:has(input[value="3"]) p {
            color: #5b21b6 !important;
            font-weight: 600 !important;
        }
        
        /* ===== VIBRANT STYLING FOR BUTTONS AND FILE UPLOADER ===== */
        /* File uploader - Vibrant green */
        section[data-testid="stFileUploader"] {
            background: linear-gradient(135deg, #ecfdf5 0%, #d1fae5 100%) !important;
            border: 2px dashed #10b981 !important;
            border-radius: 12px !important;
            padding: 1rem !important;
        }
        section[data-testid="stFileUploader"]:hover {
            background: linear-gradient(135deg, #d1fae5 0%, #a7f3d0 100%) !important;
            border-color: #059669 !important;
            box-shadow: 0 4px 12px rgba(16, 185, 129, 0.25) !important;
        }
        section[data-testid="stFileUploader"] button {
            background: linear-gradient(135deg, #10b981 0%, #059669 100%) !important;
            color: white !important;
            border: none !important;
            font-weight: 700 !important;
        }
        section[data-testid="stFileUploader"] button:hover {
            background: linear-gradient(135deg, #34d399 0%, #10b981 100%) !important;
            box-shadow: 0 4px 12px rgba(16, 185, 129, 0.4) !important;
        }
        
        /* Load Example buttons - Vibrant amber/gold */
        button:has(+ div:empty) {
            font-weight: 700 !important;
        }
        </style>
        """, unsafe_allow_html=True)
        
        # Additional vibrant button styles for Load Example buttons
        st.markdown("""
        <style>
        /* Load Single Example and Load Multi-FASTA Example buttons */
        div[data-testid="stHorizontalBlock"] button[kind="secondary"] {
            background: linear-gradient(135deg, #fbbf24 0%, #f59e0b 100%) !important;
            color: #78350f !important;
            border: none !important;
            font-weight: 700 !important;
            box-shadow: 0 2px 8px rgba(245, 158, 11, 0.3) !important;
        }
        div[data-testid="stHorizontalBlock"] button[kind="secondary"]:hover {
            background: linear-gradient(135deg, #fcd34d 0%, #fbbf24 100%) !important;
            box-shadow: 0 4px 12px rgba(245, 158, 11, 0.4) !important;
            transform: translateY(-1px) !important;
        }
        </style>
        """, unsafe_allow_html=True)
        
        input_method = st.radio(UI_TEXT['upload_input_method_prompt'],
                            [UI_TEXT['upload_method_file'], UI_TEXT['upload_method_paste'], 
                             UI_TEXT['upload_method_example'], UI_TEXT['upload_method_ncbi']],
                            horizontal=False,
                            label_visibility="collapsed",
                            key="upload_method")

        seqs, names = [], []
        disk_seq_ids = []  # Sequences saved directly to disk during file upload parsing

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
                
                    # Calculate aggregate QC stats for compact summary using ATGC-weighted GC%
                    # GC% = (G+C)/(A+T+G+C)*100  — only ATGC bases in denominator
                    total_atgc = sum(
                        p['a_count'] + p['t_count'] + p['g_count'] + p['c_count']
                        for p in preview_info['previews']
                    )
                    total_gc_bases = sum(
                        p['g_count'] + p['c_count'] for p in preview_info['previews']
                    )
                    avg_gc = (total_gc_bases / total_atgc * 100) if total_atgc > 0 else 0.0
                    # Aggregate character counts for display
                    sum_a = sum(p['a_count'] for p in preview_info['previews'])
                    sum_t = sum(p['t_count'] for p in preview_info['previews'])
                    sum_g = sum(p['g_count'] for p in preview_info['previews'])
                    sum_c = sum(p['c_count'] for p in preview_info['previews'])
                    sum_n = sum(p['n_count'] for p in preview_info['previews'])
                
                    # Compact QC Summary Strip (replaces individual cards)
                    st.markdown(f"""
                    <div style='background: linear-gradient(135deg, #10b981 0%, #059669 100%); 
                                border-radius: 8px; padding: 10px 14px; margin: 8px 0; color: white;
                                box-shadow: 0 2px 6px rgba(16, 185, 129, 0.2);'>
                        <div style='font-weight: 600; font-size: 0.9rem; margin-bottom: 6px;'>
                            {fasta_file.name}
                        </div>
                        <div style='display: flex; gap: 12px; flex-wrap: wrap; font-size: 0.8rem;'>
                            <span style='background: rgba(255,255,255,0.2); padding: 4px 8px; border-radius: 4px;'>
                                <strong>{preview_info['num_sequences']}</strong> sequences
                            </span>
                            <span style='background: rgba(255,255,255,0.2); padding: 4px 8px; border-radius: 4px;'>
                                <strong>{preview_info['total_bp']:,}</strong> bp
                            </span>
                            <span style='background: rgba(255,255,255,0.2); padding: 4px 8px; border-radius: 4px;'>
                                GC: <strong>{avg_gc:.1f}%</strong>
                            </span>
                            <span style='background: rgba(255,255,255,0.2); padding: 4px 8px; border-radius: 4px;'>
                                A:<strong>{sum_a:,}</strong> T:<strong>{sum_t:,}</strong> G:<strong>{sum_g:,}</strong> C:<strong>{sum_c:,}</strong>{f' N:<strong>{sum_n:,}</strong>' if sum_n > 0 else ''}
                            </span>
                            <span style='background: rgba(255,255,255,0.2); padding: 4px 8px; border-radius: 4px;'>
                                ✓ Valid
                            </span>
                        </div>
                    </div>
                    """, unsafe_allow_html=True)
                
                    # Now parse all sequences using chunked parsing for memory efficiency
                    seqs, names = [], []
                    has_large_sequences = False
                    use_disk = (
                        st.session_state.get('use_disk_storage') and
                        st.session_state.get('seq_storage') is not None
                    )
                
                    if preview_info['num_sequences'] > 10:
                        # Show progress bar for files with many sequences
                        progress_bar = st.progress(0)
                        status_text = st.empty()
                    
                        for idx, (name, seq) in enumerate(parse_fasta_chunked(fasta_file)):
                            if len(seq) > 10_000_000:
                                has_large_sequences = True
                            if use_disk:
                                # Save directly to disk to avoid holding large sequences in RAM
                                seq_id = st.session_state.seq_storage.save_sequence(seq, name)
                                disk_seq_ids.append(seq_id)
                                names.append(name)
                                del seq  # Free RAM immediately
                            else:
                                names.append(name)
                                seqs.append(seq)
                        
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
                            if len(seq) > 10_000_000:
                                has_large_sequences = True
                            if use_disk:
                                # Save directly to disk to avoid holding large sequences in RAM
                                seq_id = st.session_state.seq_storage.save_sequence(seq, name)
                                disk_seq_ids.append(seq_id)
                                names.append(name)
                                del seq  # Free RAM immediately
                            else:
                                names.append(name)
                                seqs.append(seq)
                
                    # Force garbage collection after loading all sequences if we had large ones
                    if has_large_sequences:
                        gc.collect()
                
                    if not seqs and not disk_seq_ids:
                        st.warning(UI_TEXT['upload_no_sequences'])

        elif input_method == UI_TEXT['upload_method_paste']:
            seq_input = st.text_area(UI_TEXT['upload_paste_prompt'], 
                                    height=120, 
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
                        # Uppercase, strip whitespace, then filter non-IUPAC characters
                        cur_seq += _NON_IUPAC_RE.sub('', line.strip().upper())
                if cur_seq:
                    seqs.append(cur_seq)
                    names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                if seqs:
                    # Compact QC summary strip for pasted sequences
                    # GC% = (G+C)/(A+T+G+C)*100 — only ATGC in denominator
                    all_a = all_t = all_g = all_c = all_n = 0
                    for s in seqs:
                        a, t, g, c = _count_bases(s)
                        all_a += a; all_t += t; all_g += g; all_c += c
                        all_n += s.count('N')
                    total_atgc = all_a + all_t + all_g + all_c
                    avg_gc = (all_g + all_c) / total_atgc * 100 if total_atgc > 0 else 0.0
                    total_bp = sum(len(s) for s in seqs)
                    st.markdown(f"""
                    <div style='background: linear-gradient(135deg, #10b981 0%, #059669 100%); 
                                border-radius: 8px; padding: 10px 14px; margin: 8px 0; color: white;
                                box-shadow: 0 2px 6px rgba(16, 185, 129, 0.2);'>
                        <div style='font-weight: 600; font-size: 0.9rem; margin-bottom: 6px;'>
                            Pasted Input
                        </div>
                        <div style='display: flex; gap: 12px; flex-wrap: wrap; font-size: 0.8rem;'>
                            <span style='background: rgba(255,255,255,0.2); padding: 4px 8px; border-radius: 4px;'>
                                <strong>{len(seqs)}</strong> sequences
                            </span>
                            <span style='background: rgba(255,255,255,0.2); padding: 4px 8px; border-radius: 4px;'>
                                <strong>{total_bp:,}</strong> bp
                            </span>
                            <span style='background: rgba(255,255,255,0.2); padding: 4px 8px; border-radius: 4px;'>
                                GC: <strong>{avg_gc:.1f}%</strong>
                            </span>
                            <span style='background: rgba(255,255,255,0.2); padding: 4px 8px; border-radius: 4px;'>
                                A:<strong>{all_a:,}</strong> T:<strong>{all_t:,}</strong> G:<strong>{all_g:,}</strong> C:<strong>{all_c:,}</strong>{f' N:<strong>{all_n:,}</strong>' if all_n > 0 else ''}
                            </span>
                            <span style='background: rgba(255,255,255,0.2); padding: 4px 8px; border-radius: 4px;'>
                                ✓ Valid
                            </span>
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
                    fasta_text = _load_example_file("example_single.fasta", EXAMPLE_FASTA)
                    parsed_fasta = parse_fasta(fasta_text)
                    seqs = list(parsed_fasta.values())
                    names = list(parsed_fasta.keys())
                    st.success(UI_TEXT['upload_example_single_success'])
            else:
                if st.button("Load Multi-FASTA Example", use_container_width=True):
                    fasta_text = _load_example_file("example_multi.fasta", EXAMPLE_MULTI_FASTA)
                    seqs, names = [], []
                    cur_seq, cur_name = "", ""
                    for line in fasta_text.splitlines():
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

        # Persist sequences to session state if any found from input.
        # Guard: never overwrite analysis results during Streamlit reruns
        # (e.g. reruns triggered by download button clicks).  State is only
        # mutated when no analysis has been completed yet.
        _analysis_locked = st.session_state.get('analysis_done', False)
        if disk_seq_ids:
            if not _analysis_locked:
                # Sequences were saved directly to disk during file upload parsing
                st.session_state.seq_ids = disk_seq_ids
                st.session_state.names = names
                st.session_state.results_storage = {}
                # Keep legacy fields empty for compatibility
                st.session_state.seqs = []
                st.session_state.results = []
        elif seqs:
            # Support both disk storage and legacy in-memory mode
            if st.session_state.get('use_disk_storage') and st.session_state.get('seq_storage'):
                if not _analysis_locked:
                    # Use disk-based storage
                    seq_ids = []
                    for seq, name in zip(seqs, names):
                        seq_id = st.session_state.seq_storage.save_sequence(seq, name)
                        seq_ids.append(seq_id)
                    
                    st.session_state.seq_ids = seq_ids
                    st.session_state.names = names
                    st.session_state.results_storage = {}
                    
                    # Keep legacy fields for compatibility but with empty lists
                    st.session_state.seqs = []
                    st.session_state.results = []
            else:
                if not _analysis_locked:
                    # Legacy in-memory mode
                    st.session_state.seqs = seqs
                    st.session_state.names = names
                    st.session_state.results = []

        # Compact sequence validation summary (single strip, no individual cards)
        # Support both disk storage and legacy in-memory mode
        has_sequences = False
        num_sequences = 0
        total_bp = 0
        
        if st.session_state.get('use_disk_storage') and st.session_state.get('seq_ids'):
            # Disk storage mode: get info from metadata
            has_sequences = True
            num_sequences = len(st.session_state.seq_ids)
            for seq_id in st.session_state.seq_ids:
                metadata = st.session_state.seq_storage.get_metadata(seq_id)
                total_bp += metadata['length']
        elif st.session_state.get('seqs'):
            # Legacy mode: use in-memory sequences
            has_sequences = True
            num_sequences = len(st.session_state.seqs)
            total_bp = sum(len(s) for s in st.session_state.seqs)
        
        if has_sequences:
            # Build cache key from actual sequence names to prevent stale-cache GC% mismatch
            # when different genomes with the same sequence count are uploaded in the same session.
            seq_names = st.session_state.get('names') or []
            # Include count and all names in hash to prevent collisions
            names_sig = hashlib.md5(
                f"{num_sequences}|{'|'.join(seq_names)}".encode()
            ).hexdigest()[:8]
            if st.session_state.get('use_disk_storage') and st.session_state.get('seq_ids'):
                cache_key = f"stats_cache_disk_{names_sig}"
            else:
                cache_key = f"stats_cache_mem_{names_sig}"

            if cache_key not in st.session_state:
                # Calculate stats for all sequences
                stats_list = []
                if st.session_state.get('use_disk_storage') and st.session_state.get('seq_ids'):
                    # Disk storage mode: calculate from metadata
                    for seq_id in st.session_state.seq_ids:
                        metadata = st.session_state.seq_storage.get_metadata(seq_id)
                        base_counts = metadata.get('base_counts', {})
                        stats_list.append({
                            'Length': metadata['length'],
                            'GC%': metadata['gc_content'],
                            'A': base_counts.get('A', 0),
                            'T': base_counts.get('T', 0),
                            'G': base_counts.get('G', 0),
                            'C': base_counts.get('C', 0),
                            'N': base_counts.get('N', 0),
                        })
                else:
                    # Legacy mode: calculate from in-memory sequences
                    for seq in st.session_state.seqs:
                        stats = get_basic_stats(seq)
                        stats_list.append(stats)
                st.session_state[cache_key] = stats_list
        
            # Use cached stats - show single compact summary strip
            cached_stats = st.session_state[cache_key]
            # Length-weighted GC% = (total G+C) / (total A+T+G+C) * 100
            total_atgc = sum(
                s.get('A', 0) + s.get('T', 0) + s.get('G', 0) + s.get('C', 0)
                for s in cached_stats
            )
            total_gc_bases = sum(s.get('G', 0) + s.get('C', 0) for s in cached_stats)
            avg_gc = (total_gc_bases / total_atgc * 100) if total_atgc > 0 else 0.0
            # Aggregate character counts
            sum_a = sum(s.get('A', 0) for s in cached_stats)
            sum_t = sum(s.get('T', 0) for s in cached_stats)
            sum_g = sum(s.get('G', 0) for s in cached_stats)
            sum_c = sum(s.get('C', 0) for s in cached_stats)
            sum_n = sum(s.get('N', 0) for s in cached_stats)
            
            st.markdown(f"""
            <div style='background: linear-gradient(135deg, #059669 0%, #047857 100%); 
                        border-radius: 8px; padding: 10px 14px; margin: 8px 0; color: white;
                        box-shadow: 0 2px 6px rgba(5, 150, 105, 0.25);'>
                <div style='font-weight: 600; font-size: 0.85rem; margin-bottom: 6px;'>
                    ✓ Ready for Analysis
                </div>
                <div style='display: flex; gap: 10px; flex-wrap: wrap; font-size: 0.75rem;'>
                    <span style='background: rgba(255,255,255,0.2); padding: 4px 8px; border-radius: 4px;'>
                        <strong>{num_sequences}</strong> sequences
                    </span>
                    <span style='background: rgba(255,255,255,0.2); padding: 4px 8px; border-radius: 4px;'>
                        <strong>{total_bp:,}</strong> bp total
                    </span>
                    <span style='background: rgba(255,255,255,0.2); padding: 4px 8px; border-radius: 4px;'>
                        GC: <strong>{avg_gc:.1f}%</strong>
                    </span>
                    <span style='background: rgba(255,255,255,0.2); padding: 4px 8px; border-radius: 4px;'>
                        A:<strong>{sum_a:,}</strong> T:<strong>{sum_t:,}</strong> G:<strong>{sum_g:,}</strong> C:<strong>{sum_c:,}</strong>{f' N:<strong>{sum_n:,}</strong>' if sum_n > 0 else ''}
                    </span>
                </div>
            </div>
            """, unsafe_allow_html=True)
    
    # ============================================================
    # RIGHT COLUMN — DETECTION SCOPE
    # All 24 classes + submotifs visible via color-encoded grid
    # ============================================================
    with right_col:
        st.markdown("""
        <div style='background: linear-gradient(135deg, #ede9fe 0%, #ddd6fe 100%); 
                    padding: 0.5rem 1rem; border-radius: 6px; margin-bottom: 0.75rem;
                    border-left: 3px solid #8b5cf6;'>
            <h4 style='margin: 0; color: #5b21b6; font-size: 0.95rem; font-weight: 600;'>
                Detection Scope: Select Non-B DNA Motifs
            </h4>
        </div>
        """, unsafe_allow_html=True)
    
        # ============================================================
        # FIXED ANALYSIS GRANULARITY - Always Submotif Level (no UI)
        # ============================================================
        # Detection is always submotif-level to preserve biological accuracy
        # and prevent biologically invalid assumptions.
        st.session_state.analysis_mode = "Submotif Level"
        st.session_state.analysis_mode_used = "Submotif Level"
    
        # ============================================================
        # 24-CLASS SUBMOTIF SELECTOR (COLOR-ENCODED, HIGH-DENSITY GRID)
        # ============================================================
        # All 24 submotifs always visible via compact grid layout.
        # Color-encoding provides implicit class identity (via CLASS_COLORS).
        # No scrolling required - all detection targets visible at once.
        # User selects individual motifs from the 24 available options.
        # ============================================================

        # ------------------------------------------------------------------
        # Use CLASS_COLORS from config.colors for color encoding
        # (single source of truth for motif class colors)
        # ------------------------------------------------------------------

        def _sanitize_key(name: str) -> str:
            return name.replace(' ', '_').replace('-', '_').replace('/', '_')

        # ------------------------------------------------------------------
        # Build flat ordered list of all subclasses (taxonomy order)
        # ------------------------------------------------------------------
        flat_submotifs = []
        for class_id in sorted(MOTIF_CLASSIFICATION.keys()):
            entry = MOTIF_CLASSIFICATION[class_id]
            class_name = entry['class']
            for subclass in entry['subclasses']:
                flat_submotifs.append((class_name, subclass))

        # Initialize session state (all ON by default)
        for class_name, subclass in flat_submotifs:
            key = f"submotif_{_sanitize_key(class_name)}_{_sanitize_key(subclass)}"
            if key not in st.session_state:
                st.session_state[key] = True

        # ------------------------------------------------------------------
        # Render compact 3-column grid for full 24-class visibility
        # (3 columns × 8 rows = 24 submotifs visible without scrolling)
        # ------------------------------------------------------------------
        rows = [flat_submotifs[i:i + GRID_COLUMNS]
                for i in range(0, len(flat_submotifs), GRID_COLUMNS)]

        # Wrap grid in container with CSS class for targeted styling
        # Add colored dot CSS styling (replaces checkbox tick with colored dot)
        # Target only submotif checkboxes by their aria-label pattern to avoid affecting other checkboxes
        
        st.markdown(f"""
        <style>
        /* ===== DETECTION SCOPE - Colored Dots instead of Checkmarks ===== */
        /* Target only submotif checkboxes (those with color annotations in aria-label) */
        
        /* Make submotif checkboxes round (dots) - target by color annotation pattern */
        label[data-baseweb="checkbox"]:has(input[aria-label*=":blue["]) > span:first-child,
        label[data-baseweb="checkbox"]:has(input[aria-label*=":orange["]) > span:first-child,
        label[data-baseweb="checkbox"]:has(input[aria-label*=":red["]) > span:first-child,
        label[data-baseweb="checkbox"]:has(input[aria-label*=":violet["]) > span:first-child,
        label[data-baseweb="checkbox"]:has(input[aria-label*=":green["]) > span:first-child,
        label[data-baseweb="checkbox"]:has(input[aria-label*=":gray["]) > span:first-child {{
            border-radius: 50% !important;
            width: 14px !important;
            height: 14px !important;
            min-width: 14px !important;
            max-width: 14px !important;
            border-width: 2px !important;
            border-style: solid !important;
            background-color: transparent !important;
        }}
        
        /* Hide the checkmark SVG inside submotif checkboxes */
        label[data-baseweb="checkbox"]:has(input[aria-label*=":blue["]) > span:first-child svg,
        label[data-baseweb="checkbox"]:has(input[aria-label*=":orange["]) > span:first-child svg,
        label[data-baseweb="checkbox"]:has(input[aria-label*=":red["]) > span:first-child svg,
        label[data-baseweb="checkbox"]:has(input[aria-label*=":violet["]) > span:first-child svg,
        label[data-baseweb="checkbox"]:has(input[aria-label*=":green["]) > span:first-child svg,
        label[data-baseweb="checkbox"]:has(input[aria-label*=":gray["]) > span:first-child svg {{
            display: none !important;
            visibility: hidden !important;
            opacity: 0 !important;
            width: 0 !important;
            height: 0 !important;
        }}
        
        /* ===== Color-specific dot styling based on label color ===== */
        /* Blue colored labels (Curved_DNA, Z-DNA) */
        label[data-baseweb="checkbox"]:has(input[aria-label*=":blue["]) > span:first-child {{
            border-color: #06b6d4 !important;
        }}
        label[data-baseweb="checkbox"]:has(input[aria-label*=":blue["]:checked) > span:first-child {{
            background-color: #06b6d4 !important;
        }}
        
        /* Orange colored labels (Slipped_DNA, A-philic_DNA) */
        label[data-baseweb="checkbox"]:has(input[aria-label*=":orange["]) > span:first-child {{
            border-color: #f59e0b !important;
        }}
        label[data-baseweb="checkbox"]:has(input[aria-label*=":orange["]:checked) > span:first-child {{
            background-color: #f59e0b !important;
        }}
        
        /* Red colored labels (Cruciform) */
        label[data-baseweb="checkbox"]:has(input[aria-label*=":red["]) > span:first-child {{
            border-color: #ef4444 !important;
        }}
        label[data-baseweb="checkbox"]:has(input[aria-label*=":red["]:checked) > span:first-child {{
            background-color: #ef4444 !important;
        }}
        
        /* Violet colored labels (R-Loop, Triplex) */
        label[data-baseweb="checkbox"]:has(input[aria-label*=":violet["]) > span:first-child {{
            border-color: #8b5cf6 !important;
        }}
        label[data-baseweb="checkbox"]:has(input[aria-label*=":violet["]:checked) > span:first-child {{
            background-color: #8b5cf6 !important;
        }}
        
        /* Green colored labels (G-Quadruplex, i-Motif) */
        label[data-baseweb="checkbox"]:has(input[aria-label*=":green["]) > span:first-child {{
            border-color: #10b981 !important;
        }}
        label[data-baseweb="checkbox"]:has(input[aria-label*=":green["]:checked) > span:first-child {{
            background-color: #10b981 !important;
        }}
        
        /* Gray colored labels (Hybrid, Non-B_DNA_Clusters) */
        label[data-baseweb="checkbox"]:has(input[aria-label*=":gray["]) > span:first-child {{
            border-color: #64748b !important;
        }}
        label[data-baseweb="checkbox"]:has(input[aria-label*=":gray["]:checked) > span:first-child {{
            background-color: #64748b !important;
        }}
        </style>
        """, unsafe_allow_html=True)
        
        for row in rows:
            cols = st.columns(GRID_COLUMNS, gap="small")
            for col, (class_name, subclass) in zip(cols, row):
                # Use gray as fallback for unmapped classes (consistent with original)
                color = CLASS_COLORS.get(class_name, "#cbd5e1")
                color_name = CLASS_COLOR_NAMES.get(class_name, "gray")
                key = f"submotif_{_sanitize_key(class_name)}_{_sanitize_key(subclass)}"
                # Use abbreviated label for display, full name in tooltip
                abbrev_label = get_abbreviated_label(subclass)

                with col:
                    st.checkbox(
                        f"**:{color_name}[{abbrev_label}]**",
                        key=key,
                        help=f"{subclass} ({class_name.replace('_', ' ')})"
                    )

        # ------------------------------------------------------------------
        # Build enabled class & subclass lists (for downstream analysis)
        # ------------------------------------------------------------------
        enabled_classes = set()
        enabled_subclasses = []

        for class_name, subclass in flat_submotifs:
            key = f"submotif_{_sanitize_key(class_name)}_{_sanitize_key(subclass)}"
            if st.session_state.get(key, True):
                enabled_classes.add(class_name)
                enabled_subclasses.append(subclass)

        st.session_state.selected_classes = list(enabled_classes)
        st.session_state.selected_subclasses = enabled_subclasses

        # Also update motif_selector_data for any code that uses it (in taxonomy order)
        st.session_state.motif_selector_data = []
        for class_id in sorted(MOTIF_CLASSIFICATION.keys()):
            entry = MOTIF_CLASSIFICATION[class_id]
            class_name = entry['class']
            for subclass in entry['subclasses']:
                key = f"submotif_{_sanitize_key(class_name)}_{_sanitize_key(subclass)}"
                st.session_state.motif_selector_data.append({
                    'Enabled': st.session_state.get(key, True),
                    'Motif Class': class_name,
                    'Submotif': subclass
                })

        # Compact summary with publication-grade terminology
        num_enabled = len(enabled_subclasses)
        total_submotifs = len(flat_submotifs)

        if enabled_classes:
            st.markdown(f"""
            <div style='background: linear-gradient(135deg, #f5f3ff 0%, #ede9fe 100%);
                        padding: 6px 10px; border-radius: 6px; margin-top: 8px;
                        border: 1px solid #c4b5fd; font-size: 0.78rem;'>
                <span style='font-weight: 600; color: #6d28d9;'>
                    {len(enabled_classes)} classes · {num_enabled}/{total_submotifs} submotifs
                </span>
            </div>
            """, unsafe_allow_html=True)
        else:
            st.warning("Select at least one submotif to enable analysis.")
        
        # ============================================================
        # RUN AND RESET BUTTONS - Inside Detection Scope Section
        # ============================================================
        # Initialize analysis_done flag if not present (idempotent run button)
        if "analysis_done" not in st.session_state:
            st.session_state.analysis_done = False
        
        st.markdown("<div style='margin-top: 0.5rem;'></div>", unsafe_allow_html=True)
        
        # Check if valid input is present (sequences AND at least one class selected)
        # Support both disk storage and legacy in-memory mode
        has_sequences_inner = bool(st.session_state.get('seqs')) or bool(st.session_state.get('seq_ids'))
        has_valid_input_inner = has_sequences_inner and bool(enabled_classes)
        
        # Create two equal columns for Run and Reset buttons (parallel layout)
        col_run_inner, col_reset_inner = st.columns([1, 1])

        # Vibrant Reset button styling — scoped to the reset column only
        # NOTE: CSS uses Streamlit's stable data-testid attributes (baseButton-secondary, column, stHorizontalBlock).
        # These are consistent across Streamlit >=1.28. Re-verify if upgrading to a major Streamlit version.
        st.markdown("""
        <style>
        /* Reset button: vibrant red gradient, scoped to last column in button row */
        div[data-testid="stHorizontalBlock"] div[data-testid="column"]:last-child button[data-testid="baseButton-secondary"] {
            background: linear-gradient(135deg, #ef4444 0%, #dc2626 100%) !important;
            color: white !important;
            border: none !important;
            font-weight: 700 !important;
            font-size: 1rem !important;
            letter-spacing: 0.02em !important;
            box-shadow: 0 2px 8px rgba(220, 38, 38, 0.35) !important;
            transition: all 0.2s ease !important;
        }
        div[data-testid="stHorizontalBlock"] div[data-testid="column"]:last-child button[data-testid="baseButton-secondary"]:hover {
            background: linear-gradient(135deg, #f87171 0%, #ef4444 100%) !important;
            box-shadow: 0 4px 16px rgba(220, 38, 38, 0.55) !important;
            transform: translateY(-1px) !important;
        }
        div[data-testid="stHorizontalBlock"] div[data-testid="column"]:last-child button[data-testid="baseButton-secondary"]:active {
            transform: translateY(0) !important;
            box-shadow: 0 2px 6px rgba(220, 38, 38, 0.4) !important;
        }
        </style>
        """, unsafe_allow_html=True)

        with col_run_inner:
            # Primary action button with clear scientific terminology
            if has_valid_input_inner and not st.session_state.analysis_done:
                run_button = st.button(
                    "🧬 Run Non-B DNA Analysis",
                    type="primary",
                    use_container_width=True,
                    key="run_motif_analysis_main",
                    help="Start analyzing uploaded sequences for Non-B DNA elements"
                )
            elif st.session_state.analysis_done:
                # Show analysis complete status
                st.markdown("""
                <div role="status" aria-label="Analysis Complete"
                     style='background: linear-gradient(135deg, #10b981 0%, #059669 100%); 
                            color: white; padding: 12px; 
                            border-radius: 8px; text-align: center; font-weight: 600;
                            font-size: 0.95rem;'>
                    ✓ Analysis Complete — View results in 'Results' tab
                </div>
                """, unsafe_allow_html=True)
                run_button = False
            else:
                # Disabled button appearance with accessibility
                st.markdown(f"""
                <div role="button" aria-disabled="true" aria-label="Run Analysis - Disabled"
                     style='background: #e5e7eb; color: #9ca3af; padding: 12px; 
                            border-radius: 8px; text-align: center; font-weight: 600;
                            font-size: 0.95rem; cursor: not-allowed;'>
                    🧬 Run Non-B DNA Analysis
                </div>
                <p style='text-align: center; color: #9ca3af; font-size: 0.75rem; margin-top: 4px;'>
                    Upload sequences and select submotifs to enable
                </p>
                """, unsafe_allow_html=True)
                run_button = False

        with col_reset_inner:
            # Reset button — clears all uploaded data, sequences, and results
            if st.button("🔄 Reset", use_container_width=True, help="Clear all uploaded files, sequences, results, and analysis state to start over", key="reset_button_inner"):
                # Clear analysis flags and computed results
                st.session_state.analysis_done = False
                st.session_state.results = []
                st.session_state.performance_metrics = None
                st.session_state.cached_visualizations = {}
                st.session_state.analysis_time = None
                # Clear loaded sequences and associated data
                st.session_state.seqs = []
                st.session_state.names = []
                st.session_state.seq_ids = []
                st.session_state.results_storage = {}
                st.session_state.summary_df = pd.DataFrame()
                st.session_state.analysis_status = "Ready"
                st.session_state.current_job_id = None
                st.session_state.selected_classes_used = []
                st.session_state.selected_subclasses_used = []
                # Free cached export bytes and per-sequence stats caches
                st.cache_data.clear()
                for _k in [k for k in list(st.session_state.keys()) if k.startswith('stats_cache_')]:
                    st.session_state.pop(_k, None)
                st.rerun()

        # ============================================================
        # PROGRESS DISPLAY AREA - Within Detection Scope
        # ============================================================
        # Create a vibrant container for progress display below the buttons
        with st.container():
            st.markdown("""
            <div style='background: linear-gradient(135deg, #f0f9ff 0%, #e0f2fe 100%); 
                        padding: 0.75rem 1rem; border-radius: 8px; margin-top: 1rem;
                        border-left: 4px solid #0ea5e9;
                        box-shadow: 0 2px 8px rgba(14, 165, 233, 0.15);'>
                <div style='color: #0369a1; font-weight: 600; font-size: 0.9rem;'>
                    📊 Analysis Progress
                </div>
            </div>
            """, unsafe_allow_html=True)
            
            # Create placeholders for progress tracking within the container
            progress_placeholder = st.empty()
            status_placeholder = st.empty()
            detailed_progress_placeholder = st.empty()
            timer_placeholder = st.empty()

    
        # ============================================================
        # ANALYSIS OPTIONS - Always ON, Hidden from UI
        # ============================================================
        # All engine-level options are enforced internally and never shown.
        # Validation always runs, parallelization enabled where applicable,
        # results are reproducible and reviewer-safe.
        # ============================================================
    
        # Enforce all analysis options as always ON (hard defaults, no UI)
        st.session_state.toggle_detailed = True
        st.session_state.toggle_validation = True
        st.session_state.toggle_parallel = True
        st.session_state.toggle_chunk_progress = True
        st.session_state.toggle_memory = True
    
        # Expose variables for use in analysis logic below
        detailed_output = True
        quality_check = True
        use_parallel_scanner = True
        show_chunk_progress = True
        show_memory_usage = True
    
        # Hardcoded default overlap handling
        nonoverlap = True
        overlap_option = "Remove overlaps within subclasses"
    
    # ============================================================
    # ANALYSIS EXECUTION - Variable Setup
    # ============================================================
    
    # Initialize analysis_done flag if not present (idempotent run button)
    if "analysis_done" not in st.session_state:
        st.session_state.analysis_done = False
    
    # Check if valid input is present (sequences AND at least one class selected)
    # Support both disk storage and legacy in-memory mode
    has_sequences = bool(st.session_state.get('seqs')) or bool(st.session_state.get('seq_ids'))
    has_valid_input = has_sequences and bool(st.session_state.get('selected_classes'))
    
    # Note: Progress placeholders are created in Detection Scope section (right column)
    # and will be used by the analysis execution logic below
    
    # ========== RUN ANALYSIS BUTTON LOGIC ========== 
    # Only run if button clicked AND not already done (idempotent)
    if run_button and not st.session_state.analysis_done:
        # Simplified validation - support both disk storage and legacy in-memory mode
        has_sequences = bool(st.session_state.get('seqs')) or bool(st.session_state.get('seq_ids'))
        if not has_sequences:
            st.error("Please upload or input sequences before running analysis.")
            st.session_state.analysis_status = "Error"
        elif not st.session_state.get('selected_classes'):
            st.error("Please select at least one motif class before running analysis.")
            st.session_state.analysis_status = "Error"
        else:
            # ============================================================
            # JOB ID GENERATION: Create unique ID for internal tracking
            # ============================================================
            # Generate job ID at the start of analysis (or reuse if already exists)
            # This prevents overwriting existing job IDs on multiple button clicks
            if not st.session_state.get('current_job_id'):
                job_id = generate_job_id()
                st.session_state.current_job_id = job_id
            else:
                job_id = st.session_state.current_job_id
            
            # Sequence length limit has been removed - the system now uses automatic chunking
            # (see NonBFinder.py CHUNK_THRESHOLD=10,000 bp) to handle sequences of any size
            st.session_state.analysis_status = "Running"
            # Clear stale cached export data so fresh results are generated
            st.cache_data.clear()
            
            # Store analysis parameters in session state for use in download section
            st.session_state.overlap_option_used = overlap_option
            st.session_state.nonoverlap_used = nonoverlap
            # analysis_mode_used is already set at the top (fixed as "Submotif Level")
            st.session_state.selected_classes_used = list(st.session_state.selected_classes)  # Store selected classes
            st.session_state.selected_subclasses_used = list(st.session_state.selected_subclasses)  # Store selected subclasses
            
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
            
            # ============================================================
            # PERFORMANCE TRACKING: Initialize tracker for statistics
            # ============================================================
            perf_tracker = PerformanceTracker()
            perf_tracker.start()
            
            # Enhanced progress tracking - timing captured exactly once at start and end
            import time
            
            # Note: progress placeholders are now created within the Detection Scope section
            # (see lines ~1037-1055) to display progress in the right column
            
            # ============================================================
            # DETERMINISTIC TIMING: Start time captured exactly once
            # ============================================================
            start_time = time.time()
            
            # Simple status message instead of toast
            status_placeholder.info("🧬 Starting Non-B DNA Analysis...")

            
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
            
            # Support both disk storage and legacy in-memory mode
            if st.session_state.get('use_disk_storage') and st.session_state.get('seq_ids'):
                # Disk storage mode: calculate from metadata
                total_bp_all_sequences = 0
                for seq_id in st.session_state.seq_ids:
                    metadata = st.session_state.seq_storage.get_metadata(seq_id)
                    total_bp_all_sequences += metadata['length']
                num_sequences = len(st.session_state.seq_ids)
            else:
                # Legacy mode: use in-memory sequences
                total_bp_all_sequences = sum(len(seq) for seq in st.session_state.seqs)
                num_sequences = len(st.session_state.seqs)
            
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
                
                # ============================================================
                # MULTI-FASTA STABILITY: Per-sequence logic isolated
                # ============================================================
                # Each sequence processed independently with clean isolation.
                # No shared state between iterations except cumulative counters.
                # Results accumulated in list, stored atomically after loop.
                # Ensures identical behavior for single and multi-FASTA inputs.
                # ============================================================
                
                # Check if we should use parallel processing (2+ sequences)
                use_parallel = num_sequences >= PARALLEL_PROCESSING_THRESHOLD
                
                if use_parallel:
                    # ============================================================
                    # PARALLEL PROCESSING MODE - Multi-FASTA Optimization
                    # ============================================================
                    status_placeholder.info(f"🚀 Using parallel processing for {num_sequences} sequences")
                    
                    # Import helper functions
                    from Utilities.parallel_analysis_helper import (
                        prepare_parallel_analysis,
                        get_optimal_workers
                    )
                    
                    # Prepare data for parallel processing
                    if st.session_state.get('use_disk_storage') and st.session_state.get('seq_ids'):
                        sequences_data, analysis_params = prepare_parallel_analysis(
                            num_sequences=num_sequences,
                            use_disk_storage=True,
                            seq_ids=st.session_state.seq_ids,
                            names=st.session_state.names,
                            seq_storage=st.session_state.seq_storage,
                            enabled_classes=list(analysis_classes) if analysis_classes else None,
                            chunk_threshold=CHUNK_ANALYSIS_THRESHOLD_BP
                        )
                    else:
                        sequences_data, analysis_params = prepare_parallel_analysis(
                            num_sequences=num_sequences,
                            use_disk_storage=False,
                            seqs=st.session_state.seqs,
                            names=st.session_state.names,
                            enabled_classes=list(analysis_classes) if analysis_classes else None,
                            chunk_threshold=CHUNK_ANALYSIS_THRESHOLD_BP
                        )
                    
                    # Progress tracking
                    progress_status = st.empty()
                    
                    # Progress callback closure - captures progress_placeholder and pbar from outer scope
                    # These must be initialized before this function is defined
                    def parallel_progress_callback(completed, total, seq_name):
                        """Update progress bar and status. Closure depends on progress_placeholder and pbar."""
                        progress = completed / total
                        with progress_placeholder.container():
                            pbar.progress(progress, text=f"Analyzed {completed}/{total} sequences")
                        progress_status.info(f"⏳ Completed: {seq_name}")
                    
                    # Run parallel analysis
                    with st.status("🧬 Running parallel sequence analysis...", expanded=True) as parallel_status:
                        st.write(f"Analyzing {num_sequences} sequences in parallel...")
                        
                        parallel_results = analyze_sequences_parallel(
                            sequences_data=sequences_data,
                            analysis_params=analysis_params,
                            max_workers=get_optimal_workers(num_sequences),
                            progress_callback=parallel_progress_callback,
                            use_processes=False  # Use threads for I/O-bound operations
                        )
                        
                        parallel_status.update(label="✅ Parallel analysis complete", state="complete")
                    
                    progress_status.empty()
                    
                    # Process parallel results
                    for result in parallel_results:
                        if result['success']:
                            i = result['index']
                            name = result['name']
                            results = result['results']
                            seq_length = result['seq_length']
                            
                            # Store results storage if using disk mode
                            if result.get('use_disk_storage') and result.get('results_storage'):
                                seq_id = result['seq_id']
                                st.session_state.results_storage[seq_id] = result['results_storage']
                            
                            # Update performance tracker
                            perf_tracker.add_sequence_result(name, seq_length, result['elapsed'], len(results))
                            
                            # Apply filters (same as sequential mode)
                            safe_results = []
                            for idx, motif in enumerate(results):
                                try:
                                    safe_results.append(ensure_subclass(motif))
                                except (KeyError, AttributeError, TypeError) as e:
                                    logger.error(f"Error processing motif at index {idx}: {e}")
                                    continue
                            results = safe_results
                            
                            # Filter by selected subclasses
                            selected_classes_set = set(st.session_state.selected_classes)
                            selected_subclasses_set = set(st.session_state.selected_subclasses)
                            
                            filtered_results = []
                            for idx, motif in enumerate(results):
                                try:
                                    motif_class = motif.get('Class', '')
                                    motif_subclass = motif.get('Subclass', '')
                                    
                                    if motif_class in selected_classes_set:
                                        if motif_class in ['Hybrid', 'Non-B_DNA_Clusters']:
                                            component_classes = motif.get('Component_Classes', [])
                                            if not component_classes or any(c in selected_classes_set for c in component_classes):
                                                filtered_results.append(motif)
                                        elif motif_subclass in selected_subclasses_set:
                                            filtered_results.append(motif)
                                except (KeyError, AttributeError, TypeError) as e:
                                    logger.error(f"Error filtering motif at index {idx}: {e}")
                                    continue
                            
                            results = filtered_results
                            all_results.append(results)
                            total_bp_processed += seq_length
                            
                            # Brief status update without popup
                            status_placeholder.success(f"✅ {name}: {len(results):,} motifs")
                        else:
                            st.error(f"❌ Failed: {result['name']} - {result.get('error')}")
                            all_results.append([])
                    
                    # Memory management
                    trigger_garbage_collection()
                    
                elif st.session_state.get('use_disk_storage') and st.session_state.get('seq_ids'):
                    # === DISK STORAGE MODE: Use ChunkAnalyzer ===
                    seq_ids = st.session_state.seq_ids
                    names = st.session_state.names
                    
                    for i, (seq_id, name) in enumerate(zip(seq_ids, names)):
                        progress = (i + 1) / num_sequences
                        
                        # Get sequence metadata
                        metadata = st.session_state.seq_storage.get_metadata(seq_id)
                        seq_length = metadata['length']
                        
                        # Calculate overall percentage
                        overall_percentage = (total_bp_processed / total_bp_all_sequences * 100) if total_bp_all_sequences > 0 else 0
                        
                        # Use st.status() for real-time disappearing progress
                        with st.status(f"🧬 Analyzing sequence {i+1}/{num_sequences}: {name} ({seq_length:,} bp)", expanded=True) as status:
                            seq_start_time = time.time()
                            
                            # Step 1: Validation
                            st.write("✓ Validating sequence...")
                            
                            # Determine whether to use chunking based on sequence size
                            if seq_length > CHUNK_ANALYSIS_THRESHOLD_BP:
                                # Step 2: Chunking
                                st.write(f"✓ Using chunk-based analysis (50KB chunks)")
                                status_placeholder.info(f"📦 Chunking large sequence: {name}")
                                
                                # Create chunk progress callback with status updates
                                chunk_progress = {'current': 0, 'total': 0}
                                last_update_pct = 0
                                
                                def chunk_callback(progress_pct):
                                    nonlocal last_update_pct
                                    chunk_progress['current'] = int(progress_pct)
                                    # Update every 20% to avoid too many updates
                                    if progress_pct - last_update_pct >= 20:
                                        st.write(f"  ⏳ Chunk progress: {progress_pct:.0f}%")
                                        last_update_pct = progress_pct
                                
                                # Step 3: Analysis
                                st.write("⏳ Running detectors on chunks...")
                                chunk_start = time.time()
                                
                                # Analyze using ChunkAnalyzer
                                analyzer = ChunkAnalyzer(
                                    st.session_state.seq_storage,
                                    use_parallel=True,
                                    max_workers=None,
                                    use_adaptive=False
                                )
                                
                                results_storage = analyzer.analyze(
                                    seq_id=seq_id,
                                    progress_callback=chunk_callback,
                                    enabled_classes=list(st.session_state.selected_classes) if st.session_state.selected_classes else None
                                )
                                
                                chunk_elapsed = time.time() - chunk_start
                                perf_tracker.add_chunk_time(chunk_elapsed)
                                
                                # Store results storage
                                st.session_state.results_storage[seq_id] = results_storage
                                
                                # Get results for filtering - Load ALL results (no limit)
                                results = list(results_storage.iter_results())
                                # Get actual total count from summary stats
                                total_motifs = results_storage.get_summary_stats()['total_count']
                                st.write(f"✓ Chunk analysis complete: {total_motifs:,} total motifs found")
                            else:
                                # Small sequence: use standard analysis
                                st.write("⏳ Running detectors...")
                                
                                # Get full sequence for standard analysis
                                seq = st.session_state.seq_storage.get_sequence_chunk(seq_id, 0, seq_length)
                                
                                # Validate sequence before analysis
                                if not seq or len(seq) == 0:
                                    st.error(f"❌ Empty sequence detected")
                                    status.update(label=f"❌ Failed: {name}", state="error")
                                    continue
                                
                                # Use standard analysis
                                analysis_start = time.time()
                                results = analyze_sequence(
                                    seq, name,
                                    enabled_classes=list(st.session_state.selected_classes) if st.session_state.selected_classes else None
                                )
                                analysis_elapsed = time.time() - analysis_start
                                
                                # Create results storage and save
                                results_storage = UniversalResultsStorage(
                                    base_dir=str(st.session_state.seq_storage.base_dir / "results"),
                                    seq_id=seq_id
                                )
                                results_storage.append_batch(results)
                                st.session_state.results_storage[seq_id] = results_storage
                                
                                st.write(f"✓ Analysis complete: {len(results):,} motifs found")
                            
                            # Calculate sequence stats
                            seq_elapsed = time.time() - seq_start_time
                            throughput = seq_length / seq_elapsed if seq_elapsed > 0 else 0
                            
                            # Update performance tracker
                            perf_tracker.add_sequence_result(name, seq_length, seq_elapsed, len(results))
                            
                            # Step 4: Filtering
                            st.write("⏳ Applying filters...")
                        
                            # Ensure all motifs have required fields - safe iteration with error handling
                            safe_results = []
                            for idx, motif in enumerate(results):
                                try:
                                    safe_results.append(ensure_subclass(motif))
                                except (KeyError, AttributeError, TypeError) as e:
                                    logger.error(f"Error processing motif at index {idx}: {e}")
                                    continue  # Skip problematic motif instead of crashing
                            results = safe_results
                            
                            # Filter results based on selected subclasses
                            selected_classes_set = set(st.session_state.selected_classes)
                            selected_subclasses_set = set(st.session_state.selected_subclasses)
                            
                            filtered_results = []
                            for idx, motif in enumerate(results):
                                try:
                                    motif_class = motif.get('Class', '')
                                    motif_subclass = motif.get('Subclass', '')
                                    
                                    if motif_class in selected_classes_set:
                                        if motif_class in ['Hybrid', 'Non-B_DNA_Clusters']:
                                            component_classes = motif.get('Component_Classes', [])
                                            if not component_classes or any(c in selected_classes_set for c in component_classes):
                                                filtered_results.append(motif)
                                        elif motif_subclass in selected_subclasses_set:
                                            filtered_results.append(motif)
                                except (KeyError, AttributeError, TypeError) as e:
                                    logger.error(f"Error filtering motif at index {idx}: {e}")
                                    continue  # Skip problematic motif instead of crashing
                            
                            results = filtered_results
                            st.write(f"✓ Filters applied: {len(results):,} motifs retained")
                            
                            # Final summary
                            st.write(f"⏱️ Time: {format_time(seq_elapsed)} | ⚡ Throughput: {throughput:,.0f} bp/s")
                            
                            # Update status to complete
                            status.update(label=f"✅ Completed: {name} ({len(results):,} motifs)", state="complete")
                            
                            # Brief status update without popup
                            status_placeholder.success(f"✅ {name}: {len(results):,} motifs found")
                        
                        all_results.append(results)
                        
                        total_bp_processed += seq_length
                        
                        # Update progress display
                        actual_percentage = (total_bp_processed / total_bp_all_sequences * 100) if total_bp_all_sequences > 0 else 0
                        with progress_placeholder.container():
                            pbar.progress(progress, text=f"Analyzed {i+1}/{num_sequences} sequences ({actual_percentage:.1f}%)")
                        
                        # Memory management
                        trigger_garbage_collection()
                    
                else:
                    # === LEGACY IN-MEMORY MODE ===
                    for i, (seq, name) in enumerate(zip(st.session_state.seqs, st.session_state.names)):
                        progress = (i + 1) / num_sequences
                        
                        # Calculate overall percentage
                        overall_percentage = (total_bp_processed / total_bp_all_sequences * 100) if total_bp_all_sequences > 0 else 0
                        
                        # Use st.status() for real-time disappearing progress
                        with st.status(f"🧬 Analyzing sequence {i+1}/{num_sequences}: {name} ({len(seq):,} bp)", expanded=True) as status:
                            seq_start_time = time.time()
                            
                            # Step 1: Validation
                            st.write("✓ Validating sequence...")
                            
                            if not seq or len(seq) == 0:
                                st.error(f"❌ Empty sequence detected")
                                status.update(label=f"❌ Failed: {name}", state="error")
                                continue
                            
                            # Step 2: Analysis
                            st.write("⏳ Running detectors...")
                            analysis_start = time.time()
                            
                            # Use standard consolidated NBDScanner analysis
                            results = analyze_sequence(
                                seq, name,
                                enabled_classes=list(st.session_state.selected_classes) if st.session_state.selected_classes else None
                            )
                            
                            analysis_elapsed = time.time() - analysis_start
                            st.write(f"✓ Analysis complete: {len(results):,} motifs found")
                            
                            # Step 3: Filtering
                            st.write("⏳ Applying filters...")
                            
                            # Ensure all motifs have required fields - safe iteration with error handling
                            safe_results = []
                            for idx, motif in enumerate(results):
                                try:
                                    safe_results.append(ensure_subclass(motif))
                                except (KeyError, AttributeError, TypeError) as e:
                                    logger.error(f"Error processing motif at index {idx}: {e}")
                                    continue  # Skip problematic motif instead of crashing
                            results = safe_results
                            
                            # Filter results based on selected subclasses
                            selected_classes_set = set(st.session_state.selected_classes)
                            selected_subclasses_set = set(st.session_state.selected_subclasses)
                            
                            filtered_results = []
                            for idx, motif in enumerate(results):
                                try:
                                    motif_class = motif.get('Class', '')
                                    motif_subclass = motif.get('Subclass', '')
                                    
                                    if motif_class in selected_classes_set:
                                        if motif_class in ['Hybrid', 'Non-B_DNA_Clusters']:
                                            component_classes = motif.get('Component_Classes', [])
                                            if not component_classes or any(c in selected_classes_set for c in component_classes):
                                                filtered_results.append(motif)
                                        elif motif_subclass in selected_subclasses_set:
                                            filtered_results.append(motif)
                                except (KeyError, AttributeError, TypeError) as e:
                                    logger.error(f"Error filtering motif at index {idx}: {e}")
                                    continue  # Skip problematic motif instead of crashing
                            
                            results = filtered_results
                            st.write(f"✓ Filters applied: {len(results):,} motifs retained")
                            
                            # Calculate stats
                            seq_elapsed = time.time() - seq_start_time
                            throughput = len(seq) / seq_elapsed if seq_elapsed > 0 else 0
                            
                            # Update performance tracker
                            perf_tracker.add_sequence_result(name, len(seq), seq_elapsed, len(results))
                            
                            # Final summary
                            st.write(f"⏱️ Time: {format_time(seq_elapsed)} | ⚡ Throughput: {throughput:,.0f} bp/s")
                            
                            # Update status to complete
                            status.update(label=f"✅ Completed: {name} ({len(results):,} motifs)", state="complete")
                            
                            # Brief status update without popup
                            status_placeholder.success(f"✅ {name}: {len(results):,} motifs found")
                        
                        all_results.append(results)
                    
                        total_bp_processed += len(seq)
                        
                        # Update progress display
                        actual_percentage = (total_bp_processed / total_bp_all_sequences * 100) if total_bp_all_sequences > 0 else 0
                        with progress_placeholder.container():
                            pbar.progress(progress, text=f"Analyzed {i+1}/{len(st.session_state.seqs)} sequences ({actual_percentage:.1f}%)")
                        
                        # Memory management
                        if len(seq) > 1_000_000:
                            trigger_garbage_collection()
                
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
                status_placeholder.info("Validating results for consistency and quality...")
                
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
                        if len(motifs) < 2:
                            continue  # Skip if less than 2 motifs
                        sorted_motifs = sorted(motifs, key=lambda m: m.get('Start', 0))
                        for j in range(len(sorted_motifs) - 1):
                            if sorted_motifs[j].get('End', 0) > sorted_motifs[j+1].get('Start', 0):
                                validation_issues.append(f"Note: Sequence {i+1}: Overlapping motifs in {subclass}")
                                break
                
                # Display validation results (ephemeral - replaces previous)
                if validation_issues:
                    warning_msg = UI_TEXT['status_validation_issues'].format(count=len(validation_issues)) + "\n"
                    for issue in validation_issues[:5]:  # Show first 5
                        warning_msg += f"\n• {issue}"
                    if len(validation_issues) > 5:
                        warning_msg += f"\n• ... and {len(validation_issues) - 5} more"
                    status_placeholder.warning(warning_msg)
                else:
                    status_placeholder.success(UI_TEXT['status_validation_passed'])
                
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
                    # Support both disk storage and legacy in-memory modes
                    if st.session_state.get('use_disk_storage') and st.session_state.get('seq_ids'):
                        # Disk storage mode: get metadata
                        seq_id = st.session_state.seq_ids[i]
                        metadata = st.session_state.seq_storage.get_metadata(seq_id)
                        name = st.session_state.names[i]
                        seq_length = metadata['length']
                        gc_content = metadata['gc_content']
                    else:
                        # Legacy in-memory mode
                        seq = st.session_state.seqs[i]
                        name = st.session_state.names[i]
                        stats = get_basic_stats(seq, results)
                        seq_length = stats['Length']
                        gc_content = stats['GC%']
                    
                    summary.append({
                        'Sequence': name,
                        'Length': seq_length,
                        'GC Content': f"{gc_content:.1f}%",
                        'Motifs Found': len(results),
                        'Unique Types': len(set(m.get('Type', 'Unknown') for m in results)),
                        'Avg Score': f"{np.mean([m.get('Score', 0) for m in results]):.3f}" if results else "0.000"
                    })
                
                # ATOMIC STORAGE: Store summary once in session state
                st.session_state.summary_df = pd.DataFrame(summary)
                
                # Display summary statistics to user
                if not st.session_state.summary_df.empty:
                    st.subheader("📊 Analysis Summary Statistics")
                    st.dataframe(
                        st.session_state.summary_df,
                        use_container_width=True,
                        hide_index=True
                    )
                
                # ============================================================
                # PRE-GENERATE ALL VISUALIZATIONS FOR CLASSES AND SUBCLASSES
                # ============================================================
                # Ephemeral info (replaces previous)
                status_placeholder.info("Generating comprehensive visualizations for all classes and subclasses...")
                
                # Cache all visualizations for each sequence
                st.session_state.cached_visualizations = {}
                
                viz_start_time = time.time()
                total_viz_count = 0
                
                # Reduce UI updates by batching - only update every N sequences or at end
                # Support both storage modes for UPDATE_INTERVAL calculation
                num_seqs = len(st.session_state.seq_ids) if st.session_state.get('use_disk_storage') and st.session_state.get('seq_ids') else len(st.session_state.seqs)
                UPDATE_INTERVAL = max(1, num_seqs // 5)  # Update 5 times max
                
                for seq_idx, results in enumerate(all_results):
                    # Support both storage modes
                    if st.session_state.get('use_disk_storage') and st.session_state.get('seq_ids'):
                        seq_id = st.session_state.seq_ids[seq_idx]
                        metadata = st.session_state.seq_storage.get_metadata(seq_id)
                        sequence_length = metadata['length']
                        name = st.session_state.names[seq_idx]
                    else:
                        seq = st.session_state.seqs[seq_idx]
                        sequence_length = len(seq)
                        name = st.session_state.names[seq_idx]
                    
                    # Continue with existing motif filtering logic
                    filtered_motifs = results
                    
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

                    # ----------------------------------------------------------
                    # Pre-generate all matplotlib figures so that:
                    #   1. viz_total_time accurately includes rendering time.
                    #   2. The Results tab renders instantly from cached bytes.
                    # ----------------------------------------------------------
                    fig_bytes = {}
                    _has_clusters = any(m.get('Class') == 'Non-B_DNA_Clusters' for m in filtered_motifs)
                    _has_hybrids  = any(m.get('Class') == 'Hybrid' for m in filtered_motifs)
                    _sm = [m for m in filtered_motifs if m.get('Class') not in ['Hybrid', 'Non-B_DNA_Clusters']]

                    _plot_tasks = [
                        ('class_track',       lambda: plot_linear_motif_track(filtered_motifs, sequence_length, title="Class Track")),
                    ]
                    if _sm:
                        _plot_tasks.append(('subclass_track', lambda: plot_linear_subclass_track(_sm, sequence_length, title="Subclass Track")))
                    _plot_tasks += [
                        ('class_dist',        lambda: plot_motif_distribution(filtered_motifs, by='Class', title="Class Distribution")),
                        ('subclass_dist',     lambda: plot_motif_distribution(filtered_motifs, by='Subclass', title="Subclass Distribution")),
                        ('length_kde',        lambda: plot_motif_length_kde(filtered_motifs, by_class=True, title="Length Distribution")),
                        ('score_violin',      lambda: plot_score_violin(filtered_motifs, by_class=True, title="Score Distribution")),
                        ('nested_pie',        lambda: plot_nested_pie_chart(filtered_motifs, title="Class → Subclass")),
                        ('struct_heatmap',    lambda: plot_structural_heatmap(filtered_motifs, sequence_length, title="Structural Potential Heatmap")),
                        ('motif_network',     lambda: plot_motif_network(filtered_motifs, title="Motif Co-occurrence Network")),
                        ('cooccurrence_mat',  lambda: plot_motif_cooccurrence_matrix(filtered_motifs, title="Co-occurrence Matrix")),
                        ('chrom_density',     lambda: plot_chromosome_density(filtered_motifs, title="Motif Density by Class")),
                        ('cluster_dist_fig',  lambda: plot_motif_clustering_distance(filtered_motifs, title="Inter-Motif Distance")),
                        ('spacer_loop',       lambda: plot_spacer_loop_variation(filtered_motifs, title="Structural Features Distribution")),
                        ('upset_fig',         lambda: plot_structural_competition_upset(filtered_motifs, title="Structural Competition")),
                    ]

                    # Density comparison uses already-computed density dicts
                    _gd = st.session_state.cached_visualizations[viz_cache_key].get('densities', {}).get('class_genomic')
                    _pd = st.session_state.cached_visualizations[viz_cache_key].get('densities', {}).get('class_positional')
                    if _gd and _pd:
                        _plot_tasks.append(('density_comparison', lambda: plot_density_comparison(_gd, _pd, title="Density Analysis")))

                    if _has_clusters or _has_hybrids:
                        _chm = [m for m in filtered_motifs if m.get('Class') in ['Hybrid', 'Non-B_DNA_Clusters']]
                        _plot_tasks.append(('cluster_track', lambda: plot_linear_motif_track(_chm, sequence_length, title="Hybrid & Cluster Track")))
                    if _has_clusters:
                        _plot_tasks.append(('cluster_size', lambda: plot_cluster_size_distribution(filtered_motifs, title="Cluster Statistics")))

                    for fig_key, fig_fn in _plot_tasks:
                        try:
                            fig_bytes[fig_key] = _fig_to_bytes(fig_fn())
                            total_viz_count += 1
                        except Exception:
                            pass

                    st.session_state.cached_visualizations[viz_cache_key]['figures'] = fig_bytes
                
                viz_total_time = time.time() - viz_start_time
                
                # Memory management: Trigger garbage collection after all visualizations
                trigger_garbage_collection()
                logger.debug(f"Triggered garbage collection after generating {total_viz_count} visualizations")
                
                # Ephemeral success with scientific time format
                status_placeholder.success(f"All visualizations prepared: {total_viz_count} components in {format_time_compact(viz_total_time)}")
                
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
                
                # ============================================================
                # HORIZONTAL METRICS STRIP - Analysis Complete Dashboard
                # ============================================================
                # Dashboard-style horizontal metrics bar for professional look
                # ============================================================
                
                timer_placeholder.markdown(f"""
                <div class='metrics-strip metrics-strip--success'>
                    <div class='metric-card'>
                        <span class='metric-card__icon'>⏱️</span>
                        <span class='metric-card__value'>{format_time_scientific(total_time)}</span>
                        <span class='metric-card__label'>Analysis</span>
                    </div>
                    <div class='metric-card'>
                        <span class='metric-card__icon'>🎨</span>
                        <span class='metric-card__value'>{format_time_scientific(viz_total_time)}</span>
                        <span class='metric-card__label'>Viz Time</span>
                    </div>
                    <div class='metric-card'>
                        <span class='metric-card__icon'>🧬</span>
                        <span class='metric-card__value'>{total_bp_processed:,}</span>
                        <span class='metric-card__label'>Base Pairs</span>
                    </div>
                    <div class='metric-card'>
                        <span class='metric-card__icon'>⚡</span>
                        <span class='metric-card__value'>{overall_speed:,.0f}</span>
                        <span class='metric-card__label'>bp/sec</span>
                    </div>
                    <div class='metric-card'>
                        <span class='metric-card__icon'>🧪</span>
                        <span class='metric-card__value'>{sum(len(r) for r in all_results)}</span>
                        <span class='metric-card__label'>Motifs</span>
                    </div>
                    <div class='metric-card'>
                        <span class='metric-card__icon'>⚠️</span>
                        <span class='metric-card__value'>{len(validation_issues)}</span>
                        <span class='metric-card__label'>Issues</span>
                    </div>
                </div>
                """, unsafe_allow_html=True)
                
                # Show comprehensive completion summary with scientific time format
                # GOLD STANDARD: Final time display after everything is done
                # Display as styled colored box
                validation_status = 'PASSED' if len(validation_issues) == 0 else f'{len(validation_issues)} issues found'
                
                # Build completion HTML with styled green box
                # Using concatenation for better readability while keeping single-line HTML for Streamlit
                box_style = "background: linear-gradient(135deg, #28a745 0%, #20c997 100%); border-radius: 12px; padding: 20px; margin: 15px 0; color: white; box-shadow: 0 4px 15px rgba(40, 167, 69, 0.3);"
                section_style = "background: rgba(255,255,255,0.15); border-radius: 8px; padding: 12px; margin: 10px 0;"
                h3_style = "margin: 0 0 15px 0; color: white; font-size: 1.4rem;"
                h4_style = "margin: 0 0 8px 0; color: white; font-size: 1.1rem;"
                ul_style = "margin: 0; padding-left: 20px; list-style-type: disc;"
                footer_style = "margin-top: 15px; font-weight: 600; font-size: 1.1rem; text-align: center;"
                
                completion_html = (
                    f'<div style="{box_style}">'
                    f'<h3 style="{h3_style}">✅ Analysis Complete! All processing stages finished successfully:</h3>'
                    f'<div style="{section_style}">'
                    f'<h4 style="{h4_style}">Detection &amp; Analysis:</h4>'
                    f'<ul style="{ul_style}">'
                    f'<li>{len(DETECTOR_PROCESSES)} detector processes completed</li>'
                    f'<li>{sum(len(r) for r in all_results)} total motifs detected across {len(st.session_state.seqs)} sequences</li>'
                    f'<li>Analysis completed in {format_time_scientific(total_time)} ({overall_speed:,.0f} bp/s)</li>'
                    f'</ul></div>'
                    f'<div style="{section_style}">'
                    f'<h4 style="{h4_style}">Quality Validation:</h4>'
                    f'<ul style="{ul_style}">'
                    f'<li>Data consistency checks: {validation_status}</li>'
                    f'<li>Non-redundancy validation: Complete</li>'
                    f'<li>Position validation: Complete</li>'
                    f'</ul></div>'
                    f'<div style="{section_style}">'
                    f'<h4 style="{h4_style}">Visualization Generation:</h4>'
                    f'<ul style="{ul_style}">'
                    f'<li>{total_viz_count} visualization components pre-generated</li>'
                    f'<li>Class-level and subclass-level analysis ready</li>'
                    f'<li>Visualizations prepared in {format_time_scientific(viz_total_time)}</li>'
                    f'</ul></div>'
                    f'<div style="{footer_style}">View detailed results in the \'Results\' tab.</div>'
                    f'</div>'
                )
                st.markdown(completion_html, unsafe_allow_html=True)
                
                # ============================================================
                # PERFORMANCE SUMMARY: Display comprehensive statistics
                # ============================================================
                perf_tracker.add_task_time('visualization', viz_total_time)
                perf_tracker.end()
                
                # Calculate motif distribution by class
                motifs_by_class = {}
                for results in all_results:
                    for motif in results:
                        motif_class = motif.get('Class', 'Unknown')
                        motifs_by_class[motif_class] = motifs_by_class.get(motif_class, 0) + 1
                
                # Display performance summary
                st.subheader("⚡ Performance Report")
                
                # Create expandable section with detailed performance stats
                with st.expander("📊 View Detailed Performance Statistics", expanded=False):
                    # Format and display comprehensive performance summary
                    perf_summary_text = format_performance_summary(perf_tracker, format_time_human)
                    st.markdown(perf_summary_text)
                    
                    # Add motif distribution
                    if motifs_by_class:
                        st.markdown("---")
                        motif_dist_text = format_motif_distribution(motifs_by_class)
                        st.markdown(motif_dist_text)
                
                # Display key metrics as cards
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric(
                        "Total Time",
                        format_time(total_time),
                        delta=None,
                        help="Total analysis time in human-readable format"
                    )
                with col2:
                    st.metric(
                        "Throughput",
                        f"{overall_speed:,.0f} bp/s",
                        delta=None,
                        help="Average processing speed"
                    )
                with col3:
                    st.metric(
                        "Total Motifs",
                        f"{sum(len(r) for r in all_results):,}",
                        delta=None,
                        help="Total motifs detected across all sequences"
                    )
                with col4:
                    mem_mb = get_memory_usage_mb()
                    st.metric(
                        "Peak Memory",
                        f"{mem_mb:.0f} MB",
                        delta=None,
                        help="Memory usage at completion"
                    )
                
                # Summary at the end (no popup)
                status_placeholder.empty()  # Clear the status
                st.success(
                    f"🎉 Analysis complete! {sum(len(r) for r in all_results):,} motifs found in {format_time(total_time)}"
                )
                
                st.session_state.analysis_status = "Complete"
                
                # Set analysis_done flag for idempotent run button
                st.session_state.analysis_done = True
                st.session_state.analysis_time = total_time
                
                # ============================================================
                # RESULT PERSISTENCE: Save results to disk (silently)
                # ============================================================
                job_id = st.session_state.get('current_job_id')
                
                # Safety check: Regenerate job ID if somehow missing
                if not job_id:
                    logger.warning("Job ID missing from session state, regenerating")
                    job_id = generate_job_id()
                    st.session_state.current_job_id = job_id
                
                # Prepare metadata
                job_metadata = {
                    'analysis_time': total_time,
                    'speed_bp_per_sec': overall_speed,
                    'detector_count': len(DETECTOR_PROCESSES),
                    'visualization_count': total_viz_count,
                    'validation_issues': len(validation_issues),
                    'analysis_mode': st.session_state.get('analysis_mode_used', 'Motif Level')
                }
                
                # Save results to disk silently (no job ID display to user)
                # Log failures but don't display to user since results are available in session
                save_success = save_job_results(
                    job_id,
                    all_results,
                    st.session_state.seqs,
                    st.session_state.names,
                    job_metadata
                )
                if not save_success:
                    logger.warning(f"Failed to persist results for job {job_id} - results available in session only")
                
            except IndexError:
                # Occurs when accessing empty lists during analysis (e.g., no patterns found, empty results)
                clear_analysis_placeholders(progress_placeholder, status_placeholder,
                                            detailed_progress_placeholder, timer_placeholder)
                
                # Log with automatic traceback for developer debugging
                logger.exception("IndexError during analysis")
                
                # Error message removed per requirements - no longer displaying to user
                # The analysis will fail silently and technical details are available in logs
                
                # Show technical details in expander for advanced users
                show_technical_details()
                
                st.session_state.analysis_status = "Error"
                
            except Exception as e:
                # Catch-all for other unexpected errors
                clear_analysis_placeholders(progress_placeholder, status_placeholder,
                                            detailed_progress_placeholder, timer_placeholder)
                
                # Log with automatic traceback for developer debugging
                logger.exception("Unexpected error during analysis")
                
                # Display error to user
                st.error(f"❌ **Analysis failed:** {str(e)}")
                
                # Show traceback in expander for debugging
                show_technical_details()
                
                st.session_state.analysis_status = "Error"

    # End of Upload & Analyze tab
    st.markdown("---")
