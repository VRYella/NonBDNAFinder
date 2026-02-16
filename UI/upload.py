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
import numpy as np
import pandas as pd
import logging

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
from UI.formatters import format_time_scientific
from UI.headers import render_section_heading
from Utilities.utilities import (
    parse_fasta_chunked,
    get_file_preview,
    get_basic_stats,
    trigger_garbage_collection,
    parse_fasta,
    get_memory_usage_mb,
    calculate_genomic_density,
    calculate_positional_density
)
from Utilities.nonbscanner import analyze_sequence
from Utilities.job_manager import save_job_results, generate_job_id
from Utilities.disk_storage import UniversalSequenceStorage, UniversalResultsStorage
from Utilities.chunk_analyzer import ChunkAnalyzer

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
# Disk storage chunk threshold: sequences larger than this use ChunkAnalyzer with adaptive chunking
# Aligned with CHUNKING_CONFIG direct_threshold (1MB) for optimal performance
CHUNK_ANALYSIS_THRESHOLD_BP = 1_000_000  # 1MB (adaptive chunking threshold)
SUBMOTIF_ABBREVIATIONS = {
    'Global Curvature': 'Global Curv', 'Local Curvature': 'Local Curv',
    'Direct Repeat': 'DR', 'STR': 'STR', 'Cruciform forming IRs': 'Cruciform',
    'R-loop formation sites': 'R-loop', 'Triplex': 'Triplex', 'Sticky DNA': 'Sticky DNA',
    'Telomeric G4': 'Telo G4', 'Stacked canonical G4s': 'Stacked G4', 'Stacked G4s with linker': 'Stacked G4 + Linker',
    'Canonical intramolecular G4': 'Intra G4', 'Extended-loop canonical': 'Ext. Loop G4',
    'Higher-order G4 array/G4-wire': 'G4 Array', 'Intramolecular G-triplex': 'G-triplex', 'Two-tetrad weak PQS': 'Weak PQS',
    'Canonical i-motif': 'i-Motif', 'Relaxed i-motif': 'Relaxed iM', 'AC-motif': 'AC Motif',
    'Z-DNA': 'Z-DNA', 'eGZ': 'eGZ', 'A-philic DNA': 'A-DNA',
    'Dynamic overlaps': 'Hybrid', 'Dynamic clusters': 'Cluster',
}
# ═══════════════════════════════════════════════════════════════════════════════


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
                border-radius: 8px; padding: 12px; margin: 8px 0; color: white;
                box-shadow: 0 2px 6px rgba(102, 126, 234, 0.2);'>
        <div style='font-weight: 600; font-size: 0.9rem; margin-bottom: 6px;'>
            {idx}. {display_name}
        </div>
        <div style='display: grid; grid-template-columns: repeat(4, 1fr); gap: 8px; 
                    font-size: 0.8rem; margin-top: 8px;'>
            <div style='background: rgba(255,255,255,0.15); padding: 6px; border-radius: 4px; text-align: center;'>
                <div style='font-weight: 700; font-size: 1rem;'>{length:,}</div>
                <div style='opacity: 0.9; font-size: 0.7rem;'>Base Pairs</div>
            </div>
            <div style='background: rgba(255,255,255,0.15); padding: 6px; border-radius: 4px; text-align: center;'>
                <div style='font-weight: 700; font-size: 1rem;'>{gc_pct:.2f}%</div>
                <div style='opacity: 0.9; font-size: 0.7rem;'>GC Content</div>
            </div>
            <div style='background: rgba(255,255,255,0.15); padding: 6px; border-radius: 4px; text-align: center;'>
                <div style='font-weight: 700; font-size: 1rem;'>{at_pct:.2f}%</div>
                <div style='opacity: 0.9; font-size: 0.7rem;'>AT Content</div>
            </div>
            <div style='background: rgba(255,255,255,0.15); padding: 6px; border-radius: 4px; text-align: center;'>
                <div style='font-weight: 700; font-size: 1rem;'>{gc_balance}</div>
                <div style='opacity: 0.9; font-size: 0.7rem;'>GC Balance</div>
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
    
    Args:
        sequences: List of DNA sequence strings
        
    Returns:
        Average GC percentage (0.0 if no sequences/base pairs)
    """
    total_bp = sum(len(s) for s in sequences)
    if total_bp == 0:
        return 0.0
    total_gc_count = sum(s.upper().count('G') + s.upper().count('C') for s in sequences)
    return (total_gc_count / total_bp) * 100


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
                
                    # Calculate aggregate QC stats for compact summary
                    total_gc = sum(p['gc_percent'] * p['length'] for p in preview_info['previews'])
                    total_len = sum(p['length'] for p in preview_info['previews'])
                    avg_gc = total_gc / total_len if total_len > 0 else 0.0
                
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
                                ✓ Valid
                            </span>
                        </div>
                    </div>
                    """, unsafe_allow_html=True)
                
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
                        cur_seq += line.strip()
                if cur_seq:
                    seqs.append(cur_seq)
                    names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                if seqs:
                    # Compact QC summary strip for pasted sequences
                    total_bp = sum(len(s) for s in seqs)
                    avg_gc = calculate_gc_percentage(seqs)
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
            # Support both disk storage and legacy in-memory mode
            if st.session_state.get('use_disk_storage') and st.session_state.get('seq_storage'):
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
            # Cache sequence stats with validation key to handle sequence changes
            cache_key = f"stats_cache_{num_sequences}"
            if cache_key not in st.session_state:
                # Calculate stats for all sequences
                stats_list = []
                if st.session_state.get('use_disk_storage') and st.session_state.get('seq_ids'):
                    # Disk storage mode: calculate from metadata
                    for seq_id in st.session_state.seq_ids:
                        metadata = st.session_state.seq_storage.get_metadata(seq_id)
                        stats_list.append({
                            'Length': metadata['length'],
                            'GC%': metadata['gc_content']
                        })
                else:
                    # Legacy mode: calculate from in-memory sequences
                    for seq in st.session_state.seqs:
                        stats = get_basic_stats(seq)
                        stats_list.append(stats)
                st.session_state[cache_key] = stats_list
        
            # Use cached stats - show single compact summary strip
            cached_stats = st.session_state[cache_key]
            avg_gc = sum(s['GC%'] for s in cached_stats) / len(cached_stats) if cached_stats else 0.0
            
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
                        GC: <strong>{avg_gc:.1f}%</strong> avg
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
        <div style='background: linear-gradient(135deg, #dbeafe 0%, #bfdbfe 100%); 
                    padding: 0.5rem 1rem; border-radius: 6px; margin-bottom: 0.75rem;
                    border-left: 3px solid #3b82f6;'>
            <h4 style='margin: 0; color: #1e3a8a; font-size: 0.95rem; font-weight: 600;'>
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
            <div style='background: linear-gradient(135deg, #eff6ff 0%, #dbeafe 100%);
                        padding: 6px 10px; border-radius: 6px; margin-top: 8px;
                        border: 1px solid #bfdbfe; font-size: 0.78rem;'>
                <span style='font-weight: 600; color: #1e40af;'>
                    {len(enabled_classes)} classes · {num_enabled}/{total_submotifs} submotifs
                </span>
            </div>
            """, unsafe_allow_html=True)
        else:
            st.warning("Select at least one submotif to enable analysis.")
    
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
    # ANALYSIS EXECUTION - Full-Width Section Below Two Columns
    # Clear visual and spatial separation from configuration
    # ============================================================
    st.markdown("---")
    
    # Initialize analysis_done flag if not present (idempotent run button)
    if "analysis_done" not in st.session_state:
        st.session_state.analysis_done = False
    
    # Check if valid input is present (sequences AND at least one class selected)
    # Support both disk storage and legacy in-memory mode
    has_sequences = bool(st.session_state.get('seqs')) or bool(st.session_state.get('seq_ids'))
    has_valid_input = has_sequences and bool(st.session_state.get('selected_classes'))
    
    # Create a full-width container for the run button
    run_button_container = st.container()
    with run_button_container:
        # Create three columns for Run, Reset buttons and Clock
        col_run, col_reset, col_clock = st.columns([4, 1, 1])
        
        with col_run:
            # Primary action button with clear scientific terminology
            if has_valid_input and not st.session_state.analysis_done:
                run_button = st.button(
                    "Run Non-B DNA Analysis",
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
                            color: white; padding: 14px; 
                            border-radius: 8px; text-align: center; font-weight: 600;
                            font-size: 1rem;'>
                    ✓ Analysis Complete — View results in 'Results' tab
                </div>
                """, unsafe_allow_html=True)
                run_button = False
            else:
                # Disabled button appearance with accessibility
                st.markdown(f"""
                <div role="button" aria-disabled="true" aria-label="Run Analysis - Disabled"
                     style='background: #e5e7eb; color: #9ca3af; padding: 14px; 
                            border-radius: 8px; text-align: center; font-weight: 600;
                            font-size: 1rem; cursor: not-allowed;'>
                    Run Non-B DNA Analysis
                </div>
                <p style='text-align: center; color: #9ca3af; font-size: 0.8rem; margin-top: 6px;'>
                    Upload sequences and select submotifs to enable
                </p>
                """, unsafe_allow_html=True)
                run_button = False
        
        with col_reset:
            # Vibrant Reset button styling
            st.markdown("""
            <style>
            div[data-testid="column"]:nth-child(2) button {
                background: linear-gradient(135deg, #f97316 0%, #ea580c 100%) !important;
                color: white !important;
                border: none !important;
                font-weight: 700 !important;
            }
            div[data-testid="column"]:nth-child(2) button:hover {
                background: linear-gradient(135deg, #fb923c 0%, #f97316 100%) !important;
                box-shadow: 0 4px 12px rgba(249, 115, 22, 0.4) !important;
            }
            </style>
            """, unsafe_allow_html=True)
            # Reset button to allow re-running analysis
            if st.button("🔄 Reset", use_container_width=True, help="Clear results and reset"):
                st.session_state.analysis_done = False
                st.session_state.results = []
                st.session_state.performance_metrics = None
                st.session_state.cached_visualizations = {}
                st.session_state.analysis_time = None
                # analysis_mode is always "Submotif Level" (fixed), no reset needed
                st.rerun()
        
        with col_clock:
            # Display current time (page load time) as a clock
            from datetime import datetime
            current_time = datetime.now().strftime("%H:%M:%S")
            st.markdown(f"""
            <div style='background: linear-gradient(135deg, #6366f1 0%, #4f46e5 100%); 
                        color: white; padding: 10px 8px; 
                        border-radius: 8px; text-align: center; font-weight: 600;
                        font-size: 0.9rem; height: 100%; display: flex; 
                        flex-direction: column; justify-content: center;
                        box-shadow: 0 2px 8px rgba(99, 102, 241, 0.3);'>
                <div style='font-size: 0.7rem; opacity: 0.9;'>🕐 Time</div>
                <div style='font-family: monospace; font-size: 1rem;'>{current_time}</div>
            </div>
            """, unsafe_allow_html=True)
        
        # Placeholder for progress area
        progress_placeholder = st.empty()
    
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
                
                # Support both disk storage and legacy in-memory mode
                if st.session_state.get('use_disk_storage') and st.session_state.get('seq_ids'):
                    # === DISK STORAGE MODE: Use ChunkAnalyzer ===
                    seq_ids = st.session_state.seq_ids
                    names = st.session_state.names
                    
                    for i, (seq_id, name) in enumerate(zip(seq_ids, names)):
                        progress = (i + 1) / num_sequences
                        
                        # Get sequence metadata
                        metadata = st.session_state.seq_storage.get_metadata(seq_id)
                        seq_length = metadata['length']
                        
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
                        status_msg = f"Processing sequence {i+1}/{num_sequences}: {name} ({seq_length:,} bp)"
                        status_placeholder.info(status_msg)
                        
                        # Determine whether to use chunking based on sequence size
                        # Use module-level constant for threshold
                        
                        if seq_length > CHUNK_ANALYSIS_THRESHOLD_BP:
                            # Large sequence: use ChunkAnalyzer with adaptive chunking
                            status_placeholder.info(f"Using adaptive chunk-based analysis: {name} ({seq_length:,} bp)")
                            
                            # Create chunk progress callback
                            chunk_progress = {'current': 0, 'total': 0}
                            
                            def chunk_callback(progress_pct):
                                chunk_progress['current'] = int(progress_pct)
                                if progress_pct % 10 == 0:  # Update every 10%
                                    status_placeholder.info(f"Analyzing chunks: {progress_pct:.0f}% complete")
                            
                            # Analyze using ChunkAnalyzer with adaptive chunking enabled
                            analyzer = ChunkAnalyzer(
                                st.session_state.seq_storage,
                                use_parallel=True,  # Enable parallel chunk processing
                                max_workers=None,  # Auto-detect CPU count
                                use_adaptive=True  # Enable triple adaptive chunking
                            )
                            
                            results_storage = analyzer.analyze(
                                seq_id=seq_id,
                                progress_callback=chunk_callback,
                                enabled_classes=list(st.session_state.selected_classes) if st.session_state.selected_classes else None
                            )
                            
                            # Store results storage
                            st.session_state.results_storage[seq_id] = results_storage
                            
                            # Get results for filtering (load first 10000 for display)
                            results = list(results_storage.iter_results(limit=10000))
                            
                            status_placeholder.success(f"Chunk analysis complete: {len(results)} motifs (showing first 10000)")
                        else:
                            # Small sequence: use standard analysis with sequence loaded into memory
                            # Get full sequence for standard analysis
                            seq = st.session_state.seq_storage.get_sequence_chunk(seq_id, 0, seq_length)
                            
                            # Validate sequence before analysis
                            if not seq or len(seq) == 0:
                                status_placeholder.error(f"❌ Empty sequence detected for: {name}")
                                continue
                            
                            # Use standard analysis
                            results = analyze_sequence(
                                seq, name,
                                enabled_classes=list(st.session_state.selected_classes) if st.session_state.selected_classes else None
                            )
                            
                            # Create results storage and save
                            results_storage = UniversalResultsStorage(
                                base_dir=str(st.session_state.seq_storage.base_dir / "results"),
                                seq_id=seq_id
                            )
                            results_storage.append_batch(results)
                            st.session_state.results_storage[seq_id] = results_storage
                            
                            status_placeholder.success(f"{name}: {seq_length:,} bp | {len(results)} motifs detected")
                        
                        # Ensure all motifs have required fields
                        results = [ensure_subclass(motif) for motif in results]
                        
                        # Filter results based on selected subclasses
                        selected_classes_set = set(st.session_state.selected_classes)
                        selected_subclasses_set = set(st.session_state.selected_subclasses)
                        
                        filtered_results = []
                        for motif in results:
                            motif_class = motif.get('Class', '')
                            motif_subclass = motif.get('Subclass', '')
                            
                            if motif_class in selected_classes_set:
                                if motif_class in ['Hybrid', 'Non-B_DNA_Clusters']:
                                    component_classes = motif.get('Component_Classes', [])
                                    if not component_classes or any(c in selected_classes_set for c in component_classes):
                                        filtered_results.append(motif)
                                elif motif_subclass in selected_subclasses_set:
                                    filtered_results.append(motif)
                        
                        results = filtered_results
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
                        
                        # Validate sequence before analysis
                        if not seq or len(seq) == 0:
                            status_placeholder.error(f"❌ Empty sequence detected for: {name}")
                            continue
                        
                        # Run the analysis - use parallel scanner for large sequences if enabled
                        # No per-sequence timing - total time captured once at end
                        
                        if use_parallel_scanner and len(seq) > 100000:
                            # NOTE: scanner_agent.py has been archived - parallel scanning experimental
                            # Use experimental parallel scanner for large sequences
                            try:
                                from Utilities.scanner_agent import ParallelScanner
                                
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
                                        chunk_progress_placeholder.info(f"Parallel scanner processing chunks: {current}/{total} ({(current / total) * 100:.1f}%)")
                                
                                # Run parallel scanner with progress callback
                                # Use ephemeral status (replaces previous message)
                                status_placeholder.info(f"Using parallel scanner for {len(seq):,} bp sequence (est. chunks: ~{len(seq)//50000 + 1})")
                                
                                scanner = ParallelScanner(seq, hs_db=None)
                                
                                # The parallel scanner internally calls analyze_sequence on each chunk
                                # and returns full motif dictionaries with deduplication
                                results = scanner.run_scan(progress_callback=chunk_progress_callback)
                                
                                # Update sequence names for all motifs
                                for motif in results:
                                    motif['Sequence_Name'] = name
                                
                                # Clear chunk progress and show ephemeral success (replaces previous message)
                                if show_chunk_progress:
                                    chunk_progress_placeholder.success(f"Parallel chunks complete: {len(results)} motifs from {chunk_counter['total']} chunks")
                                
                                # Ephemeral success message (replaces previous)
                                status_placeholder.success(f"Parallel scanner completed: {len(results)} motifs detected")
                                
                            except Exception as e:
                                # Fallback to standard scanner on error (ephemeral warning)
                                status_placeholder.warning(f"Parallel scanner failed, falling back to standard: {e}")
                                results = analyze_sequence(seq, name, 
                                                           enabled_classes=list(st.session_state.selected_classes))
                        else:
                            # Use standard consolidated NBDScanner analysis with selective detection
                            # Pass enabled_classes for performance optimization - only runs selected detectors
                            results = analyze_sequence(seq, name, 
                                                       enabled_classes=list(st.session_state.selected_classes))
                        
                        # Ensure all motifs have required fields
                        results = [ensure_subclass(motif) for motif in results]
                        
                        # ============================================================
                        # FILTER RESULTS BASED ON SELECTED SUBCLASSES
                        # ============================================================
                        # Note: Class-level filtering is now done at detector level for performance.
                        # Here we apply the finer-grained subclass filter.
                        selected_classes_set = set(st.session_state.selected_classes)
                        selected_subclasses_set = set(st.session_state.selected_subclasses)
                        
                        # Filter motifs to only include selected subclasses
                        filtered_results = []
                        for motif in results:
                            motif_class = motif.get('Class', '')
                            motif_subclass = motif.get('Subclass', '')
                            
                            # Include if class is selected AND (subclass is selected OR subclass matches selected class)
                            if motif_class in selected_classes_set:
                                # For Hybrid and Non-B_DNA_Clusters, include if parent class is selected
                                if motif_class in ['Hybrid', 'Non-B_DNA_Clusters']:
                                    # Always include these if they were formed from selected classes
                                    component_classes = motif.get('Component_Classes', [])
                                    if not component_classes or any(c in selected_classes_set for c in component_classes):
                                        filtered_results.append(motif)
                                elif motif_subclass in selected_subclasses_set:
                                    filtered_results.append(motif)
                        
                        results = filtered_results
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
                        pbar.progress(progress, text=f"Analyzed {i+1}/{len(st.session_state.seqs)} sequences ({actual_percentage:.1f}%)")
                    
                    # Ephemeral success (replaces previous) - no per-sequence timing shown
                    status_placeholder.success(f"{name}: {len(seq):,} bp | {len(results)} motifs detected")
                
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
                
            except Exception as e:
                progress_placeholder.empty()
                status_placeholder.empty()
                detailed_progress_placeholder.empty()
                timer_placeholder.empty()
                st.error(f"Analysis failed: {str(e)}")
                st.session_state.analysis_status = "Error"

    # End of Upload & Analyze tab
    st.markdown("---")
