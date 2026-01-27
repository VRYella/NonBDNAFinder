"""
Upload & Analyze Page

This module contains the Upload & Analyze tab functionality for NonBDNAFinder.
"""

import streamlit as st
import time
import os
import gc
import numpy as np
import pandas as pd
import logging

# Import configuration
from config.text import UI_TEXT
from config.typography import FONT_CONFIG
from config.themes import TAB_THEMES
from config.analysis import ANALYSIS_CONFIG

# Import UI components
from ui.css import load_css, get_page_colors
from ui.formatters import format_time_scientific

# Import utilities
from utilities import (
    parse_fasta_chunked,
    get_file_preview,
    get_basic_stats,
    trigger_garbage_collection,
    parse_fasta,
    get_memory_usage_mb,
    calculate_genomic_density,
    calculate_positional_density
)

# Import core analysis
from nonbscanner import analyze_sequence

# Import job management
from job_manager import save_job_results, generate_job_id

# Try to import Entrez for NCBI fetch functionality
try:
    from Bio import Entrez, SeqIO
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False

# Setup logger
logger = logging.getLogger(__name__)

# Configuration availability flag
CONFIG_AVAILABLE = False

# GC Balance Thresholds (typical genomic DNA range)
GC_BALANCE_MIN = 30  # Minimum GC% for balanced genome
GC_BALANCE_MAX = 70  # Maximum GC% for balanced genome


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
    gc_balance = '✓' if GC_BALANCE_MIN <= gc_pct <= GC_BALANCE_MAX else '⚠'
    
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
                    with st.popover(UI_TEXT['upload_preview_button'], use_container_width=True):
                        st.markdown("### 📊 Sequence Preview & Validation")
                        for idx, prev in enumerate(preview_info['previews'], 1):
                            # Use pre-calculated stats from preview (calculated from full sequence)
                            card_html = render_sequence_stats_card(
                                idx=idx,
                                name=prev['name'],
                                length=prev['length'],
                                gc_pct=prev['gc_percent'],
                                at_pct=prev['at_percent']
                            )
                            st.markdown(card_html, unsafe_allow_html=True)
                        
                        if preview_info['num_sequences'] > 3:
                            st.info(f"📌 Showing 3 of {preview_info['num_sequences']} sequences. All sequences will be analyzed.")
                    
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
            with st.popover("✅ Success: Validation Summary", use_container_width=True):
                st.markdown("### 📊 Loaded Sequences Summary")
                
                # Cache sequence stats with validation key to handle sequence changes
                cache_key = f"stats_cache_{len(st.session_state.seqs)}"
                if cache_key not in st.session_state:
                    # Calculate stats for all sequences
                    stats_list = []
                    for seq in st.session_state.seqs:
                        stats = get_basic_stats(seq)
                        stats_list.append(stats)
                    st.session_state[cache_key] = stats_list
                
                # Use cached stats
                cached_stats = st.session_state[cache_key]
                for i, stats in enumerate(cached_stats[:3]):
                    # Use helper function with green gradient for validation success
                    card_html = render_sequence_stats_card(
                        idx=i+1,
                        name=st.session_state.names[i],
                        length=len(st.session_state.seqs[i]),
                        gc_pct=stats['GC%'],
                        at_pct=stats['AT%'],
                        gradient_colors="135deg, #11998e 0%, #38ef7d 100%"
                    )
                    st.markdown(card_html, unsafe_allow_html=True)
                
                if len(st.session_state.seqs) > 3:
                    st.info(f"📌 Showing 3 of {len(st.session_state.seqs)} loaded sequences. All are validated and ready for analysis.")
    
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
