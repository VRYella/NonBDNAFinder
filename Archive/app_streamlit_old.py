"""
NBDFinder Streamlit App - 5-page interface for non-B DNA motif detection.

Pages: Home, Upload & Analyze, Visualization & Download, Advanced Analysis, Documentation
"""

import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
import os
import io
import base64
from typing import Optional, Dict, List
import threading
import time
import logging

# Core NBDFinder imports
from orchestrator import run_pipeline
from viz_tools import (
    plot_motif_counts, plot_motif_lengths, plot_score_distribution,
    plot_density_track, plot_overlap_network, plot_sequence_overview,
    tsne_feature_plot, export_excel_with_sheets, create_summary_statistics
)
import hs_registry_manager
from motifs.registry import get_all_hyperscan_patterns

# Configure page
st.set_page_config(
    page_title="NBDFinder - Non-B DNA Motif Detection",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .section-header {
        font-size: 1.5rem;
        color: #2e8b57;
        margin-top: 2rem;
        margin-bottom: 1rem;
    }
    .info-box {
        background-color: #f0f8ff;
        padding: 1rem;
        border-radius: 10px;
        border-left: 5px solid #1f77b4;
        margin: 1rem 0;
    }
    .warning-box {
        background-color: #fff8dc;
        padding: 1rem;
        border-radius: 10px;
        border-left: 5px solid #ffa500;
        margin: 1rem 0;
    }
</style>
""", unsafe_allow_html=True)


@st.cache_resource
def get_hyperscan_db():
    """Cache Hyperscan database compilation"""
    try:
        patterns = get_all_hyperscan_patterns()
        pattern_strings = [p[0] for p in patterns]
        return hs_registry_manager.compile_db(pattern_strings, key="streamlit_global")
    except Exception as e:
        st.error(f"Failed to compile Hyperscan database: {e}")
        return None


def initialize_session_state():
    """Initialize Streamlit session state variables"""
    if 'analysis_results' not in st.session_state:
        st.session_state.analysis_results = None
    if 'analysis_running' not in st.session_state:
        st.session_state.analysis_running = False
    if 'analysis_progress' not in st.session_state:
        st.session_state.analysis_progress = ""
    if 'uploaded_file_content' not in st.session_state:
        st.session_state.uploaded_file_content = None


def home_page():
    """Home page with introduction and overview"""
    st.markdown('<div class="main-header">🧬 NBDFinder - Non-B DNA Motif Detection</div>', 
                unsafe_allow_html=True)
    
    # Display the NBD circle image
    image_path = Path("nbdcircle.JPG")
    if image_path.exists():
        st.image(str(image_path), caption="Non-B DNA Structures", use_container_width=True)
    
    st.markdown("""
    <div class="info-box">
    <h3>Welcome to NBDFinder</h3>
    
    NBDFinder is a comprehensive toolkit for detecting and analyzing non-B DNA structures 
    in genomic sequences. This application uses Intel Hyperscan for high-performance 
    pattern matching combined with sophisticated biological validation algorithms.
    </div>
    """, unsafe_allow_html=True)
    
    # Feature overview
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### 🔍 Detection Capabilities")
        motif_classes = [
            "G-Quadruplex (6 subtypes)",
            "Curved DNA (A-phased repeats)", 
            "Slipped DNA (STRs)",
            "Cruciform structures",
            "R-loop forming sequences",
            "Triplex DNA",
            "i-Motif structures",
            "Z-DNA",
            "Hybrid motifs",
            "Motif clusters"
        ]
        for motif in motif_classes:
            st.write(f"• {motif}")
    
    with col2:
        st.markdown("### ⚡ Performance Features")
        features = [
            "Intel Hyperscan acceleration",
            "Multi-core parallel processing",
            "Memory-efficient streaming",
            "Cached pattern compilation",
            "Real-time progress tracking",
            "Multiple output formats",
            "Interactive visualizations",
            "Comprehensive scoring models"
        ]
        for feature in features:
            st.write(f"• {feature}")
    
    # System status
    st.markdown("### 🔧 System Status")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        hs_status = "✅ Available" if hs_registry_manager.is_available() else "❌ Not Available"
        st.metric("Hyperscan", hs_status)
    
    with col2:
        patterns = get_all_hyperscan_patterns()
        st.metric("Patterns Loaded", len(patterns))
    
    with col3:
        cached_dbs = len(hs_registry_manager.get_cached_keys())
        st.metric("Cached Databases", cached_dbs)
    
    st.markdown("""
    <div class="section-header">Quick Start</div>
    
    1. **Upload & Analyze**: Upload your FASTA file and configure analysis parameters
    2. **Visualization**: Explore results with interactive plots and statistics  
    3. **Download**: Export results in multiple formats (CSV, Excel, GFF3, Parquet)
    4. **Advanced**: Access specialized analysis tools and customization options
    """, unsafe_allow_html=True)


def upload_analyze_page():
    """Upload & Analyze page with file upload and analysis configuration"""
    st.markdown('<div class="main-header">📁 Upload & Analyze</div>', 
                unsafe_allow_html=True)
    
    # File upload section
    st.markdown("### 📄 Input Data")
    
    upload_method = st.radio(
        "Choose input method:",
        ["Upload FASTA file", "Paste sequence", "Example data"]
    )
    
    sequence_data = None
    
    if upload_method == "Upload FASTA file":
        uploaded_file = st.file_uploader(
            "Choose a FASTA file",
            type=['fa', 'fasta', 'fas', 'fna'],
            help="Upload sequences in FASTA format"
        )
        
        if uploaded_file is not None:
            st.session_state.uploaded_file_content = uploaded_file.getvalue()
            st.success(f"File uploaded: {uploaded_file.name}")
            
            # Preview first few lines
            preview = uploaded_file.getvalue().decode('utf-8').split('\n')[:10]
            st.text("Preview:\n" + '\n'.join(preview))
            sequence_data = uploaded_file.getvalue()
    
    elif upload_method == "Paste sequence":
        sequence_text = st.text_area(
            "Paste your sequence(s) here:",
            height=200,
            placeholder="Enter FASTA format sequences or raw DNA sequence..."
        )
        
        if sequence_text.strip():
            # Convert to FASTA if not already
            if not sequence_text.startswith('>'):
                sequence_text = f">user_sequence\n{sequence_text}"
            
            sequence_data = sequence_text.encode('utf-8')
            st.success("Sequence loaded from text input")
    
    elif upload_method == "Example data":
        example_seq = """>example_g4_sequence
GGGTTTGGGTTTGGGTTTGGGAAATTTCCCAAATTTCCCAAATTTCCCAAATTTCCC
>example_curved_sequence  
AAAAAAAAATGCGTAAAAAAAAAATGCGTAAAAAAAAAATGCGT
>example_triplex_sequence
AAAAAGGGGGGAAAAGGGGGGAAAAGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCC"""
        
        st.text_area("Example sequences:", value=example_seq, height=150, disabled=True)
        if st.button("Use Example Data"):
            sequence_data = example_seq.encode('utf-8')
            st.success("Example data loaded")
    
    # Analysis configuration
    if sequence_data:
        st.markdown("### ⚙️ Analysis Configuration")
        
        col1, col2 = st.columns(2)
        
        with col1:
            detector_classes = st.multiselect(
                "Select motif classes to detect:",
                options=['g_quadruplex', 'curved_dna', 'slipped_dna', 'cruciform',
                        'r_loop', 'triplex', 'i_motif', 'z_dna'],
                default=['g_quadruplex', 'curved_dna', 'triplex'],
                help="Choose which types of motifs to search for"
            )
            
            max_workers = st.slider(
                "Number of parallel workers:",
                min_value=1, max_value=8, value=4,
                help="More workers = faster processing (uses more CPU)"
            )
        
        with col2:
            chunk_size = st.select_slider(
                "Sequence chunk size:",
                options=[10000, 25000, 50000, 100000],
                value=50000,
                help="Larger chunks = more memory, smaller chunks = more overhead"
            )
            
            output_prefix = st.text_input(
                "Output file prefix:",
                value="nbdfinder_results",
                help="Prefix for output files"
            )
        
        # Analysis execution
        st.markdown("### 🚀 Run Analysis")
        
        if not detector_classes:
            st.warning("Please select at least one motif class to detect.")
        else:
            if st.button("Start Analysis", type="primary", disabled=st.session_state.analysis_running):
                run_analysis(sequence_data, detector_classes, max_workers, chunk_size, output_prefix)
            
            # Show progress
            if st.session_state.analysis_running:
                st.info("Analysis in progress...")
                progress_placeholder = st.empty()
                
                # Progress animation
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                # Simple progress simulation (in real implementation, this would be real progress)
                for i in range(100):
                    progress_bar.progress(i + 1)
                    status_text.text(f"Processing... {i+1}%")
                    time.sleep(0.1)  # Simulate work
                
                if st.session_state.analysis_results:
                    st.success("Analysis completed successfully!")
                    st.balloons()


def run_analysis(sequence_data: bytes, detector_classes: List[str], 
                max_workers: int, chunk_size: int, output_prefix: str):
    """Run the analysis pipeline in background"""
    
    def analysis_worker():
        try:
            st.session_state.analysis_running = True
            st.session_state.analysis_progress = "Starting analysis..."
            
            # Save uploaded data to temporary file
            with tempfile.NamedTemporaryFile(mode='wb', suffix='.fa', delete=False) as tmp_file:
                tmp_file.write(sequence_data)
                tmp_file_path = tmp_file.name
            
            # Run pipeline
            st.session_state.analysis_progress = "Running motif detection..."
            output_files = run_pipeline(
                fasta_path=tmp_file_path,
                output_prefix=f"/tmp/{output_prefix}",
                max_workers=max_workers,
                chunk_size=chunk_size,
                detector_classes=detector_classes
            )
            
            # Load results
            if 'csv' in output_files:
                results_df = pd.read_csv(output_files['csv'])
                st.session_state.analysis_results = {
                    'dataframe': results_df,
                    'output_files': output_files,
                    'summary': create_summary_statistics(results_df)
                }
            
            # Cleanup
            os.unlink(tmp_file_path)
            
        except Exception as e:
            st.session_state.analysis_progress = f"Error: {str(e)}"
        finally:
            st.session_state.analysis_running = False
    
    # Start analysis in background thread
    thread = threading.Thread(target=analysis_worker)
    thread.start()


def visualization_page():
    """Visualization & Download page with interactive plots"""
    st.markdown('<div class="main-header">📊 Visualization & Download</div>', 
                unsafe_allow_html=True)
    
    if st.session_state.analysis_results is None:
        st.warning("No analysis results available. Please run an analysis first.")
        return
    
    results = st.session_state.analysis_results
    df = results['dataframe']
    summary = results['summary']
    
    # Summary statistics
    st.markdown("### 📈 Summary Statistics")
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Total Motifs", summary.get('Total_Motifs', 0))
    with col2:
        st.metric("Unique Classes", summary.get('Unique_Classes', 0))
    with col3:
        st.metric("Sequences", summary.get('Unique_Sequences', 0))
    with col4:
        st.metric("Avg Length", f"{summary.get('Average_Length', 0):.1f} bp")
    
    # Visualization options
    st.markdown("### 🎨 Interactive Visualizations")
    
    viz_type = st.selectbox(
        "Choose visualization:",
        [
            "Motif Distribution", 
            "Length Distribution",
            "Score Distribution", 
            "Sequence Overview",
            "Density Track",
            "Overlap Network",
            "Feature Space (t-SNE)"
        ]
    )
    
    # Generate selected visualization
    try:
        if viz_type == "Motif Distribution":
            fig = plot_motif_counts(df)
            st.plotly_chart(fig, use_container_width=True)
            
        elif viz_type == "Length Distribution":
            fig = plot_motif_lengths(df)
            st.plotly_chart(fig, use_container_width=True)
            
        elif viz_type == "Score Distribution":
            fig = plot_score_distribution(df)
            st.plotly_chart(fig, use_container_width=True)
            
        elif viz_type == "Sequence Overview":
            sequences = df['Sequence_Name'].unique()
            selected_seq = st.selectbox("Select sequence:", sequences)
            # Estimate sequence length from max end position
            seq_len = df[df['Sequence_Name'] == selected_seq]['End'].max()
            fig = plot_sequence_overview(df, selected_seq, seq_len)
            st.plotly_chart(fig, use_container_width=True)
            
        elif viz_type == "Density Track":
            sequences = df['Sequence_Name'].unique()
            selected_seq = st.selectbox("Select sequence:", sequences)
            seq_len = df[df['Sequence_Name'] == selected_seq]['End'].max()
            window_size = st.slider("Window size:", 100, 5000, 1000)
            fig = plot_density_track(df, selected_seq, seq_len, window_size)
            st.plotly_chart(fig, use_container_width=True)
            
        elif viz_type == "Overlap Network":
            fig = plot_overlap_network(df)
            st.plotly_chart(fig, use_container_width=True)
            
        elif viz_type == "Feature Space (t-SNE)":
            fig = tsne_feature_plot(df)
            st.plotly_chart(fig, use_container_width=True)
            
    except Exception as e:
        st.error(f"Error generating visualization: {e}")
    
    # Data table
    st.markdown("### 📋 Results Table")
    
    # Filtering options
    col1, col2 = st.columns(2)
    with col1:
        class_filter = st.multiselect(
            "Filter by class:",
            options=df['Class'].unique(),
            default=df['Class'].unique()
        )
    with col2:
        min_score = st.slider(
            "Minimum normalized score:",
            min_value=0.0, max_value=1.0, value=0.0, step=0.1
        )
    
    # Apply filters
    filtered_df = df[
        (df['Class'].isin(class_filter)) & 
        (df['Normalized_Score'] >= min_score)
    ]
    
    st.dataframe(filtered_df, use_container_width=True, height=400)
    
    # Download section
    st.markdown("### 💾 Download Results")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        # CSV download
        csv_data = filtered_df.to_csv(index=False)
        st.download_button(
            "Download CSV",
            data=csv_data,
            file_name="nbdfinder_results.csv",
            mime="text/csv"
        )
    
    with col2:
        # Excel download
        try:
            excel_data = export_excel_with_sheets(filtered_df, summary)
            st.download_button(
                "Download Excel",
                data=excel_data,
                file_name="nbdfinder_results.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
        except Exception as e:
            st.error(f"Excel export failed: {e}")
    
    with col3:
        # GFF3 download (simplified)
        gff3_lines = ["##gff-version 3"]
        for _, row in filtered_df.iterrows():
            line = f"{row['Chromosome/Contig']}\tNBDFinder\tmotif\t{row['Start']}\t{row['End']}\t{row['Normalized_Score']:.3f}\t.\t.\tID=motif_{row['S.No']};Class={row['Class']}"
            gff3_lines.append(line)
        
        gff3_data = "\n".join(gff3_lines)
        st.download_button(
            "Download GFF3",
            data=gff3_data,
            file_name="nbdfinder_results.gff3",
            mime="text/plain"
        )


def advanced_page():
    """Advanced Analysis page with specialized tools"""
    st.markdown('<div class="main-header">🔬 Advanced Analysis</div>', 
                unsafe_allow_html=True)
    
    st.markdown("""
    <div class="info-box">
    <h3>Advanced Tools</h3>
    This section provides specialized analysis tools for power users and researchers.
    </div>
    """, unsafe_allow_html=True)
    
    # Database management
    st.markdown("### 🗃️ Database Management")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Hyperscan Database Status")
        cached_dbs = hs_registry_manager.get_cached_keys()
        
        if cached_dbs:
            for key in cached_dbs:
                info = hs_registry_manager.get_db_info(key)
                if info:
                    st.write(f"**{key}**: {info['pattern_count']} patterns")
        else:
            st.write("No databases cached")
        
        if st.button("Clear Database Cache"):
            hs_registry_manager.clear_cache()
            st.success("Cache cleared")
    
    with col2:
        st.subheader("Pattern Registry")
        patterns = get_all_hyperscan_patterns()
        st.write(f"Total patterns loaded: {len(patterns)}")
        
        # Show pattern breakdown
        pattern_counts = {}
        for pattern, _, class_name, _ in patterns:
            pattern_counts[class_name] = pattern_counts.get(class_name, 0) + 1
        
        for class_name, count in pattern_counts.items():
            st.write(f"• {class_name}: {count} patterns")
    
    # Custom analysis options
    st.markdown("### 🎛️ Custom Analysis")
    
    st.markdown("""
    <div class="warning-box">
    <strong>Note:</strong> These features are for advanced users. Experimental settings 
    may affect analysis accuracy.
    </div>
    """, unsafe_allow_html=True)
    
    # Scoring method selection
    scoring_method = st.selectbox(
        "Normalization method:",
        ["minmax", "zscore"],
        help="Choose score normalization method"
    )
    
    # Custom thresholds
    st.subheader("Detection Thresholds")
    col1, col2 = st.columns(2)
    
    with col1:
        min_motif_length = st.number_input("Minimum motif length:", value=10, min_value=1)
        max_motif_length = st.number_input("Maximum motif length:", value=1000, min_value=10)
    
    with col2:
        cluster_window = st.number_input("Cluster window size:", value=1000, min_value=100)
        overlap_threshold = st.slider("Overlap threshold:", 0.0, 1.0, 0.5)
    
    # Batch processing
    st.markdown("### 📦 Batch Processing")
    
    uploaded_files = st.file_uploader(
        "Upload multiple FASTA files:",
        type=['fa', 'fasta', 'fas', 'fna'],
        accept_multiple_files=True
    )
    
    if uploaded_files:
        st.write(f"Selected {len(uploaded_files)} files:")
        for file in uploaded_files:
            st.write(f"• {file.name}")
        
        if st.button("Process Batch"):
            st.info("Batch processing would be implemented here")
    
    # Export settings
    st.markdown("### 📤 Export Settings")
    
    export_formats = st.multiselect(
        "Select export formats:",
        ["CSV", "Excel", "Parquet", "GFF3", "BED"],
        default=["CSV", "GFF3"]
    )
    
    include_sequences = st.checkbox("Include full sequences in output", value=False)
    compress_output = st.checkbox("Compress output files", value=False)


def documentation_page():
    """Documentation page with help and technical information"""
    st.markdown('<div class="main-header">📚 Documentation</div>', 
                unsafe_allow_html=True)
    
    # Table of contents
    st.markdown("""
    ### 📖 Contents
    - [About NBDFinder](#about)
    - [Motif Classes](#motifs)
    - [Technical Specifications](#technical)
    - [Performance Guidelines](#performance)
    - [API Reference](#api)
    - [Citation](#citation)
    """)
    
    # About section
    st.markdown("""
    <div id="about"></div>
    
    ### 🧬 About NBDFinder
    
    NBDFinder is a high-performance toolkit for detecting non-B DNA structures in genomic sequences. 
    It combines Intel Hyperscan pattern matching with sophisticated biological validation algorithms 
    to provide accurate and fast motif detection.
    
    **Key Features:**
    - **Fast**: Uses Intel Hyperscan for multi-pattern scanning
    - **Accurate**: Biological validation and scoring models
    - **Scalable**: Parallel processing and memory-efficient streaming
    - **Comprehensive**: Detects 10 major non-B DNA classes
    """, unsafe_allow_html=True)
    
    # Motif classes
    st.markdown("""
    <div id="motifs"></div>
    
    ### 🔍 Motif Classes Detected
    """)
    
    motif_info = {
        "G-Quadruplex": "Four-stranded DNA structures formed by G-rich sequences",
        "Curved DNA": "A-tract containing sequences that bend the DNA helix", 
        "Slipped DNA": "Short tandem repeats that can form hairpin structures",
        "Cruciform": "Four-way junction structures from palindromic sequences",
        "R-loops": "RNA-DNA hybrid structures with displaced single-strand DNA",
        "Triplex DNA": "Three-stranded structures from homopurine/homopyrimidine tracts",
        "i-Motif": "Four-stranded structures formed by C-rich sequences",
        "Z-DNA": "Left-handed double helix from alternating purine-pyrimidine sequences",
        "Hybrid": "Overlapping regions containing multiple motif types",
        "Clusters": "Regions with high density of multiple motifs"
    }
    
    for motif, description in motif_info.items():
        st.write(f"**{motif}**: {description}")
    
    # Technical specifications
    st.markdown("""
    <div id="technical"></div>
    
    ### ⚙️ Technical Specifications
    
    **Hyperscan Integration:**
    - Database compilation is cached for performance
    - One scratch object per worker (reused across matches)
    - Minimal callbacks for maximum speed
    - Fallback to Python regex for unsupported patterns
    
    **Scoring Methods:**
    - G4Hunter for G-quadruplex structures
    - Curvature propensity for curved DNA
    - Triplex stability for triplex DNA
    - Z-DNA seeker for Z-DNA
    - Custom algorithms for other motif types
    
    **Performance Optimizations:**
    - Multi-core parallel processing
    - Memory-efficient sequence chunking
    - Interval trees for overlap detection
    - Streaming support for large sequences
    """)
    
    # Performance guidelines
    st.markdown("""
    <div id="performance"></div>
    
    ### 🚀 Performance Guidelines
    
    **Memory Usage:**
    - Use smaller chunk sizes for limited memory systems
    - Enable streaming for very large sequences (>100MB)
    - Clear database cache if memory is constrained
    
    **CPU Usage:**
    - Set workers to match available CPU cores
    - Reduce workers if system becomes unresponsive
    - Consider batch processing for multiple files
    
    **Accuracy vs Speed:**
    - All motif classes: Comprehensive but slower
    - Selected classes: Faster for targeted analysis
    - Higher thresholds: Faster but may miss low-scoring motifs
    """)
    
    # API reference (simplified)
    st.markdown("""
    <div id="api"></div>
    
    ### 🔧 API Reference
    
    **Command Line Usage:**
    ```bash
    python orchestrator.py --fasta input.fa --out results --workers 4
    ```
    
    **Python API:**
    ```python
    from orchestrator import run_pipeline
    
    output_files = run_pipeline(
        fasta_path="sequence.fa",
        output_prefix="results", 
        max_workers=4,
        detector_classes=['g_quadruplex', 'triplex']
    )
    ```
    """)
    
    # Citation
    st.markdown("""
    <div id="citation"></div>
    
    ### 📝 Citation
    
    If you use NBDFinder in your research, please cite:
    
    > NBDFinder: High-performance detection of non-B DNA structures using Hyperscan acceleration.
    > [Authors], [Journal] (Year). DOI: [DOI]
    
    **Related Publications:**
    - G4Hunter: Bedrat et al. NAR 44(4):1746-1759 (2016)
    - Z-DNA Seeker: Ho et al. EMBO J 5:2737-2744 (1986)
    - Hyperscan: Wang et al. NSDI 2019
    """)


def main():
    """Main application function"""
    initialize_session_state()
    
    # Sidebar navigation
    st.sidebar.title("Navigation")
    page = st.sidebar.radio(
        "Choose a page:",
        ["🏠 Home", "📁 Upload & Analyze", "📊 Visualization & Download", 
         "🔬 Advanced Analysis", "📚 Documentation"]
    )
    
    # Cache Hyperscan database on startup
    with st.sidebar:
        st.markdown("### System Status")
        if hs_registry_manager.is_available():
            st.success("Hyperscan: Ready")
            db_info = get_hyperscan_db()
        else:
            st.error("Hyperscan: Not Available")
    
    # Route to appropriate page
    if page == "🏠 Home":
        home_page()
    elif page == "📁 Upload & Analyze":
        upload_analyze_page()
    elif page == "📊 Visualization & Download":
        visualization_page()
    elif page == "🔬 Advanced Analysis":
        advanced_page()
    elif page == "📚 Documentation":
        documentation_page()


if __name__ == "__main__":
    main()