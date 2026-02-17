"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Streamlit Progress Integration - Real-time Analysis Feedback                 │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.2            │
│ Provides live progress bars, detector status, and performance metrics        │
└──────────────────────────────────────────────────────────────────────────────┘
"""

import time
import logging
from typing import List, Dict, Any, Optional, Callable
from collections import defaultdict

logger = logging.getLogger(__name__)

# Import Streamlit (graceful degradation if not available)
try:
    import streamlit as st
    STREAMLIT_AVAILABLE = True
except ImportError:
    STREAMLIT_AVAILABLE = False
    logger.warning("Streamlit not available. Progress tracking disabled.")


class StreamlitProgressTracker:
    """
    Real-time progress tracking for Streamlit UI during NonBDNA analysis.
    
    Provides:
    - Overall progress bar
    - Per-detector status indicators
    - Performance metrics (throughput, time estimates)
    - Live result streaming
    
    Features:
    - Minimal UI overhead (<50ms per update)
    - Automatic cleanup on completion
    - Thread-safe updates
    - Graceful degradation if Streamlit unavailable
    """
    
    def __init__(self, sequence_name: str, sequence_length: int, num_detectors: int = 9):
        """
        Initialize progress tracker.
        
        Args:
            sequence_name: Name of sequence being analyzed
            sequence_length: Length of sequence in base pairs
            num_detectors: Number of detectors to track (default: 9)
        """
        self.sequence_name = sequence_name
        self.sequence_length = sequence_length
        self.num_detectors = num_detectors
        self.start_time = time.time()
        
        # Tracking state
        self.completed_detectors = 0
        self.total_motifs = 0
        self.detector_status: Dict[str, Dict[str, Any]] = {}
        
        # Streamlit UI elements (created on first render)
        self.progress_bar = None
        self.status_text = None
        self.metrics_container = None
        self.detector_container = None
        
        self.ui_created = False
        
        if not STREAMLIT_AVAILABLE:
            logger.warning("Streamlit not available. Progress tracking will be minimal.")
    
    def create_ui(self) -> None:
        """Create Streamlit UI elements for progress tracking."""
        if not STREAMLIT_AVAILABLE or self.ui_created:
            return
        
        # Header with gradient background (no emojis, HTML-escaped for safety)
        import html
        sequence_name_safe = html.escape(self.sequence_name)
        st.markdown(f"""
        <div style='background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                    padding: 1.5rem; border-radius: 12px; margin-bottom: 1rem;
                    box-shadow: 0 4px 6px rgba(0,0,0,0.1);'>
            <h3 style='color: white; margin: 0; font-weight: 700;'>ANALYZING: {sequence_name_safe}</h3>
            <p style='color: rgba(255,255,255,0.9); margin: 0.5rem 0 0 0; font-size: 0.95rem;'>
                Sequence Length: {self.sequence_length:,} bp
            </p>
        </div>
        """, unsafe_allow_html=True)
        
        # Vibrant progress bar
        self.progress_bar = st.progress(0)
        self.status_text = st.empty()
        
        # Colorful metrics row with boxes
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            self.metric_progress = st.empty()
        with col2:
            self.metric_time = st.empty()
        with col3:
            self.metric_throughput = st.empty()
        with col4:
            self.metric_motifs = st.empty()
        
        # Detector status header with vibrant styling
        st.markdown("""
        <div style='background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); 
                    padding: 0.5rem 1rem; border-radius: 8px; margin: 1rem 0 0.5rem 0;'>
            <h4 style='color: white; margin: 0; font-weight: 600;'>DETECTOR STATUS</h4>
        </div>
        """, unsafe_allow_html=True)
        self.detector_container = st.container()
        
        self.ui_created = True
    
    def update(self, detector_name: str, status: str = "running", 
               motifs: int = 0, elapsed: float = 0.0) -> None:
        """
        Update progress for a specific detector.
        
        Args:
            detector_name: Name of detector (e.g., "g_quadruplex")
            status: Current status ("running", "complete", "error")
            motifs: Number of motifs found
            elapsed: Time elapsed for this detector
        """
        # Update internal state
        self.detector_status[detector_name] = {
            'status': status,
            'motifs': motifs,
            'elapsed': elapsed
        }
        
        if status == "complete":
            self.completed_detectors += 1
            self.total_motifs += motifs
        
        # Update UI
        self._render_ui()
    
    def _render_ui(self) -> None:
        """Render/update the Streamlit UI elements."""
        if not STREAMLIT_AVAILABLE or not self.ui_created:
            return
        
        # Calculate progress
        progress = self.completed_detectors / self.num_detectors
        elapsed_time = time.time() - self.start_time
        
        # Update progress bar
        if self.progress_bar:
            self.progress_bar.progress(progress)
        
        # Update status text with vibrant colored box (no emojis)
        if self.status_text:
            if progress < 1.0:
                status_color = "#FF9100"  # Orange for processing
                status_bg = "linear-gradient(135deg, #FFE57F 0%, #FFD54F 100%)"
                status_msg = f"PROCESSING: {self.completed_detectors}/{self.num_detectors} detectors complete"
            else:
                status_color = "#00E676"  # Green for complete
                status_bg = "linear-gradient(135deg, #B9F6CA 0%, #69F0AE 100%)"
                status_msg = "ANALYSIS COMPLETE"
            
            self.status_text.markdown(f"""
            <div style='background: {status_bg}; 
                        padding: 0.75rem; border-radius: 8px; 
                        border-left: 4px solid {status_color};
                        margin: 0.5rem 0;'>
                <strong style='color: {status_color};'>{status_msg}</strong>
            </div>
            """, unsafe_allow_html=True)
        
        # Update metrics with vibrant colored boxes
        if hasattr(self, 'metric_progress'):
            prog_pct = progress * 100
            self.metric_progress.markdown(f"""
            <div style='background: linear-gradient(135deg, #E040FB 0%, #D500F9 100%); 
                        padding: 1rem; border-radius: 8px; text-align: center;
                        box-shadow: 0 2px 4px rgba(0,0,0,0.1);'>
                <div style='color: white; font-size: 1.8rem; font-weight: 700;'>{prog_pct:.0f}%</div>
                <div style='color: rgba(255,255,255,0.9); font-size: 0.8rem;'>PROGRESS</div>
            </div>
            """, unsafe_allow_html=True)
        
        if hasattr(self, 'metric_time'):
            self.metric_time.markdown(f"""
            <div style='background: linear-gradient(135deg, #0091FF 0%, #0043A8 100%); 
                        padding: 1rem; border-radius: 8px; text-align: center;
                        box-shadow: 0 2px 4px rgba(0,0,0,0.1);'>
                <div style='color: white; font-size: 1.8rem; font-weight: 700;'>{elapsed_time:.1f}s</div>
                <div style='color: rgba(255,255,255,0.9); font-size: 0.8rem;'>ELAPSED</div>
            </div>
            """, unsafe_allow_html=True)
        
        if hasattr(self, 'metric_throughput'):
            throughput = self.sequence_length / elapsed_time if elapsed_time > 0 else 0
            self.metric_throughput.markdown(f"""
            <div style='background: linear-gradient(135deg, #00E676 0%, #00C853 100%); 
                        padding: 1rem; border-radius: 8px; text-align: center;
                        box-shadow: 0 2px 4px rgba(0,0,0,0.1);'>
                <div style='color: white; font-size: 1.8rem; font-weight: 700;'>{throughput:,.0f}</div>
                <div style='color: rgba(255,255,255,0.9); font-size: 0.8rem;'>BP/S</div>
            </div>
            """, unsafe_allow_html=True)
        
        if hasattr(self, 'metric_motifs'):
            self.metric_motifs.markdown(f"""
            <div style='background: linear-gradient(135deg, #FF6D00 0%, #FF9100 100%); 
                        padding: 1rem; border-radius: 8px; text-align: center;
                        box-shadow: 0 2px 4px rgba(0,0,0,0.1);'>
                <div style='color: white; font-size: 1.8rem; font-weight: 700;'>{self.total_motifs:,}</div>
                <div style='color: rgba(255,255,255,0.9); font-size: 0.8rem;'>MOTIFS</div>
            </div>
            """, unsafe_allow_html=True)
        
        # Update detector status
        if self.detector_container:
            with self.detector_container:
                self._render_detector_status()
    
    def _render_detector_status(self) -> None:
        """Render the detector status table with vibrant colored boxes (no emojis)."""
        if not self.detector_status:
            return
        
        # Create vibrant status boxes for each detector
        status_html = "<div style='display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 0.5rem;'>"
        
        for detector, info in sorted(self.detector_status.items()):
            status = info['status']
            
            # Choose color scheme based on status (no emojis/icons)
            if status == "complete":
                bg_color = "linear-gradient(135deg, #B9F6CA 0%, #69F0AE 100%)"
                border_color = "#00E676"
                status_text = "COMPLETE"
                text_color = "#00612E"
            elif status == "running":
                bg_color = "linear-gradient(135deg, #FFE57F 0%, #FFAB00 100%)"
                border_color = "#FF9100"
                status_text = "RUNNING"
                text_color = "#BF360C"
            elif status == "error":
                bg_color = "linear-gradient(135deg, #FFCDD2 0%, #FF5252 100%)"
                border_color = "#FF1744"
                status_text = "ERROR"
                text_color = "#B71C1C"
            else:
                bg_color = "linear-gradient(135deg, #E2E8F0 0%, #CBD5E1 100%)"
                border_color = "#64748B"
                status_text = "PENDING"
                text_color = "#334155"
            
            elapsed_str = f"{info['elapsed']:.2f}s" if info['elapsed'] > 0 else "-"
            motifs_str = f"{info['motifs']:,}" if status == "complete" else "-"
            detector_display = detector.replace('_', ' ').title()
            
            status_html += f"""
            <div style='background: {bg_color}; 
                        padding: 0.75rem; border-radius: 8px; 
                        border-left: 4px solid {border_color};
                        box-shadow: 0 2px 4px rgba(0,0,0,0.1);'>
                <div style='font-weight: 700; color: {text_color}; font-size: 0.9rem;'>{detector_display}</div>
                <div style='display: flex; justify-content: space-between; margin-top: 0.25rem;'>
                    <span style='color: {text_color}; font-size: 0.8rem; font-weight: 600;'>{status_text}</span>
                    <span style='color: {text_color}; font-size: 0.8rem;'>{elapsed_str}</span>
                </div>
                <div style='color: {text_color}; font-size: 0.75rem; margin-top: 0.1rem;'>Motifs: {motifs_str}</div>
            </div>
            """
        
        status_html += "</div>"
        st.markdown(status_html, unsafe_allow_html=True)
    
    def finalize(self) -> Dict[str, Any]:
        """
        Finalize progress tracking and return summary.
        
        Returns:
            Dictionary with analysis summary statistics
        """
        total_time = time.time() - self.start_time
        
        summary = {
            'sequence_name': self.sequence_name,
            'sequence_length': self.sequence_length,
            'total_motifs': self.total_motifs,
            'total_time': total_time,
            'throughput': self.sequence_length / total_time if total_time > 0 else 0,
            'detectors_run': self.completed_detectors,
            'detector_details': self.detector_status
        }
        
        # Show completion message with vibrant styling (no emojis)
        if STREAMLIT_AVAILABLE and self.ui_created:
            st.markdown(f"""
            <div style='background: linear-gradient(135deg, #B9F6CA 0%, #69F0AE 100%); 
                        padding: 1rem; border-radius: 8px; 
                        border-left: 4px solid #00E676;
                        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                        text-align: center;'>
                <strong style='color: #00612E; font-size: 1.1rem;'>
                    ANALYSIS COMPLETE | Found {self.total_motifs:,} motifs in {total_time:.2f}s
                </strong>
            </div>
            """, unsafe_allow_html=True)
        
        return summary


def analyze_with_streamlit_progress(sequence: str, 
                                    sequence_name: str,
                                    enabled_classes: Optional[List[str]] = None,
                                    **kwargs) -> List[Dict[str, Any]]:
    """
    Analyze sequence with real-time Streamlit progress tracking.
    
    This is a wrapper around the standard analyze_sequence() function that
    adds real-time progress updates to the Streamlit UI.
    
    Args:
        sequence: DNA sequence to analyze
        sequence_name: Name/identifier for the sequence
        enabled_classes: List of motif classes to detect (None = all)
        **kwargs: Additional arguments passed to analyze_sequence()
    
    Returns:
        List of motif dictionaries (same format as analyze_sequence())
    
    Example:
        >>> motifs = analyze_with_streamlit_progress(
        ...     sequence="ATCGATCG...",
        ...     sequence_name="chr1",
        ...     enabled_classes=["G-Quadruplex", "Z-DNA"]
        ... )
    """
    # Import here to avoid circular dependency
    try:
        from Utilities.nonbscanner import analyze_sequence, get_detector_display_names
        
        # Try optimized version first
        try:
            from Utilities.nonbscanner_optimized import NonBScannerOptimized
            scanner = NonBScannerOptimized(enable_all_detectors=True)
            use_optimized = True
            logger.info("Using optimized NonBScanner")
        except ImportError:
            scanner = None
            use_optimized = False
            logger.info("Using standard NonBScanner")
    except ImportError as e:
        logger.error(f"Failed to import NonBScanner: {e}")
        raise
    
    # Create progress tracker
    detector_names = get_detector_display_names()
    num_detectors = len(enabled_classes) if enabled_classes else len(detector_names)
    tracker = StreamlitProgressTracker(sequence_name, len(sequence), num_detectors)
    
    # Create UI if Streamlit available
    if STREAMLIT_AVAILABLE:
        tracker.create_ui()
    
    # Create progress callback
    def progress_callback(detector_name: str, completed: int, total: int, 
                         elapsed: float, motif_count: int):
        tracker.update(detector_name, status="complete", motifs=motif_count, elapsed=elapsed)
    
    # Run analysis
    if use_optimized and scanner:
        motifs = scanner.analyze_sequence(
            sequence=sequence,
            sequence_name=sequence_name,
            progress_callback=progress_callback,
            enabled_classes=enabled_classes,
            **kwargs
        )
    else:
        motifs = analyze_sequence(
            sequence=sequence,
            sequence_name=sequence_name,
            progress_callback=progress_callback,
            enabled_classes=enabled_classes,
            **kwargs
        )
    
    # Finalize and return
    summary = tracker.finalize()
    logger.info(f"Analysis complete: {summary['total_motifs']} motifs in {summary['total_time']:.2f}s")
    
    return motifs


def display_results_streaming(motifs: List[Dict[str, Any]], 
                              batch_size: int = 100) -> None:
    """
    Display analysis results with streaming updates for large result sets.
    
    For large numbers of motifs, this displays results in batches to avoid
    freezing the Streamlit UI.
    
    Args:
        motifs: List of motif dictionaries from analysis
        batch_size: Number of motifs to display per batch
    """
    if not STREAMLIT_AVAILABLE:
        logger.warning("Streamlit not available. Cannot display results.")
        return
    
    if not motifs:
        st.info("ℹ️ No motifs found in this sequence.")
        return
    
    # Summary statistics
    st.markdown("### 📊 Analysis Results")
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total Motifs", f"{len(motifs):,}")
    with col2:
        unique_classes = len(set(m.get('Class', 'Unknown') for m in motifs))
        st.metric("Unique Classes", unique_classes)
    with col3:
        if motifs:
            avg_score = sum(m.get('Score', 0) for m in motifs) / len(motifs)
            st.metric("Avg Score", f"{avg_score:.2f}")
    
    # Class distribution
    class_counts = defaultdict(int)
    for motif in motifs:
        class_counts[motif.get('Class', 'Unknown')] += 1
    
    st.markdown("#### 📈 Class Distribution")
    st.bar_chart(class_counts)
    
    # Results table (paginated for large datasets)
    st.markdown("#### 🧬 Detected Motifs")
    
    if len(motifs) > batch_size:
        st.info(f"ℹ️ Showing results in batches of {batch_size} for performance")
        
        num_batches = (len(motifs) + batch_size - 1) // batch_size
        
        for batch_idx in range(num_batches):
            start_idx = batch_idx * batch_size
            end_idx = min(start_idx + batch_size, len(motifs))
            batch = motifs[start_idx:end_idx]
            
            with st.expander(f"Batch {batch_idx + 1}/{num_batches} (motifs {start_idx+1}-{end_idx})"):
                _display_motif_batch(batch)
    else:
        _display_motif_batch(motifs)


def _display_motif_batch(motifs: List[Dict[str, Any]]) -> None:
    """Helper function to display a batch of motifs as a table."""
    if not STREAMLIT_AVAILABLE:
        return
    
    # Convert to DataFrame for display
    import pandas as pd
    
    df_data = []
    for m in motifs:
        df_data.append({
            'Class': m.get('Class', 'Unknown'),
            'Subclass': m.get('Subclass', 'Unknown'),
            'Start': m.get('Start', 0),
            'End': m.get('End', 0),
            'Length': m.get('Length', 0),
            'Score': f"{m.get('Score', 0):.2f}",
            'Sequence': m.get('Sequence', '')[:50] + ('...' if len(m.get('Sequence', '')) > 50 else '')
        })
    
    if df_data:
        df = pd.DataFrame(df_data)
        st.dataframe(df, use_container_width=True)


# ═══════════════════════════════════════════════════════════════════════════════
# COMPACT PROGRESS INDICATORS FOR NOTEBOOK/CLI
# ═══════════════════════════════════════════════════════════════════════════════

class CompactProgressTracker:
    """
    Minimal progress tracker for CLI/notebook environments.
    
    Provides simple console output without Streamlit dependency.
    """
    
    def __init__(self, sequence_name: str, sequence_length: int):
        self.sequence_name = sequence_name
        self.sequence_length = sequence_length
        self.start_time = time.time()
        self.completed = 0
        self.total = 9  # Default number of detectors
    
    def update(self, detector_name: str, motifs: int):
        """Update progress for a detector."""
        self.completed += 1
        elapsed = time.time() - self.start_time
        progress_pct = (self.completed / self.total) * 100
        
        print(f"[{progress_pct:5.1f}%] {detector_name}: {motifs} motifs ({elapsed:.1f}s)", flush=True)
    
    def finalize(self):
        """Print final summary."""
        total_time = time.time() - self.start_time
        throughput = self.sequence_length / total_time
        
        print(f"\n✓ Analysis complete in {total_time:.2f}s ({throughput:,.0f} bp/s)")


def create_compact_progress_callback(tracker: CompactProgressTracker) -> Callable:
    """
    Create a progress callback for compact tracking.
    
    Args:
        tracker: CompactProgressTracker instance
    
    Returns:
        Callback function compatible with analyze_sequence()
    """
    def callback(detector_name: str, completed: int, total: int, 
                 elapsed: float, motif_count: int):
        tracker.update(detector_name, motif_count)
    
    return callback
