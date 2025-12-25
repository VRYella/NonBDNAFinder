"""
╔══════════════════════════════════════════════════════════════════════════════╗
║               SCIENTIFIC PROGRESS PANEL FOR STREAMLIT UI                      ║
║          Nobel-Laureate-Level Interface: Biology-First Storytelling          ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: scientific_progress.py
PURPOSE: Live scientific progress visualization for Non-B DNA analysis
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.2

DESCRIPTION:
    Provides real-time progress tracking with scientific communication:
    - Motif class detection status (✓ ▶ ⏳)
    - Count and throughput metrics per detector
    - Genome coverage progress bar
    - Biology-first presentation (not developer logs)

Design Principles:
    ✓ Clarity over cleverness
    ✓ Figures before numbers
    ✓ Biology-first storytelling
    ✓ NO console spam
    ✓ NO per-sequence popups

Progress Panel Format:
┌────────────────────────────────────────────────────────────────────────────┐
│ Motif Class      │ Status │ Count │ Speed (bp/s) │ Time                    │
├────────────────────────────────────────────────────────────────────────────┤
│ G-quadruplex     │   ✓    │  142  │   8.3 Mb/s   │ 0.42s                  │
│ Cruciform DNA    │   ▶    │   61  │   7.9 Mb/s   │ Running...             │
│ Slipped DNA      │   ⏳    │   —   │      —       │ Pending                │
└────────────────────────────────────────────────────────────────────────────┘

Progress reflects genome coverage, not file load.
"""

import streamlit as st
import time
from typing import Dict, Any, Optional, Callable
import pandas as pd

# Detector display names (matches nonbscanner.py)
DETECTOR_DISPLAY_NAMES = {
    'curved_dna': 'Curved DNA',
    'slipped_dna': 'Slipped DNA',
    'cruciform': 'Cruciform',
    'r_loop': 'R-Loop',
    'triplex': 'Triplex',
    'g_quadruplex': 'G-Quadruplex',
    'i_motif': 'i-Motif',
    'z_dna': 'Z-DNA',
    'a_philic': 'A-philic DNA'
}

# Status icons for clarity
STATUS_ICONS = {
    'pending': '⏳',      # Hourglass - waiting
    'running': '▶',      # Play - in progress
    'complete': '✓',     # Checkmark - done
    'error': '✗'         # Cross - failed
}

# Color codes for status (Streamlit-compatible)
STATUS_COLORS = {
    'pending': '#FFA500',    # Orange
    'running': '#4169E1',    # Royal Blue
    'complete': '#32CD32',   # Lime Green
    'error': '#DC143C'       # Crimson
}


class StreamlitProgressPanel:
    """
    Scientific progress panel optimized for Streamlit real-time updates.
    
    Provides live status tracking with:
    - Per-detector status visualization
    - Motif count accumulation
    - Throughput metrics (bp/s)
    - Overall progress percentage
    """
    
    def __init__(self, sequence_length: int, sequence_name: str = "sequence"):
        """
        Initialize progress panel.
        
        Args:
            sequence_length: Length of sequence being analyzed (bp)
            sequence_name: Name/identifier of sequence
        """
        self.sequence_length = sequence_length
        self.sequence_name = sequence_name
        self.start_time = time.time()
        
        # Initialize detector tracking
        self.detector_status = {}
        for name in DETECTOR_DISPLAY_NAMES:
            self.detector_status[name] = {
                'status': 'pending',
                'elapsed': 0.0,
                'motifs': 0,
                'speed': 0.0
            }
        
        # Overall statistics
        self.total_motifs = 0
        self.completed_detectors = 0
        self.total_detectors = len(DETECTOR_DISPLAY_NAMES)
    
    def update_detector(self, detector_name: str, completed: bool = False,
                       motifs: int = 0, elapsed: float = 0.0,
                       error: Optional[str] = None) -> None:
        """
        Update status for a specific detector.
        
        Args:
            detector_name: Internal detector name (e.g., 'g_quadruplex')
            completed: Whether detector has finished
            motifs: Number of motifs found
            elapsed: Time elapsed in seconds
            error: Error message if detector failed
        """
        if detector_name not in self.detector_status:
            return
        
        if error:
            self.detector_status[detector_name]['status'] = 'error'
            self.detector_status[detector_name]['elapsed'] = elapsed  # Track time even for errors
            self.completed_detectors += 1  # Count errors as completed to reach 100%
        elif completed:
            self.detector_status[detector_name]['status'] = 'complete'
            self.detector_status[detector_name]['elapsed'] = elapsed
            self.detector_status[detector_name]['motifs'] = motifs
            if elapsed > 0:
                self.detector_status[detector_name]['speed'] = self.sequence_length / elapsed
            self.total_motifs += motifs
            self.completed_detectors += 1
        else:
            self.detector_status[detector_name]['status'] = 'running'
    
    def get_progress_percentage(self) -> float:
        """Get overall progress as percentage (0-100)."""
        return (self.completed_detectors / self.total_detectors) * 100
    
    def get_elapsed_time(self) -> float:
        """Get total elapsed time in seconds."""
        return time.time() - self.start_time
    
    def render_progress_table(self, container) -> None:
        """
        Render progress table in a Streamlit container.
        
        Uses native Streamlit dataframe for clean, professional display.
        
        Args:
            container: Streamlit container to render in
        """
        # Build progress data
        rows = []
        for detector_name, display_name in DETECTOR_DISPLAY_NAMES.items():
            status_info = self.detector_status[detector_name]
            status = status_info['status']
            
            # Format values
            status_icon = STATUS_ICONS.get(status, '○')
            count_str = str(status_info['motifs']) if status == 'complete' else '—'
            
            if status_info['speed'] > 0:
                # Format speed in appropriate units
                speed = status_info['speed']
                if speed >= 1_000_000:
                    speed_str = f"{speed / 1_000_000:.1f} Mb/s"
                elif speed >= 1_000:
                    speed_str = f"{speed / 1_000:.1f} kb/s"
                else:
                    speed_str = f"{speed:.0f} bp/s"
            else:
                speed_str = '—'
            
            # Show elapsed time for complete and error states (for debugging)
            if status in ('complete', 'error') and status_info['elapsed'] > 0:
                time_str = f"{status_info['elapsed']:.2f}s"
            else:
                time_str = '—'
            
            rows.append({
                'Motif Class': display_name,
                'Status': f"{status_icon} {status.capitalize()}",
                'Count': count_str,
                'Speed': speed_str,
                'Time': time_str
            })
        
        # Create DataFrame
        df = pd.DataFrame(rows)
        
        # Render in container
        with container:
            st.dataframe(
                df,
                use_container_width=True,
                hide_index=True,
                column_config={
                    'Motif Class': st.column_config.TextColumn('Motif Class', width='medium'),
                    'Status': st.column_config.TextColumn('Status', width='small'),
                    'Count': st.column_config.TextColumn('Count', width='small'),
                    'Speed': st.column_config.TextColumn('Speed', width='medium'),
                    'Time': st.column_config.TextColumn('Time', width='small')
                }
            )
    
    def render_progress_metrics(self, container) -> None:
        """
        Render summary metrics above the progress table.
        
        Args:
            container: Streamlit container to render in
        """
        elapsed = self.get_elapsed_time()
        progress_pct = self.get_progress_percentage()
        
        # Calculate overall throughput
        if elapsed > 0:
            throughput = self.sequence_length / elapsed
            if throughput >= 1_000_000:
                throughput_str = f"{throughput / 1_000_000:.2f} Mb/s"
            elif throughput >= 1_000:
                throughput_str = f"{throughput / 1_000:.1f} kb/s"
            else:
                throughput_str = f"{throughput:.0f} bp/s"
        else:
            throughput_str = "—"
        
        with container:
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric(
                    label="Progress",
                    value=f"{progress_pct:.0f}%",
                    delta=f"{self.completed_detectors}/{self.total_detectors} detectors"
                )
            
            with col2:
                st.metric(
                    label="Total Motifs",
                    value=f"{self.total_motifs:,}"
                )
            
            with col3:
                st.metric(
                    label="Throughput",
                    value=throughput_str
                )
            
            with col4:
                st.metric(
                    label="Elapsed Time",
                    value=f"{elapsed:.1f}s"
                )
    
    def render_full_panel(self, metrics_container, table_container) -> None:
        """
        Render complete progress panel (metrics + table).
        
        Args:
            metrics_container: Container for summary metrics
            table_container: Container for detector status table
        """
        self.render_progress_metrics(metrics_container)
        self.render_progress_table(table_container)


def create_streamlit_progress_callback(panel: StreamlitProgressPanel,
                                       metrics_container,
                                       table_container,
                                       update_interval: float = 0.5) -> Callable:
    """
    Create a progress callback that updates Streamlit UI in real-time.
    
    NOTE: The current implementation of analyze_sequence() only calls the callback
    AFTER each detector completes (not before it starts). This means we cannot show
    "running" status in real-time. To enable "running" status, analyze_sequence()
    would need to be modified to call the callback twice per detector:
    1. Before starting: callback(detector_name, idx, total, 0.0, 0, starting=True)
    2. After finishing: callback(detector_name, idx+1, total, elapsed, count)
    
    For now, this callback marks detectors as complete when the callback is invoked.
    
    Args:
        panel: StreamlitProgressPanel instance
        metrics_container: Streamlit container for metrics
        table_container: Streamlit container for progress table
        update_interval: Minimum time (seconds) between UI updates
        
    Returns:
        Callback function compatible with nonbscanner.analyze_sequence()
    """
    last_update = [0.0]  # Mutable container to track last update time
    
    def callback(detector_name: str, completed: int, total: int, 
                elapsed: float, motif_count: int):
        """Progress callback that updates Streamlit UI."""
        # Update panel state (mark as complete since callback is only called after completion)
        panel.update_detector(
            detector_name,
            completed=True,
            motifs=motif_count,
            elapsed=elapsed
        )
        
        # Update UI (rate-limited to avoid excessive redraws)
        current_time = time.time()
        if current_time - last_update[0] >= update_interval or completed == total:
            panel.render_full_panel(metrics_container, table_container)
            last_update[0] = current_time
    
    return callback


def create_compact_progress_callback(progress_bar_container,
                                     status_text_container,
                                     detector_order: list = None) -> tuple:
    """
    Create a compact progress callback for minimal UI updates.
    
    Uses just a progress bar and status text (no full table).
    Suitable for sidebar or narrow panels.
    
    Args:
        progress_bar_container: Streamlit container for progress bar
        status_text_container: Streamlit container for status text
        detector_order: Optional list of detector names (for ordered display)
        
    Returns:
        Tuple of (ProgressTracker, callback_function)
    """
    if detector_order is None:
        detector_order = list(DETECTOR_DISPLAY_NAMES.keys())
    
    tracker = {
        'completed': 0,
        'total': len(detector_order),
        'current_detector': None,
        'motif_counts': {}
    }
    
    def callback(detector_name: str, completed: int, total: int,
                elapsed: float, motif_count: int):
        """Compact progress callback."""
        tracker['completed'] = completed
        tracker['total'] = total
        tracker['current_detector'] = DETECTOR_DISPLAY_NAMES.get(detector_name, detector_name)
        tracker['motif_counts'][detector_name] = motif_count
        
        # Update progress bar
        progress = completed / total if total > 0 else 0
        with progress_bar_container:
            st.progress(progress, text=f"Detector {completed}/{total}")
        
        # Update status text
        with status_text_container:
            total_motifs = sum(tracker['motif_counts'].values())
            st.caption(f"**{tracker['current_detector']}**: {motif_count} motifs | **Total**: {total_motifs}")
    
    return tracker, callback


# =============================================================================
# EXAMPLE USAGE (for testing)
# =============================================================================

if __name__ == "__main__":
    st.title("Scientific Progress Panel Demo")
    
    # Create containers
    metrics_container = st.container()
    table_container = st.container()
    
    # Initialize panel
    panel = StreamlitProgressPanel(sequence_length=100000, sequence_name="test_sequence")
    
    if st.button("Run Demo Analysis"):
        # Simulate detector progress
        detectors = list(DETECTOR_DISPLAY_NAMES.keys())
        
        for i, detector in enumerate(detectors):
            # Simulate detector running
            panel.update_detector(detector, completed=False)
            panel.render_full_panel(metrics_container, table_container)
            
            time.sleep(0.3)  # Simulate work
            
            # Detector complete
            import random
            motif_count = random.randint(10, 200)
            elapsed = random.uniform(0.1, 0.5)
            
            panel.update_detector(detector, completed=True, motifs=motif_count, elapsed=elapsed)
            panel.render_full_panel(metrics_container, table_container)
        
        st.success("Analysis Complete!")
        st.balloons()
