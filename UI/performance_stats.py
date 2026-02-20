"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Performance Statistics Module - Real-time and summary performance reporting  │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.2            │
│ Provides performance tracking and formatted summary reports                  │
└──────────────────────────────────────────────────────────────────────────────┘
"""

import time
from typing import Dict, List, Any, Optional
from collections import defaultdict


class PerformanceTracker:
    """Track performance metrics during NonBDNA analysis."""
    
    def __init__(self):
        """Initialize performance tracker."""
        self.start_time = None
        self.end_time = None
        self.detector_times: Dict[str, float] = {}
        self.detector_motifs: Dict[str, int] = {}
        self.sequence_stats: List[Dict[str, Any]] = []
        self.total_bp_processed = 0
        self.total_motifs_found = 0
        self.chunk_times: List[float] = []
        self.task_times: Dict[str, float] = {}
        
    def start(self):
        """Start tracking."""
        self.start_time = time.time()
        
    def end(self):
        """End tracking."""
        self.end_time = time.time()
        
    def add_detector_result(self, detector_name: str, elapsed: float, motif_count: int):
        """Record detector performance."""
        self.detector_times[detector_name] = elapsed
        self.detector_motifs[detector_name] = motif_count
        self.total_motifs_found += motif_count
        
    def add_sequence_result(self, seq_name: str, seq_length: int, elapsed: float, motif_count: int):
        """Record sequence processing results."""
        self.sequence_stats.append({
            'name': seq_name,
            'length': seq_length,
            'elapsed': elapsed,
            'motif_count': motif_count,
            'throughput': seq_length / elapsed if elapsed > 0 else 0
        })
        self.total_bp_processed += seq_length
        
    def add_chunk_time(self, elapsed: float):
        """Record chunk processing time."""
        self.chunk_times.append(elapsed)

    def add_task_time(self, task_name: str, elapsed: float):
        """Record elapsed time for a named task (e.g., 'visualization', 'export')."""
        self.task_times[task_name] = elapsed
        
    def get_total_time(self) -> float:
        """Get total elapsed time."""
        if self.start_time and self.end_time:
            return self.end_time - self.start_time
        elif self.start_time:
            return time.time() - self.start_time
        return 0.0
        
    def get_average_throughput(self) -> float:
        """Calculate average throughput in bp/s."""
        total_time = self.get_total_time()
        if total_time > 0 and self.total_bp_processed > 0:
            return self.total_bp_processed / total_time
        return 0.0
        
    def get_summary(self) -> Dict[str, Any]:
        """Get performance summary."""
        total_time = self.get_total_time()
        
        summary = {
            'total_time': total_time,
            'total_bp_processed': self.total_bp_processed,
            'total_motifs_found': self.total_motifs_found,
            'average_throughput': self.get_average_throughput(),
            'detector_times': self.detector_times.copy(),
            'detector_motifs': self.detector_motifs.copy(),
            'sequence_count': len(self.sequence_stats),
            'chunk_count': len(self.chunk_times),
            'avg_chunk_time': sum(self.chunk_times) / len(self.chunk_times) if self.chunk_times else 0,
            'task_times': self.task_times.copy()
        }
        
        return summary


def format_performance_summary(tracker: PerformanceTracker, format_time_func) -> str:
    """
    Format a comprehensive performance summary for display.
    
    Args:
        tracker: PerformanceTracker instance with collected metrics
        format_time_func: Function to format time values (e.g., format_time_human)
        
    Returns:
        Formatted summary string with performance statistics
    """
    summary = tracker.get_summary()
    
    lines = []
    lines.append("# 📊 Performance Summary")
    lines.append("")
    
    # Overall statistics
    lines.append("## Overall Statistics")
    lines.append(f"- **Total Time**: {format_time_func(summary['total_time'])}")
    lines.append(f"- **Total Base Pairs Processed**: {summary['total_bp_processed']:,} bp")
    lines.append(f"- **Total Motifs Found**: {summary['total_motifs_found']:,}")
    lines.append(f"- **Average Throughput**: {summary['average_throughput']:,.0f} bp/s")
    
    if summary['sequence_count'] > 0:
        lines.append(f"- **Sequences Analyzed**: {summary['sequence_count']}")
    
    if summary['chunk_count'] > 0:
        lines.append(f"- **Chunks Processed**: {summary['chunk_count']}")
        lines.append(f"- **Average Chunk Time**: {format_time_func(summary['avg_chunk_time'])}")
    
    lines.append("")
    
    # Per-task statistics (visualization, export, etc.)
    if summary.get('task_times'):
        lines.append("## Task Timing Breakdown")
        lines.append("")
        for task_name, elapsed in summary['task_times'].items():
            task_display = task_name.replace('_', ' ').title()
            percentage = (elapsed / summary['total_time'] * 100) if summary['total_time'] > 0 else 0
            lines.append(f"- **{task_display}**: {format_time_func(elapsed)} ({percentage:.1f}% of total)")
        lines.append("")

    # Per-detector statistics
    if summary['detector_times']:
        lines.append("## Detector Performance")
        lines.append("")
        
        # Sort detectors by time (slowest first)
        sorted_detectors = sorted(summary['detector_times'].items(), 
                                 key=lambda x: x[1], reverse=True)
        
        for detector, elapsed in sorted_detectors:
            motif_count = summary['detector_motifs'].get(detector, 0)
            percentage = (elapsed / summary['total_time'] * 100) if summary['total_time'] > 0 else 0
            
            detector_display = detector.replace('_', ' ').title()
            lines.append(f"- **{detector_display}**")
            lines.append(f"  - Time: {format_time_func(elapsed)} ({percentage:.1f}% of total)")
            lines.append(f"  - Motifs: {motif_count:,}")
            
            if elapsed > 0 and summary['total_bp_processed'] > 0:
                detector_throughput = summary['total_bp_processed'] / elapsed
                lines.append(f"  - Throughput: {detector_throughput:,.0f} bp/s")
    
    return "\n".join(lines)


def format_motif_distribution(motifs_by_class: Dict[str, int]) -> str:
    """
    Format motif distribution by class.
    
    Args:
        motifs_by_class: Dictionary mapping class names to motif counts
        
    Returns:
        Formatted distribution string
    """
    if not motifs_by_class:
        return ""
    
    lines = []
    lines.append("## Motif Distribution by Class")
    lines.append("")
    
    total = sum(motifs_by_class.values())
    
    # Sort by count (most common first)
    sorted_classes = sorted(motifs_by_class.items(), key=lambda x: x[1], reverse=True)
    
    for class_name, count in sorted_classes:
        percentage = (count / total * 100) if total > 0 else 0
        lines.append(f"- **{class_name}**: {count:,} ({percentage:.1f}%)")
    
    return "\n".join(lines)
