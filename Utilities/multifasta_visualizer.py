"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    MULTIFASTA VISUALIZER MODULE                               ║
║     Unified Visualization and Positional Analysis for MultiFASTA Inputs      ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: multifasta_visualizer.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.4
LICENSE: MIT

DESCRIPTION:
    MultiFASTA visualization module providing:
    - Equal-length sequence detection
    - Unified summary visualizations across sequences
    - Positional occurrence analysis for equal-length datasets
    - Enhanced plotting for comparative analysis
    
    Supports two modes:
    1. MultiFASTA (different lengths): Unified summary plots only
    2. MultiFASTA (equal length): Unified plots + positional panels
    
MAIN FUNCTIONS:
    - all_sequences_equal_length(): Detect if all sequences have same length
    - compute_positional_occurrence(): Calculate motif occurrence by position
    - generate_unified_summary(): Create cross-sequence summary plots
    - generate_positional_panels(): Create position-specific occurrence plots
"""

from __future__ import annotations
from typing import Dict, List, Any, Tuple, Optional
from collections import defaultdict, Counter
import logging
import numpy as np

logger = logging.getLogger(__name__)

# Lazy imports for plotting (performance optimization)
_matplotlib_loaded = False
_plt = None
_sns = None

def _ensure_matplotlib():
    """Lazy import matplotlib/seaborn only when plotting is needed."""
    global _matplotlib_loaded, _plt, _sns
    
    if not _matplotlib_loaded:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import seaborn as sns
            _plt = plt
            _sns = sns
            _matplotlib_loaded = True
        except ImportError as e:
            raise ImportError(f"matplotlib/seaborn required for visualization: {e}")
    
    return _plt, _sns


class MultiFastaVisualizer:
    """
    Visualizer for MultiFASTA analysis with unified summaries and positional panels.
    
    Handles two cases:
    - Case A: Different sequence lengths → unified summary only
    - Case B: Equal sequence lengths → unified summary + positional panels
    """
    
    def __init__(self, annotations_by_sequence: Dict[str, List[Dict[str, Any]]]):
        """
        Initialize MultiFASTA visualizer.
        
        Args:
            annotations_by_sequence: Dict mapping FASTA_ID to list of motif annotations
                Each annotation should have: Start, End, Class, Subclass, etc.
        """
        self.annotations = annotations_by_sequence
        self.fasta_ids = list(annotations_by_sequence.keys())
        self.all_motifs = []
        
        # Flatten all motifs and ensure FASTA_ID is present
        for fasta_id, motifs in annotations_by_sequence.items():
            for motif in motifs:
                m = motif.copy()
                m['FASTA_ID'] = fasta_id
                self.all_motifs.append(m)
    
    def all_sequences_equal_length(self, sequences: Dict[str, str]) -> bool:
        """
        Check if all sequences have the same length.
        
        Args:
            sequences: Dict mapping FASTA_ID to sequence string
            
        Returns:
            True if all sequences have equal length, False otherwise
        """
        if not sequences:
            return False
        
        lengths = {len(seq) for seq in sequences.values()}
        return len(lengths) == 1
    
    def compute_positional_occurrence(
        self, 
        seq_length: int,
        by_class: bool = True,
        by_subclass: bool = True
    ) -> Dict[str, Dict[int, int]]:
        """
        Compute positional occurrence counts for equal-length sequences.
        
        For each position in the sequence (1 to seq_length), count how many
        sequences have a motif at that position.
        
        Args:
            seq_length: Length of sequences (must be equal for all)
            by_class: Compute counts per motif class
            by_subclass: Compute counts per motif subclass
            
        Returns:
            Dict with keys:
                - 'Overall': position -> count across all motifs
                - 'Class': Dict[class_name, Dict[position, count]]
                - 'Subclass': Dict[subclass_name, Dict[position, count]]
        """
        result = {
            'Overall': defaultdict(int),
            'Class': defaultdict(lambda: defaultdict(int)),
            'Subclass': defaultdict(lambda: defaultdict(int))
        }
        
        for motif in self.all_motifs:
            start = motif.get('Start', 0)
            end = motif.get('End', 0)
            motif_class = motif.get('Class', 'Unknown')
            motif_subclass = motif.get('Subclass', 'Unknown')
            
            # Count each position covered by this motif
            for pos in range(start, end + 1):
                if 1 <= pos <= seq_length:
                    result['Overall'][pos] += 1
                    if by_class:
                        result['Class'][motif_class][pos] += 1
                    if by_subclass:
                        result['Subclass'][motif_subclass][pos] += 1
        
        # Convert defaultdicts to regular dicts
        result['Overall'] = dict(result['Overall'])
        result['Class'] = {k: dict(v) for k, v in result['Class'].items()}
        result['Subclass'] = {k: dict(v) for k, v in result['Subclass'].items()}
        
        return result
    
    def generate_unified_summary(self) -> Dict[str, Any]:
        """
        Generate unified summary statistics across all sequences.
        
        Returns:
            Dict containing:
                - sequence_stats: Per-sequence statistics
                - class_distribution: Class-level aggregated stats
                - subclass_distribution: Subclass-level aggregated stats
        """
        sequence_stats = []
        
        # Per-sequence statistics
        for fasta_id in self.fasta_ids:
            motifs = self.annotations[fasta_id]
            
            # Get sequence length (from first motif or context)
            # Note: Caller should provide sequence lengths separately if needed
            total_motifs = len(motifs)
            
            # Class distribution for this sequence
            class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
            
            sequence_stats.append({
                'FASTA_ID': fasta_id,
                'Total_Motifs': total_motifs,
                'Class_Counts': dict(class_counts)
            })
        
        # Class-level aggregation
        all_classes = Counter(m.get('Class', 'Unknown') for m in self.all_motifs)
        
        # Subclass-level aggregation
        all_subclasses = Counter(m.get('Subclass', 'Unknown') for m in self.all_motifs)
        
        return {
            'sequence_stats': sequence_stats,
            'class_distribution': dict(all_classes),
            'subclass_distribution': dict(all_subclasses)
        }
    
    def generate_class_distribution_plot(self, figsize: Tuple[int, int] = (10, 6)):
        """
        Generate stacked bar plot showing motif class distribution per sequence.
        
        Args:
            figsize: Figure size (width, height)
            
        Returns:
            matplotlib Figure object
        """
        plt, sns = _ensure_matplotlib()
        
        # Prepare data for stacked bar chart
        data = []
        for fasta_id in self.fasta_ids:
            motifs = self.annotations[fasta_id]
            class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
            for cls, count in class_counts.items():
                data.append({
                    'FASTA_ID': fasta_id,
                    'Class': cls,
                    'Count': count
                })
        
        if not data:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, 'No motifs to plot', ha='center', va='center')
            return fig
        
        # Create pivot table for stacked bar
        import pandas as pd
        df = pd.DataFrame(data)
        pivot_df = df.pivot_table(index='FASTA_ID', columns='Class', values='Count', fill_value=0)
        
        # Plot stacked bar chart
        fig, ax = plt.subplots(figsize=figsize)
        pivot_df.plot(kind='bar', stacked=True, ax=ax, colormap='tab10')
        
        ax.set_xlabel('Sequence', fontsize=12, fontweight='bold')
        ax.set_ylabel('Motif Count', fontsize=12, fontweight='bold')
        ax.set_title('Motif Class Distribution per Sequence', fontsize=14, fontweight='bold')
        ax.legend(title='Class', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        return fig
    
    def generate_density_heatmap(
        self, 
        sequence_lengths: Dict[str, int],
        figsize: Tuple[int, int] = (10, 6)
    ):
        """
        Generate heatmap showing motif density (motifs per kb) per sequence and class.
        
        Args:
            sequence_lengths: Dict mapping FASTA_ID to sequence length
            figsize: Figure size (width, height)
            
        Returns:
            matplotlib Figure object
        """
        plt, sns = _ensure_matplotlib()
        
        # Prepare data
        data = []
        for fasta_id in self.fasta_ids:
            motifs = self.annotations[fasta_id]
            seq_length = sequence_lengths.get(fasta_id, 1)
            
            class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
            for cls, count in class_counts.items():
                density = (count / seq_length * 1000) if seq_length > 0 else 0
                data.append({
                    'FASTA_ID': fasta_id,
                    'Class': cls,
                    'Density': density
                })
        
        if not data:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, 'No motifs to plot', ha='center', va='center')
            return fig
        
        # Create pivot table for heatmap
        import pandas as pd
        df = pd.DataFrame(data)
        pivot_df = df.pivot_table(index='FASTA_ID', columns='Class', values='Density', fill_value=0)
        
        # Plot heatmap
        fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(pivot_df, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax, cbar_kws={'label': 'Motifs/kb'})
        
        ax.set_xlabel('Motif Class', fontsize=12, fontweight='bold')
        ax.set_ylabel('Sequence', fontsize=12, fontweight='bold')
        ax.set_title('Motif Density Heatmap (Motifs per kb)', fontsize=14, fontweight='bold')
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        plt.tight_layout()
        
        return fig
    
    def generate_positional_panels(
        self,
        seq_length: int,
        window_size: int = 1,
        smooth: bool = True,
        figsize: Tuple[int, int] = (12, 8)
    ):
        """
        Generate positional occurrence panels for equal-length sequences.
        
        Creates separate plots for:
        - Overall positional occurrence
        - Per-class positional occurrence
        
        Args:
            seq_length: Length of sequences (must be equal for all)
            window_size: Window size for smoothing (default: 1 = no smoothing)
            smooth: Whether to apply smoothing
            figsize: Figure size (width, height)
            
        Returns:
            Dict of matplotlib Figure objects:
                - 'overall': Overall positional plot
                - 'by_class': Per-class positional plot
        """
        plt, sns = _ensure_matplotlib()
        
        # Compute positional occurrence
        pos_data = self.compute_positional_occurrence(seq_length)
        
        figures = {}
        
        # Overall positional occurrence
        if pos_data['Overall']:
            fig, ax = plt.subplots(figsize=figsize)
            
            positions = sorted(pos_data['Overall'].keys())
            counts = [pos_data['Overall'][p] for p in positions]
            
            if smooth and len(positions) > 10:
                # Apply rolling mean smoothing
                window = min(window_size, len(counts) // 10)
                if window > 1:
                    counts = np.convolve(counts, np.ones(window)/window, mode='same')
            
            ax.plot(positions, counts, linewidth=2, color='#1f77b4')
            ax.fill_between(positions, counts, alpha=0.3, color='#1f77b4')
            
            ax.set_xlabel('Position (bp)', fontsize=12, fontweight='bold')
            ax.set_ylabel('Number of Sequences with Motif', fontsize=12, fontweight='bold')
            ax.set_title('Overall Positional Occurrence', fontsize=14, fontweight='bold')
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            
            figures['overall'] = fig
        
        # Per-class positional occurrence
        if pos_data['Class']:
            num_classes = len(pos_data['Class'])
            fig, axes = plt.subplots(num_classes, 1, figsize=(figsize[0], figsize[1] * num_classes / 2))
            
            if num_classes == 1:
                axes = [axes]
            
            for idx, (cls, class_data) in enumerate(sorted(pos_data['Class'].items())):
                ax = axes[idx]
                
                if class_data:
                    positions = sorted(class_data.keys())
                    counts = [class_data[p] for p in positions]
                    
                    if smooth and len(positions) > 10:
                        window = min(window_size, len(counts) // 10)
                        if window > 1:
                            counts = np.convolve(counts, np.ones(window)/window, mode='same')
                    
                    ax.plot(positions, counts, linewidth=2)
                    ax.fill_between(positions, counts, alpha=0.3)
                    
                    ax.set_ylabel('Count', fontsize=10, fontweight='bold')
                    ax.set_title(f'{cls.replace("_", " ")}', fontsize=11, fontweight='bold')
                    ax.grid(True, alpha=0.3)
            
            axes[-1].set_xlabel('Position (bp)', fontsize=12, fontweight='bold')
            plt.suptitle('Positional Occurrence by Motif Class', fontsize=14, fontweight='bold', y=1.00)
            plt.tight_layout()
            
            figures['by_class'] = fig
        
        return figures


def prepare_multifasta_excel_data(
    annotations_by_sequence: Dict[str, List[Dict[str, Any]]],
    sequence_lengths: Dict[str, int],
    equal_length: bool = False,
    seq_length: Optional[int] = None
) -> Dict[str, Any]:
    """
    Prepare data for MultiFASTA Excel export with all required sheets.
    
    Args:
        annotations_by_sequence: Dict mapping FASTA_ID to motif annotations
        sequence_lengths: Dict mapping FASTA_ID to sequence length
        equal_length: Whether all sequences have equal length
        seq_length: Sequence length if equal_length is True
        
    Returns:
        Dict with sheet data:
            - All_Motifs: All motifs with FASTA_ID column
            - Sequence_Summary: Per-sequence statistics
            - Class_Summary: Aggregated class statistics
            - Positional_Occurrence: Position-level counts (if equal_length=True)
    """
    visualizer = MultiFastaVisualizer(annotations_by_sequence)
    
    # Sheet 1: All motifs with FASTA_ID
    all_motifs = visualizer.all_motifs
    
    # Sheet 2: Sequence summary
    sequence_summary = []
    for fasta_id in visualizer.fasta_ids:
        motifs = annotations_by_sequence[fasta_id]
        seq_len = sequence_lengths.get(fasta_id, 0)
        total_motifs = len(motifs)
        motifs_per_kb = (total_motifs / seq_len * 1000) if seq_len > 0 else 0
        
        sequence_summary.append({
            'FASTA_ID': fasta_id,
            'Length': seq_len,
            'Total_Motifs': total_motifs,
            'Motifs_per_kb': round(motifs_per_kb, 2)
        })
    
    # Sheet 3: Class summary
    class_summary = []
    class_counts = Counter(m.get('Class', 'Unknown') for m in all_motifs)
    subclass_counts = Counter((m.get('Class', 'Unknown'), m.get('Subclass', 'Unknown')) 
                              for m in all_motifs)
    
    for (cls, subcls), count in sorted(subclass_counts.items()):
        class_summary.append({
            'Class': cls,
            'Subclass': subcls,
            'Total_Count': count
        })
    
    result = {
        'All_Motifs': all_motifs,
        'Sequence_Summary': sequence_summary,
        'Class_Summary': class_summary
    }
    
    # Sheet 4: Positional occurrence (only if equal length)
    if equal_length and seq_length:
        pos_data = visualizer.compute_positional_occurrence(seq_length)
        
        positional_occurrence = []
        for cls, class_data in sorted(pos_data['Class'].items()):
            for subcls, subclass_data in sorted(pos_data['Subclass'].items()):
                if subcls.startswith(cls.replace('_', ' ')) or subcls.startswith(cls):
                    for pos, count in sorted(subclass_data.items()):
                        positional_occurrence.append({
                            'Position': pos,
                            'Class': cls,
                            'Subclass': subcls,
                            'Count': count
                        })
        
        result['Positional_Occurrence'] = positional_occurrence
    
    return result
