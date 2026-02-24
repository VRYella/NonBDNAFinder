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
from typing import Dict, List, Any, Tuple, Optional, Callable
from collections import defaultdict, Counter
import logging
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
from Utilities.config.colors import UNIFIED_MOTIF_COLORS

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
        # Use numpy arrays for O(motif_count) accumulation instead of the
        # former per-position Python loop which was O(total_covered_bp).
        overall_arr = np.zeros(seq_length, dtype=np.int32)
        class_arrs: Dict[str, np.ndarray] = {}
        subclass_arrs: Dict[str, np.ndarray] = {}

        for motif in self.all_motifs:
            start = max(0, motif.get('Start', 1) - 1)   # 0-based
            end = min(motif.get('End', 0), seq_length)   # exclusive
            if end <= start:
                continue

            overall_arr[start:end] += 1

            if by_class:
                cls = motif.get('Class', 'Unknown')
                if cls not in class_arrs:
                    class_arrs[cls] = np.zeros(seq_length, dtype=np.int32)
                class_arrs[cls][start:end] += 1

            if by_subclass:
                subcls = motif.get('Subclass', 'Unknown')
                if subcls not in subclass_arrs:
                    subclass_arrs[subcls] = np.zeros(seq_length, dtype=np.int32)
                subclass_arrs[subcls][start:end] += 1

        # Convert numpy arrays to sparse {position: count} dicts (1-based positions)
        def _arr_to_dict(arr: np.ndarray) -> Dict[int, int]:
            nonzero = np.nonzero(arr)[0]
            return {int(i) + 1: int(arr[i]) for i in nonzero}

        result = {
            'Overall': _arr_to_dict(overall_arr),
            'Class': {cls: _arr_to_dict(arr) for cls, arr in class_arrs.items()},
            'Subclass': {subcls: _arr_to_dict(arr) for subcls, arr in subclass_arrs.items()},
        }
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
        
        # Plot stacked bar chart using unified motif class colors
        fig, ax = plt.subplots(figsize=figsize)
        all_classes = pivot_df.columns.tolist()
        colors = [UNIFIED_MOTIF_COLORS.get(cls, '#808080') for cls in all_classes]
        pivot_df.plot(kind='bar', stacked=True, ax=ax, color=colors)
        
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
                    
                    color = UNIFIED_MOTIF_COLORS.get(cls, '#808080')
                    ax.plot(positions, counts, linewidth=2, color=color)
                    ax.fill_between(positions, counts, alpha=0.3, color=color)
                    
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
        
        # Collect data by class and subclass
        for cls, class_data in sorted(pos_data['Class'].items()):
            # Get subclasses for this class
            subclass_data_dict = {}
            for subcls, subclass_positions in sorted(pos_data['Subclass'].items()):
                # Match subclass to class (simple heuristic: check if subclass contains class name)
                cls_normalized = cls.replace('-', '_').replace(' ', '_')
                subcls_normalized = subcls.replace('-', '_').replace(' ', '_')
                
                # Check if they're related (this is a simple check)
                # In a real scenario, we'd have a proper mapping
                if cls_normalized in subcls_normalized or cls == 'Unknown':
                    if cls not in subclass_data_dict:
                        subclass_data_dict[subcls] = subclass_positions
            
            # If no subclass data, use class data directly
            if not subclass_data_dict:
                for pos, count in sorted(class_data.items()):
                    positional_occurrence.append({
                        'Position': pos,
                        'Class': cls,
                        'Subclass': 'Unknown',
                        'Count': count
                    })
            else:
                # Add entries for each subclass
                for subcls, subclass_positions in subclass_data_dict.items():
                    for pos, count in sorted(subclass_positions.items()):
                        positional_occurrence.append({
                            'Position': pos,
                            'Class': cls,
                            'Subclass': subcls,
                            'Count': count
                        })
        
        if positional_occurrence:
            result['Positional_Occurrence'] = positional_occurrence
    
    return result


# =============================================================================
# PARALLEL VISUALIZATION GENERATION
# =============================================================================

def _generate_plot_worker(plot_spec: Dict[str, Any]) -> Dict[str, Any]:
    """
    Worker function to generate a single plot in parallel.
    
    Args:
        plot_spec: Dictionary containing:
            - plot_type: Type of plot to generate
            - plot_params: Parameters for the plot function
            - plot_name: Name/identifier for the plot
    
    Returns:
        Dict with plot result and metadata
    """
    import time
    
    try:
        start_time = time.time()
        
        plot_type = plot_spec['plot_type']
        # Create a copy to avoid mutating the input dictionary
        plot_params = plot_spec['plot_params'].copy()
        plot_name = plot_spec['plot_name']
        
        visualizer = plot_params.pop('visualizer')
        
        # Generate the appropriate plot
        if plot_type == 'class_distribution':
            fig = visualizer.generate_class_distribution_plot(**plot_params)
        elif plot_type == 'density_heatmap':
            fig = visualizer.generate_density_heatmap(**plot_params)
        elif plot_type == 'positional_panels':
            figs_dict = visualizer.generate_positional_panels(**plot_params)
            elapsed = time.time() - start_time
            return {
                'success': True,
                'plot_name': plot_name,
                'figures': figs_dict,
                'elapsed': elapsed
            }
        else:
            raise ValueError(f"Unknown plot type: {plot_type}")
        
        elapsed = time.time() - start_time
        
        return {
            'success': True,
            'plot_name': plot_name,
            'figure': fig,
            'elapsed': elapsed
        }
    
    except Exception as e:
        logger.error(f"Error generating plot {plot_name}: {e}")
        return {
            'success': False,
            'plot_name': plot_name,
            'error': str(e)
        }


def generate_visualizations_parallel(
    visualizer: 'MultiFastaVisualizer',
    plot_specs: List[Dict[str, Any]],
    max_workers: Optional[int] = None,
    progress_callback: Optional[Callable[[int, int, str], None]] = None
) -> Dict[str, Any]:
    """
    Generate multiple visualizations in parallel using ThreadPoolExecutor.
    
    Args:
        visualizer: MultiFastaVisualizer instance
        plot_specs: List of plot specifications, each containing:
            - plot_type: Type of plot ('class_distribution', 'density_heatmap', 'positional_panels')
            - plot_params: Parameters for the plot function
            - plot_name: Identifier for the plot
        max_workers: Maximum number of worker threads (default: CPU count)
        progress_callback: Optional callback function(completed, total, plot_name)
    
    Returns:
        Dict mapping plot names to figures
    """
    if not plot_specs:
        return {}
    
    # Determine number of workers
    if max_workers is None:
        max_workers = min(len(plot_specs), os.cpu_count() or 4)
    
    # Inject visualizer into each plot_spec
    for spec in plot_specs:
        if 'plot_params' not in spec:
            spec['plot_params'] = {}
        spec['plot_params']['visualizer'] = visualizer
    
    results = {}
    completed_count = 0
    total_count = len(plot_specs)
    
    logger.info(f"Starting parallel generation of {total_count} plots with {max_workers} workers")
    
    try:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all tasks
            future_to_spec = {
                executor.submit(_generate_plot_worker, spec): spec
                for spec in plot_specs
            }
            
            # Process completed tasks as they finish
            for future in as_completed(future_to_spec):
                try:
                    result = future.result()
                    completed_count += 1
                    
                    # Call progress callback if provided
                    if progress_callback:
                        progress_callback(completed_count, total_count, result['plot_name'])
                    
                    if result['success']:
                        if 'figures' in result:
                            # Multiple figures (e.g., positional panels)
                            for key, fig in result['figures'].items():
                                results[f"{result['plot_name']}_{key}"] = fig
                        else:
                            # Single figure
                            results[result['plot_name']] = result['figure']
                        
                        logger.info(f"Completed plot {completed_count}/{total_count}: {result['plot_name']} ({result['elapsed']:.2f}s)")
                    else:
                        logger.error(f"Failed plot {completed_count}/{total_count}: {result['plot_name']} - {result.get('error')}")
                
                except Exception as e:
                    logger.error(f"Error processing future: {e}")
                    completed_count += 1
    
    except Exception as e:
        logger.error(f"Error in parallel visualization: {e}")
        raise
    
    logger.info(f"Parallel visualization complete: {total_count} plots generated")
    return results
