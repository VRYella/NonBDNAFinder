"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    MULTIFASTA AGGREGATION ENGINE                              ║
║     Scalable Aggregation-First Engine for Large MultiFASTA Datasets          ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: multifasta_engine.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.5
LICENSE: MIT

DESCRIPTION:
    Aggregation-first MultiFASTA engine designed for scalability (>20 sequences).
    No full DataFrame materialization to maintain memory efficiency.
    
    Key Features:
    - Global class distribution across all sequences
    - % sequence coverage per motif class
    - Motif density histogram (motifs per kb)
    - Positional conservation analysis (equal-length only)
    - Memory-safe windowed aggregation
    - 100MB ceiling compatible
    
USAGE:
    engine = MultiFastaEngine(annotations_by_sequence, sequence_lengths)
    
    # Get global statistics
    class_dist = engine.global_class_distribution()
    coverage = engine.class_sequence_coverage()
    densities = engine.motif_density_per_sequence()
    
    # For equal-length sequences
    if all_equal_length:
        conservation = engine.positional_conservation(seq_length)
"""

from __future__ import annotations
from typing import Dict, List, Any, Optional, Callable
from collections import defaultdict
import math
import logging
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import multiprocessing
import os

logger = logging.getLogger(__name__)


class MultiFastaEngine:
    """
    Aggregation-first MultiFASTA engine.
    Designed for scalability (>20 sequences).
    No full DataFrame materialization.
    """

    def __init__(self, annotations_by_sequence: Dict[str, List[Dict[str, Any]]], 
                 sequence_lengths: Dict[str, int]):
        """
        Initialize MultiFASTA aggregation engine.
        
        Args:
            annotations_by_sequence: Dict mapping sequence ID to list of motif annotations
            sequence_lengths: Dict mapping sequence ID to sequence length
        """
        self.annotations = annotations_by_sequence
        self.sequence_lengths = sequence_lengths
        self.num_sequences = len(sequence_lengths)

    # ---------------------------------------------------------------------
    # GLOBAL CLASS DISTRIBUTION
    # ---------------------------------------------------------------------
    def global_class_distribution(self) -> Dict[str, int]:
        """
        Compute global motif class distribution across all sequences.
        
        Returns:
            Dict mapping motif class to total count across all sequences
        """
        distribution = defaultdict(int)
        for motifs in self.annotations.values():
            for motif in motifs:
                distribution[motif.get("Class")] += 1
        return dict(distribution)

    # ---------------------------------------------------------------------
    # % SEQUENCES CONTAINING EACH CLASS
    # ---------------------------------------------------------------------
    def class_sequence_coverage(self) -> Dict[str, float]:
        """
        Calculate % of sequences containing each motif class.
        
        Returns:
            Dict mapping motif class to percentage of sequences containing it
        """
        coverage = defaultdict(int)

        for seq, motifs in self.annotations.items():
            seen_classes = set()
            for motif in motifs:
                seen_classes.add(motif.get("Class"))

            for cls in seen_classes:
                coverage[cls] += 1

        # convert to percentage
        return {
            cls: (count / self.num_sequences) * 100
            for cls, count in coverage.items()
        }

    # ---------------------------------------------------------------------
    # MOTIF DENSITY PER KB (per sequence)
    # ---------------------------------------------------------------------
    def motif_density_per_sequence(self) -> List[float]:
        """
        Calculate motif density (motifs per kb) for each sequence.
        
        Returns:
            List of density values (one per sequence)
        """
        densities = []

        for seq, motifs in self.annotations.items():
            length = self.sequence_lengths.get(seq, 1)
            if length > 0:
                density = len(motifs) / (length / 1000)
            else:
                density = 0.0
            densities.append(density)

        return densities

    # ---------------------------------------------------------------------
    # POSITIONAL CONSERVATION (equal-length only)
    # ---------------------------------------------------------------------
    def positional_conservation(self, seq_length: int) -> List[int]:
        """
        Calculate positional conservation for equal-length sequences.
        Returns number of sequences containing motif at each position.
        Memory safe: uses windowed aggregation.
        
        Args:
            seq_length: Length of sequences (must be equal for all)
            
        Returns:
            List of conservation counts per window (1000bp windows)
        """
        window = 1000
        bins = math.ceil(seq_length / window)
        conservation = [0] * bins

        for seq, motifs in self.annotations.items():
            covered_bins = set()

            for motif in motifs:
                start_bin = motif["Start"] // window
                end_bin = motif["End"] // window
                for b in range(start_bin, min(end_bin + 1, bins)):
                    covered_bins.add(b)

            for b in covered_bins:
                conservation[b] += 1

        return conservation


# =============================================================================
# PARALLEL PROCESSING FUNCTIONS FOR MULTI-FASTA ANALYSIS
# =============================================================================

def _analyze_single_sequence_worker(args):
    """
    Worker function for parallel sequence analysis.
    This function is called by ProcessPoolExecutor/ThreadPoolExecutor.
    
    Args:
        args: Tuple of (seq_or_seq_id, name, index, analysis_params)
            - seq_or_seq_id: Either the sequence string or sequence ID for disk storage
            - name: Sequence name
            - index: Sequence index
            - analysis_params: Dict with analysis configuration
    
    Returns:
        Dict with analysis results and metadata
    """
    from Utilities.nonbscanner import analyze_sequence
    from Utilities.chunk_analyzer import ChunkAnalyzer
    from Utilities.disk_storage import UniversalResultsStorage
    import time
    
    seq_or_seq_id, name, index, analysis_params = args
    
    try:
        start_time = time.time()
        
        # Extract parameters
        use_disk_storage = analysis_params.get('use_disk_storage', False)
        seq_storage = analysis_params.get('seq_storage')
        enabled_classes = analysis_params.get('enabled_classes')
        chunk_threshold = analysis_params.get('chunk_threshold', 1_000_000)
        
        if use_disk_storage:
            # Disk storage mode
            seq_id = seq_or_seq_id
            metadata = seq_storage.get_metadata(seq_id)
            seq_length = metadata['length']
            
            # Determine whether to use chunking
            if seq_length > chunk_threshold:
                # Use chunk analyzer
                # NOTE: use_parallel=False to avoid nested parallelism since we're already
                # in a parallel worker. The sequence-level parallelism is sufficient.
                analyzer = ChunkAnalyzer(
                    seq_storage,
                    use_parallel=False,  # Disable nested parallelism
                    max_workers=None,
                    use_adaptive=False
                )
                
                results_storage = analyzer.analyze(
                    seq_id=seq_id,
                    progress_callback=None,  # No callback in parallel mode
                    enabled_classes=enabled_classes
                )
                
                # Get results
                results = list(results_storage.iter_results())
                total_motifs = results_storage.get_summary_stats()['total_count']
            else:
                # Standard analysis
                seq = seq_storage.get_sequence_chunk(seq_id, 0, seq_length)
                results = analyze_sequence(seq, name, enabled_classes=enabled_classes)
                
                # Create results storage and save results
                # NOTE: This path handles non-chunked disk storage scenarios
                # Test coverage: Requires test with disk storage + sequence < chunk_threshold
                results_storage = UniversalResultsStorage(
                    base_dir=str(seq_storage.base_dir / "results"),
                    seq_id=seq_id
                )
                results_storage.append_batch(results)
            
            elapsed = time.time() - start_time
            
            return {
                'success': True,
                'index': index,
                'name': name,
                'seq_id': seq_id,
                'seq_length': seq_length,
                'results': results,
                'results_storage': results_storage,
                'elapsed': elapsed,
                'use_disk_storage': True
            }
        else:
            # Legacy in-memory mode
            seq = seq_or_seq_id
            seq_length = len(seq)
            
            results = analyze_sequence(seq, name, enabled_classes=enabled_classes)
            elapsed = time.time() - start_time
            
            return {
                'success': True,
                'index': index,
                'name': name,
                'seq_length': seq_length,
                'results': results,
                'elapsed': elapsed,
                'use_disk_storage': False
            }
            
    except Exception as e:
        logger.error(f"Error analyzing sequence {name}: {e}")
        return {
            'success': False,
            'index': index,
            'name': name,
            'error': str(e)
        }


def analyze_sequences_parallel(
    sequences_data: List[tuple],
    analysis_params: Dict[str, Any],
    max_workers: Optional[int] = None,
    progress_callback: Optional[Callable[[int, int, str], None]] = None,
    use_processes: bool = False
) -> List[Dict[str, Any]]:
    """
    Analyze multiple sequences in parallel using ThreadPoolExecutor or ProcessPoolExecutor.
    
    Args:
        sequences_data: List of tuples (seq_or_seq_id, name, index)
            - For disk storage: (seq_id, name, index)
            - For in-memory: (sequence_string, name, index)
        analysis_params: Dict with analysis configuration:
            - use_disk_storage: bool
            - seq_storage: UniversalSequenceStorage (if use_disk_storage=True)
            - enabled_classes: List of class names to analyze
            - chunk_threshold: Threshold for chunking (default: 1MB)
        max_workers: Maximum number of worker threads/processes (default: CPU count)
        progress_callback: Optional callback function(completed, total, seq_name)
        use_processes: If True, use ProcessPoolExecutor; if False, use ThreadPoolExecutor
            Default is False (threads) because most detectors are I/O bound
    
    Returns:
        List of result dictionaries, ordered by sequence index
    """
    if not sequences_data:
        return []
    
    # Determine number of workers
    if max_workers is None:
        max_workers = min(len(sequences_data), os.cpu_count() or 4)
    
    # Prepare arguments for workers
    worker_args = [
        (seq_or_seq_id, name, idx, analysis_params)
        for seq_or_seq_id, name, idx in sequences_data
    ]
    
    # Choose executor based on use_processes flag
    ExecutorClass = ProcessPoolExecutor if use_processes else ThreadPoolExecutor
    
    results_dict = {}
    completed_count = 0
    total_count = len(worker_args)
    
    logger.info(f"Starting parallel analysis of {total_count} sequences with {max_workers} workers")
    
    try:
        with ExecutorClass(max_workers=max_workers) as executor:
            # Submit all tasks
            future_to_args = {
                executor.submit(_analyze_single_sequence_worker, args): args
                for args in worker_args
            }
            
            # Process completed tasks as they finish
            for future in as_completed(future_to_args):
                try:
                    result = future.result()
                    results_dict[result['index']] = result
                    completed_count += 1
                    
                    # Call progress callback if provided
                    if progress_callback:
                        progress_callback(completed_count, total_count, result['name'])
                    
                    if result['success']:
                        logger.info(f"Completed {completed_count}/{total_count}: {result['name']} "
                                  f"({result['seq_length']:,} bp, {len(result.get('results', [])):,} motifs, "
                                  f"{result['elapsed']:.2f}s)")
                    else:
                        logger.error(f"Failed {completed_count}/{total_count}: {result['name']} - {result.get('error')}")
                
                except Exception as e:
                    logger.error(f"Error processing future: {e}")
                    completed_count += 1
    
    except Exception as e:
        logger.error(f"Error in parallel execution: {e}")
        raise
    
    # Return results in original order
    ordered_results = [results_dict[i] for i in sorted(results_dict.keys())]
    
    logger.info(f"Parallel analysis complete: {total_count} sequences processed")
    return ordered_results
