"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    PARALLEL ANALYSIS HELPER MODULE                            ║
║     Helper functions for parallel multi-FASTA analysis in Streamlit UI       ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: parallel_analysis_helper.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.6
LICENSE: MIT

DESCRIPTION:
    Helper module to simplify parallel processing integration in the Streamlit UI.
    Provides high-level wrapper functions that handle both disk storage and
    in-memory modes with consistent interfaces.
"""

from __future__ import annotations
from typing import Dict, List, Any, Optional, Callable
import logging
import os

logger = logging.getLogger(__name__)


def should_use_parallel(num_sequences: int, threshold: int = 2) -> bool:
    """
    Determine if parallel processing should be used based on number of sequences.
    
    Args:
        num_sequences: Number of sequences to process
        threshold: Minimum number of sequences to trigger parallel processing
    
    Returns:
        True if parallel processing should be used
    """
    return num_sequences >= threshold


def prepare_parallel_analysis(
    num_sequences: int,
    use_disk_storage: bool,
    seq_ids: Optional[List[str]] = None,
    seqs: Optional[List[str]] = None,
    names: Optional[List[str]] = None,
    seq_storage: Optional[Any] = None,
    enabled_classes: Optional[List[str]] = None,
    chunk_threshold: int = 1_000_000
) -> tuple:
    """
    Prepare data structures for parallel analysis.
    
    Args:
        num_sequences: Number of sequences
        use_disk_storage: Whether using disk storage mode
        seq_ids: List of sequence IDs (for disk storage)
        seqs: List of sequences (for in-memory)
        names: List of sequence names
        seq_storage: UniversalSequenceStorage instance (for disk storage)
        enabled_classes: List of enabled class names
        chunk_threshold: Threshold for chunking large sequences
    
    Returns:
        Tuple of (sequences_data, analysis_params)
    """
    if use_disk_storage:
        # Disk storage mode
        if not seq_ids or not names:
            raise ValueError("seq_ids and names required for disk storage mode")
        
        analysis_params = {
            'use_disk_storage': True,
            'seq_storage': seq_storage,
            'enabled_classes': enabled_classes,
            'chunk_threshold': chunk_threshold
        }
        
        sequences_data = [(seq_id, name, i) for i, (seq_id, name) in enumerate(zip(seq_ids, names))]
    else:
        # In-memory mode
        if not seqs or not names:
            raise ValueError("seqs and names required for in-memory mode")
        
        analysis_params = {
            'use_disk_storage': False,
            'enabled_classes': enabled_classes,
            'chunk_threshold': chunk_threshold
        }
        
        sequences_data = [(seq, name, i) for i, (seq, name) in enumerate(zip(seqs, names))]
    
    return sequences_data, analysis_params


def get_optimal_workers(num_sequences: int) -> int:
    """
    Determine optimal number of worker threads/processes.
    
    Args:
        num_sequences: Number of sequences to process
    
    Returns:
        Optimal number of workers
    """
    cpu_count = os.cpu_count() or 4
    # Use minimum of sequences count and CPU count, but at least 2
    return max(2, min(num_sequences, cpu_count))
