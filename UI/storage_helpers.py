"""
Helper module for accessing sequences and results with backward compatibility.

Provides unified API that works with both disk storage and legacy in-memory modes.
"""

import streamlit as st
from typing import List, Dict, Any, Optional, Tuple
import pandas as pd
import logging

logger = logging.getLogger(__name__)


def get_sequences_info() -> Tuple[List[str], List[int], int]:
    """
    Get sequences information (names, lengths, count).
    
    Returns:
        Tuple of (names, lengths, count)
    """
    if st.session_state.get('use_disk_storage') and st.session_state.get('seq_ids'):
        # Disk storage mode
        names = st.session_state.get("names", [])
        lengths = []
        for seq_id in st.session_state.seq_ids:
            metadata = st.session_state.seq_storage.get_metadata(seq_id)
            lengths.append(metadata['length'])
        count = len(st.session_state.seq_ids)
    else:
        # Legacy mode
        names = st.session_state.get("names", [])
        seqs = st.session_state.get("seqs", [])
        lengths = [len(seq) for seq in seqs]
        count = len(seqs)
    
    return names, lengths, count


def get_sequence_length(seq_idx: int) -> int:
    """
    Get length of a specific sequence.
    
    Args:
        seq_idx: Sequence index
        
    Returns:
        Sequence length in base pairs
    """
    if st.session_state.get('use_disk_storage') and st.session_state.get('seq_ids'):
        # Disk storage mode - check bounds before access
        if seq_idx >= len(st.session_state.seq_ids):
            logger.error(f"Sequence index {seq_idx} out of range "
                        f"(have {len(st.session_state.seq_ids)} sequences)")
            return 0
        
        seq_id = st.session_state.seq_ids[seq_idx]
        metadata = st.session_state.seq_storage.get_metadata(seq_id)
        return metadata['length']
    else:
        # Legacy mode - check bounds before access
        if seq_idx >= len(st.session_state.seqs):
            logger.error(f"Sequence index {seq_idx} out of range "
                        f"(have {len(st.session_state.seqs)} sequences)")
            return 0
        
        return len(st.session_state.seqs[seq_idx])


def get_results(seq_idx: int, limit: Optional[int] = None) -> List[Dict[str, Any]]:
    """
    Get analysis results for a specific sequence.
    
    Args:
        seq_idx: Sequence index
        limit: Optional limit on number of results to return
        
    Returns:
        List of motif dictionaries
    """
    if st.session_state.get('use_disk_storage') and st.session_state.get('seq_ids'):
        # Disk storage mode - check bounds before access
        if seq_idx >= len(st.session_state.seq_ids):
            logger.error(f"Results index {seq_idx} out of range")
            return []
        
        seq_id = st.session_state.seq_ids[seq_idx]
        
        if seq_id in st.session_state.results_storage:
            results_storage = st.session_state.results_storage[seq_id]
            return list(results_storage.iter_results(limit=limit))
        else:
            return []
    else:
        # Legacy mode - already has bounds check
        results_list = st.session_state.get("results", [])
        if seq_idx < len(results_list):
            results = results_list[seq_idx]
            if limit is not None:
                return results[:limit]
            return results
        return []


def get_results_summary(seq_idx: int) -> Dict[str, Any]:
    """
    Get summary statistics for results without loading all motifs.
    
    Args:
        seq_idx: Sequence index
        
    Returns:
        Dictionary with summary statistics
    """
    if st.session_state.get('use_disk_storage') and st.session_state.get('seq_ids'):
        # Disk storage mode
        seq_id = st.session_state.seq_ids[seq_idx]
        
        if seq_id in st.session_state.results_storage:
            results_storage = st.session_state.results_storage[seq_id]
            return results_storage.get_summary_stats()
        else:
            return {
                'total_count': 0,
                'class_distribution': {},
                'subclass_distribution': {},
                'coverage_bp': 0,
                'avg_score': 0.0
            }
    else:
        # Legacy mode - compute from in-memory results
        if seq_idx < len(st.session_state.results):
            motifs = st.session_state.results[seq_idx]
            
            # Compute statistics
            class_dist = {}
            subclass_dist = {}
            coverage_bp = 0
            score_sum = 0.0
            
            for motif in motifs:
                motif_class = motif.get('Class', 'Unknown')
                motif_subclass = motif.get('Subclass', 'Unknown')
                
                class_dist[motif_class] = class_dist.get(motif_class, 0) + 1
                subclass_dist[motif_subclass] = subclass_dist.get(motif_subclass, 0) + 1
                
                coverage_bp += motif.get('Length', 0)
                score_sum += motif.get('Score', 0.0)
            
            return {
                'total_count': len(motifs),
                'class_distribution': class_dist,
                'subclass_distribution': subclass_dist,
                'coverage_bp': coverage_bp,
                'avg_score': score_sum / len(motifs) if motifs else 0.0
            }
        
        return {
            'total_count': 0,
            'class_distribution': {},
            'subclass_distribution': {},
            'coverage_bp': 0,
            'avg_score': 0.0
        }


def has_results() -> bool:
    """
    Stable, single source-of-truth result check.

    Results are considered available ONLY if:
    - analysis_done flag is True
    - results object exists in session_state

    This prevents false negatives during Streamlit reruns
    (e.g., after download button clicks).
    """
    return (
        st.session_state.get("analysis_done", False)
        and "results" in st.session_state
    )


def get_results_dataframe(seq_idx: int, limit: Optional[int] = None) -> pd.DataFrame:
    """
    Get results as pandas DataFrame.
    
    Args:
        seq_idx: Sequence index
        limit: Optional limit on number of results
        
    Returns:
        DataFrame with motif data
    """
    motifs = get_results(seq_idx, limit=limit)
    return pd.DataFrame(motifs) if motifs else pd.DataFrame()
