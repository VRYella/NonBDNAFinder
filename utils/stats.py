"""
Statistics and Memory Management Module
========================================

Functions for calculating sequence statistics and managing memory usage.
Extracted from utilities.py for modular architecture.
"""

from typing import Dict, Any, List, Optional
import gc
import numpy as np
import pandas as pd
from collections import Counter


def trigger_garbage_collection():
    """
    Trigger explicit garbage collection to free memory.
    
    Should be called at critical pipeline stages:
    - After processing large sequences
    - After generating visualizations
    - Before/after large data exports
    
    Returns:
        int: Number of objects collected
    """
    collected = gc.collect()
    return collected


def optimize_dataframe_memory(df: pd.DataFrame) -> pd.DataFrame:
    """
    Optimize pandas DataFrame memory usage by downcasting numeric types.
    
    Reduces memory footprint for large result datasets without losing precision.
    Note: float64 → float32 conversion may affect precision for some calculations.
    
    Args:
        df: Input DataFrame
        
    Returns:
        Memory-optimized DataFrame (a copy, original unchanged)
        
    Example:
        >>> df = optimize_dataframe_memory(large_motif_df)
        >>> # Overall memory usage typically reduced by 50-70%
        >>> # (Varies based on data types: integers optimize more than floats)
    """
    # Create a copy to avoid modifying original
    df_optimized = df.copy()
    
    # Collect all type changes first
    type_changes = {}
    
    for col in df_optimized.columns:
        col_type = df_optimized[col].dtype
        
        # Optimize integer columns
        if col_type == 'int64':
            c_min = df_optimized[col].min()
            c_max = df_optimized[col].max()
            if c_min >= 0:
                # Unsigned integers (boundary inclusive)
                if c_max <= 255:
                    type_changes[col] = np.uint8
                elif c_max <= 65535:
                    type_changes[col] = np.uint16
                elif c_max <= 4294967295:
                    type_changes[col] = np.uint32
            else:
                # Signed integers (boundary inclusive)
                if c_min >= np.iinfo(np.int8).min and c_max <= np.iinfo(np.int8).max:
                    type_changes[col] = np.int8
                elif c_min >= np.iinfo(np.int16).min and c_max <= np.iinfo(np.int16).max:
                    type_changes[col] = np.int16
                elif c_min >= np.iinfo(np.int32).min and c_max <= np.iinfo(np.int32).max:
                    type_changes[col] = np.int32
        
        # Optimize float columns (with precision note)
        # Note: float32 provides ~7 significant digits, sufficient for most genomic scores
        elif col_type == 'float64':
            type_changes[col] = np.float32
    
    # Apply all type changes
    for col, new_type in type_changes.items():
        df_optimized[col] = df_optimized[col].astype(new_type)
            
    return df_optimized


def get_memory_usage_mb() -> float:
    """
    Get current process memory usage in MB.
    
    Returns:
        Memory usage in megabytes
    """
    try:
        import psutil
        import os
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024
    except ImportError:
        return 0.0


def gc_content(sequence: str) -> float:
    """Calculate GC content percentage"""
    if not sequence:
        return 0.0
    seq = sequence.upper()
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq)) * 100


def at_content(sequence: str) -> float:
    """Calculate AT content percentage"""
    if not sequence:
        return 0.0
    seq = sequence.upper()
    at_count = seq.count('A') + seq.count('T')
    return (at_count / len(seq)) * 100


def calculate_tm(sequence: str) -> float:
    """
    Calculate melting temperature using nearest-neighbor method approximation.
    For short sequences (<14bp), uses Wallace rule.
    """
    if not sequence:
        return 0.0
    
    seq = sequence.upper()
    length = len(seq)
    
    if length < 14:
        # Wallace rule for short sequences
        return 2 * (seq.count('A') + seq.count('T')) + 4 * (seq.count('G') + seq.count('C'))
    else:
        # Nearest-neighbor approximation for longer sequences
        gc = gc_content(seq)
        return 64.9 + 41 * (gc - 16.4) / 100


def get_basic_stats(sequence: str, motifs: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
    """
    Calculate basic sequence statistics
    
    Args:
        sequence: DNA sequence string
        motifs: Optional list of detected motifs
        
    Returns:
        Dictionary of statistics
    """
    if not sequence:
        return {}
    
    seq = sequence.upper()
    length = len(seq)
    
    # Base composition
    base_counts = Counter(seq)
    
    stats = {
        'Length': length,
        'A': base_counts.get('A', 0),
        'T': base_counts.get('T', 0),
        'G': base_counts.get('G', 0),
        'C': base_counts.get('C', 0),
        'N': base_counts.get('N', 0),
        'GC%': round(gc_content(seq), 2),
        'AT%': round(at_content(seq), 2),
        'Tm': round(calculate_tm(seq), 1)
    }
    
    # Motif statistics if provided
    if motifs:
        stats.update(calculate_motif_statistics(motifs, length))
    
    return stats


def calculate_motif_statistics(motifs: List[Dict[str, Any]], sequence_length: int, 
                                background_motifs: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
    """
    Calculate comprehensive motif statistics with null-model awareness.
    
    Args:
        motifs: List of motif dictionaries (OBSERVED)
        sequence_length: Length of analyzed sequence
        background_motifs: Optional list of motifs from null/background model (EXPECTED)
        
    Returns:
        Dictionary of motif statistics with optional enrichment metrics
    """
    if not motifs:
        base_stats = {
            'Total_Motifs': 0,
            'Coverage%': 0.0,
            'Density': 0.0,
            'Classes_Detected': 0,
            'Subclasses_Detected': 0
        }
        
        # Add null-model structure (empty when no motifs)
        if background_motifs is not None:
            base_stats['Null_Model_Available'] = True
            base_stats['Expected_Motifs'] = len(background_motifs) if background_motifs else 0
            base_stats['Enrichment_Ratio'] = 0.0
        else:
            base_stats['Null_Model_Available'] = False
        
        return base_stats
    
    # Count by class and subclass
    class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    subclass_counts = Counter(m.get('Subclass', 'Unknown') for m in motifs)
    
    # Calculate coverage
    covered_positions = set()
    for motif in motifs:
        start = motif.get('Start', 0) - 1  # Convert 1-based to 0-based
        end = motif.get('End', 0)  # Inclusive end becomes exclusive in range()
        covered_positions.update(range(start, end))
    
    coverage_percent = (len(covered_positions) / sequence_length * 100) if sequence_length > 0 else 0
    density = len(motifs) / (sequence_length / 1000) if sequence_length > 0 else 0  # Motifs per kb
    
    stats = {
        'Total_Motifs': len(motifs),
        'Coverage%': round(coverage_percent, 2),
        'Density': round(density, 2),
        'Classes_Detected': len(class_counts),
        'Subclasses_Detected': len(subclass_counts),
        'Class_Distribution': dict(class_counts),
        'Subclass_Distribution': dict(subclass_counts)
    }
    
    # NULL-MODEL COMPARISON: Add enrichment metrics if background provided
    if background_motifs is not None:
        stats['Null_Model_Available'] = True
        stats['Expected_Motifs'] = len(background_motifs)
        
        # Calculate enrichment ratio (observed / expected)
        if len(background_motifs) > 0:
            stats['Enrichment_Ratio'] = round(len(motifs) / len(background_motifs), 2)
        else:
            stats['Enrichment_Ratio'] = float('inf') if len(motifs) > 0 else 1.0
        
        # Per-class enrichment (future enhancement hook)
        background_class_counts = Counter(m.get('Class', 'Unknown') for m in background_motifs)
        stats['Class_Enrichment'] = {}
        for cls in class_counts:
            obs = class_counts[cls]
            exp = background_class_counts.get(cls, 0)
            stats['Class_Enrichment'][cls] = round(obs / exp, 2) if exp > 0 else float('inf')
    else:
        stats['Null_Model_Available'] = False
    
    # Score statistics
    scores = [m.get('Score', 0) for m in motifs if isinstance(m.get('Score'), (int, float))]
    if scores:
        stats.update({
            'Score_Mean': round(np.mean(scores), 3),
            'Score_Std': round(np.std(scores), 3),
            'Score_Min': round(min(scores), 3),
            'Score_Max': round(max(scores), 3)
        })
    
    # Length statistics
    lengths = [m.get('Length', 0) for m in motifs if isinstance(m.get('Length'), int)]
    if lengths:
        stats.update({
            'Length_Mean': round(np.mean(lengths), 1),
            'Length_Std': round(np.std(lengths), 1),
            'Length_Min': min(lengths),
            'Length_Max': max(lengths)
        })
    
    return stats
