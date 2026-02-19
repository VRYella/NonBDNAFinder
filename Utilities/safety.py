"""
Safe list/array access utilities to prevent IndexError crashes.

Provides bounds-checked access to lists, arrays, and DataFrames with
automatic logging and graceful fallbacks.
"""

import logging
from typing import List, Any, Optional, Union, Dict
import pandas as pd
import re

logger = logging.getLogger(__name__)


def safe_list_access(lst: List, index: int, default=None, log_error: bool = True) -> Any:
    """
    Safely access list element with bounds checking.
    
    Args:
        lst: List to access
        index: Index to retrieve
        default: Value to return if index out of bounds
        log_error: Whether to log warning on bounds error
        
    Returns:
        lst[index] if valid, else default
        
    Example:
        >>> safe_list_access([1, 2, 3], 5, default=0)
        0  # Instead of IndexError
    """
    if not isinstance(lst, (list, tuple)):
        if log_error:
            logger.warning(f"safe_list_access called on non-list type: {type(lst)}")
        return default
    
    if 0 <= index < len(lst):
        return lst[index]
    
    if log_error:
        logger.warning(f"Index {index} out of bounds for list of length {len(lst)}")
    
    return default


def safe_min_max(lst: List, operation: str = 'min', default: float = 0.0) -> float:
    """
    Safely get min/max of list, handling empty case.
    
    Args:
        lst: List of numeric values
        operation: 'min' or 'max'
        default: Value to return if list is empty
        
    Returns:
        min/max of list, or default if empty
        
    Example:
        >>> safe_min_max([], 'min', default=0)
        0  # Instead of ValueError: min() arg is an empty sequence
    """
    if not lst:
        logger.warning(f"Attempted {operation} on empty list, returning default={default}")
        return default
    
    try:
        return min(lst) if operation == 'min' else max(lst)
    except (ValueError, TypeError) as e:
        logger.error(f"Error in {operation} operation: {e}")
        return default


def safe_dataframe_access(df: pd.DataFrame, index: int, default=None) -> Optional[pd.Series]:
    """
    Safely access DataFrame row with bounds checking.
    
    Args:
        df: DataFrame to access
        index: Row index to retrieve
        default: Value to return if index out of bounds
        
    Returns:
        DataFrame row or default
    """
    if not isinstance(df, pd.DataFrame):
        logger.warning(f"safe_dataframe_access called on non-DataFrame: {type(df)}")
        return default
    
    if 0 <= index < len(df):
        return df.iloc[index]
    
    logger.warning(f"DataFrame index {index} out of bounds (len={len(df)})")
    return default


def filter_valid_indices(indices: List[int], max_length: int, 
                         log_filtered: bool = True) -> List[int]:
    """
    Filter list of indices to only valid values [0, max_length).
    
    Args:
        indices: List of indices to filter
        max_length: Maximum valid index (exclusive)
        log_filtered: Whether to log filtered indices
        
    Returns:
        List of valid indices only
        
    Example:
        >>> filter_valid_indices([0, 5, 10, 15], max_length=10)
        [0, 5]  # 10 and 15 filtered out
    """
    valid = [i for i in indices if 0 <= i < max_length]
    
    if log_filtered and len(valid) < len(indices):
        filtered_count = len(indices) - len(valid)
        logger.warning(f"Filtered {filtered_count} out-of-bounds indices "
                      f"(valid range: 0-{max_length-1})")
    
    return valid


def validate_sequence_input(sequences: List[str], min_length: int = 10) -> tuple:
    """
    Validate and filter sequences before analysis.
    
    Args:
        sequences: List of DNA sequences
        min_length: Minimum sequence length to accept
        
    Returns:
        Tuple of (valid_sequences, issues_list)
    """
    issues = []
    valid_sequences = []
    
    for i, seq in enumerate(sequences):
        # Check empty
        if not seq or len(seq) == 0:
            issues.append(f"Sequence {i+1} is empty - skipping")
            continue
        
        # Check too short
        if len(seq) < min_length:
            issues.append(f"Sequence {i+1} is very short ({len(seq)}bp < {min_length}bp) - may find no motifs")
            # Still include, just warn
        
        # Check invalid characters
        if not re.match(r'^[ATCGNatcgn]+$', seq):
            issues.append(f"Sequence {i+1} contains invalid characters - skipping")
            continue
        
        valid_sequences.append(seq)
    
    if issues:
        for issue in issues:
            logger.warning(issue)
    
    return valid_sequences, issues


def generate_sequence_quality_report(sequences: List[str]) -> Dict[str, Any]:
    """
    Generate comprehensive quality report for input sequences.
    
    Uses the gold-standard preprocessing pipeline to analyze each sequence
    and provide detailed statistics.
    
    Args:
        sequences: List of DNA sequences to analyze
        
    Returns:
        Dictionary with comprehensive quality metrics:
        {
            'total_sequences': int,
            'total_length': int,
            'average_gc': float,
            'sequences': [
                {
                    'name': str,
                    'length': int,
                    'valid_bases': int,
                    'n_count': int,
                    'gc_percentage': float,
                    'at_percentage': float,
                    'invalid_chars': Dict[str, int],
                    'status': str  # "valid" | "warning" | "error"
                },
                ...
            ],
            'global_warnings': List[str],
            'global_errors': List[str]
        }
        
    Example:
        >>> sequences = ["ATCGATCG", "ATCGNNNN", "ATCG123"]
        >>> report = generate_sequence_quality_report(sequences)
        >>> report['total_sequences']
        3
        >>> len(report['global_errors'])
        1  # One sequence has invalid characters
    """
    from Utilities.sequence_preprocessor import preprocess_sequence
    
    report = {
        'total_sequences': len(sequences),
        'total_length': 0,
        'average_gc': 0.0,
        'sequences': [],
        'global_warnings': [],
        'global_errors': []
    }
    
    gc_values = []
    
    for i, seq in enumerate(sequences):
        seq_name = f"sequence_{i+1}"
        
        try:
            result = preprocess_sequence(seq)
            
            # Build sequence report
            seq_report = {
                'name': seq_name,
                'length': result.length,
                'valid_bases': result.valid_bases,
                'n_count': result.character_counts.get('N', 0),
                'gc_percentage': result.gc_percentage,
                'at_percentage': result.at_percentage,
                'invalid_chars': {char: len(positions) for char, positions in result.invalid_characters.items()},
                'status': result.validation_status
            }
            
            report['sequences'].append(seq_report)
            report['total_length'] += result.length
            
            # Collect GC values for average (only from valid sequences)
            if result.valid_bases > 0:
                gc_values.append(result.gc_percentage)
            
            # Collect global warnings and errors
            for warning in result.warnings:
                report['global_warnings'].append(f"{seq_name}: {warning}")
            
            for error in result.errors:
                report['global_errors'].append(f"{seq_name}: {error}")
                
        except Exception as e:
            # Handle any preprocessing errors
            logger.error(f"Error processing {seq_name}: {e}")
            report['global_errors'].append(f"{seq_name}: Preprocessing exception - {str(e)}")
    
    # Calculate average GC content
    if gc_values:
        report['average_gc'] = sum(gc_values) / len(gc_values)
    
    return report
