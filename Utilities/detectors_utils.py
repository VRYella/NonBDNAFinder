"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Detectors Utils - Shared utilities for Non-B DNA motif detectors             │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
│ Common functions extracted to reduce code repetition in detectors            │
└──────────────────────────────────────────────────────────────────────────────┘
"""
# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
from typing import List, Dict, Any, Callable, Optional

# ═══════════════════════════════════════════════════════════════════════════════
# TUNABLE PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════════
DEFAULT_UNKNOWN_SUBCLASS = 'unknown'
# Pre-compiled translation table for reverse complement (performance optimization)
_REVCOMP_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")
# Pre-compiled sets for O(1) membership testing
_GC_BASES = {'G', 'C', 'g', 'c'}
_AT_BASES = {'A', 'T', 'a', 't'}
# ═══════════════════════════════════════════════════════════════════════════════


def revcomp(seq: str) -> str:
    """
    Optimized reverse complement of DNA sequence using cached translation table.
    
    Args:
        seq: DNA sequence string
        
    Returns:
        Reverse complement sequence
        
    Example:
        >>> revcomp("ATCG")
        'CGAT'
    """
    return seq.translate(_REVCOMP_TABLE)[::-1]


def calc_gc_content(seq: str) -> float:
    """
    Calculate GC content percentage using genomic standard formula.
    
    Formula: GC% = (G + C) / (A + T + G + C) × 100
    
    Args:
        seq: DNA sequence string (case-insensitive)
        
    Returns:
        GC content as percentage (0-100), or 0.0 if no valid ATGC bases
        
    Scientific Basis:
        - Excludes ambiguous bases (N, R, Y, etc.) from denominator
        - Only counts definitive nucleotides (A, T, G, C)
        - Standard used in genomics (NCBI, Ensembl, UCSC Genome Browser)
        
    Example:
        >>> calc_gc_content("ATCG")
        50.0
        >>> calc_gc_content("ATCGNNNN")  # N's excluded
        50.0  # (2 / 4) * 100, not (2 / 8) * 100
    """
    if not seq:
        return 0.0
    
    # Count each valid base separately using O(1) set membership test
    a_count = sum(1 for c in seq if c in _AT_BASES and c in {'A', 'a'})
    t_count = sum(1 for c in seq if c in _AT_BASES and c in {'T', 't'})
    g_count = sum(1 for c in seq if c in _GC_BASES and c in {'G', 'g'})
    c_count = sum(1 for c in seq if c in _GC_BASES and c in {'C', 'c'})
    
    # Calculate valid base total (denominator)
    valid_bases = a_count + t_count + g_count + c_count
    
    if valid_bases == 0:
        return 0.0
    
    # Gold standard formula: (G+C) / (A+T+G+C) × 100
    return ((g_count + c_count) / valid_bases) * 100


def calc_at_content(seq: str) -> float:
    """
    Calculate AT content percentage using genomic standard formula.
    
    Formula: AT% = (A + T) / (A + T + G + C) × 100
    
    Args:
        seq: DNA sequence string (case-insensitive)
        
    Returns:
        AT content as percentage (0-100), or 0.0 if no valid ATGC bases
        
    Scientific Basis:
        - Excludes ambiguous bases (N, R, Y, etc.) from denominator
        - Only counts definitive nucleotides (A, T, G, C)
        - Consistent with GC% calculation (AT% + GC% = 100%)
        
    Example:
        >>> calc_at_content("ATCG")
        50.0
        >>> calc_at_content("ATCGNNNN")  # N's excluded
        50.0  # (2 / 4) * 100, not (2 / 8) * 100
    """
    if not seq:
        return 0.0
    
    # Count each valid base separately using O(1) set membership test
    a_count = sum(1 for c in seq if c in _AT_BASES and c in {'A', 'a'})
    t_count = sum(1 for c in seq if c in _AT_BASES and c in {'T', 't'})
    g_count = sum(1 for c in seq if c in _GC_BASES and c in {'G', 'g'})
    c_count = sum(1 for c in seq if c in _GC_BASES and c in {'C', 'c'})
    
    # Calculate valid base total (denominator)
    valid_bases = a_count + t_count + g_count + c_count
    
    if valid_bases == 0:
        return 0.0
    
    # Standard formula: (A+T) / (A+T+G+C) × 100
    return ((a_count + t_count) / valid_bases) * 100


def remove_overlaps(motifs: List[Dict[str, Any]], 
                   sort_key: Optional[Callable] = None) -> List[Dict[str, Any]]:
    """
    Remove overlapping motifs, keeping highest priority.
    
    Args:
        motifs: List of motif dictionaries with Start/End positions
        sort_key: Optional sorting function (default: by score desc, then length desc)
        
    Returns:
        List of non-overlapping motifs sorted by Start position
        
    Example:
        >>> motifs = [
        ...     {'Start': 1, 'End': 10, 'Score': 0.8},
        ...     {'Start': 5, 'End': 15, 'Score': 0.9}
        ... ]
        >>> result = remove_overlaps(motifs)
        >>> len(result)
        1
    """
    if not motifs:
        return []
    
    # Default sort: by score (desc), then length (desc)
    if sort_key is None:
        sort_key = lambda x: (-x.get('Score', 0), -(x.get('End', 0) - x.get('Start', 0)))
    
    sorted_motifs = sorted(motifs, key=sort_key)
    selected = []
    
    for motif in sorted_motifs:
        overlaps = False
        for selected_motif in selected:
            if not (motif['End'] <= selected_motif['Start'] or 
                   motif['Start'] >= selected_motif['End']):
                overlaps = True
                break
        
        if not overlaps:
            selected.append(motif)
    
    # Sort by start position for output
    selected.sort(key=lambda x: x['Start'])
    return selected


def remove_overlaps_by_subclass(motifs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Remove overlapping motifs within same subclass only.
    Allows overlaps between different subclasses.
    
    Args:
        motifs: List of motif dictionaries with Start/End/Subclass
        
    Returns:
        List of non-overlapping motifs (within subclass) sorted by Start
    """
    if not motifs:
        return []
    
    from collections import defaultdict
    
    groups = defaultdict(list)
    for motif in motifs:
        subclass = motif.get('Subclass', DEFAULT_UNKNOWN_SUBCLASS)
        groups[subclass].append(motif)
    
    non_overlapping = []
    
    # Process each subclass separately
    for subclass, group_motifs in groups.items():
        # Sort by score (desc), then length (desc)
        sorted_motifs = sorted(
            group_motifs,
            key=lambda x: (-x.get('Score', 0), -(x.get('End', 0) - x.get('Start', 0)))
        )
        
        selected = []
        for motif in sorted_motifs:
            overlaps = False
            for selected_motif in selected:
                if not (motif['End'] <= selected_motif['Start'] or 
                       motif['Start'] >= selected_motif['End']):
                    overlaps = True
                    break
            
            if not overlaps:
                selected.append(motif)
        
        non_overlapping.extend(selected)
    
    # Sort by start position for output
    non_overlapping.sort(key=lambda x: x['Start'])
    return non_overlapping


def load_patterns_with_fallback(patterns: Any, 
                                fallback: Callable[[], Dict]) -> Dict:
    """
    Load patterns with fallback to default patterns.
    
    Args:
        patterns: Imported patterns (may be None or empty dict)
        fallback: Function that returns default patterns
        
    Returns:
        Pattern dictionary (from import or fallback)
        
    Example:
        >>> def fallback():
        ...     return {"pattern1": "ATCG"}
        >>> result = load_patterns_with_fallback(None, fallback)
        >>> "pattern1" in result
        True
    """
    if patterns:
        return patterns.copy() if hasattr(patterns, 'copy') else patterns
    else:
        return fallback()
