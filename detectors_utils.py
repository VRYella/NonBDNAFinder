"""
Shared utilities for Non-B DNA motif detectors.
Extracted common functions to reduce code repetition in detectors.py

Author: Dr. Venkata Rajesh Yella
Version: 2024.1
"""

from typing import List, Dict, Any, Callable, Optional


def revcomp(seq: str) -> str:
    """
    Reverse complement of DNA sequence.
    
    Args:
        seq: DNA sequence string
        
    Returns:
        Reverse complement sequence
        
    Example:
        >>> revcomp("ATCG")
        'CGAT'
    """
    trans = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(trans)[::-1]


def calc_gc_content(seq: str) -> float:
    """
    Calculate GC content percentage.
    
    Args:
        seq: DNA sequence string
        
    Returns:
        GC content as percentage (0-100)
        
    Example:
        >>> calc_gc_content("ATCG")
        50.0
    """
    if len(seq) == 0:
        return 0.0
    return (seq.count('G') + seq.count('C')) / len(seq) * 100


def calc_at_content(seq: str) -> float:
    """
    Calculate AT content percentage.
    
    Args:
        seq: DNA sequence string
        
    Returns:
        AT content as percentage (0-100)
    """
    if len(seq) == 0:
        return 0.0
    return (seq.count('A') + seq.count('T')) / len(seq) * 100


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
    
    # Group by subclass
    groups = defaultdict(list)
    for motif in motifs:
        subclass = motif.get('Subclass', 'unknown')
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
