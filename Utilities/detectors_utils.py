"""Shared utility functions for Non-B DNA motif detectors."""
# IMPORTS
from typing import List, Dict, Any, Callable, Optional

# TUNABLE PARAMETERS
DEFAULT_UNKNOWN_SUBCLASS = 'unknown'
_REVCOMP_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_REVCOMP_TABLE)[::-1]


def _count_bases(seq: str) -> tuple:
    """Count A, T, G, C bases; returns (a, t, g, c) counts."""
    a_count = seq.count('A') + seq.count('a')
    t_count = seq.count('T') + seq.count('t')
    g_count = seq.count('G') + seq.count('g')
    c_count = seq.count('C') + seq.count('c')
    return a_count, t_count, g_count, c_count


def calc_gc_content(seq: str) -> float:
    """Return GC content as percentage (0-100), excluding ambiguous bases."""
    if not seq:
        return 0.0
    
    a_count, t_count, g_count, c_count = _count_bases(seq)
    
    valid_bases = a_count + t_count + g_count + c_count
    
    if valid_bases == 0:
        return 0.0
    
    return ((g_count + c_count) / valid_bases) * 100


def calc_at_content(seq: str) -> float:
    """Return AT content as percentage (0-100), excluding ambiguous bases."""
    if not seq:
        return 0.0
    
    a_count, t_count, g_count, c_count = _count_bases(seq)
    
    valid_bases = a_count + t_count + g_count + c_count
    
    if valid_bases == 0:
        return 0.0
    
    return ((a_count + t_count) / valid_bases) * 100


def remove_overlaps(motifs: List[Dict[str, Any]], 
                   sort_key: Optional[Callable] = None) -> List[Dict[str, Any]]:
    """Remove overlapping motifs, keeping highest-score first."""
    if not motifs:
        return []
    
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
    
    selected.sort(key=lambda x: x['Start'])
    return selected


def remove_overlaps_by_subclass(motifs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Remove overlapping motifs within each subclass independently."""
    if not motifs:
        return []
    
    from collections import defaultdict
    
    groups = defaultdict(list)
    for motif in motifs:
        subclass = motif.get('Subclass', DEFAULT_UNKNOWN_SUBCLASS)
        groups[subclass].append(motif)
    
    non_overlapping = []
    
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
    
    non_overlapping.sort(key=lambda x: x['Start'])
    return non_overlapping


def load_patterns_with_fallback(patterns: Any, 
                                fallback: Callable[[], Dict]) -> Dict:
    """Return patterns dict; uses fallback() if patterns is empty/None."""
    if patterns:
        return patterns.copy() if hasattr(patterns, 'copy') else patterns
    else:
        return fallback()
