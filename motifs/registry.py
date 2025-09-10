"""
Pattern registry for motif detection - interfaces with existing regex_registry.py
"""

from typing import List, Tuple, Dict, Any
import re
import sys
import os

# Add parent directory to path to import existing registry
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from regex_registry import ALL_PATTERNS
except ImportError:
    print("Warning: Could not import regex_registry, using empty patterns")
    ALL_PATTERNS = {}


def get_patterns_for_motif(motif_class: str) -> Dict[str, List[Tuple]]:
    """
    Get patterns for a specific motif class
    
    Args:
        motif_class: Name of motif class (e.g., 'g_quadruplex', 'triplex')
    
    Returns:
        Dictionary of subclass -> list of patterns
    """
    return ALL_PATTERNS.get(motif_class, {})


def get_all_hyperscan_patterns() -> List[Tuple[str, int, str, str]]:
    """
    Get all Hyperscan-safe patterns from the registry
    
    Returns:
        List of tuples: (pattern, pattern_id, class_name, subclass)
    """
    all_patterns = []
    pattern_id = 0
    
    for class_name, class_patterns in ALL_PATTERNS.items():
        for subclass, patterns in class_patterns.items():
            for pattern_info in patterns:
                pattern, orig_id, group_num, subclass_name = pattern_info[:4]
                
                # Skip patterns with backreferences (not Hyperscan compatible)
                if has_backreference(pattern):
                    print(f"Skipping pattern with backreference: {pattern}")
                    continue
                
                all_patterns.append((pattern, pattern_id, class_name, subclass_name))
                pattern_id += 1
    
    return all_patterns


def has_backreference(pattern: str) -> bool:
    """
    Check if a regex pattern contains backreferences (not Hyperscan compatible)
    
    Args:
        pattern: Regex pattern string
    
    Returns:
        True if pattern contains backreferences
    """
    # Check for numbered backreferences like \1, \2, etc.
    if re.search(r'\\[1-9]', pattern):
        return True
    
    # Check for named backreferences like (?P=name)
    if '(?P=' in pattern:
        return True
    
    # Check for other backreference patterns
    if re.search(r'\(\?\#', pattern):  # Comments
        return True
    
    return False


def get_hyperscan_safe_patterns() -> Dict[str, List[str]]:
    """
    Get patterns grouped by class, filtered for Hyperscan compatibility
    
    Returns:
        Dictionary of class_name -> list of pattern strings
    """
    patterns_by_class = {}
    
    for class_name, class_patterns in ALL_PATTERNS.items():
        safe_patterns = []
        for subclass, patterns in class_patterns.items():
            for pattern_info in patterns:
                pattern = pattern_info[0]
                if not has_backreference(pattern):
                    safe_patterns.append(pattern)
        
        if safe_patterns:
            patterns_by_class[class_name] = safe_patterns
    
    return patterns_by_class


def get_patterns_requiring_fallback() -> Dict[str, List[str]]:
    """
    Get patterns that require Python regex fallback (contain backreferences)
    
    Returns:
        Dictionary of class_name -> list of unsafe pattern strings
    """
    unsafe_patterns = {}
    
    for class_name, class_patterns in ALL_PATTERNS.items():
        fallback_patterns = []
        for subclass, patterns in class_patterns.items():
            for pattern_info in patterns:
                pattern = pattern_info[0]
                if has_backreference(pattern):
                    fallback_patterns.append(pattern)
        
        if fallback_patterns:
            unsafe_patterns[class_name] = fallback_patterns
    
    return unsafe_patterns


# Motif class mappings
MOTIF_CLASS_IDS = {
    'curved_dna': 1,
    'slipped_dna': 2,
    'cruciform': 3,
    'r_loop': 4,
    'triplex': 5,
    'g_quadruplex': 6,
    'i_motif': 7,
    'z_dna': 8,
    'a_philic': 9,
    'hybrid': 10,
    'cluster': 11
}

CLASS_NAMES = {v: k for k, v in MOTIF_CLASS_IDS.items()}


def get_class_id(class_name: str) -> int:
    """Get numeric class ID for a motif class name"""
    return MOTIF_CLASS_IDS.get(class_name, 0)


def get_class_name(class_id: int) -> str:
    """Get class name for a numeric class ID"""
    return CLASS_NAMES.get(class_id, 'unknown')