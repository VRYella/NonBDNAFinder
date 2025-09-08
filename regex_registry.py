"""
Central Regex Registry for NBDFinder - All Hyperscan-Safe Patterns
================================================================

This module contains all regex patterns used by NBDFinder motif detection modules,
organized by motif class for shared Hyperscan database compilation and maintenance.

Scientific References:
- G4Hunter: Bedrat et al. NAR 44(4):1746-1759 (2016)
- Z-DNA Seeker: Ho et al. EMBO J 5:2737-2744 (1986); Rich & Zhang Nature 302:209-217 (1983)
- i-motif: Zeraati et al. Nat Chem 10:631-637 (2018)
- Triplex: Frank-Kamenetskii & Mirkin Annu Rev Biochem 64:65-95 (1995)
- Cruciform: Lilley & Clegg Annu Rev Biophys Biomol Struct 22:299-328 (1993)

Pattern Structure:
Each pattern is a tuple: (regex_pattern, pattern_id, group_number, subclass_name, 
                         scoring_function, score_scale, min_runs, min_score, score_method)
"""

import re
import numpy as np

# === G-QUADRUPLEX FAMILY PATTERNS (Class 6) ===
# Based on G4Hunter algorithm and canonical G4 definitions
G_QUADRUPLEX_PATTERNS = {
    'multimeric_g4': [
        (r"(G{3,}\w{1,12}){4,}", 1, 0, "Multimeric_G4", None, 1.2, 4, 0.3, "G4Hunter_Multimer_raw")
    ],
    'bipartite_g4': [
        (r"(G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{10,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,})", 
         2, 1, "Bipartite_G4", None, 0.9, 8, 0.0, "Bipartite_raw")
    ],
    'canonical_g4': [
        (r"(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,})", 3, 1, "Canonical_G4", None, 1.0, 4, 0.5, "G4Hunter_v2_raw")
    ],
    'relaxed_g4': [
        (r"(G{3,}\w{8,12}G{3,}\w{8,12}G{3,}\w{8,12}G{3,})", 4, 1, "Relaxed_G4", None, 0.8, 4, 0.3, "G4Hunter_LongLoop_raw")
    ],
    'bulged_g4': [
        (r"(G{3,}\w{0,3}G{3,}\w{0,3}G{3,}\w{0,3}G{3,})", 5, 1, "Bulged_G4", None, 0.7, 4, 0.0, "G4Hunter_Bulge_raw")
    ],
    'imperfect_g4': [
        (r"(G{2,3}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,})", 6, 1, "Imperfect_G4", None, 1.0, 4, 0.4, "G4Hunter_Imperfect_raw"),
        (r"(G{3,}\w{1,7}G{2,3}\w{1,7}G{3,}\w{1,7}G{3,})", 7, 1, "Imperfect_G4", None, 1.0, 4, 0.4, "G4Hunter_Imperfect_raw"),
        (r"(G{3,}\w{1,7}G{3,}\w{1,7}G{2,3}\w{1,7}G{3,})", 8, 1, "Imperfect_G4", None, 1.0, 4, 0.4, "G4Hunter_Imperfect_raw"),
        (r"(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{2,3})", 9, 1, "Imperfect_G4", None, 1.0, 4, 0.4, "G4Hunter_Imperfect_raw"),
    ],
    'g_triplex': [
        (r"(G{3,}\w{1,7}G{3,}\w{1,7}G{3,})", 10, 1, "G-Triplex_intermediate", None, 1.0, 3, 0.0, "G3_raw")
    ]
}

# === I-MOTIF FAMILY PATTERNS (Class 7) ===
# Based on C-rich quadruplex structures, Zeraati et al. 2018
I_MOTIF_PATTERNS = {
    'canonical_imotif': [
        (r"C{3,}\w{1,12}C{3,}\w{1,12}C{3,}\w{1,12}C{3,}", 1, 0, "Canonical_iMotif", None, 1.0, 4, 0.0, "G4Hunter_adapted")
    ],
    'ac_motif': [
        (r"A{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}", 2, 0, "AC-motif", None, 1.0, 0, 0.0, "G4Hunter_adapted"),
        (r"C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}A{3}", 3, 0, "AC-motif", None, 1.0, 0, 0.0, "G4Hunter_adapted")
    ]
}

# === Z-DNA PATTERNS (Class 8) ===
# Based on Z-DNA seeker algorithm, Ho et al. 1986
Z_DNA_PATTERNS = {
    'z_dna_basic': [
        (r"([CG]{2}){6,}", 1, 0, "Z-DNA", None, 1.0, 0, 0.0, "Z_seeker_raw")
    ],
    'extended_gz': [
        (r"G[CG]{8,}G", 2, 0, "eGZ (Extruded-G) DNA", None, 1.0, 0, 0.0, "Z_seeker_raw")
    ]
}

# === CURVED DNA PATTERNS (Class 1) ===
# Based on A-tract and phased array detection
CURVED_DNA_PATTERNS = {
    'poly_a_tracts': [
        (r"A{7,}", 1, 0, "Global Array", None, 1.0, 0, 0.0, "Curvature_raw"),
    ],
    'poly_t_tracts': [
        (r"T{7,}", 2, 0, "Global Array", None, 1.0, 0, 0.0, "Curvature_raw"),
    ]
}

# === TRIPLEX PATTERNS (Class 5) ===
# Based on homopurine/homopyrimidine tract detection
TRIPLEX_PATTERNS = {
    'homopurine': [
        (r"[AG]{15,}", 1, 0, "Triplex", None, 1.0, 0, 0.0, "Triplex_stability_raw"),
    ],
    'homopyrimidine': [
        (r"[CT]{15,}", 2, 0, "sticky DNA", None, 1.0, 0, 0.0, "Triplex_stability_raw"),
    ]
}

# === CRUCIFORM PATTERNS (Class 3) ===
# Basic palindrome and inverted repeat detection patterns for Hyperscan pre-filtering
CRUCIFORM_PATTERNS = {
    'palindrome_candidates': [],  # Generated dynamically based on arm lengths
    'inverted_repeat_candidates': []  # Generated dynamically based on arm/loop lengths
}

# === R-LOOP PATTERNS (Class 4) ===
# R-loop forming sequences based on RLFS models
R_LOOP_PATTERNS = {
    'rlfs_m1': [
        (r"G{3,}[ATGC]{1,10}?G{3,}(?:[ATGC]{1,10}?G{3,}){1,}", 1, 0, "RLFS_m1", None, 1.0, 0, 0.0, "QmRLFS_raw")
    ],
    'rlfs_m2': [
        (r"G{4,}(?:[ATGC]{1,10}?G{4,}){1,}", 2, 0, "RLFS_m2", None, 1.0, 0, 0.0, "QmRLFS_raw")
    ]
}

# === SLIPPED DNA PATTERNS (Class 2) ===
# Note: Direct repeats and STRs require back-references (Python regex only)
# These patterns are for Hyperscan pre-filtering of repetitive regions
SLIPPED_DNA_PATTERNS = {
    'repetitive_candidates': [
        # Mononucleotide runs
        (r"A{10,}", 1, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"T{10,}", 2, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"G{10,}", 3, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"C{10,}", 4, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        # Dinucleotide repeats
        (r"(AT){8,}", 5, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"(GC){8,}", 6, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"(AG){8,}", 7, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"(CT){8,}", 8, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        # Trinucleotide repeats
        (r"(CGG){6,}", 9, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"(CAG){6,}", 10, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"(AAG){6,}", 11, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
    ]
}

# === MASTER PATTERN REGISTRY ===
ALL_PATTERNS = {
    'g_quadruplex': G_QUADRUPLEX_PATTERNS,
    'i_motif': I_MOTIF_PATTERNS,
    'z_dna': Z_DNA_PATTERNS,
    'curved_dna': CURVED_DNA_PATTERNS,
    'triplex': TRIPLEX_PATTERNS,
    'cruciform': CRUCIFORM_PATTERNS,
    'r_loop': R_LOOP_PATTERNS,
    'slipped_dna': SLIPPED_DNA_PATTERNS,
}

# === UTILITY FUNCTIONS ===

def get_patterns_for_motif(motif_class: str) -> dict:
    """
    Get all patterns for a specific motif class.
    
    Args:
        motif_class: One of 'g_quadruplex', 'i_motif', 'z_dna', 'curved_dna', 
                    'triplex', 'cruciform', 'r_loop', 'slipped_dna'
    
    Returns:
        Dictionary of pattern categories for the motif class
    """
    return ALL_PATTERNS.get(motif_class, {})

def get_all_hyperscan_patterns() -> list:
    """
    Get all patterns suitable for Hyperscan compilation.
    
    Returns:
        List of (pattern, id) tuples for Hyperscan database compilation
    """
    all_patterns = []
    pattern_id = 0
    
    for motif_class, pattern_dict in ALL_PATTERNS.items():
        for pattern_type, patterns in pattern_dict.items():
            for pattern_tuple in patterns:
                if pattern_tuple:  # Skip empty patterns
                    regex_pattern = pattern_tuple[0]
                    all_patterns.append((regex_pattern, pattern_id))
                    pattern_id += 1
    
    return all_patterns

def get_pattern_info(pattern_id: int) -> dict:
    """
    Get detailed information about a pattern by its ID.
    
    Args:
        pattern_id: Unique pattern identifier
        
    Returns:
        Dictionary with pattern details (motif_class, pattern_type, etc.)
    """
    current_id = 0
    
    for motif_class, pattern_dict in ALL_PATTERNS.items():
        for pattern_type, patterns in pattern_dict.items():
            for pattern_tuple in patterns:
                if pattern_tuple and current_id == pattern_id:
                    return {
                        'motif_class': motif_class,
                        'pattern_type': pattern_type,
                        'regex': pattern_tuple[0],
                        'original_id': pattern_tuple[1],
                        'group_number': pattern_tuple[2],
                        'subclass': pattern_tuple[3],
                        'score_scale': pattern_tuple[5] if len(pattern_tuple) > 5 else 1.0,
                        'min_runs': pattern_tuple[6] if len(pattern_tuple) > 6 else 0,
                        'min_score': pattern_tuple[7] if len(pattern_tuple) > 7 else 0.0,
                        'score_method': pattern_tuple[8] if len(pattern_tuple) > 8 else "default"
                    }
                current_id += 1
    
    return {}

def generate_cruciform_patterns(min_arm_len: int = 6, max_arm_len: int = 20, 
                               max_loop_len: int = 50) -> None:
    """
    Generate dynamic cruciform patterns for Hyperscan pre-filtering.
    
    Args:
        min_arm_len: Minimum arm length for palindromes/inverted repeats
        max_arm_len: Maximum arm length
        max_loop_len: Maximum loop length for inverted repeats
    """
    # Clear existing patterns
    CRUCIFORM_PATTERNS['palindrome_candidates'] = []
    CRUCIFORM_PATTERNS['inverted_repeat_candidates'] = []
    
    pattern_id = 1000  # Start with high ID to avoid conflicts
    
    # Generate palindrome patterns (no loop)
    for arm_len in range(min_arm_len, min(max_arm_len + 1, 21)):
        pattern = f'[ATGC]{{{2*arm_len}}}'
        CRUCIFORM_PATTERNS['palindrome_candidates'].append(
            (pattern, pattern_id, 0, "Cruciform DNA [IR]/HairPin [IR]", None, 1.0, 0, 0.0, "NN_Thermodynamics")
        )
        pattern_id += 1
    
    # Generate inverted repeat patterns (with loop)
    for arm_len in range(min_arm_len, min(max_arm_len + 1, 16)):
        for loop_len in range(1, min(21, max_loop_len + 1)):
            total_len = 2 * arm_len + loop_len
            if total_len <= 100:  # Reasonable upper limit
                pattern = f'[ATGC]{{{total_len}}}'
                CRUCIFORM_PATTERNS['inverted_repeat_candidates'].append(
                    (pattern, pattern_id, 0, "Cruciform DNA [IR]/HairPin [IR]", None, 1.0, 0, 0.0, "NN_Thermodynamics")
                )
                pattern_id += 1

def validate_patterns() -> bool:
    """
    Validate all regex patterns for syntax correctness.
    
    Returns:
        True if all patterns are valid, False otherwise
    """
    try:
        for motif_class, pattern_dict in ALL_PATTERNS.items():
            for pattern_type, patterns in pattern_dict.items():
                for pattern_tuple in patterns:
                    if pattern_tuple:
                        regex_pattern = pattern_tuple[0]
                        re.compile(regex_pattern)
        return True
    except re.error as e:
        print(f"Invalid regex pattern found: {e}")
        return False

# Initialize cruciform patterns
generate_cruciform_patterns()

# Validate all patterns on import
if not validate_patterns():
    print("Warning: Some regex patterns in registry are invalid!")

# Export commonly used pattern collections
HYPERSCAN_SAFE_PATTERNS = get_all_hyperscan_patterns()

__all__ = [
    'ALL_PATTERNS',
    'G_QUADRUPLEX_PATTERNS',
    'I_MOTIF_PATTERNS', 
    'Z_DNA_PATTERNS',
    'CURVED_DNA_PATTERNS',
    'TRIPLEX_PATTERNS',
    'CRUCIFORM_PATTERNS',
    'R_LOOP_PATTERNS',
    'SLIPPED_DNA_PATTERNS',
    'get_patterns_for_motif',
    'get_all_hyperscan_patterns',
    'get_pattern_info',
    'generate_cruciform_patterns',
    'validate_patterns',
    'HYPERSCAN_SAFE_PATTERNS'
]