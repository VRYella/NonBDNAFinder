"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    MOTIF PATTERN DEFINITIONS MODULE                           ║
║           Centralized Regular Expression Patterns for Motif Detection         ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: motif_patterns.py  
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.1
LICENSE: MIT

DESCRIPTION:
    Centralized repository of all regular expression patterns and definitions
    used for Non-B DNA motif detection. This file allows easy modification
    of patterns without needing to edit the detector classes themselves.
    
    Each pattern is defined as a tuple with the following structure:
    (regex_pattern, pattern_id, name, subclass, min_length, score_type, 
     base_score, description, reference)

PATTERN STRUCTURE:
    - regex_pattern: The regular expression to match
    - pattern_id: Unique identifier for the pattern (e.g., 'CRV_002')
    - name: Human-readable name of the pattern
    - subclass: Motif subclass classification
    - min_length: Minimum length for valid matches
    - score_type: Type of scoring method to use
    - base_score: Base score value for matches
    - description: Description of the motif
    - reference: Literature reference (e.g., 'Olson 1998')

SUPPORTED MOTIF CLASSES:
    1. Curved DNA (CurvedDNADetector)
    2. Z-DNA (ZDNADetector)
    3. A-philic DNA (APhilicDetector)
    4. Slipped DNA (SlippedDNADetector) - STRs and Direct Repeats
    5. Cruciform (CruciformDetector) - Inverted Repeats
    6. R-Loop (RLoopDetector) - QmRLFS patterns
    7. Triplex (TriplexDetector) - Sticky DNA and Mirror Repeats
    8. G-Quadruplex (GQuadruplexDetector) - 7 subclasses
    9. i-Motif (IMotifDetector) - Canonical and HUR AC-motifs

USAGE:
    from motif_patterns import CURVED_DNA_PATTERNS, ZDNA_PATTERNS
    
    # Get patterns for a specific detector
    patterns = CURVED_DNA_PATTERNS['local_curved']
"""

from typing import Dict, List, Tuple

# ============================================================================
# CURVED DNA PATTERNS
# ============================================================================

CURVED_DNA_PATTERNS: Dict[str, List[Tuple]] = {
    'local_curved': [
        # Long A-tract and T-tract patterns causing local DNA curvature
        (r'A{7,}', 'CRV_002', 'Long A-tract', 'Local Curvature', 7, 
         'curvature_score', 0.95, 'A-tract curvature', 'Olson 1998'),
        (r'T{7,}', 'CRV_003', 'Long T-tract', 'Local Curvature', 7, 
         'curvature_score', 0.95, 'T-tract curvature', 'Olson 1998'),
    ],
    # Note: Global curved patterns (A-phased and T-phased repeats) are generated 
    # programmatically in the CurvedDNADetector class
}

# ============================================================================
# Z-DNA PATTERNS
# ============================================================================

ZDNA_PATTERNS: Dict[str, List[Tuple]] = {
    'z_dna_10mers': [
        # Z-DNA detection uses a 10-mer scoring table (algorithmic, not regex)
        (r'', 'ZDN_10MER', 'Z-DNA 10-mer table', 'Z-DNA', 10, 
         'z_dna_10mer_score', 0.9, 'Z-DNA 10mer motif', 'user_table'),
    ],
    'egz_motifs': [
        # Extruded-G Z-DNA (eGZ) patterns - trinucleotide repeats
        (r'(?:CGG){3,}', 'ZDN_EGZ_CGG', 'CGG repeat (eGZ)', 'eGZ', 9, 
         'egz_score', 0.90, 'Extruded-G Z-DNA CGG repeat', 'Herbert 1997'),
        (r'(?:GGC){3,}', 'ZDN_EGZ_GGC', 'GGC repeat (eGZ)', 'eGZ', 9, 
         'egz_score', 0.90, 'Extruded-G Z-DNA GGC repeat', 'Herbert 1997'),
        (r'(?:CCG){3,}', 'ZDN_EGZ_CCG', 'CCG repeat (eGZ)', 'eGZ', 9, 
         'egz_score', 0.90, 'Extruded-G Z-DNA CCG repeat', 'Herbert 1997'),
        (r'(?:GCC){3,}', 'ZDN_EGZ_GCC', 'GCC repeat (eGZ)', 'eGZ', 9, 
         'egz_score', 0.90, 'Extruded-G Z-DNA GCC repeat', 'Herbert 1997'),
    ]
}

# ============================================================================
# A-PHILIC DNA PATTERNS
# ============================================================================

APHILIC_DNA_PATTERNS: Dict[str, List[Tuple]] = {
    'a_philic_10mers': [
        # A-philic DNA detection uses a 10-mer scoring table (algorithmic, not regex)
        (r'', 'APH_10MER', 'A-philic 10-mer table', 'A-philic DNA', 10, 
         'a_philic_10mer_score', 0.9, 'A-philic 10mer motif', 'user_table'),
    ]
}

# ============================================================================
# SLIPPED DNA PATTERNS (STRs and Direct Repeats)
# ============================================================================

SLIPPED_DNA_PATTERNS: Dict[str, List[Tuple]] = {
    'short_tandem_repeats': [
        # STRs are detected algorithmically using k-mer scanning
        # No regex patterns needed - uses optimized scanner
    ],
    'direct_repeats': [
        # Direct repeats (10-50 bp units with 0 bp spacer) detected algorithmically
        # No regex patterns needed - uses optimized scanner
    ]
}

# ============================================================================
# CRUCIFORM PATTERNS (Inverted Repeats / Palindromes)
# ============================================================================

CRUCIFORM_PATTERNS: Dict[str, List[Tuple]] = {
    'inverted_repeats': [
        # Cruciform DNA: 10-100 nt arm length with reverse complement, 0-3 nt spacer
        # Detection is algorithmic using k-mer indexing, not regex-based
        (r'palindrome_like', 'CRU_3_1', 'Potential palindrome', 'Inverted Repeats', 
         12, 'cruciform_stability', 0.95, 'DNA secondary structure', 'Lilley 2000'),
    ]
}

# ============================================================================
# R-LOOP PATTERNS (QmRLFS Models)
# ============================================================================

RLOOP_PATTERNS: Dict[str, List[Tuple]] = {
    'qmrlfs_model_1': [
        # QmRLFS Model 1: 3+ G tracts for R-loop Initiation Zone (RIZ) detection
        (r'G{3,}[ATCGU]{1,10}?G{3,}(?:[ATCGU]{1,10}?G{3,}){1,}?', 
         'QmRLFS_M1', 'QmRLFS Model 1', 'QmRLFS-m1', 25, 'qmrlfs_score', 
         0.90, 'RIZ detection with 3+ G tracts', 'Jenjaroenpun 2016'),
    ],
    'qmrlfs_model_2': [
        # QmRLFS Model 2: 4+ G tracts for higher confidence RIZ detection
        (r'G{4,}(?:[ATCGU]{1,10}?G{4,}){1,}?', 
         'QmRLFS_M2', 'QmRLFS Model 2', 'QmRLFS-m2', 30, 'qmrlfs_score', 
         0.95, 'RIZ detection with 4+ G tracts', 'Jenjaroenpun 2016'),
    ]
}

# ============================================================================
# TRIPLEX PATTERNS (Sticky DNA and Mirror Repeats)
# ============================================================================

TRIPLEX_PATTERNS: Dict[str, List[Tuple]] = {
    'triplex_forming_sequences': [
        # Sticky DNA patterns - disease-associated triplet repeats
        (r'(?:GAA){4,}', 'TRX_5_4', 'GAA repeat', 'Sticky_DNA', 12, 
         'sticky_dna_score', 0.95, 'Disease-associated repeats', 'Sakamoto 1999'),
        (r'(?:TTC){4,}', 'TRX_5_5', 'TTC repeat', 'Sticky_DNA', 12, 
         'sticky_dna_score', 0.95, 'Disease-associated repeats', 'Sakamoto 1999'),
    ],
    # Note: Mirror repeats (10-100 nt arms with 0-8 nt spacer, 90%+ purine/pyrimidine)
    # are detected using optimized k-mer scanner, not regex patterns
}

# ============================================================================
# G-QUADRUPLEX PATTERNS (7 Subclasses)
# ============================================================================

G4_PATTERNS: Dict[str, List[Tuple]] = {
    'canonical_g4': [
        # Classic G4 structure: G3+ N1-7 G3+ N1-7 G3+ N1-7 G3+
        (r'G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}', 
         'G4_0', 'Canonical G4', 'Canonical G4', 15, 'g4hunter_score', 
         0.95, 'Stable G4 structures', 'Burge 2006'),
    ],
    'relaxed_g4': [
        # Relaxed G4: G2+ with longer loops (1-12 nt)
        (r'G{2,}[ACGT]{1,12}G{2,}[ACGT]{1,12}G{2,}[ACGT]{1,12}G{2,}', 
         'G4_1', 'Relaxed G4', 'Relaxed G4', 12, 'g4hunter_score', 
         0.80, 'Potential G4 structures', 'Huppert 2005'),
    ],
    'long_loop_g4': [
        # G4 with extended first loop (8-15 nt)
        (r'G{3,}[ACGT]{8,15}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}', 
         'G4_2', 'Long-loop G4', 'Long-loop G4', 18, 'g4hunter_score', 
         0.75, 'Alternative G4 topology', 'Phan 2006'),
    ],
    'bulged_g4': [
        # G4 with bulge loops (interruptions in G-tracts)
        (r'(?:G{2,3}[ACGT]{0,3}G{1,3})[ACGT]{1,7}G{2,3}[ACGT]{1,7}G{2,3}[ACGT]{1,7}G{2,3}', 
         'G4_3', 'Bulged G4', 'Bulged G4', 20, 'g4hunter_score', 
         0.85, 'G4 with bulge loops', 'Lim 2009'),
    ],
    'multimeric_g4': [
        # Multiple G4 units in tandem (4+ G3+ tracts)
        (r'(?:G{3,}[ACGT]{1,7}){4,}G{3,}', 
         'G4_4', 'Multimeric G4', 'Multimeric G4', 25, 'g4hunter_score', 
         0.90, 'Multiple G4 units', 'Phan 2007'),
    ],
    'imperfect_g4': [
        # G4-like structures with interruptions (imperfect G-tracts)
        (r'G{2,}[ACGT]{1,10}[AG]G{1,3}[ACGT]{1,10}G{2,}[ACGT]{1,10}G{2,}', 
         'G4_5', 'Imperfect G4', 'Imperfect G4', 15, 'g4hunter_score', 
         0.65, 'G4-like with interruptions', 'Kuryavyi 2010'),
    ],
    'g_triplex': [
        # G-triplex: intermediate three-G-tract structures
        (r'G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}', 
         'G4_6', 'G-Triplex', 'G-Triplex intermediate', 12, 'g_triplex_score', 
         0.80, 'Three G-tract structures', 'Sen 1988'),
    ]
}

# ============================================================================
# I-MOTIF PATTERNS (Canonical and HUR AC-motifs)
# ============================================================================

IMOTIF_PATTERNS: Dict[str, List[Tuple]] = {
    'canonical_imotif': [
        # Classic i-motif: C3+ N1-7 C3+ N1-7 C3+ N1-7 C3+ (pH-dependent)
        (r'C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}', 
         'IM_0', 'Canonical i-motif', 'canonical_imotif', 15, 'imotif_score', 
         0.95, 'pH-dependent C-rich structure', 'Gehring 1993'),
    ],
    'hur_ac_motif': [
        # HUR AC-motifs: A-C alternating patterns with different linker lengths
        # Format: A3 N4-6 C3 N4-6 C3 N4-6 C3 (or reverse: C...A)
        
        # 4 bp linker patterns
        (r'A{3}[ACGT]{4}C{3}[ACGT]{4}C{3}[ACGT]{4}C{3}', 
         'HUR_AC_1', 'HUR AC-motif (4bp)', 'AC-motif (HUR)', 18, 'ac_motif_score', 
         0.85, 'HUR AC alternating motif', 'Hur 2021'),
        (r'C{3}[ACGT]{4}C{3}[ACGT]{4}C{3}[ACGT]{4}A{3}', 
         'HUR_AC_2', 'HUR CA-motif (4bp)', 'AC-motif (HUR)', 18, 'ac_motif_score', 
         0.85, 'HUR CA alternating motif', 'Hur 2021'),
        
        # 5 bp linker patterns
        (r'A{3}[ACGT]{5}C{3}[ACGT]{5}C{3}[ACGT]{5}C{3}', 
         'HUR_AC_3', 'HUR AC-motif (5bp)', 'AC-motif (HUR)', 21, 'ac_motif_score', 
         0.85, 'HUR AC alternating motif', 'Hur 2021'),
        (r'C{3}[ACGT]{5}C{3}[ACGT]{5}C{3}[ACGT]{5}A{3}', 
         'HUR_AC_4', 'HUR CA-motif (5bp)', 'AC-motif (HUR)', 21, 'ac_motif_score', 
         0.85, 'HUR CA alternating motif', 'Hur 2021'),
        
        # 6 bp linker patterns
        (r'A{3}[ACGT]{6}C{3}[ACGT]{6}C{3}[ACGT]{6}C{3}', 
         'HUR_AC_5', 'HUR AC-motif (6bp)', 'AC-motif (HUR)', 24, 'ac_motif_score', 
         0.85, 'HUR AC alternating motif', 'Hur 2021'),
        (r'C{3}[ACGT]{6}C{3}[ACGT]{6}C{3}[ACGT]{6}A{3}', 
         'HUR_AC_6', 'HUR CA-motif (6bp)', 'AC-motif (HUR)', 24, 'ac_motif_score', 
         0.85, 'HUR CA alternating motif', 'Hur 2021'),
    ]
}

# ============================================================================
# PATTERN REGISTRY - Maps detector classes to their patterns
# ============================================================================

PATTERN_REGISTRY: Dict[str, Dict[str, List[Tuple]]] = {
    'CurvedDNADetector': CURVED_DNA_PATTERNS,
    'ZDNADetector': ZDNA_PATTERNS,
    'APhilicDetector': APHILIC_DNA_PATTERNS,
    'SlippedDNADetector': SLIPPED_DNA_PATTERNS,
    'CruciformDetector': CRUCIFORM_PATTERNS,
    'RLoopDetector': RLOOP_PATTERNS,
    'TriplexDetector': TRIPLEX_PATTERNS,
    'GQuadruplexDetector': G4_PATTERNS,
    'IMotifDetector': IMOTIF_PATTERNS,
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_patterns_for_detector(detector_name: str) -> Dict[str, List[Tuple]]:
    """
    Get all patterns for a specific detector class.
    
    Args:
        detector_name: Name of the detector class (e.g., 'GQuadruplexDetector')
        
    Returns:
        Dictionary of pattern groups for the detector
        
    Example:
        >>> patterns = get_patterns_for_detector('GQuadruplexDetector')
        >>> canonical_patterns = patterns['canonical_g4']
    """
    return PATTERN_REGISTRY.get(detector_name, {})


def get_all_pattern_ids() -> List[str]:
    """
    Get a list of all pattern IDs across all detectors.
    
    Returns:
        List of all unique pattern IDs
        
    Example:
        >>> ids = get_all_pattern_ids()
        >>> print(ids)
        ['CRV_002', 'CRV_003', 'ZDN_10MER', 'ZDN_EGZ_CGG', ...]
    """
    pattern_ids = []
    for detector_patterns in PATTERN_REGISTRY.values():
        for pattern_group in detector_patterns.values():
            for pattern_tuple in pattern_group:
                if len(pattern_tuple) > 1 and pattern_tuple[1]:  # Has pattern_id
                    pattern_ids.append(pattern_tuple[1])
    return pattern_ids


def get_pattern_by_id(pattern_id: str) -> Tuple:
    """
    Find a pattern tuple by its ID.
    
    Args:
        pattern_id: The pattern ID to search for (e.g., 'G4_0')
        
    Returns:
        Pattern tuple if found, None otherwise
        
    Example:
        >>> pattern = get_pattern_by_id('G4_0')
        >>> print(pattern[2])  # Print pattern name
        'Canonical G4'
    """
    for detector_patterns in PATTERN_REGISTRY.values():
        for pattern_group in detector_patterns.values():
            for pattern_tuple in pattern_group:
                if len(pattern_tuple) > 1 and pattern_tuple[1] == pattern_id:
                    return pattern_tuple
    return None


def list_all_detectors() -> List[str]:
    """
    Get a list of all detector class names.
    
    Returns:
        List of detector class names
        
    Example:
        >>> detectors = list_all_detectors()
        >>> print(detectors)
        ['CurvedDNADetector', 'ZDNADetector', 'APhilicDetector', ...]
    """
    return list(PATTERN_REGISTRY.keys())


def pattern_info_summary() -> str:
    """
    Generate a summary of all patterns in the registry.
    
    Returns:
        Formatted string with pattern statistics
        
    Example:
        >>> print(pattern_info_summary())
        Pattern Registry Summary:
        ========================
        Total Detectors: 9
        Total Pattern Groups: 15
        Total Patterns: 28
        ...
    """
    total_detectors = len(PATTERN_REGISTRY)
    total_groups = sum(len(groups) for groups in PATTERN_REGISTRY.values())
    total_patterns = sum(
        len(patterns) 
        for groups in PATTERN_REGISTRY.values() 
        for patterns in groups.values()
    )
    
    summary = []
    summary.append("Pattern Registry Summary:")
    summary.append("=" * 50)
    summary.append(f"Total Detectors: {total_detectors}")
    summary.append(f"Total Pattern Groups: {total_groups}")
    summary.append(f"Total Patterns: {total_patterns}")
    summary.append("")
    
    for detector_name, pattern_groups in PATTERN_REGISTRY.items():
        num_groups = len(pattern_groups)
        num_patterns = sum(len(patterns) for patterns in pattern_groups.values())
        summary.append(f"{detector_name}: {num_groups} groups, {num_patterns} patterns")
    
    return "\n".join(summary)


# ============================================================================
# MODULE EXPORTS
# ============================================================================

__all__ = [
    # Pattern dictionaries
    'CURVED_DNA_PATTERNS',
    'ZDNA_PATTERNS',
    'APHILIC_DNA_PATTERNS',
    'SLIPPED_DNA_PATTERNS',
    'CRUCIFORM_PATTERNS',
    'RLOOP_PATTERNS',
    'TRIPLEX_PATTERNS',
    'G4_PATTERNS',
    'IMOTIF_PATTERNS',
    'PATTERN_REGISTRY',
    
    # Helper functions
    'get_patterns_for_detector',
    'get_all_pattern_ids',
    'get_pattern_by_id',
    'list_all_detectors',
    'pattern_info_summary',
]
