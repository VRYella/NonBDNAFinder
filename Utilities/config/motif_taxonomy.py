"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    CANONICAL MOTIF TAXONOMY                                  ║
║              Single Source of Truth for Class/Subclass Names                 ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: config/motif_taxonomy.py
AUTHOR: NonBDNAFinder Team
VERSION: 2024.1
LICENSE: MIT

DESCRIPTION:
    This module defines the authoritative motif classification system for
    NonBDNAFinder. ALL components must use these canonical names - no detector,
    visualization, or export module is allowed to emit free-form class or
    subclass strings.
    
    This is the ONLY place where class/subclass names are defined.

SCIENTIFIC ACCURACY:
    - Preserves A-philic DNA (not misrepresented as A-DNA)
    - Distinguishes Canonical vs Relaxed i-motifs  
    - Prevents biologically invalid subclass leakage
    - Makes results publication-safe

USAGE:
    >>> from config.motif_taxonomy import MOTIF_CLASSIFICATION, VALID_CLASSES
    >>> print(VALID_CLASSES)
    >>> print(MOTIF_CLASSIFICATION[6]['subclasses'])  # G-Quadruplex subclasses
"""

from typing import Dict, List, Set, FrozenSet, Any, Optional, Tuple

# =============================================================================
# CANONICAL MOTIF CLASSIFICATION TAXONOMY
# Single Source of Truth - Do NOT duplicate these strings elsewhere
# =============================================================================

MOTIF_CLASSIFICATION: Dict[int, Dict[str, any]] = {
    1: {
        'class': 'Curved_DNA',
        'subclasses': ['Global Curvature', 'Local Curvature']
    },
    2: {
        'class': 'Slipped_DNA',
        'subclasses': ['Direct Repeat', 'STR']
    },
    3: {
        'class': 'Cruciform',
        'subclasses': ['Cruciform forming IRs']
    },
    4: {
        'class': 'R-Loop',
        'subclasses': ['R-loop formation sites']
    },
    5: {
        'class': 'Triplex',
        'subclasses': ['Triplex', 'Sticky DNA']
    },
    6: {
        'class': 'G-Quadruplex',
        'subclasses': [
            'Telomeric G4',
            'Stacked G4',
            'Canonical intramolecular G4',
            'Extended-loop canonical',
            'Higher-order G4 array/G4-wire',
            'Intramolecular G-triplex',
            'Two-tetrad weak PQS',
            'Bulged G4'
        ]
    },
    7: {
        'class': 'i-Motif',
        'subclasses': [
            'Canonical i-motif',
            'Relaxed i-motif',
            'AC-motif'
        ]
    },
    8: {
        'class': 'Z-DNA',
        'subclasses': ['Z-DNA', 'eGZ']
    },
    9: {
        'class': 'A-philic_DNA',
        'subclasses': ['A-philic DNA']
    },
    10: {
        'class': 'Hybrid',
        'subclasses': ['Dynamic overlaps']
    },
    11: {
        'class': 'Non-B_DNA_Clusters',
        'subclasses': ['Dynamic clusters']
    }
}

# =============================================================================
# DERIVED SETS - Automatically generated from MOTIF_CLASSIFICATION
# =============================================================================

# All valid class names (immutable set)
VALID_CLASSES: FrozenSet[str] = frozenset(
    entry['class'] for entry in MOTIF_CLASSIFICATION.values()
)

# All valid subclass names (immutable set)
VALID_SUBCLASSES: FrozenSet[str] = frozenset(
    subclass
    for entry in MOTIF_CLASSIFICATION.values()
    for subclass in entry['subclasses']
)

# Mapping: Class → List of valid subclasses
CLASS_TO_SUBCLASSES: Dict[str, List[str]] = {
    entry['class']: entry['subclasses']
    for entry in MOTIF_CLASSIFICATION.values()
}

# Mapping: Subclass → Parent Class (for reverse lookup)
SUBCLASS_TO_CLASS: Dict[str, str] = {
    subclass: entry['class']
    for entry in MOTIF_CLASSIFICATION.values()
    for subclass in entry['subclasses']
}

# =============================================================================
# LEGACY COMPATIBILITY MAPPINGS
# Maps old/variant names to canonical names
# =============================================================================

# Class name normalizations (case-insensitive, underscore/space variants)
CLASS_ALIASES: Dict[str, str] = {
    # Curved DNA variants
    'curved dna': 'Curved_DNA',
    'curveddna': 'Curved_DNA',
    'curved_dna': 'Curved_DNA',
    'curved-dna': 'Curved_DNA',
    
    # Slipped DNA variants
    'slipped dna': 'Slipped_DNA',
    'slippeddna': 'Slipped_DNA',
    'slipped_dna': 'Slipped_DNA',
    'slipped-dna': 'Slipped_DNA',
    
    # Cruciform variants
    'cruciform': 'Cruciform',
    'cruciform dna': 'Cruciform',
    
    # R-Loop variants
    'r-loop': 'R-Loop',
    'rloop': 'R-Loop',
    'r loop': 'R-Loop',
    
    # Triplex variants
    'triplex': 'Triplex',
    'triplex dna': 'Triplex',
    
    # G-Quadruplex variants
    'g-quadruplex': 'G-Quadruplex',
    'gquadruplex': 'G-Quadruplex',
    'g quadruplex': 'G-Quadruplex',
    'g4': 'G-Quadruplex',
    
    # i-Motif variants
    'i-motif': 'i-Motif',
    'imotif': 'i-Motif',
    'i motif': 'i-Motif',
    
    # Z-DNA variants
    'z-dna': 'Z-DNA',
    'zdna': 'Z-DNA',
    'z dna': 'Z-DNA',
    
    # A-philic variants
    'a-philic_dna': 'A-philic_DNA',
    'a-philic dna': 'A-philic_DNA',
    'aphilic dna': 'A-philic_DNA',
    'aphilic_dna': 'A-philic_DNA',
    
    # Hybrid variants
    'hybrid': 'Hybrid',
    
    # Cluster variants
    'non-b_dna_clusters': 'Non-B_DNA_Clusters',
    'non-b dna clusters': 'Non-B_DNA_Clusters',
    'clusters': 'Non-B_DNA_Clusters',
}

# Subclass name normalizations
SUBCLASS_ALIASES: Dict[str, str] = {
    # Curved DNA subclasses
    'global curvature': 'Global Curvature',
    'local curvature': 'Local Curvature',
    
    # Slipped DNA subclasses
    'direct repeat': 'Direct Repeat',
    'str': 'STR',
    
    # Cruciform subclasses
    'cruciform forming irs': 'Cruciform forming IRs',
    'inverted repeats': 'Cruciform forming IRs',
    'inverted repeat': 'Cruciform forming IRs',
    
    # R-Loop subclasses
    'r-loop formation sites': 'R-loop formation sites',
    'qmrlfs-m1': 'R-loop formation sites',
    'qmrlfs-m2': 'R-loop formation sites',
    
    # Triplex subclasses
    'triplex': 'Triplex',
    'triplex_motif': 'Triplex',  # Legacy alias for backward compatibility
    'triplex motif': 'Triplex',  # Legacy alias for backward compatibility
    'sticky dna': 'Sticky DNA',
    'sticky_dna': 'Sticky DNA',
    
    # G-Quadruplex subclasses
    'telomeric g4': 'Telomeric G4',
    'stacked g4': 'Stacked G4',
    'stacked canonical g4s': 'Stacked G4',
    'stacked g4s with linker': 'Stacked G4',
    'canonical intramolecular g4': 'Canonical intramolecular G4',
    'extended-loop canonical': 'Extended-loop canonical',
    'higher-order g4 array/g4-wire': 'Higher-order G4 array/G4-wire',
    'intramolecular g-triplex': 'Intramolecular G-triplex',
    'two-tetrad weak pqs': 'Two-tetrad weak PQS',
    'bulged g4': 'Bulged G4',
    
    # i-Motif subclasses
    'canonical i-motif': 'Canonical i-motif',
    'canonical_imotif': 'Canonical i-motif',
    'relaxed i-motif': 'Relaxed i-motif',
    'relaxed_imotif': 'Relaxed i-motif',
    'ac-motif': 'AC-motif',
    'ac motif': 'AC-motif',
    
    # Z-DNA subclasses
    'z-dna': 'Z-DNA',
    'egz': 'eGZ',
    'egz-motif': 'eGZ',
    
    # A-philic subclasses
    'a-philic dna': 'A-philic DNA',
    
    # Hybrid/Cluster subclasses
    'dynamic overlaps': 'Dynamic overlaps',
    'dynamic clusters': 'Dynamic clusters',
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def get_all_classes() -> List[str]:
    """
    Get list of all canonical class names in alphabetical order.
    
    Note: For visualizations, use get_all_classes_taxonomy_order() instead
    to maintain biological taxonomy order.
    
    Returns:
        Alphabetically sorted list of class names
    """
    return sorted(VALID_CLASSES)


def get_all_classes_taxonomy_order() -> List[str]:
    """
    Get list of all canonical class names in taxonomy order.
    
    Returns taxonomy order as defined in MOTIF_CLASSIFICATION (by ID).
    This should be used for visualizations instead of alphabetical sorting.
    
    Returns:
        List of class names in taxonomy order
    """
    return [MOTIF_CLASSIFICATION[i]['class'] for i in sorted(MOTIF_CLASSIFICATION.keys())]


def get_all_subclasses_taxonomy_order() -> List[str]:
    """
    Get list of all canonical subclass names in taxonomy order.
    
    Returns subclasses in the order they appear in MOTIF_CLASSIFICATION.
    Within each class, subclasses are in the order defined in the taxonomy.
    
    Returns:
        List of subclass names in taxonomy order
    """
    subclasses = []
    for i in sorted(MOTIF_CLASSIFICATION.keys()):
        subclasses.extend(MOTIF_CLASSIFICATION[i]['subclasses'])
    return subclasses


def get_subclasses_for_class(class_name: str) -> List[str]:
    """
    Get list of valid subclasses for a given class.
    
    Args:
        class_name: Canonical class name
        
    Returns:
        List of subclass names for this class
        
    Raises:
        KeyError: If class_name is not valid
    """
    if class_name not in VALID_CLASSES:
        raise KeyError(f"Invalid class name: {class_name}")
    return CLASS_TO_SUBCLASSES[class_name]


def get_class_for_subclass(subclass_name: str) -> str:
    """
    Get parent class for a given subclass.
    
    Args:
        subclass_name: Canonical subclass name
        
    Returns:
        Parent class name
        
    Raises:
        KeyError: If subclass_name is not valid
    """
    if subclass_name not in VALID_SUBCLASSES:
        raise KeyError(f"Invalid subclass name: {subclass_name}")
    return SUBCLASS_TO_CLASS[subclass_name]


def is_valid_class(class_name: str) -> bool:
    """
    Check if a class name is valid.
    
    Args:
        class_name: Class name to check
        
    Returns:
        True if valid, False otherwise
    """
    return class_name in VALID_CLASSES


def is_valid_subclass(subclass_name: str) -> bool:
    """
    Check if a subclass name is valid.
    
    Args:
        subclass_name: Subclass name to check
        
    Returns:
        True if valid, False otherwise
    """
    return subclass_name in VALID_SUBCLASSES


def is_valid_pairing(class_name: str, subclass_name: str) -> bool:
    """
    Check if a class/subclass pairing is valid.
    
    Args:
        class_name: Class name
        subclass_name: Subclass name
        
    Returns:
        True if pairing is valid, False otherwise
    """
    if class_name not in VALID_CLASSES:
        return False
    if subclass_name not in VALID_SUBCLASSES:
        return False
    return subclass_name in CLASS_TO_SUBCLASSES[class_name]


def build_motif_selector_data(enabled_subclasses: Optional[Set[str]] = None) -> List[Dict[str, Any]]:
    """
    Build motif selector table data for st.data_editor.
    
    Creates a flat list of rows where each row represents a submotif,
    grouped by their parent motif class. All rows are enabled by default.
    
    Args:
        enabled_subclasses: Optional set of subclass names to enable.
                           If None, all subclasses are enabled by default.
    
    Returns:
        List of dicts with keys: 'Enabled', 'Motif Class', 'Submotif'
        
    Example:
        >>> data = build_motif_selector_data()
        >>> df = pd.DataFrame(data)
        >>> edited_df = st.data_editor(df, ...)
    """
    rows = []
    
    # Sort by class ID for consistent ordering
    for class_id in sorted(MOTIF_CLASSIFICATION.keys()):
        entry = MOTIF_CLASSIFICATION[class_id]
        motif_class = entry['class']
        subclasses = entry['subclasses']
        
        for subclass in subclasses:
            # Default: all enabled, or check if in enabled_subclasses set
            if enabled_subclasses is None:
                is_enabled = True
            else:
                is_enabled = subclass in enabled_subclasses
            
            rows.append({
                'Enabled': is_enabled,
                'Motif Class': motif_class,
                'Submotif': subclass
            })
    
    return rows


def get_enabled_from_selector_data(selector_data: List[Dict[str, Any]]) -> Tuple[List[str], List[str]]:
    """
    Extract enabled classes and subclasses from selector table data.
    
    Args:
        selector_data: List of dicts from st.data_editor with 'Enabled', 
                      'Motif Class', and 'Submotif' keys
    
    Returns:
        Tuple of (enabled_classes: List[str], enabled_subclasses: List[str])
    """
    enabled_classes = set()
    enabled_subclasses = []
    
    for row in selector_data:
        if row.get('Enabled', False):
            motif_class = row.get('Motif Class', '')
            submotif = row.get('Submotif', '')
            # Only add non-empty values
            if motif_class:
                enabled_classes.add(motif_class)
            if submotif:
                enabled_subclasses.append(submotif)
    
    return list(enabled_classes), enabled_subclasses


# =============================================================================
# MODULE METADATA
# =============================================================================

__all__ = [
    'MOTIF_CLASSIFICATION',
    'VALID_CLASSES',
    'VALID_SUBCLASSES',
    'CLASS_TO_SUBCLASSES',
    'SUBCLASS_TO_CLASS',
    'CLASS_ALIASES',
    'SUBCLASS_ALIASES',
    'get_all_classes',
    'get_all_classes_taxonomy_order',
    'get_all_subclasses_taxonomy_order',
    'get_subclasses_for_class',
    'get_class_for_subclass',
    'is_valid_class',
    'is_valid_subclass',
    'is_valid_pairing',
    'build_motif_selector_data',
    'get_enabled_from_selector_data',
]
