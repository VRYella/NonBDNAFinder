"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    EXPORT VALIDATION MODULE                                  ║
║              Hard Validation Before CSV/JSON/BED Export                      ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: export/export_validator.py
AUTHOR: NonBDNAFinder Team
VERSION: 2024.1
LICENSE: MIT

DESCRIPTION:
    Validates motif data before export to ensure all Class/Subclass fields are
    canonical and valid. Prevents corrupted data from being exported.
    
    This module enforces that:
    - All motifs have Class and Subclass fields
    - All Class names are in VALID_CLASSES
    - All Subclass names are in VALID_SUBCLASSES
    - All Class/Subclass pairings are biologically valid
    
USAGE:
    >>> from export.export_validator import validate_export_data
    >>> motifs = [{'Class': 'G-Quadruplex', 'Subclass': 'Telomeric G4', ...}]
    >>> validate_export_data(motifs)  # Raises if invalid
    >>> safe_to_export = True
"""

from typing import List, Dict, Any, Tuple
import warnings

from Utilities.config.motif_taxonomy import VALID_CLASSES, VALID_SUBCLASSES, is_valid_pairing
from Utilities.core.motif_normalizer import (
    MotifNormalizationError,
    normalize_motif_dict,
    validate_motif_dict
)


class ExportValidationError(ValueError):
    """Raised when export data fails validation"""
    pass


def validate_single_motif(
    motif: Dict[str, Any],
    index: int,
    strict: bool = True
) -> Tuple[bool, str]:
    """
    Validate a single motif for export.
    
    Args:
        motif: Motif dictionary to validate
        index: Index of motif in list (for error messages)
        strict: If True, fail on any warning. If False, only fail on errors.
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    # Check for required fields
    if 'Class' not in motif:
        return False, f"Motif at index {index} missing required 'Class' field"
    
    if 'Subclass' not in motif:
        return False, f"Motif at index {index} missing required 'Subclass' field"
    
    class_name = motif['Class']
    subclass_name = motif['Subclass']
    
    # Validate class
    if class_name not in VALID_CLASSES:
        return False, (
            f"Motif at index {index} has invalid class '{class_name}'. "
            f"Valid classes: {sorted(VALID_CLASSES)}"
        )
    
    # Validate subclass
    if subclass_name not in VALID_SUBCLASSES:
        return False, (
            f"Motif at index {index} has invalid subclass '{subclass_name}'. "
            f"Valid subclasses: {sorted(VALID_SUBCLASSES)}"
        )
    
    # Validate pairing
    if not is_valid_pairing(class_name, subclass_name):
        return False, (
            f"Motif at index {index} has invalid class/subclass pairing: "
            f"'{class_name}' / '{subclass_name}'"
        )
    
    return True, ""


def validate_export_data(
    motifs: List[Dict[str, Any]],
    auto_normalize: bool = True,
    strict: bool = True
) -> List[Dict[str, Any]]:
    """
    Validate list of motifs before export.
    
    This function must be called before any export operation (CSV, JSON, BED, Excel).
    It ensures all motifs have valid canonical Class/Subclass fields.
    
    Args:
        motifs: List of motif dictionaries to validate
        auto_normalize: If True, attempt to normalize invalid motifs before failing
        strict: If True, fail on any validation error. If False, warn and skip invalid motifs.
        
    Returns:
        List of validated (and possibly normalized) motifs
        
    Raises:
        ExportValidationError: If validation fails and strict=True
        
    Examples:
        >>> motifs = [{'Class': 'G-Quadruplex', 'Subclass': 'Telomeric G4', ...}]
        >>> validated = validate_export_data(motifs)
        >>> # Now safe to export
    """
    if not motifs:
        warnings.warn("Empty motif list provided for export validation")
        return []
    
    validated_motifs = []
    errors = []
    
    for i, motif in enumerate(motifs):
        # Try auto-normalization first if enabled
        if auto_normalize:
            try:
                normalized = normalize_motif_dict(motif, strict=False, auto_correct=True)
                motif = normalized
            except Exception as e:
                if strict:
                    errors.append(f"Failed to normalize motif at index {i}: {e}")
                    continue
                else:
                    warnings.warn(f"Failed to normalize motif at index {i}: {e}")
        
        # Validate
        is_valid, error_msg = validate_single_motif(motif, i, strict=strict)
        
        if is_valid:
            validated_motifs.append(motif)
        else:
            if strict:
                errors.append(error_msg)
            else:
                warnings.warn(error_msg)
    
    # If any errors in strict mode, fail
    if errors and strict:
        error_summary = "\n".join(errors)
        raise ExportValidationError(
            f"Export validation failed with {len(errors)} error(s):\n{error_summary}"
        )
    
    # Warn if we filtered out motifs in non-strict mode
    if not strict and len(validated_motifs) < len(motifs):
        filtered_count = len(motifs) - len(validated_motifs)
        warnings.warn(
            f"Filtered out {filtered_count} invalid motif(s) during export validation"
        )
    
    return validated_motifs


def get_export_summary(motifs: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Get summary statistics for export data.
    
    Args:
        motifs: List of motif dictionaries
        
    Returns:
        Dictionary with export statistics:
            - total_motifs: Total number of motifs
            - classes: Set of classes present
            - subclasses: Set of subclasses present
            - class_counts: Count of motifs per class
            - validation_status: Whether all motifs are valid
    """
    from collections import Counter
    
    total = len(motifs)
    classes = set()
    subclasses = set()
    class_counts = Counter()
    all_valid = True
    
    for motif in motifs:
        if 'Class' in motif:
            classes.add(motif['Class'])
            class_counts[motif['Class']] += 1
        
        if 'Subclass' in motif:
            subclasses.add(motif['Subclass'])
        
        # Check validity
        is_valid, _ = validate_single_motif(motif, 0, strict=False)
        if not is_valid:
            all_valid = False
    
    return {
        'total_motifs': total,
        'classes': sorted(classes),
        'subclasses': sorted(subclasses),
        'class_counts': dict(class_counts),
        'validation_status': 'VALID' if all_valid else 'INVALID',
        'unique_classes': len(classes),
        'unique_subclasses': len(subclasses)
    }


def check_class_completeness(motifs: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Check which canonical classes are present in the data.
    
    Useful for verifying that all 11 expected classes are detected.
    
    Args:
        motifs: List of motif dictionaries
        
    Returns:
        Dictionary with completeness information:
            - present_classes: Classes found in data
            - missing_classes: Classes not found in data
            - completeness_ratio: Fraction of classes present
    """
    present = set(m.get('Class') for m in motifs if m.get('Class'))
    present_canonical = present & VALID_CLASSES
    missing = VALID_CLASSES - present_canonical
    
    return {
        'present_classes': sorted(present_canonical),
        'missing_classes': sorted(missing),
        'completeness_ratio': len(present_canonical) / len(VALID_CLASSES),
        'total_classes_expected': len(VALID_CLASSES),
        'total_classes_found': len(present_canonical)
    }


# =============================================================================
# MODULE METADATA
# =============================================================================

__all__ = [
    'ExportValidationError',
    'validate_single_motif',
    'validate_export_data',
    'get_export_summary',
    'check_class_completeness',
]
