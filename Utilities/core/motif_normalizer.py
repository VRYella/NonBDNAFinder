"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    MOTIF NORMALIZATION & ENFORCEMENT                         ║
║           Mandatory Validation Layer for Class/Subclass Pairing              ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: core/motif_normalizer.py
AUTHOR: NonBDNAFinder Team
VERSION: 2024.1
LICENSE: MIT

DESCRIPTION:
    Enforces canonical class/subclass pairing across all detectors and exports.
    All motif data must pass through normalization before being returned to users.
    
    This module:
    - Normalizes casing and formatting
    - Blocks unknown subclasses
    - Ensures subclass belongs to class
    - Collapses internal detector labels into canonical names
    
RESPONSIBILITIES:
    1. Validate class/subclass combinations
    2. Apply legacy/alias mappings
    3. Fail fast on invalid data
    4. Provide helpful error messages

USAGE:
    >>> from core.motif_normalizer import normalize_class_subclass
    >>> class_name, subclass = normalize_class_subclass('g-quadruplex', 'Telomeric G4')
    >>> print(class_name, subclass)
    'G-Quadruplex' 'Telomeric G4'
"""

from typing import Tuple, Dict, Any, Optional
import warnings

from Utilities.config.motif_taxonomy import (
    VALID_CLASSES,
    VALID_SUBCLASSES,
    CLASS_TO_SUBCLASSES,
    SUBCLASS_TO_CLASS,
    CLASS_ALIASES,
    SUBCLASS_ALIASES,
    is_valid_pairing
)


class MotifNormalizationError(ValueError):
    """Raised when motif class/subclass normalization fails"""
    pass


def normalize_class_name(class_name: str, strict: bool = True) -> str:
    """
    Normalize a class name to canonical form.
    
    Args:
        class_name: Raw class name (may have variant spelling/casing)
        strict: If True, raise error on unknown class. If False, return as-is with warning.
        
    Returns:
        Canonical class name
        
    Raises:
        MotifNormalizationError: If class_name is invalid and strict=True
        
    Examples:
        >>> normalize_class_name('g-quadruplex')
        'G-Quadruplex'
        >>> normalize_class_name('curved dna')
        'Curved_DNA'
    """
    if not class_name:
        if strict:
            raise MotifNormalizationError("Class name cannot be empty")
        return "Unknown"
    
    # Check if already canonical
    if class_name in VALID_CLASSES:
        return class_name
    
    # Try case-insensitive lookup in aliases
    class_lower = class_name.lower().strip()
    if class_lower in CLASS_ALIASES:
        return CLASS_ALIASES[class_lower]
    
    # Try exact match in VALID_CLASSES (case-insensitive)
    for valid_class in VALID_CLASSES:
        if valid_class.lower() == class_lower:
            return valid_class
    
    # Not found
    if strict:
        raise MotifNormalizationError(
            f"Invalid class name: '{class_name}'. "
            f"Valid classes: {sorted(VALID_CLASSES)}"
        )
    else:
        warnings.warn(
            f"Unknown class name '{class_name}' - returning as-is. "
            f"Valid classes: {sorted(VALID_CLASSES)}"
        )
        return class_name


def normalize_subclass_name(subclass_name: str, strict: bool = True) -> str:
    """
    Normalize a subclass name to canonical form.
    
    Args:
        subclass_name: Raw subclass name (may have variant spelling/casing)
        strict: If True, raise error on unknown subclass. If False, return as-is with warning.
        
    Returns:
        Canonical subclass name
        
    Raises:
        MotifNormalizationError: If subclass_name is invalid and strict=True
        
    Examples:
        >>> normalize_subclass_name('telomeric g4')
        'Telomeric G4'
        >>> normalize_subclass_name('canonical_imotif')
        'Canonical i-motif'
    """
    if not subclass_name:
        if strict:
            raise MotifNormalizationError("Subclass name cannot be empty")
        return "Unknown"
    
    # Check if already canonical
    if subclass_name in VALID_SUBCLASSES:
        return subclass_name
    
    # Try case-insensitive lookup in aliases
    subclass_lower = subclass_name.lower().strip()
    if subclass_lower in SUBCLASS_ALIASES:
        return SUBCLASS_ALIASES[subclass_lower]
    
    # Try exact match in VALID_SUBCLASSES (case-insensitive)
    for valid_subclass in VALID_SUBCLASSES:
        if valid_subclass.lower() == subclass_lower:
            return valid_subclass
    
    # Not found
    if strict:
        raise MotifNormalizationError(
            f"Invalid subclass name: '{subclass_name}'. "
            f"Valid subclasses: {sorted(VALID_SUBCLASSES)}"
        )
    else:
        warnings.warn(
            f"Unknown subclass name '{subclass_name}' - returning as-is. "
            f"Valid subclasses: {sorted(VALID_SUBCLASSES)}"
        )
        return subclass_name


def normalize_class_subclass(
    class_name: str,
    subclass_name: str,
    strict: bool = True,
    auto_correct: bool = True
) -> Tuple[str, str]:
    """
    Enforce canonical class/subclass pairing.
    
    This is the main normalization function that should be called by all detectors
    before returning motif data.
    
    Args:
        class_name: Raw class name
        subclass_name: Raw subclass name
        strict: If True, raise error on invalid combinations. If False, warn and proceed.
        auto_correct: If True and subclass belongs to different class, auto-correct the class.
        
    Returns:
        Tuple of (canonical_class, canonical_subclass)
        
    Raises:
        MotifNormalizationError: If pairing is invalid and strict=True
        
    Examples:
        >>> normalize_class_subclass('G-Quadruplex', 'Telomeric G4')
        ('G-Quadruplex', 'Telomeric G4')
        
        >>> normalize_class_subclass('g-quadruplex', 'telomeric g4')
        ('G-Quadruplex', 'Telomeric G4')
        
        >>> normalize_class_subclass('Triplex', 'Telomeric G4', auto_correct=True)
        ('G-Quadruplex', 'Telomeric G4')  # Auto-corrected class
    """
    # Step 1: Normalize individual names
    canonical_class = normalize_class_name(class_name, strict=strict)
    canonical_subclass = normalize_subclass_name(subclass_name, strict=strict)
    
    # Step 2: Verify pairing is valid
    if canonical_class in VALID_CLASSES and canonical_subclass in VALID_SUBCLASSES:
        # Check if subclass belongs to this class
        if canonical_subclass not in CLASS_TO_SUBCLASSES[canonical_class]:
            # Subclass belongs to different class
            correct_class = SUBCLASS_TO_CLASS[canonical_subclass]
            
            if auto_correct:
                warnings.warn(
                    f"Class/subclass mismatch: '{canonical_subclass}' belongs to "
                    f"'{correct_class}', not '{canonical_class}'. Auto-correcting class."
                )
                canonical_class = correct_class
            elif strict:
                raise MotifNormalizationError(
                    f"Invalid class/subclass pairing: '{canonical_class}' / '{canonical_subclass}'. "
                    f"'{canonical_subclass}' belongs to '{correct_class}', not '{canonical_class}'."
                )
            else:
                warnings.warn(
                    f"Invalid class/subclass pairing: '{canonical_class}' / '{canonical_subclass}'. "
                    f"'{canonical_subclass}' belongs to '{correct_class}'."
                )
    
    return canonical_class, canonical_subclass


def normalize_motif_dict(
    motif: Dict[str, Any],
    strict: bool = False,
    auto_correct: bool = True,
    in_place: bool = False
) -> Dict[str, Any]:
    """
    Normalize Class and Subclass fields in a motif dictionary.
    
    Args:
        motif: Motif dictionary with 'Class' and 'Subclass' fields
        strict: If True, raise error on invalid data
        auto_correct: If True, auto-correct mismatched class/subclass
        in_place: If True, modify motif dict in place. If False, return new dict.
        
    Returns:
        Normalized motif dictionary
        
    Raises:
        MotifNormalizationError: If required fields missing or invalid (when strict=True)
        
    Examples:
        >>> motif = {'Class': 'g-quadruplex', 'Subclass': 'telomeric g4', 'Start': 0}
        >>> normalized = normalize_motif_dict(motif)
        >>> print(normalized['Class'], normalized['Subclass'])
        'G-Quadruplex' 'Telomeric G4'
    """
    if in_place:
        result = motif
    else:
        result = motif.copy()
    
    # Check for required fields
    if 'Class' not in result:
        if strict:
            raise MotifNormalizationError("Motif missing required 'Class' field")
        result['Class'] = "Unknown"
    
    if 'Subclass' not in result:
        if strict:
            raise MotifNormalizationError("Motif missing required 'Subclass' field")
        result['Subclass'] = "Unknown"
    
    # Normalize
    try:
        canonical_class, canonical_subclass = normalize_class_subclass(
            result['Class'],
            result['Subclass'],
            strict=strict,
            auto_correct=auto_correct
        )
        result['Class'] = canonical_class
        result['Subclass'] = canonical_subclass
    except MotifNormalizationError as e:
        if strict:
            raise
        else:
            warnings.warn(f"Motif normalization warning: {e}")
    
    return result


def validate_motif_dict(motif: Dict[str, Any], raise_on_error: bool = True) -> bool:
    """
    Validate that a motif dictionary has valid canonical Class/Subclass.
    
    Args:
        motif: Motif dictionary to validate
        raise_on_error: If True, raise exception on validation failure
        
    Returns:
        True if valid, False otherwise
        
    Raises:
        MotifNormalizationError: If validation fails and raise_on_error=True
    """
    try:
        class_name = motif.get('Class')
        subclass_name = motif.get('Subclass')
        
        if not class_name or not subclass_name:
            if raise_on_error:
                raise MotifNormalizationError("Motif missing Class or Subclass field")
            return False
        
        if class_name not in VALID_CLASSES:
            if raise_on_error:
                raise MotifNormalizationError(f"Invalid class: {class_name}")
            return False
        
        if subclass_name not in VALID_SUBCLASSES:
            if raise_on_error:
                raise MotifNormalizationError(f"Invalid subclass: {subclass_name}")
            return False
        
        if not is_valid_pairing(class_name, subclass_name):
            if raise_on_error:
                raise MotifNormalizationError(
                    f"Invalid pairing: {class_name} / {subclass_name}"
                )
            return False
        
        return True
    
    except Exception as e:
        if raise_on_error:
            raise
        return False


# =============================================================================
# MODULE METADATA
# =============================================================================

__all__ = [
    'MotifNormalizationError',
    'normalize_class_name',
    'normalize_subclass_name',
    'normalize_class_subclass',
    'normalize_motif_dict',
    'validate_motif_dict',
]
