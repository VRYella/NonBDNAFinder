"""Core modules for NonBDNAFinder"""

from .motif_normalizer import (
    MotifNormalizationError,
    normalize_class_name,
    normalize_subclass_name,
    normalize_class_subclass,
    normalize_motif_dict,
    validate_motif_dict,
)

__all__ = [
    'MotifNormalizationError',
    'normalize_class_name',
    'normalize_subclass_name',
    'normalize_class_subclass',
    'normalize_motif_dict',
    'validate_motif_dict',
]
