"""
Configuration modules for NBDScanner.

This package contains all configuration constants including:
- colors: Color palettes and schemes
- themes: Color themes and tab themes
- typography: Font configurations
- layout: Layout and spacing constants
- visualization: Visualization parameters
- analysis: Analysis configuration
- export: Export format settings
- text: UI text content
- motif_taxonomy: Canonical motif classification system
"""

from .motif_taxonomy import (
    MOTIF_CLASSIFICATION,
    VALID_CLASSES,
    VALID_SUBCLASSES,
    CLASS_TO_SUBCLASSES,
    SUBCLASS_TO_CLASS,
    get_all_classes,
    get_subclasses_for_class,
    is_valid_class,
    is_valid_subclass,
    is_valid_pairing,
)

__all__ = [
    'MOTIF_CLASSIFICATION',
    'VALID_CLASSES',
    'VALID_SUBCLASSES',
    'CLASS_TO_SUBCLASSES',
    'SUBCLASS_TO_CLASS',
    'get_all_classes',
    'get_subclasses_for_class',
    'is_valid_class',
    'is_valid_subclass',
    'is_valid_pairing',
]
