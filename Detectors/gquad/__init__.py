"""
G-Quadruplex Detection Subpackage

This subpackage provides functionality for detecting G-quadruplex DNA structures,
including G4Hunter-based scoring, structural classification, and overlap resolution.

Modules:
--------
- detector: Main GQuadruplexDetector class implementation

Classes:
--------
- GQuadruplexDetector: Primary detector for G-quadruplex motifs

Constants:
----------
- WINDOW_SIZE_DEFAULT: Default window size for G4Hunter scoring (25)
- MIN_REGION_LEN: Minimum length for valid G4 regions (8)
- CLASS_PRIORITY: Priority ordering for overlap resolution (telomeric > higher_order > stacked > canonical > bulged > extended_loop > g_triplex > weak_pqs)
"""

from .detector import GQuadruplexDetector, WINDOW_SIZE_DEFAULT, MIN_REGION_LEN, CLASS_PRIORITY

__all__ = [
    'GQuadruplexDetector',
    'WINDOW_SIZE_DEFAULT',
    'MIN_REGION_LEN',
    'CLASS_PRIORITY'
]
