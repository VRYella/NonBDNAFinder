"""
Consolidated Non-B DNA Motif Detectors Module
Dr. Venkata Rajesh Yella | 2024.1 | MIT License

This module provides backward-compatible imports.
All detector classes are re-exported for existing code.

Contains 10 detector classes for Non-B DNA structures:
BaseMotifDetector, CurvedDNADetector, ZDNADetector, APhilicDetector,
SlippedDNADetector, CruciformDetector, RLoopDetector, TriplexDetector,
GQuadruplexDetector, IMotifDetector

Performance: 5K-280K bp/s | Memory: ~5 MB/100K sequences
"""

# Import base detector
from detectors.base.base_detector import BaseMotifDetector

# Import all detector classes from submodules
from detectors.curved.detector import CurvedDNADetector
from detectors.zdna.detector import ZDNADetector
from detectors.aphilic.detector import APhilicDetector
from detectors.slipped.detector import SlippedDNADetector
from detectors.cruciform.detector import CruciformDetector
from detectors.rloop.detector import RLoopDetector
from detectors.triplex.detector import TriplexDetector
from detectors.gquad.detector import GQuadruplexDetector
from detectors.imotif.detector import IMotifDetector

__all__ = [
    "BaseMotifDetector",
    "CurvedDNADetector",
    "ZDNADetector",
    "APhilicDetector",
    "SlippedDNADetector",
    "CruciformDetector",
    "RLoopDetector",
    "TriplexDetector",
    "GQuadruplexDetector",
    "IMotifDetector",
]

__version__ = "2024.1"
__author__ = "Dr. Venkata Rajesh Yella"
