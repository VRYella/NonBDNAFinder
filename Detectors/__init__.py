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
from Detectors.base.base_detector import BaseMotifDetector

# Import all detector classes from submodules
from Detectors.curved.detector import CurvedDNADetector
from Detectors.zdna.detector import ZDNADetector
from Detectors.aphilic.detector import APhilicDetector
from Detectors.slipped.detector import SlippedDNADetector
from Detectors.cruciform.detector import CruciformDetector
from Detectors.rloop.detector import RLoopDetector
from Detectors.triplex.detector import TriplexDetector
from Detectors.gquad.detector import GQuadruplexDetector
from Detectors.imotif.detector import IMotifDetector

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
