"""
Detectors Module - Non-B DNA Motif Detectors
============================================

This module contains all specialized motif detector classes:
- CurvedDNADetector
- ZDNADetector
- APhilicDetector
- SlippedDNADetector
- CruciformDetector
- RLoopDetector
- TriplexDetector
- GQuadruplexDetector
- IMotifDetector
"""

__version__ = "2025.1"

# Import base detector
from .base import BaseMotifDetector

__all__ = [
    'BaseMotifDetector',
]
