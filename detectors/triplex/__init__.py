"""
Triplex DNA Detector Subpackage
================================

This subpackage provides detection capabilities for triplex-forming DNA structures,
including mirror repeats and sticky DNA sequences.

Classes
-------
TriplexDetector : BaseMotifDetector
    Main detector class for triplex DNA structures
"""

from .detector import TriplexDetector

__all__ = ['TriplexDetector']
