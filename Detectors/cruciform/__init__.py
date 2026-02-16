"""
Cruciform DNA Structure Detection Module
=========================================

This module provides specialized detection for cruciform DNA structures,
which are four-way junctions formed by inverted repeats that can adopt
hairpin-like conformations through intramolecular base pairing.

Classes:
--------
CruciformDetector : Primary detector class for identifying cruciform structures

Usage:
------
    from Detectors.cruciform import CruciformDetector
    
    detector = CruciformDetector()
    motifs = detector.detect_motifs(sequence, "seq_name")

Scientific Background:
----------------------
Cruciforms are non-B DNA structures characterized by:
- Two palindromic arms (inverted repeats)
- Short loop regions (typically 0-3 nucleotides)
- Intramolecular base pairing forming hairpin-like structures
- Involvement in DNA replication, recombination, and transcription

Author: Dr. Venkata Rajesh Yella
License: MIT
"""

from .detector import CruciformDetector

__all__ = ['CruciformDetector']
