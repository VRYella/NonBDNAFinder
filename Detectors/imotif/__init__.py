"""
i-Motif DNA Detection Package

This package provides functionality for detecting and analyzing i-motif DNA structures,
which are quadruplex structures formed by cytosine-rich sequences.

Main Components:
---------------
- IMotifDetector: Primary detector class for i-motif structures
- Pattern recognition for canonical i-motifs (4 C-tracts)
- HUR AC-motif detection (alternating A/C patterns)
- Validated sequence matching from literature

Usage:
------
    from Detectors.imotif import IMotifDetector
    
    detector = IMotifDetector()
    motifs = detector.detect_motifs(sequence, "seq_name")
    score = detector.calculate_score(sequence)

References:
----------
See detector.py for detailed scientific references and documentation.
"""

from .detector import IMotifDetector

__all__ = ['IMotifDetector']
