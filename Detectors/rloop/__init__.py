"""R-Loop Detector Module

This module provides detection of R-loop forming sequences using the QmRLFS
(Quantitative Model of R-Loop Forming Sequences) algorithm.

R-loops are three-stranded nucleic acid structures that form naturally during
transcription and play important roles in gene regulation and genomic stability.

The detector uses Hyperscan for accelerated pattern matching when available,
with automatic fallback to pure Python regex implementation.

Classes:
    RLoopDetector: Main detector class for R-loop forming sequences

References:
    Jenjaroenpun et al. (2016) Nucleic Acids Research 43: W527-W534
    Ginno et al. (2012) Molecular Cell 45(6): 814-825
    Skourti-Stathaki & Proudfoot (2014) Genes & Development 28(13): 1384-1396
"""

from .detector import RLoopDetector

__all__ = ['RLoopDetector']
