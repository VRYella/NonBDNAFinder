"""
Slipped DNA Detector Module
===========================

This module provides detection for slippage-prone DNA structures.

Classes:
--------
- SlippedDNADetector: Unified detector for STRs and direct repeats

Scientific Background:
---------------------
Slipped-strand DNA structures form when repetitive sequences mis-align
during DNA replication or repair, creating looped-out structures that
can lead to repeat expansion or contraction.

References:
----------
- Sinden, R.R. (1994). DNA Structure and Function. Academic Press.
- Pearson, C.E. et al. (2005). Nature Reviews Genetics 6(10):729-742.
- Mirkin, S.M. (2007). Nature 447:932-940.
"""

from .detector import SlippedDNADetector

__all__ = ['SlippedDNADetector']
