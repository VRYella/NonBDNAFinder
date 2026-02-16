"""A-philic DNA detection module.

This module provides detection of A-philic DNA regions using a 10-mer scoring table.
A-philic DNA represents sequences with unusual structural properties related to 
A-DNA formation propensity.

Main components:
- APhilicDetector: Main detector class
- TENMER_LOG2: 10-mer scoring table
"""

from .detector import APhilicDetector
from .tenmer_table import TENMER_LOG2

__all__ = ['APhilicDetector', 'TENMER_LOG2']
