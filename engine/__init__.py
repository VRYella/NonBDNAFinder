"""
Engine Module - Core Detection Logic
====================================

This module contains all backend engine components including:
- Detection orchestration
- Scoring algorithms
- Pattern processing
- Overlap merging and deduplication
- Sequence chunking
- Motif detectors
"""

__version__ = "2025.1"

# Import core engine modules
from . import scoring
from . import merging
from . import chunking
from . import sequence_ops
from . import detectors

__all__ = [
    'scoring',
    'merging',
    'chunking',
    'sequence_ops',
    'detectors',
]
