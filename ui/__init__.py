"""
UI Module - User Interface Components
=====================================

This module contains all UI-related components including:
- Layout elements and page structure
- Formatting utilities
- Metrics display components
- Progress indicators
- Input widgets
- Download buttons
"""

__version__ = "2025.1"

# Import UI modules
from . import formatting
from . import downloads
from . import layout
from . import metrics
from . import progress
from . import inputs

__all__ = [
    'formatting',
    'downloads',
    'layout',
    'metrics',
    'progress',
    'inputs',
]
