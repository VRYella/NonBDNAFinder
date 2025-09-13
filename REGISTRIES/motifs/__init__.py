"""
NBDFinder motifs package - Core data structures and registry
"""

from .base import Candidate
from .registry import get_patterns_for_motif, get_all_hyperscan_patterns

# Import get_basic_stats from utils to maintain compatibility
try:
    from BASE_CODES.utils import get_basic_stats
except ImportError:
    # Fallback if utils not available
    get_basic_stats = None

__all__ = ['Candidate', 'get_patterns_for_motif', 'get_all_hyperscan_patterns', 'get_basic_stats']