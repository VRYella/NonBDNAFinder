"""
Scanner Backends for High-Performance Non-B DNA Detection
==========================================================

This module provides accelerated scanning backends for Non-B DNA motif detection:
- Hyperscan-based scanning (fastest, requires hyperscan library)
- Numba-accelerated detectors (fallback when Hyperscan unavailable)
- Parallel worker with shared memory support
- I/O utilities for memmap and streaming FASTA

The backends automatically fallback to available implementations.
"""

from .io_utils import mmap_fasta, stream_sequence_chunks, get_overlap_size
from .parallel_worker import SharedMemoryWorker, chunk_sequence, parallel_scan

# Try to import Hyperscan backend
try:
    from .hyperscan_backend import HyperscanScanner, is_hyperscan_available
    HYPERSCAN_AVAILABLE = is_hyperscan_available()
except ImportError:
    HYPERSCAN_AVAILABLE = False
    HyperscanScanner = None
    is_hyperscan_available = lambda: False

# Try to import Numba backend
try:
    from .numba_backend import NumbaScanner, is_numba_available
    NUMBA_AVAILABLE = is_numba_available()
except ImportError:
    NUMBA_AVAILABLE = False
    NumbaScanner = None
    is_numba_available = lambda: False


def get_best_backend():
    """
    Get the best available scanning backend.
    
    Priority:
    1. Hyperscan (fastest, if available)
    2. Numba (accelerated fallback)
    3. Pure Python (always available)
    
    Returns:
        str: Backend name ('hyperscan', 'numba', or 'python')
    """
    if HYPERSCAN_AVAILABLE:
        return 'hyperscan'
    elif NUMBA_AVAILABLE:
        return 'numba'
    else:
        return 'python'


def get_scanner(backend='auto'):
    """
    Get a scanner instance for the specified backend.
    
    Args:
        backend: Backend name ('auto', 'hyperscan', 'numba', 'python')
        
    Returns:
        Scanner instance for the specified backend
    """
    if backend == 'auto':
        backend = get_best_backend()
    
    if backend == 'hyperscan' and HYPERSCAN_AVAILABLE:
        return HyperscanScanner()
    elif backend == 'numba' and NUMBA_AVAILABLE:
        return NumbaScanner()
    else:
        # Return None to indicate pure Python should be used
        return None


__all__ = [
    'mmap_fasta',
    'stream_sequence_chunks',
    'get_overlap_size',
    'SharedMemoryWorker',
    'chunk_sequence',
    'parallel_scan',
    'HyperscanScanner',
    'NumbaScanner',
    'is_hyperscan_available',
    'is_numba_available',
    'get_best_backend',
    'get_scanner',
    'HYPERSCAN_AVAILABLE',
    'NUMBA_AVAILABLE',
]
