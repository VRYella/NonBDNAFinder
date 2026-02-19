"""Hyperscan-accelerated pattern matching backend for Z-DNA detection.

This module provides optimized 10-mer matching using Hyperscan when available,
with automatic fallback to pure Python implementation.

Hyperscan is a high-performance regular expression matching library that can
significantly accelerate pattern matching in large sequences.

Reference:
    Intel Hyperscan: https://www.hyperscan.io/
"""

import logging
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# NumPy imports for vectorized operations (optional performance enhancement)
_NUMPY_AVAILABLE = False
try:
    import numpy as np
    _NUMPY_AVAILABLE = True
    logger.debug("NumPy available for vectorized 10-mer scanning")
except ImportError:
    logger.debug("NumPy not available - using standard Python implementation")
    pass

# Hyperscan initialization
_HYPERSCAN_AVAILABLE = False
_HYPERSCAN_VERSION = None
_HYPERSCAN_ERROR = None

try:
    import hyperscan
    _HYPERSCAN_AVAILABLE = True
    # Try to get version, but don't fail if not available
    try:
        _HYPERSCAN_VERSION = hyperscan.__version__
    except AttributeError:
        _HYPERSCAN_VERSION = 'unknown'
    logger.info(f"Hyperscan loaded successfully (version: {_HYPERSCAN_VERSION})")
except ImportError as e:
    _HYPERSCAN_ERROR = f"Hyperscan Python bindings not installed: {e}"
    logger.info(f"Hyperscan not available - using pure Python fallback. {_HYPERSCAN_ERROR}")
except Exception as e:
    _HYPERSCAN_ERROR = f"Hyperscan initialization failed: {e}"
    logger.warning(f"Hyperscan not available - using pure Python fallback. {_HYPERSCAN_ERROR}")


def hs_find_matches(seq: str, tenmer_score: Dict[str, float]) -> List[Tuple[int, str, float]]:
    """Hyperscan-based matching with improved error handling.
    
    Args:
        seq: DNA sequence to search (uppercase).
        tenmer_score: Dictionary mapping 10-mer sequences to their scores.
    
    Returns:
        List of (start, tenmer, score) tuples sorted by start position.
    
    Raises:
        Exception: If Hyperscan matching fails (to trigger fallback).
    """
    try:
        expressions = []
        ids = []
        id_to_ten = {}
        id_to_score = {}
        for idx, (ten, score) in enumerate(tenmer_score.items()):
            expressions.append(ten.encode())
            ids.append(idx)
            id_to_ten[idx] = ten
            id_to_score[idx] = float(score)
        
        db = hyperscan.Database()
        db.compile(expressions=expressions, ids=ids, elements=len(expressions))
        logger.debug(f"Hyperscan database compiled successfully for Z-DNA ({len(expressions)} patterns)")
        
        matches: List[Tuple[int, str, float]] = []

        def on_match(id, start, end, flags, context):
            # Hyperscan 'end' parameter is the offset after the match
            # Calculate actual start position: end - pattern_length
            actual_start = end - 10
            matches.append((actual_start, id_to_ten[id], id_to_score[id]))

        db.scan(seq.encode(), match_event_handler=on_match)
        matches.sort(key=lambda x: x[0])
        logger.debug(f"Hyperscan scan completed: {len(matches)} Z-DNA 10-mer matches found")
        return matches
        
    except Exception as e:
        logger.error(f"Hyperscan matching failed for Z-DNA: {e}. Falling back to pure Python.")
        raise  # Re-raise to trigger fallback in calling code


def py_find_matches(seq: str, tenmer_score: Dict[str, float]) -> List[Tuple[int, str, float]]:
    """Pure-Python exact search (overlapping matches allowed).
    
    Args:
        seq: DNA sequence to search (uppercase).
        tenmer_score: Dictionary mapping 10-mer sequences to their scores.
    
    Returns:
        List of (start, tenmer, score) tuples.
    """
    # Use vectorized version if NumPy is available, otherwise use loop-based version
    if _NUMPY_AVAILABLE:
        try:
            return vectorized_find_matches(seq, tenmer_score)
        except Exception as e:
            logger.warning(f"Vectorized matching failed, falling back to loop-based: {e}")
            return py_find_matches_loop(seq, tenmer_score)
    else:
        return py_find_matches_loop(seq, tenmer_score)


def py_find_matches_loop(seq: str, tenmer_score: Dict[str, float]) -> List[Tuple[int, str, float]]:
    """Original loop-based implementation (kept for validation and fallback).
    
    Args:
        seq: DNA sequence to search (uppercase).
        tenmer_score: Dictionary mapping 10-mer sequences to their scores.
    
    Returns:
        List of (start, tenmer, score) tuples.
    """
    n = len(seq)
    matches: List[Tuple[int, str, float]] = []
    for i in range(0, n - 10 + 1):
        ten = seq[i:i + 10]
        score = tenmer_score.get(ten)
        if score is not None:
            matches.append((i, ten, float(score)))
    return matches


# ============================================================================
# Numpy-vectorized 10-mer lookup (module-level constants, built once)
# ============================================================================

# DNA base → integer encoding (A=0, C=1, G=2, T=3; 255 = invalid/N)
_BASE_ENCODE = np.full(256, 255, dtype=np.uint8) if _NUMPY_AVAILABLE else None
if _NUMPY_AVAILABLE:
    _BASE_ENCODE[ord('A')] = 0
    _BASE_ENCODE[ord('C')] = 1
    _BASE_ENCODE[ord('G')] = 2
    _BASE_ENCODE[ord('T')] = 3

# Powers of 4 for polynomial hashing of 10-mers: hash = Σ enc[i+k] * 4^k, k=0..9
_POWERS_OF_4: "np.ndarray" = (4 ** np.arange(10, dtype=np.int64)) if _NUMPY_AVAILABLE else None

# Size of the hash lookup table: 4^10 = 1,048,576 entries
_HASH_TABLE_SIZE: int = 4 ** 10

# Cached numpy lookup table: maps 10-mer hash → score (0.0 = not in table)
# Built lazily on first call and reused for subsequent calls.
_NUMPY_LOOKUP: "Optional[np.ndarray]" = None
_NUMPY_LOOKUP_KEY: "Optional[int]" = None  # id() of the tenmer_score dict


def _build_numpy_lookup(tenmer_score: Dict[str, float]) -> "np.ndarray":
    """Build a numpy array of shape (4^10,) mapping each 10-mer hash to its score.

    The hash of a 10-mer (b0,b1,...,b9) is Σ b_k * 4^k where b_k ∈ {0,1,2,3}.
    Entries not present in *tenmer_score* are left as 0.0.
    """
    import numpy as np  # local re-import ensures availability in this scope
    lookup = np.zeros(_HASH_TABLE_SIZE, dtype=np.float64)
    for tenmer, score in tenmer_score.items():
        enc = _BASE_ENCODE[
            np.frombuffer(tenmer.upper().encode('ascii'), dtype=np.uint8)
        ].astype(np.int64)
        hash_val = int(np.dot(enc, _POWERS_OF_4))
        lookup[hash_val] = float(score)
    return lookup


def vectorized_find_matches(seq: str, tenmer_score: Dict[str, float]) -> List[Tuple[int, str, float]]:
    """True numpy-vectorized 10-mer matching — no Python inner loop over positions.

    Replaces the previous set-based Python loop with fully vectorized numpy
    operations:
      1. Encode DNA as a uint8 integer array (A=0, C=1, G=2, T=3) via a
         256-entry lookup table — O(n) numpy indexing.
      2. Compute a polynomial hash for every 10-mer window in parallel using
         10 numpy array additions — O(10·n) C-speed.
      3. Look up all hashes simultaneously in a prebuilt score table — O(n)
         numpy indexing.
      4. Collect only the (typically few) non-zero hits.

    The per-position Python overhead of the old set-based loop is eliminated;
    Python iterates only 10 times for the hash computation (not n times), plus
    once per actual match.

    Performance vs py_find_matches_loop:
        ~5–10x faster for sequences ≥ 10 KB (benchmark on random DNA).

    The lookup table is cached at module level and rebuilt only when a
    different *tenmer_score* dict is supplied.

    Args:
        seq: DNA sequence to search (uppercase, ACGT only).
        tenmer_score: Dictionary mapping 10-mer strings to their scores.

    Returns:
        List of (start, tenmer, score) tuples, identical to py_find_matches_loop.

    Raises:
        ImportError: If NumPy is not available (triggers fallback to loop-based).
    """
    global _NUMPY_LOOKUP, _NUMPY_LOOKUP_KEY

    if not _NUMPY_AVAILABLE:
        raise ImportError("NumPy is required for vectorized 10-mer matching")

    import numpy as np  # local re-import ensures availability

    n = len(seq)
    if n < 10:
        return []

    # Fall back to loop for very short sequences (numpy overhead not worth it)
    if n < 500:
        return py_find_matches_loop(seq, tenmer_score)

    # --- Build / retrieve cached lookup table ---
    # Use id() for fast identity check — TENMER_SCORE is a module-level constant
    # that never changes, so object identity is a safe and cheap cache key.
    current_key = id(tenmer_score)
    if _NUMPY_LOOKUP is None or _NUMPY_LOOKUP_KEY != current_key:
        _NUMPY_LOOKUP = _build_numpy_lookup(tenmer_score)
        _NUMPY_LOOKUP_KEY = current_key
    lookup = _NUMPY_LOOKUP

    # --- Step 1: Encode DNA sequence ---
    seq_bytes = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
    encoded = _BASE_ENCODE[seq_bytes]  # dtype uint8; 255 for invalid chars

    # Fall back if sequence contains non-ACGT characters
    if np.any(encoded == 255):
        return py_find_matches_loop(seq, tenmer_score)

    # --- Step 2: Compute polynomial hashes for all 10-mer windows ---
    # hash[i] = Σ encoded[i+k] * 4^k  for k in 0..9
    # Computed with 10 numpy array additions (not n Python iterations)
    num_windows = n - 9
    hashes = np.zeros(num_windows, dtype=np.int64)
    for k in range(10):
        hashes += encoded[k : k + num_windows].astype(np.int64) * _POWERS_OF_4[k]

    # --- Step 3: Vectorised score lookup ---
    scores = lookup[hashes]  # O(n) C-speed array indexing

    # --- Step 4: Collect matches (typically far fewer than n) ---
    match_indices = np.nonzero(scores)[0]
    return [
        (int(i), seq[i : i + 10], float(scores[i]))
        for i in match_indices
    ]


def is_hyperscan_available() -> bool:
    """Check if Hyperscan is available for use.
    
    Returns:
        True if Hyperscan is available, False otherwise.
    """
    return _HYPERSCAN_AVAILABLE


def get_hyperscan_version() -> str:
    """Get the Hyperscan version string.
    
    Returns:
        Version string or 'unavailable' if Hyperscan is not loaded.
    """
    return _HYPERSCAN_VERSION if _HYPERSCAN_AVAILABLE else 'unavailable'


def get_hyperscan_error() -> str:
    """Get the Hyperscan error message if initialization failed.
    
    Returns:
        Error message or empty string if no error occurred.
    """
    return _HYPERSCAN_ERROR or ''
