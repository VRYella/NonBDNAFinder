"""
Caching utilities for NBDScanner.

This module contains caching functions for:
- Genome sequences (as NumPy arrays)
- Hyperscan databases
- Analysis results
- Statistics
"""

import streamlit as st
import numpy as np
import logging

logger = logging.getLogger(__name__)

# Check for Hyperscan availability
HYPERSCAN_AVAILABLE = False
try:
    import hyperscan
    HYPERSCAN_AVAILABLE = True
except ImportError:
    pass

# Import analysis function from nonbscanner
from Utilities.nonbscanner import analyze_sequence
from Utilities.utilities import get_basic_stats


@st.cache_resource(show_spinner=False)
def cache_genome_as_numpy(sequence: str) -> np.ndarray:
    """
    Cache genome sequence as NumPy byte array for memory efficiency.
    
    This prevents reloading large genomes and reduces memory footprint
    when using Streamlit on free tier (1GB limit).
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        NumPy array of sequence bytes
    """
    return np.frombuffer(sequence.encode('utf-8'), dtype=np.uint8)


@st.cache_resource(show_spinner=False)
def cache_hyperscan_database(_patterns: list = None):
    """
    Cache compiled Hyperscan database for pattern matching.
    
    The underscore prefix on _patterns parameter is used by Streamlit
    to indicate the parameter should not be hashed for caching purposes.
    This prevents Streamlit from attempting to hash complex pattern objects.
    
    Args:
        _patterns: List of (pattern, pattern_id) tuples to compile
        
    Returns:
        Compiled Hyperscan database or None
    """
    if not HYPERSCAN_AVAILABLE or _patterns is None or len(_patterns) == 0:
        logger.debug("Hyperscan not available or no patterns provided")
        return None
    
    try:
        # Compile patterns into Hyperscan database
        logger.debug(f"Compiling Hyperscan database with {len(_patterns)} patterns...")
        
        expressions = []
        ids = []
        flags = []
        
        for pattern, pattern_id in _patterns:
            # Encode pattern with error handling
            # Note: DNA patterns should normally be ASCII (ATGC), but we provide
            # UTF-8 fallback for robustness in case patterns contain metadata or
            # special characters. A warning is logged to help identify data quality issues.
            try:
                pattern_bytes = pattern.encode('ascii')
            except UnicodeEncodeError:
                # Fall back to UTF-8 if ASCII fails
                pattern_bytes = pattern.encode('utf-8')
                logger.warning(f"Pattern {pattern_id} contains non-ASCII characters (expected ATGC). Using UTF-8 encoding.")
            
            expressions.append(pattern_bytes)
            ids.append(pattern_id)
            # Use CASELESS and DOTALL flags for DNA matching
            flags.append(hyperscan.HS_FLAG_CASELESS | hyperscan.HS_FLAG_DOTALL)
        
        db = hyperscan.Database()
        db.compile(
            expressions=expressions,
            ids=ids,
            elements=len(expressions),
            flags=flags
        )
        
        logger.info(f"Successfully compiled Hyperscan database with {len(expressions)} patterns")
        return db
        
    except Exception as e:
        logger.error(f"Hyperscan database compilation failed: {e}")
        st.warning(f"Hyperscan database compilation failed: {e}. Falling back to regex matching.")
        return None


@st.cache_data(show_spinner=False, max_entries=10, ttl=3600)
def cache_analysis_results(sequence_hash: str, sequence: str, name: str):
    """
    Cache analysis results for sequences to avoid re-computation.
    Uses sequence hash for efficient cache key lookup.
    
    Args:
        sequence_hash: Hash of sequence for cache key
        sequence: DNA sequence string
        name: Sequence name
        
    Returns:
        List of detected motifs
    """
    return analyze_sequence(sequence, name)


@st.cache_data(show_spinner=False, max_entries=20, ttl=3600)
def get_cached_stats(sequence: str, motifs_json: str):
    """
    Cache statistics calculation for sequences.
    
    Args:
        sequence: DNA sequence
        motifs_json: JSON string of motifs (for cache key)
        
    Returns:
        Dictionary of sequence statistics
    """
    import json
    motifs = json.loads(motifs_json) if motifs_json else []
    return get_basic_stats(sequence, motifs)
