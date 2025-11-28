"""
Numba-Accelerated Backend for Non-B DNA Scanning
=================================================

This module provides Numba JIT-compiled functions for accelerated
motif detection. When Hyperscan is not available, Numba provides
significant speedup over pure Python implementations.

Key optimizations:
- JIT compilation of pattern matching loops
- Vectorized sequence encoding
- Cache-friendly memory access patterns

Requirements:
    pip install numba

If Numba is not available, functions fall back to pure Python.
"""

import logging
from typing import List, Dict, Any, Optional, Tuple
import numpy as np

logger = logging.getLogger(__name__)

# Try to import numba
try:
    import numba
    from numba import jit, prange
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    numba = None
    jit = lambda *args, **kwargs: lambda f: f
    prange = range


def is_numba_available() -> bool:
    """Check if Numba library is available."""
    return NUMBA_AVAILABLE


# Nucleotide encoding
NUCLEOTIDE_ENCODING = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
NUCLEOTIDE_DECODING = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'N'}


def encode_sequence(sequence: str) -> np.ndarray:
    """
    Encode DNA sequence as numpy array for fast processing.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        Numpy array of uint8 encoded nucleotides
    """
    sequence = sequence.upper()
    encoded = np.zeros(len(sequence), dtype=np.uint8)
    
    for i, base in enumerate(sequence):
        encoded[i] = NUCLEOTIDE_ENCODING.get(base, 4)  # 4 for unknown
    
    return encoded


def decode_sequence(encoded: np.ndarray) -> str:
    """
    Decode numpy array back to DNA string.
    
    Args:
        encoded: Numpy array of encoded nucleotides
        
    Returns:
        DNA sequence string
    """
    return ''.join(NUCLEOTIDE_DECODING.get(int(b), 'N') for b in encoded)


# JIT-compiled helper functions
if NUMBA_AVAILABLE:
    @jit(nopython=True, cache=True)
    def _count_g_runs(encoded: np.ndarray) -> Tuple[int, int]:
        """
        Count G-runs (consecutive G bases) in encoded sequence.
        
        Returns:
            Tuple of (number of G-runs of length >= 2, max G-run length)
        """
        n = len(encoded)
        g_runs = 0
        max_run = 0
        current_run = 0
        
        for i in range(n):
            if encoded[i] == 2:  # G
                current_run += 1
            else:
                if current_run >= 2:
                    g_runs += 1
                    if current_run > max_run:
                        max_run = current_run
                current_run = 0
        
        # Handle final run
        if current_run >= 2:
            g_runs += 1
            if current_run > max_run:
                max_run = current_run
        
        return g_runs, max_run
    
    @jit(nopython=True, cache=True)
    def _count_c_runs(encoded: np.ndarray) -> Tuple[int, int]:
        """
        Count C-runs (consecutive C bases) in encoded sequence.
        
        Returns:
            Tuple of (number of C-runs of length >= 2, max C-run length)
        """
        n = len(encoded)
        c_runs = 0
        max_run = 0
        current_run = 0
        
        for i in range(n):
            if encoded[i] == 1:  # C
                current_run += 1
            else:
                if current_run >= 2:
                    c_runs += 1
                    if current_run > max_run:
                        max_run = current_run
                current_run = 0
        
        # Handle final run
        if current_run >= 2:
            c_runs += 1
            if current_run > max_run:
                max_run = current_run
        
        return c_runs, max_run
    
    @jit(nopython=True, cache=True)
    def _find_g_tracts(encoded: np.ndarray, min_len: int = 3) -> List[Tuple[int, int]]:
        """
        Find all G-tracts of minimum length.
        
        Returns:
            List of (start, end) positions for G-tracts
        """
        n = len(encoded)
        tracts = []
        start = -1
        
        for i in range(n):
            if encoded[i] == 2:  # G
                if start == -1:
                    start = i
            else:
                if start >= 0:
                    length = i - start
                    if length >= min_len:
                        tracts.append((start, i))
                    start = -1
        
        # Handle final tract
        if start >= 0:
            length = n - start
            if length >= min_len:
                tracts.append((start, n))
        
        return tracts
    
    @jit(nopython=True, cache=True)
    def _find_a_tracts(encoded: np.ndarray, min_len: int = 4) -> List[Tuple[int, int]]:
        """
        Find all A-tracts of minimum length.
        
        Returns:
            List of (start, end) positions for A-tracts
        """
        n = len(encoded)
        tracts = []
        start = -1
        
        for i in range(n):
            if encoded[i] == 0:  # A
                if start == -1:
                    start = i
            else:
                if start >= 0:
                    length = i - start
                    if length >= min_len:
                        tracts.append((start, i))
                    start = -1
        
        # Handle final tract
        if start >= 0:
            length = n - start
            if length >= min_len:
                tracts.append((start, n))
        
        return tracts
    
    @jit(nopython=True, cache=True)
    def _gc_content_window(encoded: np.ndarray, window_size: int = 100) -> np.ndarray:
        """
        Calculate GC content in sliding windows.
        
        Returns:
            Array of GC content (0-1) for each window
        """
        n = len(encoded)
        if n < window_size:
            # Return single value for entire sequence
            gc = 0
            for b in encoded:
                if b == 1 or b == 2:  # C or G
                    gc += 1
            return np.array([gc / n])
        
        num_windows = n - window_size + 1
        gc_values = np.zeros(num_windows)
        
        # Calculate first window
        gc = 0
        for i in range(window_size):
            if encoded[i] == 1 or encoded[i] == 2:
                gc += 1
        gc_values[0] = gc / window_size
        
        # Slide window
        for i in range(1, num_windows):
            # Remove leaving base
            if encoded[i - 1] == 1 or encoded[i - 1] == 2:
                gc -= 1
            # Add entering base
            if encoded[i + window_size - 1] == 1 or encoded[i + window_size - 1] == 2:
                gc += 1
            gc_values[i] = gc / window_size
        
        return gc_values
    
    @jit(nopython=True, cache=True)
    def _g4_hunter_score_window(encoded: np.ndarray, window_size: int = 25) -> float:
        """
        Calculate G4Hunter-like score for a window.
        
        Returns:
            Score value (can be negative for C-rich)
        """
        n = len(encoded)
        if n == 0:
            return 0.0
        
        score = 0.0
        for i in range(n):
            if encoded[i] == 2:  # G
                score += 1.0
            elif encoded[i] == 1:  # C
                score -= 1.0
        
        return score / n

else:
    # Pure Python fallbacks when Numba not available
    def _count_g_runs(encoded):
        n = len(encoded)
        g_runs = 0
        max_run = 0
        current_run = 0
        
        for i in range(n):
            if encoded[i] == 2:
                current_run += 1
            else:
                if current_run >= 2:
                    g_runs += 1
                    if current_run > max_run:
                        max_run = current_run
                current_run = 0
        
        if current_run >= 2:
            g_runs += 1
            if current_run > max_run:
                max_run = current_run
        
        return g_runs, max_run
    
    def _count_c_runs(encoded):
        n = len(encoded)
        c_runs = 0
        max_run = 0
        current_run = 0
        
        for i in range(n):
            if encoded[i] == 1:
                current_run += 1
            else:
                if current_run >= 2:
                    c_runs += 1
                    if current_run > max_run:
                        max_run = current_run
                current_run = 0
        
        if current_run >= 2:
            c_runs += 1
            if current_run > max_run:
                max_run = current_run
        
        return c_runs, max_run
    
    def _find_g_tracts(encoded, min_len=3):
        n = len(encoded)
        tracts = []
        start = -1
        
        for i in range(n):
            if encoded[i] == 2:
                if start == -1:
                    start = i
            else:
                if start >= 0:
                    length = i - start
                    if length >= min_len:
                        tracts.append((start, i))
                    start = -1
        
        if start >= 0:
            length = n - start
            if length >= min_len:
                tracts.append((start, n))
        
        return tracts
    
    def _find_a_tracts(encoded, min_len=4):
        n = len(encoded)
        tracts = []
        start = -1
        
        for i in range(n):
            if encoded[i] == 0:
                if start == -1:
                    start = i
            else:
                if start >= 0:
                    length = i - start
                    if length >= min_len:
                        tracts.append((start, i))
                    start = -1
        
        if start >= 0:
            length = n - start
            if length >= min_len:
                tracts.append((start, n))
        
        return tracts
    
    def _gc_content_window(encoded, window_size=100):
        n = len(encoded)
        if n < window_size:
            gc = sum(1 for b in encoded if b == 1 or b == 2)
            return np.array([gc / n])
        
        num_windows = n - window_size + 1
        gc_values = np.zeros(num_windows)
        
        gc = sum(1 for i in range(window_size) if encoded[i] == 1 or encoded[i] == 2)
        gc_values[0] = gc / window_size
        
        for i in range(1, num_windows):
            if encoded[i - 1] == 1 or encoded[i - 1] == 2:
                gc -= 1
            if encoded[i + window_size - 1] == 1 or encoded[i + window_size - 1] == 2:
                gc += 1
            gc_values[i] = gc / window_size
        
        return gc_values
    
    def _g4_hunter_score_window(encoded, window_size=25):
        n = len(encoded)
        if n == 0:
            return 0.0
        
        score = 0.0
        for i in range(n):
            if encoded[i] == 2:
                score += 1.0
            elif encoded[i] == 1:
                score -= 1.0
        
        return score / n


class NumbaScanner:
    """
    Numba-accelerated scanner for Non-B DNA motifs.
    
    Uses JIT-compiled functions for fast pattern detection,
    particularly effective for:
    - G-quadruplex detection (G-tract counting)
    - i-Motif detection (C-tract counting)
    - A-tract and T-tract detection
    - GC content analysis
    """
    
    def __init__(self):
        """Initialize Numba scanner."""
        self.encoded_cache = {}
    
    def encode(self, sequence: str) -> np.ndarray:
        """
        Encode sequence with caching.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            Encoded numpy array
        """
        # Use hash for cache key (avoid storing large sequences as keys)
        key = hash(sequence)
        
        if key not in self.encoded_cache:
            self.encoded_cache[key] = encode_sequence(sequence)
        
        return self.encoded_cache[key]
    
    def find_g_quadruplex_candidates(
        self,
        sequence: str,
        min_g_tract: int = 3,
        max_loop: int = 7
    ) -> List[Dict[str, Any]]:
        """
        Find potential G-quadruplex forming regions.
        
        Uses Numba-accelerated G-tract finding followed by
        pattern assembly.
        
        Args:
            sequence: DNA sequence string
            min_g_tract: Minimum G-tract length
            max_loop: Maximum loop length
            
        Returns:
            List of G4 candidate regions
        """
        encoded = self.encode(sequence)
        g_tracts = _find_g_tracts(encoded, min_g_tract)
        
        if len(g_tracts) < 4:
            return []
        
        candidates = []
        
        # Look for 4 G-tracts with loops in between
        for i in range(len(g_tracts) - 3):
            tract1 = g_tracts[i]
            tract2 = g_tracts[i + 1]
            tract3 = g_tracts[i + 2]
            tract4 = g_tracts[i + 3]
            
            # Check loop lengths
            loop1 = tract2[0] - tract1[1]
            loop2 = tract3[0] - tract2[1]
            loop3 = tract4[0] - tract3[1]
            
            if loop1 <= max_loop and loop2 <= max_loop and loop3 <= max_loop:
                start = tract1[0]
                end = tract4[1]
                
                # Calculate score based on tract lengths and loop sizes
                tract_lengths = [tract1[1] - tract1[0], tract2[1] - tract2[0],
                               tract3[1] - tract3[0], tract4[1] - tract4[0]]
                avg_tract = sum(tract_lengths) / 4
                avg_loop = (loop1 + loop2 + loop3) / 3
                
                # Score: longer tracts and shorter loops are better
                score = min(1.0, (avg_tract / 4) * (7 / (avg_loop + 1)))
                
                candidates.append({
                    'Start': start + 1,
                    'End': end,
                    'Length': end - start,
                    'Sequence': sequence[start:end],
                    'Score': round(score, 3),
                    'Class': 'G-Quadruplex',
                    'Subclass': 'Numba-detected',
                    'G_Tracts': 4,
                    'Avg_Tract_Length': avg_tract,
                    'Avg_Loop_Length': avg_loop
                })
        
        return candidates
    
    def find_imotif_candidates(
        self,
        sequence: str,
        min_c_tract: int = 3,
        max_loop: int = 7
    ) -> List[Dict[str, Any]]:
        """
        Find potential i-Motif forming regions.
        
        Args:
            sequence: DNA sequence string
            min_c_tract: Minimum C-tract length
            max_loop: Maximum loop length
            
        Returns:
            List of i-Motif candidate regions
        """
        encoded = self.encode(sequence)
        
        # Find C-tracts
        c_tracts = []
        start = -1
        for i, b in enumerate(encoded):
            if b == 1:  # C
                if start == -1:
                    start = i
            else:
                if start >= 0:
                    if i - start >= min_c_tract:
                        c_tracts.append((start, i))
                    start = -1
        if start >= 0 and len(encoded) - start >= min_c_tract:
            c_tracts.append((start, len(encoded)))
        
        if len(c_tracts) < 4:
            return []
        
        candidates = []
        
        for i in range(len(c_tracts) - 3):
            tract1 = c_tracts[i]
            tract2 = c_tracts[i + 1]
            tract3 = c_tracts[i + 2]
            tract4 = c_tracts[i + 3]
            
            loop1 = tract2[0] - tract1[1]
            loop2 = tract3[0] - tract2[1]
            loop3 = tract4[0] - tract3[1]
            
            if loop1 <= max_loop and loop2 <= max_loop and loop3 <= max_loop:
                start = tract1[0]
                end = tract4[1]
                
                tract_lengths = [tract1[1] - tract1[0], tract2[1] - tract2[0],
                               tract3[1] - tract3[0], tract4[1] - tract4[0]]
                avg_tract = sum(tract_lengths) / 4
                avg_loop = (loop1 + loop2 + loop3) / 3
                
                score = min(1.0, (avg_tract / 4) * (7 / (avg_loop + 1)))
                
                candidates.append({
                    'Start': start + 1,
                    'End': end,
                    'Length': end - start,
                    'Sequence': sequence[start:end],
                    'Score': round(score, 3),
                    'Class': 'i-Motif',
                    'Subclass': 'Numba-detected',
                    'C_Tracts': 4,
                    'Avg_Tract_Length': avg_tract,
                    'Avg_Loop_Length': avg_loop
                })
        
        return candidates
    
    def find_a_tracts(
        self,
        sequence: str,
        min_length: int = 4
    ) -> List[Dict[str, Any]]:
        """
        Find A-tracts for curved DNA detection.
        
        Args:
            sequence: DNA sequence string
            min_length: Minimum A-tract length
            
        Returns:
            List of A-tract regions
        """
        encoded = self.encode(sequence)
        tracts = _find_a_tracts(encoded, min_length)
        
        return [
            {
                'Start': start + 1,
                'End': end,
                'Length': end - start,
                'Sequence': sequence[start:end],
                'Score': min(1.0, (end - start) / 10),
                'Class': 'Curved_DNA',
                'Subclass': 'A-tract'
            }
            for start, end in tracts
        ]
    
    def calculate_gc_content(
        self,
        sequence: str,
        window_size: int = 100
    ) -> np.ndarray:
        """
        Calculate GC content in sliding windows.
        
        Args:
            sequence: DNA sequence string
            window_size: Size of sliding window
            
        Returns:
            Array of GC content values (0-1)
        """
        encoded = self.encode(sequence)
        return _gc_content_window(encoded, window_size)
    
    def scan(
        self,
        sequence: str,
        sequence_name: str = "sequence"
    ) -> List[Dict[str, Any]]:
        """
        Perform full scan for all accelerated motif types.
        
        Args:
            sequence: DNA sequence string
            sequence_name: Name for the sequence
            
        Returns:
            Combined list of all detected motifs
        """
        all_motifs = []
        
        # G-quadruplex candidates
        g4_motifs = self.find_g_quadruplex_candidates(sequence)
        for m in g4_motifs:
            m['Sequence_Name'] = sequence_name
            m['Method'] = 'numba'
        all_motifs.extend(g4_motifs)
        
        # i-Motif candidates
        im_motifs = self.find_imotif_candidates(sequence)
        for m in im_motifs:
            m['Sequence_Name'] = sequence_name
            m['Method'] = 'numba'
        all_motifs.extend(im_motifs)
        
        # A-tracts (for curved DNA)
        a_tract_motifs = self.find_a_tracts(sequence)
        for m in a_tract_motifs:
            m['Sequence_Name'] = sequence_name
            m['Method'] = 'numba'
        all_motifs.extend(a_tract_motifs)
        
        # Sort by position
        all_motifs.sort(key=lambda m: m.get('Start', 0))
        
        return all_motifs
    
    def clear_cache(self):
        """Clear encoded sequence cache."""
        self.encoded_cache.clear()
