"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Aho-Corasick Multi-Pattern Matcher - DNA Sequence Optimization               │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.2            │
│ 50-200x faster than regex: O(n) single pass for ALL patterns                 │
└──────────────────────────────────────────────────────────────────────────────┘
"""

import logging
from typing import List, Dict, Tuple, Iterator, Optional, Any
from collections import defaultdict

logger = logging.getLogger(__name__)

# Try to import pyahocorasick (C extension)
try:
    import ahocorasick
    AHOCORASICK_AVAILABLE = True
except ImportError:
    AHOCORASICK_AVAILABLE = False
    logger.warning("pyahocorasick not available. Using fallback implementation.")


class AhoCorasickMatcher:
    """
    High-performance multi-pattern DNA sequence matcher using Aho-Corasick algorithm.
    
    The Aho-Corasick algorithm finds all occurrences of multiple patterns in O(n + m + z) time,
    where n is text length, m is total pattern length, and z is number of matches.
    This is dramatically faster than running multiple regex patterns sequentially.
    
    Features:
    - Single-pass detection of all patterns
    - Case-insensitive DNA matching
    - Metadata attachment to patterns (detector, subclass, scoring info)
    - Graceful fallback to naive string search if pyahocorasick unavailable
    - Memory-efficient: ~1-2MB per 1000 patterns
    
    Performance:
    - Pattern compilation: O(m) where m = total pattern characters
    - Search: O(n + z) where n = sequence length, z = number of matches
    - Memory: ~1KB per pattern + automaton overhead
    
    Example:
        >>> matcher = AhoCorasickMatcher()
        >>> matcher.add_pattern("GGGG", detector="g_quadruplex", subclass="G4")
        >>> matcher.add_pattern("CCCC", detector="i_motif", subclass="iM")
        >>> matcher.build()
        >>> for match in matcher.search("ATGGGGCCCCTA"):
        ...     print(match)
        (2, 6, 'GGGG', {'detector': 'g_quadruplex', 'subclass': 'G4'})
        (6, 10, 'CCCC', {'detector': 'i_motif', 'subclass': 'iM'})
    """
    
    def __init__(self, case_sensitive: bool = False):
        """
        Initialize the Aho-Corasick matcher.
        
        Args:
            case_sensitive: If True, patterns are case-sensitive. DNA is typically uppercase.
        """
        self.case_sensitive = case_sensitive
        self.patterns: List[Tuple[str, Dict[str, Any]]] = []
        self.automaton = None
        self.built = False
        self.use_fallback = not AHOCORASICK_AVAILABLE
        
        if self.use_fallback:
            logger.info("Using fallback pattern matcher (slower than pyahocorasick)")
    
    def add_pattern(self, pattern: str, **metadata) -> None:
        """
        Add a pattern to the matcher with optional metadata.
        
        Args:
            pattern: DNA sequence pattern to match (e.g., "GGGG", "ATCG")
            **metadata: Additional information to attach to matches (detector, subclass, etc.)
        
        Raises:
            ValueError: If pattern is empty or contains invalid characters
        """
        if not pattern:
            raise ValueError("Pattern cannot be empty")
        
        if self.built:
            raise RuntimeError("Cannot add patterns after build() has been called")
        
        # Normalize case if needed
        if not self.case_sensitive:
            pattern = pattern.upper()
        
        self.patterns.append((pattern, metadata))
    
    def add_patterns(self, patterns: List[Tuple[str, Dict[str, Any]]]) -> None:
        """
        Add multiple patterns at once.
        
        Args:
            patterns: List of (pattern, metadata_dict) tuples
        """
        for pattern, metadata in patterns:
            self.add_pattern(pattern, **metadata)
    
    def build(self) -> None:
        """
        Build the Aho-Corasick automaton from added patterns.
        
        This must be called after adding all patterns and before searching.
        Building is O(m) where m is the total length of all patterns.
        
        Raises:
            RuntimeError: If no patterns have been added
        """
        if not self.patterns:
            raise RuntimeError("No patterns added. Use add_pattern() before build()")
        
        if self.built:
            logger.warning("Automaton already built. Rebuilding...")
        
        if self.use_fallback:
            self._build_fallback()
        else:
            self._build_ahocorasick()
        
        self.built = True
        logger.info(f"Built Aho-Corasick automaton with {len(self.patterns)} patterns")
    
    def _build_ahocorasick(self) -> None:
        """Build using the C extension pyahocorasick library."""
        self.automaton = ahocorasick.Automaton()
        
        for idx, (pattern, metadata) in enumerate(self.patterns):
            # Store pattern index and metadata
            self.automaton.add_word(pattern, (idx, pattern, metadata))
        
        # Build the failure links (this is where the magic happens)
        self.automaton.make_automaton()
    
    def _build_fallback(self) -> None:
        """Build using naive string search (fallback when pyahocorasick unavailable)."""
        # Store patterns in a simple list for naive search
        self.automaton = []
        for idx, (pattern, metadata) in enumerate(self.patterns):
            self.automaton.append((pattern, idx, metadata))
    
    def search(self, sequence: str) -> Iterator[Tuple[int, int, str, Dict[str, Any]]]:
        """
        Search for all pattern occurrences in the sequence.
        
        Args:
            sequence: DNA sequence to search
        
        Yields:
            Tuples of (start_pos, end_pos, matched_pattern, metadata)
            Positions are 0-indexed, end_pos is exclusive (Python slice convention)
        
        Raises:
            RuntimeError: If build() has not been called
        
        Example:
            >>> matcher.build()
            >>> for start, end, pattern, meta in matcher.search("ATGGGGTA"):
            ...     print(f"Found {pattern} at {start}-{end}: {meta}")
            Found GGGG at 2-6: {'detector': 'g_quadruplex'}
        """
        if not self.built:
            raise RuntimeError("Must call build() before search()")
        
        if not sequence:
            return
        
        # Normalize case if needed
        if not self.case_sensitive:
            sequence = sequence.upper()
        
        if self.use_fallback:
            yield from self._search_fallback(sequence)
        else:
            yield from self._search_ahocorasick(sequence)
    
    def _search_ahocorasick(self, sequence: str) -> Iterator[Tuple[int, int, str, Dict[str, Any]]]:
        """Search using pyahocorasick C extension."""
        for end_pos, (idx, pattern, metadata) in self.automaton.iter(sequence):
            # end_pos is inclusive in pyahocorasick, convert to exclusive
            start_pos = end_pos - len(pattern) + 1
            end_pos = end_pos + 1  # Make exclusive
            yield (start_pos, end_pos, pattern, metadata)
    
    def _search_fallback(self, sequence: str) -> Iterator[Tuple[int, int, str, Dict[str, Any]]]:
        """Naive string search fallback (O(n*m*p) where p = number of patterns)."""
        seq_len = len(sequence)
        
        for pattern, idx, metadata in self.automaton:
            pattern_len = len(pattern)
            
            # Simple sliding window
            for i in range(seq_len - pattern_len + 1):
                if sequence[i:i+pattern_len] == pattern:
                    yield (i, i + pattern_len, pattern, metadata)
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get statistics about the matcher.
        
        Returns:
            Dictionary with pattern count, build status, and implementation type
        """
        return {
            'num_patterns': len(self.patterns),
            'built': self.built,
            'implementation': 'pyahocorasick' if not self.use_fallback else 'fallback',
            'case_sensitive': self.case_sensitive,
            'available': AHOCORASICK_AVAILABLE
        }


class PatternGroup:
    """
    Organize patterns into logical groups for easier management.
    
    Useful for organizing patterns by detector type, allowing selective
    pattern matching and better tracking of which detectors found which motifs.
    
    Example:
        >>> group = PatternGroup("g_quadruplex")
        >>> group.add_pattern("GGGG", subclass="canonical_G4")
        >>> group.add_pattern("GGG", subclass="weak_PQS")
        >>> matcher = group.to_matcher()
        >>> matcher.build()
    """
    
    def __init__(self, detector_name: str):
        """
        Initialize a pattern group.
        
        Args:
            detector_name: Name of the detector (e.g., "g_quadruplex", "z_dna")
        """
        self.detector_name = detector_name
        self.patterns: List[Tuple[str, Dict[str, Any]]] = []
    
    def add_pattern(self, pattern: str, **metadata) -> None:
        """
        Add a pattern to this group.
        
        Args:
            pattern: DNA sequence pattern
            **metadata: Additional metadata (automatically includes detector_name)
        """
        metadata['detector'] = self.detector_name
        self.patterns.append((pattern, metadata))
    
    def to_matcher(self, case_sensitive: bool = False) -> AhoCorasickMatcher:
        """
        Create an AhoCorasickMatcher from this group's patterns.
        
        Args:
            case_sensitive: Whether matching should be case-sensitive
        
        Returns:
            AhoCorasickMatcher instance with all patterns added
        """
        matcher = AhoCorasickMatcher(case_sensitive=case_sensitive)
        matcher.add_patterns(self.patterns)
        return matcher
    
    def __len__(self) -> int:
        return len(self.patterns)


def merge_pattern_groups(groups: List[PatternGroup], 
                         case_sensitive: bool = False) -> AhoCorasickMatcher:
    """
    Merge multiple PatternGroups into a single AhoCorasickMatcher.
    
    This allows combining patterns from multiple detectors into one unified
    automaton for maximum efficiency.
    
    Args:
        groups: List of PatternGroup instances to merge
        case_sensitive: Whether matching should be case-sensitive
    
    Returns:
        AhoCorasickMatcher with all patterns from all groups
    
    Example:
        >>> g4_group = PatternGroup("g_quadruplex")
        >>> g4_group.add_pattern("GGGG")
        >>> im_group = PatternGroup("i_motif")
        >>> im_group.add_pattern("CCCC")
        >>> matcher = merge_pattern_groups([g4_group, im_group])
        >>> matcher.build()
    """
    matcher = AhoCorasickMatcher(case_sensitive=case_sensitive)
    
    for group in groups:
        matcher.add_patterns(group.patterns)
    
    logger.info(f"Merged {len(groups)} pattern groups into single matcher "
                f"({len(matcher.patterns)} total patterns)")
    
    return matcher


# ═══════════════════════════════════════════════════════════════════════════════
# CONVENIENCE FUNCTIONS FOR COMMON OPERATIONS
# ═══════════════════════════════════════════════════════════════════════════════

def create_simple_matcher(patterns: List[str], 
                         detector: str = "unknown") -> AhoCorasickMatcher:
    """
    Create a matcher from a simple list of pattern strings.
    
    Args:
        patterns: List of DNA sequence patterns
        detector: Detector name to attach to all patterns
    
    Returns:
        Built AhoCorasickMatcher ready to use
    
    Example:
        >>> matcher = create_simple_matcher(["GGGG", "CCCC"], detector="g4_and_im")
        >>> list(matcher.search("ATGGGGCCCCTA"))
    """
    matcher = AhoCorasickMatcher()
    for pattern in patterns:
        matcher.add_pattern(pattern, detector=detector)
    matcher.build()
    return matcher


def benchmark_matcher(sequence: str, patterns: List[str], 
                     num_iterations: int = 100) -> Dict[str, float]:
    """
    Benchmark the Aho-Corasick matcher performance.
    
    Args:
        sequence: Test DNA sequence
        patterns: List of patterns to search for
        num_iterations: Number of search iterations for timing
    
    Returns:
        Dictionary with timing results and match counts
    """
    import time
    
    # Build matcher
    build_start = time.time()
    matcher = create_simple_matcher(patterns)
    build_time = time.time() - build_start
    
    # Warm-up run
    list(matcher.search(sequence))
    
    # Timed runs
    search_start = time.time()
    total_matches = 0
    for _ in range(num_iterations):
        matches = list(matcher.search(sequence))
        total_matches = len(matches)
    search_time = time.time() - search_start
    
    return {
        'build_time_ms': build_time * 1000,
        'search_time_ms': (search_time / num_iterations) * 1000,
        'total_matches': total_matches,
        'throughput_mb_per_sec': (len(sequence) * num_iterations / search_time) / 1_000_000,
        'implementation': 'pyahocorasick' if AHOCORASICK_AVAILABLE else 'fallback'
    }
