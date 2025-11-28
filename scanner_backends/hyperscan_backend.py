"""
Hyperscan Backend for High-Performance Pattern Matching
========================================================

This module provides Hyperscan-accelerated pattern matching for Non-B DNA
detection. Hyperscan is Intel's high-performance regex matching library
that can scan hundreds of patterns simultaneously.

When Hyperscan is available, this backend provides significant speedup
over pure Python regex matching, especially for large genomes.

Requirements:
    pip install hyperscan

If Hyperscan is not available, the module gracefully falls back to
pure Python matching.
"""

import logging
from typing import List, Dict, Any, Optional, Tuple, Callable
import re

logger = logging.getLogger(__name__)

# Try to import hyperscan
try:
    import hyperscan
    HYPERSCAN_AVAILABLE = True
except ImportError:
    HYPERSCAN_AVAILABLE = False
    hyperscan = None


def is_hyperscan_available() -> bool:
    """Check if Hyperscan library is available."""
    return HYPERSCAN_AVAILABLE


class HyperscanScanner:
    """
    High-performance scanner using Intel Hyperscan.
    
    This scanner compiles patterns into a Hyperscan database for
    fast multi-pattern matching. It's particularly effective for:
    - Scanning many patterns simultaneously
    - Large sequences (>1MB)
    - Repeated scans with same patterns
    
    Usage:
        scanner = HyperscanScanner()
        scanner.compile_patterns(patterns)
        matches = scanner.scan(sequence)
    """
    
    def __init__(self):
        """Initialize Hyperscan scanner."""
        if not HYPERSCAN_AVAILABLE:
            raise ImportError(
                "Hyperscan is not available. Install with: pip install hyperscan"
            )
        
        self.db = None
        self.patterns = {}
        self.id_to_pattern = {}
        self.id_to_info = {}
    
    def compile_patterns(
        self,
        patterns: List[Tuple[str, str, Dict[str, Any]]]
    ) -> bool:
        """
        Compile patterns into Hyperscan database.
        
        Args:
            patterns: List of (pattern_regex, pattern_id, metadata_dict)
            
        Returns:
            True if compilation successful
        """
        if not patterns:
            logger.warning("No patterns to compile")
            return False
        
        try:
            expressions = []
            ids = []
            flags = []
            
            for i, (regex, pattern_id, metadata) in enumerate(patterns):
                expressions.append(regex.encode('ascii'))
                ids.append(i)
                flags.append(hyperscan.HS_FLAG_CASELESS | hyperscan.HS_FLAG_DOTALL)
                
                self.id_to_pattern[i] = regex
                self.id_to_info[i] = {
                    'pattern_id': pattern_id,
                    **metadata
                }
            
            # Compile database
            self.db = hyperscan.Database()
            self.db.compile(
                expressions=expressions,
                ids=ids,
                elements=len(expressions),
                flags=flags
            )
            
            logger.info(f"Compiled {len(patterns)} patterns into Hyperscan database")
            return True
            
        except Exception as e:
            logger.error(f"Hyperscan compilation failed: {e}")
            self.db = None
            return False
    
    def scan(
        self,
        sequence: str,
        sequence_name: str = "sequence"
    ) -> List[Dict[str, Any]]:
        """
        Scan sequence using compiled Hyperscan database.
        
        Args:
            sequence: DNA sequence string
            sequence_name: Name for the sequence
            
        Returns:
            List of match dictionaries
        """
        if self.db is None:
            logger.warning("No compiled database, returning empty results")
            return []
        
        sequence = sequence.upper()
        matches = []
        
        def on_match(pattern_id, start, end, flags, context):
            """Callback for Hyperscan matches."""
            info = self.id_to_info.get(pattern_id, {})
            
            match = {
                'Start': start + 1,  # 1-based
                'End': end,
                'Length': end - start,
                'Sequence': sequence[start:end],
                'Pattern_ID': info.get('pattern_id', f'HS_{pattern_id}'),
                'Class': info.get('class', 'Unknown'),
                'Subclass': info.get('subclass', 'Unknown'),
                'Score': info.get('score', 0.5),
                'Sequence_Name': sequence_name,
                'Method': 'hyperscan'
            }
            matches.append(match)
        
        try:
            self.db.scan(
                sequence.encode('ascii'),
                match_event_handler=on_match
            )
        except Exception as e:
            logger.error(f"Hyperscan scan failed: {e}")
        
        # Sort by position
        matches.sort(key=lambda m: m.get('Start', 0))
        
        return matches
    
    def scan_with_scoring(
        self,
        sequence: str,
        sequence_name: str = "sequence",
        score_function: Optional[Callable] = None
    ) -> List[Dict[str, Any]]:
        """
        Scan sequence and apply custom scoring to matches.
        
        Args:
            sequence: DNA sequence string
            sequence_name: Name for the sequence
            score_function: Function(match_dict) -> score
            
        Returns:
            List of scored match dictionaries
        """
        matches = self.scan(sequence, sequence_name)
        
        if score_function is not None:
            for match in matches:
                match['Score'] = score_function(match)
        
        return matches


class HyperscanPatternDatabase:
    """
    Manages Hyperscan pattern databases for different motif classes.
    
    This class handles compilation and caching of Hyperscan databases
    for different motif types, allowing efficient switching between
    pattern sets.
    """
    
    def __init__(self, cache_dir: Optional[str] = None):
        """
        Initialize pattern database.
        
        Args:
            cache_dir: Directory to cache compiled databases
        """
        if not HYPERSCAN_AVAILABLE:
            raise ImportError("Hyperscan is not available")
        
        self.cache_dir = cache_dir
        self.databases = {}
    
    def get_database(self, motif_class: str) -> Optional[HyperscanScanner]:
        """
        Get or create Hyperscan database for a motif class.
        
        Args:
            motif_class: Name of motif class (e.g., 'g_quadruplex')
            
        Returns:
            HyperscanScanner instance or None
        """
        if motif_class in self.databases:
            return self.databases[motif_class]
        
        # Get patterns for this class
        patterns = self._get_patterns_for_class(motif_class)
        
        if not patterns:
            return None
        
        # Create and compile scanner
        scanner = HyperscanScanner()
        if scanner.compile_patterns(patterns):
            self.databases[motif_class] = scanner
            return scanner
        
        return None
    
    def _get_patterns_for_class(
        self,
        motif_class: str
    ) -> List[Tuple[str, str, Dict[str, Any]]]:
        """
        Get patterns for a specific motif class.
        
        This should be connected to the pattern registry.
        """
        # Import from utilities
        try:
            from utilities import PatternRegistry
            
            all_patterns = PatternRegistry.get_all_patterns()
            
            if motif_class not in all_patterns:
                return []
            
            result = []
            class_patterns = all_patterns[motif_class]
            
            for pattern_group, patterns in class_patterns.items():
                for pattern_tuple in patterns:
                    if len(pattern_tuple) >= 4:
                        regex = pattern_tuple[0]
                        pattern_id = pattern_tuple[1]
                        name = pattern_tuple[2]
                        subclass = pattern_tuple[3]
                        
                        # Skip patterns that aren't Hyperscan-compatible
                        if not self._is_hyperscan_compatible(regex):
                            continue
                        
                        metadata = {
                            'class': motif_class,
                            'subclass': subclass,
                            'name': name,
                            'score': pattern_tuple[6] if len(pattern_tuple) > 6 else 0.5
                        }
                        
                        result.append((regex, pattern_id, metadata))
            
            return result
            
        except ImportError:
            logger.warning("Could not import PatternRegistry")
            return []
    
    def _is_hyperscan_compatible(self, pattern: str) -> bool:
        """
        Check if a regex pattern is compatible with Hyperscan.
        
        Hyperscan doesn't support some regex features like backreferences.
        """
        return not any(inc in pattern for inc in HYPERSCAN_INCOMPATIBLE_FEATURES)


# Regex features not supported by Hyperscan
HYPERSCAN_INCOMPATIBLE_FEATURES = [
    '\\1', '\\2', '\\3',  # Backreferences
    '(?=', '(?!',         # Lookahead
    '(?<=', '(?<!',       # Lookbehind
    '\\b', '\\B'          # Word boundaries
]


def compile_all_patterns() -> Dict[str, HyperscanScanner]:
    """
    Compile all available patterns into Hyperscan databases.
    
    Returns:
        Dictionary mapping motif class to compiled scanner
    """
    if not HYPERSCAN_AVAILABLE:
        logger.warning("Hyperscan not available")
        return {}
    
    db = HyperscanPatternDatabase()
    
    # List of motif classes to compile
    motif_classes = [
        'curved_dna',
        'slipped_dna',
        'cruciform',
        'r_loop',
        'triplex',
        'g_quadruplex',
        'i_motif',
        'z_dna',
        'a_philic'
    ]
    
    compiled = {}
    for cls in motif_classes:
        scanner = db.get_database(cls)
        if scanner is not None:
            compiled[cls] = scanner
    
    return compiled


def fallback_regex_scan(
    sequence: str,
    patterns: List[Tuple[str, str, Dict[str, Any]]],
    sequence_name: str = "sequence"
) -> List[Dict[str, Any]]:
    """
    Fallback regex scanning when Hyperscan is not available.
    
    Args:
        sequence: DNA sequence string
        patterns: List of (pattern_regex, pattern_id, metadata_dict)
        sequence_name: Name for the sequence
        
    Returns:
        List of match dictionaries
    """
    sequence = sequence.upper()
    matches = []
    
    for regex, pattern_id, metadata in patterns:
        try:
            compiled = re.compile(regex, re.IGNORECASE | re.ASCII)
            
            for match in compiled.finditer(sequence):
                start, end = match.span()
                
                matches.append({
                    'Start': start + 1,  # 1-based
                    'End': end,
                    'Length': end - start,
                    'Sequence': match.group(),
                    'Pattern_ID': pattern_id,
                    'Class': metadata.get('class', 'Unknown'),
                    'Subclass': metadata.get('subclass', 'Unknown'),
                    'Score': metadata.get('score', 0.5),
                    'Sequence_Name': sequence_name,
                    'Method': 'regex'
                })
                
        except re.error as e:
            logger.warning(f"Invalid regex pattern {pattern_id}: {e}")
    
    matches.sort(key=lambda m: m.get('Start', 0))
    return matches
