"""
All 10 motif detector classes - modular design with detect() and score() methods.

Each detector uses Hyperscan for fast anchoring, then applies specific validation
and scoring algorithms per the TECHNICAL_SPECIFICATIONS.md.
"""

from typing import List, Tuple, Dict, Optional
import re
import logging
import math
from abc import ABC, abstractmethod
from collections import defaultdict
import numpy as np

# Add imports near top
try:
    import hyperscan
    HYPERSCAN_AVAILABLE = True
except Exception:
    HYPERSCAN_AVAILABLE = False

try:
    from intervaltree import IntervalTree
    INTERVALTREE_AVAILABLE = True
except ImportError:
    INTERVALTREE_AVAILABLE = False
    logging.warning("intervaltree not available - using simple overlap detection")

from motifs.base import (
    Candidate, g4hunter_score, triplex_stability_score, 
    z_dna_score, curvature_score
)
from motifs.registry import get_patterns_for_motif, get_class_id
import hs_registry_manager

logger = logging.getLogger(__name__)

# Parameters for R-loop detection
params = {
    'min_perc_g_riz': 0.6,  # Minimum G percentage for R-loop candidates
    'regex_models': {}  # Will be populated with compiled regex models
}

# Utility functions for R-loop detection
def percent_g(sequence):
    """Calculate percentage of G nucleotides in sequence"""
    if not sequence:
        return 0.0
    g_count = sequence.upper().count('G')
    return g_count / len(sequence)


def build_hs_db(models_dict):
    """
    models_dict: mapping model_name -> regex string (like your `models` variable).
    Returns: tuple (db, id_to_model)
    """
    if not HYPERSCAN_AVAILABLE:
        return (None, {})

    expressions = []
    ids = []
    flags = []
    id_to_model = {}

    # assign integer ids 1..N for Hyperscan patterns
    for idx, (model_name, pattern) in enumerate(models_dict.items(), start=1):
        # Hyperscan expects bytes pattern. Use pattern as-is.
        expressions.append(pattern.encode('utf-8'))
        ids.append(idx)
        # HS_FLAG_UTF8 is fine; we don't need HS_FLAG_CASELESS for DNA upper case.
        flags.append(hyperscan.HS_FLAG_UTF8)
        id_to_model[idx] = model_name

    db = hyperscan.Database()
    try:
        db.compile(expressions=expressions, ids=ids, flags=flags)
    except Exception as e:
        # compile error: fallback
        print("Hyperscan compile error:", e)
        return (None, {})

    return (db, id_to_model)


def riz_search(seq, model_name):
    """
    Fallback regex-based search for R-loop sequences.
    Returns: list of dicts with match information
    """
    matches = []
    
    # Define the patterns for RLFS models
    patterns = {
        'm1': r'G{3,}[ATGC]{1,10}?G{3,}(?:[ATGC]{1,10}?G{3,}){1,}',
        'm2': r'G{4,}(?:[ATGC]{1,10}?G{4,}){1,}'
    }
    
    if model_name not in patterns:
        return matches
    
    pattern = patterns[model_name]
    for match in re.finditer(pattern, seq, re.IGNORECASE):
        start = match.start()
        end = match.end()
        riz_seq = seq[start:end]
        perc_g = percent_g(riz_seq)
        
        if perc_g >= params['min_perc_g_riz']:
            matches.append({
                "start": start,
                "end": end,
                "length": end - start,
                "G": riz_seq.count("G"),
                "G3s": riz_seq.count("GGG"),
                "G4s": riz_seq.count("GGGG"),
                "perc_g": perc_g,
                "model": model_name,
                "seq": riz_seq
            })
    
    return matches


def riz_search_hs(seq, model_name, hs_db, hs_id_map):
    """
    seq: DNA string (utf-8)
    model_name: 'm1' or 'm2' (string)
    hs_db: compiled Hyperscan Database (or None)
    hs_id_map: mapping int id -> model_name
    Returns: list of dicts same format as original riz_search
    """
    # fallback to regex if Hyperscan not available or db not built
    if not HYPERSCAN_AVAILABLE or hs_db is None:
        return riz_search(seq, model_name)  # call original re-based routine

    matches = []

    # local handler collects matches: Hyperscan will call this for every match
    def on_match(id, from_, to, flags, context):
        # id is integer pattern id, map to model_name
        matched_model = hs_id_map.get(id)
        # we only care matches for the requested model_name
        if matched_model != model_name:
            return 0  # continue
        # from_ and to are offsets (Python ints)
        start = int(from_)
        end = int(to)
        # extract substring (string slicing is fine; seq is str)
        riz_seq = seq[start:end]
        perc_g = percent_g(riz_seq)
        if perc_g >= params['min_perc_g_riz']:
            matches.append({
                "start": start,
                "end": end,
                "length": end - start,
                "G": riz_seq.count("G"),
                "G3s": riz_seq.count("GGG"),
                "G4s": riz_seq.count("GGGG"),
                "perc_g": perc_g,
                "model": model_name,
                "seq": riz_seq
            })
        return 0  # continue scanning

    # hyperscan requires bytes
    seq_bytes = seq.encode('ascii')  # DNA strings are ascii-safe
    try:
        hs_db.scan(seq_bytes, match_event_handler=on_match)
    except Exception as e:
        # an error at scan time — fallback to regex
        # (print once; avoid noisy logs in large runs)
        print("Hyperscan scan error, falling back to regex:", e)
        return riz_search(seq, model_name)

    # matches list might be unsorted depending on match callbacks; sort by start
    matches.sort(key=lambda x: x['start'])
    return matches


# Initialize regex models for R-loop detection
# Usage: after you create params['regex_models'], call:
# hs_db, hs_id_map = build_hs_db({k: v.pattern for k,v in params['regex_models'].items()})
# store hs_db and hs_id_map globally or pass into functions

# Initialize R-loop models  
params['regex_models'] = {
    'm1': type('Pattern', (), {'pattern': r'G{3,}[ATGC]{1,10}?G{3,}(?:[ATGC]{1,10}?G{3,}){1,}'})(),
    'm2': type('Pattern', (), {'pattern': r'G{4,}(?:[ATGC]{1,10}?G{4,}){1,}'})()
}

# Build Hyperscan database for R-loop detection
hs_db, hs_id_map = build_hs_db({k: v.pattern for k, v in params['regex_models'].items()})


class MotifBase(ABC):
    """
    Base class for all motif detectors.
    
    This abstract base class provides the common interface and functionality
    for all Non-B DNA motif detectors. Each specific detector inherits from
    this class and implements the detect() and score() methods.
    
    Attributes:
        class_name (str): Name of the motif class (e.g., 'g_quadruplex')
        class_id (int): Numeric identifier for the motif class
        patterns (list): List of regex patterns for this motif type
    
    Methods:
        detect(): Abstract method to detect motif candidates in a sequence
        score(): Abstract method to score detected candidates
        make_candidate(): Helper method to create Candidate objects
    """
    
    def __init__(self, class_name: str):
        """
        Initialize the motif detector.
        
        Args:
            class_name (str): Name of the motif class to detect
        """
        self.class_name = class_name
        self.class_id = get_class_id(class_name)
        self.patterns = get_patterns_for_motif(class_name)
        
    def make_candidate(self, seq: str, seq_name: str, contig: str, offset: int,
                      start: int, end: int, pattern_name: str, subclass: str,
                      motif_id: int) -> Candidate:
        """
        Create a candidate from match information.
        
        This helper method constructs a Candidate object from detection results,
        ensuring all required fields are properly set and derived fields are calculated.
        
        Args:
            seq (str): Full DNA sequence being analyzed
            seq_name (str): Name/identifier of the sequence
            contig (str): Chromosome or contig identifier
            offset (int): Offset for chunked sequences
            start (int): Start position of the motif (0-based)
            end (int): End position of the motif (0-based, inclusive)
            pattern_name (str): Name of the pattern that matched
            subclass (str): Subclass of the motif (e.g., 'canonical_G4')
            motif_id (int): Unique identifier for this motif instance
            
        Returns:
            Candidate: A fully initialized Candidate object
        """
        matched_seq = seq[start:end+1].encode('utf-8')
        
        return Candidate(
            sequence_name=seq_name,
            contig=contig,
            class_id=self.class_id,
            class_name=self.class_name,
            subclass=subclass,
            motif_id=motif_id,
            start=offset + start + 1,  # Convert to 1-based
            end=offset + end + 1,      # Convert to 1-based
            length=end - start + 1,
            matched_seq=matched_seq,
            pattern_name=pattern_name
        )
    
    @abstractmethod
    def detect(self, seq: str, seq_name: str, contig: str, offset: int) -> List[Candidate]:
        """
        Detect motif candidates in a DNA sequence.
        
        This abstract method must be implemented by each specific detector to
        identify potential motif sites using appropriate algorithms and patterns.
        
        Args:
            seq (str): DNA sequence to analyze (uppercase ATCG)
            seq_name (str): Identifier for the sequence
            contig (str): Chromosome or contig name
            offset (int): Position offset for chunked sequences
            
        Returns:
            List[Candidate]: List of detected motif candidates
        """
        pass
    
    @abstractmethod
    def score(self, candidates: List[Candidate]) -> List[Candidate]:
        """
        Score motif candidates using class-specific algorithms.
        
        This abstract method must be implemented by each detector to assign
        quality scores to detected candidates based on biological relevance
        and structural characteristics.
        
        Args:
            candidates (List[Candidate]): List of candidates to score
            
        Returns:
            List[Candidate]: Same candidates with scores assigned
        """
        pass


class G4Detector(MotifBase):
    """G-quadruplex family detector using G4Hunter scoring"""
    
    def __init__(self):
        super().__init__('g_quadruplex')
        self.db_key = 'g4_patterns'
        self._compile_patterns()
    
    def _compile_patterns(self):
        """Compile G4 patterns for Hyperscan"""
        all_patterns = []
        for subclass, patterns in self.patterns.items():
            for pattern_info in patterns:
                all_patterns.append(pattern_info[0])  # Just the regex string
        
        if all_patterns:
            hs_registry_manager.compile_db(all_patterns, self.db_key)
    
    def detect(self, seq: str, seq_name: str, contig: str, offset: int) -> List[Candidate]:
        """Detect G4 candidates using Hyperscan + regex validation"""
        candidates = []
        
        if not hs_registry_manager.is_available():
            return self._fallback_detect(seq, seq_name, contig, offset)
        
        try:
            # Get Hyperscan matches
            matches = hs_registry_manager.scan_block(seq, self.db_key)
            
            # Validate each match with regex to capture groups
            pattern_lookup = self._build_pattern_lookup()
            
            for pat_idx, start, end in matches:
                if pat_idx < len(pattern_lookup):
                    pattern_info = pattern_lookup[pat_idx]
                    pattern, subclass, motif_id = pattern_info
                    
                    # Use regex to capture the exact match and groups
                    regex_match = re.search(pattern, seq[start:end+10])  # Small buffer for groups
                    if regex_match:
                        actual_start = start + regex_match.start()
                        actual_end = start + regex_match.end() - 1
                        
                        candidate = self.make_candidate(
                            seq, seq_name, contig, offset,
                            actual_start, actual_end, pattern,
                            subclass, motif_id
                        )
                        candidates.append(candidate)
        
        except Exception as e:
            logger.warning(f"Hyperscan G4 detection failed: {e}, using fallback")
            return self._fallback_detect(seq, seq_name, contig, offset)
        
        return candidates
    
    def _build_pattern_lookup(self) -> List[Tuple[str, str, int]]:
        """Build lookup table for pattern index -> (pattern, subclass, motif_id)"""
        lookup = []
        for subclass, patterns in self.patterns.items():
            for pattern_info in patterns:
                pattern = pattern_info[0]
                orig_id = pattern_info[1]
                subclass_name = pattern_info[3] if len(pattern_info) > 3 else subclass
                lookup.append((pattern, subclass_name, orig_id))
        return lookup
    
    def _fallback_detect(self, seq: str, seq_name: str, contig: str, offset: int) -> List[Candidate]:
        """Fallback detection using pure Python regex"""
        candidates = []
        for subclass, patterns in self.patterns.items():
            for pattern_info in patterns:
                pattern = pattern_info[0]
                orig_id = pattern_info[1]
                subclass_name = pattern_info[3] if len(pattern_info) > 3 else subclass
                try:
                    for match in re.finditer(pattern, seq, re.IGNORECASE):
                        candidate = self.make_candidate(
                            seq, seq_name, contig, offset,
                            match.start(), match.end() - 1, pattern,
                            subclass_name, orig_id
                        )
                        candidates.append(candidate)
                except Exception as e:
                    logger.warning(f"Regex failed for pattern {pattern}: {e}")
        return candidates
    
    def score(self, candidates: List[Candidate]) -> List[Candidate]:
        """Score G4 candidates using G4Hunter algorithm"""
        for candidate in candidates:
            seq_str = candidate.matched_seq.decode('utf-8')
            candidate.raw_score = g4hunter_score(seq_str)
            candidate.scoring_method = "G4Hunter_v2"
        return candidates


class CurvedDetector(MotifBase):
    """Curved DNA detector for A-phased repeats"""
    
    def __init__(self):
        super().__init__('curved_dna')
    
    def detect(self, seq: str, seq_name: str, contig: str, offset: int) -> List[Candidate]:
        """Detect curved DNA using A-tract phasing algorithm"""
        candidates = []
        
        # Look for A-tracts with 8-12 bp spacing (phased arrays)
        a_tract_pattern = r'A{3,9}'
        a_tracts = list(re.finditer(a_tract_pattern, seq, re.IGNORECASE))
        
        # Find phased arrays (A-tracts separated by ~10bp)
        for i in range(len(a_tracts)):
            for j in range(i + 1, min(i + 6, len(a_tracts))):  # Look ahead max 5 A-tracts
                tract1 = a_tracts[i]
                tract2 = a_tracts[j]
                
                spacing = tract2.start() - tract1.end()
                if 7 <= spacing <= 13:  # Phased spacing
                    # Create candidate spanning both tracts
                    candidate = self.make_candidate(
                        seq, seq_name, contig, offset,
                        tract1.start(), tract2.end() - 1, 
                        "A-tract_phased", "Phased_Array", 1
                    )
                    candidates.append(candidate)
        
        return candidates
    
    def score(self, candidates: List[Candidate]) -> List[Candidate]:
        """Score curved DNA candidates"""
        for candidate in candidates:
            seq_str = candidate.matched_seq.decode('utf-8')
            candidate.raw_score = curvature_score(seq_str)
            candidate.scoring_method = "Curvature_propensity"
        return candidates


class SlippedDetector(MotifBase):
    """Slipped DNA detector for STRs and direct repeats"""
    
    def __init__(self):
        super().__init__('slipped_dna')
    
    def detect(self, seq: str, seq_name: str, contig: str, offset: int) -> List[Candidate]:
        """Detect slipped DNA structures (STRs, direct repeats)"""
        candidates = []
        
        # Use simple repetitive patterns since backreferences aren't Hyperscan-safe
        patterns = [
            (r'([ATGC])\1{9,}', 'Mononucleotide_STR'),  # Simple runs
            (r'([ATGC]{2})\1{7,}', 'Dinucleotide_STR'),  # Note: this uses backreference
            (r'([ATGC]{3})\1{5,}', 'Trinucleotide_STR'),  # Note: this uses backreference
        ]
        
        for pattern, subclass in patterns:
            try:
                for match in re.finditer(pattern, seq, re.IGNORECASE):
                    candidate = self.make_candidate(
                        seq, seq_name, contig, offset,
                        match.start(), match.end() - 1, pattern,
                        subclass, 1
                    )
                    candidates.append(candidate)
            except Exception as e:
                logger.warning(f"STR pattern {pattern} failed: {e}")
        
        return candidates
    
    def score(self, candidates: List[Candidate]) -> List[Candidate]:
        """Score slipped DNA candidates"""
        for candidate in candidates:
            # Simple repetitiveness score
            seq_str = candidate.matched_seq.decode('utf-8')
            candidate.raw_score = len(seq_str) / 100.0  # Length-based scoring
            candidate.scoring_method = "STR_length"
        return candidates


class CruciformDetector(MotifBase):
    """Cruciform detector for palindromes and inverted repeats"""
    
    def __init__(self):
        super().__init__('cruciform')
    
    def detect(self, seq: str, seq_name: str, contig: str, offset: int) -> List[Candidate]:
        """Detect cruciform structures using palindrome detection"""
        candidates = []
        
        # Look for palindromic sequences (simplified approach)
        min_arm = 6
        max_arm = 20
        max_loop = 10
        
        for arm_len in range(min_arm, max_arm):
            for loop_len in range(0, max_loop):
                for i in range(len(seq) - 2 * arm_len - loop_len):
                    left_arm = seq[i:i + arm_len]
                    loop = seq[i + arm_len:i + arm_len + loop_len]
                    right_arm = seq[i + arm_len + loop_len:i + 2 * arm_len + loop_len]
                    
                    # Check if right arm is reverse complement of left arm
                    if self._is_reverse_complement(left_arm, right_arm):
                        candidate = self.make_candidate(
                            seq, seq_name, contig, offset,
                            i, i + 2 * arm_len + loop_len - 1,
                            f"palindrome_{arm_len}_{loop_len}", "Palindrome", 1
                        )
                        candidates.append(candidate)
        
        return candidates
    
    def _is_reverse_complement(self, seq1: str, seq2: str) -> bool:
        """Check if seq2 is reverse complement of seq1"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                     'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
        
        if len(seq1) != len(seq2):
            return False
        
        rev_comp = ''.join(complement.get(c, c) for c in reversed(seq1))
        return rev_comp.upper() == seq2.upper()
    
    def score(self, candidates: List[Candidate]) -> List[Candidate]:
        """Score cruciform candidates"""
        for candidate in candidates:
            # Score based on arm length and symmetry
            candidate.raw_score = candidate.length / 100.0
            candidate.scoring_method = "Palindrome_stability"
        return candidates


class RLoopDetector(MotifBase):
    """R-loop detector using RLFS models with Hyperscan optimization"""
    
    def __init__(self):
        super().__init__('r_loop')
    
    def detect(self, seq: str, seq_name: str, contig: str, offset: int) -> List[Candidate]:
        """Detect R-loop forming sequences using hyperscan-optimized search"""
        candidates = []
        
        # Use hyperscan-optimized search for both m1 and m2 models
        for model_name in ['m1', 'm2']:
            matches = riz_search_hs(seq, model_name, hs_db, hs_id_map)
            
            for match in matches:
                # Map model names to subclasses
                subclass = 'RLFS_m1' if model_name == 'm1' else 'RLFS_m2'
                
                candidate = self.make_candidate(
                    seq, seq_name, contig, offset,
                    match['start'], match['end'] - 1, 
                    params['regex_models'][model_name].pattern,
                    subclass, 1
                )
                
                # Store additional R-loop specific data
                candidate.g_count = match['G']
                candidate.g3_count = match['G3s'] 
                candidate.g4_count = match['G4s']
                candidate.perc_g = match['perc_g']
                
                candidates.append(candidate)
        
        return candidates
    
    def score(self, candidates: List[Candidate]) -> List[Candidate]:
        """Score R-loop candidates using QmRLFS-inspired scoring with G-content"""
        for candidate in candidates:
            seq_str = candidate.matched_seq.decode('utf-8')
            
            # Use the pre-calculated G percentage if available, otherwise calculate
            if hasattr(candidate, 'perc_g') and candidate.perc_g is not None:
                g_content = candidate.perc_g
            else:
                g_content = percent_g(seq_str)
            
            # Enhanced scoring that considers G-richness and pattern complexity
            base_score = g_content
            
            # Bonus for higher G-tract density
            if hasattr(candidate, 'g4_count') and candidate.g4_count > 0:
                base_score += candidate.g4_count * 0.1
            elif hasattr(candidate, 'g3_count') and candidate.g3_count > 0:
                base_score += candidate.g3_count * 0.05
            
            candidate.raw_score = min(base_score, 1.0)  # Cap at 1.0
            candidate.scoring_method = "QmRLFS_hyperscan_enhanced"
            
        return candidates


class TriplexDetector(MotifBase):
    """Triplex DNA detector for homopurine/homopyrimidine tracts"""
    
    def __init__(self):
        super().__init__('triplex')
    
    def detect(self, seq: str, seq_name: str, contig: str, offset: int) -> List[Candidate]:
        """Detect triplex-forming sequences"""
        candidates = []
        
        # Homopurine and homopyrimidine tracts
        patterns = [
            (r'[AG]{15,}', 'Homopurine'),
            (r'[CT]{15,}', 'Homopyrimidine'),
        ]
        
        for pattern, subclass in patterns:
            for match in re.finditer(pattern, seq, re.IGNORECASE):
                candidate = self.make_candidate(
                    seq, seq_name, contig, offset,
                    match.start(), match.end() - 1, pattern,
                    subclass, 1
                )
                candidates.append(candidate)
        
        return candidates
    
    def score(self, candidates: List[Candidate]) -> List[Candidate]:
        """Score triplex candidates"""
        for candidate in candidates:
            seq_str = candidate.matched_seq.decode('utf-8')
            candidate.raw_score = triplex_stability_score(seq_str)
            candidate.scoring_method = "Triplex_stability"
        return candidates


class IMotifDetector(MotifBase):
    """i-Motif detector for C-rich quadruplex structures"""
    
    def __init__(self):
        super().__init__('i_motif')
    
    def detect(self, seq: str, seq_name: str, contig: str, offset: int) -> List[Candidate]:
        """Detect i-motif structures"""
        candidates = []
        
        # C-rich sequences that can form i-motifs
        pattern = r'C{3,}[ATGC]{1,12}C{3,}[ATGC]{1,12}C{3,}[ATGC]{1,12}C{3,}'
        
        for match in re.finditer(pattern, seq, re.IGNORECASE):
            candidate = self.make_candidate(
                seq, seq_name, contig, offset,
                match.start(), match.end() - 1, pattern,
                'Canonical_iMotif', 1
            )
            candidates.append(candidate)
        
        return candidates
    
    def score(self, candidates: List[Candidate]) -> List[Candidate]:
        """Score i-motif candidates using adapted G4Hunter"""
        for candidate in candidates:
            seq_str = candidate.matched_seq.decode('utf-8')
            # Adapt G4Hunter for C-rich sequences
            c_runs = seq_str.upper().count('CCC')
            g_runs = seq_str.upper().count('GGG')
            candidate.raw_score = (c_runs - g_runs) / len(seq_str)
            candidate.scoring_method = "G4Hunter_adapted_iMotif"
        return candidates


class ZDNADetector(MotifBase):
    """Z-DNA detector using Kadane algorithm for maximum subarray detection"""
    
    def __init__(self, use_kadane=True):
        super().__init__('z_dna')
        self.use_kadane = use_kadane
        
    def detect(self, seq: str, seq_name: str, contig: str, offset: int) -> List[Candidate]:
        """Detect Z-DNA forming sequences"""
        candidates = []
        
        if self.use_kadane:
            # Use the new Kadane algorithm-based detection
            try:
                from zdna_calculator import ZDNACalculatorSeq
                from constants import Params
                params = Params(threshold=5.0)
                calculator = ZDNACalculatorSeq(seq, params)
                subarrays = calculator.subarrays_above_threshold()
                
                for start, end, score, substring in subarrays:
                    candidate = self.make_candidate(
                        seq, seq_name, contig, offset,
                        start, end - 1, "kadane_zdna",
                        "Z_DNA_Kadane", 1
                    )
                    candidate.raw_score = score
                    candidate.scoring_method = "Z_DNA_Kadane_algorithm"
                    candidates.append(candidate)
                    
            except ImportError:
                # Fall back to regex if zdna_calculator module not available
                self.use_kadane = False
                
        if not self.use_kadane:
            # Original regex-based patterns
            patterns = [
                (r'([CG]{2}){6,}', 'Z_DNA_basic'),
                (r'G[CG]{8,}G', 'Extended_GZ'),
            ]
            
            for pattern, subclass in patterns:
                for match in re.finditer(pattern, seq, re.IGNORECASE):
                    candidate = self.make_candidate(
                        seq, seq_name, contig, offset,
                        match.start(), match.end() - 1, pattern,
                        subclass, 1
                    )
                    candidates.append(candidate)
        
        return candidates
    
    def score(self, candidates: List[Candidate]) -> List[Candidate]:
        """Score Z-DNA candidates"""
        for candidate in candidates:
            # Skip scoring if already scored by Kadane method
            if candidate.raw_score is not None:
                continue
                
            seq_str = candidate.matched_seq.decode('utf-8')
            candidate.raw_score = z_dna_score(seq_str)
            candidate.scoring_method = "Z_seeker_adapted"
        return candidates


class APhilicDetector(MotifBase):
    """
    Enhanced A-philic DNA detector using sophisticated tetranucleotide and trinucleotide scoring.
    
    Implements the advanced algorithm from call_Aphilic.py with:
    - Comprehensive tetra/tri propensity tables (embedded)
    - Nucleation detection (10-nt windows with all positive tetra steps + tri-window cutoff)
    - Extension with Kadane-variant to find best containing subarray
    - Hyperscan integration for performance optimization
    """
    
    def __init__(self):
        super().__init__('a_philic')
        
        # Comprehensive propensity tables from call_Aphilic.py
        # Tetra table: Step -> Log2_Odds_laplace
        self.TETRA_LOG2 = {
            "CCCC":4.389556283101704,"GGGG":4.389556283101704,"TGGG":4.167163861765255,"GGGC":3.9041294559314617,
            "CCCG":3.9041294559314617,"GCCC":3.9041294559314617,"CCCT":3.582201361044099,"GTGC":3.582201361044099,
            "AGGG":3.582201361044099,"TCCC":3.167163861765255,"CCCA":3.167163861044099,"CCTA":2.582201361044099,
            "TAGG":2.582201361044099,"CTCC":2.582201361044099,"CGGG":2.3191669552103056,"GAGG":2.1671638617652556,
            "GGGT":2.1671638617652556,"GCAC":1.9972388603229432,"CCAC":1.9041294559314614,"CCGG":1.8817616429030073,
            "GGCC":1.7077322431279585,"CCTC":1.5822013610440995,"TCCT":1.5822013610440995,"GACC":1.5822013610440993,
            "CTGT":1.5822013610440993,"CTCA":1.5822013610440993,"CCGC":1.5822013610440993,"TGCC":1.5822013610440993,
            "TAAG":1.5822013610440993,"TACC":1.3191669552103056,"TCGG":1.1671638617652556,"CTAG":1.0967745338738575,
            "GTGG":1.0967745338738575,"GTCC":0.9972388603229432,"CACG":0.9972388603229432,"GGTC":0.9972388603229432,
            "GGTA":0.9041294559314615,"GTAC":0.8452357668778931,"TACG":0.8452357668778931,"GGGA":0.8452357668778931,
            "ACGC":0.7748464389864953,"GCGG":0.7748464389864953,"CGGC":0.7342044544891495,"CGGT":0.7077322431279583,
            "ACGT":0.7077322431279583,"CGTA":0.5822013610440996,"TCTC":0.5822013610440996,"CCGA":0.5822013610440996,
            "GCCG":0.5822013610440994,"ACCG":0.5822013610440994,"TCCA":0.5822013610440991,"CAGT":0.5822013610440991,
            "TCCG":0.5822013610440991,"CACA":0.5822013610440991,"TCAG":0.5822013610440991,"CACT":0.5822013610440991,
            "TCAA":0.5822013610440991,"GGTG":0.5822013610440991,"GTAA":0.5822013610440991,"GGAG":0.5822013610440991,
            "TTGA":0.5822013610440991,"GGCT":0.5822013610440991,"GCAG":0.5822013610440991,"TTAC":0.5822013610440991,
            "TGTT":0.5822013610440991,"ACAG":0.5822013610440991,"TGTG":0.5822013610440991,"ACCC":0.5822013610440991,
            "ATCC":0.5822013610440991,"CTTA":0.5822013610440991,"ACTC":0.5822013610440991,"AGCC":0.5822013610440991,
            "AGTC":0.5822013610440991,"AGTG":0.5822013610440991,"ATAC":0.5822013610440991,"CGTG":0.41227635960178693,
            "TGCG":0.3598089397076514,"GCGC":0.3598089397076514,"GTAT":0.26027326615673674,"GTCT":0.26027326615673674,
            "GTGT":0.26027326615673674,"GCCT":0.26027326615673674,"TACA":0.26027326615673674,"GGCA":0.26027326615673674,
            "AGGC":0.26027326615673674,"CACC":0.26027326615673674,"ACAC":0.26027326615673674,"TCTG":0.26027326615673674,
            "TGAC":0.26027326615673674,"CGCA":0.1671638617652555,"GCGT":0.09677453387385747,"CATG":0.09677453387385747,
            "CAGA":-0.00276114,"ACTG":-0.00276114,"ATCA":-0.00276114,"TGCA":-0.00276114,"TGTA":-0.00276114,"CTAC":-0.00276114,
            "TGGC":-0.00276114,"GGTT":-0.00276114,"TTTA":-0.00276114,"AGTA":-0.00276114,"TAAA":-0.00276114,"GTTG":-0.00276114,
            "AGGA":-0.00276114,"CTGC":-0.00276114,"TGTC":-0.00276114,"TCAC":-0.00276114,"GATC":-0.00276114,"AACC":-0.00276114,
            "ATGG":-0.00276114,"ACCT":-0.00276114,"AGGT":-0.00276114,"TACT":-0.00276114,"TTAG":-0.00276114,"TGAA":-0.00276114,
            "AAGT":-0.00276114,"TAGT":-0.00276114,"AACT":-0.00276114,"TATT":-0.00276114,"GGAC":-0.00276114,"CAAC":-0.00276114,
            "ATGC":-0.154764233,"CCAT":-0.225153561,"CGCC":-0.225153561,"GTTC":-0.417798639,"AAGC":-0.417798639,
            "CTGA":-0.417798639,"AATC":-0.417798639,"AATA":-0.417798639,"CGAC":-0.417798639,"AAGG":-0.417798639,
            "CCTT":-0.417798639,"ACGA":-0.417798639,"TTCA":-0.417798639,"GCTC":-0.417798639,"AACA":-0.417798639,
            "GGCG":-0.417798639,"TCAT":-0.417798639,"GGAT":-0.417798639,"ATGT":-0.417798639,"ACCA":-0.417798639,
            "ATGA":-0.417798639,"CTTT":-0.417798639,"AGCA":-0.417798639,"CTAT":-0.417798639,"GTCG":-0.417798639,
            "GACA":-0.417798639,"TGAG":-0.417798639,"AGAC":-0.417798639,"TGGT":-0.417798639,"ACTT":-0.417798639,
            "ATAA":-0.417798639,"GAGC":-0.417798639,"AGTT":-0.417798639,"TAAC":-0.417798639,"TGAT":-0.417798639,
            "TGCT":-0.417798639,"GCAT":-0.533275856,"CCGT":-0.58772364,"TTGG":-0.739726734,"TTAT":-0.739726734,
            "TCGC":-0.739726734,"TAGA":-0.739726734,"CTTG":-0.739726734,"TTGT":-0.739726734,"GCTT":-0.739726734,
            "AGCG":-0.739726734,"AAAG":-0.739726734,"ACAT":-0.739726734,"CAGC":-0.739726734,"GAAC":-0.739726734,
            "CATC":-0.739726734,"CATT":-0.739726734,"GAGT":-0.739726734,"CGGA":-0.739726734,"ATCT":-0.739726734,
            "CCTG":-0.739726734,"ACTA":-0.739726734,"AGAT":-1.00276114,"AATG":-1.00276114,"CTAA":-1.00276114,
            "CGAG":-1.00276114,"CCAG":-1.00276114,"CTCT":-1.00276114,"CATA":-1.00276114,"CAAG":-1.00276114,
            "CTCG":-1.00276114,"TCGT":-1.00276114,"TTGC":-1.00276114,"AGAG":-1.00276114,"GATT":-1.00276114,
            "GATG":-1.00276114,"ATAG":-1.00276114,"GACT":-1.00276114,"GTGA":-1.00276114,"GTTA":-1.00276114,
            "TTCT":-1.00276114,"CGCG":-1.080763652,"ATCG":-1.118238357,"TATA":-1.118238357,"TTTG":-1.225153561,
            "TTCC":-1.225153561,"ACGG":-1.225153561,"AGCT":-1.225153561,"GCCA":-1.225153561,"ACAA":-1.225153561,
            "TCTT":-1.225153561,"CAGG":-1.225153561,"TCTA":-1.225153561,"AAGA":-1.225153561,"CGAT":-1.225153561,
            "GTAG":-1.225153561,"TATG":-1.225153561,"GTCA":-1.225153561,"CGCT":-1.225153561,"AGAA":-1.225153561,
            "CTGG":-1.225153561,"TTTC":-1.417798639,"ATTA":-1.417798639,"CCAA":-1.417798639,"TATC":-1.417798639,
            "GAAG":-1.417798639,"GACG":-1.417798639,"GAGA":-1.417798639,"GCTG":-1.417798639,"TGGA":-1.417798639,
            "TTAA":-1.58772364,"GTTT":-1.58772364,"CTTC":-1.58772364,"GATA":-1.58772364,"GCTA":-1.58772364,
            "TTCG":-1.739726734,"GAAA":-1.739726734,"TCGA":-1.739726734,"CAAT":-1.739726734,"TAGC":-1.739726734,
            "AAAC":-1.739726734,"TAAT":-1.877230258,"CGTC":-1.877230258,"ATTT":-2.00276114,"CGAA":-2.00276114,
            "ATTG":-2.00276114,"AACG":-2.00276114,"GCAA":-2.00276114,"AAAT":-2.118238357,"CAAA":-2.118238357,
            "GCGA":-2.118238357,"ATTC":-2.225153561,"GAAT":-2.324689235,"CGTT":-2.324689235,"GGAA":-2.324689235,
            "AAAA":-2.417798639,"ATAT":-2.417798639,"TTTT":-2.50526148,"AATT":-3.50526148
        }

        # Tri table: Step -> Log2_Odds_laplace  
        self.TRI_LOG2 = {
            "CCC":4.781079142726248,"GGG":3.9737242206686436,"CAC":1.6656019253063112,"GCC":1.557077468528142,
            "GGC":1.557077468528142,"CCG":1.4526082019721132,"GTG":1.2505644260274673,"ACC":1.2326425180302052,
            "CCT":1.1737488289766367,"GGT":1.1546400060289317,"CGG":1.0806394245851554,"AGG":1.080639424585155,
            "TAC":0.9811037510342406,"TCC":0.8582470032487074,"GTA":0.7810791427262473,"CTC":0.5952125974149134,
            "TGC":0.5660662517553969,"CCA":0.303031845921603,"GTC":0.2732845025275512,"TGG":0.1875546285016671,
            "CTA":0.1420399692492986,"TAG":0.08063942458515531,"ACG":0.04111106039851759,"GCA":0.030013351515187126,
            "GAC":0.010250096693757146,"GCG":-0.023697235,"CGT":-0.074638801,"CGC":-0.141752997,"CAT":-0.551628791,
            "TCT":-0.563216765,"ATG":-0.597432481,"TGT":-0.619800294,"GAG":-0.65632617,"CAG":-0.726715497,
            "TGA":-0.726715497,"GGA":-0.873556886,"ATC":-0.919360575,"ACA":-1.006823417,"CTG":-1.089285577,
            "GAT":-1.24128867,"AGT":-1.24128867,"ACT":-1.378792194,"TCA":-1.504323076,"TCG":-1.54385144,
            "CGA":-1.619800294,"TAT":-1.726715497,"ATA":-1.777341571,"TTA":-2.006823417,"AAG":-2.006823417,
            "TAA":-2.089285577,"CTT":-2.167288089,"AGC":-2.24128867,"AGA":-2.311677998,"AAC":-2.378792194,
            "GTT":-2.378792194,"TTG":-2.504323076,"GCT":-2.504323076,"CAA":-2.777341571,"TTC":-2.919360575,
            "GAA":-3.128813941,"TTT":-3.204762794,"AAA":-3.311677998,"ATT":-3.411213672,"AAT":-3.473949427
        }
        
        # Initialize hyperscan support
        self.hs_db = None
        self.id_to_tetra = {}
        self._compile_hyperscan_db()
    
    def _compile_hyperscan_db(self):
        """Compile Hyperscan database for positive tetranucleotides (performance optimization)"""
        if not HYPERSCAN_AVAILABLE:
            logger.info("Hyperscan not available, using fallback scanning")
            return
            
        try:
            # Get positive tetranucleotides for hyperscan patterns
            positive_tetras = [tetra for tetra, score in self.TETRA_LOG2.items() if score > 0.0]
            
            if not positive_tetras:
                return
                
            expressions = []
            ids = []
            flags = []
            
            for i, tetra in enumerate(positive_tetras):
                expressions.append(tetra.encode('utf-8'))
                ids.append(i)
                flags.append(hyperscan.HS_FLAG_UTF8)
                self.id_to_tetra[i] = tetra
            
            db = hyperscan.Database()
            db.compile(expressions=expressions, ids=ids, flags=flags)
            self.hs_db = db
            logger.info(f"Compiled Hyperscan DB for {len(positive_tetras)} A-philic tetranucleotides")
            
        except Exception as e:
            logger.warning(f"Failed to compile Hyperscan DB for A-philic detector: {e}")
            self.hs_db = None

    def _build_step_scores(self, seq: str, w4: float = 0.7, w3: float = 0.3):
        """
        Build tetra, tri, and combined step score arrays as in call_Aphilic.py
        
        Returns:
            tetra_scores: array of tetra scores (length seq_len - 3)
            tri_scores: array of tri scores (length seq_len - 2) 
            step_scores: combined scores (length seq_len - 3)
        """
        import numpy as np
        
        L = len(seq)
        if L < 4:
            return np.array([]), np.array([]), np.array([])
            
        tetra_scores = np.zeros(L - 3, dtype=float)
        tri_scores = np.zeros(L - 2, dtype=float)
        
        # Calculate tetra scores
        for i in range(L - 3):
            tetra_scores[i] = self.TETRA_LOG2.get(seq[i:i+4], 0.0)
            
        # Calculate tri scores  
        for i in range(L - 2):
            tri_scores[i] = self.TRI_LOG2.get(seq[i:i+3], 0.0)
            
        # Combined step scores (align tri at same start as tetra)
        step_scores = np.zeros(L - 3, dtype=float)
        for i in range(L - 3):
            tri_val = tri_scores[i] if i < len(tri_scores) else 0.0
            tet_val = tetra_scores[i]
            step_scores[i] = w4 * tet_val + w3 * tri_val
            
        return tetra_scores, tri_scores, step_scores

    def _find_10mer_positive_tetra_starts(self, tetra_scores):
        """
        Find 10-nt windows where all 7 tetra scores are > 0 (nucleation detection)
        """
        import numpy as np
        
        starts = []
        min_needed = 7  # 10 nt -> 7 tetra steps
        n_tetra = len(tetra_scores)
        
        for j in range(0, n_tetra - min_needed + 1):
            window = tetra_scores[j:j + min_needed]
            if np.all(window > 0.0):
                starts.append(j)
                
        return starts

    def _compute_nuc_threshold_auto(self, tri_scores, factor: float = 1.0):
        """
        Auto nucleation threshold: mean + factor * std * sqrt(window_len)
        """
        import math
        import numpy as np
        
        if len(tri_scores) == 0:
            return 0.0
            
        # Use population statistics from TRI table
        all_tri_values = list(self.TRI_LOG2.values())
        mean = float(np.mean(all_tri_values))
        std = float(np.std(all_tri_values))
        
        # Use 3 tri-steps nucleation window
        return mean * 3.0 + factor * std * math.sqrt(3.0)

    def _tri_window_sum_ok(self, tri_scores, nt_start: int, window_tri_len: int, nuc_threshold: float):
        """
        Check if any consecutive tri-window of specified length has sum >= nuc_threshold
        """
        tri_start = nt_start
        max_tri_index = len(tri_scores) - 1
        tri_indices_in_10nt = list(range(tri_start, min(tri_start + 8, max_tri_index + 1)))
        
        if len(tri_indices_in_10nt) < window_tri_len:
            return False
            
        for a in range(len(tri_indices_in_10nt) - window_tri_len + 1):
            s = sum(tri_scores[tri_indices_in_10nt[a]: tri_indices_in_10nt[a] + window_tri_len])
            if s >= nuc_threshold:
                return True
                
        return False

    def _best_subarray_containing_interval(self, step_scores, a: int, b: int):
        """
        Find best subarray containing the interval [a,b] using prefix min/max optimization
        """
        import numpy as np
        
        n = len(step_scores)
        if n == 0:
            return (0, 0, 0.0)
            
        prefix = np.empty(n + 1, dtype=float)
        prefix[0] = 0.0
        
        for i in range(n):
            prefix[i + 1] = prefix[i] + float(step_scores[i])
            
        # Find minimal prefix up to each index
        min_pref_idx = np.zeros(n + 1, dtype=int)
        min_val = prefix[0]
        min_idx = 0
        
        for i in range(n + 1):
            if prefix[i] < min_val:
                min_val = prefix[i]
                min_idx = i
            min_pref_idx[i] = min_idx
            
        # Find maximal prefix from each index
        max_pref_idx = np.zeros(n + 1, dtype=int)
        max_val = prefix[-1]
        max_idx = n
        
        for i in range(n, -1, -1):
            if prefix[i] > max_val:
                max_val = prefix[i]
                max_idx = i
            max_pref_idx[i] = max_idx
            
        # Find best L and R
        L_idx = min_pref_idx[a]
        Rp1_idx = max_pref_idx[b + 1] if (b + 1) <= n else n
        
        best_sum = prefix[Rp1_idx] - prefix[L_idx]
        bestL = L_idx
        bestR = Rp1_idx - 1
        
        return (int(bestL), int(bestR), float(best_sum))

    def _select_non_overlapping_regions(self, regions):
        """
        Greedy selection of non-overlapping regions by score (descending)
        """
        chosen = []
        regions_sorted = sorted(regions, key=lambda r: r['score'], reverse=True)
        occupied = []
        
        for r in regions_sorted:
            bad = False
            for (a, b) in occupied:
                if not (r['end_nt'] < a or r['start_nt'] > b):
                    bad = True
                    break
            if not bad:
                chosen.append(r)
                occupied.append((r['start_nt'], r['end_nt']))
                
        return sorted(chosen, key=lambda x: x['start_nt'])

    def detect(self, seq: str, seq_name: str, contig: str, offset: int) -> List[Candidate]:
        """
        Enhanced A-philic detection using the sophisticated algorithm from call_Aphilic.py
        """
        if not seq or len(seq) < 10:
            return []
            
        # Clean sequence
        seq = seq.upper().replace(" ", "").replace("\n", "")
        
        # Check for non-ATGC characters
        accept_bases = set("ATGC")
        if any(ch not in accept_bases for ch in seq):
            # For simplicity, skip sequences with ambiguous bases
            return []
            
        # Build step scores using the sophisticated algorithm
        tetra_scores, tri_scores, step_scores = self._build_step_scores(seq)
        
        if len(step_scores) == 0:
            return []
            
        # Find 10-mer starts with all positive tetra steps (nucleation)
        ten_starts = self._find_10mer_positive_tetra_starts(tetra_scores)
        
        # Auto nucleation threshold
        nuc_threshold = self._compute_nuc_threshold_auto(tri_scores, factor=1.0)
        
        # Validate seeds with tri-window criteria
        valid_seeds = []
        nuc_tri_win = 3
        
        for j in ten_starts:
            if self._tri_window_sum_ok(tri_scores, j, nuc_tri_win, nuc_threshold):
                # Seed tetra step interval = j .. j+6
                valid_seeds.append((j, j + 6))
                
        if not valid_seeds:
            return []
            
        # Extend each seed to find best containing subarray
        all_regions = []
        
        for (a, b) in valid_seeds:
            L_step, R_step, best_sum = self._best_subarray_containing_interval(step_scores, a, b)
            
            # Convert step indices to nucleotide coordinates
            start_nt = L_step
            end_nt = R_step + 3  # step covers 4 nucleotides
            n_nt = end_nt - start_nt + 1
            n_steps = R_step - L_step + 1
            mean_step = best_sum / max(1, n_steps)
            
            # Apply acceptance criteria
            if n_nt >= 10 and mean_step > 0.0:
                # Count positive tetra and tri
                pos_tetra_count = int(np.sum(tetra_scores[L_step:R_step + 1] > 0.0))
                tri_L = L_step
                tri_R = min(len(tri_scores) - 1, R_step + 1)
                pos_tri_count = int(np.sum(tri_scores[tri_L:tri_R + 1] > 0.0))
                
                region_seq = seq[L_step:end_nt + 1]
                
                all_regions.append({
                    "start_nt": start_nt + offset,
                    "end_nt": end_nt + offset,
                    "sequence": region_seq,
                    "score": float(best_sum),
                    "mean_step_score": float(mean_step),
                    "n_nt": int(n_nt),
                    "n_steps": int(n_steps),
                    "pos_tetra": pos_tetra_count,
                    "pos_tri": pos_tri_count,
                    "seed_tetra_a": a + offset,
                    "seed_tetra_b": b + offset
                })
                
        # Select non-overlapping regions
        selected_regions = self._select_non_overlapping_regions(all_regions)
        
        # Convert to Candidate objects
        final_candidates = []
        for i, region in enumerate(selected_regions):
            candidate = self.make_candidate(
                seq=seq,
                seq_name=seq_name,
                contig=contig,
                offset=0,  # offset already applied to coordinates
                start=region["start_nt"] - offset,  # adjust back for make_candidate
                end=region["end_nt"] - offset,
                pattern_name="A_philic_enhanced", 
                subclass="A-philic DNA",
                motif_id=i
            )
            
            # Store enhanced scoring information
            candidate.raw_score = region["score"]
            candidate.metadata = {
                "mean_step_score": region["mean_step_score"],
                "n_nt": region["n_nt"],
                "n_steps": region["n_steps"],
                "pos_tetra": region["pos_tetra"],
                "pos_tri": region["pos_tri"],
                "algorithm": "enhanced_nucleation_extension"
            }
            
            final_candidates.append(candidate)
            
        return final_candidates
    
    def score(self, candidates: List[Candidate]) -> List[Candidate]:
        """Enhanced scoring using the sophisticated algorithm metrics"""
        for candidate in candidates:
            # Use the pre-calculated score from the enhanced detection
            if hasattr(candidate, 'raw_score') and candidate.raw_score is not None:
                # Normalize by length for fair comparison
                normalized_score = candidate.raw_score / max(1, candidate.length)
                candidate.raw_score = max(0.01, normalized_score)
            else:
                # Fallback scoring
                candidate.raw_score = candidate.length / 100.0
                
            candidate.scoring_method = "A_philic_enhanced_nucleation_extension"
            
        return candidates


class HybridDetector(MotifBase):
    """Hybrid detector for overlapping motif structures"""
    
    def __init__(self):
        super().__init__('hybrid')
    
    def detect(self, seq: str, seq_name: str, contig: str, offset: int) -> List[Candidate]:
        """Detect hybrid structures - this is called after primary detection"""
        # This detector works on the output of other detectors
        return []
    
    def detect_from_candidates(self, candidates: List[Candidate]) -> List[Candidate]:
        """Detect hybrid structures from existing candidates"""
        hybrid_candidates = []
        
        if not INTERVALTREE_AVAILABLE:
            return self._simple_overlap_detection(candidates)
        
        # Use interval tree for efficient overlap detection
        intervals = IntervalTree()
        for i, candidate in enumerate(candidates):
            intervals[candidate.start:candidate.end+1] = candidate
        
        # Find overlapping candidates
        overlaps = defaultdict(list)
        for candidate in candidates:
            overlapping = intervals[candidate.start:candidate.end+1]
            if len(overlapping) > 1:  # Has overlaps
                overlap_classes = [interval.data.class_name for interval in overlapping if interval.data != candidate]
                if overlap_classes:
                    overlaps[candidate.motif_id] = overlap_classes
        
        # Create hybrid candidates for significant overlaps
        for candidate in candidates:
            if candidate.motif_id in overlaps:
                overlap_classes = overlaps[candidate.motif_id]
                hybrid = Candidate(
                    sequence_name=candidate.sequence_name,
                    contig=candidate.contig,
                    class_id=get_class_id('hybrid'),
                    class_name='hybrid',
                    subclass='Overlapping_motifs',
                    motif_id=candidate.motif_id,
                    start=candidate.start,
                    end=candidate.end,
                    length=candidate.length,
                    matched_seq=candidate.matched_seq,
                    pattern_name=f"hybrid_{candidate.class_name}",
                    overlap_classes=overlap_classes
                )
                hybrid_candidates.append(hybrid)
        
        return hybrid_candidates
    
    def _simple_overlap_detection(self, candidates: List[Candidate]) -> List[Candidate]:
        """Simple O(n²) overlap detection when intervaltree not available"""
        hybrid_candidates = []
        
        for i, candidate1 in enumerate(candidates):
            overlaps = []
            for j, candidate2 in enumerate(candidates):
                if i != j and self._overlaps(candidate1, candidate2):
                    overlaps.append(candidate2.class_name)
            
            if overlaps:
                hybrid = Candidate(
                    sequence_name=candidate1.sequence_name,
                    contig=candidate1.contig,
                    class_id=get_class_id('hybrid'),
                    class_name='hybrid',
                    subclass='Overlapping_motifs',
                    motif_id=candidate1.motif_id,
                    start=candidate1.start,
                    end=candidate1.end,
                    length=candidate1.length,
                    matched_seq=candidate1.matched_seq,
                    pattern_name=f"hybrid_{candidate1.class_name}",
                    overlap_classes=overlaps
                )
                hybrid_candidates.append(hybrid)
        
        return hybrid_candidates
    
    def _overlaps(self, c1: Candidate, c2: Candidate) -> bool:
        """Check if two candidates overlap"""
        return not (c1.end < c2.start or c2.end < c1.start)
    
    def score(self, candidates: List[Candidate]) -> List[Candidate]:
        """Score hybrid candidates"""
        for candidate in candidates:
            # Score based on number of overlapping classes
            overlap_count = len(candidate.overlap_classes) if candidate.overlap_classes else 0
            candidate.raw_score = overlap_count / 10.0  # Normalize
            candidate.scoring_method = "Overlap_density"
        return candidates


class ClusterDetector(MotifBase):
    """Cluster detector for motif hotspots"""
    
    def __init__(self):
        super().__init__('cluster')
    
    def detect(self, seq: str, seq_name: str, contig: str, offset: int) -> List[Candidate]:
        """Detect clusters - this is called after primary detection"""
        return []
    
    def detect_from_candidates(self, candidates: List[Candidate], window_size: int = 1000) -> List[Candidate]:
        """Detect clusters from existing candidates and report longest regions"""
        cluster_candidates = []
        
        if not candidates:
            return cluster_candidates
        
        # Sort candidates by position
        sorted_candidates = sorted(candidates, key=lambda c: c.start)
        
        # Find all potential clusters with sliding window
        potential_clusters = []
        for i in range(len(sorted_candidates)):
            window_candidates = []
            for j in range(i, len(sorted_candidates)):
                if sorted_candidates[j].start <= sorted_candidates[i].start + window_size:
                    window_candidates.append(sorted_candidates[j])
                else:
                    break
            
            # If window has multiple motifs, consider it a cluster
            if len(window_candidates) >= 3:  # Minimum cluster size
                start_pos = min(c.start for c in window_candidates)
                end_pos = max(c.end for c in window_candidates)
                
                potential_clusters.append({
                    'start': start_pos,
                    'end': end_pos,
                    'length': end_pos - start_pos + 1,
                    'motifs': window_candidates,
                    'density': len(window_candidates) / (end_pos - start_pos + 1)
                })
        
        # Merge overlapping clusters and keep longest regions
        merged_clusters = self._merge_overlapping_clusters(potential_clusters)
        
        # Create cluster candidates, prioritizing longest regions
        for i, cluster_info in enumerate(merged_clusters):
            cluster = Candidate(
                sequence_name=cluster_info['motifs'][0].sequence_name,
                contig=cluster_info['motifs'][0].contig,
                class_id=get_class_id('cluster'),
                class_name='cluster',
                subclass='Motif_hotspot',
                motif_id=i,
                start=cluster_info['start'],
                end=cluster_info['end'],
                length=cluster_info['length'],
                matched_seq=f"cluster_{len(cluster_info['motifs'])}_motifs_longest_{cluster_info['length']}bp".encode(),
                pattern_name=f"cluster_{window_size}bp"
            )
            cluster_candidates.append(cluster)
        
        return cluster_candidates
    
    def _merge_overlapping_clusters(self, clusters):
        """Merge overlapping clusters and prioritize longest regions"""
        if not clusters:
            return []
        
        # Sort by start position
        sorted_clusters = sorted(clusters, key=lambda x: x['start'])
        merged = []
        
        current_cluster = sorted_clusters[0].copy()
        
        for next_cluster in sorted_clusters[1:]:
            # Check if clusters overlap
            if next_cluster['start'] <= current_cluster['end']:
                # Merge clusters, extending to cover the maximum range
                current_cluster['end'] = max(current_cluster['end'], next_cluster['end'])
                current_cluster['length'] = current_cluster['end'] - current_cluster['start'] + 1
                current_cluster['motifs'].extend(next_cluster['motifs'])
                # Remove duplicates while preserving order
                seen = set()
                unique_motifs = []
                for motif in current_cluster['motifs']:
                    motif_key = (motif.start, motif.end, motif.class_name)
                    if motif_key not in seen:
                        seen.add(motif_key)
                        unique_motifs.append(motif)
                current_cluster['motifs'] = unique_motifs
                current_cluster['density'] = len(current_cluster['motifs']) / current_cluster['length']
            else:
                # No overlap, add current cluster to results
                merged.append(current_cluster)
                current_cluster = next_cluster.copy()
        
        # Add the last cluster
        merged.append(current_cluster)
        
        # Sort by length (longest first) to prioritize longest regions
        merged.sort(key=lambda x: x['length'], reverse=True)
        
        return merged
    
    def score(self, candidates: List[Candidate]) -> List[Candidate]:
        """Score cluster candidates with emphasis on longest regions and density"""
        for candidate in candidates:
            # Extract motif count and length from matched_seq
            seq_str = candidate.matched_seq.decode('utf-8')
            if 'cluster_' in seq_str and '_motifs' in seq_str:
                try:
                    # Parse: cluster_{count}_motifs_longest_{length}bp
                    parts = seq_str.split('_')
                    motif_count = int(parts[1])
                    
                    # Calculate density score (motifs per kb)
                    density_score = (motif_count / candidate.length) * 1000 if candidate.length > 0 else 0
                    
                    # Length bonus for longest regions (normalized by max reasonable cluster size)
                    length_score = min(candidate.length / 5000.0, 1.0)  # Cap at 5kb
                    
                    # Combined score: density + length bonus
                    candidate.raw_score = (density_score * 0.7) + (length_score * 0.3)
                    
                except (ValueError, IndexError):
                    # Fallback to length-based scoring
                    candidate.raw_score = candidate.length / 1000.0
            else:
                candidate.raw_score = candidate.length / 1000.0
            
            candidate.scoring_method = "Cluster_longest_density"
        return candidates


# Detector factory
DETECTOR_CLASSES = {
    'g_quadruplex': G4Detector,
    'curved_dna': CurvedDetector,
    'slipped_dna': SlippedDetector,
    'cruciform': CruciformDetector,
    'r_loop': RLoopDetector,
    'triplex': TriplexDetector,
    'i_motif': IMotifDetector,
    'z_dna': ZDNADetector,
    'a_philic': APhilicDetector,
    'hybrid': HybridDetector,
    'cluster': ClusterDetector,
}


def get_detector(class_name: str) -> MotifBase:
    """Get detector instance for a motif class"""
    detector_class = DETECTOR_CLASSES.get(class_name)
    if detector_class is None:
        raise ValueError(f"Unknown motif class: {class_name}")
    return detector_class()


def get_all_detectors() -> Dict[str, MotifBase]:
    """Get all detector instances"""
    return {name: cls() for name, cls in DETECTOR_CLASSES.items()}