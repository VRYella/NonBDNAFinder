"""
All 10 motif detector classes - modular design with detect() and score() methods.

Each detector uses Hyperscan for fast anchoring, then applies specific validation
and scoring algorithms per the TECHNICAL_SPECIFICATIONS.md.
"""

from typing import List, Tuple, Dict, Optional
import re
import logging
from abc import ABC, abstractmethod
from collections import defaultdict

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
    """Z-DNA detector using Z-DNA seeker algorithm"""
    
    def __init__(self):
        super().__init__('z_dna')
    
    def detect(self, seq: str, seq_name: str, contig: str, offset: int) -> List[Candidate]:
        """Detect Z-DNA forming sequences"""
        candidates = []
        
        # Z-DNA patterns: alternating purines/pyrimidines
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
            seq_str = candidate.matched_seq.decode('utf-8')
            candidate.raw_score = z_dna_score(seq_str)
            candidate.scoring_method = "Z_seeker_adapted"
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
        """Detect clusters from existing candidates"""
        cluster_candidates = []
        
        if not candidates:
            return cluster_candidates
        
        # Sort candidates by position
        sorted_candidates = sorted(candidates, key=lambda c: c.start)
        
        # Sliding window to find dense regions
        for i in range(len(sorted_candidates)):
            window_candidates = []
            for j in range(i, len(sorted_candidates)):
                if sorted_candidates[j].start <= sorted_candidates[i].start + window_size:
                    window_candidates.append(sorted_candidates[j])
                else:
                    break
            
            # If window has multiple motifs, create cluster
            if len(window_candidates) >= 3:  # Minimum cluster size
                start_pos = min(c.start for c in window_candidates)
                end_pos = max(c.end for c in window_candidates)
                
                cluster = Candidate(
                    sequence_name=window_candidates[0].sequence_name,
                    contig=window_candidates[0].contig,
                    class_id=get_class_id('cluster'),
                    class_name='cluster',
                    subclass='Motif_hotspot',
                    motif_id=i,
                    start=start_pos,
                    end=end_pos,
                    length=end_pos - start_pos + 1,
                    matched_seq=f"cluster_{len(window_candidates)}_motifs".encode(),
                    pattern_name=f"cluster_{window_size}bp"
                )
                cluster_candidates.append(cluster)
        
        return cluster_candidates
    
    def score(self, candidates: List[Candidate]) -> List[Candidate]:
        """Score cluster candidates"""
        for candidate in candidates:
            # Extract motif count from matched_seq
            seq_str = candidate.matched_seq.decode('utf-8')
            if 'cluster_' in seq_str and '_motifs' in seq_str:
                try:
                    motif_count = int(seq_str.split('_')[1])
                    candidate.raw_score = motif_count / 10.0  # Normalize
                except:
                    candidate.raw_score = candidate.length / 1000.0
            else:
                candidate.raw_score = candidate.length / 1000.0
            candidate.scoring_method = "Cluster_density"
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