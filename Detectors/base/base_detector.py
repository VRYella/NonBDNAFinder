"""Abstract base class for all Non-B DNA motif detectors."""
# IMPORTS
import re
from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple

from Utilities.detectors_utils import (
    calc_gc_content,
    calc_at_content,
    revcomp,
    remove_overlaps,
    remove_overlaps_by_subclass,
    load_patterns_with_fallback
)

# TUNABLE PARAMETERS
DEFAULT_MIN_SCORE_THRESHOLD = 0.5

# Experimentally supported maximum stable lengths per motif class (bp)
STRUCTURAL_LENGTH_CAPS = {
    "A-philic_DNA": 300,
    "Cruciform": 200,
    "Curved_DNA_local": 50,
    "Curved_DNA_global": 120,
    "G-Quadruplex": 120,
    "Z-DNA": 300,
    "Slipped_DNA_STR": 1000,
    "Slipped_DNA_Direct": 500,
    "Triplex": 150,
    "R-Loop": 2000,
    "i-Motif": 60,
}

# Disease-specific cap overrides (copy-number expansions extend stable length)
DISEASE_CAP_OVERRIDES = {
    "Slipped_DNA": {
        "Huntington": 1000,
        "Myotonic_Dystrophy": 1500,
    }
}


class BaseMotifDetector(ABC):
    """Abstract base class for all Non-B DNA motif detectors."""

    SCORE_REFERENCE = 'Override in subclass'
    
    def __init__(self):
        self.patterns = self.get_patterns()
        self.compiled_patterns = self._compile_patterns()
        self._last_raw_score = None  # Store raw score for motif dict creation
        self.audit = {
            'invoked': False,
            'windows_scanned': 0,
            'candidates_seen': 0,
            'candidates_filtered': 0,
            'reported': 0,
            'seed_hits': 0,
            'both_strands_scanned': False
        }
    
    @abstractmethod
    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return patterns specific to this motif class."""
    
    @abstractmethod  
    def get_motif_class_name(self) -> str:
        """Return the motif class name (e.g., 'Curved_DNA', 'G_Quadruplex')"""
        pass
    
    @abstractmethod
    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        """Calculate motif-specific confidence score."""
        pass
    
    def _compile_patterns(self) -> Dict[str, List[Tuple]]:
        """Compile all regex patterns once for performance."""
        compiled_patterns = {}
        
        for pattern_group, patterns in self.patterns.items():
            compiled_group = []
            for pattern_info in patterns:
                pattern, pattern_id, name, subclass = pattern_info[:4]
                try:
                    compiled_re = re.compile(pattern, re.IGNORECASE | re.ASCII)
                    compiled_group.append((compiled_re, pattern_id, name, subclass, pattern_info))
                except re.error as e:
                    print(f"Warning: Invalid pattern {pattern}: {e}")
                    continue
            compiled_patterns[pattern_group] = compiled_group
        
        return compiled_patterns
    
    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Scan sequence for all compiled patterns and return motif list."""
        self.audit['invoked'] = True
        self.audit['windows_scanned'] = 1
        self.audit['candidates_seen'] = 0
        self.audit['candidates_filtered'] = 0
        self.audit['reported'] = 0
        self.audit['seed_hits'] = 0
        
        sequence = sequence.upper().strip()
        motifs = []
        
        for pattern_group, compiled_patterns in self.compiled_patterns.items():
            for compiled_re, pattern_id, name, subclass, full_info in compiled_patterns:
                for match in compiled_re.finditer(sequence):
                    self.audit['seed_hits'] += 1
                    self.audit['candidates_seen'] += 1
                    start, end = match.span()
                    motif_seq = sequence[start:end]
                    
                    score = self.calculate_score(motif_seq, full_info)
                    
                    if self.passes_quality_threshold(motif_seq, score, full_info):
                        motifs.append({
                            'ID': f"{sequence_name}_{pattern_id}_{start+1}",
                            'Sequence_Name': sequence_name,
                            'Class': self.get_motif_class_name(),
                            'Subclass': subclass,
                            'Start': start + 1,
                            'End': end,
                            'Length': len(motif_seq),
                            'Sequence': motif_seq,
                            'Raw_Score': round(score, 6),
                            'Score': self.normalize_score(score, len(motif_seq), subclass),
                            'Strand': '+',
                            'Method': f'{self.get_motif_class_name()}_detection',
                            'Pattern_ID': pattern_id
                        })
                        self.audit['reported'] += 1
                    else:
                        self.audit['candidates_filtered'] += 1
        
        return motifs
    
    def passes_quality_threshold(self, sequence: str, score: float, pattern_info: Tuple) -> bool:
        """Apply quality thresholds - can be overridden by subclasses"""
        if len(pattern_info) > 6:
            min_threshold = pattern_info[6]
            return score >= min_threshold
        return score >= DEFAULT_MIN_SCORE_THRESHOLD

    def theoretical_min_score(self) -> float:
        """Return minimum biologically valid raw score. Override in subclass."""
        raise NotImplementedError

    def theoretical_max_score(self, sequence_length: int = None) -> float:
        """Return highest possible raw score under model. Override in subclass."""
        raise NotImplementedError

    def get_length_cap(self, subclass: str = None) -> int:
        """
        Return experimentally supported maximum stable length (bp)
        for this motif class or subclass.

        Each detector overrides this to return the structural upper bound
        beyond which the length component is capped at 1.0.
        """
        raise NotImplementedError

    def normalize_score(self, raw_score: float, motif_length: int,
                        subclass: str = None) -> float:
        """
        Convert raw detector-specific score to universal 1–3 scale using
        length-aware, structure-bounded normalization.

        Final Score = 1 + 2 × f_raw × f_length

        Where:
            f_raw    = (RawScore − RawMin) / (RawMax(L) − RawMin)
            f_length = min(1, L / L_cap)

        Score is clamped to [1.0, 3.0].

        Universal Score Interpretation:
        ┌────────────┬──────────────────┬─────────────────────────────────────┐
        │ Score      │ Interpretation   │ Biological Meaning                  │
        ├────────────┼──────────────────┼─────────────────────────────────────┤
        │ 1.0 - 1.7  │ Weak/Conditional │ Low confidence, context-dependent   │
        │ 1.7 - 2.3  │ Moderate         │ Reasonable confidence, likely valid │
        │ 2.3 - 3.0  │ Strong/High      │ High confidence, well-characterized │
        └────────────┴──────────────────┴─────────────────────────────────────┘

        Args:
            raw_score: Raw detector-specific score
            motif_length: Motif length in base pairs
            subclass: Optional motif subclass for subclass-specific length caps

        Returns:
            Normalized score in range [1.0, 3.0]
        """
        raw_min = self.theoretical_min_score()
        raw_max = self.theoretical_max_score(motif_length)

        if raw_max == raw_min:
            raw_component = 0.0
        else:
            raw_component = (raw_score - raw_min) / (raw_max - raw_min)

        raw_component = max(0.0, min(1.0, raw_component))

        length_cap = self.get_length_cap(subclass)
        length_component = min(1.0, motif_length / float(length_cap))

        combined = raw_component * length_component

        final_score = 1.0 + 2.0 * combined

        return round(min(3.0, final_score), 3)
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get detector statistics"""
        total_patterns = sum(len(patterns) for patterns in self.patterns.values())
        return {
            'motif_class': self.get_motif_class_name(),
            'total_patterns': total_patterns,
            'pattern_groups': list(self.patterns.keys()),
            'patterns_by_group': {k: len(v) for k, v in self.patterns.items()}
        }
    
    def get_audit_info(self) -> Dict[str, Any]:
        """Get detector execution audit information"""
        return self.audit.copy()
    
    def _calc_gc(self, seq: str) -> float:
        """Calculate GC content using shared utility"""
        return calc_gc_content(seq)
    
    def _calc_at(self, seq: str) -> float:
        """Calculate AT content using shared utility"""
        return calc_at_content(seq)
    
    def _revcomp(self, seq: str) -> str:
        """Reverse complement using shared utility"""
        return revcomp(seq)
    
    def _remove_overlaps(self, motifs: List[Dict[str, Any]], 
                        by_subclass: bool = False) -> List[Dict[str, Any]]:
        """Remove overlaps using shared utility"""
        if by_subclass:
            return remove_overlaps_by_subclass(motifs)
        return remove_overlaps(motifs)
    
    def _load_patterns(self, patterns: Any, fallback_func) -> Dict:
        """Load patterns with fallback using shared utility"""
        return load_patterns_with_fallback(patterns, fallback_func)
