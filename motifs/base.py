"""
Base data structures and utilities for motif detection
"""

from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
import numpy as np


@dataclass
class Candidate:
    """Base data structure for detected motif candidates"""
    sequence_name: str
    contig: str
    class_id: int
    class_name: str
    subclass: Optional[str]
    motif_id: int
    start: int    # 1-based inclusive
    end: int      # 1-based inclusive
    length: int
    matched_seq: bytes
    pattern_name: str
    raw_score: Optional[float] = None
    normalized_score: Optional[float] = None
    scoring_method: Optional[str] = None
    gc_content: Optional[float] = None
    overlap_classes: Optional[List[str]] = field(default_factory=list)

    def __post_init__(self):
        """Calculate derived fields"""
        if self.length is None:
            self.length = self.end - self.start + 1
        if self.gc_content is None and self.matched_seq:
            self.gc_content = calculate_gc_content(self.matched_seq)


def calculate_gc_content(sequence: bytes) -> float:
    """Calculate GC content of a sequence"""
    if isinstance(sequence, str):
        sequence = sequence.encode()
    
    sequence = sequence.upper()
    gc_count = sequence.count(b'G') + sequence.count(b'C')
    total_count = len(sequence)
    
    if total_count == 0:
        return 0.0
    
    return gc_count / total_count


def normalize_scores(candidates: List[Candidate], method: str = "minmax") -> List[Candidate]:
    """
    Normalize scores per motif class
    
    Args:
        candidates: List of candidates with raw scores
        method: "minmax" for 0-1 scaling, "zscore" for z-score normalization
    
    Returns:
        Candidates with normalized scores
    """
    if not candidates:
        return candidates
    
    # Group by class for normalization
    by_class: Dict[str, List[Candidate]] = {}
    for candidate in candidates:
        class_key = candidate.class_name
        if class_key not in by_class:
            by_class[class_key] = []
        by_class[class_key].append(candidate)
    
    # Normalize each class separately
    for class_name, class_candidates in by_class.items():
        scores = [c.raw_score for c in class_candidates if c.raw_score is not None]
        if not scores:
            continue
            
        scores_array = np.array(scores)
        
        if method == "minmax":
            if scores_array.max() == scores_array.min():
                normalized = np.ones_like(scores_array)
            else:
                normalized = (scores_array - scores_array.min()) / (scores_array.max() - scores_array.min())
        elif method == "zscore":
            if scores_array.std() == 0:
                normalized = np.zeros_like(scores_array)
            else:
                normalized = (scores_array - scores_array.mean()) / scores_array.std()
        else:
            raise ValueError(f"Unknown normalization method: {method}")
        
        # Update normalized scores
        score_idx = 0
        for candidate in class_candidates:
            if candidate.raw_score is not None:
                candidate.normalized_score = float(normalized[score_idx])
                score_idx += 1
    
    return candidates


def g4hunter_score(sequence: str) -> float:
    """
    Calculate G4Hunter score for a sequence
    
    Based on Bedrat et al. NAR 44(4):1746-1759 (2016)
    TODO: Implement exact formula from TECHNICAL_SPECIFICATIONS.md
    """
    if not sequence:
        return 0.0
    
    # Simple implementation - replace with exact G4Hunter algorithm
    sequence = sequence.upper()
    g_runs = sequence.count('GGG')
    c_runs = sequence.count('CCC')
    
    # Basic scoring heuristic
    score = (g_runs - c_runs) / len(sequence) if len(sequence) > 0 else 0.0
    return score


def triplex_stability_score(sequence: str) -> float:
    """
    Calculate triplex stability score
    
    TODO: Implement exact formula from TECHNICAL_SPECIFICATIONS.md
    """
    if not sequence:
        return 0.0
    
    sequence = sequence.upper()
    purine_content = (sequence.count('A') + sequence.count('G')) / len(sequence)
    pyrimidine_content = (sequence.count('C') + sequence.count('T')) / len(sequence)
    
    # Favor homopurine/homopyrimidine sequences
    return max(purine_content, pyrimidine_content)


def z_dna_score(sequence: str) -> float:
    """
    Calculate Z-DNA propensity score
    
    Based on Z-DNA seeker algorithm, Ho et al. 1986
    TODO: Implement exact formula from TECHNICAL_SPECIFICATIONS.md
    """
    if not sequence:
        return 0.0
    
    sequence = sequence.upper()
    alternating_content = 0
    
    # Look for alternating purine-pyrimidine patterns
    for i in range(len(sequence) - 1):
        if ((sequence[i] in 'AG' and sequence[i+1] in 'CT') or
            (sequence[i] in 'CT' and sequence[i+1] in 'AG')):
            alternating_content += 1
    
    return alternating_content / (len(sequence) - 1) if len(sequence) > 1 else 0.0


def curvature_score(sequence: str, window_size: int = 10) -> float:
    """
    Calculate DNA curvature propensity score
    
    TODO: Implement exact formula from TECHNICAL_SPECIFICATIONS.md
    """
    if not sequence or len(sequence) < window_size:
        return 0.0
    
    sequence = sequence.upper()
    a_tract_score = 0
    
    # Look for A-tracts
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i+window_size]
        a_content = window.count('A') / window_size
        if a_content > 0.7:  # High A content indicates curvature
            a_tract_score += a_content
    
    return a_tract_score / (len(sequence) - window_size + 1)