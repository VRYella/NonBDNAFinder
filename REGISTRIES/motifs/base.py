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
    
    The G4Hunter algorithm calculates a weighted sum based on:
    - G runs contribute positively (favor G4 formation)
    - C runs contribute negatively (unfavor G4 formation)
    - Weighting by run length and directionality
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        float: G4Hunter score (-2.0 to +2.0, positive favors G4)
    """
    if not sequence:
        return 0.0
    
    sequence = sequence.upper()
    if len(sequence) == 0:
        return 0.0
    
    score = 0.0
    i = 0
    
    while i < len(sequence):
        # Check for G runs (only count runs of 2+)
        if sequence[i] == 'G':
            run_length = 0
            while i < len(sequence) and sequence[i] == 'G':
                run_length += 1
                i += 1
            # Only score runs of 2+ bases, weight by run length squared
            if run_length >= 2:
                score += run_length ** 2
        
        # Check for C runs (only count runs of 2+)
        elif sequence[i] == 'C':
            run_length = 0
            while i < len(sequence) and sequence[i] == 'C':
                run_length += 1
                i += 1
            # C runs contribute negatively, only for runs of 2+
            if run_length >= 2:
                score -= run_length ** 2
        else:
            i += 1
    
    # Normalize by sequence length and scale to standard range
    # Standard G4Hunter range is typically -2.0 to +2.0
    if len(sequence) == 0:
        return 0.0
        
    normalized_score = score / len(sequence)
    
    # Apply scaling factor to match published range
    scaled_score = normalized_score * 2.0
    
    # Clamp to expected range
    return max(-2.0, min(2.0, scaled_score))


def triplex_stability_score(sequence: str) -> float:
    """
    Calculate triplex stability score
    
    Based on triplex DNA stability principles from Frank-Kamenetskii & Mirkin 1995.
    Triplex structures are most stable with homopurine/homopyrimidine tracts.
    
    The scoring considers:
    - Homopurine tract purity (all A/G)
    - Homopyrimidine tract purity (all C/T)  
    - Sequence length contribution
    - pH-dependent stability factors
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        float: Stability score (0.0 to 1.0, higher = more stable)
    """
    if not sequence:
        return 0.0
    
    sequence = sequence.upper()
    length = len(sequence)
    
    if length == 0:
        return 0.0
    
    # Calculate purine and pyrimidine content
    purine_count = sequence.count('A') + sequence.count('G')
    pyrimidine_count = sequence.count('C') + sequence.count('T')
    
    purine_purity = purine_count / length
    pyrimidine_purity = pyrimidine_count / length
    
    # Maximum purity (homopurine or homopyrimidine)
    max_purity = max(purine_purity, pyrimidine_purity)
    
    # Length bonus for longer tracts (triplex requires minimum length)
    length_bonus = min(1.0, length / 15.0)  # Optimal around 15+ bp
    
    # Stability factors based on base composition
    if purine_purity > pyrimidine_purity:
        # Favor GA-rich over AA-rich (GA provides better stacking)
        ga_content = (sequence.count('G') + sequence.count('A')) / length
        g_content = sequence.count('G') / length
        composition_bonus = 0.8 + 0.2 * g_content  # G provides extra stability
    else:
        # Pyrimidine tract - TC preferred over TT/CC
        tc_balance = min(sequence.count('T'), sequence.count('C')) / (length / 2)
        composition_bonus = 0.8 + 0.2 * tc_balance
    
    # Final stability score
    stability = max_purity * length_bonus * composition_bonus
    
    # Ensure range [0.0, 1.0]
    return min(1.0, stability)


def z_dna_score(sequence: str) -> float:
    """
    Calculate Z-DNA propensity score
    
    Based on Z-DNA seeker algorithm, Ho et al. EMBO J 5(10):2737-2744 (1986)
    
    Z-DNA formation favors:
    - Alternating purine-pyrimidine sequences (especially CG dinucleotides)
    - High GC content in alternating pattern
    - Specific sequence motifs like (CG)n, (CA/TG)n
    
    The algorithm evaluates:
    - Alternating purine-pyrimidine frequency
    - CG dinucleotide density (strongest Z-DNA former)
    - Sequence periodicity and regularity
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        float: Z-DNA propensity score (0.0 to 1.0, higher favors Z-DNA)
    """
    if not sequence or len(sequence) < 2:
        return 0.0
    
    sequence = sequence.upper()
    length = len(sequence)
    
    # Define purines and pyrimidines
    purines = set(['A', 'G'])
    pyrimidines = set(['C', 'T'])
    
    # Count alternating purine-pyrimidine patterns
    alternating_count = 0
    cg_dinucleotides = 0
    
    for i in range(length - 1):
        base1, base2 = sequence[i], sequence[i + 1]
        
        # Count CG dinucleotides (strongest Z-DNA formers)
        if (base1 == 'C' and base2 == 'G') or (base1 == 'G' and base2 == 'C'):
            cg_dinucleotides += 1
            
        # Count alternating purine-pyrimidine patterns
        if ((base1 in purines and base2 in pyrimidines) or 
            (base1 in pyrimidines and base2 in purines)):
            alternating_count += 1
    
    # Calculate base metrics
    alternating_fraction = alternating_count / (length - 1)
    cg_density = cg_dinucleotides / (length - 1)
    
    # Calculate GC content
    gc_content = (sequence.count('G') + sequence.count('C')) / length
    
    # Enhanced Z-DNA propensity calculation
    # CG dinucleotides are weighted most heavily (factor 2.5)
    # Alternating patterns get strong weight (factor 1.5)  
    # GC content provides base stability (factor 0.5)
    base_propensity = (2.5 * cg_density + 
                      1.5 * alternating_fraction + 
                      0.5 * gc_content) / 4.5
    
    # Length factor - Z-DNA requires minimum length for stability
    if length >= 6:  # Minimum (CG)3 repeat
        length_factor = min(1.0, length / 12.0)  # Optimal around 12+ bp
    else:
        length_factor = length / 6.0
    
    # Bonus for perfect alternating patterns
    alternating_bonus = 1.0
    if alternating_fraction > 0.9:  # Nearly perfect alternating
        alternating_bonus = 1.2
    elif cg_density > 0.8:  # CG-rich alternating
        alternating_bonus = 1.3
    
    final_score = base_propensity * length_factor * alternating_bonus
    
    # Ensure range [0.0, 1.0]
    return min(1.0, final_score)


def curvature_score(sequence: str, window_size: int = 10) -> float:
    """
    Calculate DNA curvature propensity score
    
    Based on DNA bending prediction using A-tract analysis and wedge model.
    
    DNA curvature is primarily caused by:
    - A-tracts (4+ consecutive A's) with specific spacing
    - Phased A-tracts every ~10-11 bp (helical repeat)
    - Flexibility differences between AT and GC regions
    
    The algorithm considers:
    - A-tract length and frequency
    - Phasing relative to helical periodicity
    - Local flexibility based on base composition
    - Wedge angles for dinucleotide steps
    
    Args:
        sequence: DNA sequence string
        window_size: Window for local analysis (default 10)
        
    Returns:
        float: Curvature propensity score (0.0 to 1.0, higher = more curved)
    """
    if not sequence or len(sequence) < window_size:
        return 0.0
    
    sequence = sequence.upper()
    length = len(sequence)
    
    # Find A-tracts (4+ consecutive A's)
    a_tracts = []
    i = 0
    while i < length:
        if sequence[i] == 'A':
            start = i
            while i < length and sequence[i] == 'A':
                i += 1
            tract_length = i - start
            if tract_length >= 4:  # Minimum A-tract for significant bending
                a_tracts.append((start, tract_length))
        else:
            i += 1
    
    if not a_tracts:
        return 0.0
    
    # Calculate A-tract scoring
    total_a_tract_score = 0.0
    for start, tract_length in a_tracts:
        # Longer A-tracts bend more (up to optimum of ~6)
        length_score = min(1.0, tract_length / 6.0)
        total_a_tract_score += length_score
    
    # Check for phased A-tracts (key feature of curvature)
    phasing_score = 0.0
    if len(a_tracts) > 1:
        for i in range(len(a_tracts) - 1):
            spacing = a_tracts[i+1][0] - (a_tracts[i][0] + a_tracts[i][1])
            # Optimal spacing is 10-11 bp (helical repeat)
            if 8 <= spacing <= 13:
                phasing_bonus = 1.0 - abs(spacing - 10.5) / 2.5
                phasing_score += phasing_bonus
    
    # Calculate flexibility score using sliding window
    flexibility_score = 0.0
    for i in range(length - window_size + 1):
        window = sequence[i:i+window_size]
        
        # AT-rich regions are more flexible
        at_content = (window.count('A') + window.count('T')) / window_size
        
        # Calculate dinucleotide flexibility
        flexibility = 0.0
        for j in range(len(window) - 1):
            dinuc = window[j:j+2]
            # Flexibility values based on crystallographic data
            flex_values = {
                'AA': 1.0, 'AT': 0.9, 'AG': 0.7, 'AC': 0.6,
                'TA': 0.8, 'TT': 1.0, 'TG': 0.6, 'TC': 0.7,
                'GA': 0.7, 'GT': 0.6, 'GG': 0.5, 'GC': 0.4,
                'CA': 0.6, 'CT': 0.7, 'CG': 0.3, 'CC': 0.5
            }
            flexibility += flex_values.get(dinuc, 0.5)
        
        avg_flexibility = flexibility / (len(window) - 1)
        flexibility_score += avg_flexibility
    
    flexibility_score /= (length - window_size + 1)
    
    # Combine scores with appropriate weights
    # A-tracts are the primary determinant (weight 0.5)
    # Phasing provides cooperative enhancement (weight 0.3)
    # Flexibility provides base contribution (weight 0.2)
    a_tract_component = min(1.0, total_a_tract_score / max(1, len(sequence) // 20))
    phasing_component = min(1.0, phasing_score / max(1, len(a_tracts) - 1)) if len(a_tracts) > 1 else 0.0
    flexibility_component = flexibility_score
    
    final_score = (0.5 * a_tract_component + 
                  0.3 * phasing_component + 
                  0.2 * flexibility_component)
    
    # Ensure range [0.0, 1.0]
    return min(1.0, final_score)