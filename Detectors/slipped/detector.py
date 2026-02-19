"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Slipped DNA Motif Detector                                                   │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.2            │
│ Unified slippage-prone repeat detection with mechanistic scoring             │
│ References:                                                                  │
│ - Sinden RR (1994) DNA Structure & Function                                  │
│ - Pearson CE et al. (2005) Nat Rev Genet                                     │
│ - Mirkin SM (2007) Nat Rev Genet                                             │
└──────────────────────────────────────────────────────────────────────────────┘
"""
# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
import re
import math
from typing import List, Dict, Any, Tuple
from ..base.base_detector import BaseMotifDetector
from Utilities.core.motif_normalizer import normalize_class_subclass

# Try to import Numba for JIT compilation (2-5x speedup)
try:
    from numba import jit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    # Create a no-op decorator if Numba is not available
    def jit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator

# ═══════════════════════════════════════════════════════════════════════════════
# TUNABLE PARAMETERS (Literature-grounded)
# ═══════════════════════════════════════════════════════════════════════════════
MIN_TRACT_LENGTH = 20; MIN_PURITY = 0.90; MAX_UNIT_SIZE = 100
STR_CORE_UNITS = (1, 4); STR_RELAXED_UNITS = (1, 6)
MIN_COPIES_STR_CORE = 6; MIN_COPIES_STR_RELAXED = 4
MIN_COPIES_DR = 2; DR_UNIT_MIN = 7; DR_UNIT_MAX = 50
ENTROPY_MIN = 0.5
STR_DIRECT_THRESHOLD = 7  # Unit size threshold: <7 = STR, >=7 = Direct Repeat
SCORING_MODE = "CORE"  # "CORE" or "LENIENT"
# ═══════════════════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════════════════
# NORMALIZATION PARAMETERS (Tunable)
# ═══════════════════════════════════════════════════════════════════════════════
# ┌──────────────┬─────────────┬────────────────────────────────────────┐
# │ Parameter    │ Value       │ Scientific Basis                       │
# ├──────────────┼─────────────┼────────────────────────────────────────┤
# │ RAW_MIN      │ 0.3         │ Minimal repeat instability             │
# │ RAW_MAX      │ 0.98        │ High repeat instability (STRs)         │
# │ NORM_MIN     │ 1.0         │ Universal low confidence threshold     │
# │ NORM_MAX     │ 3.0         │ Universal high confidence threshold    │
# │ METHOD       │ 'linear'    │ Linear interpolation                   │
# └──────────────┴─────────────┴────────────────────────────────────────┘
SLIPPED_RAW_SCORE_MIN = 0.3; SLIPPED_RAW_SCORE_MAX = 0.98
SLIPPED_NORMALIZED_MIN = 1.0; SLIPPED_NORMALIZED_MAX = 3.0
SLIPPED_NORMALIZATION_METHOD = 'linear'
SLIPPED_SCORE_REFERENCE = 'Schlötterer et al. 2000, Weber et al. 1989'
# ═══════════════════════════════════════════════════════════════════════════════

# Module-level cache for pre-compiled tandem-repeat regex patterns.
# Keyed by (k, min_copies); shared across all SlippedDNADetector instances.
_TR_PATTERN_CACHE: Dict[tuple, Any] = {}


# ═══════════════════════════════════════════════════════════════════════════════
# JIT-COMPILED HELPER FUNCTIONS FOR PERFORMANCE
# ═══════════════════════════════════════════════════════════════════════════════
@jit(nopython=True, cache=True)
def _calculate_entropy_jit(a_count, c_count, g_count, t_count, seq_len):
    """JIT-compiled entropy calculation for better performance.
    
    This provides 2-3x speedup for entropy calculations when Numba is available.
    """
    if seq_len == 0:
        return 0.0
    
    entropy = 0.0
    for count in [a_count, c_count, g_count, t_count]:
        if count > 0:
            prob = count / seq_len
            entropy -= prob * math.log2(prob)
    
    return entropy


@jit(nopython=True, cache=True)
def _compute_slippage_base_score_jit(tract_length, copy_number, k, purity, gc_fraction):
    """JIT-compiled slippage base score calculation for better performance.
    
    This provides 2-3x speedup for scoring calculations when Numba is available.
    """
    # Factor 1: Total tract length (log-scaled)
    L = min(1.0, math.log(max(1, tract_length)) / math.log(25) / 2.0)
    
    # Factor 2: Copy number (log-scaled)
    C = min(1.0, math.log(max(1, copy_number)) / math.log(4) / 2.0)
    
    # Factor 3: Unit-size weighting
    if 2 <= k <= 4:
        U = 1.0
    elif k == 1 or 5 <= k <= 6:
        U = 0.75
    elif k <= 20:
        U = 0.5
    else:
        U = 0.3
    
    # Factor 4: Repeat purity (quadratic penalty for impurity)
    P = purity ** 2
    
    # Factor 5: GC-based stability
    G = 0.6 + 0.4 * gc_fraction
    
    # Combine factors
    raw_score = 0.30 * L + 0.30 * C + 0.15 * U + 0.15 * P + 0.10 * G
    
    return 1.0 + 2.0 * min(1.0, raw_score)
# ═══════════════════════════════════════════════════════════════════════════════


class SlippedDNADetector(BaseMotifDetector):
    """Unified detector for slippage-prone DNA: STRs (k=1-6) and direct repeats (k≥7). Requires ≥20 bp tracts with ≥90% purity."""
    
    # Override normalization parameters
    RAW_SCORE_MIN = SLIPPED_RAW_SCORE_MIN
    RAW_SCORE_MAX = SLIPPED_RAW_SCORE_MAX
    NORMALIZED_MIN = SLIPPED_NORMALIZED_MIN
    NORMALIZED_MAX = SLIPPED_NORMALIZED_MAX
    NORMALIZATION_METHOD = SLIPPED_NORMALIZATION_METHOD
    SCORE_REFERENCE = SLIPPED_SCORE_REFERENCE
    
    MIN_TRACT_LENGTH = MIN_TRACT_LENGTH; MIN_PURITY = MIN_PURITY
    MIN_COPIES_STR_CORE = MIN_COPIES_STR_CORE; MIN_COPIES_STR_RELAXED = MIN_COPIES_STR_RELAXED
    MIN_COPIES_DR = MIN_COPIES_DR; MAX_UNIT_SIZE = MAX_UNIT_SIZE
    STR_DIRECT_THRESHOLD = STR_DIRECT_THRESHOLD; SCORING_MODE = SCORING_MODE
    
    def get_motif_class_name(self) -> str:
        return "Slipped_DNA"
    
    @staticmethod
    def compute_primitive_motif(sequence: str) -> str:
        """Compute primitive (irreducible) repeat unit. Returns shortest substring generating full sequence. E.g., CAGCAGCAGCAG → CAG, handles partial repeats."""
        n = len(sequence)
        if n == 0:
            return ""
        
        # Try all possible periods from 1 to n/2 (a repeat must occur at least twice)
        for period in range(1, n // 2 + 1):
            unit = sequence[:period]
            
            is_primitive = True
            for i in range(0, n, period):
                # Check each position against the unit
                check_len = min(period, n - i)
                if sequence[i:i+check_len] != unit[:check_len]:
                    is_primitive = False
                    break
            
            if is_primitive:
                return unit
        
        # If no period found, the sequence itself is primitive
        return sequence
    
    @staticmethod
    def compute_repeat_purity(sequence: str, unit: str) -> float:
        """Compute repeat purity: fraction of bases matching perfect repeat when unit aligned repeatedly. Returns 0-1, where 1.0 = perfect repeat."""
        if not unit or not sequence:
            return 0.0
        
        unit_len = len(unit)
        seq_len = len(sequence)
        
        if unit_len > seq_len:
            return 0.0
        
        # Count matching bases when aligning unit repeatedly
        matches = 0
        for i in range(seq_len):
            if sequence[i] == unit[i % unit_len]:
                matches += 1
        
        return matches / seq_len
    
    @staticmethod
    def calculate_entropy(sequence: str) -> float:
        """Calculate Shannon entropy of DNA sequence. Returns 0-2 bits (2 = maximum complexity).
        
        Uses JIT-compiled version for 2-3x speedup when Numba is available.
        """
        if not sequence:
            return 0.0
        
        # Count bases in a single pass
        a_count = sequence.count('A')
        c_count = sequence.count('C')
        g_count = sequence.count('G')
        t_count = sequence.count('T')
        seq_len = a_count + c_count + g_count + t_count
        
        if seq_len == 0:
            return 0.0
        
        # Use JIT-compiled version for better performance
        if NUMBA_AVAILABLE:
            return _calculate_entropy_jit(a_count, c_count, g_count, t_count, seq_len)
        else:
            # Fallback to original implementation
            from math import log2
            freq = {'A': a_count, 'C': c_count, 'G': g_count, 'T': t_count}
            entropy = -sum((count / seq_len) * log2(count / seq_len) 
                           for count in freq.values() if count > 0)
            return entropy

    def compute_slippage_energy_score(self, sequence: str, unit: str, 
                                      copy_number: float, purity: float) -> float:
        """Compute mechanistic slippage energy score (1-3 scale). Literature-aligned: unit-size instability weighting, disease motif bonus, GC-based stability. 1.0-1.5=weak, 1.5-2.5=moderate, 2.5-3.0=strong.
        
        Uses JIT-compiled base calculation for 2-3x speedup when Numba is available.
        """
        tract_length = len(sequence)
        k = len(unit)
        
        # Calculate GC fraction
        gc_count = sequence.count('G') + sequence.count('C')
        gc_fraction = gc_count / len(sequence) if len(sequence) > 0 else 0.0
        
        # Use JIT-compiled base score calculation if available
        if NUMBA_AVAILABLE:
            base_score = _compute_slippage_base_score_jit(
                tract_length, copy_number, k, purity, gc_fraction
            )
        else:
            # Fallback to original implementation
            # Factor 1: Total tract length (log-scaled)
            L = min(1.0, math.log(max(1, tract_length), 25) / 2.0)
            
            # Factor 2: Copy number (log-scaled)
            C = min(1.0, math.log(max(1, copy_number), 4) / 2.0)
            
            # Factor 3: Unit-size weighting
            if 2 <= k <= 4:
                U = 1.0
            elif k == 1 or 5 <= k <= 6:
                U = 0.75
            elif k <= 20:
                U = 0.5
            else:
                U = 0.3
            
            # Factor 4: Repeat purity (quadratic penalty for impurity)
            P = purity ** 2
            
            # Factor 5: GC-based stability
            G = 0.6 + 0.4 * gc_fraction
            
            # Combine factors
            raw_score = 0.30 * L + 0.30 * C + 0.15 * U + 0.15 * P + 0.10 * G
            base_score = 1.0 + 2.0 * min(1.0, raw_score)
        
        # Disease motif bonus (applied after base score)
        motif_bonus = 1.15 if unit in ["CAG", "CTG", "CGG", "CCG", "GAA", "TTC"] else 1.0
        
        final_score = base_score * motif_bonus
        return round(min(3.0, max(1.0, final_score)), 3)

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return empty patterns; uses optimized k-mer scanner for detection."""
        return {"short_tandem_repeats": [], "direct_repeats": []}

    def find_all_tandem_repeats(self, sequence: str) -> List[Dict[str, Any]]:
        """Unified tandem repeat finder for k=1 to MAX_UNIT_SIZE.

        Uses pre-compiled regex patterns for each unit size instead of nested
        Python while-loops, eliminating O(k×n) Python-level iteration overhead
        and avoiding spurious single-copy "candidates" that were previously
        created whenever k >= MIN_TRACT_LENGTH.

        The pattern ``(.{k})\\1{m-1,}`` requires at least *m* consecutive exact
        copies of a k-character unit, where m = max(2, ceil(MIN_TRACT_LENGTH/k)).
        This is semantically identical to the previous double-loop but runs at
        C speed via the ``re`` module.

        Compiled patterns are cached in the module-level ``_TR_PATTERN_CACHE``
        dict and shared across all ``SlippedDNADetector`` instances.
        """
        seq = sequence.upper()
        n = len(seq)
        candidates = []

        for k in range(1, min(self.MAX_UNIT_SIZE + 1, n // 2)):
            # Minimum copy count: at least 2 AND enough copies to meet MIN_TRACT_LENGTH
            min_copies = max(2, math.ceil(self.MIN_TRACT_LENGTH / k))

            key = (k, min_copies)
            if key not in _TR_PATTERN_CACHE:
                _TR_PATTERN_CACHE[key] = re.compile(
                    rf'(.{{{k}}})\1{{{min_copies - 1},}}'
                )
            pattern = _TR_PATTERN_CACHE[key]

            for m in pattern.finditer(seq):
                unit = m.group(1)
                if 'N' in unit:
                    continue
                full_seq = m.group(0)
                copies = len(full_seq) // k
                candidates.append({
                    'start': m.start(),
                    'end': m.end(),
                    'length': len(full_seq),
                    'unit': unit,
                    'unit_size': k,
                    'copies': copies,
                    'sequence': full_seq,
                })

        return candidates
    
    def apply_stringent_criteria(self, candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Apply hard gates with tiered thresholds (CORE vs LENIENT). Entropy computed on full tract, not primitive-only. Returns high-confidence candidates."""
        filtered = []
        
        for cand in candidates:
            sequence = cand['sequence']
            
            # Gate 1: Minimum tract length (already checked in detection, but reconfirm)
            if cand['length'] < self.MIN_TRACT_LENGTH:
                continue
            
            # Gate 2: Compute primitive motif from full sequence
            primitive_unit = self.compute_primitive_motif(sequence)
            purity = self.compute_repeat_purity(sequence, primitive_unit)
            
            if purity < self.MIN_PURITY:
                continue
            
            # Gate 3: Compute entropy on FULL TRACT (not just primitive motif)
            entropy = self.calculate_entropy(sequence)
            if entropy < ENTROPY_MIN:
                continue
            
            # Gate 4: Minimum copy number with tiered thresholds (CORE vs LENIENT)
            k = len(primitive_unit)
            primitive_copies = len(sequence) / k if k > 0 else 0
            
            if self.SCORING_MODE == "CORE":
                # CORE mode: stricter thresholds for disease-like regime
                if k <= 4 and primitive_copies < self.MIN_COPIES_STR_CORE:
                    continue
                if k > 4 and primitive_copies < self.MIN_COPIES_DR:
                    continue
            else:
                # LENIENT mode: relaxed thresholds for discovery
                if k <= 6 and primitive_copies < self.MIN_COPIES_STR_RELAXED:
                    continue
                if k > 6 and primitive_copies < self.MIN_COPIES_DR:
                    continue
            
            # Passed all gates - add to filtered list with enriched data
            cand['primitive_unit'] = primitive_unit
            cand['purity'] = purity
            cand['entropy'] = entropy
            cand['primitive_copies'] = primitive_copies
            filtered.append(cand)
        
        return filtered
    
    def eliminate_redundancy(self, candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Eliminate redundant calls using max-k dominance: retain longest primitive motif for overlapping candidates. One call per genomic locus."""
        if not candidates:
            return []
        
        sorted_cands = sorted(candidates, key=lambda c: (c['start'], -len(c['primitive_unit'])))
        
        non_redundant = []
        used_intervals = []  # Track (start, end) of accepted calls
        
        for cand in sorted_cands:
            start, end = cand['start'], cand['end']
            
            overlaps = False
            for (used_start, used_end) in used_intervals:
                # Two intervals overlap if neither is completely before the other
                if not (end <= used_start or start >= used_end):
                    overlaps = True
                    break
            
            if not overlaps:
                non_redundant.append(cand)
                used_intervals.append((start, end))
        
        return non_redundant
    
    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        """Mechanism-driven slipped DNA detection: 1) unified tandem repeat detection, 2) stringent criteria, 3) redundancy elimination, 4) slippage scoring. Returns non-redundant high-confidence annotations."""
        seq = sequence.upper()
        
        # Stage 1: Find all tandem repeat candidates
        candidates = self.find_all_tandem_repeats(seq)
        
        # Stage 2: Apply stringent criteria
        filtered = self.apply_stringent_criteria(candidates)
        
        # Stage 3: Eliminate redundancy (one call per locus)
        non_redundant = self.eliminate_redundancy(filtered)
        
        # Stage 4: Compute mechanistic slippage energy scores
        for cand in non_redundant:
            score = self.compute_slippage_energy_score(
                sequence=cand['sequence'],
                unit=cand['primitive_unit'],
                copy_number=cand['copies'],
                purity=cand['purity']
            )
            cand['slippage_score'] = score
        
        return non_redundant
    
    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Return publication-ready slipped DNA annotations with literature references. Guarantees: ONE call per locus, biologically defensible via tiered criteria."""
        self.audit['invoked'] = True
        self.audit['windows_scanned'] = 1
        self.audit['candidates_seen'] = 0
        self.audit['candidates_filtered'] = 0
        self.audit['reported'] = 0
        
        annotations = self.annotate_sequence(sequence)
        self.audit['candidates_seen'] = len(annotations)
        
        motifs = []
        
        for i, ann in enumerate(annotations):
            # Use primitive_copies if available, otherwise fall back to original copies
            copy_number = ann.get('primitive_copies', ann.get('copies', 0))
            
            # Determine canonical subclass based on unit size (new threshold: <7 = STR)
            unit_size = len(ann['primitive_unit'])
            canonical_subclass = 'STR' if unit_size < self.STR_DIRECT_THRESHOLD else 'Direct Repeat'
            
            # Normalize class/subclass using canonical taxonomy
            canonical_class, canonical_subclass = normalize_class_subclass(
                self.get_motif_class_name(),
                canonical_subclass,
                strict=False,
                auto_correct=True
            )
            
            # Calculate GC content
            gc_count = ann['sequence'].count('G') + ann['sequence'].count('C')
            gc_content = round(gc_count / len(ann['sequence']) * 100, 2) if len(ann['sequence']) > 0 else 0.0
            
            # Get entropy - should always be present from annotation
            entropy = ann.get('entropy', 0)  # Entropy calculated in apply_stringent_criteria
            
            # Determine disease relevance
            disease_relevance = self._get_disease_relevance(ann['primitive_unit'], copy_number, unit_size)
            
            # Get criterion explanation
            criterion = self._get_slipped_criterion(unit_size, copy_number, ann['purity'], ann.get('entropy', 0))
            
            # Describe regions involved
            regions_involved = self._describe_slipped_regions(ann['primitive_unit'], copy_number, ann['length'])
            
            motif = {
                'ID': f"{sequence_name}_SLIPPED_{ann['start']+1}",
                'Sequence_Name': sequence_name,
                'Class': canonical_class,
                'Subclass': canonical_subclass,
                'Start': ann['start'] + 1,  # 1-based coordinates
                'End': ann['end'],
                'Length': ann['length'],
                'Sequence': ann['sequence'],
                'Repeat_Unit': ann['primitive_unit'],
                'Unit_Size': len(ann['primitive_unit']),
                'Copy_Number': round(copy_number, 2),
                'Purity': round(ann['purity'], 3),
                'Entropy': round(entropy, 3),
                'GC_Content': gc_content,
                'Slippage_Score': round(ann['slippage_score'], 3),
                'Raw_Score': round(ann['slippage_score'], 3),  # Detector-specific scale
                'Score': self._normalize_score(ann['slippage_score']),  # Universal 1-3 scale
                'Strand': '+',
                'Method': 'Slipped_DNA_detection',
                'Pattern_ID': f'SLIPPED_{i+1}',
                'Arm_Length': 'N/A',  # Not applicable for tandem repeats
                'Loop_Length': 'N/A',  # Not applicable for tandem repeats
                'Type_Of_Repeat': self._classify_repeat_type(ann['primitive_unit'], unit_size),
                'Criterion': criterion,
                'Disease_Relevance': disease_relevance,
                'Regions_Involved': regions_involved,
                'References': 'Sinden 1994; Pearson 2005; Mirkin 2007'
            }
            
            motifs.append(motif)
            self.audit['reported'] += 1
        
        return motifs
    
    def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:
        """Calculate score for a sequence (mechanism-driven)."""
        # Find repeats and apply pipeline
        annotations = self.annotate_sequence(sequence)
        if annotations:
            # Return max score among annotations
            return max(ann['slippage_score'] for ann in annotations)
        return 0.0
    
    def _classify_repeat_type(self, unit: str, unit_size: int) -> str:
        """Classify repeat type based on unit composition and size"""
        if unit_size == 1:
            return f'Mononucleotide ({unit})'
        elif unit_size == 2:
            return f'Dinucleotide ({unit})'
        elif unit_size == 3:
            return f'Trinucleotide ({unit})'
        elif unit_size == 4:
            return f'Tetranucleotide ({unit})'
        elif unit_size == 5:
            return f'Pentanucleotide ({unit})'
        elif unit_size == 6:
            return f'Hexanucleotide ({unit})'
        elif unit_size < 20:
            return f'Short tandem repeat ({unit}, {unit_size}bp)'
        else:
            return f'Long direct repeat ({unit_size}bp unit)'
    
    def _get_disease_relevance(self, unit: str, copy_number: float, unit_size: int) -> str:
        """Annotate disease relevance for repeat expansions"""
        disease_notes = []
        
        # Trinucleotide repeat diseases
        if unit == 'CAG' and copy_number > 36:
            disease_notes.append('Huntington disease (n>36), Spinocerebellar ataxias')
        elif unit == 'CTG' and copy_number > 50:
            disease_notes.append('Myotonic dystrophy type 1 (n>50)')
        elif unit == 'CGG' and copy_number > 200:
            disease_notes.append('Fragile X syndrome (n>200)')
        elif unit == 'CCG' and copy_number > 200:
            disease_notes.append('Fragile X syndrome RC strand')
        elif unit == 'GAA' and copy_number > 66:
            disease_notes.append('Friedreich ataxia (n>66)')
        elif unit == 'TTC' and copy_number > 66:
            disease_notes.append('Friedreich ataxia RC strand')
        
        # Hexanucleotide repeat diseases
        elif unit == 'GGGGCC' and copy_number > 30:
            disease_notes.append('C9orf72 ALS/FTD (n>30, most common ALS/FTD mutation)')
        elif unit == 'GGCCCC' and copy_number > 30:
            disease_notes.append('C9orf72 ALS/FTD RC strand')
        
        # Tetranucleotide repeats
        elif unit == 'CCTG' and copy_number > 75:
            disease_notes.append('Myotonic dystrophy type 2 (n>75)')
        elif unit == 'CAGG' and copy_number > 75:
            disease_notes.append('Myotonic dystrophy type 2 RC strand')
        
        # Pentanucleotide repeats
        elif unit == 'ATTCT' and copy_number > 40:
            disease_notes.append('Spinocerebellar ataxia 10 (n>400-4500)')
        
        # General instability thresholds
        elif unit_size <= 6 and copy_number > 20:
            disease_notes.append(f'Expanded repeat (n={copy_number:.1f}, potential instability)')
        elif unit_size <= 6 and copy_number > 10:
            disease_notes.append(f'Intermediate repeat (n={copy_number:.1f}, monitor for expansion)')
        
        # Long direct repeats - genomic instability
        if unit_size >= 20 and copy_number >= 3:
            disease_notes.append('Long direct repeat (genomic instability, deletion/duplication risk)')
        
        return '; '.join(disease_notes) if disease_notes else 'None annotated'
    
    def _get_slipped_criterion(self, unit_size: int, copy_number: float, purity: float, entropy: float) -> str:
        """Explain the classification criterion"""
        criteria = []
        
        # Classification threshold
        if unit_size < self.STR_DIRECT_THRESHOLD:
            criteria.append(f'STR: unit size {unit_size}bp < {self.STR_DIRECT_THRESHOLD}bp threshold')
        else:
            criteria.append(f'Direct Repeat: unit size {unit_size}bp ≥ {self.STR_DIRECT_THRESHOLD}bp threshold')
        
        # Quality metrics
        criteria.append(f'tract ≥{self.MIN_TRACT_LENGTH}bp')
        criteria.append(f'purity {purity:.1%} ≥{self.MIN_PURITY:.0%}')
        criteria.append(f'entropy {entropy:.2f} ≥{ENTROPY_MIN:.1f}')
        
        # Copy number criterion
        if self.SCORING_MODE == "CORE":
            if unit_size <= 4:
                criteria.append(f'copies {copy_number:.1f} ≥{self.MIN_COPIES_STR_CORE} (CORE mode)')
            else:
                criteria.append(f'copies {copy_number:.1f} ≥{self.MIN_COPIES_DR} (CORE mode)')
        else:
            if unit_size <= 6:
                criteria.append(f'copies {copy_number:.1f} ≥{self.MIN_COPIES_STR_RELAXED} (LENIENT mode)')
            else:
                criteria.append(f'copies {copy_number:.1f} ≥{self.MIN_COPIES_DR} (LENIENT mode)')
        
        return '; '.join(criteria)
    
    def _describe_slipped_regions(self, unit: str, copy_number: float, total_length: int) -> str:
        """Describe the regions involved in slipped DNA formation"""
        return f'Tandem repeat of {unit} unit ({len(unit)}bp) × {copy_number:.1f} copies, total {total_length}bp tract'
