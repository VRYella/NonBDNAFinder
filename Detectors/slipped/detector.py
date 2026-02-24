"""Slipped DNA detector: STRs and direct repeats with mechanistic scoring."""
# IMPORTS
import re
import math
from typing import List, Dict, Any, Tuple
from ..base.base_detector import BaseMotifDetector
from Utilities.core.motif_normalizer import normalize_class_subclass
from Utilities.detectors_utils import calc_gc_content

try:
    from numba import jit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    def jit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator

# TUNABLE PARAMETERS (Literature-grounded)
MIN_TRACT_LENGTH = 20; MIN_PURITY = 0.90; MAX_UNIT_SIZE = 100
STR_CORE_UNITS = (1, 4); STR_RELAXED_UNITS = (1, 9)
MIN_COPIES_STR_CORE = 6; MIN_COPIES_STR_RELAXED = 4
MIN_COPIES_DR = 2; DR_UNIT_MIN = 10; DR_UNIT_MAX = 50
STR_DIRECT_THRESHOLD = 10  # Unit size threshold: <=9 = STR, >=10 = Direct Repeat
SCORING_MODE = "CORE"  # "CORE" or "LENIENT"

_TR_PATTERN_CACHE: Dict[tuple, Any] = {}


# JIT-COMPILED HELPER FUNCTIONS FOR PERFORMANCE
@jit(nopython=True, cache=True)
def _compute_slippage_base_score_jit(tract_length, copy_number, k, purity, gc_fraction):
    """JIT-compiled slippage base score calculation."""
    L = min(1.0, math.log(max(1, tract_length)) / math.log(25) / 2.0)
    
    C = min(1.0, math.log(max(1, copy_number)) / math.log(4) / 2.0)
    
    if 2 <= k <= 4:
        U = 1.0
    elif k == 1 or 5 <= k <= 6:
        U = 0.75
    elif k <= 20:
        U = 0.5
    else:
        U = 0.3
    
    P = purity ** 2
    
    G = 0.6 + 0.4 * gc_fraction
    
    raw_score = 0.30 * L + 0.30 * C + 0.15 * U + 0.15 * P + 0.10 * G
    
    return 1.0 + 2.0 * min(1.0, raw_score)


class SlippedDNADetector(BaseMotifDetector):
    """Unified detector for slippage-prone DNA: STRs (k=1-9) and direct repeats (k≥10). Requires ≥20 bp tracts with ≥90% purity."""
    
    SCORE_REFERENCE = 'Schlötterer et al. 2000, Weber et al. 1989'
    
    MIN_TRACT_LENGTH = MIN_TRACT_LENGTH; MIN_PURITY = MIN_PURITY
    MIN_COPIES_STR_CORE = MIN_COPIES_STR_CORE; MIN_COPIES_STR_RELAXED = MIN_COPIES_STR_RELAXED
    MIN_COPIES_DR = MIN_COPIES_DR; MAX_UNIT_SIZE = MAX_UNIT_SIZE
    STR_DIRECT_THRESHOLD = STR_DIRECT_THRESHOLD; SCORING_MODE = SCORING_MODE
    
    def get_motif_class_name(self) -> str:
        return "Slipped_DNA"

    def get_length_cap(self, subclass: str = None) -> int:
        """STR expansions stable up to ~1000 bp; direct repeats up to ~500 bp
        (Schlötterer 2000, Pearson 2005). Disease expansions may reach higher values."""
        if subclass is not None and 'STR' in str(subclass):
            return 1000
        return 500

    def theoretical_min_score(self) -> float:
        """Minimum biologically valid slippage score (already on 1–3 scale)."""
        return 1.0

    def theoretical_max_score(self, sequence_length: int = None) -> float:
        """Highest possible slippage score (already on 1–3 scale, capped at 3.0)."""
        return 3.0
    
    @staticmethod
    def compute_primitive_motif(sequence: str) -> str:
        """Return shortest substring that generates the full sequence (primitive repeat unit)."""
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
        """Return fraction of bases matching perfect repeat alignment (0-1)."""
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
    
    def compute_slippage_energy_score(self, sequence: str, unit: str, 
                                      copy_number: float, purity: float) -> float:
        """Compute mechanistic slippage energy score (1-3 scale)."""
        tract_length = len(sequence)
        k = len(unit)
        
        gc_fraction = calc_gc_content(sequence) / 100.0
        
        if NUMBA_AVAILABLE:
            base_score = _compute_slippage_base_score_jit(
                tract_length, copy_number, k, purity, gc_fraction
            )
        else:
            L = min(1.0, math.log(max(1, tract_length), 25) / 2.0)
            C = min(1.0, math.log(max(1, copy_number), 4) / 2.0)
            if 2 <= k <= 4:
                U = 1.0
            elif k == 1 or 5 <= k <= 6:
                U = 0.75
            elif k <= 20:
                U = 0.5
            else:
                U = 0.3
            P = purity ** 2
            G = 0.6 + 0.4 * gc_fraction
            raw_score = 0.30 * L + 0.30 * C + 0.15 * U + 0.15 * P + 0.10 * G
            base_score = 1.0 + 2.0 * min(1.0, raw_score)
        
        motif_bonus = 1.15 if unit in ["CAG", "CTG", "CGG", "CCG", "GAA", "TTC"] else 1.0
        
        final_score = base_score * motif_bonus
        return round(min(3.0, max(1.0, final_score)), 3)

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return empty patterns; uses optimized k-mer scanner for detection."""
        return {"short_tandem_repeats": [], "direct_repeats": []}

    def find_all_tandem_repeats(self, sequence: str) -> List[Dict[str, Any]]:
        """Find tandem repeats for k=1 to MAX_UNIT_SIZE using compiled regex patterns."""
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
        """Apply hard gates (tract length, purity, copy count)."""
        filtered = []
        
        for cand in candidates:
            sequence = cand['sequence']
            
            if cand['length'] < self.MIN_TRACT_LENGTH:
                continue
            
            primitive_unit = self.compute_primitive_motif(sequence)
            purity = self.compute_repeat_purity(sequence, primitive_unit)
            
            if purity < self.MIN_PURITY:
                continue
            
            k = len(primitive_unit)
            primitive_copies = len(sequence) / k if k > 0 else 0
            
            if self.SCORING_MODE == "CORE":
                if k <= 4 and primitive_copies < self.MIN_COPIES_STR_CORE:
                    continue
                if k > 4 and primitive_copies < self.MIN_COPIES_DR:
                    continue
            else:
                if k <= 9 and primitive_copies < self.MIN_COPIES_STR_RELAXED:
                    continue
                if k > 9 and primitive_copies < self.MIN_COPIES_DR:
                    continue
            
            cand['primitive_unit'] = primitive_unit
            cand['purity'] = purity
            cand['primitive_copies'] = primitive_copies
            filtered.append(cand)
        
        return filtered
    
    def eliminate_redundancy(self, candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Remove overlapping calls; retain longest primitive motif per locus."""
        if not candidates:
            return []
        
        sorted_cands = sorted(candidates, key=lambda c: (c['start'], -len(c['primitive_unit'])))
        
        non_redundant = []
        used_intervals = []  # Track (start, end) of accepted calls
        
        for cand in sorted_cands:
            start, end = cand['start'], cand['end']
            
            overlaps = False
            for (used_start, used_end) in used_intervals:
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
        """Return slipped DNA annotations with biological metadata."""
        self.audit['invoked'] = True
        self.audit['windows_scanned'] = 1
        self.audit['candidates_seen'] = 0
        self.audit['candidates_filtered'] = 0
        self.audit['reported'] = 0
        
        annotations = self.annotate_sequence(sequence)
        self.audit['candidates_seen'] = len(annotations)
        
        motifs = []
        
        for i, ann in enumerate(annotations):
            copy_number = ann.get('primitive_copies', ann.get('copies', 0))
            
            unit_size = len(ann['primitive_unit'])
            canonical_subclass = 'STR' if unit_size <= 9 else 'Direct Repeat'
            
            canonical_class, canonical_subclass = normalize_class_subclass(
                self.get_motif_class_name(),
                canonical_subclass,
                strict=False,
                auto_correct=True
            )
            
            gc_content = round(calc_gc_content(ann['sequence']), 2)
            
            disease_relevance = self._get_disease_relevance(ann['primitive_unit'], copy_number, unit_size)
            
            criterion = self._get_slipped_criterion(unit_size, copy_number, ann['purity'])
            
            regions_involved = self._describe_slipped_regions(ann['primitive_unit'], copy_number, ann['length'])
            
            motif = {
                'ID': f"{sequence_name}_SLIPPED_{ann['start']+1}",
                'Sequence_Name': sequence_name,
                'Class': canonical_class,
                'Subclass': canonical_subclass,
                'Start': ann['start'] + 1,
                'End': ann['end'],
                'Length': ann['length'],
                'Sequence': ann['sequence'],
                'Repeat_Unit': ann['primitive_unit'],
                'Unit_Size': len(ann['primitive_unit']),
                'Copy_Number': round(copy_number, 2),
                'Purity': round(ann['purity'], 3),
                'GC_Content': gc_content,
                'Slippage_Score': round(ann['slippage_score'], 3),
                'Raw_Score': round(ann['slippage_score'], 3),
                'Score': self.normalize_score(ann['slippage_score'], ann['length'], canonical_subclass),
                'Strand': '+',
                'Method': 'Slipped_DNA_detection',
                'Pattern_ID': f'SLIPPED_{i+1}',
                'Arm_Length': 'N/A',
                'Loop_Length': 'N/A',
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
        
        elif unit == 'GGGGCC' and copy_number > 30:
            disease_notes.append('C9orf72 ALS/FTD (n>30, most common ALS/FTD mutation)')
        elif unit == 'GGCCCC' and copy_number > 30:
            disease_notes.append('C9orf72 ALS/FTD RC strand')
        
        elif unit == 'CCTG' and copy_number > 75:
            disease_notes.append('Myotonic dystrophy type 2 (n>75)')
        elif unit == 'CAGG' and copy_number > 75:
            disease_notes.append('Myotonic dystrophy type 2 RC strand')
        
        elif unit == 'ATTCT' and copy_number > 40:
            disease_notes.append('Spinocerebellar ataxia 10 (n>400-4500)')
        
        elif unit_size <= 6 and copy_number > 20:
            disease_notes.append(f'Expanded repeat (n={copy_number:.1f}, potential instability)')
        elif unit_size <= 6 and copy_number > 10:
            disease_notes.append(f'Intermediate repeat (n={copy_number:.1f}, monitor for expansion)')
        
        if unit_size >= 20 and copy_number >= 3:
            disease_notes.append('Long direct repeat (genomic instability, deletion/duplication risk)')
        
        return '; '.join(disease_notes) if disease_notes else 'None annotated'
    
    def _get_slipped_criterion(self, unit_size: int, copy_number: float, purity: float) -> str:
        """Explain the classification criterion"""
        criteria = []
        
        if unit_size <= 9:
            criteria.append(f'STR: unit size {unit_size}bp ≤{self.STR_DIRECT_THRESHOLD - 1}bp threshold')
        else:
            criteria.append(f'Direct Repeat: unit size {unit_size}bp ≥{self.STR_DIRECT_THRESHOLD}bp threshold')
        
        criteria.append(f'tract ≥{self.MIN_TRACT_LENGTH}bp')
        criteria.append(f'purity {purity:.1%} ≥{self.MIN_PURITY:.0%}')
        
        if self.SCORING_MODE == "CORE":
            if unit_size <= 4:
                criteria.append(f'copies {copy_number:.1f} ≥{self.MIN_COPIES_STR_CORE} (CORE mode)')
            else:
                criteria.append(f'copies {copy_number:.1f} ≥{self.MIN_COPIES_DR} (CORE mode)')
        else:
            if unit_size <= 9:
                criteria.append(f'copies {copy_number:.1f} ≥{self.MIN_COPIES_STR_RELAXED} (LENIENT mode)')
            else:
                criteria.append(f'copies {copy_number:.1f} ≥{self.MIN_COPIES_DR} (LENIENT mode)')
        
        return '; '.join(criteria)
    
    def _describe_slipped_regions(self, unit: str, copy_number: float, total_length: int) -> str:
        """Describe the regions involved in slipped DNA formation"""
        return f'Tandem repeat of {unit} unit ({len(unit)}bp) × {copy_number:.1f} copies, total {total_length}bp tract'
