"""
Slipped DNA Motif Detector
==========================

Unified detector for slippage-prone DNA structures including:
- Short Tandem Repeats (STRs): k=1-9 bp unit sizes
- Direct Repeats: k≥10 bp unit sizes

Requires ≥20 bp tracts with ≥90% purity for high-confidence detection.

Scientific References:
---------------------
- Sinden, R.R. (1994). DNA Structure and Function. Academic Press.
  Describes the molecular basis of slipped-strand DNA formation and
  its role in genetic instability.

- Pearson, C.E., Nichol Edamura, K., and Cleary, J.D. (2005).
  Repeat instability: mechanisms of dynamic mutations.
  Nature Reviews Genetics 6(10):729-742.
  Comprehensive review of repeat expansion mechanisms and slippage.

- Mirkin, S.M. (2007). Expandable DNA repeats and human disease.
  Nature 447:932-940.
  Links slipped-strand structures to repeat expansion diseases.

Algorithm:
----------
1. Unified tandem repeat detection for k=1 to MAX_UNIT_SIZE
2. Stringent quality filters (length, purity, copy number, entropy)
3. Redundancy elimination (one call per genomic locus)
4. Mechanistic slippage energy scoring (1-3 scale)

Performance:
-----------
- O(n) complexity with efficient algorithmic detection
- No catastrophic regex backtracking
- Handles sequences up to 100 Kbp efficiently
"""

import re
from typing import List, Dict, Any, Tuple
from ..base.base_detector import BaseMotifDetector
from core.motif_normalizer import normalize_class_subclass


class SlippedDNADetector(BaseMotifDetector):
    """Unified detector for slippage-prone DNA: STRs (k=1-9) and direct repeats (k≥10). Requires ≥20 bp tracts with ≥90% purity. References: Sinden (1994), Pearson (2005), Mirkin (2007)."""
    
    MIN_TRACT_LENGTH = 20; MIN_PURITY = 0.90; MIN_COPIES_STR = 3; MIN_COPIES_DIRECT = 2; MAX_UNIT_SIZE = 100
    STR_DIRECT_REPEAT_THRESHOLD = 10  # Unit size threshold: <10 = STR, >=10 = Direct Repeat
    
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
        """Calculate Shannon entropy of DNA sequence. Returns 0-2 bits (2 = maximum complexity)."""
        from math import log2
        if not sequence:
            return 0.0
        
        # Calculate base frequencies in a single pass (O(n) complexity)
        freq = {}
        for base in sequence:
            if base in "ACGT":
                freq[base] = freq.get(base, 0) + 1
        
        # Normalize to probabilities
        seq_len = sum(freq.values())
        if seq_len == 0:
            return 0.0
        
        # Calculate Shannon entropy
        entropy = -sum((count / seq_len) * log2(count / seq_len) 
                       for count in freq.values() if count > 0)
        return entropy

    def compute_slippage_energy_score(self, sequence: str, unit: str, 
                                      copy_number: float, purity: float) -> float:
        """Compute mechanistic slippage energy score (1-3 scale). Integrates length, purity, unit size, copy number, GC-based stability. 1.0-1.5=weak, 1.5-2.5=moderate, 2.5-3.0=strong."""
        import math
        
        tract_length = len(sequence)
        unit_size = len(unit)
        
        # Factor 1: Total tract length (dominant, scales logarithmically)
        # Longer tracts → more stable slip-outs
        length_factor = min(1.0, math.log(max(1, tract_length), 20) / 2.0)
        
        # Factor 2: Repeat purity (interruptions penalize stability)
        # High purity required for stable slip-out
        purity_factor = purity ** 2  # Quadratic penalty for impurity
        
        # Factor 3: Repeat unit size (k)
        # Longer units → more register ambiguity
        unit_factor = min(1.0, math.log(max(1, unit_size), 2) / 3.5)
        
        # Factor 4: Copy number
        # More copies → more alignment possibilities
        copy_factor = min(1.0, math.log(max(1, copy_number), 3) / 2.0)
        
        # Factor 5: NN-based ΔG proxy (simple GC content heuristic)
        # Higher GC → more stable base stacking
        gc_count = sequence.count('G') + sequence.count('C')
        gc_fraction = gc_count / len(sequence) if len(sequence) > 0 else 0.0
        stability_factor = 0.5 + 0.5 * gc_fraction  # Range [0.5, 1.0]
        
        # Composite score (weighted sum)
        raw_score = (
            0.35 * length_factor +
            0.30 * purity_factor +
            0.15 * unit_factor +
            0.10 * copy_factor +
            0.10 * stability_factor
        )
        
        # Normalize to [1-3] scale (ΔG-inspired: 1=weak, 2=moderate, 3=strong)
        # Map [0-1] raw_score to [1-3]
        normalized_score = 1.0 + (2.0 * raw_score)
        
        return min(3.0, max(1.0, normalized_score))

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return empty patterns; uses optimized k-mer scanner for detection."""
        return {"short_tandem_repeats": [], "direct_repeats": []}

    def find_all_tandem_repeats(self, sequence: str) -> List[Dict[str, Any]]:
        """Unified tandem repeat finder for k=1 to MAX_UNIT_SIZE. Efficient algorithmic detection, no catastrophic backtracking. Returns candidate repeat dicts with raw data."""
        seq = sequence.upper()
        n = len(seq)
        candidates = []
        
        for k in range(1, min(self.MAX_UNIT_SIZE + 1, n // 2)):
            # Slide window across sequence looking for repeats
            current_pos = 0
            while current_pos < n - k:
                unit = seq[current_pos:current_pos+k]
                
                # Skip units with ambiguous bases
                if 'N' in unit:
                    current_pos += 1
                    continue
                
                # Count consecutive copies of this unit
                copies = 1
                repeat_end_pos = current_pos + k
                while repeat_end_pos + k <= n and seq[repeat_end_pos:repeat_end_pos+k] == unit:
                    copies += 1
                    repeat_end_pos += k
                
                tract_length = copies * k
                if tract_length >= self.MIN_TRACT_LENGTH:
                    # Record this candidate
                    candidates.append({
                        'start': current_pos,
                        'end': repeat_end_pos,
                        'length': tract_length,
                        'unit': unit,
                        'unit_size': k,
                        'copies': copies,
                        'sequence': seq[current_pos:repeat_end_pos]
                    })
                    
                    # Skip past this repeat to avoid overlapping detections
                    current_pos = repeat_end_pos
                else:
                    current_pos += 1
        
        return candidates
    
    def apply_stringent_criteria(self, candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Apply hard gates: tract length ≥ MIN_TRACT_LENGTH, purity ≥ MIN_PURITY, copies ≥ MIN_COPIES, entropy check. Returns high-confidence candidates."""
        filtered = []
        
        for cand in candidates:
            sequence = cand['sequence']
            copies = cand['copies']
            unit_size = cand['unit_size']
            
            # Gate 1: Minimum tract length (already checked in detection, but reconfirm)
            if cand['length'] < self.MIN_TRACT_LENGTH:
                continue
            
            # Gate 2: Compute primitive motif from full sequence (not just unit)
            # This ensures we get the true primitive even if detected unit is composite
            primitive_unit = self.compute_primitive_motif(sequence)
            purity = self.compute_repeat_purity(sequence, primitive_unit)
            
            if purity < self.MIN_PURITY:
                continue
            
            # Gate 3: Minimum copy number (recompute based on primitive unit)
            primitive_copies = len(sequence) / len(primitive_unit) if len(primitive_unit) > 0 else 0
            min_copies = self.MIN_COPIES_STR if len(primitive_unit) < self.STR_DIRECT_REPEAT_THRESHOLD else self.MIN_COPIES_DIRECT
            if primitive_copies < min_copies:
                continue
            
            # Gate 4: Entropy check (exclude low-complexity)
            entropy = self.calculate_entropy(primitive_unit)
            if entropy < 0.5:  # Very low entropy threshold for homopolymers
                continue
            
            # Passed all gates - add to filtered list with enriched data
            cand['primitive_unit'] = primitive_unit
            cand['purity'] = purity
            cand['entropy'] = entropy
            cand['primitive_copies'] = primitive_copies  # Updated copy count based on primitive
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
        """Return publication-ready slipped DNA annotations. Guarantees: ONE call per locus (no STR_1...STR_9 spam), biologically defensible via stringent criteria."""
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
            
            # Determine canonical subclass based on unit size
            unit_size = len(ann['primitive_unit'])
            canonical_subclass = 'STR' if unit_size < self.STR_DIRECT_REPEAT_THRESHOLD else 'Direct Repeat'
            
            # Normalize class/subclass using canonical taxonomy
            canonical_class, canonical_subclass = normalize_class_subclass(
                self.get_motif_class_name(),
                canonical_subclass,
                strict=False,
                auto_correct=True
            )
            
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
                'Copy_Number': copy_number,
                'Purity': round(ann['purity'], 3),
                'Slippage_Energy_Score': round(ann['slippage_score'], 3),
                'Score': round(ann['slippage_score'], 3),  # For compatibility
                'Strand': '+',
                'Method': 'Slipped_DNA_detection',
                'Pattern_ID': f'SLIPPED_{i+1}'
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
