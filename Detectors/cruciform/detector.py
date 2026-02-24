"""Cruciform DNA detector using thermodynamic seed-and-extend algorithm."""
# IMPORTS
import math
from typing import List, Dict, Any, Tuple, Optional
from collections import defaultdict

from ..base.base_detector import BaseMotifDetector
from Utilities.detectors_utils import revcomp, calc_gc_content
from Utilities.core.motif_normalizer import normalize_class_subclass

try: from motif_patterns import CRUCIFORM_PATTERNS
except ImportError: CRUCIFORM_PATTERNS = {}

# TUNABLE PARAMETERS
MIN_ARM = 8; MAX_ARM = 50; MAX_LOOP = 12; MAX_MISMATCHES = 0; MAX_SEQUENCE_LENGTH = 200000; SCORE_THRESHOLD = 0.2
SEED_SIZE = 6; DELTA_G_THRESHOLD = -5.0  # kcal/mol stability cutoff
NN_ENERGY = {"AA": -1.0, "AC": -1.44, "AG": -1.28, "AT": -0.88, "CA": -1.45, "CC": -1.84, "CG": -2.17, "CT": -1.28,
             "GA": -1.30, "GC": -2.24, "GG": -1.84, "GT": -1.44, "TA": -0.58, "TC": -1.30, "TG": -1.45, "TT": -1.0}


class CruciformDetector(BaseMotifDetector):
    """Cruciform (thermodynamic inverted repeat) DNA detector using seed-and-extend indexing."""

    SCORE_REFERENCE = 'Lilley et al. 2000, Sinden et al. 1994'

    MIN_ARM = MIN_ARM; MAX_ARM = MAX_ARM; MAX_LOOP = MAX_LOOP; MAX_MISMATCHES = MAX_MISMATCHES
    MAX_SEQUENCE_LENGTH = MAX_SEQUENCE_LENGTH; SCORE_THRESHOLD = SCORE_THRESHOLD
    SEED_SIZE = SEED_SIZE; DELTA_G_THRESHOLD = DELTA_G_THRESHOLD; NN_ENERGY = NN_ENERGY

    def get_motif_class_name(self) -> str:
        return "Cruciform"

    def get_length_cap(self, subclass: str = None) -> int:
        """Cruciform structures stable up to ~200 bp (Lilley 2000)."""
        return 200

    def theoretical_min_score(self) -> float:
        """Minimum biologically valid cruciform raw score (score threshold)."""
        return self.SCORE_THRESHOLD

    def theoretical_max_score(self, sequence_length: int = None) -> float:
        """Highest possible cruciform raw score.

        Score = min(1.0, -ΔG / 20). For a GC-rich stem of MAX_ARM,
        using most negative NN energy (GC = -2.24 kcal/mol) and minimal loop penalty.
        """
        min_nn = min(self.NN_ENERGY.values())  # most negative (e.g., -2.24)
        max_stem_dG = (self.MAX_ARM - 1) * min_nn
        loop_penalty = self._loop_penalty(1)
        max_neg_dG = max_stem_dG + loop_penalty
        return min(1.0, -max_neg_dG / 20.0)

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        return self._load_patterns(CRUCIFORM_PATTERNS, lambda: {
            'inverted_repeats': [
                (r'', 'CRU_3_1', 'Thermodynamic palindrome',
                 'Cruciform forming IRs', 12,
                 'cruciform_stability', 0.95,
                 'DNA secondary structure',
                 'Lilley 2000; SantaLucia 1998')
            ]
        })

    # -------------------------
    # Thermodynamic Functions
    # -------------------------

    def _calculate_stem_deltaG(self, stem_seq: str) -> float:
        deltaG = 0.0
        for i in range(len(stem_seq) - 1):
            dinuc = stem_seq[i:i+2]
            deltaG += self.NN_ENERGY.get(dinuc, 0.0)
        return deltaG

    def _loop_penalty(self, loop_len: int) -> float:
        if loop_len <= 0:
            return 4.0
        return 1.75 + 0.6 * math.log(loop_len)

    def _calculate_cruciform_deltaG(self, left_arm: str, loop_len: int) -> float:
        stem_dG = self._calculate_stem_deltaG(left_arm)
        loop_dG = self._loop_penalty(loop_len)
        return stem_dG + loop_dG

    def _normalize_deltaG(self, deltaG: float) -> float:
        score = min(1.0, max(0.0, -deltaG / 20.0))
        return round(score, 6)

    # -------------------------
    # Seed-and-Extend Detection
    # -------------------------

    def find_inverted_repeats(self, sequence: str,
                              min_arm: int = None,
                              max_arm: int = None,
                              max_loop: int = None,
                              max_mismatches: int = None) -> List[Dict[str, Any]]:

        seq = sequence.upper()
        n = len(seq)

        min_arm = min_arm or self.MIN_ARM
        max_arm = max_arm or self.MAX_ARM
        max_loop = max_loop or self.MAX_LOOP
        max_mismatches = max_mismatches or self.MAX_MISMATCHES

        hits: List[Dict[str, Any]] = []
        k = self.SEED_SIZE

        seed_index = defaultdict(list)
        for i in range(n - k + 1):
            seed_index[seq[i:i+k]].append(i)

        for i in range(n - k + 1):
            seed = seq[i:i+k]
            seed_rc = revcomp(seed)

            if seed_rc not in seed_index:
                continue

            for j in seed_index[seed_rc]:
                if j <= i:
                    continue

                loop_len = j - (i + k)
                if loop_len < 0 or loop_len > max_loop:
                    continue

                left_start = i
                right_start = j
                arm_len = k
                mismatches = 0

                # Extend outward
                while True:
                    if left_start == 0 or right_start + arm_len >= n:
                        break

                    next_left = seq[left_start - 1]
                    next_right = seq[right_start + arm_len]

                    if next_left != revcomp(next_right):
                        mismatches += 1
                        if mismatches > max_mismatches:
                            break

                    left_start -= 1
                    arm_len += 1

                    if arm_len >= max_arm:
                        break

                if arm_len < min_arm:
                    continue

                left_end = left_start + arm_len
                right_end = right_start + arm_len

                left_seq = seq[left_start:left_end]
                right_seq = seq[right_start:right_end]

                deltaG = self._calculate_cruciform_deltaG(left_seq, loop_len)

                if deltaG > self.DELTA_G_THRESHOLD:
                    continue

                score = self._normalize_deltaG(deltaG)

                hits.append({
                    'left_start': left_start,
                    'left_end': left_end,
                    'right_start': right_start,
                    'right_end': right_end,
                    'arm_len': arm_len,
                    'loop_len': loop_len,
                    'left_seq': left_seq,
                    'right_seq': right_seq,
                    'right_seq_rc': revcomp(right_seq),
                    'mismatches': mismatches,
                    'match_fraction': round((arm_len - mismatches) / arm_len, 4),
                    'deltaG': round(deltaG, 3),
                    'score': score
                })

        hits.sort(key=lambda h: (-h['score'], h['left_start'], -h['arm_len']))
        return hits

    # -------------------------
    # Unchanged Downstream Logic
    # -------------------------

    def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:
        hits = self.find_inverted_repeats(sequence)
        return float(sum(h['score'] for h in hits))

    def annotate_sequence(self, sequence: str, max_hits: int = 0) -> List[Dict[str, Any]]:
        hits = self.find_inverted_repeats(sequence)
        return hits[:max_hits] if max_hits and len(hits) > max_hits else hits

    def passes_quality_threshold(self, sequence: str, score: float, pattern_info: Tuple) -> bool:
        hits = self.find_inverted_repeats(sequence)
        if not hits:
            return False
        best_score = hits[0]['score']
        try:
            provided_thresh = float(pattern_info[6]) if (pattern_info and len(pattern_info) > 6) else None
        except Exception:
            provided_thresh = None
        thresh = provided_thresh if provided_thresh is not None else self.SCORE_THRESHOLD
        return best_score >= thresh

    def _remove_overlaps(self, inverted_repeats: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        if not inverted_repeats:
            return []
        sorted_repeats = sorted(inverted_repeats,
                                key=lambda x: (-x['score'], -(x['right_end'] - x['left_start'])))
        non_overlapping = []
        for repeat in sorted_repeats:
            overlaps = any(not (repeat['right_end'] <= sel['left_start']
                                or repeat['left_start'] >= sel['right_end'])
                           for sel in non_overlapping)
            if not overlaps:
                non_overlapping.append(repeat)
        non_overlapping.sort(key=lambda x: x['left_start'])
        return non_overlapping

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Detect cruciform motifs on BOTH strands (strand-agnostic)."""

        self.audit['invoked'] = True
        self.audit['windows_scanned'] = 0
        self.audit['candidates_seen'] = 0
        self.audit['candidates_filtered'] = 0
        self.audit['reported'] = 0
        self.audit['both_strands_scanned'] = True

        sequence = sequence.upper().strip()
        motifs = []

        self.audit['windows_scanned'] += 1
        inverted_repeats_fwd = self.find_inverted_repeats(sequence)
        self.audit['candidates_seen'] += len(inverted_repeats_fwd)

        filtered_fwd = [r for r in inverted_repeats_fwd if r.get('score', 0) > self.SCORE_THRESHOLD]
        self.audit['candidates_filtered'] += len(inverted_repeats_fwd) - len(filtered_fwd)

        non_overlapping_fwd = self._remove_overlaps(filtered_fwd)

        for i, repeat in enumerate(non_overlapping_fwd):
            start_pos = repeat['left_start']
            end_pos = repeat['right_end']
            full_seq = sequence[start_pos:end_pos]

            left_arm = repeat.get('left_seq', '')
            right_arm = repeat.get('right_seq', '')
            loop_seq = sequence[repeat['left_end']:repeat['right_start']]

            gc_total = calc_gc_content(full_seq)
            gc_left_arm = calc_gc_content(left_arm)
            gc_right_arm = calc_gc_content(right_arm)
            gc_loop = calc_gc_content(loop_seq)

            canonical_class, canonical_subclass = normalize_class_subclass(
                self.get_motif_class_name(),
                'Cruciform forming IRs',
                strict=False,
                auto_correct=True
            )

            motifs.append({
                'ID': f"{sequence_name}_CRU_{start_pos+1}",
                'Sequence_Name': sequence_name,
                'Class': canonical_class,
                'Subclass': canonical_subclass,
                'Start': start_pos + 1,
                'End': end_pos,
                'Length': end_pos - start_pos,
                'Sequence': full_seq,
                'Raw_Score': round(repeat['score'], 3),
                'Score': self.normalize_score(repeat['score'], end_pos - start_pos, canonical_subclass),
                'Strand': '+',
                'Method': 'Cruciform_detection',
                'Pattern_ID': f'CRU_{i+1}',
                'Left_Arm': left_arm,
                'Right_Arm': right_arm,
                'Loop_Seq': loop_seq,
                'Arm_Length': repeat.get('arm_len', 0),
                'Loop_Length': repeat.get('loop_len', 0),
                'Stem_Length': repeat.get('arm_len', 0),
                'GC_Content': round(gc_total, 2),
                'GC_Total': round(gc_total, 2),
                'GC_Left_Arm': round(gc_left_arm, 2),
                'GC_Right_Arm': round(gc_right_arm, 2),
                'GC_Loop': round(gc_loop, 2),
                'Mismatches': repeat.get('mismatches', 0),
                'Match_Fraction': repeat.get('match_fraction', 1.0),
                'DeltaG': repeat.get('deltaG', 0.0),
                'Type_Of_Repeat': 'Inverted repeat (palindromic mirror)',
                'Criterion': self._get_cruciform_criterion(repeat),
                'Disease_Relevance': self._get_cruciform_disease_relevance(repeat.get('arm_len', 0), repeat.get('deltaG', 0), gc_total),
                'Regions_Involved': f'Left arm ({repeat.get("arm_len", 0)}bp) - Loop ({repeat.get("loop_len", 0)}bp) - Right arm ({repeat.get("arm_len", 0)}bp mirror)'
            })

            self.audit['reported'] += 1

        return motifs
    
    def _get_cruciform_criterion(self, repeat: Dict[str, Any]) -> str:
        """Explain cruciform classification criterion"""
        criteria = []
        criteria.append(f'Inverted repeat: arm≥{self.MIN_ARM}bp, loop≤{self.MAX_LOOP}bp')
        criteria.append(f'arm_length={repeat.get("arm_len", 0)}bp, loop_length={repeat.get("loop_len", 0)}bp')
        criteria.append(f'ΔG={repeat.get("deltaG", 0):.2f} kcal/mol')
        criteria.append(f'Score={repeat.get("score", 0):.3f} ≥{self.SCORE_THRESHOLD}')
        return '; '.join(criteria)
    
    def _get_cruciform_disease_relevance(self, arm_len: int, deltaG: float, gc_content: float) -> str:
        """Annotate disease relevance for cruciform structures"""
        disease_notes = []
        
        if deltaG < -10:
            disease_notes.append(f'Highly stable cruciform (ΔG={deltaG:.1f}) - DNA breakage, genomic instability')
        
        if arm_len >= 30:
            disease_notes.append('Long palindrome - chromosomal translocations, deletions')
        
        if gc_content < 40:
            disease_notes.append('AT-rich palindrome - replication fork stalling, fragile sites')
        
        disease_notes.append('Cruciform formation - recombination hotspot, transcription regulation, replication origin')
        
        return '; '.join(disease_notes)
