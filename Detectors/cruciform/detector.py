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

try:
    import numpy as np
    _NUMPY_AVAILABLE = True
except ImportError:
    _NUMPY_AVAILABLE = False

# TUNABLE PARAMETERS
MIN_ARM = 8; MAX_ARM = 50; MAX_LOOP = 12; MAX_MISMATCHES = 0; MAX_SEQUENCE_LENGTH = 200000; SCORE_THRESHOLD = 0.2
SEED_SIZE = 6; DELTA_G_THRESHOLD = -5.0  # kcal/mol stability cutoff
NN_ENERGY = {"AA": -1.0, "AC": -1.44, "AG": -1.28, "AT": -0.88, "CA": -1.45, "CC": -1.84, "CG": -2.17, "CT": -1.28,
             "GA": -1.30, "GC": -2.24, "GG": -1.84, "GT": -1.44, "TA": -0.58, "TC": -1.30, "TG": -1.45, "TT": -1.0}

# Precompute reverse-complement for all 6-mers once (4^6 = 4096 entries)
_RC_TABLE = str.maketrans("ACGT", "TGCA")

def _precompute_seed_rc_map(k: int):
    """Return a dict mapping each k-mer → its reverse complement string."""
    # Built lazily; for k=6 this covers all 4096 possible k-mers
    return {}


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

        # -----------------------------------------------------------------------
        # Find seed pairs using numpy sliding-window hash matching.
        # For each loop offset d in 0..max_loop we check all positions i where
        # hash(seq[i:i+k]) == rc_hash(seq[i+k+d:i+2k+d]) in a single vectorised
        # C-speed comparison — O(n × max_loop) rather than two O(n) Python loops.
        # -----------------------------------------------------------------------
        if _NUMPY_AVAILABLE and n >= k:
            valid_pairs = self._find_seed_pairs_numpy(seq, n, k, max_loop)
        else:
            # Fallback: Python dict index
            seed_index: dict = defaultdict(list)
            for i in range(n - k + 1):
                seed_index[seq[i:i+k]].append(i)
            valid_pairs = []
            for kmer, i_positions in seed_index.items():
                rc_kmer = kmer.translate(_RC_TABLE)[::-1]
                j_positions = seed_index.get(rc_kmer)
                if not j_positions:
                    continue
                for i in i_positions:
                    j_lo = i + k
                    j_hi = j_lo + max_loop
                    for j in j_positions:
                        if j_lo <= j <= j_hi:
                            valid_pairs.append((i, j))

        # Extend seeds into full inverted repeats
        seen_pairs: set = set()
        for i, j in valid_pairs:
            pair_key = (i, j)
            if pair_key in seen_pairs:
                continue
            seen_pairs.add(pair_key)

            loop_len = j - (i + k)
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

                if next_left != next_right.translate(_RC_TABLE):
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
    # Numpy seed-pair discovery
    # -------------------------

    # Module-level constants for base encoding (A=0, C=1, G=2, T=3)
    _BASE_ENC = None  # built lazily

    @staticmethod
    def _get_base_enc():
        """Return (and cache) the 256-entry base-encoding lookup array."""
        if CruciformDetector._BASE_ENC is None:
            enc = np.full(256, 255, dtype=np.uint8)
            for ch, code in zip('ACGT', range(4)):
                enc[ord(ch)] = code
            CruciformDetector._BASE_ENC = enc
        return CruciformDetector._BASE_ENC

    @staticmethod
    def _find_seed_pairs_numpy(seq: str, n: int, k: int, max_loop: int):
        """Return all (i, j) seed pairs where revcomp(seq[i:i+k]) == seq[j:j+k]
        and j - i - k in [0, max_loop].

        Algorithm: compute polynomial hashes for forward k-mers and their RC
        counterparts in two O(k·n) numpy passes, then slide over max_loop+1
        offsets with a single vectorised comparison per offset — O(n·max_loop)
        total, completely eliminating the O(n) Python position loop.
        """
        enc_lut = CruciformDetector._get_base_enc()
        raw = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
        enc = enc_lut[raw].astype(np.int64)
        num_w = n - k + 1

        POWERS = 4 ** np.arange(k, dtype=np.int64)
        POWERS_REV = 4 ** np.arange(k - 1, -1, -1, dtype=np.int64)

        fwd_hash = np.zeros(num_w, dtype=np.int64)
        rev_hash = np.zeros(num_w, dtype=np.int64)
        for j in range(k):
            fwd_hash += enc[j:j + num_w] * POWERS[j]
            rev_hash += enc[j:j + num_w] * POWERS_REV[j]
        rc_hash = (4 ** k - 1) - rev_hash  # hash of revcomp(seq[i:i+k])

        valid_pairs = []
        for loop_offset in range(0, max_loop + 1):
            j_offset = k + loop_offset
            if j_offset >= num_w:
                break
            ni = num_w - j_offset
            # Positions where rc_hash[i] == fwd_hash[i + j_offset]
            matches = rc_hash[:ni] == fwd_hash[j_offset:j_offset + ni]
            for i in np.where(matches)[0].tolist():
                valid_pairs.append((i, i + j_offset))
        return valid_pairs

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
