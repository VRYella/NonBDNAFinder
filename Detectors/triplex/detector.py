"""Triplex DNA detector: H-DNA mirror repeats and Sticky DNA (GAA/TTC)."""

import re
import math
from typing import List, Dict, Any, Tuple
from collections import defaultdict
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

# Literature constraints
MIN_ARM = 10
MAX_ARM = 100
MAX_LOOP = 8
PURITY_THRESHOLD = 0.90
SEED_SIZE = 6
SCORE_THRESHOLD = 0.25

# Triplex scoring parameters (literature aligned - Frank-Kamenetskii 1995)
H_REF_ARM = 35          # Reference arm length for full stability
H_LOOP_ALPHA = 0.4      # Loop penalty exponent
H_WEIGHT_L = 0.35       # Arm weight
H_WEIGHT_H = 0.20       # Loop weight
H_WEIGHT_P = 0.30       # Purity weight
H_WEIGHT_I = 0.15       # Interruption weight
H_PURITY_MIN = 0.8      # Minimum purity for scoring
H_PURITY_RANGE = 0.2    # Purity scaling range (1.0 - H_PURITY_MIN)

# Sticky thresholds (Sakamoto 1999, FRDA literature)
STICKY_REPLICATION_MIN = 20
STICKY_STABLE_MIN = 40
STICKY_PATHOGENIC_MIN = 60

# Sticky DNA scoring coefficients (piecewise linear scoring)
STICKY_WEAK_SCALE = 0.015           # Score increment per copy in weak range
STICKY_REPLICATION_BASE = 1.3       # Base score at replication blockage threshold
STICKY_REPLICATION_SCALE = 0.03     # Score increment per copy in replication range
STICKY_STABLE_BASE = 2.0            # Base score at sticky threshold
STICKY_STABLE_SCALE = 0.02          # Score increment per copy in sticky range
STICKY_PATHOGENIC_BASE = 2.6        # Base score at pathogenic threshold
STICKY_PATHOGENIC_SCALE = 0.01      # Score increment per copy in pathogenic range

# Sticky DNA specific constants (GAA/TTC repeats - disease-associated)
MIN_STICKY_COPIES = 4  # Minimum 4 copies for Sticky DNA
MIN_STICKY_LENGTH = 12  # Minimum tract length (4 copies * 3 bp)

# NORMALIZATION PARAMETERS - Triplex scores are already on 1–3 scale
TRIPLEX_SCORE_REFERENCE = 'Frank-Kamenetskii et al. 1995'


# JIT-COMPILED HELPER FUNCTIONS FOR PERFORMANCE
@jit(nopython=True, cache=True)
def _compute_purity_jit(purine_count, pyr_count, length):
    """JIT-compiled purity calculation."""
    if length == 0:
        return 0.0
    
    purine_fraction = purine_count / float(length)
    pyr_fraction = pyr_count / float(length)
    
    if purine_fraction > pyr_fraction:
        return purine_fraction
    else:
        return pyr_fraction


@jit(nopython=True, cache=True)
def _score_mirror_triplex_jit(arm_len, loop_len, purity, interruptions,
                               H_REF_ARM, H_LOOP_ALPHA, H_WEIGHT_L, 
                               H_WEIGHT_H, H_WEIGHT_P, H_WEIGHT_I,
                               H_PURITY_MIN, H_PURITY_RANGE):
    """JIT-compiled mirror triplex scoring (returns 1–3 scale)."""
    if H_REF_ARM > 1:
        L = min(1.0, math.log(max(1, arm_len)) / math.log(H_REF_ARM))
    else:
        L = 1.0 if arm_len >= 1 else 0.0

    H = math.exp(-H_LOOP_ALPHA * loop_len)

    if purity >= H_PURITY_MIN and H_PURITY_RANGE > 0:
        P = max(0.0, (purity - H_PURITY_MIN) / H_PURITY_RANGE)
    else:
        P = 0.0

    I = 1.0 / (1.0 + interruptions)

    raw = (
        H_WEIGHT_L * L +
        H_WEIGHT_H * H +
        H_WEIGHT_P * P +
        H_WEIGHT_I * I
    )

    return 1.0 + 2.0 * min(1.0, raw)


class TriplexDetector(BaseMotifDetector):

    # Override normalization parameters
    SCORE_REFERENCE = TRIPLEX_SCORE_REFERENCE

    MIN_ARM = MIN_ARM
    MAX_ARM = MAX_ARM
    MAX_LOOP = MAX_LOOP
    PURITY_THRESHOLD = PURITY_THRESHOLD
    SEED_SIZE = SEED_SIZE
    SCORE_THRESHOLD = SCORE_THRESHOLD
    MIN_STICKY_COPIES = MIN_STICKY_COPIES
    MIN_STICKY_LENGTH = MIN_STICKY_LENGTH

    # Sticky DNA patterns (GAA/TTC repeats associated with Friedreich ataxia)
    STICKY_PATTERNS = [
        (re.compile(r'(?:GAA){4,}', re.IGNORECASE), 'TRX_STICKY_GAA', 'GAA repeat'),
        (re.compile(r'(?:TTC){4,}', re.IGNORECASE), 'TRX_STICKY_TTC', 'TTC repeat'),
    ]

    # ------------------------------------------------------------------ #
    # Required Interface
    # ------------------------------------------------------------------ #

    def get_motif_class_name(self) -> str:
        return "Triplex"

    def get_length_cap(self, subclass: str = None) -> int:
        """Triplex/H-DNA stable up to ~150 bp (Frank-Kamenetskii 1995)."""
        return 150

    def theoretical_min_score(self) -> float:
        """Minimum biologically valid triplex score (already on 1–3 scale)."""
        return 1.0

    def theoretical_max_score(self, sequence_length: int = None) -> float:
        """Highest possible triplex score (already on 1–3 scale, capped at 3.0)."""
        return 3.0

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        return {
            "mirror_triplex": [
                (r"", "TRX_MIRROR", "Mirror repeat triplex",
                 "Triplex", self.MIN_ARM,
                 "structural_triplex_score", self.SCORE_THRESHOLD,
                 "Seed-based mirror repeat detection",
                 "Frank-Kamenetskii 1995")
            ],
            "sticky_dna": [
                (r'(?:GAA){4,}', 'TRX_STICKY_GAA', 'GAA repeat',
                 'Sticky DNA', 12, 'sticky_dna_score', 0.95,
                 'Disease-associated repeats', 'Sakamoto 1999'),
                (r'(?:TTC){4,}', 'TRX_STICKY_TTC', 'TTC repeat',
                 'Sticky DNA', 12, 'sticky_dna_score', 0.95,
                 'Disease-associated repeats', 'Sakamoto 1999'),
            ]
        }

    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        return 1.0 if sequence else 0.0

    # ------------------------------------------------------------------ #
    # Mirror Repeat Scoring (H-DNA mechanistic scoring)
    # ------------------------------------------------------------------ #

    def _score_mirror_triplex(self, arm_len: int, loop_len: int,
                              purity: float, interruptions: int) -> float:
        """Mechanistic H-DNA scoring (Frank-Kamenetskii 1995), returns 1–3 scale."""
        if NUMBA_AVAILABLE:
            score = _score_mirror_triplex_jit(
                arm_len, loop_len, purity, interruptions,
                H_REF_ARM, H_LOOP_ALPHA, H_WEIGHT_L,
                H_WEIGHT_H, H_WEIGHT_P, H_WEIGHT_I,
                H_PURITY_MIN, H_PURITY_RANGE
            )
            return round(score, 3)
        else:
            if H_REF_ARM > 1:
                L = min(1.0, math.log(max(1, arm_len)) / math.log(H_REF_ARM))
            else:
                L = 1.0 if arm_len >= 1 else 0.0

            H = math.exp(-H_LOOP_ALPHA * loop_len)

            if purity >= H_PURITY_MIN and H_PURITY_RANGE > 0:
                P = max(0.0, (purity - H_PURITY_MIN) / H_PURITY_RANGE)
            else:
                P = 0.0

            I = 1.0 / (1.0 + interruptions)

            raw = (
                H_WEIGHT_L * L +
                H_WEIGHT_H * H +
                H_WEIGHT_P * P +
                H_WEIGHT_I * I
            )

            return round(1.0 + 2.0 * min(1.0, raw), 3)

    # ------------------------------------------------------------------ #
    # Mirror Repeat Core (C-style logic reproduction)
    # ------------------------------------------------------------------ #

    def _find_mirror_repeats(self, sequence: str) -> List[Dict[str, Any]]:

        seq = sequence.upper()
        n = len(seq)
        k = self.SEED_SIZE
        hits = []

        seed_index = defaultdict(list)

        for i in range(n - k + 1):
            seed_index[seq[i:i+k]].append(i)

        for i in range(n - k + 1):

            seed = seq[i:i+k]
            mirror = seed[::-1]

            if mirror not in seed_index:
                continue

            for j in seed_index[mirror]:

                if j <= i:
                    continue

                loop_len = j - (i + k)

                # Strict spacer constraint
                if loop_len < 0 or loop_len > self.MAX_LOOP:
                    continue

                left_start = i
                right_start = j
                arm_len = k

                while (
                    left_start > 0 and
                    right_start + arm_len < n and
                    seq[left_start - 1] == seq[right_start + arm_len] and
                    arm_len < self.MAX_ARM
                ):
                    left_start -= 1
                    arm_len += 1

                if arm_len < self.MIN_ARM:
                    continue

                left_seq = seq[left_start:left_start + arm_len]

                purine_fraction = sum(b in "AG" for b in left_seq) / arm_len
                pyr_fraction = sum(b in "CT" for b in left_seq) / arm_len
                purity = max(purine_fraction, pyr_fraction)

                if purity < self.PURITY_THRESHOLD:
                    continue

                interruptions = sum(
                    1 for b in left_seq
                    if (purine_fraction > pyr_fraction and b not in "AG")
                    or (pyr_fraction > purine_fraction and b not in "CT")
                )

                score = self._score_mirror_triplex(
                    arm_len, loop_len, purity, interruptions
                )

                hits.append({
                    "start": left_start,
                    "end": right_start + arm_len,
                    "arm_length": arm_len,
                    "loop_length": loop_len,
                    "purity": round(purity, 3),
                    "interruptions": interruptions,
                    "score": score,
                    "sequence": seq[left_start:right_start + arm_len]
                })

        hits.sort(key=lambda x: (-x["arm_length"], x["loop_length"], x["start"]))
        return hits

    # ------------------------------------------------------------------ #
    # Sticky DNA Detection (GAA/TTC repeats)
    # ------------------------------------------------------------------ #

    def _score_sticky_dna(self, copy_number: int) -> Tuple[float, Dict[str, bool]]:
        """
        Piecewise Sticky DNA scoring (Sakamoto 1999, FRDA thresholds).
        Returns score (1–3) and biological flags.
        """
        flags = {
            "Replication_Blockage_Range": False,
            "Sticky_Threshold_Range": False,
            "Pathogenic_Range": False
        }

        if copy_number < STICKY_REPLICATION_MIN:
            score = 1.0 + STICKY_WEAK_SCALE * copy_number

        elif copy_number < STICKY_STABLE_MIN:
            flags["Replication_Blockage_Range"] = True
            score = STICKY_REPLICATION_BASE + STICKY_REPLICATION_SCALE * (copy_number - STICKY_REPLICATION_MIN)

        elif copy_number < STICKY_PATHOGENIC_MIN:
            flags["Sticky_Threshold_Range"] = True
            score = STICKY_STABLE_BASE + STICKY_STABLE_SCALE * (copy_number - STICKY_STABLE_MIN)

        else:
            flags["Pathogenic_Range"] = True
            score = STICKY_PATHOGENIC_BASE + STICKY_PATHOGENIC_SCALE * (copy_number - STICKY_PATHOGENIC_MIN)

        return min(3.0, round(score, 3)), flags

    def _find_sticky_dna(self, sequence: str) -> List[Dict[str, Any]]:
        """Find Sticky DNA (GAA/TTC trinucleotide repeats)."""
        seq = sequence.upper()
        hits = []

        for pattern, pattern_id, description in self.STICKY_PATTERNS:
            for match in pattern.finditer(seq):
                start, end = match.start(), match.end()
                tract_length = end - start
                
                if 'GAA' in pattern_id:
                    unit = 'GAA'
                else:
                    unit = 'TTC'
                copy_number = tract_length // 3
                
                score, flags = self._score_sticky_dna(copy_number)
                
                hits.append({
                    "start": start,
                    "end": end,
                    "length": tract_length,
                    "score": score,
                    "sequence": seq[start:end],
                    "pattern_id": pattern_id,
                    "subclass": "Sticky DNA",
                    "repeat_unit": unit,
                    "copy_number": copy_number,
                    "description": description,
                    "flags": flags
                })

        hits.sort(key=lambda x: (-x["length"], x["start"]))
        return hits

    # ------------------------------------------------------------------ #
    # Annotation with overlap control
    # ------------------------------------------------------------------ #

    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        """
        Annotate sequence with both mirror triplex and sticky DNA motifs.
        
        Note: Sticky DNA and mirror triplex motifs are allowed to overlap since
        they represent different structural features. However, overlaps within
        the same subclass are removed.
        """
        seq = sequence.upper()
        results = []
        
        # Detect mirror repeats (Triplex subclass)
        mirror_used = [False] * len(seq)
        mirrors = self._find_mirror_repeats(seq)

        for m in mirrors:
            s = m["start"]
            e = m["end"]

            if any(mirror_used[s:e]):
                continue

            for i in range(s, e):
                mirror_used[i] = True

            results.append({
                "class_name": "Triplex",
                "subclass": "Triplex",
                "pattern_id": "TRX_MIRROR",
                "start": s,
                "end": e,
                "length": e - s,
                "score": m["score"],
                "matched_seq": m["sequence"],
                "arm_length": m.get("arm_length", 0),
                "loop_length": m.get("loop_length", 0),
                "purity": m.get("purity", 0)
            })

        # Detect sticky DNA (Sticky DNA subclass) - separate overlap tracking
        sticky_used = [False] * len(seq)
        sticky_hits = self._find_sticky_dna(seq)
        
        for hit in sticky_hits:
            s = hit["start"]
            e = hit["end"]
            
            if any(sticky_used[s:e]):
                continue
            
            for i in range(s, e):
                sticky_used[i] = True
            
            results.append({
                "class_name": "Triplex",
                "subclass": "Sticky DNA",
                "pattern_id": hit["pattern_id"],
                "start": s,
                "end": e,
                "length": hit["length"],
                "score": hit["score"],
                "matched_seq": hit["sequence"],
                "repeat_unit": hit.get("repeat_unit"),
                "copy_number": hit.get("copy_number"),
                "flags": hit.get("flags", {})
            })

        results.sort(key=lambda r: r["start"])
        return results

    # ------------------------------------------------------------------ #
    # Final Detection Interface
    # ------------------------------------------------------------------ #

    def detect_motifs(self,
                      sequence: str,
                      sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Detect triplex DNA motifs (mirror repeats and sticky DNA)."""
        sequence = sequence.upper().strip()
        motifs = []
        annotations = self.annotate_sequence(sequence)

        for annotation in annotations:
            subclass = annotation.get("subclass", "Triplex")
            
            canonical_class, canonical_subclass = normalize_class_subclass(
                self.get_motif_class_name(),
                subclass,
                strict=False,
                auto_correct=True
            )
            
            if canonical_subclass == "Sticky DNA":
                method = "Sticky_DNA_detection"
                id_prefix = "STICKY"
            else:
                method = "Triplex_seed_mirror_detection"
                id_prefix = "TRX"
            
            motif_seq = annotation["matched_seq"]
            
            gc_content = round(calc_gc_content(motif_seq), 2)

            motif = {
                "ID": f"{sequence_name}_{id_prefix}_{annotation['start']+1}",
                "Sequence_Name": sequence_name,
                "Class": canonical_class,
                "Subclass": canonical_subclass,
                "Start": annotation["start"] + 1,
                "End": annotation["end"],
                "Length": annotation["length"],
                "Sequence": motif_seq,
                "Raw_Score": annotation.get("score", 1.0),
                "Score": self.normalize_score(annotation.get("score", 1.0), annotation["length"], canonical_subclass),
                "Strand": "+",
                "Method": method,
                "Pattern_ID": annotation["pattern_id"],
                "GC_Content": gc_content
            }
            
            # Add extra fields for Sticky DNA
            if canonical_subclass == "Sticky DNA":
                if annotation.get("repeat_unit"):
                    motif["Repeat_Unit"] = annotation["repeat_unit"]
                    motif["Type_Of_Repeat"] = f"Trinucleotide ({annotation['repeat_unit']})"
                if annotation.get("copy_number"):
                    motif["Copy_Number"] = annotation["copy_number"]
                
                flags = annotation.get("flags", {})
                motif.update(flags)
                
                copy_num = annotation.get("copy_number", 0)
                motif["Arm_Length"] = "N/A"
                motif["Loop_Length"] = "N/A"
                motif["Criterion"] = self._get_sticky_criterion(copy_num, flags)
                motif["Disease_Relevance"] = self._get_sticky_disease_relevance(
                    annotation.get("repeat_unit", ""),
                    copy_num,
                    flags
                )
                motif["Regions_Involved"] = f"{annotation.get('repeat_unit', 'GAA/TTC')} trinucleotide repeat × {copy_num} copies"
                
            else:
                arm_len = annotation.get("arm_length", 0)
                loop_len = annotation.get("loop_length", 0)
                purity = annotation.get("purity", 0)
                
                if arm_len > 0:
                    motif["Arm_Length"] = arm_len
                if loop_len >= 0:
                    motif["Loop_Length"] = loop_len
                if purity > 0:
                    motif["Purity"] = round(purity, 3)
                
                motif["Type_Of_Repeat"] = "Mirror repeat (inverted)"
                motif["Criterion"] = self._get_mirror_criterion(arm_len, loop_len, purity)
                motif["Disease_Relevance"] = self._get_mirror_disease_relevance(motif_seq, arm_len)
                motif["Regions_Involved"] = self._describe_mirror_regions(arm_len, loop_len)

            motifs.append(motif)

        return motifs
    
    def _get_sticky_criterion(self, copy_number: int, flags: Dict[str, bool]) -> str:
        """Explain Sticky DNA classification criterion"""
        criteria = [f"GAA/TTC trinucleotide repeat ≥{self.MIN_STICKY_COPIES} copies"]
        criteria.append(f"n={copy_number} copies detected")
        
        if flags.get("Pathogenic_Range"):
            criteria.append("Pathogenic range (n≥60, FRDA)")
        elif flags.get("Sticky_Threshold_Range"):
            criteria.append("Sticky threshold range (40≤n<60)")
        elif flags.get("Replication_Blockage_Range"):
            criteria.append("Replication blockage range (20≤n<40)")
        else:
            criteria.append("Weak range (4≤n<20)")
        
        return "; ".join(criteria)
    
    def _get_sticky_disease_relevance(self, unit: str, copy_number: int, flags: Dict[str, bool]) -> str:
        """Annotate disease relevance for Sticky DNA"""
        disease_notes = []
        
        if unit in ['GAA', 'TTC']:
            if copy_number >= 66:
                disease_notes.append(f"Friedreich ataxia (FRDA): pathogenic n≥66, detected n={copy_number}")
            elif copy_number >= 40:
                disease_notes.append(f"Friedreich ataxia risk: intermediate n={copy_number} (normal <33)")
            elif copy_number >= 20:
                disease_notes.append(f"Replication stress: n={copy_number} may cause polymerase stalling")
            else:
                disease_notes.append(f"Weak triplex potential: n={copy_number}")
        
        if flags.get("Pathogenic_Range"):
            disease_notes.append("Forms stable H-DNA/triplex, disease-associated")
        elif flags.get("Sticky_Threshold_Range"):
            disease_notes.append("Sticky DNA threshold - stable triplex formation")
        elif flags.get("Replication_Blockage_Range"):
            disease_notes.append("May block replication fork progression")
        
        return "; ".join(disease_notes) if disease_notes else "None annotated"
    
    def _get_mirror_criterion(self, arm_len: int, loop_len: int, purity: float) -> str:
        """Explain mirror triplex classification criterion"""
        criteria = []
        criteria.append(f"Mirror repeat: arm≥{self.MIN_ARM}bp, loop≤{self.MAX_LOOP}bp")
        if arm_len > 0:
            criteria.append(f"arm_length={arm_len}bp")
        if loop_len >= 0:
            criteria.append(f"loop_length={loop_len}bp")
        if purity > 0:
            criteria.append(f"Pu/Py purity={purity:.1%} ≥{self.PURITY_THRESHOLD:.0%}")
        return "; ".join(criteria)
    
    def _get_mirror_disease_relevance(self, sequence: str, arm_len: int) -> str:
        """Annotate disease relevance for mirror triplex"""
        disease_notes = []
        
        # Check for purine-rich or pyrimidine-rich sequences
        purine_frac = sum(b in "AG" for b in sequence) / len(sequence) if len(sequence) > 0 else 0
        
        if purine_frac > 0.9 or purine_frac < 0.1:
            disease_notes.append("Strong Pu/Py bias - H-DNA formation potential")
        
        if arm_len >= 30:
            disease_notes.append("Long mirror repeat - genomic instability risk, chromosomal rearrangement")
        
        if arm_len >= 50:
            disease_notes.append("PKD1-like structure (polycystic kidney disease associated)")
        
        if not disease_notes:
            disease_notes.append("H-DNA formation potential (transcription, replication regulation)")
        
        return "; ".join(disease_notes)
    
    def _describe_mirror_regions(self, arm_len: int, loop_len: int) -> str:
        """Describe regions involved in mirror triplex formation"""
        if arm_len > 0:
            return f"Left arm ({arm_len}bp) - Loop ({loop_len}bp) - Right arm ({arm_len}bp mirror)"
        else:
            return "Mirror repeat structure with symmetrical arms"
