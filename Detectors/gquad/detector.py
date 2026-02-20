"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ G-Quadruplex Detector - Ultra-fast seeded G4 detection                       │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
│ References: Huppert 2005, Bedrat 2016                                        │
└──────────────────────────────────────────────────────────────────────────────┘
"""
# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
import re
from typing import Dict, List, Tuple, Any
from ..base.base_detector import BaseMotifDetector
from Utilities.core.motif_normalizer import normalize_class_subclass
from Utilities.detectors_utils import calc_gc_content

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
# TUNABLE PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════════
WINDOW_SIZE_DEFAULT = 25
MIN_REGION_LEN = 8
CLASS_PRIORITY = ["telomeric_g4", "stacked_canonical_g4s", "stacked_g4s_linker", "canonical_g4", "extended_loop_g4", "higher_order_g4", "g_triplex", "weak_pqs"]

# ═══════════════════════════════════════════════════════════════════════════════
# NORMALIZATION PARAMETERS (Tunable)
# ═══════════════════════════════════════════════════════════════════════════════
# ┌──────────────┬─────────────┬────────────────────────────────────────┐
# │ Parameter    │ Value       │ Scientific Basis                       │
# ├──────────────┼─────────────┼────────────────────────────────────────┤
# │ RAW_MIN      │ 0.5         │ Bedrat 2016 - minimal G4 stability     │
# │ RAW_MAX      │ 1.0         │ Bedrat 2016 - maximal G4Hunter score   │
# │ NORM_MIN     │ 1.0         │ Universal low confidence threshold     │
# │ NORM_MAX     │ 3.0         │ Universal high confidence threshold    │
# │ METHOD       │ 'g4hunter'  │ G4Hunter-specific normalization        │
# └──────────────┴─────────────┴────────────────────────────────────────┘
G4_RAW_SCORE_MIN = 0.5; G4_RAW_SCORE_MAX = 1.0
G4_NORMALIZED_MIN = 1.0; G4_NORMALIZED_MAX = 3.0
G4_NORMALIZATION_METHOD = 'g4hunter'  # G4Hunter abs-value based scoring
G4_SCORE_REFERENCE = 'Bedrat et al. 2016 (Bioinformatics)'
# ═══════════════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════════════
# JIT-COMPILED SCORING FUNCTIONS FOR PERFORMANCE
# ═══════════════════════════════════════════════════════════════════════════════
@jit(nopython=True, cache=True)
def _compute_max_window_sum_jit(vals, ws, L):
    """JIT-compiled sliding window maximum sum computation.
    
    This provides 2-5x speedup for G4 scoring when Numba is available.
    """
    if ws <= 0 or L < ws:
        return 0
    
    cur = 0
    for i in range(ws):
        cur += vals[i]
    
    max_sum = cur
    for i in range(1, L - ws + 1):
        cur += vals[i + ws - 1] - vals[i - 1]
        if cur > max_sum:
            max_sum = cur
    
    return max_sum
# ═══════════════════════════════════════════════════════════════════════════════

class GQuadruplexDetector(BaseMotifDetector):
    """Ultra-fast seeded G4 detector with priority logic retained."""

    # Override normalization parameters
    RAW_SCORE_MIN = G4_RAW_SCORE_MIN
    RAW_SCORE_MAX = G4_RAW_SCORE_MAX
    NORMALIZED_MIN = G4_NORMALIZED_MIN
    NORMALIZED_MAX = G4_NORMALIZED_MAX
    NORMALIZATION_METHOD = G4_NORMALIZATION_METHOD
    SCORE_REFERENCE = G4_SCORE_REFERENCE

    # -------------------------
    # Core Interface
    # -------------------------

    def get_motif_class_name(self) -> str:
        return "G-Quadruplex"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        return {
            'telomeric_g4': [(r'(?:TTAGGG){4,}', 'G4_TEL', 'Telomeric G4', 'Telomeric G4')],
            'stacked_canonical_g4s': [(r'(?:(?:G{3,}[ACGT]{1,7}){3}G{3,}){2,}', 'G4_STK_CAN', 'Stacked canonical G4s', 'Stacked canonical G4s')],
            'stacked_g4s_linker': [(r'(?:(?:G{3,}[ACGT]{1,7}){3}G{3,})(?:[ACGT]{0,20}(?:(?:G{3,}[ACGT]{1,7}){3}G{3,})){1,}', 'G4_STK_LNK', 'Stacked G4s with linker', 'Stacked G4s with linker')],
            'canonical_g4': [(r'G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}', 'G4_CAN', 'Canonical intramolecular G4', 'Canonical intramolecular G4')],
            'extended_loop_g4': [(r'G{3,}[ACGT]{1,12}G{3,}[ACGT]{1,12}G{3,}[ACGT]{1,12}G{3,}', 'G4_EXT', 'Extended-loop canonical', 'Extended-loop canonical')],
            'higher_order_g4': [(r'(?:G{3,}[ACGT]{1,7}){7,}', 'G4_HIGH', 'Higher-order G4 array/G4-wire', 'Higher-order G4 array/G4-wire')],
            'g_triplex': [(r'G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}', 'G4_TRX', 'Intramolecular G-triplex', 'Intramolecular G-triplex')],
            'weak_pqs': [(r'G{2,}[ACGT]{1,7}G{2,}[ACGT]{1,7}G{2,}[ACGT]{1,7}G{2,}', 'G4_WEAK', 'Two-tetrad weak PQS', 'Two-tetrad weak PQS')],
        }

    # -------------------------
    # Public API
    # -------------------------

    def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:
        annotations = self.annotate_sequence(sequence)
        return float(sum(a['score'] for a in annotations))

    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        seq = sequence.upper()
        candidates = self._seed_and_scan(seq)
        scored = [self._score_candidate(c, seq) for c in candidates]
        accepted = self._resolve_overlaps(scored)
        return accepted

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        sequence = sequence.upper().strip()
        motifs = []
        annotations = self.annotate_sequence(sequence)

        subclass_map = {
            'telomeric_g4': 'Telomeric G4',
            'stacked_canonical_g4s': 'Stacked canonical G4s',
            'stacked_g4s_linker': 'Stacked G4s with linker',
            'canonical_g4': 'Canonical intramolecular G4',
            'extended_loop_g4': 'Extended-loop canonical',
            'higher_order_g4': 'Higher-order G4 array/G4-wire',
            'g_triplex': 'Intramolecular G-triplex',
            'weak_pqs': 'Two-tetrad weak PQS'
        }

        for ann in annotations:
            subclass = subclass_map.get(ann['class_name'], ann['class_name'])
            canonical_class, canonical_subclass = normalize_class_subclass(
                self.get_motif_class_name(),
                subclass,
                strict=False,
                auto_correct=True
            )

            # Extract motif sequence and analyze structural features
            motif_seq = sequence[ann['start']:ann['end']]
            structural_features = self._extract_g4_features(motif_seq, ann['class_name'])
            
            # Store raw score and calculate normalized score
            raw_score = ann['score']
            normalized_score = self._normalize_score(raw_score)
            
            motif = {
                'ID': f"{sequence_name}_{ann['pattern_id']}_{ann['start']+1}",
                'Sequence_Name': sequence_name,
                'Class': canonical_class,
                'Subclass': canonical_subclass,
                'Start': ann['start'] + 1,
                'End': ann['end'],
                'Length': ann['end'] - ann['start'],
                'Sequence': motif_seq,
                'Raw_Score': round(raw_score, 4),      # Detector-specific scale
                'Score': normalized_score,              # Universal 1-3 scale
                'Strand': '+',
                'Method': 'Seeded_G4Hunter',
                'Pattern_ID': ann['pattern_id']
            }
            
            # Add comprehensive structural features
            motif.update(structural_features)

            motifs.append(motif)

        return motifs
    
    def _extract_g4_features(self, sequence: str, class_name: str) -> Dict[str, Any]:
        """Extract comprehensive G4 structural features including tracts, loops, GC content, etc."""
        features = {}
        
        # Calculate GC content (exclude N/ambiguous bases from denominator)
        features['GC_Content'] = round(calc_gc_content(sequence), 2)
        
        # Find G-tracts (runs of 2+ consecutive Gs)
        g_tracts = [m.group() for m in re.finditer(r'G{2,}', sequence)]
        features['Num_Tracts'] = len(g_tracts)
        
        if g_tracts:
            features['G_Tract_Lengths'] = ','.join(str(len(t)) for t in g_tracts)
            features['Min_Tract_Length'] = min(len(t) for t in g_tracts)
            features['Max_Tract_Length'] = max(len(t) for t in g_tracts)
            features['Avg_Tract_Length'] = round(sum(len(t) for t in g_tracts) / len(g_tracts), 2)
        else:
            features['G_Tract_Lengths'] = ''
            features['Min_Tract_Length'] = 0
            features['Max_Tract_Length'] = 0
            features['Avg_Tract_Length'] = 0.0
        
        # Extract loop regions (sequences between G-tracts)
        loops = []
        if len(g_tracts) > 1:
            parts = re.split(r'G{2,}', sequence)
            loops = [p for p in parts if p]  # Non-empty parts are loops
        
        if loops:
            features['Loop_Lengths'] = ','.join(str(len(l)) for l in loops)
            features['Min_Loop_Length'] = min(len(l) for l in loops)
            features['Max_Loop_Length'] = max(len(l) for l in loops)
            features['Avg_Loop_Length'] = round(sum(len(l) for l in loops) / len(loops), 2)
            features['Loop_Length'] = features['Avg_Loop_Length']  # Unified field for consistency
            features['Num_Loops'] = len(loops)
        else:
            features['Loop_Lengths'] = ''
            features['Min_Loop_Length'] = 0
            features['Max_Loop_Length'] = 0
            features['Avg_Loop_Length'] = 0.0
            features['Loop_Length'] = 0.0  # Unified field for consistency
            features['Num_Loops'] = 0
        
        # Add Arm_Length (maps to average tract length for G4s)
        features['Arm_Length'] = features['Avg_Tract_Length'] if g_tracts else 'N/A'
        
        # Determine motif type and add criterion
        features['Type_Of_Repeat'] = self._classify_g4_type(class_name, len(g_tracts), features)
        features['Criterion'] = self._get_g4_criterion(class_name, features)
        
        # Add disease relevance for telomeric G4s and specific patterns
        features['Disease_Relevance'] = self._get_g4_disease_relevance(class_name, sequence, features)
        
        # Regions involved in motif formation
        features['Regions_Involved'] = self._describe_g4_regions(g_tracts, loops)
        
        return features
    
    def _classify_g4_type(self, class_name: str, num_tracts: int, features: Dict) -> str:
        """Classify G4 structural type"""
        if class_name == 'telomeric_g4':
            return 'Telomeric tandem repeat'
        elif class_name == 'g_triplex':
            return 'Three-tetrad G-triplex'
        elif class_name in ['stacked_canonical_g4s', 'stacked_g4s_linker']:
            return 'Stacked/higher-order G4'
        elif class_name == 'higher_order_g4':
            return 'G4 array (G-wire)'
        elif num_tracts >= 4:
            return 'Four-tetrad intramolecular G4'
        elif num_tracts == 3:
            return 'Three-tetrad G-triplex'
        else:
            return 'G-rich potential quadruplex sequence (PQS)'
    
    def _get_g4_criterion(self, class_name: str, features: Dict) -> str:
        """Explain classification criterion"""
        if class_name == 'telomeric_g4':
            return 'Matches human telomeric repeat (TTAGGG)n with n>=4'
        elif class_name == 'canonical_g4':
            return f'Canonical G4: 4+ G-tracts (≥3G each), loops 1-7bp'
        elif class_name == 'extended_loop_g4':
            return f'Extended-loop G4: 4+ G-tracts (≥3G each), loops 1-12bp'
        elif class_name == 'stacked_canonical_g4s':
            return 'Multiple canonical G4 units stacked without linker'
        elif class_name == 'stacked_g4s_linker':
            return 'Multiple canonical G4 units with 0-20bp linker'
        elif class_name == 'higher_order_g4':
            return 'G4 array: 7+ G-tract repeats forming G-wire structure'
        elif class_name == 'g_triplex':
            return 'G-triplex: 3 G-tracts (≥3G each), loops 1-7bp'
        elif class_name == 'weak_pqs':
            return 'Weak PQS: 4 G-tracts (≥2G each), may form two-tetrad structure'
        else:
            return f'G4Hunter score-based: {features.get("Num_Tracts", 0)} G-tracts detected'
    
    def _get_g4_disease_relevance(self, class_name: str, sequence: str, features: Dict) -> str:
        """Annotate disease relevance for G4 motifs"""
        disease_notes = []
        
        # Telomeric G4s - aging and cancer
        if class_name == 'telomeric_g4':
            disease_notes.append('Telomeric instability (aging, cancer, ALT mechanism)')
        
        # C9orf72 ALS/FTD repeat expansion (GGGGCC)n
        if 'GGGGCC' in sequence or 'GGCCCC' in sequence:
            disease_notes.append('C9orf72 expansion (ALS/FTD, n>30 pathogenic)')
        
        # Fragile X-associated CGG repeats
        if re.search(r'(CGG){5,}', sequence):
            disease_notes.append('CGG repeat expansion (Fragile X, n>200 pathogenic)')
        
        # Myotonic dystrophy type 1 CTG expansion (RC: CAG -> forms G4 on C-rich strand)
        if re.search(r'(CAG){5,}', sequence):
            disease_notes.append('CAG/CTG expansion RC (Myotonic dystrophy, Huntington disease)')
        
        # Promoter G4s - oncogene regulation
        if features.get('Num_Tracts', 0) >= 4 and features.get('GC_Content', 0) > 70:
            disease_notes.append('Promoter-like G4 (potential oncogene regulation: MYC, BCL2, KRAS, VEGF)')
        
        # Long G4 arrays - genomic instability
        if class_name == 'higher_order_g4' or features.get('Length', 0) > 100:
            disease_notes.append('Genomic instability hotspot (DNA breakage, replication stress)')
        
        return '; '.join(disease_notes) if disease_notes else 'None annotated'
    
    def _describe_g4_regions(self, g_tracts: List[str], loops: List[str]) -> str:
        """Describe regions involved in G4 formation"""
        if not g_tracts:
            return 'No G-tracts detected'
        
        parts = []
        parts.append(f'{len(g_tracts)} G-tracts: {", ".join([f"G{len(t)}" for t in g_tracts])}')
        
        if loops:
            parts.append(f'{len(loops)} loops: {", ".join([f"{len(l)}bp" for l in loops])}')
        
        return '; '.join(parts)

    # -------------------------
    # Ultra-Fast Seeding
    # -------------------------

    def _seed_and_scan(self, seq: str) -> List[Dict[str, Any]]:
        """Seed on G3+ tracts, then local regex refinement."""
        candidates = []
        seed_positions = [m.start() for m in re.finditer(r'G{3,}', seq)]

        if not seed_positions:
            return []

        patterns = self.get_patterns()
        scanned_windows = set()

        for seed in seed_positions:
            window_start = max(0, seed - 50)
            window_end = min(len(seq), seed + 200)
            window_key = (window_start, window_end)

            if window_key in scanned_windows:
                continue
            scanned_windows.add(window_key)

            window_seq = seq[window_start:window_end]

            for class_name, pattern_list in patterns.items():
                for pat in pattern_list:
                    regex = pat[0]
                    pattern_id = pat[1] if len(pat) > 1 else f"{class_name}_pat"

                    for m in re.finditer(regex, window_seq):
                        s = window_start + m.start()
                        e = window_start + m.end()
                        if (e - s) >= MIN_REGION_LEN:
                            candidates.append({
                                'class_name': class_name,
                                'pattern_id': pattern_id,
                                'start': s,
                                'end': e
                            })

        return candidates

    # -------------------------
    # G-only G4Hunter Scoring
    # -------------------------

    def _score_candidate(self, candidate: Dict[str, Any], seq: str,
                         window_size: int = WINDOW_SIZE_DEFAULT) -> Dict[str, Any]:

        s, e = candidate['start'], candidate['end']
        region = seq[s:e]
        L = len(region)

        # G-only scoring (C ignored completely)
        if NUMBA_AVAILABLE and L > 50:
            # Use JIT-compiled version for better performance (2-5x speedup)
            import numpy as np
            vals = np.array([1 if ch == 'G' else 0 for ch in region], dtype=np.int32)
            ws = min(window_size, L)
            max_sum = _compute_max_window_sum_jit(vals, ws, L)
        else:
            # Fallback to original implementation for small sequences
            vals = [1 if ch == 'G' else 0 for ch in region]
            ws = min(window_size, L)
            max_sum = 0

            if ws > 0:
                cur = sum(vals[:ws])
                max_sum = cur
                for i in range(1, L - ws + 1):
                    cur += vals[i + ws - 1] - vals[i - 1]
                    if cur > max_sum:
                        max_sum = cur

        normalized_score = max_sum / float(ws) if ws > 0 else 0.0
        region_score = normalized_score * (L / float(ws)) if ws > 0 else 0.0

        out = candidate.copy()
        out['score'] = float(region_score)
        return out

    # -------------------------
    # Priority-based Overlap Resolution
    # -------------------------

    def _resolve_overlaps(self, scored_candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        if not scored_candidates:
            return []

        def class_prio_idx(class_name):
            try:
                return CLASS_PRIORITY.index(class_name)
            except ValueError:
                return len(CLASS_PRIORITY)

        scored_sorted = sorted(
            scored_candidates,
            key=lambda x: (
                -x['score'],
                class_prio_idx(x['class_name']),
                -(x['end'] - x['start'])
            )
        )

        accepted = []
        occupied = []

        for cand in scored_sorted:
            s, e = cand['start'], cand['end']
            conflict = any(not (e <= os or s >= oe) for (os, oe) in occupied)
            if not conflict:
                accepted.append(cand)
                occupied.append((s, e))

        accepted.sort(key=lambda x: x['start'])
        return accepted
