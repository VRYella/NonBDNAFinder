"""A-philic DNA detector using 10-mer propensity scoring."""
# IMPORTS
import logging
from typing import Any, Dict, List, Tuple

from Detectors.base.base_detector import BaseMotifDetector

from Utilities.core.motif_normalizer import normalize_class_subclass
from Utilities.detectors_utils import calc_gc_content, calc_at_content
from .tenmer_table import TENMER_LOG2

try:
    from Detectors.zdna import hyperscan_backend
    _HYPERSCAN_AVAILABLE = hyperscan_backend.is_hyperscan_available()
except ImportError:
    _HYPERSCAN_AVAILABLE = False
    hyperscan_backend = None

try:
    import numpy as np
    _NUMPY_AVAILABLE = True
except ImportError:
    _NUMPY_AVAILABLE = False
    np = None

try:
    from motif_patterns import APHILIC_DNA_PATTERNS
except ImportError: APHILIC_DNA_PATTERNS = {}

logger = logging.getLogger(__name__)

# TUNABLE PARAMETERS
MIN_SUM_LOG2 = 0.5  # Minimum sum_log2 for A-philic regions


class APhilicDetector(BaseMotifDetector):
    """A-philic DNA detector using 10-mer scoring table."""

    SCORE_REFERENCE = 'Gorin et al. 1995, Vinogradov et al. 2003'

    MIN_SUM_LOG2 = MIN_SUM_LOG2

    def get_motif_class_name(self) -> str: return "A-philic_DNA"

    def theoretical_min_score(self) -> float:
        """Minimum biologically valid A-philic raw score (MIN_SUM_LOG2)."""
        return self.MIN_SUM_LOG2

    def theoretical_max_score(self, sequence_length: int = None) -> float:
        """Highest possible sum_log2 score based on 10-mer table.

        Each position contributes max(TENMER_LOG2.values()) / 10 per base.
        """
        if sequence_length is None:
            # Fallback: assume a minimal A-philic region of ~10 bp (one 10-mer)
            sequence_length = 10
        max_log2 = max(TENMER_LOG2.values())
        return (max_log2 / 10.0) * sequence_length

    def get_length_cap(self, subclass: str = None) -> int:
        """A-philic DNA stable up to ~300 bp (Vinogradov 2003)."""
        return 300

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        return {"a_philic_10mers": [(r"", "APH_10MER", "A-philic 10-mer table", "A-philic DNA", 10, "a_philic_10mer_score", 0.9, "A-philic 10mer motif", "user_table")]}

    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        seq = sequence.upper(); merged_regions = self._find_and_merge_10mer_matches(seq)
        if not merged_regions: return 0.0
        contrib = self._build_per_base_contrib(seq)
        return float(sum(sum(contrib[s:e]) for s, e in merged_regions))

    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        """Return merged A-philic region annotations."""
        seq = sequence.upper(); matches = self._find_10mer_matches(seq)
        if not matches: return []
        merged = self._merge_matches(matches); contrib = self._build_per_base_contrib(seq); annotations = []
        for region in merged:
            s, e, region_matches = region; sum_log2 = sum(contrib[s:e]); n_10 = len(region_matches)
            mean_per10 = (sum(m[2] for m in region_matches) / n_10) if n_10 > 0 else 0.0
            annotations.append({"start": s, "end": e, "length": e - s, "sum_log2": round(sum_log2, 6),
                                "mean_log2_per10mer": round(mean_per10, 6), "n_10mers": n_10,
                                "contributing_10mers": [{"tenmer": m[1], "start": m[0], "log2": m[2]} for m in region_matches]})
        return annotations

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Detect A-philic regions using 10-mer scoring."""
        sequence = sequence.upper().strip(); motifs = []; annotations = self.annotate_sequence(sequence)
        for i, region in enumerate(annotations):
            if region.get('sum_log2', 0) > self.MIN_SUM_LOG2 and region.get('n_10mers', 0) >= 1:
                start_pos, end_pos = region['start'], region['end']
                motif_seq = sequence[start_pos:end_pos]
                
                gc_content = round(calc_gc_content(motif_seq), 2)
                
                at_content = round(calc_at_content(motif_seq), 2)
                
                contributing_10mers = region.get('contributing_10mers', [])
                tenmer_list = ','.join([c['tenmer'] for c in contributing_10mers[:10]])  # First 10
                if len(contributing_10mers) > 10:
                    tenmer_list += '...'
                
                raw_score = region['sum_log2']
                canonical_class, canonical_subclass = normalize_class_subclass(self.get_motif_class_name(), 'A-philic DNA', strict=False, auto_correct=True)
                normalized_score = self.normalize_score(raw_score, region['length'], canonical_subclass)
                motifs.append({
                    'ID': f"{sequence_name}_APHIL_{start_pos+1}", 'Sequence_Name': sequence_name, 'Class': canonical_class,
                    'Subclass': canonical_subclass, 'Start': start_pos + 1, 'End': end_pos, 'Length': region['length'],
                    'Sequence': motif_seq, 'Raw_Score': round(raw_score, 3), 'Score': normalized_score, 'Strand': '+',
                    'Method': 'A-philic_detection', 'Pattern_ID': f'APHIL_{i+1}',
                    'Contributing_10mers': region.get('n_10mers', 0),
                    'Mean_10mer_Log2': round(region.get('mean_log2_per10mer', 0), 3),
                    'GC_Content': gc_content,
                    'AT_Content': at_content,
                    'Arm_Length': 'N/A',
                    'Loop_Length': 'N/A',
                    'Type_Of_Repeat': self._classify_aphilic_type(motif_seq, at_content),
                    'Criterion': self._get_aphilic_criterion(region),
                    'Disease_Relevance': self._get_aphilic_disease_relevance(at_content, region['length']),
                    'Regions_Involved': f'{region["n_10mers"]} overlapping A-philic 10-mers: {tenmer_list}'
                })
        return motifs
    
    def _classify_aphilic_type(self, sequence: str, at_content: float) -> str:
        """Classify A-philic DNA type"""
        if at_content > 80:
            return 'High AT-content A-philic DNA (>80% AT)'
        elif at_content > 65:
            return 'Moderate AT-content A-philic DNA (65-80% AT)'
        else:
            return 'A-philic DNA'
    
    def _get_aphilic_criterion(self, region: Dict[str, Any]) -> str:
        """Explain A-philic classification criterion"""
        criteria = []
        criteria.append('10-mer propensity scoring (Vinogradov 2003)')
        criteria.append(f'Sum log2 score {region.get("sum_log2", 0):.2f} >{self.MIN_SUM_LOG2}')
        criteria.append(f'{region.get("n_10mers", 0)} contributing 10-mers')
        return '; '.join(criteria)
    
    def _get_aphilic_disease_relevance(self, at_content: float, length: int) -> str:
        """Annotate disease relevance for A-philic DNA"""
        disease_notes = []
        
        if at_content > 75:
            disease_notes.append('High AT-content - DNA flexibility, nucleosome positioning, chromatin structure')
        
        if length > 100:
            disease_notes.append('Extended A-philic region - potential regulatory element, transcription factor binding')
        
        disease_notes.append('A-philic DNA - minor groove narrowing, protein-DNA recognition, gene regulation')
        
        return '; '.join(disease_notes)

    def _find_10mer_matches(self, seq: str) -> List[Tuple[int, str, float]]:
        """Find all exact 10-mer matches. Uses Hyperscan if available, else pure-Python."""
        if _HYPERSCAN_AVAILABLE:
            try: return self._hs_find_matches(seq)
            except Exception as e:
                logger.warning(f"Hyperscan failed for A-philic, falling back: {e}"); return self._py_find_matches(seq)
        return self._py_find_matches(seq)

    def _hs_find_matches(self, seq: str) -> List[Tuple[int, str, float]]:
        try: return hyperscan_backend.hs_find_matches(seq, TENMER_LOG2)
        except Exception as e: logger.error(f"Hyperscan matching failed: {e}"); raise

    def _py_find_matches(self, seq: str) -> List[Tuple[int, str, float]]:
        try: return hyperscan_backend.py_find_matches(seq, TENMER_LOG2)
        except Exception as e: logger.warning(f"Optimized matching failed: {e}"); return self._py_find_matches_loop(seq)

    def _py_find_matches_loop(self, seq: str) -> List[Tuple[int, str, float]]:
        n = len(seq); matches: List[Tuple[int, str, float]] = []
        for i in range(0, n - 10 + 1):
            ten = seq[i:i + 10]; log2 = TENMER_LOG2.get(ten)
            if log2 is not None: matches.append((i, ten, float(log2)))
        return matches

    def _merge_matches(self, matches: List[Tuple[int, str, float]], merge_gap: int = 0) -> List[Tuple[int, int, List[Tuple[int, str, float]]]]:
        """Merge overlapping/adjacent 10-mer matches."""
        if not matches: return []
        merged = []; cur_start, cur_end = matches[0][0], matches[0][0] + 10; cur_matches = [matches[0]]
        for m in matches[1:]:
            s = m[0]; m_end = s + 10
            if s <= cur_end + merge_gap: cur_end = max(cur_end, m_end); cur_matches.append(m)
            else: merged.append((cur_start, cur_end, cur_matches)); cur_start, cur_end = s, m_end; cur_matches = [m]
        merged.append((cur_start, cur_end, cur_matches))
        return merged

    def _find_and_merge_10mer_matches(self, seq: str, merge_gap: int = 0) -> List[Tuple[int, int]]:
        matches = self._find_10mer_matches(seq); merged = self._merge_matches(matches, merge_gap=merge_gap)
        return [(s, e) for (s, e, _) in merged]

    def _build_per_base_contrib(self, seq: str) -> List[float]:
        """Build per-base contribution array from 10-mer matches."""
        n = len(seq)
        matches = self._find_10mer_matches(seq)
        
        if _NUMPY_AVAILABLE and n > 1000 and len(matches) > 0:
            try:
                contrib = np.zeros(n, dtype=np.float64)
                per_base_val = 1.0 / 10.0
                
                for (start, ten, log2) in matches:
                    per_base = float(log2) * per_base_val
                    end = min(start + 10, n)
                    contrib[start:end] += per_base
                
                return contrib.tolist()
            except Exception:
                pass
        
        contrib = [0.0] * n
        for (start, ten, log2) in matches:
            per_base = float(log2) / 10.0
            for k in range(start, min(start + 10, n)):
                contrib[k] += per_base
        return contrib
