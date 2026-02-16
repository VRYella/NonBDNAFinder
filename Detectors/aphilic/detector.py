"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ A-philic DNA Detector - 10-mer scoring table with hyperscan acceleration     │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘
"""
# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
import logging
from typing import Any, Dict, List, Tuple

try: from Detectors.base.base_detector import BaseMotifDetector
except ImportError:
    import sys; from pathlib import Path
    parent_dir = str(Path(__file__).parent.parent.parent)
    if parent_dir not in sys.path: sys.path.insert(0, parent_dir)
    from Detectors import BaseMotifDetector

from Utilities.core.motif_normalizer import normalize_class_subclass
from .tenmer_table import TENMER_LOG2

try: from Detectors.zdna import hyperscan_backend; _HYPERSCAN_AVAILABLE = hyperscan_backend.is_hyperscan_available()
except ImportError: _HYPERSCAN_AVAILABLE = False; hyperscan_backend = None

try: import numpy as np; _NUMPY_AVAILABLE = True
except ImportError: _NUMPY_AVAILABLE = False; np = None

try: from motif_patterns import APHILIC_DNA_PATTERNS
except ImportError: APHILIC_DNA_PATTERNS = {}

logger = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════════
# TUNABLE PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════════
MIN_SUM_LOG2 = 0.5  # Minimum sum_log2 for A-philic regions
# ═══════════════════════════════════════════════════════════════════════════════


class APhilicDetector(BaseMotifDetector):
    """A-philic DNA detector using 10-mer scoring table."""

    MIN_SUM_LOG2 = MIN_SUM_LOG2

    def get_motif_class_name(self) -> str: return "A-philic_DNA"

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
                
                # Calculate GC content
                gc_count = motif_seq.count('G') + motif_seq.count('C')
                gc_content = round(gc_count / len(motif_seq) * 100, 2) if len(motif_seq) > 0 else 0.0
                
                # Calculate AT content
                at_count = motif_seq.count('A') + motif_seq.count('T')
                at_content = round(at_count / len(motif_seq) * 100, 2) if len(motif_seq) > 0 else 0.0
                
                # Get contributing 10-mers summary
                contributing_10mers = region.get('contributing_10mers', [])
                tenmer_list = ','.join([c['tenmer'] for c in contributing_10mers[:10]])  # First 10
                # Add ellipsis only if there are more than 10
                if len(contributing_10mers) > 10:
                    tenmer_list += '...'
                
                canonical_class, canonical_subclass = normalize_class_subclass(self.get_motif_class_name(), 'A-philic DNA', strict=False, auto_correct=True)
                motifs.append({
                    'ID': f"{sequence_name}_APHIL_{start_pos+1}", 'Sequence_Name': sequence_name, 'Class': canonical_class,
                    'Subclass': canonical_subclass, 'Start': start_pos + 1, 'End': end_pos, 'Length': region['length'],
                    'Sequence': motif_seq, 'Score': round(region['sum_log2'], 3), 'Strand': '+',
                    'Method': 'A-philic_detection', 'Pattern_ID': f'APHIL_{i+1}',
                    'Contributing_10mers': region.get('n_10mers', 0),
                    'Mean_10mer_Log2': round(region.get('mean_log2_per10mer', 0), 3),
                    'GC_Content': gc_content,
                    'AT_Content': at_content,
                    'Arm_Length': 'N/A',  # Not applicable for A-philic regions
                    'Loop_Length': 'N/A',  # Not applicable for A-philic regions
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
        
        # High AT content - DNA flexibility, nucleosome positioning
        if at_content > 75:
            disease_notes.append('High AT-content - DNA flexibility, nucleosome positioning, chromatin structure')
        
        # Long A-philic regions
        if length > 100:
            disease_notes.append('Extended A-philic region - potential regulatory element, transcription factor binding')
        
        # General associations
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
        """Build per-base contribution array. Each 10-mer's log2 distributed as L/10 to each base.
        
        Uses NumPy vectorization for 2-3x speedup on large sequences when available.
        """
        n = len(seq)
        matches = self._find_10mer_matches(seq)
        
        # Use vectorized version for large sequences if NumPy is available
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
                # Fallback to loop-based implementation on error
                pass
        
        # Original loop-based implementation
        contrib = [0.0] * n
        for (start, ten, log2) in matches:
            per_base = float(log2) / 10.0
            for k in range(start, min(start + 10, n)):
                contrib[k] += per_base
        return contrib
