"""A-philic DNA detector implementation using 10-mer scoring table.

This module implements the APhilicDetector class for identifying A-philic DNA regions
in DNA sequences. A-philic DNA represents sequences with unusual structural properties
related to A-DNA formation propensity.

Detection method:
- Uses a 10-mer scoring table (TENMER_LOG2) to identify high-scoring regions
- Merges overlapping/adjacent 10-mer matches to avoid duplicate reporting
- Scores are computed by redistributing each 10-mer's log2 value equally across its 10 bases

Performance:
- Uses Hyperscan when available for fast pattern matching
- Falls back to pure Python implementation when Hyperscan is unavailable
"""

import logging
from typing import Any, Dict, List, Tuple

# Import base detector from the base module
try:
    from detectors.base.base_detector import BaseMotifDetector
except ImportError:
    # Fallback for when running from detectors.py directly
    import sys
    from pathlib import Path
    # Add parent directory to path to import from detectors module
    parent_dir = str(Path(__file__).parent.parent.parent)
    if parent_dir not in sys.path:
        sys.path.insert(0, parent_dir)
    from detectors import BaseMotifDetector

from core.motif_normalizer import normalize_class_subclass

from .tenmer_table import TENMER_LOG2

# Import Hyperscan backend (reusing Z-DNA backend which has the same interface)
try:
    from detectors.zdna import hyperscan_backend
    _HYPERSCAN_AVAILABLE = hyperscan_backend.is_hyperscan_available()
except ImportError:
    _HYPERSCAN_AVAILABLE = False
    hyperscan_backend = None

# Import pattern definitions
try:
    from motif_patterns import APHILIC_DNA_PATTERNS
except ImportError:
    APHILIC_DNA_PATTERNS = {}

logger = logging.getLogger(__name__)


class APhilicDetector(BaseMotifDetector):
    """Detector for A-philic DNA motifs using a 10-mer scoring table."""

    def get_motif_class_name(self) -> str:
        return "A-philic_DNA"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return synthetic pattern for A-philic 10-mer table detection."""
        return {
            "a_philic_10mers": [
                (r"", "APH_10MER", "A-philic 10-mer table", "A-philic DNA",
                 10, "a_philic_10mer_score", 0.9, "A-philic 10mer motif", "user_table"),
            ]
        }

    # -------------------------
    # Public API: calculate_score used by framework
    # -------------------------
    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        """Calculate raw sum_log2 for merged A-philic regions. Redistributes each 10-mer's avg_log2 equally over 10 bases, sums per-base values in merged regions."""
        seq = sequence.upper()
        merged_regions = self._find_and_merge_10mer_matches(seq)
        if not merged_regions:
            return 0.0
        contrib = self._build_per_base_contrib(seq)
        total_sum = 0.0
        for s, e in merged_regions:
            total_sum += sum(contrib[s:e])
        return float(total_sum)

    # -------------------------
    # Helper: return detailed annotation
    # -------------------------
    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        """Return merged A-philic region annotations. MERGING GUARANTEE: Always merges overlapping/adjacent 10-mers. Returns {start,end,length,sum_log2,mean_log2_per10mer,n_10mers,contributing_10mers} dicts."""
        seq = sequence.upper()
        
        matches = self._find_10mer_matches(seq)
        if not matches:
            return []
        
        # This is the critical step that ensures no duplicate/split reporting
        merged = self._merge_matches(matches)
        
        contrib = self._build_per_base_contrib(seq)
        
        annotations = []
        for region in merged:
            s, e, region_matches = region
            # Sum contributions across the merged region
            sum_log2 = sum(contrib[s:e])
            n_10 = len(region_matches)
            mean_per10 = (sum(m[2] for m in region_matches) / n_10) if n_10 > 0 else 0.0
            
            ann = {
                "start": s,
                "end": e,
                "length": e - s,
                "sum_log2": round(sum_log2, 6),
                "mean_log2_per10mer": round(mean_per10, 6),
                "n_10mers": n_10,
                "contributing_10mers": [{"tenmer": m[1], "start": m[0], "log2": m[2]} for m in region_matches]
            }
            annotations.append(ann)
        return annotations

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Detect A-philic regions using 10-mer scoring. ALWAYS outputs merged regions via annotate_sequence(), ensuring no duplicate reporting."""
        sequence = sequence.upper().strip()
        motifs = []
        
        # Use the annotation method to find A-philic regions.
        # This GUARANTEES that overlapping/adjacent 10-mer matches are merged.
        annotations = self.annotate_sequence(sequence)
        
        for i, region in enumerate(annotations):
            if region.get('sum_log2', 0) > 0.5 and region.get('n_10mers', 0) >= 1:
                start_pos = region['start']
                end_pos = region['end']
                
                # Normalize class/subclass using canonical taxonomy
                canonical_class, canonical_subclass = normalize_class_subclass(
                    self.get_motif_class_name(),
                    'A-philic DNA',
                    strict=False,
                    auto_correct=True
                )
                
                motifs.append({
                    'ID': f"{sequence_name}_APHIL_{start_pos+1}",
                    'Sequence_Name': sequence_name,
                    'Class': canonical_class,
                    'Subclass': canonical_subclass,
                    'Start': start_pos + 1,  # 1-based coordinates
                    'End': end_pos,
                    'Length': region['length'],
                    'Sequence': sequence[start_pos:end_pos],
                    'Score': round(region['sum_log2'], 3),
                    'Strand': '+',
                    'Method': 'A-philic_detection',
                    'Pattern_ID': f'APHIL_{i+1}'
                })
        
        return motifs

    # -------------------------
    # Core match / merge / contrib helpers
    # -------------------------
    def _find_10mer_matches(self, seq: str) -> List[Tuple[int, str, float]]:
        """Find all exact 10-mer matches. Returns (start, tenmer, avg_log2) tuples. Uses Hyperscan if available, else pure-Python. Merging happens in _merge_matches()."""
        if _HYPERSCAN_AVAILABLE:
            try:
                return self._hs_find_matches(seq)
            except Exception as e:
                logger.warning(f"Hyperscan matching failed for A-philic DNA, falling back to pure Python: {e}")
                return self._py_find_matches(seq)
        else:
            return self._py_find_matches(seq)

    def _hs_find_matches(self, seq: str) -> List[Tuple[int, str, float]]:
        """Use Hyperscan compiled database for fast 10-mer matching. Returns (start, tenmer, log2) sorted by start. Improved error handling and logging."""
        try:
            # Use the shared hyperscan backend
            matches = hyperscan_backend.hs_find_matches(seq, TENMER_LOG2)
            logger.debug(f"Hyperscan scan completed: {len(matches)} A-philic 10-mer matches found")
            return matches
            
        except Exception as e:
            logger.error(f"Hyperscan matching failed for A-philic DNA: {e}. Falling back to pure Python.")
            raise  # Re-raise to trigger fallback in _find_10mer_matches

    def _py_find_matches(self, seq: str) -> List[Tuple[int, str, float]]:
        """Pure-Python exact search with vectorization when available.
        
        Uses vectorized matching when NumPy is available for improved performance.
        Falls back to loop-based implementation otherwise.
        
        Returns:
            List of (start, tenmer, log2) tuples
        """
        # Delegate to shared vectorized implementation from hyperscan_backend
        # This provides automatic optimization when NumPy is available
        try:
            return hyperscan_backend.py_find_matches(seq, TENMER_LOG2)
        except Exception as e:
            logger.warning(f"Optimized matching failed, using loop-based fallback: {e}")
            return self._py_find_matches_loop(seq)
    
    def _py_find_matches_loop(self, seq: str) -> List[Tuple[int, str, float]]:
        """Original loop-based implementation (kept for validation and fallback)."""
        n = len(seq)
        matches: List[Tuple[int, str, float]] = []
        for i in range(0, n - 10 + 1):
            ten = seq[i:i + 10]
            log2 = TENMER_LOG2.get(ten)
            if log2 is not None:
                matches.append((i, ten, float(log2)))
        return matches

    def _merge_matches(self, matches: List[Tuple[int, str, float]],
                       merge_gap: int = 0) -> List[Tuple[int, int, List[Tuple[int, str, float]]]]:
        """Merge overlapping/adjacent 10-mer matches. Args: matches sorted by start, merge_gap (default 0=adjacent only). Returns (region_start, region_end, list_of_matches) tuples."""
        if not matches:
            return []
        
        merged = []
        cur_start, cur_end = matches[0][0], matches[0][0] + 10
        cur_matches = [matches[0]]
        
        for m in matches[1:]:
            s = m[0]
            m_end = s + 10
            
            if s <= cur_end + merge_gap:
                # Extend current region and add match
                cur_end = max(cur_end, m_end)
                cur_matches.append(m)
            else:
                # Gap too large; finalize current region and start new one
                merged.append((cur_start, cur_end, cur_matches))
                cur_start, cur_end = s, m_end
                cur_matches = [m]
        
        # Don't forget to add the last region
        merged.append((cur_start, cur_end, cur_matches))
        return merged

    def _find_and_merge_10mer_matches(self, seq: str, merge_gap: int = 0) -> List[Tuple[int, int]]:
        matches = self._find_10mer_matches(seq)
        merged = self._merge_matches(matches, merge_gap=merge_gap)
        return [(s, e) for (s, e, _) in merged]

    def _build_per_base_contrib(self, seq: str) -> List[float]:
        """Build per-base contribution array. Each 10-mer's log2 L is distributed as L/10 to each of its 10 bases. Returns list where contrib[i] = total from all 10-mers covering position i."""
        n = len(seq)
        contrib = [0.0] * n
        matches = self._find_10mer_matches(seq)
        
        # Distribute each 10-mer's log2 value across its 10 bases
        for (start, ten, log2) in matches:
            per_base = float(log2) / 10.0
            for k in range(start, min(start + 10, n)):
                contrib[k] += per_base
        
        return contrib
