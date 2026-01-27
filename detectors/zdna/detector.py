"""Z-DNA detector implementation using 10-mer scoring and eGZ-motif patterns.

This module implements the ZDNADetector class for identifying Z-DNA forming regions
in DNA sequences. Z-DNA is a left-handed double helix structure that can form in
sequences with alternating purine-pyrimidine repeats, particularly CG repeats.

Detection methods:
1. 10-mer scoring table (Ho et al. 1986) - thermodynamic prediction
2. eGZ-motifs (Herbert 1997) - extruded-G Z-DNA patterns

References:
    - Ho PS, Ellison MJ, Quigley GJ, Rich A. (1986). A computer aided thermodynamic 
      approach for predicting the formation of Z-DNA in naturally occurring sequences.
      EMBO J. 5(10):2737-44. PMID: 3780660
    - Herbert A, Alfken J, Kim YG, Mian IS, Nishikura K, Rich A. (1997). 
      A Z-DNA binding domain present in the human editing enzyme, double-stranded 
      RNA adenosine deaminase. Proc Natl Acad Sci U S A. 94(16):8421-6. PMID: 9237994
"""

import logging
import re
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

from .tenmer_table import TENMER_SCORE
from . import hyperscan_backend
from core.motif_normalizer import normalize_class_subclass

# Import pattern definitions
try:
    from motif_patterns import ZDNA_PATTERNS
except ImportError:
    ZDNA_PATTERNS = {}

logger = logging.getLogger(__name__)


class ZDNADetector(BaseMotifDetector):
    """Z-DNA (left-handed helix) detector using 10-mer scoring (Ho et al. 1986)
    and eGZ-motif patterns (Herbert 1997)."""

    # eGZ-motif detection constants
    MIN_EGZ_REPEATS = 3
    EGZ_BASE_SCORE = 0.85
    EGZ_MIN_SCORE_THRESHOLD = 0.80

    def get_motif_class_name(self) -> str:
        return "Z-DNA"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return patterns for Z-DNA detection (10-mer table & eGZ-motifs)."""
        def fallback():
            return {
                "z_dna_10mers": [
                    (r"", "ZDN_10MER", "Z-DNA 10-mer table", "Z-DNA", 10, "z_dna_10mer_score", 0.9, "Z-DNA 10mer motif", "user_table"),
                ],
                "egz_motifs": [
                    (r"(?:CGG){3,}", "ZDN_EGZ_CGG", "CGG repeat (eGZ)", "eGZ", 9, "egz_score", self.EGZ_BASE_SCORE, "Extruded-G Z-DNA CGG repeat", "Herbert 1997"),
                    (r"(?:GGC){3,}", "ZDN_EGZ_GGC", "GGC repeat (eGZ)", "eGZ", 9, "egz_score", self.EGZ_BASE_SCORE, "Extruded-G Z-DNA GGC repeat", "Herbert 1997"),
                    (r"(?:CCG){3,}", "ZDN_EGZ_CCG", "CCG repeat (eGZ)", "eGZ", 9, "egz_score", self.EGZ_BASE_SCORE, "Extruded-G Z-DNA CCG repeat", "Herbert 1997"),
                    (r"(?:GCC){3,}", "ZDN_EGZ_GCC", "GCC repeat (eGZ)", "eGZ", 9, "egz_score", self.EGZ_BASE_SCORE, "Extruded-G Z-DNA GCC repeat", "Herbert 1997"),
                ]
            }
        
        patterns = self._load_patterns(ZDNA_PATTERNS, fallback)
        if 'egz_motifs' in patterns:
            patterns['egz_motifs'] = [
                tuple([*list(p)[:6], self.EGZ_BASE_SCORE, *list(p)[7:]])
                for p in patterns['egz_motifs']
            ]
        return patterns

    # -------------------------
    # Public API
    # -------------------------
    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        """Return total sum_score across merged Z-like regions. Score computed by redistributing each 10-mer's score equally across its 10 bases."""
        seq = sequence.upper()
        merged = self._find_and_merge_10mer_matches(seq)
        if not merged:
            return 0.0
        contrib = self._build_per_base_contrib(seq)
        total = 0.0
        for s, e in merged:
            total += sum(contrib[s:e])
        return float(total)

    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        """Return list of merged Z-DNA and eGZ-motif regions. MERGING GUARANTEE: Always merges overlapping/adjacent 10-mer matches. Returns {start,end,length,sum_score,subclass,pattern_id} dicts."""
        seq = sequence.upper()
        annotations = []
        
        matches = self._find_10mer_matches(seq)
        if matches:
            # Merge overlapping/adjacent matches into regions
            merged = self._merge_matches(matches)
            
            contrib = self._build_per_base_contrib(seq)
            
            for (s, e, region_matches) in merged:
                # Sum contributions across the merged region
                sum_score = sum(contrib[s:e])
                n10 = len(region_matches)
                mean10 = (sum(m[2] for m in region_matches) / n10) if n10 > 0 else 0.0
                
                ann = {
                    "start": s,
                    "end": e,
                    "length": e - s,
                    "sum_score": round(sum_score, 6),
                    "mean_score_per10mer": round(mean10, 6),
                    "n_10mers": n10,
                    "contributing_10mers": [{"tenmer": m[1], "start": m[0], "score": m[2]} for m in region_matches],
                    "subclass": "Z-DNA",
                    "pattern_id": "ZDN_10MER"
                }
                annotations.append(ann)
        
        egz_annotations = self._find_egz_motifs(seq)
        annotations.extend(egz_annotations)
        
        annotations.sort(key=lambda x: x['start'])
        
        return annotations

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Detect Z-DNA (10-mer scoring) and eGZ-motif (regex) regions. ALWAYS outputs merged regions via annotate_sequence(), ensuring no duplicate reporting."""
        sequence = sequence.upper().strip()
        motifs = []
        
        # Use the annotation method to find Z-DNA regions.
        # This GUARANTEES that overlapping/adjacent 10-mer matches are merged.
        annotations = self.annotate_sequence(sequence)
        
        for i, region in enumerate(annotations):
            subclass = region.get('subclass', 'Z-DNA')
            
            # Handle eGZ-motif regions
            if subclass == 'eGZ':
                if region.get('sum_score', 0) >= self.EGZ_MIN_SCORE_THRESHOLD:
                    start_pos = region['start']
                    end_pos = region['end']
                    motif_seq = sequence[start_pos:end_pos]
                    
                    gc_content = self._calc_gc(motif_seq)
                    
                    # Normalize class/subclass using canonical taxonomy
                    canonical_class, canonical_subclass = normalize_class_subclass(
                        self.get_motif_class_name(),
                        'eGZ',
                        strict=False,
                        auto_correct=True
                    )
                    
                    motifs.append({
                        'ID': f"{sequence_name}_{region['pattern_id']}_{start_pos+1}",
                        'Sequence_Name': sequence_name,
                        'Class': canonical_class,
                        'Subclass': canonical_subclass,
                        'Start': start_pos + 1,  # 1-based coordinates
                        'End': end_pos,
                        'Length': region['length'],
                        'Sequence': motif_seq,
                        'Score': round(region['sum_score'], 3),
                        'Strand': '+',
                        'Method': 'Z-DNA_detection',
                        'Pattern_ID': region['pattern_id'],
                        # eGZ-specific details
                        'Repeat_Unit': region.get('repeat_unit', ''),
                        'Repeat_Count': region.get('repeat_count', 0),
                        'GC_Content': round(gc_content, 2)
                    })
            else:
                # Handle classic Z-DNA regions (from 10-mer scoring)
                if region.get('sum_score', 0) > 50.0 and region.get('n_10mers', 0) >= 1:
                    start_pos = region['start']
                    end_pos = region['end']
                    motif_seq = sequence[start_pos:end_pos]
                    
                    cg_count = motif_seq.count('CG') + motif_seq.count('GC')
                    at_count = motif_seq.count('AT') + motif_seq.count('TA')
                    
                    gc_content = self._calc_gc(motif_seq)
                    
                    alternating_cg = len(re.findall(r'(?:CG){2,}', motif_seq)) + len(re.findall(r'(?:GC){2,}', motif_seq))
                    alternating_at = len(re.findall(r'(?:AT){2,}', motif_seq)) + len(re.findall(r'(?:TA){2,}', motif_seq))
                    
                    # Normalize class/subclass using canonical taxonomy
                    canonical_class, canonical_subclass = normalize_class_subclass(
                        self.get_motif_class_name(),
                        'Z-DNA',
                        strict=False,
                        auto_correct=True
                    )
                    
                    motifs.append({
                        'ID': f"{sequence_name}_ZDNA_{start_pos+1}",
                        'Sequence_Name': sequence_name,
                        'Class': canonical_class,
                        'Subclass': canonical_subclass,
                        'Start': start_pos + 1,  # 1-based coordinates
                        'End': end_pos,
                        'Length': region['length'],
                        'Sequence': motif_seq,
                        'Score': round(region['sum_score'], 3),
                        'Strand': '+',
                        'Method': 'Z-DNA_detection',
                        'Pattern_ID': f'ZDNA_{i+1}',
                        # Component details
                        'Contributing_10mers': region.get('n_10mers', 0),
                        'Mean_10mer_Score': region.get('mean_score_per10mer', 0),
                        'CG_Dinucleotides': cg_count,
                        'AT_Dinucleotides': at_count,
                        'Alternating_CG_Regions': alternating_cg,
                        'Alternating_AT_Regions': alternating_at,
                        'GC_Content': round(gc_content, 2)
                    })
        
        return motifs

    # -------------------------
    # Core helpers
    # -------------------------
    def _find_egz_motifs(self, seq: str) -> List[Dict[str, Any]]:
        """Find eGZ-motif (Extruded-G Z-DNA) using regex. Searches (CGG/GGC/CCG/GCC)n repeats >= MIN_EGZ_REPEATS. Returns {start,end,length,sum_score,subclass,pattern_id,repeat_unit,repeat_count} dicts."""
        patterns = self.get_patterns().get('egz_motifs', [])
        annotations = []
        
        for pattern_info in patterns:
            regex, pattern_id, name, subclass, min_len, score_type, threshold, desc, ref = pattern_info
            
            for match in re.finditer(regex, seq, re.IGNORECASE):
                start, end = match.span()
                matched_seq = seq[start:end]
                
                # Determine repeat unit (first 3 characters of matched sequence)
                repeat_unit = matched_seq[:3].upper()
                repeat_count = len(matched_seq) // 3
                
                # Calculate score based on repeat count
                # Score increases with more repeats: base_score * (repeat_count / MIN_EGZ_REPEATS)
                # Minimum MIN_EGZ_REPEATS repeats gives score = base_score, more repeats give higher scores
                repeat_score = threshold * (repeat_count / float(self.MIN_EGZ_REPEATS))
                
                ann = {
                    "start": start,
                    "end": end,
                    "length": end - start,
                    "sum_score": round(repeat_score, 6),
                    "subclass": "eGZ",
                    "pattern_id": pattern_id,
                    "repeat_unit": repeat_unit,
                    "repeat_count": repeat_count,
                    "description": desc,
                    "reference": ref
                }
                annotations.append(ann)
        
        return annotations
    
    def _find_10mer_matches(self, seq: str) -> List[Tuple[int, str, float]]:
        """Find all exact 10-mer matches. Returns (start, tenmer, score) tuples. Uses Hyperscan if available, else pure-Python. Merging happens in _merge_matches()."""
        if hyperscan_backend.is_hyperscan_available():
            try:
                return hyperscan_backend.hs_find_matches(seq, TENMER_SCORE)
            except Exception as e:
                logger.warning(f"Hyperscan matching failed for Z-DNA, falling back to pure Python: {e}")
                return hyperscan_backend.py_find_matches(seq, TENMER_SCORE)
        else:
            return hyperscan_backend.py_find_matches(seq, TENMER_SCORE)

    def _merge_matches(self, matches: List[Tuple[int, str, float]],
                       merge_gap: int = 0) -> List[Tuple[int, int, List[Tuple[int, str, float]]]]:
        """Merge overlapping/adjacent 10-mer matches. Args: matches sorted by start, merge_gap (default 0=adjacent only). Returns (region_start, region_end, list_of_matches) tuples."""
        if not matches:
            return []
        
        merged = []
        cur_start = matches[0][0]
        cur_end = matches[0][0] + 10
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
        """Build per-base contribution array. Each 10-mer's score S is distributed as S/10 to each of its 10 bases. Returns list where contrib[i] = total from all 10-mers covering position i."""
        n = len(seq)
        contrib = [0.0] * n
        matches = self._find_10mer_matches(seq)
        
        # Distribute each 10-mer's score across its 10 bases
        for (start, ten, score) in matches:
            per_base = float(score) / 10.0
            for k in range(start, min(start + 10, n)):
                contrib[k] += per_base
        
        return contrib
