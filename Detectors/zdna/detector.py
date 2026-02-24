"""Z-DNA detector using 10-mer scoring (Ho 1986) and eGZ-motifs (Herbert 1997)."""
# IMPORTS
import logging
import re
from typing import Any, Dict, List, Tuple

from Detectors.base.base_detector import BaseMotifDetector

from .tenmer_table import TENMER_SCORE
from . import hyperscan_backend
from Utilities.core.motif_normalizer import normalize_class_subclass

try:
    from motif_patterns import ZDNA_PATTERNS
except ImportError:
    ZDNA_PATTERNS = {}

try:
    import numpy as np
    _NUMPY_AVAILABLE = True
except ImportError:
    _NUMPY_AVAILABLE = False

logger = logging.getLogger(__name__)

# TUNABLE PARAMETERS
MIN_EGZ_REPEATS = 3; EGZ_BASE_SCORE = 0.85; EGZ_MIN_SCORE_THRESHOLD = 0.80; MIN_Z_SCORE = 50.0


class ZDNADetector(BaseMotifDetector):
    """Z-DNA detector: 10-mer scoring (Ho 1986) + eGZ-motifs (Herbert 1997)."""

    SCORE_REFERENCE = 'Ho et al. 1986 (EMBO J) - Z-DNA propensity'

    MIN_EGZ_REPEATS = MIN_EGZ_REPEATS; EGZ_BASE_SCORE = EGZ_BASE_SCORE; EGZ_MIN_SCORE_THRESHOLD = EGZ_MIN_SCORE_THRESHOLD
    MIN_Z_SCORE = MIN_Z_SCORE

    def get_motif_class_name(self) -> str: return "Z-DNA"

    def get_length_cap(self, subclass: str = None) -> int:
        """Z-DNA structures stable up to ~300 bp (Ho 1986)."""
        return 300

    def theoretical_min_score(self) -> float:
        """Minimum biologically valid Z-DNA cumulative raw score (Ho 1986 10-mer threshold)."""
        return self.MIN_Z_SCORE

    def theoretical_max_score(self, sequence_length: int = None) -> float:
        """Highest possible Z-DNA cumulative raw score.

        Per-base max contribution = max(TENMER_SCORE) / 10.
        For a region of length L: max = (max_tenmer_score / 10) * L.
        Fallback uses a typical Z-DNA region length of 20 bp (two CpG dinucleotides × 10-mer).
        """
        max_tenmer = max(TENMER_SCORE.values())
        if sequence_length is None:
            # Fallback: typical short Z-DNA region (20 bp covers two overlapping 10-mers)
            sequence_length = 20
        return (max_tenmer / 10.0) * sequence_length

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return patterns for Z-DNA detection."""
        def fallback():
            return {
                "z_dna_10mers": [(r"", "ZDN_10MER", "Z-DNA 10-mer table", "Z-DNA", 10, "z_dna_10mer_score", 0.9, "Z-DNA 10mer motif", "user_table")],
                "egz_motifs": [
                    (r"(?:CGG){4,}", "ZDN_EGZ_CGG", "CGG repeat (eGZ)", "eGZ", 9, "egz_score", self.EGZ_BASE_SCORE, "Extruded-G Z-DNA CGG repeat", "Herbert 1997"),
                    (r"(?:GGC){4,}", "ZDN_EGZ_GGC", "GGC repeat (eGZ)", "eGZ", 9, "egz_score", self.EGZ_BASE_SCORE, "Extruded-G Z-DNA GGC repeat", "Herbert 1997"),
                    (r"(?:CCG){4,}", "ZDN_EGZ_CCG", "CCG repeat (eGZ)", "eGZ", 9, "egz_score", self.EGZ_BASE_SCORE, "Extruded-G Z-DNA CCG repeat", "Herbert 1997"),
                    (r"(?:GCC){4,}", "ZDN_EGZ_GCC", "GCC repeat (eGZ)", "eGZ", 9, "egz_score", self.EGZ_BASE_SCORE, "Extruded-G Z-DNA GCC repeat", "Herbert 1997"),
                ]
            }
        patterns = self._load_patterns(ZDNA_PATTERNS, fallback)
        if 'egz_motifs' in patterns:
            patterns['egz_motifs'] = [tuple([*list(p)[:6], self.EGZ_BASE_SCORE, *list(p)[7:]]) for p in patterns['egz_motifs']]
        return patterns

    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        """Return total sum_score across merged Z-like regions."""
        seq = sequence.upper(); merged = self._find_and_merge_10mer_matches(seq)
        if not merged: return 0.0
        contrib = self._build_per_base_contrib(seq)
        return float(sum(sum(contrib[s:e]) for s, e in merged))

    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        """Return list of merged Z-DNA and eGZ-motif regions."""
        seq = sequence.upper(); annotations = []
        matches = self._find_10mer_matches(seq)
        if matches:
            merged = self._merge_matches(matches); contrib = self._build_per_base_contrib(seq, matches=matches)
            for (s, e, region_matches) in merged:
                sum_score = sum(contrib[s:e]); n10 = len(region_matches)
                mean10 = (sum(m[2] for m in region_matches) / n10) if n10 > 0 else 0.0
                annotations.append({
                    "start": s, "end": e, "length": e - s, "sum_score": round(sum_score, 6),
                    "mean_score_per10mer": round(mean10, 6), "n_10mers": n10,
                    "contributing_10mers": [{"tenmer": m[1], "start": m[0], "score": m[2]} for m in region_matches],
                    "subclass": "Z-DNA", "pattern_id": "ZDN_10MER"
                })
        annotations.extend(self._find_egz_motifs(seq))
        annotations.sort(key=lambda x: x['start'])
        return annotations

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Detect Z-DNA (10-mer scoring) and eGZ-motif regions."""
        sequence = sequence.upper().strip(); motifs = []
        annotations = self.annotate_sequence(sequence)
        
        for i, region in enumerate(annotations):
            subclass = region.get('subclass', 'Z-DNA'); start_pos = region['start']; end_pos = region['end']
            motif_seq = sequence[start_pos:end_pos]; gc_content = self._calc_gc(motif_seq)
            
            if subclass == 'eGZ':
                if region.get('sum_score', 0) >= self.EGZ_MIN_SCORE_THRESHOLD:
                    canonical_class, canonical_subclass = normalize_class_subclass(self.get_motif_class_name(), 'eGZ', strict=False, auto_correct=True)
                    
                    repeat_unit = region.get('repeat_unit', '')
                    repeat_count = region.get('repeat_count', 0)
                    raw_score = region['sum_score']
                    normalized_score = self.normalize_score(raw_score, region['length'], canonical_subclass)
                    
                    motifs.append({
                        'ID': f"{sequence_name}_{region['pattern_id']}_{start_pos+1}", 'Sequence_Name': sequence_name,
                        'Class': canonical_class, 'Subclass': canonical_subclass,
                        'Start': start_pos + 1, 'End': end_pos, 'Length': region['length'], 'Sequence': motif_seq,
                        'Raw_Score': round(raw_score, 3),
                        'Score': normalized_score,
                        'Strand': '+', 'Method': 'Z-DNA_detection',
                        'Pattern_ID': region['pattern_id'], 
                        'Repeat_Unit': repeat_unit,
                        'Repeat_Count': repeat_count, 
                        'GC_Content': round(gc_content, 2),
                        'Arm_Length': 'N/A',
                        'Loop_Length': 'N/A',
                        'Type_Of_Repeat': f'Trinucleotide eGZ-motif ({repeat_unit})',
                        'Criterion': self._get_egz_criterion(repeat_unit, repeat_count),
                        'Disease_Relevance': self._get_egz_disease_relevance(repeat_unit, repeat_count),
                        'Regions_Involved': f'eGZ-motif: {repeat_unit} trinucleotide × {repeat_count} repeats'
                    })
            else:
                if region.get('sum_score', 0) > self.MIN_Z_SCORE and region.get('n_10mers', 0) >= 1:
                    cg_count = motif_seq.count('CG') + motif_seq.count('GC'); at_count = motif_seq.count('AT') + motif_seq.count('TA')
                    alternating_cg = len(re.findall(r'(?:CG){2,}', motif_seq)) + len(re.findall(r'(?:GC){2,}', motif_seq))
                    alternating_at = len(re.findall(r'(?:AT){2,}', motif_seq)) + len(re.findall(r'(?:TA){2,}', motif_seq))
                    
                    canonical_class, canonical_subclass = normalize_class_subclass(self.get_motif_class_name(), 'Z-DNA', strict=False, auto_correct=True)
                    
                    zdna_type = self._classify_zdna_type(motif_seq, alternating_cg, alternating_at)
                    raw_score = region['sum_score']
                    normalized_score = self.normalize_score(raw_score, region['length'], canonical_subclass)
                    
                    motifs.append({
                        'ID': f"{sequence_name}_ZDNA_{start_pos+1}", 'Sequence_Name': sequence_name,
                        'Class': canonical_class, 'Subclass': canonical_subclass,
                        'Start': start_pos + 1, 'End': end_pos, 'Length': region['length'], 'Sequence': motif_seq,
                        'Raw_Score': round(raw_score, 3),
                        'Score': normalized_score,
                        'Strand': '+', 'Method': 'Z-DNA_detection',
                        'Pattern_ID': f'ZDNA_{i+1}', 
                        'Contributing_10mers': region.get('n_10mers', 0),
                        'Mean_10mer_Score': region.get('mean_score_per10mer', 0), 
                        'CG_Dinucleotides': cg_count,
                        'AT_Dinucleotides': at_count, 
                        'Alternating_CG_Regions': alternating_cg,
                        'Alternating_AT_Regions': alternating_at, 
                        'GC_Content': round(gc_content, 2),
                        'Arm_Length': 'N/A',
                        'Loop_Length': 'N/A',
                        'Type_Of_Repeat': zdna_type,
                        'Criterion': self._get_zdna_criterion(region),
                        'Disease_Relevance': self._get_zdna_disease_relevance(motif_seq, gc_content, region['length']),
                        'Regions_Involved': self._describe_zdna_regions(cg_count, at_count, alternating_cg, alternating_at)
                    })
        return motifs
    
    def _classify_zdna_type(self, sequence: str, alternating_cg: int, alternating_at: int) -> str:
        """Classify Z-DNA structural type"""
        if alternating_cg >= 2:
            return 'CG/GC alternating purine-pyrimidine (canonical Z-DNA)'
        elif alternating_at >= 2:
            return 'AT/TA alternating purine-pyrimidine (Z-DNA)'
        else:
            return 'Mixed dinucleotide composition (Z-DNA)'
    
    def _get_egz_criterion(self, repeat_unit: str, repeat_count: int) -> str:
        """Explain eGZ classification criterion"""
        return f'eGZ-motif: {repeat_unit} trinucleotide repeat ≥{self.MIN_EGZ_REPEATS} copies (detected n={repeat_count}); Herbert 1997 extruded-G Z-DNA'
    
    def _get_zdna_criterion(self, region: Dict[str, Any]) -> str:
        """Explain Z-DNA classification criterion"""
        criteria = []
        criteria.append(f'10-mer table scoring (Ho 1986)')
        criteria.append(f'Sum score {region.get("sum_score", 0):.1f} >{self.MIN_Z_SCORE}')
        criteria.append(f'{region.get("n_10mers", 0)} contributing 10-mers')
        return '; '.join(criteria)
    
    def _get_egz_disease_relevance(self, repeat_unit: str, repeat_count: int) -> str:
        """Annotate disease relevance for eGZ motifs"""
        disease_notes = []
        
        if repeat_unit in ['CGG', 'CCG']:
            if repeat_count >= 200:
                disease_notes.append(f'Fragile X syndrome: pathogenic CGG/CCG expansion (n≥200, detected n={repeat_count})')
            elif repeat_count >= 55:
                disease_notes.append(f'Fragile X premutation (55≤n<200, detected n={repeat_count})')
            else:
                disease_notes.append(f'CGG/CCG repeat (n={repeat_count}, monitor for expansion)')
        
        disease_notes.append('Z-DNA formation - transcription regulation, genomic instability')
        
        return '; '.join(disease_notes)
    
    def _get_zdna_disease_relevance(self, sequence: str, gc_content: float, length: int) -> str:
        """Annotate disease relevance for Z-DNA"""
        disease_notes = []
        
        if re.search(r'(CG){5,}', sequence) or re.search(r'(GC){5,}', sequence):
            disease_notes.append('Long CG/GC alternating - potential methylation site, epigenetic regulation')
        
        if gc_content > 75:
            disease_notes.append('High GC Z-DNA - promoter element, gene regulation')
        
        if length > 50:
            disease_notes.append('Extended Z-DNA - chromatin structure, recombination hotspot')
        
        disease_notes.append('Z-DNA formation - immune response (ZBP1 binding), transcription, genome instability')
        
        return '; '.join(disease_notes) if disease_notes else 'Z-DNA formation potential'
    
    def _describe_zdna_regions(self, cg_count: int, at_count: int, alternating_cg: int, alternating_at: int) -> str:
        """Describe regions involved in Z-DNA formation"""
        parts = []
        
        if cg_count > 0:
            parts.append(f'{cg_count} CG/GC dinucleotides')
        if at_count > 0:
            parts.append(f'{at_count} AT/TA dinucleotides')
        if alternating_cg > 0:
            parts.append(f'{alternating_cg} alternating CG/GC regions')
        if alternating_at > 0:
            parts.append(f'{alternating_at} alternating AT/TA regions')
        
        return '; '.join(parts) if parts else 'Purine-pyrimidine alternating sequences'

    def _find_egz_motifs(self, seq: str) -> List[Dict[str, Any]]:
        """Find eGZ-motif (Extruded-G Z-DNA) using regex."""
        patterns = self.get_patterns().get('egz_motifs', []); annotations = []
        for pattern_info in patterns:
            regex, pattern_id, name, subclass, min_len, score_type, threshold, desc, ref = pattern_info
            for match in re.finditer(regex, seq, re.IGNORECASE):
                start, end = match.span(); matched_seq = seq[start:end]
                repeat_unit = matched_seq[:3].upper(); repeat_count = len(matched_seq) // 3
                repeat_score = threshold * (repeat_count / float(self.MIN_EGZ_REPEATS))
                annotations.append({
                    "start": start, "end": end, "length": end - start, "sum_score": round(repeat_score, 6),
                    "subclass": "eGZ", "pattern_id": pattern_id, "repeat_unit": repeat_unit,
                    "repeat_count": repeat_count, "description": desc, "reference": ref
                })
        return annotations

    def _find_10mer_matches(self, seq: str) -> List[Tuple[int, str, float]]:
        """Find all exact 10-mer matches."""
        if hyperscan_backend.is_hyperscan_available():
            try: return hyperscan_backend.hs_find_matches(seq, TENMER_SCORE)
            except Exception as e:
                logger.warning(f"Hyperscan failed for Z-DNA, falling back: {e}")
                return hyperscan_backend.py_find_matches(seq, TENMER_SCORE)
        return hyperscan_backend.py_find_matches(seq, TENMER_SCORE)

    def _merge_matches(self, matches: List[Tuple[int, str, float]], merge_gap: int = 0) -> List[Tuple[int, int, List[Tuple[int, str, float]]]]:
        """Merge overlapping/adjacent 10-mer matches."""
        if not matches: return []
        merged = []; cur_start = matches[0][0]; cur_end = matches[0][0] + 10; cur_matches = [matches[0]]
        for m in matches[1:]:
            s = m[0]; m_end = s + 10
            if s <= cur_end + merge_gap: cur_end = max(cur_end, m_end); cur_matches.append(m)
            else: merged.append((cur_start, cur_end, cur_matches)); cur_start, cur_end = s, m_end; cur_matches = [m]
        merged.append((cur_start, cur_end, cur_matches))
        return merged

    def _find_and_merge_10mer_matches(self, seq: str, merge_gap: int = 0) -> List[Tuple[int, int]]:
        matches = self._find_10mer_matches(seq); merged = self._merge_matches(matches, merge_gap=merge_gap)
        return [(s, e) for (s, e, _) in merged]

    def _build_per_base_contrib(self, seq: str, matches: List[Tuple] = None) -> List[float]:
        """Build per-base contribution array from 10-mer matches."""
        n = len(seq)
        if matches is None:
            matches = self._find_10mer_matches(seq)
        
        if _NUMPY_AVAILABLE and n > 1000 and len(matches) > 0:
            try:
                contrib = np.zeros(n, dtype=np.float64)
                per_base_val = 1.0 / 10.0
                
                for (start, ten, score) in matches:
                    per_base = float(score) * per_base_val
                    end = min(start + 10, n)
                    contrib[start:end] += per_base
                
                return contrib.tolist()
            except Exception:
                pass
        
        contrib = [0.0] * n
        for (start, ten, score) in matches:
            per_base = float(score) / 10.0
            for k in range(start, min(start + 10, n)):
                contrib[k] += per_base
        return contrib
