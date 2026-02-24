"""i-Motif DNA detector: canonical C-rich structures and HUR AC-motifs."""
# IMPORTS
import bisect
import re
from typing import Dict, List, Tuple, Any
from ..base.base_detector import BaseMotifDetector
from Utilities.detectors_utils import revcomp
from Utilities.core.motif_normalizer import normalize_class_subclass

try: from motif_patterns import IMOTIF_PATTERNS
except ImportError: IMOTIF_PATTERNS = {}

# TUNABLE PARAMETERS
MIN_REGION_LEN = 10; CLASS_PRIORITIES = {'canonical_imotif': 1, 'hur_ac_motif': 2}
VALIDATED_SEQS = [("IM_VAL_001", "CCCCTCCCCTCCCCTCCCC", "Validated i-motif 1", "Gehring 1993"),
                  ("IM_VAL_002", "CCCCACCCCACCCCACCCC", "Validated i-motif 2", "Leroy 1995")]

def _class_prio_idx(class_name: str) -> int: return CLASS_PRIORITIES.get(class_name, 999)

class IMotifDetector(BaseMotifDetector):
    """i-Motif DNA detector: canonical C-rich structures and HUR AC-motifs."""

    SCORE_REFERENCE = 'Gehring et al. 1993, Zeraati et al. 2018'

    def get_motif_class_name(self) -> str: return "i-Motif"

    def get_length_cap(self, subclass: str = None) -> int:
        """i-Motif structures stable up to ~60 bp (Gehring 1993, Zeraati 2018)."""
        return 60

    def theoretical_min_score(self) -> float:
        """Minimum biologically valid i-motif raw score."""
        return 0.4

    def theoretical_max_score(self, sequence_length: int = None) -> float:
        """Highest possible i-motif raw score.

        C-rich density + tract bonus, capped at 1.0.
        """
        return 1.0

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return i-motif patterns: canonical 4×C-tracts and HUR AC-motifs."""
        return self._load_patterns(IMOTIF_PATTERNS, lambda: {
            'canonical_imotif': [(r'C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}', 'IM_0', 'Canonical i-motif', 'canonical_imotif', 15, 'imotif_score', 0.95, 'pH-dependent C-rich structure', 'Gehring 1993')],
            'hur_ac_motif': [
                (r'A{3}[ACGT]{4}C{3}[ACGT]{4}C{3}[ACGT]{4}C{3}', 'HUR_AC_1', 'HUR AC-motif (4bp)', 'AC-motif (HUR)', 18, 'ac_motif_score', 0.85, 'HUR AC alternating motif', 'Hur 2021'),
                (r'C{3}[ACGT]{4}C{3}[ACGT]{4}C{3}[ACGT]{4}A{3}', 'HUR_AC_2', 'HUR CA-motif (4bp)', 'AC-motif (HUR)', 18, 'ac_motif_score', 0.85, 'HUR CA alternating motif', 'Hur 2021'),
                (r'A{3}[ACGT]{5}C{3}[ACGT]{5}C{3}[ACGT]{5}C{3}', 'HUR_AC_3', 'HUR AC-motif (5bp)', 'AC-motif (HUR)', 21, 'ac_motif_score', 0.85, 'HUR AC alternating motif', 'Hur 2021'),
                (r'C{3}[ACGT]{5}C{3}[ACGT]{5}C{3}[ACGT]{5}A{3}', 'HUR_AC_4', 'HUR CA-motif (5bp)', 'AC-motif (HUR)', 21, 'ac_motif_score', 0.85, 'HUR CA alternating motif', 'Hur 2021'),
                (r'A{3}[ACGT]{6}C{3}[ACGT]{6}C{3}[ACGT]{6}C{3}', 'HUR_AC_5', 'HUR AC-motif (6bp)', 'AC-motif (HUR)', 24, 'ac_motif_score', 0.85, 'HUR AC alternating motif', 'Hur 2021'),
                (r'C{3}[ACGT]{6}C{3}[ACGT]{6}C{3}[ACGT]{6}A{3}', 'HUR_AC_6', 'HUR CA-motif (6bp)', 'AC-motif (HUR)', 24, 'ac_motif_score', 0.85, 'HUR CA alternating motif', 'Hur 2021'),
            ]
        })

    def find_validated_matches(self, sequence: str, check_revcomp: bool = False) -> List[Dict[str, Any]]:
        seq = sequence.upper(); out = []
        for vid, vseq, desc, cite in VALIDATED_SEQS:
            idx = seq.find(vseq)
            if idx >= 0: out.append({'id': vid, 'seq': vseq, 'start': idx, 'end': idx+len(vseq), 'strand': '+', 'desc': desc, 'cite': cite})
            elif check_revcomp:
                rc = revcomp(vseq); idx2 = seq.find(rc)
                if idx2 >= 0: out.append({'id': vid, 'seq': vseq, 'start': idx2, 'end': idx2+len(vseq), 'strand': '-', 'desc': desc, 'cite': cite})
        return out

    def find_hur_ac_candidates(self, sequence: str, scan_rc: bool = True) -> List[Dict[str, Any]]:
        seq = sequence.upper(); candidates = []
        def _matches_hur_ac(target, strand):
            for nlink in (4, 5, 6):
                pat1 = r"A{3}[ACGT]{%d}C{3}[ACGT]{%d}C{3}[ACGT]{%d}C{3}" % (nlink, nlink, nlink)
                pat2 = r"C{3}[ACGT]{%d}C{3}[ACGT]{%d}C{3}[ACGT]{%d}A{3}" % (nlink, nlink, nlink)
                for pat in (pat1, pat2):
                    for m in re.finditer(pat, target):
                        s, e = m.span(); matched = m.group(0).upper()
                        candidates.append({'start': s if strand == '+' else len(seq) - e, 'end': e if strand == '+' else len(seq) - s,
                                           'strand': strand, 'linker': nlink, 'pattern': pat, 'matched_seq': matched,
                                           'loose_mode': True, 'high_confidence': (nlink == 4 or nlink == 5)})
        _matches_hur_ac(seq, '+')
        if scan_rc: _matches_hur_ac(revcomp(seq), '-')
        candidates.sort(key=lambda x: x['start']); return candidates

    def _find_regex_candidates(self, sequence: str) -> List[Dict[str, Any]]:
        seq = sequence.upper(); patterns = self.get_patterns(); out = []
        for class_name, pats in patterns.items():
            for patt in pats:
                regex = patt[0]; pid = patt[1] if len(patt) > 1 else f"{class_name}_pat"
                for m in re.finditer(regex, seq, flags=re.IGNORECASE | re.ASCII):
                    s, e = m.start(), m.end()
                    if (e - s) >= MIN_REGION_LEN: out.append({'class_name': class_name, 'pattern_id': pid, 'start': s, 'end': e, 'matched_seq': seq[s:e]})
        return out

    def _score_imotif_candidate(self, matched_seq: str) -> float:
        region = matched_seq.upper(); L = len(region)
        if L < 12: return 0.0
        c_tracts = [m.group(0) for m in re.finditer(r"C{2,}", region)]
        if len(c_tracts) < 3: return 0.0
        total_c = sum(len(t) for t in c_tracts); c_density = total_c / L
        tract_bonus = min(0.4, 0.12 * (len(c_tracts) - 2))
        return float(max(0.0, min(1.0, c_density + tract_bonus)))

    def _score_hur_ac_candidate(self, matched_seq: str, linker: int, high_confidence: bool) -> float:
        r = matched_seq.upper(); L = len(r); ac_count = r.count('A') + r.count('C')
        ac_frac = ac_count / L if L > 0 else 0.0
        a_tracts = [len(m.group(0)) for m in re.finditer(r"A{2,}", r)]; c_tracts = [len(m.group(0)) for m in re.finditer(r"C{2,}", r)]
        tract_score = 0.5 if (any(x >= 3 for x in a_tracts) and sum(1 for x in c_tracts if x >= 3) >= 3) else 0.0
        base = min(0.6, ac_frac * 0.8); linker_boost = 0.25 if high_confidence else (0.12 if linker == 6 else 0.0)
        return max(0.0, min(1.0, base + tract_score + linker_boost))

    def _resolve_overlaps_greedy(self, scored: List[Dict[str, Any]], merge_gap: int = 0) -> List[Dict[str, Any]]:
        if not scored: return []
        scored_sorted = sorted(scored, key=lambda x: (-x['score'], _class_prio_idx(x.get('class_name','')), -(x['end']-x['start'])))
        accepted: List[Dict[str, Any]] = []
        # Maintain sorted accepted-interval lists for O(log n) overlap checking.
        acc_starts: List[int] = []
        acc_ends: List[int] = []
        for cand in scored_sorted:
            s, e = cand['start'], cand['end']
            idx = bisect.bisect_left(acc_starts, e + merge_gap)
            conflict = (idx > 0 and acc_ends[idx - 1] + merge_gap > s) or \
                       (idx < len(acc_starts) and acc_starts[idx] - merge_gap < e)
            if not conflict:
                accepted.append(cand)
                ins = bisect.bisect_left(acc_starts, s)
                acc_starts.insert(ins, s)
                acc_ends.insert(ins, e)
        accepted.sort(key=lambda x: x['start']); return accepted

    def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:
        seq = sequence.upper()
        if self.find_validated_matches(seq, check_revcomp=False): return 0.99
        hur_cands = self.find_hur_ac_candidates(seq, scan_rc=True)
        hur_scored = [dict(class_name='ac_motif_hur', pattern_id=h['pattern'], start=h['start'], end=h['end'],
                           matched_seq=h['matched_seq'], linker=h['linker'], high_confidence=h['high_confidence'],
                           score=self._score_hur_ac_candidate(h['matched_seq'], h['linker'], h['high_confidence']), details=h) for h in hur_cands]
        regex_cands = self._find_regex_candidates(seq)
        regex_scored = [dict(class_name=r['class_name'], pattern_id=r['pattern_id'], start=r['start'], end=r['end'],
                             matched_seq=r['matched_seq'], score=self._score_imotif_candidate(r['matched_seq']), details={}) for r in regex_cands]
        combined = hur_scored + regex_scored; accepted = self._resolve_overlaps_greedy(combined, merge_gap=0)
        return float(sum(a['score'] * max(1, (a['end']-a['start'])/10.0) for a in accepted))

    def annotate_sequence(self, sequence: str) -> Dict[str, Any]:
        seq = sequence.upper(); res = {}
        res['validated_matches'] = self.find_validated_matches(seq, check_revcomp=True)
        hur_cands = self.find_hur_ac_candidates(seq, scan_rc=True)
        for h in hur_cands: h['score'] = self._score_hur_ac_candidate(h['matched_seq'], h['linker'], h['high_confidence'])
        res['hur_candidates'] = hur_cands
        regex_cands = self._find_regex_candidates(seq)
        for r in regex_cands: r['score'] = self._score_imotif_candidate(r['matched_seq'])
        res['regex_matches'] = regex_cands
        combined = [dict(class_name='ac_motif_hur', start=h['start'], end=h['end'], score=h['score'], details=h) for h in hur_cands]
        combined += [dict(class_name=r['class_name'], start=r['start'], end=r['end'], score=r['score'], details=r) for r in regex_cands]
        res['accepted'] = self._resolve_overlaps_greedy(combined, merge_gap=0)
        return res
    
    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Detect i-motif structures with component details and overlap resolution."""
        seq = sequence.upper(); motifs = []
        annotation_result = self.annotate_sequence(seq); accepted_motifs = annotation_result.get('accepted', [])
        subclass_map = {'canonical_imotif': 'Canonical i-motif', 'hur_ac_motif': 'AC-motif', 'ac_motif_hur': 'AC-motif'}

        for i, accepted in enumerate(accepted_motifs):
            start_pos, end_pos = accepted['start'], accepted['end']; motif_seq = seq[start_pos:end_pos]
            class_name = accepted.get('class_name', 'canonical_imotif'); subclass = subclass_map.get(class_name, 'Canonical i-motif')
            score = accepted.get('score', 0)
            canonical_class, canonical_subclass = normalize_class_subclass(self.get_motif_class_name(), subclass, strict=False, auto_correct=True)
            c_tracts = re.findall(r'C{2,}', motif_seq); loops = []
            c_tract_matches = list(re.finditer(r'C{2,}', motif_seq))
            for j in range(len(c_tract_matches) - 1):
                loop_start, loop_end = c_tract_matches[j].end(), c_tract_matches[j + 1].start()
                if loop_end > loop_start: loops.append(motif_seq[loop_start:loop_end])
            gc_total = self._calc_gc(motif_seq); gc_stems = self._calc_gc(''.join(c_tracts)) if c_tracts else 0
            motif = {'ID': f"{sequence_name}_IMOT_{start_pos+1}", 'Sequence_Name': sequence_name, 'Class': canonical_class,
                     'Subclass': canonical_subclass, 'Start': start_pos + 1, 'End': end_pos, 'Length': end_pos - start_pos,
                     'Sequence': motif_seq, 'Raw_Score': round(score, 3), 'Score': self.normalize_score(score, end_pos - start_pos, canonical_subclass), 'Strand': '+', 'Method': 'i-Motif_detection',
                     'Pattern_ID': f'IMOT_{i+1}', 'Stems': c_tracts, 'Loops': loops, 'Num_Stems': len(c_tracts),
                     'Num_Loops': len(loops), 'Stem_Lengths': [len(s) for s in c_tracts], 'Loop_Lengths': [len(l) for l in loops],
                     'GC_Content': round(gc_total, 2), 'GC_Total': round(gc_total, 2), 'GC_Stems': round(gc_stems, 2),
                     'Type_Of_Repeat': self._classify_imotif_type(canonical_subclass, len(c_tracts), motif_seq),
                     'Criterion': self._get_imotif_criterion(canonical_subclass, len(c_tracts), loops),
                     'Disease_Relevance': self._get_imotif_disease_relevance(motif_seq, gc_total, end_pos - start_pos),
                     'Regions_Involved': self._describe_imotif_regions(c_tracts, loops)}
            if c_tracts: 
                motif['Stem_Length'] = sum(len(s) for s in c_tracts) / len(c_tracts)
                motif['Arm_Length'] = motif['Stem_Length']  # Arm_Length maps to Stem_Length for i-Motif
            if loops: motif['Loop_Length'] = sum(len(l) for l in loops) / len(loops)
            motifs.append(motif)
        return motifs
    
    def _classify_imotif_type(self, subclass: str, num_stems: int, sequence: str) -> str:
        """Classify i-motif structural type"""
        if 'AC-motif' in subclass:
            return 'AC-motif (HUR) - A-tract/C-tract alternating'
        elif num_stems >= 4:
            return 'Four-stranded canonical i-motif (C-rich)'
        elif num_stems == 3:
            return 'Three-stranded i-motif (C-rich)'
        else:
            return 'C-rich i-motif-like structure'
    
    def _get_imotif_criterion(self, subclass: str, num_stems: int, loops: List[str]) -> str:
        """Explain i-motif classification criterion"""
        if 'AC-motif' in subclass:
            return f'HUR AC-motif: A3-linker-C3 pattern (Hur 2021); {num_stems} C-tracts detected'
        else:
            loop_range = f'loops 1-7bp' if loops and all(1 <= len(l) <= 7 for l in loops) else 'variable loops'
            return f'Canonical i-motif: {num_stems}+ C-tracts (≥3C each), {loop_range}; Gehring 1993'
    
    def _get_imotif_disease_relevance(self, sequence: str, gc_content: float, length: int) -> str:
        """Annotate disease relevance for i-motif"""
        disease_notes = []
        
        # i-Motifs in promoters - oncogene regulation
        if gc_content > 65:
            disease_notes.append('Promoter-like i-motif (high C-content) - oncogene regulation: BCL2, VEGF, c-MYC')
        
        # i-Motifs and telomeres
        if 'CCCTAA' in sequence or 'CCCTTA' in sequence:
            disease_notes.append('Telomeric i-motif (C-rich telomeric strand) - chromosome stability, aging')
        
        # Long i-motifs
        if length > 50:
            disease_notes.append('Extended i-motif - pH-dependent gene regulation, genome stability')
        
        # General associations
        disease_notes.append('i-Motif formation - pH sensor, transcription regulation, potential therapeutic target')
        
        return '; '.join(disease_notes)
    
    def _describe_imotif_regions(self, c_tracts: List[str], loops: List[str]) -> str:
        """Describe regions involved in i-motif formation"""
        parts = []
        
        if c_tracts:
            parts.append(f'{len(c_tracts)} C-tracts: {", ".join([f"C{len(t)}" for t in c_tracts])}')
        
        if loops:
            parts.append(f'{len(loops)} loops: {", ".join([f"{len(l)}bp" for l in loops])}')
        
        return '; '.join(parts) if parts else 'C-rich i-motif structure'
