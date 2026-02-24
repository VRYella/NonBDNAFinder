"""Curved DNA detector using A-tract and T-tract phasing patterns."""
# IMPORTS
import re
from typing import List, Dict, Any, Tuple
from ..base.base_detector import BaseMotifDetector
from Utilities.detectors_utils import revcomp
from .patterns import _generate_phased_repeat_patterns
from Utilities.core.motif_normalizer import normalize_class_subclass

# TUNABLE PARAMETERS
MIN_AT_TRACT = 3; MAX_AT_WINDOW = None; PHASING_CENTER_SPACING = 11.0
PHASING_TOL_LOW = 9.9; PHASING_TOL_HIGH = 11.1; MIN_APR_TRACTS = 3; LOCAL_LONG_TRACT = 8; SCORE_THRESHOLD = 0.1

# NORMALIZATION PARAMETERS (Tunable)
# ┌──────────────┬─────────────┬────────────────────────────────────────┐
# │ Parameter    │ Value       │ Scientific Basis                       │
# ├──────────────┼─────────────┼────────────────────────────────────────┤
# │ RAW_MIN      │ 0.1         │ Koo 1986 - minimal curvature threshold │
# │ RAW_MAX      │ 0.95        │ Koo 1986 - strong curvature (APR)      │
# │ NORM_MIN     │ 1.0         │ Universal low confidence threshold     │
# │ NORM_MAX     │ 3.0         │ Universal high confidence threshold    │
# │ METHOD       │ 'linear'    │ Linear interpolation                   │
# └──────────────┴─────────────┴────────────────────────────────────────┘
CURVED_RAW_SCORE_MIN = 0.1; CURVED_RAW_SCORE_MAX = 0.95
CURVED_NORMALIZED_MIN = 1.0; CURVED_NORMALIZED_MAX = 3.0
CURVED_NORMALIZATION_METHOD = 'linear'
CURVED_SCORE_REFERENCE = 'Koo et al. 1986, Olson et al. 1998'


class CurvedDNADetector(BaseMotifDetector):
    """Curved DNA detector using A-tract and T-tract phasing patterns."""

    # Override normalization parameters
    RAW_SCORE_MIN = CURVED_RAW_SCORE_MIN
    RAW_SCORE_MAX = CURVED_RAW_SCORE_MAX
    NORMALIZED_MIN = CURVED_NORMALIZED_MIN
    NORMALIZED_MAX = CURVED_NORMALIZED_MAX
    NORMALIZATION_METHOD = CURVED_NORMALIZATION_METHOD
    SCORE_REFERENCE = CURVED_SCORE_REFERENCE

    MIN_AT_TRACT = MIN_AT_TRACT; MAX_AT_WINDOW = MAX_AT_WINDOW; PHASING_CENTER_SPACING = PHASING_CENTER_SPACING
    PHASING_TOL_LOW = PHASING_TOL_LOW; PHASING_TOL_HIGH = PHASING_TOL_HIGH; MIN_APR_TRACTS = MIN_APR_TRACTS
    LOCAL_LONG_TRACT = LOCAL_LONG_TRACT; SCORE_THRESHOLD = SCORE_THRESHOLD

    def get_motif_class_name(self) -> str: return "Curved_DNA"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Generate curved DNA patterns programmatically."""
        patterns = {'local_curved': [
            (r'A{7,}', 'CRV_002', 'Long A-tract', 'Local Curvature', 7, 'curvature_score', 0.95, 'A-tract curvature', 'Olson 1998'),
            (r'T{7,}', 'CRV_003', 'Long T-tract', 'Local Curvature', 7, 'curvature_score', 0.95, 'T-tract curvature', 'Olson 1998'),
        ]}
        tract_range = range(3, 10)
        patterns['global_curved_a_3tract'] = _generate_phased_repeat_patterns('A', 3, tract_range, 8)
        patterns['global_curved_a_4tract'] = _generate_phased_repeat_patterns('A', 4, tract_range, 15)
        patterns['global_curved_a_5tract'] = _generate_phased_repeat_patterns('A', 5, tract_range, 22)
        patterns['global_curved_t_3tract'] = _generate_phased_repeat_patterns('T', 3, tract_range, 29)
        patterns['global_curved_t_4tract'] = _generate_phased_repeat_patterns('T', 4, tract_range, 36)
        patterns['global_curved_t_5tract'] = _generate_phased_repeat_patterns('T', 5, tract_range, 43)
        return patterns

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Detect curved DNA motifs using A/T-tract phasing patterns."""
        sequence = sequence.upper().strip(); motifs = []
        annotation = self.annotate_sequence(sequence)

        for i, apr in enumerate(annotation.get('aprs', [])):
            if apr.get('score', 0) > self.SCORE_THRESHOLD:
                center_positions = apr.get('center_positions', [])
                if not center_positions:
                    continue
                start_pos = max(0, int(min(center_positions)) - 10)
                end_pos = min(len(sequence), int(max(center_positions)) + 10)
                motif_seq = sequence[start_pos:end_pos]
                a_tracts = re.findall(r'A{3,}', motif_seq); t_tracts = re.findall(r'T{3,}', motif_seq)
                gc_total = self._calc_gc(motif_seq); at_content = self._calc_at(motif_seq)
                canonical_class, canonical_subclass = normalize_class_subclass(self.get_motif_class_name(), 'Global Curvature', strict=False, auto_correct=True)
                
                spacing_info = self._get_spacing_info(apr.get('center_positions', []))
                
                motifs.append({
                    'ID': f"{sequence_name}_CRV_APR_{start_pos+1}", 'Sequence_Name': sequence_name, 'Class': canonical_class,
                    'Subclass': canonical_subclass, 'Start': start_pos + 1, 'End': end_pos, 'Length': end_pos - start_pos,
                    'Sequence': motif_seq, 'Raw_Score': round(apr.get('score', 0), 3), 'Score': self._normalize_score(apr.get('score', 0)), 'Strand': '+', 'Method': 'Curved_DNA_detection',
                    'Pattern_ID': f'CRV_APR_{i+1}', 'A_Tracts': a_tracts, 'T_Tracts': t_tracts, 'Num_A_Tracts': len(a_tracts),
                    'Num_T_Tracts': len(t_tracts), 'A_Tract_Lengths': [len(t) for t in a_tracts],
                    'T_Tract_Lengths': [len(t) for t in t_tracts], 'GC_Content': round(gc_total, 2),
                    'AT_Content': round(at_content, 2), 'Center_Positions': apr.get('center_positions', []),
                    'Arm_Length': 'N/A',
                    'Loop_Length': 'N/A',
                    'Type_Of_Repeat': 'A/T-tract phased repeat (global APR)',
                    'Criterion': self._get_curved_criterion(apr, 'APR'),
                    'Disease_Relevance': self._get_curved_disease_relevance(at_content, end_pos - start_pos),
                    'Regions_Involved': f'{len(a_tracts)} A-tracts + {len(t_tracts)} T-tracts; {spacing_info}'
                })

        for i, tract in enumerate(annotation.get('long_tracts', [])):
            if tract.get('score', 0) > self.SCORE_THRESHOLD:
                start_pos, end_pos = tract['start'], tract['end']; motif_seq = sequence[start_pos:end_pos]
                tract_type = 'A-tract' if motif_seq.count('A') > motif_seq.count('T') else 'T-tract'
                gc_total = self._calc_gc(motif_seq); at_content = self._calc_at(motif_seq)
                canonical_class, canonical_subclass = normalize_class_subclass(self.get_motif_class_name(), 'Local Curvature', strict=False, auto_correct=True)
                motifs.append({
                    'ID': f"{sequence_name}_CRV_TRACT_{start_pos+1}", 'Sequence_Name': sequence_name, 'Class': canonical_class,
                    'Subclass': canonical_subclass, 'Start': start_pos + 1, 'End': end_pos, 'Length': end_pos - start_pos,
                    'Sequence': motif_seq, 'Raw_Score': round(tract.get('score', 0), 3), 'Score': self._normalize_score(tract.get('score', 0)), 'Strand': '+', 'Method': 'Curved_DNA_detection',
                    'Pattern_ID': f'CRV_TRACT_{i+1}', 'Tract_Type': tract_type, 'Tract_Length': end_pos - start_pos,
                    'GC_Content': round(gc_total, 2), 'AT_Content': round(at_content, 2),
                    'Arm_Length': 'N/A',
                    'Loop_Length': 'N/A',
                    'Type_Of_Repeat': f'{tract_type} homopolymer',
                    'Criterion': self._get_curved_criterion(tract, 'TRACT'),
                    'Disease_Relevance': self._get_curved_disease_relevance(at_content, end_pos - start_pos),
                    'Regions_Involved': f'Single {tract_type} ({end_pos - start_pos}bp)'
                })
        return motifs
    
    def _get_spacing_info(self, center_positions: List[float]) -> str:
        """Get spacing information for A/T-tracts"""
        if len(center_positions) < 2:
            return 'single tract'
        
        spacings = []
        for i in range(1, len(center_positions)):
            spacing = center_positions[i] - center_positions[i-1]
            spacings.append(spacing)
        
        if spacings:
            avg_spacing = sum(spacings) / len(spacings)
            return f'avg spacing {avg_spacing:.1f}bp (helical phasing ~10.5bp)'
        return 'multiple tracts'
    
    def _get_curved_criterion(self, element: Dict[str, Any], element_type: str) -> str:
        """Explain curved DNA classification criterion"""
        if element_type == 'APR':
            return f'A-phased repeat: ≥3 A/T-tracts with phasing score {element.get("score", 0):.3f} ≥{self.SCORE_THRESHOLD}; helical periodicity ~10.5bp'
        else:
            return f'Long A/T-tract: tract length ≥{self.MIN_AT_TRACT}bp, score {element.get("score", 0):.3f} ≥{self.SCORE_THRESHOLD}'
    
    def _get_curved_disease_relevance(self, at_content: float, length: int) -> str:
        """Annotate disease relevance for curved DNA"""
        disease_notes = []
        
        if at_content > 80:
            disease_notes.append('High AT-content curvature - nucleosome positioning, chromatin structure')
        
        if length > 100:
            disease_notes.append('Extended curvature - DNA packaging, gene regulation')
        
        disease_notes.append('Curved DNA - transcription factor binding, replication origin, chromatin organization')
        
        return '; '.join(disease_notes)

    def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:
        seq = sequence.upper(); ann = self.annotate_sequence(seq)
        apr_sum = sum(a['score'] for a in ann.get('aprs', []))
        local_sum = sum(l['score'] for l in ann.get('long_tracts', []))
        return float(apr_sum + local_sum)

    def find_a_tracts(self, sequence: str, minAT: int = None, max_window: int = None) -> List[Dict[str, Any]]:
        """Detect A-tract candidates using AT-window analysis. Returns dicts with start/end, maxATlen, maxTlen, a_center, call, window info."""
        seq = sequence.upper()
        if minAT is None:
            minAT = self.MIN_AT_TRACT
        if max_window is None:
            max_window = self.MAX_AT_WINDOW  # None allowed

        results: List[Dict[str, Any]] = []

        for m in re.finditer(r'[AT]{' + str(minAT) + r',}', seq):
            wstart, wend = m.start(), m.end()  # [wstart, wend)
            window_seq = seq[wstart:wend]
            window_len = wend - wstart

            maxATlen, maxATend, maxTlen = self._analyze_at_window(window_seq)
            rc_window = revcomp(window_seq)
            maxATlen_rc, maxATend_rc, maxTlen_rc = self._analyze_at_window(rc_window)

            diff_forward = maxATlen - maxTlen
            diff_rc = maxATlen_rc - maxTlen_rc
            call = False
            chosen_center = None
            chosen_maxATlen = None

            if diff_forward >= minAT or diff_rc >= minAT:
                call = True
                if diff_forward >= diff_rc:
                    chosen_maxATlen = maxATlen
                    chosen_center = (wstart + maxATend) - ((maxATlen - 1) / 2.0)
                else:
                    chosen_maxATlen = maxATlen_rc
                    rc_end_original = wstart + (window_len - 1 - maxATend_rc)
                    chosen_center = rc_end_original - ((chosen_maxATlen - 1) / 2.0)

            results.append({
                'start': wstart,
                'end': wend,
                'window_len': window_len,
                'window_seq': window_seq,
                'maxATlen': int(maxATlen),
                'maxATend': int(wstart + maxATend),
                'maxTlen': int(maxTlen),
                'maxATlen_rc': int(maxATlen_rc),
                'maxATend_rc': int(wstart + (window_len - 1 - maxATend_rc)),
                'maxTlen_rc': int(maxTlen_rc),
                'diff_forward': int(diff_forward),
                'diff_rc': int(diff_rc),
                'call': bool(call),
                'a_center': float(chosen_center) if chosen_center is not None else None,
                'chosen_maxATlen': int(chosen_maxATlen) if chosen_maxATlen is not None else None
            })

        return results

    def _analyze_at_window(self, window_seq: str) -> Tuple[int,int,int]:
        """Analyze A/T window returning (maxATlen, maxATend_index_in_window, maxTlen). Tracks Alen, Tlen, ATlen, TAlen using a sliding algorithm."""
        Alen = 0
        Tlen = 0
        ATlen = 0
        TAlen = 0
        maxATlen = 0
        maxTlen = 0
        maxATend = 0
        # we'll iterate from index 0..len(window_seq)-1
        L = len(window_seq)
        for i in range(L):
            ch = window_seq[i]
            prev = window_seq[i-1] if i>0 else None
            if ch == 'A':
                Tlen = 0
                TAlen = 0
                # if previous base was T, reset A-run counters
                if prev == 'T':
                    Alen = 1
                    ATlen = 1
                else:
                    Alen += 1
                    ATlen += 1
            elif ch == 'T':
                # if T follows A-run shorter than Alen, it's considered TAlen (T following A)
                if TAlen < Alen:
                    TAlen += 1
                    ATlen += 1
                else:
                    # T is starting a T-only run
                    Tlen += 1
                    TAlen = 0
                    ATlen = 0
                    Alen = 0
            else:
                # non-AT not expected inside this window (we only pass contiguous AT windows)
                Alen = 0
                Tlen = 0
                ATlen = 0
                TAlen = 0
            if ATlen > maxATlen:
                maxATlen = ATlen
                maxATend = i  # end index within window
            if Tlen > maxTlen:
                maxTlen = Tlen
        return int(maxATlen), int(maxATend), int(maxTlen)

    # -------------------------
    # APR grouping / phasing
    # -------------------------
    def find_aprs(self, sequence: str, min_tract: int = None, min_apr_tracts: int = None) -> List[Dict[str, Any]]:
        """Group a-tract centers into APRs (A-phased repeats). Requires min_apr_tracts centers with spacing in PHASING_TOL_LOW..PHASING_TOL_HIGH range."""
        if min_tract is None:
            min_tract = self.MIN_AT_TRACT
        if min_apr_tracts is None:
            min_apr_tracts = self.MIN_APR_TRACTS

        a_calls = [r for r in self.find_a_tracts(sequence, minAT=min_tract) if r['call'] and r['a_center'] is not None]
        centers = [r['a_center'] for r in a_calls]
        centers_sorted = sorted(centers)
        aprs: List[Dict[str, Any]] = []

        if len(centers_sorted) < min_apr_tracts:
            return aprs

        i = 0
        while i < len(centers_sorted):
            run = [centers_sorted[i]]
            j = i + 1
            while j < len(centers_sorted):
                spacing = centers_sorted[j] - centers_sorted[j-1]
                if self.PHASING_TOL_LOW <= spacing <= self.PHASING_TOL_HIGH:
                    run.append(centers_sorted[j])
                    j += 1
                else:
                    break
            if len(run) >= min_apr_tracts:
                spacings = [run[k+1] - run[k] for k in range(len(run)-1)]
                devs = [abs(sp - self.PHASING_CENTER_SPACING) for sp in spacings]
                mean_dev = sum(devs) / len(devs) if devs else 0.0
                max_dev_allowed = max(abs(self.PHASING_TOL_HIGH - self.PHASING_CENTER_SPACING),
                                      abs(self.PHASING_CENTER_SPACING - self.PHASING_TOL_LOW))
                phasing_score = max(0.0, 1.0 - (mean_dev / (max_dev_allowed if max_dev_allowed>0 else 1.0)))
                aprs.append({
                    'start_center_idx': i,
                    'end_center_idx': j-1,
                    'center_positions': run,
                    'n_tracts': len(run),
                    'spacings': spacings,
                    'mean_deviation': mean_dev,
                    'score': round(phasing_score, 6)
                })
            i = j

        return aprs

    # -------------------------
    # Local long tract finder
    # -------------------------
    def find_long_tracts(self, sequence: str, min_len: int = None) -> List[Dict[str, Any]]:
        """Find long A-tracts or T-tracts >= min_len (default LOCAL_LONG_TRACT). Returns {start,end,base,len,score} dicts."""
        if min_len is None:
            min_len = self.LOCAL_LONG_TRACT
        seq = sequence.upper()
        results = []
        for m in re.finditer(r'A{' + str(min_len) + r',}', seq):
            ln = m.end() - m.start()
            # simple local score: normalized by (len/(len+6)) to saturate
            score = float(ln) / (ln + 6.0)
            results.append({'start': m.start(), 'end': m.end(), 'base': 'A', 'len': ln, 'score': round(score, 6), 'seq': seq[m.start():m.end()]})
        for m in re.finditer(r'T{' + str(min_len) + r',}', seq):
            ln = m.end() - m.start()
            score = float(ln) / (ln + 6.0)
            results.append({'start': m.start(), 'end': m.end(), 'base': 'T', 'len': ln, 'score': round(score, 6), 'seq': seq[m.start():m.end()]})
        results.sort(key=lambda x: x['start'])
        return results

    # -------------------------
    # Scoring helpers (interpretability)
    # -------------------------
    def phasing_score(self, apr: Dict[str, Any]) -> float:
        """Return APR phasing score (already stored in apr['score'])."""
        return float(apr.get('score', 0.0))

    def local_curvature_score(self, tract: Dict[str, Any]) -> float:
        """Return local curvature score for a long tract (already stored)."""
        return float(tract.get('score', 0.0))

    # -------------------------
    # Annotate (summary)
    # -------------------------
    def annotate_sequence(self, sequence: str) -> Dict[str, Any]:
        """Returns comprehensive annotation: a_tract_windows, aprs with phasing scores, long_tracts, summary counts and combined score."""
        seq = sequence.upper()
        a_windows = self.find_a_tracts(seq, minAT=self.MIN_AT_TRACT)
        a_centers = [w for w in a_windows if w['call'] and w['a_center'] is not None]
        aprs = self.find_aprs(seq, min_tract=self.MIN_AT_TRACT, min_apr_tracts=self.MIN_APR_TRACTS)
        long_tracts = self.find_long_tracts(seq, min_len=self.LOCAL_LONG_TRACT)

        for apr in aprs:
            apr['constituent_windows'] = []
            for center in apr['center_positions']:
                # find closest a_window with that center
                best = min(a_windows, key=lambda w: abs((w['a_center'] or 1e9) - center))
                apr['constituent_windows'].append(best)

        summary = {
            'n_a_windows': len(a_windows),
            'n_a_centers': len(a_centers),
            'n_aprs': len(aprs),
            'n_long_tracts': len(long_tracts),
            'apr_score_sum': sum(self.phasing_score(a) for a in aprs),
            'long_tract_score_sum': sum(self.local_curvature_score(l) for l in long_tracts),
            'combined_score': sum(self.phasing_score(a) for a in aprs) + sum(self.local_curvature_score(l) for l in long_tracts)
        }

        return {
            'a_tract_windows': a_windows,
            'a_centers': a_centers,
            'aprs': aprs,
            'long_tracts': long_tracts,
            'summary': summary
        }
