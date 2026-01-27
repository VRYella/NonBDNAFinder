"""
Curved DNA Detector
Dr. Venkata Rajesh Yella | 2024.1 | MIT License

Detects curved DNA motifs using A-tract and T-tract phasing patterns.
"""

import re
from typing import List, Dict, Any, Tuple

from ..base.base_detector import BaseMotifDetector
from detectors_utils import revcomp
from .patterns import _generate_phased_repeat_patterns
from core.motif_normalizer import normalize_class_subclass


class CurvedDNADetector(BaseMotifDetector):
    """
    Detects curved DNA motifs using A-tract and T-tract phasing patterns.
    
    # Motif Structure:
    """
    
    def get_motif_class_name(self) -> str: return "Curved_DNA"

    # Tunable parameters
    MIN_AT_TRACT = 3; MAX_AT_WINDOW = None; PHASING_CENTER_SPACING = 11.0
    PHASING_TOL_LOW = 9.9; PHASING_TOL_HIGH = 11.1
    MIN_APR_TRACTS = 3; LOCAL_LONG_TRACT = 7

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Generate curved DNA patterns programmatically."""
        patterns = {
            'local_curved': [
                (r'A{7,}', 'CRV_002', 'Long A-tract', 'Local Curvature', 7, 
                 'curvature_score', 0.95, 'A-tract curvature', 'Olson 1998'),
                (r'T{7,}', 'CRV_003', 'Long T-tract', 'Local Curvature', 7, 
                 'curvature_score', 0.95, 'T-tract curvature', 'Olson 1998'),
            ]
        }
        
        # Generate A/T-phased repeats programmatically (3-, 4-, and 5-tract)
        tract_range = range(3, 10)
        patterns['global_curved_a_3tract'] = _generate_phased_repeat_patterns('A', 3, tract_range, 8)
        patterns['global_curved_a_4tract'] = _generate_phased_repeat_patterns('A', 4, tract_range, 15)
        patterns['global_curved_a_5tract'] = _generate_phased_repeat_patterns('A', 5, tract_range, 22)
        patterns['global_curved_t_3tract'] = _generate_phased_repeat_patterns('T', 3, tract_range, 29)
        patterns['global_curved_t_4tract'] = _generate_phased_repeat_patterns('T', 4, tract_range, 36)
        patterns['global_curved_t_5tract'] = _generate_phased_repeat_patterns('T', 5, tract_range, 43)
        
        return patterns

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Override base method to use sophisticated curved DNA detection with component details"""
        sequence = sequence.upper().strip()
        motifs = []
        
        # Use the sophisticated annotation method
        annotation = self.annotate_sequence(sequence)
        
        for i, apr in enumerate(annotation.get('aprs', [])):
            if apr.get('score', 0) > 0.1:  # Lower threshold for sensitivity
                start_pos = int(min(apr['center_positions'])) - 10  # Estimate start
                end_pos = int(max(apr['center_positions'])) + 10    # Estimate end
                start_pos = max(0, start_pos)
                end_pos = min(len(sequence), end_pos)
                
                motif_seq = sequence[start_pos:end_pos]
                
                a_tracts = re.findall(r'A{3,}', motif_seq)
                t_tracts = re.findall(r'T{3,}', motif_seq)
                
                gc_total = self._calc_gc(motif_seq); at_content = self._calc_at(motif_seq)
                
                # Normalize class/subclass using canonical taxonomy
                canonical_class, canonical_subclass = normalize_class_subclass(
                    self.get_motif_class_name(),
                    'Global Curvature',
                    strict=False,
                    auto_correct=True
                )
                
                motifs.append({
                    'ID': f"{sequence_name}_CRV_APR_{start_pos+1}",
                    'Sequence_Name': sequence_name,
                    'Class': canonical_class,
                    'Subclass': canonical_subclass,
                    'Start': start_pos + 1,  # 1-based coordinates
                    'End': end_pos,
                    'Length': end_pos - start_pos,
                    'Sequence': motif_seq,
                    'Score': round(apr.get('score', 0), 3),
                    'Strand': '+',
                    'Method': 'Curved_DNA_detection',
                    'Pattern_ID': f'CRV_APR_{i+1}',
                    # Component details
                    'A_Tracts': a_tracts,
                    'T_Tracts': t_tracts,
                    'Num_A_Tracts': len(a_tracts),
                    'Num_T_Tracts': len(t_tracts),
                    'A_Tract_Lengths': [len(t) for t in a_tracts],
                    'T_Tract_Lengths': [len(t) for t in t_tracts],
                    'GC_Content': round(gc_total, 2),
                    'AT_Content': round(at_content, 2),
                    'Center_Positions': apr.get('center_positions', [])
                })
        
        for i, tract in enumerate(annotation.get('long_tracts', [])):
            if tract.get('score', 0) > 0.1:  # Lower threshold for sensitivity
                start_pos = tract['start']
                end_pos = tract['end']
                motif_seq = sequence[start_pos:end_pos]
                
                # Identify tract type
                tract_type = 'A-tract' if motif_seq.count('A') > motif_seq.count('T') else 'T-tract'
                
                gc_total = self._calc_gc(motif_seq); at_content = self._calc_at(motif_seq)
                
                # Normalize class/subclass using canonical taxonomy
                canonical_class, canonical_subclass = normalize_class_subclass(
                    self.get_motif_class_name(),
                    'Local Curvature',
                    strict=False,
                    auto_correct=True
                )
                
                motifs.append({
                    'ID': f"{sequence_name}_CRV_TRACT_{start_pos+1}",
                    'Sequence_Name': sequence_name,
                    'Class': canonical_class,
                    'Subclass': canonical_subclass,
                    'Start': start_pos + 1,  # 1-based coordinates
                    'End': end_pos,
                    'Length': end_pos - start_pos,
                    'Sequence': motif_seq,
                    'Score': round(tract.get('score', 0), 3),
                    'Strand': '+',
                    'Method': 'Curved_DNA_detection',
                    'Pattern_ID': f'CRV_TRACT_{i+1}',
                    # Component details
                    'Tract_Type': tract_type,
                    'Tract_Length': end_pos - start_pos,
                    'GC_Content': round(gc_total, 2),
                    'AT_Content': round(at_content, 2)
                })
        
        # Remove overlaps within each subclass using base class method
        motifs = self._remove_overlaps(motifs, by_subclass=True)
        
        return motifs

    # -------------------------
    # Top-level scoring API
    # -------------------------
    def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:
        """
        Returns a combined raw score reflecting:
          - phasing_score for APRs (sum of APR phasing scores)
          - local curvature contribution (sum of local A/T tract scores)
        The sum reflects both number and quality of hits.
        """
        seq = sequence.upper()
        ann = self.annotate_sequence(seq)
        # Sum APR scores
        apr_sum = sum(a['score'] for a in ann.get('aprs', []))
        local_sum = sum(l['score'] for l in ann.get('long_tracts', []))
        return float(apr_sum + local_sum)

    # -------------------------
    # A-tract detection (core)
    # -------------------------
    def find_a_tracts(self, sequence: str, minAT: int = None, max_window: int = None) -> List[Dict[str, Any]]:
        """Detect A-tract candidates using AT-window analysis. Returns dicts with start/end, maxATlen, maxTlen, a_center, call, window info."""
        seq = sequence.upper()
        len(seq)
        if minAT is None:
            minAT = self.MIN_AT_TRACT
        if max_window is None:
            max_window = self.MAX_AT_WINDOW  # None allowed

        results: List[Dict[str, Any]] = []

        for m in re.finditer(r'[AT]{' + str(minAT) + r',}', seq):
            wstart, wend = m.start(), m.end()  # [wstart, wend)
            window_seq = seq[wstart:wend]
            window_len = wend - wstart

            # analyze forward strand window
            maxATlen, maxATend, maxTlen = self._analyze_at_window(window_seq)
            # analyze reverse complement window (to mimic C code check on reverse)
            rc_window = revcomp(window_seq)
            maxATlen_rc, maxATend_rc, maxTlen_rc = self._analyze_at_window(rc_window)

            diff_forward = maxATlen - maxTlen
            diff_rc = maxATlen_rc - maxTlen_rc
            call = False
            chosen_center = None
            chosen_maxATlen = None

            if diff_forward >= minAT or diff_rc >= minAT:
                call = True
                # choose the strand giving larger difference
                if diff_forward >= diff_rc:
                    chosen_maxATlen = maxATlen
                    # in C code: a_center = maxATend - ((maxATlen-1)/2) + 1  (1-based)
                    # we'll produce 0-based center = (wstart + maxATend - ((maxATlen-1)/2))
                    chosen_center = (wstart + maxATend) - ((maxATlen - 1) / 2.0)
                else:
                    chosen_maxATlen = maxATlen_rc
                    # maxATend_rc is position in RC sequence; convert to original coords:
                    # RC index i corresponds to original index: wstart + (window_len - 1 - i)
                    # maxATend_rc is index in window_rc (end position index)
                    # In C, they convert similarly; for simplicity compute center via rc mapping
                    # rc_end_original = wstart + (window_len - 1 - i_rc)
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
        """Analyze A/T window returning (maxATlen, maxATend_index_in_window, maxTlen). Tracks Alen, Tlen, ATlen, TAlen per C code logic."""
        Alen = 0
        Tlen = 0
        ATlen = 0
        TAlen = 0
        maxATlen = 0
        maxTlen = 0
        maxATend = 0
        # we'll iterate from index 0..len(window_seq)-1
        L = len(window_seq)
        # to mimic C code scanning with lookbacks, we iterate straightforwardly
        for i in range(L):
            ch = window_seq[i]
            prev = window_seq[i-1] if i>0 else None
            if ch == 'A':
                Tlen = 0
                TAlen = 0
                # if previous base was T, reset A-run counters per C code
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

        # find runs of centers where consecutive spacing is within tolerance
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
            # if run has enough tracts, call APR
            if len(run) >= min_apr_tracts:
                # score APR by how close spacings are to ideal spacing
                spacings = [run[k+1] - run[k] for k in range(len(run)-1)]
                # closeness = product of gaussian-like terms, but simpler: average deviation
                devs = [abs(sp - self.PHASING_CENTER_SPACING) for sp in spacings]
                # normalized closeness
                mean_dev = sum(devs) / len(devs) if devs else 0.0
                # phasing_score between 0..1: 1 when mean_dev==0, drop linearly with dev up to tolerance
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
        # A runs
        for m in re.finditer(r'A{' + str(min_len) + r',}', seq):
            ln = m.end() - m.start()
            # simple local score: normalized by (len/(len+6)) to saturate
            score = float(ln) / (ln + 6.0)
            results.append({'start': m.start(), 'end': m.end(), 'base': 'A', 'len': ln, 'score': round(score, 6), 'seq': seq[m.start():m.end()]})
        # T runs
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

        # annotate aprs with constituent windows (optional)
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
