"""Curved DNA detector using A-tract and T-tract phasing patterns."""

# IMPORTS
import re
from typing import List, Dict, Any, Tuple
from ..base.base_detector import BaseMotifDetector
from .patterns import _generate_phased_repeat_patterns
from Utilities.core.motif_normalizer import normalize_class_subclass

# TUNABLE PARAMETERS
MIN_AT_TRACT = 3; MAX_AT_WINDOW = None; PHASING_CENTER_SPACING = 11.0
PHASING_TOL_LOW = 9.9; PHASING_TOL_HIGH = 11.1; MIN_APR_TRACTS = 3
LOCAL_LONG_TRACT = 8; SCORE_THRESHOLD = 0.1


class CurvedDNADetector(BaseMotifDetector):

    SCORE_REFERENCE = 'Koo et al. 1986, Olson et al. 1998'

    MIN_AT_TRACT = MIN_AT_TRACT; MAX_AT_WINDOW = MAX_AT_WINDOW
    PHASING_CENTER_SPACING = PHASING_CENTER_SPACING
    PHASING_TOL_LOW = PHASING_TOL_LOW; PHASING_TOL_HIGH = PHASING_TOL_HIGH
    MIN_APR_TRACTS = MIN_APR_TRACTS
    LOCAL_LONG_TRACT = LOCAL_LONG_TRACT
    SCORE_THRESHOLD = SCORE_THRESHOLD

    def get_motif_class_name(self) -> str:
        return "Curved_DNA"

    def get_length_cap(self, subclass: str = None) -> int:
        """Local curvature (long A/T tract) stable up to ~50 bp;
        global curvature (APR phasing) stable up to ~120 bp (Koo 1986)."""
        if subclass == "Local Curvature":
            return 50
        return 120

    def theoretical_min_score(self) -> float:
        """Minimum biologically valid curved DNA raw score (score threshold)."""
        return self.SCORE_THRESHOLD

    def theoretical_max_score(self, sequence_length: int = None) -> float:
        """Highest possible curved DNA raw score.

        Score = ln / (ln + 7) → 1.0 as length → ∞.
        """
        return 1.0

    def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:
        ln = len(sequence)
        return float(ln) / (ln + 7.0)

    def passes_quality_threshold(self, sequence: str, score: float,
                                 pattern_info: Tuple = None) -> bool:
        return score >= self.SCORE_THRESHOLD

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        patterns = {'local_curved': [
            (r'A{8,}', 'CRV_002', 'Long A-tract', 'Local Curvature', 8,
             'curvature_score', 0.95, 'A-tract curvature', 'Olson 1998'),
            (r'T{8,}', 'CRV_003', 'Long T-tract', 'Local Curvature', 8,
             'curvature_score', 0.95, 'T-tract curvature', 'Olson 1998'),
        ]}
        return patterns

    # =========================
    # GLOBAL CURVATURE (STRICT)
    # =========================

    def find_a_tracts(self, sequence: str, minAT: int = None,
                      max_window: int = None) -> List[Dict[str, Any]]:

        seq = sequence.upper()
        if minAT is None:
            minAT = self.MIN_AT_TRACT

        tracts = []

        for m in re.finditer(r'A{' + str(minAT) + r',}', seq):
            start, end = m.start(), m.end()
            center = start + (end - start - 1) / 2.0
            tracts.append({'start': start, 'end': end, 'a_center': center})

        for m in re.finditer(r'T{' + str(minAT) + r',}', seq):
            start, end = m.start(), m.end()
            center = start + (end - start - 1) / 2.0
            tracts.append({'start': start, 'end': end, 'a_center': center})

        tracts.sort(key=lambda x: x['a_center'])
        return tracts

    def find_aprs(self, sequence: str, min_tract: int = None,
                  min_apr_tracts: int = None) -> List[Dict[str, Any]]:

        if min_apr_tracts is None:
            min_apr_tracts = self.MIN_APR_TRACTS

        tracts = self.find_a_tracts(sequence, minAT=min_tract)
        n = len(tracts)
        if n < min_apr_tracts:
            return []

        aprs = []
        run_start = 0

        for i in range(1, n):
            spacing = tracts[i]['a_center'] - tracts[i-1]['a_center']
            if not (self.PHASING_TOL_LOW <= spacing <= self.PHASING_TOL_HIGH):

                run_len = i - run_start
                if run_len >= min_apr_tracts:
                    run = tracts[run_start:i]
                    centers = [t['a_center'] for t in run]
                    spacings = [centers[j+1] - centers[j]
                                for j in range(len(centers)-1)]
                    mean_dev = sum(abs(s - self.PHASING_CENTER_SPACING)
                                   for s in spacings) / len(spacings)

                    max_dev = max(
                        abs(self.PHASING_TOL_HIGH -
                            self.PHASING_CENTER_SPACING),
                        abs(self.PHASING_CENTER_SPACING -
                            self.PHASING_TOL_LOW)
                    )

                    score = max(0.0, 1.0 - (mean_dev / max_dev))

                    aprs.append({
                        'center_positions': centers,
                        'tracts': run,
                        'n_tracts': run_len,
                        'spacings': spacings,
                        'mean_deviation': mean_dev,
                        'score': round(score, 6)
                    })

                run_start = i

        run_len = n - run_start
        if run_len >= min_apr_tracts:
            run = tracts[run_start:n]
            centers = [t['a_center'] for t in run]
            spacings = [centers[j+1] - centers[j]
                        for j in range(len(centers)-1)]
            mean_dev = sum(abs(s - self.PHASING_CENTER_SPACING)
                           for s in spacings) / len(spacings)

            max_dev = max(
                abs(self.PHASING_TOL_HIGH -
                    self.PHASING_CENTER_SPACING),
                abs(self.PHASING_CENTER_SPACING -
                    self.PHASING_TOL_LOW)
            )

            score = max(0.0, 1.0 - (mean_dev / max_dev))

            aprs.append({
                'center_positions': centers,
                'tracts': run,
                'n_tracts': run_len,
                'spacings': spacings,
                'mean_deviation': mean_dev,
                'score': round(score, 6)
            })

        return aprs

    # =========================
    # LOCAL CURVATURE (UNCHANGED)
    # =========================

    def find_long_tracts(self, sequence: str,
                         min_len: int = None) -> List[Dict[str, Any]]:

        if min_len is None:
            min_len = self.LOCAL_LONG_TRACT

        seq = sequence.upper()
        results = []

        for m in re.finditer(r'A{' + str(min_len) + r',}', seq):
            ln = m.end() - m.start()
            score = float(ln) / (ln + 7.0)
            results.append({
                'start': m.start(), 'end': m.end(),
                'base': 'A', 'len': ln,
                'score': round(score, 6),
                'seq': seq[m.start():m.end()]
            })

        for m in re.finditer(r'T{' + str(min_len) + r',}', seq):
            ln = m.end() - m.start()
            score = float(ln) / (ln + 7.0)
            results.append({
                'start': m.start(), 'end': m.end(),
                'base': 'T', 'len': ln,
                'score': round(score, 6),
                'seq': seq[m.start():m.end()]
            })

        results.sort(key=lambda x: x['start'])
        return results

    # =========================
    # ANNOTATION
    # =========================

    def detect_motifs(self, sequence: str,
                      sequence_name: str = "sequence") -> List[Dict[str, Any]]:

        sequence = sequence.upper().strip()
        motifs = []

        # Global Curvature: APR-based detection using true tract boundaries
        aprs = self.find_aprs(sequence,
                              min_tract=self.MIN_AT_TRACT,
                              min_apr_tracts=self.MIN_APR_TRACTS)

        for apr_idx, apr in enumerate(aprs):
            tract_windows = apr.get('tracts', [])
            if not tract_windows:
                continue

            start_pos = min(t['start'] for t in tract_windows)
            end_pos = max(t['end'] for t in tract_windows)

            if apr['score'] < self.SCORE_THRESHOLD:
                continue

            raw_score = apr['score']
            normalized_score = self.normalize_score(raw_score, end_pos - start_pos, 'Global Curvature')
            motif_seq = sequence[start_pos:end_pos]

            motifs.append({
                'ID': f"{sequence_name}_CRV_001_{start_pos + 1}_{apr_idx}",
                'Sequence_Name': sequence_name,
                'Class': self.get_motif_class_name(),
                'Subclass': 'Global Curvature',
                'Start': start_pos + 1,
                'End': end_pos,
                'Length': end_pos - start_pos,
                'Sequence': motif_seq,
                'Raw_Score': round(raw_score, 6),
                'Score': normalized_score,
                'Strand': '+',
                'Method': f'{self.get_motif_class_name()}_detection',
                'Pattern_ID': 'CRV_001',
                'A_Tracts': tract_windows,
                'N_Tracts': apr['n_tracts']
            })

        # Local Curvature: long single A/T tract detection via base class patterns
        local_motifs = super().detect_motifs(sequence, sequence_name)
        motifs.extend(local_motifs)

        return motifs

    def annotate_sequence(self, sequence: str) -> Dict[str, Any]:

        seq = sequence.upper()

        a_windows = self.find_a_tracts(seq, minAT=self.MIN_AT_TRACT)
        aprs = self.find_aprs(seq,
                              min_tract=self.MIN_AT_TRACT,
                              min_apr_tracts=self.MIN_APR_TRACTS)
        long_tracts = self.find_long_tracts(seq,
                                            min_len=self.LOCAL_LONG_TRACT)

        summary = {
            'n_a_windows': len(a_windows),
            'n_a_centers': len(a_windows),
            'n_aprs': len(aprs),
            'n_long_tracts': len(long_tracts),
            'apr_score_sum': sum(a['score'] for a in aprs),
            'long_tract_score_sum': sum(l['score']
                                        for l in long_tracts),
            'combined_score':
                sum(a['score'] for a in aprs) +
                sum(l['score'] for l in long_tracts)
        }

        return {
            'a_tract_windows': a_windows,
            'a_centers': a_windows,
            'aprs': aprs,
            'long_tracts': long_tracts,
            'summary': summary
        }
