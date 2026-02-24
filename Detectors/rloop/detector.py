"""R-Loop detector using QmRLFS algorithm (Jenjaroenpun 2016)."""
# IMPORTS
import re
from typing import List, Dict, Any, Tuple, Optional
from ..base.base_detector import BaseMotifDetector
from Utilities.detectors_utils import revcomp, calc_gc_content
from Utilities.core.motif_normalizer import normalize_class_subclass

try:
    import hyperscan
    HS_AVAILABLE = True
except Exception:
    HS_AVAILABLE = False

try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False

# TUNABLE PARAMETERS - QmRLFS Literature Parameters
MIN_PERC_G_RIZ = 50; NUM_LINKER = 50; WINDOW_STEP = 100
MAX_LENGTH_REZ = 2000; MIN_PERC_G_REZ = 40; QUALITY_THRESHOLD = 0.4


class RLoopDetector(BaseMotifDetector):
    """QmRLFS-finder R-loop detector (literature-faithful, accelerated)."""

    SCORE_REFERENCE = 'Aguilera et al. 2012, Jenjaroenpun et al. 2016'

    MIN_PERC_G_RIZ = MIN_PERC_G_RIZ; NUM_LINKER = NUM_LINKER; WINDOW_STEP = WINDOW_STEP
    MAX_LENGTH_REZ = MAX_LENGTH_REZ; MIN_PERC_G_REZ = MIN_PERC_G_REZ; QUALITY_THRESHOLD = QUALITY_THRESHOLD

    def __init__(self):
        super().__init__()
        self.hs_db = None
        self.hs_id_to_model = {}

        if HS_AVAILABLE:
            self._compile_hyperscan_patterns()

    # Required abstract methods

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """
        Return patterns for R-loop detection (QmRLFS algorithm).
        
        R-loops are detected using two models based on G-rich patterns:
        - Model 1: Standard G-cluster pattern for RIZ detection
        - Model 2: Extended G-tract pattern for RIZ detection
        """
        return {
            'qmrlfs_model_1': [
                (r'G{3,}[ATCG]{1,10}?G{3,}(?:[ATCG]{1,10}?G{3,}){1,}?',
                 'RLOOP_M1', 'QmRLFS Model 1', 'R-loop formation sites',
                 12, 'qmrlfs_score', 0.4, 'G-cluster RIZ pattern', 'Jenjaroenpun 2016')
            ],
            'qmrlfs_model_2': [
                (r'G{4,}(?:[ATCG]{1,10}?G{4,}){1,}?',
                 'RLOOP_M2', 'QmRLFS Model 2', 'R-loop formation sites',
                 8, 'qmrlfs_score', 0.4, 'Extended G-tract RIZ pattern', 'Jenjaroenpun 2016')
            ]
        }

    def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:
        annotations = self.annotate_sequence(sequence)
        if not annotations:
            return 0.0

        best = max(
            (ann['riz_perc_g'] / 100.0) +
            (ann.get('rez_perc_g', 0) / 100.0)
            for ann in annotations
        )

        return min(1.0, round(best, 6))

    def passes_quality_threshold(self,
                                  sequence: str,
                                  score: float,
                                  pattern_info: Tuple = None) -> bool:
        return score >= self.QUALITY_THRESHOLD


    def get_motif_class_name(self) -> str:
        return "R-Loop"

    def get_length_cap(self, subclass: str = None) -> int:
        """R-loops extend up to ~2000 bp (Aguilera 2012, Jenjaroenpun 2016)."""
        return 2000

    def theoretical_min_score(self) -> float:
        """Minimum biologically valid R-loop raw score (quality threshold)."""
        return self.QUALITY_THRESHOLD

    def theoretical_max_score(self, sequence_length: int = None) -> float:
        """Highest possible R-loop raw score.

        Score = riz_perc_g/100 + rez_perc_g/100, capped at 1.0.
        Max G% = 100% for both zones → max combined = 2.0, capped at 1.0.
        """
        return 1.0

    # Hyperscan Compilation

    def _compile_hyperscan_patterns(self):
        try:
            expressions = [
                br"G{3,}[ATCG]{1,10}?G{3,}(?:[ATCG]{1,10}?G{3,}){1,}?",
                br"G{4,}(?:[ATCG]{1,10}?G{4,}){1,}?"
            ]
            ids = [1, 2]
            flags = [hyperscan.HS_FLAG_DOTALL] * 2

            self.hs_db = hyperscan.Database()
            self.hs_db.compile(expressions=expressions,
                               ids=ids,
                               flags=flags)

            self.hs_id_to_model = {1: 'qmrlfs_model_1',
                                   2: 'qmrlfs_model_2'}
        except Exception:
            self.hs_db = None

    # RIZ Detection

    def _riz_search(self, seq: str, model: str) -> List[Dict[str, Any]]:

        results = []

        if HS_AVAILABLE and self.hs_db is not None:

            seq_bytes = seq.encode()

            def on_match(id, start, end, flags, context):
                if self.hs_id_to_model.get(id) == model:
                    riz_seq = seq[start:end]
                    if self._percent_g(riz_seq) >= self.MIN_PERC_G_RIZ:
                        results.append({
                            'start': start,
                            'end': end,
                            'sequence': riz_seq
                        })
                return 0

            self.hs_db.scan(seq_bytes, match_event_handler=on_match)

        else:
            if model == 'qmrlfs_model_1':
                pattern = re.compile(
                    r"G{3,}[ATCG]{1,10}?G{3,}(?:[ATCG]{1,10}?G{3,}){1,}?",
                    re.IGNORECASE
                )
            else:
                pattern = re.compile(
                    r"G{4,}(?:[ATCG]{1,10}?G{4,}){1,}?",
                    re.IGNORECASE
                )

            for m in pattern.finditer(seq):
                riz_seq = m.group(0)
                if self._percent_g(riz_seq) >= self.MIN_PERC_G_RIZ:
                    results.append({
                        'start': m.start(),
                        'end': m.end(),
                        'sequence': riz_seq
                    })

        return results


    def _find_rez(self,
                  seq: str,
                  riz_end: int) -> Optional[Dict[str, Any]]:

        seq_len = len(seq)
        search_start = riz_end + self.NUM_LINKER
        if search_start >= seq_len:
            return None

        if NUMPY_AVAILABLE and seq_len > 1000:
            seq_array = np.array([1 if c == 'G' else 0 for c in seq], dtype=np.int32)
            prefix_g = np.concatenate(([0], np.cumsum(seq_array)))
        else:
            prefix_g = [0] * (seq_len + 1)
            for i in range(seq_len):
                prefix_g[i + 1] = prefix_g[i] + (1 if seq[i] == 'G' else 0)

        best = None
        max_score = 0.0

        max_end = min(seq_len, riz_end + self.MAX_LENGTH_REZ)

        for start in range(search_start, max_end, self.WINDOW_STEP):

            window_seed_end = min(start + 100, seq_len)
            seed_g = prefix_g[window_seed_end] - prefix_g[start]
            seed_len = window_seed_end - start

            if seed_len == 0:
                continue

            if (seed_g / seed_len) * 100 < self.MIN_PERC_G_REZ:
                continue

            for end in range(start + 50,
                             max_end,
                             50):

                length = end - start
                g_count = prefix_g[end] - prefix_g[start]
                perc_g = (g_count / length) * 100

                if perc_g >= self.MIN_PERC_G_REZ:
                    score = perc_g * length / 100.0
                    if score > max_score:
                        max_score = score
                        best = {
                            'start': start,
                            'end': end,
                            'length': length,
                            'sequence': seq[start:end],
                            'perc_g': round(perc_g, 2)
                        }

        return best


    def _percent_g(self, seq: str) -> float:
        return round((seq.count("G") / float(len(seq))) * 100.0, 2) if seq else 0.0


    def annotate_sequence(self,
                          sequence: str,
                          models: Optional[List[str]] = None
                          ) -> List[Dict[str, Any]]:

        seq = sequence.upper()
        models = models or ['qmrlfs_model_1', 'qmrlfs_model_2']
        results = []

        for model in models:

            riz_regions = self._riz_search(seq, model)

            for riz in riz_regions:

                rez = self._find_rez(seq, riz['end'])

                result = {
                    'model': model,
                    'riz_start': riz['start'],
                    'riz_end': riz['end'],
                    'riz_length': riz['end'] - riz['start'],
                    'riz_sequence': riz['sequence'],
                    'riz_perc_g': self._percent_g(riz['sequence']),
                    'total_start': riz['start'],
                    'total_end': riz['end'],
                    'total_length': riz['end'] - riz['start']
                }

                if rez:
                    result.update({
                        'rez_start': rez['start'],
                        'rez_end': rez['end'],
                        'rez_length': rez['length'],
                        'rez_sequence': rez['sequence'],
                        'rez_perc_g': rez['perc_g'],
                        'total_end': rez['end'],
                        'total_length': rez['end'] - riz['start']
                    })

                results.append(result)

        return results


    def detect_motifs(self,
                      sequence: str,
                      sequence_name: str = "sequence"
                      ) -> List[Dict[str, Any]]:

        self.audit['invoked'] = True
        self.audit['windows_scanned'] = 2
        self.audit['candidates_seen'] = 0
        self.audit['reported'] = 0
        self.audit['both_strands_scanned'] = True

        motifs = []

        for strand, seq in [('+', sequence.upper()),
                            ('-', revcomp(sequence.upper()))]:

            annotations = self.annotate_sequence(seq)
            self.audit['candidates_seen'] += len(annotations)

            seq_len = len(sequence)

            for i, ann in enumerate(annotations):

                if strand == '-':
                    start = seq_len - ann['total_end']
                    end = seq_len - ann['total_start']
                else:
                    start = ann['total_start']
                    end = ann['total_end']

                score = (
                    ann['riz_perc_g'] / 100.0 +
                    ann.get('rez_perc_g', 0) / 100.0
                )

                motif_seq = sequence[start:end]
                gc_content = round(calc_gc_content(motif_seq), 2)
                
                riz_len = ann.get('riz_length', 0)
                rez_len = ann.get('rez_length', 0)
                total_len = end - start
                linker_len = max(0, total_len - riz_len - rez_len)
                
                g_count = motif_seq.count('G')
                c_count = motif_seq.count('C')
                gc_skew = round((g_count - c_count) / (g_count + c_count), 3) if (g_count + c_count) > 0 else 0.0

                canonical_class, canonical_subclass = normalize_class_subclass(
                    self.get_motif_class_name(),
                    'R-loop formation sites',
                    strict=False,
                    auto_correct=True
                )

                motifs.append({
                    'ID': f"{sequence_name}_RLOOP_{start+1}",
                    'Sequence_Name': sequence_name,
                    'Class': canonical_class,
                    'Subclass': canonical_subclass,
                    'Start': start + 1,
                    'End': end,
                    'Length': end - start,
                    'Sequence': motif_seq,
                    'Raw_Score': round(min(score, 1.0), 3),
                    'Score': self.normalize_score(min(score, 1.0), end - start, canonical_subclass),
                    'Strand': strand,
                    'Method': 'QmRLFS_detection',
                    'Pattern_ID': f'RLOOP_{ann["model"]}_{i+1}',
                    'Model': ann['model'],
                    'RIZ_Length': ann.get('riz_length', 0),
                    'RIZ_Perc_G': ann.get('riz_perc_g', 0),
                    'REZ_Length': ann.get('rez_length', 0),
                    'REZ_Perc_G': ann.get('rez_perc_g', 0),
                    'Linker_Length': linker_len,
                    'GC_Content': gc_content,
                    'GC_Skew': gc_skew,
                    'Arm_Length': 'N/A',
                    'Loop_Length': linker_len if linker_len > 0 else 0,
                    'Type_Of_Repeat': 'RNA-DNA hybrid (R-loop)',
                    'Criterion': self._get_rloop_criterion(ann),
                    'Disease_Relevance': self._get_rloop_disease_relevance(gc_content, total_len, ann),
                    'Regions_Involved': self._describe_rloop_regions(riz_len, rez_len, linker_len)
                })

                self.audit['reported'] += 1

        return motifs
    
    def _get_rloop_criterion(self, ann: Dict[str, Any]) -> str:
        """Explain R-loop classification criterion"""
        criteria = []
        model = ann.get('model', 'unknown')
        criteria.append(f'QmRLFS {model}')
        
        riz_perc_g = ann.get('riz_perc_g', 0)
        rez_perc_g = ann.get('rez_perc_g', 0)
        
        criteria.append(f'RIZ %G={riz_perc_g:.1f}% ≥{self.MIN_PERC_G_RIZ}%')
        if rez_perc_g > 0:
            criteria.append(f'REZ %G={rez_perc_g:.1f}% ≥{self.MIN_PERC_G_REZ}%')
        
        criteria.append(f'Score threshold ≥{self.QUALITY_THRESHOLD}')
        
        return '; '.join(criteria)
    
    def _get_rloop_disease_relevance(self, gc_content: float, length: int, ann: Dict[str, Any]) -> str:
        """Annotate disease relevance for R-loops"""
        disease_notes = []
        
        riz_perc_g = ann.get('riz_perc_g', 0)
        
        if gc_content > 70:
            disease_notes.append('High GC-content - genomic instability, DNA damage hotspot')
        
        if length > 500:
            disease_notes.append('Long R-loop (>500bp) - replication-transcription conflicts')
        
        if riz_perc_g > 65:
            disease_notes.append('Strong RIZ signal - transcription-associated R-loop')
        
        if not disease_notes:
            disease_notes.append('R-loop formation - potential role in gene regulation, DNA damage, genome instability')
        
        disease_notes.append('Associated with: neurodegeneration (ALS, Fragile X), cancer, repeat expansion diseases')
        
        return '; '.join(disease_notes)
    
    def _describe_rloop_regions(self, riz_len: int, rez_len: int, linker_len: int) -> str:
        """Describe regions involved in R-loop formation"""
        parts = []
        
        if riz_len > 0:
            parts.append(f'RIZ (RNA invasion zone, {riz_len}bp, G-rich)')
        
        if linker_len > 0:
            parts.append(f'Linker ({linker_len}bp)')
        
        if rez_len > 0:
            parts.append(f'REZ (RNA exit zone, {rez_len}bp)')
        
        return ' - '.join(parts) if parts else 'R-loop forming region'
