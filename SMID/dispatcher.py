"""SMID Dispatcher: run detectors on pre-computed windows.

The dispatcher receives a mapping of class → windows and instantiates the
appropriate detector for each window.  It does NOT modify detector internals;
it only slices the sequence and translates result coordinates back to the
original chromosome frame.

Public interface
----------------
    dispatcher = SMIDDispatcher(class_to_detector_map)
    results    = dispatcher.run(chromosome_seq, windows_by_class)
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

# Default class → detector callable mapping.
# Importing lazily to avoid import-time overhead when SMID is used as a
# library without all detectors present.
_DEFAULT_DETECTOR_MAP: Dict[str, Any] | None = None


def _build_default_detector_map() -> Dict[str, Any]:
    """Build the default class-to-detector map on first use."""
    from Detectors.aphilic.detector import APhilicDetector
    from Detectors.zdna.detector import ZDNADetector
    from Detectors.gquad.detector import GQuadruplexDetector
    from Detectors.rloop.detector import RLoopDetector
    from Detectors.imotif.detector import IMotifDetector
    from Detectors.triplex.detector import TriplexDetector

    return {
        "A-philic_DNA": APhilicDetector,
        "Z-DNA": ZDNADetector,
        "G-Quadruplex": GQuadruplexDetector,
        "R-Loop": RLoopDetector,
        "i-Motif": IMotifDetector,
        "Triplex": TriplexDetector,
    }


class SMIDDispatcher:
    """Dispatch windowed sequences to the matching detector.

    Parameters
    ----------
    class_to_detector_map : dict, optional
        Mapping of ``class_name -> DetectorClass``.  If *None*, the default
        map (all SMID-supported detectors) is used.  Pass a custom map for
        testing or partial-genome runs.
    """

    def __init__(self, class_to_detector_map: Dict[str, Any] | None = None):
        if class_to_detector_map is not None:
            self._map = class_to_detector_map
        else:
            global _DEFAULT_DETECTOR_MAP
            if _DEFAULT_DETECTOR_MAP is None:
                _DEFAULT_DETECTOR_MAP = _build_default_detector_map()
            self._map = _DEFAULT_DETECTOR_MAP

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(
        self,
        chromosome_seq: str,
        windows_by_class: Dict[str, List[Tuple[int, int]]],
        sequence_name: str = "seq",
    ) -> List[Dict[str, Any]]:
        """Run detectors over windows and return coordinate-adjusted results.

        Parameters
        ----------
        chromosome_seq : str
            Full chromosome (or contig) sequence.
        windows_by_class : dict
            ``{class_name: [(start, end), ...]}`` — 0-based, end-inclusive
            windows to scan for each class.
        sequence_name : str, optional
            Label used in motif ``ID`` and ``Sequence_Name`` fields.

        Returns
        -------
        list of dict
            Merged motif records from all detectors, with ``Start`` / ``End``
            translated back to chromosome coordinates (1-based, matching
            detector convention).
        """
        all_results: List[Dict[str, Any]] = []
        seq_upper = chromosome_seq.upper()

        for class_name, windows in windows_by_class.items():
            detector_cls = self._map.get(class_name)
            if detector_cls is None:
                continue

            for win_start, win_end in windows:
                window_seq = seq_upper[win_start: win_end + 1]
                if not window_seq:
                    continue

                all_results.extend(
                    self.run_window(window_seq, class_name, win_start, sequence_name)
                )

        return all_results

    def run_window(
        self,
        window_seq: str,
        class_name: str,
        win_start: int,
        sequence_name: str = "seq",
    ) -> List[Dict[str, Any]]:
        """Run the detector for *class_name* on a pre-extracted window sequence.

        Parameters
        ----------
        window_seq : str
            Already-extracted window sequence (0-based slice of the chromosome).
        class_name : str
            Motif class to detect.
        win_start : int
            0-based start position of the window in the chromosome; used to
            translate coordinates back to chromosome space.
        sequence_name : str, optional
            Label for result records.

        Returns
        -------
        list of dict
            Motif records with ``Start`` / ``End`` in chromosome coordinates
            (1-based, matching detector convention).
        """
        detector_cls = self._map.get(class_name)
        if detector_cls is None or not window_seq:
            return []

        detector = detector_cls()
        window_motifs = detector.detect_motifs(
            window_seq.upper(),
            sequence_name=f"{sequence_name}_w{win_start}",
        )

        for motif in window_motifs:
            motif["Start"] = motif["Start"] + win_start
            motif["End"] = motif["End"] + win_start
            motif["ID"] = (
                f"{sequence_name}_{motif.get('Pattern_ID', 'p')}"
                f"_{motif['Start']}"
            )
            motif["Sequence_Name"] = sequence_name

        return window_motifs
