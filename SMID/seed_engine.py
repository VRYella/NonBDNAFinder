"""SMID Seed Engine: multi-pattern scanner over seed database.

Uses pyahocorasick (already in requirements.txt) for exact-match 10-mer seeds
and the built-in ``re`` module for regex seeds.  Falls back to a pure-Python
linear scan when pyahocorasick is unavailable.

Public interface
----------------
    engine = SMIDSeedEngine()
    hits   = engine.scan(sequence)      # -> Dict[class_name, np.ndarray]

The returned dict maps each class name to a sorted numpy array of 0-based
start positions of seed hits.  When a seed belongs to multiple classes (e.g.
the shared G4/R-Loop seed) the *same* position array is stored under each
class name — no coordinate is duplicated.
"""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Dict, List

import numpy as np

from SMID.seed_database import SMID_SEEDS

try:
    import ahocorasick
    _AHO_AVAILABLE = True
except ImportError:
    _AHO_AVAILABLE = False


class SMIDSeedEngine:
    """Compile seed patterns once and scan sequences efficiently.

    Parameters
    ----------
    seeds : list, optional
        Override the default SMID_SEEDS for testing purposes.
    """

    def __init__(self, seeds: List[dict] | None = None):
        self._seeds = seeds if seeds is not None else SMID_SEEDS

        # Separate exact-match seeds from regex seeds
        self._exact_seeds: List[dict] = [s for s in self._seeds if not s["is_regex"]]
        self._regex_seeds: List[dict] = [s for s in self._seeds if s["is_regex"]]

        # Build Aho-Corasick automaton for exact seeds (case-insensitive)
        self._aho = self._build_aho(self._exact_seeds)

        # Pre-compile regex seeds (case-insensitive, ASCII)
        self._compiled_regex: List[tuple] = [
            (re.compile(s["pattern"], re.IGNORECASE | re.ASCII), s)
            for s in self._regex_seeds
        ]

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_aho(self, exact_seeds: List[dict]):
        """Return an Aho-Corasick automaton, or None if unavailable."""
        if not _AHO_AVAILABLE or not exact_seeds:
            return None

        A = ahocorasick.Automaton()
        for idx, seed in enumerate(exact_seeds):
            key = seed["pattern"].upper()
            # Multiple seeds can share the same k-mer (e.g. AT-rich 10-mers
            # that appear in both A-philic and another table): store a list.
            if key in A:
                A.get(key).append(idx)
            else:
                A.add_word(key, [idx])
        A.make_automaton()
        return A

    def _scan_exact_aho(self, seq_upper: str) -> Dict[str, List[int]]:
        """Return exact-match hit positions grouped by class using Aho-Corasick."""
        positions: Dict[str, List[int]] = defaultdict(list)
        kmer_len = 10  # all exact seeds are 10-mers
        for end_idx, indices in self._aho.iter(seq_upper):
            start = end_idx - kmer_len + 1
            for idx in indices:
                for cls in self._exact_seeds[idx]["classes"]:
                    positions[cls].append(start)
        return positions

    def _scan_exact_linear(self, seq_upper: str) -> Dict[str, List[int]]:
        """Fallback exact-match scan when pyahocorasick is unavailable."""
        positions: Dict[str, List[int]] = defaultdict(list)
        for seed in self._exact_seeds:
            pat = seed["pattern"].upper()
            start = 0
            while True:
                idx = seq_upper.find(pat, start)
                if idx == -1:
                    break
                for cls in seed["classes"]:
                    positions[cls].append(idx)
                start = idx + 1
        return positions

    def _scan_regex(self, seq_upper: str) -> Dict[str, List[int]]:
        """Return regex hit positions grouped by class.

        Shared-class seeds (e.g. G4/R-Loop) record the hit position once per
        class — coordinates are identical so no duplication occurs.
        """
        positions: Dict[str, List[int]] = defaultdict(list)
        for compiled_re, seed in self._compiled_regex:
            # Collect start positions for this seed in one pass
            hit_positions = [m.start() for m in compiled_re.finditer(seq_upper)]
            for cls in seed["classes"]:
                positions[cls].extend(hit_positions)
        return positions

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def scan(self, sequence: str) -> Dict[str, np.ndarray]:
        """Scan *sequence* for all seeds and return hit positions by class.

        Parameters
        ----------
        sequence : str
            DNA sequence (any case; internally uppercased).

        Returns
        -------
        dict
            ``{class_name: np.ndarray}`` of sorted 0-based start positions.
        """
        seq_upper = sequence.upper()

        # Gather positions from both backends
        if self._aho is not None:
            exact_hits = self._scan_exact_aho(seq_upper)
        else:
            exact_hits = self._scan_exact_linear(seq_upper)

        regex_hits = self._scan_regex(seq_upper)

        # Merge and convert to sorted numpy arrays
        all_classes = set(exact_hits) | set(regex_hits)
        result: Dict[str, np.ndarray] = {}
        for cls in all_classes:
            combined = exact_hits.get(cls, []) + regex_hits.get(cls, [])
            result[cls] = np.array(sorted(set(combined)), dtype=np.int64)

        return result
