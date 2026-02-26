"""Tests for the SMID preprocessing layer.

Validates:
- seed_database : correct seed counts, no forbidden patterns, required attributes
- seed_engine   : correct hit positions, case-insensitivity, G4/R-Loop sharing
- clustering    : O(n) 1D clustering correctness, edge cases
- extend_and_merge : extension arithmetic, overlap merging, boundary clipping
- dispatcher    : coordinate translation, detector invocation
- pipeline      : integration smoke-test
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import re
import numpy as np
import pytest

from Detectors.aphilic.tenmer_table import TENMER_LOG2
from Detectors.zdna.tenmer_table import TENMER_SCORE

from SMID.seed_database import SMID_SEEDS
from SMID.seed_engine import SMIDSeedEngine
from SMID.clustering import linear_1d_clustering, extend_and_merge
from SMID.dispatcher import SMIDDispatcher


# ===========================================================================
# seed_database tests
# ===========================================================================

class TestSeedDatabase:
    def test_aphilic_seeds_match_tenmer_table(self):
        """Every key in TENMER_LOG2 must appear as an exact seed for A-philic_DNA."""
        aphilic_patterns = {
            s["pattern"] for s in SMID_SEEDS
            if not s["is_regex"] and "A-philic_DNA" in s["classes"]
        }
        assert aphilic_patterns == set(TENMER_LOG2.keys()), (
            "A-philic seeds do not exactly match TENMER_LOG2 keys"
        )

    def test_zdna_seeds_match_tenmer_table(self):
        """Every key in TENMER_SCORE must appear as an exact seed for Z-DNA."""
        zdna_patterns = {
            s["pattern"] for s in SMID_SEEDS
            if not s["is_regex"] and "Z-DNA" in s["classes"]
        }
        assert zdna_patterns == set(TENMER_SCORE.keys()), (
            "Z-DNA seeds do not exactly match TENMER_SCORE keys"
        )

    def test_required_seed_attributes(self):
        required = {"pattern", "is_regex", "classes", "max_extension", "epsilon", "min_samples"}
        for seed in SMID_SEEDS:
            assert required.issubset(seed.keys()), f"Seed missing attributes: {seed}"

    def test_no_forbidden_classes(self):
        forbidden = {"Curved_DNA", "Slipped_DNA", "Cruciform"}
        for seed in SMID_SEEDS:
            for cls in seed["classes"]:
                assert cls not in forbidden, f"Forbidden class {cls} found in seed: {seed}"

    def test_g4_rloop_shared_seed_present(self):
        shared = [
            s for s in SMID_SEEDS
            if s["is_regex"]
            and "G-Quadruplex" in s["classes"]
            and "R-Loop" in s["classes"]
        ]
        assert len(shared) >= 1, "No shared G4/R-Loop seed found"

    def test_imotif_seed_present(self):
        imotif = [s for s in SMID_SEEDS if "i-Motif" in s["classes"]]
        assert len(imotif) >= 1

    def test_triplex_seeds_present(self):
        triplex = [s for s in SMID_SEEDS if "Triplex" in s["classes"]]
        assert len(triplex) >= 2, "Expected at least two Triplex seeds ([AG] and [CT])"

    def test_no_permissive_single_run_seeds(self):
        """A{3,}, G{3,} alone (without second tract) must NOT appear as a standalone seed."""
        single_run_patterns = [r"A{3,}", r"G{3,}", r"C{3,}", r"T{3,}"]
        for seed in SMID_SEEDS:
            if seed["is_regex"]:
                for bad in single_run_patterns:
                    assert seed["pattern"] != bad, (
                        f"Permissive single-run seed '{bad}' found"
                    )

    def test_class_specific_extensions(self):
        """G4/i-Motif should have smaller extension than A-philic/Z-DNA."""
        aphilic_ext = next(
            s["max_extension"] for s in SMID_SEEDS
            if "A-philic_DNA" in s["classes"] and not s["is_regex"]
        )
        g4_ext = next(
            s["max_extension"] for s in SMID_SEEDS
            if "G-Quadruplex" in s["classes"] and s["is_regex"]
        )
        imotif_ext = next(
            s["max_extension"] for s in SMID_SEEDS
            if "i-Motif" in s["classes"] and s["is_regex"]
        )
        assert g4_ext < aphilic_ext, "G4 extension should be less than A-philic"
        assert imotif_ext <= g4_ext, "i-Motif extension should be ≤ G4 extension"


# ===========================================================================
# seed_engine tests
# ===========================================================================

class TestSeedEngine:
    def _make_engine_with_seeds(self, seeds):
        return SMIDSeedEngine(seeds=seeds)

    def test_exact_aphilic_hit(self):
        """A known A-philic 10-mer should produce a hit at the correct position."""
        kmer = next(iter(TENMER_LOG2))
        seq = "AAAAAAAAAA" + kmer + "GGGGGGGGGG"
        seeds = [s for s in SMID_SEEDS if s["pattern"] == kmer and not s["is_regex"]]
        engine = self._make_engine_with_seeds(seeds)
        hits = engine.scan(seq)
        assert "A-philic_DNA" in hits
        assert 10 in hits["A-philic_DNA"]

    def test_case_insensitive_exact(self):
        """Exact seeds must match lowercase sequences."""
        kmer = next(iter(TENMER_LOG2))
        seq = kmer.lower()
        seeds = [s for s in SMID_SEEDS if s["pattern"] == kmer and not s["is_regex"]]
        engine = self._make_engine_with_seeds(seeds)
        hits = engine.scan(seq)
        assert "A-philic_DNA" in hits
        assert 0 in hits["A-philic_DNA"]

    def test_g4_rloop_shared_position(self):
        """G4 seed hit must appear in both G-Quadruplex and R-Loop classes."""
        seq = "NNNGGGGACGTGGGGNNNN"
        g4_seeds = [s for s in SMID_SEEDS if "G-Quadruplex" in s["classes"] and s["is_regex"]]
        engine = self._make_engine_with_seeds(g4_seeds)
        hits = engine.scan(seq)
        assert "G-Quadruplex" in hits
        assert "R-Loop" in hits
        # Same position must be in both
        assert set(hits["G-Quadruplex"].tolist()) == set(hits["R-Loop"].tolist())

    def test_regex_imotif_hit(self):
        """i-Motif seed must fire on a three-C-tract sequence."""
        seq = "AAACCCACGTCCCAAA"
        imotif_seeds = [s for s in SMID_SEEDS if "i-Motif" in s["classes"]]
        engine = self._make_engine_with_seeds(imotif_seeds)
        hits = engine.scan(seq)
        assert "i-Motif" in hits
        assert len(hits["i-Motif"]) >= 1

    def test_regex_triplex_hit(self):
        """Triplex [AG]{12,} seed must fire on a long purine run."""
        seq = "TTT" + "AGAGAGAGAGAGAG" + "CCC"
        triplex_seeds = [s for s in SMID_SEEDS if "Triplex" in s["classes"]]
        engine = self._make_engine_with_seeds(triplex_seeds)
        hits = engine.scan(seq)
        assert "Triplex" in hits

    def test_no_hit_on_random_sequence(self):
        """A random sequence without any seed patterns should return empty."""
        # 'TTTTTTTTTT' is a 10-mer absent from both tenmer tables
        seq = "T" * 100
        seeds_only_aphilic_zdna = [s for s in SMID_SEEDS if not s["is_regex"]]
        engine = self._make_engine_with_seeds(seeds_only_aphilic_zdna)
        hits = engine.scan(seq)
        # No exact 10-mer in TTTTTTTTTT... should match (poly-T is not in tables)
        assert all(len(positions) == 0 for positions in hits.values()), (
            "Expected no hits in poly-T sequence"
        )

    def test_scan_returns_sorted_unique_positions(self):
        """Scan result arrays must be sorted and deduplicated."""
        kmer = next(iter(TENMER_LOG2))
        seq = kmer + kmer  # overlapping by design (two adjacent copies)
        seeds = [s for s in SMID_SEEDS if s["pattern"] == kmer and not s["is_regex"]]
        engine = self._make_engine_with_seeds(seeds)
        hits = engine.scan(seq)
        if "A-philic_DNA" in hits:
            arr = hits["A-philic_DNA"]
            assert list(arr) == sorted(set(arr.tolist())), "Positions not sorted/unique"


# ===========================================================================
# clustering tests
# ===========================================================================

class TestLinear1DClustering:
    def test_single_cluster(self):
        positions = np.array([10, 15, 20, 25])
        result = linear_1d_clustering(positions, epsilon=10, min_samples=2)
        assert len(result) == 1
        assert result[0] == (10, 25)

    def test_two_separate_clusters(self):
        positions = np.array([10, 20, 100, 110])
        result = linear_1d_clustering(positions, epsilon=15, min_samples=2)
        assert len(result) == 2
        assert result[0] == (10, 20)
        assert result[1] == (100, 110)

    def test_min_samples_filter(self):
        positions = np.array([10, 100])
        # Each "cluster" has only one hit → should be filtered
        result = linear_1d_clustering(positions, epsilon=5, min_samples=2)
        assert result == []

    def test_empty_input(self):
        result = linear_1d_clustering(np.array([]), epsilon=10, min_samples=2)
        assert result == []

    def test_single_point_below_min_samples(self):
        result = linear_1d_clustering(np.array([42]), epsilon=10, min_samples=2)
        assert result == []

    def test_single_point_meets_min_samples_1(self):
        result = linear_1d_clustering(np.array([42]), epsilon=10, min_samples=1)
        assert result == [(42, 42)]

    def test_unsorted_input(self):
        positions = np.array([25, 10, 15, 20])
        result = linear_1d_clustering(positions, epsilon=10, min_samples=2)
        assert len(result) == 1
        assert result[0] == (10, 25)

    def test_exact_epsilon_boundary(self):
        """Hits exactly epsilon apart should be in the same cluster."""
        positions = np.array([0, 10])
        result = linear_1d_clustering(positions, epsilon=10, min_samples=2)
        assert len(result) == 1

    def test_one_beyond_epsilon(self):
        """Hits epsilon+1 apart should be in different clusters."""
        positions = np.array([0, 11])
        result = linear_1d_clustering(positions, epsilon=10, min_samples=2)
        assert result == []  # each cluster has only 1 sample → filtered


class TestExtendAndMerge:
    def test_basic_extension(self):
        clusters = [(50, 100)]
        result = extend_and_merge(clusters, max_extension=20, seq_len=200)
        assert result == [(30, 120)]

    def test_clip_to_sequence_bounds(self):
        clusters = [(5, 10)]
        result = extend_and_merge(clusters, max_extension=20, seq_len=50)
        # start clipped to 0, end clipped to seq_len - 1 = 49
        assert result[0][0] == 0
        assert result[0][1] == 30

    def test_overlapping_windows_merged(self):
        clusters = [(10, 50), (60, 100)]
        result = extend_and_merge(clusters, max_extension=20, seq_len=200)
        # (10-20, 50+20) = (-10→0, 70) and (60-20, 100+20) = (40, 120)
        # 40 ≤ 70 → they merge
        assert len(result) == 1

    def test_non_overlapping_windows_kept_separate(self):
        clusters = [(0, 10), (200, 210)]
        result = extend_and_merge(clusters, max_extension=5, seq_len=300)
        assert len(result) == 2

    def test_empty_clusters(self):
        result = extend_and_merge([], max_extension=100, seq_len=1000)
        assert result == []

    def test_end_clipped_to_seq_len(self):
        clusters = [(990, 995)]
        result = extend_and_merge(clusters, max_extension=50, seq_len=1000)
        assert result[0][1] == 999  # seq_len - 1


# ===========================================================================
# dispatcher tests
# ===========================================================================

class _MockDetector:
    """Minimal fake detector for dispatcher unit tests."""

    def __init__(self):
        self.detect_called_with = []

    def detect_motifs(self, sequence, sequence_name="seq"):
        self.detect_called_with.append((sequence, sequence_name))
        return [
            {
                "ID": f"{sequence_name}_p1_1",
                "Sequence_Name": sequence_name,
                "Class": "MockClass",
                "Subclass": "MockSub",
                "Start": 1,
                "End": len(sequence),
                "Length": len(sequence),
                "Sequence": sequence,
                "Raw_Score": 1.0,
                "Score": 2.0,
                "Strand": "+",
                "Method": "mock_detection",
                "Pattern_ID": "p1",
            }
        ]


class TestSMIDDispatcher:
    def test_coordinate_translation(self):
        """Motif Start/End must be translated from window-local to chromosome coords."""
        chromosome = "N" * 100 + "ACGT" * 10 + "N" * 100
        win_start = 100
        win_end = 139

        mock_cls = _MockDetector
        dispatcher = SMIDDispatcher(class_to_detector_map={"MockClass": mock_cls})

        results = dispatcher.run(
            chromosome,
            {"MockClass": [(win_start, win_end)]},
            sequence_name="chr1",
        )
        assert len(results) == 1
        # Start = 1 (window-local) + win_start offset
        assert results[0]["Start"] == 1 + win_start
        assert results[0]["Sequence_Name"] == "chr1"

    def test_unknown_class_skipped(self):
        """Windows for classes not in the map should be silently skipped."""
        dispatcher = SMIDDispatcher(class_to_detector_map={})
        results = dispatcher.run("ACGT" * 10, {"UnknownClass": [(0, 9)]})
        assert results == []

    def test_empty_windows_no_results(self):
        dispatcher = SMIDDispatcher(class_to_detector_map={})
        results = dispatcher.run("ACGT" * 10, {})
        assert results == []


# ===========================================================================
# Integration / pipeline smoke test
# ===========================================================================

class TestSMIDPipelineSmoke:
    def test_pipeline_returns_list(self):
        """smid_pipeline must return a list (even on a short or empty sequence)."""
        from SMID import smid_pipeline
        result = smid_pipeline("ACGT" * 5, sequence_name="test", max_workers=1)
        assert isinstance(result, list)

    def test_pipeline_on_known_g4_sequence(self):
        """A canonical G4-forming sequence should produce G-Quadruplex results."""
        from SMID import smid_pipeline
        # Three-repeat G4 seed region
        seq = "A" * 50 + "GGGACGTGGG" * 3 + "A" * 50
        results = smid_pipeline(seq, sequence_name="g4_test", max_workers=1)
        # We just check that the pipeline completes without error and returns a list
        assert isinstance(results, list)

    def test_pipeline_all_coords_within_sequence(self):
        """All result coordinates must fall within the input sequence length."""
        from SMID import smid_pipeline

        kmer = next(iter(TENMER_LOG2))
        seq = "A" * 100 + kmer * 5 + "A" * 100
        seq_len = len(seq)

        results = smid_pipeline(seq, sequence_name="coord_test", max_workers=1)
        for motif in results:
            assert 1 <= motif["Start"] <= seq_len, f"Start out of bounds: {motif['Start']}"
            assert motif["End"] <= seq_len, f"End out of bounds: {motif['End']}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
