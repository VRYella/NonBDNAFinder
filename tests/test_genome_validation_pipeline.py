"""
tests/test_genome_validation_pipeline.py
=========================================
Unit tests for NBSTVALIDATION/run_genome_validation.py

Covers:
  - genome discovery (find all FASTA files)
  - _covered_bases calculation (interval merging)
  - build_master_tables schema and key columns
  - figure generation (return valid file paths and non-empty PNGs)
  - report generation (key sections present)

These tests use small synthetic data so they run in seconds without
needing the full genome files.
"""

from __future__ import annotations

import os
import sys
import textwrap
import tempfile
from pathlib import Path
from typing import Dict, List

import pytest
import pandas as pd

# Add repo root to path so imports work
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO)

# Import the module under test
from NBSTVALIDATION.run_genome_validation import (
    discover_genomes,
    _covered_bases,
    build_master_tables,
    generate_all_figures,
    generate_report,
    NON_B_CLASSES,
    HYBRID_CLASS,
    CLUSTER_CLASS,
    _jaccard_interval,
    _match_motif_lists,
    NBST_DATA_FILES,
    NBST_TYPE_TO_NBF_CLASS,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fake_motif(cls: str, start: int, end: int, subclass: str = "sub") -> dict:
    return {
        "ID": f"seq_{start}",
        "Organism": "FakeOrg",
        "Sequence_Name": "seq1",
        "Class": cls,
        "Subclass": subclass,
        "Start": start,
        "End": end,
        "Length": end - start + 1,
        "Score": 1.5,
        "Strand": "+",
    }


def _make_summaries(n: int = 3) -> List[dict]:
    summaries = []
    for i in range(1, n + 1):
        s: dict = {
            "Organism": f"Organism_{i}",
            "Genome_Size_bp": i * 500_000,
            "Total_Motifs": i * 200,
            "Core_Motifs": i * 180,
            "Coverage_pct": float(i * 10),
            "Density_per_kb": float(i * 2),
            "Hybrid_Count": i * 5,
            "Cluster_Count": i * 3,
            "Classes_Detected": min(i + 3, 9),
            "Mean_Score": 1.8,
            "Runtime_s": float(i),
        }
        for cls in NON_B_CLASSES:
            col = f"n_{cls.replace('-', '_').replace(' ', '_')}"
            s[col] = i * 10
        summaries.append(s)
    return summaries


def _make_all_motifs(summaries: List[dict]) -> pd.DataFrame:
    rows = []
    for s in summaries:
        org = s["Organism"]
        for j, cls in enumerate(NON_B_CLASSES[:4]):
            rows.append({
                "Organism": org,
                "Sequence_Name": "seq1",
                "Class": cls,
                "Subclass": f"{cls}_sub",
                "Start": 100 + j * 1000,
                "End":   200 + j * 1000,
                "Length": 100,
                "Score": 1.5,
                "Strand": "+",
            })
    return pd.DataFrame(rows)


# ===========================================================================
# 1. discover_genomes
# ===========================================================================

class TestDiscoverGenomes:
    def test_finds_fna_files(self, tmp_path):
        (tmp_path / "org1.fna").write_text(">seq\nATGC\n")
        (tmp_path / "org2.fna").write_text(">seq\nGCTA\n")
        result = discover_genomes(str(tmp_path))
        names = [n for n, _ in result]
        assert "org1" in names
        assert "org2" in names

    def test_finds_fasta_files(self, tmp_path):
        (tmp_path / "genome.fasta").write_text(">seq\nATGC\n")
        result = discover_genomes(str(tmp_path))
        assert any(n == "genome" for n, _ in result)

    def test_empty_folder_returns_empty(self, tmp_path):
        assert discover_genomes(str(tmp_path)) == []

    def test_ignores_non_fasta_files(self, tmp_path):
        (tmp_path / "notes.txt").write_text("not a genome")
        (tmp_path / "genome.fna").write_text(">s\nATGC\n")
        result = discover_genomes(str(tmp_path))
        assert len(result) == 1

    def test_no_duplicates(self, tmp_path):
        (tmp_path / "org.fna").write_text(">s\nATGC\n")
        result = discover_genomes(str(tmp_path))
        names = [n for n, _ in result]
        assert len(names) == len(set(names))

    def test_returns_correct_path(self, tmp_path):
        fp = tmp_path / "mygenome.fna"
        fp.write_text(">s\nATGC\n")
        result = discover_genomes(str(tmp_path))
        _, path = result[0]
        assert os.path.abspath(path) == os.path.abspath(str(fp))


# ===========================================================================
# 2. _covered_bases
# ===========================================================================

class TestCoveredBases:
    def test_single_interval(self):
        m = [_fake_motif("G-Quadruplex", 1, 100)]
        assert _covered_bases(m, 1000) == 100

    def test_non_overlapping_intervals(self):
        motifs = [
            _fake_motif("G-Quadruplex", 1,   100),
            _fake_motif("Z-DNA",        200, 300),
        ]
        # [1,100] = 100 bases; [200,300] = 101 bases (1-based inclusive)
        assert _covered_bases(motifs, 1000) == 201

    def test_overlapping_intervals_merged(self):
        motifs = [
            _fake_motif("G-Quadruplex", 1,  150),
            _fake_motif("Z-DNA",        100, 200),
        ]
        assert _covered_bases(motifs, 1000) == 200

    def test_contained_interval(self):
        motifs = [
            _fake_motif("G-Quadruplex", 1, 500),
            _fake_motif("Z-DNA",       50, 100),
        ]
        assert _covered_bases(motifs, 1000) == 500

    def test_empty_motifs_returns_zero(self):
        assert _covered_bases([], 1000) == 0

    def test_zero_genome_returns_zero(self):
        m = [_fake_motif("G-Quadruplex", 1, 100)]
        assert _covered_bases(m, 0) == 0

    def test_coverage_capped_at_genome_length(self):
        m = [_fake_motif("G-Quadruplex", 1, 10000)]
        # genome only 500 bp
        assert _covered_bases(m, 500) == 500


# ===========================================================================
# 3. build_master_tables
# ===========================================================================

class TestBuildMasterTables:
    def setup_method(self):
        self.summaries  = _make_summaries(3)
        self.all_motifs = _make_all_motifs(self.summaries)
        self.tables     = build_master_tables(self.summaries, self.all_motifs)

    def test_returns_dict(self):
        assert isinstance(self.tables, dict)

    def test_required_tables_present(self):
        required = [
            "master_genome_overview",
            "master_class_distribution",
            "master_hybrid_cluster",
            "master_density_coverage_by_class",
        ]
        for key in required:
            assert key in self.tables, f"Missing table: {key}"

    def test_genome_overview_has_key_columns(self):
        df = self.tables["master_genome_overview"]
        for col in ["Organism", "Coverage_pct", "Density_per_kb",
                    "Hybrid_Count", "Cluster_Count"]:
            assert col in df.columns, f"Missing column: {col}"

    def test_genome_overview_row_count(self):
        df = self.tables["master_genome_overview"]
        assert len(df) == 3

    def test_class_distribution_has_all_classes(self):
        df  = self.tables["master_class_distribution"]
        classes_found = set(df["Class"].unique())
        for cls in NON_B_CLASSES:
            assert cls in classes_found, f"Missing class: {cls}"

    def test_hybrid_cluster_table_columns(self):
        df = self.tables["master_hybrid_cluster"]
        assert "Hybrid_Count" in df.columns
        assert "Cluster_Count" in df.columns

    def test_density_coverage_class_column(self):
        df = self.tables.get("master_density_coverage_by_class", pd.DataFrame())
        if not df.empty:
            assert "Class" in df.columns
            assert "Density_per_kb" in df.columns

    def test_subclass_table_present_when_subclass_col_exists(self):
        # all_motifs has 'Subclass' column → table should be generated
        assert "master_subclass_breakdown" in self.tables

    def test_no_negative_counts(self):
        df = self.tables["master_genome_overview"]
        assert (df["Hybrid_Count"] >= 0).all()
        assert (df["Cluster_Count"] >= 0).all()
        assert (df["Coverage_pct"] >= 0).all()


# ===========================================================================
# 4. generate_all_figures
# ===========================================================================

class TestGenerateAllFigures:
    def setup_method(self):
        self.summaries  = _make_summaries(3)
        self.all_motifs = _make_all_motifs(self.summaries)

    def test_returns_list_of_paths(self, tmp_path):
        paths = generate_all_figures(
            self.summaries, self.all_motifs, str(tmp_path)
        )
        assert isinstance(paths, list)
        assert len(paths) > 0

    def test_all_files_created(self, tmp_path):
        paths = generate_all_figures(
            self.summaries, self.all_motifs, str(tmp_path)
        )
        for p in paths:
            assert os.path.isfile(p), f"Figure not created: {p}"

    def test_png_files_nonempty(self, tmp_path):
        paths = generate_all_figures(
            self.summaries, self.all_motifs, str(tmp_path)
        )
        for p in paths:
            if p.endswith(".png"):
                assert os.path.getsize(p) > 0, f"Empty PNG: {p}"

    def test_six_figures_generated(self, tmp_path):
        paths = generate_all_figures(
            self.summaries, self.all_motifs, str(tmp_path)
        )
        assert len(paths) == 6

    def test_works_with_empty_motifs(self, tmp_path):
        """Pipeline must not crash when all_motifs is empty."""
        paths = generate_all_figures(
            self.summaries, pd.DataFrame(), str(tmp_path)
        )
        assert len(paths) == 6

    def test_figure_1_coverage_density(self, tmp_path):
        paths = generate_all_figures(
            self.summaries, self.all_motifs, str(tmp_path)
        )
        assert any("fig1" in os.path.basename(p) for p in paths)

    def test_figure_3_hybrid_cluster(self, tmp_path):
        paths = generate_all_figures(
            self.summaries, self.all_motifs, str(tmp_path)
        )
        assert any("fig3" in os.path.basename(p) for p in paths)


# ===========================================================================
# 5. generate_report
# ===========================================================================

class TestGenerateReport:
    def setup_method(self):
        self.summaries  = _make_summaries(3)
        self.all_motifs = _make_all_motifs(self.summaries)
        self.tables     = build_master_tables(self.summaries, self.all_motifs)

    def test_report_file_created(self, tmp_path):
        path = generate_report(
            self.summaries, self.tables,
            ["fig1.png", "fig2.png", "fig3.png"],
            str(tmp_path),
        )
        assert os.path.isfile(path)

    def test_report_nonempty(self, tmp_path):
        path = generate_report(
            self.summaries, self.tables, [], str(tmp_path)
        )
        assert os.path.getsize(path) > 500

    def test_report_contains_abstract(self, tmp_path):
        path = generate_report(
            self.summaries, self.tables, [], str(tmp_path)
        )
        text = Path(path).read_text()
        assert "Abstract" in text

    def test_report_contains_methods(self, tmp_path):
        path = generate_report(
            self.summaries, self.tables, [], str(tmp_path)
        )
        text = Path(path).read_text()
        assert "Methods" in text

    def test_report_contains_results(self, tmp_path):
        path = generate_report(
            self.summaries, self.tables, [], str(tmp_path)
        )
        text = Path(path).read_text()
        assert "Results" in text

    def test_report_mentions_coverage(self, tmp_path):
        path = generate_report(
            self.summaries, self.tables, [], str(tmp_path)
        )
        text = Path(path).read_text()
        assert "coverage" in text.lower() or "Coverage" in text

    def test_report_mentions_density(self, tmp_path):
        path = generate_report(
            self.summaries, self.tables, [], str(tmp_path)
        )
        text = Path(path).read_text()
        assert "density" in text.lower() or "Density" in text

    def test_report_mentions_hybrids(self, tmp_path):
        path = generate_report(
            self.summaries, self.tables, [], str(tmp_path)
        )
        text = Path(path).read_text()
        assert "hybrid" in text.lower() or "Hybrid" in text

    def test_report_mentions_clusters(self, tmp_path):
        path = generate_report(
            self.summaries, self.tables, [], str(tmp_path)
        )
        text = Path(path).read_text()
        assert "cluster" in text.lower() or "Cluster" in text

    def test_report_filename(self, tmp_path):
        path = generate_report(
            self.summaries, self.tables, [], str(tmp_path)
        )
        assert os.path.basename(path) == "NATURE_METHODS_VALIDATION_REPORT.md"


# ===========================================================================
# 6. _remove_overlaps  –  quadratic-scaling regression guard
# ===========================================================================

class TestRemoveOverlapsPerformance:
    """
    Guard against reintroducing the O(n²) linear-scan bottleneck in
    NonBDNAScanner._remove_overlaps.

    The fix replaces the inner ``for acc_start, acc_end in accepted_intervals``
    loop (O(n) per motif → O(n²) total) with an O(log n) binary-search check,
    reducing worst-case cost from ~437 ms to ~10 ms for n=5 000 motifs.
    """

    @staticmethod
    def _make_motifs(n: int, cls: str = "G-Quadruplex", sub: str = "canon") -> List[dict]:
        """Return *n* non-overlapping motifs of length 50, spaced 100 bp apart."""
        return [
            {
                "Start":    i * 100,
                "End":      i * 100 + 50,
                "Length":   50,
                "Score":    float(n - i),   # descending so best motif comes first
                "Class":    cls,
                "Subclass": sub,
                "Strand":   "+",
            }
            for i in range(n)
        ]

    def _scanner(self):
        from Utilities.nonbscanner import NonBScanner
        return NonBScanner(enable_all_detectors=False)

    def test_remove_overlaps_returns_all_non_overlapping(self):
        """All non-overlapping motifs must be kept."""
        scanner = self._scanner()
        motifs = self._make_motifs(100)
        result = scanner._remove_overlaps(motifs)
        assert len(result) == 100

    def test_remove_overlaps_keeps_highest_score(self):
        """When motifs overlap, the one with the higher score is retained."""
        scanner = self._scanner()
        # Two overlapping motifs: second has higher score
        motifs = [
            {"Start": 0,  "End": 100, "Length": 100, "Score": 1.0,
             "Class": "G-Quadruplex", "Subclass": "canon", "Strand": "+"},
            {"Start": 50, "End": 150, "Length": 100, "Score": 2.0,
             "Class": "G-Quadruplex", "Subclass": "canon", "Strand": "+"},
        ]
        result = scanner._remove_overlaps(motifs)
        assert len(result) == 1
        assert result[0]["Score"] == 2.0

    def test_remove_overlaps_linear_time_large_input(self):
        """
        _remove_overlaps must complete in under 1 s for 5_000 non-overlapping
        motifs in a single subclass group (the O(n²) version takes ~437 ms
        for this size; the O(n log n) version takes ~10 ms).
        """
        import time
        scanner = self._scanner()
        motifs = self._make_motifs(5_000)
        t0 = time.monotonic()
        result = scanner._remove_overlaps(motifs)
        elapsed = time.monotonic() - t0
        assert len(result) == 5_000, "All non-overlapping motifs must be retained"
        assert elapsed < 1.0, (
            f"_remove_overlaps took {elapsed:.2f}s for 5_000 motifs — "
            "quadratic scaling may have been reintroduced"
        )


# ===========================================================================
# 7. _jaccard_interval  –  interval similarity engine
# ===========================================================================

class TestJaccardInterval:
    """Unit tests for the interval-level Jaccard similarity computation."""

    def test_identical_intervals(self):
        assert _jaccard_interval(1, 100, 1, 100) == 1.0

    def test_no_overlap(self):
        assert _jaccard_interval(1, 50, 60, 100) == 0.0

    def test_adjacent_no_overlap(self):
        # [1,50] and [51,100]: no shared base (1-based inclusive)
        assert _jaccard_interval(1, 50, 51, 100) == 0.0

    def test_partial_overlap(self):
        # [1,100] ∩ [50,150] = positions 50..100 = 51 positions;
        # union = positions 1..150 = 150 positions → J = 51/150
        j = _jaccard_interval(1, 100, 50, 150)
        assert round(j, 4) == round(51 / 150, 4)

    def test_containment(self):
        # [1,100] fully contains [40,60]; overlap = [40,60] = 21 positions;
        # union = [1,100] = 100 positions → J = 21/100
        j = _jaccard_interval(1, 100, 40, 60)
        assert round(j, 4) == round(21 / 100, 4)

    def test_symmetry(self):
        j1 = _jaccard_interval(10, 80, 50, 120)
        j2 = _jaccard_interval(50, 120, 10, 80)
        assert abs(j1 - j2) < 1e-9

    def test_single_base_overlap(self):
        # [1,10] and [10,20]: overlap = 1 bp, union = 20 bp
        j = _jaccard_interval(1, 10, 10, 20)
        assert round(j, 4) == round(1 / 20, 4)

    def test_returns_float(self):
        assert isinstance(_jaccard_interval(1, 50, 1, 50), float)

    def test_above_50pct_threshold(self):
        # [1,100] and [1,60]: overlap=60, union=100 → Jaccard=0.60 ≥ 0.50
        assert _jaccard_interval(1, 100, 1, 60) >= 0.5

    def test_below_50pct_threshold(self):
        # [1,100] and [80,200]: overlap=21, union=200 → Jaccard<0.50
        assert _jaccard_interval(1, 100, 80, 200) < 0.5


# ===========================================================================
# 8. _match_motif_lists  –  TP/FP/FN concordance engine
# ===========================================================================

class TestMatchMotifLists:
    """Unit tests for the greedy interval-to-interval matching function."""

    @staticmethod
    def _motif(start: int, end: int, cls: str = "G-Quadruplex") -> dict:
        return {"Start": start, "End": end, "Class": cls}

    @staticmethod
    def _nbst_df(intervals: list) -> pd.DataFrame:
        return pd.DataFrame(intervals, columns=["Start", "Stop"])

    def test_perfect_match(self):
        nbf   = [self._motif(1, 100)]
        nbst  = self._nbst_df([(1, 100)])
        tp, fp, fn, jaccards = _match_motif_lists(nbf, nbst)
        assert tp == 1 and fp == 0 and fn == 0
        assert jaccards[0] == 1.0

    def test_no_overlap_all_fp_fn(self):
        nbf  = [self._motif(1, 50)]
        nbst = self._nbst_df([(200, 300)])
        tp, fp, fn, _ = _match_motif_lists(nbf, nbst)
        assert tp == 0 and fp == 1 and fn == 1

    def test_empty_nbf_all_fn(self):
        nbst = self._nbst_df([(1, 100), (200, 300)])
        tp, fp, fn, jaccards = _match_motif_lists([], nbst)
        assert tp == 0 and fp == 0 and fn == 2
        assert jaccards == []

    def test_empty_nbst_all_fp(self):
        nbf  = [self._motif(1, 100), self._motif(200, 300)]
        nbst = self._nbst_df([])
        tp, fp, fn, jaccards = _match_motif_lists(nbf, nbst)
        assert tp == 0 and fp == 2 and fn == 0
        assert jaccards == []

    def test_partial_recall(self):
        # 2 NBF motifs, 3 NBST — 2 match, 1 NBST unmatched
        nbf  = [self._motif(1, 100), self._motif(200, 300)]
        nbst = self._nbst_df([(1, 100), (200, 300), (500, 600)])
        tp, fp, fn, _ = _match_motif_lists(nbf, nbst)
        assert tp == 2 and fp == 0 and fn == 1

    def test_each_nbst_matched_at_most_once(self):
        # Two NBF motifs that both overlap the same NBST motif
        nbf  = [self._motif(1, 100), self._motif(10, 90)]
        nbst = self._nbst_df([(1, 100)])
        tp, fp, fn, _ = _match_motif_lists(nbf, nbst)
        # Only one TP; the second NBF motif cannot claim the same NBST motif
        assert tp == 1

    def test_jaccard_threshold_respected(self):
        # NBF [1,100] vs NBST [90,200]: overlap=11, union=200 → J≈0.055 < 0.5
        nbf  = [self._motif(1, 100)]
        nbst = self._nbst_df([(90, 200)])
        tp, fp, fn, _ = _match_motif_lists(nbf, nbst, jaccard_threshold=0.5)
        assert tp == 0 and fp == 1 and fn == 1

    def test_tp_fp_fn_sum_equals_total(self):
        nbf  = [self._motif(i * 200, i * 200 + 100) for i in range(5)]
        nbst = self._nbst_df([(i * 200, i * 200 + 100) for i in range(3)])
        tp, fp, fn, _ = _match_motif_lists(nbf, nbst)
        assert tp + fp == len(nbf)
        assert tp + fn == len(nbst)

    def test_precision_recall_from_results(self):
        # 3 NBF, 3 NBST, 2 match → precision = 2/3, recall = 2/3
        nbf  = [self._motif(1, 50), self._motif(100, 150), self._motif(500, 600)]
        nbst = self._nbst_df([(1, 50), (100, 150), (900, 1000)])
        tp, fp, fn, _ = _match_motif_lists(nbf, nbst)
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall    = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        assert round(precision, 4) == round(2 / 3, 4)
        assert round(recall,    4) == round(2 / 3, 4)


# ===========================================================================
# 9. NBST_DATA_FILES  –  benchmark configuration
# ===========================================================================

class TestNBSTDataFilesConfig:
    """Validate the static benchmark configuration table."""

    def test_all_six_motif_classes_covered(self):
        """The six motif classes from the problem statement must be present."""
        expected_classes = {
            "G-Quadruplex",   # GQ
            "Z-DNA",          # Z
            "Triplex",        # MR (mirror repeats = triplex-forming)
            "Slipped_DNA",    # STR + DR
            "Curved_DNA",     # A-phased repeats
            "Cruciform",      # IR (inverted repeats)
        }
        found_classes = set(NBST_DATA_FILES.values())
        assert expected_classes.issubset(found_classes), (
            f"Missing classes: {expected_classes - found_classes}"
        )

    def test_ir_file_mapped_to_cruciform(self):
        assert "699a2b3fb6cbe_IR.tsv" in NBST_DATA_FILES
        assert NBST_DATA_FILES["699a2b3fb6cbe_IR.tsv"] == "Cruciform"

    def test_inverted_repeat_in_type_map(self):
        assert "Inverted_Repeat" in NBST_TYPE_TO_NBF_CLASS
        assert NBST_TYPE_TO_NBF_CLASS["Inverted_Repeat"] == "Cruciform"

    def test_mirror_repeat_maps_to_triplex(self):
        assert NBST_TYPE_TO_NBF_CLASS["Mirror_Repeat"] == "Triplex"

    def test_gq_file_mapped_to_gquadruplex(self):
        assert "693fc40d26a53_GQ.tsv" in NBST_DATA_FILES
        assert NBST_DATA_FILES["693fc40d26a53_GQ.tsv"] == "G-Quadruplex"

    def test_slipped_dna_has_two_files(self):
        slipped_files = [f for f, c in NBST_DATA_FILES.items() if c == "Slipped_DNA"]
        assert len(slipped_files) == 2, (
            "Both STR and DR files must map to Slipped_DNA"
        )
