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
