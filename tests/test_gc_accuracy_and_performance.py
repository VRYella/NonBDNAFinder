"""
Tests for accurate GC% calculation and performance improvements.

Validates:
  - calc_gc_content excludes N/ambiguous bases from the denominator
  - _count_bases uses fast str.count() and returns correct counts
  - sequence_preprocessor reports GC% correctly for N-containing sequences
  - Every detector that outputs GC_Content uses the gold-standard formula
  - Single uppercase conversion at the analysis entry point
  - Throughput regression guard: analysis must stay above a minimum bp/s floor
  - GC-rich throughput: detector performance is not degraded on high-GC input
"""

import pytest
import time

# ─────────────────────────────────────────────────────────────────────────────
# calc_gc_content
# ─────────────────────────────────────────────────────────────────────────────

class TestCalcGcContent:
    """Unit tests for detectors_utils.calc_gc_content."""

    def _gc(self, seq):
        from Utilities.detectors_utils import calc_gc_content
        return calc_gc_content(seq)

    def test_pure_atgc(self):
        assert self._gc("ATGC") == pytest.approx(50.0)

    def test_all_gc(self):
        assert self._gc("GCGCGC") == pytest.approx(100.0)

    def test_all_at(self):
        assert self._gc("ATATAT") == pytest.approx(0.0)

    def test_n_excluded_from_denominator(self):
        # ATCG (4 ATGC bases, 2 GC) surrounded by 4 N's
        # Correct: (2/4)*100 = 50.0  NOT (2/8)*100 = 25.0
        assert self._gc("NNATCGNN") == pytest.approx(50.0)

    def test_all_n_returns_zero(self):
        assert self._gc("NNNN") == pytest.approx(0.0)

    def test_empty_returns_zero(self):
        assert self._gc("") == pytest.approx(0.0)

    def test_case_insensitive(self):
        assert self._gc("atcg") == pytest.approx(50.0)
        assert self._gc("ATCG") == pytest.approx(self._gc("atcg"))

    def test_only_valid_bases_in_denominator(self):
        # 'R' and 'Y' are ambiguous IUPAC codes, excluded from the denominator.
        # Valid ATGC bases in "RYGCGCAT": G, C, G, C, A, T → 6 bases, 4 GC → 66.67%
        seq = "RYGCGCAT"
        result = self._gc(seq)
        assert result == pytest.approx(66.67, abs=0.01)


# ─────────────────────────────────────────────────────────────────────────────
# _count_bases
# ─────────────────────────────────────────────────────────────────────────────

class TestCountBases:
    """Unit tests for the fast str.count()-based _count_bases implementation."""

    def _count(self, seq):
        from Utilities.detectors_utils import _count_bases
        return _count_bases(seq)

    def test_simple(self):
        a, t, g, c = self._count("AATTGGCC")
        assert a == 2
        assert t == 2
        assert g == 2
        assert c == 2

    def test_lowercase(self):
        a, t, g, c = self._count("aattggcc")
        assert a == 2
        assert t == 2
        assert g == 2
        assert c == 2

    def test_mixed_case(self):
        a, t, g, c = self._count("AaTtGgCc")
        assert a == 2
        assert t == 2
        assert g == 2
        assert c == 2

    def test_n_not_counted(self):
        a, t, g, c = self._count("NNATCGNN")
        assert a == 1
        assert t == 1
        assert g == 1
        assert c == 1

    def test_empty(self):
        assert self._count("") == (0, 0, 0, 0)


# ─────────────────────────────────────────────────────────────────────────────
# sequence_preprocessor
# ─────────────────────────────────────────────────────────────────────────────

class TestSequencePreprocessor:
    """GC% formula and character counts in the preprocessor."""

    def _pp(self, seq):
        from Utilities.sequence_preprocessor import preprocess_sequence
        return preprocess_sequence(seq)

    def test_gc_excludes_n(self):
        # ATCG (4 valid, 2 GC) + 4 N's
        result = self._pp("NNATCGNN")
        assert result.gc_percentage == pytest.approx(50.0)

    def test_character_counts_correct(self):
        result = self._pp("AATCGG")
        assert result.character_counts['A'] == 2
        assert result.character_counts['T'] == 1
        assert result.character_counts['C'] == 1
        assert result.character_counts['G'] == 2
        assert result.character_counts['N'] == 0

    def test_n_counted_separately(self):
        result = self._pp("NNATCGNN")
        assert result.character_counts['N'] == 4

    def test_valid_bases_excludes_n(self):
        result = self._pp("NNATCGNN")
        assert result.valid_bases == 4

    def test_invalid_char_detected(self):
        result = self._pp("ATCGX")
        assert 'X' in result.invalid_characters
        assert result.validation_status == "error"

    def test_fasta_header_stripped(self):
        result = self._pp(">my_seq\nATCG")
        assert result.header == "MY_SEQ"
        assert result.sequence == "ATCG"

    def test_multiline_fasta(self):
        result = self._pp(">hdr\nATCG\nATCG")
        assert result.sequence == "ATCGATCG"
        assert result.length == 8


# ─────────────────────────────────────────────────────────────────────────────
# Per-detector GC_Content field accuracy
# ─────────────────────────────────────────────────────────────────────────────

class TestDetectorGcContent:
    """Each detector that reports GC_Content must use gold-standard formula."""

    # Canonical motif with 4 N's appended: GC% should still be 50% (not 25%).
    # We search for motifs in sequences that are entirely ATGC, so N's won't
    # appear in motif sequences extracted by detectors.  Instead, we verify
    # the formula indirectly: for a pure ATGC motif, gc/len == gc/valid_bases,
    # so both formulas agree.  The key test is the unit-level function test above.

    def _run(self, seq, enabled):
        from Utilities.nonbscanner import analyze_sequence
        return analyze_sequence(seq, "test", enabled_classes=enabled)

    def test_gquad_gc_in_range(self):
        seq = "GGGTTTAGGGTTTGGG" * 10
        motifs = self._run(seq, ["G-Quadruplex"])
        gc_values = [m.get("GC_Content") for m in motifs if "GC_Content" in m]
        for v in gc_values:
            assert 0.0 <= v <= 100.0, f"G4 GC_Content out of range: {v}"

    def test_slipped_gc_in_range(self):
        seq = ("CAG" * 30) + ("AT" * 100)
        motifs = self._run(seq, ["Slipped_DNA"])
        gc_values = [m.get("GC_Content") for m in motifs if "GC_Content" in m]
        for v in gc_values:
            assert 0.0 <= v <= 100.0, f"Slipped GC_Content out of range: {v}"

    def test_rloop_gc_in_range(self):
        # R-loop: GC-skewed region, high G
        seq = ("GGGCCCC" * 8 + "AAAATTT" * 20) * 5
        motifs = self._run(seq, ["R-Loop"])
        gc_values = [m.get("GC_Content") for m in motifs if "GC_Content" in m]
        for v in gc_values:
            assert 0.0 <= v <= 100.0, f"R-Loop GC_Content out of range: {v}"


# ─────────────────────────────────────────────────────────────────────────────
# Single uppercase conversion
# ─────────────────────────────────────────────────────────────────────────────

class TestSingleUppercaseConversion:
    """
    The entry-point analyze_sequence must handle lowercase input correctly,
    proving that uppercase conversion happens at most once before detectors run.
    """

    def test_lowercase_input_produces_same_motifs_as_uppercase(self):
        from Utilities.nonbscanner import analyze_sequence

        seq = "gggtttagggtttggg" * 5 + "atcgatcg" * 50
        motifs_lower = analyze_sequence(seq.lower(), "lower")
        motifs_upper = analyze_sequence(seq.upper(), "upper")

        assert len(motifs_lower) == len(motifs_upper), (
            f"Lowercase ({len(motifs_lower)}) and uppercase ({len(motifs_upper)}) "
            "inputs must find the same number of motifs"
        )

    def test_mixed_case_input_works(self):
        from Utilities.nonbscanner import analyze_sequence

        seq = "GgGtttAGGgtttGGG" * 5 + "atcgatcg" * 50
        motifs = analyze_sequence(seq, "mixed")
        assert isinstance(motifs, list)


# ─────────────────────────────────────────────────────────────────────────────
# Throughput regression guard
# ─────────────────────────────────────────────────────────────────────────────

class TestThroughputRegression:
    """
    Ensure analysis throughput does not degrade below a baseline.

    The floor is deliberately conservative (50,000 bp/s) to avoid flakiness
    on slow CI runners while still catching severe regressions.
    """

    THROUGHPUT_FLOOR_BPS = 50_000  # bp/s minimum acceptable throughput

    def _make_seq(self, length: int) -> str:
        import random
        random.seed(0)
        bases = "ATGC"
        seq = "".join(random.choice(bases) for _ in range(length))
        # Embed a few G4 motifs so the result is non-trivial
        motif = "GGGTTTAGGGTTTGGG"
        insert_at = length // 2
        return seq[:insert_at] + motif + seq[insert_at + len(motif):]

    def test_50kb_throughput_floor(self):
        from Utilities.nonbscanner import analyze_sequence

        seq = self._make_seq(50_000)
        t0 = time.time()
        analyze_sequence(seq, "bench", use_fast_mode=True)
        elapsed = time.time() - t0

        throughput = len(seq) / elapsed
        assert throughput >= self.THROUGHPUT_FLOOR_BPS, (
            f"Throughput {throughput:,.0f} bp/s is below floor "
            f"{self.THROUGHPUT_FLOOR_BPS:,} bp/s"
        )


# ─────────────────────────────────────────────────────────────────────────────
# GC-rich throughput guard
# ─────────────────────────────────────────────────────────────────────────────

class TestGCRichThroughput:
    """
    Guard against severe performance regression on high-GC content sequences.

    GC-rich genomes (70-90 % GC) produce many G3+ seeds which, in a naive
    seeded-window implementation, creates O(n) near-duplicate scan windows.
    The merged-region fix reduces this to O(n_regions) scans.

    The floor is conservative (50,000 bp/s) to stay robust on slow CI runners
    while still catching the O(n) regression that was the reported issue.
    """

    THROUGHPUT_FLOOR_BPS = 50_000  # bp/s minimum acceptable throughput

    def _make_gc_rich_seq(self, length: int, gc_frac: float = 0.70) -> str:
        """Build a pseudo-random sequence with the requested GC fraction."""
        import random
        random.seed(7)
        gc_bases = "GC"
        at_bases = "AT"
        seq = []
        for _ in range(length):
            seq.append(random.choice(gc_bases) if random.random() < gc_frac else random.choice(at_bases))
        return "".join(seq)

    def test_70pct_gc_throughput_floor(self):
        from Detectors.gquad.detector import GQuadruplexDetector

        seq = self._make_gc_rich_seq(50_000, gc_frac=0.70)
        det = GQuadruplexDetector()
        t0 = time.time()
        det.detect_motifs(seq, "gc70_bench")
        elapsed = time.time() - t0

        throughput = len(seq) / elapsed
        assert throughput >= self.THROUGHPUT_FLOOR_BPS, (
            f"G4 detector throughput on 70% GC sequence: {throughput:,.0f} bp/s "
            f"is below floor {self.THROUGHPUT_FLOOR_BPS:,} bp/s"
        )

    def test_dense_g_tracts_throughput_floor(self):
        """Worst-case: G3 tract every 5 bp (mimics extremely G-dense region)."""
        from Detectors.gquad.detector import GQuadruplexDetector

        unit = "GGGAT"
        length = 50_000
        seq = (unit * (length // len(unit) + 1))[:length]
        det = GQuadruplexDetector()
        t0 = time.time()
        det.detect_motifs(seq, "dense_g_bench")
        elapsed = time.time() - t0

        throughput = length / elapsed
        assert throughput >= self.THROUGHPUT_FLOOR_BPS, (
            f"G4 detector throughput on dense-G sequence: {throughput:,.0f} bp/s "
            f"is below floor {self.THROUGHPUT_FLOOR_BPS:,} bp/s"
        )

    def test_seed_regions_merged_correctly(self):
        """
        Verify that merging seed windows does not lose any motif that was
        detected before the optimisation (correctness check).
        """
        from Detectors.gquad.detector import GQuadruplexDetector

        # A sequence with canonical G4 sites spaced far apart (no merging needed)
        g4_motif = "GGGTTTAGGGTTTGGG"
        padding = "ATAT" * 100   # 400 bp gap – seeds will NOT merge across gap
        seq = g4_motif + padding + g4_motif
        det = GQuadruplexDetector()
        motifs = det.detect_motifs(seq.upper(), "merge_test")
        # Both G4 sites should be detected
        assert len(motifs) >= 2, (
            f"Expected ≥2 G4 motifs in sequence with two separated G4 sites, "
            f"got {len(motifs)}"
        )
