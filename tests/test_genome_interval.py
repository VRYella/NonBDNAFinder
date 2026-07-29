"""
Unit tests for Utilities.genome_interval
-----------------------------------------
Covers:
  - GenomicInterval construction and validation
  - parse_interval_string (compact format)
  - build_interval (structured inputs)
  - validate_accession
  - validate_interval_inputs
  - absolute_position and annotate_motif coordinate calculations
"""

import pytest
import sys
import os

# Ensure the package root is on the path when running pytest from any directory
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from Utilities.genome_interval import (
    GenomicInterval,
    build_interval,
    parse_interval_string,
    validate_accession,
    validate_interval_inputs,
)


# ─── GenomicInterval construction ────────────────────────────────────────────

class TestGenomicIntervalConstruction:

    def test_basic_construction(self):
        gi = GenomicInterval(accession="NC_000913.3", start=100000, end=150000)
        assert gi.accession == "NC_000913.3"
        assert gi.start == 100000
        assert gi.end == 150000

    def test_length(self):
        gi = GenomicInterval(accession="NC_000913.3", start=100000, end=150000)
        assert gi.length == 50001

    def test_length_single_base(self):
        gi = GenomicInterval(accession="NC_000913.3", start=5, end=5)
        assert gi.length == 1

    def test_compact_property(self):
        gi = GenomicInterval(accession="NC_000913.3", start=100000, end=150000)
        assert gi.compact == "NC_000913.3:100000-150000"

    def test_start_below_one_raises(self):
        with pytest.raises(ValueError, match="start must be"):
            GenomicInterval(accession="NC_000913.3", start=0, end=100)

    def test_end_before_start_raises(self):
        with pytest.raises(ValueError, match="end.*must be"):
            GenomicInterval(accession="NC_000913.3", start=200, end=100)

    def test_optional_fields_default(self):
        gi = GenomicInterval(accession="NC_000913.3", start=1, end=100)
        assert gi.organism == ""
        assert gi.chromosome == ""
        assert gi.provider == "ncbi"

    def test_optional_fields_set(self):
        gi = GenomicInterval(
            accession="NC_000913.3",
            start=1,
            end=100,
            organism="Escherichia coli",
            chromosome="chr1",
        )
        assert gi.organism == "Escherichia coli"
        assert gi.chromosome == "chr1"


# ─── Coordinate conversion ────────────────────────────────────────────────────

class TestAbsolutePosition:

    def setup_method(self):
        # Interval starting at genomic position 100000
        self.gi = GenomicInterval(accession="NC_000913.3", start=100000, end=150000)

    def test_zero_based_first_base(self):
        # relative position 0 → absolute 100000
        assert self.gi.absolute_position(0, zero_based=True) == 100000

    def test_zero_based_arbitrary(self):
        # relative 500 → absolute 100500
        assert self.gi.absolute_position(500, zero_based=True) == 100500

    def test_one_based_first_base(self):
        # relative 1 → absolute 100000
        assert self.gi.absolute_position(1, zero_based=False) == 100000

    def test_one_based_arbitrary(self):
        # relative 501 → absolute 100500
        assert self.gi.absolute_position(501, zero_based=False) == 100500


# ─── annotate_motif ───────────────────────────────────────────────────────────

class TestAnnotateMotif:

    def setup_method(self):
        self.gi = GenomicInterval(
            accession="NC_000913.3",
            start=100000,
            end=150000,
            organism="Escherichia coli",
            chromosome="chr1",
        )

    def test_coordinate_fields_added(self):
        motif = {"Start": 500, "End": 530, "Motif": "GGGG"}
        result = self.gi.annotate_motif(motif)
        assert result["Accession"] == "NC_000913.3"
        assert result["Organism"] == "Escherichia coli"
        assert result["Chromosome"] == "chr1"
        assert result["Interval_Start"] == 100000
        assert result["Interval_End"] == 150000
        assert result["Relative_Start"] == 500
        assert result["Relative_End"] == 530

    def test_absolute_coordinates(self):
        motif = {"Start": 1, "End": 10}
        result = self.gi.annotate_motif(motif)
        # 1-based relative → absolute: interval_start + relative - 1
        assert result["Absolute_Start"] == 100000
        assert result["Absolute_End"] == 100009

    def test_returns_same_dict(self):
        motif = {"Start": 1, "End": 5}
        result = self.gi.annotate_motif(motif)
        assert result is motif  # modified in-place and returned


# ─── parse_interval_string ───────────────────────────────────────────────────

class TestParseIntervalString:

    def test_basic_parse(self):
        gi = parse_interval_string("NC_000913.3:100000-150000")
        assert gi.accession == "NC_000913.3"
        assert gi.start == 100000
        assert gi.end == 150000

    def test_no_version_accession(self):
        gi = parse_interval_string("NC_000913:100000-150000")
        assert gi.accession == "NC_000913"

    def test_whitespace_stripped(self):
        gi = parse_interval_string("  NC_000913.3:100000-150000  ")
        assert gi.accession == "NC_000913.3"

    def test_single_base_interval(self):
        gi = parse_interval_string("NC_000913.3:500-500")
        assert gi.start == 500
        assert gi.end == 500
        assert gi.length == 1

    def test_invalid_format_raises(self):
        with pytest.raises(ValueError, match="Cannot parse"):
            parse_interval_string("NC_000913.3 100000 150000")

    def test_missing_colon_raises(self):
        with pytest.raises(ValueError):
            parse_interval_string("NC_000913.3-100000-150000")

    def test_empty_string_raises(self):
        with pytest.raises(ValueError):
            parse_interval_string("")

    def test_none_raises(self):
        with pytest.raises(ValueError):
            parse_interval_string(None)  # type: ignore[arg-type]

    def test_end_before_start_raises(self):
        """parse_interval_string defers coordinate validation to GenomicInterval."""
        with pytest.raises(ValueError):
            parse_interval_string("NC_000913.3:150000-100000")

    def test_genbank_style_accession(self):
        gi = parse_interval_string("AB012345.1:1000-2000")
        assert gi.accession == "AB012345.1"
        assert gi.start == 1000
        assert gi.end == 2000


# ─── build_interval ──────────────────────────────────────────────────────────

class TestBuildInterval:

    def test_basic(self):
        gi = build_interval("NC_000913.3", 100000, 150000)
        assert gi.accession == "NC_000913.3"
        assert gi.start == 100000
        assert gi.end == 150000

    def test_organism_and_chromosome(self):
        gi = build_interval(
            "NC_000913.3", 1, 100,
            organism="E. coli", chromosome="chr1"
        )
        assert gi.organism == "E. coli"
        assert gi.chromosome == "chr1"

    def test_invalid_start_raises(self):
        with pytest.raises(ValueError):
            build_interval("NC_000913.3", 0, 100)

    def test_invalid_end_raises(self):
        with pytest.raises(ValueError):
            build_interval("NC_000913.3", 100, 50)


# ─── validate_accession ──────────────────────────────────────────────────────

class TestValidateAccession:

    @pytest.mark.parametrize("acc", [
        "NC_000913.3",
        "NM_001234.5",
        "NR_003287.2",
        "NG_012345.1",
        "NW_000107077.1",
        "NZ_CP000253.1",   # WGS master record (alphanumeric body)
        "XM_001234567.1",
        "XR_001234567.1",
        "AB012345.1",
        "AY123456",
        "AB012345",
    ])
    def test_valid_accessions(self, acc):
        assert validate_accession(acc) is True, f"Expected {acc!r} to be valid"

    @pytest.mark.parametrize("acc", [
        "",
        None,
        "hello",
        "12345",
        "NOTANACCESSION",
        "NC_",
        # validate_accession strips whitespace before matching so leading spaces
        # are accepted (matches NCBI's own handling); we test explicit rejects only
    ])
    def test_invalid_accessions(self, acc):
        assert validate_accession(acc) is False, f"Expected {acc!r} to be invalid"

    def test_leading_space_is_stripped_and_valid(self):
        """validate_accession strips leading/trailing whitespace before matching."""
        assert validate_accession(" NC_000913.3") is True


# ─── validate_interval_inputs ────────────────────────────────────────────────

class TestValidateIntervalInputs:

    def test_valid_inputs_no_errors(self):
        errors = validate_interval_inputs("NC_000913.3", 100000, 150000)
        assert errors == []

    def test_empty_accession(self):
        errors = validate_interval_inputs("", 100000, 150000)
        assert any("Accession is required" in e for e in errors)

    def test_invalid_accession_format(self):
        errors = validate_interval_inputs("NOTVALID", 100000, 150000)
        assert any("does not look like a valid" in e for e in errors)

    def test_missing_start(self):
        errors = validate_interval_inputs("NC_000913.3", None, 150000)
        assert any("Start coordinate is required" in e for e in errors)

    def test_start_below_one(self):
        errors = validate_interval_inputs("NC_000913.3", 0, 150000)
        assert any("must be ≥ 1" in e for e in errors)

    def test_missing_end(self):
        errors = validate_interval_inputs("NC_000913.3", 100000, None)
        assert any("End coordinate is required" in e for e in errors)

    def test_end_before_start(self):
        errors = validate_interval_inputs("NC_000913.3", 150000, 100000)
        assert any("End coordinate" in e and "≥ start" in e for e in errors)

    def test_interval_too_large(self):
        errors = validate_interval_inputs("NC_000913.3", 1, 60_000_000, max_length=50_000_000)
        assert any("larger than the maximum allowed" in e for e in errors)

    def test_multiple_errors_reported(self):
        errors = validate_interval_inputs("", None, None)
        assert len(errors) >= 2

    def test_equal_start_end_valid(self):
        errors = validate_interval_inputs("NC_000913.3", 500, 500)
        assert errors == []
