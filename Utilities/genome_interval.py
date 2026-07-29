"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Genome Interval Module - Genomic Coordinate Parsing and NCBI Interval Fetch  │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2026.1            │
└──────────────────────────────────────────────────────────────────────────────┘

This module provides a unified genomic coordinate system for the Non-B DNA
Finder application.  It supports three input modes:

    1. Accession  — NR_003287.2, NC_000913.3, etc.
    2. Gene       — gene symbol or Entrez Gene ID
    3. Genome Interval — NC_000913.3:100000-150000
                         or structured (accession + start + end)

Coordinate conventions
----------------------
All coordinates follow the standard 1-based, fully-closed convention
(same as GFF3 / NCBI GenBank):

    - start ≥ 1
    - end ≥ start
    - the interval [start, end] includes both endpoints

When a subsequence is fetched from NCBI the returned sequence starts at
*absolute_start* so relative coordinates (0-based internally in motif
detectors) map to absolute coordinates as:

    absolute_position = absolute_start + relative_position - 1   (1-based)

Provider architecture
---------------------
:class:`IntervalProvider` is an abstract base class.  Add new providers
(Ensembl, UCSC, local genome) by subclassing it and registering the
class in :data:`INTERVAL_PROVIDERS`.
"""

from __future__ import annotations

import logging
import re
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ─── Accession validation ─────────────────────────────────────────────────────

# RefSeq-style accession patterns (covers most common use-cases)
_ACCESSION_PATTERN = re.compile(
    r"""^
    (?:
        # RefSeq with underscore: XY_<alphanumeric> (NM_, NC_, NZ_CP..., NW_, etc.)
        # The body after _ may be digits-only or mixed alphanumeric (WGS records)
        [A-Z]{1,3}_[A-Z0-9]{5,15}(?:\.\d+)?
        |
        # GenBank/EMBL/DDBJ — 1–2 letters + 5–8 digits (+ optional version)
        # Require the letter prefix to be followed immediately by digits (no underscore)
        [A-Z]{1,2}\d{5,8}(?:\.\d+)?
    )$""",
    re.VERBOSE | re.IGNORECASE,
)

# Compact interval string: ACCESSION:START-END
_INTERVAL_COMPACT_RE = re.compile(
    r"""^
    (?P<accession>[A-Za-z0-9_\.]+?)  # accession (non-greedy)
    :                                 # colon separator
    (?P<start>\d+)                    # 1-based start
    -                                 # dash
    (?P<end>\d+)                      # 1-based end
    $""",
    re.VERBOSE,
)

# ─── Data structures ──────────────────────────────────────────────────────────


@dataclass
class GenomicInterval:
    """Fully-annotated genomic locus.

    Coordinates are 1-based, fully-closed [start, end] matching NCBI/GFF3
    convention.  The *organism* and *chromosome* fields are informational and
    may be empty when not known ahead of time.
    """

    accession: str
    """RefSeq / GenBank accession, e.g. NC_000913.3."""

    start: int
    """1-based interval start (inclusive)."""

    end: int
    """1-based interval end (inclusive)."""

    organism: str = ""
    """Organism name or taxonomy, e.g. 'Homo sapiens'."""

    chromosome: str = ""
    """Chromosome or contig label, e.g. 'chr1' or '1'."""

    provider: str = "ncbi"
    """Data provider identifier ('ncbi', 'ensembl', 'ucsc', 'local')."""

    def __post_init__(self) -> None:
        if self.start < 1:
            raise ValueError(
                f"Interval start must be ≥ 1 (got {self.start})"
            )
        if self.end < self.start:
            raise ValueError(
                f"Interval end ({self.end}) must be ≥ start ({self.start})"
            )

    # ── Convenience helpers ──────────────────────────────────────────────────

    @property
    def length(self) -> int:
        """Length in base-pairs (inclusive of both endpoints)."""
        return self.end - self.start + 1

    @property
    def compact(self) -> str:
        """Return the compact representation: ``accession:start-end``."""
        return f"{self.accession}:{self.start}-{self.end}"

    def absolute_position(self, relative_pos: int, zero_based: bool = True) -> int:
        """Convert a position relative to this interval to an absolute genome position.

        Parameters
        ----------
        relative_pos:
            Position measured from the interval start.
            If *zero_based* is True this is 0-indexed (the first base = 0).
        zero_based:
            Whether *relative_pos* is 0-indexed (True) or 1-indexed (False).

        Returns
        -------
        int
            Absolute 1-based genomic position.
        """
        if zero_based:
            return self.start + relative_pos
        return self.start + relative_pos - 1

    def annotate_motif(self, motif: dict) -> dict:
        """Inject absolute coordinate fields into a motif result dictionary.

        The motif dict is expected to contain ``Start`` and ``End`` keys
        that hold positions relative to the fetched subsequence (1-based).
        The following keys are added / updated in-place and also returned:

        - ``Organism``
        - ``Accession``
        - ``Chromosome``
        - ``Interval_Start``
        - ``Interval_End``
        - ``Relative_Start``   — original Start value (1-based within interval)
        - ``Relative_End``     — original End value (1-based within interval)
        - ``Absolute_Start``   — 1-based absolute genome position of motif start
        - ``Absolute_End``     — 1-based absolute genome position of motif end
        """
        rel_start = motif.get("Start", 0)
        rel_end = motif.get("End", 0)
        motif["Organism"] = self.organism
        motif["Accession"] = self.accession
        motif["Chromosome"] = self.chromosome
        motif["Interval_Start"] = self.start
        motif["Interval_End"] = self.end
        motif["Relative_Start"] = rel_start
        motif["Relative_End"] = rel_end
        # Positions stored in motif dicts are 0-based from the detector
        motif["Absolute_Start"] = self.absolute_position(rel_start, zero_based=False)
        motif["Absolute_End"] = self.absolute_position(rel_end, zero_based=False)
        return motif


# ─── Parsing helpers ──────────────────────────────────────────────────────────


def validate_accession(accession: str) -> bool:
    """Return True if *accession* looks like a valid RefSeq/GenBank accession.

    This performs a lightweight syntactic check only; it does **not** contact
    NCBI to confirm the accession exists.
    """
    if not accession or not isinstance(accession, str):
        return False
    return bool(_ACCESSION_PATTERN.match(accession.strip()))


def parse_interval_string(interval_str: str) -> GenomicInterval:
    """Parse a compact interval string into a :class:`GenomicInterval`.

    Supported format::

        NC_000913.3:100000-150000

    Parameters
    ----------
    interval_str:
        The compact interval string to parse.

    Returns
    -------
    GenomicInterval
        Parsed interval with the detected accession, start, and end.

    Raises
    ------
    ValueError
        If the string does not match the expected format, or the parsed
        coordinates are invalid.
    """
    if not interval_str or not isinstance(interval_str, str):
        raise ValueError("Interval string must be a non-empty string.")

    m = _INTERVAL_COMPACT_RE.match(interval_str.strip())
    if not m:
        raise ValueError(
            f"Cannot parse interval '{interval_str}'. "
            "Expected format: ACCESSION:START-END "
            "(e.g. NC_000913.3:100000-150000)"
        )

    accession = m.group("accession")
    start = int(m.group("start"))
    end = int(m.group("end"))

    # GenomicInterval.__post_init__ validates start/end constraints
    return GenomicInterval(accession=accession, start=start, end=end)


def build_interval(
    accession: str,
    start: int,
    end: int,
    organism: str = "",
    chromosome: str = "",
) -> GenomicInterval:
    """Construct a :class:`GenomicInterval` from structured fields.

    Parameters
    ----------
    accession:
        RefSeq / GenBank accession.
    start:
        1-based interval start (inclusive).
    end:
        1-based interval end (inclusive).
    organism:
        Optional organism name.
    chromosome:
        Optional chromosome / contig label.

    Returns
    -------
    GenomicInterval

    Raises
    ------
    ValueError
        If any coordinate constraint is violated.
    """
    return GenomicInterval(
        accession=accession,
        start=start,
        end=end,
        organism=organism,
        chromosome=chromosome,
    )


# ─── NCBI provider ────────────────────────────────────────────────────────────


def fetch_genome_interval(
    interval: GenomicInterval,
    email: str = "",
    api_key: Optional[str] = None,
) -> Tuple[str, str]:
    """Fetch the DNA sequence for a genomic interval from NCBI Entrez.

    Uses ``Entrez.efetch`` with the ``seq_start`` / ``seq_stop`` parameters
    so only the requested sub-sequence is transferred.

    Parameters
    ----------
    interval:
        The genomic interval to fetch.
    email:
        E-mail address to pass to NCBI (required by NCBI policy).
    api_key:
        Optional NCBI API key for higher rate limits.

    Returns
    -------
    (sequence, record_id)
        *sequence* — uppercase DNA string (U→T substitution applied).
        *record_id* — the NCBI record identifier extracted from the FASTA
                      header, or the accession when not available.

    Raises
    ------
    ImportError
        If BioPython is not installed.
    RuntimeError
        If NCBI returns no sequences or an error response.
    ValueError
        If the interval accession is empty.
    """
    try:
        from Bio import Entrez, SeqIO
    except ImportError as exc:
        raise ImportError(
            "BioPython is required for NCBI retrieval. "
            "Install it with: pip install biopython"
        ) from exc

    if not interval.accession:
        raise ValueError("Interval accession must not be empty.")

    if email:
        Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    logger.info(
        "Fetching NCBI interval %s [%d–%d] (%d bp)",
        interval.accession,
        interval.start,
        interval.end,
        interval.length,
    )

    try:
        handle = Entrez.efetch(
            db="nucleotide",
            id=interval.accession,
            rettype="fasta",
            retmode="text",
            seq_start=interval.start,
            seq_stop=interval.end,
        )
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
    except Exception as exc:
        raise RuntimeError(
            f"NCBI retrieval failed for {interval.compact}: {exc}"
        ) from exc

    if not records:
        raise RuntimeError(
            f"NCBI returned no sequences for {interval.compact}. "
            "Check that the accession and coordinates are valid."
        )

    record = records[0]
    sequence = str(record.seq).upper().replace("U", "T")
    record_id = record.id or interval.accession

    logger.info(
        "Retrieved %d bp for %s (record id: %s)",
        len(sequence),
        interval.compact,
        record_id,
    )

    return sequence, record_id


# ─── Provider architecture ────────────────────────────────────────────────────


class IntervalProvider(ABC):
    """Abstract base class for genomic interval data providers.

    Subclass this to add support for Ensembl, UCSC, local genomes, etc.
    Register new providers in :data:`INTERVAL_PROVIDERS`.
    """

    name: str = ""
    """Short identifier used in :data:`INTERVAL_PROVIDERS` and UI labels."""

    description: str = ""
    """Human-readable description shown in provider selection UI."""

    @abstractmethod
    def fetch(
        self,
        interval: GenomicInterval,
        **kwargs,
    ) -> Tuple[str, str]:
        """Fetch the DNA sequence for *interval*.

        Returns
        -------
        (sequence, record_id)
        """

    @abstractmethod
    def validate_accession(self, accession: str) -> bool:
        """Return True if *accession* is valid for this provider."""


class NCBIProvider(IntervalProvider):
    """NCBI Entrez provider — uses RefSeq / GenBank accessions."""

    name = "ncbi"
    description = "NCBI Entrez (RefSeq/GenBank)"

    def __init__(self, email: str = "", api_key: Optional[str] = None) -> None:
        self.email = email
        self.api_key = api_key

    def fetch(self, interval: GenomicInterval, **kwargs) -> Tuple[str, str]:
        return fetch_genome_interval(
            interval,
            email=kwargs.get("email", self.email),
            api_key=kwargs.get("api_key", self.api_key),
        )

    def validate_accession(self, accession: str) -> bool:
        return validate_accession(accession)


#: Registry of available interval providers.  Add new providers here.
INTERVAL_PROVIDERS: Dict[str, IntervalProvider] = {
    NCBIProvider.name: NCBIProvider(),
}


def get_provider(name: str) -> IntervalProvider:
    """Look up a provider by name.

    Raises
    ------
    KeyError
        If *name* is not registered in :data:`INTERVAL_PROVIDERS`.
    """
    if name not in INTERVAL_PROVIDERS:
        available = ", ".join(INTERVAL_PROVIDERS)
        raise KeyError(
            f"Unknown interval provider '{name}'. "
            f"Available providers: {available}"
        )
    return INTERVAL_PROVIDERS[name]


# ─── Validation helpers ───────────────────────────────────────────────────────


def validate_interval_inputs(
    accession: str,
    start: Optional[int],
    end: Optional[int],
    *,
    max_length: int = 50_000_000,
) -> List[str]:
    """Validate structured interval inputs and return a list of error messages.

    An empty list means all inputs are valid.

    Parameters
    ----------
    accession:
        The accession string to validate.
    start:
        1-based interval start coordinate (or None if not provided).
    end:
        1-based interval end coordinate (or None if not provided).
    max_length:
        Maximum allowed interval length in base-pairs (default 50 Mbp).
    """
    errors: List[str] = []

    if not accession or not accession.strip():
        errors.append("Accession is required.")
    elif not validate_accession(accession.strip()):
        errors.append(
            f"'{accession}' does not look like a valid RefSeq or GenBank accession. "
            "Expected formats: NC_000913.3, NM_001234.5, AB012345.1, etc."
        )

    if start is None:
        errors.append("Start coordinate is required.")
    elif start < 1:
        errors.append(f"Start coordinate must be ≥ 1 (got {start}).")

    if end is None:
        errors.append("End coordinate is required.")
    elif end < 1:
        errors.append(f"End coordinate must be ≥ 1 (got {end}).")

    if start is not None and end is not None and start >= 1 and end >= 1:
        if end < start:
            errors.append(
                f"End coordinate ({end}) must be ≥ start coordinate ({start})."
            )
        elif (end - start + 1) > max_length:
            errors.append(
                f"Interval is {end - start + 1:,} bp — larger than the "
                f"maximum allowed {max_length:,} bp. "
                "Consider splitting into smaller sub-intervals."
            )

    return errors
