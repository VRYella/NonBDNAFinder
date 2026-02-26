"""SMID Seed Database: Canonical seed definitions for SMID pre-screening layer.

Defines SMID_SEEDS as a structured list of seed entries.  Each entry has:
    - pattern   : exact string (10-mer) or compiled regex string
    - is_regex  : True if pattern is a regular-expression, False for exact match
    - classes   : list of motif class names this seed belongs to
    - max_extension : bp extension added on each side of a cluster window
    - epsilon   : maximum gap (bp) between hits to merge into one cluster
    - min_samples : minimum number of hits required to form a cluster

Design constraints
------------------
* A-philic and Z-DNA seeds use the original 10-mer tables verbatim.
* No permissive single-run seeds (A{3,}, G{3,}, CG, etc.).
* Curved_DNA and Slipped_DNA are excluded from SMID v1.
* G-Quadruplex and R-Loop share the same two-tract core seed; hit positions
  are recorded for both classes simultaneously.
"""

from typing import Any, Dict, List

from Detectors.aphilic.tenmer_table import TENMER_LOG2
from Detectors.zdna.tenmer_table import TENMER_SCORE

# ---------------------------------------------------------------------------
# Build seed list
# ---------------------------------------------------------------------------

SMID_SEEDS: List[Dict[str, Any]] = []

# ── A-philic ──────────────────────────────────────────────────────────────
for _kmer in TENMER_LOG2:
    SMID_SEEDS.append(
        {
            "pattern": _kmer,
            "is_regex": False,
            "classes": ["A-philic_DNA"],
            "max_extension": 300,
            "epsilon": 150,
            "min_samples": 2,
        }
    )

# ── Z-DNA ─────────────────────────────────────────────────────────────────
for _kmer in TENMER_SCORE:
    SMID_SEEDS.append(
        {
            "pattern": _kmer,
            "is_regex": False,
            "classes": ["Z-DNA"],
            "max_extension": 300,
            "epsilon": 150,
            "min_samples": 2,
        }
    )

# ── G-Quadruplex / R-Loop (shared two-tract core) ─────────────────────────
SMID_SEEDS.append(
    {
        "pattern": r"G{3,}[ACGT]{1,7}G{3,}",
        "is_regex": True,
        "classes": ["G-Quadruplex", "R-Loop"],
        "max_extension": 120,
        "epsilon": 100,
        "min_samples": 2,
    }
)

# ── i-Motif ───────────────────────────────────────────────────────────────
SMID_SEEDS.append(
    {
        "pattern": r"C{3,}[ACGT]{1,7}C{3,}",
        "is_regex": True,
        "classes": ["i-Motif"],
        "max_extension": 60,
        "epsilon": 80,
        "min_samples": 2,
    }
)

# ── Triplex ───────────────────────────────────────────────────────────────
SMID_SEEDS.append(
    {
        "pattern": r"[AG]{12,}",
        "is_regex": True,
        "classes": ["Triplex"],
        "max_extension": 150,
        "epsilon": 120,
        "min_samples": 2,
    }
)
SMID_SEEDS.append(
    {
        "pattern": r"[CT]{12,}",
        "is_regex": True,
        "classes": ["Triplex"],
        "max_extension": 150,
        "epsilon": 120,
        "min_samples": 2,
    }
)
