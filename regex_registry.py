"""
Central Regex Registry for NBDFinder - Hyperscan-safe Patterns
==============================================================
Scientific references and implementation notes included per motif class.
All regexes restricted to [ACGT], non-capturing groups, and Hyperscan-safe constructs.
Patterns requiring reverse-complement equality or backreferences are detected algorithmically.
| Motif Class      | Arm/Repeat Length | Loop Length | Detection Criterion                                 | GC Content Check | Purine/Pyrimidine Check | Periodicity/Phasing | Other Notes                            |
|------------------|-------------------|-------------|----------------------------------------------------|------------------|------------------------|---------------------|----------------------------------------|
| G-quadruplex     | ≥3 Gs, 4 runs     | 1–7 nt      | Four G-runs separated by short loops ([ACGT] only)  | Optional         | Not required           | Not required        | Bulged/Imperfect variants loosen runs  |
| G-triplex        | ≥3 Gs, 3 runs     | 1–7 nt      | Three G-runs, as lower-order motif                 | Optional         | Not required           | Not required        | Distinct from canonical G4             |
| I-motif          | ≥3 Cs, 4 runs     | 1–12 nt     | Four C-runs ([ACGT] loops)                         | Optional         | Not required           | Not required        | Strand-specific scoring recommended    |
| Z-DNA            | Dinucleotide      | N/A         | (CG){6+} or G[CG]{8+}G                             | Required         | Not required           | Not required        | Non-capturing groups for Hyperscan     |
| Triplex          | ≥15 bases         | N/A         | Homopurine ([AG]{15+})/Homopyrimidine ([CT]{15+})  | Optional         | Required               | Not required        | Mirror repeat scoring recommended      |
| Curved DNA       | ≥7 A/Ts           | N/A         | Poly-A/T runs ([A/T]{7+})                          | Not required     | Not required           | Required            | Periodicity detection (10 bp spacing)  |
| Cruciform IR     | 6–20 nt           | 0–50 nt     | Reverse complement equality (algorithmic, not regex)| Optional         | Not required           | Not required        | Use fast inverted repeat finder        |
| STR/Slipped DNA  | Mono/di/tri       | N/A         | Literal repeat patterns (e.g. (AT){8+})            | Optional         | Not required           | Not required        | Arbitrary motifs: use k-mer algorithm  |
| R-loop           | ≥3-4 Gs           | 1–10 nt     | G-rich runs, RLFS-like patterns                    | Optional         | Not required           | Not required        | Thermodynamic scoring recommended      |
"""

import re
import numpy as np

# === G-QUADRUPLEX FAMILY PATTERNS (Class 6) ===
# Canonical, bulged, imperfect, and G-triplex forms. All loop regions use [ACGT].
G_QUADRUPLEX_PATTERNS = {
    'canonical_g4': [
        # Four G-runs (>=3 Gs) with loops 1-7 nt
        (r"(?:G{3,}[ACGT]{1,7}){3}G{3,}", 1, 0, "Canonical_G4", None, 1.0, 4, 0.5, "G4Hunter_v2_raw"),
    ],
    'bulged_g4': [
        # Small bulges allowed (handled in postfilter)
        (r"(?:G{2,3}[ACGT]{0,3}){3}G{2,3}", 2, 0, "Bulged_G4", None, 0.8, 4, 0.0, "G4Hunter_Bulge_raw"),
    ],
    'imperfect_g4': [
        # G-run relaxation, loop length as in canonical
        (r"(?:G{2,3}[ACGT]{1,7}){3}G{2,3}", 3, 0, "Imperfect_G4", None, 1.0, 4, 0.4, "G4Hunter_Imperfect_raw"),
    ],
    'relaxed_g4': [
        # Long loops (8-12 nt)
        (r"(?:G{3,}[ACGT]{8,12}){3}G{3,}", 4, 0, "Relaxed_G4", None, 0.8, 4, 0.3, "G4Hunter_LongLoop_raw"),
    ],
    'multimeric_g4': [
        # Multimeric: more than 4 runs, long loops
        (r"(?:G{3,}[ACGT]{1,12}){4,}", 5, 0, "Multimeric_G4", None, 1.2, 4, 0.3, "G4Hunter_Multimer_raw"),
    ],
    'g_triplex': [
        # Three G-runs (lower-order motif, not G4)
        (r"(?:G{3,}[ACGT]{1,7}){2}G{3,}", 6, 0, "G-Triplex_intermediate", None, 1.0, 3, 0.0, "G3_raw"),
    ]
}

# === I-MOTIF FAMILY PATTERNS (Class 7) ===
# C-rich quadruplex structures, all with [ACGT] loops
I_MOTIF_PATTERNS = {
    'canonical_imotif': [
        (r"(?:C{3,}[ACGT]{1,12}){3}C{3,}", 1, 0, "Canonical_iMotif", None, 1.0, 4, 0.0, "G4Hunter_adapted"),
    ],
    'ac_motif': [
        (r"A{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}", 2, 0, "AC-motif", None, 1.0, 0, 0.0, "G4Hunter_adapted"),
        (r"C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}A{3}", 3, 0, "AC-motif", None, 1.0, 0, 0.0, "G4Hunter_adapted")
    ]
}

# === Z-DNA PATTERNS (Class 8) ===
# CG dinucleotide repeats, Hyperscan-safe grouping
Z_DNA_PATTERNS = {
    'z_dna_basic': [
        (r"(?:[CG]{2}){6,}", 1, 0, "Z-DNA", None, 1.0, 0, 0.0, "Z_seeker_raw"),
    ],
    'extended_gz': [
        (r"G[CG]{8,}G", 2, 0, "eGZ (Extruded-G) DNA", None, 1.0, 0, 0.0, "Z_seeker_raw")
    ]
}

# === CURVED DNA PATTERNS (Class 1) ===
# A-tract and T-tract runs (prefilter only)
CURVED_DNA_PATTERNS = {
    'poly_a_tracts': [
        (r"A{7,}", 1, 0, "Global Array", None, 1.0, 0, 0.0, "Curvature_raw"),
    ],
    'poly_t_tracts': [
        (r"T{7,}", 2, 0, "Global Array", None, 1.0, 0, 0.0, "Curvature_raw"),
    ]
}

# === TRIPLEX PATTERNS (Class 5) ===
# Homopurine/homopyrimidine runs; refined by algorithmic mirror repeat test
TRIPLEX_PATTERNS = {
    'homopurine': [
        (r"[AG]{15,}", 1, 0, "Triplex", None, 1.0, 0, 0.0, "Triplex_stability_raw"),
    ],
    'homopyrimidine': [
        (r"[CT]{15,}", 2, 0, "Sticky DNA", None, 1.0, 0, 0.0, "Triplex_stability_raw"),
    ]
}

# === CRUCIFORM PATTERNS (Class 3) ===
# NO regex patterns emitted; use algorithmic inverted repeat finder below
CRUCIFORM_PATTERNS = {
    'palindrome_candidates': [],
    'inverted_repeat_candidates': []
}

# === R-LOOP PATTERNS (Class 4) ===
# RLFS-like G-rich patterns
R_LOOP_PATTERNS = {
    'rlfs_m1': [
        (r"G{3,}[ACGT]{1,10}?G{3,}(?:[ACGT]{1,10}?G{3,}){1,}", 1, 0, "RLFS_m1", None, 1.0, 0, 0.0, "QmRLFS_raw")
    ],
    'rlfs_m2': [
        (r"G{4,}(?:[ACGT]{1,10}?G{4,}){1,}", 2, 0, "RLFS_m2", None, 1.0, 0, 0.0, "QmRLFS_raw")
    ]
}

# === SLIPPED DNA PATTERNS (Class 2) ===
# STRs with literal motifs only; arbitrary motif STRs require streaming algorithm
SLIPPED_DNA_PATTERNS = {
    'repetitive_candidates': [
        (r"A{10,}", 1, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"T{10,}", 2, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"G{10,}", 3, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"C{10,}", 4, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"(?:AT){8,}", 5, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"(?:GC){8,}", 6, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"(?:AG){8,}", 7, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"(?:CT){8,}", 8, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"(?:CGG){6,}", 9, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"(?:CAG){6,}", 10, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
        (r"(?:AAG){6,}", 11, 0, "STR_candidate", None, 1.0, 0, 0.0, "STR_raw"),
    ]
}

# === MASTER PATTERN REGISTRY ===
ALL_PATTERNS = {
    'g_quadruplex': G_QUADRUPLEX_PATTERNS,
    'i_motif': I_MOTIF_PATTERNS,
    'z_dna': Z_DNA_PATTERNS,
    'curved_dna': CURVED_DNA_PATTERNS,
    'triplex': TRIPLEX_PATTERNS,
    'cruciform': CRUCIFORM_PATTERNS,
    'r_loop': R_LOOP_PATTERNS,
    'slipped_dna': SLIPPED_DNA_PATTERNS,
}

# === UTILITY FUNCTIONS ===

def get_patterns_for_motif(motif_class: str) -> dict:
    """Return all pattern categories for a motif class."""
    return ALL_PATTERNS.get(motif_class, {})

def get_all_hyperscan_patterns() -> list:
    """Return list of (pattern, id) tuples for Hyperscan DB. Deduplicated and mapped."""
    all_patterns = []
    pattern_map = {}  # pattern string → (motif_class, pattern_type, original_id)
    pattern_id = 0
    seen = set()
    for motif_class, pattern_dict in ALL_PATTERNS.items():
        for pattern_type, patterns in pattern_dict.items():
            for pattern_tuple in patterns:
                if pattern_tuple:
                    regex_pattern = pattern_tuple[0]
                    if regex_pattern not in seen:
                        all_patterns.append((regex_pattern, pattern_id))
                        pattern_map[pattern_id] = {
                            'motif_class': motif_class,
                            'pattern_type': pattern_type,
                            'original_id': pattern_tuple[1],
                        }
                        seen.add(regex_pattern)
                        pattern_id += 1
    # optionally save pattern_map to disk for reproducibility
    return all_patterns

def get_pattern_info(pattern_id: int) -> dict:
    """Return pattern metadata for a given ID (stable mapping)."""
    current_id = 0
    for motif_class, pattern_dict in ALL_PATTERNS.items():
        for pattern_type, patterns in pattern_dict.items():
            for pattern_tuple in patterns:
                if pattern_tuple:
                    if current_id == pattern_id:
                        return {
                            'motif_class': motif_class,
                            'pattern_type': pattern_type,
                            'regex': pattern_tuple[0],
                            'original_id': pattern_tuple[1],
                            'group_number': pattern_tuple[2],
                            'subclass': pattern_tuple[3],
                            'score_scale': pattern_tuple[5] if len(pattern_tuple) > 5 else 1.0,
                            'min_runs': pattern_tuple[6] if len(pattern_tuple) > 6 else 0,
                            'min_score': pattern_tuple[7] if len(pattern_tuple) > 7 else 0.0,
                            'score_method': pattern_tuple[8] if len(pattern_tuple) > 8 else "default"
                        }
                    current_id += 1
    return {}

def validate_patterns() -> bool:
    """Check all regex patterns for syntax correctness (Hyperscan-safe)."""
    try:
        for motif_class, pattern_dict in ALL_PATTERNS.items():
            for pattern_type, patterns in pattern_dict.items():
                for pattern_tuple in patterns:
                    if pattern_tuple:
                        regex_pattern = pattern_tuple[0]
                        re.compile(regex_pattern)
        return True
    except re.error as e:
        print(f"Invalid regex pattern found: {e}")
        return False

# === ALGORITHMIC INVERTED REPEAT / PALINDROME DETECTOR ===
# Not expressible as regex; use this for Cruciform candidates
def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA string."""
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def find_inverted_repeats(seq: str, min_arm: int = 6, max_arm: int = 20, max_loop: int = 50):
    """
    Find inverted repeats/palindromes (for Cruciform DNA).
    Returns list of (start, end, arm_len, loop_len, left_arm, right_arm).
    """
    seq = seq.upper()
    n = len(seq)
    results = []
    for arm in range(min_arm, max_arm+1):
        for i in range(n - 2*arm + 1):
            left = seq[i:i+arm]
            if any(ch not in "ACGT" for ch in left): continue
            left_rc = reverse_complement(left)
            for loop in range(0, max_loop+1):
                j = i + arm + loop
                right = seq[j:j+arm]
                if len(right) < arm: break
                if right == left_rc:
                    results.append((i, j+arm, arm, loop, left, right))
    return results

# === EXPORT COMMONLY USED PATTERNS AND FUNCTIONS ===
HYPERSCAN_SAFE_PATTERNS = get_all_hyperscan_patterns()

__all__ = [
    'ALL_PATTERNS',
    'G_QUADRUPLEX_PATTERNS',
    'I_MOTIF_PATTERNS',
    'Z_DNA_PATTERNS',
    'CURVED_DNA_PATTERNS',
    'TRIPLEX_PATTERNS',
    'CRUCIFORM_PATTERNS',
    'R_LOOP_PATTERNS',
    'SLIPPED_DNA_PATTERNS',
    'get_patterns_for_motif',
    'get_all_hyperscan_patterns',
    'get_pattern_info',
    'validate_patterns',
    'HYPERSCAN_SAFE_PATTERNS',
    'reverse_complement',
    'find_inverted_repeats'
]

# === VALIDATION ON IMPORT ===
if not validate_patterns():
    print("Warning: Some regex patterns in registry are invalid!")

# === END OF REGISTRY ===
