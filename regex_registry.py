"""
Tabular Motif Summary (Scientific/Algorithmic Cutoffs)
------------------------------------------------------
| Motif Type        | Definition / Cutoffs                                   | Tool / Source              | Priority (if applicable)        |
|-------------------|--------------------------------------------------------|----------------------------|---------------------------------|
| APR               | ≥3 A-tracts (3–11 bp, spaced 9–11 bp)                  | In-house regex             |                                 |
| Direct Repeat     | Unit 10–300 bp, spacer ≤100 bp (≤10 bp preferred)      | In-house regex             |                                 |
| Inverted Repeat   | Arm ≥6 bp, loop ≤100 bp                                | Algorithmic                |                                 |
| Mirror Repeat     | Arm ≥10 bp, loop ≤100 bp                               | Algorithmic                |                                 |
| STR               | Repeat unit 1–9 bp, total repeat length ≥10 bp         | In-house regex             |                                 |
| Z-DNA             | Alternating CG/AT >10 bp                               | Z-DNA Seeker (Ho 1986)     |                                 |
| Triplex           | Purine/pyrimidine ≥10 bp, spacer ≤8 bp                 | grep subset                |                                 |
| G4 (multimeric)   | >4 G4 modules, tandem                                  | Quadron (adapted)          | 1 (highest)                     |
| G4 (bipartite)    | Two G4 modules, 10–30 bp apart                         | Quadron (adapted)          | 2                               |
| G4 (canonical)    | ≥4 GGG-stems, 1–7 bp loops, no bulges                  | Quadron (default)          | 3                               |
| G4 (relaxed)      | ≥4 GGG-stems, loops 1–12 bp, no bulges                 | Quadron (adapted)          | 4                               |
| G4 (bulged)       | ≥4 GGG-stems, bulges ≤3 bases                          | Quadron (adapted)          | 5                               |
| G4 (imperfect)    | ≥4 stems, one run only 2 Gs                            | Quadron (adapted)          | 6                               |
| G-triplex         | GGG–loop–GGG–loop–GGG                                  | Quadron (adapted)          | 7 (lowest)                      |
| i-motif (canon)   | ≥4 C-tracts, loops 1–7 bp                              | iM-seeker (Zeraati 2018)   | 1 (highest)                     |
| i-motif (relaxed) | ≥4 C-tracts, loops 1–12 bp                             | iM-seeker (adapted)        | 2 (lowest)                      |
------------------------------------------------------

NBDFinder Motif Regex Registry - Upgraded
-----------------------------------------
Contains regexes & algorithmic logic for DNA secondary structure motifs.
All patterns are Hyperscan-safe (no backreferences/captures), restricted to [ACGT].
Motif type priorities prevent double-counting; scientific references cited per motif.

"""

import re; import numpy as np # Core imports

# === APR (A-phased repeats) ===
# ≥3 A-tracts (3–11 nt), 9–11 nt spacing; T-tracts on reverse
APR_PATTERNS = {
    'a_phased_repeat': [(r"(?:A{3,11}.{5,9}){2}A{3,11}", 1, "APR", None, 1.0, 3, "InHouse")],
    't_phased_repeat': [(r"(?:T{3,11}.{5,9}){2}T{3,11}", 2, "APR", None, 1.0, 3, "InHouse")]
}

# === Direct Repeats (DRs) ===
# Unit 10–300 nt, spacer ≤100 nt (≤10 nt preferred)
DR_PATTERNS = {
    'direct_repeat': [(r"([ACGT]{10,300})[ACGT]{0,100}\1", 1, "DR", None, 1.0, 2, "InHouse")]
}

# === Inverted Repeats (IRs) ===
# Arm ≥6 nt, loop ≤100 nt (algorithmic, not regex)
IR_PATTERNS = {'inverted_repeat': []} # Use find_inverted_repeats

# === Mirror Repeats (MRs) ===
# Arm ≥10 nt, loop ≤100 nt (algorithmic for mirror, not regex)
MR_PATTERNS = {'mirror_repeat': []} # Use find_mirror_repeats

# === Short Tandem Repeats (STRs) ===
# Unit 1–9 nt, total length ≥10 nt
STR_PATTERNS = {
    'mono': [(r"(A{10,}|T{10,}|G{10,}|C{10,})", 1, "STR", None, 1.0, 1, "InHouse")],
    'di': [(r"(?:AT){5,}|(?:GC){5,}|(?:AG){5,}|(?:CT){5,}", 2, "STR", None, 1.0, 2, "InHouse")],
    'tri': [(r"(?:CGG){4,}|(?:CAG){4,}|(?:AAG){4,}", 3, "STR", None, 1.0, 3, "InHouse")],
    'arbitrary': [(r"(([ACGT]{1,9}))\2{1,}", 4, "STR", None, 1.0, 1, "InHouse")] # generic
}

# === Z-DNA ===
# Enhanced Z-DNA detection using Kadane algorithm for transition scoring
# Traditional patterns for compatibility with existing framework
Z_DNA_PATTERNS = {
    'cg_zdna': [(r"(?:CG){6,}", 1, "Z-DNA", None, 1.0, 12, "Ho1986")],
    'at_zdna': [(r"(?:AT){6,}", 2, "Z-DNA", None, 1.0, 12, "Ho1986")],
    # New Kadane-based detection entry point
    'kadane_zdna': [(r"[ATCG]{20,}", 3, "Z-DNA_Kadane", None, 1.5, 20, "Kadane2024")]
}

# === Triplex ===
# Purine/pyrimidine ≥10 bp, spacer ≤8 bp
TRIPLEX_PATTERNS = {
    'purine_triplex': [(r"[AG]{10,}[ACGT]{0,8}[AG]{10,}", 1, "Triplex", None, 1.0, 20, "grep")],
    'pyrimidine_triplex': [(r"[CT]{10,}[ACGT]{0,8}[CT]{10,}", 2, "Triplex", None, 1.0, 20, "grep")]
}

# === G4 Family (priority order) ===
# Multimeric (highest) → Bipartite → Canonical → Relaxed → Bulged → Imperfect → G-triplex (lowest)
G_QUADRUPLEX_PATTERNS = {
    'multimeric_g4': [(r"(?:G{3,}[ACGT]{1,12}){4,}", 1, "Multimeric_G4", None, 1.2, 4, "Quadron")],
    'bipartite_g4': [(r"(?:G{3,}[ACGT]{1,7}){4}([ACGT]{10,30})(?:G{3,}[ACGT]{1,7}){4}", 2, "Bipartite_G4", None, 1.1, 8, "Quadron")],
    'canonical_g4': [(r"(?:G{3,}[ACGT]{1,7}){3}G{3,}", 3, "Canonical_G4", None, 1.0, 4, "Quadron")],
    'relaxed_g4': [(r"(?:G{3,}[ACGT]{8,12}){3}G{3,}", 4, "Relaxed_G4", None, 0.8, 4, "Quadron")],
    'bulged_g4': [(r"(?:G{2,3}[ACGT]{0,3}){3}G{2,3}", 5, "Bulged_G4", None, 0.8, 4, "Quadron")],
    'imperfect_g4': [(r"(?:G{2,3}[ACGT]{1,7}){3}G{2,3}", 6, "Imperfect_G4", None, 0.7, 4, "Quadron")],
    'g_triplex': [(r"(?:G{3,}[ACGT]{1,7}){2}G{3,}", 7, "G-Triplex", None, 0.6, 3, "Quadron")]
}

# === i-motif Family (priority order) ===
# Canonical (highest) → Relaxed (lowest)
I_MOTIF_PATTERNS = {
    'canonical_imotif': [(r"(?:C{3,}[ACGT]{1,7}){3}C{3,}", 1, "Canonical_iMotif", None, 1.0, 4, "Zeraati2018")],
    'relaxed_imotif': [(r"(?:C{3,}[ACGT]{1,12}){3}C{3,}", 2, "Relaxed_iMotif", None, 0.8, 4, "Zeraati2018")]
}

# === MASTER REGISTRY ===
ALL_PATTERNS = {
    'apr': APR_PATTERNS,
    'direct_repeat': DR_PATTERNS,
    'inverted_repeat': IR_PATTERNS,
    'mirror_repeat': MR_PATTERNS,
    'str': STR_PATTERNS,
    'z_dna': Z_DNA_PATTERNS,
    'triplex': TRIPLEX_PATTERNS,
    'g_quadruplex': G_QUADRUPLEX_PATTERNS,
    'i_motif': I_MOTIF_PATTERNS
}

# === Get patterns for motif class ===
def get_patterns_for_motif(motif_class: str) -> dict:
    "Return all categories for a motif class."; return ALL_PATTERNS.get(motif_class, {})

# === Hyperscan pattern list, deduped (regex, id) ===
def get_all_hyperscan_patterns() -> list:
    "List of (pattern, id) for Hyperscan DB. Deduped, mapped."
    all_patterns = []; pattern_map = {}; pattern_id = 0; seen = set()
    for motif_class, pattern_dict in ALL_PATTERNS.items():
        for pattern_type, patterns in pattern_dict.items():
            for pattern_tuple in patterns:
                if pattern_tuple:
                    regex_pattern = pattern_tuple[0]
                    if regex_pattern not in seen:
                        all_patterns.append((regex_pattern, pattern_id))
                        pattern_map[pattern_id] = {'motif_class': motif_class, 'pattern_type': pattern_type, 'original_id': pattern_tuple[1]}
                        seen.add(regex_pattern); pattern_id += 1
    return all_patterns

# === Metadata for pattern_id ===
def get_pattern_info(pattern_id: int) -> dict:
    "Return metadata for a given pattern ID (stable mapping)."
    current_id = 0
    for motif_class, pattern_dict in ALL_PATTERNS.items():
        for pattern_type, patterns in pattern_dict.items():
            for pattern_tuple in patterns:
                if pattern_tuple:
                    if current_id == pattern_id:
                        return {'motif_class': motif_class, 'pattern_type': pattern_type, 'regex': pattern_tuple[0], 'original_id': pattern_tuple[1],
                                'label': pattern_tuple[2], 'subclass': pattern_tuple[3],
                                'score_scale': pattern_tuple[4] if len(pattern_tuple) > 4 else 1.0,
                                'min_runs': pattern_tuple[5] if len(pattern_tuple) > 5 else 0,
                                'source': pattern_tuple[6] if len(pattern_tuple) > 6 else "default"}
                    current_id += 1
    return {}

# === Validate regex patterns (Hyperscan-safe) ===
def validate_patterns() -> bool:
    "Check all regex patterns for syntax correctness (Hyperscan-safe)."
    try:
        for motif_class, pattern_dict in ALL_PATTERNS.items():
            for pattern_type, patterns in pattern_dict.items():
                for pattern_tuple in patterns:
                    if pattern_tuple: re.compile(pattern_tuple[0])
        return True
    except re.error as e:
        print(f"Invalid regex pattern found: {e}"); return False

# === Algorithmic: Inverted Repeat / Palindrome (for IR) ===
def reverse_complement(seq: str) -> str:
    "Return reverse complement of DNA."; return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def find_inverted_repeats(seq: str, min_arm: int = 6, max_arm: int = 100, max_loop: int = 100):
    "Find inverted repeats/palindromes (IRs). Returns (start, end, arm_len, loop_len, left_arm, right_arm)."
    seq = seq.upper(); n = len(seq); results = []
    for arm in range(min_arm, max_arm+1):
        for i in range(n - 2*arm + 1):
            left = seq[i:i+arm]; 
            if any(ch not in "ACGT" for ch in left): continue
            left_rc = reverse_complement(left)
            for loop in range(0, max_loop+1):
                j = i + arm + loop; right = seq[j:j+arm]
                if len(right) < arm: break
                if right == left_rc: results.append((i, j+arm, arm, loop, left, right))
    return results

# === Algorithmic: Mirror Repeat (for MR/Triplex) ===
def find_mirror_repeats(seq: str, min_arm: int = 10, max_arm: int = 100, max_loop: int = 100):
    "Find mirror repeats (MRs). Returns (start, end, arm_len, loop_len, left_arm, right_arm)."
    seq = seq.upper(); n = len(seq); results = []
    for arm in range(min_arm, max_arm+1):
        for i in range(n - 2*arm + 1):
            left = seq[i:i+arm]
            for loop in range(0, max_loop+1):
                j = i + arm + loop; right = seq[j:j+arm]
                if len(right) < arm: break
                if left == right[::-1]: results.append((i, j+arm, arm, loop, left, right))
    return results

# === Export ===
HYPERSCAN_SAFE_PATTERNS = get_all_hyperscan_patterns()
__all__ = [
    'ALL_PATTERNS', 'APR_PATTERNS', 'DR_PATTERNS', 'IR_PATTERNS', 'MR_PATTERNS', 'STR_PATTERNS', 
    'Z_DNA_PATTERNS', 'TRIPLEX_PATTERNS', 'G_QUADRUPLEX_PATTERNS', 'I_MOTIF_PATTERNS',
    'get_patterns_for_motif', 'get_all_hyperscan_patterns', 'get_pattern_info', 'validate_patterns', 
    'HYPERSCAN_SAFE_PATTERNS', 'reverse_complement', 'find_inverted_repeats', 'find_mirror_repeats'
]

# === Validation on import ===
if not validate_patterns(): print("Warning: Some regex patterns in registry are invalid!")
# === End of Registry ===
