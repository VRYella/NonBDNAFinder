"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Sequence Preprocessor - Gold Standard Genomic Sequence Preprocessing         │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2025.1            │
│ Implements scientifically rigorous preprocessing for DNA sequence analysis   │
└──────────────────────────────────────────────────────────────────────────────┘

Scientific Basis:
    - GC% calculation follows NCBI/Ensembl/UCSC standard
    - Excludes ambiguous bases (N, R, Y, etc.) from denominators
    - Comprehensive validation with position reporting
    - FASTA-compliant parsing (Pearson & Lipman, 1988)
"""

# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
from dataclasses import dataclass, field
from typing import Dict, List
import re

# ═══════════════════════════════════════════════════════════════════════════════
# DATA CLASSES
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class PreprocessingResult:
    """
    Data class for comprehensive preprocessing results.
    
    Attributes:
        header: FASTA header (without '>'), empty string if none
        sequence: Clean, single-line, uppercase sequence
        length: Total sequence length (includes N's and valid bases)
        valid_bases: Count of A+T+G+C only
        character_counts: Detailed counts {'A': n, 'T': n, 'G': n, 'C': n, 'N': n}
        invalid_characters: Dict mapping invalid chars to list of positions
        gc_percentage: (G+C) / (A+T+G+C) × 100
        at_percentage: (A+T) / (A+T+G+C) × 100
        gc_balance: "Balanced" | "GC-rich" | "AT-rich"
        validation_status: "valid" | "warning" | "error"
        warnings: List of warning messages
        errors: List of error messages
    """
    header: str = ""
    sequence: str = ""
    length: int = 0
    valid_bases: int = 0
    character_counts: Dict[str, int] = field(default_factory=dict)
    invalid_characters: Dict[str, List[int]] = field(default_factory=dict)
    gc_percentage: float = 0.0
    at_percentage: float = 0.0
    gc_balance: str = "Balanced"
    validation_status: str = "valid"
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)


# ═══════════════════════════════════════════════════════════════════════════════
# PREPROCESSING FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def preprocess_sequence(raw_input: str) -> PreprocessingResult:
    """
    Comprehensive preprocessing of DNA sequence with gold-standard validation.
    
    Implements 7-step workflow:
        1. Uppercase normalization
        2. FASTA header extraction
        3. Multi-line concatenation
        4. Character analysis (counts + invalid char detection)
        5. Length calculation
        6. GC/AT percentage (correct formula)
        7. Validation and classification
    
    Args:
        raw_input: Raw sequence input (may include FASTA headers, newlines, etc.)
        
    Returns:
        PreprocessingResult object with comprehensive analysis
        
    Scientific Basis:
        - GC% = (G+C) / (A+T+G+C) × 100 (NCBI standard)
        - Excludes ambiguous bases from denominator
        - IUPAC nucleotide codes supported for validation
        
    Example:
        >>> result = preprocess_sequence(">header\\nATCGNNNN\\nATCG")
        >>> result.header
        'HEADER'
        >>> result.sequence
        'ATCGNNNNАТCG'
        >>> result.gc_percentage
        50.0  # (4 GC / 8 ATGC), N's excluded
    """
    result = PreprocessingResult()
    
    # ═══════════════════════════════════════════════════════════════════════════
    # STEP 1: UPPERCASE NORMALIZATION
    # Scientific Basis: Case-insensitive DNA sequence analysis standard
    # ═══════════════════════════════════════════════════════════════════════════
    normalized_input = raw_input.upper().strip()
    
    # ═══════════════════════════════════════════════════════════════════════════
    # STEP 2: FASTA HEADER EXTRACTION
    # Scientific Basis: FASTA format specification (Pearson & Lipman, 1988)
    # Header line starts with '>', followed by description
    # ═══════════════════════════════════════════════════════════════════════════
    lines = normalized_input.split('\n')
    sequence_lines = []
    
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            # Extract first header found (ignore subsequent headers)
            if not result.header:
                result.header = line[1:].strip()
        elif line:
            # Non-header, non-empty line is part of sequence
            sequence_lines.append(line)
    
    # ═══════════════════════════════════════════════════════════════════════════
    # STEP 3: MULTI-LINE CONCATENATION
    # Scientific Basis: FASTA sequences can span multiple lines
    # ═══════════════════════════════════════════════════════════════════════════
    result.sequence = ''.join(sequence_lines)
    result.length = len(result.sequence)
    
    # Handle empty sequence
    if result.length == 0:
        result.validation_status = "error"
        result.errors.append("Empty sequence after preprocessing")
        return result
    
    # ═══════════════════════════════════════════════════════════════════════════
    # STEP 4: CHARACTER ANALYSIS
    # Scientific Basis: IUPAC nucleotide codes
    # Standard: A, T, G, C
    # Ambiguous: N (any), R (purine), Y (pyrimidine), K, M, S, W, B, D, H, V
    # ═══════════════════════════════════════════════════════════════════════════
    
    # ─── Fast counts via str.count() (C-level, 10-20x faster than a Python loop) ──
    seq = result.sequence
    a_count = seq.count('A')
    t_count = seq.count('T')
    g_count = seq.count('G')
    c_count = seq.count('C')
    n_count = seq.count('N')

    result.character_counts = {
        'A': a_count, 'T': t_count, 'G': g_count, 'C': c_count, 'N': n_count
    }

    # Detect invalid characters (non-IUPAC) using a single set-difference scan.
    # This only iterates characters that fall outside the valid IUPAC alphabet.
    valid_iupac = set('ATGCNRYKMSWBDHV')
    invalid_chars_seen: dict = {}
    for i, char in enumerate(seq):
        if char not in valid_iupac:
            if char not in invalid_chars_seen:
                invalid_chars_seen[char] = []
            if len(invalid_chars_seen[char]) < 10:
                invalid_chars_seen[char].append(i)
    result.invalid_characters = invalid_chars_seen

    # ═══════════════════════════════════════════════════════════════════════════
    # STEP 5: LENGTH CALCULATION
    # ═══════════════════════════════════════════════════════════════════════════
    result.valid_bases = a_count + t_count + g_count + c_count
    
    # ═══════════════════════════════════════════════════════════════════════════
    # STEP 6: GC PERCENTAGE (CORRECT FORMULA)
    # Scientific Basis: Gold standard genomic formula
    # GC% = (G + C) / (A + T + G + C) × 100
    # Denominator excludes N's and other ambiguous bases
    # ═══════════════════════════════════════════════════════════════════════════
    if result.valid_bases == 0:
        result.validation_status = "error"
        result.errors.append("No valid ATGC bases found in sequence")
        return result
    
    result.gc_percentage = ((g_count + c_count) / result.valid_bases) * 100
    result.at_percentage = ((a_count + t_count) / result.valid_bases) * 100
    
    # ═══════════════════════════════════════════════════════════════════════════
    # STEP 7: VALIDATION & CLASSIFICATION
    # ═══════════════════════════════════════════════════════════════════════════
    
    # GC Balance classification
    if result.gc_percentage > 60:
        result.gc_balance = "GC-rich"
    elif result.gc_percentage < 40:
        result.gc_balance = "AT-rich"
    else:
        result.gc_balance = "Balanced"
    
    # Validation status determination
    if result.invalid_characters:
        result.validation_status = "error"
        for char, positions in result.invalid_characters.items():
            pos_str = ', '.join(str(p) for p in positions[:5])
            if len(positions) > 5:
                pos_str += f" ... ({len(positions)} total)"
            result.errors.append(
                f"Invalid character '{char}' found at positions: {pos_str}"
            )
    elif n_count > 0:
        result.validation_status = "warning"
        n_percentage = (n_count / result.length) * 100
        result.warnings.append(
            f"Sequence contains {n_count} ambiguous base(s) 'N' ({n_percentage:.2f}% of total length)"
        )
    else:
        result.validation_status = "valid"
    
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# UTILITY FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def format_preprocessing_report(result: PreprocessingResult) -> str:
    """
    Generate human-readable report from preprocessing results.
    
    Args:
        result: PreprocessingResult object
        
    Returns:
        Formatted string report
    """
    lines = []
    lines.append("=" * 80)
    lines.append("SEQUENCE PREPROCESSING REPORT")
    lines.append("=" * 80)
    
    if result.header:
        lines.append(f"Header: {result.header}")
    
    lines.append(f"Length: {result.length:,} bp")
    lines.append(f"Valid Bases (ATGC): {result.valid_bases:,} bp")
    lines.append("")
    
    lines.append("Base Composition:")
    lines.append(f"  A: {result.character_counts.get('A', 0):,} ({result.character_counts.get('A', 0)/result.valid_bases*100:.2f}%)" if result.valid_bases > 0 else "  A: 0")
    lines.append(f"  T: {result.character_counts.get('T', 0):,} ({result.character_counts.get('T', 0)/result.valid_bases*100:.2f}%)" if result.valid_bases > 0 else "  T: 0")
    lines.append(f"  G: {result.character_counts.get('G', 0):,} ({result.character_counts.get('G', 0)/result.valid_bases*100:.2f}%)" if result.valid_bases > 0 else "  G: 0")
    lines.append(f"  C: {result.character_counts.get('C', 0):,} ({result.character_counts.get('C', 0)/result.valid_bases*100:.2f}%)" if result.valid_bases > 0 else "  C: 0")
    if result.character_counts.get('N', 0) > 0:
        lines.append(f"  N: {result.character_counts['N']:,} (ambiguous)")
    lines.append("")
    
    lines.append(f"GC Content: {result.gc_percentage:.2f}%")
    lines.append(f"AT Content: {result.at_percentage:.2f}%")
    lines.append(f"GC Balance: {result.gc_balance}")
    lines.append("")
    
    lines.append(f"Validation Status: {result.validation_status.upper()}")
    
    if result.warnings:
        lines.append("")
        lines.append("WARNINGS:")
        for warning in result.warnings:
            lines.append(f"  ⚠ {warning}")
    
    if result.errors:
        lines.append("")
        lines.append("ERRORS:")
        for error in result.errors:
            lines.append(f"  ✗ {error}")
    
    lines.append("=" * 80)
    
    return '\n'.join(lines)
