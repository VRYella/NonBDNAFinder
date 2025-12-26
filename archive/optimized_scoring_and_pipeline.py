"""
Optimized Motif Detection + Scoring Pipeline
============================================
Design rules (HARD):
- Hyperscan = detection only
- Exactly ONE scan per sequence
- Scoring is independent per motif
- Score range strictly [1.0, 3.0]
- ΔG used ONLY where physically valid
- Visualization is lazy
- Results cached per sequence + config
"""

import hashlib
import gc
import numpy as np
from typing import Dict, List
import matplotlib.pyplot as plt
import pandas as pd

# ============================================================
# 1. NN ΔG°37 (USER PROVIDED — EXACT)
# ============================================================

NN_DELTA_G = {
    "AA": -1.00, "TT": -1.00,
    "AT": -0.88, "TA": -0.58,
    "CA": -1.45, "TG": -1.45,
    "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28,
    "GA": -1.30, "TC": -1.30,
    "CG": -2.17,
    "GC": -2.24,
    "GG": -1.84, "CC": -1.84,
}

def nn_delta_g(seq: str) -> float:
    """
    Calculate nearest-neighbor ΔG°37 for a DNA sequence.
    
    Uses SantaLucia nearest-neighbor parameters (1998, 2004).
    
    Args:
        seq: DNA sequence (ATGC only, case-insensitive)
    
    Returns:
        Free energy in kcal/mol
        
    Raises:
        ValueError: If sequence contains non-ATGC characters or invalid dinucleotides
    """
    seq = seq.upper()
    
    # Validate sequence contains only ATGC
    valid_bases = set('ATGC')
    if not all(base in valid_bases for base in seq):
        invalid = set(seq) - valid_bases
        raise ValueError(f"Invalid bases in sequence: {invalid}. Only ATGC allowed.")
    
    if len(seq) < 2:
        return 0.0
    
    dg = 0.0
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        if dinuc not in NN_DELTA_G:
            raise ValueError(f"Invalid dinucleotide {dinuc}")
        dg += NN_DELTA_G[dinuc]
    return dg


# ============================================================
# 2. RAW ΔG MODELS (ONLY VALID MOTIFS)
# ============================================================

# Thermodynamic constants for cruciform structures
LOOP_PENALTY = 0.35  # kcal/mol per nucleotide - entropic penalty for loop formation
SUPERCOILING_FACTOR = 12.0  # kcal/mol - torsional stress relief from supercoiling

# Thermodynamic constants for slipped DNA structures
REPEAT_STABILIZATION = 0.8  # kcal/mol per repeat unit - cooperative stacking benefit
MISMATCH_PENALTY = 1.0  # kcal/mol per mismatch - hydrogen bonding disruption

# Thermodynamic constants for R-loop structures
RNA_DNA_STABILITY = -2.2  # kcal/mol per bp - RNA-DNA hybrid stability
DNA_DNA_STABILITY = -1.5  # kcal/mol per bp - DNA-DNA duplex stability
RLOOP_INITIATION = -1.5  # kcal/mol - initiation penalty for R-loop formation

def dg_cruciform(stem_seq: str, loop_len: int, sigma: float = -0.06) -> float:
    """
    Calculate ΔG for cruciform structures.
    
    Model includes:
    - Stem stability from nearest-neighbor parameters
    - Loop entropy penalty (0.35 kcal/mol per nt)
    - Supercoiling relief benefit (12.0 kcal/mol per bp under tension)
    
    Args:
        stem_seq: Stem sequence
        loop_len: Loop length
        sigma: Supercoiling density (default -0.06)
    
    Returns:
        Free energy in kcal/mol
    """
    return (
        nn_delta_g(stem_seq)
        + LOOP_PENALTY * loop_len
        - abs(sigma) * len(stem_seq) * SUPERCOILING_FACTOR
    )

def dg_slipped(stem_seq: str, repeat_units: int, mismatches: int) -> float:
    """
    Calculate ΔG for slipped DNA structures.
    
    Model includes:
    - Stem stability from nearest-neighbor parameters
    - Repeat unit stabilization (0.8 kcal/mol per unit)
    - Mismatch destabilization (1.0 kcal/mol per mismatch)
    
    Args:
        stem_seq: Stem sequence
        repeat_units: Number of repeat units
        mismatches: Number of mismatches
    
    Returns:
        Free energy in kcal/mol
    """
    return (
        nn_delta_g(stem_seq)
        - REPEAT_STABILIZATION * repeat_units
        + MISMATCH_PENALTY * mismatches
    )

def dg_rloop(length: int) -> float:
    """
    Calculate ΔG for R-loop structures.
    
    Model based on differential stability:
    - RNA-DNA hybrid (-2.2 kcal/mol per bp)
    - vs displaced DNA-DNA duplex (-1.5 kcal/mol per bp)
    - Plus initiation penalty (-1.5 kcal/mol)
    
    Net ΔG = length * (RNA-DNA stability - DNA-DNA stability) + initiation
           = length * (-2.2 - (-1.5)) + (-1.5)
           = length * (-0.7) - 1.5
    
    Args:
        length: R-loop length
    
    Returns:
        Free energy in kcal/mol
    """
    # Net stability difference: RNA-DNA hybrid is more stable than DNA-DNA
    # More negative ΔG means more stable
    net_stability_per_bp = RNA_DNA_STABILITY - DNA_DNA_STABILITY  # -2.2 - (-1.5) = -0.7
    return (length * net_stability_per_bp) + RLOOP_INITIATION


# ============================================================
# 3. PER-MOTIF NORMALIZATION (INDEPENDENT 1–3 SCALE)
# ============================================================

DG_RANGES = {
    "Cruciform": (-25.0, -5.0),
    "Slipped DNA": (-20.0, -3.0),
    "R-Loop": (-30.0, -6.0),
}

def normalize_1_3(dg_raw: float, motif: str) -> float:
    """
    Normalize ΔG to 1-3 scale independently per motif class.
    
    More stable (lower/more negative ΔG) → higher score (closer to 3.0)
    Less stable (higher/less negative ΔG) → lower score (closer to 1.0)
    
    Args:
        dg_raw: Raw ΔG value
        motif: Motif class name (must be in DG_RANGES)
    
    Returns:
        Score in range [1.0, 3.0]
        
    Raises:
        ValueError: If motif class not supported for ΔG-based scoring
    """
    if motif not in DG_RANGES:
        raise ValueError(
            f"Motif class '{motif}' not supported for ΔG-based scoring. "
            f"Supported classes: {list(DG_RANGES.keys())}"
        )
    
    lo, hi = DG_RANGES[motif]
    x = (dg_raw - lo) / (hi - lo)
    x = max(0.0, min(1.0, x))
    return round(1.0 + 2.0 * (1.0 - x), 3)


# Default scoring constants
DEFAULT_LEGACY_SCORE = 1.5  # Default score for non-ΔG motifs (middle of 1-3 range)
LARGE_SEQUENCE_THRESHOLD = 5_000_000  # 5M bp - threshold for explicit garbage collection


# ============================================================
# 4. SCORING DISPATCHER (ΔG SELECTIVE)
# ============================================================

def score_motif(m: Dict) -> Dict:
    """
    Apply appropriate scoring model to a motif.
    
    Args:
        m: Motif dictionary
    
    Returns:
        Motif dictionary with Score, ΔG_raw, and Score_Model fields
    """
    cls = m["Class"]

    if cls == "Cruciform":
        dg = dg_cruciform(m["Stem_Sequence"], m["Loop_Length"])
        m["ΔG_raw"] = dg
        m["Score"] = normalize_1_3(dg, cls)
        m["Score_Model"] = "ΔG_NN+Loop+Supercoil"
        return m

    if cls == "Slipped DNA":
        dg = dg_slipped(
            m["Stem_Sequence"],
            m["Repeat_Units"],
            m.get("Mismatches", 0)
        )
        m["ΔG_raw"] = dg
        m["Score"] = normalize_1_3(dg, cls)
        m["Score_Model"] = "ΔG_NN+Repeat"
        return m

    if cls == "R-Loop":
        dg = dg_rloop(m["Length"])
        m["ΔG_raw"] = dg
        m["Score"] = normalize_1_3(dg, cls)
        m["Score_Model"] = "ΔG_RNA-DNA"
        return m

    # NON-ΔG MOTIFS → KEEP LEGACY SCORE
    # Use existing score if present, otherwise default to middle of range
    m["Score"] = m.get("Score", DEFAULT_LEGACY_SCORE)
    m["Score_Model"] = "Heuristic_Legacy"
    return m


def score_all(motifs: List[Dict]) -> List[Dict]:
    """
    Score all motifs using appropriate models.
    
    Args:
        motifs: List of motif dictionaries
    
    Returns:
        List of scored motif dictionaries
    """
    return [score_motif(m) for m in motifs]


# ============================================================
# 5. SINGLE-PASS ANALYSIS (NO DUPLICATES)
# ============================================================

def analyze_sequence_single_pass(
    seq: str,
    name: str,
    use_parallel: bool,
    parallel_scanner=None
) -> List[Dict]:
    """
    Analyze sequence with exactly ONE scan pass.
    
    HARD RULE:
    - Either parallel scanner OR standard analyzer
    - NEVER both
    
    Args:
        seq: DNA sequence
        name: Sequence name
        use_parallel: Use parallel scanner if True
        parallel_scanner: ParallelScanner instance with run_scan() method (optional)
    
    Returns:
        List of motif dictionaries
        
    Raises:
        AttributeError: If parallel_scanner doesn't have required interface
    """
    if use_parallel and parallel_scanner is not None:
        # Validate parallel_scanner has required interface
        if not hasattr(parallel_scanner, 'run_scan'):
            raise AttributeError(
                "parallel_scanner must have a 'run_scan()' method. "
                "Expected interface: run_scan(seq: str) -> List[Dict]"
            )
        raw = parallel_scanner.run_scan(seq)
        return raw
    else:
        from nonbscanner import analyze_sequence
        return analyze_sequence(seq, name)


# ============================================================
# 6. PER-SEQUENCE CACHE KEY (CONFIG-AWARE)
# ============================================================

def seq_cache_key(seq: str, config: Dict) -> str:
    """
    Generate cache key for sequence + configuration.
    
    Args:
        seq: DNA sequence
        config: Configuration dictionary
    
    Returns:
        SHA256 hash as hex string
    """
    h = hashlib.sha256()
    h.update(seq.encode())
    h.update(str(sorted(config.items())).encode())
    return h.hexdigest()


# ============================================================
# 7. NUMPY GENOME CACHE
# ============================================================

def genome_to_numpy(seq: str) -> np.ndarray:
    """
    Convert genome sequence to NumPy array for efficient processing.
    
    Args:
        seq: DNA sequence
    
    Returns:
        NumPy uint8 array
    """
    return np.frombuffer(seq.encode("ascii"), dtype=np.uint8)


# ============================================================
# 8. LAZY VISUALIZATION (MANDATORY plt.close)
# ============================================================

def plot_when_visible(fig_fn, *args, **kwargs):
    """
    Generate plot lazily and ensure cleanup.
    
    IMPORTANT: This function closes the figure before returning it.
    The returned figure is already closed and should not be used for further operations.
    This is intentional to prevent memory leaks from matplotlib figures.
    
    Usage pattern:
        # Generate and display immediately
        fig = plot_when_visible(my_plot_function, data)
        st.pyplot(fig)  # streamlit will handle the closed figure correctly
        # Don't use fig after this point
    
    Args:
        fig_fn: Function that returns matplotlib figure
        *args: Arguments to fig_fn
        **kwargs: Keyword arguments to fig_fn
    
    Returns:
        Matplotlib figure (already closed for memory efficiency)
        Note: The figure is closed but can still be passed to display functions
    """
    fig = fig_fn(*args, **kwargs)
    plt.close(fig)
    return fig


# ============================================================
# 9. STREAMING EXCEL EXPORT
# ============================================================

def export_excel_stream(df: pd.DataFrame, path: str):
    """
    Export DataFrame to Excel with constant memory usage.
    
    Args:
        df: Pandas DataFrame to export
        path: Output file path
    """
    with pd.ExcelWriter(
        path,
        engine="xlsxwriter",
        options={"constant_memory": True}
    ) as writer:
        df.to_excel(writer, index=False, sheet_name="Motifs")


# ============================================================
# 10. MEMORY HYGIENE
# ============================================================

def cleanup_after_sequence(seq_len: int):
    """
    Perform garbage collection for large sequences.
    
    Explicitly triggers Python's garbage collector to release memory
    after processing sequences larger than LARGE_SEQUENCE_THRESHOLD.
    This helps manage memory usage in Streamlit's 1GB limit.
    
    Args:
        seq_len: Sequence length in bp
    """
    if seq_len > LARGE_SEQUENCE_THRESHOLD:
        gc.collect()
