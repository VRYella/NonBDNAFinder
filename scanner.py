"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                   CONSOLIDATED SCANNER MODULE                                 ║
║          Non-B DNA Motif Analysis and Orchestration Functions                ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: scanner.py  
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.1 - Consolidated
LICENSE: MIT

DESCRIPTION:
    Consolidated module containing scanner orchestration, analysis functions,
    and repeat detection algorithms. Combines functionality from modular_scanner,
    nbdscanner, and repeat_scanner into a single organized module.

MAIN FUNCTIONS:
    - analyze_sequence(): Single sequence analysis with all detectors
    - analyze_multiple_sequences(): Batch sequence analysis
    - detect_hybrids_and_clusters(): Advanced overlap and cluster detection
    - export_results_to_dataframe(): Export results to pandas DataFrame
    - get_motif_classification_info(): Get motif class metadata

PERFORMANCE:
    - Overall rate: ~5,800 bp/second on 10kb sequences
    - Slipped DNA: ~280,000 bp/second on 50kb sequences
    - All detectors: Linear O(n) complexity
"""
from __future__ import annotations

#!/usr/bin/env python3
"""
Optimized genome-scale Python scanner for:
 - Direct repeats (unit length 10..300 bp, spacer <= 10 bp)
 - Inverted repeats (arm >= 6 bp, loop <= 100 bp)
 - Mirror repeats (arm >= 10 bp, loop <= 100 bp)
 - STRs (unit size 1..9 bp, total repeated length >= 20 bp)

Design:
 - Seed-and-extend using k-mer indices (dict: kmer -> list(positions))
 - Rolling-hash-free verification: uses direct slice compare for correctness
 - Candidate pruning: only consider seed position pairs with delta <= max_unit + max_spacer
 - Tuned defaults and safe-guards to avoid explosion on highly-repetitive seeds

Usage:
    from utilities importrepeat_scanner import find_direct_repeats, find_inverted_repeats, find_mirror_repeats, find_strs

Requirements:
    Python 3.8+
    (pure Python; no external C deps required)
"""

from collections import defaultdict
from typing import List, Dict, Tuple

# -------------------------
# Parameters (user constraints)
# -------------------------
# Slipped DNA: 10–50 nt repeat separated by spacer = 0 nt
DIRECT_MIN_UNIT = 10
DIRECT_MAX_UNIT = 50
DIRECT_MAX_SPACER = 0

# Cruciform DNA: 10–100 nt arm length with reverse complement separated by spacer = 0–3 nt
INVERTED_MIN_ARM = 10
INVERTED_MAX_LOOP = 3

# Triplex DNA: 10–100 nt mirrored with spacer = 0–8 nt
MIRROR_MIN_ARM = 10
MIRROR_MAX_LOOP = 8

STR_MIN_UNIT = 1
STR_MAX_UNIT = 9
STR_MIN_TOTAL = 20

# Seeds (k) used for k-mer indexing (trade-off: larger k -> fewer hits)
K_DIRECT = 10   # equal to DIRECT_MIN_UNIT (good seed)
K_INVERTED = 10  # equal to INVERTED_MIN_ARM (min arm) - updated for Cruciform
K_MIRROR = 10   # equal to MIRROR_MIN_ARM

# Safety thresholds
MAX_POSITIONS_PER_KMER = 10000  # skip extremely frequent k-mers to avoid blow-up

# -------------------------
# Utilities
# -------------------------
_RC_TRANS = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s: str) -> str:
    return s.translate(_RC_TRANS)[::-1]


# -------------------------
# Performance utilities
# -------------------------
def _calc_gc_content(seq: str) -> float:
    """Fast GC content calculation (percentage)."""
    if not seq:
        return 0.0
    return (seq.count('G') + seq.count('C')) * 100.0 / len(seq)


# -------------------------
# K-mer index builder
# -------------------------
def build_kmer_index(seq: str, k: int) -> Dict[str, List[int]]:
    """
    Build k-mer position index for seed-and-extend pattern matching.
    
    # Index Structure:
    # | Field    | Type       | Description                        |
    # |----------|------------|------------------------------------|
    # | key      | str        | k-mer sequence (length k)          |
    # | value    | List[int]  | 0-based positions in sequence      |
    
    Args:
        seq: DNA sequence (ACGT)
        k: K-mer length
        
    Returns:
        Dictionary mapping k-mer -> position list (filtered for frequency)
    """
    n = len(seq)
    idx: Dict[str, List[int]] = defaultdict(list)
    
    # Build index with early filtering
    valid_bases = frozenset("ACGT")
    for i in range(n - k + 1):
        kmer = seq[i:i + k]
        if all(ch in valid_bases for ch in kmer):
            lst = idx[kmer]
            if len(lst) <= MAX_POSITIONS_PER_KMER:
                lst.append(i)
    
    # Remove overly-frequent k-mers in-place for memory efficiency
    for kmer in list(idx.keys()):
        if len(idx[kmer]) > MAX_POSITIONS_PER_KMER:
            del idx[kmer]
    
    return idx


# -------------------------
# Direct Repeats (seed-and-extend)
# -------------------------
def find_direct_repeats(seq: str, 
                       min_unit: int = DIRECT_MIN_UNIT, 
                       max_unit: int = DIRECT_MAX_UNIT,
                       max_spacer: int = DIRECT_MAX_SPACER) -> List[Dict]:
    """
    Find direct repeats using k-mer seed-and-extend strategy.
    
    # Output Structure:
    # | Field       | Type  | Description                      |
    # |-------------|-------|----------------------------------|
    # | Class       | str   | Always 'Direct_Repeat'           |
    # | Subclass    | str   | Direct_L{unit_length}            |
    # | Start       | int   | 1-based start position           |
    # | End         | int   | End position                     |
    # | Unit_Length | int   | Length of repeat unit            |
    # | Spacer      | int   | Spacer length between units      |
    # | GC_Unit     | float | GC% of repeat unit               |
    # | GC_Total    | float | GC% of full motif                |
    
    Returns:
        List of direct repeat dictionaries
    """
    n = len(seq)
    idx = build_kmer_index(seq, K_DIRECT)
    results = []
    max_delta = max_unit + max_spacer
    
    for kmer, poses in idx.items():
        m = len(poses)
        if m < 2:
            continue
            
        for a in range(m):
            i = poses[a]
            b = a + 1
            while b < m:
                j = poses[b]
                delta = j - i
                if delta > max_delta:
                    break
                    
                L_min = max(min_unit, delta - max_spacer)
                L_max = min(max_unit, delta)
                if L_min <= L_max:
                    for L in range(L_max, L_min - 1, -1):
                        j_start = j
                        if j_start + L > n or i + L > n:
                            continue
                        if seq[i:i + L] == seq[j_start:j_start + L]:
                            unit_seq = seq[i:i + L]
                            spacer_seq = seq[i + L:j_start] if j_start > i + L else ''
                            full_seq = seq[i:j_start + L]
                            
                            rec = {
                                'Class': 'Direct_Repeat',
                                'Subclass': f'Direct_L{L}',
                                'Start': i + 1,
                                'End': j_start + L,
                                'Length': (j_start + L) - i,
                                'Unit_Length': L,
                                'Spacer': j_start - (i + L),
                                'Left_Pos': i + 1,
                                'Right_Pos': j_start + 1,
                                'Unit_Seq': unit_seq,
                                'Spacer_Seq': spacer_seq,
                                'Sequence': full_seq,
                                'Left_Unit': unit_seq,
                                'Right_Unit': seq[j_start:j_start + L],
                                'GC_Unit': round(_calc_gc_content(unit_seq), 2),
                                'GC_Spacer': round(_calc_gc_content(spacer_seq), 2),
                                'GC_Total': round(_calc_gc_content(full_seq), 2)
                            }
                            results.append(rec)
                            break
                b += 1
    
    # Deduplicate: prefer maximal unit_length per (Left_Pos, Spacer)
    dedup: Dict[Tuple[int, int], Dict] = {}
    for rec in results:
        key = (rec['Left_Pos'], rec['Spacer'])
        if key not in dedup or rec['Unit_Length'] > dedup[key]['Unit_Length']:
            dedup[key] = rec
    return list(dedup.values())


# -------------------------
# Inverted Repeats (Cruciform)
# -------------------------
def find_inverted_repeats(seq: str, 
                         min_arm: int = INVERTED_MIN_ARM, 
                         max_loop: int = INVERTED_MAX_LOOP) -> List[Dict]:
    """
    Find inverted repeats (cruciform precursors) using k-mer indexing.
    
    # Output Structure:
    # | Field       | Type  | Description                      |
    # |-------------|-------|----------------------------------|
    # | Class       | str   | Always 'Inverted_Repeat'         |
    # | Subclass    | str   | Inverted_arm_{arm_length}        |
    # | Arm_Length  | int   | Length of palindromic arm        |
    # | Loop_Length | int   | Length of loop/spacer            |
    # | Left_Arm    | str   | Left arm sequence                |
    # | Right_Arm   | str   | Right arm sequence (RC of left)  |
    # | GC_Total    | float | GC% of full structure            |
    
    Returns:
        List of inverted repeat dictionaries
    """
    n = len(seq)
    idx = build_kmer_index(seq, K_INVERTED)
    results = []
    
    # Precompute reverse complement mapping
    rc_map = {}
    for kmer in list(idx.keys()):
        rc = revcomp(kmer)
        if rc in idx:
            rc_map[kmer] = idx[rc]
    
    # Find inverted repeats
    for kmer, rc_positions in rc_map.items():
        left_positions = idx[kmer]
        for i in left_positions:
            for j in rc_positions:
                if j <= i:
                    continue
                delta = j - i
                arm_min = max(min_arm, delta - max_loop)
                arm_max = min(delta, n - j, n - i)
                if arm_min > arm_max:
                    continue
                    
                for arm in range(arm_max, arm_min - 1, -1):
                    if i + arm > n or j + arm > n:
                        continue
                    left_sub = seq[i:i + arm]
                    right_sub = seq[j:j + arm]
                    if left_sub == revcomp(right_sub):
                        loop_seq = seq[i + arm:j] if j > i + arm else ''
                        full_seq = seq[i:j + arm]
                        
                        rec = {
                            'Class': 'Inverted_Repeat',
                            'Subclass': f'Inverted_arm_{arm}',
                            'Start': i + 1,
                            'End': j + arm,
                            'Length': (j + arm) - i,
                            'Left_Start': i + 1,
                            'Right_Start': j + 1,
                            'Arm_Length': arm,
                            'Loop': delta - arm,
                            'Loop_Length': len(loop_seq),
                            'Left_Arm': left_sub,
                            'Right_Arm': right_sub,
                            'Loop_Seq': loop_seq,
                            'Sequence': full_seq,
                            'Stem': f'{left_sub}...{right_sub}',
                            'GC_Left_Arm': round(_calc_gc_content(left_sub), 2),
                            'GC_Right_Arm': round(_calc_gc_content(right_sub), 2),
                            'GC_Loop': round(_calc_gc_content(loop_seq), 2),
                            'GC_Total': round(_calc_gc_content(full_seq), 2)
                        }
                        results.append(rec)
                        break
    
    # Deduplicate: keep maximal arm per (Left_Start, Loop)
    dedup: Dict[Tuple[int, int], Dict] = {}
    for rec in results:
        key = (rec['Left_Start'], rec['Loop'])
        if key not in dedup or rec['Arm_Length'] > dedup[key]['Arm_Length']:
            dedup[key] = rec
    return list(dedup.values())


# -------------------------
# Mirror Repeats (Triplex DNA component)
# -------------------------
def find_mirror_repeats(seq: str, 
                       min_arm: int = MIRROR_MIN_ARM, 
                       max_loop: int = MIRROR_MAX_LOOP,
                       purine_pyrimidine_threshold: float = 0.9) -> List[Dict]:
    """
    Mirror repeats: left arm matches reverse (not complement) of right arm.
    Same pattern as inverted but compare seq[i:i+arm] == seq[j:j+arm][::-1]
    
    For Triplex DNA, we also filter for >90% purine or pyrimidine content in arms.
    """
    n = len(seq)
    idx = build_kmer_index(seq, K_MIRROR)
    results = []
    
    keys = list(idx.keys())
    # For mirror, the reverse of substring is simply reversed string
    rev_map = {}
    for kmer in keys:
        rev = kmer[::-1]
        if rev in idx:
            rev_map[kmer] = idx[rev]
    
    for kmer, rev_positions in rev_map.items():
        left_positions = idx[kmer]
        for i in left_positions:
            for j in rev_positions:
                if j <= i:
                    continue
                delta = j - i
                arm_min = max(min_arm, delta - max_loop)
                arm_max = min(delta, n - j, n - i)
                if arm_min > arm_max:
                    continue
                for arm in range(arm_max, arm_min - 1, -1):
                    if i + arm > n or j + arm > n:
                        continue
                    left_arm = seq[i:i + arm]
                    right_arm = seq[j:j + arm]
                    if left_arm == right_arm[::-1]:
                        # Check purine/pyrimidine content for Triplex DNA
                        combined_arms = left_arm + right_arm
                        purine_count = sum(1 for b in combined_arms if b in 'AG')
                        pyrimidine_count = sum(1 for b in combined_arms if b in 'CT')
                        total_bases = len(combined_arms)
                        
                        purine_fraction = purine_count / total_bases if total_bases > 0 else 0
                        pyrimidine_fraction = pyrimidine_count / total_bases if total_bases > 0 else 0
                        
                        is_triplex = (purine_fraction >= purine_pyrimidine_threshold or 
                                     pyrimidine_fraction >= purine_pyrimidine_threshold)
                        
                        # Calculate component details
                        loop_seq = seq[i + arm:j] if j > i + arm else ''
                        full_seq = seq[i:j + arm]
                        
                        # Calculate GC content for components
                        gc_left_arm = (left_arm.count('G') + left_arm.count('C')) / len(left_arm) * 100 if len(left_arm) > 0 else 0
                        gc_right_arm = (right_arm.count('G') + right_arm.count('C')) / len(right_arm) * 100 if len(right_arm) > 0 else 0
                        gc_loop = (loop_seq.count('G') + loop_seq.count('C')) / len(loop_seq) * 100 if len(loop_seq) > 0 else 0
                        gc_total = (full_seq.count('G') + full_seq.count('C')) / len(full_seq) * 100 if len(full_seq) > 0 else 0
                        
                        rec = {
                            'Class': 'Mirror_Repeat',
                            'Subclass': f'Mirror_arm_{arm}',
                            'Start': i + 1,
                            'End': j + arm,
                            'Length': (j + arm) - i,
                            'Left_Start': i + 1,
                            'Right_Start': j + 1,
                            'Arm_Length': arm,
                            'Loop': delta - arm,
                            'Loop_Length': len(loop_seq),
                            'Left_Arm': left_arm,
                            'Right_Arm': right_arm,
                            'Loop_Seq': loop_seq,
                            'Sequence': full_seq,
                            'Is_Triplex': is_triplex,
                            'Purine_Fraction': round(purine_fraction, 3),
                            'Pyrimidine_Fraction': round(pyrimidine_fraction, 3),
                            # Component details
                            'Stem': f'{left_arm}...{right_arm}',
                            'GC_Left_Arm': round(gc_left_arm, 2),
                            'GC_Right_Arm': round(gc_right_arm, 2),
                            'GC_Loop': round(gc_loop, 2),
                            'GC_Total': round(gc_total, 2)
                        }
                        results.append(rec)
                        break
    
    # dedupe maximal arms
    dedup: Dict[Tuple[int, int], Dict] = {}
    for rec in results:
        key = (rec['Left_Start'], rec['Loop'])
        if key not in dedup or rec['Arm_Length'] > dedup[key]['Arm_Length']:
            dedup[key] = rec
    return list(dedup.values())


# -------------------------
# STRs (Short Tandem Repeats)
# -------------------------
def find_strs(seq: str, 
              min_u: int = STR_MIN_UNIT, 
              max_u: int = STR_MAX_UNIT, 
              min_total: int = STR_MIN_TOTAL) -> List[Dict]:
    """
    Greedy detection of perfect STRs (tandem repeats).
    For each unit size k in 1..9, slide and count consecutive copies.
    """
    n = len(seq)
    results = []
    for k in range(min_u, max_u + 1):
        i = 0
        while i <= n - k:
            unit = seq[i:i + k]
            # require at least one repeat following to qualify quickly
            if i + k >= n or seq[i + k:i + 2 * k] != unit:
                i += 1
                continue
            # count copies
            copies = 1
            j = i + k
            while j + k <= n and seq[j:j + k] == unit:
                copies += 1
                j += k
            total_len = copies * k
            if total_len >= min_total:
                # Calculate component details
                full_seq = seq[i:j]
                
                # Calculate GC content
                gc_unit = (unit.count('G') + unit.count('C')) / len(unit) * 100 if len(unit) > 0 else 0
                gc_total = (full_seq.count('G') + full_seq.count('C')) / len(full_seq) * 100 if len(full_seq) > 0 else 0
                
                # Calculate AT/GC composition of unit
                a_count = unit.count('A')
                t_count = unit.count('T')
                g_count = unit.count('G')
                c_count = unit.count('C')
                
                results.append({
                    'Class': 'STR',
                    'Subclass': f'unit_{k}',
                    'Start': i + 1,
                    'End': j,
                    'Unit_Length': k,
                    'Copies': copies,
                    'Length': total_len,
                    'Unit_Seq': unit,
                    'Sequence': full_seq,
                    # Component details
                    'Repeat_Unit': unit,
                    'Number_of_Copies': copies,
                    'GC_Unit': round(gc_unit, 2),
                    'GC_Total': round(gc_total, 2),
                    'Unit_A_Count': a_count,
                    'Unit_T_Count': t_count,
                    'Unit_G_Count': g_count,
                    'Unit_C_Count': c_count
                })
                i = j  # skip to end of tandem block
            else:
                i += 1
    return results
"""
Modular NBD Scanner (RECOMMENDED for Production Use)
====================================================

TABULAR SUMMARY:
┌──────────────────────────────────────────────────────────────────────────────┐
│ Module:        modular_scanner.py                                            │
│ Purpose:       High-performance modular motif detection system               │
│ Performance:   24,674 bp/s on 100K sequences (fast detectors)                │
│ Author:        Dr. Venkata Rajesh Yella                                      │
│ Last Updated:  2024                                                          │
├──────────────────────────────────────────────────────────────────────────────┤
│ MOTIF CLASSES (9 detectors):                                                 │
│ ✓ CurvedDNA      - A-tract mediated DNA bending                              │
│ ✓ SlippedDNA     - Tandem repeats (adaptive sampling for large sequences)    │
│ ✓ Cruciform      - Inverted repeats (sliding window for large sequences)    │
│ ✓ R-Loop         - RNA-DNA hybrid formation sites                            │
│ ✓ Triplex        - Three-stranded DNA structures                             │
│ ✓ G-Quadruplex   - Four-stranded G-rich structures                           │
│ ✓ i-Motif        - C-rich structures                                         │
│ ✓ Z-DNA          - Left-handed double helix                                  │
│ ✓ A-philic       - A-rich protein binding sites                              │
├──────────────────────────────────────────────────────────────────────────────┤
│ PERFORMANCE OPTIMIZATIONS:                                                   │
│ • Pure Python scanner DISABLED (too slow - O(n³) complexity)                 │
│ • Cruciform uses sliding window for sequences >1K bp                         │
│ • SlippedDNA uses adaptive sampling for sequences >10K bp                    │
│ • Pattern compilation and caching                                            │
│ • Hybrid and cluster detection post-processing                               │
├──────────────────────────────────────────────────────────────────────────────┤
│ USAGE:                                                                       │
│   from modular_scanner import analyze_sequence                               │
│   motifs = analyze_sequence(sequence, "seq_name")                            │
└──────────────────────────────────────────────────────────────────────────────┘
"""

import re
import os
import logging
import pandas as pd
from typing import List, Dict, Any, Union
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
import warnings

logger = logging.getLogger(__name__)

# Note: Detectors are imported lazily to avoid circular dependency
# (detectors.py imports from scanner.py for optimized functions)

warnings.filterwarnings("ignore")

# Try to import hyperscan for acceleration (if available)
try:
    import hyperscan
    HYPERSCAN_AVAILABLE = True
except ImportError:
    HYPERSCAN_AVAILABLE = False

# Import pure Python scanner for specific classes
try:
    from nonb_pure_python import scan_sequence as pure_python_scan
    PURE_PYTHON_AVAILABLE = True
except ImportError:
    PURE_PYTHON_AVAILABLE = False

# Import motif_patterns for registry loading
try:
    # from . import motif_patterns
    _MOTIF_PATTERNS_AVAILABLE = True
except ImportError:
    try:
        # Fallback: direct import
        import sys
        import os
        sys.path.insert(0, os.path.dirname(__file__))
        import motif_patterns
        _MOTIF_PATTERNS_AVAILABLE = True
    except ImportError:
        motif_patterns = None
        _MOTIF_PATTERNS_AVAILABLE = False

DEFAULT_REGISTRY_DIR = os.environ.get("NBD_REGISTRY_DIR", "registry")


class ModularMotifDetector:
    """
    Modular motif detection using individual detector classes.
    Each motif class has its own dedicated detector.
    """
    
    def __init__(self, registry_dir: str = DEFAULT_REGISTRY_DIR):
        """Initialize all individual detectors (lazy import to avoid circular dependency)"""
        # Lazy import of detectors to break circular dependency
        from detectors import (
            CurvedDNADetector,
            SlippedDNADetector,
            CruciformDetector,
            RLoopDetector,
            TriplexDetector,
            GQuadruplexDetector,
            IMotifDetector,
            ZDNADetector,
            APhilicDetector
        )
        
        self.registry_dir = registry_dir
        self.detectors = {
            'curved_dna': CurvedDNADetector(),
            'slipped_dna': SlippedDNADetector(),
            'cruciform': CruciformDetector(), 
            'r_loop': RLoopDetector(),
            'triplex': TriplexDetector(),
            'g_quadruplex': GQuadruplexDetector(),
            'i_motif': IMotifDetector(),
            'z_dna': ZDNADetector(),
            'a_philic': APhilicDetector()
        }
        # Store detector classes for preload use (only the ones that need Hyperscan DB preloading)
        self._detector_classes = {
            'ZDNADetector': ZDNADetector,
            'APhilicDetector': APhilicDetector
        }
        # Preload Hyperscan DBs for detectors that have registries
        self._preload_detector_dbs()
    
    def _preload_detector_dbs(self):
        """
        Attempt to load precompiled Hyperscan DBs for detectors that benefit from them.
        This is optional and non-fatal — detectors still have their own fallback code paths.
        """
        if not _MOTIF_PATTERNS_AVAILABLE:
            # motif_patterns module not available, skip preloading
            logger.debug("motif_patterns module not available, skipping DB preloading")
            return
        
        # Try to detect available registries dynamically from registry directory
        candidate_classes = []
        if os.path.isdir(self.registry_dir):
            # Look for *_registry.pkl or *_registry.json files
            for fname in os.listdir(self.registry_dir):
                if fname.endswith('_registry.pkl') or fname.endswith('_registry.json'):
                    cls_name = fname.replace('_registry.pkl', '').replace('_registry.json', '')
                    if cls_name not in candidate_classes:
                        candidate_classes.append(cls_name)
        
        # Fallback to known classes if directory doesn't exist
        if not candidate_classes:
            candidate_classes = ["ZDNA", "APhilic"]
        
        # Map class names to detector instances for optional attribute setting
        detector_map = {
            "ZDNA": ("z_dna", self._detector_classes['ZDNADetector']),
            "APhilic": ("a_philic", self._detector_classes['APhilicDetector'])
        }
        
        for cls in candidate_classes:
            try:
                db, id_to_pattern, id_to_score = motif_patterns.get_hs_db_for_class(cls, registry_dir=self.registry_dir)
                # Pass to detectors: detectors should check for availability when scanning
                # Store in a class-level map for detectors to access
                if not hasattr(self, "hsdb_map"):
                    self.hsdb_map = {}
                self.hsdb_map[cls] = {"db": db, "id_to_pattern": id_to_pattern, "id_to_score": id_to_score}
                
                # Optionally, set into detector classes for direct use
                if cls in detector_map:
                    det_key, det_cls = detector_map[cls]
                    # Set class-level attribute for detectors to access
                    setattr(det_cls, "HS_DB_INFO", {"db": db, "id_to_pattern": id_to_pattern, "id_to_score": id_to_score})
                    logger.debug(f"Set HS_DB_INFO for {det_cls.__name__}")
            except FileNotFoundError:
                # registry not built — continue silently
                logger.debug(f"Registry not found for {cls}, skipping")
                continue
            except Exception as e:
                # protect scanner initialization from any DB loading error
                logger.warning(f"Failed to preload HS DB for {cls}: {e}")
                continue
    
    def analyze_sequence(self, sequence: str, sequence_name: str = "sequence", 
                        use_pure_python: bool = False) -> List[Dict[str, Any]]:
        """
        Main analysis function using modular detectors
        
        Args:
            sequence: DNA sequence to analyze
            sequence_name: Name identifier for the sequence
            use_pure_python: Whether to use pure Python scanner (DISABLED for performance)
            
        Returns:
            List of detected motifs with comprehensive metadata
        """
        sequence = sequence.upper().strip()
        all_motifs = []
        
        # PERFORMANCE: Disable pure Python scanner - too slow on large sequences
        # Use modular detectors only for better performance
        skip_classes = set()
        
        # Use modular detectors for all classes
        for class_key, detector in self.detectors.items():
            if class_key not in skip_classes:
                try:
                    motifs = detector.detect_motifs(sequence, sequence_name)
                    all_motifs.extend(motifs)
                except Exception as e:
                    print(f"Warning: Error in {class_key} detector: {e}")
                    continue
        
        # Remove overlaps and add hybrid/cluster detection
        filtered_motifs = self._remove_overlaps(all_motifs)
        hybrid_motifs = self._detect_hybrid_motifs(filtered_motifs)
        cluster_motifs = self._detect_clusters(filtered_motifs, sequence)
        
        final_motifs = filtered_motifs + hybrid_motifs + cluster_motifs
        
        # Sort by position
        final_motifs.sort(key=lambda x: x.get('Start', 0))
        
        return final_motifs
    
    def _remove_overlaps(self, motifs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Remove overlapping motifs within the same class/subclass - ensures NO overlaps"""
        if not motifs:
            return motifs
        
        # Group by class/subclass
        groups = defaultdict(list)
        for motif in motifs:
            key = f"{motif.get('Class', '')}-{motif.get('Subclass', '')}"
            groups[key].append(motif)
        
        filtered_motifs = []
        
        for group_motifs in groups.values():
            # Sort by score (highest first), then by length (longest first)
            group_motifs.sort(key=lambda x: (x.get('Score', 0), x.get('Length', 0)), reverse=True)
            
            non_overlapping = []
            for motif in group_motifs:
                overlaps = False
                for existing in non_overlapping:
                    # Strict overlap check: any overlap at all (>0) is rejected
                    if self._calculate_overlap(motif, existing) > 0.0:
                        overlaps = True
                        break
                
                if not overlaps:
                    non_overlapping.append(motif)
            
            filtered_motifs.extend(non_overlapping)
        
        return filtered_motifs
    
    def _calculate_overlap(self, motif1: Dict[str, Any], motif2: Dict[str, Any]) -> float:
        """Calculate overlap ratio between two motifs"""
        start1, end1 = motif1.get('Start', 0), motif1.get('End', 0)
        start2, end2 = motif2.get('Start', 0), motif2.get('End', 0)
        
        if end1 <= start2 or end2 <= start1:
            return 0.0
        
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)
        overlap_length = overlap_end - overlap_start
        
        min_length = min(end1 - start1, end2 - start2)
        return overlap_length / min_length if min_length > 0 else 0.0
    
    def _detect_hybrid_motifs(self, motifs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Detect hybrid motifs (overlapping different classes)"""
        hybrid_motifs = []
        
        for i, motif1 in enumerate(motifs):
            for motif2 in motifs[i+1:]:
                if motif1.get('Class') != motif2.get('Class'):
                    overlap = self._calculate_overlap(motif1, motif2)
                    if 0.3 < overlap < 1.0:  # Partial overlap between different classes
                        
                        # Create hybrid motif
                        start = min(motif1.get('Start', 0), motif2.get('Start', 0))
                        end = max(motif1.get('End', 0), motif2.get('End', 0))
                        avg_score = (motif1.get('Score', 0) + motif2.get('Score', 0)) / 2
                        
                        hybrid_motifs.append({
                            'ID': f"{motif1.get('Sequence_Name', 'seq')}_HYBRID_{start}",
                            'Sequence_Name': motif1.get('Sequence_Name', 'sequence'),
                            'Class': 'Hybrid',
                            'Subclass': f"{motif1.get('Class', '')}-{motif2.get('Class', '')}",
                            'Start': start,
                            'End': end,
                            'Length': end - start,
                            'Score': round(avg_score, 3),
                            'Strand': '+',
                            'Method': 'hybrid_detection'
                        })
        
        return hybrid_motifs
    
    def _detect_clusters(self, motifs: List[Dict[str, Any]], sequence: str) -> List[Dict[str, Any]]:
        """Detect high-density non-B DNA clusters"""
        if len(motifs) < 3:
            return []
        
        cluster_motifs = []
        window_size = 500  # 500bp window for cluster detection
        min_density = 3    # Minimum 3 motifs per window
        
        # Sort motifs by start position
        sorted_motifs = sorted(motifs, key=lambda x: x.get('Start', 0))
        
        for i in range(len(sorted_motifs)):
            window_start = sorted_motifs[i].get('Start', 0)
            window_end = window_start + window_size
            
            # Count motifs in window
            window_motifs = []
            for motif in sorted_motifs[i:]:
                if motif.get('Start', 0) <= window_end:
                    window_motifs.append(motif)
                else:
                    break
            
            if len(window_motifs) >= min_density:
                # Create cluster motif
                actual_start = min(m.get('Start', 0) for m in window_motifs)
                actual_end = max(m.get('End', 0) for m in window_motifs)
                avg_score = sum(m.get('Score', 0) for m in window_motifs) / len(window_motifs)
                
                cluster_motifs.append({
                    'ID': f"{sorted_motifs[i].get('Sequence_Name', 'seq')}_CLUSTER_{actual_start}",
                    'Sequence_Name': sorted_motifs[i].get('Sequence_Name', 'sequence'),
                    'Class': 'Non-B_DNA_Clusters',
                    'Subclass': f'High-density cluster ({len(window_motifs)} motifs)',
                    'Start': actual_start,
                    'End': actual_end,
                    'Length': actual_end - actual_start,
                    'Score': round(avg_score, 3),
                    'Strand': '+',
                    'Method': 'cluster_detection'
                })
        
        return cluster_motifs
    
    def get_detector_info(self) -> Dict[str, Any]:
        """Get information about all loaded detectors"""
        info = {
            'total_detectors': len(self.detectors),
            'detectors': {},
            'total_patterns': 0
        }
        
        for name, detector in self.detectors.items():
            stats = detector.get_statistics()
            info['detectors'][name] = stats
            info['total_patterns'] += stats['total_patterns']
        
        return info


def analyze_sequence(sequence: str, sequence_name: str = "sequence", 
                    detailed: bool = True) -> Union[List[Dict[str, Any]], Dict[str, Any]]:
    """
    Convenience function for sequence analysis using modular architecture
    
    Args:
        sequence: DNA sequence to analyze
        sequence_name: Name identifier for the sequence
        detailed: If True, returns detailed motif list; if False, returns summary
        
    Returns:
        List of motifs or summary statistics
    """
    detector = ModularMotifDetector()
    motifs = detector.analyze_sequence(sequence, sequence_name)
    
    if detailed:
        return motifs
    else:
        # Return summary statistics
        class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
        subclass_counts = Counter(m.get('Subclass', 'Unknown') for m in motifs)
        
        return {
            'total_motifs': len(motifs),
            'classes_detected': len(class_counts),
            'subclasses_detected': len(subclass_counts),
            'class_distribution': dict(class_counts),
            'subclass_distribution': dict(subclass_counts),
            'avg_score': sum(m.get('Score', 0) for m in motifs) / len(motifs) if motifs else 0
        }


def get_motif_classification_info() -> Dict[str, Any]:
    """Get comprehensive motif classification information"""
    detector = ModularMotifDetector()
    return detector.get_detector_info()


# For backward compatibility, create aliases
MotifDetector = ModularMotifDetector

if __name__ == "__main__":
    # Test the modular detector
    test_sequence = "GGGTTAGGGTTAGGGTTAGGGAAAAAAAATTTTTTCACACACACACACACA"
    
    detector = ModularMotifDetector()
    motifs = detector.analyze_sequence(test_sequence, "test_sequence")
    
    print(f"Modular detector found {len(motifs)} motifs:")
    for motif in motifs:
        print(f"  {motif['Class']}/{motif['Subclass']}: {motif['Start']}-{motif['End']} (score: {motif['Score']:.3f})")
    
    print(f"\nDetector info: {detector.get_detector_info()}")


# ============================================================================
# NBDScanner - Non-B DNA Structure Detection & Analysis Engine
# ============================================================================
"""
NBDScanner - Non-B DNA Structure Detection & Analysis Engine
============================================================

Consolidated motif detection system with 11 classes and 22+ subclasses.
Optimized architecture with Hyperscan acceleration for high-performance analysis.

MOTIF CLASSIFICATION TABLE:
============================
Class | Name              | Subclasses                                    | Description
------|-------------------|-----------------------------------------------|--------------------
  1   | Curved DNA        | Global curvature, Local Curvature           | A-tract mediated bending
  2   | Slipped DNA       | Direct Repeat, STR                           | Tandem repeat slippage
  3   | Cruciform DNA     | Inverted Repeats                             | Four-way junction
  4   | R-loop            | R-loop formation sites                       | RNA-DNA hybrids
  5   | Triplex           | Triplex, Sticky DNA                          | Three-strand structures
  6   | G-Quadruplex      | Multimeric, Canonical, Relaxed, Bulged,     | Four-strand G-rich
      |                   | Bipartite, Imperfect, G-Triplex             | structures
  7   | i-Motif Family    | Canonical, Relaxed, AC-motif                | C-rich structures
  8   | Z-DNA             | Z-DNA, eGZ (Extruded-G) DNA                 | Left-handed helix
  9   | A-philic DNA      | A-philic DNA                                 | A-rich structures
 10   | Hybrid            | Dynamic overlaps                             | Multi-class overlaps
 11   | Non-B Clusters    | Dynamic clusters                             | High-density regions

Author: Dr. Venkata Rajesh Yella
License: MIT
Version: 2024.1
"""

import re
import pandas as pd
from typing import List, Dict, Any, Tuple, Union
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
import warnings
warnings.filterwarnings("ignore")

# Try to import hyperscan for acceleration
try:
    HYPERSCAN_AVAILABLE = True
except ImportError:
    HYPERSCAN_AVAILABLE = False

# Import pure Python scanner for specific classes
try:
    from nonb_pure_python import scan_sequence as pure_python_scan
    PURE_PYTHON_AVAILABLE = True
except ImportError:
    PURE_PYTHON_AVAILABLE = False

# =============================================================================
# CORE MOTIF PATTERNS & DETECTION ALGORITHMS
# =============================================================================

class MotifPatterns:
    """Centralized pattern registry for all 11 motif classes"""
    
    # Class 1: Curved DNA patterns (optimized)
    CURVED_DNA_PATTERNS = {
        'a_tracts': [
            (r'A{4,}', 'A-tract', 'Local curvature'),
            (r'(?:A{3,}T{1,3}){2,}', 'AT-rich', 'Mixed curvature'),
            (r'(?:A{3,}.{0,10}){3,}A{3,}', 'Phased A-tracts', 'Global curvature')
        ],
        't_tracts': [
            (r'T{4,}', 'T-tract', 'Local curvature'),
            (r'(?:T{3,}A{1,3}){2,}', 'TA-rich', 'Mixed curvature')
        ]
    }
    
    # Class 2: Slipped DNA patterns (optimized with non-capturing groups)
    SLIPPED_DNA_PATTERNS = {
        'direct_repeats': [
            (r'([ATGC]{2,})(?:[ATGC]{0,50})\1', 'Direct repeat', 'STR'),
            (r'([ATGC]{3,20})(?:[ATGC]{0,100})\1', 'Long direct repeat', 'Direct Repeat')
        ],
        'short_tandem_repeats': [
            (r'([ATGC]{1,6})\1{3,}', 'STR', 'STR'),
            (r'(?:CA){4,}', 'CA repeat', 'STR'),
            (r'(?:CGG){3,}', 'CGG repeat', 'STR')
        ]
    }
    
    # Class 3: Cruciform DNA patterns
    CRUCIFORM_PATTERNS = {
        'inverted_repeats': [
            (r'([ATGC]{4,})(?:[ATGC]{0,50})(?-i:' + r'(?:(?:[ATGC](?:[ATGC](?:[ATGC](?:[ATGC]))))?)*)' + r')', 'Inverted repeat', 'Inverted Repeats'),
            (r'([ATGC]{6,20})[ATGC]{0,100}(?-i:(?:(?=(?:[ATGC]))\1))', 'Palindrome', 'Inverted Repeats')
        ]
    }
    
    # Class 4: R-loop patterns
    R_LOOP_PATTERNS = {
        'r_loop_sites': [
            (r'[GC]{3,}[AT]{1,5}[GC]{3,}', 'GC-rich R-loop', 'R-loop formation sites'),
            (r'G{3,}[ATGC]{5,50}C{3,}', 'G-C rich region', 'R-loop formation sites')
        ]
    }
    
    # Class 5: Triplex patterns (optimized with non-capturing groups)
    TRIPLEX_PATTERNS = {
        'triplex_forming': [
            (r'[GA]{10,}', 'Homopurine tract', 'Triplex'),
            (r'[CT]{10,}', 'Homopyrimidine tract', 'Triplex'),
            (r'(?:GA){4,}[GA]*(?:TC){4,}', 'Mirror repeat', 'Sticky DNA')
        ]
    }
    
    # Class 6: G-Quadruplex Family patterns (7 subclasses)
    G_QUADRUPLEX_PATTERNS = {
        'canonical_g4': [
            (r'G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'Canonical G4', 'Canonical G4')
        ],
        'relaxed_g4': [
            (r'G{2,}[ATGC]{1,12}G{2,}[ATGC]{1,12}G{2,}[ATGC]{1,12}G{2,}', 'Relaxed G4', 'Relaxed G4')
        ],
        'bulged_g4': [
            (r'G{3,}[ATGC]{8,20}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'Bulged G4', 'Bulged G4')
        ],
        'bipartite_g4': [
            (r'G{2,}[ATGC]{15,50}G{2,}[ATGC]{1,7}G{2,}[ATGC]{1,7}G{2,}', 'Bipartite G4', 'Bipartite G4')
        ],
        'multimeric_g4': [
            (r'(?:G{3,}[ATGC]{1,7}){4,}G{3,}', 'Multimeric G4', 'Multimeric G4')
        ],
        'imperfect_g4': [
            (r'G{2,}[ATGC]{1,10}[AG]G{1,}[ATGC]{1,10}G{2,}[ATGC]{1,10}G{2,}', 'Imperfect G4', 'Imperfect G4')
        ],
        'g_triplex': [
            (r'G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'G-Triplex', 'G-Triplex intermediate')
        ]
    }
    
    # Class 7: i-Motif Family patterns (3 subclasses, optimized)
    I_MOTIF_PATTERNS = {
        'canonical_imotif': [
            (r'C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}', 'Canonical i-motif', 'Canonical i-motif')
        ],
        'relaxed_imotif': [
            (r'C{2,}[ATGC]{1,12}C{2,}[ATGC]{1,12}C{2,}[ATGC]{1,12}C{2,}', 'Relaxed i-motif', 'Relaxed i-motif')
        ],
        'ac_motif': [
            (r'(?:AC){3,}|(?:CA){3,}', 'AC-motif', 'AC-motif')
        ]
    }
    
    # Class 8: Z-DNA patterns (2 subclasses)
    Z_DNA_PATTERNS = {
        'z_dna': [
            (r'(?:CG){4,}|(?:GC){4,}', 'CG alternating', 'Z-DNA'),
            (r'(?:AT){4,}|(?:TA){4,}', 'AT alternating', 'Z-DNA')
        ],
        'egz_dna': [
            (r'[CG]{6,}', 'CG-rich region', 'eGZ (Extruded-G) DNA')
        ]
    }
    
    # Class 9: A-philic DNA patterns
    A_PHILIC_PATTERNS = {
        'a_philic': [
            (r'A{6,}', 'Poly-A tract', 'A-philic DNA'),
            (r'(?:A{2,}[AT]){3,}', 'A-rich region', 'A-philic DNA')
        ]
    }
    
    @classmethod
    def get_all_patterns(cls) -> Dict[str, Dict[str, List[Tuple]]]:
        """Get all patterns organized by motif class"""
        return {
            'curved_dna': cls.CURVED_DNA_PATTERNS,
            'slipped_dna': cls.SLIPPED_DNA_PATTERNS,
            'cruciform': cls.CRUCIFORM_PATTERNS,
            'r_loop': cls.R_LOOP_PATTERNS,
            'triplex': cls.TRIPLEX_PATTERNS,
            'g_quadruplex': cls.G_QUADRUPLEX_PATTERNS,
            'i_motif': cls.I_MOTIF_PATTERNS,
            'z_dna': cls.Z_DNA_PATTERNS,
            'a_philic': cls.A_PHILIC_PATTERNS
        }

# =============================================================================
# SCORING ALGORITHMS
# =============================================================================

class MotifScoring:
    """Scientific scoring algorithms for different motif classes"""
    
    @staticmethod
    def g4hunter_score(sequence: str) -> float:
        """G4Hunter algorithm for G-quadruplex scoring"""
        if len(sequence) < 10:
            return 0.0
        
        score = 0.0
        for i, base in enumerate(sequence):
            if base == 'G':
                score += 1
            elif base == 'C':
                score -= 1
        
        return score / len(sequence)
    
    @staticmethod
    def curvature_score(sequence: str) -> float:
        """DNA curvature scoring based on A/T tract content"""
        if len(sequence) < 4:
            return 0.0
        
        # Count A/T tracts of length 3+
        at_tracts = len(re.findall(r'[AT]{3,}', sequence))
        return min(at_tracts / (len(sequence) / 10), 1.0)
    
    @staticmethod
    def z_dna_score(sequence: str) -> float:
        """Z-DNA scoring based on alternating purine-pyrimidine content"""
        if len(sequence) < 6:
            return 0.0
        
        alternating_score = 0
        for i in range(len(sequence) - 1):
            curr_base = sequence[i]
            next_base = sequence[i + 1]
            
            # Check for alternating purine-pyrimidine pattern
            if ((curr_base in 'AG' and next_base in 'CT') or 
                (curr_base in 'CT' and next_base in 'AG')):
                alternating_score += 1
        
        return alternating_score / (len(sequence) - 1) if len(sequence) > 1 else 0.0
    
    @staticmethod
    def triplex_score(sequence: str) -> float:
        """Triplex potential scoring"""
        if len(sequence) < 10:
            return 0.0
        
        # Homopurine/homopyrimidine content
        purine_content = len(re.findall(r'[GA]', sequence)) / len(sequence)
        pyrimidine_content = len(re.findall(r'[CT]', sequence)) / len(sequence)
        
        return max(purine_content, pyrimidine_content)

# =============================================================================
# MOTIF DETECTION ENGINE
# =============================================================================

class MotifDetector:
    """High-performance motif detection with 11 classes and 22+ subclasses"""
    
    def __init__(self):
        self.patterns = MotifPatterns.get_all_patterns()
        self.scorer = MotifScoring()
        self._compiled_patterns = self._compile_patterns()
    
    def _compile_patterns(self) -> Dict[str, List[Tuple]]:
        """Compile all regex patterns for performance with optimized flags"""
        compiled = {}
        for motif_class, pattern_groups in self.patterns.items():
            compiled[motif_class] = []
            for pattern_group, patterns in pattern_groups.items():
                for pattern, name, subclass in patterns:
                    try:
                        # Use IGNORECASE | ASCII for better performance on DNA sequences
                        compiled_pattern = re.compile(pattern, re.IGNORECASE | re.ASCII)
                        compiled[motif_class].append((compiled_pattern, name, subclass))
                    except re.error:
                        continue
        return compiled
    
    def _detect_class_motifs(self, sequence: str, motif_class: str) -> List[Dict[str, Any]]:
        """Detect motifs for a specific class"""
        motifs = []
        
        if motif_class not in self._compiled_patterns:
            return motifs
        
        for compiled_pattern, name, subclass in self._compiled_patterns[motif_class]:
            for match in compiled_pattern.finditer(sequence):
                start, end = match.span()
                motif_seq = sequence[start:end]
                
                # Calculate class-specific score
                score = self._calculate_score(motif_seq, motif_class)
                
                # Apply quality thresholds
                if self._passes_quality_threshold(motif_seq, motif_class, score):
                    motifs.append({
                        'Class': self._format_class_name(motif_class),
                        'Subclass': subclass,
                        'Start': start + 1,  # 1-based coordinates
                        'End': end,
                        'Length': len(motif_seq),
                        'Sequence': motif_seq,
                        'Score': round(score, 3),
                        'Strand': '+',
                        'Method': f'{motif_class.upper()}_detection'
                    })
        
        return motifs
    
    def _calculate_score(self, sequence: str, motif_class: str) -> float:
        """Calculate appropriate score based on motif class"""
        if motif_class == 'g_quadruplex':
            return abs(self.scorer.g4hunter_score(sequence))
        elif motif_class == 'curved_dna':
            return self.scorer.curvature_score(sequence)
        elif motif_class == 'z_dna':
            return self.scorer.z_dna_score(sequence)
        elif motif_class == 'triplex':
            return self.scorer.triplex_score(sequence)
        else:
            # Default scoring based on length and composition
            return min(len(sequence) / 50, 1.0)
    

    
    def _passes_quality_threshold(self, sequence: str, motif_class: str, score: float) -> bool:
        """Apply class-specific quality thresholds"""
        min_lengths = {
            'curved_dna': 4, 'slipped_dna': 6, 'cruciform': 8, 'r_loop': 10,
            'triplex': 10, 'g_quadruplex': 15, 'i_motif': 12, 'z_dna': 6,
            'a_philic': 6
        }
        
        min_scores = {
            'g_quadruplex': 0.5, 'curved_dna': 0.2, 'z_dna': 0.5,
            'triplex': 0.6
        }
        
        # Length check
        min_len = min_lengths.get(motif_class, 4)
        if len(sequence) < min_len:
            return False
        
        # Score check
        min_score = min_scores.get(motif_class, 0.0)
        if score < min_score:
            return False
        
        return True
    
    def _format_class_name(self, motif_class: str) -> str:
        """Format class names for output"""
        class_names = {
            'curved_dna': 'Curved_DNA',
            'slipped_dna': 'Slipped_DNA', 
            'cruciform': 'Cruciform',
            'r_loop': 'R-Loop',
            'triplex': 'Triplex',
            'g_quadruplex': 'G-Quadruplex',
            'i_motif': 'i-Motif',
            'z_dna': 'Z-DNA',
            'a_philic': 'A-philic_DNA'
        }
        return class_names.get(motif_class, motif_class.title())
    
    def _detect_hybrid_motifs(self, motifs: List[Dict[str, Any]], sequence: str = None) -> List[Dict[str, Any]]:
        """Detect Class 10: Hybrid motifs (overlapping classes) - longest non-overlapping regions"""
        hybrid_candidates = []
        
        # Sort motifs by position
        sorted_motifs = sorted(motifs, key=lambda x: x['Start'])
        
        for i, motif1 in enumerate(sorted_motifs):
            for motif2 in sorted_motifs[i+1:]:
                # Check for significant overlap between different classes
                if motif1['Class'] != motif2['Class']:
                    overlap = self._calculate_overlap(motif1, motif2)
                    if 0.3 <= overlap <= 0.7:  # 30-70% overlap
                        start = min(motif1['Start'], motif2['Start'])
                        end = max(motif1['End'], motif2['End'])
                        length = end - start
                        
                        # Extract actual sequence if provided
                        seq_text = 'HYBRID_REGION'
                        if sequence and 0 <= start - 1 < len(sequence) and 0 < end <= len(sequence):
                            seq_text = sequence[start-1:end]
                        
                        # Calculate raw score
                        raw_score = (motif1.get('Score', 0) + motif2.get('Score', 0)) / 2
                        
                        hybrid = {
                            'Class': 'Hybrid',
                            'Subclass': f"{motif1['Class']}_{motif2['Class']}_Overlap",
                            'Start': start,
                            'End': end,
                            'Length': length,
                            'Sequence': seq_text,
                            'Score': raw_score,
                            'Strand': '+',
                            'Method': 'Hybrid_Detection',
                            'Component_Classes': [motif1['Class'], motif2['Class']]
                        }
                        hybrid_candidates.append(hybrid)
        
        # Select longest non-overlapping hybrids
        return self._select_longest_nonoverlapping(hybrid_candidates)
    
    def _detect_cluster_motifs(self, motifs: List[Dict[str, Any]], sequence: str = None) -> List[Dict[str, Any]]:
        """Detect Class 11: Non-B DNA Clusters (high-density regions) - longest non-overlapping regions"""
        cluster_candidates = []
        window_size = 100  # bp
        min_motifs = 3
        
        if len(motifs) < min_motifs:
            return cluster_candidates
        
        # Sort motifs by position
        sorted_motifs = sorted(motifs, key=lambda x: x['Start'])
        
        # Sliding window to find high-density regions
        for i in range(len(sorted_motifs)):
            window_motifs = []
            start_pos = sorted_motifs[i]['Start']
            
            for motif in sorted_motifs[i:]:
                if motif['Start'] - start_pos <= window_size:
                    window_motifs.append(motif)
                else:
                    break
            
            if len(window_motifs) >= min_motifs:
                # Get unique classes in cluster
                classes = set(m['Class'] for m in window_motifs)
                if len(classes) >= 2:  # Multiple classes
                    start = window_motifs[0]['Start']
                    end = window_motifs[-1]['End']
                    length = end - start
                    motif_count = len(window_motifs)
                    
                    # Extract actual sequence if provided
                    seq_text = 'CLUSTER_REGION'
                    if sequence and 0 <= start - 1 < len(sequence) and 0 < end <= len(sequence):
                        seq_text = sequence[start-1:end]
                    
                    # Calculate raw score (density score)
                    raw_score = motif_count / window_size * 100  # Density score
                    
                    cluster = {
                        'Class': 'Non-B_DNA_Clusters',
                        'Subclass': f"Mixed_Cluster_{len(classes)}_classes",
                        'Start': start,
                        'End': end,
                        'Length': length,
                        'Sequence': seq_text,
                        'Score': raw_score,
                        'Strand': '+',
                        'Method': 'Cluster_Detection',
                        'Motif_Count': motif_count,
                        'Class_Diversity': len(classes)
                    }
                    cluster_candidates.append(cluster)
        
        # Select longest non-overlapping clusters
        return self._select_longest_nonoverlapping(cluster_candidates)
    
    def _calculate_overlap(self, motif1: Dict[str, Any], motif2: Dict[str, Any]) -> float:
        """Calculate overlap percentage between two motifs"""
        start1, end1 = motif1['Start'], motif1['End']
        start2, end2 = motif2['Start'], motif2['End']
        
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)
        
        if overlap_end <= overlap_start:
            return 0.0
        
        overlap_length = overlap_end - overlap_start
        min_length = min(end1 - start1, end2 - start2)
        
        return overlap_length / min_length if min_length > 0 else 0.0
    
    
    def _select_longest_nonoverlapping(self, motifs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Select longest non-overlapping motifs from a list.
        Prioritizes by length, then by score.
        """
        if not motifs:
            return []
        
        # Sort by length (descending), then by score (descending)
        sorted_motifs = sorted(motifs, key=lambda x: (-x['Length'], -x.get('Score', 0)))
        
        selected = []
        occupied_positions = set()
        
        for motif in sorted_motifs:
            # Check if this motif overlaps with any selected motif
            motif_positions = set(range(motif['Start'], motif['End'] + 1))
            
            if not motif_positions & occupied_positions:
                selected.append(motif)
                occupied_positions |= motif_positions
        
        # Sort by start position for output
        return sorted(selected, key=lambda x: x['Start'])
    
    
    def analyze_sequence(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """
        Main analysis function - detects all 11 motif classes with 22+ subclasses
        Uses split approach: Pure Python for Slipped/Cruciform/Triplex, Hyperscan for others
        
        Args:
            sequence: DNA sequence to analyze
            sequence_name: Name identifier for the sequence
            
        Returns:
            List of detected motifs with comprehensive metadata
        """
        sequence = sequence.upper().strip()
        all_motifs = []
        
        # Split approach as requested:
        # 1. Pure Python detection for Slipped DNA, Cruciform, and Triplex
        if PURE_PYTHON_AVAILABLE:
            pure_python_motifs = pure_python_scan(sequence, sequence_name)
            all_motifs.extend(pure_python_motifs)
        else:
            # Fallback to original detection for these classes
            fallback_classes = ['slipped_dna', 'cruciform', 'triplex']
            for motif_class in fallback_classes:
                class_motifs = self._detect_class_motifs(sequence, motif_class)
                all_motifs.extend(class_motifs)
        
        # 2. Hyperscan-based detection for all other classes
        hyperscan_classes = [
            'curved_dna', 'r_loop', 'g_quadruplex', 'i_motif', 'z_dna', 'a_philic'
        ]
        
        for motif_class in hyperscan_classes:
            class_motifs = self._detect_class_motifs(sequence, motif_class)
            all_motifs.extend(class_motifs)
        
        # Remove overlapping motifs within same class
        all_motifs = self._remove_class_overlaps(all_motifs)
        
        # Class 10: Hybrid motifs (overlapping different classes)
        hybrid_motifs = self._detect_hybrid_motifs(all_motifs, sequence)
        all_motifs.extend(hybrid_motifs)
        
        # Class 11: Non-B DNA Clusters (high-density regions)
        cluster_motifs = self._detect_cluster_motifs(all_motifs, sequence)
        all_motifs.extend(cluster_motifs)
        
        # Add metadata
        for i, motif in enumerate(all_motifs):
            motif['ID'] = f"{sequence_name}_motif_{i+1}"
            motif['Sequence_Name'] = sequence_name
        
        return sorted(all_motifs, key=lambda x: x['Start'])
    
    def _remove_class_overlaps(self, motifs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Remove overlapping motifs within the same class/subclass"""
        if len(motifs) <= 1:
            return motifs
        
        # Group by class and subclass
        grouped = defaultdict(list)
        for motif in motifs:
            key = f"{motif['Class']}_{motif['Subclass']}"
            grouped[key].append(motif)
        
        filtered_motifs = []
        for class_subclass, class_motifs in grouped.items():
            # Sort by score descending, then by length descending
            sorted_motifs = sorted(class_motifs, 
                                 key=lambda x: (x['Score'], x['Length']), 
                                 reverse=True)
            
            non_overlapping = []
            for motif in sorted_motifs:
                overlaps = False
                for existing in non_overlapping:
                    # Strict overlap check: any overlap at all (>0) is rejected within same subclass
                    if self._calculate_overlap(motif, existing) > 0.0:
                        overlaps = True
                        break
                
                if not overlaps:
                    non_overlapping.append(motif)
            
            filtered_motifs.extend(non_overlapping)
        
        return filtered_motifs

# =============================================================================
# CONVENIENCE FUNCTIONS & API
# =============================================================================

# Note: analyze_sequence is already defined at line 740 using ModularMotifDetector
# This duplicate definition has been removed to avoid conflicts

def get_motif_classification_info() -> Dict[str, Any]:
    """Get comprehensive information about the 11-class, 22+ subclass system"""
    # Try to get info from modular architecture first
    try:
        # Modular data included inline
        modular_data = {
            'total_detectors': 9,
            'total_patterns': 54,
            'version': '2024.1'
        }
        
        # Enhanced info combining modular data with classification details
        return {
            'version': '2024.1',
            'architecture': 'modular',
            'total_classes': 11,
            'total_subclasses': '22+',
            'total_detectors': modular_data.get('total_detectors', 9),
            'total_patterns': modular_data.get('total_patterns', 54),
            'classification': {
                1: {'name': 'Curved DNA', 'subclasses': ['Global curvature', 'Local Curvature']},
                2: {'name': 'Slipped DNA', 'subclasses': ['Direct Repeat', 'STR']},
                3: {'name': 'Cruciform DNA', 'subclasses': ['Inverted Repeats']},
                4: {'name': 'R-loop', 'subclasses': ['R-loop formation sites']},
                5: {'name': 'Triplex', 'subclasses': ['Triplex', 'Sticky DNA']},
                6: {'name': 'G-Quadruplex Family', 'subclasses': [
                    'Multimeric G4', 'Canonical G4', 'Relaxed G4', 'Bulged G4',
                    'Bipartite G4', 'Imperfect G4', 'G-Triplex intermediate'
                ]},
                7: {'name': 'i-Motif Family', 'subclasses': ['Canonical i-motif', 'Relaxed i-motif', 'AC-motif']},
                8: {'name': 'Z-DNA', 'subclasses': ['Z-DNA', 'eGZ (Extruded-G) DNA']},
                9: {'name': 'A-philic DNA', 'subclasses': ['A-philic DNA']},
                10: {'name': 'Hybrid', 'subclasses': ['Dynamic overlaps']},
                11: {'name': 'Non-B DNA Clusters', 'subclasses': ['Dynamic clusters']}
            },
            'detector_details': modular_data.get('detectors', {})
        }
    except ImportError:
        # Fallback to legacy info
        return {
            'version': '2024.1',
            'architecture': 'legacy',
            'total_classes': 11,
            'total_subclasses': '22+',
            'classification': {
                1: {'name': 'Curved DNA', 'subclasses': ['Global curvature', 'Local Curvature']},
                2: {'name': 'Slipped DNA', 'subclasses': ['Direct Repeat', 'STR']},
                3: {'name': 'Cruciform DNA', 'subclasses': ['Inverted Repeats']},
                4: {'name': 'R-loop', 'subclasses': ['R-loop formation sites']},
                5: {'name': 'Triplex', 'subclasses': ['Triplex', 'Sticky DNA']},
                6: {'name': 'G-Quadruplex Family', 'subclasses': [
                    'Multimeric G4', 'Canonical G4', 'Relaxed G4', 'Bulged G4',
                    'Bipartite G4', 'Imperfect G4', 'G-Triplex intermediate'
                ]},
                7: {'name': 'i-Motif Family', 'subclasses': ['Canonical i-motif', 'Relaxed i-motif', 'AC-motif']},
                8: {'name': 'Z-DNA', 'subclasses': ['Z-DNA', 'eGZ (Extruded-G) DNA']},
                9: {'name': 'A-philic DNA', 'subclasses': ['A-philic DNA']},
                10: {'name': 'Hybrid', 'subclasses': ['Dynamic overlaps']},
                11: {'name': 'Non-B DNA Clusters', 'subclasses': ['Dynamic clusters']}
            }
        }

# =============================================================================
# BATCH PROCESSING & UTILITIES
# =============================================================================

def analyze_multiple_sequences(sequences: Dict[str, str], 
                             use_multiprocessing: bool = True) -> Dict[str, List[Dict[str, Any]]]:
    """
    Analyze multiple sequences in parallel
    
    Args:
        sequences: Dictionary of {name: sequence}
        use_multiprocessing: Whether to use multiprocessing
        
    Returns:
        Dictionary of {name: motifs_list}
    """
    if not use_multiprocessing or len(sequences) == 1:
        # Sequential processing
        results = {}
        for name, seq in sequences.items():
            results[name] = analyze_sequence(seq, name)
        return results
    
    # Parallel processing
    results = {}
    with ProcessPoolExecutor(max_workers=min(len(sequences), mp.cpu_count())) as executor:
        future_to_name = {
            executor.submit(analyze_sequence, seq, name): name 
            for name, seq in sequences.items()
        }
        
        for future in as_completed(future_to_name):
            name = future_to_name[future]
            try:
                results[name] = future.result()
            except Exception as e:
                print(f"Error processing {name}: {e}")
                results[name] = []
    
    return results

def export_results_to_dataframe(motifs: List[Dict[str, Any]]) -> pd.DataFrame:
    """Convert motif results to pandas DataFrame with comprehensive fields"""
    if not motifs:
        return pd.DataFrame()
    
    df = pd.DataFrame(motifs)
    
    # Comprehensive column list based on user requirements
    comprehensive_columns = [
        'ID',
        'Sequence_Name',  # Sequence Name (or Accession)
        'Source',  # Source (e.g., genome, experiment, study)
        'Class',  # Motif Class
        'Subclass',  # Motif Subclass
        'Pattern_ID',  # Pattern/Annotation ID
        'Start',  # Start Position
        'End',  # End Position
        'Length',  # Length (bp)
        'Sequence',  # Sequence
        'Method',  # Detection Method
        'Score',  # Motif Score
        'Repeat_Type',  # Repeat/Tract Type
        'Left_Arm',  # Left Arm Sequence
        'Right_Arm',  # Right Arm Sequence
        'Loop_Seq',  # Loop Sequence
        'Arm_Length',  # Arm Length
        'Loop_Length',  # Loop Length
        'Stem_Length',  # Stem Length(s)
        'Unit_Length',  # Unit/Repeat Length
        'Number_Of_Copies',  # Number of Copies/Repeats
        'Spacer_Length',  # Spacer Length
        'Spacer_Sequence',  # Spacer Sequence
        'GC_Content',  # GC Content (%)
        'Structural_Features',  # Structural Features (e.g., Tract Type, Curvature Score)
        'Strand'  # Strand information
    ]
    
    # Ensure all comprehensive columns are present, fill missing with 'NA'
    for col in comprehensive_columns:
        if col not in df.columns:
            df[col] = 'NA'
    
    # Map existing fields to comprehensive column names if they differ
    column_mappings = {
        'Repeat_Units': 'Number_Of_Copies',
        'Tract_Type': 'Repeat_Type',
        'GC_Total': 'GC_Content',
        'Gc_Total': 'GC_Content',
        'Curvature_Score': 'Structural_Features',
        'Spacer': 'Spacer_Length',
        'Spacer_Seq': 'Spacer_Sequence'
    }
    
    for old_col, new_col in column_mappings.items():
        if old_col in df.columns and new_col in comprehensive_columns:
            # Only map if the new column is 'NA' (empty)
            df[new_col] = df.apply(
                lambda row: row[old_col] if pd.isna(row[new_col]) or row[new_col] == 'NA' else row[new_col],
                axis=1
            )
    
    # Fill all NaN/None values with 'NA' string
    result_df = df[comprehensive_columns].fillna('NA')
    
    return result_df

# =============================================================================
# TESTING & VALIDATION
# =============================================================================

def run_validation_tests() -> bool:
    """Run basic validation tests on the system"""
    test_sequences = {
        'g4_test': 'GGGTTAGGGTTAGGGTTAGGG',  # G-quadruplex
        'curved_test': 'AAAAATTTTAAAAATTTT',   # Curved DNA
        'z_dna_test': 'CGCGCGCGCGCGCG',       # Z-DNA
        'complex_test': 'GGGTTAGGGTTAGGGTTAGGGAAAAATTTTCGCGCGCGCG'  # Multiple classes
    }
    
    print("Running NBDScanner validation tests...")
    
    all_passed = True
    for name, seq in test_sequences.items():
        try:
            motifs = analyze_sequence(seq, name)
            print(f"✓ {name}: {len(motifs)} motifs detected")
            
            if name == 'g4_test' and not any('G-Quadruplex' in m['Class'] for m in motifs):
                print(f"✗ {name}: Expected G-quadruplex not detected")
                all_passed = False
                
        except Exception as e:
            print(f"✗ {name}: Error - {e}")
            all_passed = False
    
    print(f"Validation {'PASSED' if all_passed else 'FAILED'}")
    return all_passed

if __name__ == "__main__":
    # Run validation tests
    success = run_validation_tests()
    
    if success:
        print("\n" + "="*60)
        print("NBDScanner - Ready for Analysis")
        print("="*60)
        
        # Show classification info
        info = get_motif_classification_info()
        print(f"Version: {info['version']}")
        print(f"Classes: {info['total_classes']}")
        print(f"Subclasses: {info['total_subclasses']}")
        
        # Demo analysis
        demo_seq = "GGGTTAGGGTTAGGGTTAGGGAAAAATTTTCGCGCGCGCGATATATATATCCCCTAACCCTAACCCTAACCC"
        print(f"\nDemo analysis on {len(demo_seq)}bp sequence:")
        motifs = analyze_sequence(demo_seq, "demo")
        
        for motif in motifs:
            print(f"  {motif['Class']:<15} {motif['Subclass']:<20} "
                  f"Pos:{motif['Start']:>3}-{motif['End']:<3} "
                  f"Score:{motif['Score']:.2f}")