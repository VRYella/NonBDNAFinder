#!/usr/bin/env python3
"""
NBST Validation Analysis - Subclass Level
==========================================
Enhanced validation comparing NonBDNAFinder subclasses against NBST benchmark data.

This script provides granular subclass-level analysis rather than coarse class-level comparison,
enabling more precise validation of detection performance.

Author: NonBDNAFinder Validation Suite
Date: 2026
Version: 2.0
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Set, Optional, Any
import warnings
warnings.filterwarnings('ignore')

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from Utilities.nonbscanner import NonBScanner, analyze_sequence
from Utilities.utilities import read_fasta_file
from Utilities.safety import filter_valid_indices
import logging

logger = logging.getLogger(__name__)

# ============================================================================
# CONFIGURATION
# ============================================================================
VALIDATION_DIR = Path(__file__).parent
FASTA_FILE = VALIDATION_DIR / "693fc40d26a53.fasta"

# NBST file mapping
NBST_FILES = {
    "curved": VALIDATION_DIR / "693fc40d26a53_curved.tsv",
    "GQ": VALIDATION_DIR / "693fc40d26a53_GQ.tsv",
    "Z": VALIDATION_DIR / "693fc40d26a53_Z.tsv",
    "STR": VALIDATION_DIR / "693fc40d26a53_Slipped_STR.tsv",
    "DR": VALIDATION_DIR / "693fc40d26a53_slipped_DR.tsv",
    "MR": VALIDATION_DIR / "693fc40d26a53_MR.tsv",
}

# Subclass to NBST mapping - defines which NBF subclasses correspond to which NBST classes
# NOTE: This is a comprehensive mapping of all possible subclasses. Not all subclasses
# will be present in every analysis run - only subclasses detected in the input sequence
# will appear in the results. Some subclasses (e.g., Telomeric G4, Stacked G4s) may be
# rare or absent depending on the sequence content.
SUBCLASS_TO_NBST_MAP = {
    # G-Quadruplex subclasses all map to NBST GQ
    "Canonical intramolecular G4": "GQ",
    "Extended-loop canonical": "GQ",
    "Intramolecular G-triplex": "GQ",
    "Two-tetrad weak PQS": "GQ",
    "Telomeric G4": "GQ",
    "Stacked canonical G4s": "GQ",
    "Stacked G4s with linker": "GQ",
    "Higher-order G4 array/G4-wire": "GQ",
    
    # Slipped DNA subclasses
    "STR": "STR",
    "Direct Repeat": "DR",
    
    # Curved DNA subclasses
    "Global Curvature": "curved",
    "Local Curvature": "curved",
    
    # Z-DNA
    "Z-DNA": "Z",
    
    # Cruciform (may correspond to Mirror Repeats)
    "Cruciform forming IRs": "MR",
    
    # i-Motif subclasses (no NBST equivalent but include for completeness)
    "Canonical i-motif": None,
    "AC-motif": None,
    
    # R-Loop (no NBST equivalent)
    "R-loop formation sites": None,
    
    # Triplex (no NBST equivalent)
    "Triplex_Motif": None,
    
    # A-philic (no NBST equivalent)
    "A-philic DNA": None,
}


# ============================================================================
# DATA LOADING FUNCTIONS
# ============================================================================

def load_fasta_sequence(fasta_path: Path) -> Optional[Tuple[str, str]]:
    """
    Load sequence from FASTA file.
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        Tuple of (sequence_name, sequence) or None if loading fails
    """
    sequences = read_fasta_file(str(fasta_path))
    if sequences:
        first_name = list(sequences.keys())[0]
        return (first_name, sequences[first_name])
    return None


def load_nbst_file(file_path: Path) -> pd.DataFrame:
    """
    Load NBST TSV file into DataFrame.
    
    Args:
        file_path: Path to NBST TSV file
        
    Returns:
        DataFrame with standardized column names
    """
    if not file_path.exists():
        print(f"Warning: NBST file not found: {file_path}")
        return pd.DataFrame()
    
    try:
        df = pd.read_csv(file_path, sep='\t')
        # Standardize column names
        df.columns = [col.strip() for col in df.columns]
        return df
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return pd.DataFrame()


def run_nonbdnafinder(sequence: str, sequence_name: str = "seq1") -> List[Dict[str, Any]]:
    """
    Run NonBDNAFinder analysis on sequence.
    
    Args:
        sequence: DNA sequence string
        sequence_name: Name for the sequence
        
    Returns:
        List of detected motif dictionaries
    """
    scanner = NonBScanner(enable_all_detectors=True)
    motifs = scanner.analyze_sequence(sequence, sequence_name)
    return motifs


# ============================================================================
# GEOMETRIC OVERLAP CALCULATIONS
# ============================================================================

def calculate_overlap(start1: int, end1: int, start2: int, end2: int) -> int:
    """
    Calculate base pair overlap between two genomic intervals.
    
    Args:
        start1, end1: First interval coordinates
        start2, end2: Second interval coordinates
        
    Returns:
        Number of overlapping base pairs (0 if no overlap)
    """
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_start < overlap_end:
        return overlap_end - overlap_start
    return 0


def calculate_jaccard(start1: int, end1: int, start2: int, end2: int) -> float:
    """
    Calculate Jaccard similarity index between two genomic intervals.
    
    Jaccard = |Intersection| / |Union|
    
    Args:
        start1, end1: First interval coordinates
        start2, end2: Second interval coordinates
        
    Returns:
        Jaccard index (0.0 to 1.0)
    """
    overlap = calculate_overlap(start1, end1, start2, end2)
    if overlap == 0:
        return 0.0
    union = max(end1, end2) - min(start1, start2)
    return overlap / union if union > 0 else 0.0


def calculate_overlap_fraction(start1: int, end1: int, start2: int, end2: int) -> Tuple[float, float]:
    """
    Calculate what fraction of each interval is covered by overlap.
    
    Args:
        start1, end1: First interval coordinates
        start2, end2: Second interval coordinates
        
    Returns:
        Tuple of (fraction_of_interval1, fraction_of_interval2)
    """
    overlap = calculate_overlap(start1, end1, start2, end2)
    len1 = end1 - start1
    len2 = end2 - start2
    
    frac1 = overlap / len1 if len1 > 0 else 0.0
    frac2 = overlap / len2 if len2 > 0 else 0.0
    
    return frac1, frac2


# ============================================================================
# SUBCLASS-LEVEL MATCHING
# ============================================================================

def match_motifs_subclass_level(
    nbf_motifs: List[Dict[str, Any]], 
    nbst_df: pd.DataFrame,
    nbf_subclass: str,
    overlap_threshold: float = 0.5
) -> Tuple[List[Dict], Set[int], Set[int]]:
    """
    Match NonBDNAFinder motifs of a specific subclass with NBST predictions.
    
    This implements the core subclass-level comparison logic.
    
    Args:
        nbf_motifs: List of NBF motif dictionaries
        nbst_df: DataFrame of NBST predictions
        nbf_subclass: Specific NBF subclass to match
        overlap_threshold: Minimum Jaccard index for match (default 0.5)
        
    Returns:
        Tuple of (matches, nbf_matched_indices, nbst_matched_indices)
    """
    matches = []
    nbf_matched = set()
    nbst_matched = set()
    
    # Filter NBF motifs for this specific subclass
    subclass_motifs = [(i, m) for i, m in enumerate(nbf_motifs) 
                       if m.get('Subclass') == nbf_subclass]
    
    for motif_idx, nbf in subclass_motifs:
        nbf_start = nbf.get('Start', 0)
        nbf_end = nbf.get('End', 0)
        
        best_match = None
        best_jaccard = 0
        
        for j, nbst in nbst_df.iterrows():
            nbst_start = nbst.get('Start', 0)
            nbst_end = nbst.get('Stop', nbst.get('End', 0))
            
            jaccard = calculate_jaccard(nbf_start, nbf_end, nbst_start, nbst_end)
            
            if jaccard > overlap_threshold and jaccard > best_jaccard:
                best_jaccard = jaccard
                best_match = j
        
        if best_match is not None:
            nbst_start = nbst_df.loc[best_match, 'Start']
            nbst_end = nbst_df.loc[best_match].get('Stop', nbst_df.loc[best_match].get('End', 0))
            overlap_bp = calculate_overlap(nbf_start, nbf_end, nbst_start, nbst_end)
            frac_nbf, frac_nbst = calculate_overlap_fraction(nbf_start, nbf_end, nbst_start, nbst_end)
            
            matches.append({
                'nbf_idx': motif_idx,
                'nbst_idx': best_match,
                'jaccard': best_jaccard,
                'overlap_bp': overlap_bp,
                'nbf_start': nbf_start,
                'nbf_end': nbf_end,
                'nbf_length': nbf_end - nbf_start,
                'nbst_start': nbst_start,
                'nbst_end': nbst_end,
                'nbst_length': nbst_end - nbst_start,
                'fraction_nbf_covered': frac_nbf,
                'fraction_nbst_covered': frac_nbst,
            })
            nbf_matched.add(motif_idx)
            nbst_matched.add(best_match)
    
    return matches, nbf_matched, nbst_matched


def analyze_subclass_results(
    nbf_motifs: List[Dict[str, Any]],
    nbst_name: str,
    nbst_df: pd.DataFrame,
    nbf_subclass: str
) -> Dict[str, Any]:
    """
    Analyze validation results for a specific NBF subclass vs NBST class.
    
    Args:
        nbf_motifs: All NBF motifs
        nbst_name: NBST class name (e.g., "GQ", "STR")
        nbst_df: NBST predictions DataFrame
        nbf_subclass: Specific NBF subclass to analyze
        
    Returns:
        Dictionary with performance metrics and match details
    """
    # Filter NBF motifs for this subclass
    subclass_motifs = [m for m in nbf_motifs if m.get('Subclass') == nbf_subclass]
    
    # Get the NBF class from the first motif of this subclass
    nbf_class = subclass_motifs[0].get('Class', 'Unknown') if subclass_motifs else 'Unknown'
    
    matches, nbf_matched, nbst_matched = match_motifs_subclass_level(
        nbf_motifs, nbst_df, nbf_subclass
    )
    
    # Calculate metrics
    tp = len(matches)  # True positives (matched)
    fp = len(subclass_motifs) - len(nbf_matched)  # False positives (NBF only)
    fn = len(nbst_df) - len(nbst_matched)  # False negatives (NBST only)
    
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    
    # Calculate match quality metrics
    if matches:
        jaccards = [m['jaccard'] for m in matches]
        overlaps = [m['overlap_bp'] for m in matches]
        mean_jaccard = np.mean(jaccards)
        median_jaccard = np.median(jaccards)
        mean_overlap = np.mean(overlaps)
        median_overlap = np.median(overlaps)
        min_jaccard = min(jaccards)
        max_jaccard = max(jaccards)
    else:
        mean_jaccard = median_jaccard = mean_overlap = median_overlap = 0.0
        min_jaccard = max_jaccard = 0.0
    
    # Get indices of unmatched motifs
    all_nbf_indices = set(i for i, m in enumerate(nbf_motifs) if m.get('Subclass') == nbf_subclass)
    nbf_only_indices = list(all_nbf_indices - nbf_matched)
    nbst_only_indices = list(set(range(len(nbst_df))) - nbst_matched)
    
    return {
        'nbf_subclass': nbf_subclass,
        'nbf_class': nbf_class,
        'nbst_class': nbst_name,
        'nbf_count': len(subclass_motifs),
        'nbst_count': len(nbst_df),
        'tp': tp,
        'fp': fp,
        'fn': fn,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'mean_jaccard': mean_jaccard,
        'median_jaccard': median_jaccard,
        'mean_overlap_bp': mean_overlap,
        'median_overlap_bp': median_overlap,
        'min_jaccard': min_jaccard,
        'max_jaccard': max_jaccard,
        'num_matches': len(matches),
        'matches': matches,
        'nbf_only_indices': nbf_only_indices,
        'nbst_only_indices': nbst_only_indices,
    }


# ============================================================================
# FALSE POSITIVE/NEGATIVE ANALYSIS
# ============================================================================

def analyze_false_positives(
    nbf_motifs: List[Dict[str, Any]],
    result: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Analyze characteristics of false positive detections.
    
    Args:
        nbf_motifs: All NBF motifs
        result: Analysis result dictionary
        
    Returns:
        Dictionary with FP analysis
    """
    fp_indices = result['nbf_only_indices']
    
    if not fp_indices:
        return {
            'count': 0,
            'mean_length': 0,
            'mean_score': 0,
            'common_characteristics': 'N/A'
        }
    
    # Filter invalid indices before access to prevent crashes
    valid_fp_indices = filter_valid_indices(fp_indices, len(nbf_motifs))
    
    if not valid_fp_indices:
        logger.warning("All FP indices were out of bounds")
        return {
            'count': 0,
            'mean_length': 0,
            'mean_score': 0,
            'common_characteristics': 'Invalid indices'
        }
    
    # Safe access with filtered indices
    fp_motifs = [nbf_motifs[i] for i in valid_fp_indices]
    lengths = [m.get('Length', 0) for m in fp_motifs]
    scores = [m.get('Score', 0) for m in fp_motifs if 'Score' in m]
    
    characteristics = []
    if lengths:
        mean_len = np.mean(lengths)
        if mean_len < 20:
            characteristics.append("Short motifs (<20bp)")
        elif mean_len > 100:
            characteristics.append("Long motifs (>100bp)")
    
    if scores:
        mean_score = np.mean(scores)
        if mean_score < 1.5:
            characteristics.append("Low scores (<1.5)")
    
    return {
        'count': len(valid_fp_indices),
        'mean_length': np.mean(lengths) if lengths else 0,
        'median_length': np.median(lengths) if lengths else 0,
        'mean_score': np.mean(scores) if scores else 0,
        'median_score': np.median(scores) if scores else 0,
        'common_characteristics': '; '.join(characteristics) if characteristics else 'Mixed characteristics',
    }


def analyze_false_negatives(
    nbst_df: pd.DataFrame,
    result: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Analyze characteristics of false negative (missed) NBST predictions.
    
    Args:
        nbst_df: NBST predictions DataFrame
        result: Analysis result dictionary
        
    Returns:
        Dictionary with FN analysis
    """
    fn_indices = result['nbst_only_indices']
    
    if not fn_indices:
        return {
            'count': 0,
            'mean_length': 0,
            'potential_causes': 'N/A'
        }
    
    # Filter invalid indices before access to prevent crashes
    valid_indices = filter_valid_indices(fn_indices, len(nbst_df))
    
    if not valid_indices:
        logger.warning("All FN indices were out of bounds")
        return {
            'count': 0,
            'mean_length': 0,
            'potential_causes': 'Invalid indices (out of bounds)'
        }
    
    fn_nbst = nbst_df.iloc[valid_indices]
    lengths = []
    for _, row in fn_nbst.iterrows():
        start = row.get('Start', 0)
        end = row.get('Stop', row.get('End', 0))
        lengths.append(end - start)
    
    causes = []
    if lengths:
        mean_len = np.mean(lengths)
        if mean_len < 15:
            causes.append("Short motifs below detection threshold")
        if mean_len > 200:
            causes.append("Very long motifs may exceed search windows")
    
    return {
        'count': len(valid_indices),
        'mean_length': np.mean(lengths) if lengths else 0,
        'median_length': np.median(lengths) if lengths else 0,
        'potential_causes': '; '.join(causes) if causes else 'Parameter/algorithm differences',
    }


# ============================================================================
# MAIN VALIDATION WORKFLOW
# ============================================================================

def main() -> Dict[str, Any]:
    """
    Execute comprehensive subclass-level validation analysis.
    
    Returns:
        Dictionary containing all analysis results
    """
    print("=" * 80)
    print("NBST VALIDATION ANALYSIS - SUBCLASS LEVEL")
    print("NonBDNAFinder Subclasses vs NBST Classes")
    print("=" * 80)
    
    # Load sequence
    print("\n[1] Loading FASTA sequence...")
    seq_record = load_fasta_sequence(FASTA_FILE)
    if seq_record is None:
        print("ERROR: Could not load FASTA file")
        return {}
    
    sequence = seq_record[1]
    seq_name = seq_record[0]
    print(f"    Loaded: {seq_name}, Length: {len(sequence):,} bp")
    
    # Run NonBDNAFinder
    print("\n[2] Running NonBDNAFinder analysis...")
    nbf_motifs = run_nonbdnafinder(sequence, seq_name)
    print(f"    Found {len(nbf_motifs)} total motifs")
    
    # Analyze subclass distribution
    print("\n[3] Analyzing subclass distribution...")
    subclass_counts = defaultdict(int)
    class_to_subclasses = defaultdict(set)
    
    for m in nbf_motifs:
        cls = m.get('Class', 'Unknown')
        subclass = m.get('Subclass', 'Unknown')
        subclass_counts[subclass] += 1
        class_to_subclasses[cls].add(subclass)
    
    print(f"    Found {len(subclass_counts)} unique subclasses across {len(class_to_subclasses)} classes")
    print("\n    Subclass distribution:")
    for subclass, count in sorted(subclass_counts.items(), key=lambda x: -x[1])[:15]:
        print(f"      - {subclass}: {count}")
    
    # Load NBST benchmark files
    print("\n[4] Loading NBST benchmark files...")
    nbst_data = {}
    total_nbst = 0
    for name, path in NBST_FILES.items():
        df = load_nbst_file(path)
        nbst_data[name] = df
        total_nbst += len(df)
        print(f"    - {name}: {len(df)} motifs")
    print(f"    Total NBST motifs: {total_nbst}")
    
    # Perform subclass-level comparisons
    print("\n[5] Performing subclass-level comparisons...")
    comparison_results = []
    
    for subclass, count in subclass_counts.items():
        # Determine which NBST class this subclass should be compared to
        nbst_class = SUBCLASS_TO_NBST_MAP.get(subclass)
        
        if nbst_class is None:
            print(f"    - {subclass}: No NBST equivalent (NBF-specific structure)")
            continue
        
        if nbst_class not in nbst_data or nbst_data[nbst_class].empty:
            print(f"    - {subclass}: NBST class '{nbst_class}' not available")
            continue
        
        # Perform comparison
        result = analyze_subclass_results(
            nbf_motifs, 
            nbst_class, 
            nbst_data[nbst_class],
            subclass
        )
        
        # Add FP/FN analysis
        result['fp_analysis'] = analyze_false_positives(nbf_motifs, result)
        result['fn_analysis'] = analyze_false_negatives(nbst_data[nbst_class], result)
        
        comparison_results.append(result)
        
        print(f"\n    {subclass} vs NBST {nbst_class}:")
        print(f"      NBF: {result['nbf_count']}, NBST: {result['nbst_count']}")
        print(f"      TP: {result['tp']}, FP: {result['fp']}, FN: {result['fn']}")
        print(f"      Precision: {result['precision']:.3f}, Recall: {result['recall']:.3f}, F1: {result['f1']:.3f}")
        print(f"      Mean Jaccard: {result['mean_jaccard']:.3f}")
    
    # Save results
    print("\n[6] Saving results...")
    
    # Save subclass comparison summary
    summary_df = pd.DataFrame([{
        'NBF_Subclass': r['nbf_subclass'],
        'NBF_Class': r['nbf_class'],
        'NBST_Class': r['nbst_class'],
        'NBF_Count': r['nbf_count'],
        'NBST_Count': r['nbst_count'],
        'True_Positives': r['tp'],
        'False_Positives': r['fp'],
        'False_Negatives': r['fn'],
        'Precision': round(r['precision'], 4),
        'Recall': round(r['recall'], 4),
        'F1_Score': round(r['f1'], 4),
        'Mean_Jaccard': round(r['mean_jaccard'], 4),
        'Median_Jaccard': round(r['median_jaccard'], 4),
        'Mean_Overlap_BP': round(r['mean_overlap_bp'], 2),
        'Num_Matches': r['num_matches'],
    } for r in comparison_results])
    
    summary_df.to_csv(VALIDATION_DIR / "subclass_comparison_summary.tsv", sep='\t', index=False)
    print(f"    Saved: subclass_comparison_summary.tsv")
    
    # Save detailed match quality table (Table 2)
    match_quality_rows = []
    for r in comparison_results:
        match_quality_rows.append({
            'NBF_Subclass': r['nbf_subclass'],
            'NBST_Class': r['nbst_class'],
            'Num_Matches': r['num_matches'],
            'Mean_Overlap_BP': round(r['mean_overlap_bp'], 2),
            'Median_Overlap_BP': round(r['median_overlap_bp'], 2),
            'Mean_Jaccard': round(r['mean_jaccard'], 4),
            'Best_Match_Jaccard': round(r['max_jaccard'], 4),
            'Worst_Match_Jaccard': round(r['min_jaccard'], 4),
        })
    
    match_quality_df = pd.DataFrame(match_quality_rows)
    match_quality_df.to_csv(VALIDATION_DIR / "table2_subclass_match_quality.csv", index=False)
    print(f"    Saved: table2_subclass_match_quality.csv")
    
    # Save false positive analysis table (Table 3)
    fp_analysis_rows = []
    for r in comparison_results:
        fp_analysis_rows.append({
            'NBF_Subclass': r['nbf_subclass'],
            'FP_Count': r['fp_analysis']['count'],
            'Mean_Length': round(r['fp_analysis']['mean_length'], 2),
            'Mean_Score': round(r['fp_analysis']['mean_score'], 3),
            'Common_Characteristics': r['fp_analysis']['common_characteristics'],
        })
    
    fp_df = pd.DataFrame(fp_analysis_rows)
    fp_df.to_csv(VALIDATION_DIR / "table3_false_positive_analysis.csv", index=False)
    print(f"    Saved: table3_false_positive_analysis.csv")
    
    # Save false negative analysis table (Table 4)
    fn_analysis_rows = []
    for r in comparison_results:
        fn_analysis_rows.append({
            'NBST_Class': r['nbst_class'],
            'NBF_Subclass_Checked': r['nbf_subclass'],
            'FN_Count': r['fn_analysis']['count'],
            'Mean_Length': round(r['fn_analysis']['mean_length'], 2),
            'Potential_Causes': r['fn_analysis']['potential_causes'],
        })
    
    fn_df = pd.DataFrame(fn_analysis_rows)
    fn_df.to_csv(VALIDATION_DIR / "table4_false_negative_analysis.csv", index=False)
    print(f"    Saved: table4_false_negative_analysis.csv")
    
    # Save Table 1: Subclass-Level Performance Metrics
    perf_metrics_df = summary_df.copy()
    perf_metrics_df.to_csv(VALIDATION_DIR / "table1_subclass_metrics.csv", index=False)
    print(f"    Saved: table1_subclass_metrics.csv")
    
    print("\n" + "=" * 80)
    print("SUBCLASS-LEVEL VALIDATION COMPLETE")
    print("=" * 80)
    
    return {
        'nbf_motifs': nbf_motifs,
        'nbst_data': nbst_data,
        'comparison_results': comparison_results,
        'subclass_counts': dict(subclass_counts),
        'summary_df': summary_df,
    }


if __name__ == "__main__":
    results = main()
