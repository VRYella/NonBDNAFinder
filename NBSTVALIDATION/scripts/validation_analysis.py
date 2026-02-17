#!/usr/bin/env python3
"""
NBST Validation Analysis Script
================================
Comprehensive validation of NonBDNAFinder against NBST benchmark data.

Author: Automated Analysis Script
Date: 2024
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from Utilities.nonbscanner import NonBScanner, analyze_sequence
from Utilities.utilities import read_fasta_file

# ============================================================================
# CONFIGURATION
# ============================================================================
# Directory structure: scripts are in NBSTVALIDATION/scripts/
BASE_DIR = Path(__file__).parent.parent  # NBSTVALIDATION directory
DATA_DIR = BASE_DIR / "data"
FASTA_FILE = DATA_DIR / "693fc40d26a53.fasta"

# NBST file mapping (based on problem statement)
NBST_FILES = {
    "curved": DATA_DIR / "693fc40d26a53_curved.tsv",      # Global Curvature
    "GQ": DATA_DIR / "693fc40d26a53_GQ.tsv",              # G-Quadruplexes  
    "Z": DATA_DIR / "693fc40d26a53_Z.tsv",                # Z-DNA
    "STR": DATA_DIR / "693fc40d26a53_Slipped_STR.tsv",    # Slipped DNA (STR)
    "DR": DATA_DIR / "693fc40d26a53_slipped_DR.tsv",      # Slipped DNA (Direct Repeat)
    "MR": DATA_DIR / "693fc40d26a53_MR.tsv",              # Mirror Repeats
}

# Class mapping from NBST to NonBDNAFinder
NBST_TO_NBF_CLASS = {
    "curved": "Curved_DNA",
    "GQ": "G-Quadruplex",
    "Z": "Z-DNA",
    "STR": "Slipped_DNA",
    "DR": "Slipped_DNA",
    "MR": "Cruciform",  # Mirror repeats may map to cruciform
}


def load_fasta_sequence(fasta_path):
    """Load sequence from FASTA file."""
    sequences = read_fasta_file(str(fasta_path))
    if sequences:
        # Returns dict {seq_name: sequence}
        first_name = list(sequences.keys())[0]
        return (first_name, sequences[first_name])
    return None


def load_nbst_file(file_path):
    """Load NBST TSV file into DataFrame."""
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


def run_nonbdnafinder(sequence, sequence_name="seq1"):
    """Run NonBDNAFinder on sequence."""
    scanner = NonBScanner(enable_all_detectors=True)
    motifs = scanner.analyze_sequence(sequence, sequence_name)
    return motifs


def calculate_overlap(start1, end1, start2, end2):
    """Calculate overlap between two intervals."""
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_start < overlap_end:
        return overlap_end - overlap_start
    return 0


def calculate_jaccard(start1, end1, start2, end2):
    """Calculate Jaccard index between two intervals."""
    overlap = calculate_overlap(start1, end1, start2, end2)
    if overlap == 0:
        return 0.0
    union = max(end1, end2) - min(start1, start2)
    return overlap / union if union > 0 else 0.0


def match_motifs(nbf_motifs, nbst_df, overlap_threshold=0.5):
    """Match NonBDNAFinder motifs with NBST predictions."""
    matches = []
    nbf_matched = set()
    nbst_matched = set()
    
    for i, nbf in enumerate(nbf_motifs):
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
            matches.append({
                'nbf_idx': i,
                'nbst_idx': best_match,
                'jaccard': best_jaccard,
                'nbf_start': nbf_start,
                'nbf_end': nbf_end,
                'nbst_start': nbst_df.loc[best_match, 'Start'],
                'nbst_end': nbst_df.loc[best_match].get('Stop', nbst_df.loc[best_match].get('End', 0)),
            })
            nbf_matched.add(i)
            nbst_matched.add(best_match)
    
    return matches, nbf_matched, nbst_matched


def analyze_class_results(nbf_motifs, nbst_name, nbst_df, nbf_class):
    """Analyze results for a specific class."""
    # Filter NBF motifs for this class
    nbf_filtered = [m for m in nbf_motifs if m.get('Class') == nbf_class]
    
    matches, nbf_matched, nbst_matched = match_motifs(nbf_filtered, nbst_df)
    
    tp = len(matches)  # True positives (matched)
    fp = len(nbf_filtered) - len(nbf_matched)  # False positives (NBF only)
    fn = len(nbst_df) - len(nbst_matched)  # False negatives (NBST only)
    
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    
    return {
        'class': nbf_class,
        'nbst_name': nbst_name,
        'nbf_count': len(nbf_filtered),
        'nbst_count': len(nbst_df),
        'tp': tp,
        'fp': fp,
        'fn': fn,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'matches': matches,
        'nbf_only_indices': [i for i in range(len(nbf_filtered)) if i not in nbf_matched],
        'nbst_only_indices': [i for i in range(len(nbst_df)) if i not in nbst_matched],
    }


def analyze_overlaps_within_tool(motifs):
    """Analyze overlaps within NonBDNAFinder results."""
    overlaps = {
        'same_class_same_subclass': [],
        'same_class_diff_subclass': [],
        'different_class': [],
    }
    
    for i, m1 in enumerate(motifs):
        for j, m2 in enumerate(motifs):
            if i >= j:
                continue
            
            start1, end1 = m1.get('Start', 0), m1.get('End', 0)
            start2, end2 = m2.get('Start', 0), m2.get('End', 0)
            
            overlap = calculate_overlap(start1, end1, start2, end2)
            if overlap > 0:
                class1, class2 = m1.get('Class'), m2.get('Class')
                sub1, sub2 = m1.get('Subclass'), m2.get('Subclass')
                
                overlap_info = {
                    'motif1_idx': i,
                    'motif2_idx': j,
                    'motif1_class': class1,
                    'motif2_class': class2,
                    'motif1_subclass': sub1,
                    'motif2_subclass': sub2,
                    'overlap_bp': overlap,
                    'start1': start1, 'end1': end1,
                    'start2': start2, 'end2': end2,
                }
                
                if class1 == class2:
                    if sub1 == sub2:
                        overlaps['same_class_same_subclass'].append(overlap_info)
                    else:
                        overlaps['same_class_diff_subclass'].append(overlap_info)
                else:
                    overlaps['different_class'].append(overlap_info)
    
    return overlaps


def generate_overlap_report(overlaps):
    """Generate a report of overlap analysis."""
    report = []
    report.append("=" * 80)
    report.append("OVERLAP ANALYSIS REPORT")
    report.append("=" * 80)
    
    report.append(f"\n1. Same Class, Same Subclass Overlaps: {len(overlaps['same_class_same_subclass'])}")
    if overlaps['same_class_same_subclass']:
        report.append("   (These should be REMOVED - indicates overlap removal may have issues)")
        for ov in overlaps['same_class_same_subclass'][:5]:
            report.append(f"   - {ov['motif1_class']}/{ov['motif1_subclass']}: {ov['start1']}-{ov['end1']} vs {ov['start2']}-{ov['end2']} ({ov['overlap_bp']}bp)")
    else:
        report.append("   (GOOD: No overlaps within same subclass)")
    
    report.append(f"\n2. Same Class, Different Subclass Overlaps: {len(overlaps['same_class_diff_subclass'])}")
    if overlaps['same_class_diff_subclass']:
        # Group by class to analyze
        by_class = defaultdict(list)
        for ov in overlaps['same_class_diff_subclass']:
            by_class[ov['motif1_class']].append(ov)
        
        for cls, cls_overlaps in by_class.items():
            report.append(f"\n   Class: {cls} ({len(cls_overlaps)} overlaps)")
            for ov in cls_overlaps[:3]:
                report.append(f"   - {ov['motif1_subclass']} vs {ov['motif2_subclass']}: {ov['start1']}-{ov['end1']} vs {ov['start2']}-{ov['end2']} ({ov['overlap_bp']}bp)")
    
    report.append(f"\n3. Different Class Overlaps: {len(overlaps['different_class'])}")
    report.append("   (These are ALLOWED - different structure types can coexist)")
    if overlaps['different_class']:
        class_pairs = defaultdict(int)
        for ov in overlaps['different_class']:
            pair = tuple(sorted([ov['motif1_class'], ov['motif2_class']]))
            class_pairs[pair] += 1
        
        report.append("\n   Class pair overlap counts:")
        for pair, count in sorted(class_pairs.items(), key=lambda x: -x[1]):
            report.append(f"   - {pair[0]} + {pair[1]}: {count} overlaps")
    
    return "\n".join(report)


def main():
    """Main validation function."""
    print("=" * 80)
    print("NBST VALIDATION ANALYSIS")
    print("NonBDNAFinder vs NBST Benchmark Comparison")
    print("=" * 80)
    
    # Load sequence
    print("\n[1] Loading FASTA sequence...")
    seq_record = load_fasta_sequence(FASTA_FILE)
    if seq_record is None:
        print("ERROR: Could not load FASTA file")
        return
    
    sequence = seq_record[1] if isinstance(seq_record, tuple) else seq_record
    seq_name = seq_record[0] if isinstance(seq_record, tuple) else "seq1"
    print(f"    Loaded: {seq_name}, Length: {len(sequence)} bp")
    
    # Run NonBDNAFinder
    print("\n[2] Running NonBDNAFinder analysis...")
    nbf_motifs = run_nonbdnafinder(sequence, seq_name)
    print(f"    Found {len(nbf_motifs)} total motifs")
    
    # Class distribution
    class_counts = defaultdict(int)
    subclass_counts = defaultdict(int)
    for m in nbf_motifs:
        class_counts[m.get('Class', 'Unknown')] += 1
        subclass_counts[f"{m.get('Class', 'Unknown')}/{m.get('Subclass', 'Unknown')}"] += 1
    
    print("\n    NonBDNAFinder Class Distribution:")
    for cls, count in sorted(class_counts.items()):
        print(f"    - {cls}: {count}")
    
    # Load NBST benchmark files
    print("\n[3] Loading NBST benchmark files...")
    nbst_data = {}
    for name, path in NBST_FILES.items():
        df = load_nbst_file(path)
        nbst_data[name] = df
        print(f"    - {name}: {len(df)} motifs")
    
    # Compare results
    print("\n[4] Comparing results...")
    comparison_results = []
    
    for nbst_name, nbst_df in nbst_data.items():
        if nbst_df.empty:
            continue
        
        nbf_class = NBST_TO_NBF_CLASS.get(nbst_name, nbst_name)
        result = analyze_class_results(nbf_motifs, nbst_name, nbst_df, nbf_class)
        comparison_results.append(result)
        
        print(f"\n    {nbst_name} ({nbf_class}):")
        print(f"      NBST: {result['nbst_count']}, NBF: {result['nbf_count']}")
        print(f"      TP: {result['tp']}, FP: {result['fp']}, FN: {result['fn']}")
        print(f"      Precision: {result['precision']:.3f}, Recall: {result['recall']:.3f}, F1: {result['f1']:.3f}")
    
    # Analyze overlaps within NBF
    print("\n[5] Analyzing overlaps within NonBDNAFinder results...")
    overlaps = analyze_overlaps_within_tool(nbf_motifs)
    
    print(f"    Same class, same subclass overlaps: {len(overlaps['same_class_same_subclass'])}")
    print(f"    Same class, different subclass overlaps: {len(overlaps['same_class_diff_subclass'])}")
    print(f"    Different class overlaps: {len(overlaps['different_class'])}")
    
    # Generate detailed report
    overlap_report = generate_overlap_report(overlaps)
    print("\n" + overlap_report)
    
    # Save results
    print("\n[6] Saving results...")
    
    # Save NBF motifs to TSV
    nbf_df = pd.DataFrame(nbf_motifs)
    nbf_df.to_csv(DATA_DIR / "nonbdnafinder_results.tsv", sep='\t', index=False)
    print(f"    Saved NonBDNAFinder results to nonbdnafinder_results.tsv")
    
    # Save comparison summary
    summary_df = pd.DataFrame(comparison_results)
    summary_df.to_csv(DATA_DIR / "comparison_summary.tsv", sep='\t', index=False)
    print(f"    Saved comparison summary to comparison_summary.tsv")
    
    # Return results for further analysis
    return {
        'nbf_motifs': nbf_motifs,
        'nbst_data': nbst_data,
        'comparison_results': comparison_results,
        'overlaps': overlaps,
        'class_counts': dict(class_counts),
        'subclass_counts': dict(subclass_counts),
    }


if __name__ == "__main__":
    results = main()
