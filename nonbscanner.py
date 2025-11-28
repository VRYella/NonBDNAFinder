"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    NONBSCANNER - MAIN API MODULE                              ║
║              Professional Non-B DNA Motif Detection Suite                     ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: nonbscanner.py  
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.1 - Professional Edition
LICENSE: MIT

DESCRIPTION:
    Primary API interface for NonBScanner - a state-of-the-art tool for detecting
    and analyzing Non-B DNA structures in genomic sequences.
    
    This module provides the main entry points for:
    - Single sequence analysis
    - Batch/multi-FASTA analysis  
    - Genome-scale analysis (100MB+)
    - Hybrid and cluster detection
    - Results export and visualization

MAIN API FUNCTIONS:
    analyze_sequence(sequence, name) -> List[motif_dict]
        Primary function for single sequence analysis
        
    analyze_fasta(fasta_content) -> Dict[name, List[motif_dict]]
        Analyze multiple sequences from FASTA format
        
    analyze_file(filename) -> Dict[name, List[motif_dict]]
        Analyze sequences from FASTA file
        
    get_motif_info() -> Dict
        Get information about detected motif classes

SUPPORTED MOTIF CLASSES:
    1. Curved DNA (A-tract mediated bending)
    2. Slipped DNA (Direct repeats, STRs)
    3. Cruciform (Inverted repeats)
    4. R-Loop (RNA-DNA hybrids)
    5. Triplex (Three-stranded structures)
    6. G-Quadruplex (Four-stranded G-rich structures, 7 subclasses)
    7. i-Motif (C-rich structures)
    8. Z-DNA (Left-handed helix)
    9. A-philic DNA (A-rich structures)
    10. Hybrid (Multi-class overlaps)
    11. Clusters (High-density regions)

PERFORMANCE:
    - Standard: ~5,800 bp/second (10kb sequences)
    - Optimized: ~24,674 bp/second (100K sequences, fast detectors)
    - Genome-scale: 100MB+ sequences supported
    - Memory efficient: ~5 MB for 100K sequences

EXAMPLE USAGE:
    >>> import nonbscanner as nbs
    >>> motifs = nbs.analyze_sequence("GGGTTAGGGTTAGGGTTAGGG", "test_seq")
    >>> print(f"Found {len(motifs)} motifs")
    >>> 
    >>> # Analyze FASTA file
    >>> results = nbs.analyze_file("sequences.fasta")
    >>> for name, motifs in results.items():
    ...     print(f"{name}: {len(motifs)} motifs")
"""

import os
import re
import math
import warnings
from typing import List, Dict, Any, Optional, Union, Tuple
from collections import defaultdict, Counter
import pandas as pd

warnings.filterwarnings("ignore")

# Import detector classes
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

# Import utilities
from utilities import (
    parse_fasta,
    read_fasta_file,
    validate_sequence,
    export_to_csv,
    export_to_bed,
    export_to_json,
    export_to_excel,
    get_basic_stats,
    calculate_motif_statistics,
    analyze_class_subclass_detection,
    print_detection_report
)

__version__ = "2024.1"
__author__ = "Dr. Venkata Rajesh Yella"

# =============================================================================
# MAIN SCANNER CLASS
# =============================================================================

class NonBScanner:
    """
    Main scanner class orchestrating all motif detectors
    """
    
    def __init__(self, enable_all_detectors: bool = True):
        """
        Initialize NonBScanner with all detector modules
        
        Args:
            enable_all_detectors: Enable all 9 detector classes (default: True)
        """
        self.detectors = {}
        
        if enable_all_detectors:
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
    
    def analyze_sequence(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """
        Detect all Non-B DNA motifs in a sequence with high performance.
        
        # Analysis Process:
        # | Step | Action                                    | Performance  |
        # |------|-------------------------------------------|--------------|
        # | 1    | Validate sequence (ACGT check)            | O(n)         |
        # | 2    | Run 9 specialized detectors in parallel   | O(n) each    |
        # | 3    | Merge results                             | O(m log m)   |
        # | 4    | Sort by position                          | O(m log m)   |
        
        # Output Motif Fields:
        # | Field         | Type  | Description                       |
        # |---------------|-------|-----------------------------------|
        # | ID            | str   | Unique motif identifier           |
        # | Sequence_Name | str   | Source sequence name              |
        # | Class         | str   | Motif class (e.g., 'G_Quadruplex')|
        # | Subclass      | str   | Motif subclass/variant            |
        # | Start         | int   | 1-based start position            |
        # | End           | int   | End position (inclusive)          |
        # | Length        | int   | Motif length in bp                |
        # | Sequence      | str   | Actual DNA sequence               |
        # | Score         | float | Detection confidence (0-1)        |
        # | Strand        | str   | '+' or '-'                        |
        # | Method        | str   | Detection method used             |
        
        Args:
            sequence: DNA sequence to analyze (ATGC characters)
            sequence_name: Identifier for the sequence
            
        Returns:
            List of motif dictionaries sorted by genomic position
            
        Performance:
            - ~5,000-8,000 bp/s on typical sequences
            - Linear O(n) complexity for most detectors
            - Optimized k-mer indexing for repeat detection
            
        Example:
            >>> scanner = NonBScanner()
            >>> motifs = scanner.analyze_sequence("GGGTTAGGGTTAGGGTTAGGG", "test")
            >>> print(f"Found {len(motifs)} motifs")
            >>> for m in motifs:
            ...     print(f"{m['Class']} at {m['Start']}-{m['End']}")
        """
        sequence = sequence.upper().strip()
        
        # Validate sequence
        is_valid, msg = validate_sequence(sequence)
        if not is_valid:
            raise ValueError(f"Invalid sequence: {msg}")
        
        all_motifs = []
        
        # Run all detectors
        for detector_name, detector in self.detectors.items():
            try:
                motifs = detector.detect_motifs(sequence, sequence_name)
                all_motifs.extend(motifs)
            except Exception as e:
                warnings.warn(f"Error in {detector_name} detector: {e}")
                continue
        
        # Remove overlaps within same class
        filtered_motifs = self._remove_overlaps(all_motifs)
        
        # Detect hybrid and cluster motifs
        hybrid_motifs = self._detect_hybrid_motifs(filtered_motifs, sequence)
        cluster_motifs = self._detect_clusters(filtered_motifs, sequence)
        
        # Apply overlap removal to hybrid and cluster motifs too
        hybrid_motifs = self._remove_overlaps(hybrid_motifs)
        cluster_motifs = self._remove_overlaps(cluster_motifs)
        
        final_motifs = filtered_motifs + hybrid_motifs + cluster_motifs
        
        # Sort by position
        final_motifs.sort(key=lambda x: x.get('Start', 0))
        
        return final_motifs
    
    def _remove_overlaps(self, motifs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Remove overlapping motifs within the same class/subclass"""
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
                    # Strict overlap check
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
    
    def _detect_hybrid_motifs(self, motifs: List[Dict[str, Any]], sequence: str) -> List[Dict[str, Any]]:
        """Detect hybrid motifs (overlapping different classes)"""
        hybrid_motifs = []
        
        for i, motif1 in enumerate(motifs):
            for motif2 in motifs[i+1:]:
                if motif1.get('Class') != motif2.get('Class'):
                    overlap = self._calculate_overlap(motif1, motif2)
                    if 0.3 < overlap < 1.0:  # Partial overlap
                        start = min(motif1.get('Start', 0), motif2.get('Start', 0))
                        end = max(motif1.get('End', 0), motif2.get('End', 0))
                        avg_score = (motif1.get('Score', 0) + motif2.get('Score', 0)) / 2
                        
                        # Extract sequence
                        seq_text = 'HYBRID_REGION'
                        if 0 <= start - 1 < len(sequence) and 0 < end <= len(sequence):
                            seq_text = sequence[start-1:end]
                        
                        hybrid_motifs.append({
                            'ID': f"{motif1.get('Sequence_Name', 'seq')}_HYBRID_{start}",
                            'Sequence_Name': motif1.get('Sequence_Name', 'sequence'),
                            'Class': 'Hybrid',
                            'Subclass': f"{motif1.get('Class', '')}_{motif2.get('Class', '')}_Overlap",
                            'Start': start,
                            'End': end,
                            'Length': end - start,
                            'Sequence': seq_text,
                            'Score': round(avg_score, 3),
                            'Strand': '+',
                            'Method': 'Hybrid_Detection',
                            'Component_Classes': [motif1.get('Class'), motif2.get('Class')]
                        })
        
        return hybrid_motifs
    
    def _detect_clusters(self, motifs: List[Dict[str, Any]], sequence: str) -> List[Dict[str, Any]]:
        """Detect high-density non-B DNA clusters"""
        if len(motifs) < 3:
            return []
        
        cluster_motifs = []
        window_size = 500  # 500bp window
        min_density = 3     # Minimum 3 motifs per window
        
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
                # Get unique classes
                classes = set(m.get('Class') for m in window_motifs)
                if len(classes) >= 2:
                    actual_start = min(m.get('Start', 0) for m in window_motifs)
                    actual_end = max(m.get('End', 0) for m in window_motifs)
                    avg_score = sum(m.get('Score', 0) for m in window_motifs) / len(window_motifs)
                    
                    # Extract sequence
                    seq_text = 'CLUSTER_REGION'
                    if 0 <= actual_start - 1 < len(sequence) and 0 < actual_end <= len(sequence):
                        seq_text = sequence[actual_start-1:actual_end]
                    
                    cluster_motifs.append({
                        'ID': f"{sorted_motifs[i].get('Sequence_Name', 'seq')}_CLUSTER_{actual_start}",
                        'Sequence_Name': sorted_motifs[i].get('Sequence_Name', 'sequence'),
                        'Class': 'Non-B_DNA_Clusters',
                        'Subclass': f'Mixed_Cluster_{len(classes)}_classes',
                        'Start': actual_start,
                        'End': actual_end,
                        'Length': actual_end - actual_start,
                        'Sequence': seq_text,
                        'Score': round(avg_score, 3),
                        'Strand': '+',
                        'Method': 'Cluster_Detection',
                        'Motif_Count': len(window_motifs),
                        'Class_Diversity': len(classes)
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


# =============================================================================
# PUBLIC API FUNCTIONS
# =============================================================================

def analyze_sequence(sequence: str, sequence_name: str = "sequence", 
                    use_fast_mode: bool = False) -> List[Dict[str, Any]]:
    """
    Analyze a single DNA sequence for all Non-B DNA motifs (high-performance API).
    
    # Detection Coverage:
    # | Class          | Subclasses | Method              | Speed (Standard) | Speed (Fast Mode) |
    # |----------------|------------|---------------------|------------------|-------------------|
    # | Curved_DNA     | 2          | APR phasing         | ~8000 bp/s       | ~8000 bp/s        |
    # | Slipped_DNA    | 2          | K-mer indexing      | ~8000 bp/s       | ~8000 bp/s        |
    # | Cruciform      | 1          | K-mer indexing      | ~8000 bp/s       | ~8000 bp/s        |
    # | R_Loop         | 3          | QmRLFS algorithm    | ~6000 bp/s       | ~6000 bp/s        |
    # | Triplex        | 2          | K-mer + regex       | ~7000 bp/s       | ~7000 bp/s        |
    # | G_Quadruplex   | 7          | G4Hunter + patterns | ~5000 bp/s       | ~5000 bp/s        |
    # | i_Motif        | 3          | Regex patterns      | ~6000 bp/s       | ~6000 bp/s        |
    # | Z_DNA          | 2          | 10-mer scoring      | ~7000 bp/s       | ~7000 bp/s        |
    # | A_Philic       | 1          | Tetranucleotide     | ~7000 bp/s       | ~7000 bp/s        |
    # Note: Individual detector speeds are similar. Fast mode achieves speedup through parallelization,
    #       running all 9 detectors simultaneously (wall-clock time reduction: ~9x on 9+ cores).
    
    # Output Structure:
    # | Field         | Type  | Always Present | Description              |
    # |---------------|-------|----------------|--------------------------|
    # | ID            | str   | Yes            | Unique identifier        |
    # | Class         | str   | Yes            | Motif class name         |
    # | Subclass      | str   | Yes            | Motif subclass           |
    # | Start         | int   | Yes            | 1-based start position   |
    # | End           | int   | Yes            | End position             |
    # | Score         | float | Yes            | Confidence score (0-1)   |
    # | Sequence      | str   | Yes            | DNA sequence             |
    # | Length        | int   | Yes            | Motif length in bp       |
    # | Sequence_Name | str   | Yes            | Source sequence name     |
    
    Args:
        sequence: DNA sequence (ATGC characters, case-insensitive)
        sequence_name: Identifier for the sequence
        use_fast_mode: Enable parallel processing for ~9x wall-clock speedup (9 parallel threads)
        
    Returns:
        List of motif dictionaries sorted by genomic position
        
    Performance:
        - Standard mode (sequential): ~5,000-8,000 bp/s per detector, total ~50s for all 9
        - Fast mode (9 parallel threads): Same bp/s per detector, but ~9x faster wall-clock time
        - Wall-clock speedup: ~9x on systems with 9+ CPU cores
        - 7.2kb sequence: ~2.0s (standard), ~0.22s (fast, 9x speedup)
        - 72kb sequence: ~20s (standard), ~2.2s (fast, 9x speedup)
        - 720kb sequence: ~200s (standard), ~22s (fast, 9x speedup)
        
    Note: Speedup is wall-clock time reduction, not throughput increase. Each detector
          runs at the same speed, but all 9 run simultaneously in parallel.
        
    Example:
        >>> import nonbscanner as nbs
        >>> # Standard mode (sequential, all detectors run one after another)
        >>> motifs = nbs.analyze_sequence("GGGTTAGGGTTAGGG", "test")
        >>> # Fast mode (parallel, all 9 detectors run simultaneously - 9x faster)
        >>> motifs_fast = nbs.analyze_sequence("GGGTTAGGGTTAGGG", "test", use_fast_mode=True)
        >>> for m in motifs:
        ...     print(f"{m['Class']}: {m['Start']}-{m['End']}, score={m['Score']}")
    """
    if use_fast_mode:
        # Use parallel processing for 9x speedup
        try:
            from parallel_scanner import analyze_sequence_parallel
            return analyze_sequence_parallel(sequence, sequence_name, use_parallel=True)
        except ImportError:
            warnings.warn("Fast mode not available (parallel_scanner module not found), falling back to standard mode")
    
    # Standard mode (sequential)
    scanner = NonBScanner()
    return scanner.analyze_sequence(sequence, sequence_name)


def analyze_fasta(fasta_content: str) -> Dict[str, List[Dict[str, Any]]]:
    """
    Analyze multiple sequences from FASTA format content
    
    Args:
        fasta_content: FASTA format string content
        
    Returns:
        Dictionary mapping sequence_name -> list of motifs
        
    Example:
        >>> fasta = ">seq1\\nGGGTTAGGGTTAGGGTTAGGG\\n>seq2\\nAAAATTTTCCCCGGGG"
        >>> results = analyze_fasta(fasta)
        >>> for name, motifs in results.items():
        ...     print(f"{name}: {len(motifs)} motifs")
    """
    sequences = parse_fasta(fasta_content)
    results = {}
    scanner = NonBScanner()
    
    for name, seq in sequences.items():
        results[name] = scanner.analyze_sequence(seq, name)
    
    return results


def analyze_file(filename: str) -> Dict[str, List[Dict[str, Any]]]:
    """
    Analyze sequences from a FASTA file
    
    Args:
        filename: Path to FASTA file
        
    Returns:
        Dictionary mapping sequence_name -> list of motifs
        
    Example:
        >>> results = analyze_file("sequences.fasta")
        >>> print(f"Analyzed {len(results)} sequences")
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    
    sequences = read_fasta_file(filename)
    results = {}
    scanner = NonBScanner()
    
    for name, seq in sequences.items():
        results[name] = scanner.analyze_sequence(seq, name)
    
    return results


def get_motif_info() -> Dict[str, Any]:
    """
    Get comprehensive information about motif classification system
    
    Returns:
        Dictionary with classification information:
            - version: NonBScanner version
            - total_classes: Number of motif classes (11)
            - total_subclasses: Number of subclasses (22+)
            - classification: Detailed class/subclass mapping
            
    Example:
        >>> info = get_motif_info()
        >>> print(f"NonBScanner v{info['version']}")
        >>> print(f"Classes: {info['total_classes']}, Subclasses: {info['total_subclasses']}")
    """
    return {
        'version': __version__,
        'author': __author__,
        'total_classes': 11,
        'total_subclasses': '22+',
        'classification': {
            1: {'name': 'Curved DNA', 'subclasses': ['Global Curvature', 'Local Curvature']},
            2: {'name': 'Slipped DNA', 'subclasses': ['Direct Repeat', 'STR']},
            3: {'name': 'Cruciform DNA', 'subclasses': ['Inverted Repeats']},
            4: {'name': 'R-loop', 'subclasses': ['R-loop formation sites', 'QmRLFS-m1', 'QmRLFS-m2']},
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


def export_results(motifs: List[Dict[str, Any]], format: str = 'csv', 
                  filename: Optional[str] = None, **kwargs) -> str:
    """
    Export motifs to various formats
    
    Args:
        motifs: List of motif dictionaries
        format: Export format ('csv', 'bed', 'json', 'excel')
        filename: Optional output filename
        **kwargs: Additional format-specific arguments
        
    Returns:
        Formatted string (and writes to file if filename provided)
        
    Example:
        >>> motifs = analyze_sequence("GGGTTAGGGTTAGGGTTAGGG")
        >>> csv_data = export_results(motifs, format='csv', filename='results.csv')
        >>> excel_data = export_results(motifs, format='excel', filename='results.xlsx')
    """
    if format.lower() == 'csv':
        return export_to_csv(motifs, filename)
    elif format.lower() == 'bed':
        sequence_name = kwargs.get('sequence_name', 'sequence')
        return export_to_bed(motifs, sequence_name, filename)
    elif format.lower() == 'json':
        pretty = kwargs.get('pretty', True)
        return export_to_json(motifs, filename, pretty)
    elif format.lower() in ['excel', 'xlsx']:
        if not filename:
            filename = 'nonbscanner_results.xlsx'
        return export_to_excel(motifs, filename)
    else:
        raise ValueError(f"Unsupported format: {format}. Use 'csv', 'bed', 'json', or 'excel'")


# =============================================================================
# BATCH PROCESSING UTILITIES
# =============================================================================

def analyze_multiple_sequences(sequences: Dict[str, str], 
                              use_multiprocessing: bool = False) -> Dict[str, List[Dict[str, Any]]]:
    """
    Analyze multiple sequences in batch
    
    Args:
        sequences: Dictionary of {name: sequence}
        use_multiprocessing: Enable parallel processing (default: False)
        
    Returns:
        Dictionary of {name: motifs_list}
    """
    results = {}
    scanner = NonBScanner()
    
    if use_multiprocessing:
        # Parallel processing
        from concurrent.futures import ProcessPoolExecutor, as_completed
        import multiprocessing as mp
        
        with ProcessPoolExecutor(max_workers=min(len(sequences), mp.cpu_count())) as executor:
            future_to_name = {
                executor.submit(scanner.analyze_sequence, seq, name): name 
                for name, seq in sequences.items()
            }
            
            for future in as_completed(future_to_name):
                name = future_to_name[future]
                try:
                    results[name] = future.result()
                except Exception as e:
                    warnings.warn(f"Error processing {name}: {e}")
                    results[name] = []
    else:
        # Sequential processing
        for name, seq in sequences.items():
            results[name] = scanner.analyze_sequence(seq, name)
    
    return results


def get_summary_statistics(results: Dict[str, List[Dict[str, Any]]]) -> pd.DataFrame:
    """
    Generate summary statistics for multiple sequence analysis
    
    Args:
        results: Dictionary of {sequence_name: motifs_list}
        
    Returns:
        Pandas DataFrame with summary statistics
    """
    summary_data = []
    
    for name, motifs in results.items():
        stats = calculate_motif_statistics(motifs, 0)  # Length not needed for summary
        
        summary_data.append({
            'Sequence_Name': name,
            'Total_Motifs': stats.get('Total_Motifs', 0),
            'Classes_Detected': stats.get('Classes_Detected', 0),
            'Subclasses_Detected': stats.get('Subclasses_Detected', 0),
            'Score_Mean': stats.get('Score_Mean', 0),
            'Length_Mean': stats.get('Length_Mean', 0)
        })
    
    return pd.DataFrame(summary_data)


# =============================================================================
# MAIN & TESTING
# =============================================================================

def main():
    """Main function for command-line usage and testing"""
    print("="*70)
    print("NonBScanner - Non-B DNA Motif Detection Suite")
    print(f"Version {__version__} by {__author__}")
    print("="*70)
    
    # Get motif info
    info = get_motif_info()
    print(f"\nDetection Capabilities:")
    print(f"  Total Classes: {info['total_classes']}")
    print(f"  Total Subclasses: {info['total_subclasses']}")
    
    # Demo analysis
    demo_seq = "GGGTTAGGGTTAGGGTTAGGGAAAAATTTTCGCGCGCGCGATATATATATCCCCTAACCCTAACCCTAACCC"
    print(f"\nDemo Analysis:")
    print(f"  Sequence: {demo_seq[:50]}... ({len(demo_seq)} bp)")
    
    motifs = analyze_sequence(demo_seq, "demo")
    
    print(f"\nResults: {len(motifs)} motifs detected")
    for motif in motifs[:5]:  # Show first 5
        print(f"  {motif['Class']:<15} {motif['Subclass']:<20} "
              f"Pos:{motif['Start']:>3}-{motif['End']:<3} "
              f"Score:{motif['Score']:.2f}")
    
    if len(motifs) > 5:
        print(f"  ... and {len(motifs) - 5} more")
    
    print("\nNonBScanner ready for analysis!")
    print("="*70)


if __name__ == "__main__":
    main()
