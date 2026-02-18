"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    MULTIFASTA AGGREGATION ENGINE                              ║
║     Scalable Aggregation-First Engine for Large MultiFASTA Datasets          ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: multifasta_engine.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.5
LICENSE: MIT

DESCRIPTION:
    Aggregation-first MultiFASTA engine designed for scalability (>20 sequences).
    No full DataFrame materialization to maintain memory efficiency.
    
    Key Features:
    - Global class distribution across all sequences
    - % sequence coverage per motif class
    - Motif density histogram (motifs per kb)
    - Positional conservation analysis (equal-length only)
    - Memory-safe windowed aggregation
    - 100MB ceiling compatible
    
USAGE:
    engine = MultiFastaEngine(annotations_by_sequence, sequence_lengths)
    
    # Get global statistics
    class_dist = engine.global_class_distribution()
    coverage = engine.class_sequence_coverage()
    densities = engine.motif_density_per_sequence()
    
    # For equal-length sequences
    if all_equal_length:
        conservation = engine.positional_conservation(seq_length)
"""

from __future__ import annotations
from typing import Dict, List, Any
from collections import defaultdict
import math


class MultiFastaEngine:
    """
    Aggregation-first MultiFASTA engine.
    Designed for scalability (>20 sequences).
    No full DataFrame materialization.
    """

    def __init__(self, annotations_by_sequence: Dict[str, List[Dict[str, Any]]], 
                 sequence_lengths: Dict[str, int]):
        """
        Initialize MultiFASTA aggregation engine.
        
        Args:
            annotations_by_sequence: Dict mapping sequence ID to list of motif annotations
            sequence_lengths: Dict mapping sequence ID to sequence length
        """
        self.annotations = annotations_by_sequence
        self.sequence_lengths = sequence_lengths
        self.num_sequences = len(sequence_lengths)

    # ---------------------------------------------------------------------
    # GLOBAL CLASS DISTRIBUTION
    # ---------------------------------------------------------------------
    def global_class_distribution(self) -> Dict[str, int]:
        """
        Compute global motif class distribution across all sequences.
        
        Returns:
            Dict mapping motif class to total count across all sequences
        """
        distribution = defaultdict(int)
        for motifs in self.annotations.values():
            for motif in motifs:
                distribution[motif.get("Class")] += 1
        return dict(distribution)

    # ---------------------------------------------------------------------
    # % SEQUENCES CONTAINING EACH CLASS
    # ---------------------------------------------------------------------
    def class_sequence_coverage(self) -> Dict[str, float]:
        """
        Calculate % of sequences containing each motif class.
        
        Returns:
            Dict mapping motif class to percentage of sequences containing it
        """
        coverage = defaultdict(int)

        for seq, motifs in self.annotations.items():
            seen_classes = set()
            for motif in motifs:
                seen_classes.add(motif.get("Class"))

            for cls in seen_classes:
                coverage[cls] += 1

        # convert to percentage
        return {
            cls: (count / self.num_sequences) * 100
            for cls, count in coverage.items()
        }

    # ---------------------------------------------------------------------
    # MOTIF DENSITY PER KB (per sequence)
    # ---------------------------------------------------------------------
    def motif_density_per_sequence(self) -> List[float]:
        """
        Calculate motif density (motifs per kb) for each sequence.
        
        Returns:
            List of density values (one per sequence)
        """
        densities = []

        for seq, motifs in self.annotations.items():
            length = self.sequence_lengths.get(seq, 1)
            density = len(motifs) / (length / 1000)
            densities.append(density)

        return densities

    # ---------------------------------------------------------------------
    # POSITIONAL CONSERVATION (equal-length only)
    # ---------------------------------------------------------------------
    def positional_conservation(self, seq_length: int) -> List[int]:
        """
        Calculate positional conservation for equal-length sequences.
        Returns number of sequences containing motif at each position.
        Memory safe: uses windowed aggregation.
        
        Args:
            seq_length: Length of sequences (must be equal for all)
            
        Returns:
            List of conservation counts per window (1000bp windows)
        """
        window = 1000
        bins = math.ceil(seq_length / window)
        conservation = [0] * bins

        for seq, motifs in self.annotations.items():
            covered_bins = set()

            for motif in motifs:
                start_bin = motif["Start"] // window
                end_bin = motif["End"] // window
                for b in range(start_bin, min(end_bin + 1, bins)):
                    covered_bins.add(b)

            for b in covered_bins:
                conservation[b] += 1

        return conservation
