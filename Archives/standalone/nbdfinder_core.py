#!/usr/bin/env python3
"""
Standalone NonBDNA Finder - Core Module
=======================================

A simplified, self-contained version of NonBDNAFinder for Jupyter notebook usage.
Includes essential functionality for detecting non-B DNA structures with minimal dependencies.

Author: Dr. Venkata Rajesh Yella
"""

import re
import os
import io
import zipfile
import tempfile
import pandas as pd
import numpy as np
from typing import List, Dict, Any, Tuple, Optional
from dataclasses import dataclass, field
from Bio import SeqIO
from collections import defaultdict, Counter
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Try to import hyperscan for optimized pattern matching
try:
    import hyperscan as hs
    HYPERSCAN_AVAILABLE = True
    logger.info("Hyperscan available for optimized pattern matching")
except ImportError:
    HYPERSCAN_AVAILABLE = False
    logger.warning("Hyperscan not available, falling back to regex")


@dataclass
class MotifCandidate:
    """Data structure for detected motif candidates"""
    sequence_name: str
    chromosome: str
    motif_class: str
    subclass: str
    motif_id: int
    start: int  # 1-based
    end: int    # 1-based
    length: int
    sequence: str
    score: float
    normalized_score: float
    scoring_method: str
    gc_content: float
    overlap_classes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for DataFrame creation"""
        return {
            'S.No': self.motif_id,
            'Sequence_Name': self.sequence_name,
            'Chromosome/Contig': self.chromosome,
            'Class': self.motif_class,
            'Subclass': self.subclass,
            'Motif_ID': f"motif_{self.motif_id}",
            'Start': self.start,
            'End': self.end,
            'Length': self.length,
            'Normalized_Score': self.normalized_score,
            'Actual_Score': self.score,
            'Scoring_Method': self.scoring_method,
            'GC_Content': self.gc_content,
            'Sequence': self.sequence,
            'Overlap_Classes': ','.join(self.overlap_classes) if self.overlap_classes else ''
        }


class NonBDNAPatterns:
    """Pattern definitions for non-B DNA structures"""
    
    def __init__(self):
        self.patterns = {
            # G-Quadruplex patterns
            'g4_canonical': r'G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,}',
            'g4_relaxed': r'G{2,}\w{1,15}G{2,}\w{1,15}G{2,}\w{1,15}G{2,}',
            'g4_bulged': r'G{2,}\w{0,3}G{2,}\w{1,12}G{2,}\w{0,3}G{2,}\w{1,12}G{2,}\w{0,3}G{2,}',
            
            # i-Motif patterns
            'i_motif_canonical': r'C{3,}\w{1,7}C{3,}\w{1,7}C{3,}\w{1,7}C{3,}',
            'i_motif_relaxed': r'C{2,}\w{1,15}C{2,}\w{1,15}C{2,}\w{1,15}C{2,}',
            
            # Z-DNA patterns
            'z_dna_alternating': r'([ATGC][CGT]){6,}',
            'z_dna_extended': r'([ATGC][CGT]){4,}',
            
            # Triplex patterns
            'triplex_purine': r'[AG]{15,}',
            'triplex_pyrimidine': r'[CT]{15,}',
            
            # Cruciform patterns (palindromes)
            'cruciform_short': r'([ATGC]{6,15})[ATGC]{0,8}(?=.*\1)',
            'cruciform_long': r'([ATGC]{10,30})[ATGC]{0,20}(?=.*\1)',
            
            # Curved DNA (A-tract)
            'curved_a_tract': r'A{5,}',
            'curved_phased': r'(A{3,}[ATGC]{5,12}){3,}',
            
            # Slipped DNA (repeats)
            'slipped_dinucleotide': r'([ATGC]{2}){6,}',
            'slipped_trinucleotide': r'([ATGC]{3}){5,}',
            'slipped_microsatellite': r'([ATGC]{1,6}){4,}',
            
            # A-philic DNA
            'a_philic_tract': r'A{4,}[ATGC]{0,5}A{4,}',
            
            # R-loop forming sequences
            'r_loop_g_rich': r'G{4,}[ATGC]{1,10}G{4,}',
        }
    
    def get_pattern_info(self, pattern_name: str) -> Tuple[str, str]:
        """Get class and subclass for a pattern"""
        class_map = {
            'g4_': ('G-Quadruplex Family', 'G4'),
            'i_motif_': ('i-motif family', 'i-motif'),
            'z_dna_': ('Z-DNA', 'Z-DNA'),
            'triplex_': ('Triplex', 'Triplex'),
            'cruciform_': ('Cruciform DNA', 'Palindrome'),
            'curved_': ('Curved DNA', 'A-tract'),
            'slipped_': ('Slipped DNA', 'STR'),
            'a_philic_': ('A-philic DNA', 'A-philic'),
            'r_loop_': ('R-loop', 'R-loop'),
        }
        
        for prefix, (cls, subcls) in class_map.items():
            if pattern_name.startswith(prefix):
                return cls, subcls
        
        return 'Unknown', 'Unknown'


class MotifScorer:
    """Scoring algorithms for different motif types"""
    
    @staticmethod
    def g4hunter_score(sequence: str) -> float:
        """G4Hunter algorithm score"""
        if not sequence:
            return 0.0
        
        score = 0.0
        for i, base in enumerate(sequence.upper()):
            if base == 'G':
                score += 1
            elif base == 'C':
                score -= 1
        
        return score / len(sequence) if len(sequence) > 0 else 0.0
    
    @staticmethod
    def gc_content(sequence: str) -> float:
        """Calculate GC content"""
        if not sequence:
            return 0.0
        
        gc_count = sequence.upper().count('G') + sequence.upper().count('C')
        return gc_count / len(sequence)
    
    @staticmethod
    def z_dna_score(sequence: str) -> float:
        """Z-DNA propensity score based on dinucleotide transitions"""
        if len(sequence) < 2:
            return 0.0
        
        # Z-DNA favoring dinucleotides with scores
        z_scores = {
            'CG': 3.0, 'GC': 3.0,
            'CA': 1.0, 'AC': 1.0,
            'GT': 1.0, 'TG': 1.0,
            'AT': 0.5, 'TA': 0.5,
        }
        
        total_score = 0.0
        count = 0
        
        for i in range(len(sequence) - 1):
            dinuc = sequence[i:i+2].upper()
            if dinuc in z_scores:
                total_score += z_scores[dinuc]
                count += 1
        
        return total_score / count if count > 0 else 0.0
    
    @staticmethod
    def triplex_score(sequence: str) -> float:
        """Triplex stability score"""
        if not sequence:
            return 0.0
        
        # Simple purine/pyrimidine tract scoring
        purine_count = sequence.upper().count('A') + sequence.upper().count('G')
        pyrimidine_count = sequence.upper().count('C') + sequence.upper().count('T')
        
        max_tract = max(purine_count, pyrimidine_count)
        return max_tract / len(sequence) if len(sequence) > 0 else 0.0


class StandaloneNonBDNAFinder:
    """
    Standalone NonBDNA motif finder for Jupyter notebooks
    """
    
    def __init__(self):
        self.patterns = NonBDNAPatterns()
        self.scorer = MotifScorer()
        self.motif_counter = 0
        
        # Compile regex patterns
        self._compile_patterns()
    
    def _compile_patterns(self):
        """Compile regex patterns for motif detection"""
        self.compiled_patterns = {}
        
        for name, pattern in self.patterns.patterns.items():
            try:
                self.compiled_patterns[name] = re.compile(pattern, re.IGNORECASE)
            except re.error as e:
                logger.warning(f"Failed to compile pattern {name}: {e}")
    
    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[MotifCandidate]:
        """
        Detect non-B DNA motifs in a sequence
        
        Args:
            sequence: DNA sequence string
            sequence_name: Name/identifier for the sequence
            
        Returns:
            List of MotifCandidate objects
        """
        motifs = []
        sequence = sequence.upper().replace('N', 'A')  # Handle ambiguous bases
        
        # Detect motifs using compiled patterns
        for pattern_name, compiled_pattern in self.compiled_patterns.items():
            matches = compiled_pattern.finditer(sequence)
            
            for match in matches:
                start = match.start() + 1  # Convert to 1-based
                end = match.end()
                matched_seq = match.group()
                length = len(matched_seq)
                
                # Skip very short matches
                if length < 10:
                    continue
                
                # Get motif class and subclass
                motif_class, subclass = self.patterns.get_pattern_info(pattern_name)
                
                # Calculate scores
                raw_score = self._calculate_score(matched_seq, motif_class)
                normalized_score = min(1.0, max(0.0, raw_score))
                
                # Calculate GC content
                gc_content = self.scorer.gc_content(matched_seq)
                
                # Create motif candidate
                motif = MotifCandidate(
                    sequence_name=sequence_name,
                    chromosome=sequence_name,
                    motif_class=motif_class,
                    subclass=subclass,
                    motif_id=self.motif_counter,
                    start=start,
                    end=end,
                    length=length,
                    sequence=matched_seq,
                    score=raw_score,
                    normalized_score=normalized_score,
                    scoring_method=self._get_scoring_method(motif_class),
                    gc_content=gc_content
                )
                
                motifs.append(motif)
                self.motif_counter += 1
        
        return motifs
    
    def _calculate_score(self, sequence: str, motif_class: str) -> float:
        """Calculate appropriate score based on motif class"""
        if 'G-Quadruplex' in motif_class:
            return abs(self.scorer.g4hunter_score(sequence))
        elif 'Z-DNA' in motif_class:
            return self.scorer.z_dna_score(sequence)
        elif 'Triplex' in motif_class:
            return self.scorer.triplex_score(sequence)
        elif 'i-motif' in motif_class:
            return abs(self.scorer.g4hunter_score(sequence))  # Similar to G4 but for C
        else:
            # Default scoring based on GC content and length
            gc_score = abs(0.5 - self.scorer.gc_content(sequence))  # Deviation from 50%
            length_score = min(1.0, len(sequence) / 100.0)
            return (gc_score + length_score) / 2.0
    
    def _get_scoring_method(self, motif_class: str) -> str:
        """Get scoring method name for motif class"""
        methods = {
            'G-Quadruplex Family': 'G4Hunter',
            'Z-DNA': 'Z-DNA_Score',
            'Triplex': 'Triplex_Score',
            'i-motif family': 'C4Hunter',
        }
        return methods.get(motif_class, 'Generic_Score')
    
    def analyze_fasta(self, fasta_content: str) -> pd.DataFrame:
        """
        Analyze FASTA content and return results as DataFrame
        
        Args:
            fasta_content: FASTA file content as string
            
        Returns:
            pandas DataFrame with motif results
        """
        all_motifs = []
        
        # Parse FASTA content
        fasta_io = io.StringIO(fasta_content)
        
        try:
            for record in SeqIO.parse(fasta_io, 'fasta'):
                sequence_name = record.id
                sequence = str(record.seq)
                
                logger.info(f"Analyzing sequence: {sequence_name} (length: {len(sequence)} bp)")
                
                # Detect motifs in this sequence
                motifs = self.detect_motifs(sequence, sequence_name)
                all_motifs.extend(motifs)
                
                logger.info(f"Found {len(motifs)} motifs in {sequence_name}")
        
        except Exception as e:
            logger.error(f"Error parsing FASTA: {e}")
            raise
        
        # Convert to DataFrame
        if not all_motifs:
            logger.warning("No motifs detected")
            return pd.DataFrame()
        
        motif_dicts = [motif.to_dict() for motif in all_motifs]
        df = pd.DataFrame(motif_dicts)
        
        # Add summary statistics
        logger.info(f"Total motifs detected: {len(df)}")
        if len(df) > 0:
            class_counts = df['Class'].value_counts()
            logger.info(f"Motif classes found: {dict(class_counts)}")
        
        return df
    
    def create_results_package(self, df: pd.DataFrame, filename_prefix: str = "nbdfinder_results") -> bytes:
        """
        Create a zip package with all results files
        
        Args:
            df: Results DataFrame
            filename_prefix: Prefix for output files
            
        Returns:
            Zip file content as bytes
        """
        zip_buffer = io.BytesIO()
        
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zipf:
            # CSV file
            csv_buffer = io.StringIO()
            df.to_csv(csv_buffer, index=False)
            zipf.writestr(f"{filename_prefix}.csv", csv_buffer.getvalue())
            
            # Excel file with multiple sheets
            excel_buffer = io.BytesIO()
            with pd.ExcelWriter(excel_buffer, engine='openpyxl') as writer:
                # All motifs sheet
                df.to_excel(writer, sheet_name='All_Motifs', index=False)
                
                # Summary statistics sheet
                if len(df) > 0:
                    summary_stats = {
                        'Total_Motifs': len(df),
                        'Unique_Classes': df['Class'].nunique(),
                        'Mean_Score': df['Normalized_Score'].mean(),
                        'Max_Score': df['Normalized_Score'].max(),
                        'Mean_Length': df['Length'].mean(),
                        'Total_Coverage': df['Length'].sum(),
                    }
                    summary_df = pd.DataFrame.from_dict(summary_stats, orient='index', columns=['Value'])
                    summary_df.to_excel(writer, sheet_name='Summary')
                
                # Class-specific sheets
                for class_name in df['Class'].unique():
                    if pd.notna(class_name):
                        class_df = df[df['Class'] == class_name]
                        sheet_name = class_name.replace(' ', '_')[:31]  # Excel sheet name limit
                        class_df.to_excel(writer, sheet_name=sheet_name, index=False)
            
            zipf.writestr(f"{filename_prefix}.xlsx", excel_buffer.getvalue())
            
            # GFF3 file
            gff3_content = self._create_gff3(df)
            zipf.writestr(f"{filename_prefix}.gff3", gff3_content)
            
            # Summary report
            report_content = self._create_summary_report(df)
            zipf.writestr(f"{filename_prefix}_summary.txt", report_content)
        
        zip_buffer.seek(0)
        return zip_buffer.getvalue()
    
    def _create_gff3(self, df: pd.DataFrame) -> str:
        """Create GFF3 format output"""
        if len(df) == 0:
            return "##gff-version 3\n# No motifs detected\n"
        
        lines = ["##gff-version 3"]
        lines.append("# NonBDNA Finder Results")
        lines.append("# Columns: seqid source type start end score strand phase attributes")
        
        for _, row in df.iterrows():
            attributes = f"ID=motif_{row['S.No']};Class={row['Class']};Subclass={row['Subclass']};Score={row['Normalized_Score']:.3f}"
            gff_line = f"{row['Chromosome/Contig']}\tNBDFinder\tmotif\t{row['Start']}\t{row['End']}\t{row['Normalized_Score']:.3f}\t.\t.\t{attributes}"
            lines.append(gff_line)
        
        return '\n'.join(lines)
    
    def _create_summary_report(self, df: pd.DataFrame) -> str:
        """Create a text summary report"""
        if len(df) == 0:
            return "NonBDNA Finder Results Summary\n=============================\n\nNo motifs detected in the input sequence(s).\n"
        
        lines = [
            "NonBDNA Finder Results Summary",
            "=" * 30,
            "",
            f"Total motifs detected: {len(df)}",
            f"Unique sequences analyzed: {df['Sequence_Name'].nunique()}",
            f"Unique motif classes: {df['Class'].nunique()}",
            "",
            "Motif Classes Found:",
            "-" * 20
        ]
        
        class_counts = df['Class'].value_counts()
        for class_name, count in class_counts.items():
            lines.append(f"{class_name}: {count} motifs")
        
        lines.extend([
            "",
            "Score Statistics:",
            "-" * 16,
            f"Mean normalized score: {df['Normalized_Score'].mean():.3f}",
            f"Maximum score: {df['Normalized_Score'].max():.3f}",
            f"Minimum score: {df['Normalized_Score'].min():.3f}",
            "",
            "Length Statistics:",
            "-" * 17,
            f"Mean motif length: {df['Length'].mean():.1f} bp",
            f"Total coverage: {df['Length'].sum()} bp",
            f"Longest motif: {df['Length'].max()} bp",
            "",
            "Analysis completed successfully!",
            f"Results package generated with {len(df)} motifs."
        ])
        
        return '\n'.join(lines)


# Convenience function for easy usage
def analyze_fasta_file(fasta_content: str) -> Tuple[pd.DataFrame, bytes]:
    """
    Analyze FASTA content and return results DataFrame and zip package
    
    Args:
        fasta_content: FASTA file content as string
        
    Returns:
        Tuple of (results_dataframe, zip_package_bytes)
    """
    finder = StandaloneNonBDNAFinder()
    df = finder.analyze_fasta(fasta_content)
    zip_package = finder.create_results_package(df)
    
    return df, zip_package


if __name__ == "__main__":
    # Example usage
    sample_fasta = """>example_sequence
GGGTTTTGGGTTTTGGGTTTTGGGAAACCCAAACCCAAACCCAAACCCATATATATATATATGCGCGCGCGCGCGC
AAAAAAAAATGCGTAAAAAAAAATGCGTAAAAAAAAATGCGT"""
    
    df, zip_data = analyze_fasta_file(sample_fasta)
    print(f"Detected {len(df)} motifs")
    print(df.head() if len(df) > 0 else "No motifs found")