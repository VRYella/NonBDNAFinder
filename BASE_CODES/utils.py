#!/usr/bin/env python3
"""
Utility functions for NonBDNAFinder
===================================

Common utility functions for sequence processing and analysis.
"""

import re
from typing import Dict, List, Any, Optional
from collections import Counter


def parse_fasta(sequence: str) -> str:
    """
    Parse a FASTA sequence by removing headers and whitespace.
    
    Args:
        sequence: Input sequence (may contain FASTA header or just DNA sequence)
        
    Returns:
        Cleaned DNA sequence string in uppercase
    """
    if not sequence:
        return ""
    
    # Remove FASTA headers (lines starting with >)
    lines = []
    for line in sequence.strip().split('\n'):
        line = line.strip()
        if line and not line.startswith('>'):
            lines.append(line)
    
    # Join all sequence lines and clean up
    cleaned_seq = ''.join(lines)
    
    # Remove any remaining whitespace and convert to uppercase
    cleaned_seq = re.sub(r'\s+', '', cleaned_seq).upper()
    
    # Validate DNA sequence (only ATCGN allowed)
    if cleaned_seq and not re.match(r'^[ATCGN]+$', cleaned_seq):
        # Remove invalid characters, keeping only valid DNA bases
        cleaned_seq = re.sub(r'[^ATCGN]', '', cleaned_seq)
    
    return cleaned_seq


def gc_content(sequence: str) -> float:
    """
    Calculate GC content percentage of a DNA sequence.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        GC content as percentage (0-100)
    """
    if not sequence:
        return 0.0
    
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total_count = len(sequence)
    
    if total_count == 0:
        return 0.0
    
    return (gc_count / total_count) * 100


def reverse_complement(sequence: str) -> str:
    """
    Calculate reverse complement of a DNA sequence.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        Reverse complement sequence
    """
    if not sequence:
        return ""
    
    complement_map = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'N': 'N', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'n': 'n'
    }
    
    sequence = sequence.upper()
    complement = ''.join(complement_map.get(base, base) for base in sequence)
    
    return complement[::-1]  # Reverse the complement


def wrap(text: str, width: int = 80) -> str:
    """
    Wrap text to specified width for display.
    
    Args:
        text: Text to wrap
        width: Maximum line width (default: 80)
        
    Returns:
        Text with line breaks inserted
    """
    if not text:
        return ""
    
    if len(text) <= width:
        return text
    
    lines = []
    for i in range(0, len(text), width):
        lines.append(text[i:i + width])
    
    return '\n'.join(lines)


def get_basic_stats(sequence: str, motifs: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
    """
    Calculate basic statistics for a DNA sequence.
    
    Args:
        sequence: DNA sequence string
        motifs: Optional list of detected motifs
        
    Returns:
        Dictionary containing sequence statistics
    """
    if not sequence:
        return {
            'Length': 0,
            'GC%': 0.0,
            'AT%': 0.0,
            'A': 0,
            'T': 0,
            'G': 0,
            'C': 0,
            'N': 0,
            'GC_Content': 0.0,
            'Motif Coverage %': 0.0
        }
    
    sequence = sequence.upper()
    length = len(sequence)
    
    # Count nucleotides
    base_counts = Counter(sequence)
    a_count = base_counts.get('A', 0)
    t_count = base_counts.get('T', 0)
    g_count = base_counts.get('G', 0)
    c_count = base_counts.get('C', 0)
    n_count = base_counts.get('N', 0)
    
    # Calculate percentages
    gc_percent = ((g_count + c_count) / length * 100) if length > 0 else 0.0
    at_percent = ((a_count + t_count) / length * 100) if length > 0 else 0.0
    
    # Calculate motif coverage if motifs provided
    motif_coverage = 0.0
    if motifs and length > 0:
        covered_positions = set()
        for motif in motifs:
            start = motif.get('Start', 0)
            end = motif.get('End', 0)
            if isinstance(start, int) and isinstance(end, int) and start > 0 and end > 0:
                # Convert to 0-based indexing for coverage calculation
                for pos in range(start - 1, end):
                    covered_positions.add(pos)
        
        motif_coverage = (len(covered_positions) / length * 100) if length > 0 else 0.0
    
    return {
        'Length': length,
        'GC%': round(gc_percent, 2),
        'AT%': round(at_percent, 2),
        'A': a_count,
        'T': t_count,
        'G': g_count,
        'C': c_count,
        'N': n_count,
        'GC_Content': round(gc_percent, 2),  # Alternative name for compatibility
        'Motif Coverage %': round(motif_coverage, 2)
    }


def validate_dna_sequence(sequence: str) -> bool:
    """
    Validate if a string is a valid DNA sequence.
    
    Args:
        sequence: String to validate
        
    Returns:
        True if valid DNA sequence, False otherwise
    """
    if not sequence:
        return False
    
    sequence = sequence.upper().strip()
    return bool(re.match(r'^[ATCGN]+$', sequence))


def format_sequence_name(name: str) -> str:
    """
    Format sequence name for display and file naming.
    
    Args:
        name: Raw sequence name
        
    Returns:
        Formatted sequence name
    """
    if not name:
        return "Unnamed_Sequence"
    
    # Remove problematic characters for file naming
    formatted = re.sub(r'[^\w\-_\.]', '_', name.strip())
    
    # Remove multiple consecutive underscores
    formatted = re.sub(r'_+', '_', formatted)
    
    # Remove leading/trailing underscores
    formatted = formatted.strip('_')
    
    return formatted if formatted else "Unnamed_Sequence"