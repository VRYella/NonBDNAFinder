#!/usr/bin/env python3
"""
Integration module to bridge the hyperscan orchestrator with the Streamlit app.

This module provides a simplified interface that matches the expected signature
of all_motifs_refactored while using the hyperscan-based pipeline.
"""

import logging
import tempfile
import os
from typing import List, Dict, Any, Optional
import pandas as pd

from orchestrator import run_pipeline
from utils import parse_fasta

logger = logging.getLogger(__name__)


def convert_candidate_to_motif_dict(candidate, sequence_name: str) -> Dict[str, Any]:
    """
    Convert a Candidate object to the dictionary format expected by the Streamlit app.
    
    Args:
        candidate: Candidate object from the orchestrator
        sequence_name: Name of the sequence
        
    Returns:
        Dictionary in the format expected by the app.py visualization system
    """
    # Map internal class names to official 11 Non-B DNA class names
    class_mapping = {
        'curved_dna': 'Curved DNA',
        'slipped_dna': 'Slipped DNA', 
        'cruciform': 'Cruciform DNA',
        'r_loop': 'R-loop',
        'triplex': 'Triplex',
        'g_quadruplex': 'G-Quadruplex Family',
        'i_motif': 'i-motif family',
        'z_dna': 'Z-DNA',
        'a_philic': 'A-philic DNA',
        'hybrid': 'Hybrid',
        'cluster': 'Non-B DNA cluster regions'
    }
    
    # Get the original class name from candidate
    original_class = candidate.get('Class', '')
    mapped_class = class_mapping.get(original_class, original_class)
    
    return {
        'Sequence Name': sequence_name,
        'Class': mapped_class,
        'Subclass': candidate.get('Subclass', ''),
        'Start': candidate.get('Start', 0),
        'End': candidate.get('End', 0),
        'Length': candidate.get('Length', 0),
        'Normalized_Score': candidate.get('Normalized_Score', 0.0),
        'Actual_Score': candidate.get('Actual_Score', 0.0),
        'Score': candidate.get('Actual_Score', 0.0),  # Alias for compatibility
        'GC Content': candidate.get('GC_Content', 0.0),
        'Sequence': candidate.get('Sequence', ''),
        'Motif_ID': candidate.get('Motif_ID', ''),
        'Scoring_Method': candidate.get('Scoring_Method', ''),
        'S.No': candidate.get('S.No', 0)
    }


def all_motifs_refactored(
    sequence: str,
    sequence_name: str = "sequence",
    nonoverlap: bool = False,
    report_hotspots: bool = True,
    calculate_conservation: bool = False
) -> List[Dict[str, Any]]:
    """
    Hyperscan-based motif detection function that replaces the stub implementation.
    
    This function maintains compatibility with the Streamlit app interface while
    using the high-performance hyperscan pipeline.
    
    Args:
        sequence: DNA sequence to analyze
        sequence_name: Name of the sequence
        nonoverlap: Whether to remove overlapping motifs (not used in hyperscan pipeline)
        report_hotspots: Whether to detect hotspot regions
        calculate_conservation: Whether to calculate conservation scores (not implemented)
        
    Returns:
        List of motif dictionaries compatible with app.py visualization
    """
    try:
        # Clean the sequence
        cleaned_sequence = parse_fasta(sequence)
        
        if not cleaned_sequence:
            logger.warning("Empty sequence provided")
            return []
        
        # Create temporary FASTA file for the orchestrator
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp_file:
            tmp_file.write(f">{sequence_name}\n{cleaned_sequence}\n")
            tmp_file_path = tmp_file.name
        
        try:
            # Run the hyperscan pipeline
            output_prefix = f"/tmp/nbdfinder_{sequence_name}"
            
            # Use all detector classes for comprehensive analysis (11 classes)
            detector_classes = [
                'g_quadruplex', 'curved_dna', 'slipped_dna', 'cruciform',
                'r_loop', 'triplex', 'i_motif', 'z_dna', 'a_philic', 'hybrid', 'cluster'
            ]
            
            output_files = run_pipeline(
                fasta_path=tmp_file_path,
                output_prefix=output_prefix,
                max_workers=2,  # Conservative for Streamlit
                chunk_size=50000,
                detector_classes=detector_classes
            )
            
            # Read the results CSV
            csv_file = output_files.get('csv')
            if not csv_file or not os.path.exists(csv_file):
                logger.warning("No CSV output file found")
                return []
            
            # Load results DataFrame
            results_df = pd.read_csv(csv_file)
            
            # Convert to the format expected by the Streamlit app
            motifs = []
            for _, row in results_df.iterrows():
                motif_dict = convert_candidate_to_motif_dict(row.to_dict(), sequence_name)
                motifs.append(motif_dict)
            
            logger.info(f"Successfully detected {len(motifs)} motifs using hyperscan pipeline")
            return motifs
            
        finally:
            # Cleanup temporary files
            try:
                os.unlink(tmp_file_path)
                # Cleanup output files
                for file_path in output_files.values():
                    if os.path.exists(file_path):
                        os.unlink(file_path)
            except Exception as e:
                logger.warning(f"Failed to cleanup temporary files: {e}")
        
    except Exception as e:
        logger.error(f"Hyperscan motif detection failed: {e}")
        # Return empty list on error to prevent app crashes
        return []


# Alias for compatibility
all_motifs = all_motifs_refactored