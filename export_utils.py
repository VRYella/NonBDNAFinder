#!/usr/bin/env python3
"""
NBDFinder Export Utilities
===========================

Utilities for exporting NBDFinder results to various genome browser formats:
- BED format for UCSC Genome Browser and IGV
- BigWig format for quantitative track data
- GFF3 format for detailed annotations

Author: Enhanced NBDFinder by Dr. Venkata Rajesh Yella
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Any, Optional
import io
import tempfile
import os
from collections import defaultdict

def export_to_bed(motifs: List[Dict[str, Any]], 
                  sequence_name: str = "sequence",
                  score_type: str = "normalized",
                  include_subclass: bool = True) -> str:
    """
    Export motifs to BED format for genome browsers.
    
    Args:
        motifs: List of motif dictionaries
        sequence_name: Name to use as chromosome/sequence identifier
        score_type: "normalized" or "actual" for score field
        include_subclass: Whether to include subclass in name field
    
    Returns:
        BED format string
    """
    bed_lines = []
    
    # BED header
    bed_lines.append(f'track name="NBDFinder_Motifs" description="Non-B DNA Motifs" itemRgb="On"')
    
    # Color mapping for motif classes
    class_colors = {
        "Curved_DNA": "255,154,162",      # Pink
        "Slipped_DNA": "255,218,193",     # Light orange  
        "Cruciform_DNA": "226,240,203",   # Light green
        "R-Loop": "255,211,182",          # Peach
        "Triplex_DNA": "181,234,215",     # Mint
        "G-Quadruplex": "162,215,216",    # Light blue
        "i-Motif": "176,196,222",         # Light steel blue
        "Z-DNA": "255,183,178",           # Light salmon
        "Hybrid": "193,161,146",          # Brown
        "Cluster": "162,200,204"          # Light cyan
    }
    
    for i, motif in enumerate(motifs):
        chrom = sequence_name
        start = int(motif.get('Start', 1)) - 1  # Convert to 0-based for BED
        end = int(motif.get('End', start + 1))
        
        # Create name field
        motif_class = motif.get('Class', 'Unknown')
        subclass = motif.get('Subclass', '')
        if include_subclass and subclass:
            name = f"{motif_class}_{subclass}_{i+1}"
        else:
            name = f"{motif_class}_{i+1}"
        
        # Score (BED format expects 0-1000)
        if score_type == "normalized":
            score = int(float(motif.get('Normalized_Score', 0)) * 1000)
        else:
            actual_score = float(motif.get('Actual_Score', 0))
            # Normalize actual score to 0-1000 range (approximate)
            score = min(1000, max(0, int(actual_score * 100)))
        
        # Strand (always + for motifs)
        strand = "+"
        
        # Color
        color = class_colors.get(motif_class, "128,128,128")
        
        # BED line: chrom start end name score strand thickStart thickEnd itemRgb
        bed_line = f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{start}\t{end}\t{color}"
        bed_lines.append(bed_line)
    
    return "\n".join(bed_lines)

def export_to_gff3(motifs: List[Dict[str, Any]], 
                   sequence_name: str = "sequence",
                   source: str = "NBDFinder") -> str:
    """
    Export motifs to GFF3 format for detailed annotations.
    
    Args:
        motifs: List of motif dictionaries
        sequence_name: Sequence/chromosome name
        source: Source program name
    
    Returns:
        GFF3 format string
    """
    gff_lines = []
    
    # GFF3 header
    gff_lines.append("##gff-version 3")
    gff_lines.append(f"##sequence-region {sequence_name} 1 {max([int(m.get('End', 0)) for m in motifs], default=1000)}")
    
    for i, motif in enumerate(motifs):
        seqid = sequence_name
        feature_type = "non_B_DNA_motif"
        start = motif.get('Start', 1)
        end = motif.get('End', start)
        
        # Score
        score = motif.get('Normalized_Score', '.')
        if score != '.':
            score = f"{float(score):.3f}"
        
        strand = "+"
        phase = "."
        
        # Attributes
        attributes = []
        attributes.append(f"ID=motif_{i+1}")
        attributes.append(f"Name={motif.get('Class', 'Unknown')}")
        attributes.append(f"motif_class={motif.get('Class', 'Unknown')}")
        attributes.append(f"subclass={motif.get('Subclass', '')}")
        attributes.append(f"motif_id={motif.get('Motif_ID', f'motif_{i+1}')}")
        attributes.append(f"actual_score={motif.get('Actual_Score', 0)}")
        attributes.append(f"scoring_method={motif.get('Scoring_Method', 'Unknown')}")
        attributes.append(f"gc_content={motif.get('GC_Content', 0)}")
        attributes.append(f"length={motif.get('Length', end-start+1)}")
        
        if motif.get('Overlap_Classes'):
            attributes.append(f"overlaps={motif.get('Overlap_Classes', '')}")
        
        attributes_str = ";".join(attributes)
        
        # GFF3 line
        gff_line = f"{seqid}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes_str}"
        gff_lines.append(gff_line)
    
    return "\n".join(gff_lines)

def create_density_bedgraph(motifs: List[Dict[str, Any]], 
                           sequence_length: int,
                           sequence_name: str = "sequence",
                           window_size: int = 100) -> str:
    """
    Create a bedGraph format for motif density visualization.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total sequence length
        sequence_name: Sequence/chromosome name
        window_size: Window size for density calculation
    
    Returns:
        bedGraph format string
    """
    bedgraph_lines = []
    
    # Header
    bedgraph_lines.append(f'track type=bedGraph name="NBDFinder_Density" description="Non-B DNA Motif Density"')
    
    # Calculate density in windows
    density_array = np.zeros(sequence_length)
    
    # Add motif coverage
    for motif in motifs:
        start = max(0, int(motif.get('Start', 1)) - 1)  # 0-based
        end = min(sequence_length, int(motif.get('End', start + 1)))
        score = float(motif.get('Normalized_Score', 1))
        
        # Add weighted coverage
        density_array[start:end] += score
    
    # Smooth with sliding window
    smoothed_density = np.convolve(density_array, np.ones(window_size)/window_size, mode='same')
    
    # Convert to bedGraph format (merge adjacent identical values)
    current_value = None
    current_start = 0
    
    for i, value in enumerate(smoothed_density):
        if value != current_value:
            if current_value is not None and current_value > 0:
                bedgraph_lines.append(f"{sequence_name}\t{current_start}\t{i}\t{current_value:.6f}")
            current_value = value
            current_start = i
    
    # Final region
    if current_value is not None and current_value > 0:
        bedgraph_lines.append(f"{sequence_name}\t{current_start}\t{sequence_length}\t{current_value:.6f}")
    
    return "\n".join(bedgraph_lines)

def export_class_specific_tracks(motifs: List[Dict[str, Any]], 
                                sequence_name: str = "sequence") -> Dict[str, str]:
    """
    Create separate BED tracks for each motif class.
    
    Args:
        motifs: List of motif dictionaries
        sequence_name: Sequence/chromosome name
    
    Returns:
        Dictionary mapping class names to BED format strings
    """
    class_motifs = defaultdict(list)
    
    # Group motifs by class
    for motif in motifs:
        motif_class = motif.get('Class', 'Unknown')
        class_motifs[motif_class].append(motif)
    
    # Generate BED for each class
    class_tracks = {}
    for motif_class, class_motif_list in class_motifs.items():
        bed_content = export_to_bed(
            class_motif_list, 
            sequence_name=sequence_name,
            include_subclass=True
        )
        # Update track name
        bed_content = bed_content.replace(
            'track name="NBDFinder_Motifs"',
            f'track name="NBDFinder_{motif_class}" description="{motif_class} Motifs"'
        )
        class_tracks[motif_class] = bed_content
    
    return class_tracks

def create_motif_browser_session(motifs: List[Dict[str, Any]], 
                                sequence_name: str = "sequence",
                                sequence_length: int = None) -> Dict[str, Any]:
    """
    Create a comprehensive browser session with multiple tracks.
    
    Args:
        motifs: List of motif dictionaries
        sequence_name: Sequence/chromosome name
        sequence_length: Total sequence length
    
    Returns:
        Dictionary with all track data for browser loading
    """
    if sequence_length is None:
        sequence_length = max([int(m.get('End', 0)) for m in motifs], default=1000)
    
    session_data = {
        "sequence_info": {
            "name": sequence_name,
            "length": sequence_length
        },
        "tracks": {}
    }
    
    # All motifs track
    session_data["tracks"]["all_motifs"] = {
        "name": "All Non-B DNA Motifs",
        "type": "bed",
        "data": export_to_bed(motifs, sequence_name)
    }
    
    # Class-specific tracks
    class_tracks = export_class_specific_tracks(motifs, sequence_name)
    for class_name, bed_data in class_tracks.items():
        session_data["tracks"][f"class_{class_name}"] = {
            "name": f"{class_name} Motifs",
            "type": "bed", 
            "data": bed_data
        }
    
    # Density track
    session_data["tracks"]["density"] = {
        "name": "Motif Density",
        "type": "bedgraph",
        "data": create_density_bedgraph(motifs, sequence_length, sequence_name)
    }
    
    # GFF3 annotation
    session_data["tracks"]["annotations"] = {
        "name": "Detailed Annotations",
        "type": "gff3",
        "data": export_to_gff3(motifs, sequence_name)
    }
    
    return session_data

def save_browser_files(session_data: Dict[str, Any], output_dir: str = ".") -> Dict[str, str]:
    """
    Save all browser-compatible files to disk.
    
    Args:
        session_data: Browser session data from create_motif_browser_session
        output_dir: Directory to save files
    
    Returns:
        Dictionary mapping file types to file paths
    """
    file_paths = {}
    
    sequence_name = session_data["sequence_info"]["name"]
    
    for track_id, track_info in session_data["tracks"].items():
        file_ext = track_info["type"]
        filename = f"{sequence_name}_{track_id}.{file_ext}"
        filepath = os.path.join(output_dir, filename)
        
        with open(filepath, 'w') as f:
            f.write(track_info["data"])
        
        file_paths[track_id] = filepath
    
    return file_paths