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


def export_detailed_motif_table(motifs: List[Dict[str, Any]], 
                               coverage_metrics: Dict[str, float] = None,
                               format_type: str = "excel") -> bytes:
    """
    Export comprehensive detailed motif table with all analysis metrics
    
    Args:
        motifs: List of motif dictionaries
        coverage_metrics: Optional coverage and density metrics
        format_type: Export format ("excel", "csv", "tsv")
    
    Returns:
        File content as bytes for download
    """
    import pandas as pd
    from io import BytesIO, StringIO
    
    if not motifs:
        if format_type == "excel":
            return BytesIO().getvalue()
        else:
            return StringIO().getvalue().encode()
    
    # Convert motifs to DataFrame
    df = pd.DataFrame(motifs)
    
    # Ensure all required columns exist
    required_columns = [
        'Sequence Name', 'Class', 'Subclass', 'Start', 'End', 'Length',
        'Normalized_Score', 'Actual_Score', 'GC Content', 'Sequence',
        'Motif_ID', 'Scoring_Method'
    ]
    
    for col in required_columns:
        if col not in df.columns:
            df[col] = 'N/A'
    
    # Add derived columns for detailed analysis
    df['Center_Position'] = (df['Start'] + df['End']) / 2
    df['Score_Rank'] = df['Normalized_Score'].rank(ascending=False, method='min')
    df['Length_Category'] = pd.cut(df['Length'], bins=[0, 20, 50, 100, 500, float('inf')], 
                                  labels=['Very Short', 'Short', 'Medium', 'Long', 'Very Long'])
    
    # Score categories
    score_quantiles = df['Normalized_Score'].quantile([0.25, 0.5, 0.75])
    df['Score_Category'] = pd.cut(df['Normalized_Score'], 
                                 bins=[0, score_quantiles[0.25], score_quantiles[0.5], 
                                       score_quantiles[0.75], 1.0],
                                 labels=['Low', 'Medium', 'High', 'Very High'],
                                 include_lowest=True)
    
    if format_type == "excel":
        output = BytesIO()
        
        with pd.ExcelWriter(output, engine='openpyxl') as writer:
            # Main detailed motifs table
            df.to_excel(writer, sheet_name='Detailed_Motif_Table', index=False)
            
            # Summary statistics
            summary_stats = {
                'Total_Motifs': len(df),
                'Unique_Classes': df['Class'].nunique(),
                'Unique_Subclasses': df['Subclass'].nunique(),
                'Mean_Length': df['Length'].mean(),
                'Median_Length': df['Length'].median(),
                'Mean_Score': df['Normalized_Score'].mean(),
                'Median_Score': df['Normalized_Score'].median(),
                'Max_Score': df['Normalized_Score'].max(),
                'Min_Score': df['Normalized_Score'].min(),
                'Score_StdDev': df['Normalized_Score'].std(),
                'Total_Sequence_Coverage': df['Length'].sum()
            }
            
            # Add coverage metrics if provided
            if coverage_metrics:
                summary_stats.update(coverage_metrics)
            
            summary_df = pd.DataFrame.from_dict(summary_stats, orient='index', columns=['Value'])
            summary_df.to_excel(writer, sheet_name='Summary_Statistics')
            
            # Class-specific sheets
            for class_name in df['Class'].unique():
                if pd.isna(class_name):
                    continue
                class_df = df[df['Class'] == class_name]
                sheet_name = str(class_name).replace('_', ' ').title()[:31]
                class_df.to_excel(writer, sheet_name=sheet_name, index=False)
            
            # Subclass distribution
            subclass_counts = df['Subclass'].value_counts().reset_index()
            subclass_counts.columns = ['Subclass', 'Count']
            subclass_counts['Percentage'] = (subclass_counts['Count'] / len(df)) * 100
            subclass_counts.to_excel(writer, sheet_name='Subclass_Distribution', index=False)
            
            # Position-based analysis (binned)
            if 'Start' in df.columns and 'End' in df.columns:
                max_pos = df['End'].max()
                bin_size = max(1000, max_pos // 100)  # Create ~100 bins
                df['Position_Bin'] = pd.cut(df['Center_Position'], bins=range(0, int(max_pos) + bin_size, bin_size))
                
                position_analysis = df.groupby('Position_Bin').agg({
                    'Class': 'count',
                    'Length': ['mean', 'sum'],
                    'Normalized_Score': 'mean'
                }).round(3)
                
                position_analysis.columns = ['Motif_Count', 'Mean_Length', 'Total_Length', 'Mean_Score']
                position_analysis = position_analysis.reset_index()
                position_analysis.to_excel(writer, sheet_name='Position_Analysis', index=False)
        
        output.seek(0)
        return output.getvalue()
    
    elif format_type == "csv":
        output = StringIO()
        df.to_csv(output, index=False)
        return output.getvalue().encode()
    
    elif format_type == "tsv":
        output = StringIO()
        df.to_csv(output, sep='\t', index=False)
        return output.getvalue().encode()
    
    else:
        raise ValueError(f"Unsupported format: {format_type}")


def create_comprehensive_analysis_report(motifs: List[Dict[str, Any]], 
                                        sequence_length: int = None,
                                        sequence_name: str = "sequence") -> bytes:
    """
    Create a comprehensive analysis report with multiple data exports
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total sequence length
        sequence_name: Name of the sequence
    
    Returns:
        ZIP file content as bytes containing multiple analysis files
    """
    import zipfile
    from io import BytesIO
    
    zip_buffer = BytesIO()
    
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        # 1. Detailed motif table (Excel)
        detailed_table = export_detailed_motif_table(motifs, format_type="excel")
        zip_file.writestr(f"{sequence_name}_detailed_motif_table.xlsx", detailed_table)
        
        # 2. BED format for genome browsers
        bed_content = export_to_bed(motifs, sequence_name)
        zip_file.writestr(f"{sequence_name}_motifs.bed", bed_content)
        
        # 3. GFF3 format for detailed annotations
        gff3_content = export_to_gff3(motifs, sequence_name)
        zip_file.writestr(f"{sequence_name}_motifs.gff3", gff3_content)
        
        # 4. Class-specific BED files
        class_tracks = export_class_specific_tracks(motifs, sequence_name)
        for class_name, bed_data in class_tracks.items():
            safe_class_name = class_name.replace('_', '-').replace(' ', '-')
            zip_file.writestr(f"{sequence_name}_{safe_class_name}_motifs.bed", bed_data)
        
        # 5. Coverage and density analysis
        if sequence_length:
            from viz_tools import calculate_coverage_density_metrics
            coverage_metrics = calculate_coverage_density_metrics(
                pd.DataFrame(motifs), sequence_length
            )
            
            # Save coverage metrics as JSON
            import json
            coverage_json = json.dumps(coverage_metrics, indent=2)
            zip_file.writestr(f"{sequence_name}_coverage_metrics.json", coverage_json)
            
            # Create density bedgraph
            density_bg = create_density_bedgraph(motifs, sequence_length, sequence_name)
            zip_file.writestr(f"{sequence_name}_density.bedgraph", density_bg)
        
        # 6. Summary report (text)
        summary_report = create_text_summary_report(motifs, sequence_length, sequence_name)
        zip_file.writestr(f"{sequence_name}_summary_report.txt", summary_report)
    
    zip_buffer.seek(0)
    return zip_buffer.getvalue()


def create_text_summary_report(motifs: List[Dict[str, Any]], 
                              sequence_length: int = None,
                              sequence_name: str = "sequence") -> str:
    """
    Create a human-readable text summary report
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total sequence length
        sequence_name: Name of the sequence
    
    Returns:
        Formatted text report as string
    """
    if not motifs:
        return f"NBDFinder Analysis Report for {sequence_name}\n\nNo motifs detected."
    
    df = pd.DataFrame(motifs)
    
    report_lines = [
        f"NBDFinder Analysis Report",
        f"=" * 50,
        f"Sequence: {sequence_name}",
        f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"",
        f"SUMMARY STATISTICS",
        f"-" * 20,
        f"Total motifs detected: {len(df)}",
        f"Unique motif classes: {df['Class'].nunique()}",
        f"Unique subclasses: {df['Subclass'].nunique()}",
    ]
    
    if sequence_length:
        coverage_metrics = calculate_coverage_density_metrics(df, sequence_length)
        report_lines.extend([
            f"Sequence length: {sequence_length:,} bp",
            f"Total bases covered: {coverage_metrics.get('bases_covered', 0):,} bp",
            f"Coverage percentage: {coverage_metrics.get('coverage_percentage', 0):.2f}%",
            f"Motif density: {coverage_metrics.get('motif_density_per_kb', 0):.2f} motifs/kb",
        ])
    
    report_lines.extend([
        f"Mean motif length: {df['Length'].mean():.1f} bp",
        f"Median motif length: {df['Length'].median():.1f} bp",
        f"Mean normalized score: {df['Normalized_Score'].mean():.3f}",
        f"",
        f"CLASS DISTRIBUTION",
        f"-" * 20,
    ])
    
    # Class distribution
    class_counts = df['Class'].value_counts()
    for class_name, count in class_counts.items():
        percentage = (count / len(df)) * 100
        report_lines.append(f"{class_name}: {count} motifs ({percentage:.1f}%)")
    
    report_lines.extend([
        f"",
        f"SUBCLASS DISTRIBUTION",
        f"-" * 20,
    ])
    
    # Subclass distribution
    subclass_counts = df['Subclass'].value_counts()
    for subclass, count in subclass_counts.head(10).items():  # Top 10
        percentage = (count / len(df)) * 100
        report_lines.append(f"{subclass}: {count} motifs ({percentage:.1f}%)")
    
    if len(subclass_counts) > 10:
        report_lines.append(f"... and {len(subclass_counts) - 10} more subclasses")
    
    report_lines.extend([
        f"",
        f"SCORE STATISTICS",
        f"-" * 20,
        f"Highest score: {df['Normalized_Score'].max():.3f}",
        f"Lowest score: {df['Normalized_Score'].min():.3f}",
        f"Score standard deviation: {df['Normalized_Score'].std():.3f}",
        f"",
        f"TOP 10 HIGHEST SCORING MOTIFS",
        f"-" * 20,
    ])
    
    # Top scoring motifs
    top_motifs = df.nlargest(10, 'Normalized_Score')
    for _, motif in top_motifs.iterrows():
        report_lines.append(
            f"{motif['Class']} ({motif['Subclass']}) at {motif['Start']}-{motif['End']}: "
            f"Score {motif['Normalized_Score']:.3f}"
        )
    
    return "\n".join(report_lines)


# Helper function to calculate coverage density metrics (if not imported from viz_tools)
def calculate_coverage_density_metrics(df: pd.DataFrame, sequence_length: int = None) -> Dict[str, float]:
    """
    Calculate comprehensive coverage and density metrics
    (Fallback implementation if viz_tools import fails)
    """
    if df.empty:
        return {}
    
    metrics = {
        'total_motifs': len(df),
        'unique_classes': df['Class'].nunique(),
        'unique_subclasses': df['Subclass'].nunique(),
        'mean_motif_length': df['Length'].mean(),
        'median_motif_length': df['Length'].median(),
        'total_motif_bases': df['Length'].sum(),
    }
    
    if 'Normalized_Score' in df.columns:
        metrics.update({
            'mean_score': df['Normalized_Score'].mean(),
            'median_score': df['Normalized_Score'].median(),
            'max_score': df['Normalized_Score'].max(),
            'min_score': df['Normalized_Score'].min(),
        })
    
    if sequence_length:
        metrics['sequence_length'] = sequence_length
        
        # Calculate coverage
        covered_positions = set()
        for _, motif in df.iterrows():
            start, end = int(motif['Start']), int(motif['End'])
            covered_positions.update(range(start, end + 1))
        
        bases_covered = len(covered_positions)
        metrics.update({
            'bases_covered': bases_covered,
            'coverage_fraction': bases_covered / sequence_length,
            'coverage_percentage': (bases_covered / sequence_length) * 100,
            'motif_density_per_kb': (len(df) / sequence_length) * 1000,
        })
    
    return metrics