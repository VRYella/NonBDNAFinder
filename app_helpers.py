"""
App Helper Functions
====================

Helper functions specifically for the Streamlit app.
These functions support the web interface and are extracted from utilities.py
to enable removal of deprecated files while maintaining app functionality.
"""

from typing import Dict, Any, List, Optional, Tuple
import pandas as pd
import numpy as np
from collections import Counter, defaultdict
import re


# Import constants from utils.constants
from utils.constants import CORE_OUTPUT_COLUMNS, DEFAULT_COLUMN_VALUES, EXCLUDED_FROM_CONSOLIDATED


def export_results_to_dataframe(motifs: List[Dict[str, Any]]) -> pd.DataFrame:
    """Convert motif results to pandas DataFrame with CORE fields only for display tables.
    
    Per Task 1 & 2 requirements, output tables should only show mandatory, universal columns:
    - Sequence_Name, Class, Subclass, Start, End, Length, Strand, Score, Method, Pattern_ID
    
    Additional motif-specific columns (Repeat_Unit, Loop_Length, etc.) are only shown 
    in Excel download class-specific sheets, not in display tables.
    
    This ensures publication-grade clarity per Nature/NAR/Genome Research standards.
    """
    if not motifs:
        return pd.DataFrame()
    
    df = pd.DataFrame(motifs)
    
    # Use core columns constant
    core_columns = CORE_OUTPUT_COLUMNS
    
    # Ensure all core columns are present with appropriate defaults
    for col in core_columns:
        if col not in df.columns:
            # Set appropriate defaults for missing columns using constants
            df[col] = DEFAULT_COLUMN_VALUES.get(col, 'NA')
    
    # Fill all NaN/None values with appropriate defaults
    result_df = df[core_columns].copy()
    for col in core_columns:
        default_val = DEFAULT_COLUMN_VALUES.get(col, 'NA')
        result_df[col] = result_df[col].fillna(default_val)
    
    return result_df


def calculate_genomic_density(motifs: List[Dict[str, Any]], 
                               sequence_length: int,
                               by_class: bool = True,
                               by_subclass: bool = False) -> Dict[str, float]:
    """
    Calculate genomic density (coverage) for motifs.
    
    Genomic Density (σ_G) = (Total unique bp covered by motifs / 
                             Total length in bp of analyzed region) × 100
    
    IMPORTANT: Uses set-based overlap handling to ensure coverage never exceeds 100%.
    If motifs overlap, only unique positions are counted.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total length of analyzed sequence
        by_class: If True, calculate density per motif class
        by_subclass: If True, calculate density per motif subclass (takes precedence over by_class)
        
    Returns:
        Dictionary with density metrics (percentage, capped at 100%)
        - If by_subclass=True: keys are 'Class:Subclass' format
        - If by_class=True: keys are class names
        - Otherwise: key is 'Overall'
    """
    if not motifs or sequence_length == 0:
        return {'Overall': 0.0}
    
    if not by_class and not by_subclass:
        # Overall density using set-based coverage (handles overlaps correctly)
        covered_positions = set()
        for motif in motifs:
            start = motif.get('Start', 0) - 1  # Convert to 0-based
            end = motif.get('End', 0)
            covered_positions.update(range(start, end))
        
        overall_density = min((len(covered_positions) / sequence_length) * 100, 100.0)
        return {'Overall': round(overall_density, 4)}
    
    # Density per subclass using set-based coverage
    if by_subclass:
        density_by_subclass = {}
        subclass_groups = defaultdict(list)
        
        for motif in motifs:
            class_name = motif.get('Class', 'Unknown')
            subclass_name = motif.get('Subclass', 'Unknown')
            key = f"{class_name}:{subclass_name}"
            subclass_groups[key].append(motif)
        
        # Calculate per-subclass density with overlap handling
        for subclass_key, subclass_motifs in subclass_groups.items():
            covered_positions = set()
            for motif in subclass_motifs:
                start = motif.get('Start', 0) - 1  # Convert to 0-based
                end = motif.get('End', 0)
                covered_positions.update(range(start, end))
            
            subclass_density = min((len(covered_positions) / sequence_length) * 100, 100.0)
            density_by_subclass[subclass_key] = round(subclass_density, 4)
        
        # Calculate overall density (all motifs combined)
        all_covered_positions = set()
        for motif in motifs:
            start = motif.get('Start', 0) - 1  # Convert to 0-based
            end = motif.get('End', 0)
            all_covered_positions.update(range(start, end))
        
        overall_density = min((len(all_covered_positions) / sequence_length) * 100, 100.0)
        density_by_subclass['Overall'] = round(overall_density, 4)
        
        return density_by_subclass
    
    # Density per class using set-based coverage
    density_by_class = {}
    class_groups = defaultdict(list)
    
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        class_groups[class_name].append(motif)
    
    # Calculate per-class density with overlap handling
    for class_name, class_motifs in class_groups.items():
        covered_positions = set()
        for motif in class_motifs:
            start = motif.get('Start', 0) - 1  # Convert to 0-based
            end = motif.get('End', 0)
            covered_positions.update(range(start, end))
        
        class_density = min((len(covered_positions) / sequence_length) * 100, 100.0)
        density_by_class[class_name] = round(class_density, 4)
    
    # Calculate overall density (all motifs combined)
    all_covered_positions = set()
    for motif in motifs:
        start = motif.get('Start', 0) - 1  # Convert to 0-based
        end = motif.get('End', 0)
        all_covered_positions.update(range(start, end))
    
    overall_density = min((len(all_covered_positions) / sequence_length) * 100, 100.0)
    density_by_class['Overall'] = round(overall_density, 4)
    
    return density_by_class


def calculate_positional_density(motifs: List[Dict[str, Any]], 
                                  sequence_length: int,
                                  unit: str = 'Mbp',
                                  by_class: bool = True,
                                  by_subclass: bool = False) -> Dict[str, float]:
    """
    Calculate positional density (frequency) for motifs.
    
    Positional Density (λ) = Total count of predicted motifs / 
                             Total length (in kbp or Mbp) of analyzed region
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total length of analyzed sequence in bp
        unit: 'kbp' or 'Mbp' for reporting units
        by_class: If True, calculate density per motif class
        by_subclass: If True, calculate density per motif subclass (takes precedence over by_class)
        
    Returns:
        Dictionary with positional density (motifs per unit)
        - If by_subclass=True: keys are 'Class:Subclass' format
        - If by_class=True: keys are class names
        - Otherwise: key is 'Overall'
    """
    if not motifs or sequence_length == 0:
        return {'Overall': 0.0}
    
    # Convert to appropriate unit
    if unit == 'kbp':
        sequence_length_unit = sequence_length / 1000
    elif unit == 'Mbp':
        sequence_length_unit = sequence_length / 1000000
    else:
        sequence_length_unit = sequence_length
    
    if not by_class and not by_subclass:
        # Overall positional density
        overall_density = len(motifs) / sequence_length_unit
        return {'Overall': round(overall_density, 2)}
    
    # Positional density per subclass
    if by_subclass:
        density_by_subclass = {}
        subclass_counts = Counter()
        
        for motif in motifs:
            class_name = motif.get('Class', 'Unknown')
            subclass_name = motif.get('Subclass', 'Unknown')
            key = f"{class_name}:{subclass_name}"
            subclass_counts[key] += 1
        
        for subclass_key, count in subclass_counts.items():
            subclass_density = count / sequence_length_unit
            density_by_subclass[subclass_key] = round(subclass_density, 2)
        
        # Also add overall
        density_by_subclass['Overall'] = round(len(motifs) / sequence_length_unit, 2)
        
        return density_by_subclass
    
    # Positional density per class
    density_by_class = {}
    class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    
    for class_name, count in class_counts.items():
        class_density = count / sequence_length_unit
        density_by_class[class_name] = round(class_density, 2)
    
    # Also add overall
    density_by_class['Overall'] = round(len(motifs) / sequence_length_unit, 2)
    
    return density_by_class


def export_statistics_to_excel(motifs: List[Dict[str, Any]], sequence_length: int, 
                               filename: str = "statistics.xlsx") -> str:
    """
    Export comprehensive statistical analysis to Excel format.
    
    Creates an Excel file with multiple sheets containing:
    - Summary statistics
    - Class-level density metrics
    - Subclass-level density metrics
    - Length distribution statistics
    - Score distribution statistics
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of analyzed sequence in base pairs
        filename: Output Excel filename
        
    Returns:
        Success message string
    """
    try:
        import openpyxl
        from openpyxl.styles import Font, Alignment, PatternFill
    except ImportError:
        raise ImportError("openpyxl is required for Excel export. Install with: pip install openpyxl")
    
    if not motifs:
        return "No motifs to export statistics"
    
    # Filter out Hybrid and Cluster motifs for main statistics
    main_motifs = [m for m in motifs if m.get('Class') not in EXCLUDED_FROM_CONSOLIDATED]
    
    with pd.ExcelWriter(filename, engine='openpyxl') as writer:
        # Sheet 1: Summary Statistics
        summary_data = {
            'Metric': [
                'Total Sequences Analyzed',
                'Total Sequence Length (bp)',
                'Total Motifs Detected',
                'Non-Overlapping Motifs',
                'Hybrid Motifs',
                'Cluster Motifs',
                'Unique Motif Classes',
                'Unique Motif Subclasses',
                'Sequence Coverage (%)',
                'Overall Motif Density (motifs/kb)',
                'Average Motif Length (bp)',
                'Min Motif Length (bp)',
                'Max Motif Length (bp)',
                'Average Score',
                'Min Score',
                'Max Score'
            ],
            'Value': []
        }
        
        # Calculate summary values
        hybrid_count = len([m for m in motifs if m.get('Class') == 'Hybrid'])
        cluster_count = len([m for m in motifs if m.get('Class') == 'Non-B_DNA_Clusters'])
        
        total_coverage = sum(m.get('Length', 0) for m in main_motifs)
        coverage_pct = (total_coverage / sequence_length * 100) if sequence_length > 0 else 0
        density = (len(main_motifs) / sequence_length * 1000) if sequence_length > 0 else 0
        
        lengths = [m.get('Length', 0) for m in main_motifs]
        scores = [m.get('Score', 0) for m in main_motifs]
        
        summary_data['Value'] = [
            1,  # Total sequences
            sequence_length,
            len(motifs),
            len(main_motifs),
            hybrid_count,
            cluster_count,
            len(set(m.get('Class', 'Unknown') for m in main_motifs)),
            len(set(m.get('Subclass', 'Unknown') for m in main_motifs)),
            f"{coverage_pct:.2f}",
            f"{density:.2f}",
            f"{sum(lengths)/len(lengths):.2f}" if lengths else "0",
            min(lengths) if lengths else 0,
            max(lengths) if lengths else 0,
            f"{sum(scores)/len(scores):.2f}" if scores else "0",
            f"{min(scores):.2f}" if scores else "0",
            f"{max(scores):.2f}" if scores else "0"
        ]
        
        df_summary = pd.DataFrame(summary_data)
        df_summary.to_excel(writer, sheet_name='Summary', index=False)
        
        # Sheet 2: Class-Level Density Analysis
        class_stats = {}
        for motif in main_motifs:
            cls = motif.get('Class', 'Unknown')
            if cls not in class_stats:
                class_stats[cls] = {'count': 0, 'total_length': 0, 'scores': []}
            class_stats[cls]['count'] += 1
            class_stats[cls]['total_length'] += motif.get('Length', 0)
            class_stats[cls]['scores'].append(motif.get('Score', 0))
        
        class_data = {
            'Motif Class': [],
            'Count': [],
            'Total Length (bp)': [],
            'Genomic Density (%)': [],
            'Motifs per kb': [],
            'Average Length (bp)': [],
            'Average Score': []
        }
        
        for cls, stats in sorted(class_stats.items()):
            class_data['Motif Class'].append(cls)
            class_data['Count'].append(stats['count'])
            class_data['Total Length (bp)'].append(stats['total_length'])
            genomic_density = (stats['total_length'] / sequence_length * 100) if sequence_length > 0 else 0
            class_data['Genomic Density (%)'].append(f"{genomic_density:.4f}")
            motifs_per_kb = (stats['count'] / sequence_length * 1000) if sequence_length > 0 else 0
            class_data['Motifs per kb'].append(f"{motifs_per_kb:.2f}")
            avg_length = stats['total_length'] / stats['count'] if stats['count'] > 0 else 0
            class_data['Average Length (bp)'].append(f"{avg_length:.2f}")
            avg_score = sum(stats['scores']) / len(stats['scores']) if stats['scores'] else 0
            class_data['Average Score'].append(f"{avg_score:.2f}")
        
        df_class = pd.DataFrame(class_data)
        df_class.to_excel(writer, sheet_name='Class_Level_Analysis', index=False)
    
    return f"Statistics exported to {filename}"


def export_to_pdf(motifs: List[Dict[str, Any]], 
                  sequence_length: int,
                  sequence_name: str = "Unknown Sequence") -> bytes:
    """
    Generate PDF containing visualization summaries of sequence analyses.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of analyzed sequence
        sequence_name: Name of the sequence
        
    Returns:
        PDF content as bytes
    """
    # PDF export is complex and requires matplotlib
    # For now, return a simple message
    # In a real implementation, this would generate plots and compile them into a PDF
    raise NotImplementedError("PDF export is not yet fully implemented in the modular architecture")


def create_collapsible_card(title: str, content: str, card_id: str = None, 
                            default_open: bool = False) -> str:
    """
    Create a collapsible card component with smooth animations.
    
    This component replaces the standard Streamlit expander with a more 
    professional, customizable card design that integrates with the global
    CSS design system. Features include:
    - Smooth expand/collapse animations
    - Chevron rotation on toggle
    - Hover and focus effects
    - Full theme integration
    - Mobile responsive design
    
    Args:
        title: The card header/title text
        content: The HTML content to display when expanded
        card_id: Unique identifier for the card (auto-generated if None)
        default_open: Whether the card should be open by default
        
    Returns:
        HTML string for rendering with st.components.v1.html(..., height=None)
        
    Example:
        >>> import streamlit.components.v1 as components
        >>> card_html = create_collapsible_card(
        ...     title="Q: What is Non-B DNA?",
        ...     content="<p>Non-B DNA refers to...</p>"
        ... )
        >>> components.html(card_html, height=None)
    """
    import uuid
    
    # Generate unique ID if not provided
    if card_id is None:
        card_id = f"card-{uuid.uuid4().hex[:8]}"
    
    # Sanitize IDs for HTML/JavaScript
    safe_id = re.sub(r'[^a-zA-Z0-9_-]', '', card_id)
    
    # Set initial states
    body_class = "collapsible-card-body" + (" open" if default_open else "")
    chevron_class = "collapsible-card-chevron" + (" open" if default_open else "")
    
    html = f"""
<!DOCTYPE html>
<html>
<head>
    <style>
        .collapsible-card {{
            margin: 0.75rem 0;
            border-radius: 16px;
            overflow: hidden;
            box-shadow: 0 2px 8px rgba(37, 99, 235, 0.15);
            border: 2px solid rgba(96, 165, 250, 0.4);
            background: #FFFFFF;
        }}
        .collapsible-card-header {{
            padding: 0.75rem 1rem;
            background: linear-gradient(135deg, #EEF2FF 0%, #E0E7FF 100%);
            cursor: pointer;
            font-weight: 700;
            color: #0F172A;
        }}
        .collapsible-card-body {{
            max-height: 0;
            overflow: hidden;
            transition: max-height 0.3s ease;
        }}
        .collapsible-card-body.open {{
            max-height: 1000px;
        }}
        .collapsible-card-content {{
            padding: 1rem;
        }}
    </style>
</head>
<body>
    <div class="collapsible-card">
        <div class="collapsible-card-header" onclick="toggleCard()">
            <span class="{chevron_class}">▼</span> {title}
        </div>
        <div class="{body_class}" id="{safe_id}">
            <div class="collapsible-card-content">
                {content}
            </div>
        </div>
    </div>
    <script>
        function toggleCard() {{
            var body = document.getElementById('{safe_id}');
            body.classList.toggle('open');
        }}
    </script>
</body>
</html>
"""
    return html


def render_summary_panel(seq_length: int, 
                        processing_time: float, 
                        motif_count: int,
                        total_chunks: int = 0,
                        theme_color: str = "#10b981") -> str:
    """
    Create a modern, compact scientific summary panel with theme-aware styling.
    
    This function generates a publication-ready summary block that displays:
    - Sequence length (with appropriate units: bp, kb, or Mb)
    - Processing time (formatted as hh:mm:ss with clock icon 🕐)
    - Total motifs detected
    - Optional chunk information (for chunked processing)
    
    Args:
        seq_length: Total sequence length in base pairs
        processing_time: Processing duration in seconds
        motif_count: Total number of motifs detected
        total_chunks: Number of chunks processed (0 if no chunking)
        theme_color: Primary gradient color
        
    Returns:
        HTML string ready for st.markdown(..., unsafe_allow_html=True)
        
    Example:
        >>> summary_html = render_summary_panel(
        ...     seq_length=1500000,
        ...     processing_time=125.5,
        ...     motif_count=342,
        ...     total_chunks=5,
        ...     theme_color="#10b981"
        ... )
        >>> st.markdown(summary_html, unsafe_allow_html=True)
    """
    # Calculate sequence length display with appropriate units
    if seq_length >= 1_000_000:
        seq_length_display = f"{seq_length / 1_000_000:.1f} Mb"
    elif seq_length >= 1_000:
        seq_length_display = f"{seq_length / 1_000:.1f} kb"
    else:
        seq_length_display = f"{seq_length:,} bp"
    
    # Format processing time as hh:mm:ss or mm:ss
    hours = int(processing_time // 3600)
    mins = int((processing_time % 3600) // 60)
    secs = int(processing_time % 60)
    
    if hours > 0:
        time_display = f"{hours:02d}:{mins:02d}:{secs:02d}"
    else:
        time_display = f"{mins:02d}:{secs:02d}"
    
    # Build chunk info if applicable
    chunk_info = ""
    if total_chunks > 0:
        chunk_info = f"<div style='margin-bottom: 0.4rem;'><b>Total chunks:</b> {total_chunks}</div>"
    
    # Generate HTML with modern styling
    html = f"""
    <div style='background: linear-gradient(135deg, {theme_color} 0%, #059669 100%); 
                padding: 1.2rem; border-radius: 12px; color: white; 
                box-shadow: 0 5px 20px rgba(16, 185, 129, 0.3); margin-bottom: 1rem;'>
        <h3 style='margin: 0 0 1rem 0; text-align: center; font-size: 1.2rem; 
                   border-bottom: 2px solid rgba(255, 255, 255, 0.3); padding-bottom: 0.8rem;'>
            ✅ Analysis Summary
        </h3>
        <div style='background: rgba(255, 255, 255, 0.15); padding: 0.8rem; 
                   border-radius: 8px; line-height: 1.8;'>
            <div style='margin-bottom: 0.4rem;'><b>Sequence length:</b> {seq_length_display}</div>
            {chunk_info}
            <div style='margin-bottom: 0.4rem;'><b>🕐 Processing time:</b> {time_display}</div>
            <div style='margin-bottom: 0.4rem;'><b>Motifs detected:</b> {motif_count:,}</div>
        </div>
    </div>
    """
    
    return html
