"""
Enhanced visualization module for comprehensive motif analysis.

Provides the create_comprehensive_information_based_visualizations function
required by app.py.
"""

import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, Tuple, Any
import logging

# Try to import plotly for interactive plots
try:
    import plotly.graph_objects as go
    import plotly.express as px
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

logger = logging.getLogger(__name__)


def sanitize_numeric_column(series: pd.Series, fill_value: float = 0.0) -> pd.Series:
    """
    Convert a pandas Series to numeric, handling mixed types.
    
    Args:
        series: Input pandas Series that may contain mixed types
        fill_value: Value to use for non-numeric entries
        
    Returns:
        Cleaned pandas Series with only numeric values
    """
    if series.empty:
        return series
    
    # Convert to numeric, coercing errors to NaN
    numeric_series = pd.to_numeric(series, errors='coerce')
    
    # Fill NaN values with fill_value
    numeric_series = numeric_series.fillna(fill_value)
    
    return numeric_series


def create_comprehensive_information_based_visualizations(
    df: pd.DataFrame, 
    sequence_length: int, 
    sequence_name: str
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """
    Create comprehensive information-based visualizations for motif analysis.
    
    Args:
        df: Motif analysis results DataFrame
        sequence_length: Length of the analyzed sequence
        sequence_name: Name of the analyzed sequence
        
    Returns:
        Tuple of (static_plots, interactive_plots, detailed_stats)
    """
    static_plots = {}
    interactive_plots = {}
    detailed_stats = {}
    
    try:
        # Create basic static plots using matplotlib
        if not df.empty:
            
            # Coverage & Density Analysis
            try:
                # Basic motif distribution plot
                fig, ax = plt.subplots(figsize=(10, 6))
                if 'Class' in df.columns:
                    class_counts = df['Class'].value_counts()
                    ax.bar(class_counts.index, class_counts.values)
                    ax.set_title(f'Motif Class Distribution - {sequence_name}')
                    ax.set_xlabel('Motif Class')
                    ax.set_ylabel('Count')
                    plt.xticks(rotation=45, ha='right')
                    plt.tight_layout()
                    static_plots['coverage_analysis'] = fig
                else:
                    # Create empty plot if no Class column
                    ax.text(0.5, 0.5, 'No motif class data available', 
                           ha='center', va='center', transform=ax.transAxes)
                    static_plots['coverage_analysis'] = fig
                    
            except Exception as e:
                logger.warning(f"Could not create coverage analysis: {e}")
            
            # Detailed coverage map
            try:
                fig, ax = plt.subplots(figsize=(12, 4))
                if 'Start' in df.columns and 'End' in df.columns:
                    for i, (_, row) in enumerate(df.iterrows()):
                        start = row.get('Start', 0)
                        end = row.get('End', 0)
                        ax.plot([start, end], [i, i], linewidth=3, alpha=0.7)
                    ax.set_xlabel('Sequence Position (bp)')
                    ax.set_ylabel('Motif Index')
                    ax.set_title(f'Motif Position Map - {sequence_name}')
                    ax.set_xlim(0, sequence_length)
                else:
                    ax.text(0.5, 0.5, 'No position data available', 
                           ha='center', va='center', transform=ax.transAxes)
                plt.tight_layout()
                static_plots['detailed_coverage_map'] = fig
                
            except Exception as e:
                logger.warning(f"Could not create detailed coverage map: {e}")
            
            # Distribution Analysis
            try:
                fig, ax = plt.subplots(figsize=(8, 6))
                if 'Length' in df.columns:
                    lengths = sanitize_numeric_column(df['Length']).dropna()
                    if len(lengths) > 0:
                        ax.hist(lengths, bins=20, alpha=0.7, edgecolor='black')
                        ax.set_xlabel('Motif Length (bp)')
                        ax.set_ylabel('Frequency')
                        ax.set_title(f'Motif Length Distribution - {sequence_name}')
                    else:
                        ax.text(0.5, 0.5, 'No length data available', 
                               ha='center', va='center', transform=ax.transAxes)
                else:
                    ax.text(0.5, 0.5, 'No length data available', 
                           ha='center', va='center', transform=ax.transAxes)
                plt.tight_layout()
                static_plots['distribution_analysis'] = fig
                
            except Exception as e:
                logger.warning(f"Could not create distribution analysis: {e}")
            
            # Sequence Analysis
            try:
                fig, ax = plt.subplots(figsize=(10, 6))
                if 'Normalized_Score' in df.columns:
                    scores = sanitize_numeric_column(df['Normalized_Score']).dropna()
                    if len(scores) > 0:
                        ax.hist(scores, bins=20, alpha=0.7, edgecolor='black')
                        ax.set_xlabel('Normalized Score')
                        ax.set_ylabel('Frequency')
                        ax.set_title(f'Score Distribution - {sequence_name}')
                    else:
                        ax.text(0.5, 0.5, 'No score data available', 
                               ha='center', va='center', transform=ax.transAxes)
                else:
                    ax.text(0.5, 0.5, 'No score data available', 
                           ha='center', va='center', transform=ax.transAxes)
                plt.tight_layout()
                static_plots['sequence_analysis'] = fig
                
            except Exception as e:
                logger.warning(f"Could not create sequence analysis: {e}")
            
            # Comparative Analysis
            try:
                fig, ax = plt.subplots(figsize=(10, 6))
                if 'Class' in df.columns and 'Length' in df.columns:
                    classes = df['Class'].unique()
                    for cls in classes:
                        cls_data = sanitize_numeric_column(df[df['Class'] == cls]['Length']).dropna()
                        if len(cls_data) > 0:
                            ax.scatter([cls] * len(cls_data), cls_data, alpha=0.6, label=cls)
                    ax.set_xlabel('Motif Class')
                    ax.set_ylabel('Length (bp)')
                    ax.set_title(f'Class vs Length Analysis - {sequence_name}')
                    plt.xticks(rotation=45, ha='right')
                    if len(classes) <= 10:  # Only show legend for reasonable number of classes
                        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                else:
                    ax.text(0.5, 0.5, 'No class/length data available', 
                           ha='center', va='center', transform=ax.transAxes)
                plt.tight_layout()
                static_plots['comparative_analysis'] = fig
                
            except Exception as e:
                logger.warning(f"Could not create comparative analysis: {e}")
            
            # Advanced Analysis
            try:
                fig, ax = plt.subplots(figsize=(10, 6))
                if 'Start' in df.columns:
                    positions = df['Start'].dropna()
                    if len(positions) > 0:
                        ax.hist(positions, bins=30, alpha=0.7, edgecolor='black')
                        ax.set_xlabel('Start Position (bp)')
                        ax.set_ylabel('Frequency')
                        ax.set_title(f'Motif Position Distribution - {sequence_name}')
                    else:
                        ax.text(0.5, 0.5, 'No position data available', 
                               ha='center', va='center', transform=ax.transAxes)
                else:
                    ax.text(0.5, 0.5, 'No position data available', 
                           ha='center', va='center', transform=ax.transAxes)
                plt.tight_layout()
                static_plots['advanced_analysis'] = fig
                
            except Exception as e:
                logger.warning(f"Could not create advanced analysis: {e}")
            
            # Generate detailed statistics
            try:
                detailed_stats = {
                    'total_motifs': len(df),
                    'sequence_length': sequence_length,
                    'sequence_name': sequence_name,
                    'coverage_percentage': (len(df) / sequence_length * 100) if sequence_length > 0 else 0
                }
                
                if 'Class' in df.columns:
                    detailed_stats['motif_classes'] = df['Class'].nunique()
                    detailed_stats['class_distribution'] = df['Class'].value_counts().to_dict()
                
                if 'Length' in df.columns:
                    lengths = sanitize_numeric_column(df['Length']).dropna()
                    if len(lengths) > 0:
                        detailed_stats['avg_motif_length'] = lengths.mean()
                        detailed_stats['max_motif_length'] = lengths.max()
                        detailed_stats['min_motif_length'] = lengths.min()
                
                if 'Normalized_Score' in df.columns:
                    scores = sanitize_numeric_column(df['Normalized_Score']).dropna()
                    if len(scores) > 0:
                        detailed_stats['avg_score'] = scores.mean()
                        detailed_stats['max_score'] = scores.max()
                        detailed_stats['min_score'] = scores.min()
                        
            except Exception as e:
                logger.warning(f"Could not create detailed statistics: {e}")
                detailed_stats = {'total_motifs': len(df)}
            
            # Create interactive plots if plotly is available
            if PLOTLY_AVAILABLE:
                try:
                    if 'Class' in df.columns:
                        # Interactive motif count plot
                        class_counts = df['Class'].value_counts()
                        fig = go.Figure(data=[go.Bar(x=class_counts.index, y=class_counts.values)])
                        fig.update_layout(title=f"Interactive Motif Distribution - {sequence_name}",
                                        xaxis_title="Motif Class", yaxis_title="Count")
                        interactive_plots['motif_distribution'] = fig
                    
                    if 'Start' in df.columns and 'End' in df.columns and 'Class' in df.columns:
                        # Interactive position plot
                        fig = go.Figure()
                        for i, (_, row) in enumerate(df.iterrows()):
                            fig.add_trace(go.Scatter(
                                x=[row['Start'], row['End']], 
                                y=[i, i],
                                mode='lines',
                                name=f"{row['Class']} ({i})",
                                showlegend=False
                            ))
                        fig.update_layout(title=f"Interactive Motif Position Map - {sequence_name}",
                                        xaxis_title="Position (bp)", yaxis_title="Motif Index")
                        interactive_plots['position_map'] = fig
                        
                except Exception as e:
                    logger.warning(f"Could not create interactive plots: {e}")
        
        else:
            # Handle empty DataFrame
            logger.info("DataFrame is empty, creating minimal visualizations")
            
            # Create basic static plots
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.5, 0.5, 'No motif data available for visualization', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'No Data - {sequence_name}')
            static_plots['coverage_analysis'] = fig
            
            # Basic stats
            detailed_stats = {
                'total_motifs': 0,
                'sequence_length': sequence_length,
                'sequence_name': sequence_name
            }
    
    except Exception as e:
        logger.error(f"Error in create_comprehensive_information_based_visualizations: {e}")
        # Return minimal empty structure
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, f'Error creating visualizations: {str(e)[:100]}...', 
               ha='center', va='center', transform=ax.transAxes)
        static_plots = {'coverage_analysis': fig}
        interactive_plots = {}
        detailed_stats = {'total_motifs': len(df) if df is not None else 0}
    
    return static_plots, interactive_plots, detailed_stats