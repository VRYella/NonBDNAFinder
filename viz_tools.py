"""
Visualization tools for NBDFinder - plotting helpers for Streamlit integration.

Provides interactive visualizations using plotly for motif analysis results.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Optional, Tuple, Any
import logging

try:
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    logging.warning("Plotly not available for visualizations")

try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False
    logging.warning("NetworkX not available for network plots")

try:
    from sklearn.manifold import TSNE
    from sklearn.preprocessing import StandardScaler
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    logging.warning("Scikit-learn not available for dimensionality reduction")

logger = logging.getLogger(__name__)


def check_dependencies():
    """Check if required dependencies are available"""
    if not PLOTLY_AVAILABLE:
        raise ImportError("plotly is required for visualizations. Install with: pip install plotly")


def plot_motif_counts(df: pd.DataFrame, title: str = "Motif Distribution") -> go.Figure:
    """
    Create bar plot of motif counts by class
    
    Args:
        df: Results DataFrame
        title: Plot title
        
    Returns:
        Plotly Figure object
    """
    check_dependencies()
    
    if df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No data to display", xref="paper", yref="paper",
                          x=0.5, y=0.5, showarrow=False)
        return fig
    
    # Count motifs by class
    counts = df.groupby('Class').size().reset_index(name='Count')
    
    # Create bar plot
    fig = px.bar(
        counts, 
        x='Class', 
        y='Count',
        title=title,
        labels={'Class': 'Motif Class', 'Count': 'Number of Motifs'},
        color='Count',
        color_continuous_scale='viridis'
    )
    
    fig.update_layout(
        xaxis_tickangle=-45,
        height=500,
        showlegend=False
    )
    
    return fig


def plot_motif_lengths(df: pd.DataFrame, title: str = "Motif Length Distribution") -> go.Figure:
    """
    Create histogram of motif lengths by class
    
    Args:
        df: Results DataFrame
        title: Plot title
        
    Returns:
        Plotly Figure object
    """
    check_dependencies()
    
    if df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No data to display", xref="paper", yref="paper",
                          x=0.5, y=0.5, showarrow=False)
        return fig
    
    fig = px.histogram(
        df,
        x='Length',
        color='Class',
        title=title,
        labels={'Length': 'Motif Length (bp)', 'count': 'Frequency'},
        marginal='box',  # Add box plot on top
        nbins=30
    )
    
    fig.update_layout(height=500)
    return fig


def plot_score_distribution(df: pd.DataFrame, score_col: str = 'Normalized_Score',
                           title: str = "Score Distribution") -> go.Figure:
    """
    Create violin plot of score distributions by class
    
    Args:
        df: Results DataFrame
        score_col: Column name for scores
        title: Plot title
        
    Returns:
        Plotly Figure object
    """
    check_dependencies()
    
    if df.empty or score_col not in df.columns:
        fig = go.Figure()
        fig.add_annotation(text="No data to display", xref="paper", yref="paper",
                          x=0.5, y=0.5, showarrow=False)
        return fig
    
    fig = px.violin(
        df,
        x='Class',
        y=score_col,
        title=title,
        box=True,  # Add box plot inside violin
        labels={'Class': 'Motif Class', score_col: 'Score'}
    )
    
    fig.update_layout(
        xaxis_tickangle=-45,
        height=500
    )
    
    return fig


def plot_density_track(df: pd.DataFrame, sequence_name: str, seq_len: int, 
                      window_size: int = 1000, title: str = "Motif Density Track") -> go.Figure:
    """
    Create density heatmap showing motif distribution along sequence
    
    Args:
        df: Results DataFrame
        sequence_name: Name of sequence to plot
        seq_len: Total sequence length
        window_size: Window size for density calculation
        title: Plot title
        
    Returns:
        Plotly Figure object
    """
    check_dependencies()
    
    # Filter for specific sequence
    seq_df = df[df['Sequence_Name'] == sequence_name]
    
    if seq_df.empty:
        fig = go.Figure()
        fig.add_annotation(text=f"No motifs found for sequence: {sequence_name}",
                          xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
        return fig
    
    # Create windows
    windows = np.arange(0, seq_len, window_size)
    classes = seq_df['Class'].unique()
    
    # Calculate density matrix
    density_matrix = []
    for class_name in classes:
        class_df = seq_df[seq_df['Class'] == class_name]
        densities = []
        
        for window_start in windows:
            window_end = window_start + window_size
            # Count motifs in window
            in_window = class_df[
                (class_df['Start'] >= window_start) & 
                (class_df['Start'] < window_end)
            ]
            density = len(in_window)
            densities.append(density)
        
        density_matrix.append(densities)
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=density_matrix,
        x=[f"{w//1000}k" for w in windows],
        y=classes,
        colorscale='viridis',
        showscale=True,
        colorbar=dict(title="Motif Count")
    ))
    
    fig.update_layout(
        title=f"{title} - {sequence_name}",
        xaxis_title="Position (kb)",
        yaxis_title="Motif Class",
        height=400
    )
    
    return fig


def plot_overlap_analysis(df: pd.DataFrame, title: str = "Overlap Analysis") -> go.Figure:
    """
    Create comprehensive overlap analysis visualization.
    
    This function generates multiple visualizations to analyze overlapping motifs:
    1. Overlap frequency by class
    2. Resolution strategy effectiveness
    3. Score distribution of overlapping vs non-overlapping motifs
    
    Args:
        df (pd.DataFrame): Results DataFrame with overlap information
        title (str): Title for the plot
        
    Returns:
        go.Figure: Plotly figure with overlap analysis
    """
    check_dependencies()
    
    if df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No data available for overlap analysis",
                          xref="paper", yref="paper", x=0.5, y=0.5,
                          showarrow=False, font=dict(size=16))
        fig.update_layout(title=title)
        return fig
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=[
            "Overlap Frequency by Class",
            "Score Distribution: Overlapping vs Non-overlapping",
            "Overlap Classes Distribution", 
            "Motif Length vs Overlap Status"
        ],
        specs=[[{"type": "bar"}, {"type": "box"}],
               [{"type": "pie"}, {"type": "scatter"}]]
    )
    
    # Prepare overlap data
    df_copy = df.copy()
    df_copy['Has_Overlap'] = df_copy.get('Overlap_Classes', '').astype(str).str.len() > 0
    df_copy['Overlap_Count'] = df_copy.get('Overlap_Classes', '').astype(str).str.count(',') + 1
    df_copy.loc[df_copy['Overlap_Classes'].astype(str) == '', 'Overlap_Count'] = 0
    
    # 1. Overlap frequency by class
    overlap_by_class = df_copy.groupby('Class')['Has_Overlap'].agg(['sum', 'count']).reset_index()
    overlap_by_class['Overlap_Rate'] = overlap_by_class['sum'] / overlap_by_class['count'] * 100
    
    fig.add_trace(
        go.Bar(
            x=overlap_by_class['Class'],
            y=overlap_by_class['Overlap_Rate'],
            name="Overlap Rate (%)",
            marker=dict(color='lightcoral')
        ),
        row=1, col=1
    )
    
    # 2. Score distribution
    overlapping_scores = df_copy[df_copy['Has_Overlap']]['Normalized_Score'].dropna()
    non_overlapping_scores = df_copy[~df_copy['Has_Overlap']]['Normalized_Score'].dropna()
    
    if len(overlapping_scores) > 0:
        fig.add_trace(
            go.Box(y=overlapping_scores, name="Overlapping", marker=dict(color='red')),
            row=1, col=2
        )
    
    if len(non_overlapping_scores) > 0:
        fig.add_trace(
            go.Box(y=non_overlapping_scores, name="Non-overlapping", marker=dict(color='blue')),
            row=1, col=2
        )
    
    # 3. Overlap classes distribution
    overlap_classes_flat = []
    for overlap_str in df_copy['Overlap_Classes'].dropna():
        if overlap_str and str(overlap_str) != 'nan':
            classes = str(overlap_str).split(',')
            overlap_classes_flat.extend([cls.strip() for cls in classes if cls.strip()])
    
    if overlap_classes_flat:
        overlap_counts = pd.Series(overlap_classes_flat).value_counts()
        fig.add_trace(
            go.Pie(
                labels=overlap_counts.index,
                values=overlap_counts.values,
                name="Overlap Classes"
            ),
            row=2, col=1
        )
    
    # 4. Length vs overlap
    fig.add_trace(
        go.Scatter(
            x=df_copy['Length'],
            y=df_copy['Normalized_Score'],
            mode='markers',
            marker=dict(
                color=df_copy['Has_Overlap'].map({True: 'red', False: 'blue'}),
                size=8,
                opacity=0.6
            ),
            name="Motifs",
            text=df_copy['Class'],
            hovertemplate="<b>%{text}</b><br>Length: %{x}<br>Score: %{y}<br>Overlapping: %{marker.color}<extra></extra>"
        ),
        row=2, col=2
    )
    
    # Update layout
    fig.update_layout(
        title=title,
        height=800,
        showlegend=True
    )
    
    # Update axis labels
    fig.update_xaxes(title_text="Motif Class", row=1, col=1)
    fig.update_yaxes(title_text="Overlap Rate (%)", row=1, col=1)
    fig.update_yaxes(title_text="Normalized Score", row=1, col=2)
    fig.update_xaxes(title_text="Motif Length", row=2, col=2)
    fig.update_yaxes(title_text="Normalized Score", row=2, col=2)
    
    return fig


def plot_resolution_effectiveness(before_df: pd.DataFrame, after_df: pd.DataFrame, 
                                title: str = "Overlap Resolution Effectiveness") -> go.Figure:
    """
    Compare motif distributions before and after overlap resolution.
    
    Args:
        before_df (pd.DataFrame): Results before overlap resolution
        after_df (pd.DataFrame): Results after overlap resolution
        title (str): Title for the plot
        
    Returns:
        go.Figure: Plotly figure showing resolution effectiveness
    """
    check_dependencies()
    
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=["Before Resolution", "After Resolution"],
        specs=[[{"type": "bar"}, {"type": "bar"}]]
    )
    
    # Count by class before and after
    before_counts = before_df['Class'].value_counts() if not before_df.empty else pd.Series()
    after_counts = after_df['Class'].value_counts() if not after_df.empty else pd.Series()
    
    # Before resolution
    if not before_counts.empty:
        fig.add_trace(
            go.Bar(
                x=before_counts.index,
                y=before_counts.values,
                name="Before",
                marker=dict(color='lightblue')
            ),
            row=1, col=1
        )
    
    # After resolution
    if not after_counts.empty:
        fig.add_trace(
            go.Bar(
                x=after_counts.index,
                y=after_counts.values,
                name="After",
                marker=dict(color='darkblue')
            ),
            row=1, col=2
        )
    
    # Add reduction percentages
    reduction_stats = []
    all_classes = set(before_counts.index) | set(after_counts.index)
    
    for cls in all_classes:
        before_count = before_counts.get(cls, 0)
        after_count = after_counts.get(cls, 0)
        reduction = (before_count - after_count) / before_count * 100 if before_count > 0 else 0
        reduction_stats.append(f"{cls}: {reduction:.1f}% reduction")
    
    # Update layout
    fig.update_layout(
        title=f"{title}<br><sub>{'<br>'.join(reduction_stats[:5])}</sub>",
        height=500,
        showlegend=False
    )
    
    fig.update_xaxes(title_text="Motif Class", row=1, col=1)
    fig.update_xaxes(title_text="Motif Class", row=1, col=2)
    fig.update_yaxes(title_text="Count", row=1, col=1)
    fig.update_yaxes(title_text="Count", row=1, col=2)
    
    return fig
    """
    Create network plot showing overlapping motifs
    
    Args:
        df: Results DataFrame
        title: Plot title
        
    Returns:
        Plotly Figure object
    """
    check_dependencies()
    
    if not NETWORKX_AVAILABLE:
        fig = go.Figure()
        fig.add_annotation(text="NetworkX required for network plots",
                          xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
        return fig
    
    # Filter for motifs with overlaps
    overlap_df = df[df['Overlap_Classes'].str.len() > 0]
    
    if overlap_df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No overlapping motifs found",
                          xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
        return fig
    
    # Build network
    G = nx.Graph()
    
    # Add nodes (motif classes)
    classes = df['Class'].unique()
    for class_name in classes:
        G.add_node(class_name)
    
    # Add edges (overlaps)
    edge_weights = {}
    for _, row in overlap_df.iterrows():
        class1 = row['Class']
        overlaps = row['Overlap_Classes'].split(',') if row['Overlap_Classes'] else []
        
        for class2 in overlaps:
            if class2.strip() and class2.strip() in classes:
                edge = tuple(sorted([class1, class2.strip()]))
                edge_weights[edge] = edge_weights.get(edge, 0) + 1
    
    # Add weighted edges
    for edge, weight in edge_weights.items():
        G.add_edge(edge[0], edge[1], weight=weight)
    
    # Layout
    pos = nx.spring_layout(G)
    
    # Create edge traces
    edge_x = []
    edge_y = []
    edge_weights_list = []
    
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
        edge_weights_list.append(G[edge[0]][edge[1]]['weight'])
    
    # Create node traces
    node_x = []
    node_y = []
    node_text = []
    node_sizes = []
    
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_text.append(node)
        # Size based on degree
        node_sizes.append(10 + G.degree(node) * 5)
    
    # Create figure
    fig = go.Figure()
    
    # Add edges
    fig.add_trace(go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=2, color='gray'),
        hoverinfo='none',
        mode='lines'
    ))
    
    # Add nodes
    fig.add_trace(go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        hoverinfo='text',
        text=node_text,
        textposition="middle center",
        marker=dict(
            size=node_sizes,
            color='lightblue',
            line=dict(width=2, color='black')
        )
    ))
    
    fig.update_layout(
        title=title,
        showlegend=False,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        height=500
    )
    
    return fig


def plot_sequence_overview(df: pd.DataFrame, sequence_name: str, seq_len: int,
                          title: str = "Sequence Overview") -> go.Figure:
    """
    Create overview plot showing all motifs on a sequence
    
    Args:
        df: Results DataFrame
        sequence_name: Name of sequence to plot
        seq_len: Total sequence length
        title: Plot title
        
    Returns:
        Plotly Figure object
    """
    check_dependencies()
    
    # Filter for specific sequence
    seq_df = df[df['Sequence_Name'] == sequence_name]
    
    if seq_df.empty:
        fig = go.Figure()
        fig.add_annotation(text=f"No motifs found for sequence: {sequence_name}",
                          xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
        return fig
    
    # Create figure
    fig = go.Figure()
    
    # Color map for classes
    classes = seq_df['Class'].unique()
    colors = px.colors.qualitative.Set1[:len(classes)]
    color_map = dict(zip(classes, colors))
    
    # Add motifs as horizontal bars
    for i, (_, row) in enumerate(seq_df.iterrows()):
        fig.add_trace(go.Scatter(
            x=[row['Start'], row['End']],
            y=[row['Class'], row['Class']],
            mode='lines',
            line=dict(color=color_map[row['Class']], width=8),
            hovertemplate=f"<b>{row['Class']}</b><br>" +
                         f"Position: {row['Start']}-{row['End']}<br>" +
                         f"Length: {row['Length']} bp<br>" +
                         f"Score: {row['Normalized_Score']:.3f}<extra></extra>",
            showlegend=False
        ))
    
    # Add sequence length reference
    fig.add_vline(x=seq_len, line_dash="dash", line_color="gray",
                  annotation_text=f"Sequence end ({seq_len} bp)")
    
    fig.update_layout(
        title=f"{title} - {sequence_name}",
        xaxis_title="Position (bp)",
        yaxis_title="Motif Class",
        height=max(300, len(classes) * 40),
        hovermode='closest'
    )
    
    return fig


def tsne_feature_plot(df: pd.DataFrame, title: str = "Motif Feature Space (t-SNE)") -> go.Figure:
    """
    Create t-SNE plot of motif features
    
    Args:
        df: Results DataFrame with numeric features
        title: Plot title
        
    Returns:
        Plotly Figure object
    """
    check_dependencies()
    
    if not SKLEARN_AVAILABLE:
        fig = go.Figure()
        fig.add_annotation(text="Scikit-learn required for t-SNE",
                          xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
        return fig
    
    if df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No data to display",
                          xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
        return fig
    
    # Select numeric features
    feature_cols = ['Length', 'Normalized_Score', 'Actual_Score', 'GC_Content']
    feature_cols = [col for col in feature_cols if col in df.columns]
    
    if len(feature_cols) < 2:
        fig = go.Figure()
        fig.add_annotation(text="Insufficient numeric features for t-SNE",
                          xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
        return fig
    
    # Prepare data
    features = df[feature_cols].fillna(0)
    
    if len(features) < 3:
        fig = go.Figure()
        fig.add_annotation(text="Need at least 3 samples for t-SNE",
                          xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
        return fig
    
    # Standardize features
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(features)
    
    # Run t-SNE
    perplexity = min(30, len(features) - 1)
    tsne = TSNE(n_components=2, perplexity=perplexity, random_state=42)
    tsne_result = tsne.fit_transform(features_scaled)
    
    # Create plot
    fig = px.scatter(
        x=tsne_result[:, 0],
        y=tsne_result[:, 1],
        color=df['Class'],
        title=title,
        labels={'x': 't-SNE 1', 'y': 't-SNE 2'},
        hover_data={'Length': df['Length'], 'Score': df['Normalized_Score']}
    )
    
    fig.update_layout(height=500)
    return fig


def export_excel_with_sheets(df: pd.DataFrame, summary_stats: Optional[Dict] = None) -> bytes:
    """
    Export DataFrame to Excel with multiple sheets and return as bytes
    
    Args:
        df: Main results DataFrame
        summary_stats: Optional summary statistics
        
    Returns:
        Excel file as bytes for Streamlit download
    """
    from io import BytesIO
    
    output = BytesIO()
    
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        # Main results sheet
        df.to_excel(writer, sheet_name='All_Motifs', index=False)
        
        # Class-specific sheets
        for class_name in df['Class'].unique():
            class_df = df[df['Class'] == class_name]
            sheet_name = class_name.replace('_', ' ').title()[:31]
            class_df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        # Summary statistics sheet
        if summary_stats:
            stats_df = pd.DataFrame.from_dict(summary_stats, orient='index', columns=['Value'])
            stats_df.to_excel(writer, sheet_name='Summary')
    
    output.seek(0)
    return output.getvalue()


def create_summary_statistics(df: pd.DataFrame) -> Dict[str, Any]:
    """
    Generate summary statistics for results
    
    Args:
        df: Results DataFrame
        
    Returns:
        Dictionary of summary statistics
    """
    if df.empty:
        return {}
    
    stats = {
        'Total_Motifs': len(df),
        'Unique_Classes': df['Class'].nunique(),
        'Unique_Sequences': df['Sequence_Name'].nunique(),
        'Average_Length': df['Length'].mean(),
        'Median_Score': df['Normalized_Score'].median(),
        'Max_Score': df['Normalized_Score'].max(),
        'Min_Score': df['Normalized_Score'].min(),
    }
    
    # Class-specific counts
    class_counts = df['Class'].value_counts()
    for class_name, count in class_counts.items():
        stats[f'{class_name}_count'] = count
    
    return stats