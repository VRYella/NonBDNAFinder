"""
Advanced Non-B DNA Motif Visualization Suite
====================================

Features: counts, distributions, locations, overlaps, interactions, 
density, networks, dimensionality reduction, interactive plots.

This module provides comprehensive visualization capabilities for 
analyzing Non-B DNA motifs detected by NBDFinder.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import networkx as nx
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import random
from matplotlib_venn import venn2, venn3
from typing import List, Dict, Any, Optional
import warnings
warnings.filterwarnings('ignore')

# Set style for matplotlib
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

def create_motif_count_plots(df: pd.DataFrame, save_plots: bool = False, output_dir: str = "."):
    """
    Create motif count and hierarchy visualizations
    
    Args:
        df: DataFrame with motif data
        save_plots: Whether to save plots to files
        output_dir: Directory to save plots
    """
    if df.empty:
        print("No data available for visualization")
        return
    
    # 1. Motif Count per Class
    plt.figure(figsize=(12, 6))
    class_counts = df['Class'].value_counts()
    sns.countplot(data=df, x='Class', order=class_counts.index)
    plt.title("Motif Count per Class", fontsize=16, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.ylabel("Count")
    plt.tight_layout()
    if save_plots:
        plt.savefig(f"{output_dir}/motif_count_per_class.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    # 2. Stacked Bar: Subclass Distribution
    if 'Subclass' in df.columns:
        pivot = df.groupby(['Class', 'Subclass']).size().reset_index(name='Count')
        pivot_table = pivot.pivot(index='Class', columns='Subclass', values='Count').fillna(0)
        
        plt.figure(figsize=(14, 8))
        pivot_table.plot(kind='bar', stacked=True, colormap='tab20', ax=plt.gca())
        plt.title("Stacked Bar: Motif Subclass Distribution", fontsize=16, fontweight='bold')
        plt.ylabel("Count")
        plt.xlabel("Class")
        plt.xticks(rotation=45, ha='right')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        if save_plots:
            plt.savefig(f"{output_dir}/subclass_distribution.png", dpi=300, bbox_inches='tight')
        plt.show()
    
    # 3. Pie/Donut Chart
    plt.figure(figsize=(10, 8))
    colors = plt.cm.Set3(np.linspace(0, 1, len(class_counts)))
    wedges, texts, autotexts = plt.pie(class_counts.values, labels=class_counts.index, 
                                      autopct='%1.1f%%', colors=colors,
                                      wedgeprops=dict(width=0.5))
    plt.title("Donut Chart: Motif Class Proportion", fontsize=16, fontweight='bold')
    if save_plots:
        plt.savefig(f"{output_dir}/class_proportion_donut.png", dpi=300, bbox_inches='tight')
    plt.show()

def create_interactive_hierarchy_plots(df: pd.DataFrame, save_plots: bool = False, output_dir: str = "."):
    """
    Create interactive sunburst and treemap plots using Plotly
    
    Args:
        df: DataFrame with motif data
        save_plots: Whether to save plots to files
        output_dir: Directory to save plots
    """
    if df.empty or 'Subclass' not in df.columns:
        print("Insufficient data for interactive hierarchy plots")
        return
    
    # Sunburst Chart
    fig_sunburst = px.sunburst(df, path=['Class', 'Subclass'], 
                              title="Sunburst: Class/Subclass Hierarchy",
                              color_discrete_sequence=px.colors.qualitative.Set3)
    fig_sunburst.update_layout(font_size=12)
    if save_plots:
        fig_sunburst.write_html(f"{output_dir}/sunburst_hierarchy.html")
    fig_sunburst.show()
    
    # Treemap
    fig_treemap = px.treemap(df, path=['Class', 'Subclass'], 
                            title="Treemap: Class/Subclass Distribution",
                            color_discrete_sequence=px.colors.qualitative.Pastel)
    fig_treemap.update_layout(font_size=12)
    if save_plots:
        fig_treemap.write_html(f"{output_dir}/treemap_hierarchy.html")
    fig_treemap.show()

def create_score_distribution_plots(df: pd.DataFrame, save_plots: bool = False, output_dir: str = "."):
    """
    Create score distribution visualizations
    
    Args:
        df: DataFrame with motif data
        save_plots: Whether to save plots to files
        output_dir: Directory to save plots
    """
    if df.empty:
        print("No data available for score distribution plots")
        return
    
    score_col = 'Normalized_Score' if 'Normalized_Score' in df.columns else 'Actual_Score'
    if score_col not in df.columns:
        print("No score column found for distribution plots")
        return
    
    # 1. Box Plot
    plt.figure(figsize=(14, 6))
    sns.boxplot(data=df, x='Class', y=score_col)
    plt.title(f"{score_col.replace('_', ' ')} Distribution by Motif Class", fontsize=16, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    if save_plots:
        plt.savefig(f"{output_dir}/score_boxplot.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    # 2. Violin Plot
    plt.figure(figsize=(14, 6))
    sns.violinplot(data=df, x='Class', y=score_col)
    plt.title(f"Violin Plot: {score_col.replace('_', ' ')} Distribution", fontsize=16, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    if save_plots:
        plt.savefig(f"{output_dir}/score_violin.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    # 3. Histogram with KDE
    plt.figure(figsize=(10, 6))
    sns.histplot(df[score_col], bins=30, kde=True, stat='density')
    plt.title(f"Histogram: {score_col.replace('_', ' ')} Distribution", fontsize=16, fontweight='bold')
    plt.xlabel(score_col.replace('_', ' '))
    plt.ylabel("Density")
    plt.tight_layout()
    if save_plots:
        plt.savefig(f"{output_dir}/score_histogram.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    # 4. Cumulative Distribution Function (CDF)
    scores_sorted = np.sort(df[score_col].dropna())
    cdf = np.arange(1, len(scores_sorted) + 1) / len(scores_sorted)
    
    plt.figure(figsize=(10, 6))
    plt.plot(scores_sorted, cdf, linewidth=2)
    plt.xlabel(score_col.replace('_', ' '))
    plt.ylabel('Cumulative Probability')
    plt.title(f"CDF of {score_col.replace('_', ' ')}", fontsize=16, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    if save_plots:
        plt.savefig(f"{output_dir}/score_cdf.png", dpi=300, bbox_inches='tight')
    plt.show()

def create_mapping_density_plots(df: pd.DataFrame, save_plots: bool = False, output_dir: str = "."):
    """
    Create motif mapping and density visualizations
    
    Args:
        df: DataFrame with motif data
        save_plots: Whether to save plots to files
        output_dir: Directory to save plots
    """
    if df.empty or 'Start' not in df.columns or 'End' not in df.columns:
        print("Insufficient data for mapping/density plots")
        return
    
    # 1. Lollipop/Track Plot
    motif_classes = df['Class'].unique()
    colors = plt.cm.Set3(np.linspace(0, 1, len(motif_classes)))
    
    plt.figure(figsize=(15, 8))
    max_pos = df['End'].max() if not df['End'].empty else 1000
    
    for i, motif_class in enumerate(motif_classes):
        hits = df[df['Class'] == motif_class]
        if not hits.empty:
            plt.hlines(i, 0, max_pos, color='gray', alpha=0.15)
            plt.scatter(hits['Start'], [i] * len(hits), 
                       label=motif_class, s=60, alpha=0.7, color=colors[i])
    
    plt.yticks(range(len(motif_classes)), motif_classes)
    plt.xlabel("Sequence Position")
    plt.ylabel("Motif Class")
    plt.title("Lollipop/Track Plot: Motif Positions", fontsize=16, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    if save_plots:
        plt.savefig(f"{output_dir}/motif_track_plot.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    # 2. Density Heatmap
    max_pos = min(df['End'].max(), 10000) if not df['End'].empty else 1000  # Limit for visualization
    density = np.zeros(int(max_pos) + 1)
    
    for _, row in df.iterrows():
        start = int(row['Start'])
        end = int(row['End'])
        if end <= max_pos:
            density[start:end + 1] += 1
    
    plt.figure(figsize=(15, 3))
    plt.imshow(density[np.newaxis, :], aspect='auto', cmap='hot', extent=[0, max_pos, 0, 1])
    plt.xlabel("Position")
    plt.ylabel("")
    plt.title("Motif Density Heatmap", fontsize=16, fontweight='bold')
    plt.colorbar(label='Motif Density')
    plt.yticks([])
    plt.tight_layout()
    if save_plots:
        plt.savefig(f"{output_dir}/density_heatmap.png", dpi=300, bbox_inches='tight')
    plt.show()

def create_overlap_interaction_plots(df: pd.DataFrame, save_plots: bool = False, output_dir: str = "."):
    """
    Create overlap and interaction visualizations
    
    Args:
        df: DataFrame with motif data
        save_plots: Whether to save plots to files
        output_dir: Directory to save plots
    """
    if df.empty:
        print("No data available for overlap/interaction plots")
        return
    
    # 1. Venn Diagram (if we have at least 2 classes)
    classes = df['Class'].unique()
    if len(classes) >= 2:
        if 'Sequence_Name' in df.columns:
            class1_seqs = set(df[df['Class'] == classes[0]]['Sequence_Name'])
            class2_seqs = set(df[df['Class'] == classes[1]]['Sequence_Name'])
            
            plt.figure(figsize=(8, 6))
            venn2([class1_seqs, class2_seqs], set_labels=(classes[0], classes[1]))
            plt.title(f"Venn: Sequence Overlap between {classes[0]} and {classes[1]}", 
                     fontsize=14, fontweight='bold')
            if save_plots:
                plt.savefig(f"{output_dir}/venn_overlap.png", dpi=300, bbox_inches='tight')
            plt.show()
    
    # 2. Network Graph of motif co-occurrences
    try:
        # Create co-occurrence matrix based on overlapping regions
        co_occurrence = np.zeros((len(classes), len(classes)))
        class_to_idx = {cls: i for i, cls in enumerate(classes)}
        
        # Simple co-occurrence based on position proximity
        for i, row1 in df.iterrows():
            for j, row2 in df.iterrows():
                if i != j and abs(row1['Start'] - row2['Start']) < 100:  # Within 100 bp
                    cls1_idx = class_to_idx[row1['Class']]
                    cls2_idx = class_to_idx[row2['Class']]
                    co_occurrence[cls1_idx, cls2_idx] += 1
        
        # Create network graph
        plt.figure(figsize=(10, 8))
        G = nx.Graph()
        
        # Add nodes
        for cls in classes:
            G.add_node(cls)
        
        # Add edges with weights
        for i, cls1 in enumerate(classes):
            for j, cls2 in enumerate(classes):
                if i < j and co_occurrence[i, j] > 0:
                    G.add_edge(cls1, cls2, weight=co_occurrence[i, j])
        
        # Draw network
        pos = nx.spring_layout(G, seed=42)
        nx.draw_networkx_nodes(G, pos, node_color='skyblue', node_size=2000, alpha=0.7)
        nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold')
        
        # Draw edges with varying thickness based on weight
        edges = G.edges()
        weights = [G[u][v]['weight'] for u, v in edges]
        nx.draw_networkx_edges(G, pos, width=[w/max(weights)*5 for w in weights], alpha=0.6)
        
        plt.title("Motif Interaction Network", fontsize=16, fontweight='bold')
        plt.axis('off')
        if save_plots:
            plt.savefig(f"{output_dir}/interaction_network.png", dpi=300, bbox_inches='tight')
        plt.show()
        
    except Exception as e:
        print(f"Could not create network plot: {e}")

def create_sequence_feature_plots(df: pd.DataFrame, save_plots: bool = False, output_dir: str = "."):
    """
    Create sequence feature analysis plots
    
    Args:
        df: DataFrame with motif data
        save_plots: Whether to save plots to files
        output_dir: Directory to save plots
    """
    if df.empty:
        print("No data available for sequence feature plots")
        return
    
    # GC Content vs Position
    if 'GC_Content' in df.columns and 'Start' in df.columns:
        plt.figure(figsize=(12, 6))
        plt.scatter(df['Start'], df['GC_Content'], alpha=0.6, s=50)
        plt.xlabel("Position")
        plt.ylabel("GC Content (%)")
        plt.title("GC Content by Motif Position", fontsize=16, fontweight='bold')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        if save_plots:
            plt.savefig(f"{output_dir}/gc_content_position.png", dpi=300, bbox_inches='tight')
        plt.show()
    
    # Length distribution by class
    if 'Length' in df.columns:
        plt.figure(figsize=(12, 6))
        sns.boxplot(data=df, x='Class', y='Length')
        plt.title("Motif Length Distribution by Class", fontsize=16, fontweight='bold')
        plt.xticks(rotation=45, ha='right')
        plt.ylabel("Length (bp)")
        plt.tight_layout()
        if save_plots:
            plt.savefig(f"{output_dir}/length_distribution.png", dpi=300, bbox_inches='tight')
        plt.show()

def create_dimensionality_reduction_plots(df: pd.DataFrame, save_plots: bool = False, output_dir: str = "."):
    """
    Create dimensionality reduction visualizations
    
    Args:
        df: DataFrame with motif data
        save_plots: Whether to save plots to files
        output_dir: Directory to save plots
    """
    if df.empty:
        print("No data available for dimensionality reduction plots")
        return
    
    # Select numeric features for analysis
    numeric_features = []
    for col in ['Normalized_Score', 'Actual_Score', 'Length', 'GC_Content']:
        if col in df.columns:
            numeric_features.append(col)
    
    if len(numeric_features) < 2:
        print("Insufficient numeric features for dimensionality reduction")
        return
    
    # Prepare data
    features_df = df[numeric_features].dropna()
    if features_df.empty or len(features_df) < 10:
        print("Insufficient data points for dimensionality reduction")
        return
    
    # Standardize features
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(features_df)
    
    # t-SNE
    try:
        tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, len(features_df)-1))
        tsne_result = tsne.fit_transform(features_scaled)
        
        # Create DataFrame with t-SNE results
        tsne_df = df.loc[features_df.index].copy()
        tsne_df['TSNE1'] = tsne_result[:, 0]
        tsne_df['TSNE2'] = tsne_result[:, 1]
        
        plt.figure(figsize=(10, 8))
        classes = tsne_df['Class'].unique()
        colors = plt.cm.tab10(np.linspace(0, 1, len(classes)))
        
        for i, cls in enumerate(classes):
            class_data = tsne_df[tsne_df['Class'] == cls]
            plt.scatter(class_data['TSNE1'], class_data['TSNE2'], 
                       label=cls, alpha=0.7, s=50, color=colors[i])
        
        plt.xlabel("t-SNE 1")
        plt.ylabel("t-SNE 2")
        plt.title("t-SNE: Motif Feature Clustering", fontsize=16, fontweight='bold')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        if save_plots:
            plt.savefig(f"{output_dir}/tsne_clustering.png", dpi=300, bbox_inches='tight')
        plt.show()
        
    except Exception as e:
        print(f"Could not create t-SNE plot: {e}")

def create_manhattan_plot(df: pd.DataFrame, save_plots: bool = False, output_dir: str = "."):
    """
    Create Manhattan-like plot for motif scores
    
    Args:
        df: DataFrame with motif data
        save_plots: Whether to save plots to files
        output_dir: Directory to save plots
    """
    if df.empty:
        print("No data available for Manhattan plot")
        return
    
    score_col = 'Normalized_Score' if 'Normalized_Score' in df.columns else 'Actual_Score'
    if score_col not in df.columns or 'Start' not in df.columns:
        print("Insufficient data for Manhattan plot")
        return
    
    plt.figure(figsize=(15, 6))
    classes = df['Class'].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(classes)))
    
    for i, cls in enumerate(classes):
        hits = df[df['Class'] == cls]
        if not hits.empty:
            plt.scatter(hits['Start'], hits[score_col], 
                       label=cls, alpha=0.6, s=40, color=colors[i])
    
    plt.xlabel("Position")
    plt.ylabel(score_col.replace('_', ' '))
    plt.title("Manhattan Plot: Motif Scores vs Position", fontsize=16, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=1)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    if save_plots:
        plt.savefig(f"{output_dir}/manhattan_plot.png", dpi=300, bbox_inches='tight')
    plt.show()

def create_comprehensive_visualization_suite(df: pd.DataFrame, save_plots: bool = False, output_dir: str = "."):
    """
    Create the complete visualization suite for Non-B DNA motifs
    
    Args:
        df: DataFrame with motif data containing columns like:
            'Class', 'Subclass', 'Start', 'End', 'Length', 'Normalized_Score', 
            'Actual_Score', 'GC_Content', 'Sequence_Name'
        save_plots: Whether to save plots to files
        output_dir: Directory to save plots
    """
    print("🧬 Creating Advanced Non-B DNA Motif Visualization Suite")
    print("=" * 60)
    
    if df.empty:
        print("❌ No data provided for visualization")
        return
    
    print(f"📊 Analyzing {len(df)} motifs across {df['Class'].nunique()} classes")
    
    # Create all visualization categories
    print("\n1️⃣ Creating motif counts & hierarchy plots...")
    create_motif_count_plots(df, save_plots, output_dir)
    
    print("\n2️⃣ Creating interactive hierarchy plots...")
    create_interactive_hierarchy_plots(df, save_plots, output_dir)
    
    print("\n3️⃣ Creating score distribution plots...")
    create_score_distribution_plots(df, save_plots, output_dir)
    
    print("\n4️⃣ Creating mapping & density plots...")
    create_mapping_density_plots(df, save_plots, output_dir)
    
    print("\n5️⃣ Creating overlap & interaction plots...")
    create_overlap_interaction_plots(df, save_plots, output_dir)
    
    print("\n6️⃣ Creating sequence feature plots...")
    create_sequence_feature_plots(df, save_plots, output_dir)
    
    print("\n7️⃣ Creating dimensionality reduction plots...")
    create_dimensionality_reduction_plots(df, save_plots, output_dir)
    
    print("\n8️⃣ Creating Manhattan plot...")
    create_manhattan_plot(df, save_plots, output_dir)
    
    print("\n✅ Visualization suite completed!")
    
    # Summary statistics
    print(f"\n📈 Summary Statistics:")
    print(f"   • Total motifs: {len(df)}")
    print(f"   • Unique classes: {df['Class'].nunique()}")
    if 'Subclass' in df.columns:
        print(f"   • Unique subclasses: {df['Subclass'].nunique()}")
    if 'Sequence_Name' in df.columns:
        print(f"   • Sequences analyzed: {df['Sequence_Name'].nunique()}")
    
    # Class distribution
    print(f"\n🧪 Class Distribution:")
    class_counts = df['Class'].value_counts()
    for cls, count in class_counts.head().items():
        print(f"   • {cls}: {count} motifs")

# Hybrid and Cluster motif definitions
def identify_hybrid_motifs(df: pd.DataFrame, proximity_threshold: int = 50) -> pd.DataFrame:
    """
    Identify hybrid motifs: 2 or more overlapping motifs
    
    Args:
        df: DataFrame with motif data
        proximity_threshold: Maximum distance between motifs to consider them overlapping
    
    Returns:
        DataFrame with hybrid motifs
    """
    hybrids = []
    
    for i, motif1 in df.iterrows():
        overlapping_classes = []
        for j, motif2 in df.iterrows():
            if i != j:
                # Check if motifs overlap or are in proximity
                if (abs(motif1['Start'] - motif2['Start']) <= proximity_threshold or
                    (motif1['Start'] <= motif2['End'] and motif2['Start'] <= motif1['End'])):
                    overlapping_classes.append(motif2['Class'])
        
        if len(set(overlapping_classes)) >= 1:  # At least 1 other class
            hybrid = motif1.copy()
            hybrid['Class'] = 'Hybrid'
            hybrid['Subclass'] = f"Hybrid_{'+'.join(sorted(set([motif1['Class']] + overlapping_classes)))}"
            hybrid['Overlap_Classes'] = ','.join(sorted(set(overlapping_classes)))
            hybrids.append(hybrid)
    
    return pd.DataFrame(hybrids)

def identify_cluster_regions(df: pd.DataFrame, window_size: int = 100, min_classes: int = 3, min_occurrences: int = 3) -> pd.DataFrame:
    """
    Identify cluster regions: regions of 100 nucleotides or longer with at least 3 different classes of motifs
    
    Args:
        df: DataFrame with motif data
        window_size: Size of the sliding window
        min_classes: Minimum number of different motif classes required
        min_occurrences: Minimum occurrences per class in the window
    
    Returns:
        DataFrame with cluster regions
    """
    if df.empty or 'Start' not in df.columns:
        return pd.DataFrame()
    
    clusters = []
    max_pos = df['End'].max() if 'End' in df.columns else df['Start'].max()
    
    # Slide window across the sequence
    for start in range(0, int(max_pos), window_size // 2):  # 50% overlap
        end = start + window_size
        
        # Get motifs in this window
        window_motifs = df[(df['Start'] >= start) & (df['Start'] <= end)]
        
        if len(window_motifs) == 0:
            continue
        
        # Count classes and their occurrences
        class_counts = window_motifs['Class'].value_counts()
        qualifying_classes = class_counts[class_counts >= min_occurrences]
        
        if len(qualifying_classes) >= min_classes:
            cluster = {
                'Class': 'Non-B DNA cluster regions',
                'Subclass': f"Cluster_{len(qualifying_classes)}classes",
                'Start': start,
                'End': end,
                'Length': window_size,
                'Motif_ID': f"cluster_{start}_{end}",
                'Overlap_Classes': ','.join(sorted(qualifying_classes.index)),
                'Normalized_Score': np.mean(window_motifs['Normalized_Score']) if 'Normalized_Score' in window_motifs.columns else 0,
                'Actual_Score': np.mean(window_motifs['Actual_Score']) if 'Actual_Score' in window_motifs.columns else 0,
                'GC_Content': np.mean(window_motifs['GC_Content']) if 'GC_Content' in window_motifs.columns else 0,
                'Sequence_Name': window_motifs['Sequence_Name'].iloc[0] if 'Sequence_Name' in window_motifs.columns else 'unknown',
                'Chromosome/Contig': window_motifs['Chromosome/Contig'].iloc[0] if 'Chromosome/Contig' in window_motifs.columns else 'unknown'
            }
            clusters.append(cluster)
    
    return pd.DataFrame(clusters)

# Export all functions
__all__ = [
    'create_motif_count_plots',
    'create_interactive_hierarchy_plots', 
    'create_score_distribution_plots',
    'create_mapping_density_plots',
    'create_overlap_interaction_plots',
    'create_sequence_feature_plots',
    'create_dimensionality_reduction_plots',
    'create_manhattan_plot',
    'create_comprehensive_visualization_suite',
    'identify_hybrid_motifs',
    'identify_cluster_regions'
]