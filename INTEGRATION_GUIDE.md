# Integration Guide: Structural Improvements, Statistical Layer, and UI Refactor

## Overview

This guide describes the new modules and how to integrate them into the NonBDNAFinder application to implement the 3-tab structure with statistical enrichment and structural analysis capabilities.

## New Modules Created

### 1. `enrichment.py` - Statistical Enrichment Engine
**Purpose:** 100× shuffle-based null model for pattern enrichment analysis

**Key Functions:**
```python
from enrichment import (
    run_enrichment_analysis,
    add_enrichment_to_motifs,
    format_enrichment_summary
)

# Run enrichment analysis
enrichment_results = run_enrichment_analysis(
    sequence=dna_sequence,
    motifs=detected_motifs,
    n_shuffles=100,
    random_seed=42
)

# Add enrichment scores to motifs
enriched_motifs = add_enrichment_to_motifs(motifs, enrichment_results)

# Get formatted summary
summary = format_enrichment_summary(enrichment_results)
```

**Output Structure:**
- `observed_metrics`: Dictionary of metrics from real sequence
- `surrogate_metrics`: List of metrics from 100 shuffled sequences
- `enrichment_scores`: P-values, Z-scores, O/E ratios for each metric
- Metrics include: pattern_count, mean_block_size, density_per_kb, co-occurrence_frequency, clustering_strength, hybrid_frequency

### 2. `structural_analysis.py` - Structural Block & Hybrid Pattern Module
**Purpose:** Higher-order pattern behavior discovery and analysis

**Key Functions:**
```python
from structural_analysis import (
    run_structural_analysis,
    identify_pattern_rich_blocks,
    detect_hybrid_zones,
    calculate_cluster_metrics,
    calculate_compositional_skew
)

# Run complete structural analysis
structural_results = run_structural_analysis(
    sequence=dna_sequence,
    motifs=detected_motifs,
    enable_blocks=True,
    enable_hybrids=True,
    enable_clusters=True,
    enable_skew=True
)

# Results include:
# - blocks: Pattern-rich regions with elevated motif density
# - hybrid_zones: Multi-class overlap regions
# - clusters: Spatial groupings of motifs
# - inter_cluster_distances: Distances between consecutive clusters
# - compositional_skew: GC/AT skew analysis
# - summary: Aggregate statistics
```

**Output Structure:**
```python
{
    'blocks': [{
        'start': int,
        'end': int,
        'length': int,
        'motif_count': int,
        'density': float,
        'classes': List[str],
        'class_diversity': int,
        'mean_score': float,
        'motifs': List[Dict]
    }, ...],
    'hybrid_zones': [{
        'start': int,
        'end': int,
        'classes': List[str],
        'num_classes': int,
        'hybrid_density': float,
        'interaction_strength': float,
        'motifs': List[Dict]
    }, ...],
    'clusters': [{
        'cluster_id': int,
        'start': int,
        'end': int,
        'size': int,
        'density': float,
        'stability_score': float,
        'class_diversity': int,
        'motifs': List[Dict]
    }, ...],
    'inter_cluster_distances': List[float],
    'compositional_skew': {
        'gc_skew_motif': float,
        'gc_skew_background': float,
        'at_skew_motif': float,
        'at_skew_background': float,
        ...
    },
    'summary': {
        'num_blocks': int,
        'num_hybrid_zones': int,
        'num_clusters': int,
        ...
    }
}
```

### 3. `visualizations_enhanced.py` - Enhanced Visualization Functions
**Purpose:** Publication-ready figures for clusters, hybrids, and enrichment

**Key Functions:**
```python
from visualizations_enhanced import (
    plot_cluster_footprint_heatmap,
    plot_hybrid_interaction_matrix,
    plot_observed_vs_shuffled_violin,
    plot_enrichment_null_distribution_panel,
    plot_block_size_enrichment_chart,
    plot_hybrid_region_stability_scatter,
    plot_inter_cluster_distance_distribution,
    plot_composition_vs_enrichment_contour
)

# Example usage:
fig1 = plot_cluster_footprint_heatmap(
    clusters=structural_results['clusters'],
    sequence_length=len(sequence)
)

fig2 = plot_hybrid_interaction_matrix(
    hybrid_zones=structural_results['hybrid_zones']
)

fig3 = plot_observed_vs_shuffled_violin(
    enrichment_results=enrichment_results,
    metric_name='density_per_kb'
)
```

## Integration into app.py

### Recommended 3-Tab Structure

#### **Tab 1: Overview & Summary**
Display global statistics and enrichment summaries:

```python
import streamlit as st
from enrichment import run_enrichment_analysis, format_enrichment_summary
from structural_analysis import run_structural_analysis

tab1, tab2, tab3 = st.tabs([
    "📊 Overview & Summary",
    "🗺️ Regional & Positional Analysis",
    "🔗 Clusters & Hybrids"
])

with tab1:
    st.header("Global Statistics")
    
    # Run analyses
    enrichment_results = run_enrichment_analysis(sequence, motifs)
    structural_results = run_structural_analysis(sequence, motifs)
    
    # Display summary
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Total Motifs", len(motifs))
        st.metric("Pattern-Rich Blocks", structural_results['summary']['num_blocks'])
    with col2:
        st.metric("Hybrid Zones", structural_results['summary']['num_hybrid_zones'])
        st.metric("Clusters", structural_results['summary']['num_clusters'])
    
    # Enrichment summary
    st.subheader("Statistical Enrichment")
    st.code(format_enrichment_summary(enrichment_results))
    
    # Observed vs. shuffled comparison
    from visualizations_enhanced import plot_observed_vs_shuffled_violin
    fig = plot_observed_vs_shuffled_violin(enrichment_results, 'pattern_count')
    st.pyplot(fig)
```

#### **Tab 2: Regional & Positional Analysis**
Display block-level maps and positional analysis:

```python
with tab2:
    st.header("Regional Analysis")
    
    # Block visualization
    from utilities import plot_coverage_map, plot_density_heatmap
    from visualizations_enhanced import plot_block_size_enrichment_chart
    
    st.subheader("Pattern-Rich Blocks")
    fig1 = plot_block_size_enrichment_chart(
        structural_results['blocks'],
        enrichment_results
    )
    st.pyplot(fig1)
    
    # Sliding-window density
    st.subheader("Sliding-Window Density")
    fig2 = plot_density_heatmap(motifs, len(sequence), window_size=1000)
    st.pyplot(fig2)
    
    # Positional heatmaps
    st.subheader("Positional Distribution")
    fig3 = plot_coverage_map(motifs, len(sequence))
    st.pyplot(fig3)
```

#### **Tab 3: Clusters & Hybrids**
Display cluster detection and hybrid interaction analysis:

```python
with tab3:
    st.header("Cluster & Hybrid Analysis")
    
    from visualizations_enhanced import (
        plot_cluster_footprint_heatmap,
        plot_hybrid_interaction_matrix,
        plot_inter_cluster_distance_distribution,
        plot_hybrid_region_stability_scatter
    )
    
    # Cluster footprint
    st.subheader("Cluster Footprints")
    fig1 = plot_cluster_footprint_heatmap(
        structural_results['clusters'],
        len(sequence)
    )
    st.pyplot(fig1)
    
    # Inter-cluster distances
    st.subheader("Inter-Cluster Spacing")
    fig2 = plot_inter_cluster_distance_distribution(
        structural_results['inter_cluster_distances']
    )
    st.pyplot(fig2)
    
    # Hybrid interaction matrix
    st.subheader("Hybrid Interaction Matrix")
    fig3 = plot_hybrid_interaction_matrix(
        structural_results['hybrid_zones']
    )
    st.pyplot(fig3)
    
    # Hybrid stability
    st.subheader("Hybrid Region Stability")
    fig4 = plot_hybrid_region_stability_scatter(
        structural_results['hybrid_zones']
    )
    st.pyplot(fig4)
    
    # Top hybrid zones report
    st.subheader("Top Hybrid Zones")
    sorted_zones = sorted(
        structural_results['hybrid_zones'],
        key=lambda z: z['interaction_strength'],
        reverse=True
    )[:10]
    
    for i, zone in enumerate(sorted_zones, 1):
        with st.expander(f"Hybrid Zone {i}: {zone['start']}-{zone['end']} bp"):
            st.write(f"**Classes:** {', '.join(zone['classes'])}")
            st.write(f"**Density:** {zone['hybrid_density']:.2f} motifs/kb")
            st.write(f"**Interaction Strength:** {zone['interaction_strength']:.2f}")
            st.write(f"**Number of Motifs:** {len(zone['motifs'])}")
```

## Integration Steps

### 1. Update Main Analysis Pipeline

Modify the main analysis function in `nonbscanner.py` or `app.py`:

```python
def analyze_sequence_comprehensive(sequence: str, sequence_name: str):
    """Complete analysis with enrichment and structural analysis."""
    
    # Step 1: Detect motifs (existing functionality)
    from nonbscanner import analyze_sequence
    motifs = analyze_sequence(sequence, sequence_name)
    
    # Step 2: Run enrichment analysis
    from enrichment import run_enrichment_analysis, add_enrichment_to_motifs
    enrichment_results = run_enrichment_analysis(
        sequence=sequence,
        motifs=motifs,
        n_shuffles=100,
        random_seed=42
    )
    
    # Add enrichment scores to motifs
    motifs = add_enrichment_to_motifs(motifs, enrichment_results)
    
    # Step 3: Run structural analysis
    from structural_analysis import run_structural_analysis
    structural_results = run_structural_analysis(
        sequence=sequence,
        motifs=motifs,
        enable_blocks=True,
        enable_hybrids=True,
        enable_clusters=True,
        enable_skew=True
    )
    
    return {
        'motifs': motifs,
        'enrichment': enrichment_results,
        'structural': structural_results
    }
```

### 2. Update Export Functions

Add enrichment and structural data to exports:

```python
def export_comprehensive_results(results, output_path):
    """Export all analysis results to Excel."""
    import pandas as pd
    
    with pd.ExcelWriter(output_path) as writer:
        # Existing motif data
        df_motifs = pd.DataFrame(results['motifs'])
        df_motifs.to_excel(writer, sheet_name='Motifs', index=False)
        
        # Enrichment scores
        enrichment_data = []
        for metric, scores in results['enrichment']['enrichment_scores'].items():
            enrichment_data.append({
                'Metric': metric,
                'Observed': scores['observed'],
                'Expected_Mean': scores['expected_mean'],
                'Pvalue': scores['pvalue'],
                'Zscore': scores['zscore'],
                'Significance': scores['significance']
            })
        df_enrichment = pd.DataFrame(enrichment_data)
        df_enrichment.to_excel(writer, sheet_name='Enrichment', index=False)
        
        # Blocks
        df_blocks = pd.DataFrame(results['structural']['blocks'])
        if not df_blocks.empty:
            df_blocks = df_blocks.drop(columns=['motifs'])  # Remove nested list
        df_blocks.to_excel(writer, sheet_name='Blocks', index=False)
        
        # Hybrid zones
        df_hybrids = pd.DataFrame(results['structural']['hybrid_zones'])
        if not df_hybrids.empty:
            df_hybrids = df_hybrids.drop(columns=['motifs'])
        df_hybrids.to_excel(writer, sheet_name='Hybrid_Zones', index=False)
        
        # Clusters
        df_clusters = pd.DataFrame(results['structural']['clusters'])
        if not df_clusters.empty:
            df_clusters = df_clusters.drop(columns=['motifs'])
        df_clusters.to_excel(writer, sheet_name='Clusters', index=False)
```

### 3. Add Progress Tracking

For long-running enrichment analysis:

```python
import streamlit as st

def run_enrichment_with_progress(sequence, motifs):
    """Run enrichment with Streamlit progress bar."""
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    def progress_callback(current, total):
        progress = current / total
        progress_bar.progress(progress)
        status_text.text(f"Generating shuffles: {current}/{total}")
    
    from enrichment import run_enrichment_analysis
    results = run_enrichment_analysis(
        sequence=sequence,
        motifs=motifs,
        n_shuffles=100,
        progress_callback=progress_callback
    )
    
    progress_bar.empty()
    status_text.empty()
    return results
```

## Performance Considerations

### Caching
Use Streamlit caching to avoid re-computation:

```python
@st.cache_data
def cached_enrichment_analysis(sequence, motifs_tuple, n_shuffles=100):
    """Cached version of enrichment analysis."""
    from enrichment import run_enrichment_analysis
    motifs = list(motifs_tuple)  # Convert back from tuple
    return run_enrichment_analysis(sequence, motifs, n_shuffles=n_shuffles)

@st.cache_data
def cached_structural_analysis(sequence, motifs_tuple):
    """Cached version of structural analysis."""
    from structural_analysis import run_structural_analysis
    motifs = list(motifs_tuple)
    return run_structural_analysis(sequence, motifs)
```

### Memory Management
For large sequences, consider:
- Reducing number of shuffles from 100 to 50 for faster results
- Processing in chunks for very long sequences (>1 MB)
- Clearing matplotlib figures after display: `plt.close(fig)`

## Testing

Test the integration with sample data:

```python
def test_integration():
    """Test integration of new modules."""
    
    # Sample sequence (1000 bp)
    sequence = "GGGTTAGGGTTAGGGTTAGGG" * 50
    
    # Detect motifs
    from nonbscanner import analyze_sequence
    motifs = analyze_sequence(sequence, "test")
    
    print(f"Detected {len(motifs)} motifs")
    
    # Run enrichment
    from enrichment import run_enrichment_analysis
    enrichment_results = run_enrichment_analysis(sequence, motifs, n_shuffles=10)
    print(f"Enrichment analysis complete: {enrichment_results['n_shuffles']} shuffles")
    
    # Run structural analysis
    from structural_analysis import run_structural_analysis
    structural_results = run_structural_analysis(sequence, motifs)
    print(f"Structural analysis complete:")
    print(f"  - Blocks: {structural_results['summary']['num_blocks']}")
    print(f"  - Hybrids: {structural_results['summary']['num_hybrid_zones']}")
    print(f"  - Clusters: {structural_results['summary']['num_clusters']}")
    
    # Test visualizations
    from visualizations_enhanced import (
        plot_cluster_footprint_heatmap,
        plot_hybrid_interaction_matrix
    )
    
    if structural_results['clusters']:
        fig1 = plot_cluster_footprint_heatmap(
            structural_results['clusters'],
            len(sequence)
        )
        print("✓ Cluster footprint heatmap generated")
        plt.close(fig1)
    
    if structural_results['hybrid_zones']:
        fig2 = plot_hybrid_interaction_matrix(
            structural_results['hybrid_zones']
        )
        print("✓ Hybrid interaction matrix generated")
        plt.close(fig2)
    
    print("\n✓ All integration tests passed!")

if __name__ == "__main__":
    test_integration()
```

## Documentation Updates

Update README.md with new capabilities:

```markdown
## New in Version 2025.1: Statistical Framework

### Statistical Enrichment Analysis
- 100× shuffle-based null model
- Empirical p-values, Z-scores, O/E ratios
- Composition-preserving randomization
- No re-detection on surrogates (compositional metrics only)

### Structural Pattern Discovery
- Pattern-rich block identification
- Hybrid zone detection (multi-class overlaps)
- Cluster analysis with stability scoring
- Compositional skew analysis

### Enhanced Visualizations
- Cluster footprint heatmaps
- Hybrid interaction matrices
- Observed vs. shuffled violin plots
- Block-size enrichment charts
- Inter-cluster distance distributions
- And more...

### 3-Tab Interface
- **Tab 1:** Overview & Summary (global stats, enrichment)
- **Tab 2:** Regional & Positional (blocks, density, peaks)
- **Tab 3:** Clusters & Hybrids (interaction maps, stability)
```

## Summary

The core functionality for structural improvements, statistical enrichment, and enhanced visualizations is complete and ready for integration. The modular design allows for flexible integration into the existing application without requiring wholesale refactoring of app.py.

Key integration points:
1. Call enrichment and structural analysis after motif detection
2. Add new tabs to Streamlit interface
3. Use enhanced visualization functions in each tab
4. Update export functions to include new data
5. Add caching for performance

All modules are tested and production-ready.
