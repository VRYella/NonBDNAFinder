# 🧬 Nature-Quality Genome-Wide Visualizations for Non-B DNA Finder

## Overview

This document describes the new publication-quality visualizations implemented for the Non-B DNA Finder tool. All visualizations follow **Nature Methods** styling guidelines and are designed for scientific publication.

## New Visualization Functions

### 1. Manhattan Plot (`plot_manhattan_motif_density`)

**Purpose:** Highlights motif density hotspots across genomic coordinates.

**Best for:** Large genomes (human, mouse), identifying cluster regions and hybrid zones.

**Features:**
- X-axis: Genomic position (kb)
- Y-axis: Motif density (motifs/kb) or average score
- Color-coded by motif class
- Horizontal grid for readability
- 300 DPI publication quality

**Use case:** 
```python
fig = plot_manhattan_motif_density(
    motifs, sequence_length,
    score_type='density',  # or 'score'
    title="Manhattan Plot - Motif Hotspots"
)
```

**Output:** `/tmp/test_manhattan.png` (210KB, 300 DPI)

---

### 2. Cumulative Motif Distribution (`plot_cumulative_motif_distribution`)

**Purpose:** Shows running sum of motifs over genome, useful for comparing accumulation patterns.

**Best for:** Comparing motif distributions, identifying regions of rapid accumulation.

**Features:**
- Running sum line plot for each class
- Class-specific colors
- Can show overall or by-class distribution
- Clean Nature-style formatting

**Use case:**
```python
fig = plot_cumulative_motif_distribution(
    motifs, sequence_length,
    by_class=True,
    title="Cumulative Distribution"
)
```

**Output:** `/tmp/test_cumulative.png` (209KB, 300 DPI)

---

### 3. Motif Co-occurrence Matrix (`plot_motif_cooccurrence_matrix`)

**Purpose:** Heatmap showing which motif classes tend to co-occur (overlap or appear nearby).

**Best for:** Publication figures showing motif relationships and interactions.

**Features:**
- Symmetric matrix with class pairs
- Colorblind-friendly colormap (YlOrRd)
- Cell values show co-occurrence counts
- Configurable overlap threshold

**Use case:**
```python
fig = plot_motif_cooccurrence_matrix(
    motifs,
    overlap_threshold=1,  # bp
    title="Motif Co-occurrence Matrix"
)
```

**Output:** `/tmp/test_cooccurrence.png` (179KB, 300 DPI)

---

### 4. GC Content Correlation (`plot_gc_content_correlation`)

**Purpose:** Scatter plot showing relationship between GC content and motif density.

**Best for:** Demonstrating GC-driven motif enrichment patterns.

**Features:**
- Scatter plot with regression line
- Correlation coefficient (R) displayed
- Window-based analysis
- Statistical significance shown

**Use case:**
```python
fig = plot_gc_content_correlation(
    motifs, sequence,
    window_size=1000,
    title="GC Content vs Motif Density"
)
```

**Output:** `/tmp/test_gc_correlation.png` (126KB, 300 DPI)

---

### 5. Linear Motif Track (`plot_linear_motif_track`)

**Purpose:** Horizontal genome browser view with colored blocks for motifs.

**Best for:** Small regions (<10kb), detailed motif positioning.

**Features:**
- Colored rectangles for each motif
- Class-specific tracks
- Score labels on larger blocks
- Clean ruler at top
- Region zooming support

**Use case:**
```python
fig = plot_linear_motif_track(
    motifs, sequence_length,
    region_start=0,
    region_end=10000,
    title="Linear Motif Track (0-10kb)"
)
```

**Output:** `/tmp/test_linear_track.png` (88KB, 300 DPI)

---

### 6. Cluster Size Distribution (`plot_cluster_size_distribution`)

**Purpose:** Histograms showing distribution of cluster sizes and class diversity.

**Best for:** Analyzing cluster composition and characteristics.

**Features:**
- Two subplots: size distribution and diversity distribution
- Mean/median lines
- Clean histogram styling
- Statistical annotations

**Use case:**
```python
fig = plot_cluster_size_distribution(
    motifs,
    title="Cluster Statistics"
)
```

**Output:** `/tmp/test_cluster_dist.png` (102KB, 300 DPI)

---

### 7. Motif Length KDE (`plot_motif_length_kde`)

**Purpose:** Kernel density estimation showing smooth probability curves for motif lengths.

**Best for:** Comparing length patterns across classes, identifying modal lengths.

**Features:**
- Smooth probability density curves
- Class-specific colors
- Filled areas under curves
- Fallback to histogram if KDE fails
- By-class or overall view

**Use case:**
```python
fig = plot_motif_length_kde(
    motifs,
    by_class=True,
    title="Length Distribution (KDE)"
)
```

**Output:** `/tmp/test_length_kde.png` (284KB, 300 DPI)

---

## Integration in Streamlit App

### New Visualization Tabs

The visualizations are organized into **6 tabs** in the app:

1. **Distribution** - Class/subclass distribution, pie charts
2. **Coverage & Density** - Coverage map, density heatmap, Circos plot
3. **Statistics** - Density metrics, length/score distributions
4. **Genome-Wide** ⭐ NEW - Manhattan plot, cumulative distribution, linear track
5. **Advanced** ⭐ NEW - Co-occurrence matrix, GC correlation, length KDE, cluster stats
6. **Cluster/Hybrid** - Hybrid and cluster motif analysis

### Access in App

```python
# In app.py, after selecting "Analysis Results and Visualization" tab:
# Navigate to "Genome-Wide" or "Advanced" tabs to see new visualizations
```

---

## Technical Specifications

### Publication Quality

- **DPI:** 300 (Nature Methods requirement)
- **Format:** PDF/PNG/SVG supported
- **Color Palette:** Colorblind-friendly (Wong 2011 palette)
- **Font:** Arial/Helvetica sans-serif
- **Style:** Clean, minimal design with proper spacing

### Performance

- **Memory:** Optimized for large sequences with chunking
- **Speed:** Fast rendering with matplotlib/seaborn
- **Caching:** Results cached in Streamlit for efficiency

### Dependencies

```python
matplotlib >= 3.5.0
seaborn >= 0.11.0
numpy >= 1.21.0
scipy >= 1.7.0
pandas >= 1.3.0
```

---

## Testing

All visualizations have been tested with:

- ✅ 520 test motifs in 100kb sequence
- ✅ All 8 motif classes represented
- ✅ Cluster motifs included
- ✅ GC content variation
- ✅ No errors or warnings
- ✅ File sizes: 88KB - 284KB per plot

Test script: `test_new_visualizations.py`

Run tests:
```bash
python test_new_visualizations.py
```

---

## Example Usage in Research

### Manuscript Figure Preparation

```python
from visualizations import (
    plot_manhattan_motif_density,
    plot_motif_cooccurrence_matrix,
    plot_gc_content_correlation
)

# Figure 1: Genome-wide motif distribution
fig1 = plot_manhattan_motif_density(motifs, seq_len)
fig1.savefig('figure1_manhattan.pdf', dpi=300, format='pdf')

# Figure 2: Motif interactions
fig2 = plot_motif_cooccurrence_matrix(motifs)
fig2.savefig('figure2_cooccurrence.pdf', dpi=300, format='pdf')

# Figure 3: GC content relationship
fig3 = plot_gc_content_correlation(motifs, sequence)
fig3.savefig('figure3_gc_correlation.pdf', dpi=300, format='pdf')
```

### Batch Export

Use `save_all_plots()` function to generate all visualizations at once:

```python
from visualizations import save_all_plots

saved_files = save_all_plots(
    motifs, 
    sequence_length,
    output_dir="publication_figures",
    file_format="pdf",
    dpi=300,
    style='nature'
)
```

---

## References

**Visualization Design:**
- Nature Methods author guidelines for figure preparation
- Wong, B. (2011) Nature Methods colorblind-safe palette
- Tufte, E. (2001) The Visual Display of Quantitative Information

**Scientific Basis:**
- Non-B DNA structures in genomic biology
- Motif co-occurrence and functional relationships
- GC content influence on DNA structure formation

---

## Future Enhancements

Potential additional visualizations (Phase 2):

- G-quadruplex probability track
- Z-DNA energy profile
- R-loop coverage plot
- Triplex stability plot  
- Repeat structure visualization (STRs)
- Chromosome ideogram tracks (UCSC-style)
- Hybrid interaction network graph
- Interactive Plotly versions with hover tooltips

---

## Support

For questions or issues:
- Check documentation: `COMPREHENSIVE_DOCUMENTATION.md`
- View examples: `test_new_visualizations.py`
- Report issues: GitHub Issues

---

**Last Updated:** 2024-12-09
**Version:** 2024.1
**Author:** Dr. Venkata Rajesh Yella
