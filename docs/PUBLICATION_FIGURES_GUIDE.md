# Publication-Ready Multi-Panel Figures Guide

## Overview

This guide describes how to create Nature-level publication-ready multi-panel figures for scientific manuscripts using the NBDScanner visualization suite.

## Table of Contents

1. [Quick Start](#quick-start)
2. [Figure Standards](#figure-standards)
3. [Multi-Panel Figures](#multi-panel-figures)
4. [Customization](#customization)
5. [Export Options](#export-options)
6. [Submission Checklist](#submission-checklist)

---

## Quick Start

### Basic Usage

```python
from publication_figures import (
    create_figure_1_overview,
    create_figure_2_statistics,
    create_figure_3_landscape
)
import nonbscanner as nbs

# Analyze your sequence
sequence = "ATCGATCG..."  # Your DNA sequence
motifs = nbs.analyze_sequence(sequence, "sequence_name")

# Generate publication figures
fig1 = create_figure_1_overview(motifs, len(sequence))
fig1.savefig('figure1.pdf', dpi=300, bbox_inches='tight')

fig2 = create_figure_2_statistics(motifs)
fig2.savefig('figure2.pdf', dpi=300, bbox_inches='tight')

fig3 = create_figure_3_landscape(motifs, len(sequence))
fig3.savefig('figure3.pdf', dpi=300, bbox_inches='tight')
```

### Command-Line Generation

```bash
# Generate all publication figures automatically
python generate_publication_figures.py --output-dir ./figures --format pdf

# Options:
#   --output-dir DIR    Output directory (default: ./publication_figures)
#   --format FORMAT     pdf, svg, or png (default: pdf)
#   --dpi DPI          Resolution for raster formats (default: 300)
```

---

## Figure Standards

### Nature Journal Requirements

All figures follow Nature Methods/Genetics submission guidelines:

| Requirement | Specification | Implementation |
|------------|---------------|----------------|
| **Resolution** | 300 DPI minimum | ✓ PUBLICATION_DPI = 300 |
| **Format** | PDF/SVG (vector) preferred | ✓ Default: PDF |
| **Fonts** | Arial/Helvetica, 5-8pt | ✓ Sans-serif, 7-9pt |
| **Column Width** | Single: 89mm (3.5"), Double: 183mm (7.2") | ✓ NATURE_SINGLE/DOUBLE_COLUMN |
| **Color Palette** | Colorblind-friendly | ✓ Wong (2011) palette |
| **Panel Labels** | Bold, uppercase (A, B, C...) | ✓ 12pt bold |
| **Line Width** | 0.5-1pt | ✓ 0.5-1.5pt |
| **File Size** | <10 MB per figure | ✓ Optimized output |

### Colorblind-Safe Palette

Based on Wong, B. (2011) *Nature Methods* 8:441

```python
MOTIF_CLASS_COLORS = {
    'Curved_DNA': '#CC79A7',         # Reddish Purple
    'Slipped_DNA': '#E69F00',        # Orange
    'Cruciform': '#56B4E9',          # Sky Blue
    'R-Loop': '#009E73',             # Bluish Green
    'Triplex': '#F0E442',            # Yellow
    'G-Quadruplex': '#0072B2',       # Blue
    'i-Motif': '#D55E00',            # Vermillion
    'Z-DNA': '#882255',              # Wine
    'A-philic_DNA': '#44AA99',       # Teal
    'Hybrid': '#999999',             # Gray
    'Non-B_DNA_Clusters': '#666666'  # Dark Gray
}
```

**Verification**: All colors tested with Coblis (Color Blindness Simulator) for deuteranopia, protanopia, and tritanopia.

---

## Multi-Panel Figures

### Figure 1: Comprehensive Overview (2×2 Layout)

**Purpose**: Provide complete overview of Non-B DNA motif analysis

**Panels**:
- **A**: Motif class distribution (bar chart)
- **B**: Class-subclass hierarchy (nested donut)
- **C**: Genomic coverage map (linear track)
- **D**: [Reserved for custom analysis]

**Dimensions**: 7.2" × 7" (double column)

```python
fig1 = create_figure_1_overview(
    motifs=motifs,
    sequence_length=len(sequence),
    title="Comprehensive Non-B DNA Motif Analysis",
    save_path="figure1_overview.pdf"
)
```

**Key Features**:
- Shows all 11 Non-B DNA classes (even if count = 0)
- Consistent colors across all panels
- Automatic count labeling on bars
- Class labels without underscores (e.g., "Slipped DNA" not "Slipped_DNA")

---

### Figure 2: Statistical Characterization (2×2 Layout)

**Purpose**: Present statistical properties of detected motifs

**Panels**:
- **A**: Score distribution by class (violin plots)
- **B**: Length distribution by class (box plots)
- **C**: Score-length correlation (scatter plot)
- **D**: Summary statistics table

**Dimensions**: 7.2" × 7" (double column)

```python
fig2 = create_figure_2_statistics(
    motifs=motifs,
    title="Statistical Characterization of Non-B DNA Motifs",
    save_path="figure2_statistics.pdf"
)
```

**Key Features**:
- Violin plots show full distribution shape
- Inner box plots indicate quartiles
- Scatter plot reveals correlations
- Summary table with n, mean score, mean length

---

### Figure 3: Genomic Landscape (3×1 Layout)

**Purpose**: Visualize spatial distribution and density patterns

**Panels**:
- **A**: Density heatmap across genome
- **B**: Manhattan plot showing hotspots
- **C**: Cumulative distribution curves

**Dimensions**: 7.2" × 9" (double column, tall)

```python
fig3 = create_figure_3_landscape(
    motifs=motifs,
    sequence_length=len(sequence),
    window_size=1000,  # Optional, auto-calculated if None
    title="Genomic Landscape and Density Analysis",
    save_path="figure3_landscape.pdf"
)
```

**Key Features**:
- Heatmap uses viridis colormap (colorblind-safe)
- Manhattan plot highlights density peaks
- Cumulative curves show temporal accumulation
- Window size auto-optimizes for sequence length

---

## Customization

### Creating Custom Multi-Panel Figures

```python
from publication_figures import create_multipanel_figure, add_panel_label
import matplotlib.pyplot as plt

# Create custom 2×3 layout
fig, axes = create_multipanel_figure(
    n_rows=2, 
    n_cols=3,
    figsize=(7.2, 6),  # Double column width
    hspace=0.4,        # Vertical spacing
    wspace=0.3         # Horizontal spacing
)

# Add your plots to each panel
add_panel_label(axes[0, 0], 'A')
# ... plot in axes[0, 0] ...

add_panel_label(axes[0, 1], 'B')
# ... plot in axes[0, 1] ...

# etc.

fig.suptitle('Custom Multi-Panel Figure', fontsize=11, fontweight='bold')
plt.tight_layout()
fig.savefig('custom_figure.pdf', dpi=300, bbox_inches='tight')
```

### Adjusting Figure Dimensions

```python
from publication_figures import (
    NATURE_SINGLE_COLUMN,      # 3.5"
    NATURE_ONE_AND_HALF,       # 4.7"
    NATURE_DOUBLE_COLUMN       # 7.2"
)

# Single column figure
fig, ax = plt.subplots(figsize=(NATURE_SINGLE_COLUMN, 3))

# One-and-a-half column figure
fig, ax = plt.subplots(figsize=(NATURE_ONE_AND_HALF, 4))

# Double column figure
fig, ax = plt.subplots(figsize=(NATURE_DOUBLE_COLUMN, 5))
```

### Custom Color Schemes

While the default colorblind-safe palette is recommended, you can customize:

```python
# Use custom colors (ensure colorblind safety!)
custom_colors = {
    'G-Quadruplex': '#1f77b4',
    'Z-DNA': '#ff7f0e',
    # ... other classes
}

# Pass to plotting functions
from visualizations import plot_motif_distribution

fig = plot_motif_distribution(
    motifs, 
    by='Class',
    # Colors will be automatically applied from MOTIF_CLASS_COLORS
)
```

---

## Export Options

### PDF (Recommended for Publication)

```python
fig.savefig(
    'figure.pdf',
    format='pdf',
    dpi=300,
    bbox_inches='tight',
    pad_inches=0.05,
    facecolor='white',
    edgecolor='none'
)
```

**Advantages**:
- Vector format (infinite scalability)
- Editable text in Adobe Illustrator/Inkscape
- Small file size
- Preferred by most journals

### SVG (Alternative Vector Format)

```python
fig.savefig(
    'figure.svg',
    format='svg',
    dpi=300,
    bbox_inches='tight'
)
```

**Advantages**:
- Web-friendly vector format
- Editable in Inkscape, Adobe Illustrator
- Good for supplementary materials

### PNG (High-Resolution Raster)

```python
fig.savefig(
    'figure.png',
    format='png',
    dpi=300,  # or 600 for extra quality
    bbox_inches='tight',
    facecolor='white'
)
```

**Advantages**:
- Universal compatibility
- Good for presentations
- Larger file sizes

**Note**: Always use 300 DPI minimum for print publication (600 DPI for line art).

---

## Submission Checklist

### Before Submitting to Journal

- [ ] **Resolution**: All figures at 300+ DPI
- [ ] **Format**: PDF or SVG (vector preferred)
- [ ] **Dimensions**: Match journal column widths
  - [ ] Single column: 89mm (3.5")
  - [ ] Double column: 183mm (7.2")
- [ ] **Fonts**: Arial/Helvetica family, readable sizes
  - [ ] Panel labels: 12pt bold
  - [ ] Axis labels: 8pt
  - [ ] Tick labels: 7pt
- [ ] **Colors**: Colorblind-safe palette verified
- [ ] **Panel Labels**: Bold uppercase letters (A, B, C...)
- [ ] **Axis Labels**: Include units (bp, kb, %, etc.)
- [ ] **Legends**: Clear, concise, properly positioned
- [ ] **File Size**: <10 MB per figure
- [ ] **Naming**: Consistent naming scheme (figure1.pdf, figure2.pdf...)

### Figure Captions

Each figure should have comprehensive caption including:

1. **Title**: Brief descriptive title
2. **Panel descriptions**: (A) Description of panel A. (B) Description of panel B...
3. **Methods**: Key analysis parameters
4. **Statistics**: Sample sizes (n=X), significance levels
5. **Scale information**: Units, time points, etc.

**Example**:

```
Figure 1. Comprehensive Non-B DNA Motif Analysis

Multi-panel overview of Non-B DNA motif detection in test genomic 
sequence (n=51 motifs). (A) Distribution of detected motifs across 
11 major Non-B DNA classes. (B) Hierarchical composition showing 
class-subclass relationships. (C) Genomic coverage map displaying 
motif positions along analyzed sequence. All analyses performed 
using NBDScanner v2024.1. Colors follow Wong (2011) colorblind-safe 
palette. Scale bars as indicated.
```

### Quality Verification

1. **Print test**: Print figure at actual size, check readability
2. **Grayscale test**: Convert to grayscale, verify distinguishability
3. **Colorblind simulation**: Use tools like Coblis to verify accessibility
4. **Zoom test**: Check details at 100%, 200%, 400% zoom
5. **File size**: Confirm <10 MB (optimize if needed)

---

## Advanced Topics

### Combining Multiple Analyses

```python
# Compare multiple samples
samples = {
    'Sample A': motifs_a,
    'Sample B': motifs_b,
    'Sample C': motifs_c
}

# Create comparison figure
fig, axes = create_multipanel_figure(3, 1, figsize=(7.2, 9))

for i, (sample_name, motifs) in enumerate(samples.items()):
    ax = axes[i, 0]
    add_panel_label(ax, chr(65 + i))  # A, B, C...
    
    # Plot distribution for this sample
    # ... custom plotting code ...
    
    ax.set_title(sample_name)

plt.tight_layout()
fig.savefig('comparison_figure.pdf', dpi=300, bbox_inches='tight')
```

### Statistical Annotations

```python
import matplotlib.pyplot as plt
from scipy import stats

# Add significance stars
def add_significance_bar(ax, x1, x2, y, p_value, height=0.05):
    """Add significance bar with stars."""
    if p_value < 0.001:
        stars = '***'
    elif p_value < 0.01:
        stars = '**'
    elif p_value < 0.05:
        stars = '*'
    else:
        stars = 'n.s.'
    
    ax.plot([x1, x1, x2, x2], [y, y+height, y+height, y], 'k-', lw=0.8)
    ax.text((x1+x2)/2, y+height, stars, ha='center', va='bottom', fontsize=8)

# Usage in your figure
ax.bar(...)  # your bars
add_significance_bar(ax, 0, 1, max_val*1.1, p_value=0.003)
```

---

## Troubleshooting

### Font Issues

**Problem**: "Font family 'Arial' not found"

**Solution**: System falls back to DejaVu Sans (equivalent). To use exact Arial:

```bash
# Install Microsoft fonts (Linux)
sudo apt-get install ttf-mscorefonts-installer

# Or use Helvetica (macOS has this)
# Or accept DejaVu Sans fallback (very similar)
```

### Figure Too Large

**Problem**: File size >10 MB

**Solutions**:
1. Use PDF/SVG (vector) instead of PNG
2. Reduce number of data points in scatter plots
3. Simplify complex elements
4. Use rasterization for dense plots: `ax.set_rasterized(True)`

### Panel Labels Overlap

**Problem**: Panel label 'A' overlaps with plot

**Solution**: Adjust label position

```python
add_panel_label(ax, 'A', x=-0.15, y=1.1)  # Move left and up
```

### Colors Not Colorblind-Safe

**Problem**: Custom colors not distinguishable

**Solution**: Use online tool to verify:
- [Coblis - Color Blindness Simulator](https://www.color-blindness.com/coblis-color-blindness-simulator/)
- [Adobe Color Accessibility Tool](https://color.adobe.com/create/color-accessibility)

---

## References

1. **Wong, B.** (2011) Points of view: Color blindness. *Nature Methods* 8:441.
   - Establishes colorblind-safe palette used in this toolkit

2. **Nature Methods** - Guide to Publication Figures
   - https://www.nature.com/nature/for-authors/formatting-guide

3. **Rougier et al.** (2014) Ten Simple Rules for Better Figures. *PLOS Comp. Biol.* 10(9):e1003833.
   - Best practices for scientific figures

4. **Cleveland & McGill** (1984) Graphical Perception. *J. Amer. Stat. Assoc.* 79(387):531-554.
   - Visual encoding principles

---

## Support

For questions, issues, or feature requests:
- **GitHub Issues**: https://github.com/VRYella/NonBDNAFinder/issues
- **Documentation**: See repository README.md
- **Examples**: Check `docs/visualization_examples/`

---

## License

This visualization toolkit is part of NBDScanner and is released under the MIT License.

---

*Last updated: December 2024*
*Version: 2024.1.5*
