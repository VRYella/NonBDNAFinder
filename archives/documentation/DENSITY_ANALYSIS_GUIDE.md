# Density and Enrichment Analysis Guide

## Overview

NonBScanner now provides comprehensive density and enrichment analysis at both **class** and **subclass** levels for Non-B DNA motifs. This guide explains how to use these features for publication-quality statistical analysis.

## Key Concepts

### 1. Genomic Density (σ_G)

**Definition**: The percentage of the analyzed sequence covered by motifs.

**Formula**: 
```
σ_G = (Total unique bp covered by motifs / Total length in bp) × 100
```

**Important Notes**:
- Uses set-based overlap handling
- Coverage never exceeds 100%, even with overlapping motifs
- Only unique positions are counted

### 2. Positional Density (λ)

**Definition**: The frequency of motifs per unit length.

**Formula**:
```
λ = Total count of predicted motifs / Total length (in kbp or Mbp)
```

**Unit Options**:
- `kbp` (kilobase pairs) - default for smaller sequences
- `Mbp` (megabase pairs) - better for large genomes

### 3. Enrichment Analysis

**Definition**: Statistical comparison of observed motif density versus background (shuffled sequences).

**Metrics**:
- **Fold Enrichment**: Observed density / Background mean
- **P-value**: Statistical significance (proportion of shuffled ≥ observed)
- **Background Statistics**: Mean and standard deviation from shuffling

## Python API Usage

### Basic Density Calculation

```python
import nonbscanner as nbs
from utilities import calculate_genomic_density, calculate_positional_density

# Analyze sequence
sequence = "GGGTTAGGGTTAGGGTTAGGG" * 100
motifs = nbs.analyze_sequence(sequence, "my_sequence")
sequence_length = len(sequence)

# Class-level density
genomic_density_class = calculate_genomic_density(
    motifs, sequence_length, 
    by_class=True, 
    by_subclass=False
)

positional_density_class = calculate_positional_density(
    motifs, sequence_length, 
    unit='kbp',
    by_class=True, 
    by_subclass=False
)

print(f"Genomic Density (Overall): {genomic_density_class['Overall']:.4f}%")
print(f"Positional Density (Overall): {positional_density_class['Overall']:.2f} motifs/kbp")

# Display per-class results
for class_name in sorted(genomic_density_class.keys()):
    if class_name != 'Overall':
        print(f"{class_name}:")
        print(f"  Genomic: {genomic_density_class[class_name]:.4f}%")
        print(f"  Positional: {positional_density_class[class_name]:.2f} motifs/kbp")
```

### Subclass-Level Analysis

```python
# Subclass-level density
genomic_density_subclass = calculate_genomic_density(
    motifs, sequence_length, 
    by_class=False,  # Turn off class-level
    by_subclass=True  # Turn on subclass-level
)

positional_density_subclass = calculate_positional_density(
    motifs, sequence_length, 
    unit='kbp',
    by_class=False, 
    by_subclass=True
)

# Display per-subclass results
# Subclass keys are in format 'Class:Subclass'
for subclass_key in sorted(genomic_density_subclass.keys()):
    if subclass_key != 'Overall':
        print(f"{subclass_key}:")
        print(f"  Genomic: {genomic_density_subclass[subclass_key]:.4f}%")
        print(f"  Positional: {positional_density_subclass[subclass_key]:.2f} motifs/kbp")
```

### Enrichment Analysis

```python
from utilities import calculate_enrichment_with_shuffling

# Class-level enrichment (100 shuffles, takes ~2-5 minutes)
enrichment_class = calculate_enrichment_with_shuffling(
    motifs, sequence, 
    n_shuffles=100,
    by_class=True,
    by_subclass=False
)

# Display enrichment results
for class_name, result in enrichment_class.items():
    if class_name != 'Overall':
        print(f"\n{class_name}:")
        print(f"  Fold Enrichment: {result['fold_enrichment']}")
        print(f"  P-value: {result['p_value']:.4f}")
        print(f"  Observed Density: {result['observed_density']:.4f}%")
        print(f"  Background Mean: {result['background_mean']:.4f}%")
        
        # Interpret significance
        if result['p_value'] < 0.05:
            print(f"  ✓ Significantly enriched (p < 0.05)")
        elif result['p_value'] < 0.1:
            print(f"  ⚠ Marginally significant (p < 0.1)")
        else:
            print(f"  - Not significantly enriched")

# Subclass-level enrichment
enrichment_subclass = calculate_enrichment_with_shuffling(
    motifs, sequence, 
    n_shuffles=100,
    by_class=False,
    by_subclass=True
)
```

## Visualization Functions

### Class-Level Density Plots

```python
from visualizations import plot_density_comparison
import matplotlib.pyplot as plt

# Create density comparison plot
fig = plot_density_comparison(
    genomic_density_class, 
    positional_density_class,
    title="Motif Density Analysis (Class Level)"
)

# Save at publication quality (300 DPI)
plt.savefig('density_class.png', dpi=300, bbox_inches='tight')
plt.close(fig)
```

### Subclass-Level Density Plots

```python
from visualizations import plot_density_comparison_by_subclass

# Create subclass density comparison
fig = plot_density_comparison_by_subclass(
    genomic_density_subclass, 
    positional_density_subclass,
    title="Motif Density Analysis (Subclass Level)"
)

plt.savefig('density_subclass.png', dpi=300, bbox_inches='tight')
plt.close(fig)
```

### Enrichment Visualization

```python
from visualizations import plot_enrichment_analysis_by_subclass

# Create enrichment analysis plot (class or subclass)
fig = plot_enrichment_analysis_by_subclass(
    enrichment_subclass,
    title="Motif Enrichment Analysis"
)

plt.savefig('enrichment.png', dpi=300, bbox_inches='tight')
plt.close(fig)
```

### Density Heatmap

```python
from visualizations import plot_subclass_density_heatmap

# Create heatmap showing subclass density along sequence
fig = plot_subclass_density_heatmap(
    motifs, sequence_length,
    window_size=1000,  # 1kb windows
    title="Subclass Density Distribution"
)

plt.savefig('density_heatmap.png', dpi=300, bbox_inches='tight')
plt.close(fig)
```

## Web Interface Usage

### Accessing Density Analysis

1. **Upload or paste your sequence** in the main interface
2. **Click "Analyze Sequence"** to detect motifs
3. **Navigate to "Statistics" tab** in the visualization section
4. **View both analysis levels** - Both class and subclass level density analysis are automatically displayed:
   - **Class Level**: Analysis grouped by major motif classes (11 classes)
   - **Subclass Level**: Detailed analysis by specific motif subtypes (22+ subclasses)

### Interpreting Results

#### Overall Metrics
- **Genomic Density**: Percentage of sequence covered by motifs
- **Motifs/kbp**: Frequency of motifs per kilobase pair

#### Per-Class/Subclass Tables
- Shows density metrics for each detected motif type at both class and subclass levels
- Sorted alphabetically for easy comparison
- Values formatted for readability (4 decimal places for %, 2 for frequency)

#### Density Comparison Charts
- **Panel A**: Genomic density (coverage %)
- **Panel B**: Positional density (motifs/kbp)
- Charts are provided for both class and subclass levels
- Color-coded by parent class (subclass plots use parent class colors)
- Values labeled on bars for precision

#### Subclass Density Heatmap
- Shows distribution of each subclass across the sequence
- Uses sliding windows (auto-calculated based on sequence length)
- Color intensity indicates motif count per window
- Useful for identifying hotspots and distribution patterns

## Publication-Quality Output

### Figure Specifications

All density and enrichment plots are designed for publication in scientific journals:

- **DPI**: 300 (Nature journal requirement)
- **Font Family**: Arial/Helvetica (sans-serif)
- **Font Sizes**: 
  - Title: 12-14pt
  - Axis labels: 11pt
  - Tick labels: 8-9pt
  - Legend: 9pt
- **Color Palette**: Colorblind-friendly (Wong 2011)
- **Style**: Clean Nature journal style (no top/right spines)
- **Format**: PDF (vector) or PNG (300 DPI raster)

### Saving Figures

```python
# Vector format (preferred for publications)
plt.savefig('figure.pdf', format='pdf', bbox_inches='tight', dpi=300)

# Raster format (for presentations)
plt.savefig('figure.png', format='png', bbox_inches='tight', dpi=300)

# SVG format (editable)
plt.savefig('figure.svg', format='svg', bbox_inches='tight')
```

## Best Practices

### 1. Choose Appropriate Analysis Level

- **Use Class-Level when**:
  - You want a broad overview of motif distribution
  - Working with shorter sequences (<10kb)
  - Comparing multiple sequences/samples
  
- **Use Subclass-Level when**:
  - You need detailed characterization
  - Studying specific motif variants
  - Publishing detailed motif analysis
  - Investigating functional differences between subtypes

### 2. Enrichment Analysis Tips

- **Shuffles**: Use 100+ shuffles for reliable p-values
- **Sequence Length**: Works best with sequences >5kb
- **Interpretation**: 
  - Fold Enrichment > 2: Strong enrichment
  - Fold Enrichment 1-2: Moderate enrichment
  - Fold Enrichment < 1: Depletion
  - P < 0.05: Statistically significant

### 3. Visualization Guidelines

- Use **class-level plots** for presentations (cleaner, easier to read)
- Use **subclass-level plots** for detailed publications
- Include **heatmaps** to show spatial distribution
- Always report **both genomic and positional density**

## Example Workflow

```python
import nonbscanner as nbs
from utilities import (
    calculate_genomic_density, 
    calculate_positional_density,
    calculate_enrichment_with_shuffling
)
from visualizations import (
    plot_density_comparison,
    plot_density_comparison_by_subclass,
    plot_enrichment_analysis_by_subclass,
    plot_subclass_density_heatmap
)
import matplotlib.pyplot as plt

# 1. Analyze sequence
sequence = open('my_sequence.fasta').read()
motifs = nbs.analyze_file('my_sequence.fasta')['my_sequence']
sequence_length = len(sequence)

print(f"Detected {len(motifs)} motifs in {sequence_length} bp sequence")

# 2. Calculate densities
genomic_class = calculate_genomic_density(motifs, sequence_length, by_class=True)
positional_class = calculate_positional_density(motifs, sequence_length, unit='kbp', by_class=True)

genomic_subclass = calculate_genomic_density(motifs, sequence_length, by_subclass=True)
positional_subclass = calculate_positional_density(motifs, sequence_length, unit='kbp', by_subclass=True)

# 3. Calculate enrichment
enrichment = calculate_enrichment_with_shuffling(
    motifs, sequence, n_shuffles=100, by_subclass=True
)

# 4. Generate figures
fig1 = plot_density_comparison(genomic_class, positional_class)
plt.savefig('density_class.pdf', dpi=300, bbox_inches='tight')
plt.close(fig1)

fig2 = plot_density_comparison_by_subclass(genomic_subclass, positional_subclass)
plt.savefig('density_subclass.pdf', dpi=300, bbox_inches='tight')
plt.close(fig2)

fig3 = plot_enrichment_analysis_by_subclass(enrichment)
plt.savefig('enrichment.pdf', dpi=300, bbox_inches='tight')
plt.close(fig3)

fig4 = plot_subclass_density_heatmap(motifs, sequence_length)
plt.savefig('heatmap.pdf', dpi=300, bbox_inches='tight')
plt.close(fig4)

print("✓ Analysis complete! Generated 4 publication-quality figures.")
```

## Troubleshooting

### Common Issues

**Q: Why is genomic density > 100%?**  
A: This should never happen with the new implementation. If you see this, please report it as a bug.

**Q: Enrichment analysis is taking too long**  
A: Reduce `n_shuffles` parameter (minimum 50 for reasonable estimates).

**Q: Subclass plot shows no data**  
A: Make sure you're using `by_subclass=True` when calculating densities.

**Q: Labels are overlapping in plots**  
A: This is handled automatically for subclass plots. If issues persist, try adjusting `figsize`.

## References

- **Genomic Density**: Based on standard coverage metrics used in genomics
- **Statistical Validation**: Permutation testing approach (Good, 2005)
- **Visualization Standards**: Nature Methods author guidelines
- **Color Palette**: Wong, B. (2011) Nature Methods colorblind-safe palette

## Support

For issues or questions:
- Check the main README.md
- Review example notebooks
- Contact: yvrajesh_bt@kluniversity.in
