# 📊 Visualization Examples

This directory contains publication-quality sample outputs from the new Nature-level genome-wide visualizations implemented for the Non-B DNA Finder tool.

## Sample Outputs

All images are generated at **300 DPI** following Nature Methods styling guidelines.

### 1. Manhattan Plot (`test_manhattan.png`)
**Size:** 210KB | **Resolution:** 2103×843px

Genome-wide visualization showing motif density hotspots across the sequence. Color-coded by motif class to identify cluster regions and hybrid zones. Ideal for large-scale genome analysis.

---

### 2. Cumulative Distribution (`test_cumulative.png`)
**Size:** 209KB | **Resolution:** 2103×843px

Running sum plot showing how different motif classes accumulate across the genome. Each line represents a motif class, making it easy to compare accumulation patterns and identify regions of rapid increase.

---

### 3. Co-occurrence Matrix (`test_cooccurrence.png`)
**Size:** 179KB | **Resolution:** 1401×993px

Heatmap showing which motif classes tend to co-occur or overlap. Cell values indicate co-occurrence counts. Excellent for publication figures demonstrating motif interactions and relationships.

---

### 4. GC Content Correlation (`test_gc_correlation.png`)
**Size:** 126KB | **Resolution:** 1401×993px

Scatter plot with regression line showing the relationship between GC content and motif density. Includes correlation coefficient (R) to quantify the relationship. Demonstrates GC-driven motif enrichment patterns.

---

### 5. Linear Motif Track (`test_linear_track.png`)
**Size:** 88KB | **Resolution:** 2103×843px

Horizontal genome browser view showing motif positions as colored blocks. Best for detailed visualization of small regions (<10kb). Each motif class has its own track with class-specific colors.

---

### 6. Cluster Size Distribution (`test_cluster_dist.png`)
**Size:** 102KB | **Resolution:** 1401×993px

Dual histogram showing:
- Left: Distribution of motif counts per cluster
- Right: Distribution of class diversity within clusters

Includes mean/median annotations for statistical reference.

---

### 7. Motif Length KDE (`test_length_kde.png`)
**Size:** 284KB | **Resolution:** 2103×843px

Kernel density estimation showing smooth probability curves for motif length distributions. Each curve represents a motif class, making it easy to compare length patterns and identify modal lengths.

---

## Test Dataset

All visualizations were generated using a synthetic test dataset:

- **Sequence Length:** 100kb
- **Total Motifs:** 520
  - Regular motifs: 500
  - Cluster motifs: 20
- **Motif Classes:** 8 (G-Quadruplex, Z-DNA, Curved_DNA, i-Motif, Triplex, R-Loop, Cruciform, Slipped_DNA)
- **GC Content:** Varying (40-60% in windows)

## Usage

These images serve as:

1. **Reference Examples** - Show expected output quality and style
2. **Documentation** - Illustrate what each visualization type produces
3. **Testing** - Verify visualization functions work correctly
4. **Publications** - Demonstrate publication-ready quality

## Regenerating Outputs

To regenerate these visualizations:

```bash
cd /path/to/NonBDNAFinder
python test_new_visualizations.py
```

Outputs will be saved to `/tmp/test_*.png`

## Citation

If you use these visualizations in your research, please cite:

```
NonBScanner - Professional Non-B DNA Motif Detection Suite
Version 2024.1.3
Author: Dr. Venkata Rajesh Yella
```

---

**Note:** All visualizations follow:
- Nature Methods author guidelines for figures
- Wong (2011) colorblind-safe palette
- 300 DPI minimum for publication
- Clean, minimal professional design
- Proper axis labels and statistical annotations

**Generated:** December 9, 2024
**Test Script:** `test_new_visualizations.py`
**Documentation:** `VISUALIZATION_GUIDE.md`
