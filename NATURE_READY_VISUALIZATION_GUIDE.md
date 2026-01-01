# Nature-Ready Visualization Standardization Guide

## Overview

This document describes the comprehensive visualization standardization implemented in NonBDNAFinder 2025.1 to meet publication requirements of top-tier journals (Nature, Science, Cell).

**Version**: 2025.1  
**Author**: Dr. Venkata Rajesh Yella  
**Date**: January 2026

---

## 🎯 Core Objectives

### 1. Eliminate Redundant Plots Automatically ✓
- **Before**: 25+ visualization functions, many showing overlapping information
- **After**: 8 core plots in 3 canonical figures, zero redundancy

### 2. Expose Only Biologically Interpretable Metrics ✓
- **Before**: 80+ columns including raw thermodynamic values
- **After**: 10 core columns + motif-specific conditionals

### 3. Guarantee Non-Overlapping, Uncluttered Figures ✓
- **Before**: Individual motif labels causing overlap
- **After**: Only top 5% hotspots labeled with collision detection

### 4. Produce Journal-Ready Panels ✓
- **Before**: Feature-rich explorer with 11 distinct colors
- **After**: Publication-first instrument with 6-color palette

### 5. Preserve Advanced Outputs ✓
- **Before**: All data in UI and exports
- **After**: Core data in UI, full data in CSV/Excel exports

---

## 📊 Redundancy Elimination Logic

### Metric De-duplication Rules

| Concept | KEEP (UI) | HIDE (Export only) |
|---------|-----------|-------------------|
| Quality | Score (1-3) | Raw ΔG, Percentile, Unnormalized |
| Position | Start, End | Intermediate calculations |
| Structure | Length, Strand | Raw counts, pre-filtered values |

**Rule**: Only normalized, cross-class comparable metrics are shown in UI.

### Plot Dominance Rules

| Biological Concept | Allowed Plot | Hidden Plots |
|-------------------|--------------|--------------|
| **Composition** | Nested donut | Bar charts, sunbursts, separate pie charts |
| **Position** | Manhattan OR Linear track | Coverage maps, stacked tracks, circos |
| **Density** | Density comparison | Multiple KDEs, heatmaps, radial plots |
| **Clustering** | Cluster distribution | Multiple histograms, cumulative plots |
| **Co-occurrence** | Co-occurrence matrix | Individual overlap counts |
| **Length** | Single KDE | Histograms, box plots per class |

**Rule**: For each biological concept, only ONE visualization in main UI.

---

## 🎨 Color Palette Standardization

### Nature-Ready Palette (6 Colors + Neutrals)

#### Primary Classes (6 Unique Colors)
```python
'Curved_DNA'     : '#CC79A7'  # Reddish Purple - Structure
'G-Quadruplex'   : '#0072B2'  # Blue - Stable structures
'Z-DNA'          : '#882255'  # Wine - Alternative helices
'Cruciform'      : '#56B4E9'  # Sky Blue - Symmetric structures
'Triplex'        : '#E69F00'  # Orange - Triple-stranded
'R-Loop'         : '#009E73'  # Bluish Green - RNA-DNA hybrids
```

#### Secondary Classes (Consolidated)
```python
'i-Motif'        : '#0072B2'  # Same as G4 (complementary)
'A-philic_DNA'   : '#CC79A7'  # Same as Curved (structural affinity)
'Slipped_DNA'    : '#E69F00'  # Same as Triplex (repeats)
```

#### Meta-Classes (Neutral)
```python
'Hybrid'                 : '#888888'  # Medium gray
'Non-B_DNA_Clusters'     : '#666666'  # Dark gray
```

### Rationale for Consolidation

| Secondary Class | Consolidated To | Scientific Rationale |
|----------------|-----------------|---------------------|
| i-Motif | G-Quadruplex | Complementary structures (C-rich vs G-rich) |
| A-philic_DNA | Curved_DNA | Both involve DNA bending/structural properties |
| Slipped_DNA | Triplex | Both involve repeat-mediated structures |

**Validation**: `test_visualization_standards.py` verifies ≤6 unique colors

---

## 📐 Figure Panel Layout (Canonical)

### Figure 1: Global Non-B DNA Landscape
**Purpose**: *What structures exist, and where?*

- **Panel A**: Motif Composition (Nested donut: Class → Subclass)
- **Panel B**: Genome-Scale Localization
  - Manhattan plot (sequences >50kb)
  - Linear track (sequences ≤50kb)
- **Panel C**: Genome Coverage (% per motif class, compact bar)

### Figure 2: Structural Clustering & Co-Occurrence  
**Purpose**: *Which structures co-localize and form regulatory hotspots?*

- **Panel D**: Cluster Size Distribution (conditional on clusters)
- **Panel E**: Motif Co-occurrence Matrix (upper triangle heatmap)

### Figure 3: Structural Constraints (Optional/Toggle)
**Purpose**: *What are the physical constraints shaping each motif class?*

- **Panel F**: Length Distributions by Class (single KDE, shared axis)
- **Additional**: Score distribution (optional)

**Default UI**: 6 panels maximum (5 if no clusters detected)

---

## 🚫 Overlap Suppression Implementation

### Label Rendering Policy

```python
RENDER_INDIVIDUAL_LABELS = False      # Never label individual motifs
RENDER_CLUSTER_CENTROIDS = True       # Only cluster centers
RENDER_TOP_QUANTILE = True            # Only top 5% density hotspots

MIN_LABEL_DISTANCE_BP = sequence_length / 20  # Minimum spacing
MAX_LABELS_PER_PLOT = 10              # Hard limit
LABEL_PRIORITY_THRESHOLD = 0.95       # Top 5% by density/score
```

### Implementation Details

#### Manhattan Plot
- Only top 5% density windows are labeled
- Minimum spacing: 5% of sequence length
- Maximum 10 labels per plot
- Collision detection: programmatic distance checking

#### Linear Track
- **Removed**: Individual motif score labels
- **Kept**: Class track labels (left side, always shown)
- **Rationale**: Prevents overlap in dense regions

### Before vs After

| Aspect | Before | After |
|--------|--------|-------|
| Labels per plot | 50-500+ | ≤10 |
| Overlap incidents | Frequent | Zero (enforced) |
| Readability | Low (cluttered) | High (clean) |

---

## 🔬 Scientific Transparency

### Transparency Badge
Every visualization section includes:

```markdown
📊 **Scientific Transparency**: Only biologically interpretable, 
non-redundant metrics are displayed. Full results including raw 
thermodynamic values and all intermediate calculations are available 
in CSV/Excel exports.
```

### Supplementary Materials Note
```markdown
💡 **Supplementary Figures**: Additional visualizations (circos plots, 
radial tracks, cumulative distributions) are available on request for 
supplementary materials.
```

### Export Integrity
- **UI Display**: 10 core columns + motif-specific
- **CSV/Excel**: ALL columns including raw values
- **BED**: Genomic coordinates only
- **JSON**: Complete metadata

**Guarantee**: No data is lost, only UI display is filtered

---

## 📏 UI Compactness Improvements

### Layout Changes

| Aspect | Before | After | Improvement |
|--------|--------|-------|-------------|
| Vertical padding | 16px | 8px | 50% reduction |
| Tab structure | 3 sprawling tabs | 3 focused figure tabs | Clear hierarchy |
| Empty spacers | Frequent | None | Cleaner flow |
| Plot arrangement | Single column | Context-dependent | Better use of space |

### Figure Tab Structure

```
📊 Figure 1: Global Landscape
   └─ Panel A, B, C (always shown)

🔗 Figure 2: Clustering & Co-occurrence
   └─ Panel D (conditional), E (always shown)

📏 Figure 3: Structural Constraints (Optional)
   └─ Panel F + additional (user toggle)
```

---

## ✅ Validation & Testing

### Test Suite: `test_visualization_standards.py`

#### Test Coverage
- [x] Color count (≤6 unique colors)
- [x] Color consolidation (secondary → primary mapping)
- [x] Plot dominance (redundant plots hidden)
- [x] Required vs optional plots
- [x] Metric filtering (core vs hidden vs conditional)
- [x] Label suppression policy
- [x] Context-dependent display logic
- [x] Validation thresholds
- [x] Integration with utilities.py
- [x] Analysis pipeline integrity

#### Test Results
```
======================================================================
NATURE-READY VISUALIZATION STANDARDS - VALIDATION TESTS
======================================================================
TEST SUMMARY: 10 passed, 0 failed
======================================================================
```

### Manual Validation Checklist

- [x] App loads without errors
- [x] Analysis produces identical results
- [x] Visualizations render correctly
- [x] Colors consolidated as expected
- [x] No text overlap observed
- [x] Export files unchanged
- [x] UI more compact and navigable

---

## 🔄 Migration Guide

### For Users

**No action required!** The changes are entirely backend and improve the UI automatically.

**What you'll notice**:
1. Cleaner, more focused visualizations
2. Only 3 main figure tabs instead of sprawling sections
3. Consistent color scheme across all plots
4. No overlapping labels
5. Optional Figure 3 can be toggled on/off

**What stays the same**:
- All data still available in exports
- Analysis results identical
- Same motif detection algorithms
- Export formats unchanged

### For Developers

**Updated modules**:
1. `visualization_standards.py` - New centralized config
2. `utilities.py` - Updated color scheme and label suppression
3. `app.py` - Refactored visualization tabs

**Key imports**:
```python
from visualization_standards import (
    NATURE_MOTIF_COLORS,
    PlotDominance,
    FigurePanel,
    MetricFilter,
    LabelPolicy,
    should_show_plot
)
```

**Breaking changes**: None (backward compatible)

---

## 📚 References

### Journal Guidelines
- **Nature**: Figure preparation guidelines (2024)
- **Science**: Graphics and multimedia standards
- **Cell**: Author instructions for figures

### Scientific Best Practices
- **Broad Institute**: Data visualization standards
- **EMBL-EBI**: Genome visualization guidelines
- **Max Planck**: Publication figure requirements

### Colorblind-Friendly Design
- Wong, B. (2011) *Points of view: Color blindness*. Nat Methods 8, 441

---

## 🎓 Training & Support

### Quick Reference

**Q**: Where did my favorite plot go?  
**A**: Likely moved to supplementary (available on request) or consolidated with a more informative plot.

**Q**: Why are colors shared between some classes?  
**A**: Scientific grouping (i-Motif with G4, etc.) reduces cognitive load and meets journal color limits.

**Q**: Can I still get all the data?  
**A**: Yes! Full data in CSV/Excel exports, only UI display is streamlined.

**Q**: How do I enable Figure 3?  
**A**: Check the box "Include Figure 3 in main report" in the third tab.

### Support Contacts
- **Technical Issues**: yvrajesh_bt@kluniversity.in
- **GitHub**: https://github.com/VRYella/NonBDNAFinder
- **Documentation**: See README.md and OUTPUT_SCHEMA.md

---

## 📝 Version History

### v2025.1 (January 2026)
- ✅ Implemented Nature-ready visualization standards
- ✅ Reduced color palette from 11 to 6 unique colors
- ✅ Eliminated redundant plots (25 → 8 core plots)
- ✅ Added label suppression with collision detection
- ✅ Created 3-figure canonical layout
- ✅ Added scientific transparency badges
- ✅ Comprehensive test suite (10/10 tests passing)
- ✅ Zero breaking changes (fully backward compatible)

---

*This implementation converts NonBDNAFinder from a feature-rich explorer into a reviewer-ready scientific instrument, aligned with best practices from leading research institutions.*
