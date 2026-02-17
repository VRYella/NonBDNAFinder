# NBST Validation Analysis - User Guide

## Overview

This directory contains comprehensive **subclass-level validation** of NonBDNAFinder against the Non-B DNA Structure Tool (NBST) benchmark. The analysis provides granular performance metrics for each NonBDNAFinder subclass (e.g., "Canonical G4", "Extended-loop G4", "Weak PQS") compared to NBST classes.

## Quick Start

### Prerequisites

```bash
# Install required Python packages
pip install pandas numpy matplotlib seaborn scipy matplotlib-venn biopython
```

### Running the Analysis

```bash
# 1. Run subclass-level validation analysis
python validation_analysis_subclass.py

# 2. Generate publication-quality figures (300 DPI)
python create_publication_figures.py
```

## Output Files

### Analysis Results

- **`subclass_comparison_summary.tsv`**: Comprehensive subclass performance metrics (TP, FP, FN, Precision, Recall, F1, Jaccard)
- **`nonbdnafinder_results.tsv`**: All NonBDNAFinder predictions (341 motifs)

### Tables (CSV format)

- **`table1_subclass_metrics.csv`**: Subclass-level performance metrics (primary results table)
- **`table2_subclass_match_quality.csv`**: Match quality statistics (mean/median Jaccard, overlap BP)
- **`table3_false_positive_analysis.csv`**: False positive characterization by subclass
- **`table4_false_negative_analysis.csv`**: False negative analysis (NBST motifs missed)

### Figures (300 DPI, colorblind-friendly)

- **`figure1_subclass_distribution.png`**: Class/subclass distribution comparison (pie charts)
- **`figure2_performance_metrics.png`**: Subclass performance bar charts (Precision, Recall, F1)
- **`figure3_genomic_coverage.png`**: Genomic coverage map with subclass differentiation
- **`figure4_detection_counts.png`**: Detection count comparison (NBST vs NBF)
- **`figure5_overlap_heatmap.png`**: Subclass-level overlap heatmap
- **`figure6_subclass_distribution.png`**: Comprehensive subclass breakdown
- **`figure7_precision_recall_curves.png`**: Precision-recall scatter plot per subclass
- **`figure8_venn_diagrams.png`**: Venn diagrams showing TP/FP/FN relationships

### Reports

- **`VALIDATION_REPORT_EXTENDED.md`**: Comprehensive validation report with:
  - Enhanced abstract with subclass-level quantification
  - Biological significance of subclass distinctions
  - Detailed results by major class with subclass breakdowns
  - Discussion of performance patterns and recommendations
  - Updated methods and references (2023-2026)

## Key Findings Summary

### High-Confidence Subclasses (Precision = 1.00)

**Best for experimental validation and functional studies:**

- **Canonical intramolecular G4**: 100% precision, perfect NBST agreement (Jaccard=1.00)
- **STR (Short Tandem Repeats)**: 100% precision, near-perfect agreement (Jaccard=0.95)
- **Direct Repeat**: 100% precision, perfect agreement (Jaccard=1.00)

**Recommendation**: Use these subclasses for CRISPR screens, ChIP validation, or mutagenesis studies where specificity is critical.

### Moderate-Performance Subclasses (F1 = 0.16-0.22)

**Suitable for hypothesis-generating studies:**

- **Extended-loop canonical G4**: F1=0.20, moderate NBST agreement
- **Intramolecular G-triplex**: F1=0.16, distinct three-stranded structures
- **Z-DNA**: F1=0.22, thermodynamic scoring approach

**Recommendation**: Use for discovery screens; validate candidates experimentally (CD spectroscopy, structural probes).

### Expanded Landscape Subclasses (High Count, Low Precision)

**Two-tetrad weak PQS** (122 motifs, 36% of all detections):
- Precision: 0.016 (very low NBST concordance)
- May represent biologically relevant transient G4s OR false positives
- **Filtering strategy**: Use score ≥1.5 and functional genomic context (promoters, enhancers)

### Divergent Subclasses (Zero Concordance)

**Curved DNA (Local/Global)**:
- No overlap with NBST curved predictions
- Represents **complementary annotations** (local bending vs. global phased curvature)
- Not competing predictions—use both for comprehensive curvature analysis

## Understanding the Results

### Performance Metrics

- **Precision**: What fraction of NBF predictions match NBST? (Higher = fewer false positives)
- **Recall**: What fraction of NBST predictions are detected by NBF? (Higher = fewer false negatives)
- **F1 Score**: Harmonic mean of precision and recall (balanced metric)
- **Jaccard Index**: Coordinate overlap quality (1.0 = perfect agreement, 0.5 = threshold)

### Interpreting Subclass Results

| Metric Pattern | Interpretation | User Action |
|----------------|----------------|-------------|
| High Precision (≥0.9), Low Recall (≤0.2) | Conservative detection, high confidence | Use for experimental validation |
| Balanced (Precision ≈ Recall, 0.2-0.5) | Moderate agreement | Suitable for discovery, validate candidates |
| Low Precision (<0.1), High Count | Expanded landscape OR noise | Apply score/context filtering |
| Zero concordance | Algorithmic divergence | Consider complementary annotations |

## Use Cases

### 1. High-Specificity Application (Functional Validation)

**Goal**: Identify high-confidence non-B DNA motifs for CRISPR editing experiments.

```bash
# Run analysis
python validation_analysis_subclass.py

# Filter results for high-precision subclasses
# Use: Canonical G4, STR, Direct Repeat
# Additional filter: Score ≥1.8, Length ≥25bp
```

**Expected Output**: Small set (~10-20 motifs) with near-zero false positive rate.

### 2. High-Sensitivity Application (Genome Annotation)

**Goal**: Comprehensive non-B DNA landscape for genome browser tracks.

```bash
# Run analysis with all subclasses
python validation_analysis_subclass.py

# Include all subclasses (Canonical + Extended-loop + Weak PQS + etc.)
# Post-process: Rank by score, annotate functional context
```

**Expected Output**: Large set (300+ motifs) requiring downstream filtering/prioritization.

### 3. Discovery Screen (Identifying Novel Structures)

**Goal**: Find potential weak G4 motifs in regulatory regions.

```bash
# Run analysis
python validation_analysis_subclass.py

# Focus on: Two-tetrad weak PQS
# Filters: Score ≥1.5, TSS ±1kb, enhancers
```

**Expected Output**: ~20-30 high-priority weak PQS candidates for experimental validation.

## Advanced Usage

### Adjusting Matching Threshold

Edit `validation_analysis_subclass.py`:

```python
# Default threshold (Jaccard ≥ 0.5)
matches, nbf_matched, nbst_matched = match_motifs_subclass_level(
    nbf_motifs, nbst_df, nbf_subclass, overlap_threshold=0.5
)

# More stringent (Jaccard ≥ 0.7 for high-quality matches only)
matches, nbf_matched, nbst_matched = match_motifs_subclass_level(
    nbf_motifs, nbst_df, nbf_subclass, overlap_threshold=0.7
)
```

### Customizing Figure Generation

Edit `create_publication_figures.py`:

```python
# Adjust DPI (default 300 for publication)
plt.rcParams['figure.dpi'] = 600  # Ultra-high resolution

# Change color palette
COLORBLIND_PALETTE = ['#your', '#custom', '#colors']

# Adjust figure size
fig, ax = plt.subplots(figsize=(20, 12))  # Larger figures
```

### Filtering Weak PQS by Score

```python
import pandas as pd

# Load results
df = pd.read_csv('nonbdnafinder_results.tsv', sep='\t')

# Filter weak PQS by score threshold
weak_pqs = df[df['Subclass'] == 'Two-tetrad weak PQS']
high_confidence_pqs = weak_pqs[weak_pqs['Score'] >= 1.5]

print(f"High-confidence weak PQS: {len(high_confidence_pqs)} / {len(weak_pqs)}")
```

## Troubleshooting

### Issue: Module not found errors

```bash
# Install missing dependencies
pip install pandas numpy matplotlib seaborn scipy matplotlib-venn biopython

# Or use requirements from parent directory
pip install -r ../requirements.txt
```

### Issue: Figures not generating

```bash
# Check matplotlib backend
python -c "import matplotlib; print(matplotlib.get_backend())"

# If GUI backend, set to Agg (non-GUI)
export MPLBACKEND=Agg
python create_publication_figures.py
```

### Issue: Low memory

```python
# For very large sequences, process in chunks
# Or increase available memory
# Analysis on 40kb sequence requires ~500MB RAM
```

## File Descriptions

### Analysis Scripts

- **`validation_analysis.py`**: Original class-level validation (PR #86)
- **`validation_analysis_subclass.py`**: **New** subclass-level validation (this PR)
- **`create_figures.py`**: Original figure generation
- **`create_publication_figures.py`**: **New** publication-quality figures with 300 DPI

### Input Data

- **`693fc40d26a53.fasta`**: Test sequence (40,523 bp)
- **`693fc40d26a53_GQ.tsv`**: NBST G-Quadruplex predictions (22 motifs)
- **`693fc40d26a53_Z.tsv`**: NBST Z-DNA predictions (6 motifs)
- **`693fc40d26a53_Slipped_STR.tsv`**: NBST STR predictions (50 motifs)
- **`693fc40d26a53_slipped_DR.tsv`**: NBST Direct Repeat predictions (8 motifs)
- **`693fc40d26a53_MR.tsv`**: NBST Mirror Repeat predictions (9 motifs)
- **`693fc40d26a53_curved.tsv`**: NBST Curved DNA predictions (1 motif)

## Citation

If you use this validation analysis in your research, please cite:

```
NonBDNAFinder Validation Suite (2026). "Comprehensive Subclass-Level Validation 
of NonBDNAFinder Against NBST Benchmark." GitHub Repository: 
https://github.com/RajeshYella/NonBDNAFinder
```

## Contact

For questions or issues with the validation analysis:
- **Open an issue**: GitHub Issues
- **Repository**: https://github.com/RajeshYella/NonBDNAFinder

## Version History

- **v2.0 (Feb 2026)**: Subclass-level validation with 8 publication-quality figures
- **v1.0 (2024)**: Initial class-level validation (PR #86)

---

**Analysis Runtime**: ~2-3 minutes for subclass analysis + ~5 minutes for figure generation  
**Total Output Size**: ~5 MB (figures, tables, reports)  
**Python Version**: 3.8+ required for type hints
