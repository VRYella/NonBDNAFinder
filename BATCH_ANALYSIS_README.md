# Batch Genome Analysis and Publication-Ready Results

Complete batch processing pipeline for analyzing multiple genomes and generating publication-ready outputs.

## Overview

This batch analysis system processes all genome files in the `Genomes/` folder and generates:
- Comprehensive motif detection results for each genome
- Cross-genome comparative analysis
- Publication-ready Excel workbook with multiple analysis sheets
- Narrative text sections for manuscripts (Materials & Methods, Results, Discussion)
- High-quality visualization figures (300 DPI, publication-ready)

## Quick Start

### 1. Run Complete Pipeline

```bash
# Run entire pipeline (recommended for first-time users)
python run_batch_pipeline.py
```

This single command runs all steps automatically:
1. Batch analysis of all genomes
2. Comparative analysis
3. Excel export
4. Narrative generation
5. Visualization generation

### 2. Run Individual Steps

You can also run each step separately:

```bash
# Step 1: Batch analysis
python batch_analysis.py --genomes-dir Genomes --output-dir batch_results

# Step 2: Comparative analysis
python comparative_analysis.py --input-dir batch_results --output-dir comparative_results

# Step 3: Publication Excel export
python publication_export.py --batch-dir batch_results --comparative-dir comparative_results --output publication_results.xlsx

# Step 4: Publication narrative
python publication_narrative.py --comparative-dir comparative_results --output publication_narrative.txt

# Step 5: Publication visualizations
python publication_visualizations.py --comparative-dir comparative_results --output-dir publication_figures
```

## Input Requirements

### Genome Files
- Place all genome files in the `Genomes/` directory
- Supported formats: `.fna`, `.fasta`, `.fa`
- Can be single-sequence or multi-sequence FASTA files

### Example Structure
```
NonBDNAFinder/
├── Genomes/
│   ├── ecoli.fna
│   ├── saureus.fna
│   ├── mtb.fasta
│   └── ...
├── batch_analysis.py
├── comparative_analysis.py
├── publication_export.py
└── run_batch_pipeline.py
```

## Output Structure

After running the pipeline, the following outputs are generated:

```
NonBDNAFinder/
├── batch_results/
│   ├── batch_summary.json          # Summary of all genomes
│   ├── ecoli_motifs.json           # Detailed motifs for E. coli
│   ├── saureus_motifs.json         # Detailed motifs for S. aureus
│   └── ...                         # One file per genome
│
├── comparative_results/
│   ├── comparative_analysis.json   # Complete comparative analysis
│   └── comparative_report.txt      # Human-readable report
│
├── publication_results.xlsx        # Unified Excel workbook
├── publication_narrative.txt       # Manuscript sections
│
└── publication_figures/
    ├── density_comparison.png
    ├── class_distribution_heatmap.png
    ├── diversity_vs_density.png
    ├── enrichment_patterns.png
    ├── genome_statistics_radar.png
    └── correlation_matrix.png
```

## Output Details

### 1. Batch Results (`batch_results/`)
- **batch_summary.json**: Overview of all analyzed genomes with sequence info and motif counts
- **[genome]_motifs.json**: Complete motif data for each genome (positions, classes, scores)

### 2. Comparative Results (`comparative_results/`)
- **comparative_analysis.json**: Full comparative analysis including:
  - Per-genome statistics (density, coverage, diversity)
  - Cross-genome class frequencies
  - Enrichment/depletion patterns
  - Feature correlations
  - Summary statistics
- **comparative_report.txt**: Human-readable text report

### 3. Publication Excel Workbook (`publication_results.xlsx`)
Multiple sheets in a single workbook:
- **Summary**: Overview table with all genomes
- **[Genome]**: Detailed sheet for each genome with class distribution
- **Comparative Analysis**: Cross-genome frequency comparison
- **Statistical Analysis**: Enrichment patterns and correlations

### 4. Publication Narrative (`publication_narrative.txt`)
Four manuscript-ready sections:
- **Materials and Methods**: Computational approach description
- **Results - Distribution Analysis**: Global and genome-specific patterns
- **Results - Comparative Genomics**: Enrichment and correlation analysis
- **Discussion and Conclusions**: Biological interpretation and future directions

### 5. Publication Figures (`publication_figures/`)
High-quality visualizations (300 DPI, PNG format):
- **density_comparison.png**: Bar plot of motif density across genomes
- **class_distribution_heatmap.png**: Heatmap showing class frequencies
- **diversity_vs_density.png**: Scatter plot of diversity vs density
- **enrichment_patterns.png**: Bar plot with enrichment/depletion patterns
- **genome_statistics_radar.png**: Multi-feature comparison radar plot
- **correlation_matrix.png**: Feature correlation heatmap

## Advanced Usage

### Custom Output Directories

```bash
python run_batch_pipeline.py --genomes-dir /path/to/genomes --output-base /path/to/output
```

### Skip Steps (use existing results)

```bash
# Skip batch analysis if already completed
python run_batch_pipeline.py --skip-batch

# Skip comparative analysis if already completed
python run_batch_pipeline.py --skip-comparative
```

### Re-run Specific Steps

```bash
# Re-generate only visualizations
python publication_visualizations.py

# Re-generate only narrative
python publication_narrative.py

# Re-generate only Excel
python publication_export.py
```

## Interpretation Guide

### Motif Density
- **Units**: Motifs per kilobase (motifs/kb)
- **Interpretation**: Higher density suggests more non-B DNA forming regions
- **Typical range**: 1-20 motifs/kb (varies by organism and genome size)

### Genomic Coverage
- **Units**: Percentage (%)
- **Interpretation**: Fraction of genome occupied by non-B DNA motifs
- **Typical range**: 0.1-5% (most genomes < 2%)

### Class Diversity (Shannon Entropy)
- **Units**: H (bits)
- **Interpretation**: Measure of structural heterogeneity
- **Range**: 0 (single class) to log₂(11) ≈ 3.46 (perfectly uniform across 11 classes)
- **Higher values**: More diverse non-B DNA landscape

### Enrichment Z-score
- **Interpretation**: 
  - Z > 2.0: Significantly enriched
  - Z < -2.0: Significantly depleted
  - |Z| < 2.0: No significant enrichment/depletion

### Coefficient of Variation (CV)
- **Interpretation**: Relative variability across genomes
- **CV > 1.0**: High variability
- **CV < 0.5**: Low variability (consistent across genomes)

## Performance Expectations

### Processing Time
- **Small genome (<1 Mb)**: 1-5 minutes
- **Medium genome (1-5 Mb)**: 5-20 minutes
- **Large genome (>5 Mb)**: 20-60 minutes

### Batch Pipeline
- **9 genomes** (as in Genomes/ folder): ~30-90 minutes total
- Progress is displayed for each genome
- Chunked processing handles large genomes efficiently

### Memory Requirements
- **Typical**: 2-4 GB RAM
- **Large genomes**: Up to 8 GB RAM
- Chunked processing keeps memory usage reasonable

## Troubleshooting

### "No genome files found"
- Check that genome files are in the `Genomes/` directory
- Verify files have `.fna`, `.fasta`, or `.fa` extensions
- Use `--genomes-dir` to specify a different directory

### "Batch summary not found"
- Run `batch_analysis.py` first before running other scripts
- Or use `--skip-batch` only after batch analysis is complete

### "Out of memory" errors
- Process genomes individually if needed
- Close other applications to free memory
- Consider chunked processing for very large genomes

### Python module errors
- Install required dependencies: `pip install -r requirements.txt`
- Ensure NumPy, Pandas, Matplotlib, Seaborn, SciPy, and OpenPyXL are installed

## Citation

If you use this batch analysis pipeline in your research, please cite:

```
NonBDNAFinder: Comprehensive Detection and Analysis of Non-B DNA Motifs
Dr. Venkata Rajesh Yella
GitHub: https://github.com/VRYella/NonBDNAFinder
Email: yvrajesh_bt@kluniversity.in
```

## Support

For questions, issues, or feature requests:
- Email: yvrajesh_bt@kluniversity.in
- GitHub Issues: https://github.com/VRYella/NonBDNAFinder/issues

## License

MIT License - see LICENSE file for details

---

**Version**: 2025.1  
**Last Updated**: January 2026  
**Status**: Production-Ready
