# 🚀 Quick Start Guide - NonBDNAFinder

## Installation

```bash
# Clone repository
git clone https://github.com/VRYella/NonBDNAFinder.git
cd NonBDNAFinder

# Install dependencies
pip install -r requirements.txt

# Verify installation
python -c "import nonbscanner; print('✅ Installation successful!')"
```

## Usage Options

### Option 1: High-Efficiency Jupyter Notebook (Recommended for Large Genomes)

**Best for**: Genomes >100 MB, batch processing, publication-ready outputs

```bash
jupyter notebook HighEfficiency_Genome_Analysis.ipynb
```

**Steps**:
1. **Box 1**: Configure `INPUT_FASTA` path to your genome file
2. **Box 2**: Run analysis (includes progress tracking)
3. **Box 3**: Generate Excel + visualizations

**Output**: 
- Multi-sheet Excel workbook
- CSV file
- 13+ publication-quality plots (300 DPI)

### Option 2: Streamlit Web Application

**Best for**: Interactive analysis, small-medium genomes, visual exploration

```bash
streamlit run app.py
```

**Features**:
- Drag-and-drop file upload
- Real-time progress tracking
- Interactive visualizations
- Multiple export formats

**Navigate to**: http://localhost:8501

### Option 3: Python API

**Best for**: Integration with other tools, custom workflows

```python
from nonbscanner import analyze_sequence
from utilities import export_to_excel, export_to_csv

# Analyze sequence
sequence = "GGGTTAGGGTTAGGGTTAGGG" * 100
motifs = analyze_sequence(sequence, "my_sequence")

# Export results
export_to_excel(motifs, "results.xlsx")
export_to_csv(motifs, "results.csv")

print(f"Found {len(motifs)} motifs")
```

## Input Format

### FASTA File
```
>Sequence_1
ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGG
CCCCCTCCCCCTCCCCCTCCCCATCGATCGCGCGCGCGATCGCACACACACAGCTGC
>Sequence_2
GGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGT
```

### Requirements
- **Format**: FASTA (.fa, .fasta, .fna, .txt)
- **Size**: No limit (automatic chunking for large files)
- **Sequences**: Single or multi-FASTA

## Output Files

### Excel Workbook
- **Consolidated Sheet**: All motifs with core columns
- **Class Sheets**: Separate sheets per motif class (G4, Z-DNA, etc.)
- **Subclass Sheets**: Detailed breakdown by subclass

### CSV File
- Simple table format
- All motifs in one file
- Core columns only

### Visualizations (Notebook Only)
1. Motif distribution by class
2. Nested pie chart (class-subclass hierarchy)
3. Coverage map
4. Density heatmap
5. Length distribution
6. Score distribution
7. Density comparison
8. Circos plot
9. Manhattan plot
10. Cumulative distribution
11. Co-occurrence matrix
12. GC content correlation
13. Length KDE

## Performance Tips

### For Large Genomes (>1 GB)
- Use the Jupyter notebook (not web app)
- Enable parallel processing: `ENABLE_PARALLEL = True`
- Adjust chunk size if needed: `CHUNK_SIZE = 1000000` (1 Mbp)
- Close other applications to free RAM
- Expected time: 1-2 hours for 1 Gbp

### For Multiple Sequences
- Use multi-FASTA format
- Process in batch mode with notebook
- Results generated per sequence

### Memory Optimization
- Default settings handle most genomes
- For very large files (>10 GB):
  - Increase chunk size to 5-10 Mbp
  - Process one chromosome at a time
  - Use systems with 16GB+ RAM

## Troubleshooting

### Issue: "ModuleNotFoundError"
```bash
# Solution: Install dependencies
pip install -r requirements.txt
```

### Issue: "Out of memory"
```python
# Solution: Increase chunk size in notebook Box 1
CHUNK_SIZE = 5000000  # 5 Mbp chunks
```

### Issue: "Analysis too slow"
```python
# Solution: Enable parallel processing
ENABLE_PARALLEL = True
```

### Issue: "Streamlit won't start"
```bash
# Solution: Check if port 8501 is available
streamlit run app.py --server.port 8502
```

## Detected Motif Classes

NonBDNAFinder detects 11 major classes:

1. **Curved DNA**: A-tract mediated curvature
2. **Slipped DNA**: Direct repeats and STRs
3. **Cruciform**: Palindromic inverted repeats
4. **R-Loop**: RNA-DNA hybrid formation sites
5. **Triplex**: Three-stranded structures
6. **G-Quadruplex**: G-rich four-stranded structures (7 variants)
7. **i-Motif**: C-rich structures
8. **Z-DNA**: Left-handed double helix
9. **A-philic DNA**: A/T-rich protein binding sites
10. **Hybrid**: Multi-class overlapping regions
11. **Clusters**: High-density motif hotspots

## Citation

```
NonBDNAFinder: Comprehensive Detection and Analysis of Non-B DNA Motifs
Dr. Venkata Rajesh Yella
GitHub: https://github.com/VRYella/NonBDNAFinder
Email: yvrajesh_bt@kluniversity.in
```

## Support

- **Documentation**: See README.md and CHANGELOG.md
- **Issues**: https://github.com/VRYella/NonBDNAFinder/issues
- **Email**: yvrajesh_bt@kluniversity.in

---

**Version**: 2024.12.12
**Last Updated**: December 12, 2024
