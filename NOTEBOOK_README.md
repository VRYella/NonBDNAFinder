# 🧬 NonBScanner High-Efficiency Genome Analysis Notebook

**Quick Start Guide for Big Genome Analysis**

---

## 📖 Overview

This Jupyter notebook (`HighEfficiency_Genome_Analysis.ipynb`) provides a streamlined, high-performance workflow for analyzing large genomic sequences with Non-B DNA motif detection.

### Key Features

- ⚡ **3 Simple Execution Boxes**: Setup → Analyze → Results
- 🚀 **Parallel Processing**: Multi-core execution for maximum speed
- 💾 **Memory Efficient**: Handles genomes 100MB+ with chunked processing
- 📊 **Consolidated Excel Output**: All motif classes in separate sheets
- 📈 **Publication-Quality Visualizations**: 7+ chart types at 300 DPI
- 🔍 **Complete Detection**: All 11 Non-B DNA motif classes

---

## 🚀 Quick Start

### 1. Install Dependencies

```bash
# Install required packages
pip install -r requirements.txt

# Or install individually
pip install numpy pandas matplotlib seaborn plotly biopython openpyxl xlsxwriter psutil
```

### 2. Launch Jupyter Notebook

```bash
# Start Jupyter Notebook
jupyter notebook HighEfficiency_Genome_Analysis.ipynb

# Or use JupyterLab
jupyter lab HighEfficiency_Genome_Analysis.ipynb
```

### 3. Configure Your Analysis

In **Box 1** (Setup), update the following parameter:

```python
# Change this to your genome FASTA file path
GENOME_FILE = "your_genome.fasta"
```

### 4. Execute the Boxes

Run each box in sequence:

- **Box 1**: Setup and configuration
- **Box 2**: Run high-efficiency analysis
- **Box 3**: Generate consolidated Excel and visualizations

---

## 📁 Input Requirements

### FASTA File Format

Your genome file should be in standard FASTA format:

```
>sequence_name_1 Optional description
ATCGATCGATCG...
>sequence_name_2 Optional description
GCTAGCTAGCTA...
```

### Example Files

A test file is included: `test_genome_example.fasta`

To use it, set in Box 1:
```python
GENOME_FILE = "test_genome_example.fasta"
```

---

## 📊 Output Files

After running all three boxes, you'll find the following outputs in the `analysis_results/` directory:

### 1. Excel File (Consolidated Results)
- **File**: `NonBScanner_Results_YYYYMMDD_HHMMSS.xlsx`
- **Sheets**:
  - `Summary`: Overall statistics
  - `All_Motifs`: Complete motif list
  - Individual sheets for each motif class (Curved_DNA, G-Quadruplex, etc.)

### 2. CSV File (Alternative Format)
- **File**: `NonBScanner_Results_YYYYMMDD_HHMMSS.csv`
- Tab-delimited format compatible with Excel and other tools

### 3. BED File (Genome Browser Format)
- **File**: `NonBScanner_Results_YYYYMMDD_HHMMSS.bed`
- Compatible with UCSC Genome Browser and IGV

### 4. Detection Report
- **File**: `Detection_Report_YYYYMMDD_HHMMSS.txt`
- Text report showing which motif classes/subclasses were detected

### 5. Visualizations
- **Directory**: `figures_YYYYMMDD_HHMMSS/`
- All plots saved as PNG files at 300 DPI
- Chart types:
  1. Comprehensive class analysis
  2. Comprehensive subclass analysis
  3. Score statistics by class
  4. Length statistics by class
  5. Motif distribution
  6. Coverage map
  7. Genome landscape track (for sequences < 50kb)

---

## ⚙️ Configuration Options

### Performance Parameters (Box 1)

```python
CHUNK_SIZE = 50000      # 50kb chunks (optimal for most systems)
OVERLAP_SIZE = 500      # 500bp overlap to catch boundary motifs
USE_PARALLEL = True     # Enable parallel processing (recommended)
```

### Visualization Parameters

```python
PLOT_DPI = 300          # Publication quality (300 DPI)
SAVE_PLOTS = True       # Save plots to files
```

### Adjust for Your System

- **Large RAM (16GB+)**: Increase `CHUNK_SIZE` to 100000 for faster processing
- **Limited RAM (4GB)**: Decrease `CHUNK_SIZE` to 25000
- **Single-core CPU**: Set `USE_PARALLEL = False`
- **Many cores (8+)**: Default settings will utilize all cores

---

## 🔬 Detected Motif Classes

The notebook detects **11 major classes** with **22+ specialized subclasses**:

1. **Curved DNA** - A-tract mediated bending
2. **Slipped DNA** - Direct repeats, STRs
3. **Cruciform** - Inverted repeats
4. **R-Loop** - RNA-DNA hybrids (QmRLFS algorithm)
5. **Triplex** - Three-stranded structures
6. **G-Quadruplex** - 7 subclasses (Canonical, Bulged, Relaxed, etc.)
7. **i-Motif** - C-rich structures
8. **Z-DNA** - Left-handed helix
9. **A-philic DNA** - A-rich structures
10. **Hybrid** - Multi-class overlaps
11. **Clusters** - High-density regions

---

## 📈 Performance Expectations

### Typical Throughput

- **Small sequences (<10kb)**: ~5,000-8,000 bp/s
- **Medium sequences (10-100kb)**: ~15,000-25,000 bp/s (with parallel processing)
- **Large sequences (>100kb)**: ~20,000-30,000 bp/s (with chunking)

### Example Analysis Times

| Genome Size | Estimated Time | Memory Usage |
|-------------|----------------|--------------|
| 10 kb       | < 1 second     | ~10 MB       |
| 100 kb      | ~4-5 seconds   | ~50 MB       |
| 1 MB        | ~40-50 seconds | ~200 MB      |
| 10 MB       | ~7-8 minutes   | ~1 GB        |
| 100 MB      | ~60-90 minutes | ~5 GB        |

*Times are approximate and vary based on motif complexity and system specifications.*

---

## 🐛 Troubleshooting

### Import Errors

If you see `ModuleNotFoundError`:

```bash
# Ensure all dependencies are installed
pip install -r requirements.txt

# Or install specific missing package
pip install <package_name>
```

### Memory Errors

If you run out of memory:

1. Reduce `CHUNK_SIZE` in Box 1:
   ```python
   CHUNK_SIZE = 25000  # Smaller chunks use less memory
   ```

2. Process sequences individually (modify Box 2 to process one at a time)

3. Close other applications to free up RAM

### File Not Found

Ensure your FASTA file path is correct in Box 1:

```python
# Use absolute path if needed
GENOME_FILE = "/full/path/to/your/genome.fasta"

# Or relative path from notebook location
GENOME_FILE = "../data/genome.fasta"
```

### Slow Performance

To improve speed:

1. Ensure `USE_PARALLEL = True` in Box 1
2. Increase `CHUNK_SIZE` if you have sufficient RAM
3. Close other applications to free up CPU cores
4. Consider using a machine with more CPU cores

---

## 📚 Additional Resources

- **Main README**: [README.md](README.md) - Full NonBScanner documentation
- **Comprehensive Documentation**: [COMPREHENSIVE_DOCUMENTATION.md](COMPREHENSIVE_DOCUMENTATION.md)
- **Performance Guide**: [docs/perf_runbook.md](docs/perf_runbook.md)
- **Visualization Guide**: [VISUALIZATION_GUIDE.md](VISUALIZATION_GUIDE.md)

---

## 🆘 Support

For issues or questions:

- **GitHub Issues**: https://github.com/VRYella/NonBScanner/issues
- **Email**: yvrajesh_bt@kluniversity.in

---

## 📜 Citation

If you use NonBScanner in your research, please cite:

```
NonBScanner: Comprehensive Detection and Analysis of Non-B DNA Motifs
Dr. Venkata Rajesh Yella
GitHub: https://github.com/VRYella/NonBScanner
```

---

**Happy Genome Analyzing! 🧬**
