# NonBDNA Finder - Standalone Version

A self-contained Jupyter notebook interface for detecting non-B DNA structures in genomic sequences.

## 🚀 Quick Start

### 1. Installation

```bash
# Clone or download the standalone folder
# Navigate to the standalone directory
cd standalone

# Install dependencies
pip install -r requirements.txt

# Optional: Install Jupyter if not already available
pip install jupyter
```

### 2. Launch the Notebook

```bash
jupyter notebook NonBDNA_Finder_Notebook.ipynb
```

### 3. Use the Interface

1. **Upload FASTA File**: Use the file upload widget in the notebook
2. **Analyze**: Click the "Analyze Sequences" button
3. **Download Results**: Get a zip package with all output formats

## 📁 What's Included

### Core Files
- `NonBDNA_Finder_Notebook.ipynb` - Main Jupyter notebook interface
- `nbdfinder_core.py` - Standalone analysis engine
- `requirements.txt` - Python dependencies

### Output Package
When you run an analysis, you'll get a zip file containing:
- **CSV file** - Spreadsheet format for data analysis
- **Excel file** - Multi-sheet workbook with class-specific data
- **GFF3 file** - Genome browser compatible format
- **Summary report** - Human-readable analysis summary

## 🧬 Detected Motif Classes

This standalone version detects 11 classes of non-B DNA structures:

1. **G-Quadruplex Family** - G-rich four-stranded structures
2. **i-motif family** - C-rich quadruplex structures  
3. **Z-DNA** - Left-handed helical structures
4. **Triplex** - Three-stranded DNA structures
5. **Cruciform DNA** - Four-way junction structures
6. **Curved DNA** - Bent DNA structures
7. **Slipped DNA** - Repetitive sequence structures
8. **A-philic DNA** - A-tract containing regions
9. **R-loop** - RNA-DNA hybrid structures
10. **Hybrid** - Overlapping motif structures
11. **Non-B DNA cluster regions** - Hotspot regions

## ⚡ Features

### Core Functionality
- **File Upload Interface** - Easy drag-and-drop FASTA upload
- **Real-time Analysis** - Process sequences with progress feedback
- **Interactive Visualizations** - Charts and plots for result exploration
- **Multi-format Output** - CSV, Excel, GFF3, and summary reports
- **Downloadable Results** - Packaged as convenient zip files

### Analysis Features
- **Pattern Recognition** - Regex-based motif detection
- **Scoring Algorithms** - Class-specific scoring methods
- **Statistical Analysis** - Summary statistics and distributions
- **Position Mapping** - Motif location visualization
- **Clustering Analysis** - Identify motif hotspots

### Performance
- **Optimized Patterns** - Efficient regex compilation
- **Hyperscan Support** - Optional acceleration (if available)
- **Memory Efficient** - Handles large sequences
- **Self-contained** - No external dependencies beyond Python packages

## 📋 Requirements

### Python Version
- Python 3.7 or higher

### Core Dependencies
- `jupyter` - Notebook interface
- `pandas` - Data manipulation
- `numpy` - Numerical computing
- `biopython` - Sequence analysis
- `matplotlib` - Basic plotting
- `openpyxl` - Excel file support
- `plotly` - Interactive visualizations
- `ipywidgets` - Notebook widgets

### Optional Dependencies
- `hyperscan` - Optimized pattern matching (Linux/Mac only)

## 🔧 Usage Examples

### Basic Usage
1. Open the notebook
2. Run the setup cell to install dependencies
3. Upload your FASTA file using the widget
4. Click "Analyze Sequences"
5. Download the results zip file

### Example Data
The notebook includes sample sequences for testing:
```
>example_sequence_with_multiple_motifs
GGGTTTTGGGTTTTGGGTTTTGGGAAACCCAAACCCAAACCCAAACCC...
```

### Advanced Analysis
The notebook also provides optional advanced features:
- Position mapping along sequences
- Clustering analysis with adjustable windows
- Custom visualization options

## 📊 Output Format

### CSV/Excel Columns
- `S.No` - Serial number
- `Sequence_Name` - Input sequence identifier
- `Chromosome/Contig` - Sequence name
- `Class` - Non-B DNA motif class
- `Subclass` - Motif subtype
- `Start` - Start position (1-based)
- `End` - End position (1-based)
- `Length` - Motif length in base pairs
- `Normalized_Score` - Score normalized to 0-1 range
- `Actual_Score` - Raw algorithm score
- `Scoring_Method` - Algorithm used for scoring
- `GC_Content` - GC percentage of motif sequence
- `Sequence` - Actual DNA sequence of motif

### GFF3 Format
Standard GFF3 format compatible with:
- UCSC Genome Browser
- IGV (Integrative Genomics Viewer)
- JBrowse
- Other genome visualization tools

## 🧪 Testing

### Test with Example Data
```python
# In the notebook, use the "Analyze Example" button
# This will run analysis on built-in test sequences
```

### Validate Installation
```python
# Run this in a notebook cell to check dependencies
import pandas as pd
import numpy as np
from Bio import SeqIO
import plotly.express as px
from nbdfinder_core import StandaloneNonBDNAFinder
print("✅ All dependencies loaded successfully!")
```

## 🔬 Scientific Background

### Algorithms Used
- **G4Hunter** - G-quadruplex detection
- **Z-DNA Scoring** - Dinucleotide transition analysis
- **Triplex Scoring** - Purine/pyrimidine tract detection
- **Pattern Matching** - Regex-based motif recognition

### Validation
This standalone version implements the same core algorithms as the full NonBDNAFinder suite, validated against:
- Known non-B DNA structures
- Literature-reported motifs
- Experimental structural data

## 🤝 Support

### Troubleshooting
- **File Upload Issues**: Ensure FASTA format with '>' headers
- **Analysis Errors**: Check sequence quality and format
- **Memory Issues**: Try analyzing smaller sequences
- **Visualization Problems**: Update plotly and jupyter

### Performance Tips
- Install hyperscan for faster analysis (Linux/Mac)
- Use smaller chunk sizes for very large genomes
- Close other applications to free memory
- Consider running on high-memory systems for genome-scale analysis

## 📝 Citation

If you use this standalone version in your research, please cite:

> NonBDNAFinder: High-performance detection of non-B DNA structures
> [Authors and publication details when available]

## 📄 License

This standalone version is released under the same license as the main NonBDNAFinder project.

## 🔄 Version History

- **v1.0** - Initial standalone release
  - Jupyter notebook interface
  - Core motif detection algorithms
  - Multi-format output support
  - Interactive visualizations

---

*For the full-featured version with additional algorithms and performance optimizations, see the main NonBDNAFinder repository.*