# NonBDNA Finder Standalone - Quick Start Guide

## 🚀 What's Been Created

A complete standalone version of NonBDNAFinder has been created in the `standalone/` folder with:

### 📓 Jupyter Notebook Interface (`NonBDNA_Finder_Notebook.ipynb`)
- **File Upload Widget**: Drag-and-drop FASTA file upload
- **Real-time Analysis**: Progress feedback during motif detection
- **Interactive Visualizations**: Plotly charts for results exploration
- **Download Package**: One-click download of complete results as ZIP

### 🐍 Command-Line Interface (`launcher.py`)
- **Batch Processing**: Process multiple files programmatically
- **Scriptable**: Integration with bioinformatics pipelines
- **Detailed Output**: Multiple file formats (CSV, Excel, GFF3, summary)

### 🧬 Core Engine (`nbdfinder_core.py`)
- **11 Motif Classes**: G4, i-motif, Z-DNA, Triplex, Cruciform, Curved, Slipped, A-philic, R-loop, Hybrid, Clusters
- **Multiple Algorithms**: G4Hunter, Z-DNA scoring, Triplex stability, Pattern matching
- **Optimized Performance**: Optional Hyperscan acceleration
- **Self-contained**: Minimal dependencies, works independently

## 📦 Installation & Usage

### Quick Start (3 steps):
```bash
# 1. Navigate to standalone folder
cd standalone/

# 2. Install dependencies
pip install -r requirements.txt

# 3. Launch notebook
jupyter notebook NonBDNA_Finder_Notebook.ipynb
```

### Alternative: Command Line
```bash
# Analyze a FASTA file
python launcher.py your_sequences.fasta

# Get detailed statistics
python launcher.py your_sequences.fasta --show-stats

# Custom output name
python launcher.py your_sequences.fasta --output my_analysis
```

## 📊 Output Package

Each analysis generates a ZIP file containing:
- **📄 CSV file** - For Excel, R, Python analysis
- **📊 Excel workbook** - Multi-sheet with class-specific data
- **🧬 GFF3 file** - For UCSC, IGV, JBrowse genome browsers
- **📋 Summary report** - Human-readable analysis overview

## ✅ Verification

Run these to verify everything works:
```bash
# Quick installation check
python verify_installation.py

# Full test suite
python test_standalone.py

# Test command-line interface
echo ">test_seq\nGGGTTTTGGGTTTTGGGTTTTGGG" > test.fasta
python launcher.py test.fasta
```

## 🎯 Requirements Met

✅ **Standalone version**: Complete self-contained package  
✅ **Upload option only**: Jupyter notebook with file upload widget  
✅ **Downloadable zip**: All results packaged as ZIP file  
✅ **Jupyter notebook**: Interactive interface with visualizations  
✅ **Separate folder**: Located in `standalone/` directory  

The standalone version provides all the core functionality of NonBDNAFinder in a user-friendly package suitable for researchers who want a simple, self-contained solution for non-B DNA motif detection.