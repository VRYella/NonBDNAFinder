# 🧬 NonBScanner - Professional Non-B DNA Motif Detection Suite

**Comprehensive, elegant, and high-performance detection of Non-B DNA structures with 11 major classes and 22+ specialized subclasses**

## 🎯 Overview

NonBScanner is a state-of-the-art bioinformatics tool for detecting and analyzing Non-B DNA motifs in genomic sequences. It combines high-performance optimized algorithms with scientific scoring methods to provide comprehensive analysis of structural DNA elements.

**New in Version 2024.1:** Professional 5-file architecture for maximum elegance and maintainability.

**✨ Latest Enhancements (2024.1.2):**
- **📏 Density Analysis**: Genomic density (coverage %) and positional density (motifs/kbp) calculations per motif class
- **✨ Enrichment Analysis**: Fold enrichment calculation with 100-iteration sequence shuffling for statistical validation
- **📊 Statistical Significance**: P-value calculation using permutation testing to validate motif enrichment
- **🎨 Enhanced Visualizations**: New density comparison and enrichment analysis plots with no text overlap
- **Consolidated Registry**: Single file (`consolidated_registry.json`) instead of 18 separate files
- **Enhanced Scientific Visualizations**: Publication-quality plotting functions
- **Comprehensive Class Analysis**: Shows all 11 classes with detection status
- **Advanced Statistics**: Mean, median, std, min/max for scores and lengths
- **300 DPI Output**: Publication-ready figures

### 📊 Detection Coverage
- **11 Major Non-B DNA Classes** with comprehensive subclass analysis  
- **22+ Specialized Subclasses** for detailed motif characterization
- **High-performance detection** (24,674 bp/s on 100K sequences)
- **Literature-validated** scoring algorithms
- **✨ Enhanced hybrid/cluster detection** with actual sequence extraction
- **✨ 25+ advanced publication-quality visualizations** (colorblind-friendly)
- **📏 Rigorous statistical analysis** with density and enrichment metrics

### ⚡ Performance Highlights
- **100,000 bp in 4 seconds** with optimized detectors
- **24,674 bp/second** processing speed
- **Memory efficient**: ~5 MB for 100K sequences
- **Production ready**: Tested on large genomic datasets

## 🏗️ Architecture - Professional 5-File Design

NonBScanner features an elegant, professional architecture with just **5 core Python files**:

```
NonBScanner/
├── nonbscanner.py      # Main API & Scanner Orchestration (~600 lines)
├── detectors.py        # All 9 Motif Detector Classes (~3,500 lines)
├── utilities.py        # Sequence I/O, Export & Statistics (~2,100 lines)
├── visualizations.py   # Complete Visualization Suite (~1,000 lines)
└── app.py              # Streamlit Web Interface (~1,800 lines)
```

**Key Features:**
- ✅ **Minimal & Clean**: Just 5 files for the entire tool
- ✅ **Professional**: Well-documented, type-annotated code
- ✅ **Modular**: Clear separation of concerns
- ✅ **Elegant**: Advanced algorithms in readable structure
- ✅ **Production-Ready**: Tested on large genomic datasets

### Supporting Files
- `scanner.py` - Low-level k-mer indexing functions (used by detectors)
- `consolidated_registry.json` - Single file with all 411 pattern definitions
- `example_motifs_multiline.fasta` - Example FASTA file with all motif types

## 🔬 Supported Motif Classes

| Class | Name | Subclasses | Key Features |
|-------|------|------------|--------------|
| **1** | Curved DNA | Global Curvature, Local Curvature | A-tract mediated curvature |
| **2** | Slipped DNA | Direct Repeat, STR | Tandem repeats, slipped structures |
| **3** | Cruciform | Palindromic Inverted Repeat | Four-way junctions |
| **4** | R-Loop | R-Loop Formation Site, QmRLFS-m1, QmRLFS-m2 | RNA-DNA hybrids, QmRLFS algorithm |
| **5** | Triplex | Mirror Repeat, Sticky DNA | Three-stranded structures |
| **6** | G-Quadruplex | Canonical, Bulged, Relaxed, Bipartite, Multimeric, Imperfect, G-Triplex | Four-stranded G-rich structures |
| **7** | i-Motif | Canonical, Extended, AC-Motif | C-rich structures |
| **8** | Z-DNA | Classic Z-DNA, eGZ | Left-handed double helix |
| **9** | Hybrid | Multi-class Overlap, Composite | Overlapping motifs |
| **10** | Cluster | Motif Hotspot, Mixed Cluster | High-density regions |

## 📚 Documentation

This README provides comprehensive documentation for using NonBScanner. The tool supports multiple interfaces:
- Python API (see Quick Start below)
- Web Interface (Streamlit app)
- Jupyter Notebook (`NonBScanner_Local.ipynb`)
- Shell Script (`generate_csv_output.sh`)

## 🚀 Quick Start

### Python API (Recommended)
```python
import nonbscanner as nbs
from utilities import analyze_class_subclass_detection, print_detection_report

# Analyze a single sequence
sequence = "GGGTTAGGGTTAGGGTTAGGG"
motifs = nbs.analyze_sequence(sequence, "my_sequence")

print(f"Found {len(motifs)} motifs:")
for motif in motifs:
    print(f"  {motif['Class']} at position {motif['Start']}-{motif['End']}")

# Analyze FASTA file
results = nbs.analyze_file("sequences.fasta")
for name, motifs in results.items():
    print(f"{name}: {len(motifs)} motifs detected")

# Combine all motifs for comprehensive analysis
all_motifs = []
for name, motifs in results.items():
    all_motifs.extend(motifs)

# Export results to Excel with separate sheets for each motif class
nbs.export_results(all_motifs, format='excel', filename='output.xlsx')

# Export to CSV (traditional format)
nbs.export_results(all_motifs, format='csv', filename='output.csv')

# Analyze class/subclass detection status
detection_report = analyze_class_subclass_detection(all_motifs)
report_text = print_detection_report(detection_report)
print(report_text)
```

### Web Interface
```bash
# Clone and setup
git clone https://github.com/VRYella/NonBScanner.git
cd NonBScanner
pip install -r requirements.txt

# Launch web interface
streamlit run app.py                    # Web interface on :8501
```

### Jupyter Notebook (Recommended for Comprehensive Analysis)
```bash
# Launch Jupyter notebook - includes Excel export and detection analysis
jupyter notebook NonBScanner_Local.ipynb

# Or with JupyterLab
jupyter lab NonBScanner_Local.ipynb
```

The notebook provides:
- **Excel output** with separate sheets for each motif class/subclass
- **Comprehensive class/subclass detection analysis** showing which motifs were not predicted
- **Publication-quality visualizations**
- **Step-by-step workflow** from FASTA input to final reports


## 📱 User Interfaces

### **Streamlit Web Application** (`http://localhost:8501`)
- Interactive sequence upload and analysis
- Comprehensive visualization suite (21+ chart types)
- Real-time analysis with progress tracking
- Export capabilities (CSV, BED, BigWig)
- Documentation and tutorial sections

### **Jupyter Notebook** (Local Interactive Analysis)
- **File**: `NonBScanner_Local.ipynb`
- Interactive cell-by-cell analysis
- Built-in visualizations and examples
- Perfect for exploratory analysis

### **Shell Script** (Batch CSV Generation)
- **File**: `generate_csv_output.sh`
- Command-line batch processing
- Generate CSV for all motif classes
- Automated summary statistics
- Run with `--help` flag for usage information

## 🛠️ Technical Features

### Performance
- **Optimized Algorithms**: Linear O(n) complexity for all major detectors
- **Seed-and-Extend K-mer Index**: Genome-scale efficient repeat detection
- **No Sequence Limits**: Handle sequences of any size without performance degradation
- **~280,000 bp/second**: Slipped DNA detector on 50kb sequences
- **~5,800 bp/second**: Overall rate including all detector types on 10kb sequences
- **Memory Efficient**: K-mer indexing with frequency limits for safety
- **Production Ready**: Rigorously tested on real genome sequences

### Architecture Highlights
- **Hybrid Detection Approach**:
  - **Optimized Python Scanner** (seed-and-extend k-mer index) for Slipped DNA, Cruciform, Triplex
  - **Hyperscan** (regex/pattern matching) for Z-DNA, G4, i-Motif, Curved DNA, A-Philic
  - **Algorithmic** (QmRLFS, GC-skew) for R-Loop detection
- See `OPTIMIZED_SCANNER_ARCHITECTURE.md` for detailed architecture documentation

### Performance Notes
- **Slipped DNA**: No size limits, O(n) complexity (was O(n²) with 50kb limit)
- **Cruciform**: No size limits, O(n) complexity (was O(n²) with 1kb limit)
- **Triplex**: No size limits, O(n) complexity with purine/pyrimidine filtering
- **All detectors**: Linear scaling validated on sequences up to 50kb+
- See `OPTIMIZED_SCANNER_ARCHITECTURE.md` for benchmarks

### Scientific Accuracy
- **Literature-Based**: Algorithms from peer-reviewed research
- **QmRLFS Integration**: Advanced R-loop detection with RIZ/REZ analysis
- **Validated Thresholds**: Biologically relevant cut-offs
- **Raw Score Reporting**: Direct algorithm outputs without normalization
- **Overlap Resolution**: Automatic removal of overlapping motifs within subclasses
- **Quality Control**: Built-in validation and error checking

### Export & Integration
- **Multiple Formats**: BED, BigWig, CSV, JSON export
- **Genome Browser**: UCSC/IGV compatible outputs
- **Batch Processing**: Multi-FASTA support

## 🔬 Scoring Algorithms

- **G4Hunter**: G-quadruplex prediction (Bedrat et al., 2016)
- **QmRLFS**: R-loop formation site detection with RIZ/REZ analysis (Jenjaroenpun & Wongsurawat, 2016)
- **Z-Seeker**: Z-DNA detection (Ho et al., 1986)
- **Enhanced Curvature**: A-tract curvature (Olson et al., 1998)
- **Instability Scoring**: Repeat instability analysis
- **Thermodynamic Models**: Cruciform stability prediction
- **Triplex Potential**: Three-strand formation scoring

## 📈 Visualization Suite

### Static Plots (Classic)
- Motif distribution analysis
- Coverage and density maps
- Length distribution analysis
- Sequence composition analysis
- Class/subclass comparisons

### Enhanced Scientific Visualizations (NEW! 🎨✨)
**Publication-quality static plots with comprehensive statistics:**

**Core Analysis Functions:**
1. **Comprehensive Class Analysis** - Shows all 11 Non-B DNA classes with detection status
   - Distribution bar chart with color coding
   - Detection status pie chart
   - Statistics table (count, avg length, avg score)
   - List of non-detected classes
   
2. **Comprehensive Subclass Analysis** - Detailed subclass breakdown organized by parent class
   - Horizontal bar chart of all subclasses
   - Color-coded by parent class
   - Summary statistics by class
   
3. **Score Statistics by Class** - Advanced statistical visualization
   - Violin plots showing distributions
   - Box plot overlays with quartiles
   - Statistical annotations (μ, σ)
   - Comprehensive stats table (mean, median, std, min, max)
   
4. **Length Statistics by Class** - Distribution analysis
   - Overlaid histograms for each class
   - Box plot comparison
   - Statistical summary table

**Previous Advanced Visualizations:**
5. **Genome Landscape Track** - Horizontal ruler with colored glyphs showing motif positions
6. **Sliding Window Heat Ribbon** - 1D heatmap with density and score overlay
7. **Ridge Plots (Joyplots)** - Stacked density ridges for length distributions
8. **Sunburst/Treemap** - Hierarchical composition visualization
9. **Hexbin with Marginals** - 2D density plot with marginal histograms
10. **UpSet Plot** - Clear intersection visualization (better than Venn diagrams)
11. **Violin + Beeswarm** - Score distributions with individual data points
12. **Cluster Hotspot Map** - Regional cluster analysis with annotations

**Design Features:**
- ✅ Export as PNG @300 DPI for publications
- ✅ Colorblind-safe palette
- ✅ Clean sans-serif typography
- ✅ Statistical annotations (mean, median, std)
- ✅ Shows ALL classes (detected and not detected)
- ✅ Publication-ready layouts

### Interactive Visualizations
- Motif browser with zoom/pan
- Class hierarchy sunburst charts
- Position-based track plots
- Statistical correlation plots
- Network analysis graphs

## 🧪 Example Analysis

Use the web interface at `http://localhost:8501` to:
- Upload FASTA sequences or paste sequence data
- Analyze G-quadruplex, Z-DNA, R-loops, and other motifs
- Visualize results with comprehensive charts
- Export findings in BED, CSV, or JSON formats

## 🔗 Hybrid and Cluster Motif Separation (NEW)

NBDScanner now separates hybrid and cluster motifs from regular Non-B DNA motifs for cleaner analysis:

### What are Hybrid and Cluster Motifs?

- **Hybrid Motifs**: Regions where different Non-B DNA classes overlap (30-70% overlap)
  - Example: `R-Loop_Cruciform_Overlap` - indicates complex genomic regions
  - Shows interaction between different structural elements

- **Cluster Motifs**: High-density regions containing multiple Non-B DNA motifs from different classes
  - Example: `Mixed_Cluster_10_classes` - hotspots of Non-B DNA activity
  - Indicates regions with exceptional structural diversity

### How It Works

1. **Main Results**: Shows only primary Non-B DNA motifs (e.g., 52 motifs)
2. **Cluster/Hybrid Tab**: Dedicated visualization tab for hybrid/cluster motifs (e.g., 71 motifs)
3. **Downloads**: Export files contain only regular motifs for cleaner downstream analysis
4. **Clear Messaging**: Info indicators guide users to the appropriate location

### Benefits

- ✅ **Cleaner Results**: Main results focus on primary motifs
- ✅ **Focused Analysis**: Analyze regular motifs without complex overlaps
- ✅ **Advanced Access**: Hybrid/cluster data still available in dedicated tab
- ✅ **Better Downloads**: Export files contain clean datasets

## 📄 Citation

If you use NonBScanner in your research, please cite:

```
NonBScanner: Comprehensive Detection and Analysis of Non-B DNA Motifs
Dr. Venkata Rajesh Yella
GitHub: https://github.com/VRYella/NonBScanner
```

## 🤝 Contributing

We welcome contributions! Please see our contributing guidelines and submit pull requests for:
- New motif detection algorithms
- Additional visualization features
- Performance improvements
- Documentation enhancements

## 📞 Contact

**Dr. Venkata Rajesh Yella**
- Email: yvrajesh_bt@kluniversity.in
- GitHub: [@VRYella](https://github.com/VRYella)

## 📜 License

This project is licensed under the MIT License - see the LICENSE file for details.

---

*Developed for the scientific community to advance our understanding of Non-B DNA structures and their biological roles.*
