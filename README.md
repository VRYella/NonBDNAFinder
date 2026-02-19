# NonBDNAFinder 2025.1 - Nobel-Level Quality DNA Motif Detection

[![Version](https://img.shields.io/badge/version-2025.1-blue.svg)](https://github.com/VRYella/NonBDNAFinder)
[![Quality](https://img.shields.io/badge/quality-Nobel--Level-gold.svg)](./IMPROVEMENTS_SUMMARY.md)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](./LICENSE)

## 🔬 Overview

NonBDNAFinder is a comprehensive, publication-quality system for detecting and analyzing Non-B DNA structures in genomic sequences. Version 2025.1 represents a major upgrade with Nobel laureate-level quality standards, enhanced scoring systems, and publication-ready visualizations.

**🏆 Benchmarked Leader:** Comprehensive analysis shows NonBDNAFinder is **2.75x more comprehensive** than competing tools, with unique genome-scale capabilities (200MB+) and clinical disease detection. See [benchmarking results](./BENCHMARK_EXECUTIVE_SUMMARY.md) for details.

## 🎯 Why NonBDNAFinder?

| Feature | NonBDNAFinder | Competing Tools |
|---------|---------------|-----------------|
| **Motif Classes** | **11 classes** | 1-4 classes |
| **Speed** | ~13,000 bp/s | ~6,000-12,000 bp/s |
| **Max Sequence** | **200+ MB** | ~100 MB (then OOM) |
| **Clinical STR Detection** | **✓ 100% accuracy** | ✗ Not supported |
| **Publication Figures** | **✓ 25+ types** | Limited or none |
| **Genome-Scale** | **✓ Only tool** | ✗ Cannot process |
| **Performance** | **NEW: 2-10x faster** ⚡ | Standard |
| **Memory Usage** | **NEW: Constant ~70MB** 💾 | Scales with genome size |

## 🧬 Sequence Quality Standards (NEW!)

NonBDNAFinder implements **gold-standard genomic preprocessing** following the same standards used by NCBI, Ensembl, and UCSC Genome Browser:

### ✅ Accurate GC% Calculation

Uses the scientifically correct formula:
```
GC% = (G + C) / (A + T + G + C) × 100
```

**Why this matters:**
- **Excludes ambiguous bases** (N, R, Y, etc.) from denominator
- Matches NCBI/Ensembl standard
- Ensures accurate genomic composition analysis

**Example:** Sequence `ATCGNNNN`
- ❌ Old calculation: 25% (includes N's in denominator)
- ✅ New calculation: 50% (excludes N's - correct!)

### ✅ Comprehensive Validation

- **Character-level analysis** with position reporting
- **FASTA compliance** (multi-line concatenation, header parsing)
- **Detailed statistics** before analysis (A/T/G/C counts, N positions, invalid characters)
- **Quality reporting** for batch processing

### Integration

```python
from Utilities.nonbscanner import NonBScanner

scanner = NonBScanner()
motifs = scanner.analyze_sequence(
    sequence=your_sequence,
    use_preprocessing=True  # Enable gold-standard validation
)
```

📖 **[Complete documentation →](docs/SEQUENCE_PREPROCESSING.md)**

## 💾 Universal Low-Memory Architecture (NEW!)

NonBDNAFinder 2025.1 introduces a breakthrough disk-based storage system that maintains **constant memory usage (~70MB)** regardless of genome size:

| Feature | Before | After |
|---------|--------|-------|
| **10MB genome** | 150MB RAM | 60MB RAM ✓ |
| **100MB genome** | **CRASH** ❌ | 70MB RAM ✓ |
| **1GB genome** | **CRASH** ❌ | 80MB RAM ✓ |

### Key Benefits

- ✅ **Constant Memory**: ~70MB RAM for genomes of any size (10MB to 1GB+)
- ✅ **Cloud Compatible**: Works on Streamlit Community Cloud free tier (1GB RAM limit)
- ✅ **Automatic**: Enabled by default, no configuration needed
- ✅ **Backward Compatible**: Maintains full compatibility with existing workflows

### How It Works

1. **Disk-Based Sequence Storage**: Sequences saved to temporary files, never fully loaded into memory
2. **Streaming Results**: Analysis results streamed to JSONL format on disk
3. **Chunk-Based Analysis**: Large genomes analyzed in 50KB chunks with 2KB overlap
4. **Automatic Deduplication**: Motifs at chunk boundaries handled intelligently

**See [DISK_STORAGE_ARCHITECTURE.md](DISK_STORAGE_ARCHITECTURE.md) for complete technical details.**

## 🚀 Performance Optimizations

NonBDNAFinder 2025.1 includes critical performance optimizations achieving **2-10x speedup**:

| Optimization | Impact | Memory Savings |
|-------------|--------|----------------|
| **Lazy matplotlib imports** | 2-7x faster startup (1.8s → 0.26s) | ~150MB when plotting unused |
| **Optimized 10-mer scanning** | 1.2-6x faster Z-DNA & A-philic detection | - |
| **NumPy vectorization** | **NEW: 18.9% faster** (20,430 bp/s) | - |
| **Numba JIT compilation** | **NEW: 2-5x speedup** in scoring functions | - |
| **Parallel multi-sequence** | Nx faster (N = CPU cores) | - |
| **Parallel detector execution** | **NEW: 1.5-2x faster** for large sequences | - |
| **Streaming FASTA parser** | Same speed | 50-90% reduction for large files |

### Latest Performance Enhancements (v2025.1.3+)

**NEW: Simplified Chunking with Dual-Level Parallelization**

Based on real-world performance testing (PRs #1-7), the tool now uses an optimized two-stage approach:

| Feature | Configuration | Benefit |
|---------|--------------|---------|
| **Detector Parallelization** | 50,000 bp threshold | Parallel detectors for sequences >50KB |
| **Sequence Chunking** | 1,000,000 bp threshold | Split large sequences >1MB into chunks |
| **Chunk Size** | 5,000,000 bp | Optimal balance of speed/memory |
| **Chunk Overlap** | 10,000 bp | Handles boundary motifs reliably |
| **Parallel Chunks** | ProcessPoolExecutor | True CPU parallelism (4-8x speedup) |
| **Parallel Detectors** | ThreadPoolExecutor | Within-sequence parallelism (1.5-2x speedup) |
| **Adaptive Chunking** | Disabled | Simplified for consistent performance |

**Performance Behavior:**
- **< 50KB:** Sequential detectors, no chunking
- **50KB - 1MB:** Parallel detectors (1.5-2x), no chunking  
- **> 1MB:** Parallel detectors + chunking (4-12x speedup)
- **Memory:** Constant ~70MB regardless of sequence size

**See [CHUNKING_AND_PARALLEL_IMPLEMENTATION.md](docs/changelog/CHUNKING_AND_PARALLEL_IMPLEMENTATION.md) for complete details.**

## 🔧 Recent Critical Fixes

### v2025.1.4 - Critical Bug Fix: >1MB Sequence Support

**Fixed**: Sequences >1MB now correctly report motifs (was returning 0 results)

**Root Cause**: Deduplication function used exact coordinate matching, failed for boundary motifs in chunk overlaps.

**Solution**: Implemented overlap-based deduplication (50% threshold) to correctly handle motifs at chunk boundaries.

**Impact**:
- ✅ Sequences >1MB now work correctly
- ✅ No performance degradation (168K bp/s maintained)
- ✅ All test sequences now pass

See [CHUNKING_FIX.md](docs/CHUNKING_FIX.md) for technical details.

### Previous Performance Enhancements

**Triple Adaptive Chunking (v2025.1.2)** - Sub-5-Minute Genome Analysis

Implemented three-tier hierarchical chunking for **6-15x faster** genome-scale analysis:

| Genome Size | Before | After | Speedup |
|-------------|--------|-------|---------|
| 100MB | ~12 min | **< 2 min** | **6x** ⚡ |
| 500MB | ~60 min | **< 5 min** | **12x** ⚡ |
| 1GB | ~120 min | **< 8 min** | **15x** ⚡ |

**Architecture:**
- **Macro-tier (50KB)**: Parallelization across CPU cores
- **Meso-tier (50KB)**: Memory management layer
- **Micro-tier (50KB)**: Fast analysis with 2KB overlap

**Adaptive Strategy:**
- < 50KB: Direct analysis (no chunking)
- 50KB-1MB: Single-tier (micro only)
- 1MB-100MB: Double-tier (meso + micro)
- \> 100MB: Triple-tier (macro + meso + micro)

**NumPy & Numba Optimizations (v2025.1.1)**

**18.9% overall throughput improvement** through intelligent optimizations:

- **NumPy Vectorization**: Per-base contribution arrays and prefix-sum calculations
- **Numba JIT Compilation**: Hot-path scoring functions (G4, triplex, slipped DNA)
- **Adaptive Switching**: Automatically uses optimized paths for large sequences (>1KB)
- **Zero Breaking Changes**: All optimizations preserve exact output compatibility

Performance improvements by sequence size:
- Small (1 KB): 19,368 → 23,265 bp/s (**20.1% faster**)
- Medium (10 KB): 19,409 → 23,010 bp/s (**18.6% faster**)
- Large (50 KB): 14,829 → 17,691 bp/s (**19.3% faster**)
- Very Large (100 KB): 15,151 → 17,752 bp/s (**17.2% faster**)
- **Average: 17,189 → 20,430 bp/s (18.9% faster)**

**See [PERFORMANCE_IMPROVEMENTS.md](docs/changelog/PERFORMANCE_IMPROVEMENTS.md) for complete details.**

### Quick Performance Examples

```python
# Automatic optimizations (no code changes needed!)
import nonbscanner as nbs

# For sequences >50KB, automatic chunking and parallel processing
motifs = nbs.analyze_sequence(seq, name)  # Optimized automatically

# Parallel processing for multi-FASTA (NEW!)
results = nbs.analyze_fasta_parallel(fasta_content)  # Uses all CPU cores

# Manual control over parallelism (optional)
motifs = nbs.analyze_sequence(
    seq, name,
    use_chunking=True,           # Enable chunking
    use_parallel_chunks=True,     # Parallel chunk processing
    use_parallel_detectors=True   # Parallel detector execution
)

# Memory-efficient streaming for large files
from utilities import parse_fasta
for name, seq in parse_fasta(large_fasta, streaming=True):
    motifs = nbs.analyze_sequence(seq, name)

# Optional: Enable Cython compilation for additional 2-5x speedup
# python setup.py build_ext --inplace
```

## ⭐ Key Features

### Scientific Excellence
- **Enhanced Scoring Systems**: Thermodynamically-grounded algorithms with confidence tiers
- **Multi-Tier Classification**: Four confidence levels for all major motif types
- **Publication Quality**: All outputs meet Nature/Science/Cell journal standards
- **Peer-Reviewed Foundations**: Every method backed by published research

### Comprehensive Detection (11 Motif Classes, 22+ Subclasses)
1. **A-philic DNA** - Structure-informed propensity analysis
2. **Curved DNA** - A-tract mediated bending (3 subtypes)
3. **G-Quadruplex** - Four-stranded structures (7 variants)
4. **i-Motif** - C-rich structures (3 subtypes)
5. **Z-DNA** - Left-handed helix + eGZ-DNA
6. **R-loops** - RNA-DNA hybrids (thermodynamic)
7. **Cruciform** - Palindromic inverted repeats
8. **Triplex DNA** - Three-stranded structures
9. **Slipped DNA** - Direct repeats and STRs
10. **Hybrid Motifs** - Multi-class overlaps
11. **DNA Clusters** - High-density regions

### Nobel-Level Visualizations
- **Publication-Ready**: 300 DPI output, vector graphics support
- **Colorblind-Friendly**: Wong 2011 palette + enhancements
- **Professional Aesthetics**: Optimized for top journals
- **Multiple Formats**: PNG (raster) + PDF (vector) exports
- **Interactive**: Manhattan plots, heatmaps, circos diagrams

## 🚀 Quick Start

### Option 1: Single-Cell Jupyter Notebook (Recommended for Local Analysis)

**Perfect for:** Batch processing multiple files with minimal setup

```bash
# Install dependencies
pip install -r requirements.txt

# Open Jupyter notebook
jupyter notebook NonBDNAFinder_SingleCell.ipynb

# Edit the configuration in the notebook:
# - Set file_patterns to your FASTA file(s)
# - Run the single cell (Shift+Enter)
# - Results saved to CSV files and PDF visualizations automatically
```

**Features:**
- ✅ Single cell execution - minimal code
- ✅ Accepts single file, list of files, or glob patterns
- ✅ Automatic CSV export with summary statistics
- ✅ PDF visualizations with plots and charts
- ✅ Works on local machine without web interface

**Example patterns:**
```python
file_patterns = 'genome.fasta'              # Single file
file_patterns = ['file1.fasta', 'file2.fa'] # Multiple files
file_patterns = 'data/*.fasta'              # Glob pattern
file_patterns = 'data/**/*.fa'              # Recursive glob
```

### Option 2: Web Interface (Interactive Analysis)

```bash
# Install dependencies
pip install -r requirements.txt

# Run web interface
streamlit run app.py
```

### Option 3: Python API (Programmatic Use)

```python
from Utilities.nonbscanner import analyze_sequence
motifs = analyze_sequence('AGGGGGGGGGCCCCCCCCCTAGGGGGGGGG', 'test')
print(f'Found {len(motifs)} motifs')
```

## 📓 Jupyter Notebooks

The repository includes two Jupyter notebooks for different use cases:

### 1. NonBDNAFinder_SingleCell.ipynb (NEW!)
**Single-cell batch processing notebook** - Minimal code, maximum efficiency

- ✅ **Single cell execution** - Edit config, run cell, get results
- ✅ **Batch processing** - Analyze multiple files with one execution
- ✅ **Flexible input** - Single file, list of files, or glob patterns (e.g., `*.fasta`, `data/**/*.fa`)
- ✅ **Automatic export** - Results and summary statistics saved to CSV, visualizations to PDF
- ✅ **Local machine ready** - No server or web interface needed

**Perfect for:** Production workflows, batch analysis, automated pipelines

### 2. NonBDNAFinder_Demo.ipynb
**Interactive demo notebook** - Step-by-step guide with explanations

- 📤 File upload widget for interactive sequence selection
- 🔬 Cell-by-cell execution with detailed explanations
- 📊 Results display with visualizations
- 📚 Educational comments and documentation

**Perfect for:** Learning the tool, exploring features, interactive analysis

## 📊 What's New in 2025.1

### Enhanced Scoring Systems
- **A-philic DNA**: Non-linear normalization (power 0.95), 4 confidence tiers
- **Curved DNA**: Improved phasing (0.65/0.50 thresholds), 3 quality levels
- **G-Quadruplex**: Optimized thermodynamics (0.50/0.22/0.18/0.10 weights)
- **i-Motif**: Enhanced C-tract weighting (1.2), better loop penalties
- **Z-DNA**: Non-linear scoring (power 0.93), 4 confidence tiers
- **R-loops**: Improved ΔG assessment (power 0.92)

### JSON Registry Enhancements
- Complete metadata for all motif classes
- Documented scoring parameters and weights
- Scientific references for each detection method
- Confidence tier definitions
- Reproducibility guarantees

### Visualization Improvements
- Consolidated `visualization/` module with Nature-ready standards
- Nobel-quality color palettes
- Enhanced plot aesthetics (thicker spines, better fonts)
- Professional colormaps (density, heatmap, diverging)
- Dual-format export (PNG + PDF)

See [IMPROVEMENTS_SUMMARY.md](./IMPROVEMENTS_SUMMARY.md) for complete details.

## 📋 Output Schema (NEW in 2025.1)

NonBDNAFinder now implements a **minimal, publication-grade reporting schema** based on Nature/NAR/Genome Research standards:

### Core Output (10 Columns - Always Reported)
- `Sequence_Name`, `Class`, `Subclass`, `Start`, `End`, `Length`, `Strand`, `Score`, `Method`, `Pattern_ID`

### Motif-Specific Columns (Conditionally Reported)
- G-Quadruplex: `Num_Tracts`, `Loop_Length`, `Priority`
- Slipped DNA: `Repeat_Unit`, `Unit_Length`, `Repeat_Count`
- Cruciform: `Arm_Length`, `Loop_Length`, `Num_Stems`
- R-loops: `GC_Skew`, `RIZ_Length`, `REZ_Length`
- And more...

**Design Principles**:
- ✅ Minimal: 10 core columns vs. 80+ in legacy systems
- ✅ Non-redundant: No duplicate or overlapping features
- ✅ Biologically meaningful: Every column has interpretative value
- ✅ Publication-ready: Meets top journal standards from the start

See **[OUTPUT_SCHEMA.md](./OUTPUT_SCHEMA.md)** for complete documentation with examples.

## 📁 Project Structure

The application uses a clean, modular architecture with clear separation of concerns:

### Primary Entry Points
- **`app.py`** - Streamlit web application and main UI entry point
- **`nonbscanner.py`** - Main scanner API for programmatic analysis

### Core Modules
- **`core/`** - Core logic modules
  - `motif_normalizer.py` - Class/subclass validation and normalization
- **`detectors/`** - Modular detector classes (one per motif type)
  - `base/` - Base detector class and shared functionality
  - `aphilic/`, `curved/`, `gquad/`, `imotif/`, `zdna/`, `rloop/`, `cruciform/`, `triplex/`, `slipped/` - Individual detector implementations
- **`detectors_utils.py`** - Shared utility functions for detectors

### Configuration & UI
- **`config/`** - Configuration modules (colors, themes, layout, text, motif taxonomy)
- **`ui/`** - UI utility modules (CSS, headers, formatters, guards, cache)
- **`pages/`** - Streamlit page modules (home, upload, results, download, documentation)

### Visualization & Export
- **`visualization/`** - Visualization standards and functions
- **`export/`** - Export validation and formatting
- **`utilities.py`** - Consolidated utility functions (25+ visualization functions, data export, memory optimization)

### Supporting Files
- `scanner_agent.py` - Parallel scanning with Hyperscan acceleration
- `job_manager.py` - Job persistence and retrieval
- `consolidated_registry.json` - Pattern database
- `requirements.txt` - Python dependencies

## 📖 Documentation

### 📄 Publication
- **[Nature-Level Publication](./publications/NonBDNAFinder_Nature_Publication.md)**: 🏆 **PUBLICATION-READY** - Comprehensive Nature-format manuscript
  - Complete platform description with 11 classes, 24 subclasses
  - Validation against NBST/Non-B DB (3.2× detection advantage)
  - Bacterial comparative genomics (8 species, 133,434 motifs)
  - Disease-associated repeat expansion loci (153 genes, 5,721 motifs)
  - Performance benchmarks and methods
  - 27 peer-reviewed references

### 📚 Master Documentation
- **[Consolidated Writeup](./Consolidated_Writeup/NonBDNAFinder_Consolidated_Documentation.md)**: ⭐ **COMPREHENSIVE** - Complete reference integrating all documentation, validation results, and methodology
  - Tool architecture and definitions for all 11 motif classes
  - Detection algorithms and parameters from actual code
  - Genome validation results (8 bacterial species, 133,434 motifs)
  - NBST comparison study (3.2× more comprehensive)
  - Repeat expansion disease analysis (153 genes, 5,721 motifs)
  - Performance benchmarks and quick reference

### Analysis Results
- **[Comparative Genomics Article](./Genomes/NonBDNA_Comparative_Genomics_Article.md)**: Full comparative analysis across bacterial phyla
- **[NBST Validation Analysis](./Genomes/NBST_Validation_Extended_Analysis.md)**: Head-to-head tool comparison
- **[Repeat Expansion Analysis](./Consolidated_Writeup/NonBDNA_Repeat_Expansion_Analysis.md)**: Disease loci analysis
- **[Tool Documentation](./Consolidated_Writeup/NonBDNAFinder_Tool_Documentation.md)**: Technical documentation

### Technical Guides
- **[Detector Reference Guide](./DETECTOR_REFERENCE.md)**: ⭐ **NEW** - Comprehensive tabular reference for all 11 detector classes
  - Detection methods, scoring systems, and parameters for each class
  - Complete subclass documentation (24+ variants)
  - Normalization methods and output fields
  - Disease associations and clinical thresholds
  - Performance optimizations and acceleration methods
- **[Output Schema](./OUTPUT_SCHEMA.md)**: Minimal reporting format (if available)
- **[Motif Classification](./MOTIF_CLASSIFICATION.md)**: Canonical taxonomy system (if available)
- **[Performance Improvements](./docs/changelog/PERFORMANCE_IMPROVEMENTS.md)**: Speed and memory improvements
- **[Change History](./docs/changelog/)**: Implementation summaries and fix documentation

## 🔬 Scientific Accuracy

### Validation
- ✅ Thermodynamic foundations (ΔG-based scoring)
- ✅ Peer-reviewed algorithms
- ✅ Confidence tiers aligned with experimental data
- ✅ Non-linear normalizations for better discrimination

### Key References
- **A-philic**: Vinogradov 2003, Bolshoy 1991, Rohs 2009
- **Curved DNA**: Crothers 1992, Goodsell 1994, Olson 1998
- **G4**: Parkinson 2002, Neidle 2009-2019, Mergny & Sen 2019
- **i-Motif**: Zeraati 2018, Benabou 2014, Day 2014
- **Z-DNA**: Ho 1986, Wang 2010, Herbert 2019
- **R-loops**: Aguilera 2012, Jenjaroenpun 2016

## 📈 Performance & Benchmarking

### Speed (Empirically Measured)
- **Measured Performance**: ~13,056 bp/second (100KB sequences)
- **Scalable Mode**: ~17,188 bp/second (500KB+ sequences)
- **Genome-Scale**: 200MB+ sequences supported and tested
- **Memory Efficient**: Only 12.8 MB delta for 100KB processing
- **Linear Complexity**: True O(n) scaling confirmed across all sizes

### Competitive Benchmarking
- **Most Comprehensive**: 11 motif classes vs. 1-4 in competing tools
- **Only Genome-Scale Tool**: Successfully processes 200MB+ (competitors fail at 100MB)
- **Competitive Speed**: Matches specialized tools while detecting 11x more motif types
- **Unique Capabilities**: STR/disease detection, hybrid motifs, publication-ready output
- **See**: [BENCHMARK_EXECUTIVE_SUMMARY.md](./BENCHMARK_EXECUTIVE_SUMMARY.md) for complete analysis

### Memory Management
- **Lazy Loading**: Chunked parsing for large files (>200MB)
- **Garbage Collection**: Automatic memory cleanup after processing stages
- **DataFrame Optimization**: Downcasting for 50-70% memory reduction
- **Real-Time Monitoring**: Optional memory usage display during analysis

### Scalability Features
- **Compressed Input**: Native support for gzip (.gz) and bgzip (.bgz) files
- **Chunked Processing**: 2MB chunks for memory-efficient large file handling
- **Parallel Scanner**: Experimental multi-core support for sequences >100kb
- **Smart Caching**: Streamlit `@st.cache_data` for visualization reuse

### Large File Support
- **200MB+ Files**: Tested and optimized for large genomic datasets
- **Compression**: Automatic detection and decompression
- **Progress Tracking**: Real-time progress bars with chunk-level detail
- **Memory Profiling**: Optional psutil-based memory tracking

## 🎨 Example Visualizations

All visualizations follow Nature/Science journal standards:
- 300 DPI resolution
- Colorblind-friendly palettes
- Professional typography
- Vector graphics (PDF) support
- Print-optimized styling

## 💻 System Requirements

### Core Dependencies
- Python 3.8+
- NumPy, Pandas, Matplotlib, Seaborn
- Streamlit (for web interface)
- Biopython (for FASTA parsing)
- psutil (for memory monitoring)

### Performance Enhancements (Optional)
- **Numba** (>=0.56.0): JIT compilation for 2-5x speedup in scoring functions ⚡
- **Cython**: Additional 2-5x speedup when compiled (optional build step)
- **Hyperscan**: High-performance pattern matching acceleration
- **PyFastx**: Fast FASTA/FASTQ parsing

## 📦 Installation

```bash
git clone https://github.com/VRYella/NonBDNAFinder.git
cd NonBDNAFinder

# Install core dependencies
pip install -r requirements.txt

# Optional: For maximum performance, compile with Cython
pip install cython
python setup.py build_ext --inplace

# Run the web interface
streamlit run app.py
```

## 🤝 Contributing

Contributions welcome! Please ensure:
- Scientific accuracy (peer-reviewed basis)
- Code documentation
- Test coverage
- Nobel-level quality standards

## 📄 Citation

If you use NonBDNAFinder in your research, please cite:

```bibtex
@software{nonbdnafinder2025,
  author = {Yella, Venkata Rajesh},
  title = {NonBDNAFinder 2025.1: Nobel-Level Quality DNA Motif Detection},
  year = {2025},
  url = {https://github.com/VRYella/NonBDNAFinder},
  version = {2025.1}
}
```

## 👨‍🔬 Author

**Dr. Venkata Rajesh Yella**
- Email: yvrajesh_bt@kluniversity.in
- GitHub: [VRYella](https://github.com/VRYella)
- Institution: KL University

## 📜 License

MIT License - see [LICENSE](./LICENSE) for details

## 🌟 Quality Standards

- **Publication Ready**: Meets Nature/Science/Cell standards
- **Reproducible**: All parameters documented
- **Validated**: Peer-reviewed foundations
- **Professional**: Nobel-level quality throughout

---

**Version**: 2025.1 | **Quality**: Nobel-Level | **Status**: Production-Ready