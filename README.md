---
title: Non-B DNA Motif Finder
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: streamlit
sdk_version: "1.49.0"
app_file: app_streamlit.py
pinned: false
---

# NBDFinder - High-Performance Non-B DNA Motif Detection

## 🚀 Core Refactor - Hyperscan-Powered Engine

NBDFinder has been completely refactored with a high-performance Hyperscan-powered core for detecting non-B DNA structures. The new architecture provides:

### ⚡ Key Features
- **Intel Hyperscan acceleration** for multi-pattern scanning
- **5 core scripts** + supporting modules for clean architecture
- **10 motif classes** with >20 subclasses detected
- **Parallel processing** with automatic chunking for large sequences
- **Multiple output formats**: CSV, Excel, Parquet, GFF3
- **Interactive Streamlit UI** with real-time visualizations
- **Comprehensive testing** with CI/CD pipeline

### 🧬 Motif Classes Detected
1. **G-Quadruplex** (6 subtypes: canonical, relaxed, bulged, bipartite, multimeric, imperfect)
2. **Curved DNA** (A-phased repeats with 8-12bp spacing)
3. **Slipped DNA** (STRs, direct repeats)
4. **Cruciform** (palindromes, inverted repeats)
5. **R-loops** (RLFS models)
6. **Triplex DNA** (homopurine/homopyrimidine tracts)
7. **i-Motif** (C-rich quadruplex structures)
8. **Z-DNA** (alternating purine-pyrimidine)
9. **Hybrid motifs** (overlapping structures)
10. **Clusters** (motif hotspots)

## 🏗️ Architecture

### Core Components (5 Scripts)
```
hs_registry_manager.py     # Hyperscan DB compilation & caching
motif_detectors.py         # 10 detector classes with detect() & score()
orchestrator.py            # Pipeline orchestration & parallel processing
viz_tools.py               # Interactive visualizations (plotly)
app_streamlit.py           # 5-page Streamlit interface
```

### Supporting Files
```
motifs/
├── __init__.py           # Package initialization
├── base.py               # Candidate dataclass & scoring functions
└── registry.py           # Pattern management & Hyperscan safety
tests/                    # Comprehensive unit tests
.github/workflows/ci.yml  # CI pipeline with multi-Python testing
```

## 🚀 Quick Start

### Command Line
```bash
# Install with all dependencies
pip install -e .

# Run on FASTA file
python orchestrator.py --fasta sequences.fa --out results --workers 4

# Outputs: results.csv, results.xlsx, results.parquet, results.gff3
```

### Streamlit Interface
```bash
streamlit run app_streamlit.py
```

### Python API
```python
from orchestrator import run_pipeline

output_files = run_pipeline(
    fasta_path="sequences.fa",
    output_prefix="results",
    max_workers=4,
    detector_classes=['g_quadruplex', 'triplex', 'curved_dna']
)
```

## ⚙️ Performance Design

### Hyperscan Integration Rules
1. **Compile Once**: Databases cached at module level
2. **Scratch Per Worker**: One `hs.Scratch(db)` per process/thread
3. **Minimal Callbacks**: Only append `(pattern_id, start, end)`
4. **Fallback Support**: Python regex for unsafe patterns
5. **Memory Efficient**: Streaming & chunking for large sequences

### Parallel Processing
- **Multi-core detection** with `ProcessPoolExecutor`
- **Automatic chunking** for sequences >50kb
- **Overlap resolution** using interval trees
- **Memory streaming** for large genomes

## 📊 Output Schema
```
S.No, Sequence_Name, Chromosome/Contig, Class, Subclass, Motif_ID, 
Start, End, Length, Normalized_Score, Actual_Score, Scoring_Method, 
GC_Content, Sequence, Overlap_Classes
```

## 🧪 Testing & CI

```bash
# Run unit tests
pytest tests/ -v

# Test with coverage
pytest tests/ --cov=. --cov-report=html

# Integration test
python orchestrator.py --fasta test_sample.fa --out test_results
```

## 📚 Scientific References

- **G4Hunter**: Bedrat et al. NAR 44(4):1746-1759 (2016)
- **Z-DNA Seeker**: Ho et al. EMBO J 5:2737-2744 (1986)
- **i-motif**: Zeraati et al. Nat Chem 10:631-637 (2018)
- **Triplex**: Frank-Kamenetskii & Mirkin Annu Rev Biochem 64:65-95 (1995)
- **Hyperscan**: Wang et al. NSDI 2019

## 🔧 Development

### Installing for Development
```bash
git clone https://github.com/VRYella/NBDFinder.git
cd NBDFinder
pip install -e .[dev]  # Install with development dependencies
```

### Running Tests
```bash
pytest tests/                    # All tests
pytest tests/test_hs_compile.py  # Hyperscan tests
pytest tests/test_g4_detection.py # G4 detection tests
```

## 📝 Citation

If you use NBDFinder in your research, please cite:

> NBDFinder: High-performance detection of non-B DNA structures using Hyperscan acceleration. 
> [Authors], [Journal] (Year). DOI: [DOI]

## 📄 License

MIT License - see LICENSE file for details.
