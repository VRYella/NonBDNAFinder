# NBDFinder

[![License](https://img.shields.io/badge/license-MIT-green.svg)](./LICENSE)

## Overview

NBDFinder is an open-source computational framework for genome-scale detection of non-B DNA–forming sequence motifs. The platform integrates classical and emerging non-canonical DNA structures within a unified analytical workflow, enabling systematic identification, subclass annotation, overlap resolution, and structural hotspot detection. The system supports both command-line/API usage and a web-based interactive interface. An interactive Jupyter notebook (`NonBDNAFinder_Analysis.ipynb`) is also included for step-by-step exploratory analysis with visualizations.

## Supported Structural Classes

NBDFinder detects 11 major structural classes encompassing 22+ subclasses. Curved DNA is identified through both global A/T-tract phasing and local long-tract curvature. Slipped DNA covers direct repeats and short tandem repeats (STRs). Cruciform DNA is detected via thermodynamic seed-and-extend analysis of inverted repeats. Triplex DNA encompasses H-DNA mirror repeats and Sticky DNA (GAA/TTC runs). R-loop–forming sequences are identified using the QmRLFS-based model. Z-DNA is scored using the Ho 1986 10-mer table, and eGZ (extruded-guanine Z-DNA) is detected via trinucleotide repeat patterns (CGG, GGC, CCG, GCC) following Herbert 1997. The G-quadruplex family includes canonical, relaxed (extended-loop), bulged, multimeric (stacked), higher-order arrays (G-wire), G-triplex, and weak PQS subclasses, all scored by the seeded G4Hunter algorithm. The i-Motif family covers canonical C-rich four-tract structures and AC-motif variants. A-philic DNA propensity regions are detected by 10-mer scoring. Hybrid regions identify multi-class overlaps, and structural clusters mark density-based hotspot regions. All motif calls include genomic coordinates, subclass labels, lengths, and class-specific scores.

## Key Capabilities

NBDFinder provides integrated multi-class detection in a single workflow, subclass-specific scoring and prioritization, overlap-aware hybrid motif identification, and density-based structural hotspot detection. The platform supports genome-scale analysis of multi-megabase sequences and produces exportable structured outputs in CSV/XLSX/BED-ready formats alongside publication-quality visualizations.

## Input

Supported input formats include FASTA (single or multi-sequence), direct nucleotide sequence input, and NCBI accession retrieval via Biopython. Sequences are automatically normalized (uppercase conversion, U→T substitution, filtering of non-ATGC characters).

## Output Schema

Core output columns are: `Sequence_Name`, `Class`, `Subclass`, `Start` (1-based inclusive), `End`, `Length`, `Strand`, `Score`, `Detection_Method`. Motif-specific metadata such as repeat unit, G-run count, and loop length are reported when applicable. Hybrid and cluster annotations include contributing motif classes, class diversity, and density score.

## Installation

```bash
git clone https://github.com/VRYella/NonBDNAFinder.git
cd NonBDNAFinder
pip install -r requirements.txt
```

Optional performance enhancements:

```bash
pip install numba
pip install cython
```

Run the web interface:

```bash
streamlit run app.py
```

## Programmatic Usage

```python
from Utilities.nonbscanner import analyze_sequence

sequence = "AGGGGGGGGGCCCCCCCCCTAGGGGGGGGG"
results = analyze_sequence(sequence, "example")
print(len(results))
```

Multi-FASTA support:

```python
from Utilities.nonbscanner import analyze_fasta_parallel

results = analyze_fasta_parallel("genome.fasta")
```

## Jupyter Notebook

The repository includes `NonBDNAFinder_Analysis.ipynb`, an interactive notebook for exploratory analysis. It provides step-by-step motif detection, annotated results tables, and built-in visualizations. It is the recommended starting point for users who prefer a guided, cell-by-cell workflow over the command-line API or web interface.

## Performance

NBDFinder operates with linear time complexity with respect to sequence length. Pattern matching is optimized using Hyperscan where available. Parallel detector execution is used for large sequences, and a constant-memory architecture supports multi-megabase genomes. Performance depends on hardware configuration and enabled optimizations.

## Web Application

The Streamlit-based web interface provides interactive motif selection, real-time execution metrics, linear motif maps, class distribution plots, hybrid and cluster visualization, and downloadable result tables.

To launch locally:

```bash
streamlit run app.py
```

## Reproducibility

All scoring parameters are documented in `Utilities/consolidated_registry.json`. Motif detection is deterministic, and the codebase is versioned and open-source.

## System Requirements

Python ≥ 3.8 with NumPy, pandas, matplotlib, seaborn, Streamlit, and Biopython. Optional: Hyperscan, Numba, Cython.

## Citation

If you use NBDFinder in your research, please cite:

```bibtex
@article{NBDFinder2025,
  author  = {Yella, Venkata Rajesh and colleagues},
  title   = {Non-B DNA Finder: A unified framework for detection of diverse non-canonical DNA structures},
  journal = {Nucleic Acids Research},
  year    = {2025}
}
```

## License

MIT License — see [LICENSE](./LICENSE) for details.