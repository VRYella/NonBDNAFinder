# NBDFinder

[![License](https://img.shields.io/badge/license-MIT-green.svg)](./LICENSE)

## Overview

NBDFinder is an open-source computational framework for genome-scale detection of non-B DNA–forming sequence motifs. The platform integrates classical and emerging non-canonical DNA structures within a unified analytical workflow, enabling systematic identification, subclass annotation, overlap resolution, and structural hotspot detection. The system supports both command-line/API usage and a web-based interactive interface. An interactive Jupyter notebook (`NonBDNAFinder_Analysis.ipynb`) is also included for step-by-step exploratory analysis with visualizations.

## Materials and Methods

### Detection Overview

NBDFinder implements nine specialised structural detectors that run independently (optionally in parallel) on the input sequence. Overlap resolution, hybrid annotation, and cluster detection are applied as post-processing steps, yielding **11 output classes** in total (9 structural + Hybrid + Non-B DNA Clusters) encompassing **24 subclasses**. Every motif call carries genomic coordinates (1-based inclusive Start/End), strand, subclass label, length, and a normalised confidence score on a 1–3 scale.

### Detector Descriptions

**Curved DNA** — Global intrinsic curvature is detected by A/T-tract phasing analysis: runs of ≥ 3 consecutive adenines (A-tracts) or thymines (T-tracts) are identified and their inter-tract centre-to-centre spacing is evaluated against the 10.5 bp helical repeat (tolerance 9.9–11.1 bp; minimum 3 in-phase tracts). Local curvature is additionally reported for long uninterrupted A-tracts or T-tracts (≥ 8 nt). Scoring follows Koo *et al.* (1986) and Olson *et al.* (1998); normalised linearly to the 1–3 scale.

**Slipped DNA** — Short tandem repeats (STRs) are detected using a k-mer index approach: mono- through tetranucleotide units (1–4 nt) with ≥ 6 copies are identified and scored by a mechanistic slippage model (Shannon entropy × copy number × GC content; Schlötterer *et al.* 2000; Weber *et al.* 1989). Longer direct repeats (unit 10–50 nt, ≥ 2 copies) are detected separately by the same engine. Subclasses reported: *STR* and *Direct Repeat*.

**Cruciform** — Inverted repeats (IRs) capable of cruciform extrusion are located by a seed-and-extend algorithm. A 6-mer seed dictionary indexes all reverse-complement seeds in the sequence; candidate IR pairs are extended up to arm lengths of 8–50 nt with zero mismatches. Thermodynamic stability is evaluated using the unified nearest-neighbor model (SantaLucia 1998): only stems with ΔG < −5.0 kcal/mol and a loop penalty-adjusted score above 0.2 are retained. Scoring reference: Lilley *et al.* (2000). Subclass reported: *Cruciform forming IRs*.

**R-Loop** — R-loop–forming sequences (RLFS) are identified using a faithful implementation of the QmRLFS-finder algorithm (Jenjaroenpun *et al.* 2016). The algorithm defines an R-loop–initiation zone (RIZ) as a G-cluster region containing overlapping G-tracts (≥ 3 Gs separated by ≤ 10 nt linkers, G-content ≥ 50 %), and an R-loop–elongation zone (REZ) extending downstream (up to 2000 nt, G-content ≥ 40 %). Two QmRLFS models are applied in parallel: Model 1 (standard G-cluster RIZ) and Model 2 (extended G-tract RIZ, ≥ 4-G tracts). Quality scores ≥ 0.4 are reported. References: Aguilera & García-Muse (2012); Jenjaroenpun *et al.* (2016). Subclass reported: *R-loop formation sites*.

**Triplex DNA** — H-DNA forming mirror repeats are detected by a seed-and-extend purity scanner: a 6-mer seed index locates candidate mirror-repeat arms (10–100 nt, ≥ 90 % purine/pyrimidine purity, loop ≤ 8 nt). Stability is scored via the arm-length, loop-penalty, purity, and interruption model of Frank-Kamenetskii *et al.* (1995); threshold 0.25. Sticky DNA (GAA/TTC trinucleotide expansions) is detected as a separate subclass using piecewise linear copy-number scoring aligned to replication-blockage (≥ 20 copies), stable-Sticky (≥ 40), and pathogenic FRDA thresholds (≥ 60; Sakamoto *et al.* 1999). Subclasses: *H-DNA*, *Sticky DNA*.

**G-Quadruplex** — A seeded G4Hunter algorithm (Bedrat *et al.* 2016) is used for all G4 subclasses. G-run seeds are located by a bisect-indexed position list; a sliding window (default 25 nt) computes the mean G4Hunter score (|G|−|C| per base) and regions above 0.5 are extended. Eight hierarchically prioritized subclasses are resolved in a single-pass overlap resolution step: *Telomeric G4* (sequence-specific TTAGGG/TTGGGG arrays), *G-wire* (higher-order G4 arrays, ≥ 4 stacked tetrads), *Stacked G4* (multi-quadruplex assemblies), *Canonical G4* (four G-tracts, loops 1–7 nt), *Bulged G4* (canonical with single-base bulge), *Extended-loop G4* (loops up to 12 nt), *G-triplex* (three G-tracts), and *Weak PQS* (two-tetrad structures). Optional Numba JIT compilation accelerates the sliding-window computation.

**i-Motif** — Canonical i-motifs are detected as four C-rich tracts (C ≥ 3, inter-tract loops 1–7 nt) using a direct regex search; C-run density and tract regularity drive the score (Gehring *et al.* 1993; Zeraati *et al.* 2018). AC-motif variants following the HUR model (Hur *et al.* 2021) are detected as alternating A₃–C₃ or C₃–A₃ patterns with 4–6 nt spacers. Subclasses: *Canonical i-Motif*, *AC-motif (HUR)*.

**Z-DNA** — Classical Z-DNA–prone regions are scored by the cumulative 10-mer propensity table of Ho *et al.* (1986 EMBO J): every overlapping 10-mer in the sequence is scored and adjacent high-scoring 10-mers are merged; minimum merged score 50. eGZ (extruded-guanine Z-DNA) motifs following Herbert (1997) are detected as runs of ≥ 4 trinucleotide repeats from the set {CGG, GGC, CCG, GCC}. Log-linear normalisation is applied to accommodate the wide dynamic range of cumulative Z-DNA scores. Optional Hyperscan acceleration is used when available. Subclasses: *Z-DNA*, *eGZ*.

**A-philic DNA** — A-philic propensity regions are identified using the 10-mer scoring table derived from Gorin *et al.* (1995) and Vinogradov *et al.* (2003). All overlapping 10-mers with positive log₂ A-philic propensity scores are located; adjacent high-scoring 10-mers are merged into contiguous A-philic regions (minimum merged sum-log₂ score 0.5). Optional Hyperscan acceleration is used when available. Subclass: *A-philic DNA*.

**Hybrid regions** — After all nine structural detectors have completed, motif pairs from distinct classes that share ≥ 50 % positional overlap (relative to the shorter motif) are consolidated into *Hybrid* annotations. Each hybrid record reports all contributing classes and class diversity.

**Non-B DNA Clusters** — A density-based scan anchored at each detected motif start identifies genomic positions where ≥ 4 structurally distinct non-B DNA motifs from ≥ 3 unique classes co-occur within a 300 nt window. Cluster records report motif count, class diversity, and the window boundaries.

### Sequence Requirements

Sequences shorter than 10 bp are automatically skipped with a warning; no error is raised and analysis continues with remaining sequences in the file. Sequences must contain only IUPAC nucleotide characters (A, T, G, C, N and standard ambiguity codes); unrecognised characters trigger a validation error.

## Key Capabilities

NBDFinder provides integrated multi-class detection in a single workflow, subclass-specific scoring and prioritization, overlap-aware hybrid motif identification, and density-based structural hotspot detection. The platform supports genome-scale analysis of multi-megabase sequences and produces exportable structured outputs in CSV/XLSX/BED-ready formats alongside publication-quality visualizations.

## Input

Supported input formats include FASTA (single or multi-sequence), direct nucleotide sequence input, and NCBI accession retrieval via Biopython. Sequences are automatically normalized (uppercase conversion, U→T substitution, filtering of non-ATGC characters).

## Output Schema

Core output columns are: `Sequence_Name`, `Class`, `Subclass`, `Start` (1-based inclusive), `End`, `Length`, `Strand`, `Score`, `Detection_Method`. Motif-specific metadata such as repeat unit, G-run count, and loop length are reported when applicable. Hybrid and cluster annotations include contributing motif classes, class diversity, and a composite score.

## Example Datasets

The `examples/` directory contains ready-to-use FASTA files for testing:

| File | Description |
|------|-------------|
| `examples/example_single.fasta` | Single sequence with mixed non-B DNA motifs |
| `examples/example_multi.fasta` | Multi-sequence FASTA covering all major motif classes |
| `examples/example.fasta` | Original test sequences (G-runs, C-runs, CGCG repeats) |

These files are available in the web interface under **Upload & Analyze → Example Data**, and can be used directly from the command line or Jupyter notebook.

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