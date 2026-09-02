# Non-B DNA Finder

## Overview

NonBDNAFinder is a scientific software platform for sequence-based detection and annotation of non-B DNA-forming motifs. The repository contains the Streamlit web application (`app.py`), a Jupyter notebook (`NonBDNAFinder_Analysis.ipynb`), example FASTA inputs, and the Python implementation of the motif-detection pipeline.

## Key Features

- Unified analysis of 9 primary motif detector families
- Derived reporting for hybrid overlaps and non-B DNA cluster hotspots
- Class- and subclass-level normalized confidence scoring on a 1.0–3.0 scale
- Interactive Streamlit workflow for upload, analysis, results review, and export
- Multi-sequence FASTA support, example datasets, and NCBI retrieval workflows
- Structured exports for downstream analysis and reporting

## Supported Non-B DNA Structures

The current implementation reports 11 output classes (9 primary detector families plus Hybrid and Non-B DNA Clusters):

- **Curved DNA** — subclasses: Global Curvature, Local Curvature
- **Slipped DNA** — subclasses: Direct Repeat, STR
- **Cruciform** — subclass: Cruciform forming IRs
- **R-Loop** — subclass: R-loop formation sites
- **Triplex** — subclasses: H-DNA, Sticky DNA
- **G-Quadruplex** — subclasses: Telomeric G4, Stacked G4, Canonical intramolecular G4, Extended-loop canonical, Higher-order G4 array/G4-wire, Intramolecular G-triplex, Two-tetrad weak PQS, Bulged G4
- **i-Motif** — subclasses: Canonical i-motif, Relaxed i-motif, AC-motif
- **Z-DNA** — subclasses: Z-DNA, eGZ
- **A-philic DNA** — subclass: A-philic DNA
- **Hybrid** — subclass: Dynamic overlaps
- **Non-B DNA Clusters** — subclass: Dynamic clusters

## Computational Workflow

1. Input FASTA or retrieve a sequence from NCBI.
2. Process and validate nucleotide sequence records.
3. Detect motif-class candidates with motif-specific algorithms.
4. Apply class-specific scoring and subclass classification.
5. Annotate hybrid overlaps and cluster hotspots.
6. Review results and visualizations.
7. Download/export the resulting tables and summaries.

## Algorithms

- **Curved DNA:** A/T-tract phasing and local tract detection
- **Slipped DNA:** repeat scanning for STR and direct-repeat architectures
- **Cruciform:** seed-and-extend inverted-repeat detection with thermodynamic filtering
- **R-Loop:** QmRLFS-related R-loop initiation/elongation logic
- **Triplex:** mirror-repeat detection plus Sticky DNA repeat handling
- **G-Quadruplex:** G4Hunter-derived scoring with subclass pattern resolution
- **i-Motif:** C-tract and AC-motif pattern logic
- **Z-DNA:** 10-mer propensity scoring plus eGZ repeat detection
- **A-philic DNA:** 10-mer A-form propensity scoring

## Scoring and Confidence

NonBDNAFinder normalizes class-specific raw scores to a **1.0–3.0** confidence range.

```text
Score = 1.0 + 2.0 × (raw − floor) / (ceiling − floor)
```

Normalized confidence values are primarily intended for relative ranking within the relevant structural class or subclass. They should not be interpreted as directly equivalent measures of structural stability across different motif classes.

A-philic DNA predictions represent sequence-derived conformational propensity and should not be interpreted as direct experimental confirmation of A-DNA formation.

## Hybrid Prediction

Hybrid overlap represents overlapping sequence-level prediction criteria and does not by itself demonstrate simultaneous physical formation of multiple DNA structures at the same locus.

For G-quadruplex/i-motif overlap, interpret the result as **co-localized structural potential on complementary strands**, not direct proof that both structures simultaneously exist at the same locus.

## Input Requirements

Supported input workflows in the current repository:

- Upload FASTA or multi-FASTA files (`.fa`, `.fasta`, `.fna`, `.txt`)
- Paste FASTA-formatted sequence records
- Load bundled example data
- Fetch from NCBI by accession, gene, or genome interval

Input records should provide sequence identifiers in FASTA headers and nucleotide sequences composed of supported IUPAC DNA characters. Results are reported with 1-based inclusive coordinates.

## Application Limits

The public web application currently supports:

- **Maximum sequences per upload:** 25
- **Maximum genome/input size:** 10 MB

For larger genomes or datasets, please use the local version of Non-B DNA Finder available through this GitHub repository:

- https://github.com/VRYella/NonBDNAFinder

## Output Formats

Current Streamlit application downloads:

- CSV
- Excel (XLSX)
- JSON
- BED
- GFF3
- PDF visualization package

## How to Use the Web Application

1. Upload FASTA.
2. Select available motif classes/subclasses.
3. Configure the available analysis options.
4. Run analysis.
5. Explore results and visualizations.
6. Download results.

## Running Locally

Clone and install:

```bash
git clone https://github.com/VRYella/NonBDNAFinder.git
cd NonBDNAFinder
pip install -r requirements.txt
```

Launch the Streamlit application locally:

```bash
streamlit run app.py
```

Additional local workflows already present in the repository:

- `NonBDNAFinder_Analysis.ipynb` for notebook-based exploration
- Python API entry points in `Utilities/nonbscanner.py`

## Interpretation and Limitations

NonBDNAFinder provides sequence-based computational predictions. These predictions do not automatically capture:

- DNA supercoiling
- transcriptional state
- chromatin state
- DNA methylation
- histone modifications

unless that contextual information is explicitly incorporated by the current implementation.

Computational prediction is not experimental confirmation.

## References

- Bedrat A, Lacroix L, Mergny JL. Re-evaluation of G-quadruplex propensity with G4Hunter. *Nucleic Acids Research* (2016).
- Jenjaroenpun P, Wongsurawat T, Sutheeworapong S, et al. QmRLFS-finder: a model, web server and stand-alone tool. *Nucleic Acids Research* (2015/2016).
- Frank-Kamenetskii MD, Mirkin SM. Triplex DNA structures. *Annual Review of Biochemistry* (1995).
- Ho PS, Frederick CA, Saal D, et al. Z-DNA propensity framework underlying the implemented 10-mer approach. *Journal of Biomolecular Structure and Dynamics* (1986).
- Zeraati M, Langley DB, et al. i-Motif DNA structures are formed in the nuclei of human cells. *Nature Chemistry* (2018).
- Crothers DM, Koo HS, Olson WK and related DNA curvature studies used for A/T-tract phasing interpretation.
- Vinogradov AE and related conformational propensity studies underlying the implemented A-philic scoring approach.

## Citation

```text
Yella VR (2025). NonBDNAFinder: Comprehensive Detection and Analysis of Non-B DNA Forming Motifs. GitHub repository: https://github.com/VRYella/NonBDNAFinder
```
