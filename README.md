# NonBDNAFinder

[![License](https://img.shields.io/badge/license-MIT-green.svg)](./LICENSE)

## Overview

NonBDNAFinder detects non-B DNA-forming sequence motifs from DNA sequence input. It provides a Streamlit web interface, a Python API, example datasets, and an analysis notebook.

The current release supports:

- 9 primary motif classes
- 11 reported output classes after hybrid and cluster annotation
- 24 subclasses
- Multiple input modes, including targeted Genome Interval retrieval from NCBI
- Export to CSV, XLSX, BED, GFF3, JSON, and PDF

## Motif Classes

- Curved DNA
- Slipped DNA
- Cruciform
- R-Loop
- Triplex DNA
- G-Quadruplex
- i-Motif
- Z-DNA
- A-philic DNA
- Hybrid regions
- Non-B DNA clusters

## Installation

```bash
git clone https://github.com/VRYella/NonBDNAFinder.git
cd NonBDNAFinder
pip install -r requirements.txt
```

Run the web app:

```bash
python -m streamlit run app.py
```

Optional packages for faster execution:

```bash
pip install numba cython
```

## Input Modes

| Mode | Description |
|------|-------------|
| Upload FASTA | Upload a `.fa`, `.fasta`, `.fna`, or `.txt` FASTA file |
| Paste FASTA | Paste one or more FASTA records |
| Example Data | Load bundled example sequences |
| NCBI Accession | Fetch by accession |
| NCBI Gene | Search and fetch by gene |
| Genome Interval | Fetch `ACCESSION:START-END` directly from NCBI |

### Genome Interval

Genome Interval mode retrieves only the requested locus instead of the full reference sequence.

Compact format:

```text
NC_000913.3:100000-150000
```

Structured fields:

```text
Accession: NC_000913.3
Start: 100000
End: 150000
```

Coordinates are 1-based and inclusive throughout the application and exports.

## Input Requirements

- Sequences shorter than 10 bp are skipped with a warning.
- FASTA input must contain valid IUPAC nucleotide characters.
- RNA input is accepted; `U` is converted to `T`.
- In the hosted web app, input is limited to 5 MB per run for uploaded FASTA, combined pasted sequence content, and Genome Interval fetch length.

## Output

Core output columns:

- `Sequence_Name`
- `Class`
- `Subclass`
- `Start`
- `End`
- `Length`
- `Strand`
- `Score`
- `Detection_Method`

Additional fields for Genome Interval input:

- `Organism`
- `Accession`
- `Chromosome`
- `Interval_Start`
- `Interval_End`
- `Relative_Start`
- `Relative_End`
- `Absolute_Start`
- `Absolute_End`

## Programmatic Usage

```python
from Utilities.nonbscanner import analyze_sequence

sequence = "AGGGGGGGGGCCCCCCCCCTAGGGGGGGGG"
results = analyze_sequence(sequence, "example")
print(len(results))
```

Genome Interval workflow:

```python
from Utilities.genome_interval import parse_interval_string, fetch_genome_interval
from Utilities.nonbscanner import analyze_sequence

interval = parse_interval_string("NC_000913.3:100000-150000")
sequence, record_id = fetch_genome_interval(interval, email="your@email.com")
results = analyze_sequence(sequence, record_id)
annotated = [interval.annotate_motif(motif) for motif in results]
```

Multi-FASTA:

```python
from Utilities.nonbscanner import analyze_fasta_parallel

results = analyze_fasta_parallel("genome.fasta")
```

## Examples

The `examples/` directory includes:

- `examples/example_single.fasta`
- `examples/example_multi.fasta`
- `examples/example.fasta`

These can also be loaded from the web interface.

## Tests

```bash
python -m pytest -q
```

## Notebook

`NonBDNAFinder_Analysis.ipynb` provides a notebook-based workflow for exploratory analysis and visualization.

## Performance

Runtime scales with sequence length. For large sequences, the application uses chunked processing and optional parallel execution.

| Sequence Size | RAM | Typical Runtime |
|--------------|-----|-----------------|
| ≤ 100 Kbp | ~200 MB | < 10 s |
| 100 Kbp – 1 Mbp | ~500 MB | 10–60 s |
| 1 Mbp – 10 Mbp | ~1 GB | 1–10 min |
| 10 Mbp – 100 Mbp | ~2 GB | 10–90 min |
| > 100 Mbp | ~4 GB+ | 90 min+ |

## Project Structure

| Path | Purpose |
|------|---------|
| `app.py` | Streamlit entry point |
| `UI/` | Web interface components |
| `Detectors/` | Motif detector implementations |
| `Utilities/` | Shared analysis and export logic |
| `examples/` | Example datasets |
| `NonBDNAFinder_Analysis.ipynb` | Notebook workflow |

## Citation

If you use NonBDNAFinder in your research, cite:

```bibtex
@article{NBDFinder2026,
  author  = {Yella, Venkata Rajesh and colleagues},
  title   = {Non-B DNA Finder: A unified framework for detection of diverse non-canonical DNA structures},
  journal = {Nucleic Acids Research},
  year    = {2026}
}
```

## License

MIT License. See [LICENSE](./LICENSE).
