# Performance Benchmarking Report - NonBDNAFinder

## Overview

This document describes the performance benchmarking system for NonBDNAFinder, which generates comprehensive reports showing the time taken by every function of the tool when analyzing sequences.

## Purpose

The performance report generator provides:
- **Individual detector timing** for each of the 9 Non-B DNA detectors
- **Memory usage** monitoring during analysis
- **Throughput calculations** (base pairs per second)
- **Motif detection statistics** including counts and coverage
- **Full pipeline benchmarking** showing end-to-end performance
- **HTML and JSON reports** for easy viewing and integration

## Test Sequence Specification

The benchmark uses a **100,000 nucleotide (100kb)** test sequence that includes:
- Random regions for baseline performance
- AT-rich and GC-rich regions
- G-quadruplex motifs
- Direct and inverted repeats
- Curved DNA (A-tracts)

This ensures realistic performance measurements across all detector types.

## Detector Functions Benchmarked

The tool benchmarks all 9 specialized Non-B DNA detectors:

### 1. **Curved DNA Detector**
   - Detects A-tract mediated DNA curvature
   - Methods: APR phasing analysis, local/global curvature
   - Output: Local and global curvature motifs

### 2. **Slipped DNA Detector**
   - Detects direct repeats and short tandem repeats (STRs)
   - Methods: K-mer indexing, repeat detection
   - Output: STR and direct repeat motifs

### 3. **Cruciform Detector**
   - Detects palindromic inverted repeats
   - Methods: K-mer indexing, inverted repeat detection
   - Output: Cruciform/inverted repeat motifs

### 4. **R-Loop Detector**
   - Detects RNA-DNA hybrid formation sites
   - Methods: QmRLFS algorithm (models 1 and 2)
   - Output: R-loop formation site motifs

### 5. **Triplex Detector**
   - Detects three-stranded DNA structures
   - Methods: Mirror repeat detection, purine/pyrimidine analysis
   - Output: Triplex and sticky DNA motifs

### 6. **G-Quadruplex Detector**
   - Detects four-stranded G-rich structures
   - Methods: G4Hunter scoring, pattern matching
   - Output: 7 G4 subclasses including canonical, telomeric, stacked, weak PQS

### 7. **i-Motif Detector**
   - Detects C-rich intramolecular structures
   - Methods: C-tract analysis, HUR AC-motif detection
   - Output: Canonical i-motif and AC-motif variants

### 8. **Z-DNA Detector**
   - Detects left-handed helix structures
   - Methods: 10-mer scoring table, eGZ pattern matching
   - Output: Z-DNA and eGZ (extruded-G) motifs

### 9. **A-philic DNA Detector**
   - Detects A-rich structural regions
   - Methods: Tetranucleotide pattern analysis
   - Output: A-philic DNA motifs

## Usage

### Basic Usage

Generate a performance report with default settings (100,000 nt sequence):

```bash
python performance_report.py
```

This creates `performance_report.html` with comprehensive timing data.

### Custom Sequence Length

Test with a different sequence length:

```bash
python performance_report.py --length 50000 --output report_50k.html
```

### Using Your Own Sequence

Benchmark with your own FASTA file:

```bash
python performance_report.py --sequence-file my_sequence.fasta --output my_report.html
```

### Generate Both HTML and JSON Reports

```bash
python performance_report.py --output report.html --json report.json
```

### All Options

```bash
python performance_report.py --help

usage: performance_report.py [-h] [--sequence-file SEQUENCE_FILE] [--length LENGTH] 
                             [--output OUTPUT] [--json JSON] [--seed SEED]

options:
  --sequence-file SEQUENCE_FILE  Input FASTA file (optional)
  --length LENGTH                Length of test sequence (default: 100000)
  --output OUTPUT                Output HTML report file (default: performance_report.html)
  --json JSON                    Output JSON report file (optional)
  --seed SEED                    Random seed for reproducibility (default: 42)
```

## Example Output

### Sample Performance Results (100,000 nt sequence)

#### Individual Detector Performance

| Detector       | Time (s) | Throughput (bp/s) | Memory (MB) | Motifs Found |
|----------------|----------|-------------------|-------------|--------------|
| i-Motif        | 0.039    | 2,561,591        | 0.16        | 82           |
| A-philic DNA   | 0.262    | 380,454          | 0.03        | 41           |
| G-Quadruplex   | 0.324    | 307,984          | 0.69        | 361          |
| Z-DNA          | 0.377    | 265,115          | 0.02        | 24           |
| Triplex        | 0.386    | 258,491          | 0.01        | 12           |
| R-Loop         | 0.539    | 185,322          | 0.12        | 34           |
| Curved DNA     | 0.799    | 124,914          | 0.12        | 31           |
| Slipped DNA    | 6.019    | 16,585           | 0.10        | 1            |
| Cruciform      | 6.295    | 15,858           | 0.00        | 0            |

#### Pipeline Performance

- **Total Pipeline Time**: 15.71 seconds
- **Pipeline Throughput**: 6,354 bp/s
- **Total Motifs Detected**: 643 motifs
- **Memory Used**: 1.24 MB
- **Sequence Coverage**: 195.35% (overlapping motifs)

### Performance Insights

1. **Fastest Detectors**:
   - i-Motif: 2.56 million bp/s
   - A-philic DNA: 380k bp/s
   - G-Quadruplex: 308k bp/s

2. **Most Productive Detectors**:
   - G-Quadruplex: 361 motifs (56% of total)
   - i-Motif: 82 motifs (13% of total)
   - A-philic DNA: 41 motifs (6% of total)

3. **Memory Efficiency**:
   - Most detectors use < 0.2 MB
   - G-Quadruplex uses 0.69 MB (highest due to complex scoring)
   - Total pipeline uses only 1.24 MB

## Report Contents

### HTML Report Sections

1. **Summary Statistics**
   - Total detector time
   - Pipeline time
   - Total motifs found
   - Average throughput

2. **Individual Detector Performance Table**
   - Time elapsed
   - Throughput
   - Memory usage
   - Motifs found
   - Coverage percentage
   - Success status

3. **Timing Comparison Chart**
   - Visual bar chart comparing detector times
   - Easy identification of performance bottlenecks

4. **Full Pipeline Performance**
   - End-to-end timing
   - Complete motif statistics
   - Class and subclass distribution

### JSON Report Structure

```json
{
  "timestamp": "2025-12-25T18:29:08.120928",
  "sequence_length": 99827,
  "detector_results": {
    "Curved DNA": {
      "elapsed_time": 0.799,
      "throughput_bps": 124914.3,
      "motif_count": 31,
      "subclasses": { ... }
    },
    ...
  },
  "pipeline_result": {
    "elapsed_time": 15.71,
    "throughput_bps": 6353.7,
    "motif_count": 643
  },
  "summary": { ... }
}
```

## Performance Optimization Tips

Based on the benchmarking results:

1. **For Large Sequences (>100kb)**:
   - Use chunking mode (automatically enabled)
   - Consider parallel processing
   - Monitor memory usage

2. **For Speed-Critical Applications**:
   - Focus on fastest detectors (i-Motif, A-philic, G4)
   - Consider disabling slower detectors if not needed
   - Use Hyperscan acceleration (if available)

3. **For Comprehensive Analysis**:
   - Full pipeline provides best balance
   - Expect ~6,000-8,000 bp/s throughput
   - Plan for ~0.1-0.2 seconds per 1000 bp

## Technical Details

### Benchmarking Methodology

1. **Sequence Generation**: Creates realistic test sequences with varied composition
2. **Timing**: Uses high-resolution `time.time()` for accurate measurements
3. **Memory Tracking**: Uses `tracemalloc` for Python memory profiling
4. **Throughput Calculation**: `throughput = sequence_length / elapsed_time`
5. **Statistics**: Collects motif counts, coverage, and subclass distributions

### System Requirements

- Python 3.8+
- Dependencies: pandas, numpy, biopython
- Recommended: 4GB+ RAM for 100kb sequences
- Disk space: Minimal (reports are typically < 50 KB)

## Interpreting Results

### What to Look For

- **High throughput** (>100,000 bp/s): Efficient detector
- **Low memory** (<1 MB): Memory-efficient implementation
- **Many motifs**: Active detector for this sequence type
- **High coverage**: Detector finding large/overlapping regions

### Common Patterns

- **i-Motif is always fastest**: Simple C-tract pattern matching
- **Slipped/Cruciform are slowest**: Complex k-mer indexing
- **G-Quadruplex finds most motifs**: Multiple subclasses detected
- **Pipeline is slower than sum**: Includes overlap resolution, hybrid detection

## Troubleshooting

### Issue: Script fails with ImportError

**Solution**: Install required dependencies:
```bash
pip install pandas numpy biopython matplotlib seaborn
```

### Issue: Out of memory on large sequences

**Solution**: Reduce sequence length or use chunking:
```bash
python performance_report.py --length 50000
```

### Issue: Report generation is slow

**Solution**: This is expected. The tool is timing the actual analysis!
- 100kb sequence: ~15-20 seconds
- 50kb sequence: ~5-8 seconds
- 10kb sequence: ~1-2 seconds

## Citations

If you use this benchmarking tool in your research, please cite:

```bibtex
@software{nonbdnafinder2025,
  author = {Yella, Venkata Rajesh},
  title = {NonBDNAFinder 2025.1: Nobel-Level Quality DNA Motif Detection},
  year = {2025},
  url = {https://github.com/VRYella/NonBDNAFinder},
  version = {2025.1}
}
```

## Support

For issues or questions:
- GitHub: https://github.com/VRYella/NonBDNAFinder
- Email: yvrajesh_bt@kluniversity.in

## License

MIT License - See LICENSE file for details

---

**Version**: 2025.1  
**Author**: Dr. Venkata Rajesh Yella  
**Institution**: KL University  
**Last Updated**: December 2025
