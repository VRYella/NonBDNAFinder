# Z-DNA Extractor & Annotator Implementation

## Overview

This implementation adds a sophisticated Z-DNA extractor & annotator (`zdna_hs.py`) that uses Intel Hyperscan for high-performance transition detection, following the Z-seeker algorithm methodology.

## Key Features

- **Hyperscan Acceleration**: Uses Intel Hyperscan for fast dinucleotide transition detection (GC, CG, GT, TG, AC, CA, AT, TA)
- **Z-seeker Algorithm**: Implements the Z-DNA seeker scoring methodology with transition-based scoring
- **Sophisticated Scoring**: Includes weighted transitions, consecutive AT bonuses, mismatch penalties, and cadence rewards
- **Kadane-like Detection**: Uses modified Kadane algorithm for optimal subarray detection with drop thresholds
- **Multi-processing**: Parallel processing with automatic task distribution
- **Fallback Support**: Pure Python implementation when Hyperscan is not available

## Usage

### Standalone Script

```bash
# Basic usage
python zdna_hs.py --fasta sequences.fa --threshold 5.0

# Advanced parameters
python zdna_hs.py --fasta sequences.fa \
    --threshold 10.0 \
    --GC_weight 3.0 \
    --AT_weight 1.0 \
    --GT_weight 2.0 \
    --AC_weight 2.0 \
    --cadence_reward 0.2 \
    --mismatch_penalty_type linear \
    --drop_threshold 50.0 \
    --method transitions \
    --n_jobs 4 \
    --output_dir zdna_results
```

### Integration with Existing Framework

The Z-DNA detector is automatically enhanced when using the main NBDFinder pipeline:

```python
from motif_detectors import ZDNADetector

# Use enhanced Hyperscan-accelerated detection
detector = ZDNADetector(use_hyperscan=True)
candidates = detector.detect(sequence, "chr1", "seq1", 0)
scored_candidates = detector.score(candidates)
```

### CLI Entry Point

When installed, the script is available as a command-line tool:

```bash
zdna-hs --fasta sequences.fa --threshold 5.0
```

## Algorithm Details

### Transition-Based Scoring

The algorithm scores dinucleotide transitions with the following weights:
- GC/CG transitions: 3.0 (strongest Z-DNA formers)
- GT/TG transitions: 2.0 
- AC/CA transitions: 2.0
- AT/TA transitions: 1.0 + consecutive bonuses

### Mismatch Penalties

Non-Z-DNA forming transitions receive penalties:
- Linear: `-3 + -1 * (mismatch_count - 1)`
- Exponential: `-3^mismatch_count` (capped at -32000)

### Consecutive AT Scoring

AT transitions receive additional bonuses when consecutive:
- First AT: +3.0
- Second AT: +1.5  
- Third+ AT: +0.7

### Subarray Detection

Uses a modified Kadane algorithm to find optimal Z-DNA regions:
- Minimum threshold for candidate regions
- Drop threshold for terminating regions when score drops significantly
- Reports (start, end, score, sequence) tuples

## Output Format

CSV output with columns:
- `Chromosome`: Sequence identifier
- `Start`: 1-based start position
- `End`: 1-based end position  
- `Z-DNA Score`: Calculated Z-DNA propensity score
- `Sequence`: DNA sequence of the region
- `totalSequenceScore`: Sum of all transition scores for the full sequence

## Performance

- **Hyperscan Mode**: Highly optimized for large sequences with Intel Hyperscan
- **Fallback Mode**: Pure Python implementation for compatibility
- **Multi-processing**: Automatic parallelization for large datasets
- **Memory Efficient**: Streaming processing for large genome files

## Requirements

- Python 3.8+
- numpy, pandas, biopython, termcolor
- hyperscan (optional, for acceleration)
- Intel Hyperscan C library (optional)

## Scientific References

Based on the Z-DNA seeker algorithm methodology:
- Ho et al. EMBO J 5(10):2737-2744 (1986)
- Z-DNA formation in alternating purine-pyrimidine sequences