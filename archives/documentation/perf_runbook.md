# Performance Runbook for NBDScanner

This document provides guidance on optimizing NBDScanner performance for large-scale genome analysis.

## Overview

NBDScanner supports multiple scanning backends and parallel processing modes to handle genomes of various sizes efficiently.

### Backend Selection

NBDScanner automatically selects the best available backend:

1. **Hyperscan** (fastest) - Intel's high-performance pattern matching library
2. **Numba** (fast) - JIT-compiled Python with NumPy acceleration
3. **Pure Python** (baseline) - Standard regex-based matching

Check available backends:
```python
from scanner_backends import get_best_backend, HYPERSCAN_AVAILABLE, NUMBA_AVAILABLE

print(f"Best backend: {get_best_backend()}")
print(f"Hyperscan available: {HYPERSCAN_AVAILABLE}")
print(f"Numba available: {NUMBA_AVAILABLE}")
```

## Performance Guidelines

### Sequence Size Recommendations

| Sequence Size | Recommended Approach | Expected Time |
|---------------|---------------------|---------------|
| < 10 KB | Standard scanning | < 5 seconds |
| 10 KB - 100 KB | Standard scanning | 5-30 seconds |
| 100 KB - 1 MB | Parallel scanning | 30 seconds - 2 minutes |
| 1 MB - 10 MB | Parallel with chunking | 2-15 minutes |
| > 10 MB | CLI with streaming output | 15+ minutes |

### CLI Usage for Large Genomes

For genomes larger than 1 MB, use the CLI with parallel processing:

```bash
# Basic parallel scan
python cli/nbdscanner_cli.py \
    --input genome.fa \
    --out results.ndjson \
    --workers 8

# Optimized for large genomes (>10 MB)
python cli/nbdscanner_cli.py \
    --input large_genome.fa \
    --out results \
    --format ndjson \
    --workers 16 \
    --chunksize 100000 \
    --streaming
```

### Chunk Size Tuning

The `--chunksize` parameter controls how the sequence is divided for parallel processing:

- **Small chunks (50 KB)**: Better load balancing, more overhead
- **Large chunks (200 KB)**: Less overhead, potential memory issues
- **Recommended**: 100 KB for most systems

```bash
# For systems with limited memory
python cli/nbdscanner_cli.py --input genome.fa --out results --chunksize 50000

# For high-memory systems
python cli/nbdscanner_cli.py --input genome.fa --out results --chunksize 200000
```

### Worker Count

The `--workers` parameter sets the number of parallel processes:

- Default: Number of CPU cores
- For I/O bound systems: 2x CPU cores
- For memory-constrained systems: CPU cores / 2

```bash
# Auto-detect (default)
python cli/nbdscanner_cli.py --input genome.fa --out results

# Explicit worker count
python cli/nbdscanner_cli.py --input genome.fa --out results --workers 4
```

## Memory Management

### Estimated Memory Usage

Memory usage scales approximately as:

```
Total Memory ≈ Sequence Size + (Workers × Chunk Size) + Motif Storage
```

For a 10 MB genome with 8 workers and 100 KB chunks:
```
≈ 10 MB + (8 × 100 KB) + ~5 MB = ~16 MB
```

### Reducing Memory Usage

1. **Use streaming output** - Results are written incrementally
2. **Reduce chunk size** - Smaller chunks use less memory
3. **Reduce workers** - Fewer parallel processes

```bash
# Memory-efficient configuration
python cli/nbdscanner_cli.py \
    --input large_genome.fa \
    --out results \
    --chunksize 50000 \
    --workers 2 \
    --streaming
```

### Shared Memory

NBDScanner uses `multiprocessing.SharedMemory` to share the sequence between workers without copying. This significantly reduces memory usage for parallel processing.

## Profiling and Benchmarking

### Quick Performance Test

```bash
# Run smoke tests
python -m pytest tests/test_perf_smoke.py -v

# Benchmark specific sequence
python -c "
import time
from nonbscanner import analyze_sequence

seq = 'ATGC' * 10000  # 40 KB
start = time.time()
motifs = analyze_sequence(seq, 'benchmark')
elapsed = time.time() - start
print(f'Time: {elapsed:.2f}s, Speed: {len(seq)/elapsed:,.0f} bp/s')
print(f'Motifs found: {len(motifs)}')
"
```

### Detailed Profiling

For detailed performance analysis:

```python
import cProfile
import pstats
from nonbscanner import analyze_sequence

# Profile the analysis
profiler = cProfile.Profile()
profiler.enable()

sequence = 'ATGC' * 25000  # 100 KB
motifs = analyze_sequence(sequence, 'profile_test')

profiler.disable()
stats = pstats.Stats(profiler)
stats.sort_stats('cumtime')
stats.print_stats(20)
```

### Memory Profiling

To profile memory usage:

```bash
pip install memory_profiler

python -m memory_profiler cli/nbdscanner_cli.py \
    --input test.fa \
    --out results \
    --no-parallel
```

## Troubleshooting

### Slow Performance

1. **Check backend availability**
   ```python
   from scanner_backends import get_best_backend
   print(get_best_backend())  # Should be 'hyperscan' or 'numba' for best performance
   ```

2. **Enable parallel processing**
   ```bash
   python cli/nbdscanner_cli.py --input genome.fa --out results --workers 4
   ```

3. **Reduce filtering**
   - Some motif classes are more expensive to detect
   - Use `--classes` to limit detection to specific types

### Memory Errors

1. **Reduce chunk size**
   ```bash
   python cli/nbdscanner_cli.py --input genome.fa --out results --chunksize 25000
   ```

2. **Reduce workers**
   ```bash
   python cli/nbdscanner_cli.py --input genome.fa --out results --workers 2
   ```

3. **Use streaming output**
   ```bash
   python cli/nbdscanner_cli.py --input genome.fa --out results --streaming
   ```

### Parallel Processing Failures

If parallel processing fails:

1. **Fall back to sequential**
   ```bash
   python cli/nbdscanner_cli.py --input genome.fa --out results --no-parallel
   ```

2. **Check for shared memory issues**
   ```bash
   # On Linux, check shared memory limits
   cat /proc/sys/kernel/shmmax
   ```

## Backend-Specific Notes

### Hyperscan

Hyperscan provides the best performance but requires the `hyperscan` package:

```bash
pip install hyperscan
```

Hyperscan works best with:
- Many patterns matching simultaneously
- Large sequences (> 1 MB)
- Pattern-based detection (not algorithmic)

### Numba

Numba provides JIT compilation for Python code:

```bash
pip install numba
```

Numba works best with:
- Numeric operations
- Loop-heavy code
- First run is slower (compilation), subsequent runs are fast

### Pure Python

Always available as fallback. Performance can be improved by:
- Using PyPy instead of CPython
- Pre-compiling regex patterns
- Reducing pattern complexity

## Example: Processing E. coli Genome

```bash
# Download E. coli genome (if needed)
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
# gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz

# Scan with optimal settings for 4.6 Mb genome
time python cli/nbdscanner_cli.py \
    --input ecoli.fa \
    --out ecoli_motifs \
    --format ndjson,csv \
    --workers 8 \
    --chunksize 100000 \
    --verbose

# Expected time: 3-10 minutes depending on system
```

## Performance Targets

| Metric | Target | Notes |
|--------|--------|-------|
| Throughput | > 5,000 bp/s | Per detector, sequential |
| Parallel speedup | ~Nx with N workers | Up to 9 detectors |
| Memory overhead | < 3x sequence size | With parallel processing |
| Startup time | < 2 seconds | Pattern compilation |
