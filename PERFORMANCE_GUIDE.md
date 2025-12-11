# NonBDNAFinder Performance Optimization Guide

## 🚀 Performance Overview

NonBDNAFinder has been optimized for maximum speed while maintaining identical output and reproducibility. This guide explains the optimizations and how to get the best performance.

## 📊 Performance Metrics

### Scanner Functions (Optimized)
- **Direct Repeats**: 2.8M bp/s (5-10x faster with Numba)
- **Inverted Repeats**: 2.9M bp/s (5-10x faster with Numba)
- **STRs**: 272K bp/s (optimized inner loops)
- **Pattern Matching**: 1.42x faster with pre-compiled patterns

### End-to-End Analysis
| Sequence Size | Throughput | Notes |
|---------------|------------|-------|
| Small (1kb)   | 64K bp/s   | Fast mode, all detectors |
| Medium (10kb) | 28K bp/s   | Mixed content, 9 detectors |
| Large (100kb) | 3.1K bp/s  | Genomic-like, comprehensive |

### Speedup Factors
- **Numba JIT**: 5-10x for k-mer operations
- **Pattern Cache**: 1.4x for regex matching
- **Parallel Mode**: 1.2-1.3x on multi-core systems
- **Overall**: 2-5x faster than baseline

## 🔧 Installation for Maximum Performance

### Basic Installation
```bash
pip install -r requirements.txt
```

### Performance Enhancements (Recommended)
```bash
# Install optional performance libraries
pip install -r requirements-performance.txt

# This includes:
# - numba>=0.58.0  (JIT compilation, 5-10x speedup)
# - regex>=2023.10 (faster pattern matching, 1.4x speedup)
```

### Verify Installation
```bash
python -c "import numba, regex; print('Performance libs installed!')"
```

## 🎯 Usage for Best Performance

### Automatic Optimization (Default)
```python
import nonbscanner as nbs

# Fast mode is now enabled by default
motifs = nbs.analyze_sequence(sequence, "my_sequence")
```

### Explicit Fast Mode
```python
# Explicitly enable parallel processing
motifs = nbs.analyze_sequence(sequence, "my_sequence", use_fast_mode=True)
```

### Large Sequences (>10kb)
```python
# Automatic chunking for sequences >10kb
motifs = nbs.analyze_sequence(large_sequence, "genome", 
                               use_chunking=True,
                               chunk_size=10000,
                               chunk_overlap=500)
```

### Batch Processing
```python
from concurrent.futures import ThreadPoolExecutor

sequences = {"seq1": seq1, "seq2": seq2, "seq3": seq3}
results = {}

with ThreadPoolExecutor(max_workers=4) as executor:
    futures = {
        executor.submit(nbs.analyze_sequence, seq, name): name
        for name, seq in sequences.items()
    }
    
    for future in as_completed(futures):
        name = futures[future]
        results[name] = future.result()
```

## 🎛️ Optimization Techniques Used

### 1. Numba JIT Compilation
**What**: Just-in-time compilation of Python code to machine code
**Where**: K-mer indexing, sequence comparison, GC calculation
**Benefit**: 5-10x speedup for numerical operations
**Fallback**: Pure Python if Numba not installed

```python
# Automatically used when available
from scanner import find_direct_repeats

# Uses Numba-optimized version if installed
repeats = find_direct_repeats(sequence)
```

### 2. Pre-compiled Pattern Cache
**What**: Regex patterns compiled once at module load
**Where**: All pattern-based detectors (G4, i-motif, Z-DNA, etc.)
**Benefit**: 1.4x speedup, eliminates compilation overhead
**Fallback**: Standard `re` module if `regex` not installed

```python
# Patterns are pre-compiled automatically
from pattern_cache import get_pattern_cache

cache = get_pattern_cache()
matches = cache.find_all('g4', 'canonical', sequence)
```

### 3. Parallel Processing
**What**: Run all 9 detectors simultaneously
**Where**: Main analysis pipeline
**Benefit**: 1.2-1.3x wall-clock speedup on multi-core systems
**When**: Enabled by default (use_fast_mode=True)

```python
# Parallel mode (default)
motifs = nbs.analyze_sequence(sequence, "test", use_fast_mode=True)

# Sequential mode (for debugging)
motifs = nbs.analyze_sequence(sequence, "test", use_fast_mode=False)
```

### 4. Chunked Processing
**What**: Divide large sequences into manageable chunks
**Where**: Sequences >10kb
**Benefit**: Better memory efficiency, enables parallelism
**When**: Auto-enabled for large sequences

```python
# Automatic for sequences >10kb
motifs = nbs.analyze_sequence(large_seq, "genome")

# Manual control
motifs = nbs.analyze_sequence(large_seq, "genome",
                               use_chunking=True,
                               chunk_size=50000,
                               chunk_overlap=500)
```

## 📈 Benchmarking Your System

### Run Comprehensive Benchmark
```bash
python benchmark_performance.py
```

This will:
- Test individual scanner functions
- Measure pattern cache performance
- Run end-to-end analysis on test sequences
- Generate performance report
- Save results to `benchmark_results.json`

### Expected Output
```
======================================================================
FULL ANALYSIS BENCHMARK
======================================================================

--- 1kb_g4 (1050 bp) ---
Standard Mode (Sequential):
  Found: 2 motifs
  Time: 0.0220s
  Speed: 47,635 bp/s
Fast Mode (Parallel):
  Found: 2 motifs
  Time: 0.0165s
  Speed: 63,760 bp/s

  Speedup: 1.34x
  Reproducible: ✓ Yes
```

### Custom Benchmark
```python
import time
import nonbscanner as nbs

sequence = load_your_sequence()

# Measure performance
start = time.time()
motifs = nbs.analyze_sequence(sequence, "test")
elapsed = time.time() - start

print(f"Found {len(motifs)} motifs")
print(f"Time: {elapsed:.4f}s")
print(f"Throughput: {len(sequence)/elapsed:,.0f} bp/s")
```

## 🎓 Performance Best Practices

### 1. **Install Performance Libraries**
Always install `numba` and `regex` for maximum speed:
```bash
pip install numba regex
```

### 2. **Use Appropriate Sequence Sizes**
- Small sequences (<1kb): Fast mode overhead minimal
- Medium sequences (1-100kb): Optimal for all optimizations
- Large sequences (>100kb): Use chunking for best results

### 3. **Parallel Processing**
- Enable fast mode (default) for multi-core systems
- Use batch processing for multiple sequences
- Consider parallel chunks for very large sequences

### 4. **Memory Management**
- Use chunking for sequences >100kb to avoid memory pressure
- Process large batches in smaller groups
- Clear results periodically in long-running analyses

### 5. **Profile Your Workload**
```python
# Enable timing information
motifs, progress = nbs.analyze_with_progress(
    sequence, "test", 
    return_progress=True,
    print_progress=True
)

print(progress.format_detector_table())
print(f"Throughput: {progress.get_throughput():,.0f} bp/s")
```

## 🔍 Troubleshooting Performance Issues

### Slow Performance
1. **Check if performance libraries are installed:**
   ```python
   try:
       import numba
       print("Numba: ✓ Installed")
   except ImportError:
       print("Numba: ✗ Not installed - install for 5-10x speedup")
   
   try:
       import regex
       print("Regex: ✓ Installed")
   except ImportError:
       print("Regex: ✗ Not installed - install for 1.4x speedup")
   ```

2. **Verify fast mode is enabled:**
   ```python
   # Should use parallel processing by default
   motifs = nbs.analyze_sequence(sequence, "test")
   ```

3. **Check system resources:**
   ```python
   import multiprocessing as mp
   print(f"CPU cores: {mp.cpu_count()}")
   ```

### High Memory Usage
1. **Enable chunking for large sequences:**
   ```python
   motifs = nbs.analyze_sequence(large_seq, "genome",
                                  use_chunking=True,
                                  chunk_size=50000)
   ```

2. **Process in batches:**
   ```python
   for i in range(0, len(sequences), batch_size):
       batch = sequences[i:i+batch_size]
       results.extend(process_batch(batch))
   ```

### Reproducibility Issues
All optimizations maintain identical output. If you see differences:
1. **Check detector versions** - ensure same code
2. **Verify sequence format** - identical input sequences
3. **Compare motif counts** - should be identical

## 📚 Additional Resources

- **Benchmark Script**: `benchmark_performance.py`
- **Pattern Cache**: `pattern_cache.py`
- **Optimized Scanner**: `scanner_optimized.py`
- **Main API**: `nonbscanner.py`

## 🤝 Contributing Performance Improvements

To contribute performance optimizations:
1. Run benchmarks before changes: `python benchmark_performance.py`
2. Implement optimization
3. Run benchmarks after changes
4. Verify reproducibility (output must be identical)
5. Document performance gains
6. Submit PR with benchmark results

## 📊 Performance Tracking

Track performance over time:
```bash
# Run benchmarks
python benchmark_performance.py

# Results saved to benchmark_results.json
cat benchmark_results.json
```

Compare with previous runs to detect regressions.

## ⚡ Summary

**Quick Start for Maximum Performance:**
```bash
# Install with performance libraries
pip install -r requirements.txt -r requirements-performance.txt

# Use default (fast) mode
python your_analysis.py

# Benchmark your system
python benchmark_performance.py
```

**Expected Performance:**
- Small sequences: 50-60K bp/s
- Medium sequences: 25-30K bp/s
- Large sequences: 3-5K bp/s
- K-mer operations: 2-3M bp/s

**Key Optimizations:**
- ✅ Numba JIT (5-10x speedup)
- ✅ Pattern Cache (1.4x speedup)
- ✅ Parallel Mode (1.2-1.3x speedup)
- ✅ Chunked Processing (memory efficient)
- ✅ Automatic fallbacks (always works)

For questions or issues, see the main README or open an issue on GitHub.
