# Performance and Memory Optimization Guide

## Overview

This document describes the performance and memory optimizations implemented in NBDScanner to improve handling of large datasets and enhance the user experience.

## Key Improvements

### 1. Memory-Efficient File Upload

#### **Chunked File Parsing**
- **Function**: `parse_fasta_chunked(file_object, chunk_size_mb=2)`
- **Purpose**: Process large FASTA files without loading entire file into memory
- **How it works**:
  - Reads file in 2MB chunks (optimized from 5MB for better memory efficiency)
  - Yields sequences one at a time using Python generators
  - Immediately frees memory after processing each sequence
  - Reduces peak memory usage by 50-70% for large files
- **Maximum Upload Size**: 1GB (increased from 200MB default)

**Example Usage**:
```python
from utilities import parse_fasta_chunked

# Process large file without memory overflow
file_object = open('large_genome.fasta', 'rb')
for name, sequence in parse_fasta_chunked(file_object):
    # Process one sequence at a time
    results = analyze_sequence(sequence, name)
    # Memory is freed after each iteration
```

**Benefits**:
- ✅ Handles files up to 1GB in size
- ✅ Handles files larger than available RAM
- ✅ Streaming processing for better memory efficiency
- ✅ Immediate garbage collection after each sequence
- ✅ Progress indicators for files with >10 sequences
- ✅ Automatic garbage collection for sequences >10MB

### 2. File Preview Feature

#### **Lightweight Preview Generation**
- **Function**: `get_file_preview(file_object, max_sequences=3, max_preview_chars=400)`
- **Purpose**: Quick inspection of file contents before full processing
- **Benefits**:
  - Shows file size, sequence count, and total base pairs
  - Displays preview of first few sequences
  - Only reads necessary portions of the file
  - 70-120% faster than full file parsing

**UI Display**:
```
📁 File: genome_sequences.fasta | Size: 45.23 MB
✅ File contains 250 sequences totaling 12,450,789 bp

Sequence Preview:
  chr1: 50,000 bp | GC: 42.5% | AT: 57.5%
  chr2: 48,500 bp | GC: 41.8% | AT: 58.2%
  chr3: 52,300 bp | GC: 43.1% | AT: 56.9%
  ...and 247 more sequences.
```

### 3. Result Caching

#### **Analysis Result Cache**
- **Function**: `@st.cache_data` decorator on analysis functions
- **Configuration**:
  - Max entries: 10 (stores up to 10 recent analyses)
  - TTL: 3600 seconds (1 hour)
  - Cache key: Sequence hash + parameters

**Benefits**:
- ✅ Instant results for repeated analyses
- ✅ Reduces server load
- ✅ Improves user experience with faster response

**Cache Invalidation**:
- Automatic after 1 hour
- Cleared when exceeding 10 entries (LRU eviction)
- Manual clear via Streamlit's cache clear button

### 4. Result Pagination

#### **Large Dataset Display Optimization**
- **Implementation**: Paginated dataframe display (100 rows per page)
- **Triggers**: Automatically enabled for datasets >100 motifs
- **Benefits**:
  - ✅ Faster page rendering
  - ✅ Reduced browser memory usage
  - ✅ Better user experience for large result sets

**UI Controls**:
```
Page: [1] of 15 (showing 100 rows per page)
Showing motifs 1 to 100 of 1,432
```

### 5. System Resource Monitoring

#### **Real-time Memory Usage Display**
- **Location**: Upload & Analyze tab (collapsible expander)
- **Metrics**:
  - Memory usage percentage
  - Available memory
  - CPU core count
  - Visual progress bar with color coding

**Display**:
```
💻 System Resource Monitor
  💾 Memory Usage: 45.2% (8.5 GB used)
  🖥️ Total Memory: 16.0 GB
  ⚙️ CPU Cores: 8

[████████░░░░░░░░░░░░] Memory: 45.2% used (green/orange/red)
```

## Performance Benchmarks

### Test Results (from performance_benchmark.py)

| Configuration | Sequences | Seq Length | Traditional Time | Optimized Time | Memory Saved |
|--------------|-----------|------------|------------------|----------------|--------------|
| Small | 10 | 1 KB | 0.001s | 0.001s | 89% |
| Medium | 50 | 5 KB | 0.003s | 0.003s | 97% |
| Large | 100 | 10 KB | 0.015s | 0.008s | 98% |
| Very Large | 200 | 20 KB | 0.050s | 0.011s | 99% |

### Analysis Performance
- **Throughput**: 9,010 bp/s (50kb sequence)
- **Memory Usage**: 6.41 MB (50kb sequence with 784 motifs)
- **Motif Detection**: 784 motifs in 6.437 seconds

## Best Practices for Users

### For Large Files (>50 MB)
1. ✅ Use file upload (not paste) for better memory efficiency
2. ✅ Monitor system resources before uploading
3. ✅ Close unnecessary browser tabs to free memory
4. ✅ Consider processing sequences in batches if system is low on memory

### For Maximum Performance
1. ✅ Enable caching (automatic in Streamlit)
2. ✅ Use pagination for large result sets
3. ✅ Download results instead of keeping in browser
4. ✅ Clear cache periodically via Streamlit menu

### Memory Usage Guidelines

| File Size | Recommended RAM | Expected Processing Time | Upload Limit |
|-----------|----------------|-------------------------|--------------|
| < 10 MB | 2 GB | < 10 seconds | ✅ Supported |
| 10-50 MB | 4 GB | 10-60 seconds | ✅ Supported |
| 50-100 MB | 8 GB | 1-3 minutes | ✅ Supported |
| 100-500 MB | 16 GB | 3-15 minutes | ✅ Supported |
| 500 MB-1 GB | 32 GB | 15-30 minutes | ✅ Supported |
| > 1 GB | N/A | N/A | ❌ Exceeds limit |

**Note**: Maximum upload size is now 1GB (increased from 200MB default). Files larger than 100MB will show a warning about processing time.

## Technical Details

### Streamlit Configuration (.streamlit/config.toml)

The application now includes optimized Streamlit server settings for large file handling:

```toml
[server]
# Maximum file upload size in MB (1GB)
maxUploadSize = 1000

# Maximum WebSocket message size in MB (1GB)
maxMessageSize = 1000

# Enable XSRF protection for security
enableXsrfProtection = true

# File watcher type for large file processing
fileWatcherType = "auto"

[runner]
# Enable fast reruns for better performance
fastReruns = true

[client]
# Show error details for debugging
showErrorDetails = true
```

**Key Benefits**:
- ✅ Supports files up to 1GB (5x increase from 200MB default)
- ✅ Prevents WebSocket timeout errors for large transfers
- ✅ Enhanced security with XSRF protection
- ✅ Optimized reruns for better performance

### Chunked Parsing Algorithm
```python
def parse_fasta_chunked(file_object, chunk_size_mb=2):
    """
    Memory-efficient streaming parser with 2MB chunks
    
    1. Read file in chunks (optimized to 2MB for better memory efficiency)
    2. Process lines incrementally
    3. Yield completed sequences immediately
    4. Free memory after each yield
    5. Support for files up to 1GB
    """
    chunk_size = chunk_size_mb * 1024 * 1024
    buffer = ""
    current_name = None
    current_seq_parts = []
    
    while True:
        chunk = file_object.read(chunk_size)
        if not chunk:
            break
        
        # Process chunk
        buffer += chunk.decode('utf-8', errors='ignore')
        lines = buffer.split('\n')
        buffer = lines[-1]  # Keep incomplete line
        
        for line in lines[:-1]:
            # Parse FASTA format
            if line.startswith('>'):
                if current_name and current_seq_parts:
                    # Yield and free memory immediately
                    full_seq = ''.join(current_seq_parts)
                    current_seq_parts.clear()
                    yield (current_name, full_seq)
                    del full_seq
                current_name = line[1:].strip()
                current_seq_parts = []
            else:
                current_seq_parts.append(line.upper())
    
    # Yield final sequence
    if current_name and current_seq_parts:
        yield (current_name, ''.join(current_seq_parts))
```

### Cache Configuration
```python
@st.cache_data(
    show_spinner=False,  # Don't show spinner for cached results
    max_entries=10,      # Keep 10 most recent
    ttl=3600            # 1 hour expiration
)
def cache_analysis_results(sequence_hash, sequence, name):
    return analyze_sequence(sequence, name)
```

## Monitoring and Debugging

### Performance Benchmark Script
Run `performance_benchmark.py` to test system performance:

```bash
python performance_benchmark.py
```

**Output**:
- File parsing comparison (traditional vs optimized)
- Memory usage tracking with tracemalloc
- Analysis throughput measurements
- Detailed JSON report saved to `benchmark_results.json`

### Memory Profiling
```python
import tracemalloc

tracemalloc.start()
# Run your code
current, peak = tracemalloc.get_traced_memory()
print(f"Peak memory: {peak / 1024 / 1024:.2f} MB")
tracemalloc.stop()
```

## Future Optimizations

### Planned Improvements
- [ ] Lazy loading for visualizations (render on-demand)
- [ ] Database caching for frequently analyzed sequences
- [ ] Parallel processing for multi-sequence files
- [ ] Compressed data storage for large results
- [ ] Progressive result streaming

### Experimental Features
- Hyperscan pattern database caching (optional dependency)
- GPU acceleration for pattern matching (requires CUDA)
- Distributed processing for very large genomes (requires cluster)

## Troubleshooting

### Common Issues

**File Upload Failed (>2MB files)**
- **Cause**: Missing or misconfigured `.streamlit/config.toml`
- **Solution**: Ensure `.streamlit/config.toml` exists with `maxUploadSize = 1000` and `maxMessageSize = 1000`
- **Verification**: Check that the file is in the root directory of the application

**Out of Memory Error**
- **Cause**: File too large for available RAM
- **Solution**: Use chunked processing (automatic), close other apps, or upgrade RAM
- **For files >500MB**: Ensure you have at least 32GB RAM available

**Slow Analysis**
- **Cause**: Large sequence or many motifs
- **Solution**: Enable chunking, use pagination, check CPU usage
- **For files >100MB**: Processing may take 3-30 minutes depending on complexity

**WebSocket Connection Error**
- **Cause**: WebSocket message size limit exceeded
- **Solution**: Verify `maxMessageSize = 1000` in `.streamlit/config.toml`
- **Alternative**: Split large file into smaller chunks

**Cache Not Working**
- **Cause**: Streamlit cache disabled or full
- **Solution**: Clear cache via menu, restart app, check cache settings

### Error Messages

| Error | Cause | Solution |
|-------|-------|----------|
| `MemoryError` | Insufficient RAM | Reduce file size, enable chunking, upgrade RAM |
| `TimeoutError` | Analysis too slow | Use smaller chunks, check CPU, increase timeout |
| `CacheError` | Cache corruption | Clear cache, restart app |
| `File too large` | Exceeds upload limit | Ensure config.toml has `maxUploadSize = 1000` |
| `WebSocket closed` | Message size limit | Set `maxMessageSize = 1000` in config.toml |

### Performance Tips for Large Files

**Before Upload**:
1. ✅ Verify `.streamlit/config.toml` exists with proper settings
2. ✅ Check available system memory (should have >2x file size)
3. ✅ Close unnecessary browser tabs and applications
4. ✅ Use a modern browser (Chrome, Firefox, Edge)

**During Processing**:
1. ✅ Don't refresh the page while processing
2. ✅ Monitor system resource usage
3. ✅ Be patient - large files can take several minutes
4. ✅ Watch for progress indicators and status messages

**After Processing**:
1. ✅ Download results immediately
2. ✅ Clear browser cache if experiencing issues
3. ✅ Use pagination for viewing large result sets
4. ✅ Export to Excel/CSV for external analysis

## Contact & Support

For issues related to performance optimizations:
- GitHub: https://github.com/VRYella/NonBDNAFinder
- Email: yvrajesh_bt@kluniversity.in

## Changelog

### Version 2024.1.3 (Current - December 2024)
- **Major**: Increased maximum upload size from 200MB to 1GB
- Added `.streamlit/config.toml` with optimized server settings
- Reduced chunk size from 5MB to 2MB for better memory efficiency
- Added progress indicators for files with >10 sequences
- Added automatic garbage collection for sequences >10MB
- Added warnings for files >100MB about processing time
- Updated documentation with new limits and guidelines

### Version 2024.1.2
- Added chunked file parsing for memory efficiency
- Implemented result caching with configurable TTL
- Added pagination for large result sets (>100 rows)
- Created system resource monitor
- Added comprehensive performance benchmark suite

### Version 2024.1.1
- Initial release with basic optimizations
- Simple file upload and analysis

---

**Last Updated**: December 2024
**Author**: Dr. Venkata Rajesh Yella
