# Production-Ready Hyperscan + Streamlit Integration Guide

## 🚀 Overview

This guide provides production-ready patterns and copy-paste code snippets for integrating Hyperscan with Streamlit apps to scan large genomic sequences efficiently without blocking the UI.

## 📋 Key Principles

### 1. **Database Compilation Strategy**
- ✅ Compile DB once per process (expensive operation)
- ✅ Cache with `st.cache_resource` in main app
- ✅ Compile once per worker via `initializer` in ProcessPoolExecutor
- ❌ Never pass Database objects to child processes (not picklable)

### 2. **Hyperscan Callback Optimization**
- ✅ Minimal work in callbacks: only append `(pattern_id, start, end)`
- ✅ Heavy validation/scoring in Python after collecting matches
- ❌ Don't do complex processing inside Hyperscan callbacks

### 3. **Parallelization Strategy**
- ✅ Parallel processing across independent units (contigs/chunks)
- ✅ Hyperscan is C-fast but single-threaded; parallelism gives linear speedups
- ✅ One `hs.Scratch(db)` per thread/process (not shareable)

### 4. **Memory Management**
- ✅ Stream results incrementally (write to file/Parquet in append mode)
- ✅ Update `st.progress` and UI regularly
- ❌ Don't hold millions of hits in RAM

### 5. **Chunking & Overlap Handling**
- ✅ Use `OVERLAP = max_pattern_length - 1` for safety
- ✅ Process chunks but deduplicate matches in overlap regions
- ✅ Alternative: Use `hs.Stream()` for seamless streaming

## 🏗️ Architecture Components

```python
# Core components
HyperscanStreamlitScanner  # Main scanner class
ScanConfig                 # Configuration dataclass
ScanResult                 # Result container
_worker_init()            # ProcessPoolExecutor initializer
_worker_scan_chunk()      # Worker function
```

## 📖 Copy-Paste Examples

### 1. Basic Streamlit App Structure

```python
import streamlit as st
import tempfile
from production_hyperscan_streamlit import (
    HyperscanStreamlitScanner, 
    ScanConfig, 
    get_cached_scanner,
    streamlit_scan_with_progress
)

# App configuration
st.set_page_config(
    page_title="Genomic Pattern Scanner",
    layout="wide",
    page_icon="🧬"
)

st.title("🧬 High-Performance Genomic Scanner")

# Pattern definitions
GENOMIC_PATTERNS = [
    r'G{3,}N{1,7}G{3,}N{1,7}G{3,}N{1,7}G{3,}',  # G-quadruplex
    r'[AT]{10,}',                                 # AT-rich regions
    r'[CG]{20,}',                                # CpG islands
    r'CA{6,}|TG{6,}',                           # Microsatellites
    r'[ACGT]{10,20}(?=[ACGT]\1)',               # Tandem repeats
]

# Sidebar configuration
st.sidebar.header("⚙️ Configuration")
chunk_size = st.sidebar.selectbox(
    "Chunk Size", 
    options=[1, 2, 5, 10], 
    index=1, 
    help="Size in MB for sequence chunking"
) * 1_000_000

max_workers = st.sidebar.slider(
    "Max Workers", 
    min_value=1, 
    max_value=8, 
    value=4,
    help="Number of parallel workers"
)

overlap = st.sidebar.slider(
    "Overlap (bp)", 
    min_value=100, 
    max_value=2000, 
    value=500,
    help="Overlap between chunks to avoid missing patterns"
)

# Pattern selection
st.header("🎯 Select Patterns")
selected_patterns = st.multiselect(
    "Choose patterns to search:",
    options=GENOMIC_PATTERNS,
    default=GENOMIC_PATTERNS[:2],
    help="Select one or more patterns to search for"
)

# Custom pattern input
custom_pattern = st.text_input(
    "Add custom pattern (regex):",
    placeholder="e.g., GGGG|CCCC",
    help="Enter a custom regex pattern"
)

if custom_pattern:
    selected_patterns.append(custom_pattern)

# File upload
st.header("📁 Upload FASTA File")
uploaded_file = st.file_uploader(
    "Choose a FASTA file",
    type=['fa', 'fasta', 'fas', 'fna'],
    help="Upload genomic sequences in FASTA format"
)

# Main scanning logic
if uploaded_file and selected_patterns:
    # Create config
    config = ScanConfig(
        chunk_size=chunk_size,
        overlap=overlap,
        max_workers=max_workers,
        stream_results=True,
        incremental_save=True,
        output_format='csv'
    )
    
    # Get cached scanner (compiled once)
    scanner = get_cached_scanner(selected_patterns, config.__dict__)
    
    # Save uploaded file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp_file:
        tmp_file.write(uploaded_file.getvalue().decode())
        tmp_path = tmp_file.name
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        if st.button("🚀 Start Scanning", type="primary", use_container_width=True):
            try:
                start_time = time.time()
                
                # Scan with progress tracking
                with st.spinner("Initializing scanner..."):
                    results = streamlit_scan_with_progress(scanner, tmp_path)
                
                scan_time = time.time() - start_time
                
                # Display results
                st.success(f"✅ Scanning completed in {scan_time:.2f} seconds!")
                
                # Statistics
                stats = results['stats']
                col1, col2, col3, col4 = st.columns(4)
                
                with col1:
                    st.metric("🎯 Total Matches", f"{stats['total_matches']:,}")
                with col2:
                    st.metric("🧬 Sequences", stats['sequences_scanned'])
                with col3:
                    st.metric("🔍 Patterns Found", stats['patterns_matched'])
                with col4:
                    st.metric("⏱️ Scan Time", f"{scan_time:.1f}s")
                
                # Results preview and download
                if 'csv' in results['files']:
                    st.subheader("📊 Results Preview")
                    
                    # Load and display sample results
                    df = pd.read_csv(results['files']['csv'])
                    
                    if len(df) > 0:
                        # Show top matches
                        st.dataframe(
                            df.head(100),
                            use_container_width=True,
                            height=300
                        )
                        
                        # Download button
                        with open(results['files']['csv'], 'r') as f:
                            st.download_button(
                                "📥 Download Full Results (CSV)",
                                f.read(),
                                file_name=f"genomic_scan_results_{int(time.time())}.csv",
                                mime="text/csv",
                                use_container_width=True
                            )
                    else:
                        st.info("No matches found with the selected patterns.")
                
            except Exception as e:
                st.error(f"❌ Scanning failed: {str(e)}")
                st.exception(e)
            
            finally:
                # Cleanup
                try:
                    os.unlink(tmp_path)
                except:
                    pass
    
    with col2:
        st.info("💡 **Tips:**\n"
                "- Larger chunk sizes = faster but more memory\n"
                "- More workers = faster on multi-core systems\n"
                "- Increase overlap for longer patterns\n"
                "- Results are saved incrementally")

# Footer
st.markdown("---")
st.markdown("Built with ❤️ using Hyperscan + Streamlit")
```

### 2. Advanced Configuration Example

```python
# Advanced scanner configuration
config = ScanConfig(
    chunk_size=5_000_000,      # 5MB chunks
    overlap=1000,              # 1kb overlap
    max_workers=6,             # 6 parallel workers
    timeout=600,               # 10 minute timeout
    stream_results=True,       # Stream for memory efficiency
    output_format='parquet',   # Better performance than CSV
    incremental_save=True,     # Save results as they come
    max_memory_mb=2000        # 2GB memory limit
)

# Pattern examples for different use cases
COMPREHENSIVE_PATTERNS = {
    'g_quadruplex': [
        r'G{3,}N{1,7}G{3,}N{1,7}G{3,}N{1,7}G{3,}',  # Canonical
        r'G{2,}N{1,12}G{2,}N{1,12}G{2,}N{1,12}G{2,}', # Relaxed
    ],
    'microsatellites': [
        r'(CA){6,}|(TG){6,}',     # Dinucleotide
        r'(CAG){4,}|(CTG){4,}',   # Trinucleotide
        r'(GATA){3,}',            # Tetranucleotide
    ],
    'cpg_islands': [
        r'[CG]{50,}',             # Simple CpG rich
        r'(?:[CG][ATCG]{0,19}){5,}[CG]',  # Spaced CpG
    ],
    'palindromes': [
        r'([ATCG]{4,})N{0,10}\1',  # Perfect palindromes
        r'([ATCG]{6,})N{1,20}(?:[ATCG](?=\1))', # Imperfect
    ]
}
```

### 3. Memory-Efficient Large Genome Processing

```python
def process_large_genome(fasta_path: str, output_dir: str):
    """
    Process very large genomes (>1GB) with memory optimization
    """
    
    # Conservative settings for large genomes
    config = ScanConfig(
        chunk_size=1_000_000,    # 1MB chunks (smaller for memory)
        overlap=500,
        max_workers=2,           # Fewer workers to reduce memory
        stream_results=True,
        incremental_save=True,
        max_memory_mb=500       # Strict memory limit
    )
    
    patterns = [
        r'G{3,}N{1,7}G{3,}N{1,7}G{3,}N{1,7}G{3,}',  # G-quadruplex only
    ]
    
    scanner = HyperscanStreamlitScanner(patterns, config)
    
    # Progress tracking
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    def progress_callback(progress: float, message: str):
        progress_bar.progress(progress)
        status_text.text(f"🔄 {message}")
        
        # Force Streamlit to update
        time.sleep(0.01)
    
    # Process with streaming output
    results = scanner.scan_fasta_file(
        fasta_path, 
        output_prefix=f"{output_dir}/large_genome_scan",
        progress_callback=progress_callback
    )
    
    progress_bar.empty()
    status_text.empty()
    
    return results
```

### 4. Real-time Results Display

```python
def create_realtime_results_display():
    """
    Create a real-time updating results display
    """
    
    # Containers for dynamic updates
    results_container = st.container()
    stats_container = st.container()
    
    # Initialize display
    with stats_container:
        col1, col2, col3 = st.columns(3)
        matches_metric = col1.empty()
        sequences_metric = col2.empty()
        time_metric = col3.empty()
    
    with results_container:
        results_table = st.empty()
    
    def update_display(current_results: List[dict], scan_time: float):
        """Update the display with current results"""
        
        # Update metrics
        matches_metric.metric("🎯 Matches Found", len(current_results))
        sequences_metric.metric("🧬 Sequences", len(set(r['sequence_name'] for r in current_results)))
        time_metric.metric("⏱️ Elapsed", f"{scan_time:.1f}s")
        
        # Update results table (show latest 50)
        if current_results:
            df = pd.DataFrame(current_results[-50:])  # Latest 50 results
            results_table.dataframe(
                df[['sequence_name', 'start', 'end', 'pattern', 'length']],
                use_container_width=True,
                height=200
            )
    
    return update_display

# Usage in scanning loop
update_display = create_realtime_results_display()
all_results = []
start_time = time.time()

# Simulate streaming results
for batch in result_batches:
    all_results.extend(batch)
    current_time = time.time() - start_time
    update_display(all_results, current_time)
    time.sleep(0.1)  # Allow UI to update
```

### 5. Error Handling and Validation

```python
def robust_scanner_with_validation(patterns: List[str], fasta_file) -> dict:
    """
    Robust scanner with comprehensive error handling and validation
    """
    
    # Input validation
    if not patterns:
        st.error("❌ No patterns provided")
        return {}
    
    if not fasta_file:
        st.error("❌ No FASTA file uploaded")
        return {}
    
    # Validate patterns
    valid_patterns = []
    invalid_patterns = []
    
    for pattern in patterns:
        try:
            # Test compile pattern
            import re
            re.compile(pattern)
            valid_patterns.append(pattern)
        except re.error as e:
            invalid_patterns.append((pattern, str(e)))
    
    if invalid_patterns:
        st.warning(f"⚠️ Invalid patterns found:")
        for pattern, error in invalid_patterns:
            st.text(f"  • {pattern}: {error}")
    
    if not valid_patterns:
        st.error("❌ No valid patterns found")
        return {}
    
    # File size check
    file_size = len(fasta_file.getvalue())
    if file_size > 100_000_000:  # 100MB limit
        st.warning(f"⚠️ Large file detected ({file_size/1e6:.1f} MB). This may take a while.")
        
        # Adjust config for large files
        config = ScanConfig(
            chunk_size=1_000_000,
            max_workers=min(2, multiprocessing.cpu_count()),
            timeout=1800  # 30 minutes
        )
    else:
        config = ScanConfig()
    
    try:
        # Create scanner with validation
        scanner = HyperscanStreamlitScanner(valid_patterns, config)
        
        # Save uploaded file safely
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp_file:
            content = fasta_file.getvalue().decode('utf-8', errors='replace')
            tmp_file.write(content)
            tmp_path = tmp_file.name
        
        # Scan with timeout
        start_time = time.time()
        
        with st.spinner("🔍 Scanning in progress..."):
            results = scanner.scan_fasta_file(tmp_path)
        
        scan_time = time.time() - start_time
        
        # Cleanup
        os.unlink(tmp_path)
        
        # Success metrics
        st.success(f"✅ Scan completed successfully in {scan_time:.2f} seconds")
        
        return {
            'results': results,
            'scan_time': scan_time,
            'patterns_used': valid_patterns,
            'file_size': file_size
        }
        
    except Exception as e:
        st.error(f"❌ Scanning failed: {str(e)}")
        
        # Detailed error information
        with st.expander("🔍 Error Details"):
            st.exception(e)
        
        return {}
```

### 6. Performance Monitoring

```python
def add_performance_monitoring():
    """
    Add performance monitoring to your Streamlit app
    """
    
    # Performance metrics container
    perf_container = st.sidebar.container()
    
    with perf_container:
        st.subheader("📊 Performance")
        
        # System info
        cpu_count = multiprocessing.cpu_count()
        st.text(f"💻 CPU Cores: {cpu_count}")
        
        # Memory usage (if psutil available)
        try:
            import psutil
            memory = psutil.virtual_memory()
            st.text(f"🧠 RAM: {memory.percent:.1f}% used")
            st.text(f"💾 Available: {memory.available / 1e9:.1f} GB")
        except ImportError:
            st.text("📊 Install psutil for memory monitoring")
        
        # Performance recommendations
        st.info("💡 **Performance Tips:**\n"
                "- Use 50-80% of available CPU cores\n"
                "- Monitor memory usage during scans\n"
                "- Larger chunks = less overhead\n"
                "- Stream results for large datasets")

# Add to your main app
add_performance_monitoring()
```

## 🧪 Testing Your Implementation

### Basic Test Script

```python
# test_your_implementation.py
import tempfile
from production_hyperscan_streamlit import HyperscanStreamlitScanner, ScanConfig

def test_basic_functionality():
    """Basic test to ensure everything works"""
    
    # Test patterns
    patterns = ['GGGG', 'CCCC', 'AAAA']
    
    # Test sequence
    test_seq = "ATCGGGGGCCCCAAAATAGC" * 100
    
    # Create scanner
    scanner = HyperscanStreamlitScanner(patterns)
    
    # Scan
    results = scanner.scan_sequence(test_seq, "test")
    
    print(f"✅ Found {len(results)} matches")
    print(f"✅ Patterns found: {set(r.pattern_id for r in results)}")
    
    return len(results) > 0

def test_file_processing():
    """Test FASTA file processing"""
    
    # Create test file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(">test_seq\n")
        f.write("ATCGGGGGCCCCAAAATAGC" * 200)
        test_file = f.name
    
    # Test processing
    patterns = ['GGGG', 'CCCC']
    scanner = HyperscanStreamlitScanner(patterns)
    
    results = scanner.scan_fasta_file(test_file)
    
    print(f"✅ File scan completed")
    print(f"✅ Stats: {results['stats']}")
    
    # Cleanup
    os.unlink(test_file)
    
    return results['stats']['total_matches'] > 0

if __name__ == "__main__":
    print("🧪 Testing Hyperscan-Streamlit Integration...")
    
    success1 = test_basic_functionality()
    success2 = test_file_processing()
    
    if success1 and success2:
        print("🎉 All tests passed! Ready for production.")
    else:
        print("❌ Some tests failed. Check your setup.")
```

## 🚀 Deployment Checklist

### Pre-deployment Validation

- [ ] Hyperscan properly installed (`pip install hyperscan`)
- [ ] All dependencies available (`streamlit`, `pandas`, `biopython`)
- [ ] Test with sample data
- [ ] Performance benchmarks completed
- [ ] Error handling tested
- [ ] Memory limits configured appropriately

### Production Configuration

```python
# production_config.py
PRODUCTION_CONFIG = ScanConfig(
    chunk_size=2_000_000,      # 2MB chunks
    overlap=500,               # Conservative overlap
    max_workers=4,             # Conservative worker count
    timeout=600,               # 10 minute timeout
    stream_results=True,       # Always stream in production
    incremental_save=True,     # Always save incrementally
    max_memory_mb=1000        # 1GB memory limit
)
```

### Monitoring and Logging

```python
import logging

# Configure logging for production
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('hyperscan_scanner.log'),
        logging.StreamHandler()
    ]
)
```

## 📚 Additional Resources

- **Hyperscan Documentation**: https://hyperscan.io/
- **Streamlit Documentation**: https://docs.streamlit.io/
- **Pattern Testing Tool**: Use `re.compile()` to test patterns before deployment
- **Performance Profiling**: Use `cProfile` to identify bottlenecks

## 🤝 Support

For issues or questions:
1. Check the test suite for examples
2. Review error logs for specific problems
3. Verify Hyperscan installation and compatibility
4. Test with smaller datasets first

---

**Happy Scanning! 🧬🚀**