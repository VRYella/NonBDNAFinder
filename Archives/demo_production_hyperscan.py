#!/usr/bin/env python3
"""
Demo script showing the production-ready Hyperscan-Streamlit integration in action.

This script demonstrates:
- Basic scanning functionality
- Chunking with large sequences
- File processing
- Performance characteristics
"""

import time
import tempfile
import os
from production_hyperscan_streamlit import (
    HyperscanStreamlitScanner, 
    ScanConfig, 
    HYPERSCAN_AVAILABLE
)

def demo_basic_scanning():
    """Demonstrate basic scanning functionality"""
    print("🧪 Demo 1: Basic Scanning")
    print("=" * 50)
    
    # Sample genomic patterns
    patterns = [
        r'G{4,}',           # G-quadruplex motif
        r'C{4,}',           # C-rich regions
        r'A{6,}',           # A-rich regions
        r'T{6,}',           # T-rich regions
        r'[AT]{10,}',       # AT-rich regions
    ]
    
    # Test sequence with known patterns
    test_sequence = (
        "ATCGATCG" + "GGGG" + "ATCG" + "CCCC" + "ATCG" +
        "AAAAAA" + "ATCG" + "TTTTTT" + "ATCG" +
        "ATATATATATATATATAT" + "GCTAGCTA"
    ) * 50  # Repeat to make it interesting
    
    print(f"📊 Patterns: {len(patterns)}")
    print(f"📏 Sequence length: {len(test_sequence):,} bp")
    
    # Create scanner
    start_time = time.time()
    scanner = HyperscanStreamlitScanner(patterns)
    compile_time = time.time() - start_time
    
    print(f"⚡ Database compiled in {compile_time:.3f}s")
    
    # Scan sequence
    start_time = time.time()
    results = scanner.scan_sequence(test_sequence, "demo_sequence")
    scan_time = time.time() - start_time
    
    print(f"🔍 Scanning completed in {scan_time:.3f}s")
    print(f"🎯 Found {len(results)} total matches")
    
    # Show pattern breakdown
    pattern_counts = {}
    for result in results:
        pattern = result.pattern
        pattern_counts[pattern] = pattern_counts.get(pattern, 0) + 1
    
    print("\n📈 Pattern breakdown:")
    for pattern, count in pattern_counts.items():
        print(f"  • {pattern}: {count} matches")
    
    print()


def demo_chunking():
    """Demonstrate chunking with large sequences"""
    print("🧪 Demo 2: Chunking Large Sequences")
    print("=" * 50)
    
    # Create a large sequence
    base_sequence = "ATCGATCG" + "GGGG" + "ATCG" + "CCCC" + "ATCGCGTA"
    large_sequence = base_sequence * 10000  # ~640KB sequence
    
    patterns = [r'GGGG', r'CCCC']
    
    print(f"📏 Large sequence length: {len(large_sequence):,} bp")
    
    # Configure for chunking
    config = ScanConfig(
        chunk_size=50000,  # 50KB chunks
        overlap=100,       # 100bp overlap
        max_workers=2      # 2 workers for demo
    )
    
    scanner = HyperscanStreamlitScanner(patterns, config)
    
    # Track progress
    progress_updates = []
    def progress_callback(progress, message):
        progress_updates.append((progress, message))
        print(f"📊 Progress: {progress*100:.1f}% - {message}")
    
    print("\n🚀 Starting chunked scan...")
    start_time = time.time()
    
    results = scanner.scan_sequence(
        large_sequence, 
        "large_demo_sequence",
        progress_callback
    )
    
    scan_time = time.time() - start_time
    
    print(f"\n✅ Chunked scanning completed!")
    print(f"⏱️  Total time: {scan_time:.3f}s")
    print(f"🎯 Found {len(results)} matches")
    print(f"📊 Progress updates: {len(progress_updates)}")
    
    # Show chunking effectiveness
    chunks = scanner._create_chunks(large_sequence)
    print(f"🔄 Sequence split into {len(chunks)} chunks")
    print(f"📦 Chunk size: {config.chunk_size:,} bp")
    print(f"🔗 Overlap: {config.overlap} bp")
    print()


def demo_file_processing():
    """Demonstrate FASTA file processing"""
    print("🧪 Demo 3: FASTA File Processing")
    print("=" * 50)
    
    # Create a test FASTA file
    sequences = [
        ("seq1", "ATCGATCG" + "GGGG" * 10 + "ATCG" + "CCCC" * 5 + "ATCG"),
        ("seq2", "GCTAGCTA" + "AAAA" * 8 + "GCTA" + "TTTT" * 6 + "GCTA"),
        ("seq3", "ATCGATCG" + "GGGG" + "CCCC" + "AAAA" + "TTTT" + "ATCG" * 100),
    ]
    
    # Write FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        for seq_id, sequence in sequences:
            f.write(f">{seq_id}\n{sequence}\n")
        fasta_path = f.name
    
    print(f"📁 Created test FASTA with {len(sequences)} sequences")
    
    # Set up scanner
    patterns = [r'G{4,}', r'C{4,}', r'A{4,}', r'T{4,}']
    config = ScanConfig(
        chunk_size=1000,   # Small chunks for demo
        max_workers=2,
        stream_results=True,
        incremental_save=True
    )
    
    scanner = HyperscanStreamlitScanner(patterns, config)
    
    # Process file
    print("\n🚀 Processing FASTA file...")
    start_time = time.time()
    
    results = scanner.scan_fasta_file(fasta_path, output_prefix="/tmp/demo_scan")
    
    process_time = time.time() - start_time
    
    print(f"\n✅ File processing completed!")
    print(f"⏱️  Total time: {process_time:.3f}s")
    print(f"📊 Statistics:")
    
    stats = results['stats']
    for key, value in stats.items():
        print(f"  • {key}: {value}")
    
    # Show output files
    print(f"\n📂 Output files:")
    for file_type, file_path in results.get('files', {}).items():
        if file_path and os.path.exists(file_path):
            file_size = os.path.getsize(file_path)
            print(f"  • {file_type}: {file_path} ({file_size} bytes)")
    
    # Cleanup
    try:
        os.unlink(fasta_path)
        for file_path in results.get('files', {}).values():
            if file_path and os.path.exists(file_path):
                os.unlink(file_path)
    except Exception:
        pass
    
    print()


def demo_performance_comparison():
    """Demonstrate performance characteristics"""
    print("🧪 Demo 4: Performance Characteristics")
    print("=" * 50)
    
    # Test different sequence sizes
    base_seq = "ATCGATCG" + "GGGG" + "CCCC" + "AAAA" + "TTTT" + "ATCGCGTA"
    patterns = [r'GGGG', r'CCCC', r'AAAA', r'TTTT']
    
    sizes = [1000, 10000, 100000]  # Different sequence sizes
    
    scanner = HyperscanStreamlitScanner(patterns)
    
    print("📊 Performance test results:")
    print("Size (bp)    | Compile (s) | Scan (s) | Matches | Rate (KB/s)")
    print("-" * 65)
    
    for size in sizes:
        test_seq = (base_seq * (size // len(base_seq) + 1))[:size]
        
        # Measure compilation (first time)
        start_time = time.time()
        temp_scanner = HyperscanStreamlitScanner(patterns)
        compile_time = time.time() - start_time
        
        # Measure scanning
        start_time = time.time()
        results = scanner.scan_sequence(test_seq, f"perf_test_{size}")
        scan_time = time.time() - start_time
        
        # Calculate rate
        rate = (size / 1000) / scan_time if scan_time > 0 else 0
        
        print(f"{size:8,} | {compile_time:10.4f} | {scan_time:7.4f} | {len(results):6,} | {rate:9.1f}")
    
    print()


def run_all_demos():
    """Run all demonstration scripts"""
    if not HYPERSCAN_AVAILABLE:
        print("❌ Hyperscan not available. Please install with: pip install hyperscan")
        return
    
    print("🚀 Production Hyperscan-Streamlit Integration Demo")
    print("=" * 60)
    print("This demo shows the key features of the production-ready scanner:")
    print("• Efficient pattern compilation and scanning")
    print("• Chunking for large sequences")
    print("• File processing with multiple sequences")
    print("• Performance characteristics")
    print("=" * 60)
    print()
    
    try:
        demo_basic_scanning()
        demo_chunking()
        demo_file_processing()
        demo_performance_comparison()
        
        print("🎉 All demos completed successfully!")
        print("\n💡 Next steps:")
        print("1. Try the Streamlit app: streamlit run production_hyperscan_streamlit.py")
        print("2. Check the comprehensive guide: HYPERSCAN_STREAMLIT_GUIDE.md")
        print("3. Run the test suite: pytest test_production_hyperscan.py")
        
    except Exception as e:
        print(f"❌ Demo failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    run_all_demos()