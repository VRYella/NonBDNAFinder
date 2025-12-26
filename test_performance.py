"""
Test suite for performance and memory optimizations.

This module validates:
- Compressed file support
- Memory management utilities
- DataFrame optimization
- Large file handling
"""

import os
import gzip
import tempfile
import pandas as pd
import numpy as np
from utilities import (
    open_compressed_file,
    parse_fasta_chunked_compressed,
    trigger_garbage_collection,
    optimize_dataframe_memory,
    get_memory_usage_mb
)


def test_compressed_file_detection():
    """Test automatic detection of compressed files."""
    # Create a temporary gzipped FASTA file
    test_fasta = ">test_seq\nACGTACGTACGTACGT\n"
    
    with tempfile.NamedTemporaryFile(mode='wb', suffix='.fa.gz', delete=False) as f:
        temp_file = f.name
        f.write(gzip.compress(test_fasta.encode('utf-8')))
    
    try:
        # Test automatic decompression
        with open_compressed_file(temp_file) as f:
            content = f.read()
            assert ">test_seq" in content
            assert "ACGTACGTACGTACGT" in content
        print("✓ Compressed file detection test passed")
    finally:
        os.unlink(temp_file)


def test_chunked_compressed_parsing():
    """Test parsing of compressed FASTA files."""
    test_sequences = {
        "seq1": "ACGTACGTACGTACGT",
        "seq2": "GGGGCCCCAAAATTTT",
        "seq3": "TATATATATATATAT"
    }
    
    # Create multi-FASTA content
    fasta_content = ""
    for name, seq in test_sequences.items():
        fasta_content += f">{name}\n{seq}\n"
    
    # Create compressed file
    with tempfile.NamedTemporaryFile(mode='wb', suffix='.fa.gz', delete=False) as f:
        temp_file = f.name
        f.write(gzip.compress(fasta_content.encode('utf-8')))
    
    try:
        # Parse compressed file
        parsed_seqs = {}
        for name, seq in parse_fasta_chunked_compressed(temp_file):
            parsed_seqs[name] = seq
        
        # Validate
        assert len(parsed_seqs) == 3
        for name, expected_seq in test_sequences.items():
            assert name in parsed_seqs
            assert parsed_seqs[name] == expected_seq
        
        print("✓ Chunked compressed parsing test passed")
    finally:
        os.unlink(temp_file)


def test_memory_optimization():
    """Test DataFrame memory optimization."""
    # Create a large DataFrame with various types
    data = {
        'id': np.arange(10000, dtype='int64'),
        'score': np.random.rand(10000).astype('float64'),
        'category': np.random.randint(0, 10, 10000, dtype='int64')
    }
    df = pd.DataFrame(data)
    
    # Get initial memory usage
    initial_memory = df.memory_usage(deep=True).sum()
    
    # Optimize
    df_optimized = optimize_dataframe_memory(df)
    
    # Get optimized memory usage
    optimized_memory = df_optimized.memory_usage(deep=True).sum()
    
    # Should see significant reduction
    reduction_pct = (1 - optimized_memory / initial_memory) * 100
    
    print(f"✓ DataFrame optimization test passed")
    print(f"  Memory reduced by {reduction_pct:.1f}%")
    print(f"  Initial: {initial_memory / 1024:.1f} KB")
    print(f"  Optimized: {optimized_memory / 1024:.1f} KB")
    
    assert optimized_memory < initial_memory, "Optimization should reduce memory"


def test_garbage_collection():
    """Test garbage collection trigger."""
    initial_mem = get_memory_usage_mb()
    
    # Create some temporary objects
    large_list = [i for i in range(1000000)]
    del large_list
    
    # Trigger collection
    collected = trigger_garbage_collection()
    
    final_mem = get_memory_usage_mb()
    
    print(f"✓ Garbage collection test passed")
    print(f"  Objects collected: {collected}")
    print(f"  Memory before: {initial_mem:.1f} MB")
    print(f"  Memory after: {final_mem:.1f} MB")


def test_memory_monitoring():
    """Test memory usage monitoring."""
    mem_mb = get_memory_usage_mb()
    
    assert mem_mb > 0, "Memory usage should be positive"
    assert mem_mb < 10000, "Memory usage should be reasonable"
    
    print(f"✓ Memory monitoring test passed")
    print(f"  Current memory usage: {mem_mb:.1f} MB")


def test_large_sequence_handling():
    """Test handling of large sequences with memory management."""
    # Create a large test sequence (10 MB)
    large_seq = "ACGT" * 2500000  # 10 MB sequence
    
    test_fasta = f">large_seq\n{large_seq}\n"
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        temp_file = f.name
        f.write(test_fasta)
    
    try:
        initial_mem = get_memory_usage_mb()
        
        # Parse using chunked method
        parsed_seqs = {}
        for name, seq in parse_fasta_chunked_compressed(temp_file, chunk_size_mb=1):
            parsed_seqs[name] = seq
            # Validate sequence was parsed correctly
            assert len(seq) == len(large_seq)
        
        # Cleanup
        del parsed_seqs
        trigger_garbage_collection()
        
        final_mem = get_memory_usage_mb()
        
        print(f"✓ Large sequence handling test passed")
        print(f"  Sequence size: {len(large_seq) / 1024 / 1024:.1f} MB")
        print(f"  Memory delta: {final_mem - initial_mem:.1f} MB")
    finally:
        os.unlink(temp_file)


def run_all_tests():
    """Run all performance tests."""
    print("\n" + "="*70)
    print("Running Performance & Memory Optimization Tests")
    print("="*70 + "\n")
    
    tests = [
        test_compressed_file_detection,
        test_chunked_compressed_parsing,
        test_memory_optimization,
        test_garbage_collection,
        test_memory_monitoring,
        test_large_sequence_handling
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            print(f"\n{test.__name__}:")
            test()
            passed += 1
        except Exception as e:
            print(f"✗ {test.__name__} failed: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "="*70)
    print(f"Test Results: {passed} passed, {failed} failed")
    print("="*70 + "\n")
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    exit(0 if success else 1)
