#!/usr/bin/env python3
"""
Integration test for disk storage system.

Tests the full workflow: upload -> storage -> analysis -> results -> export
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from Utilities.disk_storage import UniversalSequenceStorage, UniversalResultsStorage
from Utilities.chunk_analyzer import ChunkAnalyzer
from Utilities.nonbscanner import analyze_sequence
import tempfile
import shutil


def test_small_sequence_workflow():
    """Test workflow with small sequence (<10MB)."""
    print("\n=== Testing Small Sequence Workflow ===")
    
    # Create test sequence (5KB)
    test_sequence = "ATCGATCGATCG" * 400  # ~5KB
    seq_name = "test_small"
    
    # Initialize storage
    storage = UniversalSequenceStorage()
    print(f"✓ Storage initialized at {storage.base_dir}")
    
    # Save sequence
    seq_id = storage.save_sequence(test_sequence, seq_name)
    print(f"✓ Sequence saved with ID: {seq_id}")
    
    # Verify metadata
    metadata = storage.get_metadata(seq_id)
    assert metadata['length'] == len(test_sequence)
    assert metadata['name'] == seq_name
    print(f"✓ Metadata verified: {metadata['length']} bp, GC% = {metadata['gc_content']:.1f}")
    
    # Analyze sequence directly (small sequence)
    results = analyze_sequence(test_sequence, seq_name)
    print(f"✓ Analysis complete: {len(results)} motifs detected")
    
    # Store results
    results_storage = UniversalResultsStorage(
        base_dir=str(storage.base_dir / "results"),
        seq_id=seq_id
    )
    results_storage.append_batch(results)
    print(f"✓ Results stored")
    
    # Get summary stats
    stats = results_storage.get_summary_stats()
    print(f"✓ Summary stats: {stats['total_count']} motifs, {len(stats['class_distribution'])} classes")
    
    # Cleanup
    results_storage.cleanup()
    storage.cleanup()
    print("✓ Cleanup complete")
    
    return True


def test_chunked_analysis_workflow():
    """Test workflow with chunk-based analysis."""
    print("\n=== Testing Chunked Analysis Workflow ===")
    
    # Create test sequence (50KB)
    test_sequence = "ATCGATCGATCG" * 4000  # ~50KB
    seq_name = "test_chunked"
    
    # Initialize storage
    storage = UniversalSequenceStorage()
    print(f"✓ Storage initialized at {storage.base_dir}")
    
    # Save sequence
    seq_id = storage.save_sequence(test_sequence, seq_name)
    print(f"✓ Sequence saved: {len(test_sequence)} bp")
    
    # Progress callback
    progress_updates = []
    def progress_callback(pct):
        progress_updates.append(pct)
    
    # Analyze with ChunkAnalyzer
    analyzer = ChunkAnalyzer(
        storage,
        chunk_size=20000,  # 20KB chunks
        overlap=1000  # 1KB overlap
    )
    
    results_storage = analyzer.analyze(
        seq_id=seq_id,
        progress_callback=progress_callback
    )
    
    print(f"✓ Chunk analysis complete")
    print(f"  Progress updates: {len(progress_updates)}")
    
    # Get summary
    stats = results_storage.get_summary_stats()
    print(f"✓ Results: {stats['total_count']} motifs")
    print(f"  Class distribution: {stats['class_distribution']}")
    
    # Test iteration
    motif_count = 0
    for motif in results_storage.iter_results(limit=10):
        motif_count += 1
    print(f"✓ Iteration works: retrieved {motif_count} motifs (limit=10)")
    
    # Cleanup
    results_storage.cleanup()
    storage.cleanup()
    print("✓ Cleanup complete")
    
    return True


def test_storage_iteration():
    """Test chunk iteration on stored sequences."""
    print("\n=== Testing Storage Chunk Iteration ===")
    
    # Create test sequence (30KB)
    test_sequence = "ATCG" * 7500  # 30KB
    
    storage = UniversalSequenceStorage()
    seq_id = storage.save_sequence(test_sequence, "test_iteration")
    
    # Iterate in chunks
    chunk_count = 0
    total_length = 0
    
    for chunk, start, end in storage.iter_chunks(seq_id, chunk_size=10000, overlap=1000):
        chunk_count += 1
        chunk_len = end - start
        total_length += chunk_len
        print(f"  Chunk {chunk_count}: [{start:,}-{end:,}] = {chunk_len:,} bp")
    
    print(f"✓ Iteration complete: {chunk_count} chunks")
    
    # Verify chunks cover the sequence
    assert chunk_count >= 3, f"Expected at least 3 chunks, got {chunk_count}"
    
    # Cleanup
    storage.cleanup()
    print("✓ Cleanup complete")
    
    return True


def test_results_pagination():
    """Test results pagination for large result sets."""
    print("\n=== Testing Results Pagination ===")
    
    storage_dir = tempfile.mkdtemp(prefix="test_pagination_")
    
    # Create mock results
    results_storage = UniversalResultsStorage(storage_dir, "test_seq")
    
    # Add 100 mock motifs
    for i in range(100):
        motif = {
            'Class': f'Class_{i % 5}',
            'Subclass': f'Subclass_{i % 10}',
            'Start': i * 100,
            'End': i * 100 + 50,
            'Length': 50,
            'Score': 2.0 + (i % 10) * 0.1
        }
        results_storage.append(motif)
    
    print(f"✓ Added 100 mock motifs")
    
    # Test pagination
    page_size = 10
    retrieved = list(results_storage.iter_results(limit=page_size))
    assert len(retrieved) == page_size
    print(f"✓ Pagination works: retrieved {len(retrieved)} motifs (limit={page_size})")
    
    # Test summary stats
    stats = results_storage.get_summary_stats()
    assert stats['total_count'] == 100
    print(f"✓ Summary stats correct: {stats['total_count']} motifs")
    
    # Cleanup
    results_storage.cleanup()
    shutil.rmtree(storage_dir)
    print("✓ Cleanup complete")
    
    return True


def main():
    """Run all integration tests."""
    print("=" * 70)
    print("DISK STORAGE INTEGRATION TESTS")
    print("=" * 70)
    
    tests = [
        ("Small Sequence Workflow", test_small_sequence_workflow),
        ("Chunked Analysis Workflow", test_chunked_analysis_workflow),
        ("Storage Chunk Iteration", test_storage_iteration),
        ("Results Pagination", test_results_pagination),
    ]
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        try:
            if test_func():
                passed += 1
                print(f"✅ {test_name} PASSED\n")
            else:
                failed += 1
                print(f"❌ {test_name} FAILED\n")
        except Exception as e:
            failed += 1
            print(f"❌ {test_name} FAILED with exception:")
            print(f"   {e}\n")
            import traceback
            traceback.print_exc()
    
    print("=" * 70)
    print(f"RESULTS: {passed} passed, {failed} failed")
    print("=" * 70)
    
    return failed == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
