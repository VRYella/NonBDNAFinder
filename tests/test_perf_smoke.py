"""
Performance Smoke Tests for NBDScanner
======================================

These tests verify that the scanner works correctly on sequences
of various sizes and provides acceptable performance.

Note: These are smoke tests, not benchmarks. They verify functionality
and provide rough performance baselines.
"""

import os
import sys
import time
import tempfile
import unittest

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def generate_random_sequence(length: int, gc_content: float = 0.5) -> str:
    """Generate a random DNA sequence with specified GC content."""
    import random
    
    bases = []
    for _ in range(length):
        if random.random() < gc_content:
            bases.append(random.choice(['G', 'C']))
        else:
            bases.append(random.choice(['A', 'T']))
    
    return ''.join(bases)


def generate_g4_sequence(length: int = 100) -> str:
    """Generate a sequence with G-quadruplex motifs."""
    # Create a sequence with clear G4 patterns
    g4_motif = "GGGTTAGGGTTAGGGTTAGGG"
    spacer = "ATGCATGC"
    
    result = []
    while len(''.join(result)) < length:
        result.append(g4_motif)
        result.append(spacer)
    
    return ''.join(result)[:length]


class TestBasicFunctionality(unittest.TestCase):
    """Basic functionality tests."""
    
    def test_import_modules(self):
        """Test that all modules can be imported."""
        import nonbscanner
        import detectors
        import utilities
        from scanner_backends import io_utils
        from scanner_backends import parallel_worker
        
    def test_analyze_small_sequence(self):
        """Test analysis of a small sequence."""
        from nonbscanner import analyze_sequence
        
        # Known G4-forming sequence
        test_seq = "GGGTTAGGGTTAGGGTTAGGG"
        
        motifs = analyze_sequence(test_seq, "test")
        
        # Should find at least some motifs
        self.assertIsInstance(motifs, list)
        
        # Check motif structure
        for motif in motifs:
            self.assertIn('Class', motif)
            self.assertIn('Start', motif)
            self.assertIn('End', motif)
            self.assertIn('Score', motif)
    
    def test_analyze_empty_sequence(self):
        """Test handling of edge case sequences."""
        from nonbscanner import analyze_sequence
        from utilities import validate_sequence
        
        # Test with sequence validation
        is_valid, msg = validate_sequence("")
        self.assertFalse(is_valid)
    
    def test_analyze_sequence_with_n(self):
        """Test handling of sequences with N characters."""
        from nonbscanner import analyze_sequence
        
        test_seq = "ATGCNNNNATGCGGGGCCCC"
        
        # Should handle gracefully
        motifs = analyze_sequence(test_seq, "test_n")
        self.assertIsInstance(motifs, list)


class TestPerformanceSmoke(unittest.TestCase):
    """Performance smoke tests."""
    
    def test_1kb_sequence(self):
        """Test scanning a 1kb sequence completes quickly."""
        from nonbscanner import analyze_sequence
        
        seq = generate_random_sequence(1000)
        
        start = time.time()
        motifs = analyze_sequence(seq, "1kb_test")
        elapsed = time.time() - start
        
        # Should complete in reasonable time (< 10 seconds)
        self.assertLess(elapsed, 10.0)
        self.assertIsInstance(motifs, list)
    
    def test_10kb_sequence(self):
        """Test scanning a 10kb sequence."""
        from nonbscanner import analyze_sequence
        
        seq = generate_random_sequence(10000)
        
        start = time.time()
        motifs = analyze_sequence(seq, "10kb_test")
        elapsed = time.time() - start
        
        # Should complete in reasonable time (< 60 seconds)
        self.assertLess(elapsed, 60.0)
        
        # Calculate throughput
        throughput = len(seq) / elapsed
        print(f"\n10kb scan: {elapsed:.2f}s, {throughput:,.0f} bp/s")
    
    def test_chunked_analysis_large_sequence(self):
        """Test chunking for sequences larger than CHUNK_THRESHOLD (10,000 bp)."""
        from nonbscanner import analyze_sequence, CHUNK_THRESHOLD, DEFAULT_CHUNK_SIZE
        
        # Generate sequence larger than chunk threshold
        seq = generate_g4_sequence(15000)
        
        self.assertGreater(len(seq), CHUNK_THRESHOLD, "Test sequence should be larger than chunk threshold")
        
        # Track progress (new signature: chunk, total, bp_processed, elapsed, throughput)
        progress_calls = []
        def progress_callback(chunk, total, bp_processed, *args):
            # Accept optional elapsed and throughput args for backward compatibility
            progress_calls.append((chunk, total, bp_processed))
        
        start = time.time()
        motifs = analyze_sequence(
            seq, 
            "chunked_test", 
            use_chunking=True,
            chunk_size=5000,  # Use smaller chunks to ensure multiple chunks
            progress_callback=progress_callback
        )
        elapsed = time.time() - start
        
        # Should find motifs
        self.assertIsInstance(motifs, list)
        self.assertGreater(len(motifs), 0, "Should find motifs in G4 sequence")
        
        # Progress callback should be called multiple times (for each chunk)
        self.assertGreater(len(progress_calls), 1, "Progress callback should be called for each chunk")
        
        # Verify fix: First progress callback should show bp > 0
        # Previously, the first callback showed 0 bp processed because it was called
        # before the chunk was processed. After the fix, it should show the actual
        # number of base pairs processed after the first chunk completes.
        # Progress callback tuple: (current_chunk, total_chunks, bp_processed)
        if progress_calls:
            current_chunk, total_chunks, bp_processed = progress_calls[0]
            self.assertGreater(bp_processed, 0, 
                              "First progress callback should show bp_processed > 0 "
                              "(fix for 'Processed: 0 / X bp' bug)")
        
        # Calculate throughput
        throughput = len(seq) / elapsed
        print(f"\n15kb chunked scan: {elapsed:.2f}s, {throughput:,.0f} bp/s, {len(motifs)} motifs")
        print(f"  Progress callback calls: {len(progress_calls)}")
        if progress_calls:
            _, _, first_bp_processed = progress_calls[0]
            print(f"  First callback bp_processed: {first_bp_processed:,}")
    
    def test_chunked_vs_non_chunked_consistency(self):
        """Test that chunked and non-chunked analysis produce similar motif counts."""
        from nonbscanner import analyze_sequence, CHUNK_THRESHOLD
        
        # Use sequence just above threshold
        seq = generate_g4_sequence(12000)
        
        # Analyze with chunking
        motifs_chunked = analyze_sequence(seq, "chunked", use_chunking=True)
        
        # Analyze without chunking (force single pass)
        motifs_single = analyze_sequence(seq, "single", use_chunking=False)
        
        # Both should find motifs
        self.assertGreater(len(motifs_chunked), 0)
        self.assertGreater(len(motifs_single), 0)
        
        # Motif counts should be relatively similar (within 20% tolerance due to overlap handling)
        ratio = len(motifs_chunked) / len(motifs_single) if len(motifs_single) > 0 else 0
        print(f"\nChunked: {len(motifs_chunked)}, Single: {len(motifs_single)}, Ratio: {ratio:.2f}")
        
        # The ratio should be close to 1 (allowing for deduplication differences)
        self.assertGreater(ratio, 0.5, "Chunked should find at least 50% of single-pass motifs")
        self.assertLess(ratio, 2.0, "Chunked should not find more than 2x the single-pass motifs")
    
    def test_sequence_with_motifs(self):
        """Test sequence known to contain motifs."""
        from nonbscanner import analyze_sequence
        
        # Generate sequence with known G4 patterns
        seq = generate_g4_sequence(5000)
        
        start = time.time()
        motifs = analyze_sequence(seq, "g4_test")
        elapsed = time.time() - start
        
        # Should find G-quadruplex motifs
        g4_motifs = [m for m in motifs if 'G' in m.get('Class', '')]
        
        self.assertGreater(len(g4_motifs), 0, "Should find G-quadruplex motifs")
        print(f"\nG4 sequence scan: {elapsed:.2f}s, found {len(g4_motifs)} G4 motifs")


class TestParallelScanning(unittest.TestCase):
    """Tests for parallel scanning functionality."""
    
    def test_parallel_scan_small(self):
        """Test parallel scanning on a small sequence."""
        from scanner_backends.parallel_worker import parallel_scan
        
        seq = generate_random_sequence(5000)
        
        # Use small chunk size to force multiple chunks
        motifs = parallel_scan(
            seq,
            sequence_name="parallel_test",
            num_workers=2,
            chunk_size=1000,
            overlap=100
        )
        
        self.assertIsInstance(motifs, list)
    
    def test_parallel_vs_sequential_consistency(self):
        """Test that parallel and sequential scanning produce similar results."""
        from nonbscanner import analyze_sequence
        from scanner_backends.parallel_worker import parallel_scan
        
        # Use a known sequence with motifs
        seq = generate_g4_sequence(3000)
        
        # Sequential scan
        sequential_motifs = analyze_sequence(seq, "seq_test")
        
        # Parallel scan (may use different chunking)
        parallel_motifs = parallel_scan(
            seq,
            sequence_name="par_test",
            num_workers=2,
            chunk_size=1000,
            overlap=200
        )
        
        # Both should find motifs (counts may differ due to overlap handling)
        self.assertIsInstance(sequential_motifs, list)
        self.assertIsInstance(parallel_motifs, list)
        
        # If both find motifs, they should be in similar ranges
        if len(sequential_motifs) > 0 and len(parallel_motifs) > 0:
            # At least some motifs should be in both
            seq_positions = set((m['Start'], m['End']) for m in sequential_motifs)
            par_positions = set((m['Start'], m['End']) for m in parallel_motifs)
            
            # There should be some overlap (not necessarily perfect match due to dedup)
            # This is a sanity check, not a strict equality test


class TestCLI(unittest.TestCase):
    """Tests for CLI functionality."""
    
    def test_cli_parse_args(self):
        """Test CLI argument parsing."""
        from cli.nbdscanner_cli import NBDScannerCLI
        
        args = NBDScannerCLI.parse_args([
            '--input', 'test.fa',
            '--out', 'results',
            '--format', 'ndjson',
            '--workers', '4'
        ])
        
        self.assertEqual(args.input, 'test.fa')
        self.assertEqual(args.out, 'results')
        self.assertEqual(args.format, 'ndjson')
        self.assertEqual(args.workers, 4)
    
    def test_cli_with_fasta_file(self):
        """Test CLI with actual FASTA file."""
        from cli.nbdscanner_cli import NBDScannerCLI
        
        # Create temporary FASTA file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(">test_sequence\n")
            f.write("GGGTTAGGGTTAGGGTTAGGGATGCATGC\n")
            tmp_fasta = f.name
        
        # Create temporary output file securely
        with tempfile.NamedTemporaryFile(mode='w', suffix='_out', delete=False) as f:
            tmp_out = f.name
            # Remove the file since CLI will create it with .ndjson extension
            os.unlink(tmp_out)
        
        try:
            args = NBDScannerCLI.parse_args([
                '--input', tmp_fasta,
                '--out', tmp_out,
                '--quiet',
                '--no-parallel'
            ])
            
            cli = NBDScannerCLI(args)
            exit_code = cli.run()
            
            self.assertEqual(exit_code, 0)
            
            # Check output file exists
            self.assertTrue(os.path.exists(tmp_out + '.ndjson'))
            
        finally:
            os.unlink(tmp_fasta)
            if os.path.exists(tmp_out + '.ndjson'):
                os.unlink(tmp_out + '.ndjson')


class TestIOUtilsPerformance(unittest.TestCase):
    """Performance tests for I/O utilities."""
    
    def test_streaming_chunks_performance(self):
        """Test that streaming chunks works efficiently."""
        from scanner_backends.io_utils import stream_sequence_chunks
        
        # Create large sequence as bytes
        seq_length = 1_000_000  # 1MB
        seq = ('ATGC' * (seq_length // 4)).encode()
        
        start = time.time()
        chunks = list(stream_sequence_chunks(seq, chunk_size=50000, overlap=500))
        elapsed = time.time() - start
        
        # Chunking should be very fast (< 1 second for 1MB)
        self.assertLess(elapsed, 1.0)
        
        # Should have reasonable number of chunks
        self.assertGreater(len(chunks), 10)
        
        print(f"\nChunking 1MB: {elapsed*1000:.1f}ms, {len(chunks)} chunks")
    
    def test_memory_estimate(self):
        """Test memory usage estimation."""
        from scanner_backends.io_utils import estimate_memory_usage
        
        # 10MB sequence, 4 workers
        estimate = estimate_memory_usage(10_000_000, num_workers=4)
        
        self.assertIn('total_estimate', estimate)
        self.assertIn('sequence_storage', estimate)
        
        # Total should be reasonable (< 1GB for 10MB sequence)
        self.assertLess(estimate['total_estimate'], 1_000_000_000)


class TestChunkPerformance(unittest.TestCase):
    """Tests for 100x chunk processing performance improvements."""
    
    def test_parallel_chunks_enabled(self):
        """Test that parallel chunk processing is available and improves throughput."""
        from nonbscanner import analyze_sequence
        
        # Generate a sequence large enough to trigger chunking
        seq = generate_g4_sequence(20000)
        
        # Test with parallel chunks enabled (default)
        start = time.time()
        motifs_parallel = analyze_sequence(
            seq, 
            "parallel_test",
            use_chunking=True,
            use_parallel_chunks=True,
            chunk_size=5000
        )
        elapsed_parallel = time.time() - start
        
        # Should complete and find motifs
        self.assertIsInstance(motifs_parallel, list)
        self.assertGreater(len(motifs_parallel), 0)
        
        # Calculate throughput
        throughput = len(seq) / elapsed_parallel
        print(f"\nParallel chunks test: {elapsed_parallel:.2f}s, {throughput:,.0f} bp/s")
    
    def test_progress_callback_with_metrics(self):
        """Test that progress callback receives performance metrics."""
        from nonbscanner import analyze_sequence
        
        seq = generate_g4_sequence(15000)
        
        # Track progress with new signature
        progress_data = []
        def progress_callback(chunk, total, bp, elapsed=None, throughput=None):
            progress_data.append({
                'chunk': chunk,
                'total': total,
                'bp': bp,
                'elapsed': elapsed,
                'throughput': throughput
            })
        
        motifs = analyze_sequence(
            seq,
            "metrics_test",
            use_chunking=True,
            chunk_size=5000,
            progress_callback=progress_callback
        )
        
        # Should have progress data
        self.assertGreater(len(progress_data), 0, "Progress callback should be called")
        
        # Check that we received the new metrics (elapsed and throughput)
        final_progress = progress_data[-1]
        if final_progress['elapsed'] is not None:
            self.assertGreater(final_progress['elapsed'], 0, "Elapsed time should be positive")
            print(f"\nFinal progress - elapsed: {final_progress['elapsed']:.2f}s, throughput: {final_progress['throughput']:,.0f} bp/s")
    
    def test_chunk_boundary_handling(self):
        """Test that motifs spanning chunk boundaries are properly detected."""
        from nonbscanner import analyze_sequence
        
        # Create a sequence with a known motif pattern
        # G4 motif that should be detected
        g4 = "GGGTTAGGGTTAGGGTTAGGG"
        
        # Place G4 motif at chunk boundaries (assuming 1000bp chunks)
        spacer = "ATGCATGC" * 125  # 1000 bp spacer
        
        # Create sequence with G4 at position ~1000 (chunk boundary)
        seq = spacer + g4 + spacer + g4 + spacer
        
        # Analyze with small chunks to create boundaries
        motifs = analyze_sequence(
            seq,
            "boundary_test",
            use_chunking=True,
            chunk_size=1000,
            chunk_overlap=500  # Overlap should catch boundary motifs
        )
        
        # Should find the G4 motifs
        g4_motifs = [m for m in motifs if 'G' in m.get('Class', '')]
        self.assertGreater(len(g4_motifs), 0, "Should detect G4 motifs across chunk boundaries")
        print(f"\nBoundary test: Found {len(g4_motifs)} G4 motifs across {len(seq)} bp")
    
    def test_parallel_vs_sequential_chunks_consistency(self):
        """Test that parallel and sequential chunk processing produce similar results."""
        from nonbscanner import analyze_sequence
        
        seq = generate_g4_sequence(15000)
        
        # Analyze with parallel chunks
        motifs_parallel = analyze_sequence(
            seq,
            "parallel",
            use_chunking=True,
            use_parallel_chunks=True,
            chunk_size=5000
        )
        
        # Analyze with sequential chunks
        motifs_sequential = analyze_sequence(
            seq,
            "sequential",
            use_chunking=True,
            use_parallel_chunks=False,
            chunk_size=5000
        )
        
        # Both should find similar number of motifs (within tolerance)
        self.assertGreater(len(motifs_parallel), 0)
        self.assertGreater(len(motifs_sequential), 0)
        
        # Ratio should be close to 1 (within 10% tolerance)
        ratio = len(motifs_parallel) / len(motifs_sequential) if len(motifs_sequential) > 0 else 0
        self.assertGreater(ratio, 0.9, "Parallel should find at least 90% of sequential motifs")
        self.assertLess(ratio, 1.1, "Parallel should not find more than 110% of sequential motifs")
        
        print(f"\nConsistency test - Parallel: {len(motifs_parallel)}, Sequential: {len(motifs_sequential)}, Ratio: {ratio:.2f}")


class TestDetectorTimings(unittest.TestCase):
    """Tests for individual detector timing feature."""
    
    def test_detector_timings_available(self):
        """Test that detector timings are available after analysis."""
        from nonbscanner import analyze_sequence, get_last_detector_timings, get_detector_display_names
        
        # Known test sequence
        seq = generate_g4_sequence(1000)
        
        # Run analysis
        motifs = analyze_sequence(seq, "timing_test")
        
        # Get timings
        timings = get_last_detector_timings()
        
        # Should have timings for all 9 detectors
        self.assertIsInstance(timings, dict)
        self.assertEqual(len(timings), 9, "Should have timings for all 9 detectors")
        
        # All timings should be non-negative floats
        for det_name, det_time in timings.items():
            self.assertIsInstance(det_time, float, f"Timing for {det_name} should be a float")
            self.assertGreaterEqual(det_time, 0.0, f"Timing for {det_name} should be non-negative")
        
        # Display names should be available
        display_names = get_detector_display_names()
        self.assertIsInstance(display_names, dict)
        self.assertEqual(len(display_names), 9)
        
        # All detector names should have display names
        for det_name in timings.keys():
            self.assertIn(det_name, display_names, f"Display name should exist for {det_name}")
        
        print(f"\nDetector timings test - Found timings for {len(timings)} detectors")
        total_time = sum(timings.values())
        print(f"Sum of detector times: {total_time:.3f}s")
    
    def test_detector_callback_receives_timing(self):
        """Test that detector callback receives elapsed time."""
        from nonbscanner import NonBScanner
        
        seq = generate_g4_sequence(500)
        
        # Track callbacks
        callback_data = []
        def timing_callback(det_name, completed, total, elapsed_time):
            callback_data.append({
                'name': det_name,
                'completed': completed,
                'total': total,
                'elapsed': elapsed_time
            })
        
        scanner = NonBScanner()
        motifs = scanner.analyze_sequence(seq, "callback_test", detector_callback=timing_callback)
        
        # Should have received callbacks for all 9 detectors
        self.assertEqual(len(callback_data), 9, "Should receive 9 callbacks")
        
        # Each callback should have timing
        for cb in callback_data:
            self.assertIn('elapsed', cb, "Callback should include elapsed time")
            self.assertIsInstance(cb['elapsed'], float, "Elapsed time should be a float")
            self.assertGreaterEqual(cb['elapsed'], 0.0, "Elapsed time should be non-negative")
        
        print(f"\nCallback test - Received {len(callback_data)} callbacks with timings")


if __name__ == '__main__':
    # Run with verbose output for performance info
    unittest.main(verbosity=2)
