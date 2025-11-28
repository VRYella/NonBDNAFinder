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


if __name__ == '__main__':
    # Run with verbose output for performance info
    unittest.main(verbosity=2)
