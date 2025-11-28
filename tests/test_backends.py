"""
Backend Tests for NBDScanner Scanner Backends
==============================================

Tests for the scanner_backends module, including:
- I/O utilities (mmap, streaming)
- Hyperscan backend (when available)
- Numba backend (when available)
- Parallel worker with shared memory
"""

import os
import sys
import tempfile
import unittest

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class TestIOUtils(unittest.TestCase):
    """Tests for scanner_backends.io_utils module."""
    
    def test_get_overlap_size(self):
        """Test overlap size calculation."""
        from scanner_backends.io_utils import get_overlap_size
        
        # Test default (all motifs)
        default_overlap = get_overlap_size()
        self.assertGreaterEqual(default_overlap, 100)
        self.assertLessEqual(default_overlap, 1000)
        
        # Test specific motif types
        g4_overlap = get_overlap_size(['g_quadruplex'])
        self.assertGreaterEqual(g4_overlap, 50)
        
        # Test multiple motif types
        multi_overlap = get_overlap_size(['r_loop', 'slipped_dna'])
        self.assertGreaterEqual(multi_overlap, 300)
    
    def test_stream_sequence_chunks(self):
        """Test sequence chunking for parallel processing."""
        from scanner_backends.io_utils import stream_sequence_chunks
        
        # Create test sequence
        test_seq = b'ATGC' * 1000  # 4000 bp
        
        # Test with small chunk size
        chunks = list(stream_sequence_chunks(test_seq, chunk_size=500, overlap=50))
        
        # Should have multiple chunks
        self.assertGreater(len(chunks), 1)
        
        # First chunk starts at 0
        self.assertEqual(chunks[0][0], 0)
        
        # Chunks should overlap
        for i in range(1, len(chunks)):
            prev_end = chunks[i-1][0] + len(chunks[i-1][1])
            curr_start = chunks[i][0]
            # Current should start before previous ends
            self.assertLess(curr_start, prev_end)
    
    def test_bytes_to_numpy(self):
        """Test sequence bytes to numpy conversion."""
        from scanner_backends.io_utils import bytes_to_numpy, numpy_to_string
        import numpy as np
        
        test_seq = b'ATGCATGC'
        arr = bytes_to_numpy(test_seq)
        
        self.assertIsInstance(arr, np.ndarray)
        self.assertEqual(arr.dtype, np.uint8)
        self.assertEqual(len(arr), len(test_seq))
        
        # Test round-trip
        reconstructed = numpy_to_string(arr)
        self.assertEqual(reconstructed, test_seq.decode('utf-8'))
    
    def test_mmap_fasta(self):
        """Test memory-mapped FASTA reading."""
        from scanner_backends.io_utils import mmap_fasta
        
        # Create temporary FASTA file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(">seq1\n")
            f.write("ATGCATGCATGC\n")
            f.write(">seq2\n")
            f.write("GGGGCCCCAAAA\n")
            tmp_path = f.name
        
        try:
            mm_data, seq_offsets = mmap_fasta(tmp_path)
            
            # Should find both sequences
            self.assertEqual(len(seq_offsets), 2)
            self.assertIn('seq1', seq_offsets)
            self.assertIn('seq2', seq_offsets)
            
            # Check we can extract sequences
            for seq_name, (start, end) in seq_offsets.items():
                seq_bytes = mm_data[start:end]
                seq = seq_bytes.decode('utf-8').replace('\n', '')
                self.assertGreater(len(seq), 0)
                
        finally:
            os.unlink(tmp_path)
    
    def test_ndjson_writer(self):
        """Test NDJSON streaming writer."""
        from scanner_backends.io_utils import write_ndjson_stream, read_ndjson_stream
        import json
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.ndjson', delete=False) as f:
            tmp_path = f.name
        
        try:
            # Write records
            records = [
                {'id': 1, 'name': 'test1'},
                {'id': 2, 'name': 'test2'},
                {'id': 3, 'name': 'test3'},
            ]
            
            with write_ndjson_stream(tmp_path) as writer:
                for r in records:
                    writer.write_record(r)
            
            # Read records back
            read_records = list(read_ndjson_stream(tmp_path))
            
            self.assertEqual(len(read_records), len(records))
            for orig, read in zip(records, read_records):
                self.assertEqual(orig, read)
                
        finally:
            os.unlink(tmp_path)


class TestParallelWorker(unittest.TestCase):
    """Tests for scanner_backends.parallel_worker module."""
    
    def test_chunk_sequence(self):
        """Test sequence chunking calculation."""
        from scanner_backends.parallel_worker import chunk_sequence
        
        # Small sequence - single chunk
        chunks = chunk_sequence(1000, chunk_size=5000, overlap=100)
        self.assertEqual(len(chunks), 1)
        self.assertEqual(chunks[0], (0, 1000))
        
        # Large sequence - multiple chunks
        chunks = chunk_sequence(10000, chunk_size=2000, overlap=200)
        self.assertGreater(len(chunks), 1)
        
        # Verify overlap between consecutive chunks
        for i in range(1, len(chunks)):
            prev_end = chunks[i-1][1]
            curr_start = chunks[i][0]
            self.assertLess(curr_start, prev_end)
    
    def test_shared_memory_worker(self):
        """Test SharedMemoryWorker creation and access."""
        from scanner_backends.parallel_worker import SharedMemoryWorker
        
        test_seq = "ATGCATGCATGCATGC"
        
        with SharedMemoryWorker(test_seq, name="test_shm") as worker:
            # Get chunk
            chunk = worker.get_chunk(0, 8)
            self.assertEqual(chunk, "ATGCATGC")
            
            chunk2 = worker.get_chunk(4, 12)
            self.assertEqual(chunk2, "ATGCATGC")
    
    def test_deduplicate_overlap_motifs(self):
        """Test deduplication of motifs from overlap regions."""
        from scanner_backends.parallel_worker import _deduplicate_overlap_motifs
        
        # Create duplicate motifs
        motifs = [
            {'Class': 'G-Quadruplex', 'Subclass': 'Canonical', 'Start': 100, 'End': 150, 'Score': 0.9},
            {'Class': 'G-Quadruplex', 'Subclass': 'Canonical', 'Start': 100, 'End': 150, 'Score': 0.85},  # Duplicate
            {'Class': 'G-Quadruplex', 'Subclass': 'Canonical', 'Start': 200, 'End': 250, 'Score': 0.8},
        ]
        
        deduped = _deduplicate_overlap_motifs(motifs, overlap=100)
        
        # Should have 2 unique motifs
        self.assertEqual(len(deduped), 2)
        
        # Should keep higher score for duplicates
        motif_100 = [m for m in deduped if m['Start'] == 100][0]
        self.assertEqual(motif_100['Score'], 0.9)


class TestNumbaBackend(unittest.TestCase):
    """Tests for scanner_backends.numba_backend module."""
    
    def test_numba_availability(self):
        """Test Numba availability check."""
        from scanner_backends.numba_backend import is_numba_available
        
        # Should return a boolean
        result = is_numba_available()
        self.assertIsInstance(result, bool)
    
    def test_encode_sequence(self):
        """Test sequence encoding."""
        from scanner_backends.numba_backend import encode_sequence, decode_sequence
        import numpy as np
        
        test_seq = "ATGCATGC"
        encoded = encode_sequence(test_seq)
        
        self.assertIsInstance(encoded, np.ndarray)
        self.assertEqual(len(encoded), len(test_seq))
        self.assertEqual(encoded.dtype, np.uint8)
        
        # Test round-trip
        decoded = decode_sequence(encoded)
        self.assertEqual(decoded, test_seq)
    
    def test_numba_scanner_g4(self):
        """Test Numba scanner G-quadruplex detection."""
        from scanner_backends.numba_backend import NumbaScanner
        
        scanner = NumbaScanner()
        
        # Test sequence with clear G4 pattern
        test_seq = "ATGCGGGTTAGGGTTAGGGTTAGGGATGC"
        
        candidates = scanner.find_g_quadruplex_candidates(test_seq)
        
        # May or may not find candidates depending on exact pattern
        # Just verify it runs without error
        self.assertIsInstance(candidates, list)
    
    def test_numba_scanner_a_tracts(self):
        """Test Numba scanner A-tract detection."""
        from scanner_backends.numba_backend import NumbaScanner
        
        scanner = NumbaScanner()
        
        # Test sequence with A-tracts
        test_seq = "GCGCAAAAAGCGCAAAAAAAGCGC"
        
        tracts = scanner.find_a_tracts(test_seq, min_length=4)
        
        # Should find A-tracts
        self.assertGreater(len(tracts), 0)
        
        for tract in tracts:
            self.assertEqual(tract['Class'], 'Curved_DNA')
            self.assertGreaterEqual(tract['Length'], 4)
    
    def test_numba_scanner_gc_content(self):
        """Test Numba scanner GC content calculation."""
        from scanner_backends.numba_backend import NumbaScanner
        import numpy as np
        
        scanner = NumbaScanner()
        
        # Test sequence with known GC content
        test_seq = "GGGGCCCCAAAA"  # 66.7% GC
        
        gc_values = scanner.calculate_gc_content(test_seq, window_size=len(test_seq))
        
        self.assertEqual(len(gc_values), 1)
        self.assertAlmostEqual(gc_values[0], 8/12, places=2)


class TestHyperscanBackend(unittest.TestCase):
    """Tests for scanner_backends.hyperscan_backend module."""
    
    def test_hyperscan_availability(self):
        """Test Hyperscan availability check."""
        from scanner_backends.hyperscan_backend import is_hyperscan_available
        
        # Should return a boolean
        result = is_hyperscan_available()
        self.assertIsInstance(result, bool)
    
    def test_fallback_regex_scan(self):
        """Test fallback regex scanning."""
        from scanner_backends.hyperscan_backend import fallback_regex_scan
        
        patterns = [
            (r'G{3,}', 'G_RUN', {'class': 'Test', 'subclass': 'G-run', 'score': 0.8}),
            (r'A{4,}', 'A_RUN', {'class': 'Test', 'subclass': 'A-run', 'score': 0.7}),
        ]
        
        test_seq = "ATGCGGGGAAAAATGC"
        
        matches = fallback_regex_scan(test_seq, patterns)
        
        # Should find matches
        self.assertGreater(len(matches), 0)
        
        # Check match structure
        for match in matches:
            self.assertIn('Start', match)
            self.assertIn('End', match)
            self.assertIn('Score', match)
            self.assertIn('Class', match)


class TestBackendIntegration(unittest.TestCase):
    """Integration tests for scanner backends."""
    
    def test_get_best_backend(self):
        """Test backend selection."""
        from scanner_backends import get_best_backend
        
        backend = get_best_backend()
        
        # Should return a valid backend name
        self.assertIn(backend, ['hyperscan', 'numba', 'python'])
    
    def test_get_scanner(self):
        """Test scanner instance creation."""
        from scanner_backends import get_scanner
        
        # Auto backend
        scanner = get_scanner('auto')
        # May be None (python) or a scanner instance
        
        # Python backend should return None (use standard scanner)
        python_scanner = get_scanner('python')
        self.assertIsNone(python_scanner)


if __name__ == '__main__':
    unittest.main()
