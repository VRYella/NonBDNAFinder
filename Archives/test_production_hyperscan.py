#!/usr/bin/env python3
"""
Comprehensive tests for production_hyperscan_streamlit module.

Tests cover:
- Basic scanner functionality
- Chunking and overlap handling
- Parallel processing
- Memory management
- Progress tracking
- Error handling
"""

import pytest
import tempfile
import os
import time
from unittest.mock import Mock, patch
from typing import List

# Test imports
from production_hyperscan_streamlit import (
    HyperscanStreamlitScanner,
    ScanConfig,
    ScanResult,
    _worker_init,
    _worker_scan_chunk,
    HYPERSCAN_AVAILABLE
)


class TestScanResult:
    """Test ScanResult dataclass"""
    
    def test_scan_result_creation(self):
        """Test ScanResult creation and methods"""
        result = ScanResult(
            pattern_id=1,
            start=100,
            end=150,
            pattern="GGGG",
            sequence_name="test_seq",
            score=0.8
        )
        
        assert result.pattern_id == 1
        assert result.start == 100
        assert result.end == 150
        assert result.pattern == "GGGG"
        assert result.sequence_name == "test_seq"
        assert result.score == 0.8
    
    def test_to_dict(self):
        """Test ScanResult to_dict method"""
        result = ScanResult(
            pattern_id=1,
            start=100,
            end=150,
            pattern="GGGG"
        )
        
        expected = {
            'pattern_id': 1,
            'start': 100,
            'end': 150,
            'length': 50,
            'pattern': 'GGGG',
            'sequence_name': '',
            'score': 0.0
        }
        
        assert result.to_dict() == expected


class TestScanConfig:
    """Test ScanConfig dataclass"""
    
    def test_default_config(self):
        """Test default configuration values"""
        config = ScanConfig()
        
        assert config.chunk_size == 2_000_000
        assert config.overlap == 500
        assert config.max_workers >= 1
        assert config.timeout == 300
        assert config.stream_results == True
        assert config.output_format == 'csv'
        assert config.incremental_save == True
        assert config.max_memory_mb == 1000
    
    def test_custom_config(self):
        """Test custom configuration"""
        config = ScanConfig(
            chunk_size=1000000,
            overlap=200,
            max_workers=2,
            timeout=600
        )
        
        assert config.chunk_size == 1000000
        assert config.overlap == 200
        assert config.max_workers == 2
        assert config.timeout == 600


@pytest.mark.skipif(not HYPERSCAN_AVAILABLE, reason="Hyperscan not available")
class TestHyperscanStreamlitScanner:
    """Test main scanner class"""
    
    @pytest.fixture
    def sample_patterns(self):
        """Sample patterns for testing"""
        return [
            r'GGGG',
            r'CCCC',
            r'AAAA',
            r'TTTT'
        ]
    
    @pytest.fixture
    def test_sequence(self):
        """Test DNA sequence with known patterns"""
        return "ATCGATCGGGGGCCCCTTTAAAAAAAGCTAGCT" * 10  # ~320bp
    
    @pytest.fixture
    def large_sequence(self):
        """Large test sequence for chunking tests"""
        base = "ATCGATCGGGGGCCCCTTTAAAAAAAGCTAGCT"
        return base * 100000  # ~3.2MB sequence
    
    @pytest.fixture
    def scanner(self, sample_patterns):
        """Create scanner for testing"""
        return HyperscanStreamlitScanner(sample_patterns)
    
    def test_scanner_initialization(self, sample_patterns):
        """Test scanner initialization"""
        scanner = HyperscanStreamlitScanner(sample_patterns)
        
        assert scanner.patterns == sample_patterns
        assert scanner._db is not None
        assert len(scanner._pattern_map) == len(sample_patterns)
    
    def test_scanner_initialization_no_patterns(self):
        """Test scanner initialization with no patterns"""
        with pytest.raises(ValueError, match="No patterns provided"):
            HyperscanStreamlitScanner([])
    
    def test_scan_sequence_basic(self, scanner, test_sequence):
        """Test basic sequence scanning"""
        results = scanner.scan_sequence(test_sequence, "test_seq")
        
        assert isinstance(results, list)
        assert all(isinstance(r, ScanResult) for r in results)
        assert all(r.sequence_name == "test_seq" for r in results)
        
        # Should find multiple GGGG, CCCC, AAAA patterns
        pattern_ids = [r.pattern_id for r in results]
        assert len(set(pattern_ids)) > 1  # Multiple different patterns found
    
    def test_scan_empty_sequence(self, scanner):
        """Test scanning empty sequence"""
        results = scanner.scan_sequence("", "empty")
        assert results == []
    
    def test_chunking_functionality(self, scanner):
        """Test sequence chunking"""
        sequence = "A" * 1000000  # 1MB sequence
        chunks = scanner._create_chunks(sequence)
        
        assert len(chunks) > 1  # Should be chunked
        
        # Check overlap
        if len(chunks) > 1:
            chunk1_end = chunks[0][1] + len(chunks[0][0])
            chunk2_start = chunks[1][1]
            overlap = chunk1_end - chunk2_start
            assert overlap == scanner.config.overlap
    
    def test_deduplicate_results(self, scanner):
        """Test result deduplication"""
        # Create duplicate results
        results = [
            ScanResult(1, 100, 150, "GGGG", "seq1"),
            ScanResult(1, 100, 150, "GGGG", "seq1"),  # Duplicate
            ScanResult(2, 200, 250, "CCCC", "seq1"),
            ScanResult(1, 100, 150, "GGGG", "seq2"),  # Different sequence
        ]
        
        deduplicated = scanner._deduplicate_results(results)
        
        assert len(deduplicated) == 3  # One duplicate removed
        
        # Check all remaining results are unique
        seen = set()
        for result in deduplicated:
            key = (result.sequence_name, result.start, result.end, result.pattern_id)
            assert key not in seen
            seen.add(key)
    
    def test_scan_chunked_sequence(self, scanner, large_sequence):
        """Test scanning of large sequence with chunking"""
        # Use smaller chunk size for testing
        scanner.config.chunk_size = 100000
        scanner.config.overlap = 1000
        
        progress_calls = []
        def progress_callback(progress, message):
            progress_calls.append((progress, message))
        
        results = scanner._scan_chunked_sequence(
            large_sequence, 
            "large_seq", 
            progress_callback
        )
        
        assert isinstance(results, list)
        assert len(progress_calls) > 0  # Progress was tracked
        assert all(0 <= call[0] <= 1 for call in progress_calls)  # Valid progress values
    
    def create_test_fasta_file(self, sequences):
        """Helper to create test FASTA file"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            for seq_id, sequence in sequences:
                f.write(f">{seq_id}\n{sequence}\n")
            return f.name
    
    def test_scan_fasta_file_single_sequence(self, scanner, test_sequence):
        """Test scanning FASTA file with single sequence"""
        fasta_path = self.create_test_fasta_file([("seq1", test_sequence)])
        
        try:
            results = scanner.scan_fasta_file(fasta_path)
            
            assert 'files' in results
            assert 'stats' in results
            assert results['stats']['sequences_scanned'] == 1
            assert results['stats']['total_matches'] > 0
            
        finally:
            os.unlink(fasta_path)
    
    def test_scan_fasta_file_multiple_sequences(self, scanner, test_sequence):
        """Test scanning FASTA file with multiple sequences"""
        sequences = [
            ("seq1", test_sequence),
            ("seq2", test_sequence[:100]),
            ("seq3", test_sequence[50:])
        ]
        fasta_path = self.create_test_fasta_file(sequences)
        
        try:
            # Use only 1 worker to test sequential processing
            scanner.config.max_workers = 1
            results = scanner.scan_fasta_file(fasta_path)
            
            assert results['stats']['sequences_scanned'] == 3
            
        finally:
            os.unlink(fasta_path)
    
    def test_memory_management(self, scanner):
        """Test memory management with incremental saves"""
        # Create scanner with small memory limit
        scanner.config.incremental_save = True
        scanner.config.max_memory_mb = 1  # Very small limit
        
        # This should trigger incremental saves
        large_results = [
            ScanResult(i % 4, i * 10, (i + 1) * 10, f"pattern_{i % 4}")
            for i in range(20000)  # Many results
        ]
        
        output_prefix = "/tmp/test_memory"
        scanner._save_incremental_results(large_results, output_prefix)
        
        # Check that file was created
        csv_path = f"{output_prefix}_incremental.csv"
        assert os.path.exists(csv_path)
        
        # Cleanup
        try:
            os.unlink(csv_path)
        except Exception:
            pass


class TestWorkerFunctions:
    """Test multiprocessing worker functions"""
    
    @pytest.mark.skipif(not HYPERSCAN_AVAILABLE, reason="Hyperscan not available")
    def test_worker_init(self):
        """Test worker initialization"""
        patterns = ['GGGG', 'CCCC']
        
        # Test initialization
        _worker_init(patterns)
        
        # Check global variables are set
        from production_hyperscan_streamlit import _worker_db, _worker_patterns
        assert _worker_db is not None
        assert _worker_patterns == patterns
    
    @pytest.mark.skipif(not HYPERSCAN_AVAILABLE, reason="Hyperscan not available")
    def test_worker_scan_chunk(self):
        """Test worker chunk scanning"""
        patterns = ['GGGG', 'CCCC']
        _worker_init(patterns)
        
        work_item = {
            'sequence': 'ATCGGGGGCCCCTAG',
            'sequence_name': 'test_chunk_1',
            'offset': 0,
            'original_name': 'test_seq'
        }
        
        results = _worker_scan_chunk(work_item)
        
        assert isinstance(results, list)
        assert all(isinstance(r, ScanResult) for r in results)
        assert all(r.sequence_name == 'test_seq' for r in results)  # Uses original name


class TestErrorHandling:
    """Test error handling and edge cases"""
    
    def test_scanner_without_hyperscan(self):
        """Test scanner behavior when Hyperscan is not available"""
        with patch('production_hyperscan_streamlit.HYPERSCAN_AVAILABLE', False):
            with pytest.raises(RuntimeError, match="Hyperscan is not available"):
                HyperscanStreamlitScanner(['GGGG'])
    
    @pytest.mark.skipif(not HYPERSCAN_AVAILABLE, reason="Hyperscan not available")
    def test_invalid_fasta_file(self, sample_patterns):
        """Test handling of invalid FASTA file"""
        scanner = HyperscanStreamlitScanner(sample_patterns)
        
        with pytest.raises(FileNotFoundError):
            scanner.scan_fasta_file("/nonexistent/file.fa")
    
    @pytest.mark.skipif(not HYPERSCAN_AVAILABLE, reason="Hyperscan not available")
    def test_empty_fasta_file(self, sample_patterns):
        """Test handling of empty FASTA file"""
        scanner = HyperscanStreamlitScanner(sample_patterns)
        
        # Create empty file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write("")  # Empty file
            empty_path = f.name
        
        try:
            with pytest.raises(ValueError, match="No sequences found"):
                scanner.scan_fasta_file(empty_path)
        finally:
            os.unlink(empty_path)


class TestStreamlitIntegration:
    """Test Streamlit-specific functionality"""
    
    def test_streamlit_imports(self):
        """Test that Streamlit imports are handled gracefully"""
        # This should not raise an error even if Streamlit is not available
        from production_hyperscan_streamlit import STREAMLIT_AVAILABLE
        assert isinstance(STREAMLIT_AVAILABLE, bool)
    
    @pytest.mark.skipif(not HYPERSCAN_AVAILABLE, reason="Hyperscan not available")
    def test_progress_callback(self, sample_patterns):
        """Test progress callback functionality"""
        scanner = HyperscanStreamlitScanner(sample_patterns)
        
        progress_calls = []
        def progress_callback(progress, message):
            progress_calls.append((progress, message))
            assert 0 <= progress <= 1
            assert isinstance(message, str)
        
        # Test with small sequence that will use chunking
        scanner.config.chunk_size = 50  # Very small chunks
        sequence = "A" * 200  # Force chunking
        
        results = scanner.scan_sequence(sequence, "test", progress_callback)
        
        # Should have received progress updates
        assert len(progress_calls) > 0


class TestPerformance:
    """Performance-related tests"""
    
    @pytest.mark.skipif(not HYPERSCAN_AVAILABLE, reason="Hyperscan not available")
    def test_large_sequence_performance(self, sample_patterns):
        """Test performance with large sequences"""
        scanner = HyperscanStreamlitScanner(sample_patterns)
        
        # Create large sequence
        large_seq = "ATCGATCGGGGGCCCCTTTAAAAAAAGCTAGCT" * 10000  # ~320KB
        
        start_time = time.time()
        results = scanner.scan_sequence(large_seq, "performance_test")
        end_time = time.time()
        
        # Should complete reasonably quickly
        assert end_time - start_time < 30  # Less than 30 seconds
        assert isinstance(results, list)
    
    @pytest.mark.skipif(not HYPERSCAN_AVAILABLE, reason="Hyperscan not available")
    def test_many_patterns_performance(self):
        """Test performance with many patterns"""
        # Create many patterns
        patterns = [f"{'A' * i}{'T' * (10-i)}" for i in range(1, 10)]
        patterns.extend([f"{'G' * i}{'C' * (10-i)}" for i in range(1, 10)])
        
        start_time = time.time()
        scanner = HyperscanStreamlitScanner(patterns)
        compile_time = time.time() - start_time
        
        # Database compilation should be reasonable
        assert compile_time < 10  # Less than 10 seconds
        
        # Test scanning
        test_seq = "ATCGATCG" * 1000
        start_time = time.time()
        results = scanner.scan_sequence(test_seq, "many_patterns_test")
        scan_time = time.time() - start_time
        
        assert scan_time < 5  # Less than 5 seconds
        assert isinstance(results, list)


if __name__ == "__main__":
    # Run tests if executed directly
    pytest.main([__file__, "-v"])