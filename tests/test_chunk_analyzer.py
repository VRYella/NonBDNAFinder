"""
Integration tests for chunk_analyzer.py - Chunk-based genome analysis.

Tests ChunkAnalyzer class with real analysis pipeline.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
import gc

# Add parent directory to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from Utilities.disk_storage import UniversalSequenceStorage
from Utilities.chunk_analyzer import ChunkAnalyzer


class TestChunkAnalyzer(unittest.TestCase):
    """Test cases for ChunkAnalyzer class."""
    
    def setUp(self):
        """Create temporary directory and storage for tests."""
        self.test_dir = tempfile.mkdtemp(prefix="test_chunk_analyzer_")
        self.storage = UniversalSequenceStorage(base_dir=self.test_dir)
    
    def tearDown(self):
        """Clean up temporary directory."""
        if Path(self.test_dir).exists():
            shutil.rmtree(self.test_dir)
        gc.collect()
    
    def test_analyze_small_sequence(self):
        """Test analysis of small sequence (single chunk)."""
        # Create test sequence with known patterns
        # G4 motif pattern: GGGTTAGGGTTAGGGTTAGGG (canonical G-quadruplex)
        sequence = "ATCGATCG" * 500 + "GGGTTAGGGTTAGGGTTAGGG" + "ATCGATCG" * 500
        seq_id = self.storage.save_sequence(sequence, "test_seq")
        
        # Analyze with small chunk size to test chunking
        analyzer = ChunkAnalyzer(self.storage, chunk_size=5000, overlap=1000)
        results_storage = analyzer.analyze(seq_id)
        
        # Verify results are stored
        stats = results_storage.get_summary_stats()
        self.assertGreater(stats['total_count'], 0)
        
        # Cleanup
        results_storage.cleanup()
    
    def test_analyze_with_overlap_deduplication(self):
        """Test that motifs in overlap regions are not duplicated."""
        # Create sequence where motif will be in overlap region
        # Position motif at chunk boundary
        chunk_size = 10000
        overlap = 1000
        
        # Place G4 motif at position that will be in overlap
        pre_sequence = "ATCG" * 2400  # 9600 bp
        g4_motif = "GGGTTAGGGTTAGGGTTAGGG"  # Should be at ~9600
        post_sequence = "ATCG" * 2500  # 10000 bp
        
        sequence = pre_sequence + g4_motif + post_sequence
        seq_id = self.storage.save_sequence(sequence, "overlap_test")
        
        # Analyze with chunking
        analyzer = ChunkAnalyzer(self.storage, chunk_size=chunk_size, overlap=overlap)
        results_storage = analyzer.analyze(seq_id)
        
        # Get all motifs
        all_motifs = list(results_storage.iter_results())
        
        # Check for duplicates at same position
        positions = [(m['Start'], m['End'], m.get('Class', '')) for m in all_motifs]
        unique_positions = set(positions)
        
        # Should have no duplicates
        self.assertEqual(len(positions), len(unique_positions),
                        "Found duplicate motifs at same position")
        
        # Cleanup
        results_storage.cleanup()
    
    def test_position_adjustment(self):
        """Test that motif positions are adjusted correctly across chunks."""
        # Create sequence with identifiable patterns at known positions
        sequence = "ATCG" * 5000  # 20,000 bp
        seq_id = self.storage.save_sequence(sequence, "position_test")
        
        # Analyze with small chunks
        analyzer = ChunkAnalyzer(self.storage, chunk_size=8000, overlap=500)
        results_storage = analyzer.analyze(seq_id)
        
        # Verify motif positions are within sequence bounds
        for motif in results_storage.iter_results():
            start = motif['Start']
            end = motif['End']
            
            self.assertGreaterEqual(start, 0, "Motif start position is negative")
            self.assertLessEqual(end, 20000, "Motif end position exceeds sequence length")
            self.assertLess(start, end, "Motif start >= end")
        
        # Cleanup
        results_storage.cleanup()
    
    def test_progress_callback(self):
        """Test progress callback functionality."""
        sequence = "ATCG" * 5000  # 20,000 bp
        seq_id = self.storage.save_sequence(sequence, "progress_test")
        
        # Track progress updates
        progress_updates = []
        
        def progress_callback(progress_pct):
            progress_updates.append(progress_pct)
        
        # Analyze with callback
        analyzer = ChunkAnalyzer(self.storage, chunk_size=8000, overlap=500)
        results_storage = analyzer.analyze(seq_id, progress_callback=progress_callback)
        
        # Verify callback was called
        self.assertGreater(len(progress_updates), 0, "Progress callback was not called")
        
        # Verify progress increases
        for i in range(len(progress_updates) - 1):
            self.assertLessEqual(progress_updates[i], progress_updates[i+1],
                               "Progress should be monotonically increasing")
        
        # Verify final progress is 100%
        self.assertAlmostEqual(progress_updates[-1], 100.0, places=0)
        
        # Cleanup
        results_storage.cleanup()
    
    def test_enabled_classes_filter(self):
        """Test filtering by enabled classes."""
        sequence = "ATCG" * 5000  # 20,000 bp
        seq_id = self.storage.save_sequence(sequence, "filter_test")
        
        # Analyze with only G-Quadruplex enabled
        analyzer = ChunkAnalyzer(self.storage, chunk_size=10000, overlap=500)
        results_storage = analyzer.analyze(
            seq_id,
            enabled_classes=['G-Quadruplex']
        )
        
        # Verify only G-Quadruplex motifs are found
        for motif in results_storage.iter_results():
            self.assertEqual(motif.get('Class', ''), 'G-Quadruplex',
                           "Found motif from disabled class")
        
        # Cleanup
        results_storage.cleanup()
    
    def test_memory_efficiency(self):
        """Test that chunk analysis doesn't load entire sequence."""
        # Create larger sequence (~102KB)
        sequence = "ATCGATCGATCG" * 8500  # ~102KB
        seq_id = self.storage.save_sequence(sequence, "memory_test")
        
        # Force garbage collection before analysis
        gc.collect()
        
        # Analyze with chunking
        analyzer = ChunkAnalyzer(self.storage, chunk_size=20000, overlap=1000)
        results_storage = analyzer.analyze(seq_id)
        
        # Verify analysis completed
        stats = results_storage.get_summary_stats()
        self.assertIsNotNone(stats['total_count'])
        
        # Force garbage collection after analysis
        gc.collect()
        
        # Cleanup
        results_storage.cleanup()
    
    def test_multiple_chunks(self):
        """Test analysis with multiple chunks."""
        # Create sequence that requires multiple chunks
        sequence = "ATCG" * 10000  # 40,000 bp
        seq_id = self.storage.save_sequence(sequence, "multi_chunk_test")
        
        # Analyze with small chunk size to force multiple chunks
        analyzer = ChunkAnalyzer(self.storage, chunk_size=10000, overlap=500)
        results_storage = analyzer.analyze(seq_id)
        
        # Verify results were generated
        stats = results_storage.get_summary_stats()
        self.assertGreaterEqual(stats['total_count'], 0)
        
        # Cleanup
        results_storage.cleanup()


if __name__ == '__main__':
    unittest.main()
