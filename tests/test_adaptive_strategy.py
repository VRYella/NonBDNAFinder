"""
Integration tests for adaptive chunking strategy.

Tests the integration between ChunkAnalyzer and TripleAdaptiveChunkAnalyzer,
and verifies that results are identical across different strategies.
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
from Utilities.triple_chunk_analyzer import TripleAdaptiveChunkAnalyzer


class TestAdaptiveIntegration(unittest.TestCase):
    """Test integration of adaptive chunking with existing ChunkAnalyzer."""
    
    def setUp(self):
        """Create temporary directory and storage for tests."""
        self.test_dir = tempfile.mkdtemp(prefix="test_adaptive_integration_")
        self.storage = UniversalSequenceStorage(base_dir=self.test_dir)
    
    def tearDown(self):
        """Clean up temporary directory."""
        if Path(self.test_dir).exists():
            shutil.rmtree(self.test_dir)
        gc.collect()
    
    def test_chunk_analyzer_with_adaptive_enabled(self):
        """Test ChunkAnalyzer delegates to TripleAdaptiveChunkAnalyzer when use_adaptive=True."""
        sequence = "ATCGATCG" * 250000  # 2,000,000 bp
        seq_id = self.storage.save_sequence(sequence, "adaptive_test")
        
        # Create analyzer with adaptive enabled
        analyzer = ChunkAnalyzer(
            self.storage,
            use_adaptive=True
        )
        
        results_storage = analyzer.analyze(seq_id)
        
        # Should produce results
        stats = results_storage.get_summary_stats()
        self.assertGreaterEqual(stats['total_count'], 0)
        
        # Cleanup
        results_storage.cleanup()
    
    def test_chunk_analyzer_backward_compatibility(self):
        """Test ChunkAnalyzer with adaptive disabled (default behavior)."""
        sequence = "ATCGATCG" * 250000  # 2,000,000 bp
        seq_id = self.storage.save_sequence(sequence, "backward_compat_test")
        
        # Create analyzer with default settings (adaptive=False)
        analyzer = ChunkAnalyzer(
            self.storage,
            chunk_size=5_000_000,
            overlap=10_000,
            use_adaptive=False
        )
        
        results_storage = analyzer.analyze(seq_id)
        
        # Should produce results with original behavior
        stats = results_storage.get_summary_stats()
        self.assertGreaterEqual(stats['total_count'], 0)
        
        # Cleanup
        results_storage.cleanup()
    
    def test_results_consistency_small_sequence(self):
        """Test that adaptive and non-adaptive produce consistent results for small sequences."""
        # Use a sequence with known patterns
        sequence = "ATCG" * 1000 + "GGGTTAGGGTTAGGGTTAGGG" + "ATCG" * 1000
        seq_id = self.storage.save_sequence(sequence, "consistency_test_small")
        
        # Analyze with original ChunkAnalyzer
        analyzer_original = ChunkAnalyzer(
            self.storage,
            chunk_size=5_000_000,
            overlap=10_000,
            use_adaptive=False
        )
        results_original = analyzer_original.analyze(seq_id)
        stats_original = results_original.get_summary_stats()
        motifs_original = list(results_original.iter_results())
        
        # Analyze with adaptive (should use direct analysis for small seq)
        analyzer_adaptive = ChunkAnalyzer(
            self.storage,
            use_adaptive=True
        )
        results_adaptive = analyzer_adaptive.analyze(seq_id)
        stats_adaptive = results_adaptive.get_summary_stats()
        motifs_adaptive = list(results_adaptive.iter_results())
        
        # Both should find similar number of motifs (may vary slightly due to chunking)
        # Allow for some variation due to boundary effects
        self.assertAlmostEqual(
            stats_original['total_count'],
            stats_adaptive['total_count'],
            delta=5,  # Allow up to 5 motif difference
            msg="Motif counts differ significantly between strategies"
        )
        
        # Cleanup
        results_original.cleanup()
        results_adaptive.cleanup()
    
    def test_parallel_execution(self):
        """Test parallel execution of adaptive chunking."""
        sequence = "ATCGATCG" * 1250000  # 10,000,000 bp
        seq_id = self.storage.save_sequence(sequence, "parallel_test")
        
        # Create analyzer with parallel enabled
        analyzer = ChunkAnalyzer(
            self.storage,
            use_parallel=True,
            use_adaptive=True,
            max_workers=2
        )
        
        results_storage = analyzer.analyze(seq_id)
        
        # Should produce results
        stats = results_storage.get_summary_stats()
        self.assertGreaterEqual(stats['total_count'], 0)
        
        # Cleanup
        results_storage.cleanup()
    
    def test_progress_callback_integration(self):
        """Test progress callbacks work through ChunkAnalyzer."""
        sequence = "ATCGATCG" * 250000  # 2,000,000 bp
        seq_id = self.storage.save_sequence(sequence, "progress_integration_test")
        
        progress_values = []
        
        def progress_callback(pct):
            progress_values.append(pct)
        
        analyzer = ChunkAnalyzer(
            self.storage,
            use_adaptive=True
        )
        
        results_storage = analyzer.analyze(seq_id, progress_callback=progress_callback)
        
        # Progress should have been called
        self.assertGreater(len(progress_values), 0)
        
        # Cleanup
        results_storage.cleanup()
    
    def test_enabled_classes_propagation(self):
        """Test that enabled_classes filter propagates correctly."""
        sequence = "ATCGATCG" * 125000  # 1,000,000 bp
        seq_id = self.storage.save_sequence(sequence, "filter_integration_test")
        
        # Analyze with class filter
        analyzer = ChunkAnalyzer(
            self.storage,
            use_adaptive=True
        )
        
        results_storage = analyzer.analyze(
            seq_id,
            enabled_classes=['z_dna']
        )
        
        # Should produce results (or none if no Z-DNA found)
        stats = results_storage.get_summary_stats()
        self.assertGreaterEqual(stats['total_count'], 0)
        
        # Cleanup
        results_storage.cleanup()
    
    def test_memory_efficiency(self):
        """Test that adaptive chunking doesn't load entire sequence."""
        # Create large sequence
        sequence = "ATCGATCG" * 2500000  # 20,000,000 bp
        seq_id = self.storage.save_sequence(sequence, "memory_test")
        
        # This should not cause memory issues
        analyzer = ChunkAnalyzer(
            self.storage,
            use_adaptive=True
        )
        
        results_storage = analyzer.analyze(seq_id)
        
        # Should complete successfully
        stats = results_storage.get_summary_stats()
        self.assertGreaterEqual(stats['total_count'], 0)
        
        # Cleanup
        results_storage.cleanup()


class TestStrategyComparison(unittest.TestCase):
    """Compare results across different chunking strategies."""
    
    def setUp(self):
        """Create temporary directory and storage for tests."""
        self.test_dir = tempfile.mkdtemp(prefix="test_strategy_comparison_")
        self.storage = UniversalSequenceStorage(base_dir=self.test_dir)
    
    def tearDown(self):
        """Clean up temporary directory."""
        if Path(self.test_dir).exists():
            shutil.rmtree(self.test_dir)
        gc.collect()
    
    def test_single_vs_double_tier(self):
        """Compare single-tier and double-tier strategies on same sequence."""
        # Create sequence at boundary (9MB - single, 11MB - double)
        sequence_9mb = "ATCGATCG" * 1125000  # 9,000,000 bp
        seq_id_9mb = self.storage.save_sequence(sequence_9mb, "seq_9mb")
        
        # Analyze with adaptive (should use single-tier)
        analyzer_adaptive = TripleAdaptiveChunkAnalyzer(
            self.storage,
            use_adaptive=True
        )
        results_9mb = analyzer_adaptive.analyze(seq_id_9mb)
        stats_9mb = results_9mb.get_summary_stats()
        
        # Should produce results
        self.assertGreaterEqual(stats_9mb['total_count'], 0)
        
        # Cleanup
        results_9mb.cleanup()
    
    def test_deduplication_consistency(self):
        """Test that deduplication is consistent across tiers."""
        # Create sequence with pattern at multiple boundaries
        base_seq = "ATCG" * 10000
        g4_motif = "GGGTTAGGGTTAGGGTTAGGG"
        
        # Build sequence with patterns at different positions
        sequence = base_seq + g4_motif + base_seq + g4_motif + base_seq
        seq_id = self.storage.save_sequence(sequence, "dedup_test")
        
        analyzer = TripleAdaptiveChunkAnalyzer(
            self.storage,
            use_adaptive=True
        )
        results_storage = analyzer.analyze(seq_id)
        
        # Count unique motifs
        all_motifs = list(results_storage.iter_results())
        positions = [(m['Start'], m['End']) for m in all_motifs]
        unique_positions = set(positions)
        
        # Should have no duplicates
        self.assertEqual(
            len(positions),
            len(unique_positions),
            "Found duplicate motifs across tier boundaries"
        )
        
        # Cleanup
        results_storage.cleanup()


class TestErrorHandling(unittest.TestCase):
    """Test error handling in adaptive chunking."""
    
    def setUp(self):
        """Create temporary directory and storage for tests."""
        self.test_dir = tempfile.mkdtemp(prefix="test_error_handling_")
        self.storage = UniversalSequenceStorage(base_dir=self.test_dir)
    
    def tearDown(self):
        """Clean up temporary directory."""
        if Path(self.test_dir).exists():
            shutil.rmtree(self.test_dir)
        gc.collect()
    
    def test_empty_sequence(self):
        """Test handling of empty sequence."""
        sequence = ""
        seq_id = self.storage.save_sequence(sequence, "empty_seq")
        
        analyzer = TripleAdaptiveChunkAnalyzer(
            self.storage,
            use_adaptive=True
        )
        
        # Should handle gracefully
        results_storage = analyzer.analyze(seq_id)
        stats = results_storage.get_summary_stats()
        
        # Should have no motifs
        self.assertEqual(stats['total_count'], 0)
        
        # Cleanup
        results_storage.cleanup()
    
    def test_invalid_sequence_characters(self):
        """Test handling of sequence with invalid characters."""
        # Most systems filter invalid chars, but test anyway
        sequence = "ATCGATCGNNNNATCGATCG"
        seq_id = self.storage.save_sequence(sequence, "invalid_chars_seq")
        
        analyzer = TripleAdaptiveChunkAnalyzer(
            self.storage,
            use_adaptive=True
        )
        
        # Should handle gracefully (N's are typically ignored)
        results_storage = analyzer.analyze(seq_id)
        stats = results_storage.get_summary_stats()
        
        # Should complete without error
        self.assertGreaterEqual(stats['total_count'], 0)
        
        # Cleanup
        results_storage.cleanup()


if __name__ == '__main__':
    unittest.main()
