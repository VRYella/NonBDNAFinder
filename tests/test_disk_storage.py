"""
Unit tests for disk_storage.py - Universal disk-based storage system.

Tests UniversalSequenceStorage and UniversalResultsStorage classes.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
import json

# Add parent directory to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from Utilities.disk_storage import UniversalSequenceStorage, UniversalResultsStorage


class TestUniversalSequenceStorage(unittest.TestCase):
    """Test cases for UniversalSequenceStorage class."""
    
    def setUp(self):
        """Create temporary directory for tests."""
        self.test_dir = tempfile.mkdtemp(prefix="test_seq_storage_")
        self.storage = UniversalSequenceStorage(base_dir=self.test_dir)
    
    def tearDown(self):
        """Clean up temporary directory."""
        if Path(self.test_dir).exists():
            shutil.rmtree(self.test_dir)
    
    def test_save_sequence(self):
        """Test saving a sequence to disk."""
        sequence = "ATCGATCGATCG" * 100
        name = "test_sequence"
        
        seq_id = self.storage.save_sequence(sequence, name)
        
        # Verify seq_id is returned
        self.assertIsNotNone(seq_id)
        self.assertIsInstance(seq_id, str)
        
        # Verify metadata is stored
        metadata = self.storage.get_metadata(seq_id)
        self.assertEqual(metadata['name'], name)
        self.assertEqual(metadata['length'], len(sequence))
        self.assertGreater(metadata['gc_content'], 0)
    
    def test_get_sequence_chunk(self):
        """Test retrieving specific sequence chunks."""
        sequence = "ATCGATCGATCG" * 100  # 1200 bp
        seq_id = self.storage.save_sequence(sequence, "test_seq")
        
        # Get first 100 bp
        chunk = self.storage.get_sequence_chunk(seq_id, 0, 100)
        self.assertEqual(len(chunk), 100)
        self.assertEqual(chunk, sequence[0:100])
        
        # Get middle chunk
        chunk = self.storage.get_sequence_chunk(seq_id, 500, 600)
        self.assertEqual(len(chunk), 100)
        self.assertEqual(chunk, sequence[500:600])
    
    def test_iter_chunks(self):
        """Test chunk iteration with overlap."""
        sequence = "ATCG" * 10000  # 40,000 bp
        seq_id = self.storage.save_sequence(sequence, "test_seq")
        
        chunk_size = 15000
        overlap = 1000
        
        chunks = list(self.storage.iter_chunks(seq_id, chunk_size, overlap))
        
        # Should have 3 chunks: [0-15000], [14000-29000], [28000-40000]
        self.assertEqual(len(chunks), 3)
        
        # Verify first chunk
        chunk1, start1, end1 = chunks[0]
        self.assertEqual(start1, 0)
        self.assertEqual(end1, 15000)
        self.assertEqual(len(chunk1), 15000)
        
        # Verify overlap between chunks
        chunk2, start2, end2 = chunks[1]
        self.assertEqual(start2, 14000)  # 15000 - 1000
        
        # Verify last chunk
        chunk3, start3, end3 = chunks[2]
        self.assertEqual(end3, 40000)
    
    def test_gc_content_calculation(self):
        """Test GC content calculation."""
        # 50% GC content
        sequence = "ATCGATCGATCG" * 100
        seq_id = self.storage.save_sequence(sequence, "test_seq")
        
        metadata = self.storage.get_metadata(seq_id)
        gc_content = metadata['gc_content']
        
        # Should be 50% (ATCG has 2 GC out of 4)
        self.assertAlmostEqual(gc_content, 50.0, places=1)
    
    def test_list_sequences(self):
        """Test listing all stored sequences."""
        seq1 = self.storage.save_sequence("ATCG" * 100, "seq1")
        seq2 = self.storage.save_sequence("GCTA" * 100, "seq2")
        
        sequences = self.storage.list_sequences()
        self.assertEqual(len(sequences), 2)
        
        names = [seq['name'] for seq in sequences]
        self.assertIn("seq1", names)
        self.assertIn("seq2", names)
    
    def test_cleanup_specific_sequence(self):
        """Test deleting a specific sequence."""
        seq1 = self.storage.save_sequence("ATCG" * 100, "seq1")
        seq2 = self.storage.save_sequence("GCTA" * 100, "seq2")
        
        # Delete seq1
        self.storage.cleanup(seq1)
        
        # Verify seq1 is gone but seq2 remains
        sequences = self.storage.list_sequences()
        self.assertEqual(len(sequences), 1)
        self.assertEqual(sequences[0]['name'], "seq2")
    
    def test_cleanup_all(self):
        """Test deleting all sequences."""
        self.storage.save_sequence("ATCG" * 100, "seq1")
        self.storage.save_sequence("GCTA" * 100, "seq2")
        
        # Clean up all
        self.storage.cleanup()
        
        # Verify directory is deleted
        self.assertFalse(Path(self.test_dir).exists())


class TestUniversalResultsStorage(unittest.TestCase):
    """Test cases for UniversalResultsStorage class."""
    
    def setUp(self):
        """Create temporary directory for tests."""
        self.test_dir = tempfile.mkdtemp(prefix="test_results_storage_")
        self.results = UniversalResultsStorage(self.test_dir, "test_seq")
    
    def tearDown(self):
        """Clean up temporary directory."""
        if Path(self.test_dir).exists():
            shutil.rmtree(self.test_dir)
    
    def test_append_motif(self):
        """Test appending single motif."""
        motif = {
            'Class': 'G-Quadruplex',
            'Subclass': 'Canonical',
            'Start': 100,
            'End': 150,
            'Score': 2.5
        }
        
        self.results.append(motif)
        
        # Verify motif is stored
        motifs = list(self.results.iter_results())
        self.assertEqual(len(motifs), 1)
        self.assertEqual(motifs[0]['Class'], 'G-Quadruplex')
    
    def test_append_batch(self):
        """Test appending multiple motifs efficiently."""
        motifs = [
            {'Class': 'G-Quadruplex', 'Start': i * 100, 'End': i * 100 + 50, 'Score': 2.0}
            for i in range(100)
        ]
        
        self.results.append_batch(motifs)
        
        # Verify all motifs are stored
        stored_motifs = list(self.results.iter_results())
        self.assertEqual(len(stored_motifs), 100)
    
    def test_iter_results_with_limit(self):
        """Test iterating with limit."""
        motifs = [
            {'Class': 'G-Quadruplex', 'Start': i * 100, 'End': i * 100 + 50}
            for i in range(100)
        ]
        self.results.append_batch(motifs)
        
        # Get only first 10
        limited_motifs = list(self.results.iter_results(limit=10))
        self.assertEqual(len(limited_motifs), 10)
    
    def test_get_summary_stats(self):
        """Test summary statistics calculation."""
        motifs = [
            {'Class': 'G-Quadruplex', 'Subclass': 'Canonical', 
             'Start': 100, 'End': 150, 'Length': 50, 'Score': 2.5},
            {'Class': 'G-Quadruplex', 'Subclass': 'Canonical',
             'Start': 200, 'End': 250, 'Length': 50, 'Score': 2.0},
            {'Class': 'Z-DNA', 'Subclass': 'Standard',
             'Start': 300, 'End': 320, 'Length': 20, 'Score': 1.5},
        ]
        self.results.append_batch(motifs)
        
        stats = self.results.get_summary_stats()
        
        # Verify statistics
        self.assertEqual(stats['total_count'], 3)
        self.assertEqual(stats['class_distribution']['G-Quadruplex'], 2)
        self.assertEqual(stats['class_distribution']['Z-DNA'], 1)
        self.assertEqual(stats['coverage_bp'], 120)  # 50 + 50 + 20
        self.assertAlmostEqual(stats['avg_score'], 2.0, places=1)
    
    def test_to_dataframe(self):
        """Test conversion to pandas DataFrame."""
        motifs = [
            {'Class': 'G-Quadruplex', 'Start': i * 100, 'End': i * 100 + 50}
            for i in range(10)
        ]
        self.results.append_batch(motifs)
        
        df = self.results.to_dataframe()
        
        # Verify DataFrame
        self.assertEqual(len(df), 10)
        self.assertIn('Class', df.columns)
        self.assertIn('Start', df.columns)
    
    def test_to_dataframe_with_limit(self):
        """Test DataFrame conversion with limit."""
        motifs = [
            {'Class': 'G-Quadruplex', 'Start': i * 100, 'End': i * 100 + 50}
            for i in range(100)
        ]
        self.results.append_batch(motifs)
        
        df = self.results.to_dataframe(limit=20)
        
        # Verify only 20 rows
        self.assertEqual(len(df), 20)
    
    def test_cleanup(self):
        """Test cleanup of results files."""
        self.results.append({'Class': 'G-Quadruplex', 'Start': 100, 'End': 150})
        
        # Verify files exist
        self.assertTrue(self.results.results_file.exists())
        
        # Clean up
        self.results.cleanup()
        
        # Verify files are deleted
        self.assertFalse(self.results.results_file.exists())


if __name__ == '__main__':
    unittest.main()
