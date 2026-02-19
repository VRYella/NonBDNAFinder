"""
Unit tests for disk_storage.py FASTA sanitization feature.

Tests that UniversalSequenceStorage correctly sanitizes sequences with
FASTA formatting (newlines, whitespace) to enable correct byte-offset based chunking.

This addresses the issue where >1MB sequences with FASTA formatting would
return 0 motifs due to incorrect chunk reading.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from Utilities.disk_storage import UniversalSequenceStorage


class TestFastaSanitization(unittest.TestCase):
    """Test cases for FASTA sanitization in UniversalSequenceStorage."""
    
    def setUp(self):
        """Create temporary directory for tests."""
        self.test_dir = tempfile.mkdtemp(prefix="test_fasta_sanitization_")
        self.storage = UniversalSequenceStorage(base_dir=self.test_dir)
    
    def tearDown(self):
        """Clean up temporary directory."""
        if Path(self.test_dir).exists():
            shutil.rmtree(self.test_dir)
    
    def test_save_sequence_strips_newlines(self):
        """Test that save_sequence removes newlines from FASTA-formatted sequences."""
        # Create a sequence with FASTA newlines (80 chars per line)
        base_seq = "ATCGATCGATCGATCG" * 100  # 1600 bp
        fasta_seq = ""
        line_length = 80
        for i in range(0, len(base_seq), line_length):
            fasta_seq += base_seq[i:i+line_length] + "\n"
        
        # Input has newlines
        self.assertIn('\n', fasta_seq)
        self.assertEqual(len(fasta_seq), len(base_seq) + len(base_seq) // line_length)
        
        # Save sequence
        seq_id = self.storage.save_sequence(fasta_seq, "test_fasta")
        
        # Verify metadata reflects sanitized length (without newlines)
        metadata = self.storage.get_metadata(seq_id)
        self.assertEqual(metadata['length'], len(base_seq))
        
        # Verify stored sequence has no newlines
        stored_chunk = self.storage.get_sequence_chunk(seq_id, 0, len(base_seq))
        self.assertNotIn('\n', stored_chunk)
        self.assertEqual(len(stored_chunk), len(base_seq))
        self.assertEqual(stored_chunk, base_seq.upper())  # Should also be uppercase
    
    def test_save_sequence_strips_all_whitespace(self):
        """Test that save_sequence removes all types of whitespace."""
        # Create sequence with various whitespace
        base_seq = "ATCGATCG"
        whitespace_seq = "ATCG\nATCG\t\r \n"  # newlines, tabs, carriage returns, spaces
        
        seq_id = self.storage.save_sequence(whitespace_seq, "test_whitespace")
        
        # Verify all whitespace is removed
        metadata = self.storage.get_metadata(seq_id)
        self.assertEqual(metadata['length'], len(base_seq))
        
        stored_seq = self.storage.get_sequence_chunk(seq_id, 0, len(base_seq))
        self.assertEqual(stored_seq, base_seq.upper())
        self.assertFalse(any(c.isspace() for c in stored_seq))
    
    def test_save_sequence_converts_to_uppercase(self):
        """Test that save_sequence converts sequences to uppercase."""
        seq = "atcg"
        seq_id = self.storage.save_sequence(seq, "test_case")
        
        stored_seq = self.storage.get_sequence_chunk(seq_id, 0, 4)
        self.assertEqual(stored_seq, "ATCG")
    
    def test_iter_chunks_with_sanitized_sequence(self):
        """Test that iter_chunks works correctly with sanitized sequences."""
        # Create a sequence with FASTA formatting
        base_seq = "ATCG" * 5000  # 20KB
        fasta_seq = ""
        for i in range(0, len(base_seq), 80):
            fasta_seq += base_seq[i:i+80] + "\n"
        
        seq_id = self.storage.save_sequence(fasta_seq, "test_chunking")
        
        # Iterate chunks - should work without errors
        chunk_count = 0
        total_bases_seen = 0
        for chunk_seq, start, end in self.storage.iter_chunks(seq_id, chunk_size=8000, overlap=500):
            chunk_count += 1
            
            # Verify no whitespace in chunks
            self.assertFalse(any(c.isspace() for c in chunk_seq))
            
            # Verify chunk length matches expected
            expected_length = end - start
            self.assertEqual(len(chunk_seq), expected_length)
            
            # Verify chunk matches original sequence at correct positions
            self.assertEqual(chunk_seq, base_seq[start:end].upper())
            
            total_bases_seen += len(chunk_seq)
        
        # Verify we processed the expected number of chunks
        self.assertGreater(chunk_count, 1)  # Should have multiple chunks
        self.assertGreaterEqual(total_bases_seen, len(base_seq))  # Due to overlaps
    
    def test_iter_chunks_validates_clean_storage(self):
        """Test that iter_chunks validates chunks don't contain whitespace."""
        # Manually create a sequence file with newlines (bypassing save_sequence)
        # This simulates old data that wasn't sanitized
        seq_id = "manually_created"
        seq_file = self.storage.base_dir / f"{seq_id}.seq"
        
        dirty_seq = "ATCG\nATCG\nATCG\n"
        with open(seq_file, 'w') as f:
            f.write(dirty_seq)
        
        # Manually add metadata
        self.storage.metadata[seq_id] = {
            'seq_id': seq_id,
            'name': 'dirty_seq',
            'length': len(dirty_seq),  # Including newlines (wrong!)
            'gc_content': 50.0,
            'file_path': str(seq_file),
            'created_at': '2024-01-01T00:00:00'
        }
        
        # iter_chunks should raise ValueError when it detects whitespace
        with self.assertRaises(ValueError) as context:
            list(self.storage.iter_chunks(seq_id, chunk_size=10, overlap=2))
        
        # Verify error message is helpful
        error_msg = str(context.exception)
        self.assertIn("whitespace", error_msg.lower())
        self.assertIn("newline", error_msg.lower())
        self.assertIn("re-upload", error_msg.lower())
    
    def test_large_sequence_sanitization(self):
        """Test sanitization with a sequence > 1MB."""
        # Create a 1.5MB sequence with FASTA formatting
        base_seq = "ATCGATCGATCGATCG" * 100_000  # 1.6MB
        
        # Add FASTA newlines every 80 characters
        fasta_seq = ""
        line_length = 80
        for i in range(0, len(base_seq), line_length):
            fasta_seq += base_seq[i:i+line_length] + "\n"
        
        # Input size should be larger due to newlines
        self.assertGreater(len(fasta_seq), len(base_seq))
        
        # Save and verify
        seq_id = self.storage.save_sequence(fasta_seq, "large_fasta")
        
        metadata = self.storage.get_metadata(seq_id)
        self.assertEqual(metadata['length'], len(base_seq))
        
        # Verify chunks are clean
        chunk_count = 0
        for chunk_seq, start, end in self.storage.iter_chunks(seq_id, chunk_size=50_000, overlap=2_000):
            chunk_count += 1
            self.assertFalse(any(c.isspace() for c in chunk_seq))
            
            # Spot check: verify chunk matches original at correct position
            if chunk_count == 1:  # First chunk
                self.assertEqual(chunk_seq, base_seq[start:end].upper())
        
        self.assertGreater(chunk_count, 10)  # Should have many chunks for 1.6MB
    
    def test_mixed_case_and_whitespace(self):
        """Test sequence with mixed case and various whitespace."""
        input_seq = "AtCg\n GaTc \t\r\nCgAt"
        expected_seq = "ATCGGATCCGAT"
        
        seq_id = self.storage.save_sequence(input_seq, "mixed")
        
        metadata = self.storage.get_metadata(seq_id)
        self.assertEqual(metadata['length'], len(expected_seq))
        
        stored = self.storage.get_sequence_chunk(seq_id, 0, len(expected_seq))
        self.assertEqual(stored, expected_seq)
    
    def test_empty_lines_in_fasta(self):
        """Test FASTA with empty lines (common in some FASTA files)."""
        input_seq = "ATCG\n\nGATC\n\n\nCGAT"
        expected_seq = "ATCGGATCCGAT"
        
        seq_id = self.storage.save_sequence(input_seq, "empty_lines")
        
        stored = self.storage.get_sequence_chunk(seq_id, 0, len(expected_seq))
        self.assertEqual(stored, expected_seq)


if __name__ == '__main__':
    unittest.main()
