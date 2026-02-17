"""
Test chunk size consistency across the codebase.

Ensures that all chunking operations use 50,000 bp as specified in requirements.
"""

import unittest
import sys
import os
import re

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class TestChunkSizeConsistency(unittest.TestCase):
    """Test that chunk size is consistently set to 50,000 bp everywhere."""
    
    def test_analysis_config_chunk_sizes(self):
        """Test that chunk sizes in analysis config are 50,000."""
        config_file = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'Utilities', 'config', 'analysis.py'
        )
        
        with open(config_file, 'r') as f:
            content = f.read()
        
        # Check default_chunk_size (handles both 50000 and 50_000)
        match = re.search(r"'default_chunk_size':\s*([\d_]+)", content)
        self.assertIsNotNone(match, "default_chunk_size not found in analysis.py")
        value = int(match.group(1).replace('_', ''))
        self.assertEqual(value, 50_000, 
                        "default_chunk_size should be 50,000")
        
        # Check micro_chunk_size
        match = re.search(r"'micro_chunk_size':\s*([\d_]+)", content)
        self.assertIsNotNone(match, "micro_chunk_size not found in analysis.py")
        value = int(match.group(1).replace('_', ''))
        self.assertEqual(value, 50_000,
                        "micro_chunk_size should be 50,000")
        
        # Check meso_chunk_size
        match = re.search(r"'meso_chunk_size':\s*([\d_]+)", content)
        self.assertIsNotNone(match, "meso_chunk_size not found in analysis.py")
        value = int(match.group(1).replace('_', ''))
        self.assertEqual(value, 50_000,
                        "meso_chunk_size should be 50,000")
        
        # Check macro_chunk_size
        match = re.search(r"'macro_chunk_size':\s*([\d_]+)", content)
        self.assertIsNotNone(match, "macro_chunk_size not found in analysis.py")
        value = int(match.group(1).replace('_', ''))
        self.assertEqual(value, 50_000,
                        "macro_chunk_size should be 50,000")
    
    def test_nonbscanner_default_chunk_size(self):
        """Test that nonbscanner DEFAULT_CHUNK_SIZE is 50,000."""
        scanner_file = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'Utilities', 'nonbscanner.py'
        )
        
        with open(scanner_file, 'r') as f:
            content = f.read()
        
        # Find DEFAULT_CHUNK_SIZE
        match = re.search(r'DEFAULT_CHUNK_SIZE\s*=\s*(\d+)', content)
        self.assertIsNotNone(match, "DEFAULT_CHUNK_SIZE not found in nonbscanner.py")
        self.assertEqual(int(match.group(1)), 50_000,
                        "DEFAULT_CHUNK_SIZE should be 50,000")
    
    def test_chunk_analyzer_default(self):
        """Test that ChunkAnalyzer default is 50,000."""
        analyzer_file = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'Utilities', 'chunk_analyzer.py'
        )
        
        with open(analyzer_file, 'r') as f:
            content = f.read()
        
        # Find chunk_size parameter default (handles both 50000 and 50_000)
        match = re.search(r'chunk_size:\s*int\s*=\s*([\d_]+)', content)
        self.assertIsNotNone(match, "chunk_size default not found in chunk_analyzer.py")
        value = int(match.group(1).replace('_', ''))
        self.assertEqual(value, 50_000,
                        "ChunkAnalyzer default chunk_size should be 50,000")
    
    def test_disk_storage_default(self):
        """Test that disk_storage default is 50,000."""
        storage_file = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'Utilities', 'disk_storage.py'
        )
        
        with open(storage_file, 'r') as f:
            content = f.read()
        
        # Find chunk_size parameter default (handles both 50000 and 50_000)
        match = re.search(r'chunk_size:\s*int\s*=\s*([\d_]+)', content)
        self.assertIsNotNone(match, "chunk_size default not found in disk_storage.py")
        value = int(match.group(1).replace('_', ''))
        self.assertEqual(value, 50_000,
                        "disk_storage default chunk_size should be 50,000")


if __name__ == '__main__':
    unittest.main()
