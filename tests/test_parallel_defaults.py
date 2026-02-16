"""
Test to verify parallel processing defaults are enabled.

Verifies that ChunkAnalyzer defaults to parallel=True and adaptive=False.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from Utilities.disk_storage import UniversalSequenceStorage
from Utilities.chunk_analyzer import ChunkAnalyzer


class TestParallelDefaults(unittest.TestCase):
    """Test that parallel processing is enabled by default."""
    
    def setUp(self):
        """Create temporary directory for tests."""
        self.test_dir = tempfile.mkdtemp(prefix="test_parallel_defaults_")
        self.storage = UniversalSequenceStorage(base_dir=self.test_dir)
    
    def tearDown(self):
        """Clean up temporary directory."""
        if Path(self.test_dir).exists():
            shutil.rmtree(self.test_dir)
    
    def test_parallel_default_is_true(self):
        """Test that use_parallel defaults to True."""
        analyzer = ChunkAnalyzer(self.storage)
        self.assertTrue(analyzer.use_parallel, 
                       "use_parallel should default to True")
    
    def test_adaptive_default_is_false(self):
        """Test that use_adaptive defaults to False."""
        analyzer = ChunkAnalyzer(self.storage)
        self.assertFalse(analyzer.use_adaptive, 
                       "use_adaptive should default to False")
    
    def test_parallel_can_be_disabled(self):
        """Test that parallel processing can still be explicitly disabled."""
        analyzer = ChunkAnalyzer(self.storage, use_parallel=False)
        self.assertFalse(analyzer.use_parallel,
                        "use_parallel should be False when explicitly set")
    
    def test_adaptive_can_be_enabled(self):
        """Test that adaptive chunking can still be explicitly enabled."""
        analyzer = ChunkAnalyzer(self.storage, use_adaptive=True)
        self.assertTrue(analyzer.use_adaptive,
                        "use_adaptive should be True when explicitly set")
    
    def test_max_workers_is_set(self):
        """Test that max_workers is properly configured."""
        analyzer = ChunkAnalyzer(self.storage)
        self.assertIsNotNone(analyzer.max_workers,
                           "max_workers should be set")
        self.assertGreater(analyzer.max_workers, 0,
                          "max_workers should be greater than 0")
    
    def test_custom_max_workers(self):
        """Test that custom max_workers can be specified."""
        custom_workers = 2
        analyzer = ChunkAnalyzer(self.storage, max_workers=custom_workers)
        self.assertEqual(analyzer.max_workers, custom_workers,
                        f"max_workers should be {custom_workers}")


if __name__ == '__main__':
    unittest.main()
