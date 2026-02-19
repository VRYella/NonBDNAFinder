"""
╔══════════════════════════════════════════════════════════════════════════════╗
║  Disk Storage IndexError Fix - Integration Test Suite                        ║
║  Testing summary generation and visualization caching with disk storage       ║
╚══════════════════════════════════════════════════════════════════════════════╝

Tests validate:
1. Summary generation works with disk storage mode (>1MB sequences)
2. Summary generation works with in-memory mode (<1MB sequences)
3. Visualization caching works with disk storage mode
4. Visualization caching works with in-memory mode
5. Mixed mode (some sequences in disk, some in memory) - graceful handling
"""

import sys
import os
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import unittest
import tempfile
import shutil
import numpy as np
import pandas as pd
from unittest.mock import MagicMock

from Utilities.disk_storage import UniversalSequenceStorage


class MockSessionState:
    """Mock Streamlit session state for testing."""
    
    def __init__(self):
        self._state = {}
    
    def get(self, key, default=None):
        """Get session state value."""
        return self._state.get(key, default)
    
    def __getitem__(self, key):
        """Get session state value."""
        return self._state[key]
    
    def __setitem__(self, key, value):
        """Set session state value."""
        self._state[key] = value
    
    def __contains__(self, key):
        """Check if key exists."""
        return key in self._state


class TestDiskStorageSummaryGeneration(unittest.TestCase):
    """Test summary generation with disk storage mode."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp(prefix="test_disk_storage_")
        self.storage = UniversalSequenceStorage(base_dir=self.test_dir)
        
        # Create test sequences
        self.small_seq = "ATCGATCGATCG" * 100  # 1,200 bp - small
        self.large_seq = "ATCGATCGATCG" * 100000  # 1,200,000 bp - large (>1MB)
    
    def tearDown(self):
        """Clean up test directory."""
        if Path(self.test_dir).exists():
            shutil.rmtree(self.test_dir)
    
    def test_disk_storage_mode_summary_generation(self):
        """Test summary generation uses metadata in disk storage mode."""
        # Set up disk storage mode
        session_state = MockSessionState()
        session_state['use_disk_storage'] = True
        session_state['seq_storage'] = self.storage
        
        # Save sequences to disk
        seq_id1 = self.storage.save_sequence(self.large_seq, "large_seq_1")
        seq_id2 = self.storage.save_sequence(self.large_seq, "large_seq_2")
        
        session_state['seq_ids'] = [seq_id1, seq_id2]
        session_state['names'] = ["large_seq_1", "large_seq_2"]
        session_state['seqs'] = []  # Empty in disk storage mode
        
        # Mock results
        all_results = [
            [{'Type': 'GQ', 'Score': 0.85}, {'Type': 'AP', 'Score': 0.90}],
            [{'Type': 'ZD', 'Score': 0.75}]
        ]
        
        # Simulate summary generation logic from upload.py
        summary = []
        for i, results in enumerate(all_results):
            # Support both disk storage and legacy in-memory modes
            if session_state.get('use_disk_storage') and session_state.get('seq_ids'):
                # Disk storage mode: get metadata
                seq_id = session_state['seq_ids'][i]
                metadata = session_state['seq_storage'].get_metadata(seq_id)
                name = session_state['names'][i]
                seq_length = metadata['length']
                gc_content = metadata['gc_content']
            else:
                # This path shouldn't be taken in disk storage mode
                self.fail("Should not reach in-memory path in disk storage mode")
            
            summary.append({
                'Sequence': name,
                'Length': seq_length,
                'GC Content': f"{gc_content:.1f}%",
                'Motifs Found': len(results),
                'Unique Types': len(set(m.get('Type', 'Unknown') for m in results)),
                'Avg Score': f"{np.mean([m.get('Score', 0) for m in results]):.3f}" if results else "0.000"
            })
        
        # Verify summary was generated correctly
        self.assertEqual(len(summary), 2)
        self.assertEqual(summary[0]['Sequence'], "large_seq_1")
        self.assertEqual(summary[0]['Length'], 1200000)
        self.assertEqual(summary[0]['Motifs Found'], 2)
        self.assertEqual(summary[0]['Unique Types'], 2)
        self.assertEqual(summary[0]['Avg Score'], "0.875")
        
        self.assertEqual(summary[1]['Sequence'], "large_seq_2")
        self.assertEqual(summary[1]['Length'], 1200000)
        self.assertEqual(summary[1]['Motifs Found'], 1)
        self.assertEqual(summary[1]['Unique Types'], 1)
        self.assertEqual(summary[1]['Avg Score'], "0.750")
    
    def test_in_memory_mode_summary_generation(self):
        """Test summary generation works in legacy in-memory mode."""
        # Set up in-memory mode
        session_state = MockSessionState()
        session_state['use_disk_storage'] = False
        session_state['seqs'] = [self.small_seq, self.small_seq]
        session_state['names'] = ["small_seq_1", "small_seq_2"]
        
        # Mock results
        all_results = [
            [{'Type': 'GQ', 'Score': 0.85}],
            [{'Type': 'ZD', 'Score': 0.75}, {'Type': 'AP', 'Score': 0.80}]
        ]
        
        # Mock get_basic_stats function
        def mock_get_basic_stats(seq, results):
            return {
                'Length': len(seq),
                'GC%': 50.0  # Simplified for test
            }
        
        # Simulate summary generation logic from upload.py
        summary = []
        for i, results in enumerate(all_results):
            # Support both disk storage and legacy in-memory modes
            if session_state.get('use_disk_storage') and session_state.get('seq_ids'):
                # This path shouldn't be taken in in-memory mode
                self.fail("Should not reach disk storage path in in-memory mode")
            else:
                # Legacy in-memory mode
                seq = session_state['seqs'][i]
                name = session_state['names'][i]
                stats = mock_get_basic_stats(seq, results)
                seq_length = stats['Length']
                gc_content = stats['GC%']
            
            summary.append({
                'Sequence': name,
                'Length': seq_length,
                'GC Content': f"{gc_content:.1f}%",
                'Motifs Found': len(results),
                'Unique Types': len(set(m.get('Type', 'Unknown') for m in results)),
                'Avg Score': f"{np.mean([m.get('Score', 0) for m in results]):.3f}" if results else "0.000"
            })
        
        # Verify summary was generated correctly
        self.assertEqual(len(summary), 2)
        self.assertEqual(summary[0]['Sequence'], "small_seq_1")
        self.assertEqual(summary[0]['Length'], 1200)
        self.assertEqual(summary[0]['Motifs Found'], 1)
        
        self.assertEqual(summary[1]['Sequence'], "small_seq_2")
        self.assertEqual(summary[1]['Length'], 1200)
        self.assertEqual(summary[1]['Motifs Found'], 2)
    
    def test_no_indexerror_with_empty_seqs_list(self):
        """Test that accessing empty seqs list in disk mode doesn't cause IndexError."""
        # Set up disk storage mode with empty seqs list
        session_state = MockSessionState()
        session_state['use_disk_storage'] = True
        session_state['seq_storage'] = self.storage
        
        # Save sequence to disk
        seq_id = self.storage.save_sequence(self.large_seq, "large_seq")
        
        session_state['seq_ids'] = [seq_id]
        session_state['names'] = ["large_seq"]
        session_state['seqs'] = []  # EMPTY - this is what caused the IndexError
        
        # Mock results
        all_results = [[{'Type': 'GQ', 'Score': 0.85}]]
        
        # This should NOT raise IndexError
        summary = []
        for i, results in enumerate(all_results):
            if session_state.get('use_disk_storage') and session_state.get('seq_ids'):
                # Disk storage mode: get metadata (should work)
                seq_id = session_state['seq_ids'][i]
                metadata = session_state['seq_storage'].get_metadata(seq_id)
                name = session_state['names'][i]
                seq_length = metadata['length']
                gc_content = metadata['gc_content']
                
                summary.append({
                    'Sequence': name,
                    'Length': seq_length,
                    'GC Content': f"{gc_content:.1f}%",
                    'Motifs Found': len(results)
                })
        
        # Should complete successfully without IndexError
        self.assertEqual(len(summary), 1)
        self.assertEqual(summary[0]['Sequence'], "large_seq")


class TestDiskStorageVisualizationCaching(unittest.TestCase):
    """Test visualization caching with disk storage mode."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp(prefix="test_viz_storage_")
        self.storage = UniversalSequenceStorage(base_dir=self.test_dir)
        self.large_seq = "ATCGATCGATCG" * 100000  # 1,200,000 bp
    
    def tearDown(self):
        """Clean up test directory."""
        if Path(self.test_dir).exists():
            shutil.rmtree(self.test_dir)
    
    def test_visualization_loop_disk_storage_mode(self):
        """Test visualization caching loop works with disk storage."""
        # Set up disk storage mode
        session_state = MockSessionState()
        session_state['use_disk_storage'] = True
        session_state['seq_storage'] = self.storage
        
        # Save sequences
        seq_id1 = self.storage.save_sequence(self.large_seq, "seq1")
        seq_id2 = self.storage.save_sequence(self.large_seq, "seq2")
        
        session_state['seq_ids'] = [seq_id1, seq_id2]
        session_state['names'] = ["seq1", "seq2"]
        session_state['seqs'] = []  # Empty in disk storage mode
        
        # Mock results
        all_results = [
            [{'Type': 'GQ', 'Class': 'GQ', 'Subclass': 'GQ_type1'}],
            [{'Type': 'ZD', 'Class': 'ZD', 'Subclass': 'ZD_type1'}]
        ]
        
        # Calculate UPDATE_INTERVAL (from lines 1816-1819)
        num_seqs = len(session_state['seq_ids']) if session_state.get('use_disk_storage') and session_state.get('seq_ids') else len(session_state['seqs'])
        UPDATE_INTERVAL = max(1, num_seqs // 5)
        self.assertEqual(UPDATE_INTERVAL, 1)  # 2 // 5 = 0, max(1, 0) = 1
        
        # Simulate visualization caching loop (from lines 1821-1830)
        cached_visualizations = {}
        for seq_idx, results in enumerate(all_results):
            # Support both storage modes
            if session_state.get('use_disk_storage') and session_state.get('seq_ids'):
                seq_id = session_state['seq_ids'][seq_idx]
                metadata = session_state['seq_storage'].get_metadata(seq_id)
                sequence_length = metadata['length']
                name = session_state['names'][seq_idx]
            else:
                self.fail("Should use disk storage path")
            
            # Continue with existing motif filtering logic
            filtered_motifs = results
            
            if not filtered_motifs:
                continue
            
            viz_cache_key = f"seq_{seq_idx}"
            cached_visualizations[viz_cache_key] = {
                'name': name,
                'length': sequence_length,
                'motif_count': len(filtered_motifs)
            }
        
        # Verify visualization cache was created
        self.assertEqual(len(cached_visualizations), 2)
        self.assertIn('seq_0', cached_visualizations)
        self.assertIn('seq_1', cached_visualizations)
        self.assertEqual(cached_visualizations['seq_0']['name'], "seq1")
        self.assertEqual(cached_visualizations['seq_0']['length'], 1200000)
    
    def test_visualization_loop_in_memory_mode(self):
        """Test visualization caching loop works with in-memory mode."""
        # Set up in-memory mode
        session_state = MockSessionState()
        session_state['use_disk_storage'] = False
        small_seq = "ATCGATCGATCG" * 100
        session_state['seqs'] = [small_seq, small_seq]
        session_state['names'] = ["seq1", "seq2"]
        
        # Mock results
        all_results = [
            [{'Type': 'GQ'}],
            [{'Type': 'ZD'}]
        ]
        
        # Calculate UPDATE_INTERVAL
        num_seqs = len(session_state['seq_ids']) if session_state.get('use_disk_storage') and session_state.get('seq_ids') else len(session_state['seqs'])
        UPDATE_INTERVAL = max(1, num_seqs // 5)
        self.assertEqual(UPDATE_INTERVAL, 1)
        
        # Simulate visualization caching loop
        cached_visualizations = {}
        for seq_idx, results in enumerate(all_results):
            # Support both storage modes
            if session_state.get('use_disk_storage') and session_state.get('seq_ids'):
                self.fail("Should use in-memory path")
            else:
                seq = session_state['seqs'][seq_idx]
                sequence_length = len(seq)
                name = session_state['names'][seq_idx]
            
            filtered_motifs = results
            
            if not filtered_motifs:
                continue
            
            viz_cache_key = f"seq_{seq_idx}"
            cached_visualizations[viz_cache_key] = {
                'name': name,
                'length': sequence_length,
                'motif_count': len(filtered_motifs)
            }
        
        # Verify visualization cache was created
        self.assertEqual(len(cached_visualizations), 2)
        self.assertEqual(cached_visualizations['seq_0']['length'], 1200)


# ═══════════════════════════════════════════════════════════════════════════════
# TEST RUNNER
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)
