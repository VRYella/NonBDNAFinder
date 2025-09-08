"""
Unit tests for Hyperscan database compilation and pattern management.
"""

import pytest
import hs_registry_manager
from motifs.registry import get_all_hyperscan_patterns, get_hyperscan_safe_patterns


class TestHyperscanCompilation:
    """Test Hyperscan database compilation"""
    
    def test_hyperscan_available(self):
        """Test that Hyperscan is available"""
        assert hs_registry_manager.is_available()
    
    def test_pattern_loading(self):
        """Test pattern loading from registry"""
        patterns = get_all_hyperscan_patterns()
        assert len(patterns) > 0
        
        # Check pattern format
        for pattern, pat_id, class_name, subclass in patterns:
            assert isinstance(pattern, str)
            assert isinstance(pat_id, int)
            assert isinstance(class_name, str)
            assert isinstance(subclass, str)
    
    def test_safe_patterns_only(self):
        """Test that only Hyperscan-safe patterns are included"""
        safe_patterns = get_hyperscan_safe_patterns()
        
        for class_name, patterns in safe_patterns.items():
            for pattern in patterns:
                # Check for backreferences (should not be present)
                assert '\\1' not in pattern
                assert '\\2' not in pattern
                assert '(?P=' not in pattern
    
    def test_database_compilation(self):
        """Test Hyperscan database compilation"""
        test_patterns = [
            'GGGTTTGGGTTTGGGTTTGGG',
            'AAAAAAA',
            'CCCAAACCCAAACCCAAACCC'
        ]
        
        # Compile database
        db_info = hs_registry_manager.compile_db(test_patterns, key="test_db")
        
        assert db_info is not None
        assert 'db' in db_info
        assert 'patterns' in db_info
        assert db_info['count'] == len(test_patterns)
    
    def test_database_caching(self):
        """Test that databases are properly cached"""
        test_patterns = ['GGGG', 'AAAA', 'TTTT']
        
        # First compilation
        db_info1 = hs_registry_manager.compile_db(test_patterns, key="cache_test")
        
        # Second compilation (should use cache)
        db_info2 = hs_registry_manager.compile_db(test_patterns, key="cache_test")
        
        # Should be the same object (cached)
        assert db_info1 is db_info2
        
        # Check cache keys
        cached_keys = hs_registry_manager.get_cached_keys()
        assert "cache_test" in cached_keys
    
    def test_scanning_basic(self):
        """Test basic sequence scanning"""
        test_patterns = ['GGG', 'AAA', 'TTT']
        test_sequence = "AAAGGGAAATTTGGGAAA"
        
        # Compile and scan
        hs_registry_manager.compile_db(test_patterns, key="scan_test")
        matches = hs_registry_manager.scan_block(test_sequence, key="scan_test")
        
        assert len(matches) > 0
        
        # Check match format
        for pat_id, start, end in matches:
            assert isinstance(pat_id, int)
            assert isinstance(start, int)
            assert isinstance(end, int)
            assert start < end
    
    def test_cache_clearing(self):
        """Test cache clearing functionality"""
        test_patterns = ['ATCG']
        
        # Compile database
        hs_registry_manager.compile_db(test_patterns, key="clear_test")
        assert "clear_test" in hs_registry_manager.get_cached_keys()
        
        # Clear specific cache
        hs_registry_manager.clear_cache("clear_test")
        assert "clear_test" not in hs_registry_manager.get_cached_keys()
    
    def test_invalid_patterns(self):
        """Test handling of invalid patterns"""
        with pytest.raises(ValueError):
            hs_registry_manager.compile_db([], key="empty_test")
    
    def test_database_info(self):
        """Test database information retrieval"""
        test_patterns = ['ATCG', 'GCTA']
        hs_registry_manager.compile_db(test_patterns, key="info_test")
        
        info = hs_registry_manager.get_db_info("info_test")
        assert info is not None
        assert info['pattern_count'] == 2
        assert info['key'] == "info_test"
        
        # Test non-existent database
        info = hs_registry_manager.get_db_info("nonexistent")
        assert info is None