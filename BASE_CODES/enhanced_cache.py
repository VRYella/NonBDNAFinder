#!/usr/bin/env python3
"""
Enhanced Caching System for NBDFinder
=====================================

Advanced caching system for improved performance with repeated genome queries.
Implements multi-level caching for sequences, analysis results, and conservation data.

Features:
- LRU cache for analysis results with configurable size limits
- Sequence-based caching with hash-based lookup
- Conservation analysis caching
- Cache statistics and management
- Automatic cache cleanup and memory management

Author: Enhanced NBDFinder by Dr. Venkata Rajesh Yella
"""

import hashlib
import json
import time
import pickle
import os
import threading
from typing import Dict, Any, Optional, List, Tuple
from collections import OrderedDict
from dataclasses import dataclass
from datetime import datetime, timedelta
import tempfile
import shutil

@dataclass
class CacheEntry:
    """Cache entry with metadata"""
    data: Any
    created_at: datetime
    access_count: int
    last_accessed: datetime
    size_bytes: int
    cache_type: str

class LRUCache:
    """
    Thread-safe LRU Cache with size limits and automatic cleanup.
    """
    
    def __init__(self, max_size: int = 100, max_memory_mb: int = 512):
        self.max_size = max_size
        self.max_memory_bytes = max_memory_mb * 1024 * 1024
        self.cache = OrderedDict()
        self.lock = threading.RLock()
        self.stats = {
            'hits': 0,
            'misses': 0,
            'evictions': 0,
            'total_requests': 0,
            'memory_used': 0
        }
    
    def _calculate_size(self, obj: Any) -> int:
        """Estimate object size in bytes"""
        try:
            return len(pickle.dumps(obj))
        except:
            # Fallback estimation
            if isinstance(obj, str):
                return len(obj.encode('utf-8'))
            elif isinstance(obj, (list, tuple)):
                return sum(self._calculate_size(item) for item in obj[:10])  # Sample first 10
            elif isinstance(obj, dict):
                return sum(self._calculate_size(k) + self._calculate_size(v) 
                          for k, v in list(obj.items())[:10])  # Sample first 10
            else:
                return 1024  # Default estimate
    
    def _cleanup_if_needed(self):
        """Remove old entries if cache is too large"""
        # Remove entries until we're under limits
        while (len(self.cache) > self.max_size or 
               self.stats['memory_used'] > self.max_memory_bytes):
            if not self.cache:
                break
            
            # Remove least recently used
            oldest_key = next(iter(self.cache))
            removed_entry = self.cache.pop(oldest_key)
            self.stats['memory_used'] -= removed_entry.size_bytes
            self.stats['evictions'] += 1
    
    def get(self, key: str) -> Optional[Any]:
        """Get value from cache"""
        with self.lock:
            self.stats['total_requests'] += 1
            
            if key in self.cache:
                # Move to end (most recently used)
                entry = self.cache.pop(key)
                entry.access_count += 1
                entry.last_accessed = datetime.now()
                self.cache[key] = entry
                
                self.stats['hits'] += 1
                return entry.data
            else:
                self.stats['misses'] += 1
                return None
    
    def put(self, key: str, value: Any, cache_type: str = "analysis"):
        """Put value in cache"""
        with self.lock:
            size_bytes = self._calculate_size(value)
            now = datetime.now()
            
            # Remove if already exists
            if key in self.cache:
                old_entry = self.cache.pop(key)
                self.stats['memory_used'] -= old_entry.size_bytes
            
            # Create new entry
            entry = CacheEntry(
                data=value,
                created_at=now,
                access_count=1,
                last_accessed=now,
                size_bytes=size_bytes,
                cache_type=cache_type
            )
            
            self.cache[key] = entry
            self.stats['memory_used'] += size_bytes
            
            # Cleanup if needed
            self._cleanup_if_needed()
    
    def remove(self, key: str) -> bool:
        """Remove specific key from cache"""
        with self.lock:
            if key in self.cache:
                entry = self.cache.pop(key)
                self.stats['memory_used'] -= entry.size_bytes
                return True
            return False
    
    def clear(self):
        """Clear all cache entries"""
        with self.lock:
            self.cache.clear()
            self.stats = {
                'hits': 0,
                'misses': 0,
                'evictions': 0,
                'total_requests': 0,
                'memory_used': 0
            }
    
    def get_stats(self) -> Dict[str, Any]:
        """Get cache statistics"""
        with self.lock:
            hit_rate = (self.stats['hits'] / max(self.stats['total_requests'], 1)) * 100
            
            return {
                'cache_size': len(self.cache),
                'max_size': self.max_size,
                'memory_used_mb': self.stats['memory_used'] / (1024 * 1024),
                'max_memory_mb': self.max_memory_bytes / (1024 * 1024),
                'hit_rate_percent': round(hit_rate, 2),
                'total_requests': self.stats['total_requests'],
                'hits': self.stats['hits'],
                'misses': self.stats['misses'],
                'evictions': self.stats['evictions']
            }

class EnhancedCacheManager:
    """
    Comprehensive cache manager for NBDFinder with multiple cache types.
    """
    
    def __init__(self, 
                 analysis_cache_size: int = 100,
                 sequence_cache_size: int = 50,
                 conservation_cache_size: int = 30):
        
        # Multiple cache instances for different data types
        self.analysis_cache = LRUCache(analysis_cache_size, 256)  # Analysis results
        self.sequence_cache = LRUCache(sequence_cache_size, 128)   # Processed sequences
        self.conservation_cache = LRUCache(conservation_cache_size, 128)  # Conservation data
        
        # Persistent cache directory
        self.cache_dir = os.path.join(tempfile.gettempdir(), "nbdfinder_cache")
        os.makedirs(self.cache_dir, exist_ok=True)
        
        # Cache configuration
        self.enable_persistent = True
        self.persistent_ttl_hours = 24
    
    def _get_sequence_hash(self, sequence: str, params: Dict[str, Any] = None) -> str:
        """Generate hash for sequence and analysis parameters"""
        # Normalize sequence
        normalized_seq = sequence.upper().replace(' ', '').replace('\n', '')
        
        # Include analysis parameters in hash
        if params:
            param_str = json.dumps(params, sort_keys=True)
            hash_input = f"{normalized_seq}:{param_str}"
        else:
            hash_input = normalized_seq
        
        return hashlib.sha256(hash_input.encode()).hexdigest()[:16]
    
    def get_analysis_result(self, sequence: str, params: Dict[str, Any] = None) -> Optional[List[Dict[str, Any]]]:
        """Get cached analysis result"""
        cache_key = self._get_sequence_hash(sequence, params)
        
        # Try memory cache first
        result = self.analysis_cache.get(cache_key)
        if result is not None:
            return result
        
        # Try persistent cache
        if self.enable_persistent:
            result = self._load_persistent_cache(cache_key, "analysis")
            if result is not None:
                # Load back into memory cache
                self.analysis_cache.put(cache_key, result, "analysis")
                return result
        
        return None
    
    def store_analysis_result(self, sequence: str, result: List[Dict[str, Any]], params: Dict[str, Any] = None):
        """Store analysis result in cache"""
        cache_key = self._get_sequence_hash(sequence, params)
        
        # Store in memory cache
        self.analysis_cache.put(cache_key, result, "analysis")
        
        # Store in persistent cache
        if self.enable_persistent:
            self._save_persistent_cache(cache_key, result, "analysis")
    
    def get_conservation_data(self, sequence: str) -> Optional[Dict[str, Any]]:
        """Get cached conservation analysis data"""
        cache_key = self._get_sequence_hash(sequence)
        return self.conservation_cache.get(cache_key)
    
    def store_conservation_data(self, sequence: str, data: Dict[str, Any]):
        """Store conservation analysis data"""
        cache_key = self._get_sequence_hash(sequence)
        self.conservation_cache.put(cache_key, data, "conservation")
    
    def _save_persistent_cache(self, key: str, data: Any, cache_type: str):
        """Save data to persistent cache file"""
        try:
            cache_file = os.path.join(self.cache_dir, f"{cache_type}_{key}.pkl")
            metadata_file = os.path.join(self.cache_dir, f"{cache_type}_{key}.meta")
            
            # Save data
            with open(cache_file, 'wb') as f:
                pickle.dump(data, f)
            
            # Save metadata
            metadata = {
                'created_at': datetime.now().isoformat(),
                'cache_type': cache_type,
                'ttl_hours': self.persistent_ttl_hours
            }
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f)
                
        except Exception as e:
            print(f"Warning: Could not save persistent cache: {e}")
    
    def _load_persistent_cache(self, key: str, cache_type: str) -> Optional[Any]:
        """Load data from persistent cache file"""
        try:
            cache_file = os.path.join(self.cache_dir, f"{cache_type}_{key}.pkl")
            metadata_file = os.path.join(self.cache_dir, f"{cache_type}_{key}.meta")
            
            if not os.path.exists(cache_file) or not os.path.exists(metadata_file):
                return None
            
            # Check if cache is still valid
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
            
            created_at = datetime.fromisoformat(metadata['created_at'])
            ttl = timedelta(hours=metadata.get('ttl_hours', 24))
            
            if datetime.now() - created_at > ttl:
                # Cache expired, remove files
                os.remove(cache_file)
                os.remove(metadata_file)
                return None
            
            # Load data
            with open(cache_file, 'rb') as f:
                return pickle.load(f)
                
        except Exception as e:
            print(f"Warning: Could not load persistent cache: {e}")
            return None
    
    def cleanup_persistent_cache(self):
        """Remove expired persistent cache files"""
        if not os.path.exists(self.cache_dir):
            return
        
        removed_count = 0
        for filename in os.listdir(self.cache_dir):
            if filename.endswith('.meta'):
                try:
                    metadata_file = os.path.join(self.cache_dir, filename)
                    with open(metadata_file, 'r') as f:
                        metadata = json.load(f)
                    
                    created_at = datetime.fromisoformat(metadata['created_at'])
                    ttl = timedelta(hours=metadata.get('ttl_hours', 24))
                    
                    if datetime.now() - created_at > ttl:
                        # Remove both data and metadata files
                        cache_key = filename.replace('.meta', '')
                        cache_file = os.path.join(self.cache_dir, f"{cache_key}.pkl")
                        
                        if os.path.exists(cache_file):
                            os.remove(cache_file)
                        os.remove(metadata_file)
                        removed_count += 1
                        
                except Exception:
                    pass  # Skip corrupted files
        
        if removed_count > 0:
            print(f"Cleaned up {removed_count} expired cache files")
    
    def get_comprehensive_stats(self) -> Dict[str, Any]:
        """Get comprehensive statistics for all caches"""
        stats = {
            'analysis_cache': self.analysis_cache.get_stats(),
            'sequence_cache': self.sequence_cache.get_stats(),
            'conservation_cache': self.conservation_cache.get_stats(),
            'persistent_cache': self._get_persistent_cache_stats()
        }
        
        # Calculate totals
        total_memory_mb = sum(
            cache_stats['memory_used_mb'] 
            for cache_stats in [stats['analysis_cache'], stats['sequence_cache'], stats['conservation_cache']]
        )
        total_requests = sum(
            cache_stats['total_requests'] 
            for cache_stats in [stats['analysis_cache'], stats['sequence_cache'], stats['conservation_cache']]
        )
        total_hits = sum(
            cache_stats['hits'] 
            for cache_stats in [stats['analysis_cache'], stats['sequence_cache'], stats['conservation_cache']]
        )
        
        stats['totals'] = {
            'memory_used_mb': round(total_memory_mb, 2),
            'hit_rate_percent': round((total_hits / max(total_requests, 1)) * 100, 2),
            'total_requests': total_requests,
            'total_hits': total_hits
        }
        
        return stats
    
    def _get_persistent_cache_stats(self) -> Dict[str, Any]:
        """Get statistics for persistent cache"""
        if not os.path.exists(self.cache_dir):
            return {'files': 0, 'size_mb': 0}
        
        file_count = 0
        total_size = 0
        
        for filename in os.listdir(self.cache_dir):
            if filename.endswith('.pkl'):
                file_count += 1
                filepath = os.path.join(self.cache_dir, filename)
                total_size += os.path.getsize(filepath)
        
        return {
            'files': file_count,
            'size_mb': round(total_size / (1024 * 1024), 2)
        }
    
    def clear_all_caches(self):
        """Clear all cache data"""
        self.analysis_cache.clear()
        self.sequence_cache.clear()
        self.conservation_cache.clear()
        
        # Clear persistent cache
        if os.path.exists(self.cache_dir):
            shutil.rmtree(self.cache_dir)
            os.makedirs(self.cache_dir, exist_ok=True)

# Global cache manager instance
enhanced_cache_manager = EnhancedCacheManager()

def get_cache_manager() -> EnhancedCacheManager:
    """Get the global cache manager instance"""
    return enhanced_cache_manager

def clear_all_caches():
    """Clear all NBDFinder caches"""
    enhanced_cache_manager.clear_all_caches()
    
    # Also clear Hyperscan cache
    try:
        from motifs.hyperscan_manager import clear_hyperscan_cache
        clear_hyperscan_cache()
    except ImportError:
        pass

def get_cache_stats() -> Dict[str, Any]:
    """Get comprehensive cache statistics"""
    return enhanced_cache_manager.get_comprehensive_stats()