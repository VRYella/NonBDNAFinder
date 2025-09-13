"""
Hyperscan registry manager: compile & cache DB(s), expose scan APIs.
- compile_db(patterns, key="global")
- scan_block(seq: str, key="global") -> List[ (pat_index,start0,end0) ]
- scan_stream(seq_iter, key="global") -> generator of same tuples
"""

from typing import List, Tuple, Dict, Generator, Optional
import logging

try:
    import hyperscan as hs
    HS_AVAILABLE = True
except Exception as e:
    HS_AVAILABLE = False
    logging.warning(f"Hyperscan binding not available: {e}")

# module-level cache
_DB_CACHE: Dict[str, Dict] = {}

logger = logging.getLogger(__name__)


def compile_db(patterns: List[str], key: str = "global") -> Dict:
    """
    Compile expressions (bytes or str) into a Hyperscan Database and cache under key.
    
    Args:
        patterns: list of regex strings (Hyperscan-safe)
        key: cache key for the compiled database
        
    Returns:
        cache entry { db, patterns, metadata }
    """
    if not HS_AVAILABLE:
        raise RuntimeError("Hyperscan binding not available")
    
    if key in _DB_CACHE:
        logger.debug(f"Using cached DB for key: {key}")
        return _DB_CACHE[key]
    
    if not patterns:
        raise ValueError("No patterns provided for compilation")
    
    logger.info(f"Compiling Hyperscan DB with {len(patterns)} patterns for key: {key}")
    
    # Convert patterns to bytes and prepare compilation
    exprs = []
    ids = []
    flags = []
    
    for i, pattern in enumerate(patterns):
        if isinstance(pattern, str):
            expr = pattern.encode('utf-8')
        else:
            expr = pattern
        
        exprs.append(expr)
        ids.append(i)
        flags.append(hs.HS_FLAG_UTF8 | hs.HS_FLAG_CASELESS | hs.HS_FLAG_DOTALL)
    
    try:
        # Create and compile database
        db = hs.Database()
        db.compile(expressions=exprs, ids=ids, flags=flags)
        
        cache_entry = {
            "db": db,
            "patterns": patterns,
            "count": len(patterns),
            "key": key
        }
        
        _DB_CACHE[key] = cache_entry
        logger.info(f"Successfully compiled DB for key: {key}")
        return cache_entry
        
    except Exception as e:
        logger.error(f"Failed to compile Hyperscan DB for key {key}: {e}")
        raise RuntimeError(f"Hyperscan compilation failed: {e}")


def scan_block(seq: str, key: str = "global") -> List[Tuple[int, int, int]]:
    """
    Scan full sequence once; returns list of (pat_index, start0, end0).
    
    Args:
        seq: sequence to scan
        key: database cache key
        
    Returns:
        List of matches as (pattern_index, start_pos, end_pos) tuples
    """
    if not HS_AVAILABLE:
        raise RuntimeError("Hyperscan binding not available")
    
    cache = _DB_CACHE.get(key)
    if not cache:
        raise RuntimeError(f"DB {key} not compiled. Call compile_db() first.")
    
    db = cache["db"]
    
    # Create scratch space (reuse per worker, not per match)
    scratch = hs.Scratch(db)
    matches = []
    
    def on_match(pat_id, start, end, flags, ctx):
        """Minimal callback - only append match info"""
        matches.append((int(pat_id), int(start), int(end)))
    
    try:
        # Scan the sequence
        if isinstance(seq, str):
            seq_bytes = seq.encode('utf-8')
        else:
            seq_bytes = seq
            
        db.scan(seq_bytes, match_event_handler=on_match, scratch=scratch)
        
    except Exception as e:
        logger.error(f"Hyperscan scanning failed for key {key}: {e}")
        raise RuntimeError(f"Scanning failed: {e}")
    
    logger.debug(f"Found {len(matches)} matches for key: {key}")
    return matches


def scan_stream(seq_iter: Generator[bytes, None, None], key: str = "global") -> Generator[Tuple[int, int, int], None, None]:
    """
    Use HS stream API (if available) to feed chunks and yield matches w/ absolute offsets.
    
    Args:
        seq_iter: generator yielding byte chunks
        key: database cache key
        
    Yields:
        Matches as (pattern_index, absolute_start, absolute_end) tuples
    """
    if not HS_AVAILABLE:
        raise RuntimeError("Hyperscan binding not available")
    
    cache = _DB_CACHE.get(key)
    if not cache:
        raise RuntimeError(f"DB {key} not compiled. Call compile_db() first.")
    
    db = cache["db"]
    scratch = hs.Scratch(db)
    
    try:
        # Create streaming context
        stream = hs.Stream(db)
        offset = 0
        
        def on_match(pat_id, start, end, flags, ctx):
            """Stream callback with offset adjustment"""
            return (int(pat_id), offset + int(start), offset + int(end))
        
        # Process chunks
        for chunk in seq_iter:
            if not isinstance(chunk, bytes):
                chunk = chunk.encode('utf-8')
                
            matches = []
            
            def chunk_callback(pat_id, start, end, flags, ctx):
                matches.append((int(pat_id), offset + int(start), offset + int(end)))
            
            stream.scan(chunk, match_event_handler=chunk_callback, scratch=scratch)
            
            # Yield all matches from this chunk
            for match in matches:
                yield match
            
            offset += len(chunk)
        
        # Close stream and get final matches
        final_matches = []
        
        def final_callback(pat_id, start, end, flags, ctx):
            final_matches.append((int(pat_id), int(start), int(end)))
        
        stream.close(match_event_handler=final_callback, scratch=scratch)
        
        # Yield final matches
        for match in final_matches:
            yield match
            
    except Exception as e:
        logger.error(f"Hyperscan stream scanning failed for key {key}: {e}")
        raise RuntimeError(f"Stream scanning failed: {e}")


def get_db_info(key: str = "global") -> Optional[Dict]:
    """
    Get information about a compiled database
    
    Args:
        key: database cache key
        
    Returns:
        Database info dict or None if not found
    """
    cache = _DB_CACHE.get(key)
    if not cache:
        return None
    
    return {
        "key": cache["key"],
        "pattern_count": cache["count"],
        "patterns": cache["patterns"][:5] + ["..."] if len(cache["patterns"]) > 5 else cache["patterns"]
    }


def clear_cache(key: Optional[str] = None):
    """
    Clear database cache
    
    Args:
        key: specific key to clear, or None to clear all
    """
    global _DB_CACHE
    
    if key is None:
        _DB_CACHE.clear()
        logger.info("Cleared all cached databases")
    elif key in _DB_CACHE:
        del _DB_CACHE[key]
        logger.info(f"Cleared cached database for key: {key}")
    else:
        logger.warning(f"No cached database found for key: {key}")


def is_available() -> bool:
    """Check if Hyperscan is available"""
    return HS_AVAILABLE


def get_cached_keys() -> List[str]:
    """Get list of cached database keys"""
    return list(_DB_CACHE.keys())


# Performance notes for documentation:
"""
HYPERSCAN DESIGN RULES (enforce in code review):

1. COMPILE ONCE: Database compilation is expensive. Always cache compiled DBs
   using the module-level _DB_CACHE. Each compile_db() call should check cache first.

2. SCRATCH PER WORKER: Create one hs.Scratch(db) per worker/thread and reuse it.
   Never create scratch objects inside tight loops or per-match callbacks.

3. MINIMAL CALLBACKS: Hyperscan callbacks should only append (pat_id, start, end).
   Heavy Python processing (validation, scoring) must be done after scanning.

4. MEMORY MANAGEMENT: For large sequences, prefer scan_stream() over scan_block()
   to process data in chunks and avoid loading entire sequences into memory.

5. PATTERN GROUPING: Group patterns by motif class for better cache locality.
   Use separate DBs for different analysis types to optimize memory usage.
"""