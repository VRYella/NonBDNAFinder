"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                  HYPERSCAN PARALLEL SCANNER AGENT                             ║
║          Memory-Efficient Parallel Motif Location Detection                  ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: scanner_agent.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.2 - Parallel Hyperscan Architecture
LICENSE: MIT

DESCRIPTION:
    Dedicated scanning agent responsible only for finding raw motif locations.
    Uses Hyperscan (when available) and parallel processing for optimal performance.
    Designed to work within Streamlit's 1GB memory limit.

COORDINATE SYSTEM:
    This module uses 1-based INCLUSIVE coordinates for input/output to match the
    existing NonBScanner convention:
    - Start: 1-based position (first base = 1)
    - End: 1-based position (INCLUSIVE - last base in motif)
    - Length: End - Start + 1
    
    Internally, Python range() uses 0-based half-open intervals:
    - range(Start-1, End) covers positions Start through End inclusive
    
    Example: Start=1, End=20 represents 20 bases (positions 1-20 inclusive)
             Converted to range(0, 20) for Python indexing

ARCHITECTURE:
    - Chunks genome into overlapping segments (1000nt overlap)
    - Parallel processing using multiprocessing.Pool
    - Deduplication of overlapping matches
    - No scoring logic (delegated to existing ScoringEngine)

PERFORMANCE:
    - Chunk size: 50,000 bp (configurable)
    - Overlap: 1,000 bp (handles motifs at boundaries)
    - Worker processes: CPU count
    - Progress reporting: Per-chunk completion

USAGE:
    from scanner_agent import ParallelScanner
    
    scanner = ParallelScanner(genome_sequence, hs_database)
    raw_motifs = scanner.run_scan()
"""

import multiprocessing as mp
from typing import List, Dict, Tuple, Optional, Any, Callable
import numpy as np
from collections import defaultdict

# Try to import Hyperscan (optional dependency)
try:
    import hyperscan
    HYPERSCAN_AVAILABLE = True
except ImportError:
    HYPERSCAN_AVAILABLE = False
    hyperscan = None

# Configuration constants
CHUNK_SIZE = 50000  # 50kb chunks for memory efficiency
OVERLAP_SIZE = 1000  # 1kb overlap to catch motifs at chunk boundaries
MIN_MOTIF_LENGTH = 4  # Minimum motif length to consider


def hs_worker_task(args: Tuple[int, np.ndarray, Any]) -> List[Tuple[int, int, int]]:
    """
    Worker function for parallel Hyperscan scanning.
    
    This is the core worker that:
    1. Takes a chunk of the genome and its global offset
    2. Runs Hyperscan database in Block Mode
    3. Converts local coordinates to global coordinates
    4. Returns raw (start, end, pattern_id) tuples
    
    Args:
        args: Tuple of (offset, genome_chunk, hs_database)
            - offset: Global starting position of this chunk
            - genome_chunk: NumPy byte array of sequence chunk
            - hs_database: Compiled Hyperscan database (or None for fallback)
    
    Returns:
        List of (global_start, global_end, pattern_id) tuples
    """
    offset, genome_chunk, hs_db = args
    local_results = []
    
    # Convert numpy array to bytes for Hyperscan
    chunk_bytes = genome_chunk.tobytes()
    
    if HYPERSCAN_AVAILABLE and hs_db is not None:
        # Hyperscan match handler (callback function)
        def match_handler(id, start, end, flags, context):
            """Callback for Hyperscan matches"""
            # Convert local match coordinates to global coordinates
            global_start = offset + start
            global_end = offset + end
            local_results.append((global_start, global_end, id))
            return 0  # Continue scanning
        
        try:
            # Execute Hyperscan in Block Mode
            hs_db.scan(chunk_bytes, match_handler)
        except Exception as e:
            # If Hyperscan fails, log but don't crash
            print(f"Warning: Hyperscan scan failed at offset {offset}: {e}")
            return local_results
    else:
        # Fallback: Use regex-based detection from existing scanner
        # This maintains compatibility when Hyperscan is not available
        try:
            # Use absolute import to avoid sys.path manipulation
            from scanner import analyze_sequence
            
            # Decode chunk to string
            chunk_str = chunk_bytes.decode('utf-8', errors='ignore')
            
            # Run detection on chunk
            chunk_motifs = analyze_sequence(chunk_str, f"chunk_{offset}")
            
            # Convert to (start, end, pattern_id) format
            # COORDINATE SYSTEM NOTE:
            # - Input motifs use 1-based INCLUSIVE coordinates (Start=1, End=20 = 20 bases)
            # - We convert to 0-based half-open for internal tuple format
            # - Conversion: Start-1 (to 0-based), End (already exclusive in half-open)
            for i, motif in enumerate(chunk_motifs):
                # Start: 1-based inclusive → 0-based (subtract 1)
                # End: 1-based inclusive → 0-based exclusive (no change, becomes half-open)
                global_start = offset + motif.get('Start', 0) - 1
                global_end = offset + motif.get('End', 0)  # Becomes 0-based half-open
                pattern_id = hash(f"{motif.get('Class', '')}_{motif.get('Subclass', '')}")
                local_results.append((global_start, global_end, pattern_id))
        except Exception as e:
            # Silently skip if detection fails - worker shouldn't crash the whole job
            pass
    
    return local_results


class ParallelScanner:
    """
    Parallel scanning agent for high-performance motif detection.
    
    This class orchestrates the parallel scanning process:
    - Chunks the genome with overlap
    - Distributes work to multiple CPU cores
    - Aggregates and deduplicates results
    - Reports progress dynamically
    
    Attributes:
        genome: Full genome sequence as NumPy byte array
        hs_db: Compiled Hyperscan database (optional)
        chunk_size: Size of each chunk in base pairs
        overlap_size: Overlap between chunks in base pairs
        num_workers: Number of parallel worker processes
    """
    
    def __init__(self, 
                 genome: str, 
                 hs_db: Optional[Any] = None,
                 chunk_size: int = CHUNK_SIZE,
                 overlap_size: int = OVERLAP_SIZE,
                 num_workers: Optional[int] = None):
        """
        Initialize the parallel scanner.
        
        Args:
            genome: DNA sequence string (will be converted to NumPy array)
            hs_db: Compiled Hyperscan database (optional, uses fallback if None)
            chunk_size: Size of each chunk (default: 50kb)
            overlap_size: Overlap between chunks (default: 1kb)
            num_workers: Number of worker processes (default: CPU count)
        """
        # Convert genome to NumPy byte array for efficient chunking
        self.genome_array = np.frombuffer(genome.encode('utf-8'), dtype=np.uint8)
        self.genome_length = len(self.genome_array)
        self.hs_db = hs_db
        self.chunk_size = chunk_size
        self.overlap_size = overlap_size
        self.num_workers = num_workers or mp.cpu_count()
        
        # Calculate number of chunks
        self.num_chunks = self._calculate_num_chunks()
    
    def _calculate_num_chunks(self) -> int:
        """Calculate the number of chunks needed to cover the genome."""
        effective_chunk_size = self.chunk_size - self.overlap_size
        if effective_chunk_size <= 0:
            raise ValueError("Chunk size must be larger than overlap size")
        
        num_chunks = (self.genome_length + effective_chunk_size - 1) // effective_chunk_size
        return max(1, num_chunks)
    
    def _create_tasks(self) -> List[Tuple[int, np.ndarray, Any]]:
        """
        Create task list for parallel processing.
        
        Each task is a tuple of (offset, chunk_array, hs_db).
        Chunks overlap by OVERLAP_SIZE to ensure motifs at boundaries are found.
        
        Returns:
            List of (offset, genome_chunk, hs_db) tuples for workers
        """
        tasks = []
        effective_chunk_size = self.chunk_size - self.overlap_size
        
        for i in range(self.num_chunks):
            # Calculate chunk boundaries
            chunk_start = i * effective_chunk_size
            chunk_end = min(chunk_start + self.chunk_size, self.genome_length)
            
            # Extract chunk (includes overlap with next chunk)
            chunk_array = self.genome_array[chunk_start:chunk_end]
            
            # Create task tuple
            tasks.append((chunk_start, chunk_array, self.hs_db))
        
        return tasks
    
    def run_scan(self, progress_callback: Optional[Callable[[int, int], None]] = None) -> List[Tuple[int, int, int]]:
        """
        Execute parallel scanning across all chunks.
        
        This method:
        1. Creates chunking tasks
        2. Distributes to worker pool
        3. Collects results
        4. Deduplicates overlapping matches
        5. Returns raw motif locations
        
        Args:
            progress_callback: Optional callback function(current, total) for progress
        
        Returns:
            List of unique (start, end, pattern_id) tuples
        """
        tasks = self._create_tasks()
        all_raw_results = []
        
        # Use multiprocessing.Pool for parallel execution
        with mp.Pool(processes=self.num_workers) as pool:
            # Process tasks and collect results
            for i, result in enumerate(pool.imap(hs_worker_task, tasks)):
                all_raw_results.append(result)
                
                # Report progress if callback provided
                if progress_callback:
                    progress_callback(i + 1, len(tasks))
        
        # Deduplicate results from overlapping regions
        unique_motifs = self._deduplicate(all_raw_results)
        
        return unique_motifs
    
    def _deduplicate(self, results_list: List[List[Tuple[int, int, int]]]) -> List[Tuple[int, int, int]]:
        """
        Deduplicate motifs found in overlapping regions.
        
        Uses a set-based approach where duplicates are identified by 
        identical (start, end, pattern_id) tuples.
        
        Args:
            results_list: List of result lists from each worker
        
        Returns:
            Deduplicated list of (start, end, pattern_id) tuples
        """
        # Flatten the list of lists
        all_motifs = []
        for result_chunk in results_list:
            all_motifs.extend(result_chunk)
        
        # Use set for deduplication (handles exact duplicates)
        # For motifs in overlap regions, they'll have identical coordinates
        unique_motifs_set = set(all_motifs)
        
        # Convert back to sorted list
        unique_motifs = sorted(list(unique_motifs_set), key=lambda x: (x[0], x[1]))
        
        return unique_motifs
    
    def get_statistics(self) -> Dict[str, Any]:
        """
        Get scanning statistics.
        
        Returns:
            Dictionary with scan configuration and statistics
        """
        return {
            'genome_length': self.genome_length,
            'chunk_size': self.chunk_size,
            'overlap_size': self.overlap_size,
            'num_chunks': self.num_chunks,
            'num_workers': self.num_workers,
            'hyperscan_available': HYPERSCAN_AVAILABLE,
            'using_hyperscan': HYPERSCAN_AVAILABLE and self.hs_db is not None
        }


# Convenience function for single-call scanning
def scan_genome_parallel(genome: str, 
                         hs_db: Optional[Any] = None,
                         progress_callback: Optional[Callable[[int, int], None]] = None) -> List[Tuple[int, int, int]]:
    """
    Convenience function for parallel genome scanning.
    
    Args:
        genome: DNA sequence string
        hs_db: Compiled Hyperscan database (optional)
        progress_callback: Optional callback function(current, total)
    
    Returns:
        List of (start, end, pattern_id) tuples
    """
    scanner = ParallelScanner(genome, hs_db)
    return scanner.run_scan(progress_callback=progress_callback)


# Example usage and testing
if __name__ == "__main__":
    import time
    
    # Generate test sequence
    test_genome = "GGGTTAGGGTTAGGGTTAGGG" * 1000 + "AAAAATTTTAAAAATTTT" * 500
    
    print("=" * 70)
    print("Parallel Scanner Agent - Test Run")
    print("=" * 70)
    print(f"Genome length: {len(test_genome):,} bp")
    print(f"Hyperscan available: {HYPERSCAN_AVAILABLE}")
    
    # Progress callback
    def progress_cb(current, total):
        percent = (current / total) * 100
        print(f"Progress: {current}/{total} chunks ({percent:.1f}%)", end='\r')
    
    # Run scan
    scanner = ParallelScanner(test_genome, hs_db=None)
    stats = scanner.get_statistics()
    
    print(f"\nScan configuration:")
    for key, value in stats.items():
        print(f"  {key}: {value}")
    
    print(f"\nStarting parallel scan...")
    start_time = time.time()
    
    motifs = scanner.run_scan(progress_callback=progress_cb)
    
    scan_time = time.time() - start_time
    
    print(f"\n\nResults:")
    print(f"  Motifs found: {len(motifs)}")
    print(f"  Scan time: {scan_time:.3f}s")
    print(f"  Speed: {len(test_genome)/scan_time:,.0f} bp/s")
    
    # Show first few motifs
    if motifs:
        print(f"\nFirst 5 motifs:")
        for i, (start, end, pattern_id) in enumerate(motifs[:5]):
            print(f"  {i+1}. Start: {start}, End: {end}, Pattern ID: {pattern_id}")
