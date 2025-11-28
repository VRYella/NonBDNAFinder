"""
Parallel Worker with Shared Memory for Non-B DNA Scanning
==========================================================

This module provides parallel scanning capabilities using multiprocessing
with shared memory. Key features:

- SharedMemory-based sequence sharing (avoids pickling large sequences)
- Overlapping chunk processing for boundary-aware scanning
- Worker pool management with result aggregation

Usage:
    from scanner_backends.parallel_worker import parallel_scan
    
    motifs = parallel_scan(sequence, num_workers=4, chunk_size=100000)
"""

import os
import sys
import multiprocessing as mp
from multiprocessing import shared_memory
from typing import List, Dict, Any, Optional, Tuple, Callable
import numpy as np
from functools import partial


# Default configuration
DEFAULT_CHUNK_SIZE = 100_000  # 100KB chunks
DEFAULT_OVERLAP = 500  # 500bp overlap for motif spanning
DEFAULT_NUM_WORKERS = None  # None means use all available CPUs


def chunk_sequence(
    sequence_length: int,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    overlap: int = DEFAULT_OVERLAP
) -> List[Tuple[int, int]]:
    """
    Calculate chunk boundaries for parallel processing.
    
    Args:
        sequence_length: Total length of sequence in bp
        chunk_size: Size of each chunk in bp
        overlap: Overlap between chunks to avoid missing boundary motifs
        
    Returns:
        List of (start, end) tuples defining chunk boundaries
    """
    if sequence_length <= chunk_size:
        return [(0, sequence_length)]
    
    chunks = []
    step = chunk_size - overlap
    pos = 0
    
    while pos < sequence_length:
        chunk_end = min(pos + chunk_size, sequence_length)
        chunks.append((pos, chunk_end))
        pos += step
        
        if pos >= sequence_length:
            break
    
    return chunks


class SharedMemoryWorker:
    """
    Worker class for scanning sequence chunks using shared memory.
    
    This class manages a shared memory block containing the sequence,
    allowing multiple workers to access it without copying.
    """
    
    def __init__(
        self,
        sequence: str,
        name: Optional[str] = None
    ):
        """
        Initialize shared memory worker with sequence.
        
        Args:
            sequence: DNA sequence string
            name: Optional name for the shared memory block
        """
        self.sequence = sequence.upper()
        self.sequence_bytes = self.sequence.encode('utf-8')
        self.sequence_length = len(self.sequence)
        
        # Create shared memory block
        self.shm_name = name or f"nbdscan_{os.getpid()}"
        self.shm = None
        self._create_shared_memory()
    
    def _create_shared_memory(self):
        """Create shared memory block and copy sequence into it."""
        try:
            # Create shared memory block
            self.shm = shared_memory.SharedMemory(
                name=self.shm_name,
                create=True,
                size=len(self.sequence_bytes)
            )
            
            # Copy sequence into shared memory
            self.shm.buf[:len(self.sequence_bytes)] = self.sequence_bytes
            
        except FileExistsError:
            # If already exists, try to unlink and recreate
            try:
                old_shm = shared_memory.SharedMemory(name=self.shm_name, create=False)
                old_shm.close()
                old_shm.unlink()
            except (FileNotFoundError, PermissionError) as e:
                # Ignore if already cleaned up or no permission
                pass
            except Exception as e:
                # Log unexpected exceptions but continue
                import logging
                logging.getLogger(__name__).debug(f"Shared memory cleanup warning: {e}")
            
            self.shm = shared_memory.SharedMemory(
                name=self.shm_name,
                create=True,
                size=len(self.sequence_bytes)
            )
            self.shm.buf[:len(self.sequence_bytes)] = self.sequence_bytes
    
    def get_chunk(self, start: int, end: int) -> str:
        """
        Get a chunk of sequence from shared memory.
        
        Args:
            start: Start position (0-based)
            end: End position (exclusive)
            
        Returns:
            Sequence chunk as string
        """
        chunk_bytes = bytes(self.shm.buf[start:end])
        return chunk_bytes.decode('utf-8')
    
    def cleanup(self):
        """Release shared memory resources."""
        if self.shm is not None:
            try:
                self.shm.close()
                self.shm.unlink()
            except Exception:
                pass
            self.shm = None
    
    def __del__(self):
        """Destructor to clean up shared memory."""
        self.cleanup()
    
    def __enter__(self):
        return self
    
    def __exit__(self, *args):
        self.cleanup()


def _worker_scan_chunk(
    args: Tuple[str, int, int, str],
    scan_function: Callable = None
) -> List[Dict[str, Any]]:
    """
    Worker function to scan a single chunk.
    
    This function is called by worker processes and accesses the
    shared memory block to get the sequence chunk.
    
    Args:
        args: Tuple of (shm_name, start, end, sequence_name)
        scan_function: Function to use for scanning (default: analyze_sequence)
        
    Returns:
        List of motifs with adjusted positions
    """
    shm_name, start, end, sequence_name = args
    
    try:
        # Attach to shared memory
        shm = shared_memory.SharedMemory(name=shm_name, create=False)
        
        # Extract chunk
        chunk_bytes = bytes(shm.buf[start:end])
        chunk_str = chunk_bytes.decode('utf-8')
        
        # Close shared memory reference (don't unlink - main process owns it)
        shm.close()
        
        # Scan chunk - use provided function or import default
        # Note: We import inside the function to avoid issues with multiprocessing
        # and to allow custom scan functions to be provided
        if scan_function is None:
            from nonbscanner import analyze_sequence
            scan_function = analyze_sequence
        
        motifs = scan_function(chunk_str, f"{sequence_name}_chunk_{start}")
        
        # Adjust positions to global coordinates
        for motif in motifs:
            if 'Start' in motif:
                motif['Start'] += start
            if 'End' in motif:
                motif['End'] += start
            motif['Sequence_Name'] = sequence_name
            motif['Chunk_Start'] = start
        
        return motifs
        
    except Exception as e:
        # Return error info instead of crashing worker
        return [{'error': str(e), 'chunk_start': start, 'chunk_end': end}]


def parallel_scan(
    sequence: str,
    sequence_name: str = "sequence",
    num_workers: Optional[int] = DEFAULT_NUM_WORKERS,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    overlap: int = DEFAULT_OVERLAP,
    scan_function: Callable = None,
    progress_callback: Callable = None
) -> List[Dict[str, Any]]:
    """
    Scan sequence in parallel using multiple worker processes.
    
    This function:
    1. Creates shared memory with the sequence
    2. Divides sequence into overlapping chunks
    3. Dispatches chunks to worker processes
    4. Merges results and removes duplicates from overlap regions
    
    Args:
        sequence: DNA sequence string
        sequence_name: Name for the sequence
        num_workers: Number of parallel workers (default: CPU count)
        chunk_size: Size of each chunk in bp
        overlap: Overlap between chunks
        scan_function: Custom scanning function (default: analyze_sequence)
        progress_callback: Function to call with (completed, total) counts
        
    Returns:
        List of motifs sorted by position
        
    Example:
        >>> motifs = parallel_scan(sequence, num_workers=4)
        >>> print(f"Found {len(motifs)} motifs")
    """
    if num_workers is None:
        num_workers = mp.cpu_count()
    
    sequence_length = len(sequence)
    
    # For small sequences, just use single process
    if sequence_length < chunk_size * 2 or num_workers == 1:
        if scan_function is None:
            from nonbscanner import analyze_sequence
            scan_function = analyze_sequence
        return scan_function(sequence, sequence_name)
    
    # Calculate chunk boundaries
    chunks = chunk_sequence(sequence_length, chunk_size, overlap)
    
    if progress_callback:
        progress_callback(0, len(chunks))
    
    # Create shared memory worker
    with SharedMemoryWorker(sequence) as worker:
        shm_name = worker.shm_name
        
        # Prepare worker arguments
        work_args = [
            (shm_name, start, end, sequence_name)
            for start, end in chunks
        ]
        
        # Process chunks in parallel
        all_motifs = []
        
        # Use multiprocessing pool
        with mp.Pool(processes=num_workers) as pool:
            # Create partial function with scan_function
            worker_func = partial(_worker_scan_chunk, scan_function=scan_function)
            
            # Map work to pool
            results = pool.map(worker_func, work_args)
            
            # Collect results
            for i, chunk_motifs in enumerate(results):
                if progress_callback:
                    progress_callback(i + 1, len(chunks))
                
                # Filter out error entries and collect valid motifs
                for motif in chunk_motifs:
                    if 'error' not in motif:
                        all_motifs.append(motif)
    
    # Remove duplicates from overlap regions
    all_motifs = _deduplicate_overlap_motifs(all_motifs, overlap)
    
    # Sort by position
    all_motifs.sort(key=lambda m: m.get('Start', 0))
    
    # Reassign IDs
    for i, motif in enumerate(all_motifs):
        motif['ID'] = f"{sequence_name}_motif_{i+1}"
    
    return all_motifs


def _deduplicate_overlap_motifs(
    motifs: List[Dict[str, Any]],
    overlap: int
) -> List[Dict[str, Any]]:
    """
    Remove duplicate motifs from overlap regions.
    
    When chunks overlap, the same motif may be detected in both chunks.
    This function removes duplicates by keeping the motif with higher score.
    
    Args:
        motifs: List of motifs potentially with duplicates
        overlap: Overlap size used in chunking
        
    Returns:
        Deduplicated list of motifs
    """
    if not motifs:
        return []
    
    # Group motifs by (Class, Subclass, Start, End)
    seen = {}
    
    for motif in motifs:
        key = (
            motif.get('Class', ''),
            motif.get('Subclass', ''),
            motif.get('Start', 0),
            motif.get('End', 0)
        )
        
        if key in seen:
            # Keep the one with higher score
            if motif.get('Score', 0) > seen[key].get('Score', 0):
                seen[key] = motif
        else:
            seen[key] = motif
    
    return list(seen.values())


def get_system_info() -> Dict[str, Any]:
    """
    Get system information relevant for parallel processing.
    
    Returns:
        Dictionary with CPU count, memory info, etc.
    """
    info = {
        'cpu_count': mp.cpu_count(),
        'platform': sys.platform,
    }
    
    try:
        import psutil
        mem = psutil.virtual_memory()
        info['total_memory_gb'] = round(mem.total / (1024**3), 2)
        info['available_memory_gb'] = round(mem.available / (1024**3), 2)
    except ImportError:
        pass
    
    return info


def suggest_num_workers(sequence_length: int) -> int:
    """
    Suggest optimal number of workers based on sequence length and system.
    
    Args:
        sequence_length: Length of sequence in bp
        
    Returns:
        Suggested number of workers
    """
    cpu_count = mp.cpu_count()
    
    # For very small sequences, use fewer workers
    if sequence_length < 50_000:
        return 1
    elif sequence_length < 200_000:
        return min(2, cpu_count)
    elif sequence_length < 1_000_000:
        return min(4, cpu_count)
    else:
        # Use most CPUs for large sequences, leave 1 for main process
        return max(1, cpu_count - 1)
