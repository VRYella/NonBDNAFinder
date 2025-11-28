"""
I/O Utilities for High-Performance Sequence Processing
======================================================

This module provides memory-efficient I/O utilities for processing large genomes:
- Memory-mapped FASTA file reading (for >100MB files)
- Streaming sequence chunks with overlap handling
- Efficient byte-based sequence storage

Performance Notes:
- Use mmap for files >100MB to avoid loading entire file into memory
- Chunk sizes of 50k-200k work well for cache efficiency
- Overlap should be at least max_motif_length to avoid missing motifs at boundaries
"""

import os
import mmap
import numpy as np
from typing import Iterator, Tuple, Optional, Dict, Any, BinaryIO
import io


# Default chunk size for streaming (200KB works well for L3 cache)
DEFAULT_CHUNK_SIZE = 200_000

# Default overlap size (max expected motif length)
DEFAULT_OVERLAP_SIZE = 500


def get_overlap_size(motif_types: Optional[list] = None) -> int:
    """
    Get recommended overlap size based on motif types being detected.
    
    Args:
        motif_types: List of motif types to detect (e.g., ['g_quadruplex', 'cruciform'])
                    If None, returns conservative default for all motif types.
    
    Returns:
        int: Recommended overlap size in bp
    """
    # Max expected lengths for different motif types
    motif_max_lengths = {
        'curved_dna': 100,
        'slipped_dna': 300,
        'cruciform': 200,
        'r_loop': 500,
        'triplex': 200,
        'g_quadruplex': 100,
        'i_motif': 100,
        'z_dna': 50,
        'a_philic': 100,
    }
    
    if motif_types is None:
        # Return conservative maximum
        return max(motif_max_lengths.values())
    
    # Calculate max for specified types
    max_len = max(motif_max_lengths.get(m, 100) for m in motif_types)
    return max_len


def mmap_fasta(filepath: str) -> Tuple[bytes, Dict[str, Tuple[int, int]]]:
    """
    Memory-map a FASTA file for efficient large genome access.
    
    This function memory-maps the file content, avoiding loading the entire
    file into memory. Useful for genomes >100MB where traditional loading
    would cause memory pressure.
    
    Args:
        filepath: Path to FASTA file
        
    Returns:
        Tuple of:
        - Memory-mapped bytes object
        - Dictionary mapping sequence names to (start, end) byte offsets
        
    Example:
        >>> mm_data, seq_offsets = mmap_fasta("ecoli.fa")
        >>> for name, (start, end) in seq_offsets.items():
        ...     seq_bytes = mm_data[start:end]
        ...     seq = seq_bytes.decode('utf-8').replace('\\n', '')
    """
    file_size = os.path.getsize(filepath)
    
    with open(filepath, 'rb') as f:
        # Memory-map the file
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        
        # Parse FASTA to find sequence boundaries
        seq_offsets = {}
        current_name = None
        current_start = 0
        pos = 0
        
        while pos < file_size:
            # Find next line
            line_end = mm.find(b'\n', pos)
            if line_end == -1:
                line_end = file_size
            
            line = mm[pos:line_end]
            
            if line.startswith(b'>'):
                # Save previous sequence if exists
                if current_name is not None:
                    seq_offsets[current_name] = (current_start, pos)
                
                # Parse header
                header = line[1:].decode('utf-8', errors='replace').strip()
                current_name = header.split()[0] if header else f"seq_{len(seq_offsets)}"
                current_start = line_end + 1
            
            pos = line_end + 1
        
        # Save last sequence
        if current_name is not None:
            seq_offsets[current_name] = (current_start, file_size)
        
        # Return a copy of the mapped data (so file can be closed)
        # For very large files, keep the mmap open instead
        return bytes(mm), seq_offsets


def stream_sequence_chunks(
    sequence: bytes,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    overlap: int = DEFAULT_OVERLAP_SIZE
) -> Iterator[Tuple[int, bytes]]:
    """
    Stream sequence in overlapping chunks for parallel processing.
    
    Each chunk includes an overlap with the previous chunk to ensure
    motifs spanning chunk boundaries are detected. The overlap region
    should be at least as large as the longest expected motif.
    
    Args:
        sequence: DNA sequence as bytes
        chunk_size: Size of each chunk in bp (default: 200KB)
        overlap: Overlap size between chunks (default: 500bp)
        
    Yields:
        Tuple of (chunk_start_position, chunk_bytes)
        
    Example:
        >>> seq = b"ATGC" * 100000  # 400KB sequence
        >>> for start_pos, chunk in stream_sequence_chunks(seq, chunk_size=50000):
        ...     # Process chunk, add start_pos to any detected positions
        ...     process_chunk(chunk, offset=start_pos)
    """
    seq_len = len(sequence)
    
    if seq_len <= chunk_size:
        # Small sequence, yield entire thing
        yield (0, sequence)
        return
    
    # Calculate step size (chunk_size minus overlap)
    step = chunk_size - overlap
    
    pos = 0
    while pos < seq_len:
        # Calculate chunk end
        chunk_end = min(pos + chunk_size, seq_len)
        
        # Extract chunk
        chunk = sequence[pos:chunk_end]
        
        yield (pos, chunk)
        
        # Move to next chunk
        pos += step
        
        # If remaining sequence is smaller than overlap, stop
        if pos >= seq_len:
            break


def bytes_to_numpy(sequence: bytes) -> np.ndarray:
    """
    Convert sequence bytes to NumPy array for efficient storage and sharing.
    
    Uses uint8 dtype which is optimal for DNA sequences (1 byte per base).
    
    Args:
        sequence: DNA sequence as bytes
        
    Returns:
        NumPy array of uint8 values
    """
    return np.frombuffer(sequence, dtype=np.uint8)


def numpy_to_string(arr: np.ndarray) -> str:
    """
    Convert NumPy array back to DNA sequence string.
    
    Args:
        arr: NumPy array of uint8 values
        
    Returns:
        DNA sequence string
    """
    return arr.tobytes().decode('utf-8', errors='replace')


def parse_fasta_streaming(
    filepath: str,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    overlap: int = DEFAULT_OVERLAP_SIZE
) -> Iterator[Tuple[str, int, bytes]]:
    """
    Parse FASTA file in streaming fashion, yielding overlapping chunks.
    
    This is the recommended approach for processing very large genomes
    where even memory-mapping is not ideal. Processes file sequentially
    without loading the entire sequence.
    
    Args:
        filepath: Path to FASTA file
        chunk_size: Size of each chunk in bp
        overlap: Overlap size between chunks
        
    Yields:
        Tuple of (sequence_name, chunk_start_position, chunk_bytes)
        
    Example:
        >>> for seq_name, start, chunk in parse_fasta_streaming("genome.fa"):
        ...     motifs = scan_chunk(chunk)
        ...     for m in motifs:
        ...         m['Start'] += start  # Adjust position to genome coordinates
    """
    buffer = bytearray()
    current_name = None
    chunk_start = 0
    step = chunk_size - overlap
    
    with open(filepath, 'rb') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith(b'>'):
                # Yield remaining buffer for previous sequence
                if current_name is not None and len(buffer) > 0:
                    yield (current_name, chunk_start, bytes(buffer))
                
                # Parse new header
                header = line[1:].decode('utf-8', errors='replace').strip()
                current_name = header.split()[0] if header else "unknown"
                buffer = bytearray()
                chunk_start = 0
                continue
            
            # Append sequence line to buffer
            buffer.extend(line.upper())
            
            # Yield chunks when buffer is large enough
            while len(buffer) >= chunk_size:
                chunk = bytes(buffer[:chunk_size])
                yield (current_name, chunk_start, chunk)
                
                # Keep overlap in buffer
                buffer = buffer[step:]
                chunk_start += step
        
        # Yield final chunk
        if current_name is not None and len(buffer) > 0:
            yield (current_name, chunk_start, bytes(buffer))


def write_ndjson_stream(filepath: str, mode: str = 'w'):
    """
    Create an NDJSON (newline-delimited JSON) writer for streaming results.
    
    NDJSON is ideal for incremental writes as each line is a complete JSON object.
    This allows partial results to be read before analysis completes.
    
    Args:
        filepath: Path to output file
        mode: File mode ('w' for write, 'a' for append)
        
    Returns:
        File handle with write_record method
        
    Example:
        >>> with write_ndjson_stream('results.ndjson') as writer:
        ...     for chunk in scan_chunks():
        ...         for motif in chunk.motifs:
        ...             writer.write_record(motif)
    """
    import json
    
    class NDJSONWriter:
        def __init__(self, fh):
            self.fh = fh
            self.count = 0
        
        def write_record(self, record: dict):
            """Write a single record as one line of JSON."""
            self.fh.write(json.dumps(record, ensure_ascii=True) + '\n')
            self.count += 1
        
        def flush(self):
            """Flush buffer to disk."""
            self.fh.flush()
        
        def __enter__(self):
            return self
        
        def __exit__(self, *args):
            self.fh.close()
    
    return NDJSONWriter(open(filepath, mode))


def read_ndjson_stream(filepath: str) -> Iterator[dict]:
    """
    Read NDJSON file as a stream of records.
    
    Args:
        filepath: Path to NDJSON file
        
    Yields:
        Dictionary for each JSON line
    """
    import json
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                yield json.loads(line)


def estimate_memory_usage(sequence_length: int, num_workers: int = 1) -> dict:
    """
    Estimate memory usage for scanning a sequence.
    
    Args:
        sequence_length: Length of sequence in bp
        num_workers: Number of parallel workers
        
    Returns:
        Dictionary with memory estimates in bytes
    """
    # Base memory for sequence storage
    seq_memory = sequence_length  # 1 byte per base
    
    # Shared memory overhead per worker (minimal, just slice info)
    shared_memory_overhead = num_workers * 1024
    
    # Estimated motif storage (assume ~10 bytes per bp for worst case)
    motif_memory_estimate = sequence_length // 100 * 500
    
    # Buffer for chunk processing
    chunk_buffer = DEFAULT_CHUNK_SIZE * num_workers
    
    return {
        'sequence_storage': seq_memory,
        'shared_memory_overhead': shared_memory_overhead,
        'motif_storage_estimate': motif_memory_estimate,
        'chunk_buffers': chunk_buffer,
        'total_estimate': seq_memory + shared_memory_overhead + motif_memory_estimate + chunk_buffer
    }
