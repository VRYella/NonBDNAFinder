"""
FASTA Parsing Module
====================

Functions for parsing FASTA format sequences.
Extracted from utilities.py for focused file format handling.
"""

from typing import Dict, Tuple
import re


def parse_fasta(fasta_content: str) -> Dict[str, str]:
    """
    Parse FASTA format content into dictionary.
    
    Args:
        fasta_content: FASTA format string
        
    Returns:
        Dictionary mapping sequence_name -> sequence
        
    Example:
        >>> content = ">seq1\\nATGC\\n>seq2\\nGCTA"
        >>> parse_fasta(content)
        {'seq1': 'ATGC', 'seq2': 'GCTA'}
    """
    sequences = {}
    current_name = None
    current_seq = []
    
    for line in fasta_content.split('\n'):
        line = line.strip()
        
        if line.startswith('>'):
            # Save previous sequence if exists
            if current_name is not None:
                sequences[current_name] = ''.join(current_seq)
            
            # Start new sequence
            current_name = line[1:].split()[0]  # Take first word after '>'
            current_seq = []
        elif line and current_name is not None:
            current_seq.append(line)
    
    # Save last sequence
    if current_name is not None:
        sequences[current_name] = ''.join(current_seq)
    
    return sequences


def read_fasta_file(filename: str) -> Dict[str, str]:
    """
    Read and parse FASTA file.
    
    Args:
        filename: Path to FASTA file
        
    Returns:
        Dictionary mapping sequence_name -> sequence
        
    Example:
        >>> sequences = read_fasta_file("sequences.fasta")
        >>> len(sequences)
        2
    """
    with open(filename, 'r') as f:
        content = f.read()
    
    return parse_fasta(content)


def format_fasta(sequences: Dict[str, str], line_width: int = 80) -> str:
    """
    Format sequences as FASTA with line wrapping.
    
    Args:
        sequences: Dictionary of name -> sequence
        line_width: Maximum characters per line
        
    Returns:
        FASTA formatted string
        
    Example:
        >>> seqs = {'seq1': 'ATGC' * 30}
        >>> fasta = format_fasta(seqs, line_width=60)
        >>> fasta.startswith('>seq1\\n')
        True
    """
    lines = []
    
    for name, seq in sequences.items():
        lines.append(f'>{name}')
        
        # Wrap sequence at line_width
        for i in range(0, len(seq), line_width):
            lines.append(seq[i:i + line_width])
    
    return '\n'.join(lines)


def parse_fasta_chunked(file_object, chunk_size_mb: int = 2):
    """
    Memory-efficient FASTA parser using chunked reading for large files.
    Yields (name, sequence) tuples one at a time to avoid loading entire file in memory.
    
    This version is optimized to use less memory than loading the entire file at once
    by processing sequences one at a time and yielding them immediately.
    
    Args:
        file_object: File-like object (from st.file_uploader or open())
        chunk_size_mb: Size of read chunks in MB (default: 2MB - optimized for memory efficiency)
        
    Yields:
        Tuple of (sequence_name, sequence_string)
    """
    chunk_size = chunk_size_mb * 1024 * 1024  # Convert to bytes
    current_name = None
    current_seq_parts = []  # Use list to accumulate, join only when yielding
    
    # Read file in smaller chunks to avoid memory overflow
    buffer = ""
    while True:
        chunk = file_object.read(chunk_size)
        if not chunk:
            break
        
        # Decode if bytes
        if isinstance(chunk, bytes):
            chunk = chunk.decode('utf-8', errors='ignore')
        
        buffer += chunk
        lines = buffer.split('\n')
        
        # Keep the last incomplete line in buffer
        buffer = lines[-1]
        lines = lines[:-1]
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                # Yield previous sequence
                if current_name and current_seq_parts:
                    # Join and yield, then clear parts to free memory immediately
                    full_seq = ''.join(current_seq_parts)
                    current_seq_parts.clear()  # Free memory immediately
                    yield (current_name, full_seq)
                    del full_seq  # Help garbage collector
                
                # Start new sequence
                current_name = line[1:].strip()
                if not current_name:
                    current_name = f"sequence_{id(line)}"
                current_seq_parts = []
            else:
                # Add to current sequence
                current_seq_parts.append(line.upper())
    
    # Process remaining buffer
    if buffer.strip() and not buffer.startswith('>'):
        current_seq_parts.append(buffer.strip().upper())
    
    # Yield last sequence
    if current_name and current_seq_parts:
        full_seq = ''.join(current_seq_parts)
        current_seq_parts.clear()
        yield (current_name, full_seq)
        del full_seq


def open_compressed_file(file_path_or_object):
    """
    Open a file that may be compressed (gzip, bgzip) or uncompressed.
    
    Automatically detects compression by file extension or magic bytes.
    Supports both file paths and file-like objects.
    
    Args:
        file_path_or_object: File path string or file-like object
        
    Returns:
        File-like object ready for reading (text mode)
        
    Example:
        >>> with open_compressed_file('genome.fa.gz') as f:
        ...     sequences = parse_fasta_chunked(f)
    """
    import gzip
    
    # Check if it's a file path or file object
    if isinstance(file_path_or_object, str):
        # File path - check extension
        if file_path_or_object.endswith('.gz') or file_path_or_object.endswith('.bgz'):
            return gzip.open(file_path_or_object, 'rt', encoding='utf-8')
        else:
            return open(file_path_or_object, 'r', encoding='utf-8')
    else:
        # File object - check magic bytes
        file_object = file_path_or_object
        file_object.seek(0)
        
        # Read first 2 bytes to check for gzip magic number (1f 8b)
        magic = file_object.read(2)
        file_object.seek(0)
        
        if isinstance(magic, bytes) and magic == b'\x1f\x8b':
            # Gzip compressed
            return gzip.open(file_object, 'rt', encoding='utf-8')
        else:
            # Not compressed, return as-is
            return file_object


def parse_fasta_chunked_compressed(file_path_or_object, chunk_size_mb: int = 2):
    """
    Parse FASTA file with automatic compression detection.
    
    Combines parse_fasta_chunked with open_compressed_file for seamless
    handling of both compressed (.gz, .bgz) and uncompressed FASTA files.
    
    Args:
        file_path_or_object: File path string or file-like object
        chunk_size_mb: Size of read chunks in MB (default: 2MB)
        
    Yields:
        Tuple of (sequence_name, sequence_string)
        
    Example:
        >>> # Works with compressed files
        >>> for name, seq in parse_fasta_chunked_compressed('large_genome.fa.gz'):
        ...     print(f"{name}: {len(seq)} bp")
        >>> 
        >>> # Also works with regular files
        >>> for name, seq in parse_fasta_chunked_compressed('sequences.fasta'):
        ...     print(f"{name}: {len(seq)} bp")
    """
    with open_compressed_file(file_path_or_object) as f:
        yield from parse_fasta_chunked(f, chunk_size_mb=chunk_size_mb)


def get_file_preview(file_object, max_sequences: int = 3, max_preview_chars: int = 400):
    """
    Get a preview of FASTA file content without loading entire file.
    
    Args:
        file_object: File-like object from st.file_uploader
        max_sequences: Maximum number of sequences to preview (default: 3)
        max_preview_chars: Maximum characters to show per sequence (default: 400)
        
    Returns:
        Dict with keys: 'num_sequences', 'total_bp', 'previews' (list of dicts)
    """
    file_object.seek(0)  # Reset to beginning
    
    previews = []
    total_bp = 0
    num_sequences = 0
    
    for name, seq in parse_fasta_chunked(file_object):
        num_sequences += 1
        seq_len = len(seq)
        total_bp += seq_len
        
        if len(previews) < max_sequences:
            preview_seq = seq[:max_preview_chars]
            if len(seq) > max_preview_chars:
                preview_seq += "..."
            
            previews.append({
                'name': name,
                'length': seq_len,
                'preview': preview_seq
            })
    
    file_object.seek(0)  # Reset for subsequent use
    
    return {
        'num_sequences': num_sequences,
        'total_bp': total_bp,
        'previews': previews
    }


def wrap(sequence: str, width: int = 80) -> str:
    """
    Wrap sequence to specified width
    
    Args:
        sequence: DNA sequence string
        width: Line width for wrapping
        
    Returns:
        Wrapped sequence string
    """
    if not sequence or width <= 0:
        return sequence
    
    return '\n'.join(sequence[i:i+width] for i in range(0, len(sequence), width))
