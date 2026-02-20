"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Universal Disk-Based Storage System - Memory-Constant Genome Analysis        │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Implements disk-based storage for sequences and results to maintain constant
    memory usage. Enables analysis of genomes up to 100MB in size on Streamlit
    Community Cloud free tier (1GB RAM limit) with disk-based storage.

ARCHITECTURE:
    - UniversalSequenceStorage: Saves sequences to disk, provides chunk-based iteration
    - UniversalResultsStorage: Streams results to disk in JSONL format, supports pagination
    
MEMORY GUARANTEES:
    - Sequence storage: Never loads full sequence into memory
    - Results storage: Never loads all results into memory
    - Chunk iteration: 50Kbp chunks with 2Kbp overlap (balanced for performance/accuracy)
    - Memory overhead: ~50-70MB total including active chunk and results buffering
    - Summary stats: Computed incrementally without full load
"""

import os
import json
import tempfile
import shutil
import logging
from typing import Dict, Any, List, Optional, Iterator, Tuple
from pathlib import Path
import hashlib
from datetime import datetime
from Utilities.detectors_utils import calc_gc_content, _count_bases

logger = logging.getLogger(__name__)


class UniversalSequenceStorage:
    """
    Disk-based sequence storage with chunk-based iteration.
    
    Saves all genome sequences to disk immediately upon upload to avoid
    memory bottlenecks. Provides chunk-based iteration for analysis with
    configurable chunk size and overlap.
    
    Features:
        - Automatic temp directory management
        - Chunk-based iteration (default: 50Kbp chunks, 2Kbp overlap for performance/accuracy balance)
        - Metadata caching (length, GC%, etc.)
        - Memory-safe operations (never loads full sequence)
        
    Usage:
        storage = UniversalSequenceStorage()
        seq_id = storage.save_sequence(sequence, "chr1")
        
        # Iterate in chunks
        for chunk_data in storage.iter_chunks(seq_id, chunk_size=50_000):
            seq_chunk, start, end = chunk_data
            # Process chunk...
            
        # Get metadata without loading sequence
        metadata = storage.get_metadata(seq_id)
        print(f"Length: {metadata['length']}, GC%: {metadata['gc_content']:.2f}")
        
        # Cleanup
        storage.cleanup()
    """
    
    def __init__(self, base_dir: Optional[str] = None):
        """
        Initialize storage system.
        
        Args:
            base_dir: Optional base directory for storage. If None, uses system temp.
        """
        if base_dir is None:
            self.base_dir = Path(tempfile.mkdtemp(prefix="nonbdna_sequences_"))
        else:
            self.base_dir = Path(base_dir)
            self.base_dir.mkdir(parents=True, exist_ok=True)
        
        self.metadata_file = self.base_dir / "metadata.json"
        self.metadata: Dict[str, Dict[str, Any]] = {}
        
        # Load existing metadata if present
        if self.metadata_file.exists():
            with open(self.metadata_file, 'r') as f:
                self.metadata = json.load(f)
        
        logger.info(f"UniversalSequenceStorage initialized at {self.base_dir}")
    
    def _generate_seq_id(self, name: str) -> str:
        """Generate unique sequence ID from name and timestamp."""
        timestamp = datetime.now().isoformat()
        hash_input = f"{name}_{timestamp}".encode('utf-8')
        return hashlib.md5(hash_input).hexdigest()[:12]
    
    def _calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content percentage using standardized method."""
        # Use the same method as detectors_utils.calc_gc_content for consistency
        return calc_gc_content(sequence)
    
    def save_sequence(self, sequence: str, name: str) -> str:
        """
        Save sequence to disk and return unique ID.
        
        Sanitizes the sequence by removing all whitespace and converting to uppercase
        to ensure consistent byte-offset based chunking.
        
        Args:
            sequence: DNA sequence string (may contain whitespace/newlines)
            name: Sequence identifier (e.g., "chr1", "genome")
            
        Returns:
            seq_id: Unique identifier for retrieving this sequence
        """
        seq_id = self._generate_seq_id(name)
        seq_file = self.base_dir / f"{seq_id}.seq"
        
        # Sanitize sequence: remove all whitespace (newlines, spaces, tabs) and convert to uppercase
        # This ensures byte-offset based chunking works correctly
        original_length = len(sequence)
        sanitized_sequence = ''.join(sequence.split()).upper()
        sanitized_length = len(sanitized_sequence)
        
        # Log if sanitization removed characters
        if sanitized_length != original_length:
            logger.info(f"Sanitized sequence '{name}': removed {original_length - sanitized_length} "
                       f"whitespace characters ({original_length:,} -> {sanitized_length:,} bp)")
        
        # Write sanitized sequence to disk (no newlines, no whitespace)
        with open(seq_file, 'w') as f:
            f.write(sanitized_sequence)
        
        # Calculate and store metadata using sanitized sequence
        a_count, t_count, g_count, c_count = _count_bases(sanitized_sequence)
        n_count = sanitized_sequence.count('N')
        metadata = {
            'seq_id': seq_id,
            'name': name,
            'length': sanitized_length,  # Length of actual stored sequence
            'gc_content': self._calculate_gc_content(sanitized_sequence),
            'base_counts': {
                'A': a_count,
                'T': t_count,
                'G': g_count,
                'C': c_count,
                'N': n_count,
            },
            'file_path': str(seq_file),
            'created_at': datetime.now().isoformat()
        }
        
        self.metadata[seq_id] = metadata
        self._save_metadata()
        
        logger.info(f"Saved sequence '{name}' (length: {sanitized_length:,} bp) as {seq_id}")
        return seq_id
    
    def _save_metadata(self):
        """Persist metadata to disk."""
        with open(self.metadata_file, 'w') as f:
            json.dump(self.metadata, f, indent=2)
    
    def get_sequence_chunk(self, seq_id: str, start: int, end: int) -> str:
        """
        Retrieve specific sequence chunk by position.
        
        Args:
            seq_id: Sequence identifier
            start: Start position (0-based)
            end: End position (exclusive)
            
        Returns:
            Sequence substring from start to end
        """
        if seq_id not in self.metadata:
            raise KeyError(f"Sequence ID '{seq_id}' not found")
        
        seq_file = Path(self.metadata[seq_id]['file_path'])
        
        with open(seq_file, 'r') as f:
            f.seek(start)
            return f.read(end - start)
    
    def iter_chunks(
        self,
        seq_id: str,
        chunk_size: int = 50_000,       # Recommended: 50Kbp chunks
        overlap: int = 2_000            # 2Kbp overlap (optimized)
    ) -> Iterator[Tuple[str, int, int]]:
        """
        Iterate over sequence in overlapping chunks.
        
        Validates that chunks contain only expected base characters (no whitespace/newlines).
        If whitespace is detected, raises an error instructing to re-save the sequence.
        
        Args:
            seq_id: Sequence identifier
            chunk_size: Size of each chunk in base pairs (default: 50Kbp, recommended for optimal performance)
            overlap: Overlap between chunks in base pairs (default: 2Kbp, balanced for performance and accuracy)
            
        Yields:
            Tuple of (chunk_sequence, chunk_start, chunk_end)
            
        Raises:
            ValueError: If chunk contains unexpected whitespace/newlines
            
        Example:
            for seq_chunk, start, end in storage.iter_chunks(seq_id):
                motifs = analyze_sequence(seq_chunk)
                # Process motifs at positions [start:end]
        """
        if seq_id not in self.metadata:
            raise KeyError(f"Sequence ID '{seq_id}' not found")
        
        seq_file = Path(self.metadata[seq_id]['file_path'])
        seq_length = self.metadata[seq_id]['length']
        
        chunk_num = 0
        start = 0
        
        with open(seq_file, 'r') as f:
            while start < seq_length:
                end = min(start + chunk_size, seq_length)
                
                f.seek(start)
                chunk = f.read(end - start)
                
                chunk_num += 1
                
                # Defensive validation: check for unexpected whitespace/newlines
                # This catches cases where sequences were stored with FASTA formatting
                if any(c.isspace() for c in chunk):
                    logger.error(f"Chunk {chunk_num} [{start}-{end}] contains whitespace/newlines!")
                    logger.error(f"First 100 chars: {repr(chunk[:100])}")
                    raise ValueError(
                        f"Sequence file '{seq_file}' contains whitespace/newlines in chunk {chunk_num} "
                        f"[{start}-{end}]. This indicates the sequence was saved with FASTA formatting. "
                        f"Please re-upload the sequence to allow automatic sanitization, or use the "
                        f"updated save_sequence() method which now strips whitespace automatically."
                    )
                
                logger.debug(f"Yielding chunk {chunk_num} [{start}-{end}], length={len(chunk)}")
                yield (chunk, start, end)
                
                # If we've reached the end, stop
                if end >= seq_length:
                    break
                
                # Move to next chunk with overlap
                start = end - overlap
        
        logger.info(f"Completed iteration of {chunk_num} chunks for sequence {seq_id}")
    
    def get_metadata(self, seq_id: str) -> Dict[str, Any]:
        """
        Get sequence metadata without loading sequence.
        
        Args:
            seq_id: Sequence identifier
            
        Returns:
            Dictionary with metadata (name, length, gc_content, etc.)
        """
        if seq_id not in self.metadata:
            raise KeyError(f"Sequence ID '{seq_id}' not found")
        
        return self.metadata[seq_id].copy()
    
    def list_sequences(self) -> List[Dict[str, Any]]:
        """
        List all stored sequences with metadata.
        
        Returns:
            List of metadata dictionaries for all sequences
        """
        return list(self.metadata.values())
    
    def cleanup(self, seq_id: Optional[str] = None):
        """
        Clean up storage files.
        
        Args:
            seq_id: If provided, delete only this sequence. If None, delete all.
        """
        if seq_id is not None:
            # Delete specific sequence
            if seq_id in self.metadata:
                seq_file = Path(self.metadata[seq_id]['file_path'])
                if seq_file.exists():
                    seq_file.unlink()
                del self.metadata[seq_id]
                self._save_metadata()
                logger.info(f"Deleted sequence {seq_id}")
        else:
            # Delete entire storage directory
            if self.base_dir.exists():
                shutil.rmtree(self.base_dir)
                logger.info(f"Deleted storage directory {self.base_dir}")


class UniversalResultsStorage:
    """
    Disk-based results storage with streaming and pagination.
    
    Streams analysis results to disk in JSONL format (one JSON per line).
    Never loads all results into memory. Supports iterator-based access,
    pagination, and cached summary statistics.
    
    Features:
        - Streaming append (no memory bottleneck)
        - Iterator-based access (lazy loading)
        - Cached summary statistics
        - Pagination support for display
        - JSONL format (one motif per line)
        
    Usage:
        results = UniversalResultsStorage("results_dir", "seq1")
        
        # Append results as they're generated
        for motif in detected_motifs:
            results.append(motif)
        
        # Get summary without loading all results
        stats = results.get_summary_stats()
        print(f"Total motifs: {stats['total_count']}")
        
        # Iterate for display (lazy loading)
        for motif in results.iter_results(limit=100):
            print(motif['Class'], motif['Start'])
        
        # Convert to DataFrame for visualization
        df = results.to_dataframe(limit=1000)
    """
    
    def __init__(self, base_dir: str, seq_id: str):
        """
        Initialize results storage for a sequence.
        
        Args:
            base_dir: Base directory for results storage
            seq_id: Sequence identifier
        """
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)
        
        self.seq_id = seq_id
        self.results_file = self.base_dir / f"{seq_id}_results.jsonl"
        self.stats_file = self.base_dir / f"{seq_id}_stats.json"
        
        # Cached statistics
        self._stats: Optional[Dict[str, Any]] = None
        self._motif_count = 0
        
        # Load existing stats if present
        if self.stats_file.exists():
            with open(self.stats_file, 'r') as f:
                self._stats = json.load(f)
                self._motif_count = self._stats.get('total_count', 0)
        
        logger.info(f"UniversalResultsStorage initialized for {seq_id}")
    
    def append(self, motif: Dict[str, Any]):
        """
        Append single motif to results file.
        
        Args:
            motif: Motif dictionary with fields (Class, Start, End, Score, etc.)
        """
        # Write motif as JSONL (one line per motif)
        with open(self.results_file, 'a') as f:
            f.write(json.dumps(motif) + '\n')
        
        self._motif_count += 1
        
        # Invalidate cached stats
        self._stats = None
    
    def append_batch(self, motifs: List[Dict[str, Any]]):
        """
        Append multiple motifs efficiently.
        
        Args:
            motifs: List of motif dictionaries
        """
        with open(self.results_file, 'a') as f:
            for motif in motifs:
                f.write(json.dumps(motif) + '\n')
        
        self._motif_count += len(motifs)
        self._stats = None
    
    def iter_results(self, limit: Optional[int] = None) -> Iterator[Dict[str, Any]]:
        """
        Iterate over results (lazy loading).
        
        Args:
            limit: Maximum number of results to return (None = all)
            
        Yields:
            Motif dictionaries
        """
        if not self.results_file.exists():
            return
        
        count = 0
        with open(self.results_file, 'r') as f:
            for line in f:
                if limit is not None and count >= limit:
                    break
                
                motif = json.loads(line.strip())
                yield motif
                count += 1
    
    def get_summary_stats(self) -> Dict[str, Any]:
        """
        Get summary statistics (computed without full load).
        
        Returns:
            Dictionary with statistics:
                - total_count: Total number of motifs
                - class_distribution: Count per class
                - subclass_distribution: Count per subclass
                - coverage_bp: Total bases covered
                - avg_score: Average motif score
                - score_range: (min, max) score
        """
        if self._stats is not None:
            return self._stats
        
        # Compute statistics by streaming through file
        stats = {
            'total_count': 0,
            'class_distribution': {},
            'subclass_distribution': {},
            'coverage_bp': 0,
            'score_sum': 0.0,
            'min_score': float('inf'),
            'max_score': float('-inf'),
            'positions': []
        }
        
        if not self.results_file.exists():
            self._stats = stats
            return stats
        
        with open(self.results_file, 'r') as f:
            for line in f:
                motif = json.loads(line.strip())
                
                stats['total_count'] += 1
                
                # Class distribution
                motif_class = motif.get('Class', 'Unknown')
                stats['class_distribution'][motif_class] = \
                    stats['class_distribution'].get(motif_class, 0) + 1
                
                # Subclass distribution
                motif_subclass = motif.get('Subclass', 'Unknown')
                stats['subclass_distribution'][motif_subclass] = \
                    stats['subclass_distribution'].get(motif_subclass, 0) + 1
                
                # Coverage
                length = motif.get('Length', motif.get('End', 0) - motif.get('Start', 0))
                stats['coverage_bp'] += length
                
                # Score statistics
                score = motif.get('Score', 0.0)
                stats['score_sum'] += score
                stats['min_score'] = min(stats['min_score'], score)
                stats['max_score'] = max(stats['max_score'], score)
        
        # Calculate averages
        if stats['total_count'] > 0:
            stats['avg_score'] = stats['score_sum'] / stats['total_count']
            stats['score_range'] = (stats['min_score'], stats['max_score'])
        else:
            stats['avg_score'] = 0.0
            stats['score_range'] = (0.0, 0.0)
        
        # Remove temporary fields
        del stats['score_sum']
        del stats['min_score']
        del stats['max_score']
        
        # Cache and persist
        self._stats = stats
        with open(self.stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        
        return stats
    
    def to_dataframe(self, limit: Optional[int] = None):
        """
        Convert results to pandas DataFrame.
        
        Args:
            limit: Maximum number of results to load (None = all)
            
        Returns:
            pandas DataFrame with motif data
        """
        import pandas as pd
        
        motifs = list(self.iter_results(limit=limit))
        return pd.DataFrame(motifs)
    
    def get_results_file_path(self) -> str:
        """Get path to JSONL results file (for direct export)."""
        return str(self.results_file)
    
    def cleanup(self):
        """Delete results files."""
        if self.results_file.exists():
            self.results_file.unlink()
        if self.stats_file.exists():
            self.stats_file.unlink()
        logger.info(f"Deleted results for {self.seq_id}")
