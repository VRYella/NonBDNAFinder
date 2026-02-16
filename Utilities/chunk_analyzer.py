"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Chunk-Based Analyzer - Memory-Efficient Genome Analysis                      │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Implements chunk-based genome analysis with deduplication at boundaries.
    Analyzes genomes in 5MB chunks with 10KB overlap to handle motifs that
    span chunk boundaries without duplication.

ARCHITECTURE:
    - ChunkAnalyzer: Orchestrates chunk-based analysis workflow
    - Deduplication: Removes duplicate motifs at chunk boundaries
    - Progress tracking: Callbacks for UI updates
    - Memory management: Aggressive garbage collection between chunks

MEMORY GUARANTEES:
    - Processes one chunk at a time
    - Maximum memory: chunk_size + results buffer (~50MB)
    - Aggressive GC between chunks
"""

import gc
import logging
from typing import Dict, Any, List, Optional, Callable, Set, Tuple
from collections import defaultdict

logger = logging.getLogger(__name__)


class ChunkAnalyzer:
    """
    Memory-efficient chunk-based analyzer for large genomes.
    
    Analyzes sequences in overlapping chunks to handle motifs at boundaries.
    Automatically deduplicates motifs found in overlap regions.
    
    Features:
        - Configurable chunk size (default: 5MB)
        - Configurable overlap (default: 10KB)
        - Automatic boundary deduplication
        - Progress callbacks for UI integration
        - Aggressive garbage collection
        
    Usage:
        from Utilities.disk_storage import UniversalSequenceStorage, UniversalResultsStorage
        from Utilities.nonbscanner import analyze_sequence
        
        storage = UniversalSequenceStorage()
        seq_id = storage.save_sequence(sequence, "chr1")
        
        analyzer = ChunkAnalyzer(storage)
        results_storage = analyzer.analyze(
            seq_id=seq_id,
            progress_callback=lambda p: print(f"Progress: {p}%")
        )
        
        stats = results_storage.get_summary_stats()
        print(f"Found {stats['total_count']} motifs")
    """
    
    def __init__(
        self,
        sequence_storage,
        chunk_size: int = 5_000_000,
        overlap: int = 10_000
    ):
        """
        Initialize chunk analyzer.
        
        Args:
            sequence_storage: UniversalSequenceStorage instance
            chunk_size: Size of each chunk in base pairs (default: 5MB)
            overlap: Overlap between chunks in base pairs (default: 10KB)
        """
        self.sequence_storage = sequence_storage
        self.chunk_size = chunk_size
        self.overlap = overlap
        
        logger.info(f"ChunkAnalyzer initialized (chunk_size={chunk_size:,}, overlap={overlap:,})")
    
    def _create_motif_key(self, motif: Dict[str, Any]) -> Tuple:
        """
        Create unique key for motif deduplication.
        
        Args:
            motif: Motif dictionary
            
        Returns:
            Tuple representing motif identity (Class, Subclass, Start, End)
        """
        return (
            motif.get('Class', ''),
            motif.get('Subclass', ''),
            motif.get('Start', 0),
            motif.get('End', 0)
        )
    
    def _is_in_overlap_region(
        self,
        motif: Dict[str, Any],
        chunk_start: int,
        chunk_end: int,
        overlap: int
    ) -> bool:
        """
        Check if motif is in the overlap region between chunks.
        
        Args:
            motif: Motif dictionary
            chunk_start: Start position of current chunk
            chunk_end: End position of current chunk
            overlap: Overlap size
            
        Returns:
            True if motif overlaps with next chunk's overlap region
        """
        motif_start = motif.get('Start', 0)
        motif_end = motif.get('End', 0)
        
        # Check if motif is in the last 'overlap' bp of the chunk
        overlap_region_start = chunk_end - overlap
        
        # Motif overlaps if it starts or ends in the overlap region
        return (motif_start >= overlap_region_start or 
                motif_end > overlap_region_start)
    
    def _adjust_motif_positions(
        self,
        motifs: List[Dict[str, Any]],
        chunk_start: int
    ) -> List[Dict[str, Any]]:
        """
        Adjust motif positions from chunk-local to genome-global coordinates.
        
        Args:
            motifs: List of motifs with chunk-local positions
            chunk_start: Start position of chunk in genome
            
        Returns:
            List of motifs with adjusted positions
        """
        adjusted_motifs = []
        
        for motif in motifs:
            adjusted_motif = motif.copy()
            adjusted_motif['Start'] = motif['Start'] + chunk_start
            adjusted_motif['End'] = motif['End'] + chunk_start
            
            # Update ID to reflect global position
            if 'ID' in adjusted_motif:
                parts = adjusted_motif['ID'].split('_')
                if len(parts) >= 3:
                    # Format: seqname_CLASS_position
                    parts[-1] = str(adjusted_motif['Start'])
                    adjusted_motif['ID'] = '_'.join(parts)
            
            adjusted_motifs.append(adjusted_motif)
        
        return adjusted_motifs
    
    def _deduplicate_motifs(
        self,
        all_motifs: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Remove duplicate motifs (same position, class, subclass).
        
        Args:
            all_motifs: List of all motifs (may contain duplicates)
            
        Returns:
            Deduplicated list of motifs
        """
        seen_keys: Set[Tuple] = set()
        unique_motifs = []
        
        for motif in all_motifs:
            key = self._create_motif_key(motif)
            
            if key not in seen_keys:
                seen_keys.add(key)
                unique_motifs.append(motif)
        
        return unique_motifs
    
    def analyze(
        self,
        seq_id: str,
        progress_callback: Optional[Callable[[float], None]] = None,
        enabled_classes: Optional[List[str]] = None
    ):
        """
        Analyze sequence in chunks with deduplication.
        
        Args:
            seq_id: Sequence identifier from UniversalSequenceStorage
            progress_callback: Optional callback function(progress_pct: float)
            enabled_classes: Optional list of motif classes to analyze
            
        Returns:
            UniversalResultsStorage instance with results
        """
        from Utilities.disk_storage import UniversalResultsStorage
        from Utilities.nonbscanner import analyze_sequence
        
        # Get sequence metadata
        metadata = self.sequence_storage.get_metadata(seq_id)
        seq_name = metadata['name']
        seq_length = metadata['length']
        
        logger.info(f"Starting chunk analysis for {seq_name} (length: {seq_length:,} bp)")
        
        # Initialize results storage
        results_storage = UniversalResultsStorage(
            base_dir=str(self.sequence_storage.base_dir / "results"),
            seq_id=seq_id
        )
        
        # Track motifs in overlap regions for deduplication
        overlap_motifs: Set[Tuple] = set()
        
        # Calculate total chunks for progress
        total_chunks = 0
        for _ in self.sequence_storage.iter_chunks(seq_id, self.chunk_size, self.overlap):
            total_chunks += 1
        
        logger.info(f"Processing {total_chunks} chunks")
        
        # Process each chunk
        chunk_num = 0
        for chunk_seq, chunk_start, chunk_end in self.sequence_storage.iter_chunks(
            seq_id, self.chunk_size, self.overlap
        ):
            chunk_num += 1
            
            logger.info(f"Processing chunk {chunk_num}/{total_chunks} "
                       f"[{chunk_start:,}-{chunk_end:,}]")
            
            # Analyze chunk
            chunk_motifs = analyze_sequence(
                sequence=chunk_seq,
                sequence_name=f"{seq_name}_chunk{chunk_num}",
                use_fast_mode=True,
                enabled_classes=enabled_classes
            )
            
            # Adjust positions to global coordinates
            adjusted_motifs = self._adjust_motif_positions(chunk_motifs, chunk_start)
            
            # Filter out duplicates from overlap region
            unique_motifs = []
            for motif in adjusted_motifs:
                motif_key = self._create_motif_key(motif)
                
                # Skip if we've already seen this motif in a previous chunk
                if motif_key in overlap_motifs:
                    continue
                
                unique_motifs.append(motif)
                
                # If motif is in overlap region, track it for next chunk
                if self._is_in_overlap_region(motif, chunk_start, chunk_end, self.overlap):
                    overlap_motifs.add(motif_key)
            
            # Save results to disk
            results_storage.append_batch(unique_motifs)
            
            logger.info(f"Chunk {chunk_num}: {len(chunk_motifs)} motifs detected, "
                       f"{len(unique_motifs)} unique")
            
            # Update progress
            if progress_callback:
                progress_pct = (chunk_num / total_chunks) * 100
                progress_callback(progress_pct)
            
            # Aggressive garbage collection
            del chunk_seq
            del chunk_motifs
            del adjusted_motifs
            del unique_motifs
            gc.collect()
        
        # Final statistics
        stats = results_storage.get_summary_stats()
        logger.info(f"Analysis complete: {stats['total_count']} total motifs detected")
        
        return results_storage
