"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Chunk-Based Analyzer - Memory-Efficient Genome Analysis                      │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Implements chunk-based genome analysis with deduplication at boundaries.
    Analyzes genomes in 50Kbp chunks with 2Kbp overlap to handle motifs that
    span chunk boundaries without duplication.

ARCHITECTURE:
    - ChunkAnalyzer: Orchestrates chunk-based analysis workflow
    - Deduplication: Removes duplicate motifs at chunk boundaries
    - Progress tracking: Callbacks for UI updates
    - Memory management: Aggressive garbage collection between chunks

MEMORY GUARANTEES:
    - Processes one chunk at a time
    - Maximum memory: ~50-70MB for processing (chunk + results buffer)
    - Aggressive GC between chunks
    - Supports file uploads up to 100MB
"""

import gc
import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, Any, List, Optional, Callable, Set, Tuple
from collections import defaultdict

logger = logging.getLogger(__name__)


class ChunkAnalyzer:
    """
    Memory-efficient chunk-based analyzer for large genomes.
    
    Analyzes sequences in overlapping chunks to handle motifs at boundaries.
    Automatically deduplicates motifs found in overlap regions.
    
    Features:
        - Configurable chunk size (default: 50Kbp)
        - Configurable overlap (default: 2Kbp, optimized for performance)
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
        chunk_size: int = 50_000,       # Recommended: 50Kbp chunks for consistency
        overlap: int = 2_000,           # 2Kbp overlap (optimized)
        use_parallel: bool = True,
        max_workers: Optional[int] = None,
        use_adaptive: bool = False
    ):
        """
        Initialize chunk analyzer.
        
        Args:
            sequence_storage: UniversalSequenceStorage instance
            chunk_size: Size of each chunk in base pairs (default: 50Kbp, recommended for optimal performance)
            overlap: Overlap between chunks in base pairs (default: 2Kbp, optimized for performance)
            use_parallel: Enable parallel processing of chunks (default: True)
            max_workers: Maximum number of parallel workers (default: CPU count - 1)
            use_adaptive: Enable triple adaptive chunking strategy (default: False)
        """
        self.sequence_storage = sequence_storage
        self.chunk_size = chunk_size
        self.overlap = overlap
        self.use_parallel = use_parallel
        self.max_workers = max_workers or max(1, multiprocessing.cpu_count() - 1)
        self.use_adaptive = use_adaptive
        
        logger.info(f"ChunkAnalyzer initialized (chunk_size={chunk_size:,}, overlap={overlap:,}, "
                   f"parallel={use_parallel}, workers={self.max_workers}, adaptive={use_adaptive})")
        
        # Log multiprocessing configuration for diagnostics
        try:
            logger.info(f"Multiprocessing start method: {multiprocessing.get_start_method()}")
            logger.info(f"Available CPU cores: {multiprocessing.cpu_count()}")
        except (NotImplementedError, RuntimeError) as e:
            logger.warning(f"Multiprocessing diagnostics failed: {e}")
    
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
    
    @staticmethod
    def _process_chunk_worker(args: Tuple[str, str, int, int, Optional[List[str]]]) -> Tuple[int, List[Dict[str, Any]]]:
        """
        Static worker method for parallel chunk processing.
        
        Args:
            args: Tuple of (chunk_seq, seq_name_chunk, chunk_start, chunk_end, enabled_classes)
            
        Returns:
            Tuple of (chunk_start, list of motifs with adjusted positions)
        """
        from Utilities.nonbscanner import analyze_sequence
        
        chunk_seq, seq_name_chunk, chunk_start, chunk_end, enabled_classes = args
        
        # Analyze chunk
        chunk_motifs = analyze_sequence(
            sequence=chunk_seq,
            sequence_name=seq_name_chunk,
            use_fast_mode=True,
            enabled_classes=enabled_classes
        )
        
        # Adjust positions to global coordinates
        adjusted_motifs = []
        for motif in chunk_motifs:
            adjusted_motif = motif.copy()
            adjusted_motif['Start'] = motif['Start'] + chunk_start
            adjusted_motif['End'] = motif['End'] + chunk_start
            
            # Update ID to reflect global position
            if 'ID' in adjusted_motif:
                parts = adjusted_motif['ID'].split('_')
                if len(parts) >= 3:
                    parts[-1] = str(adjusted_motif['Start'])
                    adjusted_motif['ID'] = '_'.join(parts)
            
            adjusted_motifs.append(adjusted_motif)
        
        return (chunk_start, adjusted_motifs)
    
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
        # Use triple adaptive chunking if enabled
        if self.use_adaptive:
            from Utilities.triple_chunk_analyzer import TripleAdaptiveChunkAnalyzer
            
            logger.info("Using triple adaptive chunking strategy")
            adaptive_analyzer = TripleAdaptiveChunkAnalyzer(
                sequence_storage=self.sequence_storage,
                max_workers=self.max_workers,
                use_adaptive=True
            )
            return adaptive_analyzer.analyze(seq_id, progress_callback, enabled_classes)
        
        # Original single-tier chunking logic
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
        
        logger.info(f"Processing {total_chunks} chunks (parallel={self.use_parallel}, "
                   f"chunk_size={self.chunk_size:,}, overlap={self.overlap:,})")
        
        # Choose processing mode based on configuration
        _run_parallel = self.use_parallel and total_chunks > 1

        if _run_parallel:
            # Parallel processing mode
            logger.info(f"Using parallel processing with {self.max_workers} workers")
            logger.info(f"Processing {total_chunks} chunks in parallel batches")
            
            # Collect all chunk data
            chunk_data = []
            for chunk_seq, chunk_start, chunk_end in self.sequence_storage.iter_chunks(
                seq_id, self.chunk_size, self.overlap
            ):
                chunk_num = len(chunk_data) + 1
                logger.debug(f"Prepared chunk {chunk_num} [{chunk_start:,}-{chunk_end:,}] ({len(chunk_seq):,} bp)")
                chunk_data.append((
                    chunk_seq,
                    f"{seq_name}_chunk{chunk_num}",
                    chunk_start,
                    chunk_end,
                    enabled_classes
                ))
            
            # Process chunks in parallel, with fallback to sequential on any failure
            completed_chunks = 0
            all_chunk_results = {}
            
            try:
                with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                    # Submit all chunks
                    future_to_chunk = {
                        executor.submit(self._process_chunk_worker, data): idx
                        for idx, data in enumerate(chunk_data)
                    }
                    
                    # Collect results as they complete
                    for future in as_completed(future_to_chunk):
                        chunk_idx = future_to_chunk[future]
                        try:
                            chunk_start, adjusted_motifs = future.result()
                            all_chunk_results[chunk_idx] = (chunk_start, adjusted_motifs)
                            completed_chunks += 1
                            
                            logger.info(f"Chunk {completed_chunks}/{total_chunks} completed: "
                                       f"detected {len(adjusted_motifs)} motifs at position {chunk_start:,}")
                            
                            # Update progress
                            if progress_callback:
                                progress_pct = (completed_chunks / total_chunks) * 100
                                progress_callback(progress_pct)
                                logger.info(f"Overall progress: {progress_pct:.1f}% ({completed_chunks}/{total_chunks} chunks)")
                        except Exception as e:
                            logger.error(f"Error processing chunk {chunk_idx}: {e}")
                            chunk_tuple = chunk_data[chunk_idx] if chunk_idx < len(chunk_data) else None
                            chunk_pos = chunk_tuple[2] if chunk_tuple and len(chunk_tuple) > 2 else "unknown"
                            raise RuntimeError(
                                f"Failed to process chunk {chunk_idx + 1}/{total_chunks} "
                                f"at position {chunk_pos}: {e}"
                            ) from e
                
                # Process results in order and deduplicate
                logger.info(f"Deduplicating motifs across {len(all_chunk_results)} chunk boundaries")
                for chunk_idx in sorted(all_chunk_results.keys()):
                    chunk_start, adjusted_motifs = all_chunk_results[chunk_idx]
                    
                    # Validate chunk_data access with bounds checking
                    if chunk_idx >= len(chunk_data):
                        logger.error(f"Missing chunk_data for chunk {chunk_idx}")
                        continue
                    
                    chunk_tuple = chunk_data[chunk_idx]
                    if len(chunk_tuple) < 4:
                        logger.error(f"Invalid chunk_data tuple for chunk {chunk_idx}: "
                                    f"expected 4+ elements, got {len(chunk_tuple)}")
                        continue
                    
                    # Get chunk_end from original data (index 3 in the tuple)
                    chunk_end = chunk_tuple[3]
                    
                    # Filter out duplicates from overlap region
                    unique_motifs = []
                    duplicates_removed = 0
                    for motif in adjusted_motifs:
                        motif_key = self._create_motif_key(motif)
                        
                        if motif_key not in overlap_motifs:
                            unique_motifs.append(motif)
                            
                            # If motif is in overlap region, track it for next chunk
                            if self._is_in_overlap_region(motif, chunk_start, chunk_end, self.overlap):
                                overlap_motifs.add(motif_key)
                        else:
                            duplicates_removed += 1
                    
                    # Save results to disk
                    results_storage.append_batch(unique_motifs)
                    
                    logger.info(f"Chunk {chunk_idx + 1}: {len(adjusted_motifs)} motifs detected, "
                               f"{len(unique_motifs)} unique after deduplication "
                               f"({duplicates_removed} duplicates removed)")
                    
                    # Free memory for processed chunk results
                    del all_chunk_results[chunk_idx]
                    del adjusted_motifs
                    del unique_motifs
                    gc.collect()
                
                # Final cleanup
                del all_chunk_results
                del chunk_data
                gc.collect()

            except Exception as e:
                logger.warning(
                    f"Parallel chunk analysis failed ({type(e).__name__}: {e}), "
                    f"falling back to sequential processing"
                )
                # Reset state for clean sequential fallback.
                # all_chunk_results and overlap_motifs may be partially populated;
                # reset them so the sequential path starts fresh.
                # chunk_data is freed here; sequential path uses iter_chunks() directly.
                _run_parallel = False
                all_chunk_results = {}
                overlap_motifs = set()
                del chunk_data  # sequential path uses iter_chunks(), not chunk_data
                gc.collect()

        if not _run_parallel:
            # Sequential processing mode (also used as fallback when parallel fails)
            chunk_num = 0
            for chunk_seq, chunk_start, chunk_end in self.sequence_storage.iter_chunks(
                seq_id, self.chunk_size, self.overlap
            ):
                chunk_num += 1
                
                logger.info(f"Processing chunk {chunk_num}/{total_chunks} "
                           f"[{chunk_start:,}-{chunk_end:,}] ({len(chunk_seq):,} bp)")
                
                # Analyze chunk
                chunk_motifs = analyze_sequence(
                    sequence=chunk_seq,
                    sequence_name=f"{seq_name}_chunk{chunk_num}",
                    use_fast_mode=True,
                    enabled_classes=enabled_classes
                )
                
                logger.info(f"Chunk {chunk_num}: detected {len(chunk_motifs)} raw motifs")
                
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
                           f"{len(unique_motifs)} unique after deduplication")
                
                # Update progress (including percentage)
                if progress_callback:
                    progress_pct = (chunk_num / total_chunks) * 100
                    progress_callback(progress_pct)
                    logger.info(f"Overall progress: {progress_pct:.1f}% ({chunk_num}/{total_chunks} chunks)")
                
                # Aggressive garbage collection
                del chunk_seq
                del chunk_motifs
                del adjusted_motifs
                del unique_motifs
                gc.collect()
        
        # Final statistics
        stats = results_storage.get_summary_stats()
        logger.info(f"Analysis complete: {stats['total_count']} total motifs detected")
        
        # Final memory cleanup
        del overlap_motifs
        gc.collect()
        
        return results_storage
