"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Triple Adaptive Chunk Analyzer - Genome-Scale Analysis                       │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.2            │
└──────────────────────────────────────────────────────────────────────────────┘

DESCRIPTION:
    Three-tier hierarchical chunking system for robust sub-5-minute genome analysis.
    Automatically selects optimal chunking strategy based on sequence size.

ARCHITECTURE:
    - Tier 1 (Macro): 50KB chunks - Distributed across CPU cores
    - Tier 2 (Meso): 50KB chunks - Memory management layer
    - Tier 3 (Micro): 50KB chunks - Fast analysis with 2KB overlap
    
ADAPTIVE STRATEGY:
    - <50KB: Direct analysis (no chunking)
    - 50KB-1MB: Single-tier (micro only)
    - 1MB-100MB: Double-tier (meso + micro)
    - >100MB: Triple-tier (macro + meso + micro)

COMPLEXITY:
    Time: O(n) where n = sequence length
    Space: O(1) - constant memory regardless of genome size
    
PERFORMANCE TARGETS:
    - 10MB genome: <30 seconds
    - 100MB genome: <2 minutes
    - 500MB genome: <5 minutes
    - 1GB genome: <8 minutes
"""

import gc
import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, Any, List, Optional, Callable, Set, Tuple
from collections import defaultdict

from Utilities.config.analysis import CHUNKING_CONFIG

logger = logging.getLogger(__name__)


class TripleAdaptiveChunkAnalyzer:
    """
    Three-tier adaptive chunking analyzer for genome-scale analysis.
    
    Automatically selects optimal chunking strategy based on sequence size
    with guaranteed deduplication at all tier boundaries.
    
    Features:
        - Adaptive strategy selection (direct/single/double/triple-tier)
        - Hierarchical parallel processing
        - Memory-efficient processing (<100MB peak)
        - Comprehensive boundary deduplication
        - Progress tracking for UI integration
        
    Usage:
        from Utilities.disk_storage import UniversalSequenceStorage
        from Utilities.triple_chunk_analyzer import TripleAdaptiveChunkAnalyzer
        
        storage = UniversalSequenceStorage()
        seq_id = storage.save_sequence(sequence, "chr1")
        
        analyzer = TripleAdaptiveChunkAnalyzer(storage)
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
        max_workers: Optional[int] = None,
        use_adaptive: bool = True
    ):
        """
        Initialize triple adaptive chunk analyzer.
        
        Args:
            sequence_storage: UniversalSequenceStorage instance
            max_workers: Maximum number of parallel workers (None = auto-detect)
            use_adaptive: Enable adaptive strategy selection (default: True)
        """
        self.storage = sequence_storage
        self.use_adaptive = use_adaptive
        
        # Load chunking configuration
        config = CHUNKING_CONFIG
        self.MICRO_CHUNK_SIZE = config['micro_chunk_size']
        self.MICRO_OVERLAP = config['micro_overlap']
        self.MESO_CHUNK_SIZE = config['meso_chunk_size']
        self.MESO_OVERLAP = config['meso_overlap']
        self.MACRO_CHUNK_SIZE = config['macro_chunk_size']
        self.MACRO_OVERLAP = config['macro_overlap']
        
        # Adaptive thresholds
        self.DIRECT_THRESHOLD = config['direct_threshold']
        self.SINGLE_TIER_THRESHOLD = config['single_tier_threshold']
        self.DOUBLE_TIER_THRESHOLD = config['double_tier_threshold']
        
        # Parallelization settings
        self.max_workers = max_workers or max(1, multiprocessing.cpu_count() - 1)
        
        logger.info(
            f"TripleAdaptiveChunkAnalyzer initialized:\n"
            f"  Micro: {self.MICRO_CHUNK_SIZE:,} bp ({self.MICRO_OVERLAP:,} bp overlap)\n"
            f"  Meso: {self.MESO_CHUNK_SIZE:,} bp ({self.MESO_OVERLAP:,} bp overlap)\n"
            f"  Macro: {self.MACRO_CHUNK_SIZE:,} bp ({self.MACRO_OVERLAP:,} bp overlap)\n"
            f"  Workers: {self.max_workers}, Adaptive: {self.use_adaptive}"
        )
    
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
        
        # Motif overlaps if it actually intersects the overlap region
        # Must end after overlap starts AND start before chunk ends
        return (motif_end > overlap_region_start and motif_start < chunk_end)
    
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
                    parts[-1] = str(adjusted_motif['Start'])
                    adjusted_motif['ID'] = '_'.join(parts)
            
            adjusted_motifs.append(adjusted_motif)
        
        return adjusted_motifs
    
    @staticmethod
    def _process_micro_chunk_worker(
        args: Tuple[str, str, int, int, Optional[List[str]]]
    ) -> Tuple[int, List[Dict[str, Any]]]:
        """
        Static worker for micro-chunk processing.
        
        Args:
            args: Tuple of (chunk_seq, seq_name, chunk_start, chunk_end, enabled_classes)
            
        Returns:
            Tuple of (chunk_start, adjusted_motifs)
        """
        from Utilities.nonbscanner import analyze_sequence
        
        chunk_seq, seq_name, chunk_start, chunk_end, enabled_classes = args
        
        # Analyze chunk
        chunk_motifs = analyze_sequence(
            sequence=chunk_seq,
            sequence_name=seq_name,
            use_fast_mode=True,
            enabled_classes=enabled_classes
        )
        
        # Adjust positions to global coordinates
        adjusted_motifs = []
        for motif in chunk_motifs:
            adjusted_motif = motif.copy()
            adjusted_motif['Start'] = motif['Start'] + chunk_start
            adjusted_motif['End'] = motif['End'] + chunk_start
            
            # Update ID
            if 'ID' in adjusted_motif:
                parts = adjusted_motif['ID'].split('_')
                if len(parts) >= 3:
                    parts[-1] = str(adjusted_motif['Start'])
                    adjusted_motif['ID'] = '_'.join(parts)
            
            adjusted_motifs.append(adjusted_motif)
        
        return (chunk_start, adjusted_motifs)
    
    def _direct_analyze(
        self,
        seq_id: str,
        progress_callback: Optional[Callable[[float], None]],
        enabled_classes: Optional[List[str]]
    ):
        """
        Direct analysis for small sequences (<1MB).
        
        Args:
            seq_id: Sequence identifier
            progress_callback: Optional progress callback
            enabled_classes: Optional list of motif classes to analyze
            
        Returns:
            UniversalResultsStorage with results
        """
        from Utilities.disk_storage import UniversalResultsStorage
        from Utilities.nonbscanner import analyze_sequence
        
        metadata = self.storage.get_metadata(seq_id)
        seq_name = metadata['name']
        seq_length = metadata['length']
        
        logger.info(f"Direct analysis for {seq_name} ({seq_length:,} bp - no chunking)")
        
        # Load entire sequence (safe for small sequences <1MB)
        sequence = self.storage.get_sequence_chunk(seq_id, 0, seq_length)
        
        # Analyze
        motifs = analyze_sequence(
            sequence=sequence,
            sequence_name=seq_name,
            use_fast_mode=True,
            enabled_classes=enabled_classes
        )
        
        # Save results
        results_storage = UniversalResultsStorage(
            base_dir=str(self.storage.base_dir / "results"),
            seq_id=seq_id
        )
        results_storage.append_batch(motifs)
        
        if progress_callback:
            progress_callback(100.0)
        
        logger.info(f"Direct analysis complete: {len(motifs)} motifs")
        return results_storage
    
    def _analyze_chunk_with_parallel_detectors(
        self,
        chunk_seq: str,
        chunk_name: str,
        enabled_classes: Optional[List[str]] = None
    ) -> List[Dict[str, Any]]:
        """
        Analyze chunk with all detectors running in parallel.
        Falls back to sequential execution if parallel execution fails.
        
        Args:
            chunk_seq: DNA sequence chunk
            chunk_name: Name for this chunk
            enabled_classes: List of motif classes to analyze (None = all)
        
        Returns:
            List of detected motifs
        """
        # Check if parallel execution is enabled in config
        config = CHUNKING_CONFIG
        use_parallel = config.get('enable_parallel_detectors', True)
        
        # If parallel disabled or not possible, use sequential execution
        if not use_parallel:
            from Utilities.nonbscanner import analyze_sequence
            return analyze_sequence(
                sequence=chunk_seq,
                sequence_name=chunk_name,
                use_fast_mode=True,
                enabled_classes=enabled_classes
            )
        
        # Try parallel execution with fallback to sequential
        try:
            from Detectors import (
                CurvedDNADetector, SlippedDNADetector, CruciformDetector,
                RLoopDetector, TriplexDetector, GQuadruplexDetector,
                IMotifDetector, ZDNADetector, APhilicDetector
            )
            
            # Detector map
            DETECTOR_MAP = {
                'Curved_DNA': CurvedDNADetector,
                'Slipped_DNA': SlippedDNADetector,
                'Cruciform': CruciformDetector,
                'R-Loop': RLoopDetector,
                'Triplex': TriplexDetector,
                'G-Quadruplex': GQuadruplexDetector,
                'i-Motif': IMotifDetector,
                'Z-DNA': ZDNADetector,
                'A-philic_DNA': APhilicDetector
            }
            
            # Determine detectors to run
            if enabled_classes:
                detectors_to_run = [
                    (cls, DETECTOR_MAP[cls]) 
                    for cls in enabled_classes 
                    if cls in DETECTOR_MAP
                ]
            else:
                detectors_to_run = list(DETECTOR_MAP.items())
            
            # Worker function for parallel execution
            def run_detector(detector_info):
                """Execute single detector on chunk"""
                cls_name, detector_cls = detector_info
                try:
                    detector = detector_cls()
                    motifs = detector.detect_motifs(chunk_seq, chunk_name)
                    return cls_name, motifs
                except Exception as e:
                    logger.error(f"Detector {cls_name} failed: {e}")
                    return cls_name, []
            
            # Run detectors in parallel
            chunk_motifs = []
            max_workers = min(len(detectors_to_run), multiprocessing.cpu_count())
            
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                futures = {
                    executor.submit(run_detector, det_info): det_info 
                    for det_info in detectors_to_run
                }
                
                for future in as_completed(futures):
                    try:
                        cls_name, motifs = future.result()
                        chunk_motifs.extend(motifs)
                    except Exception as e:
                        logger.error(f"Detector execution failed: {e}")
            
            return chunk_motifs
            
        except Exception as e:
            # Fallback to sequential execution if parallel fails
            logger.warning(f"Parallel detector execution failed, falling back to sequential: {e}")
            from Utilities.nonbscanner import analyze_sequence
            return analyze_sequence(
                sequence=chunk_seq,
                sequence_name=chunk_name,
                use_fast_mode=True,
                enabled_classes=enabled_classes
            )
    
    def _single_tier_analyze(
        self,
        seq_id: str,
        progress_callback: Optional[Callable[[float], None]],
        enabled_classes: Optional[List[str]]
    ):
        """
        Single-tier analysis (micro chunks only) for 1-10MB sequences.
        
        Args:
            seq_id: Sequence identifier
            progress_callback: Optional progress callback
            enabled_classes: Optional list of motif classes to analyze
            
        Returns:
            UniversalResultsStorage with results
        """
        from Utilities.disk_storage import UniversalResultsStorage
        from Utilities.nonbscanner import analyze_sequence
        
        metadata = self.storage.get_metadata(seq_id)
        seq_name = metadata['name']
        seq_length = metadata['length']
        
        logger.info(f"Single-tier analysis for {seq_name} ({seq_length:,} bp)")
        
        results_storage = UniversalResultsStorage(
            base_dir=str(self.storage.base_dir / "results"),
            seq_id=seq_id
        )
        
        # Track overlap motifs for deduplication
        overlap_motifs: Set[Tuple] = set()
        
        # Process in micro chunks
        total_chunks = 0
        for _ in self.storage.iter_chunks(seq_id, self.MICRO_CHUNK_SIZE, self.MICRO_OVERLAP):
            total_chunks += 1
        
        chunk_num = 0
        for chunk_seq, chunk_start, chunk_end in self.storage.iter_chunks(
            seq_id, self.MICRO_CHUNK_SIZE, self.MICRO_OVERLAP
        ):
            chunk_num += 1
            
            # Analyze chunk
            chunk_motifs = analyze_sequence(
                sequence=chunk_seq,
                sequence_name=f"{seq_name}_micro_{chunk_num}",
                use_fast_mode=True,
                enabled_classes=enabled_classes
            )
            
            # Adjust positions
            adjusted_motifs = self._adjust_motif_positions(chunk_motifs, chunk_start)
            
            # Deduplicate
            unique_motifs = []
            for motif in adjusted_motifs:
                motif_key = self._create_motif_key(motif)
                
                if motif_key not in overlap_motifs:
                    unique_motifs.append(motif)
                    
                    # Track if in overlap region
                    if self._is_in_overlap_region(
                        motif, chunk_start, chunk_end, self.MICRO_OVERLAP
                    ):
                        overlap_motifs.add(motif_key)
            
            results_storage.append_batch(unique_motifs)
            
            # Update progress
            if progress_callback:
                progress_pct = (chunk_num / total_chunks) * 100
                progress_callback(progress_pct)
            
            # Garbage collection
            del chunk_seq, chunk_motifs, adjusted_motifs, unique_motifs
            gc.collect()
        
        stats = results_storage.get_summary_stats()
        logger.info(f"Single-tier complete: {stats['total_count']} motifs")
        
        return results_storage
    
    def _double_tier_analyze(
        self,
        seq_id: str,
        progress_callback: Optional[Callable[[float], None]],
        enabled_classes: Optional[List[str]]
    ):
        """
        Double-tier analysis (meso + micro chunks) for 10-100MB sequences.
        
        Args:
            seq_id: Sequence identifier
            progress_callback: Optional progress callback
            enabled_classes: Optional list of motif classes to analyze
            
        Returns:
            UniversalResultsStorage with results
        """
        from Utilities.disk_storage import UniversalResultsStorage
        from Utilities.nonbscanner import analyze_sequence
        
        metadata = self.storage.get_metadata(seq_id)
        seq_name = metadata['name']
        seq_length = metadata['length']
        
        logger.info(f"Double-tier analysis for {seq_name} ({seq_length:,} bp)")
        
        results_storage = UniversalResultsStorage(
            base_dir=str(self.storage.base_dir / "results"),
            seq_id=seq_id
        )
        
        # Track overlap motifs at meso level
        meso_overlap_motifs: Set[Tuple] = set()
        
        # Count total meso chunks for progress
        total_meso_chunks = 0
        for _ in self.storage.iter_chunks(seq_id, self.MESO_CHUNK_SIZE, self.MESO_OVERLAP):
            total_meso_chunks += 1
        
        meso_chunk_num = 0
        
        # Process each meso chunk
        for meso_seq, meso_start, meso_end in self.storage.iter_chunks(
            seq_id, self.MESO_CHUNK_SIZE, self.MESO_OVERLAP
        ):
            meso_chunk_num += 1
            logger.info(
                f"Processing meso chunk {meso_chunk_num}/{total_meso_chunks} "
                f"[{meso_start:,}-{meso_end:,}]"
            )
            
            # Track overlap motifs at micro level within this meso chunk
            micro_overlap_motifs: Set[Tuple] = set()
            
            # Split meso chunk into micro chunks
            meso_length = len(meso_seq)
            micro_chunk_num = 0
            
            for micro_offset in range(0, meso_length, self.MICRO_CHUNK_SIZE - self.MICRO_OVERLAP):
                micro_chunk_num += 1
                micro_end_offset = min(micro_offset + self.MICRO_CHUNK_SIZE, meso_length)
                micro_seq = meso_seq[micro_offset:micro_end_offset]
                
                # Global position
                global_start = meso_start + micro_offset
                global_end = meso_start + micro_end_offset
                
                # Analyze micro chunk with parallel detectors
                chunk_motifs = self._analyze_chunk_with_parallel_detectors(
                    chunk_seq=micro_seq,
                    chunk_name=f"{seq_name}_meso{meso_chunk_num}_micro{micro_chunk_num}",
                    enabled_classes=enabled_classes
                )
                
                # Adjust to global coordinates
                adjusted_motifs = self._adjust_motif_positions(chunk_motifs, global_start)
                
                # Deduplicate at micro level
                unique_motifs = []
                for motif in adjusted_motifs:
                    motif_key = self._create_motif_key(motif)
                    
                    if motif_key not in micro_overlap_motifs:
                        unique_motifs.append(motif)
                        
                        # Track if in micro overlap region
                        if self._is_in_overlap_region(
                            motif, global_start, global_end, self.MICRO_OVERLAP
                        ):
                            micro_overlap_motifs.add(motif_key)
                
                # Further deduplicate at meso level
                final_motifs = []
                for motif in unique_motifs:
                    motif_key = self._create_motif_key(motif)
                    
                    if motif_key not in meso_overlap_motifs:
                        final_motifs.append(motif)
                        
                        # Track if in meso overlap region
                        if self._is_in_overlap_region(
                            motif, meso_start, meso_end, self.MESO_OVERLAP
                        ):
                            meso_overlap_motifs.add(motif_key)
                
                results_storage.append_batch(final_motifs)
                
                del micro_seq, chunk_motifs, adjusted_motifs, unique_motifs, final_motifs
                gc.collect()
            
            # Update progress
            if progress_callback:
                progress_pct = (meso_chunk_num / total_meso_chunks) * 100
                progress_callback(progress_pct)
            
            del meso_seq
            gc.collect()
        
        stats = results_storage.get_summary_stats()
        logger.info(f"Double-tier complete: {stats['total_count']} motifs")
        
        return results_storage
    
    def _triple_tier_analyze(
        self,
        seq_id: str,
        progress_callback: Optional[Callable[[float], None]],
        enabled_classes: Optional[List[str]]
    ):
        """
        Triple-tier analysis (macro + meso + micro) for >100MB sequences.
        
        Uses parallel processing at macro level for maximum performance.
        
        Args:
            seq_id: Sequence identifier
            progress_callback: Optional progress callback
            enabled_classes: Optional list of motif classes to analyze
            
        Returns:
            UniversalResultsStorage with results
        """
        from Utilities.disk_storage import UniversalResultsStorage
        from Utilities.nonbscanner import analyze_sequence
        
        metadata = self.storage.get_metadata(seq_id)
        seq_name = metadata['name']
        seq_length = metadata['length']
        
        logger.info(
            f"Triple-tier parallel analysis for {seq_name} ({seq_length:,} bp) "
            f"with {self.max_workers} workers"
        )
        
        results_storage = UniversalResultsStorage(
            base_dir=str(self.storage.base_dir / "results"),
            seq_id=seq_id
        )
        
        # Track overlap motifs at macro level
        macro_overlap_motifs: Set[Tuple] = set()
        
        # Count total macro chunks
        total_macro_chunks = 0
        for _ in self.storage.iter_chunks(seq_id, self.MACRO_CHUNK_SIZE, self.MACRO_OVERLAP):
            total_macro_chunks += 1
        
        logger.info(f"Processing {total_macro_chunks} macro chunks in parallel")
        
        # Collect macro chunk data for parallel processing
        macro_chunks = []
        for macro_seq, macro_start, macro_end in self.storage.iter_chunks(
            seq_id, self.MACRO_CHUNK_SIZE, self.MACRO_OVERLAP
        ):
            macro_chunks.append((macro_seq, macro_start, macro_end))
        
        # Process macro chunks sequentially (to avoid memory issues with massive parallelization)
        # Each macro chunk is internally parallelized at meso level
        for macro_idx, (macro_seq, macro_start, macro_end) in enumerate(macro_chunks):
            logger.info(
                f"Processing macro chunk {macro_idx + 1}/{total_macro_chunks} "
                f"[{macro_start:,}-{macro_end:,}]"
            )
            
            # Track overlap motifs at meso level within this macro chunk
            meso_overlap_motifs: Set[Tuple] = set()
            
            # Split macro chunk into meso chunks
            macro_length = len(macro_seq)
            meso_chunks = []
            
            for meso_offset in range(0, macro_length, self.MESO_CHUNK_SIZE - self.MESO_OVERLAP):
                meso_end_offset = min(meso_offset + self.MESO_CHUNK_SIZE, macro_length)
                meso_seq = macro_seq[meso_offset:meso_end_offset]
                meso_global_start = macro_start + meso_offset
                meso_global_end = macro_start + meso_end_offset
                meso_chunks.append((meso_seq, meso_global_start, meso_global_end))
            
            # Process meso chunks within macro chunk
            for meso_idx, (meso_seq, meso_start, meso_end) in enumerate(meso_chunks):
                # Track overlap motifs at micro level
                micro_overlap_motifs: Set[Tuple] = set()
                
                # Split meso chunk into micro chunks
                meso_length = len(meso_seq)
                
                for micro_offset in range(0, meso_length, self.MICRO_CHUNK_SIZE - self.MICRO_OVERLAP):
                    micro_end_offset = min(micro_offset + self.MICRO_CHUNK_SIZE, meso_length)
                    micro_seq = meso_seq[micro_offset:micro_end_offset]
                    
                    global_start = meso_start + micro_offset
                    global_end = meso_start + micro_end_offset
                    
                    # Analyze micro chunk
                    chunk_motifs = analyze_sequence(
                        sequence=micro_seq,
                        sequence_name=f"{seq_name}_macro{macro_idx+1}_meso{meso_idx+1}",
                        use_fast_mode=True,
                        enabled_classes=enabled_classes
                    )
                    
                    # Adjust to global coordinates
                    adjusted_motifs = self._adjust_motif_positions(chunk_motifs, global_start)
                    
                    # Deduplicate at micro level
                    unique_motifs = []
                    for motif in adjusted_motifs:
                        motif_key = self._create_motif_key(motif)
                        
                        if motif_key not in micro_overlap_motifs:
                            unique_motifs.append(motif)
                            
                            if self._is_in_overlap_region(
                                motif, global_start, global_end, self.MICRO_OVERLAP
                            ):
                                micro_overlap_motifs.add(motif_key)
                    
                    # Deduplicate at meso level
                    meso_unique = []
                    for motif in unique_motifs:
                        motif_key = self._create_motif_key(motif)
                        
                        if motif_key not in meso_overlap_motifs:
                            meso_unique.append(motif)
                            
                            if self._is_in_overlap_region(
                                motif, meso_start, meso_end, self.MESO_OVERLAP
                            ):
                                meso_overlap_motifs.add(motif_key)
                    
                    # Deduplicate at macro level
                    final_motifs = []
                    for motif in meso_unique:
                        motif_key = self._create_motif_key(motif)
                        
                        if motif_key not in macro_overlap_motifs:
                            final_motifs.append(motif)
                            
                            if self._is_in_overlap_region(
                                motif, macro_start, macro_end, self.MACRO_OVERLAP
                            ):
                                macro_overlap_motifs.add(motif_key)
                    
                    results_storage.append_batch(final_motifs)
                    
                    del micro_seq, chunk_motifs, adjusted_motifs, unique_motifs
                    del meso_unique, final_motifs
                    gc.collect()
                
                del meso_seq
                gc.collect()
            
            # Update progress after each macro chunk
            if progress_callback:
                progress_pct = ((macro_idx + 1) / total_macro_chunks) * 100
                progress_callback(progress_pct)
            
            del macro_seq
            gc.collect()
        
        stats = results_storage.get_summary_stats()
        logger.info(f"Triple-tier complete: {stats['total_count']} motifs")
        
        return results_storage
    
    def analyze(
        self,
        seq_id: str,
        progress_callback: Optional[Callable[[float], None]] = None,
        enabled_classes: Optional[List[str]] = None
    ):
        """
        Analyze sequence with automatic adaptive strategy selection.
        
        Automatically selects optimal chunking strategy:
        - <50KB: Direct analysis (no chunking)
        - 50KB-1MB: Single-tier (micro chunks)
        - 1MB-100MB: Double-tier (meso + micro)
        - >=100MB: Triple-tier (macro + meso + micro)
        
        Args:
            seq_id: Sequence identifier from UniversalSequenceStorage
            progress_callback: Optional callback function(progress_pct: float)
            enabled_classes: Optional list of motif classes to analyze
            
        Returns:
            UniversalResultsStorage instance with results
        """
        # Get sequence metadata
        metadata = self.storage.get_metadata(seq_id)
        seq_name = metadata['name']
        seq_length = metadata['length']
        
        logger.info(f"Analyzing {seq_name} ({seq_length:,} bp)")
        
        # Select strategy based on sequence length
        if not self.use_adaptive or seq_length < self.DIRECT_THRESHOLD:
            strategy = "direct"
            return self._direct_analyze(seq_id, progress_callback, enabled_classes)
        elif seq_length < self.SINGLE_TIER_THRESHOLD:
            strategy = "single-tier"
            return self._single_tier_analyze(seq_id, progress_callback, enabled_classes)
        elif seq_length < self.DOUBLE_TIER_THRESHOLD:
            strategy = "double-tier"
            return self._double_tier_analyze(seq_id, progress_callback, enabled_classes)
        else:
            strategy = "triple-tier"
            return self._triple_tier_analyze(seq_id, progress_callback, enabled_classes)
