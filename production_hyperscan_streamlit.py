#!/usr/bin/env python3
"""
Production-ready Hyperscan-Streamlit integration for heavy-duty genomic sequence scanning.

This module provides a complete solution for scanning large genomic sequences with Hyperscan
in Streamlit applications without blocking the UI while maintaining sane resource usage.

Key Features:
- Efficient DB compilation and caching (compile once per process)
- Chunking/streaming with proper match carryover handling
- ProcessPoolExecutor with worker initialization for parallel processing
- Progress tracking and incremental output
- Memory-efficient streaming for large genomes
- Copy-paste ready code snippets

Architecture follows these principles:
1. Compile DB once per process (expensive) - cached with st.cache_resource
2. Minimal work in Hyperscan callbacks - only append (pattern_id, start, end)
3. Parallel processing across independent units (contigs/chunks)
4. Chunk with overlap to avoid missing cross-boundary matches
5. One hs.Scratch(db) per thread/process (not thread-shareable)
6. Stream results incrementally - don't hold millions of hits in RAM
"""

import os
import sys
import time
import tempfile
import logging
from typing import List, Dict, Tuple, Optional, Generator, Union, Any
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from io import StringIO
import multiprocessing

try:
    import streamlit as st
    STREAMLIT_AVAILABLE = True
except ImportError:
    STREAMLIT_AVAILABLE = False

import pandas as pd
import numpy as np

try:
    import hyperscan as hs
    HYPERSCAN_AVAILABLE = True
except ImportError:
    HYPERSCAN_AVAILABLE = False

from Bio import SeqIO


logger = logging.getLogger(__name__)


# Configuration constants
DEFAULT_CHUNK_SIZE = 2_000_000  # 2MB chunks
DEFAULT_OVERLAP = 500  # Safe margin for most biological patterns
DEFAULT_MAX_WORKERS = min(4, multiprocessing.cpu_count())  # Conservative for Streamlit
DEFAULT_TIMEOUT = 300  # 5 minutes default timeout


@dataclass
class ScanResult:
    """Container for scan results with metadata"""
    pattern_id: int
    start: int
    end: int
    pattern: str
    sequence_name: str = ""
    score: float = 0.0
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'pattern_id': self.pattern_id,
            'start': self.start,
            'end': self.end,
            'length': self.end - self.start,
            'pattern': self.pattern,
            'sequence_name': self.sequence_name,
            'score': self.score
        }


@dataclass
class ScanConfig:
    """Configuration for scanning operations"""
    chunk_size: int = DEFAULT_CHUNK_SIZE
    overlap: int = DEFAULT_OVERLAP
    max_workers: int = DEFAULT_MAX_WORKERS
    timeout: int = DEFAULT_TIMEOUT
    stream_results: bool = True
    output_format: str = 'csv'  # csv, parquet, json
    incremental_save: bool = True
    max_memory_mb: int = 1000  # Limit memory usage


class HyperscanStreamlitScanner:
    """
    Production-ready Hyperscan scanner optimized for Streamlit applications.
    
    This class implements all the best practices for high-performance genomic scanning:
    - Efficient DB compilation and caching
    - Worker process initialization
    - Chunking with overlap handling
    - Progress tracking
    - Memory management
    """
    
    def __init__(self, patterns: List[str], config: Optional[ScanConfig] = None):
        """
        Initialize scanner with patterns and configuration.
        
        Args:
            patterns: List of regex patterns for Hyperscan
            config: Scanning configuration
        """
        if not HYPERSCAN_AVAILABLE:
            raise RuntimeError("Hyperscan is not available. Install with: pip install hyperscan")
        
        self.patterns = patterns
        self.config = config or ScanConfig()
        self._db = None
        self._pattern_map = {i: pattern for i, pattern in enumerate(patterns)}
        
        # Compile database immediately
        self._compile_database()
    
    def _compile_database(self) -> None:
        """Compile Hyperscan database with patterns"""
        logger.info(f"Compiling Hyperscan database with {len(self.patterns)} patterns")
        
        if not self.patterns:
            raise ValueError("No patterns provided for compilation")
        
        # Convert patterns to bytes and prepare compilation
        exprs = []
        ids = []
        flags = []
        
        for i, pattern in enumerate(self.patterns):
            expr = pattern.encode('utf-8') if isinstance(pattern, str) else pattern
            exprs.append(expr)
            ids.append(i)
            flags.append(hs.HS_FLAG_UTF8 | hs.HS_FLAG_CASELESS | hs.HS_FLAG_DOTALL)
        
        try:
            self._db = hs.Database()
            self._db.compile(expressions=exprs, ids=ids, flags=flags)
            logger.info("Successfully compiled Hyperscan database")
        except Exception as e:
            logger.error(f"Failed to compile Hyperscan database: {e}")
            raise RuntimeError(f"Database compilation failed: {e}")
    
    def scan_sequence(
        self, 
        sequence: str, 
        sequence_name: str = "sequence",
        progress_callback: Optional[callable] = None
    ) -> List[ScanResult]:
        """
        Scan a single sequence with progress tracking.
        
        Args:
            sequence: DNA sequence to scan
            sequence_name: Name/ID of the sequence
            progress_callback: Optional callback for progress updates
            
        Returns:
            List of scan results
        """
        if not sequence:
            return []
        
        # For smaller sequences, scan directly
        if len(sequence) <= self.config.chunk_size:
            return self._scan_chunk(sequence, sequence_name, 0)
        
        # For larger sequences, use chunking
        return self._scan_chunked_sequence(sequence, sequence_name, progress_callback)
    
    def _scan_chunk(
        self, 
        chunk: str, 
        sequence_name: str, 
        offset: int
    ) -> List[ScanResult]:
        """Scan a single chunk of sequence"""
        if not self._db:
            raise RuntimeError("Database not compiled")
        
        scratch = hs.Scratch(self._db)
        matches = []
        
        def on_match(pat_id, start, end, flags, ctx):
            """Minimal callback - only append match info"""
            matches.append((int(pat_id), int(start) + offset, int(end) + offset))
        
        try:
            chunk_bytes = chunk.encode('utf-8')
            self._db.scan(chunk_bytes, match_event_handler=on_match, scratch=scratch)
        except Exception as e:
            logger.error(f"Scanning failed for chunk: {e}")
            raise
        
        # Convert to ScanResult objects
        results = []
        for pat_id, start, end in matches:
            pattern = self._pattern_map.get(pat_id, f"pattern_{pat_id}")
            result = ScanResult(
                pattern_id=pat_id,
                start=start,
                end=end,
                pattern=pattern,
                sequence_name=sequence_name
            )
            results.append(result)
        
        return results
    
    def _scan_chunked_sequence(
        self, 
        sequence: str, 
        sequence_name: str,
        progress_callback: Optional[callable] = None
    ) -> List[ScanResult]:
        """Scan large sequence using chunking with overlap"""
        chunks = self._create_chunks(sequence)
        all_results = []
        
        logger.info(f"Scanning sequence {sequence_name} with {len(chunks)} chunks")
        
        for i, (chunk, offset) in enumerate(chunks):
            if progress_callback:
                progress = (i + 1) / len(chunks)
                progress_callback(progress, f"Scanning chunk {i+1}/{len(chunks)}")
            
            chunk_results = self._scan_chunk(chunk, sequence_name, offset)
            all_results.extend(chunk_results)
        
        # Remove duplicate matches from overlapping regions
        return self._deduplicate_results(all_results)
    
    def _create_chunks(self, sequence: str) -> List[Tuple[str, int]]:
        """Create overlapping chunks from sequence"""
        chunks = []
        seq_len = len(sequence)
        chunk_size = self.config.chunk_size
        overlap = self.config.overlap
        
        offset = 0
        while offset < seq_len:
            end = min(offset + chunk_size, seq_len)
            chunk = sequence[offset:end]
            chunks.append((chunk, offset))
            
            if end >= seq_len:
                break
            
            # Move forward by chunk_size - overlap
            offset += chunk_size - overlap
        
        return chunks
    
    def _deduplicate_results(self, results: List[ScanResult]) -> List[ScanResult]:
        """Remove duplicate results from overlapping chunks"""
        if not results:
            return results
        
        # Sort by position
        results.sort(key=lambda r: (r.sequence_name, r.start, r.end))
        
        # Remove exact duplicates
        unique_results = []
        seen = set()
        
        for result in results:
            key = (result.sequence_name, result.start, result.end, result.pattern_id)
            if key not in seen:
                unique_results.append(result)
                seen.add(key)
        
        logger.debug(f"Deduplicated {len(results)} -> {len(unique_results)} results")
        return unique_results
    
    def scan_fasta_file(
        self, 
        fasta_path: str,
        output_prefix: Optional[str] = None,
        progress_callback: Optional[callable] = None
    ) -> Dict[str, Any]:
        """
        Scan FASTA file with parallel processing and incremental output.
        
        Args:
            fasta_path: Path to FASTA file
            output_prefix: Prefix for output files
            progress_callback: Optional callback for progress updates
            
        Returns:
            Dictionary with output file paths and statistics
        """
        # Load sequences
        sequences = []
        with open(fasta_path, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences.append((record.id, str(record.seq)))
        
        if not sequences:
            raise ValueError("No sequences found in FASTA file")
        
        logger.info(f"Loaded {len(sequences)} sequences from {fasta_path}")
        
        # Use parallel processing for multiple sequences
        if len(sequences) > 1 and self.config.max_workers > 1:
            return self._scan_parallel(sequences, output_prefix, progress_callback)
        else:
            return self._scan_sequential(sequences, output_prefix, progress_callback)
    
    def _scan_parallel(
        self, 
        sequences: List[Tuple[str, str]], 
        output_prefix: Optional[str],
        progress_callback: Optional[callable] = None
    ) -> Dict[str, Any]:
        """Scan sequences in parallel using ProcessPoolExecutor"""
        if output_prefix is None:
            output_prefix = f"/tmp/hyperscan_scan_{int(time.time())}"
        
        # Prepare work items
        work_items = []
        for seq_id, sequence in sequences:
            # For very large sequences, chunk them
            if len(sequence) > self.config.chunk_size * 2:
                chunks = self._create_chunks(sequence)
                for i, (chunk, offset) in enumerate(chunks):
                    work_items.append({
                        'sequence': chunk,
                        'sequence_name': f"{seq_id}_chunk_{i}",
                        'offset': offset,
                        'original_name': seq_id
                    })
            else:
                work_items.append({
                    'sequence': sequence,
                    'sequence_name': seq_id,
                    'offset': 0,
                    'original_name': seq_id
                })
        
        logger.info(f"Processing {len(work_items)} work items with {self.config.max_workers} workers")
        
        all_results = []
        completed = 0
        
        with ProcessPoolExecutor(
            max_workers=self.config.max_workers,
            initializer=_worker_init,
            initargs=(self.patterns,)
        ) as executor:
            # Submit all tasks
            future_to_item = {
                executor.submit(_worker_scan_chunk, item): item 
                for item in work_items
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_item, timeout=self.config.timeout):
                try:
                    results = future.result()
                    all_results.extend(results)
                    completed += 1
                    
                    if progress_callback:
                        progress = completed / len(work_items)
                        progress_callback(progress, f"Completed {completed}/{len(work_items)} chunks")
                    
                    # Incremental save to manage memory
                    if self.config.incremental_save and len(all_results) > 10000:
                        self._save_incremental_results(all_results, output_prefix)
                        all_results = []
                        
                except Exception as e:
                    logger.error(f"Worker task failed: {e}")
        
        # Final save
        if all_results:
            self._save_incremental_results(all_results, output_prefix)
        
        # Generate output files
        return self._generate_output_files(output_prefix)
    
    def _scan_sequential(
        self, 
        sequences: List[Tuple[str, str]], 
        output_prefix: Optional[str],
        progress_callback: Optional[callable] = None
    ) -> Dict[str, Any]:
        """Scan sequences sequentially"""
        if output_prefix is None:
            output_prefix = f"/tmp/hyperscan_scan_{int(time.time())}"
        
        all_results = []
        
        for i, (seq_id, sequence) in enumerate(sequences):
            if progress_callback:
                progress = (i + 1) / len(sequences)
                progress_callback(progress, f"Scanning sequence {i+1}/{len(sequences)}")
            
            results = self.scan_sequence(sequence, seq_id)
            all_results.extend(results)
            
            # Incremental save for large result sets
            if self.config.incremental_save and len(all_results) > 10000:
                self._save_incremental_results(all_results, output_prefix)
                all_results = []
        
        # Final save
        if all_results:
            self._save_incremental_results(all_results, output_prefix)
        
        return self._generate_output_files(output_prefix)
    
    def _save_incremental_results(self, results: List[ScanResult], output_prefix: str) -> None:
        """Save results incrementally to manage memory"""
        if not results:
            return
        
        # Convert to DataFrame
        df_data = [result.to_dict() for result in results]
        df = pd.DataFrame(df_data)
        
        # Append to CSV file
        csv_path = f"{output_prefix}_incremental.csv"
        df.to_csv(csv_path, mode='a', header=not os.path.exists(csv_path), index=False)
        
        logger.debug(f"Saved {len(results)} results to {csv_path}")
    
    def _generate_output_files(self, output_prefix: str) -> Dict[str, Any]:
        """Generate final output files from incremental saves"""
        csv_path = f"{output_prefix}_incremental.csv"
        
        if not os.path.exists(csv_path):
            logger.warning("No incremental results file found")
            return {'csv': None, 'stats': {'total_matches': 0}}
        
        # Read incremental results
        df = pd.read_csv(csv_path)
        
        # Generate different output formats
        output_files = {}
        
        # CSV (already exists)
        final_csv = f"{output_prefix}.csv"
        df.to_csv(final_csv, index=False)
        output_files['csv'] = final_csv
        
        # Parquet for better performance
        if self.config.output_format in ['parquet', 'all']:
            parquet_path = f"{output_prefix}.parquet"
            df.to_parquet(parquet_path, index=False)
            output_files['parquet'] = parquet_path
        
        # Excel for compatibility
        if self.config.output_format in ['excel', 'all']:
            excel_path = f"{output_prefix}.xlsx"
            df.to_excel(excel_path, index=False)
            output_files['excel'] = excel_path
        
        # Statistics
        stats = {
            'total_matches': len(df),
            'sequences_scanned': df['sequence_name'].nunique(),
            'patterns_matched': df['pattern_id'].nunique(),
            'total_length': df['length'].sum() if 'length' in df.columns else 0
        }
        
        # Cleanup incremental file
        try:
            os.unlink(csv_path)
        except Exception:
            pass
        
        return {'files': output_files, 'stats': stats}


# Worker functions for multiprocessing (must be at module level)
_worker_db = None
_worker_patterns = None


def _worker_init(patterns: List[str]) -> None:
    """Initialize worker process with compiled Hyperscan database"""
    global _worker_db, _worker_patterns
    
    try:
        import hyperscan as hs
        
        _worker_patterns = patterns
        
        # Compile database in worker process
        exprs = []
        ids = []
        flags = []
        
        for i, pattern in enumerate(patterns):
            expr = pattern.encode('utf-8') if isinstance(pattern, str) else pattern
            exprs.append(expr)
            ids.append(i)
            flags.append(hs.HS_FLAG_UTF8 | hs.HS_FLAG_CASELESS | hs.HS_FLAG_DOTALL)
        
        _worker_db = hs.Database()
        _worker_db.compile(expressions=exprs, ids=ids, flags=flags)
        
        logger.info(f"Worker initialized with {len(patterns)} patterns")
        
    except Exception as e:
        logger.error(f"Worker initialization failed: {e}")
        raise


def _worker_scan_chunk(work_item: Dict[str, Any]) -> List[ScanResult]:
    """Worker function to scan a chunk of sequence"""
    global _worker_db, _worker_patterns
    
    if _worker_db is None:
        raise RuntimeError("Worker not properly initialized")
    
    sequence = work_item['sequence']
    sequence_name = work_item['sequence_name']
    offset = work_item['offset']
    original_name = work_item['original_name']
    
    scratch = hs.Scratch(_worker_db)
    matches = []
    
    def on_match(pat_id, start, end, flags, ctx):
        """Minimal callback - only append match info"""
        matches.append((int(pat_id), int(start) + offset, int(end) + offset))
    
    try:
        sequence_bytes = sequence.encode('utf-8')
        _worker_db.scan(sequence_bytes, match_event_handler=on_match, scratch=scratch)
    except Exception as e:
        logger.error(f"Worker scanning failed: {e}")
        return []
    
    # Convert to ScanResult objects
    results = []
    pattern_map = {i: pattern for i, pattern in enumerate(_worker_patterns)}
    
    for pat_id, start, end in matches:
        pattern = pattern_map.get(pat_id, f"pattern_{pat_id}")
        result = ScanResult(
            pattern_id=pat_id,
            start=start,
            end=end,
            pattern=pattern,
            sequence_name=original_name  # Use original sequence name
        )
        results.append(result)
    
    return results


# Streamlit integration utilities
if STREAMLIT_AVAILABLE:
    
    @st.cache_resource
    def get_cached_scanner(patterns: List[str], config_dict: Dict[str, Any]) -> HyperscanStreamlitScanner:
        """Cache compiled Hyperscan scanner in Streamlit"""
        config = ScanConfig(**config_dict)
        return HyperscanStreamlitScanner(patterns, config)
    
    
    def streamlit_scan_with_progress(
        scanner: HyperscanStreamlitScanner,
        fasta_path: str,
        output_prefix: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Scan FASTA file with Streamlit progress tracking.
        
        Args:
            scanner: Initialized HyperscanStreamlitScanner
            fasta_path: Path to FASTA file
            output_prefix: Prefix for output files
            
        Returns:
            Dictionary with results and output files
        """
        # Create progress bars
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        def progress_callback(progress: float, message: str):
            progress_bar.progress(progress)
            status_text.text(message)
        
        # Scan with progress tracking
        results = scanner.scan_fasta_file(fasta_path, output_prefix, progress_callback)
        
        # Clear progress indicators
        progress_bar.empty()
        status_text.empty()
        
        return results


# Copy-paste ready code snippets
EXAMPLE_PATTERNS = [
    r'G{3,}N{1,7}G{3,}N{1,7}G{3,}N{1,7}G{3,}',  # G-quadruplex
    r'A{4,}T{4,}|T{4,}A{4,}',  # AT-rich
    r'[CG]{20,}',  # CpG island
    r'CA{6,}|TG{6,}',  # Microsatellite
]


def example_streamlit_app():
    """Example Streamlit app using the production-ready scanner"""
    if not STREAMLIT_AVAILABLE:
        print("Streamlit not available. Install with: pip install streamlit")
        return
    
    st.title("🧬 Production Hyperscan Genomic Scanner")
    st.markdown("High-performance pattern scanning for large genomic sequences")
    
    # Configuration
    st.sidebar.header("Configuration")
    chunk_size = st.sidebar.slider("Chunk Size (MB)", 1, 10, 2) * 1_000_000
    overlap = st.sidebar.slider("Overlap (bp)", 100, 1000, 500)
    max_workers = st.sidebar.slider("Max Workers", 1, 8, 4)
    
    # Pattern input
    st.header("Patterns to Search")
    patterns_text = st.text_area(
        "Enter patterns (one per line):",
        value="\n".join(EXAMPLE_PATTERNS),
        height=150
    )
    patterns = [p.strip() for p in patterns_text.split('\n') if p.strip()]
    
    # File upload
    st.header("Upload FASTA File")
    uploaded_file = st.file_uploader("Choose a FASTA file", type=['fa', 'fasta', 'fas'])
    
    if uploaded_file and patterns:
        # Save uploaded file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp_file:
            tmp_file.write(uploaded_file.getvalue().decode())
            tmp_path = tmp_file.name
        
        # Configure scanner
        config = ScanConfig(
            chunk_size=chunk_size,
            overlap=overlap,
            max_workers=max_workers,
            stream_results=True,
            incremental_save=True
        )
        
        # Get cached scanner
        scanner = get_cached_scanner(patterns, config.__dict__)
        
        # Scan button
        if st.button("Start Scanning", type="primary"):
            try:
                with st.spinner("Scanning sequences..."):
                    results = streamlit_scan_with_progress(scanner, tmp_path)
                
                # Display results
                st.success(f"Scanning complete! Found {results['stats']['total_matches']} matches")
                
                # Show statistics
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Total Matches", results['stats']['total_matches'])
                with col2:
                    st.metric("Sequences Scanned", results['stats']['sequences_scanned'])
                with col3:
                    st.metric("Patterns Matched", results['stats']['patterns_matched'])
                
                # Download results
                if 'csv' in results['files']:
                    with open(results['files']['csv'], 'r') as f:
                        st.download_button(
                            "Download Results (CSV)",
                            f.read(),
                            file_name="scan_results.csv",
                            mime="text/csv"
                        )
                
            except Exception as e:
                st.error(f"Scanning failed: {e}")
            finally:
                # Cleanup
                try:
                    os.unlink(tmp_path)
                except Exception:
                    pass


if __name__ == "__main__":
    # Run example if executed directly
    if STREAMLIT_AVAILABLE:
        example_streamlit_app()
    else:
        print("This module provides production-ready Hyperscan-Streamlit integration.")
        print("Install Streamlit to run the example app: pip install streamlit")
        print("Then run: streamlit run production_hyperscan_streamlit.py")