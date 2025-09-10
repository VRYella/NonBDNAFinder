"""
Pipeline orchestrator - main entry point for NBDFinder analysis.

Coordinates: FASTA loading -> pattern compilation -> parallel detection -> 
scoring -> hybrid/cluster detection -> normalization -> output generation.
"""

import os
import sys
import argparse
import logging
from typing import List, Dict, Optional
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import multiprocessing
from io import StringIO

import pandas as pd
from Bio import SeqIO

from motifs.base import Candidate, normalize_scores
from motifs.registry import get_all_hyperscan_patterns, get_hyperscan_safe_patterns
import hs_registry_manager
from motif_detectors import get_all_detectors, get_detector
from overlap_resolution import EnhancedOverlapResolver, OverlapConfig, OverlapStrategy

logger = logging.getLogger(__name__)


def setup_logging(level: str = "INFO"):
    """Setup logging configuration"""
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )


def load_fasta_sequences(fasta_path: str) -> List[tuple]:
    """
    Load sequences from FASTA file
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        List of (seq_id, sequence) tuples
    """
    sequences = []
    
    try:
        with open(fasta_path, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences.append((record.id, str(record.seq)))
        
        logger.info(f"Loaded {len(sequences)} sequences from {fasta_path}")
        return sequences
        
    except Exception as e:
        logger.error(f"Failed to load FASTA file {fasta_path}: {e}")
        raise


def compile_hyperscan_databases() -> None:
    """Compile Hyperscan databases for all motif classes"""
    logger.info("Compiling Hyperscan databases...")
    
    if not hs_registry_manager.is_available():
        logger.warning("Hyperscan not available, will use fallback regex")
        return
    
    # Compile global database with all patterns
    all_patterns = get_all_hyperscan_patterns()
    if all_patterns:
        pattern_strings = [p[0] for p in all_patterns]
        hs_registry_manager.compile_db(pattern_strings, key="global")
        logger.info(f"Compiled global database with {len(pattern_strings)} patterns")
    
    # Compile per-class databases for better performance
    class_patterns = get_hyperscan_safe_patterns()
    for class_name, patterns in class_patterns.items():
        if patterns:
            hs_registry_manager.compile_db(patterns, key=f"{class_name}_patterns")
            logger.info(f"Compiled {class_name} database with {len(patterns)} patterns")


def detect_motifs_in_sequence(args: tuple) -> List[Candidate]:
    """
    Worker function for detecting motifs in a single sequence
    
    Args:
        args: Tuple of (seq_id, sequence, detector_classes, offset)
        
    Returns:
        List of detected candidates
    """
    seq_id, sequence, detector_classes, offset = args
    all_candidates = []
    
    # Process each detector class
    for class_name in detector_classes:
        try:
            detector = get_detector(class_name)
            
            # Skip hybrid and cluster detectors in primary detection
            if class_name in ['hybrid', 'cluster']:
                continue
            
            # Detect candidates
            candidates = detector.detect(sequence, seq_id, seq_id, offset)
            
            # Score candidates
            scored_candidates = detector.score(candidates)
            all_candidates.extend(scored_candidates)
            
            logger.debug(f"Found {len(candidates)} {class_name} candidates in {seq_id}")
            
        except Exception as e:
            logger.error(f"Detection failed for {class_name} in {seq_id}: {e}")
    
    return all_candidates


def chunk_sequence(sequence: str, chunk_size: int = 50000, overlap: int = 5000) -> List[tuple]:
    """
    Split long sequence into overlapping chunks
    
    Args:
        sequence: DNA sequence string
        chunk_size: Size of each chunk
        overlap: Overlap between chunks
        
    Returns:
        List of (chunk_sequence, offset) tuples
    """
    chunks = []
    seq_len = len(sequence)
    
    if seq_len <= chunk_size:
        return [(sequence, 0)]
    
    offset = 0
    while offset < seq_len:
        end = min(offset + chunk_size, seq_len)
        chunk = sequence[offset:end]
        chunks.append((chunk, offset))
        
        if end >= seq_len:
            break
            
        offset += chunk_size - overlap
    
    logger.debug(f"Split sequence of length {seq_len} into {len(chunks)} chunks")
    return chunks


def merge_overlapping_candidates(candidates: List[Candidate], overlap_size: int = 5000) -> List[Candidate]:
    """
    Remove duplicate candidates from overlapping chunks and resolve overlaps.
    
    This function now performs two operations:
    1. Remove exact duplicates from chunked processing
    2. Apply enhanced overlap resolution within same subclass
    
    Args:
        candidates: List of candidates from chunked processing
        overlap_size: Size of overlap between chunks (for duplicate detection)
        
    Returns:
        Deduplicated and overlap-resolved candidates
    """
    if not candidates:
        return candidates
    
    # Sort by position
    sorted_candidates = sorted(candidates, key=lambda c: (c.sequence_name, c.start))
    
    # Step 1: Remove exact duplicates from chunked processing
    unique_candidates = []
    seen_positions = set()
    
    for candidate in sorted_candidates:
        position_key = (candidate.sequence_name, candidate.start, candidate.end, candidate.class_name)
        if position_key not in seen_positions:
            unique_candidates.append(candidate)
            seen_positions.add(position_key)
    
    logger.debug(f"Removed duplicates: {len(candidates)} -> {len(unique_candidates)} candidates")
    
    # Step 2: Apply enhanced overlap resolution
    overlap_config = OverlapConfig(
        strategy=OverlapStrategy.HIGHEST_SCORE,
        min_overlap_percent=0.1,
        same_class_only=True,  # Resolve overlaps within same subclass
        preserve_hybrid=True
    )
    
    resolver = EnhancedOverlapResolver(overlap_config)
    resolved_candidates = resolver.resolve_overlaps(unique_candidates)
    
    logger.debug(f"Overlap resolution: {len(unique_candidates)} -> {len(resolved_candidates)} candidates")
    return resolved_candidates


def detect_hybrid_and_cluster_motifs(candidates: List[Candidate]) -> List[Candidate]:
    """
    Detect hybrid and cluster motifs from primary candidates
    
    Args:
        candidates: List of primary motif candidates
        
    Returns:
        List including hybrid and cluster candidates
    """
    all_candidates = candidates.copy()
    
    # Detect hybrid motifs
    try:
        hybrid_detector = get_detector('hybrid')
        hybrid_candidates = hybrid_detector.detect_from_candidates(candidates)
        scored_hybrids = hybrid_detector.score(hybrid_candidates)
        all_candidates.extend(scored_hybrids)
        logger.info(f"Found {len(hybrid_candidates)} hybrid motifs")
    except Exception as e:
        logger.error(f"Hybrid detection failed: {e}")
    
    # Detect cluster motifs
    try:
        cluster_detector = get_detector('cluster')
        cluster_candidates = cluster_detector.detect_from_candidates(candidates)
        scored_clusters = cluster_detector.score(cluster_candidates)
        all_candidates.extend(scored_clusters)
        logger.info(f"Found {len(cluster_candidates)} cluster motifs")
    except Exception as e:
        logger.error(f"Cluster detection failed: {e}")
    
    return all_candidates


def candidates_to_dataframe(candidates: List[Candidate]) -> pd.DataFrame:
    """
    Convert candidates to pandas DataFrame with specified schema
    
    Args:
        candidates: List of Candidate objects
        
    Returns:
        DataFrame with standardized columns
    """
    if not candidates:
        return pd.DataFrame(columns=[
            'S.No', 'Sequence_Name', 'Chromosome/Contig', 'Class', 'Subclass',
            'Motif_ID', 'Start', 'End', 'Length', 'Normalized_Score', 
            'Actual_Score', 'Scoring_Method', 'GC_Content', 'Sequence', 'Overlap_Classes'
        ])
    
    data = []
    for i, candidate in enumerate(candidates, 1):
        row = {
            'S.No': i,
            'Sequence_Name': candidate.sequence_name,
            'Chromosome/Contig': candidate.contig,
            'Class': candidate.class_name,
            'Subclass': candidate.subclass or '',
            'Motif_ID': candidate.motif_id,
            'Start': candidate.start,
            'End': candidate.end,
            'Length': candidate.length,
            'Normalized_Score': candidate.normalized_score or 0.0,
            'Actual_Score': candidate.raw_score or 0.0,
            'Scoring_Method': candidate.scoring_method or '',
            'GC_Content': candidate.gc_content or 0.0,
            'Sequence': candidate.matched_seq.decode('utf-8') if isinstance(candidate.matched_seq, bytes) else str(candidate.matched_seq),
            'Overlap_Classes': ','.join(candidate.overlap_classes) if candidate.overlap_classes else ''
        }
        data.append(row)
    
    return pd.DataFrame(data)


def export_results(df: pd.DataFrame, output_prefix: str = "results") -> Dict[str, str]:
    """
    Export results in multiple formats
    
    Args:
        df: Results DataFrame
        output_prefix: Prefix for output files
        
    Returns:
        Dictionary of format -> filename
    """
    output_files = {}
    
    # CSV export
    csv_file = f"{output_prefix}.csv"
    df.to_csv(csv_file, index=False)
    output_files['csv'] = csv_file
    logger.info(f"Exported CSV: {csv_file}")
    
    # Excel export
    try:
        excel_file = f"{output_prefix}.xlsx"
        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='All_Motifs', index=False)
            
            # Create separate sheets by class
            for class_name in df['Class'].unique():
                class_df = df[df['Class'] == class_name]
                sheet_name = class_name.replace('_', ' ').title()[:31]  # Excel sheet name limit
                class_df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        output_files['excel'] = excel_file
        logger.info(f"Exported Excel: {excel_file}")
    except Exception as e:
        logger.warning(f"Excel export failed: {e}")
    
    # Parquet export
    try:
        parquet_file = f"{output_prefix}.parquet"
        df.to_parquet(parquet_file, index=False)
        output_files['parquet'] = parquet_file
        logger.info(f"Exported Parquet: {parquet_file}")
    except Exception as e:
        logger.warning(f"Parquet export failed: {e}")
    
    # GFF3 export
    try:
        gff_file = f"{output_prefix}.gff3"
        export_gff3(df, gff_file)
        output_files['gff3'] = gff_file
        logger.info(f"Exported GFF3: {gff_file}")
    except Exception as e:
        logger.warning(f"GFF3 export failed: {e}")
    
    return output_files


def export_gff3(df: pd.DataFrame, output_file: str):
    """
    Export results in GFF3 format
    
    Args:
        df: Results DataFrame
        output_file: Output GFF3 filename
    """
    with open(output_file, 'w') as f:
        f.write("##gff-version 3\n")
        
        for _, row in df.iterrows():
            # GFF3 format: seqname\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes
            attributes = f"ID=motif_{row['S.No']};Class={row['Class']};Subclass={row['Subclass']};Score={row['Normalized_Score']:.3f};Method={row['Scoring_Method']}"
            
            gff_line = f"{row['Chromosome/Contig']}\tNBDFinder\tmotif\t{row['Start']}\t{row['End']}\t{row['Normalized_Score']:.3f}\t.\t.\t{attributes}\n"
            f.write(gff_line)


def run_pipeline(fasta_path: str, output_prefix: str = "results", 
                max_workers: Optional[int] = None, chunk_size: int = 50000,
                detector_classes: Optional[List[str]] = None) -> Dict[str, str]:
    """
    Main pipeline function
    
    Args:
        fasta_path: Path to input FASTA file
        output_prefix: Prefix for output files
        max_workers: Number of parallel workers (default: CPU count)
        chunk_size: Size for sequence chunking
        detector_classes: List of detector classes to run (default: all)
        
    Returns:
        Dictionary of output files created
    """
    logger.info(f"Starting NBDFinder pipeline on {fasta_path}")
    
    # Setup
    if max_workers is None:
        max_workers = multiprocessing.cpu_count()
    
    if detector_classes is None:
        detector_classes = ['g_quadruplex', 'curved_dna', 'slipped_dna', 'cruciform',
                           'r_loop', 'triplex', 'i_motif', 'z_dna', 'a_philic', 'hybrid', 'cluster']
    
    # Load sequences
    sequences = load_fasta_sequences(fasta_path)
    if not sequences:
        raise ValueError("No sequences found in FASTA file")
    
    # Compile Hyperscan databases
    compile_hyperscan_databases()
    
    # Prepare work chunks
    work_chunks = []
    for seq_id, sequence in sequences:
        if len(sequence) > chunk_size:
            # Split long sequences
            chunks = chunk_sequence(sequence, chunk_size)
            for chunk_seq, offset in chunks:
                work_chunks.append((seq_id, chunk_seq, detector_classes, offset))
        else:
            work_chunks.append((seq_id, sequence, detector_classes, 0))
    
    logger.info(f"Processing {len(work_chunks)} work chunks with {max_workers} workers")
    
    # Parallel detection
    all_candidates = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_chunk = {
            executor.submit(detect_motifs_in_sequence, chunk): chunk 
            for chunk in work_chunks
        }
        
        for future in as_completed(future_to_chunk):
            try:
                candidates = future.result()
                all_candidates.extend(candidates)
            except Exception as e:
                chunk = future_to_chunk[future]
                logger.error(f"Worker failed for chunk {chunk[0]}: {e}")
    
    logger.info(f"Primary detection found {len(all_candidates)} candidates")
    
    # Merge overlapping candidates
    unique_candidates = merge_overlapping_candidates(all_candidates)
    
    # Detect hybrid and cluster motifs
    final_candidates = detect_hybrid_and_cluster_motifs(unique_candidates)
    
    # Normalize scores
    normalized_candidates = normalize_scores(final_candidates)
    
    logger.info(f"Final result: {len(normalized_candidates)} candidates")
    
    # Convert to DataFrame
    results_df = candidates_to_dataframe(normalized_candidates)
    
    # Export results
    output_files = export_results(results_df, output_prefix)
    
    logger.info("Pipeline completed successfully")
    return output_files


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(description="NBDFinder - Non-B DNA motif detection")
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument("--out", default="results", help="Output prefix (default: results)")
    parser.add_argument("--workers", type=int, help="Number of parallel workers")
    parser.add_argument("--chunk-size", type=int, default=50000, help="Sequence chunk size")
    parser.add_argument("--classes", nargs='+', help="Motif classes to detect")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.log_level)
    
    # Validate input
    if not os.path.exists(args.fasta):
        logger.error(f"FASTA file not found: {args.fasta}")
        sys.exit(1)
    
    # Run pipeline
    try:
        output_files = run_pipeline(
            fasta_path=args.fasta,
            output_prefix=args.out,
            max_workers=args.workers,
            chunk_size=args.chunk_size,
            detector_classes=args.classes
        )
        
        print("Output files created:")
        for format_name, filename in output_files.items():
            print(f"  {format_name.upper()}: {filename}")
            
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()