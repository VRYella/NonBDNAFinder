#!/usr/bin/env python3
"""
Batch Genome Analysis Script
=============================

Analyzes all genome files in the Genomes folder and generates comprehensive results.
This script processes multiple FASTA/FNA files in batch mode with simplified progress display.

Usage:
    python batch_analysis.py [--genomes-dir GENOMES_DIR] [--output-dir OUTPUT_DIR]

Features:
    - Automatic genome file discovery (.fna/.fasta)
    - Simplified progress display (time + chunks only)
    - Consolidated results storage
    - Per-genome detailed analysis
    - Summary statistics across all genomes
"""

import os
import sys
import time
import argparse
import glob
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any, Tuple
import json

# Add current directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import analysis modules
from engine.detection import analyze_sequence, NonBScanner
from utils.fasta import parse_fasta, read_fasta_file
from utils.stats import get_basic_stats


class BatchAnalysisProgress:
    """Simple progress tracker for batch analysis"""
    
    def __init__(self, total_genomes: int):
        self.total_genomes = total_genomes
        self.current_genome = 0
        self.start_time = time.time()
        self.genome_times = []
    
    def start_genome(self, genome_name: str):
        """Start processing a genome"""
        self.current_genome += 1
        self.genome_start = time.time()
        elapsed = time.time() - self.start_time
        print(f"\n{'='*80}")
        print(f"[{self.current_genome}/{self.total_genomes}] Processing: {genome_name}")
        print(f"Elapsed time: {self._format_time(elapsed)}")
        print(f"{'='*80}")
    
    def update_chunks(self, current_chunk: int, total_chunks: int, chunk_time: float):
        """Update chunk progress"""
        print(f"  Chunk {current_chunk}/{total_chunks} completed in {chunk_time:.2f}s")
    
    def finish_genome(self, genome_name: str, motif_count: int):
        """Finish processing a genome"""
        genome_time = time.time() - self.genome_start
        self.genome_times.append(genome_time)
        print(f"\n  ✓ {genome_name} completed: {motif_count} motifs in {self._format_time(genome_time)}")
        
        # Show ETA for remaining genomes
        if self.current_genome < self.total_genomes:
            avg_time = sum(self.genome_times) / len(self.genome_times)
            remaining = self.total_genomes - self.current_genome
            eta = avg_time * remaining
            print(f"  → ETA for remaining {remaining} genomes: {self._format_time(eta)}")
    
    def _format_time(self, seconds: float) -> str:
        """Format time in human-readable format"""
        if seconds < 60:
            return f"{seconds:.1f}s"
        elif seconds < 3600:
            mins = int(seconds // 60)
            secs = int(seconds % 60)
            return f"{mins}m {secs}s"
        else:
            hours = int(seconds // 3600)
            mins = int((seconds % 3600) // 60)
            return f"{hours}h {mins}m"
    
    def summary(self):
        """Print batch analysis summary"""
        total_time = time.time() - self.start_time
        print(f"\n{'='*80}")
        print(f"BATCH ANALYSIS COMPLETE")
        print(f"{'='*80}")
        print(f"Total genomes processed: {self.total_genomes}")
        print(f"Total time: {self._format_time(total_time)}")
        print(f"Average time per genome: {self._format_time(sum(self.genome_times) / len(self.genome_times))}")
        print(f"{'='*80}\n")


def find_genome_files(genomes_dir: str) -> List[str]:
    """
    Find all genome files in the specified directory.
    
    Args:
        genomes_dir: Directory containing genome files
        
    Returns:
        List of genome file paths
    """
    patterns = ['*.fna', '*.fasta', '*.fa']
    genome_files = []
    
    for pattern in patterns:
        genome_files.extend(glob.glob(os.path.join(genomes_dir, pattern)))
    
    # Sort by filename for consistent order
    genome_files.sort()
    
    return genome_files


def analyze_genome_file(genome_path: str, progress: BatchAnalysisProgress) -> Dict[str, Any]:
    """
    Analyze a single genome file.
    
    Args:
        genome_path: Path to genome file
        progress: Progress tracker
        
    Returns:
        Dictionary containing analysis results
    """
    genome_name = os.path.basename(genome_path)
    progress.start_genome(genome_name)
    
    # Read genome file
    try:
        sequences = read_fasta_file(genome_path)
        if not sequences:
            print(f"  ✗ ERROR: No sequences found in {genome_name}")
            return None
    except Exception as e:
        print(f"  ✗ ERROR reading {genome_name}: {str(e)}")
        return None
    
    # Process all sequences in the file
    all_motifs = []
    sequence_info = []
    
    for seq_name, seq_data in sequences:
        sequence = seq_data.upper()
        print(f"  → Analyzing sequence: {seq_name} ({len(sequence):,} bp)")
        
        # Analyze sequence
        try:
            motifs = analyze_sequence(sequence, seq_name)
            all_motifs.extend(motifs)
            
            # Calculate statistics
            stats = get_basic_stats(sequence, motifs)
            sequence_info.append({
                'name': seq_name,
                'length': len(sequence),
                'motif_count': len(motifs),
                'stats': stats
            })
            
            print(f"    Found {len(motifs)} motifs")
            
        except Exception as e:
            print(f"  ✗ ERROR analyzing {seq_name}: {str(e)}")
            continue
    
    # Finish genome processing
    progress.finish_genome(genome_name, len(all_motifs))
    
    return {
        'genome_file': genome_name,
        'genome_path': genome_path,
        'sequences': sequence_info,
        'total_motifs': len(all_motifs),
        'motifs': all_motifs,
        'timestamp': datetime.now().isoformat()
    }


def save_results(results: List[Dict[str, Any]], output_dir: str):
    """
    Save batch analysis results.
    
    Args:
        results: List of genome analysis results
        output_dir: Output directory for results
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Save individual genome results
    for result in results:
        if result is None:
            continue
        
        genome_name = result['genome_file']
        safe_name = genome_name.replace('.fna', '').replace('.fasta', '').replace(' ', '_')
        
        # Save motifs as JSON
        motifs_file = os.path.join(output_dir, f"{safe_name}_motifs.json")
        with open(motifs_file, 'w') as f:
            json.dump(result['motifs'], f, indent=2)
        
        print(f"  Saved: {motifs_file}")
    
    # Save batch summary
    summary = {
        'analysis_date': datetime.now().isoformat(),
        'total_genomes': len([r for r in results if r is not None]),
        'genomes': []
    }
    
    for result in results:
        if result is None:
            continue
        
        summary['genomes'].append({
            'genome_file': result['genome_file'],
            'sequences': result['sequences'],
            'total_motifs': result['total_motifs']
        })
    
    summary_file = os.path.join(output_dir, 'batch_summary.json')
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n  Saved batch summary: {summary_file}")


def main():
    """Main batch analysis function"""
    parser = argparse.ArgumentParser(
        description='Batch analysis of genome files in Genomes folder'
    )
    parser.add_argument(
        '--genomes-dir',
        default='Genomes',
        help='Directory containing genome files (default: Genomes)'
    )
    parser.add_argument(
        '--output-dir',
        default='batch_results',
        help='Output directory for results (default: batch_results)'
    )
    
    args = parser.parse_args()
    
    # Find genome files
    print("="*80)
    print("BATCH GENOME ANALYSIS")
    print("="*80)
    print(f"Genomes directory: {args.genomes_dir}")
    print(f"Output directory: {args.output_dir}")
    print()
    
    genome_files = find_genome_files(args.genomes_dir)
    
    if not genome_files:
        print(f"ERROR: No genome files found in {args.genomes_dir}")
        print("Looking for: *.fna, *.fasta, *.fa")
        return 1
    
    print(f"Found {len(genome_files)} genome file(s):")
    for i, gf in enumerate(genome_files, 1):
        print(f"  {i}. {os.path.basename(gf)}")
    
    # Initialize progress tracker
    progress = BatchAnalysisProgress(len(genome_files))
    
    # Process each genome
    results = []
    for genome_path in genome_files:
        result = analyze_genome_file(genome_path, progress)
        results.append(result)
    
    # Save all results
    print(f"\n{'='*80}")
    print("SAVING RESULTS")
    print(f"{'='*80}")
    save_results(results, args.output_dir)
    
    # Print summary
    progress.summary()
    
    print(f"All results saved to: {args.output_dir}/")
    print("Next steps:")
    print("  1. Run comparative analysis: python comparative_analysis.py")
    print("  2. Generate publication outputs: python publication_export.py")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
