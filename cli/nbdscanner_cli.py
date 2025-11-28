#!/usr/bin/env python3
"""
NBDScanner CLI - Headless Non-B DNA Motif Detection
====================================================

Command-line interface for scanning large genomes for Non-B DNA motifs.
Designed for long-running analysis jobs with:
- Streaming output (NDJSON for incremental results)
- Parallel processing with shared memory
- Progress reporting
- Memory-efficient FASTA processing

Usage:
    python -m cli.nbdscanner_cli --input genome.fa --out results.ndjson
    python cli/nbdscanner_cli.py --input genome.fa --out results.ndjson --workers 4
    
Examples:
    # Basic scan with auto-detected backend
    python cli/nbdscanner_cli.py --input ecoli.fa --out ecoli_motifs.ndjson
    
    # Parallel scan with 8 workers
    python cli/nbdscanner_cli.py --input human_chr1.fa --out chr1.ndjson --workers 8
    
    # Specify chunk size for large genomes
    python cli/nbdscanner_cli.py --input genome.fa --out results.ndjson --chunksize 50000
    
    # Export to multiple formats
    python cli/nbdscanner_cli.py --input seq.fa --out results --format ndjson,csv,bed

Author: Dr. Venkata Rajesh Yella
License: MIT
"""

import argparse
import json
import os
import sys
import time
from typing import Dict, List, Any, Optional
from datetime import datetime

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class NBDScannerCLI:
    """
    Command-line interface for NBDScanner.
    
    Provides a headless scanning mode suitable for:
    - Long-running genome scans
    - Batch processing
    - Integration with pipelines
    - Server/daemon mode
    """
    
    def __init__(self, args: Optional[argparse.Namespace] = None):
        """
        Initialize CLI scanner.
        
        Args:
            args: Parsed command-line arguments
        """
        self.args = args or self.parse_args()
        self.start_time = None
        self.total_bp_scanned = 0
        self.total_motifs_found = 0
    
    @staticmethod
    def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
        """
        Parse command-line arguments.
        
        Args:
            argv: Argument list (default: sys.argv)
            
        Returns:
            Parsed arguments namespace
        """
        parser = argparse.ArgumentParser(
            description='NBDScanner CLI - Non-B DNA Motif Detection',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='''
Examples:
  %(prog)s --input ecoli.fa --out results.ndjson
  %(prog)s --input genome.fa --out results.ndjson --workers 8 --chunksize 100000
  %(prog)s --input sequences.fa --out output --format ndjson,csv,bed
            '''
        )
        
        # Input/Output
        parser.add_argument(
            '--input', '-i',
            required=True,
            help='Input FASTA file path'
        )
        parser.add_argument(
            '--out', '-o',
            required=True,
            help='Output file path (without extension for multiple formats)'
        )
        parser.add_argument(
            '--format', '-f',
            default='ndjson',
            help='Output format(s): ndjson, json, csv, bed (comma-separated for multiple)'
        )
        
        # Performance
        parser.add_argument(
            '--backend', '-b',
            choices=['auto', 'hyperscan', 'numba', 'python'],
            default='auto',
            help='Scanning backend (default: auto-detect best available)'
        )
        parser.add_argument(
            '--workers', '-w',
            type=int,
            default=None,
            help='Number of parallel workers (default: CPU count)'
        )
        parser.add_argument(
            '--chunksize', '-c',
            type=int,
            default=100000,
            help='Chunk size in bp for parallel processing (default: 100000)'
        )
        
        # Options
        parser.add_argument(
            '--streaming',
            action='store_true',
            help='Enable streaming mode (write results incrementally)'
        )
        parser.add_argument(
            '--quiet', '-q',
            action='store_true',
            help='Suppress progress output'
        )
        parser.add_argument(
            '--verbose', '-v',
            action='store_true',
            help='Enable verbose output'
        )
        parser.add_argument(
            '--no-parallel',
            action='store_true',
            help='Disable parallel processing'
        )
        
        # Filtering
        parser.add_argument(
            '--min-score',
            type=float,
            default=0.0,
            help='Minimum score threshold for reporting motifs'
        )
        parser.add_argument(
            '--min-length',
            type=int,
            default=0,
            help='Minimum motif length to report'
        )
        parser.add_argument(
            '--classes',
            type=str,
            default=None,
            help='Comma-separated list of motif classes to detect (default: all)'
        )
        
        return parser.parse_args(argv)
    
    def log(self, message: str, level: str = 'info'):
        """
        Log message with timestamp.
        
        Args:
            message: Message to log
            level: Log level (info, warn, error, debug)
        """
        if self.args.quiet and level == 'info':
            return
        if not self.args.verbose and level == 'debug':
            return
        
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        prefix = {'info': '→', 'warn': '⚠', 'error': '✗', 'debug': '…'}
        print(f"[{timestamp}] {prefix.get(level, '•')} {message}", file=sys.stderr)
    
    def progress_callback(self, completed: int, total: int):
        """
        Progress callback for parallel scanning.
        
        Args:
            completed: Number of completed chunks
            total: Total number of chunks
        """
        if not self.args.quiet:
            percent = (completed / total) * 100
            elapsed = time.time() - self.start_time
            eta = (elapsed / completed * (total - completed)) if completed > 0 else 0
            print(
                f"\rProgress: {completed}/{total} chunks ({percent:.1f}%) "
                f"- Elapsed: {elapsed:.1f}s - ETA: {eta:.1f}s",
                end='', file=sys.stderr
            )
            if completed == total:
                print(file=sys.stderr)  # Newline at end
    
    def scan_file(self) -> Dict[str, List[Dict[str, Any]]]:
        """
        Scan input FASTA file.
        
        Returns:
            Dictionary mapping sequence names to motif lists
        """
        from utilities import read_fasta_file
        from nonbscanner import analyze_sequence
        
        input_path = self.args.input
        
        if not os.path.exists(input_path):
            self.log(f"Input file not found: {input_path}", 'error')
            sys.exit(1)
        
        self.log(f"Reading FASTA file: {input_path}")
        sequences = read_fasta_file(input_path)
        
        if not sequences:
            self.log("No sequences found in input file", 'error')
            sys.exit(1)
        
        self.log(f"Found {len(sequences)} sequence(s)")
        
        results = {}
        
        # Determine if parallel processing should be used
        use_parallel = not self.args.no_parallel
        
        for seq_name, sequence in sequences.items():
            seq_len = len(sequence)
            self.log(f"Scanning {seq_name} ({seq_len:,} bp)")
            
            self.start_time = time.time()
            
            if use_parallel and seq_len > self.args.chunksize * 2:
                # Use parallel scanning for large sequences
                try:
                    from scanner_backends.parallel_worker import parallel_scan
                    motifs = parallel_scan(
                        sequence,
                        sequence_name=seq_name,
                        num_workers=self.args.workers,
                        chunk_size=self.args.chunksize,
                        progress_callback=self.progress_callback
                    )
                except Exception as e:
                    self.log(f"Parallel scan failed, falling back to sequential: {e}", 'warn')
                    motifs = analyze_sequence(sequence, seq_name)
            else:
                # Use sequential scanning for small sequences
                motifs = analyze_sequence(sequence, seq_name)
            
            elapsed = time.time() - self.start_time
            
            # Apply filters
            motifs = self.filter_motifs(motifs)
            
            results[seq_name] = motifs
            
            self.total_bp_scanned += seq_len
            self.total_motifs_found += len(motifs)
            
            speed = seq_len / elapsed if elapsed > 0 else 0
            self.log(
                f"Found {len(motifs)} motifs in {elapsed:.2f}s "
                f"({speed:,.0f} bp/s)"
            )
        
        return results
    
    def filter_motifs(self, motifs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Apply filtering criteria to motifs.
        
        Args:
            motifs: List of detected motifs
            
        Returns:
            Filtered list of motifs
        """
        filtered = motifs
        
        # Filter by score
        if self.args.min_score > 0:
            filtered = [m for m in filtered if m.get('Score', 0) >= self.args.min_score]
        
        # Filter by length
        if self.args.min_length > 0:
            filtered = [m for m in filtered if m.get('Length', 0) >= self.args.min_length]
        
        # Filter by class
        if self.args.classes:
            allowed_classes = [c.strip() for c in self.args.classes.split(',')]
            filtered = [m for m in filtered if m.get('Class', '') in allowed_classes]
        
        return filtered
    
    def write_results(self, results: Dict[str, List[Dict[str, Any]]]):
        """
        Write results to output file(s).
        
        Args:
            results: Dictionary mapping sequence names to motif lists
        """
        formats = [f.strip().lower() for f in self.args.format.split(',')]
        output_base = self.args.out
        
        # Flatten results for export
        all_motifs = []
        for seq_name, motifs in results.items():
            all_motifs.extend(motifs)
        
        for fmt in formats:
            if fmt == 'ndjson':
                self.write_ndjson(all_motifs, output_base + '.ndjson')
            elif fmt == 'json':
                self.write_json(all_motifs, output_base + '.json')
            elif fmt == 'csv':
                self.write_csv(all_motifs, output_base + '.csv')
            elif fmt == 'bed':
                self.write_bed(all_motifs, output_base + '.bed')
            else:
                self.log(f"Unknown format: {fmt}", 'warn')
    
    def write_ndjson(self, motifs: List[Dict[str, Any]], filepath: str):
        """Write motifs as NDJSON (newline-delimited JSON)."""
        self.log(f"Writing NDJSON to {filepath}")
        with open(filepath, 'w') as f:
            for motif in motifs:
                f.write(json.dumps(motif, ensure_ascii=True) + '\n')
    
    def write_json(self, motifs: List[Dict[str, Any]], filepath: str):
        """Write motifs as JSON array."""
        self.log(f"Writing JSON to {filepath}")
        with open(filepath, 'w') as f:
            json.dump({
                'version': '2024.1',
                'timestamp': datetime.now().isoformat(),
                'total_motifs': len(motifs),
                'motifs': motifs
            }, f, indent=2)
    
    def write_csv(self, motifs: List[Dict[str, Any]], filepath: str):
        """Write motifs as CSV."""
        self.log(f"Writing CSV to {filepath}")
        from utilities import export_to_csv
        export_to_csv(motifs, filepath)
    
    def write_bed(self, motifs: List[Dict[str, Any]], filepath: str):
        """Write motifs as BED format."""
        self.log(f"Writing BED to {filepath}")
        from utilities import export_to_bed
        # Get sequence name from first motif
        seq_name = motifs[0].get('Sequence_Name', 'sequence') if motifs else 'sequence'
        export_to_bed(motifs, seq_name, filepath)
    
    def print_summary(self):
        """Print summary statistics."""
        if self.args.quiet:
            return
        
        print("\n" + "=" * 60, file=sys.stderr)
        print("SCAN SUMMARY", file=sys.stderr)
        print("=" * 60, file=sys.stderr)
        print(f"Total sequences scanned: {self.total_bp_scanned:,} bp", file=sys.stderr)
        print(f"Total motifs found: {self.total_motifs_found:,}", file=sys.stderr)
        
        if self.total_bp_scanned > 0:
            density = self.total_motifs_found / (self.total_bp_scanned / 1000)
            print(f"Motif density: {density:.2f} per kb", file=sys.stderr)
        
        print("=" * 60, file=sys.stderr)
    
    def run(self) -> int:
        """
        Run the CLI scanner.
        
        Returns:
            Exit code (0 for success)
        """
        self.log("NBDScanner CLI starting")
        self.log(f"Backend: {self.args.backend}")
        
        # Show backend info
        try:
            from scanner_backends import get_best_backend, HYPERSCAN_AVAILABLE, NUMBA_AVAILABLE
            actual_backend = get_best_backend()
            self.log(f"Using backend: {actual_backend}", 'debug')
            self.log(f"Hyperscan available: {HYPERSCAN_AVAILABLE}", 'debug')
            self.log(f"Numba available: {NUMBA_AVAILABLE}", 'debug')
        except ImportError:
            self.log("Scanner backends not available, using standard scanner", 'debug')
        
        try:
            # Scan input file
            results = self.scan_file()
            
            # Write results
            self.write_results(results)
            
            # Print summary
            self.print_summary()
            
            self.log("Scan complete!")
            return 0
            
        except KeyboardInterrupt:
            self.log("Scan interrupted by user", 'warn')
            return 130
        except Exception as e:
            self.log(f"Error: {e}", 'error')
            if self.args.verbose:
                import traceback
                traceback.print_exc()
            return 1


def main(argv: Optional[List[str]] = None) -> int:
    """
    Main entry point for CLI.
    
    Args:
        argv: Command-line arguments (default: sys.argv)
        
    Returns:
        Exit code
    """
    args = NBDScannerCLI.parse_args(argv)
    cli = NBDScannerCLI(args)
    return cli.run()


if __name__ == '__main__':
    sys.exit(main())
