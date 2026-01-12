#!/usr/bin/env python3
"""
Master Batch Analysis Pipeline
===============================

Runs the complete batch genome analysis pipeline:
1. Batch analysis of all genomes
2. Comparative analysis
3. Publication Excel export
4. Publication narrative generation
5. Publication visualizations

Usage:
    python run_batch_pipeline.py [--genomes-dir DIR] [--output-base DIR]
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path
from datetime import datetime


def run_command(cmd: list, description: str) -> int:
    """
    Run a command and display progress.
    
    Args:
        cmd: Command as list of strings
        description: Description of the step
        
    Returns:
        Return code (0 for success)
    """
    print("\n" + "="*80)
    print(f"STEP: {description}")
    print("="*80)
    print(f"Running: {' '.join(cmd)}\n")
    
    result = subprocess.run(cmd)
    
    if result.returncode != 0:
        print(f"\n✗ ERROR: {description} failed with code {result.returncode}")
        return result.returncode
    
    print(f"\n✓ {description} completed successfully")
    return 0


def main():
    """Main pipeline function"""
    parser = argparse.ArgumentParser(
        description='Run complete batch genome analysis pipeline'
    )
    parser.add_argument(
        '--genomes-dir',
        default='Genomes',
        help='Directory containing genome files (default: Genomes)'
    )
    parser.add_argument(
        '--output-base',
        default='.',
        help='Base output directory (default: current directory)'
    )
    parser.add_argument(
        '--skip-batch',
        action='store_true',
        help='Skip batch analysis (use existing results)'
    )
    parser.add_argument(
        '--skip-comparative',
        action='store_true',
        help='Skip comparative analysis (use existing results)'
    )
    
    args = parser.parse_args()
    
    # Define output directories
    batch_dir = os.path.join(args.output_base, 'batch_results')
    comparative_dir = os.path.join(args.output_base, 'comparative_results')
    figures_dir = os.path.join(args.output_base, 'publication_figures')
    excel_file = os.path.join(args.output_base, 'publication_results.xlsx')
    narrative_file = os.path.join(args.output_base, 'publication_narrative.txt')
    
    print("="*80)
    print("BATCH GENOME ANALYSIS PIPELINE")
    print("="*80)
    print(f"Genomes directory: {args.genomes_dir}")
    print(f"Output base: {args.output_base}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    
    start_time = datetime.now()
    
    # Step 1: Batch Analysis
    if not args.skip_batch:
        if run_command(
            ['python', 'batch_analysis.py', 
             '--genomes-dir', args.genomes_dir,
             '--output-dir', batch_dir],
            'Batch Genome Analysis'
        ) != 0:
            return 1
    else:
        print("\n⊳ Skipping batch analysis (using existing results)")
    
    # Step 2: Comparative Analysis
    if not args.skip_comparative:
        if run_command(
            ['python', 'comparative_analysis.py',
             '--input-dir', batch_dir,
             '--output-dir', comparative_dir],
            'Comparative Analysis'
        ) != 0:
            return 1
    else:
        print("\n⊳ Skipping comparative analysis (using existing results)")
    
    # Step 3: Publication Excel Export
    if run_command(
        ['python', 'publication_export.py',
         '--batch-dir', batch_dir,
         '--comparative-dir', comparative_dir,
         '--output', excel_file],
        'Publication Excel Export'
    ) != 0:
        return 1
    
    # Step 4: Publication Narrative
    if run_command(
        ['python', 'publication_narrative.py',
         '--comparative-dir', comparative_dir,
         '--output', narrative_file],
        'Publication Narrative Generation'
    ) != 0:
        return 1
    
    # Step 5: Publication Visualizations
    if run_command(
        ['python', 'publication_visualizations.py',
         '--comparative-dir', comparative_dir,
         '--output-dir', figures_dir],
        'Publication Visualizations'
    ) != 0:
        return 1
    
    # Final summary
    end_time = datetime.now()
    duration = end_time - start_time
    
    print("\n" + "="*80)
    print("PIPELINE COMPLETE!")
    print("="*80)
    print(f"Started:  {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Finished: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Duration: {duration}")
    print()
    print("Generated outputs:")
    print(f"  • Batch results:       {batch_dir}/")
    print(f"  • Comparative results: {comparative_dir}/")
    print(f"  • Excel workbook:      {excel_file}")
    print(f"  • Narrative text:      {narrative_file}")
    print(f"  • Figures:             {figures_dir}/")
    print("="*80)
    
    print("\nNext steps:")
    print("  1. Review publication_results.xlsx for comprehensive data")
    print("  2. Read publication_narrative.txt for manuscript sections")
    print("  3. Use figures in publication_figures/ for manuscript figures")
    print()
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
