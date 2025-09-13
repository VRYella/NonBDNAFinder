#!/usr/bin/env python3
"""
Simple command-line launcher for NonBDNA Finder Standalone
Alternative to Jupyter notebook for command-line users
"""

import os
import sys
import argparse
from pathlib import Path

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from nbdfinder_core import analyze_fasta_file
    import pandas as pd
except ImportError as e:
    print(f"❌ Import error: {e}")
    print("Please run: pip install -r requirements.txt")
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="NonBDNA Finder Standalone - Command Line Interface",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python launcher.py input.fasta
  python launcher.py input.fasta --output my_results
  python launcher.py input.fasta --output results --show-stats

Output files:
  - results.csv          - Spreadsheet format
  - results.xlsx         - Excel workbook
  - results.gff3         - Genome browser format
  - results_summary.txt  - Human-readable summary
  - results.zip          - Complete package
        """
    )
    
    parser.add_argument(
        'fasta_file',
        help='Input FASTA file'
    )
    
    parser.add_argument(
        '--output', '-o',
        default='nbdfinder_results',
        help='Output file prefix (default: nbdfinder_results)'
    )
    
    parser.add_argument(
        '--show-stats',
        action='store_true',
        help='Display detailed statistics'
    )
    
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Suppress progress messages'
    )

    args = parser.parse_args()
    
    # Check input file
    if not os.path.exists(args.fasta_file):
        print(f"❌ Error: Input file '{args.fasta_file}' not found")
        sys.exit(1)
    
    if not args.quiet:
        print("🧬 NonBDNA Finder Standalone")
        print("=" * 40)
        print(f"📁 Input file: {args.fasta_file}")
        print(f"📤 Output prefix: {args.output}")
        print("")
    
    try:
        # Read FASTA file
        if not args.quiet:
            print("📖 Reading FASTA file...")
        
        with open(args.fasta_file, 'r') as f:
            fasta_content = f.read()
        
        if not fasta_content.strip().startswith('>'):
            print("❌ Error: Input file does not appear to be in FASTA format")
            sys.exit(1)
        
        # Run analysis
        if not args.quiet:
            print("🔬 Running motif analysis...")
        
        df, zip_data = analyze_fasta_file(fasta_content)
        
        if len(df) == 0:
            print("🔍 No non-B DNA motifs detected in the input sequence(s)")
            print("💡 This could be normal for sequences without detectable motifs")
            return
        
        # Save individual files
        if not args.quiet:
            print(f"💾 Saving results...")
        
        # CSV
        csv_file = f"{args.output}.csv"
        df.to_csv(csv_file, index=False)
        
        # Excel
        excel_file = f"{args.output}.xlsx"
        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='All_Motifs', index=False)
            
            # Summary statistics
            if len(df) > 0:
                summary_stats = {
                    'Total_Motifs': len(df),
                    'Unique_Classes': df['Class'].nunique(),
                    'Mean_Score': df['Normalized_Score'].mean(),
                    'Max_Score': df['Normalized_Score'].max(),
                    'Mean_Length': df['Length'].mean(),
                    'Total_Coverage': df['Length'].sum(),
                }
                summary_df = pd.DataFrame.from_dict(summary_stats, orient='index', columns=['Value'])
                summary_df.to_excel(writer, sheet_name='Summary')
            
            # Class-specific sheets
            for class_name in df['Class'].unique():
                if pd.notna(class_name):
                    class_df = df[df['Class'] == class_name]
                    sheet_name = class_name.replace(' ', '_')[:31]
                    class_df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        # GFF3
        gff3_file = f"{args.output}.gff3"
        gff3_lines = ["##gff-version 3"]
        for _, row in df.iterrows():
            attributes = f"ID=motif_{row['S.No']};Class={row['Class']};Subclass={row['Subclass']};Score={row['Normalized_Score']:.3f}"
            gff_line = f"{row['Chromosome/Contig']}\tNBDFinder\tmotif\t{row['Start']}\t{row['End']}\t{row['Normalized_Score']:.3f}\t.\t.\t{attributes}"
            gff3_lines.append(gff_line)
        
        with open(gff3_file, 'w') as f:
            f.write('\n'.join(gff3_lines))
        
        # Summary report
        summary_file = f"{args.output}_summary.txt"
        lines = [
            "NonBDNA Finder Results Summary",
            "=" * 30,
            "",
            f"Total motifs detected: {len(df)}",
            f"Unique sequences analyzed: {df['Sequence_Name'].nunique()}",
            f"Unique motif classes: {df['Class'].nunique()}",
            "",
            "Motif Classes Found:",
            "-" * 20
        ]
        
        class_counts = df['Class'].value_counts()
        for class_name, count in class_counts.items():
            lines.append(f"{class_name}: {count} motifs")
        
        lines.extend([
            "",
            "Score Statistics:",
            "-" * 16,
            f"Mean normalized score: {df['Normalized_Score'].mean():.3f}",
            f"Maximum score: {df['Normalized_Score'].max():.3f}",
            f"Minimum score: {df['Normalized_Score'].min():.3f}",
            "",
            "Length Statistics:",
            "-" * 17,
            f"Mean motif length: {df['Length'].mean():.1f} bp",
            f"Total coverage: {df['Length'].sum()} bp",
            f"Longest motif: {df['Length'].max()} bp",
        ])
        
        with open(summary_file, 'w') as f:
            f.write('\n'.join(lines))
        
        # Save zip package
        zip_file = f"{args.output}.zip"
        with open(zip_file, 'wb') as f:
            f.write(zip_data)
        
        # Report results
        if not args.quiet:
            print("✅ Analysis completed successfully!")
            print("")
            print("📊 RESULTS SUMMARY:")
            print(f"   Total motifs: {len(df)}")
            print(f"   Motif classes: {df['Class'].nunique()}")
            print(f"   Mean score: {df['Normalized_Score'].mean():.3f}")
            print("")
            print("📁 Output files created:")
            print(f"   📄 {csv_file}")
            print(f"   📊 {excel_file}")
            print(f"   🧬 {gff3_file}")
            print(f"   📋 {summary_file}")
            print(f"   📦 {zip_file}")
        
        if args.show_stats:
            print("")
            print("📈 DETAILED STATISTICS:")
            print("-" * 30)
            print("\nClass distribution:")
            for class_name, count in class_counts.items():
                print(f"  {class_name}: {count} motifs")
            
            print(f"\nScore statistics:")
            print(f"  Mean: {df['Normalized_Score'].mean():.3f}")
            print(f"  Std:  {df['Normalized_Score'].std():.3f}")
            print(f"  Min:  {df['Normalized_Score'].min():.3f}")
            print(f"  Max:  {df['Normalized_Score'].max():.3f}")
            
            print(f"\nLength statistics:")
            print(f"  Mean: {df['Length'].mean():.1f} bp")
            print(f"  Std:  {df['Length'].std():.1f} bp")
            print(f"  Min:  {df['Length'].min()} bp")
            print(f"  Max:  {df['Length'].max()} bp")
    
    except Exception as e:
        print(f"❌ Error during analysis: {e}")
        if not args.quiet:
            import traceback
            print("\nDebug information:")
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()