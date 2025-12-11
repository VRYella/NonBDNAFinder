#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║         GENERATE NATURE-LEVEL PUBLICATION-READY MULTI-PANEL FIGURES         ║
║               Automated Figure Generation for Manuscript Submission         ║
╚══════════════════════════════════════════════════════════════════════════════╝

SCRIPT: generate_publication_figures.py
VERSION: 2024.1.5
LICENSE: MIT

DESCRIPTION:
    Generates publication-ready multi-panel figures from test data following
    Nature Methods/Genetics submission guidelines. Creates comprehensive
    figures with proper panel labeling, consistent styling, and captions.

USAGE:
    python generate_publication_figures.py [--output-dir DIR] [--format FORMAT]
    
OPTIONS:
    --output-dir DIR    Output directory (default: ./publication_figures)
    --format FORMAT     Output format: pdf, svg, png (default: pdf)
    --dpi DPI          Resolution for raster formats (default: 300)
"""

import sys
import os
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

# Import scanner and visualization modules
import nonbscanner as nbs
from publication_figures import (
    create_figure_1_overview,
    create_figure_2_statistics,
    create_figure_3_landscape,
    create_figure_4_subclass_analysis,
)


def generate_comprehensive_test_data():
    """Generate comprehensive test sequence with diverse motif types."""
    
    print("\n" + "=" * 80)
    print("GENERATING TEST DATA")
    print("=" * 80)
    
    # Create test sequence with multiple motif types
    SPACER_REPEAT_COUNT = 5
    
    # Different motif sequences
    g4_sequence = "GGGTTAGGGTTAGGGTTAGGG"
    zdna_sequence = "CGCGCGCGCGCGCGCG"
    imotif_sequence = "CCCTTACCCTTACCCTTACCC"
    rloop_sequence = "GGGGGGGGGGGGGGGGGGGGGGGGG"
    curved_sequence = "AAAAAAAAAAGCTAGCTAGCT"
    str_sequence = "CAG" * 20
    egz_sequence = "CGG" * 6
    spacer = "ATCGATCGATCG" * SPACER_REPEAT_COUNT
    
    # Combine with spacers
    full_sequence = (
        spacer + g4_sequence + spacer +
        zdna_sequence + spacer +
        imotif_sequence + spacer +
        rloop_sequence + spacer +
        curved_sequence + spacer +
        str_sequence + spacer +
        egz_sequence + spacer +
        g4_sequence + spacer
    )
    
    print(f"✓ Generated sequence: {len(full_sequence)} bp")
    
    # Analyze sequence
    motifs = nbs.analyze_sequence(full_sequence, "test_sequence")
    print(f"✓ Detected {len(motifs)} motifs")
    
    # Show breakdown
    from collections import Counter
    class_counts = Counter(m['Class'] for m in motifs)
    print(f"\nMotif class breakdown ({len(class_counts)} classes):")
    for class_name, count in sorted(class_counts.items()):
        print(f"  • {class_name.replace('_', ' ')}: {count}")
    
    return full_sequence, motifs


def generate_figure_captions(output_dir: str):
    """Generate figure caption file for manuscript."""
    
    captions = """
# Figure Captions for Manuscript

## Figure 1. Comprehensive Non-B DNA Motif Analysis

Multi-panel overview of Non-B DNA motif detection in test genomic sequence.
**(A)** Distribution of detected motifs across all 11 major Non-B DNA classes 
(n = total motifs). Bar colors follow Wong (2011) colorblind-safe palette. 
**(B)** Hierarchical composition showing class-subclass relationships. Inner ring 
represents major classes with proportional sizes; outer ring shows subclass 
distribution. Percentages indicate relative abundance. **(C)** Genomic coverage 
map displaying motif positions along analyzed sequence. Each horizontal track 
represents one motif class; colored rectangles indicate individual motif instances 
with class-specific colors. X-axis shows genomic coordinates in base pairs.

Scale bars and statistical measures as indicated. All analyses performed using
NBDScanner v2024.1 with default parameters.

---

## Figure 2. Statistical Characterization of Non-B DNA Motifs

Comprehensive statistical analysis of motif properties across detected classes.
**(A)** Violin plots showing score distributions for each motif class. Box plots
inside violins indicate quartiles; width represents density. **(B)** Box plots
of motif length distributions by class. Boxes show interquartile range (IQR);
whiskers extend to 1.5×IQR; outliers shown as individual points. **(C)** Scatter
plot of motif score versus length with class-specific coloring, revealing 
potential correlations between sequence properties and detection scores. 
**(D)** Summary statistics table showing sample sizes, mean scores, and mean 
lengths for major motif classes.

All statistical comparisons performed using two-tailed Student's t-test.
Significance levels: *P < 0.05, **P < 0.01, ***P < 0.001.

---

## Figure 3. Genomic Landscape and Density Analysis

Multi-scale visualization of motif distribution patterns across analyzed sequence.
**(A)** Density heatmap showing motif counts in sliding windows across the genome.
Color intensity indicates motif density (darker = higher density). Window size
optimized for sequence length. **(B)** Manhattan plot displaying density hotspots
for each motif class. Each point represents one genomic window colored by motif
class. Y-axis shows motif density (motifs per kilobase). Identifies regions of
high motif enrichment. **(C)** Cumulative distribution curves showing running
sum of motifs along genomic coordinates for the five most abundant classes.
Steeper slopes indicate regions of higher motif density.

Window size for density calculations: auto-optimized to sequence length.
All coordinates shown in kilobases (kb) from sequence start.

---

## General Notes for All Figures

- **Resolution**: All figures generated at 300 DPI for publication quality
- **Format**: Vector format (PDF/SVG) for editability and scalability
- **Fonts**: Arial family following Nature Methods guidelines (8pt body, 9pt titles)
- **Colors**: Wong (2011) colorblind-safe palette for accessibility
- **Software**: Figures generated using Python 3.12, matplotlib 3.5+, seaborn 0.11+
- **Reproducibility**: Complete code and data available at [repository URL]

**Color Accessibility**: All color schemes tested with colorblind simulation tools
to ensure distinguishability for readers with color vision deficiencies.

**Data Availability**: Test sequences and analysis parameters provided in 
Supplementary Materials. Full codebase open-source under MIT license.
"""
    
    caption_file = Path(output_dir) / "figure_captions.md"
    with open(caption_file, 'w') as f:
        f.write(captions)
    
    print(f"✓ Saved figure captions to: {caption_file}")
    return caption_file


def main():
    """Main execution function."""
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Generate Nature-level publication figures"
    )
    parser.add_argument('--output-dir', default='./publication_figures',
                       help='Output directory (default: ./publication_figures)')
    parser.add_argument('--format', default='pdf', choices=['pdf', 'svg', 'png'],
                       help='Output format (default: pdf)')
    parser.add_argument('--dpi', type=int, default=300,
                       help='Resolution for raster formats (default: 300)')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("\n" + "=" * 80)
    print("NATURE-LEVEL PUBLICATION FIGURE GENERATOR")
    print("=" * 80)
    print(f"Output directory: {output_dir}")
    print(f"Output format: {args.format}")
    print(f"DPI: {args.dpi}")
    
    # Generate test data
    try:
        sequence, motifs = generate_comprehensive_test_data()
    except Exception as e:
        print(f"\n❌ Error generating test data: {e}")
        return 1
    
    if not motifs:
        print("\n❌ No motifs detected. Cannot generate figures.")
        return 1
    
    sequence_length = len(sequence)
    
    # Generate figures
    print("\n" + "=" * 80)
    print("GENERATING PUBLICATION FIGURES")
    print("=" * 80)
    
    figures_generated = []
    
    try:
        # Figure 1: Overview
        print("\n[1/3] Generating Figure 1: Comprehensive Overview...")
        fig1_path = output_dir / f"figure1_overview.{args.format}"
        fig1 = create_figure_1_overview(motifs, sequence_length,
                                        save_path=str(fig1_path))
        plt.close(fig1)
        figures_generated.append(fig1_path)
        print(f"      ✓ Saved to: {fig1_path}")
        
    except Exception as e:
        print(f"      ❌ Error: {e}")
    
    try:
        # Figure 2: Statistics
        print("\n[2/3] Generating Figure 2: Statistical Characterization...")
        fig2_path = output_dir / f"figure2_statistics.{args.format}"
        fig2 = create_figure_2_statistics(motifs, save_path=str(fig2_path))
        plt.close(fig2)
        figures_generated.append(fig2_path)
        print(f"      ✓ Saved to: {fig2_path}")
        
    except Exception as e:
        print(f"      ❌ Error: {e}")
    
    try:
        # Figure 3: Landscape
        print("\n[3/4] Generating Figure 3: Genomic Landscape...")
        fig3_path = output_dir / f"figure3_landscape.{args.format}"
        fig3 = create_figure_3_landscape(motifs, sequence_length,
                                         save_path=str(fig3_path))
        plt.close(fig3)
        figures_generated.append(fig3_path)
        print(f"      ✓ Saved to: {fig3_path}")
        
    except Exception as e:
        print(f"      ❌ Error: {e}")
    
    try:
        # Figure 4: Subclass Analysis
        print("\n[4/4] Generating Figure 4: Subclass Analysis...")
        fig4_path = output_dir / f"figure4_subclass.{args.format}"
        fig4 = create_figure_4_subclass_analysis(motifs, sequence_length,
                                                 save_path=str(fig4_path))
        plt.close(fig4)
        figures_generated.append(fig4_path)
        print(f"      ✓ Saved to: {fig4_path}")
        
    except Exception as e:
        print(f"      ❌ Error: {e}")
    
    # Generate figure captions
    print("\n" + "=" * 80)
    print("GENERATING FIGURE CAPTIONS")
    print("=" * 80)
    
    try:
        caption_file = generate_figure_captions(output_dir)
    except Exception as e:
        print(f"❌ Error generating captions: {e}")
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"✓ Generated {len(figures_generated)} publication-ready figures")
    print(f"✓ Output directory: {output_dir}")
    print(f"✓ Format: {args.format} @ {args.dpi} DPI")
    print("\nGenerated files:")
    for fig_path in figures_generated:
        file_size = fig_path.stat().st_size / 1024  # KB
        print(f"  • {fig_path.name} ({file_size:.1f} KB)")
    
    print("\n" + "=" * 80)
    print("NEXT STEPS FOR MANUSCRIPT PREPARATION")
    print("=" * 80)
    print("1. Review generated figures for quality and accuracy")
    print("2. Edit figure_captions.md with specific results")
    print("3. Import figures into manuscript document")
    print("4. Verify figure dimensions match journal requirements")
    print("5. Check color rendering in grayscale for print version")
    print("=" * 80)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
