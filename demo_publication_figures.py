#!/usr/bin/env python3
"""
Quick demonstration of publication-ready subplot figures.

This script shows how to use the enhanced multi-panel figure functions
for creating publication-quality visualizations.

Usage:
    python demo_publication_figures.py
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

# Import necessary modules
import nonbscanner as nbs
from publication_figures import (
    create_figure_1_overview,
    create_figure_2_statistics,
    create_figure_3_landscape,
    create_figure_4_subclass_analysis,
    create_figure_5_quality_metrics,
    create_figure_6_comparative_analysis,
    create_supplementary_figure_s1,
    generate_example_figures,
)


def create_demo_sequence():
    """Create a demo sequence with multiple motif types."""
    print("Creating demo sequence...")
    
    # Build sequence with various motifs
    g4 = "GGGTTAGGGTTAGGGTTAGGG"
    zdna = "CGCGCGCGCGCGCGCG"
    imotif = "CCCTTACCCTTACCCTTACCC"
    rloop = "GGGGGGGGGGGGGGGGGGGGGGGGG"
    curved = "AAAAAAAAAAGCTAGCTAGCT"
    str_seq = "CAG" * 20
    spacer = "ATCGATCGATCG" * 5
    
    sequence = (
        spacer + g4 + spacer +
        zdna + spacer +
        imotif + spacer +
        rloop + spacer +
        curved + spacer +
        str_seq + spacer
    )
    
    print(f"✓ Created sequence: {len(sequence)} bp")
    return sequence


def demo_individual_figures(motifs, sequence_length, output_dir="/tmp/demo_figures"):
    """Demonstrate creating individual figures."""
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    print("\n" + "=" * 70)
    print("DEMO: Creating Individual Publication Figures")
    print("=" * 70)
    
    # Example 1: Create just Figure 1
    print("\n1. Creating Figure 1 (Comprehensive Overview)...")
    fig1 = create_figure_1_overview(
        motifs, 
        sequence_length,
        title="Non-B DNA Motif Analysis - Figure 1",
        save_path=f"{output_dir}/demo_figure1.pdf"
    )
    plt.close(fig1)
    print("   ✓ Saved to demo_figure1.pdf")
    
    # Example 2: Create Figure 5 (Quality Metrics)
    print("\n2. Creating Figure 5 (Quality Metrics)...")
    fig5 = create_figure_5_quality_metrics(
        motifs, 
        sequence_length,
        title="Motif Quality Analysis",
        save_path=f"{output_dir}/demo_figure5.pdf"
    )
    plt.close(fig5)
    print("   ✓ Saved to demo_figure5.pdf")
    
    # Example 3: Create Supplementary Figure
    print("\n3. Creating Supplementary Figure S1...")
    figS1 = create_supplementary_figure_s1(
        motifs, 
        sequence_length,
        save_path=f"{output_dir}/demo_figureS1.pdf"
    )
    plt.close(figS1)
    print("   ✓ Saved to demo_figureS1.pdf")
    
    print(f"\n✓ Individual figures saved to: {output_dir}")
    print("=" * 70)


def demo_batch_generation(motifs, sequence_length, output_dir="/tmp/demo_batch"):
    """Demonstrate batch generation of all figures."""
    print("\n" + "=" * 70)
    print("DEMO: Batch Generation of All Figures")
    print("=" * 70)
    
    # Generate all figures at once
    figures = generate_example_figures(
        motifs, 
        sequence_length, 
        output_dir=output_dir,
        include_supplementary=True
    )
    
    print(f"\n✓ Generated {len(figures)} figures")
    print(f"✓ All figures saved to: {output_dir}")
    print("=" * 70)


def main():
    """Main demo function."""
    print("\n" + "=" * 70)
    print("PUBLICATION-READY SUBPLOT FIGURES - DEMONSTRATION")
    print("=" * 70)
    
    # Create demo data
    sequence = create_demo_sequence()
    
    # Analyze sequence
    print("\nAnalyzing sequence for Non-B DNA motifs...")
    motifs = nbs.analyze_sequence(sequence, "demo_sequence")
    print(f"✓ Detected {len(motifs)} motifs")
    
    # Show motif breakdown
    from collections import Counter
    class_counts = Counter(m['Class'] for m in motifs)
    print(f"\nMotif classes detected ({len(class_counts)}):")
    for class_name, count in sorted(class_counts.items()):
        print(f"  • {class_name.replace('_', ' ')}: {count}")
    
    sequence_length = len(sequence)
    
    # Demo 1: Individual figures
    demo_individual_figures(motifs, sequence_length)
    
    # Demo 2: Batch generation
    demo_batch_generation(motifs, sequence_length)
    
    print("\n" + "=" * 70)
    print("DEMO COMPLETE!")
    print("=" * 70)
    print("\nKey Features Demonstrated:")
    print("  ✓ Multi-panel subplot layouts (2x2, 3x1, 3x2)")
    print("  ✓ Publication-quality styling (Nature/Science standards)")
    print("  ✓ Proper panel labeling (A, B, C, D...)")
    print("  ✓ Comprehensive scientific storytelling")
    print("  ✓ 300 DPI output for print publication")
    print("  ✓ Colorblind-safe color palettes")
    print("  ✓ Individual and batch generation modes")
    print("\nFigure Types Created:")
    print("  • Figure 1: Comprehensive Overview (2x2 panels)")
    print("  • Figure 2: Statistical Characterization (2x2 panels)")
    print("  • Figure 3: Genomic Landscape (3x1 panels)")
    print("  • Figure 4: Subclass Analysis (2x2 panels)")
    print("  • Figure 5: Quality Metrics (2x2 panels) [NEW]")
    print("  • Figure 6: Comparative Analysis (2x2 panels) [NEW]")
    print("  • Figure S1: Comprehensive Supplementary (3x2 panels) [NEW]")
    print("\nNext Steps:")
    print("  1. Review generated figures in /tmp/demo_figures/")
    print("  2. Import publication_figures module in your own scripts")
    print("  3. Customize titles, colors, and layouts as needed")
    print("  4. Export to PDF/SVG for manuscript submission")
    print("=" * 70)


if __name__ == '__main__':
    main()
