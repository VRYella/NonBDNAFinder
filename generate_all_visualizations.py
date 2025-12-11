#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║          COMPREHENSIVE VISUALIZATION GENERATOR WITH SCREENSHOTS              ║
║     Generate All Visualizations and Save Screenshots for Documentation      ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: generate_all_visualizations.py
VERSION: 2024.1
LICENSE: MIT

DESCRIPTION:
    Generates comprehensive test data with all detected motif classes and
    creates ALL available visualizations with screenshots for documentation.
    
OUTPUT:
    - PNG screenshots of all visualizations
    - Summary report of generated visualizations
    - Test data used for generation
"""

import sys
import os
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

# Import scanner and visualization functions
import nonbscanner as nbs
from visualizations import (
    plot_motif_distribution,
    plot_coverage_map,
    plot_density_heatmap,
    plot_score_distribution,
    plot_length_distribution,
    plot_nested_pie_chart,
    plot_class_analysis_comprehensive,
    plot_subclass_analysis_comprehensive,
    plot_score_statistics_by_class,
    plot_length_statistics_by_class,
    plot_manhattan_motif_density,
    plot_cumulative_motif_distribution,
    plot_motif_cooccurrence_matrix,
    plot_linear_motif_track,
    plot_cluster_size_distribution,
    plot_motif_length_kde,
    plot_genome_landscape_track,
    plot_sliding_window_heat_ribbon,
    plot_density_comparison_by_subclass,
    plot_enrichment_analysis_by_subclass,
    plot_subclass_density_heatmap,
    save_all_plots,
)
from utilities import calculate_genomic_density, calculate_positional_density


# ============================================================================
# TEST SEQUENCE GENERATION
# ============================================================================

def generate_comprehensive_test_sequence() -> str:
    """Generate a comprehensive test sequence with various motif types"""
    
    # G-Quadruplex canonical
    g4_sequence = "GGGTTAGGGTTAGGGTTAGGG"
    
    # Z-DNA
    zdna_sequence = "CGCGCGCGCGCGCGCG"
    
    # i-Motif
    imotif_sequence = "CCCTTACCCTTACCCTTACCC"
    
    # R-Loop (G-rich)
    rloop_sequence = "GGGGGGGGGGGGGGGGGGGGGGGGG"
    
    # A-tract curvature
    curved_sequence = "AAAAAAAAAAGCTAGCTAGCT"
    
    # STR (CAG repeat)
    str_sequence = "CAG" * 20
    
    # eGZ motif
    egz_sequence = "CGG" * 6
    
    # Spacer sequences
    spacer = "ATCGATCGATCG" * 5
    
    # Combine all motifs with spacers
    full_sequence = (
        spacer + g4_sequence + spacer +
        zdna_sequence + spacer +
        imotif_sequence + spacer +
        rloop_sequence + spacer +
        curved_sequence + spacer +
        str_sequence + spacer +
        egz_sequence + spacer +
        g4_sequence + spacer  # Repeat for more diversity
    )
    
    return full_sequence


# ============================================================================
# VISUALIZATION GENERATOR
# ============================================================================

class VisualizationGenerator:
    """Generate all visualizations and save screenshots"""
    
    def __init__(self, output_dir: str = "/tmp/visualizations"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.generated_plots = []
        self.sequence = None
        self.motifs = None
    
    def generate_test_data(self):
        """Generate comprehensive test data"""
        print("\n" + "=" * 80)
        print("GENERATING TEST DATA")
        print("=" * 80)
        
        self.sequence = generate_comprehensive_test_sequence()
        print(f"Generated sequence length: {len(self.sequence)} bp")
        
        # Analyze sequence
        self.motifs = nbs.analyze_sequence(self.sequence, "test_sequence")
        print(f"Detected motifs: {len(self.motifs)}")
        
        # Show motif breakdown
        from collections import Counter
        class_counts = Counter(m['Class'] for m in self.motifs)
        print("\nMotif class breakdown:")
        for class_name, count in sorted(class_counts.items()):
            print(f"  - {class_name}: {count}")
        
        subclass_counts = Counter(f"{m['Class']}:{m['Subclass']}" for m in self.motifs)
        print(f"\nTotal unique subclasses detected: {len(subclass_counts)}")
    
    def save_plot(self, fig, name: str, description: str):
        """Save a plot as PNG screenshot"""
        filepath = self.output_dir / f"{name}.png"
        
        try:
            fig.savefig(filepath, dpi=300, bbox_inches='tight', 
                       facecolor='white', edgecolor='none')
            plt.close(fig)
            
            self.generated_plots.append({
                'name': name,
                'filepath': str(filepath),
                'description': description
            })
            print(f"  ✅ {name}.png")
        except Exception as e:
            print(f"  ❌ {name}.png - Error: {e}")
            plt.close(fig)
    
    def generate_all_visualizations(self):
        """Generate all available visualizations"""
        print("\n" + "=" * 80)
        print("GENERATING VISUALIZATIONS")
        print("=" * 80)
        
        if not self.motifs:
            print("❌ No test data available. Run generate_test_data() first.")
            return
        
        seq_len = len(self.sequence)
        
        # 1. Distribution plots
        print("\n1. Distribution Plots:")
        fig = plot_motif_distribution(self.motifs, by='Class')
        self.save_plot(fig, "01_motif_distribution_by_class", 
                      "Bar chart showing count of each motif class")
        
        fig = plot_motif_distribution(self.motifs, by='Subclass')
        self.save_plot(fig, "02_motif_distribution_by_subclass",
                      "Bar chart showing count of each motif subclass")
        
        # 2. Coverage and positional plots
        print("\n2. Coverage and Positional Plots:")
        fig = plot_coverage_map(self.motifs, seq_len)
        self.save_plot(fig, "03_coverage_map",
                      "Linear track showing motif positions along sequence")
        
        fig = plot_density_heatmap(self.motifs, seq_len)
        self.save_plot(fig, "04_density_heatmap",
                      "Heatmap showing motif density in windows across sequence")
        
        # 3. Statistical plots
        print("\n3. Statistical Plots:")
        fig = plot_score_distribution(self.motifs, by_class=True)
        self.save_plot(fig, "05_score_distribution",
                      "Box plot showing score distributions by class")
        
        fig = plot_length_distribution(self.motifs, by_class=True)
        self.save_plot(fig, "06_length_distribution",
                      "Violin plot showing length distributions by class")
        
        # 4. Hierarchical plots
        print("\n4. Hierarchical Plots:")
        fig = plot_nested_pie_chart(self.motifs)
        self.save_plot(fig, "07_nested_donut_chart",
                      "Nested donut chart showing class-subclass hierarchy")
        
        # 5. Comprehensive analysis plots
        print("\n5. Comprehensive Analysis Plots:")
        fig = plot_class_analysis_comprehensive(self.motifs)
        self.save_plot(fig, "08_class_analysis_comprehensive",
                      "Comprehensive class analysis with detection status")
        
        fig = plot_subclass_analysis_comprehensive(self.motifs)
        self.save_plot(fig, "09_subclass_analysis_comprehensive",
                      "Comprehensive subclass analysis organized by class")
        
        fig = plot_score_statistics_by_class(self.motifs)
        self.save_plot(fig, "10_score_statistics_by_class",
                      "Detailed score statistics with violin plots and tables")
        
        fig = plot_length_statistics_by_class(self.motifs)
        self.save_plot(fig, "11_length_statistics_by_class",
                      "Detailed length statistics with histograms and tables")
        
        # 6. Advanced genome-wide visualizations
        print("\n6. Advanced Genome-Wide Visualizations:")
        fig = plot_manhattan_motif_density(self.motifs, seq_len)
        self.save_plot(fig, "12_manhattan_motif_density",
                      "Manhattan plot showing motif density hotspots")
        
        fig = plot_cumulative_motif_distribution(self.motifs, seq_len)
        self.save_plot(fig, "13_cumulative_distribution",
                      "Cumulative motif count across genomic positions")
        
        fig = plot_motif_cooccurrence_matrix(self.motifs)
        self.save_plot(fig, "14_cooccurrence_matrix",
                      "Heatmap showing co-occurrence between motif classes")
        
        fig = plot_linear_motif_track(self.motifs, seq_len)
        self.save_plot(fig, "15_linear_motif_track",
                      "Horizontal track with colored blocks for each motif")
        
        fig = plot_motif_length_kde(self.motifs, by_class=True)
        self.save_plot(fig, "16_length_kde",
                      "Kernel density estimation of motif length distributions")
        
        # 7. Specialized visualizations
        print("\n7. Specialized Visualizations:")
        fig = plot_genome_landscape_track(self.motifs, seq_len)
        self.save_plot(fig, "17_genome_landscape_track",
                      "Genome landscape with density and motif tracks")
        
        fig = plot_sliding_window_heat_ribbon(self.motifs, seq_len)
        self.save_plot(fig, "18_sliding_window_heat_ribbon",
                      "1D heatmap ribbon with density profile")
        
        # 8. Density and enrichment analysis
        print("\n8. Density and Enrichment Analysis:")
        
        # Calculate densities for subclass analysis
        genomic_density = calculate_genomic_density(self.motifs, seq_len, by_subclass=True)
        positional_density = calculate_positional_density(self.motifs, seq_len, by_subclass=True)
        
        fig = plot_density_comparison_by_subclass(genomic_density, positional_density)
        self.save_plot(fig, "19_density_comparison_subclass",
                      "Comparison of genomic and positional density by subclass")
        
        fig = plot_subclass_density_heatmap(self.motifs, seq_len)
        self.save_plot(fig, "20_subclass_density_heatmap",
                      "Heatmap showing density of each subclass across sequence")
        
        print(f"\n✅ Generated {len(self.generated_plots)} visualizations")
    
    def generate_summary_report(self) -> str:
        """Generate markdown summary report"""
        report = []
        report.append("# Visualization Test Report\n")
        report.append(f"\n**Generated:** {len(self.generated_plots)} visualizations\n")
        report.append(f"**Output Directory:** `{self.output_dir}`\n")
        
        if self.sequence and self.motifs:
            report.append(f"\n## Test Data Summary\n")
            report.append(f"- **Sequence Length:** {len(self.sequence)} bp\n")
            report.append(f"- **Total Motifs:** {len(self.motifs)}\n")
            
            from collections import Counter
            class_counts = Counter(m['Class'] for m in self.motifs)
            report.append(f"- **Unique Classes:** {len(class_counts)}\n")
            
            subclass_counts = Counter(f"{m['Class']}:{m['Subclass']}" for m in self.motifs)
            report.append(f"- **Unique Subclasses:** {len(subclass_counts)}\n")
        
        report.append("\n## Generated Visualizations\n")
        for i, plot_info in enumerate(self.generated_plots, 1):
            report.append(f"\n### {i}. {plot_info['name']}\n")
            report.append(f"**Description:** {plot_info['description']}\n")
            report.append(f"**File:** `{plot_info['filepath']}`\n")
            report.append(f"\n![{plot_info['name']}]({plot_info['filepath']})\n")
        
        return ''.join(report)


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Generate all visualizations with screenshots"""
    print("\n" + "=" * 80)
    print("COMPREHENSIVE VISUALIZATION GENERATOR")
    print("=" * 80)
    
    # Create generator
    generator = VisualizationGenerator()
    
    # Generate test data
    generator.generate_test_data()
    
    # Generate all visualizations
    generator.generate_all_visualizations()
    
    # Generate and save summary report
    report = generator.generate_summary_report()
    report_file = generator.output_dir / "visualization_report.md"
    with open(report_file, 'w') as f:
        f.write(report)
    
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"✅ Generated {len(generator.generated_plots)} visualizations")
    print(f"📁 Output directory: {generator.output_dir}")
    print(f"📄 Report saved to: {report_file}")
    print("=" * 80)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
