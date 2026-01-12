#!/usr/bin/env python3
"""
Publication Narrative Generator
================================

Generates publication-quality narrative text for manuscripts based on
comparative genome analysis results. Creates four main sections suitable
for scientific journals.

Sections:
1. Materials and Methods (computational approach)
2. Results - Distribution Analysis
3. Results - Comparative Genomics
4. Discussion and Conclusions

Usage:
    python publication_narrative.py [--comparative-dir DIR] [--output OUTPUT]
"""

import os
import sys
import json
import argparse
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, List

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


class NarrativeGenerator:
    """Generator for publication-quality narrative text"""
    
    def __init__(self, comparative_results: Dict[str, Any]):
        self.results = comparative_results
        self.summary = comparative_results['summary']
        self.stats = comparative_results['genome_statistics']
        self.enrichment = comparative_results['enrichment_patterns']
        self.correlations = comparative_results['correlations']
    
    def generate_materials_methods(self) -> str:
        """Generate Materials and Methods section"""
        text = []
        
        text.append("MATERIALS AND METHODS")
        text.append("=" * 80)
        text.append("")
        text.append("Genome Data and Processing")
        text.append("-" * 80)
        text.append(f"We analyzed {self.summary['total_genomes']} complete genome sequences")
        text.append(f"totaling {self.summary['total_sequence_length']:,} base pairs. All sequences")
        text.append("were obtained in FASTA format and processed using the NonBDNAFinder")
        text.append("computational pipeline (version 2025.1).")
        text.append("")
        
        text.append("Non-B DNA Motif Detection")
        text.append("-" * 80)
        text.append("Non-B DNA motif detection was performed using NonBDNAFinder, a comprehensive")
        text.append("computational tool that identifies 11 major structural classes including")
        text.append("G-quadruplexes, Z-DNA, i-motifs, curved DNA, A-philic DNA, cruciform structures,")
        text.append("triplex DNA, R-loops, slipped DNA structures, hybrid motifs, and non-B DNA")
        text.append("clusters. The detection pipeline employs pattern-matching algorithms combined")
        text.append("with thermodynamic scoring models to ensure high sensitivity and specificity.")
        text.append("")
        
        text.append("Each detected motif was assigned a confidence score on a 1-3 scale based on")
        text.append("structural propensity calculations. Motifs were classified into specific")
        text.append("subclasses based on structural features such as loop length, tract arrangement,")
        text.append("and thermodynamic stability. Overlapping motifs within the same subclass were")
        text.append("merged to avoid redundancy, while overlaps between different classes were")
        text.append("flagged as hybrid structures for further analysis.")
        text.append("")
        
        text.append("Statistical Analysis")
        text.append("-" * 80)
        text.append("For each genome, we calculated motif density (motifs per kilobase), genomic")
        text.append("coverage (percentage of sequence occupied by motifs), and class diversity")
        text.append("(Shannon entropy). Cross-genome comparisons were performed using standardized")
        text.append("frequency normalization. Enrichment analysis identified motif classes showing")
        text.append("significantly elevated or depleted frequencies (Z-score > 2.0 or < -2.0)")
        text.append("relative to the cross-genome mean. Pearson correlation analysis was used to")
        text.append("assess relationships between genome features (length, GC content, motif density).")
        text.append("")
        
        return "\n".join(text)
    
    def generate_distribution_results(self) -> str:
        """Generate Results - Distribution Analysis section"""
        text = []
        
        text.append("RESULTS: DISTRIBUTION ANALYSIS")
        text.append("=" * 80)
        text.append("")
        
        # Overall statistics
        text.append("Global Non-B DNA Landscape")
        text.append("-" * 80)
        text.append(f"Across all {self.summary['total_genomes']} genomes analyzed, we detected")
        text.append(f"{self.summary['total_motifs_detected']:,} non-B DNA motifs spanning")
        text.append(f"{self.summary['total_sequence_length']:,} bp of sequence. The average")
        text.append(f"motif density was {self.summary['average_density_per_kb']:.2f} motifs per")
        text.append(f"kilobase, with genomic coverage averaging {self.summary['average_coverage_percent']:.2f}%.")
        text.append(f"Class diversity, measured by Shannon entropy, averaged")
        text.append(f"{self.summary['average_class_diversity']:.3f}, indicating moderate structural")
        text.append("heterogeneity across the analyzed genomes.")
        text.append("")
        
        # Genome-specific patterns
        text.append("Genome-Specific Patterns")
        text.append("-" * 80)
        
        # Highest and lowest density
        highest = self.summary['highest_density_genome']
        lowest = self.summary['lowest_density_genome']
        
        text.append(f"Motif density varied substantially across genomes. {highest['name']}")
        text.append(f"exhibited the highest density ({highest['density']:.2f} motifs/kb), while")
        text.append(f"{lowest['name']} showed the lowest density ({lowest['density']:.2f} motifs/kb),")
        text.append(f"representing a {highest['density']/lowest['density']:.1f}-fold difference.")
        text.append("")
        
        # Diversity patterns
        highest_div = self.summary['highest_diversity_genome']
        lowest_div = self.summary['lowest_diversity_genome']
        
        text.append(f"Class diversity also showed marked variation. {highest_div['name']}")
        text.append(f"displayed the highest diversity (H={highest_div['diversity']:.3f}), suggesting")
        text.append(f"a broad distribution of structural classes. In contrast, {lowest_div['name']}")
        text.append(f"showed the lowest diversity (H={lowest_div['diversity']:.3f}), indicating")
        text.append("a more homogeneous non-B DNA landscape dominated by fewer motif types.")
        text.append("")
        
        # Class distribution summary
        text.append("Motif Class Distribution")
        text.append("-" * 80)
        
        # Collect all class counts across genomes
        all_class_counts = {}
        for genome_name, stats in self.stats.items():
            for class_name, count in stats['class_counts'].items():
                all_class_counts[class_name] = all_class_counts.get(class_name, 0) + count
        
        # Sort by frequency
        sorted_classes = sorted(all_class_counts.items(), key=lambda x: x[1], reverse=True)
        
        text.append("The most frequently detected motif classes across all genomes were:")
        for i, (class_name, count) in enumerate(sorted_classes[:5], 1):
            percentage = (count / self.summary['total_motifs_detected']) * 100
            text.append(f"  {i}. {class_name}: {count:,} motifs ({percentage:.1f}%)")
        
        text.append("")
        
        return "\n".join(text)
    
    def generate_comparative_results(self) -> str:
        """Generate Results - Comparative Genomics section"""
        text = []
        
        text.append("RESULTS: COMPARATIVE GENOMICS")
        text.append("=" * 80)
        text.append("")
        
        # Enrichment patterns
        text.append("Motif Class Enrichment Patterns")
        text.append("-" * 80)
        
        # Find classes with significant enrichment/depletion
        significant_enrichment = []
        for class_name, enrichment in self.enrichment.items():
            if enrichment['enriched_genomes'] or enrichment['depleted_genomes']:
                significant_enrichment.append((class_name, enrichment))
        
        if significant_enrichment:
            text.append("Several motif classes showed significant enrichment or depletion in")
            text.append("specific genomes (Z-score > |2.0|):")
            text.append("")
            
            for class_name, enrichment in significant_enrichment[:5]:  # Top 5
                if enrichment['enriched_genomes']:
                    enriched_list = ', '.join([g['genome'] for g in enrichment['enriched_genomes']])
                    text.append(f"• {class_name}: Enriched in {enriched_list}")
                
                if enrichment['depleted_genomes']:
                    depleted_list = ', '.join([g['genome'] for g in enrichment['depleted_genomes']])
                    text.append(f"• {class_name}: Depleted in {depleted_list}")
            
            text.append("")
        else:
            text.append("No motif classes showed significant enrichment or depletion patterns")
            text.append("across the analyzed genomes, suggesting relatively uniform distribution.")
            text.append("")
        
        # Coefficient of variation analysis
        text.append("Variability Analysis")
        text.append("-" * 80)
        
        # Find classes with highest and lowest variation
        cv_values = [(cn, e['coefficient_variation']) for cn, e in self.enrichment.items()]
        cv_values.sort(key=lambda x: x[1], reverse=True)
        
        if cv_values:
            highest_cv = cv_values[0]
            lowest_cv = cv_values[-1]
            
            text.append(f"Cross-genome variability, measured by coefficient of variation (CV),")
            text.append(f"revealed distinct patterns. {highest_cv[0]} showed the highest variability")
            text.append(f"(CV={highest_cv[1]:.3f}), indicating highly genome-specific distribution.")
            text.append(f"Conversely, {lowest_cv[0]} displayed the lowest variability (CV={lowest_cv[1]:.3f}),")
            text.append("suggesting more uniform presence across genomes.")
            text.append("")
        
        # Correlations
        if self.correlations:
            text.append("Feature Correlations")
            text.append("-" * 80)
            
            # Find significant correlations
            sig_correlations = [(name, data) for name, data in self.correlations.items() 
                               if data['significant']]
            
            if sig_correlations:
                text.append("Several genome features showed significant correlations (p < 0.05):")
                text.append("")
                
                for corr_name, corr_data in sig_correlations[:5]:  # Top 5
                    feat1, feat2 = corr_name.replace('_vs_', ' vs ').replace('_', ' ').split(' vs ')
                    corr = corr_data['correlation']
                    pval = corr_data['p_value']
                    
                    direction = "positive" if corr > 0 else "negative"
                    strength = "strong" if abs(corr) > 0.7 else "moderate" if abs(corr) > 0.4 else "weak"
                    
                    text.append(f"• {feat1.title()} vs {feat2.title()}: {strength} {direction}")
                    text.append(f"  correlation (r={corr:.3f}, p={pval:.4f})")
                
                text.append("")
            else:
                text.append("No significant correlations were detected between genome features,")
                text.append("suggesting independent variation of motif density, diversity, and coverage.")
                text.append("")
        
        return "\n".join(text)
    
    def generate_discussion(self) -> str:
        """Generate Discussion and Conclusions section"""
        text = []
        
        text.append("DISCUSSION AND CONCLUSIONS")
        text.append("=" * 80)
        text.append("")
        
        text.append("Summary of Findings")
        text.append("-" * 80)
        text.append(f"Our comprehensive analysis of {self.summary['total_genomes']} genomes revealed")
        text.append("substantial variation in non-B DNA motif distribution, density, and class")
        text.append("composition. These findings highlight the genome-specific nature of non-B DNA")
        text.append("landscapes and their potential functional implications.")
        text.append("")
        
        text.append("Biological Significance")
        text.append("-" * 80)
        text.append("The observed differences in motif density and diversity likely reflect")
        text.append("evolutionary pressures, genomic architecture, and functional requirements")
        text.append("specific to each organism. High-density regions may correspond to regulatory")
        text.append("hotspots, chromatin organization domains, or recombination-prone sequences.")
        text.append("The enrichment of specific motif classes in certain genomes suggests")
        text.append("organism-specific structural requirements or functional constraints.")
        text.append("")
        
        text.append("Methodological Considerations")
        text.append("-" * 80)
        text.append("This analysis employed computational prediction methods that, while highly")
        text.append("sensitive and specific, represent potential structural propensities rather")
        text.append("than experimentally confirmed structures. In vivo formation depends on")
        text.append("cellular context including chromatin state, protein binding, and DNA")
        text.append("supercoiling. Nevertheless, computational approaches enable genome-scale")
        text.append("analyses infeasible through experimental methods alone.")
        text.append("")
        
        text.append("Future Directions")
        text.append("-" * 80)
        text.append("Future studies should integrate these structural predictions with")
        text.append("experimental data (ChIP-seq, DNA structure mapping) to validate functional")
        text.append("relevance. Cross-species comparative analyses could reveal evolutionary")
        text.append("conservation patterns and identify functionally constrained non-B DNA motifs.")
        text.append("Additionally, correlating motif distributions with gene expression, chromatin")
        text.append("accessibility, and replication timing data would provide mechanistic insights")
        text.append("into the biological roles of these alternative DNA structures.")
        text.append("")
        
        text.append("Conclusions")
        text.append("-" * 80)
        text.append(f"This study presents a comprehensive catalog of non-B DNA motifs across")
        text.append(f"{self.summary['total_genomes']} genomes, identifying {self.summary['total_motifs_detected']:,}")
        text.append("potential structural elements. The substantial genome-to-genome variation")
        text.append("in motif density, class composition, and enrichment patterns underscores")
        text.append("the biological significance of these alternative DNA structures and their")
        text.append("potential roles in genome function, regulation, and evolution.")
        text.append("")
        
        return "\n".join(text)
    
    def generate_complete_narrative(self) -> str:
        """Generate complete publication narrative"""
        sections = []
        
        # Title
        sections.append("=" * 80)
        sections.append("PUBLICATION NARRATIVE: COMPARATIVE NON-B DNA ANALYSIS")
        sections.append("=" * 80)
        sections.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        sections.append("=" * 80)
        sections.append("")
        
        # Add all sections
        sections.append(self.generate_materials_methods())
        sections.append("\n")
        sections.append(self.generate_distribution_results())
        sections.append("\n")
        sections.append(self.generate_comparative_results())
        sections.append("\n")
        sections.append(self.generate_discussion())
        
        # Footer
        sections.append("")
        sections.append("=" * 80)
        sections.append("END OF NARRATIVE")
        sections.append("=" * 80)
        
        return "\n".join(sections)


def main():
    """Main function for narrative generation"""
    parser = argparse.ArgumentParser(
        description='Generate publication-quality narrative text'
    )
    parser.add_argument(
        '--comparative-dir',
        default='comparative_results',
        help='Comparative results directory (default: comparative_results)'
    )
    parser.add_argument(
        '--output',
        default='publication_narrative.txt',
        help='Output text file (default: publication_narrative.txt)'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("PUBLICATION NARRATIVE GENERATOR")
    print("="*80)
    print(f"Input directory: {args.comparative_dir}")
    print(f"Output file: {args.output}\n")
    
    # Load comparative results
    print("Loading comparative results...")
    comparative_file = os.path.join(args.comparative_dir, 'comparative_analysis.json')
    
    if not os.path.exists(comparative_file):
        print(f"ERROR: Comparative results not found: {comparative_file}")
        print("Please run comparative_analysis.py first.")
        return 1
    
    with open(comparative_file, 'r') as f:
        comparative_results = json.load(f)
    
    # Generate narrative
    print("Generating narrative sections...")
    generator = NarrativeGenerator(comparative_results)
    narrative = generator.generate_complete_narrative()
    
    # Save narrative
    print(f"Saving narrative: {args.output}")
    with open(args.output, 'w') as f:
        f.write(narrative)
    
    print("\n" + "="*80)
    print("NARRATIVE GENERATION COMPLETE")
    print("="*80)
    print(f"Output file: {args.output}")
    print("\nGenerated sections:")
    print("  • Materials and Methods (computational approach)")
    print("  • Results - Distribution Analysis")
    print("  • Results - Comparative Genomics")
    print("  • Discussion and Conclusions")
    print("="*80)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
