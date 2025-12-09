"""
Test script for new Nature-quality visualizations.
Generates sample outputs with test datasets to demonstrate the new plotting functions.
"""

import matplotlib.pyplot as plt
import numpy as np
from visualizations import (
    plot_manhattan_motif_density,
    plot_cumulative_motif_distribution,
    plot_motif_cooccurrence_matrix,
    plot_gc_content_correlation,
    plot_linear_motif_track,
    plot_cluster_size_distribution,
    plot_motif_length_kde,
    set_scientific_style
)

# Set scientific style globally
set_scientific_style('nature')

def generate_test_motifs(sequence_length=100000, n_motifs=500):
    """
    Generate realistic test motif data for visualization testing.
    
    Creates a synthetic dataset with a mix of regular motifs and cluster motifs
    to simulate real Non-B DNA detection results.
    
    Args:
        sequence_length (int): Length of the sequence in bp (default: 100,000)
        n_motifs (int): Number of regular motifs to generate (default: 500)
        
    Returns:
        list: List of motif dictionaries with Class, Subclass, Start, End, Length, Score keys
        
    Note:
        Also generates 20 cluster motifs with Motif_Count and Class_Diversity attributes.
    """
    motif_classes = [
        'G-Quadruplex', 'Z-DNA', 'Curved_DNA', 'i-Motif', 
        'Triplex', 'R-Loop', 'Cruciform', 'Slipped_DNA'
    ]
    
    motifs = []
    
    # Generate regular motifs
    for i in range(n_motifs):
        motif_class = np.random.choice(motif_classes)
        length = np.random.randint(15, 100)
        start = np.random.randint(1, sequence_length - length)
        end = start + length
        score = np.random.uniform(0.5, 1.0)
        
        motifs.append({
            'Class': motif_class,
            'Subclass': f'{motif_class}_subtype',
            'Start': start,
            'End': end,
            'Length': length,
            'Score': score
        })
    
    # Generate some cluster motifs
    for i in range(20):
        start = np.random.randint(1, sequence_length - 500)
        end = start + np.random.randint(200, 500)
        motif_count = np.random.randint(3, 10)
        diversity = np.random.randint(2, 5)
        
        motifs.append({
            'Class': 'Non-B_DNA_Clusters',
            'Subclass': 'High-density cluster',
            'Start': start,
            'End': end,
            'Length': end - start,
            'Score': 0.8,
            'Motif_Count': motif_count,
            'Class_Diversity': diversity
        })
    
    return motifs

def generate_test_sequence(length=100000):
    """
    Generate a random DNA sequence with realistic GC content variation.
    
    Creates a synthetic DNA sequence with oscillating GC content to simulate
    natural genomic variation and test GC correlation analysis.
    
    Args:
        length (int): Length of the sequence in bp (default: 100,000)
        
    Returns:
        str: DNA sequence string containing A, T, G, C nucleotides
        
    Note:
        GC content varies sinusoidally from 40% to 60% across the sequence
        to create realistic local variation patterns.
    """
    sequence = []
    for i in range(0, length, 1000):
        # Vary GC content in windows
        gc_pct = 0.4 + 0.2 * np.sin(i / 5000)  # Oscillating GC content
        gc_bases = int(1000 * gc_pct)
        at_bases = 1000 - gc_bases
        
        window_seq = ['G', 'C'] * (gc_bases // 2) + ['A', 'T'] * (at_bases // 2)
        np.random.shuffle(window_seq)
        sequence.extend(window_seq)
    
    return ''.join(sequence[:length])

def test_all_visualizations():
    """Test all new visualization functions and save outputs"""
    print("🧬 Testing Nature-Quality Visualizations for Non-B DNA Finder")
    print("=" * 70)
    
    # Generate test data
    print("\n📊 Generating test data...")
    sequence_length = 100000  # 100kb test sequence
    motifs = generate_test_motifs(sequence_length, n_motifs=500)
    sequence = generate_test_sequence(sequence_length)
    print(f"✓ Generated {len(motifs)} test motifs in {sequence_length/1000:.1f}kb sequence")
    
    # Test 1: Manhattan Plot
    print("\n1️⃣  Testing Manhattan Plot...")
    try:
        fig = plot_manhattan_motif_density(
            motifs, sequence_length,
            title="Manhattan Plot - Test Dataset"
        )
        fig.savefig('/tmp/test_manhattan.png', dpi=300, bbox_inches='tight')
        plt.close(fig)
        print("   ✓ Manhattan plot saved to /tmp/test_manhattan.png")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    # Test 2: Cumulative Distribution
    print("\n2️⃣  Testing Cumulative Distribution...")
    try:
        fig = plot_cumulative_motif_distribution(
            motifs, sequence_length,
            title="Cumulative Distribution - Test Dataset",
            by_class=True
        )
        fig.savefig('/tmp/test_cumulative.png', dpi=300, bbox_inches='tight')
        plt.close(fig)
        print("   ✓ Cumulative plot saved to /tmp/test_cumulative.png")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    # Test 3: Co-occurrence Matrix
    print("\n3️⃣  Testing Co-occurrence Matrix...")
    try:
        fig = plot_motif_cooccurrence_matrix(
            motifs,
            title="Motif Co-occurrence Matrix - Test Dataset"
        )
        fig.savefig('/tmp/test_cooccurrence.png', dpi=300, bbox_inches='tight')
        plt.close(fig)
        print("   ✓ Co-occurrence matrix saved to /tmp/test_cooccurrence.png")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    # Test 4: GC Content Correlation
    print("\n4️⃣  Testing GC Content Correlation...")
    try:
        fig = plot_gc_content_correlation(
            motifs, sequence,
            title="GC Content vs Motif Density - Test Dataset"
        )
        fig.savefig('/tmp/test_gc_correlation.png', dpi=300, bbox_inches='tight')
        plt.close(fig)
        print("   ✓ GC correlation plot saved to /tmp/test_gc_correlation.png")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    # Test 5: Linear Motif Track (small region)
    print("\n5️⃣  Testing Linear Motif Track...")
    try:
        fig = plot_linear_motif_track(
            motifs, sequence_length,
            region_start=0, region_end=10000,
            title="Linear Motif Track - Test Dataset (0-10kb)"
        )
        fig.savefig('/tmp/test_linear_track.png', dpi=300, bbox_inches='tight')
        plt.close(fig)
        print("   ✓ Linear track saved to /tmp/test_linear_track.png")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    # Test 6: Cluster Size Distribution
    print("\n6️⃣  Testing Cluster Size Distribution...")
    try:
        fig = plot_cluster_size_distribution(
            motifs,
            title="Cluster Statistics - Test Dataset"
        )
        fig.savefig('/tmp/test_cluster_dist.png', dpi=300, bbox_inches='tight')
        plt.close(fig)
        print("   ✓ Cluster distribution saved to /tmp/test_cluster_dist.png")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    # Test 7: Motif Length KDE
    print("\n7️⃣  Testing Motif Length KDE...")
    try:
        fig = plot_motif_length_kde(
            motifs,
            by_class=True,
            title="Length Distribution (KDE) - Test Dataset"
        )
        fig.savefig('/tmp/test_length_kde.png', dpi=300, bbox_inches='tight')
        plt.close(fig)
        print("   ✓ Length KDE saved to /tmp/test_length_kde.png")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    print("\n" + "=" * 70)
    print("✅ All visualization tests complete!")
    print("\n📁 Output files saved to /tmp/:")
    print("   - test_manhattan.png")
    print("   - test_cumulative.png")
    print("   - test_cooccurrence.png")
    print("   - test_gc_correlation.png")
    print("   - test_linear_track.png")
    print("   - test_cluster_dist.png")
    print("   - test_length_kde.png")
    print("\n🎨 All plots are publication-quality (300 DPI) following Nature Methods guidelines.")

if __name__ == "__main__":
    test_all_visualizations()
