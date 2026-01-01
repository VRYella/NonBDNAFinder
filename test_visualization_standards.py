"""
Test suite for Nature-ready visualization standards.

Validates:
1. Color palette (max 6 colors)
2. Plot redundancy elimination
3. Label suppression
4. Figure panel layout
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from visualization_standards import (
    NATURE_MOTIF_COLORS, PlotDominance, FigurePanel, MetricFilter,
    LabelPolicy, ValidationThresholds, should_show_plot
)
from utilities import MOTIF_CLASS_COLORS
import matplotlib.pyplot as plt


def test_color_count():
    """Test that we use at most 6 unique colors for main motif classes."""
    # Get unique colors (excluding neutrals)
    unique_colors = set()
    for cls, color in NATURE_MOTIF_COLORS.items():
        if cls not in ['Hybrid', 'Non-B_DNA_Clusters']:
            unique_colors.add(color)
    
    print(f"✓ Color palette test:")
    print(f"  Unique colors (excluding neutrals): {len(unique_colors)}")
    print(f"  Colors: {unique_colors}")
    
    assert len(unique_colors) <= 6, f"Too many colors: {len(unique_colors)} > 6"
    print(f"  ✓ PASS: Color limit enforced (≤6 colors)")


def test_color_consolidation():
    """Test that secondary classes share colors with primary classes."""
    print(f"\n✓ Color consolidation test:")
    
    # i-Motif should match G-Quadruplex
    assert NATURE_MOTIF_COLORS['i-Motif'] == NATURE_MOTIF_COLORS['G-Quadruplex'], \
        "i-Motif should share color with G-Quadruplex"
    print(f"  ✓ i-Motif → G-Quadruplex (complementary structures)")
    
    # A-philic should match Curved DNA
    assert NATURE_MOTIF_COLORS['A-philic_DNA'] == NATURE_MOTIF_COLORS['Curved_DNA'], \
        "A-philic should share color with Curved DNA"
    print(f"  ✓ A-philic_DNA → Curved_DNA (structural affinity)")
    
    # Slipped should match Triplex
    assert NATURE_MOTIF_COLORS['Slipped_DNA'] == NATURE_MOTIF_COLORS['Triplex'], \
        "Slipped should share color with Triplex"
    print(f"  ✓ Slipped_DNA → Triplex (repeat structures)")
    
    print(f"  ✓ PASS: Color consolidation working")


def test_plot_dominance():
    """Test that redundant plots are correctly hidden."""
    print(f"\n✓ Plot dominance test:")
    
    # These should be hidden
    hidden_plots = [
        'plot_motif_distribution',
        'plot_coverage_map',
        'plot_density_heatmap',
        'plot_circos_motif_density',
        'plot_cumulative_motif_distribution'
    ]
    
    for plot_name in hidden_plots:
        assert not PlotDominance.is_plot_allowed(plot_name), \
            f"{plot_name} should be hidden"
        print(f"  ✓ {plot_name} correctly hidden")
    
    print(f"  ✓ PASS: Redundant plots hidden")


def test_required_plots():
    """Test that required plots are identified correctly."""
    print(f"\n✓ Required plots test:")
    
    required = FigurePanel.get_required_plots()
    print(f"  Required plots: {required}")
    
    # Check that essential plots are included
    essential = ['nested_pie_chart', 'manhattan_plot', 'density_comparison', 'cooccurrence_matrix']
    for plot in essential:
        assert plot in required, f"{plot} should be in required plots"
        print(f"  ✓ {plot} is required")
    
    print(f"  ✓ PASS: Required plots correctly identified")


def test_optional_plots():
    """Test that optional plots are correctly marked."""
    print(f"\n✓ Optional plots test:")
    
    optional = FigurePanel.get_optional_plots()
    print(f"  Optional plots: {optional}")
    
    # length_kde should be optional
    assert 'length_kde' in optional, "length_kde should be optional"
    print(f"  ✓ length_kde is optional")
    
    print(f"  ✓ PASS: Optional plots correctly identified")


def test_metric_filtering():
    """Test that metrics are correctly filtered."""
    print(f"\n✓ Metric filtering test:")
    
    # Core metrics should always show
    assert MetricFilter.should_show_metric('Score'), "Score should be shown"
    assert MetricFilter.should_show_metric('Start'), "Start should be shown"
    assert MetricFilter.should_show_metric('Class'), "Class should be shown"
    print(f"  ✓ Core metrics shown")
    
    # Hidden metrics should never show
    assert not MetricFilter.should_show_metric('Raw_DeltaG'), "Raw_DeltaG should be hidden"
    assert not MetricFilter.should_show_metric('Percentile'), "Percentile should be hidden"
    print(f"  ✓ Raw metrics hidden")
    
    # Conditional metrics
    assert MetricFilter.should_show_metric('Num_Tracts', 'G-Quadruplex'), \
        "Num_Tracts should show for G4"
    assert not MetricFilter.should_show_metric('Num_Tracts', 'Curved_DNA'), \
        "Num_Tracts should not show for Curved DNA"
    print(f"  ✓ Conditional metrics work")
    
    print(f"  ✓ PASS: Metric filtering working")


def test_label_policy():
    """Test label suppression policy."""
    print(f"\n✓ Label policy test:")
    
    assert not LabelPolicy.RENDER_INDIVIDUAL_LABELS, \
        "Individual labels should be suppressed"
    print(f"  ✓ Individual labels suppressed")
    
    assert LabelPolicy.RENDER_CLUSTER_CENTROIDS, \
        "Cluster centroids should be labeled"
    print(f"  ✓ Cluster centroids labeled")
    
    assert LabelPolicy.USE_COLLISION_DETECTION, \
        "Collision detection should be enabled"
    print(f"  ✓ Collision detection enabled")
    
    print(f"  ✓ PASS: Label policy correct")


def test_should_show_plot():
    """Test context-dependent plot display logic."""
    print(f"\n✓ Plot display logic test:")
    
    # Large sequence should show Manhattan, not linear
    assert should_show_plot('plot_manhattan_motif_density', 100000, False), \
        "Manhattan should show for large sequences"
    assert not should_show_plot('plot_linear_motif_track', 100000, False), \
        "Linear track should not show for large sequences"
    print(f"  ✓ Large sequence: Manhattan plot")
    
    # Small sequence should show linear, not Manhattan
    assert not should_show_plot('plot_manhattan_motif_density', 10000, False), \
        "Manhattan should not show for small sequences"
    assert should_show_plot('plot_linear_motif_track', 10000, False), \
        "Linear track should show for small sequences"
    print(f"  ✓ Small sequence: Linear track")
    
    # Cluster plot only shows with clusters
    assert not should_show_plot('plot_cluster_size_distribution', 50000, False), \
        "Cluster plot should not show without clusters"
    assert should_show_plot('plot_cluster_size_distribution', 50000, True), \
        "Cluster plot should show with clusters"
    print(f"  ✓ Cluster plot conditional on clusters")
    
    print(f"  ✓ PASS: Context-dependent display working")


def test_validation_thresholds():
    """Test validation threshold functions."""
    print(f"\n✓ Validation thresholds test:")
    
    # Test color count validation
    valid_colors = ['#CC79A7', '#0072B2', '#882255', '#56B4E9', '#E69F00', '#009E73']
    assert ValidationThresholds.validate_color_count(valid_colors), \
        "6 colors should be valid"
    print(f"  ✓ 6 colors valid")
    
    invalid_colors = valid_colors + ['#AAAAAA', '#BBBBBB']
    assert not ValidationThresholds.validate_color_count(invalid_colors), \
        "8 colors should be invalid"
    print(f"  ✓ >6 colors invalid")
    
    # Test figure quality validation
    assert ValidationThresholds.validate_figure_quality(300, ['png', 'pdf']), \
        "300 DPI with png and pdf should be valid"
    print(f"  ✓ 300 DPI valid")
    
    assert not ValidationThresholds.validate_figure_quality(150, ['png', 'pdf']), \
        "150 DPI should be invalid"
    print(f"  ✓ Low DPI invalid")
    
    print(f"  ✓ PASS: Validation thresholds working")


def test_utilities_color_update():
    """Test that utilities.py uses the new color scheme."""
    print(f"\n✓ Utilities color update test:")
    
    # Check that utilities uses consolidated colors
    assert MOTIF_CLASS_COLORS['i-Motif'] == MOTIF_CLASS_COLORS['G-Quadruplex'], \
        "utilities.py should use consolidated i-Motif color"
    print(f"  ✓ utilities.py uses consolidated colors")
    
    # Count unique colors in utilities
    unique_colors = set()
    for cls, color in MOTIF_CLASS_COLORS.items():
        if cls not in ['Hybrid', 'Non-B_DNA_Clusters']:
            unique_colors.add(color)
    
    assert len(unique_colors) <= 6, \
        f"utilities.py should have ≤6 colors, has {len(unique_colors)}"
    print(f"  ✓ utilities.py respects 6-color limit")
    
    print(f"  ✓ PASS: utilities.py updated correctly")


def run_all_tests():
    """Run all validation tests."""
    print("=" * 70)
    print("NATURE-READY VISUALIZATION STANDARDS - VALIDATION TESTS")
    print("=" * 70)
    
    tests = [
        test_color_count,
        test_color_consolidation,
        test_plot_dominance,
        test_required_plots,
        test_optional_plots,
        test_metric_filtering,
        test_label_policy,
        test_should_show_plot,
        test_validation_thresholds,
        test_utilities_color_update
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            passed += 1
        except AssertionError as e:
            print(f"\n✗ FAILED: {test.__name__}")
            print(f"  Error: {e}")
            failed += 1
        except Exception as e:
            print(f"\n✗ ERROR: {test.__name__}")
            print(f"  Error: {e}")
            failed += 1
    
    print("\n" + "=" * 70)
    print(f"TEST SUMMARY: {passed} passed, {failed} failed")
    print("=" * 70)
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
