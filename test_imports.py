#!/usr/bin/env python3
"""
Test script to verify all imports in app.py work correctly.
This test ensures the ImportError is resolved.
"""

import sys
import traceback

def test_utilities_imports():
    """Test imports from utilities module."""
    print("Testing utilities imports...")
    try:
        from utilities import (
            parse_fasta, parse_fasta_chunked, get_file_preview, wrap, 
            get_basic_stats, export_to_bed, export_to_csv,
            export_to_json, export_to_excel, calculate_genomic_density, 
            calculate_positional_density, export_results_to_dataframe
        )
        print("✓ All utilities imports successful")
        return True
    except ImportError as e:
        print(f"✗ Utilities import failed: {e}")
        traceback.print_exc()
        return False

def test_nonbscanner_imports():
    """Test imports from nonbscanner module."""
    print("\nTesting nonbscanner imports...")
    try:
        from nonbscanner import (
            analyze_sequence, get_motif_info as get_motif_classification_info
        )
        print("✓ All nonbscanner imports successful")
        return True
    except ImportError as e:
        print(f"✗ Nonbscanner import failed: {e}")
        traceback.print_exc()
        return False

def test_visualizations_imports():
    """Test imports from visualizations module."""
    print("\nTesting visualizations imports...")
    try:
        from visualizations import (
            plot_motif_distribution, plot_coverage_map, plot_density_heatmap,
            plot_length_distribution, plot_score_distribution, plot_nested_pie_chart, 
            MOTIF_CLASS_COLORS, plot_density_comparison,
            plot_circos_motif_density,
            plot_density_comparison_by_subclass, plot_enrichment_analysis_by_subclass,
            plot_subclass_density_heatmap,
            # New Nature-quality visualizations
            plot_manhattan_motif_density, plot_cumulative_motif_distribution,
            plot_motif_cooccurrence_matrix, plot_gc_content_correlation,
            plot_linear_motif_track, plot_cluster_size_distribution,
            plot_motif_length_kde
        )
        print("✓ All visualizations imports successful")
        print(f"  - Standard visualization functions: 8")
        print(f"  - Nature-quality visualization functions: 7")
        print(f"  - Color constants: 1")
        return True
    except ImportError as e:
        print(f"✗ Visualizations import failed: {e}")
        traceback.print_exc()
        return False

def test_app_imports():
    """Test that app.py can be imported without errors."""
    print("\nTesting app.py full import...")
    try:
        import app
        print("✓ app.py imported successfully")
        return True
    except ImportError as e:
        print(f"✗ app.py import failed: {e}")
        traceback.print_exc()
        return False

def main():
    """Run all import tests."""
    print("=" * 70)
    print("NonBDNAFinder Import Verification Test")
    print("=" * 70)
    
    results = []
    results.append(("Utilities", test_utilities_imports()))
    results.append(("Nonbscanner", test_nonbscanner_imports()))
    results.append(("Visualizations", test_visualizations_imports()))
    results.append(("App.py", test_app_imports()))
    
    print("\n" + "=" * 70)
    print("Test Summary:")
    print("=" * 70)
    
    all_passed = True
    for module, passed in results:
        status = "PASS ✓" if passed else "FAIL ✗"
        print(f"{module:20} : {status}")
        if not passed:
            all_passed = False
    
    print("=" * 70)
    if all_passed:
        print("✓ All imports successful! The ImportError is resolved.")
        return 0
    else:
        print("✗ Some imports failed. Please check the errors above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
