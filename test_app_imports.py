#!/usr/bin/env python3
"""
Test script to verify that all modules required by app.py can be imported successfully.

This test prevents regression of the ImportError issue on Streamlit Cloud.
"""

import sys
import warnings
warnings.filterwarnings('ignore')


def test_app_imports():
    """Test that all required modules for app.py can be imported."""
    print("Testing app.py imports...")
    print("-" * 80)
    
    try:
        # Test utilities imports (the problematic import from line 37)
        from utilities import (
            parse_fasta, parse_fasta_chunked, get_file_preview, wrap, get_basic_stats, 
            export_to_bed, export_to_csv, export_to_json, export_to_excel, 
            calculate_genomic_density, calculate_positional_density,
            export_results_to_dataframe, CORE_OUTPUT_COLUMNS
        )
        print("✅ utilities imports successful!")
        
        # Test nonbscanner imports
        from nonbscanner import (
            analyze_sequence, get_motif_info as get_motif_classification_info
        )
        print("✅ nonbscanner imports successful!")
        
        # Test visualizations imports
        from visualizations import (
            plot_motif_distribution, plot_coverage_map, plot_density_heatmap,
            plot_length_distribution, plot_score_distribution, plot_nested_pie_chart, 
            MOTIF_CLASS_COLORS, plot_density_comparison,
            plot_circos_motif_density,
            plot_density_comparison_by_subclass, plot_enrichment_analysis_by_subclass,
        )
        print("✅ visualizations imports successful!")
        
        print("\n" + "=" * 80)
        print("✅ All app.py imports working correctly!")
        print("=" * 80)
        return True
        
    except ImportError as e:
        print(f"❌ IMPORT ERROR: {e}")
        print("\nThe app will fail with this error on Streamlit Cloud.")
        import traceback
        traceback.print_exc()
        return False
        
    except Exception as e:
        print(f"❌ UNEXPECTED ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_required_dependencies():
    """Test that all required dependencies are installed."""
    print("\nTesting required dependencies...")
    print("-" * 80)
    
    required_deps = [
        'streamlit',
        'numpy',
        'pandas',
        'scipy',
        'sklearn',
        'matplotlib',
        'seaborn',
        'plotly',
        'Bio',  # biopython
        'openpyxl',
        'xlsxwriter',
        'psutil',
    ]
    
    all_present = True
    for dep in required_deps:
        try:
            __import__(dep)
            print(f"✅ {dep}")
        except ImportError:
            print(f"❌ {dep} - NOT INSTALLED")
            all_present = False
    
    # Check hyperscan separately (optional)
    try:
        __import__('hyperscan')
        print(f"✅ hyperscan (optional)")
    except ImportError:
        print(f"⚠️  hyperscan (optional) - not installed")
    
    if all_present:
        print("\n✅ All required dependencies are installed!")
        return True
    else:
        print("\n❌ Some required dependencies are missing!")
        return False


if __name__ == "__main__":
    print("\n" + "=" * 80)
    print("APP.PY IMPORTS TEST")
    print("=" * 80 + "\n")
    
    deps_success = test_required_dependencies()
    import_success = test_app_imports()
    
    if deps_success and import_success:
        print("\n" + "=" * 80)
        print("✅ ALL TESTS PASSED - app.py should work on Streamlit Cloud")
        print("=" * 80)
        sys.exit(0)
    else:
        print("\n" + "=" * 80)
        print("❌ SOME TESTS FAILED - app.py will fail on Streamlit Cloud")
        print("=" * 80)
        sys.exit(1)
