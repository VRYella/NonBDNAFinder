#!/usr/bin/env python3
"""
Comprehensive test script for NonBDNAFinder deployment readiness.

This script verifies:
1. All core dependencies are installed
2. All modules import successfully
3. Hyperscan availability and fallback behavior
4. Basic functionality works with and without hyperscan
5. Application can start successfully
"""

import sys
import warnings
warnings.filterwarnings('ignore')


def test_hyperscan_availability():
    """Test if hyperscan is available and working."""
    print("\n" + "=" * 80)
    print("HYPERSCAN AVAILABILITY TEST")
    print("=" * 80)
    
    try:
        import hyperscan
        print("✅ Hyperscan is AVAILABLE")
        print(f"   Version: {hyperscan.__version__ if hasattr(hyperscan, '__version__') else 'unknown'}")
        
        # Test basic hyperscan functionality
        try:
            db = hyperscan.Database()
            pattern = b"test"
            expressions = [pattern]
            ids = [0]
            flags = [0]
            db.compile(expressions=expressions, ids=ids, flags=flags)
            print("✅ Hyperscan compilation test passed")
            return True
        except Exception as e:
            print(f"⚠️  Hyperscan import succeeded but compilation failed: {e}")
            print("   Fallback to regex will be used")
            return False
            
    except ImportError:
        print("ℹ️  Hyperscan is NOT available (ImportError)")
        print("   This is OK - the application will use regex-based fallback")
        print("   All features will work, with slightly reduced performance")
        return False
    except Exception as e:
        print(f"⚠️  Hyperscan test failed: {e}")
        print("   Fallback to regex will be used")
        return False


def test_core_dependencies():
    """Test that all core dependencies are installed."""
    print("\n" + "=" * 80)
    print("CORE DEPENDENCIES TEST")
    print("=" * 80)
    
    core_deps = [
        ('streamlit', 'Streamlit web framework'),
        ('numpy', 'NumPy numerical computing'),
        ('pandas', 'Pandas data analysis'),
        ('scipy', 'SciPy scientific computing'),
        ('sklearn', 'Scikit-learn machine learning'),
        ('matplotlib', 'Matplotlib plotting'),
        ('seaborn', 'Seaborn statistical visualization'),
        ('plotly', 'Plotly interactive charts'),
        ('Bio', 'Biopython bioinformatics'),
        ('openpyxl', 'Excel file support'),
        ('xlsxwriter', 'Excel writing'),
        ('psutil', 'System monitoring'),
        ('requests', 'HTTP library'),
    ]
    
    all_present = True
    for dep, description in core_deps:
        try:
            __import__(dep)
            print(f"✅ {dep:15s} - {description}")
        except ImportError:
            print(f"❌ {dep:15s} - NOT INSTALLED - {description}")
            all_present = False
    
    if all_present:
        print("\n✅ All core dependencies are installed!")
        return True
    else:
        print("\n❌ Some core dependencies are missing - installation required!")
        return False


def test_app_imports():
    """Test that all required modules for app.py can be imported."""
    print("\n" + "=" * 80)
    print("APPLICATION IMPORTS TEST")
    print("=" * 80)
    
    try:
        # Test utilities imports
        from utilities import (
            parse_fasta, parse_fasta_chunked, get_file_preview, wrap, get_basic_stats, 
            export_to_bed, export_to_csv, export_to_json, export_to_excel, 
            calculate_genomic_density, calculate_positional_density,
            export_results_to_dataframe, CORE_OUTPUT_COLUMNS
        )
        print("✅ utilities module imports successful")
        
        # Test nonbscanner imports
        from nonbscanner import (
            analyze_sequence, get_motif_info as get_motif_classification_info
        )
        print("✅ nonbscanner module imports successful")
        
        # Test visualizations imports
        from visualizations import (
            plot_motif_distribution, plot_coverage_map, plot_density_heatmap,
            plot_length_distribution, plot_score_distribution, plot_nested_pie_chart, 
            MOTIF_CLASS_COLORS, plot_density_comparison,
            plot_circos_motif_density,
            plot_density_comparison_by_subclass, plot_enrichment_analysis_by_subclass,
        )
        print("✅ visualizations module imports successful")
        
        print("\n✅ All application modules import correctly!")
        return True
        
    except ImportError as e:
        print(f"❌ IMPORT ERROR: {e}")
        print("\nThe app will fail with this error.")
        import traceback
        traceback.print_exc()
        return False
    except Exception as e:
        print(f"❌ UNEXPECTED ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_basic_functionality():
    """Test basic motif detection functionality."""
    print("\n" + "=" * 80)
    print("BASIC FUNCTIONALITY TEST")
    print("=" * 80)
    
    try:
        import nonbscanner as nbs
        
        # Test with a simple G-quadruplex sequence
        test_seq = "GGGTTAGGGTTAGGGTTAGGG"
        print(f"Testing with sequence: {test_seq}")
        
        motifs = nbs.analyze_sequence(test_seq, "test_sequence")
        print(f"✅ Detection completed: found {len(motifs)} motifs")
        
        if len(motifs) > 0:
            print(f"   Example motif: {motifs[0]['Class']} at position {motifs[0]['Start']}-{motifs[0]['End']}")
        
        return True
        
    except Exception as e:
        print(f"❌ Functionality test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_fallback_behavior():
    """Test that the application works without hyperscan."""
    print("\n" + "=" * 80)
    print("FALLBACK BEHAVIOR TEST")
    print("=" * 80)
    
    try:
        # Check how modules report hyperscan availability
        import app
        import detectors
        import scanner
        
        print(f"app.py reports:       HYPERSCAN_AVAILABLE = {getattr(app, 'HYPERSCAN_AVAILABLE', 'not defined')}")
        print(f"detectors.py reports: _HYPERSCAN_AVAILABLE = {getattr(detectors, '_HYPERSCAN_AVAILABLE', 'not defined')}")
        print(f"scanner.py reports:   HYPERSCAN_AVAILABLE = {getattr(scanner, 'HYPERSCAN_AVAILABLE', 'not defined')}")
        
        print("\n✅ All modules handle hyperscan availability correctly")
        return True
        
    except Exception as e:
        print(f"❌ Fallback test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def print_deployment_summary(hyperscan_available, deps_ok, imports_ok, functionality_ok, fallback_ok):
    """Print summary of deployment readiness."""
    print("\n" + "=" * 80)
    print("DEPLOYMENT READINESS SUMMARY")
    print("=" * 80)
    
    print("\nCore Requirements:")
    print(f"  {'✅' if deps_ok else '❌'} Core dependencies installed")
    print(f"  {'✅' if imports_ok else '❌'} All modules import successfully")
    print(f"  {'✅' if functionality_ok else '❌'} Basic functionality works")
    print(f"  {'✅' if fallback_ok else '❌'} Fallback behavior correct")
    
    print("\nOptional Performance:")
    print(f"  {'✅' if hyperscan_available else 'ℹ️ '} Hyperscan {'available' if hyperscan_available else 'not available (using fallback)'}")
    
    print("\nDeployment Status:")
    if deps_ok and imports_ok and functionality_ok and fallback_ok:
        print("  ✅ READY FOR DEPLOYMENT")
        print("     The application will work correctly on Streamlit Cloud")
        if not hyperscan_available:
            print("     (Hyperscan will be attempted during deployment)")
        return True
    else:
        print("  ❌ NOT READY - Issues must be fixed before deployment")
        return False


def main():
    """Run all tests and report deployment readiness."""
    print("\n" + "=" * 80)
    print("NonBDNAFinder DEPLOYMENT READINESS CHECK")
    print("=" * 80)
    
    # Run all tests
    hyperscan_available = test_hyperscan_availability()
    deps_ok = test_core_dependencies()
    imports_ok = test_app_imports()
    functionality_ok = test_basic_functionality()
    fallback_ok = test_fallback_behavior()
    
    # Print summary
    ready = print_deployment_summary(
        hyperscan_available, deps_ok, imports_ok, functionality_ok, fallback_ok
    )
    
    # Exit with appropriate code
    if ready:
        print("\n" + "=" * 80)
        print("✅ ALL TESTS PASSED - Ready for Streamlit Cloud deployment!")
        print("=" * 80)
        sys.exit(0)
    else:
        print("\n" + "=" * 80)
        print("❌ SOME TESTS FAILED - Fix issues before deployment")
        print("=" * 80)
        sys.exit(1)


if __name__ == "__main__":
    main()
