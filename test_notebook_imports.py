#!/usr/bin/env python3
"""
Test script to verify that all visualization functions required by 
HighEfficiency_Genome_Analysis.ipynb can be imported successfully.

This test prevents regression of the ImportError issue.
"""

import sys
import warnings
warnings.filterwarnings('ignore')


def test_imports():
    """Test that all required visualization functions can be imported."""
    print("Testing Jupyter notebook imports...")
    print("-" * 80)
    
    try:
        # Import exactly as the notebook does
        from visualizations import (
            plot_comprehensive_class_analysis,
            plot_comprehensive_subclass_analysis,
            plot_score_statistics_by_class,
            plot_length_statistics_by_class,
            plot_motif_distribution,
            plot_coverage_map,
            plot_genome_landscape_track,
            plot_sliding_window_heat_ribbon
        )
        
        print("✅ All imports successful!\n")
        
        # Verify functions are callable
        functions = [
            plot_comprehensive_class_analysis,
            plot_comprehensive_subclass_analysis,
            plot_score_statistics_by_class,
            plot_length_statistics_by_class,
            plot_motif_distribution,
            plot_coverage_map,
            plot_genome_landscape_track,
            plot_sliding_window_heat_ribbon
        ]
        
        print("Imported functions:")
        for i, func in enumerate(functions, 1):
            print(f"  {i}. {func.__name__}")
            assert callable(func), f"{func.__name__} is not callable"
        
        print("\n" + "=" * 80)
        print("✅ All tests passed! Notebook imports are working correctly.")
        print("=" * 80)
        return True
        
    except ImportError as e:
        print(f"❌ IMPORT ERROR: {e}")
        print("\nThe notebook will fail with this error.")
        return False
        
    except AssertionError as e:
        print(f"❌ ASSERTION ERROR: {e}")
        return False
        
    except Exception as e:
        print(f"❌ UNEXPECTED ERROR: {e}")
        return False


def test_function_execution():
    """Test that functions can be executed with sample data."""
    print("\nTesting function execution with sample data...")
    print("-" * 80)
    
    try:
        from visualizations import (
            plot_comprehensive_class_analysis,
            plot_genome_landscape_track,
            plot_sliding_window_heat_ribbon
        )
        import matplotlib.pyplot as plt
        
        # Sample test data
        test_motifs = [
            {
                'Class': 'G-Quadruplex',
                'Subclass': 'Canonical G4',
                'Start': 100,
                'End': 150,
                'Length': 50,
                'Score': 0.95
            }
        ]
        sequence_length = 1000
        
        # Test new functions
        fig1 = plot_genome_landscape_track(test_motifs, sequence_length)
        plt.close(fig1)
        print("✅ plot_genome_landscape_track executes successfully")
        
        fig2 = plot_sliding_window_heat_ribbon(test_motifs, sequence_length)
        plt.close(fig2)
        print("✅ plot_sliding_window_heat_ribbon executes successfully")
        
        fig3 = plot_comprehensive_class_analysis(test_motifs)
        plt.close(fig3)
        print("✅ plot_comprehensive_class_analysis executes successfully")
        
        print("\n✅ Function execution tests passed!")
        return True
        
    except Exception as e:
        print(f"❌ Function execution failed: {e}")
        return False


if __name__ == "__main__":
    print("\n" + "=" * 80)
    print("VISUALIZATION IMPORTS TEST")
    print("=" * 80 + "\n")
    
    import_success = test_imports()
    execution_success = test_function_execution()
    
    if import_success and execution_success:
        print("\n" + "=" * 80)
        print("✅ ALL TESTS PASSED")
        print("=" * 80)
        sys.exit(0)
    else:
        print("\n" + "=" * 80)
        print("❌ SOME TESTS FAILED")
        print("=" * 80)
        sys.exit(1)
