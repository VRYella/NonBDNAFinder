#!/usr/bin/env python3
"""
Simple test to verify unified application structure works correctly.

Tests:
1. All core modules import successfully
2. Basic motif detection works
3. Export functions are available
"""

def test_imports():
    """Test that all core modules import correctly"""
    print("Testing core module imports...")
    
    try:
        import utilities
        print("✓ utilities module imported")
    except ImportError as e:
        print(f"✗ Failed to import utilities: {e}")
        return False
    
    try:
        import nonbscanner
        print("✓ nonbscanner module imported")
    except ImportError as e:
        print(f"✗ Failed to import nonbscanner: {e}")
        return False
    
    try:
        import detectors
        print("✓ detectors module imported")
    except ImportError as e:
        print(f"✗ Failed to import detectors: {e}")
        return False
    
    return True


def test_basic_detection():
    """Test basic motif detection functionality"""
    print("\nTesting basic motif detection...")
    
    try:
        from nonbscanner import analyze_sequence
        
        # Simple test sequence with known G4 motif
        test_seq = "GGGTTAGGGTTAGGGTTAGGG"
        motifs = analyze_sequence(test_seq, "test_sequence")
        
        print(f"✓ Analysis completed: found {len(motifs)} motifs")
        
        if len(motifs) > 0:
            print(f"  - First motif: {motifs[0]['Class']} at position {motifs[0]['Start']}")
            return True
        else:
            print("  - No motifs found (unexpected for test sequence)")
            return False
            
    except Exception as e:
        print(f"✗ Detection failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_export_functions():
    """Test that export functions are available"""
    print("\nTesting export function availability...")
    
    try:
        from utilities import (
            export_to_csv, export_to_bed, export_to_json, 
            export_to_excel, plot_motif_distribution
        )
        print("✓ All export functions available")
        return True
    except ImportError as e:
        print(f"✗ Failed to import export functions: {e}")
        return False


def main():
    """Run all tests"""
    print("="*60)
    print("Testing Unified NonBDNAFinder Application (4 Core Files)")
    print("="*60)
    
    results = []
    
    # Test 1: Imports
    results.append(("Module Imports", test_imports()))
    
    # Test 2: Basic detection (only if imports work)
    if results[0][1]:
        results.append(("Basic Detection", test_basic_detection()))
        results.append(("Export Functions", test_export_functions()))
    
    # Summary
    print("\n" + "="*60)
    print("Test Summary:")
    print("="*60)
    
    for test_name, passed in results:
        status = "PASS" if passed else "FAIL"
        symbol = "✓" if passed else "✗"
        print(f"{symbol} {test_name}: {status}")
    
    all_passed = all(result[1] for result in results)
    
    print("\n" + "="*60)
    if all_passed:
        print("✓ All tests PASSED - Unified application is functional!")
    else:
        print("✗ Some tests FAILED - Please check the errors above")
    print("="*60)
    
    return 0 if all_passed else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
