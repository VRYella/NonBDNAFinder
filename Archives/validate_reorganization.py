#!/usr/bin/env python3
"""
Basic validation test for the reorganized NonBDNAFinder structure.
Tests that core functionality works after reorganization.
"""

import sys
import os

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def test_imports():
    """Test that all main modules can be imported"""
    print("Testing imports...")
    
    try:
        from BASE_CODES.utils import parse_fasta, gc_content
        print("✓ BASE_CODES.utils imports successfully")
    except Exception as e:
        print(f"✗ BASE_CODES.utils import failed: {e}")
        return False
    
    try:
        from REGISTRIES.motifs import get_basic_stats
        print("✓ REGISTRIES.motifs imports successfully")
    except Exception as e:
        print(f"✗ REGISTRIES.motifs import failed: {e}")
        return False
    
    try:
        from HYPERSCAN.hyperscan_integration import all_motifs_refactored
        print("✓ HYPERSCAN.hyperscan_integration imports successfully")
    except Exception as e:
        print(f"✗ HYPERSCAN.hyperscan_integration import failed: {e}")
        return False
    
    return True

def test_basic_functionality():
    """Test basic motif detection functionality"""
    print("\nTesting basic functionality...")
    
    try:
        from BASE_CODES.utils import parse_fasta
        from HYPERSCAN.hyperscan_integration import all_motifs_refactored
        
        # Test sequence with known motifs
        test_seq = "ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC"
        
        # Parse the sequence
        parsed_seq = parse_fasta(test_seq)
        print(f"✓ Sequence parsing works: {len(parsed_seq)} bp")
        
        # Run motif detection
        motifs = all_motifs_refactored(parsed_seq, sequence_name="test", nonoverlap=False)
        print(f"✓ Motif detection works: {len(motifs)} motifs found")
        
        return True
        
    except Exception as e:
        print(f"✗ Basic functionality test failed: {e}")
        return False

def test_app_import():
    """Test that the main app can be imported"""
    print("\nTesting main app import...")
    
    try:
        import app
        print("✓ Main app imports successfully")
        return True
    except Exception as e:
        print(f"✗ Main app import failed: {e}")
        return False

def main():
    """Run all validation tests"""
    print("="*60)
    print("NONBDNA FINDER REORGANIZATION VALIDATION TEST")
    print("="*60)
    
    tests = [
        test_imports,
        test_basic_functionality,
        test_app_import
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"✗ Test {test.__name__} failed with exception: {e}")
    
    print("\n" + "="*60)
    print(f"RESULTS: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! Reorganization successful!")
        return True
    else:
        print("❌ Some tests failed. Please check the issues above.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)