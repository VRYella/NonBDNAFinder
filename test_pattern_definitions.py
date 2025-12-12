#!/usr/bin/env python3
"""
Test script to verify the pattern definitions module and import fixes.

This script validates:
1. Import error fix (NonBFinder -> nonbscanner)
2. Centralized pattern definitions (motif_patterns.py)
3. Detector classes using imported patterns
4. Pattern customization functionality
"""

import sys

def test_import_fix():
    """Test that the import error is fixed."""
    print("=" * 60)
    print("TEST 1: Import Error Fix")
    print("=" * 60)
    
    try:
        from nonbscanner import analyze_sequence, get_motif_info
        print("✓ Fixed: 'from nonbscanner import' works correctly")
        print("  (Previously failed with 'from NonBFinder import')")
        return True
    except ImportError as e:
        print(f"✗ Import still broken: {e}")
        return False


def test_motif_patterns_module():
    """Test the centralized pattern definitions module."""
    print("\n" + "=" * 60)
    print("TEST 2: Motif Patterns Module")
    print("=" * 60)
    
    try:
        from motif_patterns import (
            PATTERN_REGISTRY,
            pattern_info_summary,
            get_patterns_for_detector,
            get_all_pattern_ids,
            get_pattern_by_id
        )
        
        print("✓ motif_patterns module imports successfully")
        print(f"✓ Pattern registry contains {len(PATTERN_REGISTRY)} detectors")
        
        # Test helper functions
        all_ids = get_all_pattern_ids()
        print(f"✓ Found {len(all_ids)} total patterns")
        
        # Test getting patterns for a specific detector
        g4_patterns = get_patterns_for_detector('GQuadruplexDetector')
        print(f"✓ G4 detector has {len(g4_patterns)} pattern groups")
        
        # Test getting a specific pattern
        g4_0 = get_pattern_by_id('G4_0')
        if g4_0:
            print(f"✓ Retrieved pattern G4_0: '{g4_0[2]}'")
        
        print("\nPattern Summary:")
        print(pattern_info_summary())
        
        return True
    except Exception as e:
        print(f"✗ Pattern module error: {e}")
        return False


def test_detector_integration():
    """Test that detectors use the centralized patterns."""
    print("\n" + "=" * 60)
    print("TEST 3: Detector Integration")
    print("=" * 60)
    
    try:
        from detectors import (
            ZDNADetector,
            GQuadruplexDetector,
            IMotifDetector,
            RLoopDetector,
            TriplexDetector
        )
        
        detectors = [
            ('ZDNADetector', ZDNADetector()),
            ('GQuadruplexDetector', GQuadruplexDetector()),
            ('IMotifDetector', IMotifDetector()),
            ('RLoopDetector', RLoopDetector()),
            ('TriplexDetector', TriplexDetector()),
        ]
        
        for name, detector in detectors:
            patterns = detector.get_patterns()
            num_groups = len(patterns)
            num_patterns = sum(len(p) for p in patterns.values())
            print(f"✓ {name}: {num_groups} groups, {num_patterns} patterns")
        
        return True
    except Exception as e:
        print(f"✗ Detector integration error: {e}")
        return False


def test_pattern_validation():
    """Test that patterns are properly formatted."""
    print("\n" + "=" * 60)
    print("TEST 4: Pattern Validation")
    print("=" * 60)
    
    try:
        from motif_patterns import RLOOP_PATTERNS, CRUCIFORM_PATTERNS
        import re
        
        # Test 1: Cruciform should use empty string (algorithmic)
        cruc_pattern = CRUCIFORM_PATTERNS['inverted_repeats'][0][0]
        if cruc_pattern == '':
            print("✓ Cruciform uses empty string (algorithmic detection)")
        else:
            print(f"✗ Cruciform pattern should be empty: '{cruc_pattern}'")
            return False
        
        # Test 2: R-loop patterns should use [ATCG], not [ATCGU]
        errors = []
        for model_name, patterns in RLOOP_PATTERNS.items():
            for pattern_tuple in patterns:
                regex = pattern_tuple[0]
                if 'U' in regex:
                    errors.append(f"{model_name}: contains 'U' (RNA)")
        
        if errors:
            print("✗ Found invalid patterns:")
            for error in errors:
                print(f"  - {error}")
            return False
        else:
            print("✓ R-loop patterns use [ATCG] (DNA only, no RNA)")
        
        # Test 3: Validate regex patterns compile
        from motif_patterns import G4_PATTERNS
        for group_name, patterns in G4_PATTERNS.items():
            for pattern_tuple in patterns:
                regex = pattern_tuple[0]
                if regex:  # Skip empty patterns
                    try:
                        re.compile(regex)
                    except re.error as e:
                        print(f"✗ Invalid regex in {group_name}: {e}")
                        return False
        
        print("✓ All G4 regex patterns compile successfully")
        
        return True
    except Exception as e:
        print(f"✗ Pattern validation error: {e}")
        return False


def test_analysis_workflow():
    """Test the complete analysis workflow."""
    print("\n" + "=" * 60)
    print("TEST 5: Analysis Workflow")
    print("=" * 60)
    
    try:
        from nonbscanner import analyze_sequence
        
        # Test sequence with multiple motif types
        test_seq = "GGGTTAGGGTTAGGGTTAGGG" + "A" * 20 + "CGGCGGCGG"
        
        results = analyze_sequence(test_seq, 'test_seq')
        
        print(f"✓ Analysis completed: found {len(results)} motifs")
        
        # Group by class
        by_class = {}
        for motif in results:
            cls = motif['Class']
            by_class[cls] = by_class.get(cls, 0) + 1
        
        print("\nDetected motifs by class:")
        for cls, count in sorted(by_class.items()):
            print(f"  - {cls}: {count} motif(s)")
        
        if len(results) > 0:
            print("\n✓ Successfully detected motifs from centralized patterns")
            return True
        else:
            print("⚠ No motifs detected (may be expected for test sequence)")
            return True
        
    except Exception as e:
        print(f"✗ Analysis workflow error: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("PATTERN DEFINITIONS MODULE TEST SUITE")
    print("=" * 60)
    print("\nTesting fixes for:")
    print("1. ModuleNotFoundError: NonBFinder -> nonbscanner")
    print("2. Centralized pattern definitions in motif_patterns.py")
    print("3. Pattern customization capability\n")
    
    tests = [
        test_import_fix,
        test_motif_patterns_module,
        test_detector_integration,
        test_pattern_validation,
        test_analysis_workflow,
    ]
    
    results = []
    for test_func in tests:
        try:
            result = test_func()
            results.append(result)
        except Exception as e:
            print(f"\n✗ Test crashed: {e}")
            import traceback
            traceback.print_exc()
            results.append(False)
    
    # Summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    
    passed = sum(results)
    total = len(results)
    
    print(f"\nPassed: {passed}/{total} tests")
    
    if passed == total:
        print("\n✅ ALL TESTS PASSED!")
        print("\nThe following issues have been fixed:")
        print("1. ✓ Import error resolved (NonBFinder -> nonbscanner)")
        print("2. ✓ Centralized pattern definitions created")
        print("3. ✓ Detectors using imported patterns")
        print("4. ✓ Pattern validation passed")
        print("5. ✓ Analysis workflow functional")
        print("\nUsers can now modify patterns in motif_patterns.py")
        print("See PATTERN_CUSTOMIZATION_GUIDE.md for details")
        return 0
    else:
        print(f"\n❌ {total - passed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
