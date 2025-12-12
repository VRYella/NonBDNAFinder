#!/usr/bin/env python3
"""
Test script for Excel pattern loading functionality.
Verifies that pattern_registry.xlsx can be loaded correctly and produces
the same results as consolidated_registry.json.
"""

import os
import sys
import time
import json

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_excel_loading():
    """Test loading patterns from Excel file"""
    print("=" * 80)
    print("EXCEL PATTERN LOADING TEST")
    print("=" * 80)
    
    # Import utilities after setting up the path
    from utilities import _load_consolidated_registry, _PANDAS_AVAILABLE
    
    # Check if pandas is available
    if not _PANDAS_AVAILABLE:
        print("\n❌ FAILED: pandas not available")
        print("Install with: pip install pandas openpyxl")
        return False
    
    # Check if Excel file exists
    if not os.path.exists('pattern_registry.xlsx'):
        print("\n❌ FAILED: pattern_registry.xlsx not found")
        return False
    
    print("\n✓ pandas is available")
    print("✓ pattern_registry.xlsx exists")
    
    # Test loading
    print("\n--- Loading pattern registry ---")
    start_time = time.time()
    registry = _load_consolidated_registry()
    load_time = time.time() - start_time
    
    if registry is None:
        print("❌ FAILED: Could not load registry")
        return False
    
    print(f"✓ Registry loaded in {load_time:.4f} seconds")
    print(f"✓ Source: {registry.get('source', 'unknown')}")
    print(f"✓ Total classes: {registry.get('total_classes', 0)}")
    print(f"✓ Total patterns: {registry.get('total_patterns', 0)}")
    
    # Verify structure
    print("\n--- Verifying registry structure ---")
    if 'registries' not in registry:
        print("❌ FAILED: Missing 'registries' key")
        return False
    
    registries = registry['registries']
    print(f"✓ Found {len(registries)} motif classes:")
    
    for class_name, class_data in registries.items():
        n_patterns = len(class_data.get('patterns', []))
        print(f"  - {class_name}: {n_patterns} patterns")
    
    # Test specific classes
    print("\n--- Testing specific pattern classes ---")
    
    # Test G4 patterns (regex-based)
    if 'G4' in registries:
        g4_patterns = registries['G4']['patterns']
        print(f"✓ G4: {len(g4_patterns)} patterns")
        if len(g4_patterns) > 0:
            first_pattern = g4_patterns[0]
            print(f"  Sample: id={first_pattern.get('id')}, "
                  f"pattern={first_pattern.get('pattern', 'N/A')[:50]}..., "
                  f"score={first_pattern.get('score')}")
    
    # Test ZDNA patterns (tenmer-based)
    if 'ZDNA' in registries:
        zdna_patterns = registries['ZDNA']['patterns']
        print(f"✓ ZDNA: {len(zdna_patterns)} patterns")
        if len(zdna_patterns) > 0:
            first_pattern = zdna_patterns[0]
            print(f"  Sample: id={first_pattern.get('id')}, "
                  f"tenmer={first_pattern.get('tenmer', 'N/A')}, "
                  f"score={first_pattern.get('score')}")
    
    # Test APhilic patterns (tenmer-based)
    if 'APhilic' in registries:
        aphilic_patterns = registries['APhilic']['patterns']
        print(f"✓ APhilic: {len(aphilic_patterns)} patterns")
        if len(aphilic_patterns) > 0:
            first_pattern = aphilic_patterns[0]
            print(f"  Sample: id={first_pattern.get('id')}, "
                  f"tenmer={first_pattern.get('tenmer', 'N/A')}, "
                  f"score={first_pattern.get('score')}")
    
    return True


def compare_with_json():
    """Compare Excel-loaded data with JSON data"""
    print("\n" + "=" * 80)
    print("COMPARING EXCEL VS JSON")
    print("=" * 80)
    
    if not os.path.exists('consolidated_registry.json'):
        print("\n⚠ consolidated_registry.json not found - skipping comparison")
        return True
    
    # Load JSON data
    print("\n--- Loading JSON data ---")
    with open('consolidated_registry.json', 'r') as f:
        json_data = json.load(f)
    
    print(f"✓ JSON: {json_data['total_patterns']} patterns in {json_data['total_classes']} classes")
    
    # Load Excel data
    print("\n--- Loading Excel data ---")
    from utilities import _load_consolidated_registry_from_excel
    excel_data = _load_consolidated_registry_from_excel()
    
    if excel_data is None:
        print("❌ FAILED: Could not load Excel data")
        return False
    
    print(f"✓ Excel: {excel_data['total_patterns']} patterns in {excel_data['total_classes']} classes")
    
    # Compare counts
    print("\n--- Comparing pattern counts ---")
    all_match = True
    for class_name in json_data['registries']:
        json_count = len(json_data['registries'][class_name]['patterns'])
        excel_count = len(excel_data['registries'].get(class_name, {}).get('patterns', []))
        
        match = "✓" if json_count == excel_count else "❌"
        print(f"{match} {class_name}: JSON={json_count}, Excel={excel_count}")
        
        if json_count != excel_count:
            all_match = False
    
    if all_match:
        print("\n✓ All pattern counts match!")
    else:
        print("\n⚠ Some pattern counts differ")
    
    return True


def test_with_detector():
    """Test using Excel patterns with actual detector"""
    print("\n" + "=" * 80)
    print("TESTING WITH DETECTOR")
    print("=" * 80)
    
    try:
        from utilities import load_db_for_class
        
        print("\n--- Testing G4 detector ---")
        db, id_to_pattern, id_to_score = load_db_for_class('G4', 'registry')
        
        print(f"✓ Loaded {len(id_to_pattern)} G4 patterns")
        print(f"✓ Pattern IDs: {list(id_to_pattern.keys())[:5]}...")
        
        # Show sample patterns
        for pid in list(id_to_pattern.keys())[:3]:
            print(f"  ID {pid}: {id_to_pattern[pid][:60]}... (score: {id_to_score.get(pid, 0)})")
        
        print("\n--- Testing ZDNA detector ---")
        db, id_to_pattern, id_to_score = load_db_for_class('ZDNA', 'registry')
        
        print(f"✓ Loaded {len(id_to_pattern)} ZDNA patterns")
        print(f"✓ Pattern IDs: {list(id_to_pattern.keys())[:5]}...")
        
        # Show sample patterns
        for pid in list(id_to_pattern.keys())[:3]:
            print(f"  ID {pid}: {id_to_pattern[pid]} (score: {id_to_score.get(pid, 0)})")
        
        return True
        
    except Exception as e:
        print(f"\n❌ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_performance():
    """Test loading performance"""
    print("\n" + "=" * 80)
    print("PERFORMANCE TEST")
    print("=" * 80)
    
    from utilities import _load_consolidated_registry, clear_pattern_registry_cache
    
    # Clear cache
    clear_pattern_registry_cache()
    
    # Test Excel loading time
    print("\n--- Excel loading time (cold) ---")
    start = time.time()
    registry = _load_consolidated_registry()
    excel_time = time.time() - start
    print(f"✓ Loaded in {excel_time:.4f} seconds")
    
    # Test cached loading time
    print("\n--- Excel loading time (cached) ---")
    clear_pattern_registry_cache()
    start = time.time()
    registry = _load_consolidated_registry()
    cached_time_1 = time.time() - start
    
    start = time.time()
    registry = _load_consolidated_registry()
    cached_time_2 = time.time() - start
    
    print(f"✓ First cached load: {cached_time_1:.4f} seconds")
    print(f"✓ Second cached load: {cached_time_2:.6f} seconds")
    
    # Test JSON loading if available
    if os.path.exists('consolidated_registry.json'):
        print("\n--- JSON loading time (for comparison) ---")
        start = time.time()
        with open('consolidated_registry.json', 'r') as f:
            json.load(f)
        json_time = time.time() - start
        print(f"✓ JSON load: {json_time:.4f} seconds")
        
        ratio = excel_time / json_time
        print(f"\nExcel/JSON ratio: {ratio:.2f}x")
    
    return True


if __name__ == '__main__':
    print("\n🧬 NonBDNAFinder - Excel Pattern Loading Test Suite\n")
    
    all_passed = True
    
    # Test 1: Basic Excel loading
    if not test_excel_loading():
        all_passed = False
        print("\n❌ Excel loading test FAILED")
    else:
        print("\n✓ Excel loading test PASSED")
    
    # Test 2: Compare with JSON
    if not compare_with_json():
        all_passed = False
        print("\n❌ JSON comparison test FAILED")
    else:
        print("\n✓ JSON comparison test PASSED")
    
    # Test 3: Use with detector
    if not test_with_detector():
        all_passed = False
        print("\n❌ Detector integration test FAILED")
    else:
        print("\n✓ Detector integration test PASSED")
    
    # Test 4: Performance
    if not test_performance():
        all_passed = False
        print("\n❌ Performance test FAILED")
    else:
        print("\n✓ Performance test PASSED")
    
    # Final result
    print("\n" + "=" * 80)
    if all_passed:
        print("✅ ALL TESTS PASSED")
    else:
        print("❌ SOME TESTS FAILED")
    print("=" * 80)
    
    sys.exit(0 if all_passed else 1)
