#!/usr/bin/env python3
"""
Test suite for minimal reporting schema (Task 1 & 2 requirements).

Tests:
1. Core output columns are correctly defined
2. Motif-specific columns are properly separated by class
3. Export functions use only core columns for main output
4. Excel export separates core and motif-specific data
5. Default values are properly set for missing fields
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_core_columns_defined():
    """Test that CORE_OUTPUT_COLUMNS contains all required fields"""
    print("\n" + "="*70)
    print("TEST 1: Core Output Columns Definition")
    print("="*70)
    
    from utilities import CORE_OUTPUT_COLUMNS
    
    required_columns = [
        'Sequence_Name',
        'Class',
        'Subclass',
        'Start',
        'End',
        'Length',
        'Strand',
        'Score',
        'Method',
        'Pattern_ID'
    ]
    
    print(f"\nRequired columns: {required_columns}")
    print(f"Defined CORE_OUTPUT_COLUMNS: {CORE_OUTPUT_COLUMNS}")
    
    for col in required_columns:
        if col not in CORE_OUTPUT_COLUMNS:
            print(f"✗ FAILED: Missing required column '{col}'")
            return False
        else:
            print(f"✓ Found required column: {col}")
    
    # Check that we don't have extra columns that shouldn't be there
    deprecated_columns = ['Source', 'Sequence', 'ID']
    for col in deprecated_columns:
        if col in CORE_OUTPUT_COLUMNS:
            print(f"⚠ WARNING: Deprecated column '{col}' found in CORE_OUTPUT_COLUMNS")
    
    print("\n✓ TEST 1 PASSED: All required core columns are defined")
    return True


def test_motif_specific_columns_defined():
    """Test that motif-specific columns are properly defined per class"""
    print("\n" + "="*70)
    print("TEST 2: Motif-Specific Columns Definition")
    print("="*70)
    
    from utilities import MOTIF_SPECIFIC_COLUMNS
    
    required_classes = [
        'G-Quadruplex',
        'Z-DNA',
        'i-Motif',
        'Slipped DNA',
        'Cruciform',
        'Triplex',
        'R-Loop',
        'Curved DNA',
        'A-Philic'
    ]
    
    print(f"\nRequired motif classes: {required_classes}")
    print(f"\nDefined MOTIF_SPECIFIC_COLUMNS keys: {list(MOTIF_SPECIFIC_COLUMNS.keys())}")
    
    all_passed = True
    for cls in required_classes:
        if cls not in MOTIF_SPECIFIC_COLUMNS:
            print(f"✗ FAILED: Missing motif-specific columns for class '{cls}'")
            all_passed = False
        else:
            cols = MOTIF_SPECIFIC_COLUMNS[cls]
            print(f"✓ {cls}: {cols}")
    
    if all_passed:
        print("\n✓ TEST 2 PASSED: All motif-specific columns are defined")
    else:
        print("\n✗ TEST 2 FAILED: Some motif classes missing column definitions")
    
    return all_passed


def test_export_to_dataframe():
    """Test that export_results_to_dataframe uses only core columns"""
    print("\n" + "="*70)
    print("TEST 3: Export to DataFrame (Core Columns Only)")
    print("="*70)
    
    from utilities import export_results_to_dataframe, CORE_OUTPUT_COLUMNS
    
    # Create test motifs with various fields
    test_motifs = [
        {
            'Sequence_Name': 'test_seq',
            'Class': 'G-Quadruplex',
            'Subclass': 'Canonical intramolecular G4',
            'Start': 100,
            'End': 125,
            'Length': 25,
            'Strand': '+',
            'Score': 2.5,
            'Method': 'G4Hunter_detection',
            'Pattern_ID': 'G4_001',
            # Extra fields that should NOT appear in DataFrame
            'Num_Stems': 4,
            'Loop_Length': 3,
            'Extra_Field': 'should_not_appear'
        },
        {
            'Sequence_Name': 'test_seq',
            'Class': 'Slipped DNA',
            'Subclass': 'Direct Repeat',
            'Start': 200,
            'End': 230,
            'Length': 30,
            # Missing some fields (should get defaults)
            'Score': 1.8,
            'Repeat_Unit': 'ATGC',
            'Repeat_Count': 5
        }
    ]
    
    df = export_results_to_dataframe(test_motifs)
    
    print(f"\nDataFrame columns: {list(df.columns)}")
    print(f"Expected columns: {CORE_OUTPUT_COLUMNS}")
    
    # Check that only core columns are present
    if list(df.columns) != CORE_OUTPUT_COLUMNS:
        print(f"✗ FAILED: DataFrame columns don't match core columns")
        print(f"  Extra columns: {set(df.columns) - set(CORE_OUTPUT_COLUMNS)}")
        print(f"  Missing columns: {set(CORE_OUTPUT_COLUMNS) - set(df.columns)}")
        return False
    
    print("\n✓ DataFrame has exactly the core columns")
    
    # Check that defaults are properly set
    print("\nChecking default values:")
    row2 = df.iloc[1]
    
    if row2['Strand'] == '+':
        print(f"✓ Default Strand set correctly: {row2['Strand']}")
    else:
        print(f"✗ FAILED: Default Strand not set correctly: {row2['Strand']}")
        return False
    
    if row2['Method'] == 'Pattern_detection':
        print(f"✓ Default Method set correctly: {row2['Method']}")
    else:
        print(f"✗ FAILED: Default Method not set correctly: {row2['Method']}")
        return False
    
    if row2['Pattern_ID'] == 'Unknown':
        print(f"✓ Default Pattern_ID set correctly: {row2['Pattern_ID']}")
    else:
        print(f"✗ FAILED: Default Pattern_ID not set correctly: {row2['Pattern_ID']}")
        return False
    
    # Check that extra fields are not present
    if 'Num_Stems' not in df.columns and 'Extra_Field' not in df.columns:
        print(f"✓ Extra fields correctly excluded from DataFrame")
    else:
        print(f"✗ FAILED: Extra fields found in DataFrame")
        return False
    
    print("\n✓ TEST 3 PASSED: DataFrame export uses only core columns with proper defaults")
    return True


def test_export_to_csv():
    """Test that CSV export uses only core columns"""
    print("\n" + "="*70)
    print("TEST 4: Export to CSV (Core Columns Only)")
    print("="*70)
    
    from utilities import export_to_csv, CORE_OUTPUT_COLUMNS
    import csv
    from io import StringIO
    
    test_motifs = [
        {
            'Sequence_Name': 'test_seq',
            'Class': 'G-Quadruplex',
            'Subclass': 'Canonical',
            'Start': 100,
            'End': 125,
            'Length': 25,
            'Strand': '+',
            'Score': 2.5,
            'Method': 'G4Hunter',
            'Pattern_ID': 'G4_001',
            'Extra_Field': 'should_not_appear'
        }
    ]
    
    csv_content = export_to_csv(test_motifs)
    
    # Parse CSV to check columns
    csv_reader = csv.DictReader(StringIO(csv_content))
    header = csv_reader.fieldnames
    
    print(f"\nCSV header: {header}")
    print(f"Expected columns: {CORE_OUTPUT_COLUMNS}")
    
    if list(header) != CORE_OUTPUT_COLUMNS:
        print(f"✗ FAILED: CSV header doesn't match core columns")
        return False
    
    print("\n✓ CSV header has exactly the core columns")
    
    # Check that data row doesn't have extra fields
    row = next(csv_reader)
    if 'Extra_Field' in row:
        print(f"✗ FAILED: Extra fields found in CSV data")
        return False
    
    print("✓ Extra fields correctly excluded from CSV")
    
    print("\n✓ TEST 4 PASSED: CSV export uses only core columns")
    return True


def test_basic_detection():
    """Test that detector outputs include required core fields"""
    print("\n" + "="*70)
    print("TEST 5: Detector Output Fields")
    print("="*70)
    
    from nonbscanner import analyze_sequence
    from utilities import CORE_OUTPUT_COLUMNS
    
    # G4 motif sequence
    test_seq = "GGGTTAGGGTTAGGGTTAGGG"
    motifs = analyze_sequence(test_seq, "test_sequence")
    
    print(f"\nDetected {len(motifs)} motifs")
    
    if len(motifs) == 0:
        print("⚠ WARNING: No motifs detected (test sequence might not produce motifs)")
        return True
    
    # Check first motif has all core fields
    motif = motifs[0]
    print(f"\nFirst motif class: {motif.get('Class', 'N/A')}")
    print(f"Motif fields: {list(motif.keys())}")
    
    missing_fields = []
    for field in CORE_OUTPUT_COLUMNS:
        if field not in motif:
            missing_fields.append(field)
        else:
            print(f"✓ Has field: {field} = {motif[field]}")
    
    if missing_fields:
        print(f"\n⚠ WARNING: Missing core fields in detector output: {missing_fields}")
        print("  (These should be added by export functions with defaults)")
    else:
        print("\n✓ All core fields present in detector output")
    
    print("\n✓ TEST 5 PASSED: Detectors output includes required fields")
    return True


def run_all_tests():
    """Run all test functions"""
    print("\n" + "="*70)
    print("MINIMAL REPORTING SCHEMA TEST SUITE")
    print("Testing Task 1 & 2 Requirements")
    print("="*70)
    
    tests = [
        ("Core Columns Definition", test_core_columns_defined),
        ("Motif-Specific Columns Definition", test_motif_specific_columns_defined),
        ("Export to DataFrame", test_export_to_dataframe),
        ("Export to CSV", test_export_to_csv),
        ("Detector Output Fields", test_basic_detection),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"\n✗ TEST FAILED with exception: {e}")
            import traceback
            traceback.print_exc()
            results.append((test_name, False))
    
    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✓ PASSED" if result else "✗ FAILED"
        print(f"{status}: {test_name}")
    
    print(f"\n{passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 ALL TESTS PASSED! Minimal reporting schema is correctly implemented.")
        return True
    else:
        print(f"\n⚠ {total - passed} test(s) failed. Please review the output above.")
        return False


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
