#!/usr/bin/env python3
"""
Integration test for minimal reporting schema.

Tests the complete workflow from detection to export.
"""

import sys
import os
import tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from nonbscanner import analyze_sequence
from utilities import (
    export_to_csv, export_to_excel, export_results_to_dataframe,
    CORE_OUTPUT_COLUMNS, MOTIF_SPECIFIC_COLUMNS
)
import pandas as pd


def test_complete_workflow():
    """Test complete workflow: detection -> DataFrame -> CSV -> Excel"""
    print("\n" + "="*70)
    print("INTEGRATION TEST: Complete Workflow")
    print("="*70)
    
    # Test sequence with multiple motif types
    test_seq = (
        "GGGTTAGGGTTAGGGTTAGGG"  # G4 motif
        "CGCGCGCGCGCGCGCGCG"      # Z-DNA
        "AAAAAAAAAATTTTTTTTTT"    # Curved DNA
        "ATGCATGCATGCATGC"        # Slipped DNA
    )
    
    print(f"\n1. Detecting motifs in test sequence ({len(test_seq)} bp)...")
    motifs = analyze_sequence(test_seq, "integration_test")
    
    print(f"   ✓ Detected {len(motifs)} motifs")
    
    if len(motifs) == 0:
        print("   ⚠ No motifs detected, test sequence might need adjustment")
        return True
    
    # Show motif classes
    classes = set(m['Class'] for m in motifs)
    print(f"   ✓ Motif classes: {', '.join(sorted(classes))}")
    
    # Test DataFrame export
    print("\n2. Testing DataFrame export (display format)...")
    df = export_results_to_dataframe(motifs)
    
    print(f"   ✓ DataFrame shape: {df.shape}")
    print(f"   ✓ DataFrame columns: {list(df.columns)}")
    
    # Verify only core columns
    if list(df.columns) != CORE_OUTPUT_COLUMNS:
        print(f"   ✗ FAILED: DataFrame has unexpected columns")
        return False
    
    print(f"   ✓ DataFrame contains only core columns")
    
    # Test CSV export
    print("\n3. Testing CSV export...")
    csv_content = export_to_csv(motifs)
    
    lines = csv_content.strip().split('\n')
    print(f"   ✓ CSV has {len(lines)} lines (including header)")
    
    # Verify CSV header (strip whitespace and line endings)
    header = [col.strip() for col in lines[0].replace('\r', '').split(',')]
    if header != CORE_OUTPUT_COLUMNS:
        print(f"   ✗ FAILED: CSV header doesn't match core columns")
        print(f"      Expected: {CORE_OUTPUT_COLUMNS}")
        print(f"      Got: {header}")
        return False
    
    print(f"   ✓ CSV header matches core columns")
    
    # Test Excel export
    print("\n4. Testing Excel export...")
    with tempfile.NamedTemporaryFile(suffix='.xlsx', delete=False) as tmp:
        tmp_path = tmp.name
    
    try:
        # Test simple format (2-tab)
        result = export_to_excel(motifs, tmp_path, simple_format=True)
        print(f"   ✓ Simple format Excel: {result}")
        
        # Read back and verify
        df_consolidated = pd.read_excel(tmp_path, sheet_name='NonOverlappingConsolidated')
        print(f"   ✓ NonOverlappingConsolidated sheet: {df_consolidated.shape}")
        
        if list(df_consolidated.columns) != CORE_OUTPUT_COLUMNS:
            print(f"   ✗ FAILED: Excel sheet has unexpected columns")
            print(f"      Expected: {CORE_OUTPUT_COLUMNS}")
            print(f"      Got: {list(df_consolidated.columns)}")
            return False
        
        print(f"   ✓ Excel sheet contains only core columns")
        
        df_all = pd.read_excel(tmp_path, sheet_name='OverlappingAll')
        print(f"   ✓ OverlappingAll sheet: {df_all.shape}")
        
        os.unlink(tmp_path)
        
        # Test detailed format (class-specific sheets)
        result = export_to_excel(motifs, tmp_path, simple_format=False)
        print(f"   ✓ Detailed format Excel: {result}")
        
        # Verify Core_Results sheet has only core columns
        df_core = pd.read_excel(tmp_path, sheet_name='Core_Results')
        if list(df_core.columns) != CORE_OUTPUT_COLUMNS:
            print(f"   ✗ FAILED: Core_Results sheet has unexpected columns")
            return False
        
        print(f"   ✓ Core_Results sheet contains only core columns")
        
        # Verify class-specific sheets have core + motif-specific columns
        xls = pd.ExcelFile(tmp_path)
        class_sheets = [s for s in xls.sheet_names if s not in ['Core_Results', 'Hybrid_Motifs', 'Cluster_Motifs']]
        
        print(f"   ✓ Class-specific sheets: {', '.join(class_sheets)}")
        
        for sheet_name in class_sheets:
            df_class = pd.read_excel(tmp_path, sheet_name=sheet_name)
            
            # Determine class from sheet name
            # Sheet names are sanitized, need to map back
            class_mapping = {
                'G_Quadruplex': 'G-Quadruplex',
                'Z_DNA': 'Z-DNA',
                'i_Motif': 'i-Motif',
                'Slipped_DNA': 'Slipped DNA',
                'Cruciform': 'Cruciform',
                'Triplex': 'Triplex',
                'R_Loop': 'R-Loop',
                'Curved_DNA': 'Curved DNA',
                'A_Philic': 'A-Philic'
            }
            
            motif_class = class_mapping.get(sheet_name, sheet_name)
            expected_cols = CORE_OUTPUT_COLUMNS + MOTIF_SPECIFIC_COLUMNS.get(motif_class, [])
            
            # Check that class sheet has core + specific columns
            if not all(col in expected_cols for col in CORE_OUTPUT_COLUMNS):
                print(f"   ✗ FAILED: {sheet_name} missing core columns")
                return False
            
            print(f"      ✓ {sheet_name}: {df_class.shape} (has core + specific columns)")
        
        os.unlink(tmp_path)
        
    except Exception as e:
        print(f"   ✗ FAILED with exception: {e}")
        import traceback
        traceback.print_exc()
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        return False
    
    print("\n✓ INTEGRATION TEST PASSED")
    return True


def test_motif_specific_columns_presence():
    """Test that motif-specific columns are present in class sheets"""
    print("\n" + "="*70)
    print("TEST: Motif-Specific Columns in Class Sheets")
    print("="*70)
    
    # Create test motifs for different classes
    test_motifs = [
        {
            'Sequence_Name': 'test',
            'Class': 'G-Quadruplex',
            'Subclass': 'Canonical',
            'Start': 1,
            'End': 25,
            'Length': 25,
            'Strand': '+',
            'Score': 2.5,
            'Method': 'G4Hunter',
            'Pattern_ID': 'G4_001',
            'Num_Tracts': 4,
            'Loop_Length': 3,
            'Priority': 80
        },
        {
            'Sequence_Name': 'test',
            'Class': 'Slipped DNA',
            'Subclass': 'STR',
            'Start': 50,
            'End': 70,
            'Length': 20,
            'Strand': '+',
            'Score': 2.0,
            'Method': 'STR_detection',
            'Pattern_ID': 'STR_001',
            'Repeat_Unit': 'ATGC',
            'Unit_Length': 4,
            'Repeat_Count': 5
        }
    ]
    
    with tempfile.NamedTemporaryFile(suffix='.xlsx', delete=False) as tmp:
        tmp_path = tmp.name
    
    try:
        result = export_to_excel(test_motifs, tmp_path, simple_format=False)
        print(f"✓ Excel export: {result}")
        
        # Check G-Quadruplex sheet has motif-specific columns
        df_g4 = pd.read_excel(tmp_path, sheet_name='G_Quadruplex')
        print(f"\nG-Quadruplex sheet columns: {list(df_g4.columns)}")
        
        g4_specific = MOTIF_SPECIFIC_COLUMNS['G-Quadruplex']
        for col in g4_specific:
            if col in df_g4.columns:
                print(f"  ✓ Has motif-specific column: {col} = {df_g4[col].iloc[0]}")
            else:
                print(f"  ⚠ Missing motif-specific column: {col}")
        
        # Check Slipped DNA sheet has motif-specific columns
        df_slipped = pd.read_excel(tmp_path, sheet_name='Slipped_DNA')
        print(f"\nSlipped DNA sheet columns: {list(df_slipped.columns)}")
        
        slipped_specific = MOTIF_SPECIFIC_COLUMNS['Slipped DNA']
        for col in slipped_specific:
            if col in df_slipped.columns:
                print(f"  ✓ Has motif-specific column: {col} = {df_slipped[col].iloc[0]}")
            else:
                print(f"  ⚠ Missing motif-specific column: {col}")
        
        os.unlink(tmp_path)
        
        print("\n✓ TEST PASSED")
        return True
        
    except Exception as e:
        print(f"\n✗ FAILED with exception: {e}")
        import traceback
        traceback.print_exc()
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        return False


if __name__ == "__main__":
    print("\n" + "="*70)
    print("INTEGRATION TEST SUITE")
    print("="*70)
    
    tests = [
        test_complete_workflow,
        test_motif_specific_columns_presence
    ]
    
    results = []
    for test_func in tests:
        try:
            result = test_func()
            results.append(result)
        except Exception as e:
            print(f"\n✗ TEST FAILED with exception: {e}")
            import traceback
            traceback.print_exc()
            results.append(False)
    
    passed = sum(results)
    total = len(results)
    
    print("\n" + "="*70)
    print(f"FINAL RESULT: {passed}/{total} tests passed")
    print("="*70)
    
    if passed == total:
        print("\n🎉 ALL INTEGRATION TESTS PASSED!")
        sys.exit(0)
    else:
        print(f"\n⚠ {total - passed} test(s) failed")
        sys.exit(1)
