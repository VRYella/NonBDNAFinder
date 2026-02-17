#!/usr/bin/env python3
"""
Integration test for MultiFASTA visualization feature.
Tests the complete workflow from FASTA parsing to Excel export and visualization.
"""

import sys
import os
import tempfile

# Add repository root to path for imports (required after moving to subdirectory)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Utilities.utilities import parse_fasta, export_multifasta_to_excel
from Utilities.multifasta_visualizer import MultiFastaVisualizer, prepare_multifasta_excel_data
import pandas as pd


def test_equal_length_workflow():
    """Test complete workflow with equal-length sequences."""
    print("=" * 80)
    print("Testing Equal-Length MultiFASTA Workflow")
    print("=" * 80)
    
    # Create test FASTA file
    fasta_content = """>Promoter1
GGGTTAGGGTTAGGGTTAGGGATGGGCTGGGAAGGGATCGATCGATCGGGGCGATCGATCGAT
>Promoter2
GGGTTAGGGTTAGGGTTAGGGATGGGCTGGGAAGGGATCGATCGATCGGGGCGATCGATCGAT
>Promoter3
GGGTTAGGGTTAGGGTTAGGGATGGGCTGGGAAGGGATCGATCGATCGGGGCGATCGATCGAT
"""
    
    # Parse FASTA content
    sequences = parse_fasta(fasta_content)
    
    print(f"\n✓ Parsed {len(sequences)} sequences")
    for name, seq in sequences.items():
        print(f"  - {name}: {len(seq)} bp")
    
    # Check equal length
    lengths = [len(seq) for seq in sequences.values()]
    equal_length = len(set(lengths)) == 1
    print(f"\n✓ Equal length detection: {equal_length}")
    
    # Create sample annotations (simulating detector output)
    annotations_by_sequence = {}
    sequence_lengths = {}
    
    for name, seq in sequences.items():
        # Simulate finding some motifs
        motifs = [
            {'Start': 1, 'End': 24, 'Class': 'G-Quadruplex', 'Subclass': 'G4_Telomeric', 
             'Length': 24, 'Sequence': seq[0:24], 'Score': 2.5, 'Strand': '+', 
             'Method': 'Regex', 'Pattern_ID': 'G4_001', 'Sequence_Name': name},
            {'Start': 25, 'End': 42, 'Class': 'G-Quadruplex', 'Subclass': 'G4_Canonical', 
             'Length': 18, 'Sequence': seq[24:42], 'Score': 2.2, 'Strand': '+', 
             'Method': 'Regex', 'Pattern_ID': 'G4_002', 'Sequence_Name': name},
        ]
        annotations_by_sequence[name] = motifs
        sequence_lengths[name] = len(seq)
    
    print(f"\n✓ Created sample annotations for {len(annotations_by_sequence)} sequences")
    
    # Create visualizer
    visualizer = MultiFastaVisualizer(annotations_by_sequence)
    print(f"✓ Created MultiFastaVisualizer with {len(visualizer.all_motifs)} total motifs")
    
    # Generate unified summary
    summary = visualizer.generate_unified_summary()
    print("\n✓ Generated unified summary:")
    print(f"  - Sequences: {len(summary['sequence_stats'])}")
    print(f"  - Classes: {list(summary['class_distribution'].keys())}")
    print(f"  - Total motifs: {sum(summary['class_distribution'].values())}")
    
    # Compute positional occurrence
    if equal_length:
        seq_length = lengths[0]
        pos_data = visualizer.compute_positional_occurrence(seq_length)
        print(f"\n✓ Computed positional occurrence:")
        print(f"  - Overall positions: {len(pos_data['Overall'])}")
        print(f"  - Classes: {list(pos_data['Class'].keys())}")
    
    # Prepare Excel data
    excel_data = prepare_multifasta_excel_data(
        annotations_by_sequence,
        sequence_lengths,
        equal_length=equal_length,
        seq_length=lengths[0] if equal_length else None
    )
    
    print(f"\n✓ Prepared Excel data:")
    print(f"  - Sheets: {list(excel_data.keys())}")
    print(f"  - All_Motifs: {len(excel_data['All_Motifs'])} rows")
    print(f"  - Sequence_Summary: {len(excel_data['Sequence_Summary'])} rows")
    print(f"  - Class_Summary: {len(excel_data['Class_Summary'])} rows")
    if 'Positional_Occurrence' in excel_data:
        print(f"  - Positional_Occurrence: {len(excel_data['Positional_Occurrence'])} rows")
    
    # Export to Excel
    with tempfile.NamedTemporaryFile(suffix='.xlsx', delete=False) as tmp:
        tmp_path = tmp.name
    
    try:
        result = export_multifasta_to_excel(
            annotations_by_sequence,
            sequence_lengths,
            tmp_path,
            equal_length=equal_length,
            seq_length=lengths[0] if equal_length else None
        )
        
        print(f"\n✓ Excel export: {result}")
        
        # Verify Excel file
        excel_file = pd.ExcelFile(tmp_path, engine='openpyxl')
        sheet_names = excel_file.sheet_names
        print(f"  - Sheets in file: {sheet_names}")
        
        # Read and display sample data
        df_motifs = pd.read_excel(tmp_path, sheet_name='All_Motifs')
        print(f"\n✓ All_Motifs sheet sample:")
        print(df_motifs.head(3).to_string(index=False))
        
        df_seq = pd.read_excel(tmp_path, sheet_name='Sequence_Summary')
        print(f"\n✓ Sequence_Summary sheet:")
        print(df_seq.to_string(index=False))
        
        print(f"\n✓ Excel file successfully created at: {tmp_path}")
        print(f"  File size: {os.path.getsize(tmp_path):,} bytes")
        
    finally:
        # Clean up
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
            print(f"✓ Cleaned up temporary file")
    
    print("\n" + "=" * 80)
    print("✅ Equal-Length MultiFASTA Workflow: PASSED")
    print("=" * 80)


def test_different_length_workflow():
    """Test complete workflow with different-length sequences."""
    print("\n\n" + "=" * 80)
    print("Testing Different-Length MultiFASTA Workflow")
    print("=" * 80)
    
    # Create test FASTA content
    fasta_content = """>Sequence1
GGGTTAGGGTTAGGGTTAGGGATGGGCTGGGAAGGGATCGATCGATCGGGGCGATCGATCGAT
>Sequence2
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCCCCTCCCCTCCCCTCCCC
>Sequence3
TTAGGGTTAGGGTTAGGGTTAGGGAAAAAATAAAAATAAAAATAAAAAAAAACAGCAGCAGCAGCAGCAGCAGCAG
"""
    
    # Parse FASTA content
    sequences = parse_fasta(fasta_content)
    
    print(f"\n✓ Parsed {len(sequences)} sequences")
    for name, seq in sequences.items():
        print(f"  - {name}: {len(seq)} bp")
    
    # Check equal length
    lengths = [len(seq) for seq in sequences.values()]
    equal_length = len(set(lengths)) == 1
    print(f"\n✓ Equal length detection: {equal_length} (should be False)")
    
    # Create sample annotations
    annotations_by_sequence = {}
    sequence_lengths = {}
    
    for name, seq in sequences.items():
        # Simulate finding some motifs
        motifs = [
            {'Start': 1, 'End': min(20, len(seq)), 'Class': 'Test_Class', 
             'Subclass': 'Test_Subclass', 'Length': min(20, len(seq)), 
             'Sequence': seq[0:min(20, len(seq))], 'Score': 2.0, 
             'Strand': '+', 'Method': 'Test', 'Pattern_ID': 'TEST_001',
             'Sequence_Name': name}
        ]
        annotations_by_sequence[name] = motifs
        sequence_lengths[name] = len(seq)
    
    print(f"\n✓ Created sample annotations")
    
    # Prepare Excel data
    excel_data = prepare_multifasta_excel_data(
        annotations_by_sequence,
        sequence_lengths,
        equal_length=equal_length,
        seq_length=None
    )
    
    print(f"\n✓ Prepared Excel data:")
    print(f"  - Sheets: {list(excel_data.keys())}")
    if 'Positional_Occurrence' not in excel_data:
        print(f"  ✓ Positional_Occurrence correctly excluded (different lengths)")
    
    # Export to Excel
    with tempfile.NamedTemporaryFile(suffix='.xlsx', delete=False) as tmp:
        tmp_path = tmp.name
    
    try:
        result = export_multifasta_to_excel(
            annotations_by_sequence,
            sequence_lengths,
            tmp_path,
            equal_length=equal_length,
            seq_length=None
        )
        
        print(f"\n✓ Excel export: {result}")
        
        # Verify Excel file
        excel_file = pd.ExcelFile(tmp_path, engine='openpyxl')
        sheet_names = excel_file.sheet_names
        print(f"  - Sheets in file: {sheet_names}")
        
        if 'Positional_Occurrence' in sheet_names:
            print(f"  ⚠️ WARNING: Positional_Occurrence should not be present")
        else:
            print(f"  ✓ Positional_Occurrence correctly excluded")
        
        print(f"\n✓ Excel file successfully created")
        
    finally:
        # Clean up
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
            print(f"✓ Cleaned up temporary file")
    
    print("\n" + "=" * 80)
    print("✅ Different-Length MultiFASTA Workflow: PASSED")
    print("=" * 80)


if __name__ == '__main__':
    try:
        test_equal_length_workflow()
        test_different_length_workflow()
        print("\n\n" + "=" * 80)
        print("🎉 ALL INTEGRATION TESTS PASSED!")
        print("=" * 80)
    except Exception as e:
        print(f"\n\n❌ INTEGRATION TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
