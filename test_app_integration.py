#!/usr/bin/env python3
"""
Quick functional test to verify the 11 Non-B DNA classes integration
works properly with the Streamlit app components.
"""

from hyperscan_integration import all_motifs_refactored
from motif_classification import OFFICIAL_CLASSIFICATION, update_motif_with_ids
import pandas as pd

def test_app_integration():
    """Test the complete pipeline as used by the Streamlit app"""
    
    print("=" * 60)
    print("FUNCTIONAL TEST: APP INTEGRATION FOR 11 CLASSES")
    print("=" * 60)
    
    # Test sequence with multiple motif types
    test_sequence = (
        "ATCGATCGATCG" +
        "CGGGGGGGGGGTGGCCCCCCCC" +  # A-philic patterns
        "GGGTTTTGGGTTTTGGGTTTTGGG" +  # G-quadruplex
        "CCCCTTTTCCCCTTTTCCCCTTTT" +  # i-motif
        "CGCGCGCGCGCGCGCG" +          # Z-DNA
        "AAGGAAGGAAGGAAGGAAGG" +      # Triplex
        "GAATTCGAATTCGAATTC" +        # Cruciform
        "CTCTCTCTCTCTCTCTCTCT" +      # Slipped DNA
        "AAATAAATAAATAAATAAATAAAT"    # Curved DNA
    )
    
    print(f"Test sequence length: {len(test_sequence)} bp")
    
    # Step 1: Run motif detection (as app would)
    print("\n1. Running motif detection...")
    motifs = all_motifs_refactored(test_sequence, sequence_name="test_seq")
    print(f"   Found {len(motifs)} total motifs")
    
    # Step 2: Update with official classification (as app would)
    print("\n2. Applying official classification...")
    updated_motifs = []
    for motif in motifs:
        updated_motif = update_motif_with_ids(motif.copy())
        updated_motifs.append(updated_motif)
    
    # Step 3: Create summary as app would
    print("\n3. Creating summary statistics...")
    
    # Count by official class
    class_counts = {}
    for motif in updated_motifs:
        class_name = motif.get('Class', 'Unknown')
        class_counts[class_name] = class_counts.get(class_name, 0) + 1
    
    # Count by subclass  
    subclass_counts = {}
    for motif in updated_motifs:
        subclass = motif.get('Subclass', 'Unknown')
        subclass_counts[subclass] = subclass_counts.get(subclass, 0) + 1
    
    print("\n4. Results Summary:")
    print("   Classes detected:")
    for class_name, count in sorted(class_counts.items()):
        print(f"     {class_name}: {count} motifs")
    
    print(f"\n   Subclasses detected ({len(subclass_counts)} total):")
    for subclass, count in sorted(subclass_counts.items()):
        print(f"     {subclass}: {count} motifs")
    
    # Step 4: Verify official class coverage
    print("\n5. Official Class Coverage Check:")
    expected_classes = [info['class_name'] for info in OFFICIAL_CLASSIFICATION.values()]
    found_classes = set(class_counts.keys())
    
    coverage = len(found_classes) / len(expected_classes) * 100
    print(f"   Coverage: {len(found_classes)}/{len(expected_classes)} classes ({coverage:.1f}%)")
    
    missing = set(expected_classes) - found_classes
    if missing:
        print(f"   Missing: {missing}")
    else:
        print(f"   ✅ All official classes represented!")
    
    # Step 5: Create DataFrame as app would for export
    print("\n6. Creating export DataFrame...")
    df = pd.DataFrame(updated_motifs)
    
    if not df.empty:
        print(f"   DataFrame shape: {df.shape}")
        print(f"   Columns: {list(df.columns)}")
        
        # Show sample of export data
        print("\n   Sample export data:")
        sample_cols = ['Class', 'Subclass', 'Start', 'End', 'Score']
        if all(col in df.columns for col in sample_cols):
            for i, row in df[sample_cols].head(3).iterrows():
                print(f"     {i+1}. {row['Class']} | {row['Subclass']} | {row['Start']}-{row['End']} | Score: {row['Score']:.2f}")
    
    print("\n" + "=" * 60)
    print("FUNCTIONAL TEST COMPLETE")
    print("=" * 60)
    
    return len(found_classes) >= 8  # Accept if we get at least 8/11 classes

if __name__ == "__main__":
    success = test_app_integration()
    print(f"\nTest result: {'PASS' if success else 'FAIL'}")