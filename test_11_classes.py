#!/usr/bin/env python3
"""
Test script to validate all 11 Non-B DNA classes are working properly.

This script tests each of the 11 official Non-B DNA classes according to:
Class 1:  Curved DNA
Class 2:  Slipped DNA  
Class 3:  Cruciform DNA
Class 4:  R-loop
Class 5:  Triplex
Class 6:  G-Quadruplex Family
Class 7:  i-motif family
Class 8:  Z-DNA
Class 9:  A-philic DNA   
Class 10: Hybrid               
Class 11: Non-B DNA clusters
"""

from hyperscan_integration import all_motifs_refactored
from motif_classification import OFFICIAL_CLASSIFICATION
import sys

def test_sequence_for_classes(sequence, description, expected_classes=None):
    """Test a sequence and report which classes are found"""
    print(f"\n=== Testing {description} ===")
    print(f"Sequence length: {len(sequence)} bp")
    print(f"Sequence: {sequence[:100]}{'...' if len(sequence) > 100 else ''}")
    
    try:
        results = all_motifs_refactored(sequence)
        print(f"Found {len(results)} motifs")
        
        # Count by class
        class_counts = {}
        for motif in results:
            class_name = motif.get('Class', 'Unknown')
            class_counts[class_name] = class_counts.get(class_name, 0) + 1
        
        print("Classes detected:")
        for class_name, count in sorted(class_counts.items()):
            print(f"  {class_name}: {count} motifs")
        
        if expected_classes:
            found_expected = set(class_counts.keys()) & set(expected_classes)
            missing_expected = set(expected_classes) - set(class_counts.keys())
            if missing_expected:
                print(f"  Missing expected classes: {missing_expected}")
            if found_expected:
                print(f"  Found expected classes: {found_expected}")
        
        return class_counts
        
    except Exception as e:
        print(f"Error analyzing sequence: {e}")
        return {}

def main():
    """Run comprehensive tests for all 11 Non-B DNA classes"""
    print("=" * 60)
    print("COMPREHENSIVE TEST OF 11 NON-B DNA CLASSES")
    print("=" * 60)
    
    # Print official classification
    print("\nOfficial 11 Non-B DNA Classes:")
    for class_id, info in OFFICIAL_CLASSIFICATION.items():
        print(f"Class {class_id}: {info['class_name']} ({len(info['subclasses'])} subclasses)")
    
    all_found_classes = set()
    
    # Test sequence 1: Rich in various motifs  
    seq1 = (
        "ATCGATCGATCG"          # Basic DNA
        "AAAATTTTATTTAAATTTAAA"  # Curved DNA potential (A-tracts)
        "GGTTAGGGTTAGGGTTAGGG"   # G-quadruplex potential
        "CCCCCTCCCCCTCCCCCTCCCC" # i-motif potential
        "CGCGCGCGCGCGCGCGCGCG"   # Z-DNA potential  
        "CGGGGGGGGGGTGGCCCCCCCC" # A-philic specific patterns
        "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"  # R-loop potential (long G-tract)
        "AAGGAAGGAAGGAAGGAAGG"   # Triplex potential
        "CGGCGGCGGCGGCGGCGGCGG"   # eGZ potential
        "CTCTCTCTCTCTCTCTCTCT"    # Slipped DNA potential
        "AATTCCGGAATTCCGG"        # Cruciform potential
    )
    
    results1 = test_sequence_for_classes(
        seq1, 
        "Multi-motif rich sequence",
        ["Curved DNA", "G-Quadruplex Family", "i-motif family", "Z-DNA", 
         "A-philic DNA", "Triplex", "Slipped DNA", "Cruciform DNA", "R-loop"]
    )
    all_found_classes.update(results1.keys())
    
    # Test sequence 2: Focus on specific classes
    seq2 = (
        "AAAAAAAAAAAAAAAAAAAAAA"    # A-tracts for Curved DNA
        "TTTTTTTTTTTTTTTTTTTTT"     # T-tracts 
        "CGGGGGGGGGGGGGGGGGGG"      # Long G-tract for R-loop
        "CCCCCCCCCCCCCCCCCCCC"      # Long C-tract for i-motif
        "ATATATATATATATATATATAT"    # AT repeats for Z-DNA
        "GAGAGAGAGAGAGAGAGAGA"      # GA repeats for slipped
        "GAATTCGAATTCGAATTC"        # Palindromes for cruciform
    )
    
    results2 = test_sequence_for_classes(
        seq2,
        "Specific class targeting sequence", 
        ["Curved DNA", "R-loop", "i-motif family", "Z-DNA", "Slipped DNA", "Cruciform DNA"]
    )
    all_found_classes.update(results2.keys())
    
    # Test sequence 3: Comprehensive coverage
    seq3 = (
        "ATCGATCGATCG" +
        "CGGGGGGGGGGTGGCCCCCCCC" +  # A-philic DNA (high-scoring patterns)
        "GGGTTTTGGGTTTTGGGTTTTGGG" +  # G-quadruplex
        "CCCCTTTTCCCCTTTTCCCCTTTT" +  # i-motif
        "CGCGCGCGCGCGCGCG" +          # Z-DNA
        "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" +     # R-loop (longer G-tract) 
        "GAATTCGAATTCGAATTC" +        # Cruciform
        "AAAGGAAAGGAAAGGAAAGG" +      # Triplex
        "CTCTCTCTCTCTCTCTCTCT" +      # Slipped DNA
        "AAATAAATAAATAAATAAATAAAT"    # Curved DNA
    )
    
    results3 = test_sequence_for_classes(
        seq3,
        "Comprehensive coverage sequence",
        ["A-philic DNA", "G-Quadruplex Family", "i-motif family", "Z-DNA", 
         "R-loop", "Cruciform DNA", "Triplex", "Slipped DNA", "Curved DNA"]
    )
    all_found_classes.update(results3.keys())
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    expected_classes = [info['class_name'] for info in OFFICIAL_CLASSIFICATION.values()]
    
    print(f"\nExpected 11 classes: {expected_classes}")
    print(f"\nFound classes across all tests: {sorted(all_found_classes)}")
    
    missing_classes = set(expected_classes) - all_found_classes
    if missing_classes:
        print(f"\n⚠️  Missing classes: {missing_classes}")
    else:
        print(f"\n✅ All 11 classes detected successfully!")
    
    found_extra = all_found_classes - set(expected_classes)
    if found_extra:
        print(f"\n📝 Additional classes found: {found_extra}")
    
    print(f"\nTotal classes found: {len(all_found_classes)}/11")
    
    # Detailed breakdown
    if missing_classes:
        print(f"\n🔍 Classes that need attention:")
        for missing in missing_classes:
            class_id = next(k for k, v in OFFICIAL_CLASSIFICATION.items() if v['class_name'] == missing)
            subclasses = OFFICIAL_CLASSIFICATION[class_id]['subclasses'] 
            print(f"  Class {class_id}: {missing} (subclasses: {subclasses})")
    
    return len(missing_classes) == 0

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)