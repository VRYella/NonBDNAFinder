#!/usr/bin/env python3
"""
Test script to verify eGZ motif scoring matches expected behavior
and validates motif classification system integration.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from motif_detectors import ZDNADetector
from motif_classification import get_motif_id, CURRENT_TO_OFFICIAL
import re


def test_scoring_patterns():
    """Test scoring patterns and verify Kadane algorithm raw scores."""
    print("Testing scoring patterns and Kadane raw scores...")
    
    detector = ZDNADetector(use_kadane=True)
    
    # Test sequences designed to have different scores
    test_cases = [
        ("CGGCGGCGGCGGCGGCGGCGGCGG", "High CGG content - should be eGZ"),
        ("CGCGCGCGCGCGCGCGCGCGCGCG", "High CG content - should be Z-DNA"),
        ("CGGCGCGGCGCGGCGCGGCGCGG", "Mixed CGG/CG - should be eGZ"),
        ("ATATATATATATATATATATAT", "AT repeats - should be Z-DNA"),
    ]
    
    print("\nSequence Classification and Scoring:")
    print("-" * 80)
    
    for seq, description in test_cases:
        candidates = detector.detect(seq, "test", "test", 0)
        
        print(f"\nSequence: {seq}")
        print(f"Description: {description}")
        
        if candidates:
            best = max(candidates, key=lambda c: c.raw_score if c.raw_score else 0)
            
            # Count CGG vs CG patterns
            cgg_count = len(re.findall(r'CGG', seq))
            cg_count = len(re.findall(r'CG', seq))
            
            print(f"CGG repeats: {cgg_count}, CG dinucleotides: {cg_count}")
            print(f"Detected as: {best.subclass}")
            print(f"Raw Kadane score: {best.raw_score:.2f}")
            print(f"Scoring method: {best.scoring_method}")
            
            # Verify scoring method contains "Kadane" 
            if "Kadane" in best.scoring_method:
                print("✓ Using Kadane algorithm for raw scoring")
            else:
                print("✗ Not using Kadane algorithm")
                
        else:
            print("No candidates detected")
    
    print("\n" + "="*80)


def test_classification_system():
    """Test integration with motif classification system."""
    print("Testing motif classification system integration...")
    
    # Test that eGZ maps to correct official classification
    official_name = CURRENT_TO_OFFICIAL.get("eGZ")
    motif_id = get_motif_id("eGZ (Extruded-G) DNA")
    
    print(f"\neGZ subclass mapping:")
    print(f"  Current name: eGZ")
    print(f"  Official name: {official_name}")
    print(f"  Motif ID: {motif_id}")
    
    if official_name == "eGZ (Extruded-G) DNA":
        print("✓ eGZ classification mapping is correct")
    else:
        print("✗ eGZ classification mapping is incorrect")
    
    print("\n" + "="*80)


def test_raw_score_preservation():
    """Test that Kadane algorithm raw scores are preserved correctly."""
    print("Testing raw score preservation...")
    
    detector = ZDNADetector(use_kadane=True)
    
    # Test with a sequence that should get eGZ classification
    seq = "CGGCGGCGGCGGCGGCGG"
    candidates = detector.detect(seq, "test", "test", 0)
    
    print(f"\nTest sequence: {seq}")
    
    if candidates:
        candidate = candidates[0]
        print(f"Raw score from Kadane: {candidate.raw_score}")
        print(f"Scoring method: {candidate.scoring_method}")
        
        # Verify the raw score is a float (not normalized)
        if isinstance(candidate.raw_score, (int, float)) and candidate.raw_score > 0:
            print("✓ Raw score is preserved as numeric value")
        else:
            print("✗ Raw score not preserved correctly")
            
        # Verify the scoring method indicates it's from Kadane
        if "Kadane" in candidate.scoring_method:
            print("✓ Scoring method correctly indicates Kadane algorithm")
        else:
            print("✗ Scoring method does not indicate Kadane algorithm")
            
    else:
        print("No candidates detected for test sequence")
        
    print("\n" + "="*80)


def main():
    """Run all verification tests."""
    print("=" * 80)
    print("eGZ MOTIF SYSTEM VERIFICATION")
    print("=" * 80)
    
    test_scoring_patterns()
    test_classification_system()
    test_raw_score_preservation()
    
    print("System verification completed!")
    print("=" * 80)


if __name__ == "__main__":
    main()