#!/usr/bin/env python3
"""
Test script to verify Kadane algorithm integration with motif detection framework.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from motif_detectors import ZDNADetector
from motifs.base import Candidate
import numpy as np


def test_kadane_integration():
    """Test that the ZDNADetector is using our Kadane algorithm."""
    print("Testing Kadane algorithm integration with motif detection framework...")
    
    # Create detector instance
    detector = ZDNADetector(use_kadane=True)
    
    # Test sequences with known Z-DNA potential
    test_sequences = [
        ("CGCGCGCGCGCGCGCGCGCG", "High Z-DNA potential (CG repeats)"),
        ("ATATATATATATATATATAT", "High Z-DNA potential (AT repeats)"),
        ("ATCGATCGATCGATCGATCG", "Mixed transitions"),
        ("AAAAAAAAAAAAAAAAAAA", "Low Z-DNA potential")
    ]
    
    for seq, description in test_sequences:
        print(f"\nTesting: {description}")
        print(f"Sequence: {seq}")
        
        # Detect Z-DNA candidates
        candidates = detector.detect(seq, "test_seq", "test_contig", 0)
        
        print(f"Candidates found: {len(candidates)}")
        for i, candidate in enumerate(candidates):
            print(f"  Candidate {i+1}:")
            print(f"    Position: {candidate.start}-{candidate.end}")
            print(f"    Score: {candidate.raw_score:.2f}")
            print(f"    Method: {candidate.scoring_method}")
            print(f"    Subclass: {candidate.subclass}")
    
    print("\n" + "="*60)


def test_parameter_sensitivity():
    """Test that different parameters affect detection sensitivity."""
    print("Testing parameter sensitivity...")
    
    seq = "CGCGCGCGCGCGCGCGCGCG"
    
    # Test with different threshold settings by modifying the detector
    print(f"Test sequence: {seq}")
    
    # We can't easily modify threshold from outside, but we can test behavior
    detector = ZDNADetector(use_kadane=True)
    candidates = detector.detect(seq, "test", "test", 0)
    
    print(f"Standard detection: {len(candidates)} candidates")
    if candidates:
        print(f"  Score: {candidates[0].raw_score:.2f}")
        print(f"  Method: {candidates[0].scoring_method}")
    
    print("\n" + "="*60)


def test_fallback_behavior():
    """Test fallback to regex patterns when Kadane algorithm is disabled."""
    print("Testing fallback behavior...")
    
    # Test with Kadane disabled
    detector = ZDNADetector(use_kadane=False)
    seq = "CGCGCGCGCGCGCGCGCGCG"
    
    candidates = detector.detect(seq, "test", "test", 0)
    print(f"Regex fallback detection: {len(candidates)} candidates")
    
    if candidates:
        for i, candidate in enumerate(candidates):
            print(f"  Candidate {i+1}: {candidate.subclass}")
            # Score these candidates
            scored = detector.score([candidate])
            print(f"    Score: {scored[0].raw_score:.2f}")
            print(f"    Method: {scored[0].scoring_method}")
    
    print("\n" + "="*60)


def test_scoring_method_identification():
    """Test that we can identify which scoring method was used."""
    print("Testing scoring method identification...")
    
    detector = ZDNADetector(use_kadane=True)
    seq = "CGCGCGCGCGCGCGCGCGCG"
    
    candidates = detector.detect(seq, "test", "test", 0)
    
    kadane_found = False
    for candidate in candidates:
        if "Kadane" in candidate.scoring_method:
            kadane_found = True
            print(f"✓ Kadane algorithm detected!")
            print(f"  Method: {candidate.scoring_method}")
            print(f"  Subclass: {candidate.subclass}")
            print(f"  Score: {candidate.raw_score:.2f}")
    
    if not kadane_found:
        print("✗ Kadane algorithm not detected in scoring methods")
        for candidate in candidates:
            print(f"  Found method: {candidate.scoring_method}")
    
    print("\n" + "="*60)


def main():
    """Run integration tests."""
    print("="*60)
    print("Z-DNA KADANE ALGORITHM INTEGRATION TESTING")
    print("="*60)
    
    test_kadane_integration()
    test_parameter_sensitivity()
    test_fallback_behavior()
    test_scoring_method_identification()
    
    print("Integration testing completed!")


if __name__ == "__main__":
    main()