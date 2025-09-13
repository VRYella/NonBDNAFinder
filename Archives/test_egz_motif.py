#!/usr/bin/env python3
"""
Test script for eGZ (Extruded-G) DNA motif detection.

This script validates that eGZ motifs formed by CGG trinucleotide repeats
are properly detected and scored using the Kadane algorithm.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from motif_detectors import ZDNADetector
from motifs.base import Candidate
import numpy as np


def test_egz_detection():
    """Test eGZ motif detection using example sequences from problem statement."""
    print("Testing eGZ motif detection...")
    
    # Example sequences and expected scores from problem statement
    test_sequences = [
        ("TGCGTGCGCGCGCGCG", 87),
        ("GCGCCCGCGCGCGCGC", 82),
        ("GCGCGCGCGCGT", 71),
        ("CGCGCGCGCGC", 70),
        ("GCGCGCGCGCG", 70),
        ("GCGCGTGCGCGC", 65),
        ("GCGCGCCCGTACGCGC", 64),
        ("GCACGCACACGCGCGT", 64),
        ("CGCACGCGCACGCA", 62),
        ("CGCGCGCGCACA", 59),
        ("TGTGCGCGCGCACATG", 58),
    ]
    
    detector = ZDNADetector(use_kadane=True)
    
    print("\nTesting sequences:")
    print("Sequence\t\t\tExpected\tActual\t\tMethod\t\tSubclass")
    print("-" * 90)
    
    for seq, expected_score in test_sequences:
        candidates = detector.detect(seq, "test", "test", 0)
        
        if candidates:
            # Take the highest scoring candidate
            best_candidate = max(candidates, key=lambda c: c.raw_score if c.raw_score else 0)
            actual_score = best_candidate.raw_score if best_candidate.raw_score else 0
            method = best_candidate.scoring_method
            subclass = best_candidate.subclass
            
            print(f"{seq:<20}\t{expected_score:<8}\t{actual_score:<8.1f}\t{method:<20}\t{subclass}")
        else:
            print(f"{seq:<20}\t{expected_score:<8}\t{'No detection':<8}\t{'N/A':<20}\t{'N/A'}")
    
    print("\n" + "="*90)


def test_cgg_classification():
    """Test that CGG-rich sequences are classified as eGZ."""
    print("Testing CGG sequence classification...")
    
    detector = ZDNADetector(use_kadane=True)
    
    # Test sequences with different CGG content
    test_cases = [
        ("CGGCGGCGGCGGCGGCGG", "Should be eGZ (pure CGG repeats)"),
        ("CGCGCGCGCGCGCGCGCG", "Should be Z-DNA (CG repeats, not CGG)"),
        ("CGGCGCGGCGCGGCGCGG", "Should be eGZ (mixed but CGG-rich)"),
        ("ATATATATATATATATA", "Should be Z-DNA (AT repeats)"),
    ]
    
    for seq, description in test_cases:
        candidates = detector.detect(seq, "test", "test", 0)
        
        print(f"\nSequence: {seq}")
        print(f"Description: {description}")
        
        if candidates:
            for i, candidate in enumerate(candidates):
                print(f"  Candidate {i+1}:")
                print(f"    Subclass: {candidate.subclass}")
                print(f"    Score: {candidate.raw_score:.2f}")
                print(f"    Method: {candidate.scoring_method}")
        else:
            print("  No candidates detected")
    
    print("\n" + "="*90)


def test_kadane_scoring():
    """Test that eGZ motifs use Kadane algorithm scoring."""
    print("Testing Kadane algorithm scoring for eGZ...")
    
    detector = ZDNADetector(use_kadane=True)
    seq = "CGGCGGCGGCGGCGGCGG"  # Pure CGG repeats
    
    candidates = detector.detect(seq, "test", "test", 0)
    
    print(f"Test sequence: {seq}")
    print(f"Candidates found: {len(candidates)}")
    
    for i, candidate in enumerate(candidates):
        print(f"\nCandidate {i+1}:")
        print(f"  Position: {candidate.start}-{candidate.end}")
        print(f"  Raw score: {candidate.raw_score}")
        print(f"  Scoring method: {candidate.scoring_method}")
        print(f"  Subclass: {candidate.subclass}")
        print(f"  Pattern: {candidate.pattern_name}")
        
        # Verify it's using Kadane scoring
        if "Kadane" in candidate.scoring_method:
            print("  ✓ Using Kadane algorithm")
        else:
            print("  ✗ Not using Kadane algorithm")
    
    print("\n" + "="*90)


def test_fallback_behavior():
    """Test fallback behavior when Kadane is disabled."""
    print("Testing fallback behavior for eGZ...")
    
    detector = ZDNADetector(use_kadane=False)
    seq = "CGGCGGCGGCGGCGGCGG"
    
    candidates = detector.detect(seq, "test", "test", 0)
    
    print(f"Fallback detection for: {seq}")
    print(f"Candidates found: {len(candidates)}")
    
    for i, candidate in enumerate(candidates):
        print(f"\nCandidate {i+1}:")
        print(f"  Subclass: {candidate.subclass}")
        print(f"  Pattern: {candidate.pattern_name}")
        
        # Score the candidate using fallback scoring
        scored = detector.score([candidate])
        if scored:
            print(f"  Score: {scored[0].raw_score:.2f}")
            print(f"  Scoring method: {scored[0].scoring_method}")
    
    print("\n" + "="*90)


def main():
    """Run all eGZ motif tests."""
    print("=" * 90)
    print("eGZ (EXTRUDED-G) DNA MOTIF TESTING")
    print("=" * 90)
    
    test_egz_detection()
    test_cgg_classification()
    test_kadane_scoring()
    test_fallback_behavior()
    
    print("All eGZ motif tests completed!")
    print("=" * 90)


if __name__ == "__main__":
    main()