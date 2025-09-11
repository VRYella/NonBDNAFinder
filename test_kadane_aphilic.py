#!/usr/bin/env python3
"""
Test the enhanced A-philic detector with Kadane's algorithm
"""

import sys
import numpy as np
from motif_detectors import APhilicDetector

def test_enhanced_aphilic_detector():
    """Test the enhanced A-philic detector"""
    print("=== Testing Enhanced A-philic Detector with Kadane's Algorithm ===")
    
    # Create detector instance
    detector = APhilicDetector()
    
    # Test sequence with A-philic patterns
    test_sequence = "ATGCCCCCGGGGAAAATTTTCCCCGGGGAAAAATGCCCAAATTTGGGTACCCCGGGGTTTAAACCCGGGAATTCCCGGGAAATTTCCCGGGTTTAAACCCGGGTTTT"
    
    print(f"Test sequence length: {len(test_sequence)}")
    print(f"Test sequence: {test_sequence}")
    
    # Run detection
    candidates = detector.detect(test_sequence, "test_seq", "test_contig", 0)
    
    print(f"\nFound {len(candidates)} A-philic regions:")
    
    for i, candidate in enumerate(candidates):
        print(f"\nRegion {i+1}:")
        print(f"  Position: {candidate.start}-{candidate.end}")
        print(f"  Length: {candidate.length}")
        print(f"  Raw score: {candidate.raw_score:.6f}")
        print(f"  Sequence: {candidate.matched_seq.decode('utf-8')}")
        if hasattr(candidate, 'metadata') and candidate.metadata:
            print(f"  Metadata: {candidate.metadata}")
    
    # Score the candidates
    scored_candidates = detector.score(candidates)
    
    print(f"\nAfter scoring:")
    for i, candidate in enumerate(scored_candidates):
        print(f"Region {i+1}: Score = {candidate.raw_score:.6f}, Method = {candidate.scoring_method}")
    
    return len(candidates) > 0

if __name__ == "__main__":
    success = test_enhanced_aphilic_detector()
    print(f"\nTest {'PASSED' if success else 'FAILED'}")
    sys.exit(0 if success else 1)