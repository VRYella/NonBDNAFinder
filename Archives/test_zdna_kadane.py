#!/usr/bin/env python3
"""
Test script for the new Kadane algorithm-based Z-DNA detection.

This script validates that the Z-DNA calculator and Kadane algorithm
implementation are working correctly for various test cases.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from constants import Params
from zdna_calculator import ZDNACalculatorSeq
import numpy as np


def test_basic_functionality():
    """Test basic functionality of the Z-DNA calculator."""
    print("Testing basic Z-DNA calculator functionality...")
    
    # Test simple CG repeats
    seq = "CGCGCGCGCGCG"
    params = Params()
    calculator = ZDNACalculatorSeq(seq, params)
    
    # Test transitions scoring
    scoring = calculator.zdna_calculator_transitions()
    print(f"CG sequence scoring: {scoring}")
    print(f"Total score: {np.sum(scoring)}")
    
    # Test subarray detection
    subarrays = calculator.subarrays_above_threshold()
    print(f"Subarrays found: {len(subarrays)}")
    for subarray in subarrays:
        print(f"  Start: {subarray[0]}, End: {subarray[1]}, Score: {subarray[2]:.2f}, Sequence: {subarray[3]}")
    
    print()


def test_different_transitions():
    """Test different dinucleotide transitions."""
    print("Testing different dinucleotide transitions...")
    
    params = Params(threshold=5.0)
    
    test_sequences = [
        ("CGCGCGCGCG", "CG repeats"),
        ("ATATATATATAT", "AT repeats"),
        ("GTGTGTGTGT", "GT repeats"),
        ("ACACACACAC", "AC repeats"),
        ("ATCGATCGATCG", "Mixed transitions"),
        ("AAAAAAAAAA", "No transitions (A only)")
    ]
    
    for seq, desc in test_sequences:
        calculator = ZDNACalculatorSeq(seq, params)
        scoring = calculator.zdna_calculator_transitions()
        subarrays = calculator.subarrays_above_threshold()
        total_score = np.sum(scoring)
        
        print(f"{desc}: {seq}")
        print(f"  Total score: {total_score:.2f}")
        print(f"  Subarrays: {len(subarrays)}")
        if subarrays:
            for subarray in subarrays:
                print(f"    [{subarray[0]}:{subarray[1]}] Score: {subarray[2]:.2f}")
        print()


def test_kadane_algorithm():
    """Test the Kadane algorithm behavior with specific scoring patterns."""
    print("Testing Kadane algorithm behavior...")
    
    # Create a sequence with alternating high and low scoring regions
    seq = "AAACGCGCGCGAAACGCGCGCGAAA"  # A's should get penalties, CG's should get high scores
    params = Params(threshold=8.0, drop_threshold=20.0)
    
    calculator = ZDNACalculatorSeq(seq, params)
    scoring = calculator.zdna_calculator_transitions()
    subarrays = calculator.subarrays_above_threshold()
    
    print(f"Test sequence: {seq}")
    print(f"Scoring array: {[f'{s:.1f}' for s in scoring]}")
    print(f"Subarrays found: {len(subarrays)}")
    
    for i, subarray in enumerate(subarrays):
        print(f"  Subarray {i+1}: [{subarray[0]}:{subarray[1]}] Score: {subarray[2]:.2f} Seq: '{subarray[3]}'")
    
    print()


def test_consecutive_at_scoring():
    """Test consecutive AT scoring functionality."""
    print("Testing consecutive AT scoring...")
    
    params = Params(
        threshold=3.0,
        consecutive_AT_scoring=(2.0, 1.0, 0.5),  # Decreasing bonus for consecutive ATs
        AT_weight=1.0
    )
    
    # Test with multiple consecutive ATs
    seq = "ATATATATATAT"  # Should get decreasing bonuses for consecutive ATs
    calculator = ZDNACalculatorSeq(seq, params)
    scoring = calculator.zdna_calculator_transitions()
    
    print(f"AT repeat sequence: {seq}")
    print(f"Scoring with consecutive AT bonuses: {[f'{s:.1f}' for s in scoring]}")
    print(f"Expected pattern: decreasing scores due to consecutive AT penalty")
    
    subarrays = calculator.subarrays_above_threshold()
    print(f"Subarrays: {len(subarrays)}")
    for subarray in subarrays:
        print(f"  Score: {subarray[2]:.2f}, Sequence: {subarray[3]}")
    
    print()


def test_parameter_variations():
    """Test different parameter configurations."""
    print("Testing parameter variations...")
    
    seq = "CGCGCGCGCGCG"
    base_params = Params()
    
    # Test different thresholds
    thresholds = [5.0, 10.0, 20.0]
    for threshold in thresholds:
        params = Params(threshold=threshold)
        calculator = ZDNACalculatorSeq(seq, params)
        subarrays = calculator.subarrays_above_threshold()
        print(f"Threshold {threshold}: {len(subarrays)} subarrays")
    
    # Test different weights
    weights = [1.0, 3.0, 5.0]
    for gc_weight in weights:
        params = Params(GC_weight=gc_weight, threshold=5.0)
        calculator = ZDNACalculatorSeq(seq, params)
        scoring = calculator.zdna_calculator_transitions()
        total_score = np.sum(scoring)
        print(f"GC weight {gc_weight}: Total score {total_score:.2f}")
    
    print()


def main():
    """Run all Z-DNA tests."""
    print("="*60)
    print("Z-DNA KADANE ALGORITHM TESTING")
    print("="*60)
    
    test_basic_functionality()
    test_different_transitions()
    test_kadane_algorithm()
    test_consecutive_at_scoring()
    test_parameter_variations()
    
    print("="*60)
    print("All tests completed!")
    print("="*60)


if __name__ == "__main__":
    main()