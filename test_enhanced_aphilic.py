#!/usr/bin/env python3
"""
Test script to validate the enhanced A-philic detection with Kadane's algorithm
"""

import subprocess
import sys
from pathlib import Path
import csv

def test_enhanced_aphilic():
    """Test the enhanced A-philic implementation"""
    print("=== Testing Enhanced A-philic Detection with Kadane's Algorithm ===")
    
    # Test parameters
    test_sequence = "ATGCCCCCGGGGAAAATTTTCCCCGGGGAAAAATGCCCAAATTTGGGTACCCCGGGGTTTAAACCCGGGAATTCCCGGGAAATTTCCCGGGTTTAAACCCGGGTTTT"
    
    # Create test FASTA
    with open("test_enhanced.fasta", "w") as f:
        f.write(">test_enhanced_seq\n")
        f.write(test_sequence + "\n")
    
    # Run enhanced implementation
    cmd = [
        "python3", "enhanced_aphilic_hs.py",
        "--fasta", "test_enhanced.fasta",
        "--tetra_table", "test_tetra_propensities.csv",
        "--tri_table", "test_tri_propensities.csv",
        "--outdir", "/tmp/test_enhanced_validation",
        "--threshold", "20.0",  # Test with different threshold
        "--tri_nucleation_min", "2"  # Test with stricter nucleation
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    print("Enhanced implementation output:")
    print(result.stdout)
    if result.stderr:
        print("Errors:", result.stderr)
    
    # Read and display results
    results_file = Path("/tmp/test_enhanced_validation/test_enhanced_enhanced_aphilic_regions.csv")
    if results_file.exists():
        print("\n=== Enhanced Results ===")
        with open(results_file) as f:
            reader = csv.DictReader(f)
            for i, row in enumerate(reader):
                print(f"Region {i+1}:")
                for key, value in row.items():
                    print(f"  {key}: {value}")
                print()
    else:
        print("No results file found")
    
    # Compare with original implementation
    print("\n=== Comparing with Original Implementation ===")
    cmd_orig = [
        "python3", "check_a_philic_stringent.py",
        "--tetra", "test_tetra_propensities.csv",
        "--tri", "test_tri_propensities.csv",
        "--seq", test_sequence,
        "--out", "/tmp/original_comparison.tsv"
    ]
    
    result_orig = subprocess.run(cmd_orig, capture_output=True, text=True)
    print("Original implementation output:")
    print(result_orig.stdout)
    
    # Cleanup
    Path("test_enhanced.fasta").unlink(missing_ok=True)
    
    return result.returncode == 0

if __name__ == "__main__":
    success = test_enhanced_aphilic()
    print(f"\nTest {'PASSED' if success else 'FAILED'}")
    sys.exit(0 if success else 1)