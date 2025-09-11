#!/usr/bin/env python3
"""
Comprehensive comparison of A-philic detection methods:
1. Original stringent method (check_a_philic_stringent.py)
2. Enhanced method with Kadane's algorithm (enhanced_aphilic_hs.py)
3. Enhanced detector class (motif_detectors.py)
"""

import subprocess
import sys
from pathlib import Path
import csv
import tempfile
from motif_detectors import APhilicDetector

def run_original_method(sequence, tetra_file, tri_file):
    """Run the original stringent A-philic detection"""
    print("=== Original Stringent Method ===")
    
    cmd = [
        "python3", "check_a_philic_stringent.py",
        "--tetra", str(tetra_file),
        "--tri", str(tri_file),
        "--seq", sequence,
        "--out", "/tmp/original_results.tsv"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    print("Output:", result.stdout)
    
    return result.stdout

def run_enhanced_standalone(sequence, tetra_file, tri_file):
    """Run the enhanced standalone script"""
    print("\n=== Enhanced Standalone Method (Kadane's Algorithm) ===")
    
    # Create temporary FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(">test_sequence\n")
        f.write(sequence + "\n")
        fasta_file = f.name
    
    cmd = [
        "python3", "enhanced_aphilic_hs.py",
        "--fasta", fasta_file,
        "--tetra_table", str(tetra_file),
        "--tri_table", str(tri_file),
        "--outdir", "/tmp/enhanced_results",
        "--threshold", "10.0"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    print("Output:", result.stdout)
    
    # Read results
    results_file = Path("/tmp/enhanced_results") / "test_sequence_enhanced_aphilic_regions.csv"
    if results_file.exists():
        print("\nEnhanced Results:")
        with open(results_file) as f:
            reader = csv.DictReader(f)
            for i, row in enumerate(reader):
                print(f"  Region {i+1}: {row['Start']}-{row['End']} (len={row['length']}, score={row['combined_score']})")
    
    # Cleanup
    Path(fasta_file).unlink(missing_ok=True)
    
    return result.stdout

def run_enhanced_detector_class(sequence):
    """Run the enhanced detector class"""
    print("\n=== Enhanced Detector Class (Integrated) ===")
    
    detector = APhilicDetector()
    candidates = detector.detect(sequence, "test_seq", "test_contig", 0)
    
    print(f"Found {len(candidates)} regions:")
    for i, candidate in enumerate(candidates):
        print(f"  Region {i+1}: {candidate.start}-{candidate.end} (len={candidate.length}, score={candidate.raw_score:.3f})")
        if hasattr(candidate, 'metadata') and candidate.metadata:
            metadata = candidate.metadata
            print(f"    Tetra sum: {metadata.get('tetra_sum', 0):.3f}, Tri sum: {metadata.get('tri_sum', 0):.3f}")
    
    return candidates

def main():
    """Compare all three methods"""
    print("=== Comprehensive A-philic Detection Comparison ===")
    
    # Test sequence with known A-philic patterns
    test_sequence = "ATGCCCCCGGGGAAAATTTTCCCCGGGGAAAAATGCCCAAATTTGGGTACCCCGGGGTTTAAACCCGGGAATTCCCGGGAAATTTCCCGGGTTTAAACCCGGGTTTT"
    
    print(f"Test sequence length: {len(test_sequence)}")
    print(f"Sequence: {test_sequence[:50]}...{test_sequence[-20:]}")
    
    # Use test propensity tables
    tetra_file = Path("test_tetra_propensities.csv")
    tri_file = Path("test_tri_propensities.csv")
    
    if not tetra_file.exists() or not tri_file.exists():
        print("Error: Test propensity files not found")
        return False
    
    # Run all three methods
    try:
        original_output = run_original_method(test_sequence, tetra_file, tri_file)
        enhanced_output = run_enhanced_standalone(test_sequence, tetra_file, tri_file)
        detector_candidates = run_enhanced_detector_class(test_sequence)
        
        # Summary comparison
        print("\n" + "="*60)
        print("COMPARISON SUMMARY")
        print("="*60)
        
        print("\n1. Original Method:")
        if "chrom\tstart\tend" in original_output:
            lines = [l for l in original_output.split('\n') if l and not l.startswith('chrom')]
            print(f"   Found {len(lines)} regions using fixed 10-mer windows")
        else:
            print("   No regions found or different output format")
        
        print("\n2. Enhanced Standalone:")
        print("   Uses Kadane's algorithm for optimal subarray detection")
        print("   Integrates Hyperscan for k-mer matching")
        
        print("\n3. Enhanced Detector Class:")
        print(f"   Found {len(detector_candidates)} regions using Kadane's algorithm")
        print("   Integrated with existing motif detection framework")
        
        print("\nKey Improvements:")
        print("- Dynamic region boundaries instead of fixed 10-mer windows")
        print("- Optimal subarray detection maximizing combined scores")
        print("- Enhanced nucleation filtering")
        print("- Better integration with existing framework")
        
        return True
        
    except Exception as e:
        print(f"Error during comparison: {e}")
        return False

if __name__ == "__main__":
    success = main()
    print(f"\nComparison {'COMPLETED SUCCESSFULLY' if success else 'FAILED'}")
    sys.exit(0 if success else 1)