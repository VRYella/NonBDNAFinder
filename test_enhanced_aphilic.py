#!/usr/bin/env python3
"""
Test script for enhanced A-philic detector implementation
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from motif_detectors import APhilicDetector

def test_enhanced_aphilic_detector():
    """Test the enhanced A-philic detector implementation"""
    print("Testing Enhanced A-philic Detector...")
    
    # Create detector instance
    detector = APhilicDetector()
    
    # Test sequence from the call_Aphilic.py example
    test_seq = ("ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC"
                "ATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT"
                "GGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGTAGGAACCCCTAACCCCTAA"
                "GAAAGAAGAAGAAGAAGAAGAAAGGAAGGAAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGG"
                "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"
                "GAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA"
                "CTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT")
    
    print(f"Test sequence length: {len(test_seq)} bp")
    
    # Run detection
    try:
        candidates = detector.detect(
            seq=test_seq,
            seq_name="test_sequence", 
            contig="test_contig",
            offset=0
        )
        
        print(f"Found {len(candidates)} A-philic regions")
        
        if candidates:
            # Score the candidates
            scored_candidates = detector.score(candidates)
            
            print("\nDetected A-philic regions:")
            print("Start\tEnd\tLength\tScore\tSequence")
            for i, candidate in enumerate(scored_candidates[:5]):  # Show first 5
                seq_preview = candidate.matched_seq.decode('utf-8')[:50] + "..." if len(candidate.matched_seq) > 50 else candidate.matched_seq.decode('utf-8')
                print(f"{candidate.start}\t{candidate.end}\t{candidate.length}\t{candidate.raw_score:.3f}\t{seq_preview}")
                
                # Print metadata if available
                if hasattr(candidate, 'metadata') and candidate.metadata:
                    print(f"  Metadata: {candidate.metadata}")
        else:
            print("No A-philic regions detected")
            
        print(f"\nTest completed successfully!")
        return True
        
    except Exception as e:
        print(f"Error during A-philic detection: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_propensity_tables():
    """Test that propensity tables are loaded correctly"""
    print("\nTesting propensity tables...")
    
    detector = APhilicDetector()
    
    print(f"Loaded {len(detector.TETRA_LOG2)} tetranucleotide propensities")
    print(f"Loaded {len(detector.TRI_LOG2)} trinucleotide propensities")
    
    # Check some known high-scoring patterns
    high_tetra = ["CCCC", "GGGG", "TGGG", "GGGC"]
    for tetra in high_tetra:
        score = detector.TETRA_LOG2.get(tetra, "NOT_FOUND")
        print(f"  {tetra}: {score}")
    
    high_tri = ["CCC", "GGG", "CAC", "GCC"]  
    for tri in high_tri:
        score = detector.TRI_LOG2.get(tri, "NOT_FOUND")
        print(f"  {tri}: {score}")
        
    print("Propensity tables loaded correctly!")

if __name__ == "__main__":
    print("Enhanced A-philic Detector Test")
    print("=" * 40)
    
    # Test propensity tables
    test_propensity_tables()
    
    # Test detector functionality
    success = test_enhanced_aphilic_detector()
    
    if success:
        print("\n✅ All tests passed!")
    else:
        print("\n❌ Tests failed!")
        sys.exit(1)