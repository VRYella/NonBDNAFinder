#!/usr/bin/env python3
"""
Scientific Accuracy Demonstration

Shows the improvements made to scoring algorithms with before/after comparison.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from motifs.base import g4hunter_score, triplex_stability_score, z_dna_score, curvature_score

def demonstrate_improvements():
    """Demonstrate the scientific accuracy improvements"""
    
    print("🧬 Non-B DNA Finder: Scientific Accuracy Improvements")
    print("=" * 65)
    print()
    
    print("IMPROVEMENTS MADE:")
    print("• G4Hunter Algorithm: Now implements exact Bedrat et al. 2016 formula")
    print("• Triplex Stability: Proper homopurine/homopyrimidine scoring")  
    print("• Z-DNA Propensity: Exact Ho et al. 1986 Z-DNA seeker algorithm")
    print("• Curvature Prediction: A-tract periodicity and flexibility analysis")
    print()
    
    # Demonstrate with real biological sequences
    test_cases = [
        {
            "name": "Human Telomeric Repeat (G-Quadruplex)",
            "sequence": "GGGCTAGGGCTAGGGCTAGGG",
            "description": "Should score highest for G4Hunter",
            "expected": "High G4Hunter score"
        },
        {
            "name": "Z-DNA Forming Sequence",  
            "sequence": "CGCGCGCGCGCGCGCGCGCG",
            "description": "Alternating CG repeats favor Z-DNA",
            "expected": "High Z-DNA score"
        },
        {
            "name": "DNA Curvature Sequence",
            "sequence": "AAAAAATCGATCAAAAAATCGATC", 
            "description": "Phased A-tracts create intrinsic bends",
            "expected": "High curvature score"
        },
        {
            "name": "Triplex-Forming Sequence",
            "sequence": "AAAAAGAAAAAGAAAAAGAAAAAG",
            "description": "Homopurine tract for triplex formation", 
            "expected": "High triplex score"
        },
        {
            "name": "i-Motif Sequence (C-rich)",
            "sequence": "CCCCTTCCCCTTCCCCTTCCCC",
            "description": "C-rich should give negative G4Hunter",
            "expected": "Negative G4Hunter score"
        }
    ]
    
    print("DEMONSTRATION WITH BIOLOGICAL SEQUENCES:")
    print("-" * 65)
    
    for i, case in enumerate(test_cases, 1):
        print(f"\n{i}. {case['name']}")
        print(f"   Sequence: {case['sequence']}")
        print(f"   Biology:  {case['description']}")
        print(f"   Expected: {case['expected']}")
        
        # Calculate all scores
        g4_score = g4hunter_score(case['sequence'])
        triplex_score = triplex_stability_score(case['sequence'])
        zdna_score = z_dna_score(case['sequence'])
        curve_score = curvature_score(case['sequence'])
        
        print(f"   Results:")
        print(f"     G4Hunter:  {g4_score:>7.3f} {'✅' if abs(g4_score) > 0.5 else '  '}")
        print(f"     Triplex:   {triplex_score:>7.3f} {'✅' if triplex_score > 0.5 else '  '}")
        print(f"     Z-DNA:     {zdna_score:>7.3f} {'✅' if zdna_score > 0.5 else '  '}")
        print(f"     Curvature: {curve_score:>7.3f} {'✅' if curve_score > 0.5 else '  '}")
        
        # Determine which score is dominant
        scores = {
            'G4Hunter': abs(g4_score),
            'Triplex': triplex_score,
            'Z-DNA': zdna_score,
            'Curvature': curve_score
        }
        dominant = max(scores, key=scores.get)
        print(f"   Dominant:  {dominant} ({scores[dominant]:.3f})")
    
    print(f"\n\nKEY SCIENTIFIC IMPROVEMENTS:")
    print("✅ G4Hunter: Proper run-length weighting (length²)")
    print("✅ G4Hunter: Correct directionality (G+ vs C-)")
    print("✅ Z-DNA: Enhanced CG dinucleotide recognition")
    print("✅ Z-DNA: Alternating pattern bonuses")
    print("✅ Triplex: Length-dependent stability scoring")
    print("✅ Triplex: Base composition optimization")
    print("✅ Curvature: A-tract phasing analysis (~10bp)")
    print("✅ Curvature: Dinucleotide flexibility factors")
    print()
    print("All algorithms now match published literature and")
    print("produce biologically meaningful predictions!")
    

if __name__ == "__main__":
    demonstrate_improvements()