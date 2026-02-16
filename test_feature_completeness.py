#!/usr/bin/env python3
"""
Test script to verify all detectors output comprehensive features.
Verifies that all required fields are present in detector output:
- Type_Of_Repeat
- Criterion
- Disease_Relevance
- Regions_Involved
- GC_Content
- Arm_Length
- Loop_Length
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import all detectors
from Detectors.gquad.detector import GQuadruplexDetector
from Detectors.zdna.detector import ZDNADetector
from Detectors.cruciform.detector import CruciformDetector
from Detectors.triplex.detector import TriplexDetector
from Detectors.imotif.detector import IMotifDetector
from Detectors.rloop.detector import RLoopDetector
from Detectors.slipped.detector import SlippedDNADetector
from Detectors.curved.detector import CurvedDNADetector
from Detectors.aphilic.detector import APhilicDetector

# Test sequences
TEST_SEQUENCES = {
    'G-Quadruplex': ('GGGAGGGAGGGAGGG', GQuadruplexDetector),
    'Slipped_DNA': ('CAGCAGCAGCAGCAGCAGCAGCAG', SlippedDNADetector),
    'i-Motif': ('CCCCACCCCACCCCACCCC', IMotifDetector),
    'Z-DNA': ('CGCGCGCGCGCGCGCG', ZDNADetector),
    'Curved': ('AAAAAAATTTTTTTAAAAAA', CurvedDNADetector),
    'Cruciform': ('ATCGATCGATCGGGGCGATCGATCGAT', CruciformDetector),
    'Triplex': ('GAAGAAGAAGAAGAAGAA', TriplexDetector),
    'A-philic': ('AACGAACGAACGAACGAACG', APhilicDetector),
}

# Required universal fields
REQUIRED_FIELDS = [
    'Type_Of_Repeat',
    'Criterion',
    'Disease_Relevance',
    'Regions_Involved',
    'GC_Content',
    'Arm_Length',
    'Loop_Length'
]

# Core fields (should always be present)
CORE_FIELDS = [
    'ID', 'Sequence_Name', 'Class', 'Subclass',
    'Start', 'End', 'Length', 'Sequence',
    'Score', 'Strand', 'Method', 'Pattern_ID'
]

print("=" * 80)
print("Testing Feature Completeness for All Detectors")
print("=" * 80)
print()

all_passed = True
results = {}

for detector_name, (sequence, DetectorClass) in TEST_SEQUENCES.items():
    print(f"\n{'='*60}")
    print(f"Testing {detector_name}")
    print(f"{'='*60}")
    print(f"Sequence: {sequence}")
    print(f"Length: {len(sequence)} bp\n")
    
    try:
        detector = DetectorClass()
        motifs = detector.detect_motifs(sequence, f'{detector_name}_test')
        
        if not motifs:
            print(f"⚠️  No motifs detected (this might be expected for some test sequences)")
            results[detector_name] = {'status': 'no_motifs', 'missing': []}
            continue
        
        print(f"✓ Detected {len(motifs)} motif(s)")
        
        # Check first motif for completeness
        motif = motifs[0]
        print(f"\nMotif Details:")
        print(f"  Class: {motif.get('Class', 'N/A')}")
        print(f"  Subclass: {motif.get('Subclass', 'N/A')}")
        print(f"  Position: {motif.get('Start', 'N/A')}-{motif.get('End', 'N/A')}")
        
        # Check core fields
        missing_core = [f for f in CORE_FIELDS if f not in motif]
        if missing_core:
            print(f"\n❌ Missing core fields: {', '.join(missing_core)}")
            all_passed = False
        else:
            print(f"\n✓ All core fields present")
        
        # Check required universal fields
        missing_required = []
        present_required = []
        
        for field in REQUIRED_FIELDS:
            if field in motif:
                value = motif[field]
                present_required.append(field)
                # Show a preview of the value
                if isinstance(value, str) and len(str(value)) > 50:
                    print(f"  ✓ {field}: {str(value)[:50]}...")
                else:
                    print(f"  ✓ {field}: {value}")
            else:
                missing_required.append(field)
                print(f"  ❌ {field}: MISSING")
        
        results[detector_name] = {
            'status': 'tested',
            'missing': missing_required,
            'present': present_required,
            'total_fields': len(motif.keys())
        }
        
        if missing_required:
            print(f"\n❌ Missing required fields: {', '.join(missing_required)}")
            all_passed = False
        else:
            print(f"\n✓ All required universal fields present!")
        
        print(f"\nTotal fields in motif: {len(motif.keys())}")
        
    except Exception as e:
        print(f"❌ Error testing {detector_name}: {e}")
        import traceback
        traceback.print_exc()
        all_passed = False
        results[detector_name] = {'status': 'error', 'error': str(e)}

# Summary
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

for detector_name, result in results.items():
    status = result['status']
    if status == 'tested':
        if result['missing']:
            print(f"❌ {detector_name}: Missing {len(result['missing'])} fields - {', '.join(result['missing'])}")
        else:
            print(f"✓ {detector_name}: All required fields present ({result['total_fields']} total fields)")
    elif status == 'no_motifs':
        print(f"⚠️  {detector_name}: No motifs detected")
    elif status == 'error':
        print(f"❌ {detector_name}: Error - {result.get('error', 'Unknown')}")

print("\n" + "=" * 80)
if all_passed:
    print("✓ ALL TESTS PASSED!")
    print("All detectors output comprehensive features for meticulous analysis.")
else:
    print("❌ SOME TESTS FAILED")
    print("Please review the errors above.")
print("=" * 80)
print()

sys.exit(0 if all_passed else 1)
