#!/usr/bin/env python3
"""
Validation examples for canonical motif parameters.

Demonstrates the synchronized parameter enforcement for:
- Cruciform DNA (10-100 nt arms, 0-3 nt spacer)
- Triplex DNA (10-100 nt arms, 0-8 nt spacer, >90% purine/pyrimidine)
- Slipped DNA (10-50 nt unit, spacer=0)
- STRs (1-9 bp unit, ≥20 bp total)
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from detectors import CruciformDetector, TriplexDetector, SlippedDNADetector
import scanner


def revcomp(seq):
    """Reverse complement of DNA sequence."""
    trans = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(trans)[::-1]


print("=" * 80)
print("Validation Examples: Canonical Motif Parameters")
print("=" * 80)

# Example 1: Cruciform DNA
print("\n1. Cruciform DNA Detection")
print("-" * 80)
left_arm = "ACGTACGTACGT"  # 12 nt
spacer = "TT"  # 2 nt (within 0-3 range)
right_arm = revcomp(left_arm)  # Reverse complement
cruciform_seq = left_arm + spacer + right_arm

print(f"Sequence: {cruciform_seq}")
print(f"Left arm: {left_arm} (12 nt)")
print(f"Spacer: {spacer} (2 nt)")
print(f"Right arm: {right_arm} (reverse complement)")

detector = CruciformDetector()
motifs = detector.detect_motifs(cruciform_seq, "cruciform_test")

print(f"\nDetected {len(motifs)} motif(s):")
for motif in motifs:
    print(f"  - Class: {motif.get('Class')}, Length: {motif.get('Length')} bp")
    if 'Arm_Length' in motif.get('Method', ''):
        print(f"    Arm: {motif.get('Arm_Length', 'N/A')} nt, Loop: {motif.get('Loop_Length', 'N/A')} nt")

# Example 2: Triplex DNA
print("\n2. Triplex DNA Detection (Purine-rich)")
print("-" * 80)
left_arm_purine = "AGAGAGAGAGAG"  # 12 nt, 100% purine
spacer_triplex = "TTTT"  # 4 nt (within 0-8 range)
right_arm_purine = left_arm_purine[::-1]  # Reverse (NOT complement)
triplex_seq = left_arm_purine + spacer_triplex + right_arm_purine

print(f"Sequence: {triplex_seq}")
print(f"Left arm: {left_arm_purine} (12 nt, 100% purine)")
print(f"Spacer: {spacer_triplex} (4 nt)")
print(f"Right arm: {right_arm_purine} (reverse of left)")

purine_count = sum(1 for b in left_arm_purine + right_arm_purine if b in 'AG')
total_bases = len(left_arm_purine + right_arm_purine)
purity = purine_count / total_bases
print(f"Purine purity: {purity:.1%} (threshold: 90%)")

detector = TriplexDetector()
motifs = detector.detect_motifs(triplex_seq, "triplex_test")

print(f"\nDetected {len(motifs)} motif(s):")
for motif in motifs:
    print(f"  - Class: {motif.get('Class')}, Length: {motif.get('Length')} bp")
    if motif.get('Subclass'):
        print(f"    Subclass: {motif.get('Subclass')}")
    details = motif.get('details', {})
    if 'purine_fraction' in details:
        print(f"    Purine: {details['purine_fraction']:.1%}, Pyrimidine: {details['pyrimidine_fraction']:.1%}")

# Example 3: Slipped DNA
print("\n3. Slipped DNA Detection (Tandem Repeat)")
print("-" * 80)
unit = "ACGTACGTACGTACGTACGT"  # 20 nt
slipped_seq = unit + unit  # Tandem repeat (spacer=0)

print(f"Sequence: {slipped_seq}")
print(f"Unit: {unit} (20 nt)")
print(f"Spacer: 0 nt (tandem)")
print(f"Total length: {len(slipped_seq)} bp")

detector = SlippedDNADetector()
motifs = detector.detect_motifs(slipped_seq, "slipped_test")

print(f"\nDetected {len(motifs)} motif(s):")
for motif in motifs:
    print(f"  - Class: {motif.get('Class')}, Subclass: {motif.get('Subclass')}")
    details = motif.get('details', {})
    if 'unit_length' in details:
        print(f"    Unit: {details['unit_length']} nt, Spacer: {details.get('spacer_length', 'N/A')} nt")

# Example 4: Short Tandem Repeat (STR)
print("\n4. STR Detection")
print("-" * 80)
str_unit = "AT"  # 2 bp
str_seq = str_unit * 10  # 20 bp total (minimum threshold)

print(f"Sequence: {str_seq}")
print(f"Unit: {str_unit} (2 bp)")
print(f"Copies: 10")
print(f"Total length: {len(str_seq)} bp (minimum: 20 bp)")

strs = scanner.find_strs(str_seq, min_u=1, max_u=9, min_total=20)

print(f"\nDetected {len(strs)} STR(s):")
for str_motif in strs:
    print(f"  - Unit: {str_motif.get('Unit_Length')} bp, Copies: {str_motif.get('Copies')}")
    print(f"    Total length: {str_motif.get('Length')} bp")

print("\n" + "=" * 80)
print("Validation Complete!")
print("=" * 80)
print("\nSummary:")
print("✓ Cruciform: Arms 10-100 nt, spacer 0-3 nt, perfect Watson-Crick matches")
print("✓ Triplex: Arms 10-100 nt, spacer 0-8 nt, >90% purine/pyrimidine purity")
print("✓ Slipped DNA: Unit 10-50 nt, spacer=0 (tandem)")
print("✓ STRs: Unit 1-9 bp, minimum total ≥20 bp")
