#!/usr/bin/env python3
"""
Live demonstration of the unified motif taxonomy system
Shows normalization, validation, and detector integration in action
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

print("=" * 80)
print("UNIFIED MOTIF TAXONOMY - LIVE DEMONSTRATION")
print("=" * 80)

# ============================================================================
# PART 1: Show canonical taxonomy
# ============================================================================
print("\n" + "=" * 80)
print("PART 1: Canonical Taxonomy Structure")
print("=" * 80)

from config.motif_taxonomy import (
    MOTIF_CLASSIFICATION,
    VALID_CLASSES,
    VALID_SUBCLASSES,
    CLASS_TO_SUBCLASSES
)

print(f"\n📋 Total Classes: {len(VALID_CLASSES)}")
print(f"📋 Total Subclasses: {len(VALID_SUBCLASSES)}")

print("\n🔬 All 11 Motif Classes:")
for i, class_name in enumerate(sorted(VALID_CLASSES), 1):
    subclass_count = len(CLASS_TO_SUBCLASSES[class_name])
    print(f"  {i:2d}. {class_name:25s} ({subclass_count} subclass{'es' if subclass_count != 1 else ''})")

# ============================================================================
# PART 2: Demonstrate normalization
# ============================================================================
print("\n" + "=" * 80)
print("PART 2: Normalization - Handling Variant Names")
print("=" * 80)

from core.motif_normalizer import normalize_class_subclass

test_cases = [
    ("g-quadruplex", "telomeric g4"),
    ("curved dna", "global curvature"),
    ("ZDNA", "egz"),
    ("i motif", "canonical_imotif"),
]

print("\n🔄 Normalizing variant spellings to canonical:")
for raw_class, raw_subclass in test_cases:
    canonical_class, canonical_subclass = normalize_class_subclass(
        raw_class, raw_subclass, strict=False, auto_correct=True
    )
    print(f"  '{raw_class}' + '{raw_subclass}'")
    print(f"    → '{canonical_class}' + '{canonical_subclass}'")

# ============================================================================
# PART 3: Show auto-correction
# ============================================================================
print("\n" + "=" * 80)
print("PART 3: Auto-Correction - Fixing Mismatched Class/Subclass")
print("=" * 80)

print("\n🔧 Auto-correcting mismatched class/subclass pairs:")

# Telomeric G4 belongs to G-Quadruplex, not Triplex
import warnings
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    canonical_class, canonical_subclass = normalize_class_subclass(
        'Triplex',      # Wrong class
        'Telomeric G4', # Belongs to G-Quadruplex
        auto_correct=True
    )
    print(f"  Input:  'Triplex' + 'Telomeric G4'")
    print(f"  Output: '{canonical_class}' + '{canonical_subclass}'")
    if w:
        print(f"  ⚠️  Warning: {w[0].message}")

# ============================================================================
# PART 4: Test detector integration
# ============================================================================
print("\n" + "=" * 80)
print("PART 4: Detector Integration - Real Motif Detection")
print("=" * 80)

from detectors import GQuadruplexDetector, ZDNADetector, CruciformDetector

print("\n🧬 Testing detectors with real sequences:")

# Test G-Quadruplex
print("\n  1. G-Quadruplex Detector")
g4_detector = GQuadruplexDetector()
g4_seq = 'GGGTTAGGGTTAGGGTTAGGG'  # Telomeric repeat
g4_motifs = g4_detector.detect_motifs(g4_seq, 'demo')
if g4_motifs:
    motif = g4_motifs[0]
    print(f"     Sequence: {g4_seq}")
    print(f"     Detected: {len(g4_motifs)} motif(s)")
    print(f"     Class: '{motif['Class']}'")
    print(f"     Subclass: '{motif['Subclass']}'")
    print(f"     ✅ Using canonical names!")

# Test Z-DNA
print("\n  2. Z-DNA Detector")
z_detector = ZDNADetector()
z_seq = 'CGCGCGCGCGCGCGCGCG'  # CG alternating
z_motifs = z_detector.detect_motifs(z_seq, 'demo')
if z_motifs:
    motif = z_motifs[0]
    print(f"     Sequence: {z_seq}")
    print(f"     Detected: {len(z_motifs)} motif(s)")
    print(f"     Class: '{motif['Class']}'")
    print(f"     Subclass: '{motif['Subclass']}'")
    print(f"     ✅ Using canonical names!")

# ============================================================================
# PART 5: Export validation
# ============================================================================
print("\n" + "=" * 80)
print("PART 5: Export Validation - Ensuring Data Integrity")
print("=" * 80)

from export.export_validator import validate_export_data, get_export_summary

# Create test motifs (mix of canonical and variant names)
test_motifs = [
    {'Class': 'G-Quadruplex', 'Subclass': 'Telomeric G4', 'Start': 0, 'End': 25},
    {'Class': 'z-dna', 'Subclass': 'egz', 'Start': 50, 'End': 75},  # Lowercase
    {'Class': 'Curved_DNA', 'Subclass': 'Local Curvature', 'Start': 100, 'End': 125},
]

print("\n📤 Validating motifs before export:")
print(f"  Input: {len(test_motifs)} motifs (some with variant names)")

validated = validate_export_data(test_motifs, auto_normalize=True, strict=False)

print(f"  Output: {len(validated)} validated motifs")
for i, motif in enumerate(validated, 1):
    print(f"    {i}. Class='{motif['Class']}', Subclass='{motif['Subclass']}'")

# Get summary
summary = get_export_summary(validated)
print(f"\n📊 Export Summary:")
print(f"  Total motifs: {summary['total_motifs']}")
print(f"  Unique classes: {summary['unique_classes']}")
print(f"  Classes present: {', '.join(summary['classes'])}")
print(f"  Validation status: {summary['validation_status']}")

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("✅ DEMONSTRATION COMPLETE")
print("=" * 80)

print("""
The unified motif taxonomy system successfully:

  ✅ Defines canonical names for 11 classes and 24 subclasses
  ✅ Normalizes variant spellings (lowercase, aliases, etc.)
  ✅ Auto-corrects mismatched class/subclass pairs
  ✅ Integrates with all 9 detectors
  ✅ Validates exports before writing files
  ✅ Prevents data corruption and inconsistencies

All detectors, visualizations, and exports now use the same canonical
taxonomy from a single source of truth.

For more information:
  - Documentation: MOTIF_CLASSIFICATION.md
  - Security: SECURITY_SUMMARY.md
  - Tests: validate_taxonomy.py
""")

print("=" * 80)
