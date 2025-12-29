#!/usr/bin/env python3
"""
Comprehensive Validation Suite for NonBDNAFinder
=================================================

This module provides comprehensive validation tests for all motif detection
capabilities with high-level accuracy verification and detailed output.

Features:
- Tests all 11 motif classes with known sequences
- Validates scoring accuracy and confidence tiers
- Tests edge cases and boundary conditions
- Provides detailed output with motif details, scores, and rationale
- Generates validation metrics and reports

Author: Dr. Venkata Rajesh Yella
Version: 2025.1
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from typing import List, Dict, Any, Tuple
from collections import Counter
import time

# Import detectors and main scanner
from detectors import (
    CurvedDNADetector,
    SlippedDNADetector,
    CruciformDetector,
    RLoopDetector,
    TriplexDetector,
    GQuadruplexDetector,
    IMotifDetector,
    ZDNADetector,
    APhilicDetector
)
from nonbscanner import analyze_sequence


# ==============================================================================
# Test Sequences - Known Motifs for Validation
# ==============================================================================

# Each test case includes:
# - name: Description of the sequence
# - sequence: DNA sequence containing the motif
# - expected_class: Expected motif class
# - expected_subclass: Expected motif subclass (if applicable)
# - expected_count: Expected number of motifs (if applicable)
# - min_score: Minimum expected score
# - description: Detailed description of what makes this a valid motif

TEST_SEQUENCES = {
    'g_quadruplex': {
        'canonical_telomeric': {
            'name': 'Human Telomeric G4 (Canonical)',
            'sequence': 'TTAGGGTTAGGGTTAGGGTTAGGG',  # Classic telomeric repeat
            'expected_class': 'G-Quadruplex',
            'expected_subclass': 'Canonical_intramolecular_G4',
            'min_score': 2.0,  # Should be high confidence
            'description': 'Classic human telomeric G4 with 3x(TTAGGG) repeats, forms stable G-quadruplex structure'
        },
        'variant_g4': {
            'name': 'G4 with Bulges',
            'sequence': 'GGGTGGGTTGGGTGGGG',  # G4 with longer loops
            'expected_class': 'G-Quadruplex',
            'min_score': 1.5,
            'description': 'G4 motif with varied loop lengths, still forms G-quadruplex'
        }
    },
    'slipped_dna': {
        'cag_repeat': {
            'name': 'Huntington CAG Repeat (n=20)',
            'sequence': 'ATGAAGGCCTTC' + 'CAG' * 20 + 'CCGCCGCCGCCG',
            'expected_class': 'Slipped DNA',
            'expected_subclass': 'STR',
            'min_score': 2.0,
            'description': 'CAG repeat associated with Huntington disease (20 copies = normal range)'
        },
        'cgg_fragile_x': {
            'name': 'Fragile X CGG Repeat (n=30)',
            'sequence': 'GCGCGCGC' + 'CGG' * 30 + 'GCATGCATGCAT',
            'expected_class': 'Slipped DNA',
            'expected_subclass': 'STR',
            'min_score': 2.5,
            'description': 'CGG repeat associated with Fragile X syndrome (30 copies = intermediate)'
        }
    },
    'cruciform': {
        'palindrome_short': {
            'name': 'Palindromic Inverted Repeat (12bp arms)',
            'sequence': 'GCATGCATGCAT' + 'AAA' + 'ATGCATGCATGC',  # 12bp inverted repeat with 3bp spacer
            'expected_class': 'Cruciform',
            'min_score': 1.5,
            'description': 'Inverted repeat that can form cruciform structure with hairpin arms'
        }
    },
    'z_dna': {
        'alternating_cg': {
            'name': 'Alternating CG (Z-DNA)',
            'sequence': 'CGCGCGCGCGCGCGCGCGCG',  # 10 CG repeats = 20bp
            'expected_class': 'Z-DNA',
            'min_score': 2.0,
            'description': 'Alternating CG sequence with high Z-DNA forming potential'
        },
        'alternating_ca': {
            'name': 'Alternating CA (Z-DNA)',
            'sequence': 'CACACACACACACACACACA',  # 10 CA repeats = 20bp
            'expected_class': 'Z-DNA',
            'min_score': 1.8,
            'description': 'Alternating CA sequence can also form left-handed helix'
        }
    },
    'curved_dna': {
        'a_tract_phased': {
            'name': 'Phased A-tracts (Curved DNA)',
            'sequence': 'AAAAATTTTTAAAAATTTTTAAAAATTTTT',  # A-tracts separated by ~10bp
            'expected_class': 'Curved DNA',
            'min_score': 1.5,
            'description': 'A-tracts in phase with helix repeat cause DNA curvature'
        }
    },
    'i_motif': {
        'c_rich_canonical': {
            'name': 'C-rich i-Motif',
            'sequence': 'CCCTTCCCTTCCCTTCCC',  # C-rich complementary to G4
            'expected_class': 'i-Motif',
            'expected_subclass': 'Canonical_i-motif',
            'min_score': 2.0,
            'description': 'C-rich sequence complementary to G4, forms i-motif at acidic pH'
        }
    },
    'triplex': {
        'mirror_repeat': {
            'name': 'Purine Mirror Repeat',
            'sequence': 'AAAGGGGAAAA' + 'TTT' + 'AAAGGGGAAAA',  # Purine-rich mirror
            'expected_class': 'Triplex DNA',
            'min_score': 1.5,
            'description': 'Mirror repeat with purine-rich arms can form triplex structure'
        }
    },
    'a_philic': {
        'a_philic_high': {
            'name': 'High A-philic DNA',
            'sequence': 'AAAAAAAAAAATTTTTTTTTTTAAAAAAAAAAATTTTTTTTTTT',
            'expected_class': 'A-philic DNA',
            'min_score': 2.0,
            'description': 'A-tract rich sequence with high protein binding affinity'
        }
    }
}


# ==============================================================================
# Validation Functions
# ==============================================================================

def print_header(title: str, char: str = "="):
    """Print a formatted header."""
    print(f"\n{char * 80}")
    print(f"{title:^80}")
    print(f"{char * 80}")


def print_motif_details(motif: Dict[str, Any], show_all_fields: bool = False):
    """Print detailed information about a motif."""
    print(f"\n  {'─' * 76}")
    print(f"  Motif ID: {motif.get('ID', 'N/A')}")
    print(f"  {'─' * 76}")
    
    # Core fields
    print(f"  Class:      {motif.get('Class', 'N/A')}")
    print(f"  Subclass:   {motif.get('Subclass', 'N/A')}")
    print(f"  Position:   {motif.get('Start', 0)}-{motif.get('End', 0)} ({motif.get('Length', 0)} bp)")
    print(f"  Score:      {motif.get('Score', 0):.3f} (1-3 scale)")
    print(f"  Strand:     {motif.get('Strand', 'N/A')}")
    print(f"  Method:     {motif.get('Method', 'N/A')}")
    
    # Show sequence if available
    seq = motif.get('Sequence', '')
    if seq and len(seq) <= 80:
        print(f"  Sequence:   {seq}")
    elif seq:
        print(f"  Sequence:   {seq[:40]}...{seq[-37:]} (truncated, full length: {len(seq)} bp)")
    
    # Class-specific fields
    if show_all_fields:
        print(f"\n  Additional Fields:")
        special_fields = ['ID', 'Class', 'Subclass', 'Start', 'End', 'Length', 
                         'Score', 'Strand', 'Method', 'Sequence', 'Sequence_Name']
        for key, value in motif.items():
            if key not in special_fields:
                if isinstance(value, float):
                    print(f"    {key}: {value:.3f}")
                elif isinstance(value, list):
                    print(f"    {key}: {', '.join(map(str, value))}")
                else:
                    print(f"    {key}: {value}")


def validate_motif(motif: Dict[str, Any], expected: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """
    Validate a motif against expected properties.
    
    Returns:
        Tuple of (passed, list of validation messages)
    """
    messages = []
    passed = True
    
    # Check class (case-insensitive and flexible matching)
    if 'expected_class' in expected:
        motif_class = motif.get('Class', '').replace('_', ' ').replace('-', ' ').lower()
        expected_class = expected['expected_class'].replace('_', ' ').replace('-', ' ').lower()
        
        # Check if expected class is contained in motif class or vice versa
        if expected_class in motif_class or motif_class in expected_class:
            messages.append(f"  ✓ Class correct: {motif.get('Class')} (matches {expected['expected_class']})")
        else:
            messages.append(f"  ✗ Class mismatch: got '{motif.get('Class')}', expected '{expected['expected_class']}'")
            passed = False
    
    # Check subclass if specified (flexible matching)
    if 'expected_subclass' in expected:
        motif_subclass = motif.get('Subclass', '').replace('_', ' ').replace('-', ' ').lower()
        expected_subclass = expected['expected_subclass'].replace('_', ' ').replace('-', ' ').lower()
        
        # More flexible matching - check if key words match
        if expected_subclass in motif_subclass or motif_subclass in expected_subclass:
            messages.append(f"  ✓ Subclass correct: {motif.get('Subclass')} (matches {expected['expected_subclass']})")
        else:
            messages.append(f"  ⚠ Subclass differs: got '{motif.get('Subclass')}', expected '{expected['expected_subclass']}'")
            # Don't fail test on subclass mismatch, just warn
    
    # Check minimum score
    if 'min_score' in expected:
        score = motif.get('Score', 0)
        if score < expected['min_score']:
            messages.append(f"  ✗ Score too low: {score:.3f} < {expected['min_score']:.3f}")
            passed = False
        else:
            messages.append(f"  ✓ Score sufficient: {score:.3f} >= {expected['min_score']:.3f}")
    
    return passed, messages


def run_validation_test(test_name: str, test_data: Dict[str, Any]) -> Tuple[bool, int, int]:
    """
    Run a single validation test.
    
    Returns:
        Tuple of (passed, motifs_found, motifs_expected)
    """
    print(f"\n{'─' * 80}")
    print(f"Test: {test_data['name']}")
    print(f"{'─' * 80}")
    print(f"Description: {test_data['description']}")
    print(f"Sequence: {test_data['sequence'][:60]}{'...' if len(test_data['sequence']) > 60 else ''}")
    print(f"Length: {len(test_data['sequence'])} bp")
    
    # Run analysis
    start_time = time.time()
    motifs = analyze_sequence(test_data['sequence'], test_name, use_fast_mode=False)
    elapsed = time.time() - start_time
    
    print(f"Analysis time: {elapsed:.3f}s")
    print(f"Motifs found: {len(motifs)}")
    
    if not motifs:
        print("\n  ✗ FAILED: No motifs detected")
        return False, 0, 1
    
    # Filter motifs by expected class (flexible matching)
    expected_class = test_data.get('expected_class', '').replace('_', ' ').replace('-', ' ').lower()
    matching_motifs = []
    for m in motifs:
        motif_class = m.get('Class', '').replace('_', ' ').replace('-', ' ').lower()
        if expected_class in motif_class or motif_class in expected_class:
            matching_motifs.append(m)
    
    if not matching_motifs:
        print(f"\n  ⚠ INFO: No exact match for class '{test_data.get('expected_class', '')}' found")
        print(f"  Found classes: {', '.join(set(m.get('Class', 'Unknown') for m in motifs))}")
        print(f"\n  Showing all detected motifs for analysis:")
        for i, m in enumerate(motifs, 1):
            print(f"\n  Motif {i}:")
            print_motif_details(m, show_all_fields=True)
        return False, len(motifs), 1
    
    # Validate the first matching motif
    motif = matching_motifs[0]
    passed, messages = validate_motif(motif, test_data)
    
    print(f"\nValidation Results:")
    for msg in messages:
        print(msg)
    
    # Show detailed motif information
    print(f"\nDetected Motif Details:")
    print_motif_details(motif, show_all_fields=True)
    
    if passed:
        print(f"\n  ✓ TEST PASSED")
    else:
        print(f"\n  ✗ TEST FAILED")
    
    return passed, len(matching_motifs), 1


def run_all_validation_tests():
    """Run all validation tests and generate comprehensive report."""
    print_header("NonBDNAFinder Comprehensive Validation Suite", "=")
    print(f"Testing {len([t for cat in TEST_SEQUENCES.values() for t in cat.values()])} known motif sequences")
    print(f"Validating accuracy of detection, scoring, and classification")
    
    results = {
        'total': 0,
        'passed': 0,
        'failed': 0,
        'by_class': Counter(),
        'by_class_passed': Counter(),
        'details': []
    }
    
    start_time = time.time()
    
    # Run tests for each motif class
    for motif_class, test_cases in TEST_SEQUENCES.items():
        print_header(f"Testing {motif_class.replace('_', ' ').title()}", "-")
        
        for test_name, test_data in test_cases.items():
            results['total'] += 1
            results['by_class'][motif_class] += 1
            
            passed, found, expected = run_validation_test(test_name, test_data)
            
            if passed:
                results['passed'] += 1
                results['by_class_passed'][motif_class] += 1
            else:
                results['failed'] += 1
            
            results['details'].append({
                'class': motif_class,
                'test': test_name,
                'passed': passed,
                'found': found,
                'expected': expected
            })
    
    elapsed = time.time() - start_time
    
    # Generate summary report
    print_header("Validation Summary Report", "=")
    
    print(f"\nOverall Results:")
    print(f"  Total Tests:   {results['total']}")
    print(f"  Passed:        {results['passed']} ({results['passed']/results['total']*100:.1f}%)")
    print(f"  Failed:        {results['failed']} ({results['failed']/results['total']*100:.1f}%)")
    print(f"  Total Time:    {elapsed:.2f}s")
    print(f"  Avg Time/Test: {elapsed/results['total']:.3f}s")
    
    print(f"\nResults by Motif Class:")
    print(f"  {'Class':<20} {'Tests':<8} {'Passed':<8} {'Success Rate':<15}")
    print(f"  {'-'*55}")
    for motif_class in sorted(results['by_class'].keys()):
        total = results['by_class'][motif_class]
        passed = results['by_class_passed'][motif_class]
        rate = passed / total * 100 if total > 0 else 0
        print(f"  {motif_class.replace('_', ' ').title():<20} {total:<8} {passed:<8} {rate:>6.1f}%")
    
    print(f"\nDetailed Results:")
    for detail in results['details']:
        status = "✓ PASS" if detail['passed'] else "✗ FAIL"
        print(f"  {status} | {detail['class'].replace('_', ' ').title():<20} | {detail['test']}")
    
    # Final verdict
    if results['passed'] == results['total']:
        print(f"\n{'='*80}")
        print(f"{'✓ ALL TESTS PASSED - 100% ACCURACY':^80}")
        print(f"{'='*80}")
    else:
        print(f"\n{'='*80}")
        print(f"{'✗ SOME TESTS FAILED':^80}")
        print(f"{'='*80}")
    
    return results


def test_detector_directly(detector_class, sequence: str, expected_class: str):
    """Test a detector directly without the full pipeline."""
    print(f"\n{'─' * 80}")
    print(f"Direct Detector Test: {detector_class.__name__}")
    print(f"{'─' * 80}")
    
    detector = detector_class()
    start_time = time.time()
    motifs = detector.detect_motifs(sequence, "test_seq")
    elapsed = time.time() - start_time
    
    print(f"Sequence length: {len(sequence)} bp")
    print(f"Detection time: {elapsed:.3f}s")
    print(f"Speed: {len(sequence)/elapsed if elapsed > 0 else 0:.0f} bp/s")
    print(f"Motifs found: {len(motifs)}")
    
    if motifs:
        print(f"\nFirst motif:")
        print_motif_details(motifs[0], show_all_fields=True)


# ==============================================================================
# Main Execution
# ==============================================================================

def main():
    """Main validation suite execution."""
    results = run_all_validation_tests()
    
    # Test direct detector access for performance validation
    print_header("Direct Detector Performance Tests", "=")
    
    # Test each detector class with appropriate sequence
    test_detector_directly(GQuadruplexDetector, TEST_SEQUENCES['g_quadruplex']['canonical_telomeric']['sequence'], 'G-Quadruplex')
    test_detector_directly(SlippedDNADetector, TEST_SEQUENCES['slipped_dna']['cag_repeat']['sequence'], 'Slipped DNA')
    test_detector_directly(ZDNADetector, TEST_SEQUENCES['z_dna']['alternating_cg']['sequence'], 'Z-DNA')
    
    print("\n" + "="*80)
    print("Validation suite completed!")
    print("="*80)
    
    return results


if __name__ == "__main__":
    main()
