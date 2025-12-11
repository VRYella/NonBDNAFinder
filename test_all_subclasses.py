#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║           COMPREHENSIVE SUBCLASS DETECTION TEST SUITE                        ║
║     Rigorous Testing of All Non-B DNA Detector Subclasses                   ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: test_all_subclasses.py
VERSION: 2024.1
LICENSE: MIT

DESCRIPTION:
    Comprehensive test suite that validates detection of ALL subclasses across
    all 9 detector classes. Uses known sequences for each subclass type to
    ensure accurate detection.

TEST COVERAGE:
    ✓ All 9 detector classes
    ✓ All 22+ subclasses
    ✓ Known sequences for each subclass
    ✓ Verification of proper subclass assignment
    ✓ Performance metrics
    ✓ Statistical reporting
"""

import sys
import time
from collections import defaultdict
from typing import Dict, List, Any, Tuple

# Import all detector classes
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


# ============================================================================
# TEST SEQUENCES FOR EACH SUBCLASS
# ============================================================================

TEST_SEQUENCES = {
    # Curved DNA subclasses
    'Curved_DNA:Global Curvature': {
        'sequences': [
            'AAAAATGCATAAAAATGCATAAAAATGCA',  # A-phased repeats
            'AAAACGAAAACGAAAACG',  # A3 phased
            'AAAAAACGAAAAACGAAAAACG',  # A5 phased
        ],
        'min_expected': 1
    },
    'Curved_DNA:Local Curvature': {
        'sequences': [
            'AAAAAAAAAAGCTAGCTAGCT',  # Long A-tract (10 A's)
            'TTTTTTTTTGCTAGCTAGCT',  # Long T-tract
            'GCAAAAAAAATCG',  # A-tract 8bp
        ],
        'min_expected': 1
    },
    
    # Slipped DNA subclasses
    'Slipped_DNA:Direct Repeat': {
        'sequences': [
            'ATGATGATGATGATGATG',  # ATG direct repeat
            'GCTAGCTAGCTAGCT',  # GCTA direct repeat
            'TTATTATTATTATTA',  # TTA repeat
        ],
        'min_expected': 1
    },
    'Slipped_DNA:STR': {
        'sequences': [
            'CAG' * 20,  # Huntington's CAG repeat
            'CGG' * 15,  # Fragile X CGG repeat
            'CTG' * 18,  # Myotonic dystrophy repeat
        ],
        'min_expected': 1
    },
    
    # Cruciform subclasses
    'Cruciform:Inverted Repeats': {
        'sequences': [
            'ATCGATGCTAGCATCGAT',  # Palindromic inverted repeat
            'GCGCTAGCTACGCGC',  # Inverted repeat
            'ATGCTAGCTACGCAT',  # Perfect palindrome
        ],
        'min_expected': 1
    },
    
    # R-Loop subclasses
    'R-Loop:QmRLFS-m1': {
        'sequences': [
            'GGGGGGGGGGGGGGGGGGGGG' + 'ATCG' * 10,  # High G% RIZ
            'G' * 50 + 'ATCG' * 20,  # G-rich region
        ],
        'min_expected': 1
    },
    'R-Loop:QmRLFS-m2': {
        'sequences': [
            'G' * 30 + 'ATCG' * 10 + 'G' * 30,  # RIZ + REZ pattern
            'GGG' * 15 + 'ATCG' * 5 + 'GGG' * 15,
        ],
        'min_expected': 1
    },
    
    # Triplex subclasses
    'Triplex:Homopurine mirror repeat': {
        'sequences': [
            'A' * 15 + 'TT' + 'A' * 15,  # Homopurine mirror
            'G' * 12 + 'AAA' + 'G' * 12,  # G-rich mirror
            'AGAGAGAGAGAG' + 'TT' + 'AGAGAGAGAGAG',  # Purine-rich mirror
        ],
        'min_expected': 1
    },
    'Triplex:Homopyrimidine mirror repeat': {
        'sequences': [
            'C' * 15 + 'AA' + 'C' * 15,  # Homopyrimidine mirror
            'T' * 12 + 'GGG' + 'T' * 12,  # T-rich mirror
            'CTCTCTCTCTCT' + 'AA' + 'CTCTCTCTCTCT',  # Pyrimidine-rich mirror
        ],
        'min_expected': 1
    },
    'Triplex:Sticky_DNA': {
        'sequences': [
            'GAA' * 10,  # GAA repeat
            'TTC' * 12,  # TTC repeat
            'GAAGAAGAAGAAGAA',  # Extended GAA
        ],
        'min_expected': 1
    },
    
    # G-Quadruplex subclasses
    'G-Quadruplex:Canonical G4': {
        'sequences': [
            'GGGTTAGGGTTAGGGTTAGGG',  # Classic G4
            'GGGAGGGAGGGAGGG',  # Tight G4
            'GGGTTGGGTTGGGTTGGG',  # Medium loop G4
        ],
        'min_expected': 1
    },
    'G-Quadruplex:Relaxed G4': {
        'sequences': [
            'GGATTGGATTGGATTGG',  # GG-based G4
            'GGATATATGGATTTTGGAAGG',  # Longer loops
        ],
        'min_expected': 1
    },
    'G-Quadruplex:Long-loop G4': {
        'sequences': [
            'GGGATATATATATAGGGTTGGGTAGGG',  # Long first loop
            'GGGATTTTTTTTGGGTTGGGTAGGG',  # Extended loop G4
        ],
        'min_expected': 1
    },
    'G-Quadruplex:Bulged G4': {
        'sequences': [
            'GGAGATGGGTTGGGTAGGG',  # Bulge in first tract
            'GGGATTGGAGTTGGGTAGGG',  # Bulge in second tract
        ],
        'min_expected': 1
    },
    'G-Quadruplex:Multimeric G4': {
        'sequences': [
            'GGGTTGGGTTGGGTTGGGTTTGGGTTGGGTTGGGTTGGG',  # Multiple G4 units
            'GGGTGGGGTGGGGTGGGG' * 2,  # Repeated G4 pattern
        ],
        'min_expected': 1
    },
    'G-Quadruplex:Imperfect G4': {
        'sequences': [
            'GGGATTGGAATTGGGTTGGG',  # Imperfect tracts
            'GGATTGGATTTGGATTGG',  # Interruptions
        ],
        'min_expected': 1
    },
    'G-Quadruplex:G-Triplex intermediate': {
        'sequences': [
            'GGGTTGGGTTGGG',  # 3-tract G structure
            'GGGATGGGTAGGG',  # G-triplex
        ],
        'min_expected': 1
    },
    
    # i-Motif subclasses
    'i-Motif:Canonical i-Motif': {
        'sequences': [
            'CCCTTACCCTTACCCTTACCC',  # Classic i-motif
            'CCCACCCACCCACCC',  # Tight i-motif
        ],
        'min_expected': 1
    },
    'i-Motif:HUR AC-Motif': {
        'sequences': [
            'AAATTCCCTTTCCCTTTCCC',  # AC-motif pattern
            'CCCTTTTCCCTTTTCCCTTTTAAA',  # CA-motif pattern
        ],
        'min_expected': 1
    },
    
    # Z-DNA subclasses
    'Z-DNA:Z-DNA': {
        'sequences': [
            'CGCGCGCGCGCGCGCG',  # Alternating CG
            'GTATATATATAC',  # Alternating GT/AC
        ],
        'min_expected': 1
    },
    'Z-DNA:eGZ': {
        'sequences': [
            'CGG' * 5,  # CGG repeat (eGZ)
            'GGC' * 5,  # GGC repeat (eGZ)
            'CCG' * 5,  # CCG repeat (eGZ)
        ],
        'min_expected': 1
    },
    
    # A-philic DNA
    'A-philic_DNA:A-philic DNA': {
        'sequences': [
            'AAAATTTTAAAATTTT',  # A-philic pattern
            'TAAAATAAAAATAAA',  # T/A-philic
        ],
        'min_expected': 1
    },
}


# ============================================================================
# TEST RUNNER
# ============================================================================

class SubclassDetectionTester:
    """Comprehensive tester for all detector subclasses"""
    
    def __init__(self):
        self.detectors = {
            'Curved_DNA': CurvedDNADetector(),
            'Slipped_DNA': SlippedDNADetector(),
            'Cruciform': CruciformDetector(),
            'R-Loop': RLoopDetector(),
            'Triplex': TriplexDetector(),
            'G-Quadruplex': GQuadruplexDetector(),
            'i-Motif': IMotifDetector(),
            'Z-DNA': ZDNADetector(),
            'A-philic_DNA': APhilicDetector(),
        }
        
        self.results = defaultdict(dict)
        self.performance = {}
    
    def test_all_subclasses(self) -> Dict[str, Any]:
        """Test all subclasses and collect results"""
        print("\n" + "=" * 80)
        print("COMPREHENSIVE SUBCLASS DETECTION TEST")
        print("=" * 80)
        
        total_tests = 0
        passed_tests = 0
        failed_tests = 0
        
        for subclass_key, test_data in TEST_SEQUENCES.items():
            class_name, subclass_name = subclass_key.split(':', 1)
            
            if class_name not in self.detectors:
                print(f"\n⚠️  Skipping {subclass_key} - detector not found")
                continue
            
            print(f"\n{'─' * 80}")
            print(f"Testing: {subclass_key}")
            print(f"{'─' * 80}")
            
            detector = self.detectors[class_name]
            sequences = test_data['sequences']
            min_expected = test_data['min_expected']
            
            detected = False
            detection_count = 0
            test_times = []
            
            for i, sequence in enumerate(sequences, 1):
                start_time = time.time()
                motifs = detector.detect_motifs(sequence, f"test_seq_{i}")
                elapsed = time.time() - start_time
                test_times.append(elapsed)
                
                # Check if correct subclass was detected
                subclass_motifs = [m for m in motifs 
                                  if m['Class'] == class_name and 
                                  subclass_name in m['Subclass']]
                
                if subclass_motifs:
                    detected = True
                    detection_count += len(subclass_motifs)
                    print(f"  ✅ Sequence {i}: Detected {len(subclass_motifs)} motifs")
                    # Show first motif details
                    if subclass_motifs:
                        m = subclass_motifs[0]
                        print(f"      Position: {m['Start']}-{m['End']}, "
                              f"Score: {m['Score']:.3f}, "
                              f"Subclass: {m['Subclass']}")
                else:
                    print(f"  ❌ Sequence {i}: No detection")
                    # Show what was detected instead
                    if motifs:
                        print(f"      Found {len(motifs)} motifs of other types:")
                        for m in motifs[:2]:  # Show first 2
                            print(f"        - {m['Class']}: {m['Subclass']}")
            
            # Evaluate test result
            total_tests += 1
            if detected and detection_count >= min_expected:
                passed_tests += 1
                status = "✅ PASS"
            else:
                failed_tests += 1
                status = "❌ FAIL"
            
            avg_time = sum(test_times) / len(test_times) if test_times else 0
            
            print(f"\n  {status}: {detection_count} detections (expected >= {min_expected})")
            print(f"  Average time: {avg_time*1000:.2f} ms")
            
            # Store results
            self.results[subclass_key] = {
                'detected': detected,
                'count': detection_count,
                'expected': min_expected,
                'passed': detected and detection_count >= min_expected,
                'avg_time': avg_time
            }
        
        # Print summary
        print("\n" + "=" * 80)
        print("TEST SUMMARY")
        print("=" * 80)
        print(f"Total subclasses tested: {total_tests}")
        print(f"Passed: {passed_tests} ({passed_tests/total_tests*100:.1f}%)")
        print(f"Failed: {failed_tests} ({failed_tests/total_tests*100:.1f}%)")
        print("=" * 80)
        
        if failed_tests > 0:
            print("\nFailed subclasses:")
            for subclass_key, result in self.results.items():
                if not result['passed']:
                    print(f"  ❌ {subclass_key}: {result['count']} detected "
                          f"(expected >= {result['expected']})")
        
        print("\n" + "=" * 80)
        
        return {
            'total': total_tests,
            'passed': passed_tests,
            'failed': failed_tests,
            'success_rate': passed_tests / total_tests if total_tests > 0 else 0,
            'details': dict(self.results)
        }
    
    # Report date format constant
    DATE_FORMAT = '%Y-%m-%d %H:%M:%S'
    
    def generate_report(self) -> str:
        """Generate detailed markdown report"""
        report = []
        report.append("# Subclass Detection Test Report\n")
        report.append(f"**Test Date:** {time.strftime(self.DATE_FORMAT)}\n")
        report.append("\n## Summary\n")
        
        total = len(self.results)
        passed = sum(1 for r in self.results.values() if r['passed'])
        
        report.append(f"- **Total Subclasses Tested:** {total}\n")
        report.append(f"- **Passed:** {passed} ({passed/total*100:.1f}%)\n")
        report.append(f"- **Failed:** {total - passed}\n")
        
        report.append("\n## Detailed Results\n")
        report.append("\n| Subclass | Status | Detections | Expected | Avg Time (ms) |\n")
        report.append("|----------|--------|------------|----------|---------------|\n")
        
        for subclass_key in sorted(self.results.keys()):
            result = self.results[subclass_key]
            status = "✅ PASS" if result['passed'] else "❌ FAIL"
            report.append(f"| {subclass_key} | {status} | "
                         f"{result['count']} | {result['expected']} | "
                         f"{result['avg_time']*1000:.2f} |\n")
        
        return ''.join(report)


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Run comprehensive subclass detection tests"""
    tester = SubclassDetectionTester()
    
    # Run tests
    results = tester.test_all_subclasses()
    
    # Generate and save report
    report = tester.generate_report()
    
    # Save report to file
    report_file = '/tmp/subclass_detection_report.md'
    with open(report_file, 'w') as f:
        f.write(report)
    
    print(f"\n📄 Detailed report saved to: {report_file}")
    
    # Exit with appropriate code
    if results['failed'] == 0:
        print("\n✅ ALL SUBCLASS TESTS PASSED!")
        return 0
    else:
        print(f"\n❌ {results['failed']} SUBCLASS TESTS FAILED")
        return 1


if __name__ == '__main__':
    sys.exit(main())
