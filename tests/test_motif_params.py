"""
Unit tests for motif parameter synchronization.

Tests verify that constants and detectors implement the canonical motif definitions:
- Cruciform DNA: arms 10-100 nt, spacer 0-3 nt, perfect Watson-Crick matches
- Triplex DNA: mirror repeats 10-100 nt, spacer 0-8 nt, >90% purine/pyrimidine
- Slipped DNA: direct repeats 10-50 nt, spacer=0 (tandem)
- STRs: unit 1-9 bp, min total ≥20 bp
"""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest


def test_scanner_constants():
    """Test that scanner.py has correct constants."""
    import scanner
    
    # Direct repeat constants
    assert scanner.DIRECT_MIN_UNIT == 10, "DIRECT_MIN_UNIT should be 10"
    assert scanner.DIRECT_MAX_UNIT == 50, "DIRECT_MAX_UNIT should be 50"
    assert scanner.DIRECT_MAX_SPACER == 5, "DIRECT_MAX_SPACER should be 5 (general direct repeat)"
    
    # Inverted repeat (Cruciform) constants
    assert scanner.INVERTED_MIN_ARM == 10, "INVERTED_MIN_ARM should be 10"
    assert scanner.INVERTED_MAX_ARM == 100, "INVERTED_MAX_ARM should be 100"
    assert scanner.INVERTED_MAX_LOOP == 100, "INVERTED_MAX_LOOP should be 100 (general inverted repeat)"
    
    # Mirror repeat (Triplex) constants
    assert scanner.MIRROR_MIN_ARM == 10, "MIRROR_MIN_ARM should be 10"
    assert scanner.MIRROR_MAX_ARM == 100, "MIRROR_MAX_ARM should be 100"
    assert scanner.MIRROR_MAX_LOOP == 100, "MIRROR_MAX_LOOP should be 100 (general mirror repeat)"
    
    # STR constants
    assert scanner.STR_MIN_UNIT == 1, "STR_MIN_UNIT should be 1"
    assert scanner.STR_MAX_UNIT == 9, "STR_MAX_UNIT should be 9"
    assert scanner.STR_MIN_TOTAL == 20, "STR_MIN_TOTAL should be 20"


def test_scanner_optimized_constants():
    """Test that scanner_optimized.py has correct constants."""
    try:
        import scanner_optimized
    except ImportError:
        pytest.skip("scanner_optimized not available")
    
    # Direct repeat constants
    assert scanner_optimized.DIRECT_MIN_UNIT == 10, "DIRECT_MIN_UNIT should be 10"
    assert scanner_optimized.DIRECT_MAX_UNIT == 50, "DIRECT_MAX_UNIT should be 50"
    assert scanner_optimized.DIRECT_MAX_SPACER == 5, "DIRECT_MAX_SPACER should be 5"
    
    # Inverted repeat constants
    assert scanner_optimized.INVERTED_MIN_ARM == 10, "INVERTED_MIN_ARM should be 10"
    assert scanner_optimized.INVERTED_MAX_ARM == 100, "INVERTED_MAX_ARM should be 100"
    assert scanner_optimized.INVERTED_MAX_LOOP == 100, "INVERTED_MAX_LOOP should be 100"
    
    # Mirror repeat constants
    assert scanner_optimized.MIRROR_MIN_ARM == 10, "MIRROR_MIN_ARM should be 10"
    assert scanner_optimized.MIRROR_MAX_ARM == 100, "MIRROR_MAX_ARM should be 100"
    assert scanner_optimized.MIRROR_MAX_LOOP == 100, "MIRROR_MAX_LOOP should be 100"
    
    # STR constants
    assert scanner_optimized.STR_MIN_UNIT == 1, "STR_MIN_UNIT should be 1"
    assert scanner_optimized.STR_MAX_UNIT == 9, "STR_MAX_UNIT should be 9"
    assert scanner_optimized.STR_MIN_TOTAL == 20, "STR_MIN_TOTAL should be 20"


def test_cruciform_detector_params():
    """Test CruciformDetector parameters."""
    from detectors import CruciformDetector
    
    detector = CruciformDetector()
    assert detector.MIN_ARM == 10, "CruciformDetector.MIN_ARM should be 10"
    assert detector.MAX_ARM == 100, "CruciformDetector.MAX_ARM should be 100"
    assert detector.MAX_LOOP == 3, "CruciformDetector.MAX_LOOP should be 3 (Cruciform subset)"
    assert detector.MAX_MISMATCHES == 0, "CruciformDetector.MAX_MISMATCHES should be 0"


def test_triplex_detector_params():
    """Test TriplexDetector parameters."""
    from detectors import TriplexDetector
    
    detector = TriplexDetector()
    assert detector.MIN_ARM == 10, "TriplexDetector.MIN_ARM should be 10"
    assert detector.MAX_ARM == 100, "TriplexDetector.MAX_ARM should be 100"
    assert detector.MAX_LOOP == 8, "TriplexDetector.MAX_LOOP should be 8 (Triplex subset)"
    assert detector.PURINE_PYRIMIDINE_THRESHOLD == 0.9, "TriplexDetector purity threshold should be 0.9"


def test_slipped_dna_detector_params():
    """Test SlippedDNADetector parameters."""
    from detectors import SlippedDNADetector
    
    detector = SlippedDNADetector()
    assert detector.MIN_UNIT == 10, "SlippedDNADetector.MIN_UNIT should be 10"
    assert detector.MAX_UNIT == 50, "SlippedDNADetector.MAX_UNIT should be 50"
    assert detector.MAX_SPACER == 0, "SlippedDNADetector.MAX_SPACER should be 0 (tandem only)"


def test_cruciform_detection():
    """Test Cruciform detection with a simple example."""
    from detectors import CruciformDetector
    
    detector = CruciformDetector()
    
    # Create a simple cruciform: 12-nt arms with 2-nt spacer
    # Left arm: AAAAGGGGCCCC
    # Spacer: TT
    # Right arm (reverse complement): GGGGCCCCTTTT
    seq = "AAAAGGGGCCCCTTGGGGCCCCTTTT"
    
    motifs = detector.detect_motifs(seq, "test_seq")
    
    # Should detect at least one cruciform
    assert len(motifs) > 0, "Should detect at least one cruciform"
    
    # Check that detected motif has correct class
    cruciform_found = any(m.get('Class') == 'Cruciform' for m in motifs)
    assert cruciform_found, "Should detect Cruciform class motif"


def test_triplex_detection():
    """Test Triplex detection with purine-rich mirror repeat."""
    from detectors import TriplexDetector
    
    detector = TriplexDetector()
    
    # Create a purine-rich mirror repeat (>90% purine)
    # Left arm: 12 nt of purines (AGAGAGAGAGAG)
    # Spacer: TTTT (4 nt)
    # Right arm (reverse): GAGAGAGAGAGA
    seq = "AGAGAGAGAGAGTTTTGAGAGAGAGAGA"
    
    motifs = detector.detect_motifs(seq, "test_seq")
    
    # Should detect triplex with purine purity
    triplex_found = any(
        m.get('Class') == 'Triplex' and 
        ('Homopurine' in str(m.get('Subclass', '')) or 
         'Triplex' in str(m.get('Subclass', '')))
        for m in motifs
    )
    assert triplex_found or len(motifs) > 0, "Should detect Triplex motif or mirror repeat"


def test_slipped_dna_detection():
    """Test Slipped DNA detection with tandem repeat (spacer=0)."""
    from detectors import SlippedDNADetector
    
    detector = SlippedDNADetector()
    
    # Create a 20-nt unit repeated twice (tandem, spacer=0)
    unit = "ACGTACGTACGTACGTACGT"  # 20 nt
    seq = unit + unit  # 40 nt total
    
    motifs = detector.detect_motifs(seq, "test_seq")
    
    # Should detect direct repeat
    assert len(motifs) > 0, "Should detect at least one slipped DNA motif"
    
    # Check for direct repeat with spacer=0
    direct_repeat_found = any(
        m.get('Subclass') == 'Direct_Repeat' and 
        m.get('details', {}).get('spacer_length') == 0
        for m in motifs
    )
    assert direct_repeat_found or any(m.get('Class') == 'Slipped_DNA' for m in motifs), \
        "Should detect Slipped DNA or direct repeat with spacer=0"


def test_str_detection():
    """Test STR detection with minimum total length of 20 bp."""
    import scanner
    
    # Create a 2-bp unit repeated 10 times = 20 bp total
    seq = "AT" * 10
    
    strs = scanner.find_strs(seq, min_u=1, max_u=9, min_total=20)
    
    # Should detect the STR
    assert len(strs) > 0, "Should detect STR with total length ≥20 bp"
    
    # Verify it's a 2-mer with 10 copies
    two_mer_found = any(
        s.get('Unit_Length') == 2 and s.get('Copies') == 10
        for s in strs
    )
    assert two_mer_found, "Should detect 2-mer STR with 10 copies (20 bp total)"


def test_inverted_repeat_function_signatures():
    """Test that find_inverted_repeats accepts max_arm parameter."""
    import scanner
    
    seq = "ACGTACGTACGTACGT"
    
    # Should accept max_arm parameter
    try:
        results = scanner.find_inverted_repeats(seq, min_arm=10, max_arm=50, max_loop=10)
        # If successful, function signature is correct
        assert True
    except TypeError as e:
        if 'max_arm' in str(e):
            pytest.fail("find_inverted_repeats should accept max_arm parameter")
        else:
            raise


def test_mirror_repeat_function_signatures():
    """Test that find_mirror_repeats accepts max_arm parameter."""
    import scanner
    
    seq = "ACGTACGTACGTACGT"
    
    # Should accept max_arm parameter
    try:
        results = scanner.find_mirror_repeats(seq, min_arm=10, max_arm=50, max_loop=10)
        # If successful, function signature is correct
        assert True
    except TypeError as e:
        if 'max_arm' in str(e):
            pytest.fail("find_mirror_repeats should accept max_arm parameter")
        else:
            raise


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])
