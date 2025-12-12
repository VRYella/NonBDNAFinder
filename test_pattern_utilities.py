#!/usr/bin/env python3
"""
Test script for pattern export and Hyperscan DB builder utilities.

Tests:
1. Export patterns to TSV (default behavior)
2. Export patterns to TSV with R-loop
3. Export patterns to TSV with empty patterns
4. Verify TSV file structure and content
5. Test Hyperscan DB builder (if hyperscan is available)
"""

import os
import sys
import csv
from pathlib import Path

# Add current directory to path
sys.path.insert(0, os.path.dirname(__file__))

from export_patterns_full import (
    extract_all_patterns,
    test_python_re_compile,
    test_hyperscan_compatible,
    export_to_tsv
)

try:
    from build_hyperscan_db import (
        load_patterns_from_tsv,
        filter_hyperscan_compatible,
        test_hyperscan_compatible as hs_test
    )
    BUILD_AVAILABLE = True
except Exception as e:
    BUILD_AVAILABLE = False
    print(f"Warning: build_hyperscan_db tests disabled: {e}")


def test_pattern_extraction():
    """Test pattern extraction from registry."""
    print("\n" + "=" * 70)
    print("TEST 1: Pattern Extraction")
    print("=" * 70)
    
    # Test default extraction (no R-loop, no empty)
    patterns = extract_all_patterns(include_rloop=False, include_empty=False)
    print(f"✓ Extracted {len(patterns)} patterns (default)")
    assert len(patterns) > 0, "Should extract at least some patterns"
    
    # Test with R-loop
    patterns_with_rloop = extract_all_patterns(include_rloop=True, include_empty=False)
    print(f"✓ Extracted {len(patterns_with_rloop)} patterns (with R-loop)")
    assert len(patterns_with_rloop) >= len(patterns), "Should have more patterns with R-loop"
    
    # Test with empty patterns
    patterns_with_empty = extract_all_patterns(include_rloop=False, include_empty=True)
    print(f"✓ Extracted {len(patterns_with_empty)} patterns (with empty)")
    assert len(patterns_with_empty) >= len(patterns), "Should have more patterns with empty"
    
    # Verify pattern structure
    first_pattern = patterns[0]
    required_fields = [
        'Class', 'Subclass', 'Pattern_Group', 'Pattern', 'Pattern_ID',
        'Name', 'Min_Length', 'Score_Type', 'Base_Score',
        'Python_re_compiles', 'Hyperscan_compatible',
        'Description', 'Reference'
    ]
    for field in required_fields:
        assert field in first_pattern, f"Pattern missing required field: {field}"
    print(f"✓ All required fields present in patterns")
    
    return patterns


def test_pattern_compilation():
    """Test pattern compilation checks."""
    print("\n" + "=" * 70)
    print("TEST 2: Pattern Compilation")
    print("=" * 70)
    
    # Test valid patterns
    valid_patterns = [
        r'A{7,}',
        r'(?:CGG){3,}',
        r'G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}',
    ]
    
    for pattern in valid_patterns:
        result = test_python_re_compile(pattern)
        print(f"✓ Pattern compiles: {pattern[:50]}... = {result}")
        assert result, f"Valid pattern should compile: {pattern}"
    
    # Test invalid patterns
    invalid_patterns = [
        '',  # Empty
        r'(?P<invalid)',  # Incomplete named group
    ]
    
    for pattern in invalid_patterns:
        result = test_python_re_compile(pattern)
        print(f"✓ Invalid pattern detected: {repr(pattern)} = {result}")
        assert not result, f"Invalid pattern should not compile: {pattern}"


def test_tsv_export():
    """Test TSV export functionality."""
    print("\n" + "=" * 70)
    print("TEST 3: TSV Export")
    print("=" * 70)
    
    patterns = extract_all_patterns(include_rloop=False, include_empty=False)
    
    # Export to temporary file
    output_path = Path('/tmp/test_export.tsv')
    export_to_tsv(patterns, output_path)
    
    assert output_path.exists(), "TSV file should be created"
    print(f"✓ TSV file created: {output_path}")
    
    # Verify file content
    with open(output_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)
    
    assert len(rows) == len(patterns), "TSV should contain all patterns"
    print(f"✓ TSV contains {len(rows)} rows")
    
    # Verify first row
    first_row = rows[0]
    assert 'Class' in first_row, "TSV should have Class column"
    assert 'Pattern' in first_row, "TSV should have Pattern column"
    assert 'Pattern_ID' in first_row, "TSV should have Pattern_ID column"
    print(f"✓ TSV structure verified")
    
    return output_path


def test_hyperscan_compatibility():
    """Test Hyperscan compatibility checks."""
    print("\n" + "=" * 70)
    print("TEST 4: Hyperscan Compatibility")
    print("=" * 70)
    
    patterns = extract_all_patterns(include_rloop=False, include_empty=False)
    
    compatible_count = 0
    incompatible_count = 0
    
    for pattern in patterns:
        if pattern['Hyperscan_compatible'] == 'Yes':
            compatible_count += 1
        else:
            incompatible_count += 1
    
    print(f"✓ Compatible patterns: {compatible_count}")
    print(f"✓ Incompatible patterns: {incompatible_count}")
    assert compatible_count > 0, "Should have at least some compatible patterns"


def test_tsv_loading():
    """Test loading patterns from TSV."""
    if not BUILD_AVAILABLE:
        print("\n" + "=" * 70)
        print("TEST 5: TSV Loading (SKIPPED - build module not available)")
        print("=" * 70)
        return
    
    print("\n" + "=" * 70)
    print("TEST 5: TSV Loading")
    print("=" * 70)
    
    # Export patterns first
    patterns = extract_all_patterns(include_rloop=False, include_empty=False)
    output_path = Path('/tmp/test_load.tsv')
    export_to_tsv(patterns, output_path)
    
    # Load patterns back
    loaded_patterns = load_patterns_from_tsv(output_path)
    
    assert len(loaded_patterns) == len(patterns), "Should load all patterns"
    print(f"✓ Loaded {len(loaded_patterns)} patterns from TSV")
    
    # Verify structure
    first_loaded = loaded_patterns[0]
    assert 'Pattern' in first_loaded, "Loaded pattern should have Pattern field"
    assert 'Pattern_ID' in first_loaded, "Loaded pattern should have Pattern_ID field"
    print(f"✓ Loaded pattern structure verified")


def test_filter_compatible():
    """Test filtering Hyperscan-compatible patterns."""
    if not BUILD_AVAILABLE:
        print("\n" + "=" * 70)
        print("TEST 6: Filter Compatible (SKIPPED - build module not available)")
        print("=" * 70)
        return
    
    print("\n" + "=" * 70)
    print("TEST 6: Filter Compatible")
    print("=" * 70)
    
    patterns = extract_all_patterns(include_rloop=False, include_empty=False)
    
    compatible, incompatible = filter_hyperscan_compatible(patterns)
    
    print(f"✓ Compatible patterns: {len(compatible)}")
    print(f"✓ Incompatible patterns: {len(incompatible)}")
    
    assert len(compatible) + len(incompatible) == len(patterns), "Should categorize all patterns"


def main():
    """Run all tests."""
    print("=" * 70)
    print("PATTERN UTILITIES TEST SUITE")
    print("=" * 70)
    
    try:
        test_pattern_extraction()
        test_pattern_compilation()
        test_tsv_export()
        test_hyperscan_compatibility()
        test_tsv_loading()
        test_filter_compatible()
        
        print("\n" + "=" * 70)
        print("✓ ALL TESTS PASSED")
        print("=" * 70)
        
        return 0
        
    except Exception as e:
        print("\n" + "=" * 70)
        print(f"✗ TEST FAILED: {e}")
        print("=" * 70)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
