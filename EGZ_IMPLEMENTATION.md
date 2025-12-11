# eGZ-motif (Extruded-G Z-DNA) Implementation

## Overview
This document describes the implementation of eGZ-motif detection in the ZDNADetector class.

## What is eGZ-motif?
eGZ-motif (Extruded-G Z-DNA) refers to a specific type of Z-DNA structure characterized by long trinucleotide repeats that form left-handed helices with guanine extrusion. The primary patterns are:

- **(CGG)n** repeats
- **(GGC)n** repeats  
- **(CCG)n** repeats
- **(GCC)n** repeats

These patterns are associated with fragile sites and neurological disorders when expanded (e.g., Fragile X syndrome with CGG repeats).

## Implementation Details

### Files Modified
- **detectors.py**: Updated `ZDNADetector` class to include eGZ-motif detection

### Changes Made

#### 1. Added eGZ Patterns to `get_patterns()` method
```python
"egz_motifs": [
    (r"(?:CGG){3,}", "ZDN_EGZ_CGG", "CGG repeat (eGZ)", "eGZ", 9, "egz_score", 0.85, ...),
    (r"(?:GGC){3,}", "ZDN_EGZ_GGC", "GGC repeat (eGZ)", "eGZ", 9, "egz_score", 0.85, ...),
    (r"(?:CCG){3,}", "ZDN_EGZ_CCG", "CCG repeat (eGZ)", "eGZ", 9, "egz_score", 0.85, ...),
    (r"(?:GCC){3,}", "ZDN_EGZ_GCC", "GCC repeat (eGZ)", "eGZ", 9, "egz_score", 0.85, ...),
]
```

#### 2. Added `_find_egz_motifs()` helper method
This method uses regex patterns to find long trinucleotide repeats:
- Minimum 3 repeats (9 bp) required for detection
- Scores increase with repeat count: `base_score * (repeat_count / 3.0)`
- Returns annotations with repeat unit and count information

#### 3. Updated `annotate_sequence()` method
Now detects both:
- Classic Z-DNA regions (from 10-mer scoring) - labeled as 'Z-DNA' subclass
- eGZ-motif regions (from regex patterns) - labeled as 'eGZ' subclass

#### 4. Updated `detect_motifs()` method
Handles both subclasses differently:
- **Classic Z-DNA**: Reports 10-mer contributions, CG/AT dinucleotides, alternating patterns
- **eGZ-motifs**: Reports repeat unit, repeat count, and GC content

## Detection Examples

### Example 1: CGG Repeat
```python
sequence = "ATGCGGCGGCGGCGGCGGATGC"  # Contains (CGG)5
motifs = detector.detect_motifs(sequence)
# Detects: eGZ-motif with CGG repeat, count=5, score=1.417
```

### Example 2: Mixed Z-DNA
```python
sequence = "GCGCGCGCGCGCATGCGGCGGCGGATGC"  # Both types
motifs = detector.detect_motifs(sequence)
# Detects:
#   1. Classic Z-DNA from GCGCGCGCGCGC
#   2. eGZ-motif from CGGCGGCGG
```

## Scoring Scheme

### Classic Z-DNA
- Based on 10-mer scoring table from Ho et al. 1986
- Score threshold: > 50.0 for detection
- Redistributes scores across overlapping 10-mers

### eGZ-motifs
- Based on repeat count: `base_score * (repeat_count / 3.0)`
- Base score: 0.85 (from pattern threshold)
- Minimum 3 repeats required (score ≈ 0.85)
- Score increases linearly with more repeats
- Example: 6 repeats → score ≈ 1.70

## Output Format

### eGZ-motif Fields
```python
{
    'ID': 'seq_ZDN_EGZ_CGG_45',
    'Class': 'Z-DNA',
    'Subclass': 'eGZ',
    'Start': 45,                # 1-based
    'End': 60,
    'Length': 15,
    'Sequence': 'CGGCGGCGGCGGCGG',
    'Score': 1.417,
    'Repeat_Unit': 'CGG',
    'Repeat_Count': 5,
    'GC_Content': 100.0
}
```

## Testing

### Test File: test_egz_motif.py
Comprehensive test suite covering:
- Detection of all 4 repeat types (CGG, GGC, CCG, GCC)
- Minimum repeat threshold (3 repeats required)
- Long repeat sequences (10+ repeats)
- Score increase with repeat count
- Case insensitivity
- Output format validation
- Mixed detection (both Z-DNA and eGZ)

### Test Results
All 10 tests pass successfully:
```
test_both_zdna_and_egz_detection ... ok
test_case_insensitivity ... ok
test_ccg_repeat_detection ... ok
test_cgg_repeat_detection ... ok
test_gcc_repeat_detection ... ok
test_ggc_repeat_detection ... ok
test_long_cgg_repeat ... ok
test_minimum_repeat_threshold ... ok
test_motif_output_format ... ok
test_score_increases_with_repeats ... ok
```

## References
- Herbert et al. 1997 - Extruded-G Z-DNA
- Ho et al. 1986 - Z-DNA 10-mer scoring table
- Fragile X syndrome - CGG repeat expansion

## Compatibility
- Maintains backward compatibility with existing Z-DNA detection
- No changes to existing API
- Classic Z-DNA detection still uses 10-mer scoring
- eGZ detection is additive, doesn't replace existing functionality

## Performance
- eGZ regex matching is O(n) where n is sequence length
- Minimal performance impact on existing detection
- Tested with sequences up to 50kb with no issues
