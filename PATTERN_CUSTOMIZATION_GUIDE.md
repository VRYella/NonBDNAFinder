# Pattern Customization Guide

## Overview

This guide explains how to customize the regular expression patterns used for Non-B DNA motif detection in the NBDScanner application.

## Where to Find the Patterns

All regex patterns and motif definitions are centralized in the **`motif_patterns.py`** file. This file contains pattern definitions for all 9 detector classes:

1. **CurvedDNADetector** - A-tract mediated DNA curvature
2. **ZDNADetector** - Z-DNA and eGZ motifs
3. **APhilicDetector** - A-philic DNA tetranucleotides
4. **SlippedDNADetector** - STRs and Direct Repeats
5. **CruciformDetector** - Inverted Repeats (Palindromes)
6. **RLoopDetector** - QmRLFS patterns
7. **TriplexDetector** - Sticky DNA and Mirror Repeats
8. **GQuadruplexDetector** - 7 G4 subclasses
9. **IMotifDetector** - Canonical and HUR AC-motifs

## Pattern Structure

Each pattern is defined as a tuple with 9 elements:

```python
(regex_pattern, pattern_id, name, subclass, min_length, score_type, 
 base_score, description, reference)
```

### Element Descriptions:

- **regex_pattern**: The regular expression to match (e.g., `r'G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}'`)
- **pattern_id**: Unique identifier (e.g., `'G4_0'`)
- **name**: Human-readable name (e.g., `'Canonical G4'`)
- **subclass**: Motif subclass (e.g., `'Canonical G4'`)
- **min_length**: Minimum length for valid matches (e.g., `15`)
- **score_type**: Scoring method name (e.g., `'g4hunter_score'`)
- **base_score**: Base score value (e.g., `0.95`)
- **description**: Description of the motif (e.g., `'Stable G4 structures'`)
- **reference**: Literature citation (e.g., `'Burge 2006'`)

## How to Modify Patterns

### Example 1: Modifying an Existing G4 Pattern

Let's say you want to change the canonical G4 pattern to require 4+ guanines instead of 3+:

**Before:**
```python
'canonical_g4': [
    (r'G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}', 
     'G4_0', 'Canonical G4', 'Canonical G4', 15, 'g4hunter_score', 
     0.95, 'Stable G4 structures', 'Burge 2006'),
],
```

**After:**
```python
'canonical_g4': [
    (r'G{4,}[ACGT]{1,7}G{4,}[ACGT]{1,7}G{4,}[ACGT]{1,7}G{4,}', 
     'G4_0', 'Canonical G4', 'Canonical G4', 16, 'g4hunter_score', 
     0.95, 'Stable G4 structures with 4+ G-runs', 'Modified from Burge 2006'),
],
```

### Example 2: Adding a New Pattern

To add a new G4 variant pattern, add it to the appropriate pattern group:

```python
'canonical_g4': [
    # Existing pattern
    (r'G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}', 
     'G4_0', 'Canonical G4', 'Canonical G4', 15, 'g4hunter_score', 
     0.95, 'Stable G4 structures', 'Burge 2006'),
    
    # NEW: Add your custom pattern here
    (r'G{5,}[ACGT]{1,5}G{5,}[ACGT]{1,5}G{5,}[ACGT]{1,5}G{5,}', 
     'G4_CUSTOM', 'Ultra-stable G4', 'Canonical G4', 20, 'g4hunter_score', 
     0.98, 'Very stable G4 with long G-runs', 'Custom 2024'),
],
```

### Example 3: Modifying Z-DNA Patterns

To change the eGZ motif patterns:

```python
'egz_motifs': [
    # Change minimum repeat count from 3 to 4
    (r'(?:CGG){4,}', 'ZDN_EGZ_CGG', 'CGG repeat (eGZ)', 'eGZ', 12, 
     'egz_score', 0.90, 'Extruded-G Z-DNA CGG repeat', 'Herbert 1997'),
    # ... other patterns
],
```

### Example 4: Adjusting i-Motif Loop Lengths

To modify HUR AC-motif linker lengths:

```python
'hur_ac_motif': [
    # Change from 4bp to 3bp linker
    (r'A{3}[ACGT]{3}C{3}[ACGT]{3}C{3}[ACGT]{3}C{3}', 
     'HUR_AC_CUSTOM', 'HUR AC-motif (3bp)', 'AC-motif (HUR)', 15, 'ac_motif_score', 
     0.85, 'HUR AC alternating motif with shorter linker', 'Custom 2024'),
],
```

## Pattern Groups by Detector

### G-Quadruplex Patterns (GQuadruplexDetector)

Located in `G4_PATTERNS` dictionary with these groups:
- `canonical_g4` - Classic G4 structure
- `relaxed_g4` - Looser G4 requirements
- `long_loop_g4` - Extended first loop
- `bulged_g4` - G4 with bulge loops
- `multimeric_g4` - Multiple G4 units
- `imperfect_g4` - G4-like with interruptions
- `g_triplex` - Three G-tract structures

### Z-DNA Patterns (ZDNADetector)

Located in `ZDNA_PATTERNS` dictionary:
- `z_dna_10mers` - Algorithmic detection (not regex-based)
- `egz_motifs` - Extruded-G Z-DNA repeats (CGG, GGC, CCG, GCC)

### i-Motif Patterns (IMotifDetector)

Located in `IMOTIF_PATTERNS` dictionary:
- `canonical_imotif` - Classic C-rich i-motif
- `hur_ac_motif` - HUR AC alternating patterns (6 variants)

### R-Loop Patterns (RLoopDetector)

Located in `RLOOP_PATTERNS` dictionary:
- `qmrlfs_model_1` - 3+ G tracts
- `qmrlfs_model_2` - 4+ G tracts (higher confidence)

### Triplex Patterns (TriplexDetector)

Located in `TRIPLEX_PATTERNS` dictionary:
- `triplex_forming_sequences` - GAA and TTC repeats

## Testing Your Changes

After modifying patterns, test them:

```bash
# Test pattern module loads correctly
python3 -c "import motif_patterns; print(motif_patterns.pattern_info_summary())"

# Test detectors use new patterns
python3 -c "from detectors import GQuadruplexDetector; g4 = GQuadruplexDetector(); print(g4.get_patterns())"

# Test full analysis pipeline
python3 -c "from nonbscanner import analyze_sequence; result = analyze_sequence('GGGTTAGGGTTAGGGTTAGGG', 'test'); print(f'Found {len(result)} motifs')"
```

## Helper Functions

The `motif_patterns.py` module provides several helper functions:

```python
from motif_patterns import (
    get_patterns_for_detector,  # Get patterns for a specific detector
    get_all_pattern_ids,        # List all pattern IDs
    get_pattern_by_id,          # Find pattern by ID
    list_all_detectors,         # List all detector names
    pattern_info_summary        # Get statistics summary
)

# Example usage
patterns = get_patterns_for_detector('GQuadruplexDetector')
canonical_g4 = patterns['canonical_g4']
print(canonical_g4[0])  # Print first canonical G4 pattern
```

## Regular Expression Syntax Guide

Common regex patterns used in NBDScanner:

- `G{3,}` - 3 or more consecutive G nucleotides
- `[ACGT]{1,7}` - 1 to 7 nucleotides (any combination)
- `(?:CGG){3,}` - Non-capturing group: 3+ CGG repeats
- `[AG]` - Either A or G (purine)
- `[CT]` - Either C or T (pyrimidine)
- `[ACGT]{1,10}?` - Non-greedy matching of 1-10 nucleotides

## Best Practices

1. **Always test** your changes before using in production
2. **Document your changes** - update the description and reference fields
3. **Use unique pattern IDs** - avoid conflicts with existing patterns
4. **Keep base scores consistent** - typically 0.65-0.98 range
5. **Consider performance** - complex regex can be slow on large sequences
6. **Backup original patterns** - keep a copy before making changes

## Algorithmic Detectors

Note that some detectors don't use regex patterns:
- **APhilicDetector** - Uses 10-mer scoring table
- **SlippedDNADetector** - Uses k-mer scanning for STRs and direct repeats
- **CruciformDetector** - Uses k-mer indexing for inverted repeats
- **TriplexDetector** - Uses k-mer scanner for mirror repeats (in addition to sticky DNA patterns)

These detectors have empty or placeholder pattern definitions in `motif_patterns.py`.

## Troubleshooting

### Pattern not detected?

1. Check regex syntax is valid
2. Verify min_length is appropriate
3. Ensure base_score threshold is not too high
4. Test pattern in isolation with a test sequence

### Import errors?

```python
# Test import
python3 -c "from motif_patterns import G4_PATTERNS; print('Success!')"

# If fails, check:
# 1. File is in correct directory
# 2. No syntax errors in motif_patterns.py
# 3. Python path includes the directory
```

### Pattern changes not reflected?

1. Restart Python interpreter (patterns are cached)
2. Check detector's `get_patterns()` method
3. Verify fallback patterns aren't being used

## Support

For questions or issues with pattern customization:
- Check the main README.md for general documentation
- Review detector-specific documentation in detectors.py
- Open an issue on GitHub for bugs or feature requests

---

**Last Updated:** December 2024
**Version:** 2024.1
