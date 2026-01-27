# Motif Classification System

## Overview

NonBDNAFinder uses a **canonical motif taxonomy** to ensure consistent classification across all detectors, visualizations, and exports. This document describes the unified classification system and how to use it.

## Canonical Taxonomy

The authoritative motif classification is defined in `config/motif_taxonomy.py`. This is the **single source of truth** for all class and subclass names.

### 11 Motif Classes

| ID | Class Name | Subclasses |
|----|------------|------------|
| 1  | `Curved_DNA` | Global Curvature, Local Curvature |
| 2  | `Slipped_DNA` | Direct Repeat, STR |
| 3  | `Cruciform` | Cruciform forming IRs |
| 4  | `R-Loop` | R-loop formation sites |
| 5  | `Triplex` | Triplex, Sticky DNA |
| 6  | `G-Quadruplex` | Telomeric G4, Stacked canonical G4s, Stacked G4s with linker, Canonical intramolecular G4, Extended-loop canonical, Higher-order G4 array/G4-wire, Intramolecular G-triplex, Two-tetrad weak PQS |
| 7  | `i-Motif` | Canonical i-motif, Relaxed i-motif, AC-motif |
| 8  | `Z-DNA` | Z-DNA, eGZ |
| 9  | `A-philic_DNA` | A-philic DNA |
| 10 | `Hybrid` | Dynamic overlaps |
| 11 | `Non-B_DNA_Clusters` | Dynamic clusters |

**Total: 11 Classes, 24 Subclasses**

## Architecture

### 1. Canonical Taxonomy (`config/motif_taxonomy.py`)

The master definition of all valid class/subclass names.

```python
from config.motif_taxonomy import (
    MOTIF_CLASSIFICATION,  # Full taxonomy dict
    VALID_CLASSES,         # Set of all valid class names
    VALID_SUBCLASSES,      # Set of all valid subclass names
    CLASS_TO_SUBCLASSES,   # Map: class → subclasses
    SUBCLASS_TO_CLASS,     # Map: subclass → class
)
```

**Key Functions:**
- `get_all_classes()` - Get sorted list of class names
- `get_subclasses_for_class(class_name)` - Get subclasses for a class
- `is_valid_class(class_name)` - Check if class is valid
- `is_valid_subclass(subclass_name)` - Check if subclass is valid
- `is_valid_pairing(class_name, subclass_name)` - Check if pairing is valid

### 2. Normalization Layer (`core/motif_normalizer.py`)

Enforces canonical naming and auto-corrects variants.

```python
from core.motif_normalizer import normalize_class_subclass

# Normalize variant spellings to canonical
canonical_class, canonical_subclass = normalize_class_subclass(
    'g-quadruplex',  # Input (lowercase variant)
    'telomeric g4',  # Input (lowercase variant)
    strict=False,
    auto_correct=True
)
# Returns: ('G-Quadruplex', 'Telomeric G4')
```

**Key Functions:**
- `normalize_class_name(name)` - Normalize class name
- `normalize_subclass_name(name)` - Normalize subclass name
- `normalize_class_subclass(class, subclass)` - Normalize both and validate pairing
- `normalize_motif_dict(motif)` - Normalize Class/Subclass in motif dict
- `validate_motif_dict(motif)` - Validate motif has valid Class/Subclass

**Features:**
- **Case-insensitive**: `'g-quadruplex'` → `'G-Quadruplex'`
- **Alias support**: `'gquadruplex'` → `'G-Quadruplex'`
- **Auto-correction**: If subclass belongs to different class, auto-corrects class
- **Strict mode**: Can raise errors on invalid data

### 3. Export Validation (`export/export_validator.py`)

Validates data before export to CSV/JSON/BED.

```python
from export.export_validator import validate_export_data

# Validate before export
motifs = [...]
validated_motifs = validate_export_data(
    motifs,
    auto_normalize=True,  # Auto-normalize variants
    strict=False          # Warn instead of error
)
```

**Key Functions:**
- `validate_single_motif(motif, index)` - Validate one motif
- `validate_export_data(motifs)` - Validate list of motifs
- `get_export_summary(motifs)` - Get statistics about motifs
- `check_class_completeness(motifs)` - Check which classes are present

## Usage in Detectors

All detectors MUST use normalization before returning results:

```python
from core.motif_normalizer import normalize_class_subclass

class MyDetector(BaseMotifDetector):
    def detect_motifs(self, sequence, sequence_name):
        motifs = []
        
        # ... detection logic ...
        
        # Normalize before adding to results
        canonical_class, canonical_subclass = normalize_class_subclass(
            self.get_motif_class_name(),
            raw_subclass_name,
            strict=False,
            auto_correct=True
        )
        
        motif = {
            'Class': canonical_class,
            'Subclass': canonical_subclass,
            # ... other fields ...
        }
        motifs.append(motif)
        
        return motifs
```

## Usage in Export Functions

Export functions automatically validate data:

```python
from utilities import export_to_csv, export_to_json, export_to_bed

# Export functions automatically validate and normalize
csv_data = export_to_csv(motifs)  # Validates before export
json_data = export_to_json(motifs)  # Validates before export
bed_data = export_to_bed(motifs, sequence_name)  # Validates before export
```

## Benefits

### 1. **Consistency**
- No more duplicate bars in plots due to name variants
- CSV/JSON exports have standardized fields
- Visualizations group correctly by class/subclass

### 2. **Scientific Accuracy**
- Preserves `A-philic_DNA` (not misrepresented as A-DNA)
- Distinguishes `Canonical i-motif` vs `Relaxed i-motif`
- Prevents invalid subclass combinations

### 3. **Maintainability**
- Single source of truth in `config/motif_taxonomy.py`
- Changes propagate automatically to all components
- Easy to add new classes/subclasses

### 4. **Validation**
- Export functions fail fast on invalid data
- Auto-normalization fixes common variants
- Clear error messages guide corrections

## Testing

Run validation tests:

```bash
python validate_taxonomy.py
```

This tests:
- Taxonomy structure (11 classes, 24 subclasses)
- Normalization functions
- Export validation
- Detector integration

## Examples

### Example 1: Validate Motif Data

```python
from core.motif_normalizer import validate_motif_dict

motif = {
    'Class': 'G-Quadruplex',
    'Subclass': 'Telomeric G4',
    'Start': 0,
    'End': 25
}

is_valid = validate_motif_dict(motif, raise_on_error=False)
# Returns: True
```

### Example 2: Auto-Correct Class Mismatch

```python
from core.motif_normalizer import normalize_class_subclass

# Telomeric G4 belongs to G-Quadruplex, not Triplex
canonical_class, canonical_subclass = normalize_class_subclass(
    'Triplex',      # Wrong class
    'Telomeric G4', # Belongs to G-Quadruplex
    auto_correct=True
)
# Returns: ('G-Quadruplex', 'Telomeric G4')
# Automatically corrected class!
```

### Example 3: Get Subclasses for a Class

```python
from config.motif_taxonomy import get_subclasses_for_class

g4_subclasses = get_subclasses_for_class('G-Quadruplex')
# Returns: ['Telomeric G4', 'Stacked canonical G4s', ...]
print(f"G-Quadruplex has {len(g4_subclasses)} subclasses")
# Output: G-Quadruplex has 8 subclasses
```

### Example 4: Export with Validation

```python
from export.export_validator import validate_export_data

motifs = [
    {'Class': 'g-quadruplex', 'Subclass': 'telomeric g4', ...},  # Lowercase
    {'Class': 'Z-DNA', 'Subclass': 'eGZ', ...},                  # Canonical
]

# Auto-normalize before export
validated = validate_export_data(motifs, auto_normalize=True)
# All motifs now have canonical Class/Subclass names
```

## Migration Guide

If you have existing code that hardcodes class/subclass names:

### Before (❌ Don't do this):
```python
motif = {
    'Class': 'G_Quadruplex',        # Non-canonical
    'Subclass': 'canonical_imotif'  # Wrong class!
}
```

### After (✅ Do this):
```python
from core.motif_normalizer import normalize_class_subclass

canonical_class, canonical_subclass = normalize_class_subclass(
    'G-Quadruplex',
    'Canonical intramolecular G4',
    strict=False,
    auto_correct=True
)

motif = {
    'Class': canonical_class,
    'Subclass': canonical_subclass,
}
```

## Contributing

When adding a new detector or modifying an existing one:

1. **Never hardcode** class/subclass strings outside `config/motif_taxonomy.py`
2. **Always use** `normalize_class_subclass()` before returning results
3. **Test** your detector with `validate_taxonomy.py`
4. **Verify** exports pass validation

## Files

| File | Purpose |
|------|---------|
| `config/motif_taxonomy.py` | Canonical taxonomy definition |
| `core/motif_normalizer.py` | Normalization and validation |
| `export/export_validator.py` | Export validation |
| `validate_taxonomy.py` | Test suite |
| `MOTIF_CLASSIFICATION.md` | This documentation |

## References

- [NonBDNAFinder Documentation](README.md)
- [Detector Development Guide](DETECTORS_REFACTORING_SUMMARY.md)
- [API Documentation](nonbscanner.py)
