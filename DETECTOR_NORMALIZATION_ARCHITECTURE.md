# Detector Self-Normalization Architecture

## Overview

NonBDNAFinder 2024.2+ uses a **self-normalization architecture** where each detector class embeds its own normalization logic, eliminating the need for centralized score transformation. This design improves encapsulation, makes detectors self-contained, and enables per-detector tuning of normalization parameters.

## Architecture Design

### Before (Pre-2024.2)
```
[Detector] → Raw Score → [utilities.normalize_motif_scores()] → Normalized Score
```

**Issues:**
- Separation of concerns: scoring in detectors, normalization in utilities
- Scattered logic: users must look in two places to understand final scores
- Less tunable: normalization parameters buried in SCORE_NORMALIZATION_PARAMS dictionary

### After (2024.2+)
```
[Detector] → Raw Score → [Detector._normalize_score()] → Normalized Score
```

**Benefits:**
- ✅ Self-contained detectors (scoring + normalization in one place)
- ✅ Tunable parameters at top of each detector file
- ✅ Clear separation: Raw_Score (detector-specific) and Score (universal 1-3)
- ✅ Easier per-detector customization

## Universal Score Scale

All detectors normalize scores to a **universal 1-3 scale** for cross-motif comparability:

### Score Interpretation Bands

| Score Range | Interpretation   | Biological Meaning                  |
|-------------|------------------|-------------------------------------|
| 1.0 - 1.7   | Weak/Conditional | Low confidence, context-dependent   |
| 1.7 - 2.3   | Moderate         | Reasonable confidence, likely valid |
| 2.3 - 3.0   | Strong/High      | High confidence, well-characterized |

### Cross-Motif Comparability

- A G-Quadruplex score of 2.5 is **directly comparable** to a Curved DNA score of 2.5
- Higher scores **always** indicate stronger/more confident predictions
- Scores are normalized **after** detection (never influence candidate finding)

## Detector Implementation

### Normalization Parameters

Each detector defines normalization parameters as class constants:

```python
# ═══════════════════════════════════════════════════════════════════════════════
# NORMALIZATION PARAMETERS (Tunable)
# ═══════════════════════════════════════════════════════════════════════════════
# ┌──────────────┬─────────────┬────────────────────────────────────────┐
# │ Parameter    │ Value       │ Scientific Basis                       │
# ├──────────────┼─────────────┼────────────────────────────────────────┤
# │ RAW_MIN      │ 0.5         │ Bedrat 2016 - minimal G4 stability     │
# │ RAW_MAX      │ 1.0         │ Bedrat 2016 - maximal G4Hunter score   │
# │ NORM_MIN     │ 1.0         │ Universal low confidence threshold     │
# │ NORM_MAX     │ 3.0         │ Universal high confidence threshold    │
# │ METHOD       │ 'g4hunter'  │ G4Hunter-specific normalization        │
# └──────────────┴─────────────┴────────────────────────────────────────┘
G4_RAW_SCORE_MIN = 0.5; G4_RAW_SCORE_MAX = 1.0
G4_NORMALIZED_MIN = 1.0; G4_NORMALIZED_MAX = 3.0
G4_NORMALIZATION_METHOD = 'g4hunter'
G4_SCORE_REFERENCE = 'Bedrat et al. 2016 (Bioinformatics)'
# ═══════════════════════════════════════════════════════════════════════════════
```

### Class Implementation

```python
class GQuadruplexDetector(BaseMotifDetector):
    """Ultra-fast seeded G4 detector with self-normalization."""
    
    # Override normalization parameters
    RAW_SCORE_MIN = G4_RAW_SCORE_MIN
    RAW_SCORE_MAX = G4_RAW_SCORE_MAX
    NORMALIZED_MIN = G4_NORMALIZED_MIN
    NORMALIZED_MAX = G4_NORMALIZED_MAX
    NORMALIZATION_METHOD = G4_NORMALIZATION_METHOD
    SCORE_REFERENCE = G4_SCORE_REFERENCE
    
    def detect_motifs(self, sequence: str, sequence_name: str) -> List[Dict]:
        # ... detection logic ...
        
        # Self-normalize scores
        raw_score = calculated_score
        normalized_score = self._normalize_score(raw_score)
        
        motif = {
            'Raw_Score': round(raw_score, 3),      # Detector-specific scale
            'Score': normalized_score,              # Universal 1-3 scale
            # ... other fields ...
        }
        return motifs
```

## Normalization Methods

The base `_normalize_score()` method supports multiple normalization strategies:

### 1. Linear Normalization (Default)

Used by most detectors (Curved DNA, Slipped DNA, Cruciform, Triplex, i-Motif, R-Loop, A-philic):

```python
normalized = 1.0 + (raw - RAW_MIN) * (3.0 - 1.0) / (RAW_MAX - RAW_MIN)
```

**Characteristics:**
- Simple linear interpolation
- Preserves relative ordering
- Raw scores typically in [0.0, 1.0] range

### 2. G4Hunter Method

Used by G-Quadruplex detector for G4Hunter-specific scoring:

```python
effective_score = abs(raw_score)  # Use absolute value
if effective_score < 0.5:
    return 1.0
elif effective_score >= 1.0:
    return min(3.0, 2.0 + (effective_score - 0.5) * 2.0)
else:
    return 1.0 + (effective_score - 0.5) * 2.0
```

**Characteristics:**
- Handles both positive and negative G4Hunter scores
- Non-linear scaling for G4 formation potential
- Based on Bedrat et al. 2016 G4Hunter algorithm

### 3. Log-Scale Normalization

Used by Z-DNA detector for cumulative 10-mer scores:

```python
log_raw = math.log10(max(1, raw_score))
log_min = math.log10(max(1, RAW_MIN))
log_max = math.log10(RAW_MAX)
normalized = 1.0 + (log_raw - log_min) * 2.0 / (log_max - log_min)
```

**Characteristics:**
- Log-linear scaling for cumulative scores
- Handles wide range of raw scores (50-2000)
- Better distribution for Z-DNA 10-mer accumulation

## Detector-Specific Parameters

| Detector     | RAW_MIN | RAW_MAX | Method           | Reference                       |
|--------------|---------|---------|------------------|---------------------------------|
| G-Quadruplex | 0.5     | 1.0     | g4hunter         | Bedrat et al. 2016              |
| Z-DNA        | 50.0    | 2000.0  | log              | Ho et al. 1986                  |
| Curved DNA   | 0.1     | 0.95    | linear           | Koo et al. 1986                 |
| Slipped DNA  | 0.3     | 0.98    | linear           | Schlötterer et al. 2000         |
| Cruciform    | 0.5     | 0.95    | linear           | Lilley 2000, Sinden 1994        |
| Triplex      | 0.5     | 0.95    | linear           | Frank-Kamenetskii 1995          |
| i-Motif      | 0.4     | 0.98    | linear           | Gehring et al. 1993             |
| R-Loop       | 0.4     | 0.95    | linear           | Aguilera 2012                   |
| A-philic DNA | 0.5     | 0.95    | linear           | Gorin et al. 1995               |

## Customizing Normalization

To tune normalization for a specific detector, modify the class constants at the top of the detector file:

```python
# Example: Make G-Quadruplex detector more stringent
G4_RAW_SCORE_MIN = 0.7  # Increase minimum threshold
G4_RAW_SCORE_MAX = 1.2  # Extend maximum range
```

Changes take effect immediately without modifying any other code.

## Output Format

All detector methods now return motifs with both scores:

```python
motif = {
    'ID': 'seq_1_123',
    'Class': 'G-Quadruplex',
    'Subclass': 'Canonical intramolecular G4',
    'Start': 123,
    'End': 145,
    'Sequence': 'GGGATGGGCTGGGAAGGG',
    'Raw_Score': 0.85,      # Detector-specific scale (0-1 for G4)
    'Score': 2.5,           # Universal 1-3 scale
    'Strand': '+',
    # ... other fields ...
}
```

## Backward Compatibility

The deprecated centralized functions remain available with warnings:

```python
from Utilities.utilities import normalize_motif_scores

# DeprecationWarning: normalize_motif_scores() is deprecated.
# Detectors now self-normalize scores.
motifs = normalize_motif_scores(motifs)  # Returns motifs unchanged
```

## Testing

Comprehensive unit tests validate normalization:

```bash
# Run normalization tests
python -m unittest tests.test_normalization -v

# Verify all tests pass
# - 10/10 normalization tests pass
# - 74/74 existing detector tests pass
```

## Migration Guide

### For Users

No changes required! The new architecture is transparent to end users:

```python
from Utilities.nonbscanner import NonBScanner

scanner = NonBScanner()
results = scanner.analyze_sequence(sequence, 'test')

# Results automatically include normalized scores
for motif in results:
    print(f"Score: {motif['Score']}")  # Universal 1-3 scale
    print(f"Raw: {motif['Raw_Score']}")  # Detector-specific
```

### For Developers

If you've created custom detectors, update them to use self-normalization:

1. Add normalization parameters to your detector file
2. Override class constants in your detector class
3. Update `detect_motifs()` to call `self._normalize_score()`
4. Return both `Raw_Score` and `Score` in motif dictionaries

See any of the 9 built-in detectors for reference implementations.

## Scientific References

- **G-Quadruplex**: Bedrat et al. 2016 (Bioinformatics) - G4Hunter scoring
- **Z-DNA**: Ho et al. 1986 (EMBO J) - Z-DNA propensity scores
- **Curved DNA**: Koo et al. 1986 (Biochemistry) - DNA bending models
- **Slipped DNA**: Schlötterer & Tautz 2000 (Nat Rev Genet) - Microsatellite instability
- **Cruciform**: Lilley 1980 (PNAS) - Inverted repeat stability
- **Triplex**: Frank-Kamenetskii & Mirkin 1995 (Annu Rev Biochem) - Triple helix formation
- **i-Motif**: Gehring et al. 1993, Zeraati et al. 2018 - i-Motif structures
- **R-Loop**: Aguilera & García-Muse 2012 - R-loop formation
- **A-philic DNA**: Gorin et al. 1995, Vinogradov 2003 - A-philic sequences
