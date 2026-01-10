# Quick Start: Modular Architecture

This repository is in the process of transitioning to a modular architecture. The foundation has been established with example implementations.

## Current Structure

```
NonBDNAFinder/
├── app.py                          # Main Streamlit UI (to be refactored)
├── nonbscanner.py                  # Main API (to be refactored)
├── detectors.py                    # All detectors (to be split)
├── utilities.py                    # All utilities (to be split)
├── ui/                             # ✅ UI Components (new)
│   └── __init__.py
├── engine/                         # ✅ Detection Engine (new)
│   ├── __init__.py
│   ├── sequence_ops.py             # ✅ Sequence operations (90 lines)
│   └── detectors/
│       └── __init__.py
└── utils/                          # ✅ Utilities (new)
    ├── __init__.py
    ├── fasta.py                    # ✅ FASTA parsing (120 lines)
    ├── validation.py               # ✅ Validation functions (260 lines)
    └── plotting/
        └── __init__.py
```

## Example Modules Implemented

### 1. engine/sequence_ops.py
```python
from engine.sequence_ops import reverse_complement, validate_sequence, gc_content

# Reverse complement
rc = reverse_complement("ATGC")  # Returns: "GCAT"

# Validate sequence
is_valid, msg = validate_sequence("ATGC")  # Returns: (True, "")

# Calculate GC content
gc = gc_content("ATGC")  # Returns: 50.0
```

### 2. utils/fasta.py
```python
from utils.fasta import parse_fasta, read_fasta_file, format_fasta

# Parse FASTA content
sequences = parse_fasta(">seq1\nATGC\n>seq2\nGCTA")

# Read FASTA file
sequences = read_fasta_file("sequences.fasta")

# Format as FASTA
fasta = format_fasta({'seq1': 'ATGC' * 30}, line_width=60)
```

### 3. utils/validation.py
```python
from utils.validation import (
    validate_sequence,
    validate_motif,
    validate_score,
    validate_coordinates,
    validate_fasta_format
)

# Validate DNA sequence
is_valid, msg = validate_sequence("ATGC")

# Validate motif dictionary
motif = {'Class': 'G-Quadruplex', 'Start': 100, 'End': 120, 'Score': 0.9}
is_valid, msg = validate_motif(motif)

# Validate score range
is_valid, msg = validate_score(2.5, min_score=0.0, max_score=3.0)

# Validate genomic coordinates
is_valid, msg = validate_coordinates(start=100, end=200, sequence_length=1000)
```

## Using the Migration Script

```bash
# See what would be extracted (dry run)
python migrate_to_modules.py --dry-run

# Extract specific module
python migrate_to_modules.py --module validation

# Extract all with backup
python migrate_to_modules.py --backup
```

## Design Principles

1. **Single Responsibility**: Each module does one thing well
2. **Size Control**: Target ~200 lines per module
3. **Type Safety**: Type hints for all public functions
4. **Documentation**: Docstrings with examples
5. **Testing**: Each module has corresponding test file

## Module Template

```python
"""
Module Name
===========

Brief description of module purpose.
Extracted from [source file] for [reason].
"""

from typing import Tuple, Optional


def function_name(param: str) -> Tuple[bool, str]:
    """
    Brief description.
    
    Args:
        param: Description
        
    Returns:
        Tuple of (success: bool, message: str)
        
    Example:
        >>> function_name("test")
        (True, "success")
    """
    # Implementation
    pass
```

## Full Documentation

- **Architecture Guide**: See `MODULAR_ARCHITECTURE_GUIDE.md`
- **Implementation Summary**: See `IMPLEMENTATION_SUMMARY.md`
- **Migration Script**: See `migrate_to_modules.py`

## Next Steps

To complete the modular architecture:

1. **Engine Modules**: Extract detection logic from `nonbscanner.py`
2. **Detector Classes**: Split `detectors.py` into individual files
3. **Utility Modules**: Split `utilities.py` by functionality
4. **UI Components**: Extract components from `app.py`
5. **Integration**: Update all import statements
6. **Testing**: Add comprehensive tests

See `IMPLEMENTATION_SUMMARY.md` for detailed timeline.

## Benefits

- ✅ **Maintainability**: Find code faster
- ✅ **Testability**: Test modules independently
- ✅ **Scalability**: Add features easily
- ✅ **Collaboration**: Multiple developers can work in parallel
- ✅ **Performance**: Lazy loading, better caching

## Current Progress

- ✅ Directory structure created
- ✅ Example modules implemented (3)
- ✅ Migration tooling ready
- ✅ Documentation complete
- ⏳ Full migration remaining (~90%)

**Estimated completion**: 6-7 days full-time work

---

For questions or assistance, refer to the detailed guides in the repository root.
