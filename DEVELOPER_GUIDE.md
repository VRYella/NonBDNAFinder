# NonBDNAFinder Modular Architecture - Developer Guide

## Quick Start

The NonBDNAFinder codebase is being refactored into a modular architecture for improved maintainability and scalability.

### Current Structure

```
NonBDNAFinder/
├── engine/                 # Core detection engine
│   ├── __init__.py
│   ├── scoring.py         # Score normalization ✓
│   ├── merging.py         # Overlap removal ✓
│   ├── chunking.py        # Sequence chunking ✓
│   ├── sequence_ops.py    # Sequence operations ✓
│   └── detectors/         # Detector classes
│       ├── __init__.py
│       └── base.py        # Base detector ✓
├── utils/                 # Utility functions
│   ├── __init__.py
│   ├── export.py          # Data export ✓
│   ├── constants.py       # Constants ✓
│   ├── fasta.py           # FASTA parsing ✓
│   ├── validation.py      # Validation ✓
│   └── plotting/          # Visualization
│       └── __init__.py
├── ui/                    # UI components
│   └── __init__.py
├── app.py                 # Streamlit app (monolithic)
├── detectors.py           # All detectors (monolithic)
├── nonbscanner.py         # Core API (monolithic)
└── utilities.py           # All utilities (monolithic)
```

### Using the Modules

#### Import Examples

```python
# Score normalization
from engine.scoring import normalize_motif_scores

motifs = [{'Class': 'G-Quadruplex', 'Score': 0.8}]
normalized = normalize_motif_scores(motifs)

# Overlap removal
from engine.merging import remove_overlaps

filtered_motifs = remove_overlaps(motifs)

# Data export
from utils.export import export_to_csv, export_to_excel

csv_content = export_to_csv(motifs)
export_to_excel(motifs, "results.xlsx")

# Sequence operations
from engine.sequence_ops import reverse_complement, gc_content

rc = reverse_complement("ATGC")
gc = gc_content("ATGCATGC")

# FASTA parsing
from utils.fasta import parse_fasta

sequences = parse_fasta(">seq1\nATGC\n>seq2\nGCTA")
```

#### Creating a Custom Detector

```python
from engine.detectors.base import BaseMotifDetector
from typing import Dict, List, Tuple

class MyDetector(BaseMotifDetector):
    def get_motif_class_name(self) -> str:
        return "My_Motif"
    
    def get_patterns(self) -> Dict[str, List[Tuple]]:
        return {
            'my_patterns': [
                (r'ATGC{3,}', 'MY_001', 'Pattern1', 'Subclass1',
                 7, 'simple_score', 0.5, 'Description', 'Reference')
            ]
        }
    
    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        # Implement scoring logic
        return 0.8

# Use the detector
detector = MyDetector()
motifs = detector.detect_motifs(sequence, "my_sequence")
```

### Development Workflow

#### 1. Verify Module Imports

```bash
cd /home/runner/work/NonBDNAFinder/NonBDNAFinder
python3 -c "from engine import scoring; print('✓ Scoring module OK')"
python3 -c "from engine import merging; print('✓ Merging module OK')"
python3 -c "from utils import export; print('✓ Export module OK')"
```

#### 2. Check Module Structure

```bash
# List all created modules
find engine utils ui -name "*.py" -type f | grep -v __pycache__ | sort
```

#### 3. Use Migration Script

```bash
# Preview extraction
python migrate_to_modules.py --dry-run --module validation

# Extract module
python migrate_to_modules.py --module validation

# Extract with backup
python migrate_to_modules.py --backup --module export
```

### Module Guidelines

#### Size Targets

- **Target**: ~200 lines per module
- **Maximum**: ~300 lines before considering split
- **Minimum**: ~50 lines (smaller modules should be combined)

#### Documentation Standards

```python
"""
Module Name
===========

Brief description of module purpose.
Extracted from [source_file] for focused [purpose] following modular architecture.
"""

from typing import List, Dict

def function_name(param: str) -> str:
    """
    Brief description.
    
    Args:
        param: Parameter description
        
    Returns:
        Return value description
        
    Example:
        >>> function_name("test")
        'result'
    """
    pass
```

#### Import Organization

```python
# Standard library
import os
import sys

# Third-party
import numpy as np
import pandas as pd

# Local modules
from engine.scoring import normalize_score
from utils.constants import CORE_OUTPUT_COLUMNS
```

### Testing Modules

#### Manual Testing

```python
# Test scoring module
from engine.scoring import normalize_score_to_1_3

score = normalize_score_to_1_3(0.8, 'G-Quadruplex')
assert 1.0 <= score <= 3.0
print(f"✓ Score normalization works: {score}")

# Test export module
from utils.export import export_to_csv

motifs = [{'Class': 'Test', 'Score': 1.5, 'Start': 1, 'End': 10}]
csv = export_to_csv(motifs)
assert 'Class' in csv
print("✓ CSV export works")
```

#### Integration Testing

```python
# Test module integration
from engine.scoring import normalize_motif_scores
from engine.merging import remove_overlaps
from utils.export import export_to_csv

# Create test motifs
motifs = [
    {'Class': 'A', 'Score': 0.8, 'Start': 1, 'End': 10},
    {'Class': 'A', 'Score': 0.6, 'Start': 5, 'End': 15}
]

# Process through modules
normalized = normalize_motif_scores(motifs)
filtered = remove_overlaps(normalized)
csv = export_to_csv(filtered)

print(f"✓ Pipeline works: {len(motifs)} → {len(filtered)} motifs")
```

### Common Tasks

#### Adding a New Module

1. Create module file in appropriate directory
2. Add module docstring and functions
3. Update `__init__.py` to export module
4. Test imports
5. Commit with descriptive message

#### Extracting Functions

1. Identify functions in monolithic file
2. Check dependencies (imports, constants)
3. Create new module with complete implementations
4. Update `__init__.py`
5. Verify no breaking changes

#### Updating Documentation

1. Update `MODULE_STATUS.md` with progress
2. Update module docstrings
3. Add usage examples
4. Commit documentation changes

### Troubleshooting

#### Import Errors

```python
# If you get "ModuleNotFoundError"
import sys
sys.path.insert(0, '/home/runner/work/NonBDNAFinder/NonBDNAFinder')

# Then try import again
from engine.scoring import normalize_score
```

#### Circular Imports

- Avoid importing from parent modules
- Use lazy imports if needed
- Restructure dependencies

#### Missing Dependencies

```python
# Check what a module needs
import engine.scoring
print(engine.scoring.__file__)

# View module contents
help(engine.scoring)
```

### Status and Progress

See `MODULE_STATUS.md` for detailed progress tracking.

**Current Status**: 29% Complete (9 of 31 modules)

- ✅ Foundation established
- ✅ Core engine modules (scoring, merging, chunking)
- ✅ Utility modules (export, constants, FASTA, validation)
- ⏳ Individual detectors (pending)
- ⏳ Plotting modules (pending)
- ⏳ UI modules (pending)

### Resources

- **Architecture Guide**: `MODULAR_ARCHITECTURE_GUIDE.md`
- **Status Tracking**: `MODULE_STATUS.md`
- **Implementation Summary**: `IMPLEMENTATION_SUMMARY.md`
- **Migration Script**: `migrate_to_modules.py`

### Contributing

When contributing new modules:

1. Follow existing module patterns
2. Include comprehensive docstrings
3. Add type hints
4. Keep modules focused and small
5. Update documentation
6. Test thoroughly

### Questions?

Refer to the documentation files or examine existing modules for patterns and examples.
