# Migration Tools

This directory contains automated scripts used for extracting code from monolithic files into modular components.

## Scripts

### extract_detectors.py
Extracts all 9 detector classes from `detectors.py` into individual module files in `engine/detectors/`.

**Usage:**
```bash
python tools/extract_detectors.py
```

**Output:**
- `engine/detectors/curved_dna.py` (559 lines)
- `engine/detectors/z_dna.py` (626 lines)
- `engine/detectors/a_philic.py` (536 lines)
- `engine/detectors/slipped_dna.py` (544 lines)
- `engine/detectors/cruciform.py` (511 lines)
- `engine/detectors/r_loop.py` (516 lines)
- `engine/detectors/triplex.py` (485 lines)
- `engine/detectors/g_quadruplex.py` (421 lines)
- `engine/detectors/i_motif.py` (280 lines)

### extract_utilities.py
Extracts utility functions from `utilities.py` and `app.py` into categorized modules.

**Usage:**
```bash
python tools/extract_utilities.py
```

**Output:**
- `utils/registry.py` - Pattern registry management
- `utils/plotting/distributions.py` - Distribution visualizations
- `utils/plotting/coverage.py` - Coverage maps
- `utils/plotting/density.py` - Density heatmaps
- `utils/plotting/statistical.py` - Statistical plots
- `utils/plotting/genomic.py` - Genomic visualizations
- `ui/formatting.py` - Text formatting helpers
- `ui/downloads.py` - Download functionality

### migrate_to_modules.py
Original migration script for simpler module extractions.

**Usage:**
```bash
# Dry run
python tools/migrate_to_modules.py --dry-run

# Extract specific module
python tools/migrate_to_modules.py --module validation

# Extract with backup
python tools/migrate_to_modules.py --backup --module export
```

## Notes

- These scripts have already been run and the modules extracted
- They are kept for documentation and potential future use
- The extraction process is idempotent - running again will overwrite existing files
- All extracted modules have been tested and verified working

## Module Verification

After extraction, verify imports work:

```python
# Test detector imports
from engine.detectors import CurvedDNADetector, ZDNADetector

# Test utility imports
from utils import registry, export, validation
from utils.plotting import distributions, density

# Test UI imports
from ui import formatting, downloads
```
