# Archive Directory

This directory contains deprecated and development files that are no longer part of the core application, as well as a reference copy of the parallel scanner.

## Archived Files

### Python Modules (Superseded)
- `Newdetector.py` - Old A-philic detector (functionality merged into `detectors.py`)
- `optimized_scoring_and_pipeline.py` - Experimental optimization code
- `scanner_agent.py` - Reference copy (active version in root directory for parallel scanning)
- `scientific_progress.py` - Progress UI code (integrated into `app.py`)
- `visualization_enhancements.py` - Duplicate visualization functions (in `utilities.py`)

### Test Files
- `test_performance.py` - Performance benchmarking scripts
- `test_slipped_dna.py` - Slipped DNA detector tests
- `validate_improvements.py` - Development validation scripts

### Documentation Files
- `FINAL_SUMMARY.md` - Development summary
- `IMPLEMENTATION_SUMMARY.md` - Implementation notes
- `IMPROVEMENTS_SUMMARY.md` - Improvement tracking
- `SLIPPED_DNA_REFINEMENTS.md` - Slipped DNA refinement notes
- `VALIDATION_REPORT.md` - Validation results

## Core Application Files

The application has been unified into 4 main scripts with 1 optional parallel scanner:

1. **`app.py`** - Streamlit web application
2. **`utilities.py`** - Utility functions, export, and visualization
3. **`nonbscanner.py`** - Main scanner API and analysis orchestration
4. **`detectors.py`** - All 9 detector classes consolidated
5. **`scanner_agent.py`** - Optional parallel scanner for large sequences (>100kb)

**Note:** `scanner_agent.py` has been restored from archive to root directory to enable parallel scanning functionality for sequences larger than 100kb. The archive contains a reference copy.

See the main README.md for current architecture documentation.
