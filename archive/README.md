# Archive Directory

This directory contains deprecated and development files that are no longer part of the core application.

## Archived Files

### Python Modules (Superseded)
- `Newdetector.py` - Old A-philic detector (functionality merged into `detectors.py`)
- `optimized_scoring_and_pipeline.py` - Experimental optimization code
- `scanner_agent.py` - Parallel scanner experiment (functionality integrated)
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

The application has been unified into 4 main scripts:

1. **`app.py`** - Streamlit web application
2. **`utilities.py`** - Utility functions, export, and visualization
3. **`nonbscanner.py`** - Main scanner API and analysis orchestration
4. **`detectors.py`** - All 9 detector classes consolidated

See the main README.md for current architecture documentation.
