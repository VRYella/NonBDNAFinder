# Archives

This directory contains files that are not required for the current architecture of NBDScanner but are preserved for historical reference.

## Contents

### `/documentation`
Contains various documentation files that were created during development:
- Implementation summaries and guides
- Performance optimization documents
- UI improvement notes
- Refactoring recommendations
- Visualization guides

These documents provide valuable historical context but are not needed for running the tool.

### `/tests`
Contains test files:
- `test_detectors.py` - Detector module tests
- `test_imports.py` - Import verification tests
- `test_new_visualizations.py` - Visualization tests

### `/benchmarks`
Contains performance benchmark scripts and results:
- `performance_benchmark.py` - Benchmark suite
- `benchmark_results.json` - Historical benchmark data

### Root Files
- `code_quality_improvements.py` - Code quality analysis scripts

## Why These Files Are Archived

These files have been moved to archives to:
1. Keep the root directory clean and focused on core functionality
2. Reduce cognitive load when navigating the codebase
3. Separate historical documentation from active documentation (README.md)
4. Make the current architecture more apparent

## Accessing Archived Files

All files remain in git history and can be accessed:
- In this `archives/` directory
- Through git history: `git log --follow <filename>`
- On GitHub in the commit history

## Active Documentation

For current documentation, please see:
- `/README.md` - Main project documentation
- Code comments in source files
- Docstrings in Python modules
