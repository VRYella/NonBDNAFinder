# Pattern Export and Hyperscan DB Builder Utilities

This directory contains utilities for exporting Non-B DNA patterns and building Hyperscan databases for high-performance pattern matching.

## Utilities

### 1. `export_patterns_full.py` - Pattern Exporter

Exports all regex patterns from `motif_patterns.PATTERN_REGISTRY` to TSV and Excel formats.

#### Features
- Exports to TSV and/or Excel (.xlsx) formats
- Tests each pattern for Python `re` compilation
- Tests each pattern for Hyperscan compatibility
- Excludes R-loop patterns by default (optional)
- Excludes algorithmic/empty patterns by default (optional)
- Includes comprehensive metadata columns

#### Output Columns
- `Class` - Detector class name
- `Subclass` - Motif subclass
- `Pattern_Group` - Pattern group within detector
- `Pattern` - The regex pattern
- `Pattern_ID` - Unique identifier
- `Name` - Human-readable name
- `Min_Length` - Minimum match length
- `Score_Type` - Scoring method
- `Base_Score` - Base score value
- `Python_re_compiles` - Whether pattern compiles in Python re
- `Hyperscan_compatible` - Whether pattern is Hyperscan-compatible
- `Description` - Motif description
- `Reference` - Literature reference

#### Usage Examples

```bash
# Export to both TSV and Excel (default, skip R-loop and empty patterns)
python export_patterns_full.py --out patterns_export

# Include R-loop patterns
python export_patterns_full.py --out patterns_export --include-rloop

# Include algorithmic/empty patterns
python export_patterns_full.py --out patterns_export --include-empty

# Export only TSV format
python export_patterns_full.py --out patterns_export --format tsv

# Export only Excel format
python export_patterns_full.py --out patterns_export --format excel

# Full export with all patterns
python export_patterns_full.py --out patterns_complete --include-rloop --include-empty
```

#### Requirements
- Python 3.7+
- `openpyxl` (for Excel export) - `pip install openpyxl`

---

### 2. `build_hyperscan_db.py` - Hyperscan Database Builder

Builds Hyperscan pattern databases for ultra-fast motif detection.

#### Features
- Loads patterns from TSV/Excel files or directly from registry
- Automatically filters to Hyperscan-compatible patterns
- Builds and serializes Hyperscan databases
- Supports building databases for specific detector classes
- Reports incompatible patterns

#### Usage Examples

```bash
# Build from all patterns in registry (exclude R-loop by default)
python build_hyperscan_db.py --out all_patterns.hsdb

# Build from specific detector class
python build_hyperscan_db.py --class GQuadruplexDetector --out g4.hsdb

# Build from exported TSV file
python build_hyperscan_db.py --patterns patterns.tsv --out patterns.hsdb

# Build from exported Excel file
python build_hyperscan_db.py --patterns patterns.xlsx --out patterns.hsdb

# Include R-loop patterns
python build_hyperscan_db.py --out all_patterns.hsdb --include-rloop
```

#### Requirements
- Python 3.7+
- `hyperscan` - `pip install hyperscan`
- `openpyxl` (for reading Excel files) - `pip install openpyxl`

#### Output
Creates a serialized Hyperscan database file (`.hsdb`) that can be quickly loaded for pattern matching operations.

---

## Complete Workflow Example

```bash
# Step 1: Export all patterns with metadata
python export_patterns_full.py --out nonb_patterns --include-rloop

# This creates:
#   nonb_patterns.tsv   - Tab-separated file for easy viewing/editing
#   nonb_patterns.xlsx  - Excel file with formatting

# Step 2: Build Hyperscan database from patterns
python build_hyperscan_db.py --patterns nonb_patterns.tsv --out nonb_patterns.hsdb

# This creates:
#   nonb_patterns.hsdb  - Compiled Hyperscan database for fast matching

# Step 3: Use the database in your code
# (The utilities.py module already has functions to load .hsdb files)
```

---

## Testing

A comprehensive test suite is available in `test_pattern_utilities.py`:

```bash
# Run all tests
python test_pattern_utilities.py
```

Tests cover:
- Pattern extraction from registry
- Pattern compilation validation
- TSV export functionality
- Excel export functionality
- Hyperscan compatibility checking
- Pattern loading from files
- Filter operations

---

## Performance Notes

### Pattern Export
- Exporting ~20-25 patterns takes a few seconds
- Excel export is slightly slower than TSV
- Pattern compilation tests are fast (< 1ms per pattern)

### Hyperscan Database Building
- Building a database with ~20 patterns takes < 1 second
- Database size is typically 5-20 KB for 20-30 patterns
- Loading a pre-compiled database is ~10-100x faster than compiling from scratch

---

## Integration with NBDScanner

These utilities integrate seamlessly with the main NBDScanner system:

1. **Pattern Registry**: Both utilities read from `motif_patterns.PATTERN_REGISTRY`
2. **Utilities Module**: The `utilities.py` module already supports loading `.hsdb` files
3. **Scanner Backends**: Hyperscan databases can be used with scanner backends for performance

---

## Troubleshooting

### Issue: "openpyxl not available"
**Solution**: Install with `pip install openpyxl`

### Issue: "hyperscan package not available"
**Solution**: Install with `pip install hyperscan`
Note: Hyperscan requires a compatible system (Linux x86_64). For other systems, TSV export still works.

### Issue: "No patterns found"
**Solution**: Make sure you're running from the NonBDNAFinder directory where `motif_patterns.py` exists.

### Issue: Pattern marked as "Hyperscan incompatible"
**Explanation**: Some regex features aren't supported by Hyperscan:
- Backreferences (`\1`, `\2`)
- Lookahead/lookbehind assertions
- Named capture groups

The utilities automatically detect and report these patterns.

---

## License

MIT License - See repository LICENSE file for details.

## Author

Dr. Venkata Rajesh Yella

## Version

2024.1
