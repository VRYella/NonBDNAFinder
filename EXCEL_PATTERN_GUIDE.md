# Excel Pattern Loading Guide

## Overview

NonBDNAFinder now supports loading pattern data from an Excel file (`pattern_registry2.xlsx`) instead of JSON. This provides a more user-friendly way to view, edit, and manage pattern definitions with updated normalized scores (1-3 range for most patterns).

## Features

✅ **Excel-First Loading**: Automatically loads from `pattern_registry2.xlsx` if available  
✅ **JSON Fallback**: Falls back to `consolidated_registry.json` if Excel not found  
✅ **Performance**: Fast loading with caching  
✅ **Easy Editing**: Edit patterns in Excel with formulas and formatting  
✅ **Complete Coverage**: 170,008 patterns across 9 motif classes with updated scores  
✅ **Normalized Scores**: All patterns (except A-philic and Z-DNA) use 1-3 normalized scoring

## Excel File Structure

### Sheets

The Excel file contains 10 sheets:

1. **Summary** - Overview of all motif classes
   - Class name
   - Number of patterns
   - Generation date
   - Pattern type

2. **APhilic** - A-philic DNA patterns (208 patterns)
   - id, tenmer, score (NOT normalized - range 1.00-2.70)

3. **Cruciform** - Cruciform DNA patterns (9,595 patterns)
   - id, pattern, subclass, score, description, min_arm, max_arm, max_loop, arm_len, loop_len, purity
   - Scores normalized 1.20-3.00

4. **CurvedDNA** - Curved DNA patterns (44 patterns)
   - id, pattern, subclass, score
   - Scores normalized 1.60-3.00

5. **G4** - G-Quadruplex patterns (140,969 patterns)
   - id, pattern, subclass, score
   - Scores normalized 0.60-3.00

6. **IMotif** - i-Motif patterns (454 patterns)
   - id, pattern, subclass, score
   - Scores normalized 1.30-2.90

7. **RLoop** - R-Loop patterns (5 patterns)
   - id, pattern, subclass, score
   - Scores normalized 0.75-0.95

8. **SlippedDNA** - Slipped DNA patterns (127 patterns)
   - id, pattern, subclass, score, description
   - Scores normalized 1.01-3.00

9. **Triplex** - Triplex DNA patterns (18,476 patterns)
   - id, pattern, subclass, score
   - Scores normalized 1.20-3.00

10. **ZDNA** - Z-DNA patterns (130 patterns)
    - id, tenmer, score, pattern, subclass, description, reference
    - Scores NOT normalized - range 0.90-63.00

### Pattern Types

**Tenmer-based patterns** (10-mer sequences):
- APhilic: 208 patterns (scores 1.00-2.70, not normalized)
- ZDNA: 130 patterns (scores 0.90-63.00, not normalized)

**Regex-based patterns with normalized scores (1-3)**:
- CurvedDNA: 44 patterns
- G4: 140,969 patterns
- IMotif: 454 patterns
- RLoop: 5 patterns
- SlippedDNA: 127 patterns
- Triplex: 18,476 patterns
- Cruciform: 9,595 patterns

**Note**: A-philic and Z-DNA patterns explicitly do NOT use 1-3 normalization. They maintain their original scoring schemes for scientific accuracy.

## Usage

### Automatic Loading

The system automatically tries to load `pattern_registry2.xlsx` first, then falls back to JSON:

```python
from utilities import load_db_for_class

# Automatically uses Excel if available
db, id_to_pattern, id_to_score = load_db_for_class('G4', 'registry')
```

### Manual Excel Loading

```python
from utilities import _load_consolidated_registry_from_excel

# Force load from Excel
registry = _load_consolidated_registry_from_excel()
if registry:
    print(f"Loaded {registry['total_patterns']} patterns")
    print(f"Classes: {list(registry['registries'].keys())}")
```

### Editing Patterns

1. Open `pattern_registry2.xlsx` in Excel or LibreOffice
2. Navigate to the sheet for the motif class you want to edit
3. Modify patterns, scores, or add new rows
4. Save the file
5. Restart your application (cache will be cleared)

**Important**: When adding new patterns:
- Ensure `id` column has unique sequential integers
- For regex patterns, use the `pattern` column
- For 10-mer patterns, use the `tenmer` column
- Always provide a `score` value (normalized 1-3 except for A-philic and Z-DNA)

## Performance

### Loading Times

| Method | Cold Load | Cached Load |
|--------|-----------|-------------|
| Excel | 0.035s | 0.000s |
| JSON | 0.0005s | 0.000s |

**Note**: Excel is ~65x slower than JSON for cold loads, but caching ensures this only happens once per session.

### Recommendations

- Use Excel for development and pattern editing
- Keep JSON file as backup for production deployments
- The automatic fallback ensures compatibility

## Testing

Run the comprehensive test suite:

```bash
python test_excel_pattern_loading.py
```

Tests verify:
- ✓ Excel file loading
- ✓ Pattern count accuracy
- ✓ JSON/Excel equivalence
- ✓ Detector integration
- ✓ Performance metrics

## Migration Guide

### From JSON to Excel

The Excel file is already created from the JSON data. To regenerate:

```python
import json
import pandas as pd

# Load JSON
with open('consolidated_registry.json', 'r') as f:
    data = json.load(f)

# Create Excel writer
writer = pd.ExcelWriter('pattern_registry2.xlsx', engine='openpyxl')

# Create summary
summary_data = []
for class_name, registry in data['registries'].items():
    summary_data.append({
        'Class': class_name,
        'Patterns': registry['n_patterns'],
        'Generated_At': registry.get('generated_at', ''),
        'Pattern_Type': registry.get('meta', {}).get('pattern_type', '')
    })

pd.DataFrame(summary_data).to_excel(writer, sheet_name='Summary', index=False)

# Create class sheets
for class_name, registry in data['registries'].items():
    df = pd.DataFrame(registry['patterns'])
    df.to_excel(writer, sheet_name=class_name, index=False)

writer.close()
```

### From Excel to JSON

To export Excel back to JSON format:

```python
from utilities import _load_consolidated_registry_from_excel
import json

# Load from Excel
registry = _load_consolidated_registry_from_excel()

# Save as JSON
with open('consolidated_registry_new.json', 'w') as f:
    json.dump(registry, f, indent=2)
```

## Requirements

### Python Packages

```bash
pip install pandas openpyxl
```

These are optional - the system falls back to JSON if pandas is not available.

### File Structure

```
NonBDNAFinder/
├── pattern_registry2.xlsx         # Primary pattern source (Excel, 170K+ patterns)
├── consolidated_registry.json     # Backup pattern source (JSON)
├── utilities.py                   # Loading logic
└── test_excel_pattern_loading.py  # Test suite
```

## Troubleshooting

### Excel file not loading

1. Check if `pandas` and `openpyxl` are installed:
   ```bash
   pip install pandas openpyxl
   ```

2. Verify file exists:
   ```python
   import os
   print(os.path.exists('pattern_registry2.xlsx'))
   ```

3. Check for file corruption:
   ```python
   import pandas as pd
   pd.read_excel('pattern_registry2.xlsx', sheet_name='Summary')
   ```

### Pattern counts don't match

Run the comparison test:
```bash
python test_excel_pattern_loading.py
```

Look for mismatches in the "COMPARING EXCEL VS JSON" section.

### Slow loading

First load from Excel is slower than JSON (~0.035s vs 0.0005s). This is normal and only happens once per session due to caching.

To improve:
- Use JSON in production
- Keep Excel for development only
- Pre-load patterns at startup

## FAQ

**Q: Do I need to delete the JSON file?**  
A: No. The system uses Excel first, then falls back to JSON if Excel is unavailable.

**Q: Can I use Excel on Windows/Mac/Linux?**  
A: Yes. Any tool that can edit .xlsx files works (Excel, LibreOffice, Google Sheets).

**Q: Will my edits be preserved?**  
A: Yes. Just save the Excel file. Clear the Python cache by restarting your application.

**Q: Can I add new motif classes?**  
A: Yes. Add a new sheet with the same column structure as existing classes.

**Q: What about performance in production?**  
A: Use the JSON file in production. The Excel feature is mainly for development and pattern editing.

## Example: Adding a New Pattern

### G4 Pattern

1. Open `pattern_registry2.xlsx`
2. Go to the `G4` sheet
3. Add a new row at the end:
   ```
   id: <next_id>
   pattern: G{4,}[ACGT]{1,5}G{4,}[ACGT]{1,5}G{4,}[ACGT]{1,5}G{4,}
   subclass: strict_g4
   score: 2.5
   ```
4. Save the file
5. Restart your Python session
6. Test: `load_db_for_class('G4', 'registry')`

**Note**: Use normalized scores (1-3) for G4 patterns.

### ZDNA Pattern (10-mer)

1. Open `pattern_registry2.xlsx`
2. Go to the `ZDNA` sheet
3. Add a new row at the end:
   ```
   id: 130
   tenmer: CGCGCGCGCG
   score: 63.0
   ```
4. Save the file
5. Restart your Python session
6. Test: `load_db_for_class('ZDNA', 'registry')`

**Note**: Z-DNA scores are NOT normalized (can be > 3).

## Summary

The Excel pattern loading system provides:

- ✅ **User-Friendly**: Edit patterns in familiar spreadsheet format
- ✅ **Backward Compatible**: Falls back to JSON automatically
- ✅ **Well-Tested**: Comprehensive test suite included
- ✅ **Production-Ready**: Performance optimized with caching
- ✅ **Flexible**: Supports both regex and 10-mer patterns

For most users, simply using the system as-is will work perfectly. The Excel file is automatically detected and used when available.
