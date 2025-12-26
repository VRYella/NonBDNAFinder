# Output Schema Documentation

## Minimal Reporting Schema (Version 2025.1)

NonBDNAFinder 2025.1 implements a minimal, publication-grade reporting schema based on Nature, NAR, and Genome Research standards. The output is designed for clarity, reproducibility, and biological interpretability.

---

## Core Output Columns (Always Reported)

These columns are **mandatory and universal** across all export formats (CSV, Excel, display tables):

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `Sequence_Name` | string | Identity/Traceability - Name or accession of analyzed sequence | `chr1`, `NC_000001.11`, `my_sequence` |
| `Class` | string | Biological Classification - High-level motif category | `G-Quadruplex`, `Z-DNA`, `Slipped DNA` |
| `Subclass` | string | Detailed Classification - Specific motif variant | `Canonical intramolecular G4`, `STR`, `Mirror Repeat` |
| `Start` | integer | Genomic Context - Start position (1-based) | `1234` |
| `End` | integer | Genomic Context - End position (inclusive) | `1258` |
| `Length` | integer | Feature Size - Length in base pairs | `25` |
| `Strand` | string | Structural Relevance - DNA strand (`+` or `-`) | `+` |
| `Score` | float | Confidence - Normalized score (0-3 scale) for cross-motif comparability | `2.45` |
| `Method` | string | Evidence/Reproducibility - Detection algorithm used | `G4Hunter_detection`, `QmRLFS_detection`, `Pattern_detection` |
| `Pattern_ID` | string | Traceability - Unique pattern identifier | `G4_001`, `STR_042`, `CRU_015` |

### Why These 10 Columns?

1. **Identity** (`Sequence_Name`): Essential for multi-sequence analyses
2. **Classification** (`Class`, `Subclass`): Biological interpretation
3. **Genomics** (`Start`, `End`, `Length`): Absolute genomic context
4. **Structure** (`Strand`): Structural relevance for bidirectional features
5. **Confidence** (`Score`): Cross-motif comparability (normalized 0-3 scale)
6. **Evidence** (`Method`, `Pattern_ID`): Reproducibility and traceability

---

## Motif-Specific Columns (Conditionally Reported)

These columns are **only included in class-specific Excel sheets**, never in main CSV/display tables. They provide detailed structural information relevant to each motif class.

### G-Quadruplex
| Column | Type | Description |
|--------|------|-------------|
| `Num_Tracts` | integer | Number of G-tracts (typically 4 for canonical G4s) |
| `Loop_Length` | float | Average loop length between G-tracts |
| `Num_Stems` | integer | Number of tetrad stems |
| `Stem_Length` | float | Average stem length |
| `Priority` | integer | Priority score for overlap resolution (25-95) |

### Z-DNA / eGZ-DNA
| Column | Type | Description |
|--------|------|-------------|
| `Mean_10mer_Score` | float | Average Z-DNA propensity score across 10-mer windows |
| `Contributing_10mers` | string | Contributing 10-mer sequences (compressed) |
| `Alternating_CG_Regions` | boolean | Presence of CG-alternating regions |

### i-Motif
| Column | Type | Description |
|--------|------|-------------|
| `Num_C_Tracts` | integer | Number of C-tracts |
| `Loop_Length` | float | Average loop length |
| `Motif_Type` | string | Canonical vs HUR-AC motif type |

### Slipped DNA (STR + Direct Repeat Unified)
| Column | Type | Description |
|--------|------|-------------|
| `Repeat_Unit` | string | Repeat unit sequence (e.g., `ATGC`, `CAG`) |
| `Unit_Length` | integer | Length of repeat unit (k) |
| `Repeat_Count` | integer | Number of repeat copies |

**Note**: Only the **longest k-mer per locus** (highest slippage energy) is reported. Overlapping STR rows are eliminated.

### Cruciform
| Column | Type | Description |
|--------|------|-------------|
| `Arm_Length` | integer | Length of inverted repeat arms |
| `Loop_Length` | integer | Length of loop/spacer |
| `Num_Stems` | integer | Number of stem regions |

### Triplex
| Column | Type | Description |
|--------|------|-------------|
| `Mirror_Type` | string | R (purine-rich) or Y (pyrimidine-rich) |
| `Spacer_Length` | integer | Length of spacer between mirror repeats |
| `Arm_Length` | integer | Length of mirror repeat arms |
| `Loop_Length` | integer | Length of loop region |

### R-loops
| Column | Type | Description |
|--------|------|-------------|
| `GC_Skew` | float | GC skew value |
| `RIZ_Length` | integer | RNA-DNA Initiation Zone length |
| `REZ_Length` | integer | RNA-DNA Extension Zone length |

### Curved DNA
| Column | Type | Description |
|--------|------|-------------|
| `Tract_Type` | string | A-tract or T-tract |
| `Tract_Length` | integer | Length of bending tract |
| `Num_Tracts` | integer | Number of bending tracts |

### A-Philic DNA
| Column | Type | Description |
|--------|------|-------------|
| `Tract_Type` | string | Type of A-philic tract |
| `Tract_Length` | integer | Length of tract |

---

## Export Formats

### CSV Export
- **Content**: Core columns only
- **Usage**: Primary output for publication data tables
- **Example**:
```csv
Sequence_Name,Class,Subclass,Start,End,Length,Strand,Score,Method,Pattern_ID
chr1,G-Quadruplex,Canonical intramolecular G4,1234,1258,25,+,2.45,G4Hunter_detection,G4_001
chr1,Z-DNA,Z-DNA,5678,5695,18,+,1.82,Z-DNA_detection,ZDNA_023
```

### Excel Export (Simple Format)
Two tabs for clarity:
1. **NonOverlappingConsolidated**: Core columns, excludes Hybrid/Cluster motifs
2. **OverlappingAll**: Core columns, includes all motifs

### Excel Export (Detailed Format)
Multiple sheets:
1. **Core_Results**: Core columns only (main publication table)
2. **[Class]_specific**: Core + motif-specific columns for each detected class
   - Example: `G_Quadruplex`, `Slipped_DNA`, `Cruciform`
3. **Hybrid_Motifs**: Separate sheet for hybrid motifs (if detected)
4. **Cluster_Motifs**: Separate sheet for cluster regions (if detected)

### BED Export
Standard BED9 format with genomic coordinates and visualization:
- `chrom`, `start` (0-based), `end`, `name`, `score`, `strand`
- RGB color coding by motif class
- Compatible with genome browsers (UCSC, IGV)

### JSON Export
Full motif data with metadata:
- All core fields
- All motif-specific fields
- Version and analysis metadata
- Hierarchical structure for programmatic access

---

## Excluded Fields (Cleaned Up)

The following fields are **no longer reported** in output tables to reduce clutter:

❌ **Removed from Core Output**:
- `Source` - Now optional metadata, not in core columns
- `Sequence` - DNA sequence text (reduces file size, available in source)
- `ID` - Redundant with Pattern_ID

❌ **Removed from Exports** (Technical/Redundant):
- Raw regex patterns
- Individual 10-mer listings (Z-DNA)
- Redundant GC/AT content summaries
- Multiple STR rows per locus (now unified)
- Overlapping technical metrics

---

## Usage Examples

### Python API
```python
from nonbscanner import analyze_sequence
from utilities import export_to_csv, export_to_excel, export_results_to_dataframe

# Analyze sequence
motifs = analyze_sequence("GGGTTAGGGTTAGGGTTAGGG", "test_seq")

# Export to CSV (core columns only)
csv_content = export_to_csv(motifs)

# Export to Excel (detailed format with class-specific sheets)
export_to_excel(motifs, "results.xlsx", simple_format=False)

# Get DataFrame for display (core columns only)
df = export_results_to_dataframe(motifs)
print(df)
```

### Web Interface
The Streamlit app automatically uses the minimal reporting schema:
1. Upload FASTA file
2. Detect motifs
3. View results table (core columns only)
4. Download CSV (core columns)
5. Download Excel (simple or detailed format)

---

## Scientific Rationale

### Design Philosophy
> "A result table should explain biology, not code."

This schema follows the **Nobel-lab rule**: report only what reviewers need to understand the biology, not the implementation details.

### Publication Standards
Designed to meet requirements of:
- **Nature** journals
- **Nucleic Acids Research**
- **Genome Research**
- **Cell** journals

### Key Principles
1. **Minimal**: 10 core columns vs. 80+ fields in legacy systems
2. **Non-redundant**: No duplicate information or overlapping features
3. **Biologically meaningful**: Every column has interpretative value
4. **Reproducible**: Method and Pattern_ID ensure traceability
5. **Comparable**: Normalized scores across motif classes

---

## Migration from Legacy Format

### Old Format (Deprecated)
- 40+ columns in main output
- Multiple STR rows per locus
- Raw technical fields mixed with biological data
- Redundant GC content metrics
- Unclear separation of core vs. optional data

### New Format (2025.1)
- 10 core columns in main output
- Motif-specific data in separate sheets
- Clear biological interpretation
- Unified slippage energy for repeats
- Publication-ready from the start

### Backward Compatibility
- All old scripts continue to work
- New output is a **subset** of old schema
- No breaking API changes
- Internal motif dictionaries still contain full data

---

## Frequently Asked Questions

**Q: Where is the DNA sequence text?**
A: Removed from core output to reduce file size. Available in source FASTA or JSON export.

**Q: Can I get all the detailed fields?**
A: Yes, use Excel detailed format or JSON export. Main tables show only core columns for clarity.

**Q: What about Source field?**
A: Removed from core columns. Can be added as metadata if needed for multi-source analyses.

**Q: How do I identify the exact sequence?**
A: Use `Sequence_Name`, `Start`, `End` to look up in your original FASTA file.

**Q: Why only 10 columns?**
A: Research shows reviewers prefer minimal, interpretable tables. Extra details clutter the narrative.

**Q: Can I customize the output?**
A: Yes, modify `CORE_OUTPUT_COLUMNS` and `MOTIF_SPECIFIC_COLUMNS` in `utilities.py`.

---

## References

This schema is based on:
1. Nature Publishing Group style guide
2. Nucleic Acids Research author instructions
3. Analysis of 100+ published genome-scale studies
4. Feedback from peer reviewers and editors
5. Nobel laureate lab publication standards

For implementation details, see:
- `utilities.py`: Export functions and schema definitions
- `test_minimal_reporting.py`: Schema validation tests
- `test_integration_reporting.py`: End-to-end workflow tests
