# NonBScanner - Complete Documentation

## 📋 Table of Contents

1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Analysis Pipeline](#analysis-pipeline)
4. [API Reference](#api-reference)
5. [Detection Coverage](#detection-coverage)
6. [Parameters Reference](#parameters-reference)
7. [Output Format](#output-format)
8. [Performance](#performance)
9. [Usage Examples](#usage-examples)

---

## Overview

NonBScanner is a state-of-the-art bioinformatics tool for detecting and analyzing Non-B DNA motifs in genomic sequences. It provides comprehensive detection of **11 major motif classes** and **22+ specialized subclasses**.

### Key Features

| Feature | Description |
|---------|-------------|
| 🧬 **11 Motif Classes** | Comprehensive Non-B DNA structure detection |
| ⚡ **High Performance** | 24,674 bp/second processing speed |
| 🔬 **Scientific Scoring** | Literature-validated algorithms |
| 📊 **Rich Visualizations** | 25+ publication-quality plots |
| 📁 **Multiple Exports** | CSV, BED, JSON, Excel formats |
| 🖥️ **Multiple Interfaces** | Python API, Web UI, CLI |

---

## Architecture

### File Structure Diagram

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                          NonBScanner Architecture                            │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  ┌─────────────────────────────────────────────────────────────────────────┐ │
│  │                         User Interfaces                                 │ │
│  │  ┌─────────────┐   ┌─────────────────┐   ┌────────────────┐            │ │
│  │  │  Python API │   │ Web Interface   │   │  CLI (Headless)│            │ │
│  │  │ (import nbs)│   │ (app.py)        │   │ (cli/*.py)     │            │ │
│  │  └──────┬──────┘   └────────┬────────┘   └───────┬────────┘            │ │
│  └─────────┼──────────────────┬┼────────────────────┼────────────────────┬┘ │
│            │                  ││                    │                    │   │
│            ▼                  ▼▼                    ▼                    │   │
│  ┌─────────────────────────────────────────────────────────────────────────┐ │
│  │                      Core Analysis Layer                                │ │
│  │  ┌─────────────────────────────────────────────────────────────────────┐│ │
│  │  │                    nonbscanner.py (~1000 lines)                     ││ │
│  │  │  ┌─────────────────┐  ┌────────────────┐  ┌──────────────────┐      ││ │
│  │  │  │ NonBScanner     │  │ analyze_*()    │  │ Progress         │      ││ │
│  │  │  │ (Main Class)    │  │ (Public API)   │  │ Tracking         │      ││ │
│  │  │  └────────┬────────┘  └───────┬────────┘  └─────────┬────────┘      ││ │
│  │  └───────────┼──────────────────┬┼────────────────────┬┼───────────────┘│ │
│  │              │                  ││                    ││                │ │
│  │              ▼                  ▼▼                    ▼▼                │ │
│  │  ┌─────────────────────────────────────────────────────────────────────┐│ │
│  │  │              detectors.py (~3500 lines)                             ││ │
│  │  │  ┌──────────────────────────────────────────────────────────────────┤│ │
│  │  │  │  9 Specialized Detector Classes                                  ││ │
│  │  │  │  ┌─────────────┐ ┌─────────────┐ ┌─────────────┐                 ││ │
│  │  │  │  │ CurvedDNA   │ │ SlippedDNA  │ │ Cruciform   │                 ││ │
│  │  │  │  └─────────────┘ └─────────────┘ └─────────────┘                 ││ │
│  │  │  │  ┌─────────────┐ ┌─────────────┐ ┌─────────────┐                 ││ │
│  │  │  │  │ R-Loop      │ │ Triplex     │ │ G-Quadruplex│                 ││ │
│  │  │  │  └─────────────┘ └─────────────┘ └─────────────┘                 ││ │
│  │  │  │  ┌─────────────┐ ┌─────────────┐ ┌─────────────┐                 ││ │
│  │  │  │  │ i-Motif     │ │ Z-DNA       │ │ A-philic    │                 ││ │
│  │  │  │  └─────────────┘ └─────────────┘ └─────────────┘                 ││ │
│  │  │  └──────────────────────────────────────────────────────────────────┘│ │
│  │  └─────────────────────────────────────────────────────────────────────┘│ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                    │                                         │
│                                    ▼                                         │
│  ┌─────────────────────────────────────────────────────────────────────────┐ │
│  │                        Support Modules                                  │ │
│  │  ┌────────────────┐  ┌────────────────┐  ┌────────────────────────────┐│ │
│  │  │ utilities.py   │  │visualizations.py│ │ scanner.py                 ││ │
│  │  │ (~2100 lines)  │  │ (~1000 lines)  │  │ (k-mer indexing)           ││ │
│  │  │ • I/O          │  │ • 25+ plots    │  │ • Direct repeats           ││ │
│  │  │ • Export       │  │ • Statistics   │  │ • Inverted repeats         ││ │
│  │  │ • Statistics   │  │ • Scientific   │  │ • Mirror repeats           ││ │
│  │  └────────────────┘  └────────────────┘  └────────────────────────────┘│ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                    │                                         │
│                                    ▼                                         │
│  ┌─────────────────────────────────────────────────────────────────────────┐ │
│  │                     Pattern Registry                                    │ │
│  │  ┌─────────────────────────────────────────────────────────────────────┐│ │
│  │  │           consolidated_registry.json (411 patterns)                 ││ │
│  │  │  • A-philic patterns    • Z-DNA patterns                            ││ │
│  │  │  • Curvature patterns   • G-Quadruplex patterns                     ││ │
│  │  └─────────────────────────────────────────────────────────────────────┘│ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

## Analysis Pipeline

### High-Level Analysis Flow

```
┌────────────────────────────────────────────────────────────────────────────┐
│                         ANALYSIS PIPELINE FLOWCHART                        │
└────────────────────────────────────────────────────────────────────────────┘

                              ┌──────────────┐
                              │   INPUT      │
                              │  Sequence    │
                              │  (FASTA)     │
                              └──────┬───────┘
                                     │
                                     ▼
                    ┌────────────────────────────────┐
                    │      VALIDATION STAGE          │
                    │  • Check ATGC characters       │
                    │  • Remove invalid characters   │
                    │  • Determine processing mode   │
                    └────────────────┬───────────────┘
                                     │
                          ┌──────────┴──────────┐
                          │   Sequence Size?    │
                          └──────────┬──────────┘
                                     │
              ┌──────────────────────┼──────────────────────┐
              │                      │                      │
         <10kb bp                10kb-100kb            >100kb bp
              │                      │                      │
              ▼                      ▼                      ▼
    ┌─────────────────┐   ┌──────────────────┐   ┌─────────────────┐
    │ Standard Mode   │   │  Chunked Mode    │   │ Parallel Mode   │
    │ (Sequential)    │   │  (10kb chunks)   │   │ (Multi-threaded)│
    └────────┬────────┘   └────────┬─────────┘   └────────┬────────┘
              │                      │                      │
              └──────────────────────┼──────────────────────┘
                                     │
                                     ▼
                    ┌────────────────────────────────┐
                    │      DETECTION STAGE           │
                    │  Run 9 specialized detectors   │
                    └────────────────┬───────────────┘
                                     │
     ┌───────────────────────────────┼───────────────────────────────┐
     │           │           │       │       │       │       │       │
     ▼           ▼           ▼       ▼       ▼       ▼       ▼       ▼
┌─────────┐┌─────────┐┌─────────┐┌─────────┐┌─────────┐┌─────────┐┌─────────┐
│ Curved  ││ Slipped ││Cruciform││ R-Loop  ││ Triplex ││   G4    ││ i-Motif │
│  DNA    ││   DNA   ││         ││         ││         ││         ││         │
└────┬────┘└────┬────┘└────┬────┘└────┬────┘└────┬────┘└────┬────┘└────┬────┘
     │          │          │          │          │          │          │
     └──────────┴──────────┴──────────┼──────────┴──────────┴──────────┘
                                      │
                              ┌───────┴───────┐
                              │  + Z-DNA      │
                              │  + A-philic   │
                              └───────┬───────┘
                                      │
                                      ▼
                    ┌────────────────────────────────┐
                    │     OVERLAP RESOLUTION         │
                    │  • Remove duplicate motifs     │
                    │  • Merge overlapping regions   │
                    │  • Score-based prioritization  │
                    └────────────────┬───────────────┘
                                     │
                                     ▼
                    ┌────────────────────────────────┐
                    │    HYBRID & CLUSTER DETECTION  │
                    │  • Detect overlapping classes  │
                    │  • Find high-density regions   │
                    │  • Calculate combined scores   │
                    └────────────────┬───────────────┘
                                     │
                                     ▼
                    ┌────────────────────────────────┐
                    │      SCORE NORMALIZATION       │
                    │  • Normalize to 1-3 scale      │
                    │  • Class-specific calibration  │
                    └────────────────┬───────────────┘
                                     │
                                     ▼
                    ┌────────────────────────────────┐
                    │         OUTPUT STAGE           │
                    │  • Sort by position            │
                    │  • Generate statistics         │
                    │  • Export to file formats      │
                    └────────────────────────────────┘
```

### Sequence of Operations (Detailed)

```
┌───────────────────────────────────────────────────────────────────────────────┐
│                    DETAILED SEQUENCE OF OPERATIONS                            │
├───────────────────────────────────────────────────────────────────────────────┤
│                                                                               │
│  Step │ Operation              │ Complexity │ Description                    │
│ ──────┼────────────────────────┼────────────┼────────────────────────────────│
│   1   │ Input Parsing          │ O(n)       │ Parse FASTA, extract sequences │
│   2   │ Sequence Validation    │ O(n)       │ Validate ATGC characters       │
│   3   │ Chunking Decision      │ O(1)       │ Determine processing strategy  │
│   4   │ Detector Initialization│ O(1)       │ Create 9 detector instances    │
│ ──────┼────────────────────────┼────────────┼────────────────────────────────│
│   5a  │ CurvedDNA Detection    │ O(n)       │ A-tract phasing analysis       │
│   5b  │ SlippedDNA Detection   │ O(n)       │ K-mer indexed repeat search    │
│   5c  │ Cruciform Detection    │ O(n)       │ Inverted repeat search         │
│   5d  │ R-Loop Detection       │ O(n)       │ QmRLFS algorithm               │
│   5e  │ Triplex Detection      │ O(n)       │ Mirror repeat & purine runs    │
│   5f  │ G-Quadruplex Detection │ O(n)       │ G4Hunter + pattern matching    │
│   5g  │ i-Motif Detection      │ O(n)       │ C-run pattern detection        │
│   5h  │ Z-DNA Detection        │ O(n)       │ CG/CA repeat scoring           │
│   5i  │ A-philic Detection     │ O(n)       │ 10-mer tetranucleotide search  │
│ ──────┼────────────────────────┼────────────┼────────────────────────────────│
│   6   │ Result Aggregation     │ O(m)       │ Collect all motif hits         │
│   7   │ Overlap Removal        │ O(m log m) │ Remove redundant motifs        │
│   8   │ Hybrid Detection       │ O(m²)      │ Find overlapping class pairs   │
│   9   │ Cluster Detection      │ O(m)       │ Sliding window density scan    │
│  10   │ Score Normalization    │ O(m)       │ Scale scores to 1-3 range      │
│  11   │ Position Sorting       │ O(m log m) │ Sort by genomic position       │
│  12   │ Output Generation      │ O(m)       │ Format results for export      │
│                                                                               │
│  Legend: n = sequence length (bp), m = number of motifs detected             │
└───────────────────────────────────────────────────────────────────────────────┘
```

### Progress Tracking Visualization

```
┌───────────────────────────────────────────────────────────────────────────────┐
│                         PROGRESS TRACKING DISPLAY                             │
├───────────────────────────────────────────────────────────────────────────────┤
│                                                                               │
│  Analysis Progress                                                            │
│  ════════════════                                                             │
│                                                                               │
│  [■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■░░░░░░░░░░░░░░░░░░░░] 60%                     │
│                                                                               │
│  ┌─────────────────────────────────────────────────────────────────────────┐ │
│  │ Detector          │ Status    │ Time    │ Motifs │ Speed (bp/s)        │ │
│  ├───────────────────┼───────────┼─────────┼────────┼─────────────────────┤ │
│  │ ✓ Curved DNA      │ Complete  │ 0.124s  │ 12     │ 8,065               │ │
│  │ ✓ Slipped DNA     │ Complete  │ 0.089s  │ 8      │ 11,236              │ │
│  │ ✓ Cruciform       │ Complete  │ 0.156s  │ 5      │ 6,410               │ │
│  │ ✓ R-Loop          │ Complete  │ 0.201s  │ 15     │ 4,975               │ │
│  │ ✓ Triplex         │ Complete  │ 0.134s  │ 7      │ 7,463               │ │
│  │ ► G-Quadruplex    │ Running   │ 0.089s  │ ...    │ ...                 │ │
│  │ ○ i-Motif         │ Pending   │ -       │ -      │ -                   │ │
│  │ ○ Z-DNA           │ Pending   │ -       │ -      │ -                   │ │
│  │ ○ A-philic DNA    │ Pending   │ -       │ -      │ -                   │ │
│  └───────────────────┴───────────┴─────────┴────────┴─────────────────────┘ │
│                                                                               │
│  Summary: 6/9 detectors complete | 47 motifs found | 1.2 seconds elapsed     │
│                                                                               │
└───────────────────────────────────────────────────────────────────────────────┘
```

---

## API Reference

### Main Functions

#### `analyze_sequence()`

Analyze a single DNA sequence for all Non-B DNA motifs.

```python
def analyze_sequence(
    sequence: str,
    sequence_name: str = "sequence",
    use_fast_mode: bool = False,
    use_chunking: bool = None,
    chunk_size: int = 10000,
    chunk_overlap: int = 500,
    progress_callback: Optional[Callable] = None,
    use_parallel_chunks: bool = True
) -> List[Dict[str, Any]]
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sequence` | str | (required) | DNA sequence (ATGC characters) |
| `sequence_name` | str | "sequence" | Identifier for the sequence |
| `use_fast_mode` | bool | False | Enable parallel processing (~9x speedup) |
| `use_chunking` | bool | None | Enable chunking (auto-enabled for >10kb) |
| `chunk_size` | int | 10000 | Size of each chunk in bp |
| `chunk_overlap` | int | 500 | Overlap between chunks |
| `progress_callback` | Callable | None | Callback for progress tracking |
| `use_parallel_chunks` | bool | True | Enable parallel chunk processing |

**Returns:** `List[Dict[str, Any]]` - List of motif dictionaries sorted by position.

**Example:**
```python
import nonbscanner as nbs

# Basic usage
motifs = nbs.analyze_sequence("GGGTTAGGGTTAGGGTTAGGG", "my_sequence")

# With progress tracking
def progress(chunk, total, bp, elapsed, throughput):
    print(f"Chunk {chunk}/{total}: {bp:,} bp @ {throughput:,.0f} bp/s")

motifs = nbs.analyze_sequence(large_seq, "genome", 
                               use_chunking=True, 
                               progress_callback=progress)
```

---

#### `analyze_fasta()`

Analyze multiple sequences from FASTA format content.

```python
def analyze_fasta(fasta_content: str) -> Dict[str, List[Dict[str, Any]]]
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `fasta_content` | str | FASTA format string |

**Returns:** `Dict[str, List[Dict]]` - Dictionary mapping sequence names to motif lists.

---

#### `analyze_file()`

Analyze sequences from a FASTA file.

```python
def analyze_file(filename: str) -> Dict[str, List[Dict[str, Any]]]
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `filename` | str | Path to FASTA file |

**Returns:** `Dict[str, List[Dict]]` - Dictionary mapping sequence names to motif lists.

---

#### `export_results()`

Export motifs to various formats.

```python
def export_results(
    motifs: List[Dict[str, Any]],
    format: str = 'csv',
    filename: Optional[str] = None,
    **kwargs
) -> str
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `motifs` | List[Dict] | (required) | List of motif dictionaries |
| `format` | str | 'csv' | Export format: 'csv', 'bed', 'json', 'excel' |
| `filename` | str | None | Optional output filename |

---

#### `get_motif_info()`

Get comprehensive information about motif classification system.

```python
def get_motif_info() -> Dict[str, Any]
```

**Returns:** Dictionary with version, author, and classification details.

---

### Progress Callback Signature

The `progress_callback` parameter accepts a function with this signature:

```python
def callback(
    chunk_num: int,      # Current chunk number
    total_chunks: int,   # Total number of chunks
    bp_processed: int,   # Base pairs processed so far
    elapsed: float,      # Elapsed time in seconds
    throughput: float    # Processing speed in bp/s
) -> None
```

---

## Detection Coverage

### Motif Class Summary

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                         MOTIF DETECTION COVERAGE                             │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Class #  │ Name           │ Subclasses │ Detection Method    │ Speed      │
│ ──────────┼────────────────┼────────────┼─────────────────────┼────────────│
│     1     │ Curved DNA     │ 2          │ A-tract phasing     │ ~8,000 bp/s│
│     2     │ Slipped DNA    │ 2          │ K-mer indexing      │ ~280k bp/s │
│     3     │ Cruciform      │ 1          │ Inverted repeat     │ ~8,000 bp/s│
│     4     │ R-Loop         │ 3          │ QmRLFS algorithm    │ ~6,000 bp/s│
│     5     │ Triplex        │ 2          │ Mirror + purine     │ ~7,000 bp/s│
│     6     │ G-Quadruplex   │ 7          │ G4Hunter + patterns │ ~5,000 bp/s│
│     7     │ i-Motif        │ 3          │ C-run patterns      │ ~6,000 bp/s│
│     8     │ Z-DNA          │ 2          │ CG/CA scoring       │ ~7,000 bp/s│
│     9     │ A-philic DNA   │ 1          │ 10-mer patterns     │ ~7,000 bp/s│
│    10     │ Hybrid         │ Dynamic    │ Overlap detection   │ O(m²)      │
│    11     │ Clusters       │ Dynamic    │ Window density      │ O(m)       │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### Subclass Details

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                         DETAILED SUBCLASS BREAKDOWN                          │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  CURVED DNA                                                                  │
│  ├─ Global Curvature     : A-tract induced overall bending                  │
│  └─ Local Curvature      : Regional DNA flexibility                         │
│                                                                              │
│  SLIPPED DNA                                                                 │
│  ├─ Direct Repeat (DR)   : Tandem repeats (10-50 bp units)                  │
│  └─ Short Tandem Repeat  : STRs (1-9 bp units, ≥10 bp total)               │
│                                                                              │
│  CRUCIFORM                                                                   │
│  └─ Inverted Repeat      : Palindromic sequences (10-100 bp arms)           │
│                                                                              │
│  R-LOOP                                                                      │
│  ├─ R-loop Formation     : RNA-DNA hybrid regions                           │
│  ├─ QmRLFS-m1            : Quantitative R-loop model 1                      │
│  └─ QmRLFS-m2            : Quantitative R-loop model 2                      │
│                                                                              │
│  TRIPLEX                                                                     │
│  ├─ Mirror Repeat        : H-DNA forming sequences                          │
│  └─ Sticky DNA           : GAA/TTC triplet repeats                          │
│                                                                              │
│  G-QUADRUPLEX                                                                │
│  ├─ Canonical G4         : Standard G3+N1-7G3+N1-7G3+N1-7G3+ pattern        │
│  ├─ Bulged G4            : G-quadruplex with bulge insertions               │
│  ├─ Relaxed G4           : Longer loop variants (N1-12)                     │
│  ├─ Bipartite G4         : Two-block G4 structures                          │
│  ├─ Multimeric G4        : Intermolecular G4 formations                     │
│  ├─ Imperfect G4         : G4 with single G interruptions                   │
│  └─ G-Triplex            : Three G-tract intermediate structure             │
│                                                                              │
│  i-MOTIF                                                                     │
│  ├─ Canonical i-Motif    : C-rich intercalated structure                    │
│  ├─ Extended i-Motif     : Longer loop i-motif variants                     │
│  └─ AC-Motif             : Adenine-cytosine alternating pattern             │
│                                                                              │
│  Z-DNA                                                                       │
│  ├─ Classic Z-DNA        : (CG)n or (CA)n left-handed helix                 │
│  └─ eGZ (Extruded-G)     : Extended Z-DNA with G extrusion                  │
│                                                                              │
│  A-PHILIC DNA                                                                │
│  └─ A-philic DNA         : Adenine-rich tetranucleotide patterns            │
│                                                                              │
│  HYBRID                                                                      │
│  └─ Multi-class Overlap  : 30-99% overlap between different classes         │
│                                                                              │
│  CLUSTERS                                                                    │
│  └─ Non-B DNA Cluster    : ≥3 motifs, ≥2 classes within 500bp window       │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

## Parameters Reference

### Chunking Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `CHUNK_THRESHOLD` | 10,000 bp | Auto-chunking threshold |
| `DEFAULT_CHUNK_SIZE` | 10,000 bp | Standard chunk size |
| `DEFAULT_CHUNK_OVERLAP` | 500 bp | Overlap to prevent boundary misses |

### Repeat Detection Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `DIRECT_MIN_UNIT` | 10 bp | Minimum direct repeat unit |
| `DIRECT_MAX_UNIT` | 50 bp | Maximum direct repeat unit |
| `DIRECT_MAX_SPACER` | 0 bp | Maximum spacer between repeats |
| `INVERTED_MIN_ARM` | 10 bp | Minimum inverted repeat arm |
| `INVERTED_MAX_LOOP` | 3 bp | Maximum loop in cruciform |
| `MIRROR_MIN_ARM` | 10 bp | Minimum mirror repeat arm |
| `MIRROR_MAX_LOOP` | 8 bp | Maximum triplex spacer |
| `STR_MIN_UNIT` | 1 bp | Minimum STR unit size |
| `STR_MAX_UNIT` | 9 bp | Maximum STR unit size |
| `STR_MIN_TOTAL` | 10 bp | Minimum total STR length |

### K-mer Indexing Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `K_DIRECT` | 10 | K-mer size for direct repeats |
| `K_INVERTED` | 10 | K-mer size for inverted repeats |
| `K_MIRROR` | 10 | K-mer size for mirror repeats |
| `MAX_POSITIONS_PER_KMER` | 10,000 | Safety limit for frequent k-mers |

---

## Output Format

### Motif Dictionary Structure

Each detected motif is returned as a dictionary with these fields:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                          MOTIF OUTPUT STRUCTURE                              │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Field         │ Type   │ Always │ Description                              │
│ ───────────────┼────────┼────────┼──────────────────────────────────────────│
│  ID            │ str    │  Yes   │ Unique identifier (seq_CLASS_pos)        │
│  Sequence_Name │ str    │  Yes   │ Source sequence identifier               │
│  Class         │ str    │  Yes   │ Motif class (e.g., 'G-Quadruplex')       │
│  Subclass      │ str    │  Yes   │ Motif subclass (e.g., 'Canonical_G4')    │
│  Start         │ int    │  Yes   │ 1-based start position                   │
│  End           │ int    │  Yes   │ End position (inclusive)                 │
│  Length        │ int    │  Yes   │ Motif length in base pairs               │
│  Sequence      │ str    │  Yes   │ Actual DNA sequence                      │
│  Score         │ float  │  Yes   │ Normalized confidence (1-3 scale)        │
│  Strand        │ str    │  Yes   │ '+' or '-' strand                        │
│  Method        │ str    │  Yes   │ Detection algorithm used                 │
│ ───────────────┼────────┼────────┼──────────────────────────────────────────│
│  Component_Classes │ list │ Hybrid│ List of overlapping classes             │
│  Motif_Count   │ int    │ Cluster│ Number of motifs in cluster              │
│  Class_Diversity│ int   │ Cluster│ Number of distinct classes               │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### Score Interpretation

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                         SCORE INTERPRETATION GUIDE                           │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Score Range  │ Interpretation         │ Biological Significance            │
│ ──────────────┼────────────────────────┼─────────────────────────────────────│
│  1.0 - 1.5    │ Minimal                │ Weak motif formation potential      │
│  1.5 - 2.0    │ Low-Moderate           │ Below-average formation potential   │
│  2.0 - 2.5    │ Moderate               │ Typical formation potential         │
│  2.5 - 2.8    │ High                   │ Above-average formation potential   │
│  2.8 - 3.0    │ Highest                │ Strong formation potential          │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

## Performance

### Benchmarks

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                         PERFORMANCE BENCHMARKS                               │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Sequence Size  │ Standard Mode │ Fast Mode (9x) │ Memory Usage             │
│ ────────────────┼───────────────┼────────────────┼──────────────────────────│
│  1 KB           │ < 0.5s        │ < 0.1s         │ ~2 MB                    │
│  10 KB          │ ~2s           │ ~0.3s          │ ~3 MB                    │
│  100 KB         │ ~15s          │ ~2s            │ ~5 MB                    │
│  1 MB           │ ~150s         │ ~18s           │ ~15 MB                   │
│  10 MB          │ ~25 min       │ ~3 min         │ ~100 MB                  │
│                                                                              │
│  Processing Speed:                                                           │
│  • Standard mode: ~5,800 bp/second (all detectors)                          │
│  • Fast mode: ~52,000 bp/second (wall-clock time)                           │
│  • Slipped DNA detector: ~280,000 bp/second                                 │
│  • Memory efficiency: ~5 MB for 100K sequences                              │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### Performance Optimization Tips

| Scenario | Recommendation |
|----------|----------------|
| Large genome (>1MB) | Use `use_fast_mode=True` and `use_chunking=True` |
| Many small sequences | Use `analyze_multiple_sequences()` with multiprocessing |
| Memory constrained | Reduce `chunk_size`, use streaming output |
| Maximum speed | Use CLI with `--workers 8+` and Hyperscan backend |
| Interactive analysis | Use standard mode for real-time progress tracking |

---

## Usage Examples

### Basic Analysis

```python
import nonbscanner as nbs

# Simple sequence analysis
sequence = "GGGTTAGGGTTAGGGTTAGGGAAAAATTTTCGCGCGCGCG"
motifs = nbs.analyze_sequence(sequence, "example")

print(f"Found {len(motifs)} motifs:")
for motif in motifs:
    print(f"  {motif['Class']}: {motif['Start']}-{motif['End']} (score: {motif['Score']:.2f})")
```

### Batch Analysis with Progress

```python
import nonbscanner as nbs

# Progress callback
def on_progress(chunk, total, bp, elapsed, throughput):
    pct = (chunk / total) * 100
    print(f"\r[{'█' * int(pct/2)}{' ' * (50-int(pct/2))}] {pct:.1f}% - {throughput:,.0f} bp/s", end='')

# Analyze large sequence
results = nbs.analyze_sequence(
    large_sequence,
    "chr1",
    use_chunking=True,
    progress_callback=on_progress
)
print(f"\nComplete! Found {len(results)} motifs")
```

### Export Results

```python
import nonbscanner as nbs

# Analyze and export
motifs = nbs.analyze_file("genome.fasta")
all_motifs = [m for mlist in motifs.values() for m in mlist]

# Export to multiple formats
nbs.export_results(all_motifs, format='csv', filename='results.csv')
nbs.export_results(all_motifs, format='bed', filename='results.bed')
nbs.export_results(all_motifs, format='excel', filename='results.xlsx')
nbs.export_results(all_motifs, format='json', filename='results.json')
```

### Web Interface

```bash
# Launch Streamlit web interface
streamlit run app.py

# Access at http://localhost:8501
```

### Command Line Interface

```bash
# Basic CLI usage
python cli/nbdscanner_cli.py --input genome.fa --out results.ndjson

# Parallel processing
python cli/nbdscanner_cli.py --input genome.fa --out results.ndjson --workers 8

# Multiple output formats
python cli/nbdscanner_cli.py --input genome.fa --out results --format ndjson,csv,bed
```

---

## Additional Resources

- **README.md**: Quick start guide and overview
- **docs/perf_runbook.md**: Performance optimization guide
- **consolidated_registry.json**: Complete pattern definitions

---

*NonBScanner - Professional Non-B DNA Motif Detection Suite*  
*Version 2024.1 | MIT License*  
*Author: Dr. Venkata Rajesh Yella*
