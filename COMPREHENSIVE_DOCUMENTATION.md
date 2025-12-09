# NonBScanner - Comprehensive Parameters, Flow Diagrams, and API Reference
# =========================================================================

**Version:** 2024.1 - Professional Edition  
**Author:** Dr. Venkata Rajesh Yella  
**License:** MIT

---

## Table of Contents

1. [System Overview](#1-system-overview)
2. [Complete Parameter Reference](#2-complete-parameter-reference)
3. [Code Flow Diagrams](#3-code-flow-diagrams)
4. [API Function Parameters](#4-api-function-parameters)
5. [Detector Class Parameters](#5-detector-class-parameters)
6. [Export and Output Formats](#6-export-and-output-formats)
7. [Performance Tuning Parameters](#7-performance-tuning-parameters)

---

## 1. System Overview

### 1.1 Architecture Summary

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                          SYSTEM ARCHITECTURE                                 │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Module             │ Lines │ Purpose                                       │
│ ────────────────────┼───────┼───────────────────────────────────────────────│
│  nonbscanner.py     │ 1577  │ Main API & Scanner Orchestration              │
│  detectors.py       │ 4055  │ All 9 Motif Detector Classes                  │
│  utilities.py       │ 3110  │ Sequence I/O, Export & Statistics             │
│  visualizations.py  │ 2235  │ Complete Visualization Suite                  │
│  scanner.py         │ 1728  │ K-mer Indexing Functions                      │
│  app.py             │ 3093  │ Streamlit Web Interface                       │
│  parallel_scanner.py│  505  │ Parallel Processing Engine                    │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 1.2 Supported Motif Classes

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                         MOTIF DETECTION COVERAGE                             │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Class # │ Name             │ Subclasses │ Detection Method                 │
│ ─────────┼──────────────────┼────────────┼──────────────────────────────────│
│     1    │ Curved DNA       │     2      │ A-tract phasing analysis         │
│     2    │ Slipped DNA      │     2      │ K-mer indexed repeat search      │
│     3    │ Cruciform        │     1      │ Inverted repeat detection        │
│     4    │ R-Loop           │     3      │ QmRLFS algorithm                 │
│     5    │ Triplex          │     2      │ Mirror repeat & purine runs      │
│     6    │ G-Quadruplex     │     7      │ G4Hunter + pattern matching      │
│     7    │ i-Motif          │     3      │ C-run pattern detection          │
│     8    │ Z-DNA            │     2      │ CG/CA repeat scoring             │
│     9    │ A-philic DNA     │     1      │ 10-mer tetranucleotide search    │
│    10    │ Hybrid           │  Dynamic   │ Overlap detection                │
│    11    │ Clusters         │  Dynamic   │ Density-based windowing          │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

## 2. Complete Parameter Reference

### 2.1 Global Constants

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                        GLOBAL SYSTEM PARAMETERS                              │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Parameter Name         │ Value  │ Unit │ Description                       │
│ ────────────────────────┼────────┼──────┼───────────────────────────────────│
│  CHUNK_THRESHOLD        │ 10,000 │  bp  │ Auto-chunking activation size     │
│  DEFAULT_CHUNK_SIZE     │ 10,000 │  bp  │ Standard chunk size               │
│  DEFAULT_CHUNK_OVERLAP  │    500 │  bp  │ Overlap to prevent boundary miss  │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 2.2 Repeat Detection Parameters

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                      REPEAT DETECTION PARAMETERS                             │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Parameter Name       │ Value │ Unit │ Purpose                              │
│ ──────────────────────┼───────┼──────┼──────────────────────────────────────│
│  DIRECT_MIN_UNIT      │   10  │  bp  │ Minimum direct repeat unit size      │
│  DIRECT_MAX_UNIT      │   50  │  bp  │ Maximum direct repeat unit size      │
│  DIRECT_MAX_SPACER    │    0  │  bp  │ Maximum spacer between repeats       │
│  INVERTED_MIN_ARM     │   10  │  bp  │ Minimum inverted repeat arm length   │
│  INVERTED_MAX_LOOP    │    3  │  bp  │ Maximum loop in cruciform            │
│  MIRROR_MIN_ARM       │   10  │  bp  │ Minimum mirror repeat arm length     │
│  MIRROR_MAX_LOOP      │    8  │  bp  │ Maximum triplex spacer               │
│  STR_MIN_UNIT         │    1  │  bp  │ Minimum STR unit size                │
│  STR_MAX_UNIT         │    9  │  bp  │ Maximum STR unit size                │
│  STR_MIN_TOTAL        │   20  │  bp  │ Minimum total STR length             │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 2.3 K-mer Indexing Parameters

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                       K-MER INDEXING PARAMETERS                              │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Parameter Name          │ Value  │ Unit │ Description                      │
│ ─────────────────────────┼────────┼──────┼──────────────────────────────────│
│  K_DIRECT                │   10   │  bp  │ K-mer size for direct repeats    │
│  K_INVERTED              │   10   │  bp  │ K-mer size for inverted repeats  │
│  K_MIRROR                │   10   │  bp  │ K-mer size for mirror repeats    │
│  MAX_POSITIONS_PER_KMER  │ 10,000 │ cnt  │ Safety limit for frequent k-mers │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 2.4 Hybrid and Cluster Detection Parameters

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                   HYBRID & CLUSTER DETECTION PARAMETERS                      │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Parameter Name           │ Value  │ Unit │ Description                     │
│ ──────────────────────────┼────────┼──────┼─────────────────────────────────│
│  HYBRID_MIN_OVERLAP       │  0.3   │  %   │ Minimum overlap for hybrid (30%)│
│  HYBRID_MAX_OVERLAP       │  1.0   │  %   │ Maximum overlap for hybrid (99%)│
│  CLUSTER_WINDOW_SIZE      │  500   │  bp  │ Sliding window size for clusters│
│  CLUSTER_MIN_MOTIFS       │    3   │ cnt  │ Minimum motifs to form cluster  │
│  CLUSTER_MIN_CLASSES      │    2   │ cnt  │ Minimum classes for mixed type  │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 2.5 Score Normalization Parameters

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                      SCORE NORMALIZATION (1-3 SCALE)                         │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Motif Class      │ Min Score │ Max Score │ Normalization Method           │
│ ──────────────────┼───────────┼───────────┼────────────────────────────────│
│  Curved DNA       │    0.1    │    0.5    │ Linear mapping                 │
│  Slipped DNA      │    10     │    100    │ Log scale                      │
│  Cruciform        │    10     │    100    │ Log scale                      │
│  R-Loop           │   -50     │    50     │ Centered linear                │
│  Triplex          │    10     │    100    │ Log scale                      │
│  G-Quadruplex     │   -2.0    │    4.0    │ G4Hunter score                 │
│  i-Motif          │    0.5    │    3.0    │ Linear mapping                 │
│  Z-DNA            │    0.0    │    10.0   │ Sigmoid curve                  │
│  A-philic         │    1      │    5      │ Count-based                    │
│                                                                              │
│  Normalized Range │ 1.0 (minimal) to 3.0 (highest formation potential)     │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

## 3. Code Flow Diagrams

### 3.1 Main Analysis Pipeline

```
┌────────────────────────────────────────────────────────────────────────────┐
│                       ANALYSIS PIPELINE FLOWCHART                          │
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

### 3.2 Detector Execution Sequence

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

### 3.3 K-mer Indexing Flow (Repeat Detection)

```
┌────────────────────────────────────────────────────────────────────────────┐
│                      K-MER INDEXING WORKFLOW                               │
└────────────────────────────────────────────────────────────────────────────┘

                           ┌─────────────────┐
                           │  Input Sequence │
                           │   (length n)    │
                           └────────┬────────┘
                                    │
                                    ▼
                  ┌──────────────────────────────────┐
                  │  Extract All K-mers (k=10)       │
                  │  • Generate k-mer dictionary     │
                  │  • Store positions for each k-mer│
                  └────────────┬─────────────────────┘
                               │
                               ▼
                  ┌──────────────────────────────────┐
                  │  Filter Frequent K-mers          │
                  │  • Skip if >10,000 occurrences   │
                  │  • Prevents memory explosion     │
                  └────────────┬─────────────────────┘
                               │
          ┌────────────────────┼────────────────────┐
          │                    │                    │
          ▼                    ▼                    ▼
┌──────────────────┐ ┌──────────────────┐ ┌──────────────────┐
│ Direct Repeats   │ │ Inverted Repeats │ │ Mirror Repeats   │
│ • Min unit: 10bp │ │ • Min arm: 10bp  │ │ • Min arm: 10bp  │
│ • Max unit: 50bp │ │ • Max loop: 3bp  │ │ • Max loop: 8bp  │
│ • Max spacer: 0  │ │ • Palindromic    │ │ • Triplex motifs │
└─────────┬────────┘ └─────────┬────────┘ └─────────┬────────┘
          │                    │                    │
          └────────────────────┼────────────────────┘
                               │
                               ▼
                  ┌──────────────────────────────────┐
                  │  Validate and Score Matches      │
                  │  • Check length constraints      │
                  │  • Calculate stability scores    │
                  │  • Filter by quality threshold   │
                  └────────────┬─────────────────────┘
                               │
                               ▼
                  ┌──────────────────────────────────┐
                  │  Return Motif List               │
                  │  • Sorted by position            │
                  │  • With score and metadata       │
                  └──────────────────────────────────┘
```

### 3.4 Hybrid Detection Algorithm

```
┌────────────────────────────────────────────────────────────────────────────┐
│                     HYBRID MOTIF DETECTION FLOW                            │
└────────────────────────────────────────────────────────────────────────────┘

                           ┌─────────────────┐
                           │  Input: Motifs  │
                           │  (sorted by pos)│
                           └────────┬────────┘
                                    │
                                    ▼
                  ┌──────────────────────────────────┐
                  │  Sort by Start Position          │
                  │  • Enables sweep line algorithm  │
                  │  • O(n log n) complexity         │
                  └────────────┬─────────────────────┘
                               │
                               ▼
                  ┌──────────────────────────────────┐
                  │  For Each Motif Pair             │
                  │  • Compare adjacent motifs only  │
                  │  • Early termination if no overlap│
                  └────────────┬─────────────────────┘
                               │
                        ┌──────┴──────┐
                        │   Classes   │
                        │  Different? │
                        └──────┬──────┘
                               │
                          YES  │  NO → Skip
                               ▼
                  ┌──────────────────────────────────┐
                  │  Calculate Overlap Percentage    │
                  │  overlap = intersection / union  │
                  └────────────┬─────────────────────┘
                               │
                     ┌─────────┴─────────┐
                     │  30% < overlap    │
                     │      < 99% ?      │
                     └─────────┬─────────┘
                               │
                          YES  │  NO → Skip
                               ▼
                  ┌──────────────────────────────────┐
                  │  Create Hybrid Motif             │
                  │  • Class: 'Hybrid'               │
                  │  • Subclass: 'ClassA_ClassB'     │
                  │  • Score: average of both        │
                  │  • Sequence: merged region       │
                  └────────────┬─────────────────────┘
                               │
                               ▼
                  ┌──────────────────────────────────┐
                  │  Deduplicate Hybrids             │
                  │  • Track seen (pos, classes)     │
                  │  • Avoid reporting same overlap  │
                  └────────────┬─────────────────────┘
                               │
                               ▼
                  ┌──────────────────────────────────┐
                  │  Return Hybrid Motif List        │
                  └──────────────────────────────────┘
```

### 3.5 Cluster Detection Algorithm

```
┌────────────────────────────────────────────────────────────────────────────┐
│                     CLUSTER DETECTION WORKFLOW                             │
└────────────────────────────────────────────────────────────────────────────┘

                           ┌─────────────────┐
                           │  Input: Motifs  │
                           │  (min 3 needed) │
                           └────────┬────────┘
                                    │
                                    ▼
                  ┌──────────────────────────────────┐
                  │  Initialize Sliding Window       │
                  │  • Window size: 500 bp           │
                  │  • Min motifs: 3                 │
                  │  • Min classes: 2                │
                  └────────────┬─────────────────────┘
                               │
                               ▼
                  ┌──────────────────────────────────┐
                  │  For Each Window Position        │
                  │  • Slide from start to end       │
                  └────────────┬─────────────────────┘
                               │
                               ▼
                  ┌──────────────────────────────────┐
                  │  Count Motifs in Window          │
                  │  • Overlapping with window       │
                  │  • Track distinct classes        │
                  └────────────┬─────────────────────┘
                               │
                     ┌─────────┴─────────┐
                     │  ≥3 motifs AND    │
                     │  ≥2 classes?      │
                     └─────────┬─────────┘
                               │
                          YES  │  NO → Next window
                               ▼
                  ┌──────────────────────────────────┐
                  │  Create Cluster Motif            │
                  │  • Class: 'Non-B_DNA_Clusters'   │
                  │  • Subclass: 'Mixed_Cluster_N'   │
                  │  • Start/End: window bounds      │
                  │  • Count motifs and classes      │
                  └────────────┬─────────────────────┘
                               │
                               ▼
                  ┌──────────────────────────────────┐
                  │  Calculate Cluster Score         │
                  │  • Average of all motif scores   │
                  │  • Weighted by motif count       │
                  └────────────┬─────────────────────┘
                               │
                               ▼
                  ┌──────────────────────────────────┐
                  │  Deduplicate Clusters            │
                  │  • Remove overlapping clusters   │
                  │  • Keep highest density          │
                  └────────────┬─────────────────────┘
                               │
                               ▼
                  ┌──────────────────────────────────┐
                  │  Return Cluster List             │
                  └──────────────────────────────────┘
```


---

## 4. API Function Parameters

### 4.1 analyze_sequence() - Main Analysis Function

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                     analyze_sequence() PARAMETERS                            │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Parameter          │ Type     │ Default   │ Description                    │
│ ────────────────────┼──────────┼───────────┼────────────────────────────────│
│  sequence           │ str      │ Required  │ DNA sequence (ATGC)            │
│  sequence_name      │ str      │ "sequence"│ Identifier for sequence        │
│  use_fast_mode      │ bool     │ False     │ Enable parallel processing     │
│  use_chunking       │ bool     │ None      │ Enable chunking (auto if >10kb)│
│  chunk_size         │ int      │ 10000     │ Size of each chunk (bp)        │
│  chunk_overlap      │ int      │ 500       │ Overlap between chunks (bp)    │
│  progress_callback  │ Callable │ None      │ Progress tracking function     │
│  use_parallel_chunks│ bool     │ True      │ Parallel chunk processing      │
│                                                                              │
│  Returns: List[Dict[str, Any]] - Motif dictionaries sorted by position      │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

**Progress Callback Signature:**
```python
def callback(
    chunk_num: int,        # Current chunk number (1 to total_chunks)
    total_chunks: int,     # Total number of chunks
    bp_processed: int,     # Base pairs processed so far
    elapsed: float,        # Elapsed time in seconds
    throughput: float      # Processing speed in bp/s
) -> None
```

### 4.2 analyze_file() - File Input Function

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                       analyze_file() PARAMETERS                              │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Parameter   │ Type │ Default  │ Description                                │
│ ─────────────┼──────┼──────────┼────────────────────────────────────────────│
│  filename    │ str  │ Required │ Path to FASTA file                         │
│                                                                              │
│  Returns: Dict[str, List[Dict]] - Mapping sequence names to motif lists     │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 4.3 analyze_fasta() - Direct FASTA Content Analysis

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                      analyze_fasta() PARAMETERS                              │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Parameter      │ Type │ Default  │ Description                             │
│ ────────────────┼──────┼──────────┼─────────────────────────────────────────│
│  fasta_content  │ str  │ Required │ FASTA format string                     │
│                                                                              │
│  Returns: Dict[str, List[Dict]] - Mapping sequence names to motif lists     │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 4.4 export_results() - Data Export Function

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                      export_results() PARAMETERS                             │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Parameter  │ Type       │ Default  │ Description                           │
│ ────────────┼────────────┼──────────┼───────────────────────────────────────│
│  motifs     │ List[Dict] │ Required │ List of motif dictionaries            │
│  format     │ str        │ 'csv'    │ Output format: csv, bed, json, excel  │
│  filename   │ str        │ None     │ Optional output filename              │
│  **kwargs   │ Dict       │ {}       │ Format-specific options               │
│                                                                              │
│  Supported Formats:                                                          │
│  • 'csv'    - Comma-separated values                                         │
│  • 'bed'    - UCSC Genome Browser format                                     │
│  • 'json'   - JavaScript Object Notation                                     │
│  • 'excel'  - Microsoft Excel (.xlsx) with separate sheets per class        │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

## 5. Detector Class Parameters

### 5.1 CurvedDNADetector Parameters

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                        CURVED DNA DETECTOR                                   │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Parameter Name       │ Value     │ Description                             │
│ ──────────────────────┼───────────┼─────────────────────────────────────────│
│  A_TRACT_MIN_LENGTH   │ 4 bp      │ Minimum A-tract length                  │
│  A_TRACT_MAX_LENGTH   │ 6 bp      │ Maximum A-tract length                  │
│  PHASING_PERIOD       │ ~10-11 bp │ Expected helical turn period            │
│  MIN_CURVATURE_SCORE  │ 0.1       │ Minimum curvature score threshold       │
│  WINDOW_SIZE          │ 100 bp    │ Sliding window for local curvature      │
│                                                                              │
│  Detection Method: A-tract phasing analysis                                  │
│  Complexity: O(n)                                                            │
│  Speed: ~8,000 bp/s                                                          │
│                                                                              │
│  Subclasses:                                                                 │
│  • Local Curvature   - Regional DNA flexibility                             │
│  • Global Curvature  - A-tract induced overall bending                      │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 5.2 SlippedDNADetector Parameters

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                        SLIPPED DNA DETECTOR                                  │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Parameter Name       │ Value      │ Description                            │
│ ──────────────────────┼────────────┼────────────────────────────────────────│
│  DIRECT_MIN_UNIT      │ 10 bp      │ Minimum direct repeat unit size        │
│  DIRECT_MAX_UNIT      │ 50 bp      │ Maximum direct repeat unit size        │
│  DIRECT_MAX_SPACER    │ 0 bp       │ Maximum spacer between repeats         │
│  STR_MIN_UNIT         │ 1 bp       │ Minimum STR unit size                  │
│  STR_MAX_UNIT         │ 9 bp       │ Maximum STR unit size                  │
│  STR_MIN_TOTAL        │ 20 bp      │ Minimum total STR length               │
│  K_MER_SIZE           │ 10 bp      │ K-mer size for indexing                │
│                                                                              │
│  Detection Method: K-mer indexed repeat search                               │
│  Complexity: O(n)                                                            │
│  Speed: ~280,000 bp/s (fastest detector)                                     │
│                                                                              │
│  Subclasses:                                                                 │
│  • Direct Repeat (DR)  - Tandem repeats (10-50 bp units)                    │
│  • Short Tandem Repeat - STRs (1-9 bp units, ≥20 bp total)                  │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 5.3 CruciformDetector Parameters

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                        CRUCIFORM DETECTOR                                    │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Parameter Name       │ Value    │ Description                              │
│ ──────────────────────┼──────────┼──────────────────────────────────────────│
│  INVERTED_MIN_ARM     │ 10 bp    │ Minimum inverted repeat arm length       │
│  INVERTED_MAX_ARM     │ 100 bp   │ Maximum inverted repeat arm length       │
│  INVERTED_MAX_LOOP    │ 3 bp     │ Maximum loop size in cruciform           │
│  K_MER_SIZE           │ 10 bp    │ K-mer size for indexing                  │
│  MIN_STABILITY_SCORE  │ 10       │ Minimum thermodynamic stability          │
│                                                                              │
│  Detection Method: Inverted repeat detection via k-mer indexing              │
│  Complexity: O(n)                                                            │
│  Speed: ~8,000 bp/s                                                          │
│                                                                              │
│  Subclasses:                                                                 │
│  • Palindromic Inverted Repeat - Four-way junction structures               │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 5.4 RLoopDetector Parameters

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                           R-LOOP DETECTOR                                    │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Parameter Name       │ Value    │ Description                              │
│ ──────────────────────┼──────────┼──────────────────────────────────────────│
│  WINDOW_SIZE          │ 50 bp    │ Analysis window size                     │
│  MIN_GC_SKEW          │ 0.3      │ Minimum GC-skew threshold                │
│  QMRLFS_THRESHOLD_M1  │ 0.0      │ QmRLFS model 1 threshold                 │
│  QMRLFS_THRESHOLD_M2  │ 0.0      │ QmRLFS model 2 threshold                 │
│  RIZ_MIN_LENGTH       │ 15 bp    │ Minimum RIZ region length                │
│  REZ_MIN_LENGTH       │ 15 bp    │ Minimum REZ region length                │
│                                                                              │
│  Detection Methods:                                                          │
│  • QmRLFS Algorithm - Quantitative R-loop formation scoring                  │
│  • GC-skew Analysis - Strand asymmetry detection                             │
│  Complexity: O(n)                                                            │
│  Speed: ~6,000 bp/s                                                          │
│                                                                              │
│  Subclasses:                                                                 │
│  • R-loop Formation   - RNA-DNA hybrid regions                               │
│  • QmRLFS-m1          - Quantitative model 1                                 │
│  • QmRLFS-m2          - Quantitative model 2                                 │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 5.5 TriplexDetector Parameters

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                          TRIPLEX DETECTOR                                    │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Parameter Name       │ Value    │ Description                              │
│ ──────────────────────┼──────────┼──────────────────────────────────────────│
│  MIRROR_MIN_ARM       │ 10 bp    │ Minimum mirror repeat arm length         │
│  MIRROR_MAX_ARM       │ 100 bp   │ Maximum mirror repeat arm length         │
│  MIRROR_MAX_LOOP      │ 8 bp     │ Maximum spacer length                    │
│  PURINE_RUN_MIN       │ 15 bp    │ Minimum purine run length                │
│  PYRIMIDINE_RUN_MIN   │ 15 bp    │ Minimum pyrimidine run length            │
│  GAA_TTC_MIN_REPEATS  │ 3        │ Minimum GAA/TTC triplet repeats          │
│  K_MER_SIZE           │ 10 bp    │ K-mer size for mirror repeats            │
│                                                                              │
│  Detection Methods:                                                          │
│  • Mirror repeat detection via k-mer indexing                                │
│  • Purine/pyrimidine run identification                                      │
│  • GAA/TTC sticky DNA pattern matching                                       │
│  Complexity: O(n)                                                            │
│  Speed: ~7,000 bp/s                                                          │
│                                                                              │
│  Subclasses:                                                                 │
│  • Mirror Repeat - H-DNA forming sequences                                   │
│  • Sticky DNA    - GAA/TTC triplet repeats                                   │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 5.6 GQuadruplexDetector Parameters

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                       G-QUADRUPLEX DETECTOR                                  │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Parameter Name       │ Value      │ Description                            │
│ ──────────────────────┼────────────┼────────────────────────────────────────│
│  G4HUNTER_WINDOW      │ 25 bp      │ G4Hunter sliding window size           │
│  G4HUNTER_THRESHOLD   │ 1.0        │ Minimum G4Hunter score                 │
│  MIN_G_RUNS           │ 4          │ Minimum number of G-runs               │
│  MIN_G_PER_RUN        │ 2          │ Minimum Gs per run                     │
│  MAX_G_PER_RUN        │ 4          │ Maximum Gs per run                     │
│  MIN_LOOP_LENGTH      │ 1 bp       │ Minimum loop length                    │
│  MAX_LOOP_LENGTH      │ 7 bp       │ Maximum loop length (canonical)        │
│  MAX_LOOP_RELAXED     │ 12 bp      │ Maximum loop length (relaxed)          │
│  BULGE_MAX_SIZE       │ 3 bp       │ Maximum bulge size                     │
│                                                                              │
│  Detection Methods:                                                          │
│  • G4Hunter algorithm - Sliding window G-richness scoring                    │
│  • Pattern matching   - Regex for G4 structure variants                      │
│  Complexity: O(n)                                                            │
│  Speed: ~5,000 bp/s                                                          │
│                                                                              │
│  Subclasses (7 types):                                                       │
│  • Canonical G4   - Standard G3+N1-7G3+N1-7G3+N1-7G3+ pattern               │
│  • Bulged G4      - G-quadruplex with bulge insertions                      │
│  • Relaxed G4     - Longer loop variants (N1-12)                            │
│  • Bipartite G4   - Two-block G4 structures                                 │
│  • Multimeric G4  - Intermolecular G4 formations                            │
│  • Imperfect G4   - G4 with single G interruptions                          │
│  • G-Triplex      - Three G-tract intermediate structure                    │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 5.7 IMotifDetector Parameters

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                          I-MOTIF DETECTOR                                    │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Parameter Name       │ Value    │ Description                              │
│ ──────────────────────┼──────────┼──────────────────────────────────────────│
│  MIN_C_RUNS           │ 4        │ Minimum number of C-runs                 │
│  MIN_C_PER_RUN        │ 2        │ Minimum Cs per run                       │
│  MAX_C_PER_RUN        │ 4        │ Maximum Cs per run                       │
│  MIN_LOOP_LENGTH      │ 1 bp     │ Minimum loop length                      │
│  MAX_LOOP_LENGTH      │ 7 bp     │ Maximum loop length (canonical)          │
│  MAX_LOOP_EXTENDED    │ 15 bp    │ Maximum loop length (extended)           │
│  MIN_AC_ALTERNATING   │ 8 bp     │ Minimum AC alternating pattern           │
│                                                                              │
│  Detection Methods:                                                          │
│  • C-run pattern matching - Regex for C-rich structures                      │
│  • AC alternating pattern detection                                          │
│  Complexity: O(n)                                                            │
│  Speed: ~6,000 bp/s                                                          │
│                                                                              │
│  Subclasses:                                                                 │
│  • Canonical i-Motif - C-rich intercalated structure                         │
│  • Extended i-Motif  - Longer loop i-motif variants                          │
│  • AC-Motif          - Adenine-cytosine alternating pattern                  │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 5.8 ZDNADetector Parameters

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                           Z-DNA DETECTOR                                     │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Parameter Name       │ Value    │ Description                              │
│ ──────────────────────┼──────────┼──────────────────────────────────────────│
│  WINDOW_SIZE          │ 10 bp    │ Sliding window size for scoring          │
│  Z_SCORE_THRESHOLD    │ 0.5      │ Minimum Z-DNA score                      │
│  MIN_CG_DINUC         │ 3        │ Minimum CG dinucleotides                 │
│  MIN_CA_DINUC         │ 3        │ Minimum CA dinucleotides                 │
│  MIN_TOTAL_LENGTH     │ 12 bp    │ Minimum Z-DNA region length              │
│  CG_WEIGHT            │ 1.0      │ Scoring weight for CG                    │
│  CA_WEIGHT            │ 0.8      │ Scoring weight for CA                    │
│                                                                              │
│  Detection Method: 10-mer sliding window CG/CA repeat scoring                │
│  Complexity: O(n)                                                            │
│  Speed: ~7,000 bp/s                                                          │
│                                                                              │
│  Subclasses:                                                                 │
│  • Classic Z-DNA - (CG)n or (CA)n left-handed helix                          │
│  • eGZ           - Extended Z-DNA with G extrusion                           │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 5.9 APhilicDetector Parameters

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                         A-PHILIC DNA DETECTOR                                │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Parameter Name       │ Value     │ Description                             │
│ ──────────────────────┼───────────┼─────────────────────────────────────────│
│  TETRANUC_PATTERNS    │ 411 total │ Loaded from consolidated_registry.json  │
│  MIN_PATTERN_MATCHES  │ 1         │ Minimum pattern matches required        │
│  MIN_MOTIF_LENGTH     │ 10 bp     │ Minimum motif length                    │
│  A_CONTENT_THRESHOLD  │ 60%       │ Minimum adenine content                 │
│                                                                              │
│  Detection Method: 10-mer tetranucleotide pattern matching                   │
│  Complexity: O(n)                                                            │
│  Speed: ~7,000 bp/s                                                          │
│                                                                              │
│  Subclasses:                                                                 │
│  • A-philic DNA - Adenine-rich tetranucleotide patterns                      │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```



---

## 6. Export and Output Formats

### 6.1 Motif Output Structure

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                          MOTIF OUTPUT STRUCTURE                              │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Field          │ Type   │ Always │ Description                             │
│ ────────────────┼────────┼────────┼─────────────────────────────────────────│
│  ID             │ str    │  Yes   │ Unique identifier (seq_CLASS_pos)       │
│  Sequence_Name  │ str    │  Yes   │ Source sequence identifier              │
│  Class          │ str    │  Yes   │ Motif class (e.g., 'G-Quadruplex')      │
│  Subclass       │ str    │  Yes   │ Motif subclass (e.g., 'Canonical_G4')   │
│  Start          │ int    │  Yes   │ 1-based start position                  │
│  End            │ int    │  Yes   │ End position (inclusive)                │
│  Length         │ int    │  Yes   │ Motif length in base pairs              │
│  Sequence       │ str    │  Yes   │ Actual DNA sequence                     │
│  Score          │ float  │  Yes   │ Normalized confidence (1-3 scale)       │
│  Strand         │ str    │  Yes   │ '+' or '-' strand                       │
│  Method         │ str    │  Yes   │ Detection algorithm used                │
│ ────────────────┼────────┼────────┼─────────────────────────────────────────│
│  Component_     │ list   │ Hybrid │ List of overlapping classes             │
│   Classes       │        │  only  │                                         │
│  Motif_Count    │ int    │Cluster │ Number of motifs in cluster             │
│  Class_Diversity│ int    │Cluster │ Number of distinct classes              │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 6.2 CSV Export Format

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                            CSV OUTPUT FORMAT                                 │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Column Order (left to right):                                               │
│  1. ID              11. Sequence                                             │
│  2. Sequence_Name   12. Length                                               │
│  3. Class           13. Start                                                │
│  4. Subclass        14. End                                                  │
│  5. Score           15. Strand                                               │
│  6. Method          16. Component_Classes (if hybrid)                        │
│  7-10. Reserved     17. Motif_Count (if cluster)                             │
│                                                                              │
│  Delimiter: Comma (,)                                                        │
│  Header Row: Yes (column names)                                              │
│  Quote Character: " (double quote for fields with commas)                    │
│  Encoding: UTF-8                                                             │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 6.3 BED Export Format

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                            BED OUTPUT FORMAT                                 │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Standard BED6 Format (6 columns):                                           │
│  ───────────────────────────────────────────────────────────────────────────│
│  Column 1: chrom      - Sequence name                                        │
│  Column 2: chromStart - Start position (0-based)                             │
│  Column 3: chromEnd   - End position (exclusive)                             │
│  Column 4: name       - Motif ID (Class_Subclass_Start)                      │
│  Column 5: score      - Normalized score * 333 (0-1000 scale)                │
│  Column 6: strand     - '+' or '-'                                           │
│                                                                              │
│  Compatible with: UCSC Genome Browser, IGV, bedtools                         │
│  Note: Start position converted from 1-based to 0-based for BED standard     │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 6.4 JSON Export Format

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                           JSON OUTPUT FORMAT                                 │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Structure: Array of motif objects                                           │
│  ───────────────────────────────────────────────────────────────────────────│
│  [                                                                           │
│    {                                                                         │
│      "ID": "seq_G-Quadruplex_123",                                           │
│      "Sequence_Name": "chr1",                                                │
│      "Class": "G-Quadruplex",                                                │
│      "Subclass": "Canonical_G4",                                             │
│      "Start": 123,                                                           │
│      "End": 145,                                                             │
│      "Length": 22,                                                           │
│      "Sequence": "GGGTTAGGGTTAGGGTTAGGG",                                    │
│      "Score": 2.75,                                                          │
│      "Strand": "+",                                                          │
│      "Method": "G4Hunter"                                                    │
│    },                                                                        │
│    ...                                                                       │
│  ]                                                                           │
│                                                                              │
│  Encoding: UTF-8                                                             │
│  Indentation: 2 spaces (pretty-printed)                                      │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 6.5 Excel Export Format

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                          EXCEL OUTPUT FORMAT                                 │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Format: .xlsx (Microsoft Excel 2007+)                                       │
│  ───────────────────────────────────────────────────────────────────────────│
│  Sheet Structure:                                                            │
│  • Sheet 1: "Summary"     - Overview statistics                              │
│  • Sheet 2: "All_Motifs"  - Complete motif list                              │
│  • Sheet 3+: Individual sheets per motif class                               │
│                                                                              │
│  Per-Class Sheets (one per detected class):                                  │
│  • Curved_DNA                                                                │
│  • Slipped_DNA                                                               │
│  • Cruciform                                                                 │
│  • R_Loop                                                                    │
│  • Triplex                                                                   │
│  • G_Quadruplex                                                              │
│  • i_Motif                                                                   │
│  • Z_DNA                                                                     │
│  • A_Philic_DNA                                                              │
│  • Hybrid                                                                    │
│  • Clusters                                                                  │
│                                                                              │
│  Formatting:                                                                 │
│  • Header row: Bold, frozen                                                  │
│  • Score column: Number format (2 decimal places)                            │
│  • Auto-sized columns                                                        │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

## 7. Performance Tuning Parameters

### 7.1 Performance Mode Comparison

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                       PERFORMANCE MODE COMPARISON                            │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Mode        │ Size Range │ Throughput  │ Wall-Clock│ CPU Usage│ Memory    │
│ ─────────────┼────────────┼─────────────┼───────────┼──────────┼───────────│
│  Standard    │ <10kb      │ ~6,000 bp/s │ Baseline  │ 1 core   │ ~5 MB     │
│  Chunked     │ 10kb-100kb │ ~6,000 bp/s │ Baseline  │ 1 core   │ ~10 MB    │
│  Fast Mode   │ Any        │ ~6,000 bp/s │ 9x faster │ 9 cores  │ ~15 MB    │
│  Parallel    │ >100kb     │ ~6,000 bp/s │ 9x faster │ Multi    │ ~50 MB    │
│   Chunks     │            │ per chunk   │           │          │           │
│                                                                              │
│  Note: Fast mode provides wall-clock speedup by running detectors in         │
│  parallel, not throughput increase per detector.                             │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 7.2 Performance Benchmarks

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                         PERFORMANCE BENCHMARKS                               │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Sequence  │ Standard Mode │ Fast Mode    │ Memory  │ Speedup              │
│  Size      │ (Sequential)  │ (9 threads)  │ Usage   │ Factor               │
│ ───────────┼───────────────┼──────────────┼─────────┼──────────────────────│
│  1 KB      │ < 0.5s        │ < 0.1s       │ ~2 MB   │ 5x                   │
│  10 KB     │ ~2s           │ ~0.3s        │ ~3 MB   │ 6.7x                 │
│  100 KB    │ ~15s          │ ~2s          │ ~5 MB   │ 7.5x                 │
│  1 MB      │ ~150s         │ ~18s         │ ~15 MB  │ 8.3x                 │
│  10 MB     │ ~25 min       │ ~3 min       │ ~100 MB │ 8.3x                 │
│                                                                              │
│  Individual Detector Speeds (bp/s):                                          │
│  • Slipped DNA:  ~280,000 bp/s (fastest)                                     │
│  • Curved DNA:   ~8,000 bp/s                                                 │
│  • Cruciform:    ~8,000 bp/s                                                 │
│  • Triplex:      ~7,000 bp/s                                                 │
│  • Z-DNA:        ~7,000 bp/s                                                 │
│  • A-philic:     ~7,000 bp/s                                                 │
│  • i-Motif:      ~6,000 bp/s                                                 │
│  • R-Loop:       ~6,000 bp/s                                                 │
│  • G4:           ~5,000 bp/s                                                 │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 7.3 Optimization Guidelines

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                         OPTIMIZATION GUIDELINES                              │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Scenario                │ Recommendation                                    │
│ ─────────────────────────┼───────────────────────────────────────────────────│
│  Small sequences (<10kb) │ Use standard mode (analyze_sequence)              │
│  Medium (10kb-100kb)     │ Enable chunking (use_chunking=True)               │
│  Large genomes (>1MB)    │ Use fast mode + chunking                          │
│  Many small sequences    │ Use parallel file processing                      │
│  Memory constrained      │ Reduce chunk_size to 5000 bp                      │
│  Maximum speed           │ use_fast_mode=True + use_parallel_chunks=True     │
│  Real-time progress      │ Provide progress_callback function                │
│  Interactive analysis    │ Use standard mode for immediate feedback          │
│  Batch processing        │ Use CLI with --workers parameter                  │
│  Web application         │ Use chunking + progress for large uploads         │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 7.4 Memory Management

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                          MEMORY MANAGEMENT                                   │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Component           │ Memory Usage  │ Scaling Factor                       │
│ ─────────────────────┼───────────────┼──────────────────────────────────────│
│  Base overhead       │ ~2 MB         │ Constant                             │
│  Sequence storage    │ ~1 byte/bp    │ Linear with sequence length          │
│  K-mer index         │ ~40 bytes/kmer│ Depends on sequence complexity       │
│  Motif storage       │ ~200 bytes/   │ Linear with motif count              │
│                      │  motif        │                                      │
│  Pattern compilation │ ~1 MB         │ One-time cost                        │
│  Thread overhead     │ ~0.5 MB/      │ Only in fast/parallel mode           │
│                      │  thread       │                                      │
│                                                                              │
│  Estimated Total Memory:                                                     │
│  • 1 KB sequence:     ~2 MB                                                  │
│  • 10 KB sequence:    ~3 MB                                                  │
│  • 100 KB sequence:   ~5 MB                                                  │
│  • 1 MB sequence:     ~15 MB                                                 │
│  • 10 MB sequence:    ~100 MB                                                │
│                                                                              │
│  Memory Reduction Strategies:                                                │
│  • Use smaller chunk_size (reduces peak memory)                              │
│  • Process sequences one at a time (avoid loading all into memory)           │
│  • Use streaming mode for very large files                                   │
│  • Disable parallel mode if memory-constrained                               │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

## Appendix: Complete Parameter Summary Table

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                    ALL PARAMETERS - QUICK REFERENCE                          │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Category: Global System Parameters                                          │
│ ──────────────────────────────────────────────────────────────────────────── │
│  CHUNK_THRESHOLD              │ 10,000  │ bp  │ Auto-chunking activation    │
│  DEFAULT_CHUNK_SIZE           │ 10,000  │ bp  │ Standard chunk size         │
│  DEFAULT_CHUNK_OVERLAP        │    500  │ bp  │ Overlap size                │
│                                                                              │
│  Category: Repeat Detection (Slipped DNA, Cruciform, Triplex)                │
│ ──────────────────────────────────────────────────────────────────────────── │
│  DIRECT_MIN_UNIT              │     10  │ bp  │ Min direct repeat unit      │
│  DIRECT_MAX_UNIT              │     50  │ bp  │ Max direct repeat unit      │
│  DIRECT_MAX_SPACER            │      0  │ bp  │ Max spacer                  │
│  INVERTED_MIN_ARM             │     10  │ bp  │ Min inverted arm            │
│  INVERTED_MAX_LOOP            │      3  │ bp  │ Max cruciform loop          │
│  MIRROR_MIN_ARM               │     10  │ bp  │ Min mirror arm              │
│  MIRROR_MAX_LOOP              │      8  │ bp  │ Max triplex spacer          │
│  STR_MIN_UNIT                 │      1  │ bp  │ Min STR unit                │
│  STR_MAX_UNIT                 │      9  │ bp  │ Max STR unit                │
│  STR_MIN_TOTAL                │     20  │ bp  │ Min total STR length        │
│                                                                              │
│  Category: K-mer Indexing                                                    │
│ ──────────────────────────────────────────────────────────────────────────── │
│  K_DIRECT                     │     10  │ bp  │ K-mer for direct repeats    │
│  K_INVERTED                   │     10  │ bp  │ K-mer for inverted repeats  │
│  K_MIRROR                     │     10  │ bp  │ K-mer for mirror repeats    │
│  MAX_POSITIONS_PER_KMER       │ 10,000  │ cnt │ Frequency safety limit      │
│                                                                              │
│  Category: Hybrid & Cluster Detection                                        │
│ ──────────────────────────────────────────────────────────────────────────── │
│  HYBRID_MIN_OVERLAP           │    0.3  │ %   │ Min overlap (30%)           │
│  HYBRID_MAX_OVERLAP           │    1.0  │ %   │ Max overlap (99%)           │
│  CLUSTER_WINDOW_SIZE          │    500  │ bp  │ Window size                 │
│  CLUSTER_MIN_MOTIFS           │      3  │ cnt │ Min motifs in cluster       │
│  CLUSTER_MIN_CLASSES          │      2  │ cnt │ Min classes for mixed       │
│                                                                              │
│  Category: Score Normalization (all normalized to 1-3 scale)                 │
│ ──────────────────────────────────────────────────────────────────────────── │
│  Curved DNA min/max           │ 0.1/0.5 │     │ Linear mapping              │
│  Slipped DNA min/max          │ 10/100  │     │ Log scale                   │
│  Cruciform min/max            │ 10/100  │     │ Log scale                   │
│  R-Loop min/max               │ -50/50  │     │ Centered linear             │
│  Triplex min/max              │ 10/100  │     │ Log scale                   │
│  G-Quadruplex min/max         │ -2/4    │     │ G4Hunter score              │
│  i-Motif min/max              │ 0.5/3   │     │ Linear mapping              │
│  Z-DNA min/max                │ 0/10    │     │ Sigmoid curve               │
│  A-philic min/max             │ 1/5     │     │ Count-based                 │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

## Document Information

**Document Title:** NonBScanner Comprehensive Parameters, Flow Diagrams, and API Reference  
**Version:** 2024.1  
**Last Updated:** December 2024  
**Author:** Dr. Venkata Rajesh Yella  
**License:** MIT  
**Repository:** https://github.com/VRYella/NonBScanner  

**Purpose:** This document provides a complete reference for all parameters, algorithms, and 
workflows in the NonBScanner system. All diagrams are text-based (ASCII) for easy viewing in 
any text editor or terminal.

**Usage:** This document is intended for:
- Developers integrating NonBScanner into their workflows
- Researchers understanding the detection algorithms
- System administrators tuning performance parameters
- Users seeking detailed API documentation

**Related Documents:**
- `README.md` - Quick start guide
- `docs/DOCUMENTATION.md` - User-focused documentation
- `docs/perf_runbook.md` - Performance optimization guide
- Source code comments - Implementation details

---

*End of Comprehensive Documentation*
