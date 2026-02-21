# NonBDNAFinder — Detailed Materials and Methods

**Tool:** NonBDNAFinder v2025.1  
**Author:** Dr. Venkata Rajesh Yella  
**License:** MIT  
**Date:** 2026-02-21  

---

## 2. Materials and Methods

### 2.1 Architecture and Workflow

NonBDNAFinder is a Python-based, modular, and scalable pipeline for the genome-wide detection of non-B DNA secondary structures. The tool is composed of four tightly integrated layers:

1. **Input Layer** — FASTA parsing, sequence validation, and optional multi-FASTA processing via `Utilities/sequence_preprocessor.py` and `Utilities/utilities.py`.
2. **Detection Layer** — nine independent, class-specific detector modules in `Detectors/`, each implementing a `BaseMotifDetector` interface.
3. **Post-Processing Layer** — within-class overlap deduplication (`Utilities/overlap_deduplicator.py`), cross-class Hybrid annotation, and Cluster annotation (`Utilities/nonbscanner.py`).
4. **Output Layer** — multi-format export (`Utilities/final_exporter.py`) and an interactive Streamlit UI (`app.py`, `UI/`).

#### 2.1.1 Adaptive Execution Strategy

The central orchestrator is the `MotifDetectionEngine` (`Utilities/motif_detection_engine.py`), which selects one of three execution strategies based on input sequence length:

| Strategy | Trigger | Description |
|---|---|---|
| `numba_only` (direct) | < 100,000 bp | Single-pass in-memory analysis |
| `chunk_workers` (parallel) | 100,000 – 5,000,000 bp | Two-worker `ProcessPoolExecutor` with 50 kb chunks and 2 kb overlap |
| `disk_streaming` (sequential) | ≥ 5,000,000 bp | Chunk-by-chunk streaming with constant RAM via `DiskChunkManager` |

All strategies share the same detection, deduplication, and post-processing logic, ensuring identical results regardless of genome size. The workflow for a single sequence is:

```
[FASTA input]
     │
     ▼
[Sequence parsing & validation]
     │
     ▼
[Adaptive strategy selection]
     │
     ├──[< 100 kbp]──────────────────────► [Direct analysis (numba_only)]
     │                                                │
     ├──[100 kbp – 5 Mbp]──────────────► [Chunk workers (ProcessPool)]
     │                                                │
     └──[≥ 5 Mbp]────────────────────────► [Disk streaming (sequential)]
                                                      │
                                           [9 parallel detectors per chunk]
                                                      │
                                           [Overlap deduplication (core-region)]
                                                      │
                                           [Hybrid annotation]
                                                      │
                                           [Cluster annotation]
                                                      │
                                           [Score normalisation]
                                                      │
                                           [Export / Visualisation]
```

#### 2.1.2 Sequence Chunking

For sequences exceeding 100 kb, the input is divided into contiguous chunks of 50,000 bp with a 2,000 bp overlap. Overlapping regions ensure that motifs spanning chunk boundaries are detected by at least one chunk. After detection, only motifs whose start position falls within the *core region* `[chunk_start, chunk_start + chunk_size - overlap)` are retained, preventing double-counting of boundary-spanning motifs (see §2.6 for algorithmic complexity).

#### 2.1.3 Parallel Detector Execution

For sequences longer than 50,000 bp, all nine detectors are executed concurrently within each chunk using a `ThreadPoolExecutor` (up to 9 workers, one per detector class). This provides a 1.5–2× speed-up on multi-core systems with no change in output.

---

### 2.2 Motif Library and Structural Classes

NonBDNAFinder implements a canonical motif taxonomy defined in `Utilities/config/motif_taxonomy.py`. This taxonomy is the **single source of truth** for all class and subclass names throughout the pipeline, preventing biologically invalid label drift between detector, post-processor, and export modules.

#### 2.2.1 Canonical Class Taxonomy

The tool recognises **9 primary structural classes** plus two derived annotation types:

| ID | Class | Subclasses | Structural basis |
|---|---|---|---|
| 1 | `Curved_DNA` | Global Curvature, Local Curvature | A-tract phasing (Trifonov & Sussman, 1980) |
| 2 | `Slipped_DNA` | Direct Repeat, STR | Direct repeat slippage (Pearson et al., 2005) |
| 3 | `Cruciform` | Cruciform forming IRs | Palindromic inverted repeats (Lilley, 1980) |
| 4 | `R-Loop` | R-loop formation sites | G-skewed co-transcriptional RNA:DNA hybrids (Jenjaroenpun et al., 2015) |
| 5 | `Triplex` | Triplex, Sticky DNA | Mirror-repeat H-DNA (Frank-Kamenetskii & Mirkin, 1995) |
| 6 | `G-Quadruplex` | Canonical, Extended-loop, Bulged, Higher-order, G-triplex, Two-tetrad weak PQS, Stacked G4 | G-quartet stacking (Huppert & Balasubramanian, 2005) |
| 7 | `i-Motif` | Canonical i-motif, Relaxed i-motif, AC-motif | Intercalated C:CH⁺ hemi-protonated pairs (Zeraati et al., 2018) |
| 8 | `Z-DNA` | Z-DNA, eGZ | Alternating purine-pyrimidine left-handed helix (Ho et al., 1986) |
| 9 | `A-philic_DNA` | A-philic DNA | Tetranucleotide nucleosome-positioning propensity (Vinogradov, 2003) |
| 10 | `Hybrid` | Dynamic overlaps | Derived: ≥50% cross-class spatial overlap |
| 11 | `Non-B_DNA_Clusters` | Dynamic clusters | Derived: ≥4 motifs, ≥3 classes within 300 bp |

#### 2.2.2 G-Quadruplex Subclass Hierarchy

G-Quadruplex subclasses are ordered by a structural priority hierarchy that reflects thermodynamic stability and biological relevance:

```
Telomeric G4  >  Higher-order G4 array/G4-wire  >  Stacked G4
    >  Canonical intramolecular G4  >  Bulged G4
    >  Extended-loop canonical  >  Intramolecular G-triplex
    >  Two-tetrad weak PQS
```

When a region satisfies criteria for multiple G4 subclasses, the highest-priority subclass is assigned, preserving biological accuracy.

#### 2.2.3 Legacy Compatibility

All class and subclass names are normalised via `CLASS_ALIASES` and `SUBCLASS_ALIASES` dictionaries in `motif_taxonomy.py`, allowing the tool to correctly handle legacy names (e.g., `"curved dna"` → `"Curved_DNA"`, `"g4"` → `"G-Quadruplex"`) from external inputs or older pipeline versions.

---

### 2.3 Motif Detection Algorithms

Each of the nine detector classes inherits from `BaseMotifDetector` (`Detectors/base/base_detector.py`) and implements two required methods: `get_motif_class_name()` and `detect(sequence)`. Scores are normalised to a universal `[1.0, 3.0]` confidence scale by `Utilities/core/motif_normalizer.py`.

#### 2.3.1 G-Quadruplex Detector (`Detectors/gquad/detector.py`)

**Algorithm:** Seeded G4Hunter scoring (Bedrat et al., 2016) extended with regex-based multi-subclass pattern matching.

**Core scoring function:** A JIT-compiled (Numba) sliding-window maximum-sum function (`_compute_max_window_sum_jit`) identifies regions where the G-richness score (G → +1, C → −1, other → 0) over a 25-bp window exceeds the threshold. Eight subclass-specific regex patterns are then applied in priority order.

**Parameters:**
- Window size: 25 bp (default)
- G4Hunter threshold: ≥ 1.2 (canonical) to ≥ 0.5 (weak PQS)
- Score normalization: G4Hunter abs-value scaling (raw 0.5–1.0 → normalised 1.0–3.0, Bedrat et al. 2016)

**Subclass patterns:**
- *Canonical G4*: `(G{3,}[ACGT]{1,7}){3,}G{3,}`
- *Extended-loop canonical*: loop ≤ 12 nt (Chambers et al., 2015)
- *Bulged G4*: G-tract with ≤1 non-G substitution (Mukundan & Bhattacharyya, 2011)
- *Higher-order G4-wire*: ≥3 tandem G4 units
- *Intramolecular G-triplex*: three G-tracts only
- *Two-tetrad weak PQS*: G{2} tracts, G4Hunter ≥ 0.5
- *Stacked G4*: adjacent G4 units with short linker

#### 2.3.2 Z-DNA Detector (`Detectors/zdna/detector.py`)

**Algorithm:** Cumulative 10-mer Z-DNA propensity scoring (Ho et al., 1986) supplemented by eGZ-motif detection (Herbert, 1997).

The detector uses a pre-computed 10-mer score table (`tenmer_table.py`) derived from Ho et al.'s dinucleotide stacking free energies. A sliding-window scan accumulates propensity scores; regions with cumulative score ≥ 50.0 are reported as Z-DNA. Extended GZ (eGZ) dinucleotide repeats (`(GC)_n` or `(CG)_n` with ≥3 repeats) are detected separately via regex.

**Score normalisation:** log-linear scaling (raw 50–2000 → normalised 1.0–3.0).

#### 2.3.3 Curved DNA Detector (`Detectors/curved/detector.py`)

**Algorithm:** A-tract and T-tract phasing detection based on the wedge model of intrinsic DNA curvature (Koo et al., 1986; Olson et al., 1998).

*Local Curvature* is detected by identifying long A-tracts (A{7,}) or T-tracts (T{7,}) that create localised helical bending.

*Global Curvature* is detected by matching phased A-tract (or T-tract) repeat patterns where tracts of length 3–9 are separated by inter-tract spacers consistent with the ~10.5 bp helical repeat (phasing tolerance: 9.9–11.1 bp). Patterns are programmatically generated for tract lengths 3–5 A/T bases in both orientations (`_generate_phased_repeat_patterns`).

**Parameters:** minimum 3 A/T residues per tract, minimum 3 phased tracts for global curvature, phasing centre spacing 11.0 bp ± 1.1 bp. Score normalisation: linear (raw 0.1–0.95 → normalised 1.0–3.0, Koo et al. 1986).

#### 2.3.4 R-Loop Detector (`Detectors/rloop/detector.py`)

**Algorithm:** QmRLFS-finder model (Jenjaroenpun et al., 2015; 2016), which identifies R-loop-forming sequences (RLFS) based on G-richness initiation zones (G-RIZ) and R-loop extension zones (G-REZ).

Detection proceeds as:
1. Scan for G-RIZ windows (≥50% G, length ≥ 20 bp) at step size 100 bp.
2. Extend each G-RIZ into a candidate RLFS until the G-REZ minimum criterion drops below threshold (G-REZ: ≥40% G, maximum length 2,000 bp).
3. Quality filtering: composite quality score ≥ 0.4.

**Score normalisation:** linear (raw 0.4–0.95 → normalised 1.0–3.0).

#### 2.3.5 Cruciform Detector (`Detectors/cruciform/detector.py`)

**Algorithm:** Inverted repeat (IR) detection. The detector identifies palindromic sequences capable of extruding cruciform structures under negative superhelical tension.

IR pairs are identified using reverse complement matching with configurable arm length (≥10 bp arm) and spacer tolerance. Stem stability is scored by arm length × GC-weighted base-pair count.

#### 2.3.6 Triplex Detector (`Detectors/triplex/detector.py`)

**Algorithm:** Mirror repeat detection for intramolecular triplex (H-DNA) forming sequences. Purine-only or pyrimidine-only homopurine·homopyrimidine tracts (purity ≥ 90%) of ≥ 10 bp, capable of forming Hoogsteen or reverse-Hoogsteen base pairs with the adjacent duplex region, are identified.

Sticky DNA (intermolecular triplex formed by two H-DNA-forming sequences in *trans*) is detected as a derived subclass when two proximal mirror-repeat loci are identified.

#### 2.3.7 i-Motif Detector (`Detectors/imotif/detector.py`)

**Algorithm:** C-tract pattern matching. i-Motifs (intercalated cytosine structures) form from C-rich sequences on the strand complementary to G4-forming sequences (Day et al., 2014; Zeraati et al., 2018).

Detection uses regex patterns targeting cytosine runs: canonical i-motif `(C{3,}[ACGT]{1,7}){3,}C{3,}`, relaxed i-motif (shorter C-tracts or longer loops), and AC-motif (alternating AC repeats). Score is derived from C-run count and GC content.

#### 2.3.8 Slipped DNA Detector (`Detectors/slipped/detector.py`)

**Algorithm:** Direct repeat (DR) and short tandem repeat (STR) detection. Slippage mutagenesis arises when polymerase slips on repetitive sequences during replication (Pearson et al., 2005; Wells, 2007).

DRs are detected as exact repeats of a unit (6–50 bp unit) with ≥ 2 copies and ≤ 5 bp spacer. STRs (microsatellites) are detected as tandem repetitions of di-/tri-/tetra-nucleotide units (≥ 3 copies).

#### 2.3.9 A-philic DNA Detector (`Detectors/aphilic/detector.py`)

**Algorithm:** Tetranucleotide log₂-odds scoring. A-philic DNA refers to sequences with high nucleosome positioning potential, characterised by a distinct tetranucleotide composition (Vinogradov, 2003; Rohs et al., 2009).

The detector applies a sliding window over the sequence and computes a propensity score from a pre-compiled 10-mer table (`tenmer_table.py`) of A-phased nucleotide frequencies. Regions exceeding the propensity threshold are annotated as A-philic DNA.

---

### 2.4 Hybrid and Cluster Annotation

Post-detection, NonBDNAFinder applies two additional annotation passes that identify higher-order structural complexity: Hybrid regions and Non-B DNA Clusters. Both are implemented in `Utilities/nonbscanner.py`.

#### 2.4.1 Hybrid Regions

A **Hybrid region** is defined as a genomic locus where two motifs from **different** structural classes overlap spatially by 50–99% (Jaccard overlap). These regions represent sites of structural ambivalence where the local sequence satisfies formation criteria for two distinct non-B DNA conformations simultaneously.

**Algorithm:**
1. For each pair of motifs (class A, class B) where class A ≠ class B, compute the Jaccard overlap:  
   `overlap = |intersection| / |union|`
2. If `0.50 ≤ overlap < 1.00`, annotate as a Hybrid with type label `"ClassA_ClassB_Overlap"`.
3. Hybrids are stored as derived records and excluded from the core motif counts used for density/coverage statistics.

**Parameters:**
- `HYBRID_MIN_OVERLAP = 0.50` (50% minimum spatial overlap)
- `HYBRID_MAX_OVERLAP = 0.99` (excludes identical-coordinate duplicates)

#### 2.4.2 Non-B DNA Clusters

A **Non-B DNA Cluster** is a 300 bp sliding-window region harbouring at least 4 motifs from at least 3 distinct structural classes. These dense hotspots represent the most structurally complex genomic loci and may correspond to regulatory super-enhancers, replication origins, or recombination hotspots (Besnard et al., 2024; Georgakopoulos-Soares et al., 2024).

**Algorithm:**
1. Sort all core motifs by start position.
2. Apply a sliding window of 300 bp.
3. For each window: count distinct motif classes and total motifs within the window.
4. If distinct classes ≥ 3 **and** total motifs ≥ 4, record a Cluster at the window midpoint.
5. Clusters are stored as derived records with metadata: motif count, class count, class list.

**Parameters:**
- `CLUSTER_WINDOW_SIZE = 300` bp
- `CLUSTER_MIN_MOTIFS = 4` motifs
- `CLUSTER_MIN_CLASSES = 3` distinct classes

#### 2.4.3 Biological Interpretation of Hybrid Types

The 30+ distinct Hybrid type-pairs detected across genomes include structurally informative co-occupancy events:

- **Cruciform–G-Quadruplex**: GC-rich palindromic regions supporting both IR extrusion and G4 folding; enriched at oncogene promoters (Kouzine et al., 2017).
- **G-Quadruplex–R-Loop**: G4 formation stabilises negative supercoiling that also promotes R-loop formation; linked to transcriptional pausing (Skourti-Stathaki & Proudfoot, 2014).
- **Z-DNA–G-Quadruplex**: purine-pyrimidine alternation adjacent to G-runs; co-occurs at promoters under torsional stress (Rich & Zhang, 2003).
- **Curved DNA–Cruciform**: palindromic A-tracts generating both intrinsic curvature and cruciform-prone inverted repeats; fragile site determinants (Pearson et al., 1996).

---

### 2.5 Statistical Validation and Benchmarking

NonBDNAFinder was benchmarked against the Non-B DNA Structure Tool (NBST; Buske et al., 2012) and the non-B gfa suite (Cer et al., 2012) using a 40,523 bp test sequence (`NBSTVALIDATION/data/693fc40d26a53.fasta`) for which NBST reference annotations were available for G-Quadruplex, Z-DNA, Direct Repeat, STR, Mirror Repeat, and Curved DNA.

#### 2.5.1 Performance Metrics

For each subclass comparison the following metrics were computed:

| Metric | Formula | Interpretation |
|---|---|---|
| True Positives (TP) | Predicted ∩ Reference (Jaccard ≥ 0.5) | Correct detections |
| False Positives (FP) | Predicted − Reference | Over-predictions |
| False Negatives (FN) | Reference − Predicted | Missed detections |
| Precision | TP / (TP + FP) | Prediction specificity |
| Recall | TP / (TP + FN) | Detection sensitivity |
| F1 Score | 2·Precision·Recall / (Precision + Recall) | Harmonic balance metric |
| Jaccard Index | \|Intersection\| / \|Union\| (coordinate-level) | Coordinate agreement quality |

#### 2.5.2 Matching Criterion

Two motifs are considered a match when their genomic coordinates overlap with a Jaccard Index ≥ 0.5 (50% minimum coordinate overlap). This threshold is standard for genomic interval benchmarking and is consistent with BEDTools intersect default behaviour (Quinlan & Hall, 2010).

#### 2.5.3 Benchmark Results Summary

Key results for high-confidence subclasses:

| Subclass | Precision | Recall | F1 | Jaccard (mean) |
|---|---|---|---|---|
| Canonical intramolecular G4 | 1.00 | 1.00 | 1.00 | 1.00 |
| Direct Repeat | 1.00 | 1.00 | 1.00 | 1.00 |
| STR (microsatellite) | 1.00 | 0.94 | 0.97 | 0.95 |
| Extended-loop canonical G4 | 0.80 | 0.13 | 0.22 | 0.62 |
| Intramolecular G-triplex | 0.75 | 0.09 | 0.16 | 0.58 |
| Z-DNA | 0.67 | 0.14 | 0.23 | 0.54 |

Curved DNA subclasses (Local/Global Curvature) showed zero concordance with NBST curved predictions, reflecting algorithmic divergence: NonBDNAFinder implements A-tract phasing (Koo et al., 1986) whereas NBST uses the DNase I bendability model (Brukner et al., 1995). These represent **complementary annotations** rather than competing predictions.

#### 2.5.4 Multi-Genome Statistical Analysis

For the 9-genome validation suite, the following statistical summaries are reported per genome (tables in `NBSTVALIDATION/tables/`):

- `master_genome_overview.csv` / `master_genome_overview_extended.csv`: genome size, GC%, core motifs, coverage%, density/kb, hybrid count, cluster count, classes detected, mean score, runtime (s)
- `master_class_distribution.csv`: per-genome counts for each of 9 structural classes
- `master_subclass_stats.csv` / `master_subclass_breakdown.csv`: subclass-level counts and fractions
- `master_density_coverage_by_class.csv`: per-class density and coverage per genome
- `master_hybrid_types.csv` / `master_hybrid_types_per_genome.csv`: Hybrid type-pair frequencies
- `master_hybrid_cluster.csv`: combined hybrid + cluster statistics per genome
- `master_cluster_composition.csv` / `master_cluster_complexity.csv`: cluster class-richness distribution

Pearson correlations between GC content and non-B DNA metrics were computed in Python (SciPy `pearsonr`): coverage r = 0.84, density r = 0.11, hybrid count r = 0.83, cluster count r = 0.85.

---

### 2.6 Performance Optimization

NonBDNAFinder implements a multi-layer performance stack to achieve sub-linear scaling with genome size on commodity hardware.

#### 2.6.1 JIT Compilation (Numba)

Performance-critical inner loops in the G-Quadruplex detector (sliding-window maximum-sum) are decorated with `@jit(nopython=True, cache=True)` via the Numba library, providing 2–5× speedup over pure Python on first call and subsequent instant execution from cache. Numba is an optional dependency; the code gracefully falls back to pure Python when Numba is unavailable (e.g., on Streamlit Cloud).

#### 2.6.2 Aho-Corasick Multi-Pattern Matching

The A-philic DNA detector and Z-DNA eGZ detector use the `pyahocorasick` library (Aho & Corasick, 1975) for simultaneous multi-pattern matching. This reduces the asymptotic complexity from O(P×N) (P independent regex scans of length N) to O(N + P + M) (single pass with M matches), yielding 50–200× speedup over sequential regex for large pattern sets.

#### 2.6.3 Interval Tree Overlap Detection

Hybrid annotation uses `intervaltree` (an augmented balanced BST) for interval overlap queries, reducing the Hybrid search from O(N²) naive pairwise comparison to O(N log N + K) (N motifs, K output Hybrid pairs).

#### 2.6.4 Core-Region Deduplication

Chunk-boundary overlap deduplication is implemented as an O(N) reverse scan (`OverlapDeduplicator`) that retains only motifs with start positions within the authoritative core region of each chunk. This avoids costly sort-based deduplication at final assembly and keeps post-processing strictly linear in the number of detected motifs.

#### 2.6.5 Disk-First Architecture for Large Genomes

For sequences ≥ 5 Mbp, the `disk_streaming` strategy writes each chunk's results to temporary disk files via `DiskChunkManager` (`Utilities/disk_chunk_manager.py`) and `UniversalSequenceStorage` (`Utilities/disk_storage.py`). This bounds RAM consumption to a constant amount regardless of genome size (one chunk in memory at a time, ~50 MB for a 50 kb chunk with all 9 detectors). The full result table is only materialised on disk when the user explicitly requests export (`FinalExporter.to_dataframe`).

#### 2.6.6 Cython Extensions (Optional)

`setup.py` provides optional Cython compilation of the Z-DNA 10-mer scanning hot path (`Detectors/zdna/hyperscan_backend.py`), yielding a further 2–5× speedup for Z-DNA detection on production deployments. Cython is entirely optional; the build gracefully falls back to pure Python.

#### 2.6.7 Observed Throughput

On commodity hardware (8-core CPU, 16 GB RAM), NonBDNAFinder achieves:

| Genome size | Strategy | Throughput | Wall-clock time |
|---|---|---|---|
| 174 kbp (*C. ruddii*) | chunk_workers | ~40 kbp/s | 4.4 s |
| 4.4 Mbp (*M. tuberculosis*) | chunk_workers | ~63 kbp/s | 70 s |
| 4.6 Mbp (*E. coli*) | chunk_workers | ~84 kbp/s | 55 s |

The mean throughput across the 9-genome validation suite is approximately 62 kbp/s.

---

### 2.7 Output Generation and Visualization

#### 2.7.1 Export Formats

NonBDNAFinder supports five export formats for detected motifs, implemented in `Utilities/final_exporter.py` and `Utilities/export/`:

| Format | Extension | Contents |
|---|---|---|
| Tab-delimited text | `.tsv` | All motif fields: Sequence_ID, Class, Subclass, Start, End, Length, Score, Strand |
| Comma-separated values | `.csv` | Same as TSV; compatible with Excel and R |
| Browser Extensible Data | `.bed` | BED6 format for genome browser upload |
| JSON | `.json` | Structured motif records for downstream API use |
| Excel Workbook | `.xlsx` | Multi-sheet workbook: motifs, summary statistics, class distribution |

#### 2.7.2 Output Fields

Each motif record contains the following fields:

| Field | Description |
|---|---|
| `Sequence_ID` | FASTA sequence identifier |
| `Class` | Canonical motif class name (e.g., `G-Quadruplex`) |
| `Subclass` | Structural subclass (e.g., `Canonical intramolecular G4`) |
| `Start` | 0-based start coordinate (bp) |
| `End` | 0-based exclusive end coordinate (bp) |
| `Length` | Motif length (bp) |
| `Score` | Normalised confidence score [1.0–3.0] |
| `Strand` | Strand of detection (`+` / `−` / `.`) |
| `GC_Content` | GC% of motif sequence (where applicable) |

#### 2.7.3 Visualization Suite

The `VisualizationPipeline` (`Utilities/visualization_pipeline.py`) renders publication-quality figures from pre-binned summary arrays generated by `VisualizationAccumulator`. Figures are never rendered from the full motif table; all plotting inputs are fixed-size summary arrays (≤ 200 bins), ensuring constant memory usage.

**Available plot types:**

| Plot | Function | Description |
|---|---|---|
| Density histogram | `plot_density_histogram()` | Motif density per 10 kb genomic window |
| Length distribution | `plot_length_histogram()` | Motif length frequency per class |
| Class distribution | `plot_class_distribution()` | Bar chart of motif counts by structural class |
| Co-occurrence heatmap | `plot_cooccurrence_heatmap()` | Matrix of pairwise class co-detection frequency |
| Stacked class distribution | `stacked_bar_class_subclass.py` | Stacked bar of class/subclass per genome |
| Coverage plot | — | Genome-wide coverage track |

**Multi-genome validation figures** (300 DPI, generated by `NBSTVALIDATION/run_genome_validation.py` and `NBSTVALIDATION/extend_report.py`):

| Figure | File | Content |
|---|---|---|
| Fig. 1 | `fig1_coverage_density.png` | Coverage vs genome size; density ranking |
| Fig. 2 | `fig2_class_distribution_stacked.png` | Stacked class distribution per genome |
| Fig. 3 | `fig3_hybrid_cluster.png` | Hybrid and cluster counts per genome |
| Fig. 4 | `fig4_density_heatmap.png` | Per-class density heatmap (log scale) |
| Fig. 5 | `fig5_cooccurrence.png` | Class co-occurrence matrix |
| Fig. 6 | `fig6_genome_size_vs_hybrids_clusters.png` | Genome size vs hybrid/cluster scaling |
| Fig. 7 | `fig7_gc_vs_motif_landscape.png` | GC% vs coverage and density |
| Fig. 8 | `fig8_hybrid_type_breakdown.png` | Hybrid type-pair frequency |
| Fig. 9 | `fig9_subclass_top15.png` | Top-15 subclass frequencies |
| Fig. 10 | `fig10_cluster_composition.png` / `fig10_cluster_complexity.png` | Cluster class-richness distribution |
| Fig. 11 | `fig11_subclass_diversity.png` | Per-genome subclass diversity |

#### 2.7.4 Interactive Streamlit UI

The Streamlit front-end (`app.py`, `UI/`) provides an interactive web application with:

- **Upload tab**: drag-and-drop FASTA upload (single or multi-FASTA), sequence preview and validation.
- **Analysis tab**: configurable motif class and subclass selection (via `motif_taxonomy.py` selector), adaptive chunk size and overlap sliders, real-time progress tracking.
- **Results tab**: interactive data tables, downloadable exports, in-browser visualizations.
- **Documentation tab**: embedded algorithm descriptions, parameter reference, peer-reviewed citations (`UI/documentation.py`).

---

### 2.8 Implementation and Availability

#### 2.8.1 Programming Language and Dependencies

NonBDNAFinder is implemented entirely in **Python 3.8+**. Core dependencies are:

| Package | Version | Purpose |
|---|---|---|
| `streamlit` | ≥ 1.28.0 | Interactive web UI |
| `numpy` | ≥ 1.21.0 | Numerical arrays; sliding window operations |
| `pandas` | ≥ 1.3.0 | Tabular data management and export |
| `scipy` | ≥ 1.7.0 | Statistical testing (Pearson correlation) |
| `matplotlib` | ≥ 3.5.0 | Static figure rendering |
| `seaborn` | ≥ 0.11.0 | Statistical plot styling |
| `plotly` | ≥ 5.17.0 | Interactive visualizations in Streamlit |
| `biopython` | ≥ 1.79 | FASTA parsing and sequence utilities |
| `openpyxl` / `xlsxwriter` | ≥ 3.0.0 | Excel export |
| `psutil` | ≥ 5.8.0 | System resource monitoring |
| `pyahocorasick` | ≥ 2.0.0 | Aho-Corasick multi-pattern matching |
| `intervaltree` | ≥ 3.1.0 | Efficient interval overlap detection |

**Optional performance dependencies:**

| Package | Purpose |
|---|---|
| `numba` | JIT compilation of G4 scoring hot path (2–5× speedup) |
| `hyperscan` | High-performance regex engine (optional Z-DNA acceleration) |
| `Cython` | AOT compilation of Z-DNA 10-mer backend (2–5× speedup) |

#### 2.8.2 Software Architecture Pattern

NonBDNAFinder follows a **modular strategy pattern**:

- All detectors implement a common `BaseMotifDetector` interface, enabling addition of new structural classes without modifying existing code.
- The `MotifDetectionEngine` selects execution strategy at runtime based on input characteristics.
- All class/subclass names are centralised in `motif_taxonomy.py`, preventing label inconsistency.
- Visualization is decoupled from detection: the `VisualizationAccumulator` accumulates fixed-size summary statistics during detection; plotting is deferred to `VisualizationPipeline`.

#### 2.8.3 Availability

- **Source code:** [https://github.com/VRYella/NonBDNAFinder](https://github.com/VRYella/NonBDNAFinder)
- **License:** MIT (open source; unrestricted academic and commercial use)
- **Version:** 2025.1
- **Live web application:** Deployable on Streamlit Cloud; all cloud-incompatible optimisations (Numba, Cython) are optional and gracefully disabled.

#### 2.8.4 Installation

```bash
# Clone repository
git clone https://github.com/VRYella/NonBDNAFinder.git
cd NonBDNAFinder

# Install Python dependencies
pip install -r requirements.txt

# Optional: compile Cython extensions (2-5x speedup)
python setup.py build_ext --inplace

# Launch interactive Streamlit application
streamlit run app.py
```

#### 2.8.5 Citation

> Yella VR (2025). NonBDNAFinder: Comprehensive Detection and Analysis of Non-B DNA Forming Motifs. GitHub: https://github.com/VRYella/NonBDNAFinder

---

### 2.9 Genome Datasets

#### 2.9.1 Multi-Genome Validation Suite

Nine complete bacterial genome sequences were used for validation and benchmarking of NonBDNAFinder. These were selected to span the full range of prokaryotic genome sizes (174 kbp – 4.6 Mbp) and GC content (17.6% – 76.2%), ensuring that performance is characterised across diverse compositional and structural contexts.

**Table 2.9.1 — Validation genome characteristics.**

| Organism | Assembly accession | Size (bp) | GC (%) | Lifestyle | Runtime (s) |
|---|---|---|---|---|---|
| *Candidatus* Carsonella ruddii | — | 174,014 | 17.63 | Obligate endosymbiont (insect) | 4.36–4.58 |
| *Buchnera aphidicola* | — | 452,078 | 18.28 | Obligate endosymbiont (aphid) | 11.51–11.94 |
| *Helicobacter pylori* | — | 1,674,010 | 38.79 | Gastric pathogen | 24.61–25.54 |
| *Streptococcus pneumoniae* | — | 2,110,968 | 39.73 | Respiratory pathogen | 25.42–26.65 |
| *Staphylococcus aureus* | — | 2,821,361 | 32.87 | Gram-positive pathogen | 36.83–37.04 |
| *Miltoncostaea marina* | — | 3,370,274 | 76.16 | Marine actinobacterium | 76.39–83.38 |
| *Cellulomonas shaoxiangyii* | — | 3,909,366 | 75.30 | Soil actinobacterium | 81.60–89.55 |
| *Mycobacterium tuberculosis* H37Rv | — | 4,411,532 | 65.61 | Intracellular pathogen | 65.62–70.10 |
| *Escherichia coli* K-12 | — | 4,641,652 | 50.79 | Model gram-negative | 53.45–55.22 |

Genome files are stored as uncompressed FASTA in `NBSTVALIDATION/` (`.fna` / `.fasta` format). All analyses were performed on complete, unmasked genomic sequences.

#### 2.9.2 Dataset Selection Rationale

The nine genomes were chosen to represent:

1. **AT-rich endosymbionts** (*Candidatus* Carsonella ruddii, *Buchnera aphidicola*; GC ≈ 18%): genomes dominated by Curved DNA, Cruciform, and Slipped DNA motifs driven by high A-tract frequency.
2. **Moderate-GC pathogens** (*H. pylori*, *S. pneumoniae*, *S. aureus*, *E. coli*; GC 33–51%): diverse non-B DNA profiles with representation of all 9 structural classes.
3. **High-GC soil and clinical bacteria** (*Miltoncostaea marina*, *Cellulomonas shaoxiangyii*, *M. tuberculosis*; GC 65–76%): genomes with extremely high G-Quadruplex and Z-DNA density driven by G-run and purine-pyrimidine alternation frequency.

#### 2.9.3 Single-Sequence Benchmark

An additional 40,523 bp test sequence (`NBSTVALIDATION/data/693fc40d26a53.fasta`) was used for subclass-level benchmarking against NBST reference annotations. NBST output files (`.tsv`) for G-Quadruplex, Z-DNA, Direct Repeat, STR, Mirror Repeat, and Curved DNA are provided in `NBSTVALIDATION/data/` for reproducibility.

#### 2.9.4 Genome-Level Results Summary

Across the 9 validation genomes:

| Metric | Value |
|---|---|
| Total genomic content | 23,565,255 bp (23.57 Mbp) |
| Total non-B DNA motifs detected (core) | ~144,060 – 161,872 |
| Mean genomic coverage | 25.93 – 25.96% |
| Mean motif density | 6.92 motifs/kb |
| Most frequent class | G-Quadruplex (78,434–78,483 motifs; ~48.5%) |
| Total Hybrid regions | 5,562 – 5,613 |
| Total Non-B DNA Clusters | 12,140 – 12,159 |
| GC range | 17.6% – 76.2% |
| Total pipeline runtime (9 genomes) | 380 – 404 s (6.3–6.7 min) |

Detailed per-genome and per-class statistics are tabulated in `NBSTVALIDATION/tables/master_genome_overview_extended.csv`.

---

## References

- **Aho AV, Corasick MJ** (1975). Efficient string matching: an aid to bibliographic search. *Commun ACM* 18(6):333–340.
- **Bedrat A, Lacroix L, Mergny JL** (2016). Re-evaluation of G-quadruplex propensity with G4Hunter. *Nucleic Acids Res* 44(4):1746–1759. doi:10.1093/nar/gkw006
- **Besnard E, et al.** (2024). Non-B DNA at replication origins. *Molecular Cell* 84:201–215.
- **Brukner I, Sanchez R, Suck D, Pongor S** (1995). Sequence-dependent bending propensity of DNA. *EMBO J* 14(8):1812–1818.
- **Buske FA, et al.** (2012). Intrinsic DNA curvature: an online tool and its biological implications. *Nucleic Acids Res* 40(W1):W568–W572.
- **Cer RZ, et al.** (2012). Non-B DB v2.0: a database of predicted non-B DNA-forming motifs. *Nucleic Acids Res* 41(D1):D94–D100.
- **Chambers VS, et al.** (2015). High-throughput sequencing of DNA G-quadruplex structures in the human genome. *Nat Biotechnol* 33(8):877–881.
- **Day HA, Pavlou P, Waller ZAE** (2014). i-Motif DNA: structure, stability and targeting with ligands. *Bioorg Med Chem* 22(16):4407–4418.
- **Frank-Kamenetskii MD, Mirkin SM** (1995). Triplex DNA structures. *Annu Rev Biochem* 64:65–95.
- **Georgakopoulos-Soares I, et al.** (2024). Genome-wide profiling of non-B DNA structures and their regulatory roles. *Nature Reviews Genetics*.
- **Herbert A** (1997). Left-handed Z-DNA: structure and function. *Trends Biochem Sci* 22(11):427–433.
- **Ho PS, et al.** (1986). The interactions of ruthenium hexaammine with Z-DNA. *J Biomol Struct Dyn* 4(3):521–534.
- **Huppert JL, Balasubramanian S** (2005). Prevalence of quadruplexes in the human genome. *Nucleic Acids Res* 33(9):2908–2916.
- **Jenjaroenpun P, et al.** (2015). QmRLFS-finder: a model, web server and stand-alone tool for predicting R-loop-forming sequences. *Nucleic Acids Res* 43(W1):W527–W534.
- **Koo HS, Wu HM, Crothers DM** (1986). DNA bending at adenine-thymine tracts. *Nature* 320(6062):501–506.
- **Kouzine F, et al.** (2017). Permanganate/S1 nuclease footprinting reveals non-B DNA structures with regulatory potential across a mammalian genome. *Cell Systems* 4(3):344–356.
- **Lilley DM** (1980). The inverted repeat as a recognisable structural feature in supercoiled DNA molecules. *Proc Natl Acad Sci USA* 77(11):6468–6472.
- **Mukundan VT, Bhattacharyya D** (2011). Thermodynamic and fluorescence studies on bulged G-quadruplexes. *J Am Chem Soc* 133(15):5615–5617.
- **Olson WK, et al.** (1998). New parameters for DNA sequence-dependent curvature. *J Mol Biol* 313:229–237.
- **Pearson CE, Nichol Edamura K, Cleary JD** (2005). Repeat instability: mechanisms of dynamic mutations. *Nat Rev Genet* 6(10):729–742.
- **Pearson CE, et al.** (1996). Inverted repeats, stem-loops, and cruciforms. *J Cell Biochem* 63(1):1–22.
- **Quinlan AR, Hall IM** (2010). BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics* 26(6):841–842.
- **Rich A, Zhang S** (2003). Z-DNA: the long road to biological function. *Nat Rev Genet* 4(8):566–572.
- **Rohs R, et al.** (2009). The role of DNA shape in protein-DNA recognition. *Nature* 461(7268):1248–1253.
- **Skourti-Stathaki K, Proudfoot NJ** (2014). A double-edged sword: R loops as threats to genome integrity and powerful regulators of gene expression. *Genes Dev* 28(13):1384–1396.
- **Trifonov EN, Sussman JL** (1980). The pitch of chromatin DNA is reflected in its nucleotide sequence. *Proc Natl Acad Sci USA* 77(7):3816–3820.
- **Vinogradov AE** (2003). DNA helix: the importance of being GC-rich. *Nucleic Acids Res* 31(7):1838–1844.
- **Wells RD** (2007). Non-B DNA conformations, mutagenesis and disease. *Trends Biochem Sci* 32(6):271–278.
- **Zeraati M, et al.** (2018). I-motif DNA structures are formed in the nuclei of human cells. *Nat Chem* 10(6):631–637.
- **Zhao J, et al.** (2024). Non-B DNA structures: a comprehensive overview of their roles in human disease. *Nucleic Acids Res* 52(1):1–22.
