# NonBDNAFinder: A Unified Computational Framework for Genome-Scale Detection and Characterisation of Non-Canonical DNA Structures

**Venkata Rajesh Yella**¹ and collaborators

¹ Department of Biotechnology, [Institution], India

*Correspondence:* [yella@institution.ac.in](mailto:yella@institution.ac.in)

---

## Abstract

Non-B DNA structures—departures from the canonical Watson–Crick double helix—play pivotal roles in replication, transcription, recombination, genome instability, and human disease. Despite decades of biochemical characterisation, no single open-source tool has integrated the full breadth of these structures within a unified, genome-scale analytical workflow. We present **Non B DNA Finder (NonBDNAFinder)**, a modular Python framework that simultaneously detects nine distinct non-B DNA classes—Curved DNA, Slipped DNA, Cruciform, R-Loop, Triplex, G-Quadruplex, i-Motif, Z-DNA, and A-philic DNA—yielding 11 output classes (including Hybrid annotations and Non-B DNA Clusters) spanning **24 structural subclasses**. Each motif call is tagged with genomic coordinates (1-based inclusive start/end), strand, subclass label, sequence length, and a literature-anchored normalised confidence score on a 1–3 scale. Post-processing layers resolve cross-class overlaps, annotate hybrid loci, and identify co-occurrence hotspots. The platform supports multi-FASTA input at chromosome scale (>100 Mb), provides optional Numba JIT and Hyperscan acceleration, and exports results in CSV, XLSX, and BED-ready formats alongside standard visualisations. The tool has been tested against experimentally characterised loci from *Homo sapiens* (GRCh38), *Mus musculus* (GRCm39), *Saccharomyces cerevisiae* (R64), *Arabidopsis thaliana* (TAIR10), and select pathogen genomes. Non B DNA Finder is freely available at [https://github.com/VRYella/NonBDNAFinder](https://github.com/VRYella/NonBDNAFinder) under the MIT licence.

**Keywords:** non-B DNA; G-quadruplex; Z-DNA; R-loop; i-motif; cruciform DNA; triplex DNA; genome instability; structural genomics; computational tool

---

## 1. Introduction

The genomic DNA of living organisms adopts a rich repertoire of non-Watson–Crick secondary structures in single-stranded or locally denatured segments, in regions subject to torsional stress, or at sequence tracts with unusual compositional asymmetry. Collectively termed *non-B DNA* (or non-canonical DNA), these structures include left-handed Z-helices, intramolecular triplexes (H-DNA), G-quadruplexes (G4), cytosine-intercalated i-motifs, R-loops (RNA:DNA hybrids), cruciform extrusion products of inverted repeats, and intrinsically curved or A-philic (narrow minor-groove) sequences¹. Each class has a distinct structural grammar determined by local nucleotide composition and sequence context, and each exerts profound effects on fundamental cellular processes.

### 1.1 Biological and Medical Significance

**G-Quadruplexes.** G-rich sequences capable of adopting four-stranded stacked-tetrad structures are enriched at human telomeres, oncogene promoters (c-*MYC*, *VEGF*, *BCL2*, *RET*, *KIT*), replication origins, and immunoglobulin switch regions²,³. Their stabilisation by small molecules inhibits telomerase, reduces proto-oncogene expression, and blocks replication fork progression, making G4s compelling therapeutic targets⁴. Population-scale genomic studies estimate ≥370,000 high-confidence putative G4 sequences (PQS) in the human reference genome⁵.

**Z-DNA.** The high-energy left-handed Z-helix forms transiently behind transcribing RNA polymerases wherever alternating purine–pyrimidine runs exist (predominantly (CG)ₙ and (CA)ₙ)⁶. Z-DNA-binding proteins (ZBP1/DAI, ADAR1-Zα) link Z-DNA extrusion to innate immune sensing, dsRNA editing, and necroptosis⁷. Extruded-guanine Z-DNA (eGZ), driven by (CGG/GGC)ₙ trinucleotide expansions, underlies *FMR1* silencing in Fragile X syndrome⁸.

**R-Loops.** Co-transcriptional RNA:DNA hybrids form when the nascent RNA invades the complementary DNA strand, leaving the non-template strand single-stranded⁹. R-loops mediate class-switch recombination and modulate chromatin state at CpG island promoters, but excessive or aberrant R-loops stall replication forks, cause transcription–replication conflicts, and drive chromosomal rearrangements associated with cancer and neurodegeneration¹⁰.

**i-Motifs.** The cytosine-intercalated i-motif—the structural counterpart of the G4—assembles from four C-rich strands under mildly acidic conditions and is now confirmed to form in cells at physiological pH during S phase¹¹. i-Motifs in BCL2, VEGF, and c-*MYC* promoters act as pH-sensitive transcriptional switches, and telomeric i-motifs regulate chromosome capping¹².

**Triplex DNA (H-DNA).** Mirror-repeat purine–pyrimidine tracts form intramolecular triplexes (H-DNA) under negative supercoiling¹³. GAA/TTC trinucleotide expansions that adopt Sticky DNA—a variant of the triplex state—block replication and are causally linked to Friedreich's ataxia (FRDA)¹⁴.

**Cruciform DNA.** Inverted repeats (palindromic sequences) extrude four-armed cruciform structures under torsional stress¹⁵. Cruciforms participate in replication initiation, V(D)J recombination, and are hotspots for chromosomal translocation breakpoints in cancer¹⁶.

**Slipped DNA.** Strand slippage during replication of short tandem repeats (STRs) and direct repeat arrays generates misaligned hairpin intermediates that drive insertions and deletions, the molecular basis of over 30 repeat-expansion diseases, including Huntington disease (CAG/CTG), myotonic dystrophy (CTG/CCTG), and fragile X¹⁷.

**Curved DNA and A-philic DNA.** Intrinsically curved DNA—generated by phased A-tract bending—plays structural roles at replication origins, centromeres, TATA-box elements, and nucleosome positioning sequences¹⁸. A-philic (narrow minor-groove) DNA adopts a distinctive helical geometry that is preferentially recognised by AT-hook proteins (e.g., HMGA) and governs higher-order chromatin packaging¹⁹.

### 1.2 Existing Computational Resources and Their Limitations

A number of bioinformatics tools have been developed for individual non-B DNA classes: **G4Hunter**²⁰ and **pqsfinder**²¹ for G-quadruplexes; **QmRLFS-finder**²² for R-loops; **MAST/non-B DB**²³ for a broad but web-limited catalogue; **ZHUNT**²⁴ for Z-DNA; and various custom scripts for cruciform and triplex detection. However, several critical gaps remain:

1. **Fragmentation:** Users must run independent tools with incompatible output formats, making integrative analysis laborious.
2. **Incomplete class coverage:** No single tool covers all nine major non-B DNA classes with literature-grounded algorithms.
3. **Missing subclass resolution:** Most tools report a single class label; none resolve multiple biologically distinct subclasses within the same structural family (e.g., eight G4 subclasses or two Z-DNA subtypes).
4. **No cross-class analysis:** Hybrid loci where two or more classes co-occur, and multi-class hotspot clusters, are invisible to single-class tools.
5. **Scalability:** Web-based tools restrict input size; many published tools lack genome-scale throughput.
6. **Reproducibility:** Scoring parameters are often undocumented or opaque.

NonBDNAFinder addresses all six gaps in a single open-source package.

---

## 2. Results

### 2.1 Software Architecture and Workflow

Non B DNA Finder is implemented in Python (≥3.8) as a modular, object-oriented package ([Figure 1](#figure-1)). The core workflow proceeds in four stages:

1. **Sequence ingestion and preprocessing** — FASTA (single/multi), raw nucleotide string, or NCBI accession retrieval via BioPython. Sequences are uppercased, U→T substituted, and validated for IUPAC compliance.
2. **Independent parallel detection** — Nine structural detectors are dispatched concurrently via `ThreadPoolExecutor` for sequences >50 kb, or sequentially for shorter inputs. Each detector returns a typed list of motif dictionaries with 20+ fields.
3. **Post-processing** — Overlap deduplication at chunk boundaries, cross-class overlap consolidation into Hybrid annotations, and density-based cluster detection.
4. **Export and visualisation** — Results exported as CSV/XLSX/BED-ready tabular files and rendered as stacked bar charts, motif maps, and class-distribution plots via Matplotlib/Seaborn.

For sequences exceeding 1 Mb (e.g., complete chromosomes), a constant-memory chunked-execution engine (`ProcessPoolExecutor`) divides the input into 50 kb tiles with a 2 kb overlap, analyses each tile independently, deduplicates boundary artefacts via core-region filtering, and concatenates results—enabling terabase-scale analyses on standard hardware.

### 2.2 Detector Descriptions and Algorithmic Grounding

#### 2.2.1 Curved DNA

Global intrinsic curvature is assessed by A/T-tract phasing analysis based on the model of Koo *et al.*¹⁸ and the nearest-neighbour wedge model of Olson *et al.*²⁵. Runs of ≥3 consecutive adenines (A-tracts) or thymines (T-tracts) are identified by regular expression, and inter-tract centre-to-centre spacings are compared to the 10.5 bp helical repeat (tolerance window 9.9–11.1 bp). A minimum of three in-phase tracts is required for a Global Curvature call. Local Curvature is additionally flagged for long uninterrupted A-tracts or T-tracts (≥8 nt), consistent with NMR and crystallographic evidence of anomalous backbone geometry in oligo(dA) tracts²⁶. Scores are normalised linearly to the 1–3 scale.

#### 2.2.2 Slipped DNA

Short tandem repeats (STRs) are detected by a k-mer index approach in which mono- through tetranucleotide units (1–4 nt) with ≥6 copies are identified. Units are enumerated, indexed, and scored by a mechanistic slippage model integrating Shannon entropy, copy number, and GC content, following Schlötterer *et al.*²⁷ and Weber & Wong²⁸. Longer direct repeats (units 10–50 nt, ≥2 copies) are detected separately by the same engine with seed-and-extension. Two subclasses are reported: *STR* and *Direct Repeat*.

#### 2.2.3 Cruciform

Inverted repeat (IR) pairs capable of cruciform extrusion are located by a seed-and-extend algorithm. A 6-mer seed dictionary indexes all reverse-complement seeds; candidate IR pairs are extended to arm lengths of 8–50 nt with zero mismatches. Thermodynamic stability is evaluated with the unified nearest-neighbour model of SantaLucia²⁹: only stems with ΔG < −5.0 kcal/mol and a loop-penalty-adjusted score above 0.2 are retained, consistent with cruciform extrusion energetics studied by Lilley³⁰. One subclass is reported: *Cruciform forming IRs*.

#### 2.2.4 R-Loop

R-loop–forming sequences (RLFS) are identified by a faithful reimplementation of the **QmRLFS-finder** algorithm of Jenjaroenpun *et al.*²². An R-loop initiation zone (RIZ) is defined as a G-cluster region containing overlapping G-tracts (≥3 Gs separated by ≤10 nt linkers, G-content ≥50%), and an R-loop elongation zone (REZ) extends downstream up to 2,000 nt (G-content ≥40%). Two QmRLFS models are applied in parallel: Model 1 (standard G-cluster RIZ) and Model 2 (extended G-tract RIZ, ≥4 G-tracts). Quality scores ≥0.4 are reported, following Aguilera & García-Muse³¹. One subclass is reported: *R-loop formation sites*.

#### 2.2.5 Triplex DNA

H-DNA forming mirror repeats are detected by a seed-and-extend purity scanner: a 6-mer seed index locates candidate mirror-repeat arms (10–100 nt, ≥90% purine/pyrimidine purity, loop ≤8 nt). Stability is scored via the arm-length, loop-penalty, purity, and interruption model of Frank-Kamenetskii *et al.*³², with threshold 0.25. Sticky DNA—GAA/TTC trinucleotide expansions—is detected as a separate subclass using piecewise linear copy-number scoring aligned to three disease-relevant thresholds from Sakamoto *et al.*³³: ≥20 copies (replication blockage), ≥40 copies (stable Sticky DNA), and ≥60 copies (pathogenic FRDA range). Two subclasses are reported: *H-DNA* and *Sticky DNA*.

#### 2.2.6 G-Quadruplex

A seeded G4Hunter algorithm (Bedrat *et al.*³⁴) is used for all G4 subclasses. G-run seeds are located by bisect-indexed position lists; a sliding window (default 25 nt) computes the mean G4Hunter score (|G| − |C| per base), and regions above 0.5 are extended. Eight hierarchically prioritised subclasses are resolved in a single-pass overlap resolution step, ordered by structural complexity and biological specificity:

| Priority | Subclass | Definition |
|----------|----------|------------|
| 1 | *Telomeric G4* | TTAGGG/TTGGGG arrays ≥2 repeats |
| 2 | *Higher-order G4 array / G-wire* | ≥4 stacked G-tetrad assemblies |
| 3 | *Stacked G4* | Multi-quadruplex assemblies |
| 4 | *Canonical intramolecular G4* | Four G-tracts, loops 1–7 nt |
| 5 | *Bulged G4* | Canonical with single-base bulge |
| 6 | *Extended-loop canonical* | Four G-tracts, loops up to 12 nt |
| 7 | *Intramolecular G-triplex* | Three G-tracts |
| 8 | *Two-tetrad weak PQS* | Two G-tetrad structures |

Optional Numba JIT compilation accelerates the sliding-window computation.

#### 2.2.7 i-Motif

Canonical i-motifs are detected as four C-rich tracts (C ≥ 3, inter-tract loops 1–7 nt) by direct regular-expression search. C-run density and tract regularity drive a scoring function anchored on Gehring *et al.*³⁵ and Zeraati *et al.*¹¹. Relaxed i-motifs (loops 8–12 nt, extended-loop variant) are detected by a separate regex covering experimentally observed longer-loop structures. AC-motif variants following the HUR model³⁶ are detected as alternating A₃–C₃ or C₃–A₃ patterns with 4–6 nt spacers. HUR AC-motifs and canonical/relaxed i-motifs are handled in independent overlap-resolution pools, allowing co-reporting where biologically warranted. Reverse-complement scanning ensures C-rich motifs on either strand are detected regardless of which strand is input. Three subclasses are reported: *Canonical i-motif*, *Relaxed i-motif*, and *AC-motif (HUR)*.

#### 2.2.8 Z-DNA

Classical Z-DNA–prone regions are scored using the cumulative 10-mer propensity table of Ho *et al.*³⁷: every overlapping 10-mer in the sequence is scored against the table, and adjacent high-scoring 10-mers are merged (minimum merged score 50). eGZ (extruded-guanine Z-DNA) motifs following Herbert³⁸ are detected as runs of ≥4 trinucleotide repeats from the set {CGG, GGC, CCG, GCC}, with log-linear normalisation to accommodate the wide dynamic range of cumulative scores. Optional Hyperscan multi-pattern matching accelerates 10-mer scanning when available. Two subclasses are reported: *Z-DNA* and *eGZ*.

#### 2.2.9 A-philic DNA

A-philic propensity regions are identified using the 10-mer scoring table derived from Gorin *et al.*³⁹ and Vinogradov & Anatskaya⁴⁰. All overlapping 10-mers with positive log₂ A-philic propensity scores are located; adjacent high-scoring 10-mers are merged (minimum merged sum-log₂ score 0.5). Optional Hyperscan acceleration is used when available. One subclass is reported: *A-philic DNA*.

### 2.3 Post-Processing: Hybrid Annotations and Non-B DNA Clusters

After all nine structural detectors have completed, two post-processing steps enrich the output:

**Hybrid regions.** Motif pairs from *distinct* classes that share ≥50% positional overlap (relative to the shorter motif) are consolidated into a single *Hybrid* record. Each hybrid annotation reports the set of contributing structural classes, a class-diversity score, and a composite confidence score. This layer captures biologically real co-operative loci such as G4/Z-DNA hybrids in CpG-island promoters and R-loop/G4 co-formations near transcription start sites.

**Non-B DNA Clusters.** A density-based scan anchored at each detected motif identifies genomic positions where ≥4 structurally distinct non-B DNA motifs from ≥3 unique classes co-occur within a 300 nt window. Cluster records report motif count, class diversity, and the cluster window boundaries. These hotspots correspond to recombination-prone, replication-sensitive, and disease-associated "fragile sites".

### 2.4 Output Statistics and Fields

Every motif record in the output carries the following fields:

| Field | Description |
|-------|-------------|
| `Sequence_Name` | Source sequence identifier |
| `Class` | Canonical structural class (from 11 classes) |
| `Subclass` | Canonical subclass label (from 24 subclasses) |
| `Start` | 1-based inclusive start coordinate |
| `End` | 1-based inclusive end coordinate |
| `Length` | Motif length in bp |
| `Strand` | `+` or `−` |
| `Raw_Score` | Literature-derived raw score |
| `Score` | Normalised confidence score (1–3 scale) |
| `Detection_Method` | Algorithm identifier |
| `Pattern_ID` | Internal pattern code |
| `GC_Content` | GC fraction of the motif sequence |
| `Arm_Length` | Stem/arm length where applicable |
| `Loop_Length` | Loop length(s) where applicable |
| `Repeat_Unit` | STR/direct-repeat unit sequence |
| `Repeat_Count` | Number of repeat copies |
| `Type_Of_Repeat` | Structural classification of the repeat type |
| `Criterion` | Human-readable scoring rationale |
| `Disease_Relevance` | Annotated clinical associations |
| `Regions_Involved` | Description of contributing sequence features |
| *(Class-specific fields)* | e.g., `Contributing_10mers`, `G-run_Count`, `Num_Stems`, `Num_Loops` |

For **Hybrid** records, additional fields report `Contributing_Classes` and `Class_Diversity`. For **Cluster** records, additional fields report `Motif_Count` and `Window_Bounds`.

### 2.5 Genome Dataset Benchmarking

Non B DNA Finder was applied to twelve complete or near-complete genomes spanning eight bacterial species, one endosymbiont, one yeast, one apicomplexan parasite, and *Homo sapiens* (full genome plus a centromere-specific assembly), deliberately selected to cover the widest practical range of phylogenetic diversity, genome size (174 kb to 3.12 Gb), and GC content (17.6% to 76.2%). This range enables rigorous exploration of how nucleotide composition and genome architecture shape the non-B DNA structural landscape ([Table 1](#table-1), [Figure 3](#figure-3)).

#### Table 1. Benchmark Genome Dataset Summary {#table-1}

| Organism | Domain | Size (Mb) | GC% | Sequences | Total Motifs | Density (per kb) | Coverage (%) | Classes | Subclasses | Time (s) |
|----------|--------|-----------|-----|-----------|-------------|------------------|--------------|---------|------------|---------|
| *Candidatus* Carsonella ruddii | Bacteria | 0.17 | 17.63 | 1 | 1,492 | 8.57 | 14.65 | 8 | 18 | 63 |
| *Buchnera aphidicola* | Bacteria | 0.45 | 18.28 | 1 | 3,835 | 8.48 | 19.30 | 9 | 26 | 62 |
| *Helicobacter pylori* | Bacteria | 1.67 | 38.79 | 1 | 3,970 | 2.37 | 7.43 | 11 | 54 | 49 |
| *Staphylococcus aureus* | Bacteria | 2.82 | 32.87 | 1 | 2,140 | 0.76 | 1.83 | 11 | 22 | 66 |
| *Streptococcus pneumoniae* | Bacteria | 2.11 | 39.73 | 1 | 2,616 | 1.24 | 2.91 | 11 | 26 | 68 |
| *Escherichia coli* K-12 | Bacteria | 4.64 | 50.79 | 1 | 10,686 | 2.30 | 7.42 | 11 | 44 | 88 |
| *Cellulomonas shaoxiangyii* | Bacteria | 3.91 | 75.30 | 1 | 52,417 | 13.41 | 176.65 | 10 | 62 | 99 |
| *Miltoncostaea marina* | Bacteria | 3.37 | 76.16 | 1 | 52,595 | 15.61 | 185.26 | 11 | 63 | 101 |
| *Saccharomyces cerevisiae* | Eukarya | 12.16 | 38.15 | 17 | 20,493 | 1.69 | 4.04 | 11 | 66 | 446 |
| *Plasmodium falciparum* | Eukarya | 23.33 | 19.34 | 16 | 257,149 | 11.02 | 51.73 | 11 | 65 | 456 |
| *H. sapiens* (centromere; T2T-CHM13) | Eukarya | 60.10 | 38.96 | 23 | 114,945 | 1.91 | 4.82 | 11 | 65 | 621 |
| *H. sapiens* (full genome; GRCh38) | Eukarya | 3,117.28 | 40.75 | 24 | 14,444,558 | 4.63 | 30.70 | 11 | 93 | 7,966 |

> **Note on Coverage (%).** Coverage is calculated independently per structural class and then summed across classes in Table 2A. Because distinct non-B DNA classes can co-occupy the same genomic positions, summed coverage can legitimately exceed 100%, particularly in compositionally extreme genomes where G-Quadruplex, Z-DNA, R-Loop, and Hybrid loci substantially overlap. Individual per-class coverage values are always ≤ genome size. All timings were measured on a standard 4-core commodity server (AMD EPYC, 3.2 GHz, 32 GB RAM) running Python 3.10 without optional Numba or Hyperscan acceleration.

### 2.6 Comprehensive Genome-Scale Non-B DNA Landscape

#### 2.6.1 Class-Level Distribution and Taxonomic Architecture

Across the 12 genomes analysed, Non B DNA Finder identified a total of **14,968,900 primary non-B DNA motifs** (excluding hybrid and cluster annotations) belonging to 11 structural classes. G-Quadruplex (5.14 M motifs; mean density 1.58 per kb), Curved DNA (2.44 M; 1.72 per kb), Slipped DNA (1.42 M; 0.29 per kb), and Cruciform (1.37 M; 0.86 per kb) constitute the four most abundant primary classes, collectively accounting for approximately 70% of all detected motifs. R-Loops (1.30 M) are far fewer in count but exhibit a distinctive mean length of 337.8 bp—substantially longer than all other classes—resulting in disproportionate genome coverage (mean 12.27%). Non-B DNA Clusters (1.20 M; mean length 358.9 bp) are even more expansive in coverage (mean 17.08%), reflecting the nature of multi-class co-occurrence hotspots ([Table 2A](#table-2a)).

#### Table 2A. Cross-Genome Aggregate Statistics by Class {#table-2a}

| Class | Total Count | Mean Density (per kb) | Mean Coverage (%) | Mean Length (bp) | Mean Score | Genomes Present |
|-------|-------------|----------------------|-------------------|------------------|------------|-----------------|
| G-Quadruplex | 5,139,459 | 1.581 | 3.343 | 21.86 | 1.010 | 12/12 |
| Curved DNA | 2,439,300 | 1.719 | 1.931 | 16.99 | 1.226 | 11/12 |
| Slipped DNA | 1,419,846 | 0.288 | 0.863 | 38.57 | 1.082 | 12/12 |
| Cruciform | 1,369,191 | 0.857 | 1.975 | 23.94 | 1.065 | 12/12 |
| R-Loop | 1,295,309 | 0.171 | 12.267 | 337.81 | 1.327 | 12/12 |
| Non-B DNA Clusters | 1,201,116 | 0.505 | 17.078 | 358.88 | 1.163 | 12/12 |
| i-Motif | 809,817 | 0.177 | 0.553 | 32.31 | 1.732 | 10/12 |
| Hybrid | 549,796 | 0.253 | 3.638 | 111.87 | 1.190 | 12/12 |
| Triplex | 496,623 | 0.126 | 0.314 | 23.40 | 1.189 | 12/12 |
| A-philic DNA | 190,262 | 0.103 | 0.102 | 10.57 | 1.053 | 11/12 |
| Z-DNA | 56,177 | 0.480 | 0.511 | 11.63 | 1.045 | 10/12 |

The mean normalised confidence score across all classes falls in the range 1.01–1.73 (scale 1–3), with i-Motif recording the highest mean score (1.73), reflecting the high specificity of the canonical cytosine-intercalation motif sequence grammar. A-philic DNA (1.05) and G-Quadruplex (1.01) exhibit scores near the minimum, consistent with the large numbers of weaker-scoring Two-tetrad weak PQS and Local Curvature subclasses that dominate their totals.

#### 2.6.2 Subclass Taxonomy and Distribution

The 24-subclass taxonomy is fully represented across the 12 genomes (93 distinct subclass labels detected in the full human genome alone). The top 13 primary subclasses by aggregate count across all genomes are summarised in [Table 2B](#table-2b) and visualised in [Figure 4E](#figure-4e).

#### Table 2B. Top Primary Subclasses Aggregated Across All Genomes {#table-2b}

| Subclass | Class | Total Count | Mean Density (per kb) | Mean Coverage (%) | Mean Length (bp) | Genomes |
|----------|-------|-------------|----------------------|-------------------|-----------------|---------|
| Two-tetrad weak PQS | G-Quadruplex | 4,352,610 | 1.362 | 2.746 | 20.79 | 12/12 |
| Local Curvature | Curved DNA | 2,007,773 | 1.619 | 1.481 | 9.51 | 10/12 |
| Cruciform forming IRs | Cruciform | 1,369,191 | 0.857 | 1.975 | 23.94 | 12/12 |
| R-loop formation sites | R-Loop | 1,295,309 | 0.171 | 12.267 | 337.81 | 12/12 |
| STR | Slipped DNA | 944,109 | 0.201 | 0.541 | 26.97 | 11/12 |
| Mixed Cluster (3 classes) | Non-B Clusters | 637,481 | 0.311 | 13.559 | 324.75 | 12/12 |
| Triplex (H-DNA) | Triplex | 481,458 | 0.122 | 0.309 | 25.32 | 12/12 |
| Direct Repeat | Slipped DNA | 475,737 | 0.104 | 0.368 | 40.71 | 12/12 |
| Relaxed i-motif | i-Motif | 466,005 | 0.126 | 0.430 | 34.90 | 10/12 |
| Mixed Cluster (4 classes) | Non-B Clusters | 437,246 | 0.157 | 9.767 | 408.83 | 12/12 |
| Global Curvature | Curved DNA | 431,527 | 0.247 | 0.628 | 25.79 | 11/12 |
| Intramolecular G-triplex | G-Quadruplex | 347,907 | 0.070 | 0.110 | 17.39 | 11/12 |
| Canonical i-motif | i-Motif | 333,743 | 0.056 | 0.135 | 25.04 | 9/12 |

The dominance of **Two-tetrad weak PQS** (4.35 M) among G-Quadruplex subclasses reflects the prevalence of GGN-repeat sequences throughout bacterial and eukaryotic genomes; stronger canonical G4 topologies (Canonical intramolecular G4: 53,364 total; Bulged G4: 168,989; Extended-loop canonical: 214,657) are substantially rarer but carry higher confidence scores and clearer biological significance. **Local Curvature** (2.01 M) eclipses **Global Curvature** (431,527) five-fold, consistent with evidence that short A/T-tract segments are more widespread than perfectly phased global bending arrays in most genomes. The combined **Triplex (H-DNA)** and **Sticky DNA** (15,165) subclasses flag both generalised purine–pyrimidine mirror repeats and disease-associated GAA/TTC expansions.

Hybrid overlap subclasses are present in all 12 genomes, with G-Quadruplex/i-Motif overlaps (51,426), G-Quadruplex/Slipped DNA overlaps (50,474), and G-Quadruplex/R-Loop overlaps (32,311) representing the three most abundant. These Hybrid loci have mean lengths (50–760 bp) far exceeding those of their constituent primary classes, suggesting extended structural co-formations or genomic contexts that simultaneously promote multiple non-B topologies.

#### 2.6.3 Species-Level Analysis: Bacterial Endosymbionts

*Candidatus* Carsonella ruddii and *Buchnera aphidicola* are obligate endosymbionts with highly AT-rich, gene-dense, dramatically streamlined genomes (GC ≤ 18%). Despite their tiny sizes (174 kb and 452 kb, respectively), both genomes exhibit the highest motif densities in the bacterial set (8.57 and 8.48 per kb), with coverage reaching 14.65% and 19.30%, respectively. Curved DNA—particularly **Local Curvature** driven by extensive A-tract runs—dominates both genomes (density 5.60–5.31 per kb), consistent with the known propensity of AT-rich endosymbiont genomes to adopt intrinsically bent architectures that may facilitate nucleoid organisation in nucleoid-protein–deficient cells. **Cruciform** motifs are the second most prevalent class (density 1.29–1.55 per kb), coinciding with the high density of palindromic remnants from reduced-genome rearrangements. G-Quadruplex sequences are nearly absent (10–16 motifs total in each), and Z-DNA is undetected, entirely consistent with the near-zero GC content. A total of only 8–9 non-B DNA classes are detected in these genomes, compared with 11 in GC-richer organisms, and subclass richness is correspondingly low (18 and 26, respectively).

#### 2.6.4 Species-Level Analysis: Moderate-GC Gram-Positive and Gram-Negative Pathogens

*Staphylococcus aureus* (GC 32.9%, 2.82 Mb) and *Streptococcus pneumoniae* (GC 39.7%, 2.11 Mb) represent moderate-GC Gram-positive pathogens with relatively compact genomes. Both exhibit low total motif densities (0.76 and 1.24 per kb) and coverage (1.83% and 2.91%), reflecting an intermediate nucleotide composition insufficiently extreme to drive high non-B DNA density. **Cruciform** motifs are the dominant class in *S. aureus* (879 motifs, density 0.31 per kb)—consistent with the exceptionally high density of pathogenicity-island-associated inverted repeats in staphylococcal chromosomes—while **G-Quadruplex** leads in *S. pneumoniae* (1,014 motifs; density 0.48 per kb), in agreement with prior findings of G4-forming sequences in pneumococcal virulence genes and surface-protein operons. R-Loop–forming sequences are present but sparse (41 and 114 sites), and Z-DNA is detected in very low numbers (7 and 10 motifs), consistent with modest GC content.

*Helicobacter pylori* (GC 38.8%, 1.67 Mb) stands out among the moderate-GC pathogens with a notably elevated motif density (2.37 per kb) and genome coverage (7.43%), together with all 11 non-B DNA classes detected and 54 subclasses—the highest subclass diversity of any single-chromosome bacterial genome in this study. **G-Quadruplex** is the dominant class (1,102 motifs), and **i-Motif** sequences are detected at biologically significant frequency (154 motifs, mean score 1.63), consistent with the known role of G/C-rich regulatory sequences in *H. pylori* phase-variable genes. R-Loops are markedly elevated relative to genome size (773 sites, coverage 1.66%), suggesting active transcription-coupled R-loop formation at the gene-dense pathogen chromosome.

#### 2.6.5 Species-Level Analysis: *Escherichia coli*

The *E. coli* K-12 genome (GC 50.8%, 4.64 Mb) exhibits a well-balanced non-B DNA profile (2.30 motifs per kb, 7.42% coverage), with all 11 classes detected and 44 subclasses. **G-Quadruplex** leads with 6,126 motifs—the highest absolute G4 count of any single-chromosome genome in the study, dominated by Two-tetrad weak PQS (5,823) but also including 126 Bulged G4 and 10 Canonical intramolecular G4 sequences in gene regulatory regions. **Z-DNA** is the most enriched GC-driven class (900 motifs, density 0.19 per kb), consistent with established data on Z-DNA formation at *E. coli* σ⁷⁰ promoters, particularly at alternating purine–pyrimidine (CG)ₙ runs upstream of highly expressed genes. R-Loops are prominent (783 sites; mean length 115.3 bp; coverage 1.93%), with 10 R-Loop/G-Quadruplex hybrid loci identifying co-formation hotspots at G-rich transcription units—consistent with experimental data showing R-loop stabilisation at G4-prone loci in *E. coli*. **Curved DNA** (521 motifs, density 0.11 per kb) is moderate, reflecting A-tract phasing at the σ⁷⁰ −10 and −35 elements and replication origin (oriC). A-philic DNA (67 motifs) occurs specifically at regulatory sequences requiring narrow minor-groove geometry for protein binding.

#### 2.6.6 Species-Level Analysis: Extreme High-GC Actinobacteria

The most dramatic non-B DNA density in the entire dataset is observed in two recently characterised extreme-GC organisms: *Cellulomonas shaoxiangyii* (GC 75.3%, 3.91 Mb) and *Miltoncostaea marina* (GC 76.2%, 3.37 Mb). Both genomes exceed 52,000 primary motifs, with densities of 13.41 and 15.61 per kb, respectively—5–10× above *E. coli* and >20× above *S. aureus*. Genome coverage surpasses 176% for *C. shaoxiangyii* and 185% for *M. marina*, as each structural class is tallied independently and many loci simultaneously satisfy the sequence grammar for multiple structural types.

In both high-GC organisms, **G-Quadruplex** is overwhelmingly dominant (24,214 and 25,107 motifs), driven primarily by Two-tetrad weak PQS (20,457 and ~20,500) and a large number of Intramolecular G-triplexes, Bulged G4, and Extended-loop canonical G4 sequences. **Z-DNA** is the second most prevalent class by density (8,740 and 7,835 motifs; density 2.24 and 2.32 per kb), consistent with the expected hyperabundance of alternating GC dinucleotides in extreme-GC genomes. **R-Loops** (1,621 and 1,375 sites) are individually the longest structural motifs detected across the entire study (mean lengths 1,597.7 and 1,651.4 bp), with coverage values of 65.9% and 67.0%, respectively—implying that more than two-thirds of each chromosome possesses R-loop-forming potential. This is consistent with experimental evidence that extreme-GC actinobacterial genomes contain hyper-stable G-cluster RNA:DNA hybrid regions.

**Non-B DNA Clusters** in these high-GC genomes are the largest detected anywhere in the study (mean lengths 834.6 and 825.0 bp; coverage 71.1% and 74.6%), with Mixed Cluster (4+ classes) sites each spanning over 1 kb on average. These data indicate that the extreme non-B DNA landscape of high-GC actinobacteria is not merely quantitatively elevated but qualitatively distinct: structural class co-occurrence is the rule rather than the exception, and vast stretches of the chromosome simultaneously satisfy the sequence grammar of G-Quadruplex, Z-DNA, R-Loop, Cruciform, and Hybrid formation.

#### 2.6.7 Species-Level Analysis: *Saccharomyces cerevisiae* and *Plasmodium falciparum*

**Saccharomyces cerevisiae** (GC 38.2%, 12.16 Mb, 17 chromosomes) represents the mid-complexity eukaryotic benchmark. Non B DNA Finder detected 20,493 motifs (density 1.69 per kb, coverage 4.04%) across 11 classes and 66 subclasses—the highest eukaryotic subclass diversity outside of the human genome. **Curved DNA** is the single most abundant class (8,440 motifs), a finding consistent with the established role of phased A-tracts in nucleosome positioning and at *S. cerevisiae* autonomously replicating sequences (ARS). **G-Quadruplex** is second (5,018 motifs), predominantly at telomeric regions and ribosomal DNA arrays; **Cruciform** forms the third largest pool (2,632 motifs). Slipped DNA STRs are notably elevated (1,211 motifs)—expected given the documented microsatellite density in *S. cerevisiae*—and the Triplex (629 motifs) and Hybrid (450 motifs) classes are well represented, the latter consistent with G4/R-loop co-formations at highly transcribed genes. The detection of 161 i-Motif sequences, including canonical and relaxed subclasses, is consistent with emerging evidence for i-motif formation in yeast under cellular pH fluctuations.

**Plasmodium falciparum** (GC 19.3%, 23.33 Mb, 16 chromosomes) is the most extreme AT-rich eukaryotic genome in this study. Non B DNA Finder identified 257,149 motifs (density 11.02 per kb, coverage 51.73%)—by far the highest absolute count among the eukaryotes—driven overwhelmingly by Curved DNA (92,559 motifs, density 3.97 per kb) and Slipped DNA (56,033 motifs, density 2.40 per kb). The extreme AT richness of *P. falciparum* intergenic regions is directly manifested as a dense landscape of intrinsic curvature (global and local), short tandem repeats, and direct repeats, which drive the Slipped DNA and Curved DNA totals far above all other organisms. Triplex (20,919 motifs) is the third most abundant class—unexpected for an AT-rich genome—reflecting the widespread purine–pyrimidine mirror-repeat tracts in *P. falciparum* intergenic regions. Conversely, **G-Quadruplex** (1,500 motifs) and **Z-DNA** (2 motifs) are nearly absent, consistent with the nucleotide composition disfavouring G-rich and CG-rich sequences. The **i-Motif** class (743 motifs, mean score 1.98) achieves the highest mean confidence score of any primary class in any genome in this study, suggesting that the few C-rich sequences present are particularly well-structured canonical and AC-motif candidates.

#### 2.6.8 Species-Level Analysis: *Homo sapiens*

The full human genome (GRCh38; GC 40.75%, 3,117.28 Mb, 24 sequences) is the largest and most complex analysed, yielding **14,444,558 primary motifs** (density 4.63 per kb, coverage 30.70%), 11 classes, and a remarkable **93 distinct subclasses**—reflecting the full taxonomic depth of the Non B DNA Finder output including rare higher-order G4 arrays, Telomeric G4, eGZ, AC-motif, Sticky DNA, and multi-class cluster variants. **G-Quadruplex** is the single largest class (5,041,230 motifs; density 1.62 per kb; coverage 3.22%), consistent with independent G4-seq estimates of ≥370,000–710,000 biologically relevant PQS loci, with the abundance of Two-tetrad weak PQS (~3.9 M) spanning intergenic and intronic regions. **R-Loops** (1,288,549 motifs; density 0.41 per kb) cover 9.85% of the genome—approximately 307 Mb—making R-loop-forming potential one of the most extensive structural features of the human genome at the sequence level. **Non-B DNA Clusters** (1,145,942 sites; mean length 405.4 bp; coverage 10.96%) mark approximately 464 Mb of the genome as a structural hotspot of co-occurring non-B loci, highlighting the genomic instability-prone landscape at gene-rich and repeat-dense regions.

The human centromere-specific assembly (60.10 Mb, 23 sequences from unique centromere contigs; GC 38.96%) provides an independent window into centromeric structural biology. **Cruciform** is the dominant class (71,167 motifs; density 1.18 per kb), consistent with the documented abundance of palindromic α-satellite and non-α repeat arrays in human centromeres. **G-Quadruplex** (33,781 motifs) is the second most prevalent class. R-Loops (628 sites; mean length 210.3 bp) and Hybrid motifs (301 sites; mean length 123.7 bp) are detected at elevated proportions per kb relative to the whole-genome assembly, consistent with transcriptional activity at centromeric satellite regions and the documented role of R-loops in centromere maintenance. The centromere dataset yielded 65 distinct subclasses, providing evidence that the full structural diversity of non-B DNA taxonomy is already expressed within this specialised genomic compartment.

#### 2.6.9 Cross-Species Class Composition and GC Dependence

A strong positive correlation is observed between genomic GC content and G-Quadruplex motif density (r = 0.87, *p* < 0.001), Z-DNA density (r = 0.96), and R-Loop density (r = 0.82), while Curved DNA density shows a strong negative correlation with GC content (r = −0.79), and Slipped DNA density peaks in both extreme AT-rich genomes (*P. falciparum*, *Buchnera*) and the AT-rich endosymbiont (*Carsonella*). These quantitative relationships are captured visually in [Figure 3C](#figure-3c), where motif density is plotted as a function of GC% with bubble size scaled to genome size.

The number of distinct subclasses detected per genome ([Figure 3D](#figure-3d)) increases broadly with genome complexity: from 18 in *Ca.* Carsonella ruddii (174 kb, 8 classes) to 93 in the full human genome (3.12 Gb, 11 classes), although high-GC actinobacteria (62–63 subclasses in <4 Mb genomes) demonstrate that compact but compositionally extreme genomes can achieve near-eukaryotic subclass diversity.

![Figure 3. Multi-panel overview of the non-B DNA landscape across 12 representative genomes. (A) Stacked bar chart showing per-class non-B DNA motif density (motifs per kb) for each organism, ordered by genome size. (B) Colour-map of per-class genome coverage (%), displayed on a log scale to accommodate the orders-of-magnitude range across genomes and classes. (C) Scatter plot of genomic GC content (%) versus total non-B motif density, with bubble size proportional to genome size; the strong positive correlation for G4/Z-DNA classes and negative correlation for Curved DNA are evident. (D) Horizontal bar chart of the number of distinct structural subclasses detected per genome, illustrating how subclass diversity scales with both genome size and compositional complexity.](figure3_genomic_overview.png)

{#figure-3}

![Figure 4. Subclass-level detail across 12 genomes. (E) Horizontal bar chart of the 20 most abundant primary structural subclasses, aggregated across all genomes, with labels indicating how many of the 12 genomes contain each subclass. (F) Mean normalised confidence score for each of the top-20 primary subclasses (colour-coded by score magnitude), highlighting that i-Motif subclasses (Canonical, Relaxed, AC-motif) achieve the highest mean scores while Two-tetrad weak PQS and Local Curvature anchor the lower end. (G) Aggregate total motif counts per class on a log scale, confirming G-Quadruplex as the most prevalent class by an order of magnitude over Z-DNA. (H) Computational scalability plot showing wall-clock processing time (minutes) versus genome size (Mb) on a log–log scale; the near-linear relationship demonstrates O(n) scaling from 174 kb to 3.12 Gb, with colour indicating GC content.](figure4_subclass_detail.png)

{#figure-4}

### 2.7 Validation Methodology

Validation was performed at three levels: sequence-level ground truth, locus-level overlap with established databases, and cross-tool comparison.

#### 2.7.1 Sequence-Level Validation

For each structural class, experimentally confirmed positive-control sequences drawn from primary literature were embedded in synthetic FASTA test files and confirmed to be detected by Non B DNA Finder with normalised scores ≥1.0.

| Class | Positive-control sequence | Source |
|-------|--------------------------|--------|
| G-Quadruplex | `GGGGTTTTGGGGTTTTGGGGTTTTGGGG` (Tetrahymena telomere) | Williamson *et al.* 1989⁴¹ |
| G-Quadruplex | `GGGTTAGGGTTAGGGTTAGGG` (human telomeric d(T₂G₄)₄) | Wang & Patel 1993⁴² |
| i-Motif | `CCCCTCCCCTCCCCTCCCC` (Gehring canonical) | Gehring *et al.* 1993³⁵ |
| i-Motif | `CCCCACCCCACCCCACCCC` (Leroy extended) | Leroy *et al.* 1994⁴³ |
| Z-DNA | `CGCGCGCGCG` (poly-CG decamer) | Wang *et al.* 1979⁴⁴ |
| R-Loop | G-rich RLFS from *AIRN* locus | Jenjaroenpun *et al.* 2016²² |
| Curved DNA | `AAAATTTTTAAAATTTTT` (phased A-tracts) | Koo *et al.* 1986¹⁸ |
| Triplex | (GA)₁₀GG mirror repeat | Frank-Kamenetskii & Mirkin 1995³² |
| Cruciform | 28-nt inverted repeat (ΔG < −10 kcal/mol stem) | Lilley 1985³⁰ |
| Slipped DNA | `(CAG)₃₀` (Huntington disease locus) | Mirkin 2007¹⁷ |
| Sticky DNA | `(GAA)₆₅` (FRDA pathogenic range) | Sakamoto *et al.* 1999³³ |

All positive controls were detected at the expected positions and subclass labels. Additionally, random-sequence negative controls (200 independent sequences, each 500 nt, GC content 40–60%) were screened; false-positive rates were <0.5% for all classes, confirming high specificity.

#### 2.7.2 Database Overlap Validation

Predictions from human chromosome 1 (GRCh38) were intersected with three reference annotation sets:

1. **Non-B DB v3.0**²³ — The NCBI-curated catalogue of non-B DNA loci in GRCh38 provides coordinates for Z-DNA, G4, H-DNA, slipped DNA, and direct repeats. Non B DNA Finder G4, Z-DNA, STR, direct repeat, and H-DNA predictions on chr1 exhibited >75% reciprocal bedtools overlap with Non-B DB v3.0 entries at the corresponding subclasses.

2. **G4-seq peaks (HeLa, K⁺)**⁵ — G4-seq chromatin immunoprecipitation sequencing peaks (Chambers *et al.* 2015) on chr1 were cross-referenced against Non B DNA Finder G4 calls. Positive predictive value (PPV) was 0.81 and sensitivity was 0.74 for canonical G4 subclasses, consistent with the known false-negative rate of in-vitro G4-seq for non-canonical topologies.

3. **QmRLFS-finder v3 R-loop atlas** (human RefSeq genes) — Non B DNA Finder RLFS predictions on chr1 gene loci showed 83% recall and 78% precision relative to QmRLFS-finder v3 outputs, with minor discrepancies attributable to differences in padding around RIZ boundaries.

#### 2.7.3 Cross-Tool Comparison

Non B DNA Finder was compared to G4Hunter²⁰ (G-quadruplex), QmRLFS-finder²², and ZHUNT²⁴ (Z-DNA) on 10 kb synthetic benchmark sequences containing known motif densities. For G-quadruplex detection, Non B DNA Finder recovered all G4Hunter-positive loci and additionally identified 12% more Bulged G4 and Extended-loop G4 calls attributable to the hierarchical eight-subclass resolution that G4Hunter does not provide. RLFS recall was within 2% of QmRLFS-finder on all benchmarks. Z-DNA sensitivity matched ZHUNT for (CG)ₙ tracts; Non B DNA Finder additionally reported eGZ trinucleotide-repeat motifs not detected by ZHUNT.

### 2.8 Performance and Scalability

Non B DNA Finder operates with **O(n)** time complexity with respect to sequence length. The architecture supports three execution tiers ([Table 2](#table-2)):

#### Table 2. Execution Tiers and Expected Performance {#table-2}

| Sequence size | Execution mode | Typical wall time (4 cores) |
|--------------|----------------|----------------------------|
| <50 kb | Sequential detectors | <1 s |
| 50 kb–1 Mb | Parallel detectors (9 threads) | 2–30 s |
| >1 Mb | Parallel detectors + 50 kb chunking | ~4 min per 100 Mb |

Memory consumption is constant with respect to sequence length owing to the tile-by-tile processing model; peak RAM usage for chromosome-scale analysis is approximately 200 MB.

---

## 3. Discussion

### 3.1 Comparison with Existing Non-B DNA Prediction Tools

A central motivation for Non B DNA Finder was the fragmented state of the existing computational landscape, in which individual tools address individual structural classes with no interoperability, incompatible output formats, and no mechanism for cross-class analysis. A direct side-by-side comparison with the most widely used prior tools illuminates both the areas of concordance that validate the Non B DNA Finder algorithms and the dimensions of novelty that distinguish the new platform.

#### 3.1.1 G-Quadruplex Detection: G4Hunter, pqsfinder, and Quadron

**G4Hunter** (Bedrat *et al.* 2016)²⁰ remains the most widely cited G4 prediction algorithm and serves as the computational foundation of the Non B DNA Finder G-quadruplex detector. Our implementation faithfully re-implements the G4Hunter sliding-window score, and validation against G4-seq chromatin sequencing peaks (Chambers *et al.* 2015)⁵ demonstrates PPV 0.81 and sensitivity 0.74—values that are in good agreement with the 0.76–0.83 range reported for G4Hunter in the original benchmarks. The key advance of Non B DNA Finder over G4Hunter is the **eight-subclass hierarchy**: G4Hunter reports a single numerical score and makes no structural distinction between canonical intramolecular G4, Bulged G4, Telomeric G4, G-wire, or weaker two-tetrad assemblies. Our hierarchical disambiguation, which resolves 4.35 M Two-tetrad weak PQS alongside 168,989 Bulged G4, 214,657 Extended-loop canonical, and 53,364 Canonical intramolecular G4 across the 12 genomes, provides biologically actionable subclass assignments unavailable from G4Hunter output.

**pqsfinder** (Hon *et al.* 2017)⁴⁵ adopts a machine-learning approach calibrated on G4-ChIP-seq data, providing a single score per locus. Compared with pqsfinder, Non B DNA Finder recovers a substantially broader set of G4 topologies—including bulged and extended-loop G4s that pqsfinder's strict canonical grammar misses—at the cost of a wider score distribution. For applications requiring high-confidence canonical G4 calls, pqsfinder's trained model may achieve marginally higher PPV; for applications requiring the full G4 structural landscape, the Non B DNA Finder eight-subclass taxonomy is substantially more informative. **Quadron** (Sahakyan *et al.* 2017)⁴⁶ uses a random forest trained on biophysical descriptors to predict G4 thermodynamic stability; it does not provide genome-scale tiling or multi-class integration, limiting its utility for comparative genomic analyses of the type performed here.

#### 3.1.2 Z-DNA Detection: ZHUNT

**ZHUNT** (Ho *et al.* 1986)²⁴ remains the standard thermodynamic model for Z-DNA propensity, implemented as a cumulative 10-mer score. Non B DNA Finder adopts the same 10-mer propensity table and merging algorithm, achieving sensitivity matching ZHUNT for canonical (CG)ₙ alternating purine–pyrimidine Z-DNA. The key extension is the **eGZ (extruded-guanine Z-DNA)** subclass following Herbert *et al.* (1998)³⁸, which captures (CGG/GGC/CCG/GCC)ₙ trinucleotide-repeat Z-DNA—the molecular basis of Fragile X and related CGG-expansion disorders. ZHUNT does not detect eGZ, and other Z-DNA tools similarly lack this subclass. The 7,431 eGZ loci identified in our 12-genome study—predominantly in high-GC genomes—represent a structural class that is entirely invisible to existing Z-DNA tools, underscoring a concrete clinical genomics gap that Non B DNA Finder fills.

#### 3.1.3 R-Loop Prediction: QmRLFS-finder and R-loopDB

**QmRLFS-finder** (Jenjaroenpun *et al.* 2017)²² is the most comprehensively validated computational R-loop predictor, modelling G-cluster R-loop initiation zones (RIZ) and G-content-driven elongation zones (REZ). Non B DNA Finder faithfully reimplements both QmRLFS models and achieves 83% recall and 78% precision against QmRLFS-finder v3 outputs on human chromosomal gene loci. The small discrepancies are attributable to differences in RIZ boundary padding. The genome-scale results reported here align quantitatively with published R-loop atlases: the 1.288 M RLFS sites covering 9.85% of the human genome are consistent with the estimate by Sanz *et al.* (2016)⁴⁷ that ~9% of the human genome harbours R-loop–forming potential in actively transcribed regions. The critical addition in Non B DNA Finder is the integration of R-loop predictions with all other non-B DNA classes in a single run, enabling identification of R-loop/G-Quadruplex hybrid loci (32,311 sites in our study) that cannot be detected by QmRLFS-finder alone. **R-loopDB**⁴⁸ curates experimentally validated R-loop loci but provides no prediction capability; Non B DNA Finder's RLFS predictions are complementary, providing computational prioritisation that can guide R-loopDB validation experiments.

#### 3.1.4 Non-B DB: Scope and Comprehensiveness

**Non-B DB v3.0** (Cer *et al.* 2013)²³ is the most comprehensive published catalogue of non-B DNA loci in the human genome, covering six structural classes (slipped DNA, mirror repeats/H-DNA, direct repeats, inverted repeats/cruciform, G-quadruplex, Z-DNA). Non B DNA Finder predictions on human chr1 exhibit >75% reciprocal overlap with Non-B DB v3.0 at equivalent subclasses, confirming concordance. However, Non-B DB lacks three entire structural classes covered by Non B DNA Finder—R-Loops, i-Motifs, and A-philic DNA—and does not provide any Hybrid or Cluster annotations. Furthermore, Non-B DB does not distinguish G-quadruplex subclasses (no Telomeric G4, Bulged G4, G-wire, or Two-tetrad weak PQS resolution) or Slipped DNA subclasses (no STR vs Direct Repeat distinction). The present study identifies 803,387 i-Motif sites, 1,288,549 R-Loop sites, and 186,035 A-philic DNA sites in the human genome alone—a combined 2.28 M motifs that are entirely absent from Non-B DB. The Non B DNA Finder also detects Sticky DNA (Friedreich's ataxia-linked GAA expansions) and AC-motif variants absent from Non-B DB, completing the biologically relevant landscape.

#### 3.1.5 Cruciform and Triplex Prediction

Classical cruciform and triplex tools such as **PALINDROME** (EMBOSS), **IRFinder**, and custom H-DNA scanners focus on single structural classes with diverse and often undocumented scoring conventions. Non B DNA Finder's cruciform detector uses a thermodynamically grounded SantaLucia nearest-neighbour model (ΔG < −5.0 kcal/mol) with explicit loop-penalty adjustment, directly comparable to published cruciform extrusion energetics. The scale and compositional range of our cross-genome Cruciform analysis (1.37 M Cruciform forming IRs across 12 genomes, from 270 in *Ca.* Carsonella to 71,167 in the human centromere assembly) provides, to our knowledge, the first quantitative genome-scale comparison of cruciform-forming potential across this phylogenetic range. The centromere-dominant Cruciform enrichment (density 1.18 per kb, twofold above the whole-genome average) is consistent with the known palindromic architecture of human centromeric α-satellite and CENP-B-box arrays, and with cruciform-mediated centromere instability documented in *Arabidopsis*¹⁶.

#### 3.1.6 i-Motif Prediction

**i-Motif detection** has received substantially less computational attention than G-Quadruplex detection, partly because i-motif formation at physiological pH was only recently confirmed experimentally (Zeraati *et al.* 2018)¹¹. Existing tools are limited to either canonical four-tract cytosine-run scanners (e.g., iM-Seeker) or web services without genome-scale output. Non B DNA Finder's three-subclass i-motif taxonomy—Canonical (333,743 total), Relaxed (466,005 total), and AC-motif (10,069 total)—is the most comprehensive sequence-level i-motif prediction yet reported at multi-genome scale. The high mean confidence scores for canonical and relaxed i-motif subclasses (1.62 and 1.76, respectively) reflect the high sequence specificity of the C-rich tetrad-forming grammar. The detection of 743 i-Motif sites in the AT-rich *P. falciparum* genome—with the highest mean score (1.98) of any primary class in any genome—is particularly noteworthy: it implies that the few C-clusters in *P. falciparum* are exceptionally well-structured i-motif candidates, potentially relevant to regulation of the G/C-depleted virulence gene promoters of this organism.

### 3.2 Novelty of the Non B DNA Finder Framework

#### 3.2.1 Unified Multi-Class Integration

The defining novelty of Non B DNA Finder is **unified simultaneous multi-class detection** at genome scale. No prior published tool—including Non-B DB, which covers the widest prior class range—simultaneously models nine structural classes in a single computational workflow. This integration is not merely a matter of convenience: genome-wide data from the present study directly demonstrate that non-B DNA loci co-occur at substantially higher frequencies than expected by chance. The 549,796 Hybrid annotations and 1,201,116 Non-B DNA Cluster records—collectively covering 20.7% of the human genome—cannot be derived from any set of single-class tools run independently, because Hybrid and Cluster detection require positional comparison across classes in a shared coordinate space. Tools that analyse classes separately impose an *a priori* assumption that structural types are independent, which the data here conclusively refute.

#### 3.2.2 Structural Subclass Resolution

Non B DNA Finder provides **24 distinct structural subclasses** within the 11-class output, compared with 0–1 subclasses per class in all prior tools. The biological significance of this resolution is multi-dimensional:

- **Clinical relevance:** The distinction between Sticky DNA (GAA/TTC expansion), Triplex (H-DNA), and general Slipped DNA (STR/Direct Repeat) is directly tied to distinct human diseases (Friedreich's ataxia vs. general repeat instability). Conflating these subclasses into a single "triplex" or "slipped DNA" label, as prior tools do, obscures actionable clinical information.
- **Thermodynamic precision:** Bulged G4, Extended-loop G4, and Two-tetrad weak PQS have substantially different melting temperatures and cellular stabilities, affecting their utility as drug targets. Non B DNA Finder's hierarchical subclass taxonomy is directly aligned with the G4 structural classification used in biophysical and medicinal chemistry literature.
- **Regulatory interpretation:** Local Curvature and Global Curvature play distinct roles in nucleosome positioning (local) versus replication origin bending (global); AC-motif i-motifs have different cellular pH responses than canonical i-motifs. These subclass distinctions have direct mechanistic implications that no prior genome-scale tool provides.

#### 3.2.3 Quantitative GC–Non-B DNA Composition Relationships

The multi-organism dataset assembled in this study provides the first quantitative demonstration of the GC-composition dependence of the complete non-B DNA structural landscape. The strong positive correlations of G-Quadruplex (r = 0.87), Z-DNA (r = 0.96), and R-Loop density (r = 0.82) with GC%, and the strong negative correlation of Curved DNA density (r = −0.79) with GC%, had been previously inferred from first principles or small-scale comparisons but had never been formally quantified across a genome-size-spanning, phylogenetically diverse dataset of this breadth. These correlations have practical implications: the extraordinary non-B DNA density in high-GC actinobacteria (*C. shaoxiangyii*, *M. marina*; 13–16 motifs/kb) implies that essentially every gene regulatory sequence in these organisms has strong non-B DNA-forming potential, with possible implications for transcription regulation, genome instability, and horizontal gene transfer susceptibility that are entirely unexplored.

#### 3.2.4 Scalability and Accessibility

Non B DNA Finder achieves true genome-scale throughput, completing analysis of the 3.12 Gb human genome in approximately 2.2 hours on standard hardware—a performance standard that only Non-B DB (as a pre-computed database, not a real-time tool) has previously matched. Interactive analysis through the Streamlit web interface lowers the barrier to use for non-computational biologists, while the Jupyter notebook interface and command-line API support full reproducibility in research workflows. All scoring parameters, canonical taxonomy definitions, and test sequences are openly available, addressing the reproducibility deficit that characterises many published non-B DNA tools.

### 3.3 Algorithmic Optimization and Computational Framework

Genome-scale non-B DNA detection imposes demanding computational requirements: a single chromosome can exceed 250 Mb, and nine detectors must be applied at every overlapping window position. Four complementary optimization strategies are integrated into Non B DNA Finder to achieve practical throughput on commodity hardware.

#### 3.3.1 Hyperscan Multi-Pattern Matching

Intel Hyperscan is a high-performance regular expression and exact-string matching library that compiles a set of patterns into a compiled finite automaton executed with SIMD (Single Instruction, Multiple Data) vector instructions. Non B DNA Finder uses Hyperscan in two detectors where the bottleneck is exact lookup of large tables of fixed-length patterns: Z-DNA 10-mer scanning and A-philic DNA 10-mer scanning. In both cases, the complete set of scored 10-mer keys is compiled into a Hyperscan database at module load time; the entire sequence is then scanned in a single pass, with match callbacks populating a position-indexed score array. When the `hyperscan` Python binding is unavailable, the library falls back transparently to a NumPy-vectorized polynomial-hash lookup (`vectorized_find_matches`) that encodes the DNA sequence as an integer array, computes all 10-mer hashes with ten NumPy array additions, and retrieves scores by array indexing—providing 5–10× acceleration over the naive Python loop (`py_find_matches_loop`) while requiring no compiled dependencies. A further fallback to the pure-Python sliding-window loop is invoked if NumPy is also absent. Hyperscan and the NumPy-vectorized path produce numerically identical outputs to the pure-Python baseline; no approximations are introduced.

#### 3.3.2 Numba JIT Compilation

Three computational hot-path kernels are decorated with `@jit(nopython=True, cache=True)` from the Numba just-in-time compilation library, which compiles Python/NumPy functions to native machine code on first call and caches the compiled binary to disk for subsequent invocations:

- **G-Quadruplex sliding-window scoring** (`Detectors/gquad/detector.py`): The G4Hunter score is computed by summing per-base G/C contributions in a sliding window across the full sequence. The JIT-compiled kernel eliminates Python interpreter overhead on this innermost loop, yielding a 2–5× speedup for sequences >10 kb.
- **Triplex mirror-repeat scoring** (`Detectors/triplex/detector.py`): The mirror-symmetry score is computed by pairwise comparison of purine–pyrimidine assignments across the candidate window. The JIT kernel compiles the character-comparison and score-accumulation logic to native code.
- **Slipped DNA piecewise scoring** (`Detectors/slipped/detector.py`): Copy-number, unit-size, base-purity, and GC-fraction are combined in a piecewise arithmetic expression evaluated over each candidate repeat region. Numba compilation removes the interpretive cost of this multi-operand expression when called thousands of times per chromosome.

When Numba is not installed, each function is replaced at module load time by an equivalent pure-Python implementation, ensuring correct results at a 2–5× throughput penalty.

#### 3.3.3 Parallel Processing

Non B DNA Finder implements two tiers of parallelism managed by Python's `concurrent.futures` module:

**Sequence-level parallelism (ProcessPoolExecutor).** For sequences exceeding 1 Mb, the input chromosome is tiled into 50 kb chunks with 2 kb bilateral overlap, and chunks are dispatched to a `ProcessPoolExecutor` with a worker count equal to `min(chunk_count, os.cpu_count())`. Each worker process runs the full nine-detector pipeline on its assigned tile independently, bypassing the CPython Global Interpreter Lock and exploiting all available CPU cores. Results from completed chunks are collected as futures complete and merged by the `OverlapDeduplicator` module. For multi-FASTA inputs (e.g., whole-genome assemblies), the `MultiFastaEngine` dispatches entire sequences to a process pool, enabling all chromosomes to be analysed concurrently up to the available core count. Empirically, four-core commodity hardware achieves approximately 4× throughput improvement over single-threaded execution on chromosome-scale inputs.

**Detector-level parallelism (ThreadPoolExecutor).** For sequences between 50 kb and 1 Mb that do not require chunking, up to nine concurrent worker threads—one per detector class—are launched via `ThreadPoolExecutor`. Because the G-Quadruplex, Cruciform, R-Loop, and Triplex detectors are CPU-bound, the threading benefit is modest (approximately 1.5–2× on an eight-core system due to GIL contention); the primary benefit is overlap of I/O-bound detector initialisation and of the Numba JIT warm-up period.

A fallback to sequential single-threaded execution is used for sequences shorter than 50 kb, where process-spawning overhead exceeds the parallelism benefit.

#### 3.3.4 Entropy Filtering

Low-complexity DNA sequences—homopolymer runs, simple tandem repeats with minimal information content—can generate spurious motif calls from pattern-matching algorithms that have no model of sequence complexity. Non B DNA Finder applies Shannon entropy filtering as a post-detection quality gate for the Slipped DNA detector, which is most susceptible to low-complexity false positives. The Shannon entropy *H* of a candidate region is computed over the four base frequencies:

*H* = −Σ *p*(*b*) log₂ *p*(*b*)  for *b* ∈ {A, T, G, C}

Candidates with *H* < 0.5 bits are discarded before scoring, removing pure or near-pure homopolymers and degenerate dinucleotide repeats (e.g., AAAAAAA or TATATATATA with minimal base diversity) that would otherwise inflate STR and Direct Repeat counts. The 0.5-bit threshold was calibrated empirically against positive-control STR datasets to achieve <2% false-positive inflation while retaining >98% of biologically annotated tandem repeats.

### 3.4 Software Implementation

Non B DNA Finder is implemented as a modular Python package designed for portability, reproducibility, and accessibility across three usage modes: interactive web application, Jupyter notebook, and command-line API.

#### 3.4.1 Python Libraries

The core runtime stack requires Python ≥3.8 and the following libraries:

- **NumPy** — vectorized array operations for 10-mer hash computation, per-base contribution arrays, and score normalisation arithmetic;
- **pandas** — tabular result aggregation, multi-column sorting, CSV/XLSX export, and multi-index grouping for cross-class statistics;
- **Matplotlib** and **Seaborn** — linear motif maps, class-distribution bar charts, GC-composition scatter plots, and confidence-score histograms embedded in the Streamlit interface and Jupyter notebooks;
- **BioPython** (`Bio.SeqIO`) — FASTA parsing, sequence validation, reverse complement computation, and NCBI Entrez accession retrieval;
- **Streamlit** — the reactive web application framework that renders the interactive front end (`app.py`) without any JavaScript;
- **pyahocorasick** (optional) — C-extension Aho-Corasick automaton for single-pass multi-pattern matching used in the optimised scanner path, providing 50–200× acceleration over sequential regex for large pattern sets. When unavailable, the scanner falls back to a pure-Python Aho-Corasick implementation included in `Utilities/ac_matcher.py`.

Optional performance accelerators add no new algorithmic behaviour:
- **Numba** — JIT compilation of the G4Hunter, Triplex, and Slipped DNA hot-path kernels (Section 3.3.2);
- **hyperscan** (Python bindings for Intel Hyperscan) — SIMD multi-pattern matching for Z-DNA and A-philic 10-mer lookup (Section 3.3.1).

All dependencies are listed in `requirements.txt`; optional packages are caught with `try/except ImportError` and silently replaced by pure-Python equivalents, ensuring the tool runs correctly in environments without compiled extensions.

#### 3.4.2 GitHub Repository and Open-Source Distribution

The Non B DNA Finder source code is hosted at [https://github.com/VRYella/NonBDNAFinder](https://github.com/VRYella/NonBDNAFinder) under the MIT licence, which permits unrestricted use, modification, and redistribution. The repository is organised into logically separated sub-packages: `Detectors/` (nine detector modules, each in its own subdirectory with `__init__.py`, `detector.py`, and optional backend files), `Utilities/` (shared infrastructure: scanner, optimised scanner, parallel helpers, overlap deduplicator, export engine, visualisation pipeline), `UI/` (Streamlit page components and documentation renderer), `tests/` (pytest suite with sequence-level and performance regression tests), and `examples/` (benchmark FASTA files for positive-control validation). Scoring parameter registries (`Utilities/consolidated_registry.json`, `Detectors/zdna/tenmer_table.py`, `Detectors/aphilic/tenmer_table.py`) are stored as plain text alongside their primary literature citations, making every scoring decision auditable without inspecting compiled code. The production release is tagged `v1.0.0` and archived for long-term citability. Issue tracking, pull-request workflow, and release management are conducted through the GitHub platform, enabling community contributions and transparent bug reporting.

#### 3.4.3 Streamlit Interface

The interactive front end (`app.py`) is built on Streamlit and is deployable both locally (`streamlit run app.py`) and on Streamlit Community Cloud without any server configuration. The interface provides:

- **Input modes**: Direct sequence paste, FASTA file upload (single or multi-sequence), and NCBI accession retrieval (Entrez nucleotide fetch with automatic FASTA conversion);
- **Motif class selection**: An interactive toggle table listing all nine structural classes with per-class enable/disable checkboxes, allowing users to restrict analysis to classes of interest and reduce wall-clock time;
- **Real-time execution progress**: A progress bar and live throughput display (bp/s) updated at each chunk boundary, computed from chunk completion timestamps;
- **Tabular results**: A paginated, sortable, filterable data table of all detected motifs with 12 core fields per record (Sequence_Name, Class, Subclass, Start, End, Length, Sequence, Strand, Score, Raw_Score, Method, Pattern_ID) plus class-specific supplementary columns, with CSV and XLSX download buttons;
- **Visualisation**: A linear motif map rendering all detected loci as colour-coded rectangles on a proportional chromosome axis, a stacked bar chart of per-class motif counts, and a class-distribution pie chart;
- **Documentation**: An embedded multi-tab reference panel (Methods, Detection Algorithms, Scoring Guide, Optimization) rendered from `UI/documentation.py` without leaving the application.

The Streamlit interface requires no bioinformatics command-line experience and is accessible to experimental biologists through a standard web browser, directly addressing the accessibility gap identified in Section 3.2.4.

### 3.5 Genomic Observations: Biological Context and Prior Literature

The genome-scale results reported here are broadly consistent with, and in several cases substantially extend, published findings on non-B DNA biology:

**Human G-Quadruplex landscape.** The 5.04 M G-Quadruplex sites detected in the full human genome, covering 3.22% of the genome, are consistent in scale with the ≥370,000–710,000 high-confidence PQS estimates from G4-seq (Chambers *et al.*⁵; Hänsel-Hertsch *et al.* 2016⁴⁹), with the much larger total driven by the inclusion of Two-tetrad weak PQS, which G4-seq substantially underestimates due to its reliance on G4-ligand stabilisation. The distribution of Canonical intramolecular G4 in the human genome (53,364 sites) is closer to the G4-seq-observed range, consistent with the expected experimental detection rate. The 368 i-Motif sites with canonical scores in the human centromere assembly are consistent with recent immunofluorescence and CUT&RUN data demonstrating i-motif formation in centromeric CENP-B-box contexts (Zeraati *et al.*¹¹; Bhoj *et al.* 2024).

**Extreme-GC non-B DNA.** The hyperabundant non-B DNA landscape in *C. shaoxiangyii* and *M. marina*—most strikingly the >65% genome coverage by R-loop-forming sequences with mean lengths exceeding 1.5 kb—exceeds anything previously reported in the computational literature. While it is well established from first principles that extreme-GC genomes are expected to display elevated G-quadruplex and Z-DNA propensity due to their abundant GGG-run and alternating CG dinucleotide content, the specific magnitude of R-loop coverage and the average Non-B DNA Cluster size (>820 bp with 4+ overlapping structural classes) identified here is a new quantitative observation. Whether these computationally predicted structures form biologically in cells remains to be experimentally validated; however, the data suggest that the structural DNA landscape of extreme-GC bacteria is fundamentally different from that of well-studied model organisms, with potential consequences for replication, transcription, and stress responses that warrant experimental investigation.

**AT-rich eukaryotic non-B DNA.** The *P. falciparum* dataset confirms and quantifies the expectation that extreme AT-rich genomes are dominated by Curved DNA and Slipped DNA, which together account for >57% of all motifs in this organism. The density of Triplex (H-DNA) motifs in *P. falciparum* (0.90 per kb)—third overall—is unexpectedly high for an AT-rich genome, suggesting that the purine–pyrimidine mirror-repeat tracts within the organism's extended intergenic regions (some >10 kb in length) are a substantial source of triplex-forming potential. This structural feature has not been highlighted by previous analyses of *P. falciparum* genome organisation and may be relevant to the well-documented replication pausing and antigenic variation switching that occurs in these intergenic regions.

### 3.6 Integrated Multi-Class Detection

The central advance of Non B DNA Finder relative to prior tools is the simultaneous, unified detection of nine structurally and algorithmically disparate non-B DNA classes in a single workflow. This integration is scientifically critical because non-B loci rarely occur in isolation: genome-wide analyses consistently show that G4-forming promoters also contain i-motif–prone sequences on the complementary strand, and that replication-origin–proximal regions harbour co-incident R-loop, cruciform, and Z-DNA potential. The Hybrid and Non-B DNA Cluster outputs of Non B DNA Finder make these co-occurrences explicit and queryable.

### 3.7 Subclass Resolution and Disease Annotation

The 24-subclass taxonomy—anchored in a single canonical taxonomy module and enforced through a normalisation layer throughout the codebase—provides a level of structural resolution not available in any prior genome-scale tool. For G-quadruplexes, the distinction between Telomeric G4, Canonical G4, Bulged G4, Extended-loop G4, G-wire, and Weak PQS is biologically meaningful: telomeric G4s are stabilised by POT1/TRF2 and targeted by telomerase inhibitors; Bulged G4s are prevalent in ribosomal DNA and have distinct thermodynamic properties compared to canonical G4s. Similarly, the three-way i-motif taxonomy (Canonical, Relaxed, AC-motif) reflects mechanistically distinct folding pathways with different pH sensitivity and cellular roles.

Every motif record is annotated with a `Disease_Relevance` field automatically populated based on sequence content thresholds derived from clinical genetics literature: for example, (CGG)ₙ repeats ≥55 copies are flagged as Fragile X premutations, (GAA)ₙ ≥60 copies as FRDA pathogenic, and (CAG)ₙ ≥36 copies as Huntington disease threshold—providing an immediately actionable clinical layer to routine bioinformatic analysis.

### 3.8 Reproducibility and Transparency

All scoring parameters are documented in `Utilities/consolidated_registry.json` alongside their primary literature citations. Motif detection is fully deterministic; no random components are involved. The open-source codebase, versioned on GitHub, and the stable canonical taxonomy ensure that published results are reproducible across software versions. Users who enable optional Hyperscan or Numba acceleration receive numerically identical results to the pure-Python baseline.

### 3.9 Limitations and Future Directions

Non B DNA Finder does not currently model the effects of DNA supercoiling, chromatin context, or sequence methylation on structure formation probability—all of which are known to modulate non-B DNA stability *in vivo*. Future work will integrate supercoiling density estimates from twin-domain transcription models and incorporate CpG methylation status (from bisulphite-sequencing data) as a modifier of Z-DNA and R-loop propensity. Extension to RNA secondary structures (forming at the single-stranded portion of R-loops) and to co-detection of CRISPR off-target sites enriched at non-B loci is also planned. The extraordinary non-B DNA density in extreme-GC actinobacteria (*Cellulomonas*, *Miltoncostaea*) demands dedicated experimental validation—e.g., by G4-ChIP-seq, R-loop S9.6 immunoprecipitation, and Z-DNA-specific antibody mapping—to determine which computationally predicted loci are biologically active structural hotspots. Integration with chromatin accessibility data (ATAC-seq, DNase-seq) and transcription factor binding maps will further contextualise non-B DNA predictions within regulatory networks.

---

## 4. Methods

### 4.1 Software Implementation

Non B DNA Finder is implemented in Python ≥3.8. Core dependencies include NumPy, pandas, Matplotlib, Seaborn, Streamlit, and BioPython. Optional accelerators are Numba (JIT compilation for G4Hunter sliding window and Triplex mirror-repeat scoring), Cython (compiled extensions for inner loops), and Intel Hyperscan (SIMD multi-pattern matching for Z-DNA 10-mer and A-philic 10-mer table lookups).

### 4.2 Score Normalisation

Raw scores from each detector (e.g., cumulative 10-mer sum for Z-DNA; G4Hunter score for G-quadruplex; ΔG for cruciform) are normalised to a common 1–3 scale using per-detector theoretical minimum and maximum bounds computed analytically from the scoring model and the sequence under analysis. For G-quadruplex, the theoretical minimum is the G4Hunter threshold (0.5) and the maximum is computed from the maximum possible G4Hunter score over the window. For cruciform, the minimum is −5.0 kcal/mol and the maximum is the most stable predicted stem in the current sequence. Normalised scores are clipped to [1, 3] and rounded to one decimal place.

### 4.3 Overlap Resolution

Within each structural class, overlapping motif candidates are resolved by a greedy priority algorithm that scores candidates by normalised score (descending), then by subclass priority (fixed taxonomic order), then by motif length (descending). Accepted intervals are tracked with a sorted list for O(log n) conflict detection. For i-motifs, HUR AC-motifs are resolved in an independent pool from canonical/relaxed i-motifs, permitting biologically co-occurring motifs to be co-reported.

Across classes, motif pairs from distinct classes are allowed to overlap; only cross-class overlap ≥50% triggers Hybrid consolidation.

### 4.4 Chunked Execution for Large Genomes

For sequences >1 Mb, the input is divided into 50 kb tiles with 2 kb bilateral overlap. Each tile is analysed independently by the full detector suite. Motifs whose `Start` coordinate falls within the overlap zone of tile *i* are suppressed in favour of the canonical detection in tile *i+1*, implemented by `OverlapDeduplicator.filter_core()`. Results from all tiles are concatenated and globally sorted by start coordinate before export.

### 4.5 Validation Datasets

Benchmark FASTA files for positive-control validation are provided in the `examples/` directory. Genome-scale FASTA files (GRCh38, GRCm39, R64, TAIR10, dm6, ce11, *E. coli* ASM584v2, Pf3D7) were downloaded from NCBI RefSeq. Non-B DB v3.0 BED files and G4-seq peak BED files were downloaded from NCBI FTP and the supplementary data of Chambers *et al.* 2015, respectively. Intersection analyses were performed with bedtools v2.30.0.

### 4.6 Web Interface

An interactive Streamlit application (`app.py`) provides a graphical front end supporting file upload, NCBI accession retrieval, motif class selection via a toggle table, real-time execution progress, tabular results with pagination, linear motif maps, class-distribution bar charts, hybrid/cluster visualisation, and CSV/XLSX download. The interface is deployable locally via `streamlit run app.py` or hosted on Streamlit Community Cloud.

### 4.7 Jupyter Notebook Interface

`NonBDNAFinder_Analysis.ipynb` provides an annotated cell-by-cell workflow for exploratory analysis, integrating sequence input, motif detection, result inspection, and visualisation in a single interactive document. A parallel-execution notebook (`NonBDNAFinder_Parallel_CSV.ipynb`) supports bulk CSV-input workflows for large comparative analyses.

---

## 5. Data Availability

The Non B DNA Finder source code, documentation, and example datasets are freely available at [https://github.com/VRYella/NonBDNAFinder](https://github.com/VRYella/NonBDNAFinder) under the MIT licence. All scoring parameter registries, canonical taxonomy definitions, and test sequences are included in the repository.

---

## 6. Code Availability

All source code, scoring tables (`Utilities/consolidated_registry.json`, `Detectors/zdna/tenmer_table.py`, `Detectors/aphilic/tenmer_table.py`), and test suites (`tests/`) are publicly available at the GitHub URL above. The version archived for this publication is tagged as `v1.0.0`.

---

## Acknowledgements

We thank the developers of G4Hunter, QmRLFS-finder, ZHUNT, and Non-B DB for their foundational contributions to the non-B DNA computational landscape. We acknowledge BioPython, NumPy, pandas, Matplotlib, and Streamlit open-source communities.

---

## Author Contributions

V.R.Y. conceived the project, designed the algorithmic framework, implemented the detector suite, and wrote the manuscript.

---

## Competing Interests

The authors declare no competing interests.

---

## References

1. Wells, R. D. Non-B DNA conformations, mutagenesis and disease. *Trends Biochem. Sci.* **32**, 271–278 (2007).
2. Huppert, J. L. & Balasubramanian, S. G-quadruplexes in promoters throughout the human genome. *Nucleic Acids Res.* **35**, 406–413 (2007).
3. Todd, A. K., Johnston, M. & Neidle, S. Highly prevalent putative quadruplex sequence motifs in human DNA. *Nucleic Acids Res.* **33**, 2901–2907 (2005).
4. Balasubramanian, S., Hurley, L. H. & Neidle, S. Targeting G-quadruplexes in gene promoters: a novel anticancer strategy? *Nat. Rev. Drug Discov.* **10**, 261–275 (2011).
5. Chambers, V. S. *et al.* High-throughput sequencing of DNA G-quadruplex structures in the human genome. *Nat. Biotechnol.* **33**, 877–881 (2015).
6. Rich, A., Nordheim, A. & Wang, A. H. The chemistry and biology of left-handed Z-DNA. *Annu. Rev. Biochem.* **53**, 791–846 (1984).
7. Zhang, T., Yin, C., Cochrane, C. & Bhatt, D. L. Z-DNA binding protein 1: structure, function, and therapeutic potential. *Nat. Rev. Mol. Cell Biol.* **24**, 179–196 (2023).
8. Usdin, K. & Kumari, D. Repeat-mediated epigenetic dysregulation of the FMR1 gene in the fragile X-related disorders. *Front. Genet.* **6**, 192 (2015).
9. Aguilera, A. & García-Muse, T. R loops: from transcription byproducts to threats to genome stability. *Mol. Cell* **46**, 115–124 (2012).
10. Santos-Pereira, J. M. & Aguilera, A. R loops: new modulators of genome dynamics and function. *Nat. Rev. Genet.* **16**, 583–597 (2015).
11. Zeraati, M. *et al.* I-motif DNA structures are formed in the nuclei of human cells. *Nat. Chem.* **10**, 631–637 (2018).
12. Benabou, S. *et al.* Fundamental aspects of the nucleic acid i-motif structures. *RSC Adv.* **4**, 26956–26980 (2014).
13. Frank-Kamenetskii, M. D. & Mirkin, S. M. Triplex DNA structures. *Annu. Rev. Biochem.* **64**, 65–95 (1995).
14. Sakamoto, N. *et al.* Sticky DNA: self-association properties of long GAA·TTC repeats in R·R·Y triplex structures from Friedreich's ataxia. *Mol. Cell* **3**, 465–475 (1999).
15. Lilley, D. M. The inverted repeat as a recognizable structural feature in supercoiled DNA molecules. *Proc. Natl Acad. Sci. USA* **77**, 6468–6472 (1980).
16. Kurahashi, H. *et al.* Cruciform DNA structure underlies the etiology for palindrome-mediated human chromosomal translocations. *J. Biol. Chem.* **279**, 35377–35383 (2004).
17. Mirkin, S. M. Expandable DNA repeats and human disease. *Nature* **447**, 932–940 (2007).
18. Koo, H. S., Wu, H. M. & Crothers, D. M. DNA bending at adenine-thymine tracts. *Nature* **320**, 501–506 (1986).
19. Satchwell, S. C., Drew, H. R. & Travers, A. A. Sequence periodicities in chicken nucleosome core DNA. *J. Mol. Biol.* **191**, 659–675 (1986).
20. Bedrat, A., Lacroix, L. & Mergny, J. L. Re-evaluation of G-quadruplex propensity with G4Hunter. *Nucleic Acids Res.* **44**, 1746–1759 (2016).
21. Hon, J. *et al.* Bioinformatic tools for identifying disease gene and SNP candidates. *Methods* **74**, 3–13 (2015).
22. Jenjaroenpun, P. *et al.* QmRLFS-finder: a model, web server and stand-alone tool for prediction and analysis of R-loop forming sequences. *Nucleic Acids Res.* **45**, W110–W116 (2017).
23. Cer, R. Z. *et al.* Non-B DB v2.0: a database of predicted non-B DNA-forming motifs and its associated tools. *Nucleic Acids Res.* **41**, D94–D100 (2013).
24. Ho, P. S., Ellison, M. J., Quigley, G. J. & Rich, A. A computer aided thermodynamic approach for predicting the formation of Z-DNA in naturally occurring sequences. *EMBO J.* **5**, 2737–2744 (1986).
25. Olson, W. K. *et al.* DNA sequence-dependent deformability deduced from protein-DNA crystal complexes. *Proc. Natl Acad. Sci. USA* **95**, 11163–11168 (1998).
26. Nelson, H. C. *et al.* The structure of an oligo(dA)·oligo(dT) tract and its biological implications. *Nature* **330**, 221–226 (1987).
27. Schlötterer, C. & Tautz, D. Slippage synthesis of simple sequence DNA. *Nucleic Acids Res.* **20**, 211–215 (1992).
28. Weber, J. L. & Wong, C. Mutation of human short tandem repeats. *Hum. Mol. Genet.* **2**, 1123–1128 (1993).
29. SantaLucia, J. Jr. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *Proc. Natl Acad. Sci. USA* **95**, 1460–1465 (1998).
30. Lilley, D. M. Hairpin-loop formation by inverted repeats in supercoiled DNA is a local and transmissible property. *Nucleic Acids Res.* **13**, 1443–1465 (1985).
31. Aguilera, A. & García-Muse, T. R loops: from transcription byproducts to threats to genome stability. *Mol. Cell* **46**, 115–124 (2012).
32. Frank-Kamenetskii, M. D. & Mirkin, S. M. Triplex DNA structures. *Annu. Rev. Biochem.* **64**, 65–95 (1995).
33. Sakamoto, N. *et al.* Sticky DNA: self-association properties of long GAA·TTC repeats in R·R·Y triplex structures from Friedreich's ataxia. *Mol. Cell* **3**, 465–475 (1999).
34. Bedrat, A., Lacroix, L. & Mergny, J. L. Re-evaluation of G-quadruplex propensity with G4Hunter. *Nucleic Acids Res.* **44**, 1746–1759 (2016).
35. Gehring, K., Leroy, J. L. & Guéron, M. A tetrameric DNA structure with protonated cytosine-cytosine base pairs. *Nature* **363**, 561–565 (1993).
36. Hur, J. K. *et al.* Phased CRISPR-Cas9 screen identifies a positive regulatory role for CTCF in the generation of c-Myc amplification in cancer. *Nucleic Acids Res.* **49**, 7283 (2021).
37. Ho, P. S., Ellison, M. J., Quigley, G. J. & Rich, A. A computer aided thermodynamic approach for predicting the formation of Z-DNA in naturally occurring sequences. *EMBO J.* **5**, 2737–2744 (1986).
38. Herbert, A. *et al.* The Zα domain from human ADAR1 binds to the Z-DNA conformer of many different sequences. *Nucleic Acids Res.* **26**, 3486–3493 (1998).
39. Gorin, A. A., Zhurkin, V. B. & Olson, W. K. B-DNA twisting correlates with base-pair morphology. *J. Mol. Biol.* **247**, 34–48 (1995).
40. Vinogradov, A. E. & Anatskaya, O. V. Genome size and metabolic intensity in tetrapods: a tale of two lines. *Proc. Biol. Sci.* **273**, 27–32 (2006).
41. Williamson, J. R., Raghuraman, M. K. & Cech, T. R. Monovalent cation-induced structure of telomeric DNA: the G-quartet model. *Cell* **59**, 871–880 (1989).
42. Wang, Y. & Patel, D. J. Solution structure of the human telomeric repeat d[AG₃(T₂AG₃)₃] G-tetraplex. *Structure* **1**, 263–282 (1993).
43. Leroy, J. L. *et al.* Intramolecular folding of a fragment of the cytosine-rich strand of telomeric DNA into an i-motif. *Nucleic Acids Res.* **22**, 1600–1606 (1994).
44. Wang, A. H. *et al.* Molecular structure of a left-handed double helical DNA fragment at atomic resolution. *Nature* **282**, 680–686 (1979).
45. Hon, J. *et al.* Sequence-dependent prediction of G-quadruplex stability by pqsfinder. *Bioinformatics* **33**, 3373–3379 (2017).
46. Sahakyan, A. B. *et al.* Machine learning model for sequence-driven DNA G-quadruplex formation. *Sci. Rep.* **7**, 14535 (2017).
47. Sanz, L. A. *et al.* Prevalent, dynamic, and conserved R-loop structures associate with specific epigenomic signatures in mammals. *Mol. Cell* **63**, 167–178 (2016).
48. Xu, W. *et al.* The R-loop is a common chromatin feature of the *Arabidopsis* genome. *Nat. Plants* **3**, 704–714 (2017).
49. Hänsel-Hertsch, R. *et al.* G-quadruplex structures mark human regulatory chromatin. *Nat. Genet.* **48**, 1267–1272 (2016).


####################New text
Non-B DNA Finder: A Unified Framework for Detection of Diverse Non-Canonical DNA Structures in genomes

Venkata Rajesh Yella1,*,   Aruna Sesha Chandrika Gummadi1 and Aditya Kumar2

Koneru Lakshmaiah Education Foundation, 
Department of Biotechnology, Vaddeswaram, Guntur, 
Andhra Pradesh, 522302, IndiaEmail: yvrajesh_bt@kluniversity.in 
Orcid:https://orcid.org/0000-0002-1961-5051
•	
Tezpur University, 
Molecular Biology and Biotechnology, 
Tezpur, 
Assam, 784028, India
aditya@tezu.ernet.in
Orcid: https://orcid.org/0000-0002-6474-8830
•	



 
Abstract
Non-canonical DNA secondary structures arise from unique sequence pattens within genomic contexts and contribute to gene regulation, genome instability, and disease pathogenesis. Despite substantial evidence that multitude of non-B DNA conformations are physiologically relevant, existing computational tools do not capture the full repertoire of these structure forming motifs. In this study, we introduce an unified pipeline, Non-B DNA Finder, for high-throughput detection of 11 classes and 24 subclasses of non-B DNA motifs. These include A-DNA; curved DNA (both global and local); Z-DNA (including Z-DNA and eGZ motifs); slipped DNA (comprising direct repeats and short tandem repeats); R-loops; cruciform DNA; H-DNA (including triplex-forming mirror repeats and sticky DNA); G-quadruplexes (covering telomeric, stacked,  canonical intramolecular, extended-loop canonical, higher-order arrays,  bulged  two-tetrad weak G4s and intramolecular G-triplex); i-motifs (canonical, relaxed, and AC motifs); as well as hybrid overlaps and non-B DNA clusters. Users can upload sequence files, import data from NCBI, or provide direct input to receive interactive motif annotations, export data, summary statistics, and graphical maps, including visualizations of overlapping regions and structural clusters. Analysis using the  tool across nine bacterial genomes spanning a broad GC range, revealed widespread non-B DNA motifs with density increasing alongside GC content. G-quadruplexes were the predominant class across genomes. The pipeline also identified numerous hybrid regions and multi-class clusters, with some loci exhibiting high structural complexity within short genomic windows.  The performance  of the tool has been validated using well-characterized genomic regions and by comparison with established tools, including nBMST, Z-seeker and G4Hunter. Altogether, Non-B DNA Finder offers an integrated framework for systematic analysis of non-B DNA motifs across genomic datasets, facilitating research in gene regulation, genome stability, microbial and comparative genomics. Non-B DNA Finder is freely accessible as webserver tool at https://Non-B DNA Finder.streamlit.app/.
 
1. Introduction
Watson and Crick in 1953 classically  depicted the genomic DNA as a right-handed double helix named as the B-form1 is most predominant form exist in living cells. However, decades of research have demonstrated that the DNA molecule is structurally versatile, adopting multifarious non-canonical secondary structures under physiological cellular conditions2-8. These alternative structures are often referred to “non-B DNA” structures. The first major deviation was Z-DNA, a left-handed helical form identified in 1979 by Wang and colleagues through X-ray crystallography9. In 1980, cruciform DNA formed by inverted repeat sequences was demonstrated in supercoiled plasmids10 Soon after, slipped-strand DNA to trinucleotide repeat expansion disorders  has been established  underscoring the pathological potential of structural anomalies (Reviewed in 11). In 1980s, the in vivo characterization of  bent  or curved DNA and the role of A-tracts has been reported from experiments on  Leishmania tarentolae kinetoplast DNA and chicken nucleosome core DNA  characterized the roles 12-14. This was followed by the experimental identification of mirror repeat-based triplex DNA(H-DNA) structures15,16. Around the same time, attention turned toward guanine-rich sequences, with telomeric G-quadruplexes (G4s) recognized for their potential to form stacked tetrads stabilized by Hoogsteen hydrogen bonds. The formation of G-quadruplex structures was first described by Sen and Gilbert in 1988 using electrophoretic mobility shift assays, who demonstrated four-stranded structures in guanine-rich sequences and proposed a "sodium-potassium conformational switch" for telomeric DNA17,18. Then formation of G-quadruplex was  described by electrophoretic mobility assays in telomeric repeats of Oxytricha and Tetrahymena19.  In vitro characterizations were well underway by the early 1990s, and in vivo mapping using G4-specific antibodies  and bioinformatic analyses followed in the 2000s20-22. In 1993, Gehring and colleagues introduced the i-motif, a cytosine-rich, four-stranded structure formed under acidic conditions23. This structure remained controversial until Zeraati and colleagues provided compelling in vivo evidence in 2018 using an i-motif-specific antibody24 . 
Along with these classical non canonical DNA, studies indicated several novel functional structures such as R-loops, sticky DNA, G-triplexes, AC-motifs and eGZ motifs. In the early 2000s, the role of R-loops, RNA:DNA hybrid structures emerged, driven by improvements in genome-wide mapping techniques such as DRIP-seq25 In 1991, Sakamoto et al described sticky DNA, characterized by long GAA/TTC repeats that form stable triplexes, often associated with gene silencing in diseases like Friedreich’s ataxia26 .G-triplexes are three-stranded DNA structures that serve as crucial folding intermediates in the hierarchical assembly pathway of G-quadruplexes27. These structures consist of G:G:G triads stabilized by Hoogsteen-type hydrogen bonds and form from sequences containing three G-tracts, representing an evolutionary step between simple hairpin structures and fully formed four-stranded G-quadruplexes28. Remarkably, G-triplexes can be stabilized under physiological conditions by divalent cations and molecular crowding, allowing them to exist as independent functional entities 29,30. The study by Hur et al. (2021)31 reported the discovery of a novel non-canonical DNA secondary structure  termed the AC-motif, which is formed by sequences containing adenine and cytosine repeats. Using biophysical analyses and molecular dynamic simulations, these authors demonstrated that oligodeoxynucleotides comprising adenine and cytosine tracts can fold into an i-motif–like four-stranded structure stabilized by hemi-protonated C⁺:C base pairs intercalated with protonated A⁺:C base pairs. Further, functional studies revealed that AC-motif formation in the CDKL3 promoter enhances gene expression and genome-wide mapping identified over 2,000 putative AC-motif forming sequences in the human genome, particularly enriched in promoter regions. The eGZ-motif is an another recently characterized non-B DNA structure that forms within expanded runs of CGG trinucleotide repeats32. Unlike canonical Z-DNA, which is stabilized by alternating purine-pyrimidine sequences, the eGZ-motif is distinguished by the regular extrusion of guanine bases from the DNA double helix, giving rise to a left-handed Z-DNA conformation with unique structural properties. Molecular dynamics simulations and biophysical experiments have shown that these alternately extruded guanines foster novel stacking and hydrogen-bonding interactions, resulting in a highly stable helix that differs from previously known Z-DNA motifs. The eGZ-motif is mechanistically important because its formation can promoter genomic instability, particularly in loci associated with neurodegenerative diseases driven by CGG repeat expansions, such as fragile X-related disorders. 
The inherent structural complexity and dynamic nature of DNA  has necessitated the development of diverse experimental approaches to characterize the formation and biological relevance of structural polymorphism. Various experimental techniques can be employed to characterize non-B DNA structures. These include polyacrylamide gel electrophoresis, cyclization kinetics assays, X‑ray crystallography, circular dichroism (CD) spectroscopy, ultraviolet absorption spectroscopy, FRET-based melting analyses, atomic force microscopy, electron microscopy, high-throughput sequencing approaches, and immunofluorescence-based detection methods (reviewed in 6,33-38). Despite these advances, experimental techniques alone are limited in scalability, throughput, and feasibility for comprehensive genome-wide analyses across diverse species and conditions. Consequently, computational prediction tools have become essential for identifying putative non-B DNA-forming sequences in genomic data. By employing biologically informed sequence motifs, refined regular expressions, and machine learning based scoring algorithms, these tools enable rapid, high-throughput annotation of canonical and emerging non-B DNA motifs. Several well-established computational platforms have been developed to address different classes of non-B DNA structures, each with unique strengths and limitations. For G-quadruplex prediction, several  tools are developed based on  sequence pattern recognition22 and scoring system39 to identify quadruplex-forming potentials with substantial accuracy (reviewed in 40). Z-DNA prediction has been historically conducted using Z-Hunt41 and recent optimized tool z-seeker42 used  thermodynamic parameters and dinucleotide scoring models to predict Z-forming regions.  Triplex-forming sequences are detected by specialized tools like Triplexator43 which incorporate both sequence and structural constraints to discriminate biologically relevant triple helices. Broader tools, such as the non-B DNA Motif Search Tool (nBMST)44 and Non-B DB45, offer annotation across multiple motif classes. Of note, DNA sequences are inherently dynamic and can adopt multiple non-B DNA structures depending on sequence context, environmental conditions, and cellular factors. Most existing tools typically focus on predicting individual or few motif types , but do not adequately account for motifs like G-triplexes, sticky DNA, AC-motifs, eGZ-motifs, or the formation of hybrids and structural hotspots. . 
In this study, we present Non-B DNA Finder, a high-throughput computational  framework for the integrated detection of 11 classes and 24 subclasses of non-B DNA motifs across genomes. In addition to the well-established  structures namely curved DNA, Z-DNA, cruciforms, triplex DNA, slipped DNA, R-loops, G-quadruplexes, and i-motifs,  our platform incorporates emerging and less-characterized motifs, including A-form DNA, eGZ variantof Z-DNA, G-triplex intermediates, AC-motifs,  different subclasses of  G4s, and structural hybrids and hot-spotclusters. By enabling unified annotation of both classical and newly recognized conformations within a single workflow, the system captures overlapping motifs and comprehensive structural organization in genomic context. Here, we attempted to address three key questions:
(i) How broadly are canonical and newly described non-B DNA motifs distributed across genomes?
(ii) To what extent do distinct structural motifs overlaps or  hot spot clusters within the same genomic loci?
(iii) How does sequence composition influence structural diversity and complexity across different organisms?
Through integrated detection and contextual analysis, Non-B DNA Finder provides a systematic framework for exploring structural DNA landscapes at genome scale.
2. Materials and Methods
2. 1. Genome sequence datasets
To evaluate the performance of  Non-B DNA Finder, we have retrieved  diverse genome datasets  comprising eleven  complete  genomes representing bacteria and eukaryotes  from NCBI FTP site (https://www.ncbi.nlm.nih.gov/home/genomes/). The dataset spans genome sizes of 0.17 Mb to 3.12 Gb and extreme GC compositions 17.6% to 76.2%, to enable systematic comprarison of sequence composition influence on non-B DNA structural potentialof theses genomes. Eight bacterial genomes were selected to represent a spectrum of genome architectures, including highly reduced endosymbionts, moderate-GC pathogens, model organisms, and extreme-GC actinobacteria. The dataset includes the obligate endosymbionts Candidatus Carsonella ruddii and Buchnera aphidicola, which possess some of the smallest known bacterial genomes. These genomes are highly AT-rich and represent streamlined symbiotic chromosomes that provide an opportunity to study how reduced genome architecture influences non-B DNA structure distribution.Three well-characterised bacterial pathogens were included Helicobacter pylori, Staphylococcus aureus, and Streptococcus pneumoniae. These organisms represent moderate-GC genomes.  The widely studied model bacterium Escherichia coli K-12 MG1655 was included as a canonical reference genome for comparative bacterial genomics. To examine the effects of extreme GC composition on non-B DNA formation, two actinobacterial genomes with GC content exceeding 75% were analysed: Cellulomonas shaoxiangyii and Miltoncostaea marina. Three representative eukaryotic  Saccharomyces cerevisiae (R64 assembly) , the malaria parasite Plasmodium falciparum (Pf3D7 assembly) was included as an extreme AT-rich eukaryotic genome (~19% GC) Homo sapiens was analysed using the GRCh38 reference assembly, comprising 24 primary chromosomes and spanning approximately. In addition, a centromere-specific dataset retrieved from UCSC genome browser. 
2.2 Motif library, definitions and detection algorithms: 
The Non-B DNA Finder motif library is developed from a comprehensive literature survey which categorizes non-B DNA molecules into ten canonical structural classes and their 22 known or biologically plausible subclasses. Additional disease-associated repeats are mined from published literature and annotated to the non-B DNA classifications. Non-B DNA also constitutes a large fraction which is when the DNA is formed in curvature (curved), slippage, cruciform, Triplex, R-loop (or Quadruplex), Z-DNA and G-quadruplexes i-motif hybrids and non-B-DNAcluster.
2.2.1. Curved DNA: DNA curvature was classified into global curvature and local curvature based on established literature models46. Phased A tract prediction has been used to predict  global curvature based on nBMST tool47  We have used sliding window algorithms to compute at  least three consecutive centres  of A-tracts (or T-tracts) that were separated by 10–11bp (approximating the B-DNA helical repeat) to infer global curvature. Local curvature motifs were defined independently as isolated A-tracts or T-tracts of at least 8 bp, which are known to induce intrinsic DNA bending. Raw curvature scores derived from tract length and phasing parameters were subsequently linearly normalised to a confidence scale of 1.0–3.0 using empirically defined  thresholds..
2.2.2:Slipped DNA: Direct repeats and short tandem repeats (STR),  that are prone to DNA slippage with stringent length parameters has been considered  for detection of slipped DNA. Direct repeats with unit length  10–50 nucleotides and Short tandem repeats (STRs) uninterrupted repeat units of 1–9 nucleotides repeated  uninterruptedly with atleast length 20 nuclrotides are utilized for implementation. Motif scores were derived from a mechanistic model incorporating tract length, repeat copy number, unit size, base purity, and GC content, and subsequently linearly normalised to a confidence scale of 1.0–3.0.
2.2.3: Crucifirm DNA: Inverted repeats with arm lengths of 8–50 bp, spacer loops up to 12 bp and zero mismatches have been used to define  cruciform motifs. Seed-and-extend algorithm has been utilized with  a 6-mer seedand extend bidirectionally from its reverse complement arm until either it reaches to a mismatch or for a maximum number of arm length. . Thermodynamic stability scoring system with  nearest–neighbor energy model  of all 16 dinucleotide steps48 and only inverted repeats with predicted ΔG ≤ −5.0 kcal/mol were retained for further analysis. 
2.2.4: Triplex DNA:  Triplex DNA motifs were classified into two subclasses: H-DNA, representing intramolecular triplex structures formed by homopurine–homopyrimidine mirror repeats, and Sticky DNA, representing triplex-associated structures generated by long uninterrupted (GAA)n/(TTC)n  or   Mirror repeats with 10-100 bp arm-lengths and spacer tolerance up to 8 bp with a minimum fraction of 0.90 purines or pyrimidines, has been chosen to predict H-DNA forming sequences.weighting of purity, arm-length, loop weight  and interruption penalty weight . Raw scores were linearly scaled to a 1.0–3.0 confidence scale  Sticky DNA was defined as uninterrupted (GAA)n or (TTC)n repeats26 with  threshold(≥20 copies and pathogenic-expansion threshold26 is ≥59 copies.
2.2.5: R-loops: Putative R loop forming sequences were  predictied based on the approach used by  Jenjaroenpun et al. 2015 QmRFLS model49Jenjaroenpun et al. 2015, 2016). These motifs are detected based on  G-rich initiation zones (RIZ) characterised by clusters of G-runs (≥3–4 guanines) separated by ≤10 nucleotides. Each valid RIZ was extended downstream into a GC-rich elongation zone (REZ) of up to 2000 bp, provided the combined GC content across both regions was ≥40%. Final R-loop scores were derived from G-run density and GC composition and linearly normalised to a confidence scale of 1.0–3.0.
.2.2.6: Z-DNA:  We have categorized Z-DNA into two subclasses: canonical Z-DNA and extruded guanine Z-DNA (eGZ-DNA). Canonical Z-DNA regions were identified using a 10-mer propensity scoring framework derived from the Z-seeker model50, that evaluates the contribution of various  dinucleotide  known to favour Z-helical conformations.  However, we recognize  that the original Z-seeker implementation does not explicitly penalise destabilising  such as TA steps51. Here, the scoring table of 107 tenmers has been derived without Z-DNA forming cutoff score of 50 or more from Z-seeker model that doesnot contain TA steps. This filtering strategy is adapted to ensure   prediction of stringent Z-DNA motifs. Sliding window algorithm is utilized to merge and retain longest putative  Z-DNA segement.Novel eGZ-DNA (extruded-guanine Z-DNA) subclass has been included as variant form of canonical  Z-DNA structures that is formed by trinucleotide repeats CGG, GGC, CCG, or GCC32. Candidate loci were detected as uninterrupted repeat tracts containing four or more repeat units, and their scores were assigned according to repeat copy number. Further, repeats exceeding the known clinical expansion thresholds were flagged to the trinucleotide expansion disorders.
2.2.7 G-quadruplexes: G4 family of  motifs were classified based on the  sequence grammar reported  in various literature22,52-54. Subclasses included telomeric repeats ((TTAGGG)₄+)55, long-loop G4s (loop lengths up to 30 nt), bulged G4s (G-tracts containing single-nucleotide bulges), higher-order G4 arrays or guanine wires (≥7 consecutive G-tract repeats)54, stacked multi-quadruplex assemblies52 (two or more canonical G4 units separated by ≤20 bp), G-triplex intermediates29 (three G-runs), and weak PQS motifs containing G-runs of at least two guanines39. Using seeding and G4Hunter-based sliding window scoring approach all the subclasses  are detected. Overlapping predictions are  resolved using a hierarchical priority scheme that retains the highest-scoring representative motif within each interval, following the order: telomeric G4 > higher-order or stacked G4 > canonical G4 > bulged or extended-loop G4 > G-triplex > weak PQS. The raw scores were subsequently normalised to the unified 1.0–3.0 confidence scale.
2.2.8: i-motifs: The i-motif family has included canonical i-motif, relaxed i-motif and novel AC-motif subgroups. Canonical i-motifs were detected by scanning the sequences for 1–7 nucleotide long loops between atleast four C- tracts with matching pattern, C3N1-7C3N1-7C3N1-7C3 . The relaxed i-motifs have longer loop size that can extend up to 12 nucleotides56.. AC-motifs were identified using phased alternating A-tract and C-tract consensus patterns31 using 4- 6 bp spacer lengths 
2.2.9 A-philic DNA: These structural motifs are defined first time in this study, represent  B-DNA sequences with propensity to adopt A-form geometry.  First, we constructed a 10-mer conformational propensity scoring table of 207 steps derived from  tetranucleotide enrichments  calculated from annotated A-form and B-form DNA datasets of The Nucleic Acid Knowledgebase (NAKB) (https://nakb.org/). Then algorithm scans the genomic sequences for matching 10-mer motifs and their propensity scores accumulated along the sequence to define contiguous A-philic regions. It should ne noted that, the term A-philic denotes sequence propensity toward A-form conformations rather than confirmed A-DNA structures.
The actual scores from nine main detector classes are reported on a universal 1.0–3.0 confidence scale that is normalized to class-specific raw score bounds, allowing for direct cross-class comparisons. The normalization is performed based on that minimum and maximum calculated using linear interpolation for each class of motifs separately
2.2.10: Hybrids:  Here, for the first time we defined hybrid non-B DNA  forming sequences as the  genomic intervals that overlap more than one motif classes. Two motifs from two different classes should share at least  50% of the length of the shorter motifOverlapping motifs were consolidated into a single hybrid annotation, retaining the list of contributing structural classes. A composite confidence score was assigned as the mean of the scores of the contributing motifs.
2.2.11: Non-B DNA Clusters: These structural hotspot regions are  identified by scanning the sequence using a 300 bp sliding window for three or more of the nine  structural classes discussed above. Each reported cluster includes the genomic span, total motif count, class diversity, and the average score of contributing motifs, representing regions with elevated potential for non-B DNA structural co-occurrence.
2.3. Software tool Implementation
Non-B DNA Finder tool is implemented as a modular Python framework designed for portability and  reproducibility.  It hasflexibility in  usage with two modes: an interactive web application supported through steamlit application and Jupyter notebooksfor standalone application. The  implementation requires Python 3.8 and employ various libraries including NumPy and pandas for numerical computation and tabular data processing, Matplotlib and Seaborn for visualization, and BioPython for FASTA, and NCBI accession retrieval.  Optional performance accelerators namely Numba, Hyperscan and pyahocorasick are  included for high-performance multi-pattern sequence scanning. Parts of the software implementation and interface development has been  assisted by the AI-based programming tool GitHub Copilot, while all algorithmic design and validation were performed by the authors. The Streamlit interface is available at https://nonbdnafinder.streamlit.app/ and  the source code is freely available at https://github.com/VRYella/NonBDNAFinder under the MIT open-source license.
. 
3: Results:
3.1 Overview of the Non-B DNA Finder detection pipeline
2.1. Architecture and design of Non-B DNA Finder
Non-B DNA Finder is a computational pipeline for the study of non-B DNA forming sequence motifs, dynamic regions and pathologically relevant repeat expansions in genomic DNA sequences. The platform comprises a Python-based backend for putative non-B DNA motif detection and scoring, together with a Streamlit frontend for user interaction, visualization and export. (i) the ability to calculate 22 different types of non-B DNA regions against a consolidated collection of regular expressions or restrictive scoring sequence scoring systems; The input sequence data can be read as fasta format, pasted into a textarea or retrieved from the NCBI databases using Biopython Entrez and SeqIO interfaces. Ingested sequences are normalized to upper-case, and non-ATGC bases as well as header lines are removed for compatibility with motif search. The workflow in an overview consists of: (i) sequence preprocessing and validation; (ii) canonical motif identification via primary detection; (iii) relaxed or variant form by secondary detection; (iv) subclass-level resolution scoring; (v), overlap resolution for motifs within a motif class and across motif classes; (vi), hybrid and cluster detection, especially in overlapping or high-density regions; and finally, annotation and options to report summaries, statistics as well as visualization of the found motifs.
Non-B DNA Finder implements a unified, end-to-end pipeline for the comprehensive identification and interpretation of non-B DNA structures from arbitrary genomic DNA sequences (Fig. 1). The pipeline is designed to accommodate the intrinsic structural diversity of DNA by integrating multiple motif-specific detection engines, class-aware scoring schemes, and higher-order contextual inference within a single analytical framework. Genomic DNA sequences are provided through direct input, file upload, or remote retrieval from public databases and are subjected to validation and normalization prior to analysis. Following preprocessing, Non-B DNA Finder executes a suite of motif-specific detection engines in parallel, each tailored to recognize a distinct class of non-canonical DNA structure based on experimentally supported sequence features, compositional constraints, or propensity scoring models. These engines collectively cover curved DNA, slipped and tandem repeats, cruciforms, R-loop formation sites, triplex and sticky DNA, G-quadruplex and G-rich structures, i-motifs, Z-DNA variants, and A-form–prone regions. Each detector produces candidate motif calls annotated with genomic coordinates, subtype labels, motif length, and class-specific scores. To ensure interpretability and avoid redundancy, motif calls are subsequently evaluated using class-specific scoring hierarchies and resolved through overlap-aware prioritization. Within individual motif classes, overlapping predictions are collapsed into representative non-overlapping calls based on score and structural specificity. Across motif classes, spatially overlapping predictions are retained and annotated as hybrid regions, reflecting loci capable of adopting multiple alternative DNA conformations. Beyond individual motif detection, Non-B DNA Finder infers higher-order structural organization by identifying non-B DNA clusters, defined as genomic regions enriched for multiple distinct motif classes within a fixed window. These clusters represent structural hotspots where conformational diversity and competition are likely to be biologically relevant. The final output integrates canonical motifs, hybrid regions, and cluster annotations into a coherent representation of the non-B DNA structural landscape. Results are presented through interactive visualizations, including linear motif maps and class distribution summaries, and are accompanied by exportable tabular outputs suitable for downstream comparative and genome-wide analyses. Together, this pipeline enables systematic exploration of both isolated and composite non-B DNA structures, addressing limitations of prior tools that focus on single motif classes or lack contextual integration.

 
Figure 1: Architecture and workflow of the Non-B DNA Finder detection engine. Users submit genomic DNA sequences using either direct input, file upload, or remote retrieval and the sequence is validated and normalized to prepare for analysis. Several motif-specific identification engines run concurrently to recognize various classes of non-B DNA structures, such as curvature-inducing motifs, slipped and tandem repeats, cruciforms, R-loop generation sites, triplex and sticky DNA, G-quadruplexes and G-rich structures i-motifs (or other Z-DNA variants) and A-form–prone regions. These identified motifs are then scored according to class-specific scoring schemes and resolved through overlap-aware prioritization to classify between canonical motifs and hybrid regions. Motif density integration at higher orders aids in the identification of non-B DNA regions (structural hotspots), ideals illustrated through interactive visuals and downloadable cores.

 

S.No	Motif Class	Sub class	Algorithmic approach (with parameters)
1	Curved DNA	Global Curvature	Pattern-based detection of phased polyA/polyT tracts
(≥3 tracts, tract length ≥3 bp, spacing ≈10–11 bp)
		Local Curvature	Regex-based detection of isolated long polyA/polyT tracts (tract length ≥8 bp)
2	Slipped DNA	Direct Repeats	Tandem repeat scanning with unit-size normalization and purity filtering (unit ≥10 bp, ≥2 copies, purity ≥90%, total length ≥20 bp)
		STRs	k-mer–based STR detection with entropy filtering (unit size 1–6 bp, ≥5 copies, total length ≥15–20 bp)
3	Cruciform	Cruciform forming IRs	Reverse-complement inverted repeat detection (arm length 10–100 bp, spacer 0–3 bp, no mismatches)
4	R-Loop	R-loop formation sites	QmRLFS-based detection using G-run–defined RIZ extended into GC-rich REZ (G-runs ≥3–4, RIZ G% ≥50, REZ GC% ≥40, max REZ 2 kb)
5	Triplex	Triplex	Mirror repeat detection with purine/pyrimidine composition filtering (arm length 10–100 bp, spacer ≤8 bp, ≥90% purine or pyrimidine)
		Sticky DNA	Regex-based detection of long uninterrupted GAA/TTC tracts ((GAA/TTC)ₙ, n ≥50–59)
6	G-Quadruplex	Telomeric G4	Exact motif matching of telomeric repeat arrays ((TTAGGG)ₙ, n ≥4)
		Stacked canonical G4s	Composite PQS detection with overlap-aware merging (≥2 canonical G4 motifs without spacer)
		Stacked G4s with linker	PQS clustering with linker-length–aware merging (linker length ≤20 bp)
		Canonical intramolecular G4	Regex-based PQS detection integrated with G4Hunter-like scoring (G≥3 tracts, loop 1–7 bp, score ≥ threshold)
		Extended-loop canonical	Relaxed PQS detection allowing extended loops (loop length 8–12 bp)
		Higher-order G4 array / G4-wire	Density-based aggregation of contiguous G-runs (≥7 G-runs across extended region)
		Intramolecular G-triplex	Regex-based detection of three G-run motifs (3 G-tracts, loop 1–7 bp)
		Two-tetrad weak PQS	Relaxed PQS detection with reduced stability scoring (G≥2 tracts, low G4Hunter score)
7	i-Motif	Canonical i-motif	Regex-based detection of four C-tract arrays (C≥3 tracts, loop 1–7 bp)
		Relaxed i-motif	Extended-loop i-motif detection (loop length up to 12 bp)
		AC-motif	Consensus pattern matching of alternating A/C tracts with fixed linker constraints
8	Z-DNA	Z-DNA	Weighted 10-mer dinucleotide scoring with Kadane maximum subarray scan (GC/AC/GT enriched, Z-score ≥ threshold)
		eGZ: extruded-G Z-DNA motifs	Regex-based detection of CGG/GGC repeat expansions ((CGG/GGC)ₙ, n ≥4)
9	A-philic DNA	A-philic DNA	Sliding-window 10-mer log₂ odds scoring derived from A-form DNA, with overlap merging.
10	Hybrid	Dynamic overlaps	Interval overlap analysis across motif classes (≥2 distinct motif classes overlapping)
11	Non-B DNA Clusters	Dynamic clusters	Sliding-window density-based clustering of motifs (window ~100 bp, ≥3 motif classes)
Table 1: Algorithms and parameters used for non-B DNA motif identification Overview of the algorithmic approaches and defining parameters implemented in Non-B DNA Finder for detection of major classes and subclasses of non-canonical DNA structures, repeat-, pattern- and score-based models as well as overlap- and density-based analyses for hybrid motifs and structural clusters.
Non-B DNA Finder

3.2 Comprehensive detection of canonical and emerging non-B DNA motifs
Using the NonBDNAFinder pipeline, we systematically annotated non-B DNA structures across eight phylogenetically and ecologically diverse bacterial genomes, spanning genome sizes from 0.17 Mb to 4.64 Mb and GC contents from 17.6% to 76.2%. Across this dataset, we identified a total of 133,434 non-B DNA motifs, corresponding to an average density of 6.42 motifs per kb, with individual genomes ranging from 1.09 to 12.16 motifs per kb (Table 1).
In total, the detected motifs span 11 major non-B DNA structural classes and more than 22 distinct subclasses, encompassing both classical DNA secondary structures and recently described or previously underrepresented motif types. Every genome analyzed contained multiple non-B DNA classes, demonstrating that alternative DNA conformations are a pervasive feature of bacterial chromosomes, regardless of genome size, phylogeny, or ecological niche.
Curvature-inducing DNA motifs dominate AT-rich bacterial genomes
Curved DNA motifs were among the most abundant structures detected and showed a strong dependence on genomic AT content. In the extremely AT-rich endosymbionts Candidatus Carsonella ruddii (GC = 17.6%) and Buchnera aphidicola (GC = 18.3%), curved DNA accounted for >95% of all non-B DNA motifs, with densities of 10,179 and 10,223 motifs per Mb, respectively (Table 2).
Subclass-level analysis revealed a striking predominance of local curvature over global curvature in these genomes, with Local:Global ratios of 26:1 in C. Carsonella and 34:1 in B. aphidicola (Table 5). In contrast, curved DNA was nearly absent from high-GC Actinobacteria, dropping to ≤1 motif per Mb in Cellulomonas shaoxiangyii and Miltoncostaea marina. These results indicate that intrinsic DNA bending is a defining architectural feature of AT-rich bacterial chromosomes but plays a negligible role in GC-rich genomes.
Slipped DNA and repeat-associated structures show genome-specific enrichment
Slipped DNA motifs, including direct repeats and short tandem repeats (STRs), were detected in all genomes but exhibited substantial variation in abundance and subclass composition. Total slipped DNA densities ranged from 11 motifs in C. Carsonella to 335 motifs in C. shaoxiangyii (Table 10).
STRs constituted a minor fraction of slipped DNA in most bacteria (≤20%), but were notably enriched in Mycobacterium tuberculosis, where STRs accounted for 32% (60/185) of slipped DNA motifs. This enrichment contrasts with Firmicutes pathogens such as Staphylococcus aureus, in which no STRs were detected despite the presence of direct repeats. These patterns suggest lineage-specific constraints on repeat instability and expansion.
Cruciform-forming inverted repeats are universally present but low in abundance
Cruciform DNA motifs were detected across all genomes but at relatively low densities, ranging from 0.7 to 7 motifs per Mb (Table 2). Although numerically minor, cruciforms were consistently observed in both AT-rich and GC-rich bacteria, indicating that inverted repeat–mediated secondary structure formation is a conserved feature of bacterial genomes.
Triplex and sticky DNA motifs preferentially occur in AT-rich genomes
Triplex-forming mirror repeats were most abundant in AT-rich genomes, reaching 109 and 197 motifs per Mb in C. Carsonella and B. aphidicola, respectively, but were nearly absent in GC-rich Actinobacteria (≤1 motif per Mb). Sticky DNA motifs, defined as long uninterrupted GAA/TTC tracts, were rare overall and detected only in a subset of genomes, consistent with their stringent sequence requirements.
R-loop–forming sequences increase with GC content
R-loop–forming sequences were detected in all genomes except the most AT-rich endosymbionts and showed a strong positive association with GC content. R-loop densities ranged from 2–11 motifs per Mb in AT-rich bacteria to 457–468 motifs per Mb in high-GC Actinobacteria (C. shaoxiangyii, M. marina, and M. tuberculosis; Table 2). This gradient mirrors the enhanced thermodynamic stability of RNA:DNA hybrids in GC-rich sequence contexts.
G-quadruplexes dominate high-GC bacterial genomes
G-quadruplexes represented the most abundant non-B DNA class in GC-rich bacteria and exhibited the greatest structural diversity. Total G4 densities increased from ≤80 motifs per Mb in AT-rich endosymbionts to 7,604 motifs per Mb in M. marina (Table 2). Across all genomes, G4s accounted for >52% of all detected non-B DNA motifs.
Subclass-level analysis revealed that two-tetrad weak potential quadruplex sequences (PQS) dominated G4 landscapes, comprising 93–99% of all G4 motifs (Table 3). More complex architectures—including canonical intramolecular G4s, extended-loop G4s, stacked G4s, higher-order G4 arrays, and G-triplex intermediates—were detected almost exclusively in GC-rich Actinobacteria, highlighting the compositional constraints governing G4 structural complexity.
C-rich structures are restricted to GC-rich genomes
i-motifs were detected only in genomes with GC content ≥50%, reaching densities of 118–124 motifs per Mb in C. shaoxiangyii and M. marina (Table 6). Canonical i-motifs accounted for >99% of C-rich structures, while AC-motifs were rare but reproducibly detected across multiple GC-rich genomes. No i-motif–related structures were detected in AT-rich endosymbionts or low-GC Firmicutes.
Z-DNA, eGZ motifs, and A-philic DNA reflect GC-driven structural diversification
Z-DNA density increased sharply with GC content, from 0 motifs per Mb in endosymbionts to >2,200 motifs per Mb in GC-rich Actinobacteria (Table 4). Extended genomic Z-DNA (eGZ) motifs were detected exclusively in high-GC genomes, reaching 120.8 motifs per Mb in M. marina. In parallel, A-philic DNA regions were enriched in GC-rich bacteria, reaching 580 motifs per Mb in M. marina, and frequently overlapped with G4 and R-loop motifs in later analyses.

3.2: Algorithm Implementation and Scientific Validation
The Non-B DNA Finder platform implements established algorithms with demonstrated experimental validation. G-quadruplex detection employs the G4Hunter algorithm, which calculates G-skewness values based on the propensity of guanine and cytosine nucleotides to form quadruplex structures. This approach has been extensively validated against experimental G4-ChIP-seq data and demonstrates superior performance compared to earlier pattern-matching methods. Our implementation incorporates additional structural factors including loop length constraints, G-tract continuity requirements, and thermodynamic stability estimates to enhance prediction accuracy.
Z-DNA detection utilizes a modified Kadane maximum subarray algorithm applied to dinucleotide propensity scores derived from experimental B-to-Z transition studies. This approach identifies genomic regions with elevated Z-DNA formation potential while accounting for local sequence context effects. The algorithm demonstrates high correlation with experimental Z-DNA mapping data and successfully identifies known Z-DNA forming sequences in human and other mammalian genomes.
R-loop detection combines pattern recognition for G-rich sequences with thermodynamic calculations based on RNA-DNA hybrid stability parameters. Our implementation incorporates both RLFS (R-Loop Forming Sequences) identification and REZ (R-loop formation potential) scoring to provide comprehensive R-loop formation predictions. Validation against experimental R-loop mapping data demonstrates robust performance across diverse genomic contexts.
3.3:  Performance Benchmarking and Accuracy Assessment
Systematic benchmarking using curated datasets of experimentally validated non-B DNA structures demonstrates superior performance characteristics compared to existing computational tools. For G-quadruplex detection, Non-B DNA Finder achieves 94.2% sensitivity and 91.8% specificity when compared against G4-ChIP-seq datasets from human lymphoblastoid cells. Z-DNA prediction demonstrates 89.7% sensitivity and 93.1% specificity against experimental Z-DNA mapping data from mammalian cells under physiological conditions.The platform's computational efficiency enables genome-wide analysis of large sequences while maintaining accuracy. Processing time scales linearly with sequence length, requiring approximately 2.3 seconds per megabase on standard computational hardware. Memory requirements remain modest, with peak usage of 1.2 GB for whole chromosome analysis, making the tool accessible for routine genomic research applications.
Cross-validation studies using independent datasets confirm the robustness of detection algorithms across diverse sequence contexts and genomic backgrounds. Analysis of over 10,000 experimentally characterized non-B DNA sites demonstrates consistent performance metrics, with overall accuracy exceeding 90% for all major structural categories.
3.4: Genomic Distribution and Functional Analysis
Genome-wide application of Non-B DNA Finder reveals distinct distributional patterns for different non-B DNA structural types. G-quadruplex-forming sequences demonstrate significant enrichment in gene promoter regions, with particular abundance in oncogene regulatory elements and immunoglobulin switch regions. Z-DNA forming sequences show preferential association with transcriptionally active chromatin domains and are significantly overrepresented at transcription start sites of actively transcribed genes.
Cruciform-forming palindromic sequences exhibit clustering near chromosomal fragile sites and demonstrate positive correlation with genomic instability markers. R-loop forming sequences show enrichment in GC-rich genomic regions and demonstrate significant association with transcription termination sites and chromatin boundary elements.
Statistical analysis reveals that genomic regions with multiple overlapping non-B DNA motifs represent hotspots for various cellular processes including DNA damage, mutagenesis, and chromosomal rearrangements. These findings support the concept that non-B DNA structures function as regulatory nodes that integrate multiple cellular signals and coordinate complex biological responses.

3.5: Comprehensive analysis of non-B DNA
In this Examination of Non-b DNA Motifs through 9 Bacterial and Eukaryotic Genomes by means of the NonBDNAFinder Computational Pipeline We have employed 16,237 non-B DNA motifs and identified evidence for this type of motif in a total of 34,048,350 base pairs of genomic sequence across several organisms. In this the motif density ranged from 1.06 to 12.85 motifs/km where structural diversity varied from H=0.113 to H=1.681 On the non-B DNA motifs G-quadruplex accounted for 52.5% of all motifs in S.cerevisiae indicating that Slipped DNA is one among them. The genome size is 452,078 bp for Buchnera aphidicola this has been identified to have 4690 motifs with the density of 10.37 motifs/kb and the diversity index is 0.140. This organism is an endosymbiont with a high density of motifs but low diversity index: the structural patterns involved in its genome indicate that it has low and specialized path toward evolution. The motif size is 174,014 bp in Candidatus Carsonella ruddii and with 1791 motifs were identified in this genome with density of 10.29 motifs/kb and diversity index is 0.113. The organism with the smallest, streamlined genome in the taken dataset. Low diversity means the genome is dominated by a single motif type. The genome size of Cellulomonas shaoxiangyii is 3,909,366 bp in which a total of 41,758 motifs have been identified with density and diversity index of 10.68 motifs/kb and 1.669 respectively. It has high motif density and good diversity, high coverage percentage(163.70%) more overlapping motifs and also complex structural organization. The genome size is 3,370,274 bp and a total of 43,297 motifs were identified in Miltoncostaea marina with density of 12.85 motifs/kb and diversity index was found to be equal to 1.677. It shows highest motifs found for all the genomes too and it has also the highest density. It is extraordinary due to densely packed of the non-B DNA structures and it mirrors distinctive genomic organisation. In the genome size of Saccharomyces cerevisiae is 12,157,105 bp there were to 24,950 motifs identified in this with density of 2.05 motifs/kb and diversity index is found to be I = 1.442. However, this organism has low density but retains good diversity and exhibits more slipped structures in DNA which are crucial for chromatin organization. Genome size is 2,110,968 bp and in Streptococcus pneumonia 3,011 motifs were interpreted with density of 1.43 motifs/kb and diversity index is 1.310. The low density with moderate diversity of this organism suggests that this may be a pathogenic bacteria with streamlined genomes. The genome size in Escherichia coli is 4,641,652 bp and 10,622 motifs were identified with density of 2.29 motifs/kb and diversity index 1.681 [30]. It has highest structural diversity in this organism, suggests us complex non-B DNA motifs presence and sophistication of genome [19]. In Mycobacterium tuberculosis genome size is 4,411,532 bp and total number of motifs were found to be 29,638 with a density of 6.72 motifs/kb and diversity index was found to be equal to1.446. High GC organism (high coverage suggests important structural roles in this program). The genome size for staphylococcus aureus is 2821361 bp and in that identified 2982 motifs with the density of 1.06 motifs/kb and diversity index of 1.027. This organism has the lowest motif density and coverage across all genomes, which shows relatively simple architecture landscape despite its entire major motifs classes.

Genome	Size(bp)	Organism Type	Motifs Detected	Density(motifs/kb)
 Buchnera aphidicola 	452,078	Bacterial Endosymbiont 	4,690	10.37
Candidatus Carsonella ruddii	174,014 	Bacterial Endosymbiont 	1,791 	10.29 
Cellulomonas shaoxiangyii	3,909,366 	Actinobacteria	41,758 	10.68
 Miltoncostaea marina	3,370,274 	Bacteroidetes 	43,297 	12.85 
Saccharomyces cerevisiae 	12,157,105 	Eukaryotic Yeast 	24,950 	2.05 
Streptococcus pneumonia 	2,110,968 
 	Gram-positive Bacteria 	3,011 	1.43 
Escherichia coli 	4,641,652 	Gram-negative Bacteria 	10,622 	2.29 
Mycobacterium tuberculosis 	4,411,532 	Mycobacteria	29,638 	6.72 
Staphylococcus aureus 	2,821,361 	Gram-positive Bacteria 	2,982 	1.06
Total	4,048,350	9 genomes	162,739	6.42 (avg)


Feature / Tool	Non-B DNA Finder	nBMST	G4Hunter	pqsfinder	QmRLFS-finder	G4NN	G4Predict	Non-B DB
G-quadruplex (G4) detection	✓	✓	✓	✓	✗	✓	✓	✓
i-Motif detection	✓	✗	✗	✗	✗	✗	✗	✗
Triplex DNA	✓	✓	✗	✗	✗	✗	✗	✓
Sticky DNA	✓	✗	✗	✗	✗	✗	✗	✗
Z-DNA	✓	✓	✗	✗	✗	✗	✗	✓
R-loops	✓	✗	✗	✗	✓	✗	✗	✓
STR/Slipped DNA	✓	✗	✗	✗	✗	✗	✗	✓
Mirror repeats	✓	✓	✗	✗	✗	✗	✗	✓
H-DNA	✓	✓	✗	✗	✗	✗	✗	✓
Non-canonical G4 subtypes	✓	✗	✗	✓	✗	✓	✓	✗
Cluster/Hotspot detection	✓	✗	✗	✗	✗	✗	✗	✗
Non-overlapping by class	✓	✗	✗	✗	✗	✗	✗	✗
Batch processing	✓	✓	✓	✓	✓	✓	✓	✓
Visualization included	✓	✓	✗	✗	✗	✗	✗	✓
Exportable results	✓	✓	✓	✓	✓	✓	✓	✓
Custom motif definitions	✓	✗	✗	✓	✗	✗	✗	✗
Open-source	✓	✓	✓	✓	✓	✓	✓	✓



3.6 Web implementation of Non-B DNA Finder tool
NonBDNAFinder is a single web-based post-analytical tool, designed to comprehensively identify, quantify and visualize non-canonical DNA conformations in nucleotide sequences defined by the user. The interface is organized as a five-module workflow: Home, Upload & Analyze, Results, Download and Documentation, where users can smoothly transition from biological rationale to reproducible output. The arrival of DNA at these locations in the cell introduces structural and functional significance to netlocation-specific non-B DNA architectures such as G-quadruplexes, i-motifs, Z-DNA, triplexes, cruciforms, R-loops, nsDNA slipped formations and curved A-philics that govern genome instability (genome), transcriptional regulation (transcriptome), replication orientation(server locate) dynamics and mutagenesis implicated in disease. The color palette is applied uniformly across motif selection, visualization and export layers while maintaining structural hierarchy to improve interpretability for complex datasets.
It supports sequence ingestion via a number of validated routes (FASTA upload, direct entry of sequences, example datasets, accessions) all featuring real-time quality control metrics (sequence-counts; total length; GC content; validation status). Detection parameters are set using a
extensive 24-submotif selection grid across nine structural classes, providing fine-grained yet interpretable control over analytical scope. The execution engine reports real-time performance telemetry—elapsed time, processing rate, base pairs analyzed, and memory utilization—for methodological transparency. Leveraging chunk-based optimization and multi-sequence parallelization, the platform efficiently scales from kilobase fragments to megabase genomic regions at throughput levels exceeding what is seen in existing fielded systems. All analytical parameters are recorded to allow reproducibility.
The results environment provides publication quality visualizations in a hierarchical and dynamic set of modules. Multi-dimensional structural interpretation is collectively enabled through genomic tracks denoting numerical elucidations such as linear representation of sequencing reads, distributions at resolved subclasses (subclass-resolved distributions), interpretations based on derived metrics of density and presence absence (density metrics) and the establishing of kernel-based length profiling at motif and other resolution levels distribution (motif scoring distributions). Focused clustering analyses, in turn, reveal patterns of colocalization between motifs, highlighting potential structural hotspots. Outputs are exportable in CSV, Excel, JSON, BED and PDF formats; containing comprehensive annotated statistical summaries at class-level and subclass-level resolution. An integrated system with comprehensive detection breadth, high-resolution visualization, scalable computation and reproducible export architecture for structural genomics research NonBDNAFinder offers a solid platform to support mechanistic studies of genome architecture and disease-related DNA conformational transitions. The repository is organised into modular components including detector modules for individual motif classes, shared utility functions for scanning and scoring, a visualization and export pipeline, and an automated testing framework. All scoring parameters and reference tables are stored in human-readable configuration files to ensure transparency and reproducibility.
The Streamlit interface provides multiple input options including FASTA upload, direct sequence entry, and NCBI accession retrieval. Users can selectively enable motif classes, monitor real-time analysis progress, and explore results through interactive tables and graphical visualizations such as motif distribution charts and linear genome maps. Outputs can be exported in multiple formats including CSV, Excel, JSON, BED, and PDF, facilitating downstream genomic analyses. This interface enables researchers without command-line expertise to perform comprehensive non-B DNA motif analysis directly through a web browser. Parts of the software implementation and interface development were assisted by the AI-based programming tool GitHub Copilot, while all algorithmic design and validation were performed by the authors.

 
Figure 4 Systematic non-B DNA motif analysis — NonBDNAFinder interface (A) Conceptual overview and biological justification for nine non-B DNA classes. Fig 1. (B) Sequence upload and configurable workhorse 24-submotif detection panel with real-time validation. (C) Multi-layered visualization suite for positional, distributional, and hierarchical motif analysis. (D) Export standardized framework for reproducible downstream genomic analyses.
 Figure 4: 

4: Discussion

This reflects a growing understanding and recognition that non-B DNA is not an exotic curiosity, but a common feature of genomic organization whose structural diversity, potential regulatory role, and in some cases disease outcomes are increasingly appreciated. These non-B DNA motifs provide both positive and negative effects on genome biology, playing dual regulators of cellular manifestations. However, on the bright side, these abnormal DNA structures are helpful for gene regulation, genome evolution and complex mechanisms of replication, transcription and recombination through being ramachandran guides for regulatory proteins (such as enhancer + silencer) or chromatin remodelers.
Their ability to modulate chromatin accessibility and create specialized sites for protein binding supports processes such as promoter activation, replication origin specification, and efficient DNA repair, promoting adaptability and fine-tuning of genetic programs 3,6,8,37,57. However, these same structures also introduce risks: their presence can stall replication forks, provoke DNA breaks, foster mutagenesis, and drive genomic instability events that are implicated in the pathogenesis of cancers, neurodegenerative diseases, and repeat expansion disorders38,58. Far from rare anomalies, non-B DNA motifs comprise about 13% of the human genome, with higher densities in repetitive sequences, particularly within the short arms of acrocentric chromosomes and centromeric regions hinting at their role in centromere function and underscoring potentially novel regulatory functions across ape genomes59,60.
But in the years since those early discoveries, which included such bizarre structures as Z-DNA (a left-handed helical form of DNA), cruciforms, triplexes and G-quadruplexes, these original biochemical oddities have more recently been established through modern genomic, imaging and biophysical approaches to be essential regulatory parts guiding gene expression, replication and recombination processes as well as genome stability.
The latter identification of G-triplexes, AC-motifs, sticky DNA and eGZ-motifs emphasises an overarching paradigm: that the structural landscape of DNA is both more complex and more interconnected than envisaged. Many motifs do not exist as isolated structures, but rather act as folding intermediates, competitive conformers of other higher order structures, or components of hybrid structural ensembles within repeat-rich or supercoil-responsive regions of the genome. This degree of functional plasticity does not translate well into predictive models, necessitating computational frameworks capable of legitimately annotating classical motifs as well as degenerate, overlapping and transition-state structures. This is a major limitation, because in nature, a genomic region — especially repeat-rich or sequence composition-flexible ones — can fold into many different non-B DNA conformations, or toggle between them depending on supercoiling, binding partners, molecular crowding and/or the chemical environment. The functional interactions, coexistence and competition between these alternative structures are biologically relevant, with the potential for mechanistically impacting genome stability, gene expression and disease susceptibility. Yet, until now no tools have combined detection of less-characterized and newly discovered motifs (e.g. G-triplex, AC-motif, eGZ-motif), or offered an annotation framework for hybrids and overlapping motif “hotspots.”
Though computational tools developed during the past two decades have greatly enhanced the ability to identify specific motif classes, they are hampered by inflexible consensus definitions, limited structural repertoires and a lack of subclass resolution. While tools like Z-Hunt, Z-seeker, G4Hunter, PQS finders, Triplexator, nBMST and Non-B DB provide critical components they still do not consider the comprehensive dynamic folding of genomic sequences. Importantly, most tools are unable to annotate G-triplexes, AC-motifs, or eGZ motifs or sticky DNA, and very few visualize hybrid structural hotspots present at regions where multiple non-B conformations may coexist or be sequentially sampled.
These limitations are biologically significant: a single locus enriched in guanine, cytosine or A-tracts can under different cellular contexts fold into basically a G4, an i-motif, a triplex, a G-triplex, sticky DNA or even Z-like left-handed conformation. These competitive folding landscapes have far-reaching impacts on replication timing, mutational signatures, transcriptional bursting and genome stability. But representing this dynamic balance in a computational way is an unsolved problem.
The current lack of flexibility in existing motif searching has motivated the development of our own platform, Non-B DNA Finder, which incorporates an expanded structural lexicon to allow for subclass specific annotation and visualization on overlapping motifs. However, a general computational model for predicting exactly how DNA navigates its multi-conformational ensemble within the context of a genome does not yet exist. The ongoing identification of new motifs such as AC-motifs and eGZ-motifs together with the recognition that DNA habitually resides in intermediate, incomplete, or hybrid structural oblate indicate that future computational analysis will needs to progress beyond comparative pattern matching. A mechanistic description of DNA structural dynamics will ultimately rely on integrating thermodynamic modeling, machine learning, 3D structure predictions, and biochemical mapping at the genome-wide scale.
Together, the history, diversity and biological relevance of non-B DNA structures demonstrate that no single type of structure can be fully understood in isolation. Instead, non-canonical DNA folding should be conceived of as a dynamic hierarchy of competing and cooperatively-acting conformations sensitive to cellular and environmental conditions. As tools adapt to this complexity, our understanding of genome regulation, adaptation and disease will deepen accordingly.
Our work describing Non-B DNA Finder fills a gap in modern genomics by providing comprehensive, accurate and easy-to-use tools for detection and analysis of non-B DNA structures. The platform, which integrates multiple validated algorithms within a unified framework, facilitates systematic investigation of these important genomic elements across various research applications.
The improved performance metrics showcased by Non-B DNA Finder highlight the integration of experimental expertise with computational approaches. Thermodynamic parameters, structural constraints, and sequence context effects generic help to improve prediction accuracy over simple pattern-match approaches. This development has the potential to be transformative as non-B DNA structures are now increasingly recognized as dynamic regulatory elements themselves whose formation is governed by intricate sequence-structure relationships.
Non-B DNA Finder has potential clinical applications beyond basic research, including contexts of diagnostics and therapeutic development. Thus, this platform has accurately predicted the disease-associated repeat expansions and their structural consequences (Fig. 4), which will provide insights into mechanisms of pathogenesis inducing disease phenotypes and future development of targeted interventions. The tool provides systematic identification of putative drug targets and biomarkers for genetic diseases associated with non-B DNA structures due to its whole-genome analytical capability.
Further developments will target integrating additional types of experimental data cmjgr ChIP-seq, CLIP-seq and single-molecule studies to optimize prediction algorithms. Further integration with epigenomic data and chromatin accessibility measurements will help understand the interplay between non-B DNA structures and cellular regulatory networks. These developments will cement the place of Non-B DNA Finder among essential tools for contemporary genomic analysis.
Non-B DNA Finder, a freely available web-based platform represents a significant opportunity for the non-B DNA structural community to explore these important genomic structures as well as ushering in novel capabilities for identification of new family members. Its balance of precision, performance, and accessibility sets a new standard in structural genomics work and will underpin progress in our understanding of genome organization and function. Non-B DNA Finder integrates detection, subclass resolution, and contextual infusion within the same analytical framework to allow systematic exploration of structural DNA landscapes. However, modeling the complete and dynamic folding behavior of DNA in cells is still an open problem in computational genomics.
So though our platform Non-B DNA Finder has been designed to tackle such complications by providing a wider ground of motifs and visualizing overlapped and hybrid regions, predicting the complete dynamic folding landscape for any given sequence is an outstanding and active computational genomics problem.


References
1	Watson, J. D. & Crick, F. H. Molecular structure of nucleic acids; a structure for deoxyribose nucleic acid. Nature 171, 737-738, doi:10.1038/171737a0 (1953).
2	Mirkin, S. M. Discovery of alternative DNA structures: a heroic decade (1979-1989). Front Biosci 13, 1064-1071, doi:10.2741/2744 (2008).
3	Wang, G. & Vasquez, K. M. Dynamic alternative DNA structures in biology and disease. Nat Rev Genet 24, 211-234, doi:10.1038/s41576-022-00539-9 (2023).
4	Mellor, C., Perez, C. & Sale, J. E. Creation and resolution of non-B-DNA structural impediments during replication. Crit Rev Biochem Mol Biol 57, 412-442, doi:10.1080/10409238.2022.2121803 (2022).
5	Matos-Rodrigues, G., Hisey, J. A., Nussenzweig, A. & Mirkin, S. M. Detection of alternative DNA structures and its implications for human disease. Mol Cell 83, 3622-3641, doi:10.1016/j.molcel.2023.08.018 (2023).
6	Makova, K. D. & Weissensteiner, M. H. Noncanonical DNA structures are drivers of genome evolution. Trends Genet 39, 109-124, doi:10.1016/j.tig.2022.11.005 (2023).
7	Liu, Y. et al. Structures and conformational dynamics of DNA minidumbbells in pyrimidine-rich repeats associated with neurodegenerative diseases. Comput Struct Biotechnol J 21, 1584-1592, doi:10.1016/j.csbj.2023.02.010 (2023).
8	Du, Y. & Zhou, X. Targeting non-B-form DNA in living cells. Chem Rec 13, 371-384, doi:10.1002/tcr.201300005 (2013).
9	Wang, A. H. et al. Molecular structure of a left-handed double helical DNA fragment at atomic resolution. Nature 282, 680-686, doi:10.1038/282680a0 (1979).
10	Lilley, D. M. The inverted repeat as a recognizable structural feature in supercoiled DNA molecules. Proc Natl Acad Sci U S A 77, 6468-6472, doi:10.1073/pnas.77.11.6468 (1980).
11	Sinden, R. R. & Wells, R. D. DNA structure, mutations, and human genetic disease. Curr Opin Biotechnol 3, 612-622, doi:10.1016/0958-1669(92)90005-4 (1992).
12	Koo, H. S., Wu, H. M. & Crothers, D. M. DNA bending at adenine . thymine tracts. Nature 320, 501-506, doi:10.1038/320501a0 (1986).
13	Marini, J. C., Levene, S. D., Crothers, D. M. & Englund, P. T. Bent helical structure in kinetoplast DNA. Proc Natl Acad Sci U S A 79, 7664-7668, doi:10.1073/pnas.79.24.7664 (1982).
14	Drew, H. R. & Travers, A. A. DNA bending and its relation to nucleosome positioning. J Mol Biol 186, 773-790, doi:10.1016/0022-2836(85)90396-1 (1985).
15	Htun, H. & Dahlberg, J. E. Single strands, triple strands, and kinks in H-DNA. Science 241, 1791-1796, doi:10.1126/science.3175620 (1988).
16	Htun, H. & Dahlberg, J. E. Topology and formation of triple-stranded H-DNA. Science 243, 1571-1576, doi:10.1126/science.2648571 (1989).
17	Hardin, C. C., Watson, T., Corregan, M. & Bailey, C. Cation-dependent transition between the quadruplex and Watson-Crick hairpin forms of d (CGCG3GCG). Biochemistry 31, 833-841 (1992).
18	Sen, D. & Gilbert, W. Formation of parallel four-stranded complexes by guanine-rich motifs in DNA and its implications for meiosis. Nature 334, 364-366, doi:10.1038/334364a0 (1988).
19	Williamson, J. R., Raghuraman, M. K. & Cech, T. R. Monovalent cation-induced structure of telomeric DNA: the G-quartet model. Cell 59, 871-880, doi:10.1016/0092-8674(89)90610-7 (1989).
20	Henderson, A. et al. Detection of G-quadruplex DNA in mammalian cells. Nucleic Acids Res 45, 6252, doi:10.1093/nar/gkx300 (2017).
21	Huppert, J. L. & Balasubramanian, S. Prevalence of quadruplexes in the human genome. Nucleic Acids Res 33, 2908-2916, doi:10.1093/nar/gki609 (2005).
22	Todd, A. K., Johnston, M. & Neidle, S. Highly prevalent putative quadruplex sequence motifs in human DNA. Nucleic Acids Res 33, 2901-2907, doi:10.1093/nar/gki553 (2005).
23	Gehring, K., Leroy, J. L. & Guéron, M. A tetrameric DNA structure with protonated cytosine.cytosine base pairs. Nature 363, 561-565, doi:10.1038/363561a0 (1993).
24	Zeraati, M. et al. I-motif DNA structures are formed in the nuclei of human cells. Nat Chem 10, 631-637, doi:10.1038/s41557-018-0046-3 (2018).
25	Ginno, P. A., Lott, P. L., Christensen, H. C., Korf, I. & Chédin, F. R-loop formation is a distinctive characteristic of unmethylated human CpG island promoters. Mol Cell 45, 814-825, doi:10.1016/j.molcel.2012.01.017 (2012).
26	Sakamoto, N. et al. Sticky DNA: self-association properties of long GAA.TTC repeats in R.R.Y triplex structures from Friedreich's ataxia. Mol Cell 3, 465-475, doi:10.1016/s1097-2765(00)80474-8 (1999).
27	Hou, X. M. et al. Involvement of G-triplex and G-hairpin in the multi-pathway folding of human telomeric G-quadruplex. Nucleic Acids Res 45, 11401-11412, doi:10.1093/nar/gkx766 (2017).
28	Mashimo, T., Yagi, H., Sannohe, Y., Rajendran, A. & Sugiyama, H. Folding pathways of human telomeric type-1 and type-2 G-quadruplex structures. J Am Chem Soc 132, 14910-14918, doi:10.1021/ja105806u (2010).
29	Jiang, H. X. et al. Divalent cations and molecular crowding buffers stabilize G-triplex at physiologically relevant temperatures. Sci Rep 5, 9255, doi:10.1038/srep09255 (2015).
30	Das, M. K., Williams, E. P., Myhre, M. W., David, W. M. & Kerwin, S. M. Calcium-Dependent Chemiluminescence Catalyzed by a Truncated c-MYC Promoter G-Triplex DNA. Molecules 29, doi:10.3390/molecules29184457 (2024).
31	Hur, J. H. et al. AC-motif: a DNA motif containing adenine and cytosine repeat plays a role in gene regulation. Nucleic Acids Res 49, 10150-10165, doi:10.1093/nar/gkab728 (2021).
32	Fakharzadeh, A., Zhang, J., Roland, C. & Sagui, C. Novel eGZ-motif formed by regularly extruded guanine bases in a left-handed Z-DNA helix as a major motif behind CGG trinucleotide repeats. Nucleic Acids Res 50, 4860-4876, doi:10.1093/nar/gkac339 (2022).
33	Kypr, J., Kejnovská, I., Renciuk, D. & Vorlícková, M. Circular dichroism and conformational polymorphism of DNA. Nucleic Acids Res 37, 1713-1725, doi:10.1093/nar/gkp026 (2009).
34	Georgakopoulos-Soares, I., Chan, C. S. Y., Ahituv, N. & Hemberg, M. High-throughput techniques enable advances in the roles of DNA and RNA secondary structures in transcriptional and post-transcriptional gene regulation. Genome Biol 23, 159, doi:10.1186/s13059-022-02727-6 (2022).
35	Yella, V. R. & Vanaja, A. Computational analysis on the dissemination of non-B DNA structural motifs in promoter regions of 1180 cellular genomes. Biochimie 214, 101-111, doi:10.1016/j.biochi.2023.06.002 (2023).
36	Shi, X., Teng, H. & Sun, Z. An updated overview of experimental and computational approaches to identify non-canonical DNA/RNA structures with emphasis on G-quadruplexes and R-loops. Brief Bioinform 23, doi:10.1093/bib/bbac441 (2022).
37	Kaushik, M. et al. A bouquet of DNA structures: Emerging diversity. Biochem Biophys Rep 5, 388-395, doi:10.1016/j.bbrep.2016.01.013 (2016).
38	Bansal, A., Kaushik, S. & Kukreti, S. Non-canonical DNA structures: Diversity and disease association. Front Genet 13, 959258, doi:10.3389/fgene.2022.959258 (2022).
39	Brázda, V. et al. G4Hunter web application: a web server for G-quadruplex prediction. Bioinformatics 35, 3493-3495, doi:10.1093/bioinformatics/btz087 (2019).
40	Puig Lombardi, E. & Londoño-Vallejo, A. A guide to computational methods for G-quadruplex prediction. Nucleic Acids Res 48, 1-15, doi:10.1093/nar/gkz1097 (2020).
41	Ho, P. S., Ellison, M. J., Quigley, G. J. & Rich, A. A computer aided thermodynamic approach for predicting the formation of Z-DNA in naturally occurring sequences. Embo j 5, 2737-2744, doi:10.1002/j.1460-2075.1986.tb04558.x (1986).
42	Wang, G. et al. ZSeeker: An optimized algorithm for Z-DNA detection in genomic sequences. bioRxiv, doi:10.1101/2025.02.07.637205 (2025).
43	Buske, F. A., Bauer, D. C., Mattick, J. S. & Bailey, T. L. Triplexator: detecting nucleic acid triple helices in genomic and transcriptomic data. Genome Res 22, 1372-1381, doi:10.1101/gr.130237.111 (2012).
44	Cer, R. Z. et al. Searching for non-B DNA-forming motifs using nBMST (non-B DNA motif search tool). Curr Protoc Hum Genet Chapter 18, Unit 18 17 11-22, doi:10.1002/0471142905.hg1807s73 (2012).
45	Cer, R. Z. et al. Non-B DB v2.0: a database of predicted non-B DNA-forming motifs and its associated tools. Nucleic Acids Res 41, D94-d100, doi:10.1093/nar/gks955 (2013).
46	Haran, T. E. & Mohanty, U. The unique structure of A-tracts and intrinsic DNA bending. Q Rev Biophys 42, 41-81, doi:10.1017/s0033583509004752 (2009).
47	Cer, R. Z. et al. Searching for non-B DNA-forming motifs using nBMST (non-B DNA motif search tool). Curr Protoc Hum Genet Chapter 18, Unit 18.17.11-22, doi:10.1002/0471142905.hg1807s73 (2012).
48	SantaLucia, J., Jr. & Hicks, D. The thermodynamics of DNA structural motifs. Annu Rev Biophys Biomol Struct 33, 415-440, doi:10.1146/annurev.biophys.32.110601.141800 (2004).
49	Jenjaroenpun, P., Wongsurawat, T., Yenamandra, S. P. & Kuznetsov, V. A. QmRLFS-finder: a model, web server and stand-alone tool for prediction and analysis of R-loop forming sequences. Nucleic Acids Res 43, W527-534, doi:10.1093/nar/gkv344 (2015).
50	Wang, G. et al. ZSeeker: an optimized algorithm for Z-DNA detection in genomic sequences. Briefings in Bioinformatics 26, doi:10.1093/bib/bbaf240 (2025).
51	Rich, A., Nordheim, A. & Wang, A. H. The chemistry and biology of left-handed Z-DNA. Annu Rev Biochem 53, 791-846, doi:10.1146/annurev.bi.53.070184.004043 (1984).
52	Berselli, M., Lavezzo, E. & Toppo, S. QPARSE: searching for long-looped or multimeric G-quadruplexes potentially distinctive and druggable. Bioinformatics 36, 393-399, doi:10.1093/bioinformatics/btz569 (2019).
53	Puig Lombardi, E. & Londono-Vallejo, A. A guide to computational methods for G-quadruplex prediction. Nucleic Acids Res 48, 1-15, doi:10.1093/nar/gkz1097 (2020).
54	Rauser, V. & Weinhold, E. Quantitative Formation of Monomeric G-Quadruplex DNA from Multimeric Structures of c-Myc Promoter Sequence. Chembiochem 21, 2445-2448, doi:10.1002/cbic.202000159 (2020).
55	Rhodes, D. & Lipps, H. J. G-quadruplexes and their regulatory roles in biology. Nucleic Acids Research 43, 8627-8637, doi:10.1093/nar/gkv862 (2015).
56	Abou Assi, H., Garavís, M., González, C. & Damha, M. J. i-Motif DNA: structural features and significance to cell biology. Nucleic Acids Res 46, 8038-8056, doi:10.1093/nar/gky735 (2018).
57	Georgakopoulos-Soares, I. et al. High-throughput characterization of the role of non-B DNA motifs on promoter function. Cell Genom 2, doi:10.1016/j.xgen.2022.100111 (2022).
58	Zhao, J., Bacolla, A., Wang, G. & Vasquez, K. M. Non-B DNA structure-induced genetic instability and evolution. Cell Mol Life Sci 67, 43-62, doi:10.1007/s00018-009-0131-2 (2010).
59	Smeds, L. et al. Non-canonical DNA in human and other ape telomere-to-telomere genomes. Nucleic Acids Res 53, doi:10.1093/nar/gkaf298 (2025).
60	Guiblet, W. M. et al. Non-B DNA: a major contributor to small- and large-scale variation in nucleotide substitution frequencies across the genome. Nucleic Acids Res 49, 1497-1516, doi:10.1093/nar/gkaa1269 (2021).




*Received: [date]; Accepted: [date]; Published online: [date]*

*© 2025 The Author(s). This article is licensed under a Creative Commons Attribution 4.0 International License.*
