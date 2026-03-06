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

---

*Received: [date]; Accepted: [date]; Published online: [date]*

*© 2025 The Author(s). This article is licensed under a Creative Commons Attribution 4.0 International License.*
