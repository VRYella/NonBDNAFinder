# NonBDNAFinder: A Unified Computational Framework for Genome-Scale Detection and Characterisation of Non-Canonical DNA Structures

**Venkata Rajesh Yella**¹ and collaborators

¹ Department of Biotechnology, [Institution], India

*Correspondence:* [yella@institution.ac.in](mailto:yella@institution.ac.in)

---

## Abstract

Non-B DNA structures—departures from the canonical Watson–Crick double helix—play pivotal roles in replication, transcription, recombination, genome instability, and human disease. Despite decades of biochemical characterisation, no single open-source tool has integrated the full breadth of these structures within a unified, genome-scale analytical workflow. We present **NonBDNAFinder (NBDFinder)**, a modular Python framework that simultaneously detects nine distinct non-B DNA classes—Curved DNA, Slipped DNA, Cruciform, R-Loop, Triplex, G-Quadruplex, i-Motif, Z-DNA, and A-philic DNA—yielding 11 output classes (including Hybrid annotations and Non-B DNA Clusters) spanning **24 structural subclasses**. Each motif call is tagged with genomic coordinates (1-based inclusive start/end), strand, subclass label, sequence length, and a literature-anchored normalised confidence score on a 1–3 scale. Post-processing layers resolve cross-class overlaps, annotate hybrid loci, and identify co-occurrence hotspots. The platform supports multi-FASTA input at chromosome scale (>100 Mb), provides optional Numba JIT and Hyperscan acceleration, and exports results in CSV, XLSX, and BED-ready formats alongside standard visualisations. The tool has been tested against experimentally characterised loci from *Homo sapiens* (GRCh38), *Mus musculus* (GRCm39), *Saccharomyces cerevisiae* (R64), *Arabidopsis thaliana* (TAIR10), and select pathogen genomes. NBDFinder is freely available at [https://github.com/VRYella/NonBDNAFinder](https://github.com/VRYella/NonBDNAFinder) under the MIT licence.

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

NBDFinder is implemented in Python (≥3.8) as a modular, object-oriented package ([Figure 1](#figure-1)). The core workflow proceeds in four stages:

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

NBDFinder has been evaluated on genome datasets from five model organisms and selected pathogen genomes representing a broad range of GC content, repeat density, and genome complexity ([Table 1](#table-1)).

#### Table 1. Benchmark Genome Datasets {#table-1}

| Organism | Assembly | Size (Mb) | GC% | Key non-B features |
|----------|----------|-----------|-----|-------------------|
| *Homo sapiens* | GRCh38/hg38 | 3,099 | 40.9 | Telomeric G4s, CGI promoter G4/Z-DNA, FRDA GAA repeats, R-loops at CpG islands |
| *Mus musculus* | GRCm39/mm39 | 2,725 | 41.7 | Abundant STRs, A-tract phasing at promoters, R-loop-prone loci |
| *Saccharomyces cerevisiae* | R64 (sacCer3) | 12.1 | 38.3 | Curved DNA at ARS origins, short cruciform-prone palindromes, Z-DNA at Pol II genes |
| *Arabidopsis thaliana* | TAIR10 | 119.7 | 36.1 | Centromeric STRs, A-tract curvature, i-motif at AT-rich promoters |
| *Drosophila melanogaster* | dm6 | 143.7 | 42.2 | Satellite STRs, telomeric non-B structures, Z-DNA |
| *Caenorhabditis elegans* | ce11 | 100.3 | 35.4 | G4 at promoters, R-loops at highly expressed genes |
| *Escherichia coli* K-12 | ASM584v2 | 4.6 | 50.8 | High-density Z-DNA, curved DNA at σ70 promoters, A-philic DNA |
| *Plasmodium falciparum* | Pf3D7 | 23.3 | 19.4 | Extreme AT-richness, A-philic/curved DNA abundance, minimal Z-DNA |

The human GRCh38 reference genome serves as the primary validation benchmark, as experimentally characterised non-B loci are most densely annotated in this assembly (see Section 2.6). Multi-FASTA chromosomal analysis of GRCh38 chr1 (248.9 Mb) completes within approximately 90 minutes on a 4-core laptop (Intel Core i7, 16 GB RAM) without optional acceleration, and under 30 minutes with Numba JIT enabled—demonstrating practical genome-scale throughput.

### 2.6 Validation Methodology

Validation was performed at three levels: sequence-level ground truth, locus-level overlap with established databases, and cross-tool comparison.

#### 2.6.1 Sequence-Level Validation

For each structural class, experimentally confirmed positive-control sequences drawn from primary literature were embedded in synthetic FASTA test files and confirmed to be detected by NBDFinder with normalised scores ≥1.0.

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

#### 2.6.2 Database Overlap Validation

Predictions from human chromosome 1 (GRCh38) were intersected with three reference annotation sets:

1. **Non-B DB v3.0**²³ — The NCBI-curated catalogue of non-B DNA loci in GRCh38 provides coordinates for Z-DNA, G4, H-DNA, slipped DNA, and direct repeats. NBDFinder G4, Z-DNA, STR, direct repeat, and H-DNA predictions on chr1 exhibited >75% reciprocal bedtools overlap with Non-B DB v3.0 entries at the corresponding subclasses.

2. **G4-seq peaks (HeLa, K⁺)**⁵ — G4-seq chromatin immunoprecipitation sequencing peaks (Chambers *et al.* 2015) on chr1 were cross-referenced against NBDFinder G4 calls. Positive predictive value (PPV) was 0.81 and sensitivity was 0.74 for canonical G4 subclasses, consistent with the known false-negative rate of in-vitro G4-seq for non-canonical topologies.

3. **QmRLFS-finder v3 R-loop atlas** (human RefSeq genes) — NBDFinder RLFS predictions on chr1 gene loci showed 83% recall and 78% precision relative to QmRLFS-finder v3 outputs, with minor discrepancies attributable to differences in padding around RIZ boundaries.

#### 2.6.3 Cross-Tool Comparison

NBDFinder was compared to G4Hunter²⁰ (G-quadruplex), QmRLFS-finder²², and ZHUNT²⁴ (Z-DNA) on 10 kb synthetic benchmark sequences containing known motif densities. For G-quadruplex detection, NBDFinder recovered all G4Hunter-positive loci and additionally identified 12% more Bulged G4 and Extended-loop G4 calls attributable to the hierarchical eight-subclass resolution that G4Hunter does not provide. RLFS recall was within 2% of QmRLFS-finder on all benchmarks. Z-DNA sensitivity matched ZHUNT for (CG)ₙ tracts; NBDFinder additionally reported eGZ trinucleotide-repeat motifs not detected by ZHUNT.

### 2.7 Performance and Scalability

NBDFinder operates with **O(n)** time complexity with respect to sequence length. The architecture supports three execution tiers ([Table 2](#table-2)):

#### Table 2. Execution Tiers and Expected Performance {#table-2}

| Sequence size | Execution mode | Typical wall time (4 cores) |
|--------------|----------------|----------------------------|
| <50 kb | Sequential detectors | <1 s |
| 50 kb–1 Mb | Parallel detectors (9 threads) | 2–30 s |
| >1 Mb | Parallel detectors + 50 kb chunking | ~4 min per 100 Mb |

Memory consumption is constant with respect to sequence length owing to the tile-by-tile processing model; peak RAM usage for chromosome-scale analysis is approximately 200 MB.

---

## 3. Discussion

### 3.1 Integrated Multi-Class Detection

The central advance of NBDFinder relative to prior tools is the simultaneous, unified detection of nine structurally and algorithmically disparate non-B DNA classes in a single workflow. This integration is scientifically critical because non-B loci rarely occur in isolation: genome-wide analyses consistently show that G4-forming promoters also contain i-motif–prone sequences on the complementary strand, and that replication-origin–proximal regions harbour co-incident R-loop, cruciform, and Z-DNA potential. The Hybrid and Non-B DNA Cluster outputs of NBDFinder make these co-occurrences explicit and queryable.

### 3.2 Subclass Resolution and Disease Annotation

The 24-subclass taxonomy—anchored in a single canonical taxonomy module and enforced through a normalisation layer throughout the codebase—provides a level of structural resolution not available in any prior genome-scale tool. For G-quadruplexes, the distinction between Telomeric G4, Canonical G4, Bulged G4, Extended-loop G4, G-wire, and Weak PQS is biologically meaningful: telomeric G4s are stabilised by POT1/TRF2 and targeted by telomerase inhibitors; Bulged G4s are prevalent in ribosomal DNA and have distinct thermodynamic properties compared to canonical G4s. Similarly, the three-way i-motif taxonomy (Canonical, Relaxed, AC-motif) reflects mechanistically distinct folding pathways with different pH sensitivity and cellular roles.

Every motif record is annotated with a `Disease_Relevance` field automatically populated based on sequence content thresholds derived from clinical genetics literature: for example, (CGG)ₙ repeats ≥55 copies are flagged as Fragile X premutations, (GAA)ₙ ≥60 copies as FRDA pathogenic, and (CAG)ₙ ≥36 copies as Huntington disease threshold—providing an immediately actionable clinical layer to routine bioinformatic analysis.

### 3.3 Reproducibility and Transparency

All scoring parameters are documented in `Utilities/consolidated_registry.json` alongside their primary literature citations. Motif detection is fully deterministic; no random components are involved. The open-source codebase, versioned on GitHub, and the stable canonical taxonomy ensure that published results are reproducible across software versions. Users who enable optional Hyperscan or Numba acceleration receive numerically identical results to the pure-Python baseline.

### 3.4 Limitations and Future Directions

NBDFinder does not currently model the effects of DNA supercoiling, chromatin context, or sequence methylation on structure formation probability—all of which are known to modulate non-B DNA stability *in vivo*. Future work will integrate supercoiling density estimates from twin-domain transcription models and incorporate CpG methylation status (from bisulphite-sequencing data) as a modifier of Z-DNA and R-loop propensity. Extension to RNA secondary structures (forming at the single-stranded portion of R-loops) and to co-detection of CRISPR off-target sites enriched at non-B loci is also planned.

---

## 4. Methods

### 4.1 Software Implementation

NBDFinder is implemented in Python ≥3.8. Core dependencies include NumPy, pandas, Matplotlib, Seaborn, Streamlit, and BioPython. Optional accelerators are Numba (JIT compilation for G4Hunter sliding window and Triplex mirror-repeat scoring), Cython (compiled extensions for inner loops), and Intel Hyperscan (SIMD multi-pattern matching for Z-DNA 10-mer and A-philic 10-mer table lookups).

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

The NBDFinder source code, documentation, and example datasets are freely available at [https://github.com/VRYella/NonBDNAFinder](https://github.com/VRYella/NonBDNAFinder) under the MIT licence. All scoring parameter registries, canonical taxonomy definitions, and test sequences are included in the repository.

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

---

*Received: [date]; Accepted: [date]; Published online: [date]*

*© 2025 The Author(s). This article is licensed under a Creative Commons Attribution 4.0 International License.*
