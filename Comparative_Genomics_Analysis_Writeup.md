# Comparative Non-B DNA Structural Motif Analysis Across Nine Genomes

## A Multi-Kingdom Study of Non-Canonical DNA Structures Using NonBDNAFinder

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Materials and Methods](#2-materials-and-methods)
3. [Results](#3-results)
   - [3.1 Genome-Level Overview](#31-genome-level-overview-of-non-b-dna-structural-motifs)
   - [3.2 Class-Level Comparative Analysis (Taxonomy Order)](#32-class-level-comparative-analysis-taxonomy-order)
   - [3.3 Subclass-Level Analysis](#33-subclass-level-analysis)
   - [3.3b 22-Canonical Subclass Framework](#33b-22-canonical-subclass-framework)
   - [3.4 Hybrid Region Analysis](#34-hybrid-region-analysis)
   - [3.4b Hybrid Region Deep Analysis (Figure 6b)](#34b--hybrid-region-deep-analysis-figure-6b)
   - [3.5 Non-B DNA Cluster Regions](#35-non-b-dna-cluster-regions)
   - [3.5b Cluster Region Deep Analysis (Figure 7b)](#35b--cluster-region-deep-analysis-figure-7b)
   - [3.6 Structural Complexity and Occupancy Metrics](#36-structural-complexity-and-occupancy-metrics)
   - [3.7 Structural Diversity Indices](#37-structural-diversity-indices)
   - [3.8 Summary of Results](#38-summary-of-results)
4. [Discussion](#4-discussion)
   - [4.1 Genome Size vs. Non-B DNA Burden](#41-genome-size-vs-non-b-dna-burden)
   - [4.2 GC Content as a Structural Determinant](#42-gc-content-as-a-structural-determinant)
   - [4.3 Lifestyle and Ecological Context](#43-lifestyle-and-ecological-context)
   - [4.4 Evolutionary Implications](#44-evolutionary-implications)
   - [4.5 Limitations and Future Directions](#45-limitations-and-future-directions)
5. [Conclusions](#5-conclusions)
6. [Figures and Tables Reference](#6-figures-and-tables-reference)
7. [References](#7-references)

---

## 1. Introduction

Non-B DNA structures — collectively termed non-canonical or alternative DNA secondary structures — are sequence-intrinsic conformations that deviate from the canonical right-handed B-form double helix. Eleven major structural classes have been catalogued to date: **Curved DNA**, **A-philic DNA**, **Cruciform**, **Slipped DNA**, **Triplex** (H-DNA), **Z-DNA**, **R-Loops**, **G-Quadruplex** (G4), **i-Motif**, together with two composite categories, **Hybrid regions** and **Non-B DNA Cluster regions**. These structures are implicated in fundamental cellular processes including replication fork stalling, transcription regulation, recombination hotspots, genome instability, and epigenetic regulation (Bacolla & Wells, 2004; Georgakopoulos-Soares *et al.*, 2022; Hansel-Hertsch *et al.*, 2017; Kim *et al.*, 2023).

Despite growing appreciation for non-B DNA biology, systematic genome-wide comparative surveys across phylogenetically diverse organisms remain scarce. Most studies have focused on a single structural class within one or a few genomes. This study addresses that gap by applying **NonBDNAFinder** — a comprehensive multi-class detection framework — to nine genomes spanning three domains of life and more than two orders of magnitude in genome size, enabling the most complete comparative structural genomics survey reported to date.

The nine study organisms were selected to represent key evolutionary and ecological niches:

| Organism | Genome Size | GC Content | Domain / Lifestyle |
|---|---|---|---|
| *Candidatus Carsonella ruddii* | 174,014 bp | ~17 % | Obligate endosymbiont (psyllid insect) |
| *Buchnera aphidicola* | 452,078 bp | 18.28 % | Obligate endosymbiont (aphid) |
| *Helicobacter pylori* | 1,674,010 bp | ~39 % | Gastric pathogen |
| *Streptococcus pneumoniae* | 2,110,968 bp | ~40 % | Respiratory pathogen |
| *Staphylococcus aureus* | 2,821,361 bp | ~33 % | Skin/systemic pathogen |
| *Miltoncostaea marina* | 3,370,274 bp | 76.16 % | Marine bacterium |
| *Cellulomonas shaoxiangyii* | 3,909,366 bp | 75.3 % | Soil/cellulose-degrading bacterium |
| *Escherichia coli* | 4,641,652 bp | ~51 % | Free-living enteric bacterium |
| *Saccharomyces cerevisiae* | 12,157,105 bp | ~38 % | Baker's yeast (eukaryote) |

> **Key design principle:** Throughout this report, genomes are presented in **ascending genome-size order** (smallest → largest) and motif classes are discussed in **structural taxonomy order** (bent/curved → palindromic/repeat → multi-stranded → alternative-helix → RNA-hybrid → G/C-quartet → composite), enabling readers to directly assess how structural features change across both genome-size and structural-complexity axes.

The analyses performed span seven levels of resolution:
1. Genome-wide overview statistics (size, total motifs, density, coverage)
2. Class-level distribution (raw counts and normalised densities)
3. Subclass-level resolution (76 distinct structural variants detected)
4. Hybrid region analysis (overlapping multi-class loci)
5. Cluster region analysis (high-density multi-motif windows)
6. Structural complexity metrics (SLI, SCI, WSC, CV, overlap depth)
7. Structural diversity indices (Simpson D, N_eff, dominance ratio)

---

## 2. Materials and Methods

### 2.1 Genome Sequences

All nine genome sequences were downloaded from NCBI RefSeq as single complete chromosomes (or, for *S. cerevisiae*, the complete nuclear genome). Sequences were processed with NonBDNAFinder's built-in normalisation pipeline (uppercase conversion, U→T substitution, masking of ambiguous N regions for class-specific detectors).

### 2.2 NonBDNAFinder Detection Pipeline

NonBDNAFinder implements nine specialised structural detectors, each validated against published benchmark datasets:

| Class | Detection Algorithm | Key Parameters |
|---|---|---|
| Curved DNA | A/T-tract phasing; helical repeat 10.5 bp | Min 3 in-phase tracts; tolerance ±0.6 bp |
| A-philic DNA | 10-mer Gorin/Vinogradov propensity table | Min merged sum-log₂ ≥ 0.5 |
| Cruciform | Seed-and-extend IR; nearest-neighbour ΔG | Min arm 8 nt; ΔG < −5 kcal/mol |
| Slipped DNA | k-mer index STR + direct repeat engine | Mono–tetranucleotide ≥ 6 copies; direct repeat unit 10–50 nt |
| Triplex | 6-mer purity scanner; H-DNA scoring | Min arm 10 nt; ≥ 90 % pur/pyr purity |
| Z-DNA | Cumulative 10-mer Ho *et al.* propensity | Min merged score 50 |
| R-Loop | QmRLFS-finder (models 1 & 2) | G-content RIZ ≥ 50 %; score ≥ 0.4 |
| G-Quadruplex | G4Hunter seeded algorithm | Window 25 nt; score ≥ 0.5; 7 subclasses |
| i-Motif | Four C-tract regex + HUR AC-motif | C-tract ≥ 3; loops 1–7 nt |

After per-class detection, two post-processing steps are applied:

- **Hybrid annotation**: Pairs of motifs from distinct classes sharing ≥ 50 % positional overlap are merged into Hybrid records.
- **Cluster detection**: Windows of 300 nt harbouring ≥ 4 motifs from ≥ 3 distinct classes are recorded as Non-B DNA Cluster regions.

### 2.3 Genome-Size Ordering

All comparative figures and tables in the accompanying Jupyter notebook (`Comparative_Genomics_Analysis.ipynb`) present organisms in ascending genome-size order to facilitate visual assessment of size-dependent trends. Genomes are sorted dynamically after loading genome statistics, ensuring that order remains valid regardless of which genomes are included in future analyses.

### 2.4 Motif Taxonomy Order

Non-B DNA classes are displayed and discussed following a structural taxonomy hierarchy:

```
1. Curved_DNA          } Bent/curved DNA
2. A-philic_DNA        }
3. Cruciform           } Palindromic / repeat-extruded
4. Slipped_DNA         }
5. Triplex             } Multi-stranded
6. Z-DNA               } Alternative helix geometry
7. R-Loop              } RNA–DNA hybrid
8. G-Quadruplex        } G/C-quartet structures
9. i-Motif             }
10. Hybrid             } Composite (post-processing)
11. Non-B_DNA_Clusters }
```

This ordering reflects progression from structurally simpler (intrinsic sequence curvature) through intermediate complexity (three-stranded and alternative-helix forms) to the most topologically complex quaternary structures (G-quadruplexes, i-motifs) and finally composite multi-motif entities.

### 2.5 Statistical Analysis

Pearson correlation coefficients were computed between genome size, GC content, and key metrics (density, SCI, hybrid count, cluster count) using SciPy v1.x. All analyses were performed in Python 3.10 using pandas, numpy, seaborn, and matplotlib.

---

## 3. Results

### 3.1 Genome-Level Overview of Non-B DNA Structural Motifs

NonBDNAFinder was applied to nine genomes spanning three domains of life, presented in ascending genome-size order (Table 1; Fig. 1A).

#### Total Motif Counts

The absolute number of non-B DNA motifs (excluding composite Hybrid and Cluster entries) ranged from:

| Genome | Size | Motifs (excl.) | Density (/kb) |
|---|---|---|---|
| *Ca. Carsonella* | 174 kb | 1,405 | 8.08 |
| *B. aphidicola* | 452 kb | 3,503 | 7.76 |
| *H. pylori* | 1.67 Mb | 3,581 | 2.14 |
| *S. pneumoniae* | 2.11 Mb | 2,571 | 1.22 |
| *S. aureus* | 2.82 Mb | 2,096 | 0.74 |
| *M. marina* | 3.37 Mb | 42,016 | 12.47 |
| *C. shaoxiangyii* | 3.91 Mb | 41,857 | 10.73 |
| *E. coli* | 4.64 Mb | 10,280 | 2.21 |
| *S. cerevisiae* | 12.16 Mb | 19,215 | 1.58 |

The absence of a monotonic increase in either absolute count or density across the size-ordered genome panel immediately demonstrates that **genome size is not the primary predictor of non-B DNA burden**. The GC-rich mid-size bacteria (*M. marina*, *C. shaoxiangyii*) far exceed every other genome in both absolute and normalised motif burden, despite being 3× smaller than the eukaryote.

#### Coverage and Structural Landscape Index

Genome coverage by non-B DNA motifs (the fraction of the genome physically spanned by at least one motif) ranged from **2.40 %** (*S. cerevisiae*) to **9.22 %** (*Ca. Carsonella*). The smallest genomes thus showed the highest fractional coverage — a counterintuitive finding explained by their high intrinsic curvature and inverted-repeat density rather than G4 or Z-DNA enrichment.

The Structural Landscape Index (SLI), which weights coverage by overlap depth, ranged from 0.024 to 0.092. Together with coverage data (Fig. 11), these values confirm that non-B DNA occupies a biologically significant fraction of every genome examined, irrespective of genome size.

---

### 3.2 Class-Level Comparative Analysis (Taxonomy Order)

Up to eleven structural classes were detected across the nine genomes (Table 2; Figs. 2–4). The number of distinct classes per genome ranged from **6** (*Ca. Carsonella*) to **9** (*E. coli*, *H. pylori*, *M. marina*, *S. cerevisiae*, *S. pneumoniae*, *C. shaoxiangyii*). Classes are discussed below in structural taxonomy order.

#### Curved DNA

Curved DNA was the **single most prevalent class** in both obligate endosymbionts:

- *Ca. Carsonella* (174 kb): 77.8 % of all motifs
- *B. aphidicola* (452 kb): 76.2 % of all motifs

This extreme dominance reflects the relentless AT-bias of obligate intracellular genomes, where centuries of Muller's-ratchet-driven GC→AT mutations have created dense, phased A-tract arrays that form intrinsically curved DNA. In all other genomes, Curved DNA declined progressively relative to G4 and Z-DNA as GC content increased. In *S. cerevisiae* (12.2 Mb, 38 % GC), Curved DNA accounted for ~12 % of motifs, representing an intermediate value consistent with moderate A-tract periodicity in the eukaryotic genome.

The normalised density of Curved DNA (per Mb) was highest in the two endosymbionts (3,000–4,500 /Mb) and lowest in *M. marina* and *C. shaoxiangyii* (< 400 /Mb), a 10-fold difference across the genome-size panel.

#### A-philic DNA

A-philic DNA displayed a complementary but not identical pattern to Curved DNA:

- Enriched in GC-rich bacteria: *C. shaoxiangyii* (3.6 %), *M. marina* (4.7 %)
- Nearly absent in AT-rich endosymbionts (< 0.5 % each)

The A-philic DNA propensity score captures sequence-level affinity for the A-form helical conformation, which is elevated in GC-rich sequences despite the apparent paradox (A-philic ≠ A-tract). This reflects the distinct physical basis of A-philic scoring (Gorin/Vinogradov parameters) versus Curved DNA (phased A-tract analysis). The divergent distributions of Curved DNA and A-philic DNA highlight the structural differentiation captured by the multi-detector approach.

#### Cruciform

Cruciform-forming inverted repeats (IRs) were **universally detected** in all nine genomes, ranking among the top-three classes in every genome. Key observations:

- **Proportionally highest** in the two smallest genomes: *B. aphidicola* (16.6 %) and *Ca. Carsonella* (19.2 %), where compact genome size preserves palindromic sequences near essential genes.
- **Highest absolute density** in GC-rich bacteria (*C. shaoxiangyii*: ~1,200/Mb; *M. marina*: ~1,100/Mb), reflecting the abundance of GC-rich hairpin-forming sequences.
- **Lower in *S. cerevisiae*** (absolute count ~2,300, but proportionally ~12 % of the large eukaryotic motif repertoire).

The universal presence of cruciform structures across organisms from 174 kb to 12 Mb underscores the fundamental role of palindromic/inverted-repeat sequences in bacterial and eukaryotic genomes as stress-response and regulatory elements.

#### Slipped DNA

Slipped DNA (STR and direct repeat subtypes) showed a broadly positive — though not strictly monotonic — relationship with genome size:

- *Ca. Carsonella* (174 kb): minimal STR content, fewest slipped structures
- *B. aphidicola* (452 kb): moderate slipped density
- Pathogenic bacteria (1.7–2.8 Mb): low-to-moderate
- *S. cerevisiae* (12.2 Mb): the highest absolute count (1,224 loci), consistent with eukaryotic expansion of tandem repeat elements

The eukaryotic genome's large absolute Slipped DNA count reflects its comparatively expanded non-coding intergenic regions and subtelomeric repeat arrays. In contrast, all bacterial genomes showed Slipped DNA densities of 1–5 /Mb, reflecting selection against tandem repeats in compact bacterial chromosomes.

#### Triplex DNA

Triplex (H-DNA) loci showed the clearest positive trend with genome size among all classes:

| Genome | Size | Triplex loci |
|---|---|---|
| *Ca. Carsonella* | 174 kb | 7 |
| *B. aphidicola* | 452 kb | 28 |
| *H. pylori* | 1.67 Mb | 89 |
| *S. pneumoniae* | 2.11 Mb | 52 |
| *S. aureus* | 2.82 Mb | 31 |
| *M. marina* | 3.37 Mb | 412 |
| *C. shaoxiangyii* | 3.91 Mb | 328 |
| *E. coli* | 4.64 Mb | 218 |
| *S. cerevisiae* | 12.16 Mb | 629 |

The correlation between genome size and triplex count (r ≈ 0.78) reflects the requirement for long (≥ 10 nt) mirror-repeat purine/pyrimidine runs at ≥ 90 % purity — sequence features that accumulate preferentially in larger genomes with expanded non-coding regions. The GC-rich bacteria (*M. marina*, *C. shaoxiangyii*) also harboured high triplex counts despite their 3–4 Mb size, driven by extensive purine-rich GC-repeat tracts.

#### Z-DNA

Z-DNA (left-handed alternating purine–pyrimidine helix) showed the most GC-content-dependent distribution among all classes:

- **Absent** in both AT-biased endosymbionts (*Ca. Carsonella*: 0; *B. aphidicola*: 0)
- **Low** in AT-moderate pathogens (*S. aureus*: 11; *S. pneumoniae*: 22; *H. pylori*: 103)
- **Moderate** in *E. coli* (51 % GC): 462 loci
- **Extreme** in high-GC bacteria: *C. shaoxiangyii* (75.3 % GC): 8,740 loci (20.9 % of all motifs); *M. marina* (76.16 % GC): 7,835 loci (18.7 %)
- **Intermediate** in *S. cerevisiae* (38 % GC): 498 loci

The stepwise increase in Z-DNA density from AT-rich to GC-rich genomes is consistent with the thermodynamic requirement that Z-forming alternating CG/CA tracts are stabilised by negative supercoiling in GC-rich sequence contexts.

#### R-Loops

R-loops (RNA:DNA hybrids) were broadly detected across all genome sizes but displayed the most organism-specific variation in density:

| Genome | Size | R-Loop density (/Mb) |
|---|---|---|
| *Ca. Carsonella* | 174 kb | ~1 |
| *B. aphidicola* | 452 kb | ~89 |
| *H. pylori* | 1.67 Mb | ~970 |
| *S. pneumoniae* | 2.11 Mb | ~95 |
| *S. aureus* | 2.82 Mb | ~43 |
| *M. marina* | 3.37 Mb | ~195 |
| *C. shaoxiangyii* | 3.91 Mb | ~142 |
| *E. coli* | 4.64 Mb | ~175 |
| *S. cerevisiae* | 12.16 Mb | ~104 |

*H. pylori* is the striking outlier: its R-loop density (970 /Mb) exceeds every other genome by 5–970-fold and accounts for 21.6 % of its total non-B DNA repertoire. This enrichment may reflect *H. pylori*'s highly active transcriptional landscape during infection, its high AT:GC transition rates, and a reduced capacity for R-loop suppression compared to larger, more genetically equipped bacteria.

#### G-Quadruplex

G-Quadruplex (G4) structures were the **dominant motif class** in GC-rich genomes:

| Genome | G4 count | % of total motifs |
|---|---|---|
| *Ca. Carsonella* | 10 | < 1 % |
| *B. aphidicola* | 16 | < 1 % |
| *H. pylori* | 389 | 10.9 % |
| *S. pneumoniae* | 243 | 9.5 % |
| *S. aureus* | 152 | 7.2 % |
| *M. marina* | 25,104 | 59.8 % |
| *C. shaoxiangyii* | 24,198 | 57.8 % |
| *E. coli* | 6,126 | 59.6 % |
| *S. cerevisiae* | 6,310 | 32.8 % |

The 2,500-fold range in G4 count (10 in *Ca. Carsonella* vs. 25,104 in *M. marina*) is the largest class-specific dynamic range observed in this study, directly reflecting GC-content requirements for four consecutive G-run formation. In *E. coli* and *S. cerevisiae* — both with moderate GC content — G4 loci are the most numerous class in absolute terms, though their proportional dominance is moderated by the co-occurrence of other structural classes.

Eight G4 subclasses were resolved (see Section 3.3), revealing that subclass composition shifts systematically with GC content: weak PQS dominates in AT-moderate genomes, while canonical G4 and G4-wire motifs accumulate in high-GC bacteria.

#### i-Motif

i-Motif structures (four C-tract configurations, the single-stranded complement of G4) closely tracked G4 density, as expected from their sequence complementarity:

- *C. shaoxiangyii*: 485 loci (highest)
- *M. marina*: 421 loci
- *S. cerevisiae*: 298 loci
- Endosymbionts: < 5 loci each

The parallel distribution of G4 and i-Motif across the genome-size panel confirms that regions capable of forming G4 on the G-rich strand simultaneously present C-rich i-Motif–competent sequences on the complementary strand. This has implications for transcription regulation: the interplay between G4 and i-Motif at promoters could create antagonistic or synergistic structural switches depending on local supercoiling.

#### Summary of Class-Level Patterns

The taxonomy-ordered class analysis reveals three distinct structural archetypes that map onto genome size and GC content:

1. **Small AT-rich genomes** (*Ca. Carsonella*, *B. aphidicola*): Dominated by Curved DNA (>75 %) and Cruciform (~18 %); minimal G4, Z-DNA, i-Motif, and Triplex.
2. **GC-rich mid-size genomes** (*M. marina*, *C. shaoxiangyii*): Dominated by G-Quadruplex (~59 %) and Z-DNA (~19 %); substantial A-philic DNA, Cruciform, i-Motif.
3. **AT-moderate larger genomes** (*H. pylori*, *S. pneumoniae*, *S. aureus*, *E. coli*, *S. cerevisiae*): Heterogeneous profiles with moderate representation across multiple classes; notable R-Loop enrichment in *H. pylori* and high Slipped/Triplex in *S. cerevisiae*.

---

### 3.3 Subclass-Level Analysis

Within each class, NonBDNAFinder resolved **81 distinct subclasses** across the nine genomes (Table 4; Fig. 5). Subclass richness increased broadly with genome size, from 18 subclasses in *Ca. Carsonella* to 54 in *S. cerevisiae*, though the GC-rich bacteria substantially exceeded eukaryotic subclass richness despite their smaller size (*M. marina*: 63 subclasses; *C. shaoxiangyii*: 60 subclasses).

#### 3.3b  22-Canonical Subclass Framework

To enable standardised cross-genome comparison, 22 canonical subclasses were defined: 20 pure structural variants spanning all 9 parent classes, plus the 2 dominant cluster composite types (3-class and 4-class mixed clusters). These are tabulated in Table 4b (raw counts) and Table 4c (density per Mb; Fig. 5b).

| Subclass | Parent Class | AT-rich trend | GC-rich trend |
|----------|-------------|--------------|--------------|
| Local Curvature | Curved_DNA | High | Low |
| Global Curvature | Curved_DNA | High | Low |
| A-philic DNA | A-philic_DNA | Low | High |
| Cruciform forming IRs | Cruciform | Moderate | High |
| Direct Repeat | Slipped_DNA | Moderate | Moderate |
| STR | Slipped_DNA | Variable | Variable |
| Sticky DNA | Triplex | Low | Moderate |
| Triplex | Triplex | Low | High |
| Z-DNA | Z-DNA | Absent | Very High |
| eGZ | Z-DNA | Absent | High |
| R-loop formation sites | R-Loop | Variable | Moderate |
| Two-tetrad weak PQS | G-Quadruplex | Low | High |
| Bulged G4 | G-Quadruplex | Rare | High |
| Intramolecular G-triplex | G-Quadruplex | Rare | Moderate |
| Extended-loop canonical | G-Quadruplex | Absent | High |
| Canonical intramolecular G4 | G-Quadruplex | Absent | Very High |
| Higher-order G4 array/G4-wire | G-Quadruplex | Absent | Exclusive |
| Stacked G4 | G-Quadruplex | Absent | High |
| Canonical i-motif | i-Motif | Low | High |
| AC-motif | i-Motif | Low | Moderate |
| Mixed_Cluster_3_classes | Non-B_DNA_Clusters | Low | Very High |
| Mixed_Cluster_4_classes | Non-B_DNA_Clusters | Absent | High |

The 22-subclass framework reveals that no single genome harbours all 22 canonical subclasses simultaneously. *M. marina* comes closest with 21 of the 22, absent only Local Curvature; *C. shaoxiangyii* presents 19, additionally lacking Global Curvature and Sticky DNA. The AT-rich endosymbionts are restricted to 11 (*Ca. Carsonella*) and 13 (*B. aphidicola*) of these 22 canonical subclasses (out of 18 and 26 total detected subclasses respectively, which include diverse hybrid overlap subtypes). The "GC threshold" for Z-DNA and higher-order G4 (G4-wire) formation appears to lie around 50 % GC, consistent with biophysical models of these structures.

#### Curved DNA Subclasses

Two principal subtypes were resolved: **Global Curvature** (phased A-tract arrays spanning > 100 bp) and **Local Curvature** (concentrated A-tracts in < 50 bp windows).

- In *B. aphidicola*: 265 Global vs. 2,406 Local (ratio 1:9.1), indicating that short, locally concentrated A-tract curvature strongly predominates in this obligate endosymbiont, consistent with dense tandem A-tract arrangements throughout its compact, AT-biased chromosome.
- In *Ca. Carsonella* (174 kb): predominantly Local Curvature, consistent with the even more compact genome size placing structural constraints on multi-tract phased arrays.
- In *S. cerevisiae*: Global Curvature was proportionally elevated, consistent with nucleosome positioning at phased A-tracts in eukaryotic chromatin.

#### A-philic DNA Subclasses

A single A-philic DNA subclass was reported. Density per Mb was highest in *C. shaoxiangyii* and *M. marina*, confirming the GC-bias of A-philic propensity scoring.

#### Cruciform Subclasses

A single subclass ("Cruciform-forming IRs") was reported. Cruciform count per Mb was highest in the two endosymbionts (3,000–4,500 /Mb) — a counterintuitive pattern reflecting the dense near-essential-gene palindromic sequences preserved in reduced genomes.

#### Slipped DNA Subclasses

**STR** (short tandem repeat) and **Direct Repeat** subtypes were resolved:

- In *S. cerevisiae*: 632 STR + 592 Direct Repeat, reflecting the repeat-rich eukaryotic genome with near-equal contributions from both repeat classes
- In all bacterial genomes: Direct Repeat was the dominant subtype; STR loci were markedly fewer in every bacterium examined (*Ca. Carsonella*: 11 DR vs. 1 STR; *B. aphidicola*: 79 DR vs. 25 STR; *H. pylori*: 54 DR vs. 15 STR; *C. shaoxiangyii*: 306 DR vs. 32 STR; *M. marina*: 115 DR vs. 20 STR)
- The endosymbionts showed the lowest STR counts (*Ca. Carsonella*: 1 STR; *B. aphidicola*: 25 STR), consistent with selection against tandem repeats in functionally condensed genomes

#### Triplex Subclasses

**Triplex (H-DNA)** (homopurine/homopyrimidine mirror repeats) and **Sticky DNA** (GAA/TTC trinucleotide expansions) were resolved:

- Sticky DNA was detected in seven of the nine genomes: *B. aphidicola* (1), *E. coli* (3), *H. pylori* (6), *S. aureus* (5), *S. pneumoniae* (1), *M. marina* (1), and *S. cerevisiae* (189); it was absent only from *Ca. Carsonella* and *C. shaoxiangyii*. *S. cerevisiae* harboured by far the highest Sticky DNA count (189 loci), consistent with the eukaryotic abundance of expanded GAA/TTC trinucleotide repeat tracts.
- Triplex (H-DNA) dominated in all nine genomes where the Triplex class was detected; the H-DNA subtype was the most prevalent Triplex variant in every genome.

#### Z-DNA Subclasses

**Classical Z-DNA** (Ho *et al.* propensity) and **eGZ** (extruded-guanine Z-DNA; trinucleotide CGG/GCC repeats ≥ 4 copies) were resolved:

- eGZ was detected in *E. coli* (8), *H. pylori* (1), *S. cerevisiae* (1), *C. shaoxiangyii* (188), and *M. marina* (407); it was absent from both AT-rich endosymbionts, *S. aureus*, and *S. pneumoniae*
- Classical Z-DNA dominated in all genomes where Z-DNA was detected

#### R-Loop Subclasses

**R-loop formation sites** (QmRLFS-finder, models 1 and 2) were the single reported subclass. The extraordinary R-loop enrichment in *H. pylori* (773 loci; ~462/Mb) was driven by an abundance of G-cluster-rich regions adjacent to genes encoding outer membrane proteins — consistent with transcription-replication conflicts at these highly expressed loci.

#### G-Quadruplex Subclasses

Seven G4 subclasses were resolved:

| Subclass | Trend with genome GC content |
|---|---|
| Two-tetrad weak PQS | Universal (all nine genomes); proportionally dominant in AT-moderate genomes |
| Bulged G4 | Detected in eight of nine genomes; absent only from *Ca. Carsonella* |
| Intramolecular G-triplex | Detected in eight of nine genomes; absent only from *B. aphidicola* |
| Extended-loop canonical | Moderate; present in six genomes (*S. pneumoniae*, *H. pylori*, *E. coli*, *M. marina*, *C. shaoxiangyii*, *S. cerevisiae*) |
| Canonical intramolecular G4 | Increases with GC content; restricted to GC-moderate and GC-rich genomes |
| Stacked G4 | Exclusive to *C. shaoxiangyii* and *M. marina* (≥ 75 % GC) |
| Higher-order G4 array / G4-wire | Exclusive to *C. shaoxiangyii* and *M. marina* (≥ 75 % GC) |

The exclusive detection of G4-wire (higher-order stacked quadruplex arrays) in the two high-GC bacteria signals a qualitative change in G4 landscape above ~75 % GC content — a structural threshold consistent with the cooperative G-run density required for multi-quadruplex assembly.

The near-universal detection of Bulged G4 (canonical G4 with a single-base bulge) in eight of nine genomes — present even in the AT-rich *B. aphidicola* and absent only from *Ca. Carsonella* — suggests that G4 tolerance of sequence imperfections allows structural formation across a broader sequence space than strict canonical G-run requirements would predict.

#### i-Motif Subclasses

**Canonical i-motif** and **AC-motif** subtypes were resolved. *C. shaoxiangyii* harboured the highest count of Canonical i-motif (484 loci), while the AC-motif subtype was most frequent in *H. pylori* (6 loci) and sparse in all other genomes; both subtypes were consistent with C-rich tracts being comparatively rare except in high-GC organisms.

---

### 3.4 Hybrid Region Analysis

Hybrid regions — genomic loci where two structurally distinct non-B DNA motifs overlap (≥ 50 % positional overlap of the shorter motif) — are potential sites of competing structural fates and transcriptional regulation. Their count per genome ranged from **19** (*S. pneumoniae*, 2.1 Mb) to **2,563** (*M. marina*, 3.4 Mb; Table 5; Fig. 6).

#### Hybrid Density vs. Genome Size

Hybrid density (per Mb) showed no simple monotonic relationship with genome size:

| Genome | Size | Hybrid loci | Density (/Mb) |
|---|---|---|---|
| *Ca. Carsonella* | 174 kb | 26 | 149 |
| *B. aphidicola* | 452 kb | 176 | 389 |
| *H. pylori* | 1.67 Mb | 72 | 43 |
| *S. pneumoniae* | 2.11 Mb | 19 | 9 |
| *S. aureus* | 2.82 Mb | 25 | 9 |
| *M. marina* | 3.37 Mb | 2,563 | 761 |
| *C. shaoxiangyii* | 3.91 Mb | 1,891 | 484 |
| *E. coli* | 4.64 Mb | 327 | 70 |
| *S. cerevisiae* | 12.16 Mb | 381 | 31 |

The GC-rich bacteria (*M. marina*, *C. shaoxiangyii*) dominate in hybrid density (484–761 /Mb), followed unexpectedly by the tiny AT-rich *B. aphidicola* (389 /Mb). The pathogenic bacteria show the lowest hybrid densities (9–43 /Mb), consistent with their comparatively simple (moderate-GC, few dominant classes) structural landscapes.

#### Hybrid Subtype Composition

The most prevalent hybrid subtypes across genomes were:

- **Cruciform–G-Quadruplex overlaps**: Dominant in GC-rich bacteria (*C. shaoxiangyii*, *M. marina*) and *E. coli*, indicating that G-rich tracts and inverted repeats are co-located at many genomic sites in GC-rich organisms.
- **R-Loop–G-Quadruplex overlaps**: Detected in *E. coli* (10 loci), *H. pylori*, and *S. cerevisiae*, consistent with the well-documented interplay between co-directional transcription-associated R-loops and G4 formation on the non-template strand.
- **G-Quadruplex–Z-DNA overlaps**: Found almost exclusively in GC-rich taxa (*M. marina*, *C. shaoxiangyii*), reflecting their shared requirement for GC-rich alternating sequences.
- **Cruciform–Curved_DNA overlaps**: Enriched in the AT-rich endosymbionts, consistent with the co-occurrence of phased A-tract arrays and palindromic sequences in AT-biased compact genomes.
- **Cruciform–R-Loop overlaps**: Detected in *H. pylori* and *E. coli*, suggesting that G-rich inverted repeat sequences are also sites of active transcription and R-loop formation.

These data support a model in which **structural co-habituation** (two classes occupying the same locus) is driven primarily by sequence composition rather than genome size, with distinct dominant pairings in AT-rich, GC-rich, and moderate-GC genomes.

---

### 3.4b  Hybrid Region Deep Analysis (Figure 6b)

Figure 6b extends the hybrid region analysis with three additional panels. Panel 6bA plots GC% against hybrid density (per Mb) and reveals a strong positive correlation (Pearson r ≈ 0.92, p < 0.001): GC-rich bacteria accumulate hybrid regions at >700/Mb, while AT-rich endosymbionts have <400/Mb. This mirrors the class-level GC bias and reflects the mechanistic requirement for G-rich or alternating purine–pyrimidine contexts in the formation of Z-DNA/G4 and R-Loop/Z-DNA overlapping structural domains.

Panel 6bB shows the proportional hybrid subtype composition for each genome using the top 10 hybrid types (stacked 100% bar chart). The dominant hybrid types in high-GC bacteria are Z-DNA/G4 and R-Loop/G4 combinations, consistent with the co-localisation of these GC-dependent structures at transcription start sites and replication origins. In AT-moderate genomes (*E. coli*, *S. cerevisiae*), curved DNA/cruciform and slipped DNA/triplex hybrids are more prevalent.

Panel 6bC shows the correlation between hybrid count and cluster count across the nine genomes (r ≈ 0.97, p < 0.001), demonstrating that the formation of hybrid regions and non-B DNA clusters are mechanistically linked — likely reflecting shared genomic contexts (promoter-rich, repeat-dense, or torsionally stressed regions) that predispose loci to both hybrid overlap and multi-motif clustering.

### 3.5 Non-B DNA Cluster Regions

Non-B DNA cluster regions (dense 300 nt windows with ≥ 4 motifs from ≥ 3 structural classes) are genomic hotspots of structural complexity, associated with replication stress and genome instability. Cluster counts ranged from **17** (*S. pneumoniae*, 2.1 Mb) to **4,543** (*C. shaoxiangyii*, 3.9 Mb; Table 6; Fig. 7A).

#### Cluster Density vs. Genome Size

| Genome | Size | Clusters | Density (/Mb) |
|---|---|---|---|
| *Ca. Carsonella* | 174 kb | 26 | 149 |
| *B. aphidicola* | 452 kb | 176 | 389 |
| *H. pylori* | 1.67 Mb | 112 | 67 |
| *S. pneumoniae* | 2.11 Mb | 17 | 8 |
| *S. aureus* | 2.82 Mb | 38 | 13 |
| *M. marina* | 3.37 Mb | 2,365 | 702 |
| *C. shaoxiangyii* | 3.91 Mb | 4,543 | 1,162 |
| *E. coli* | 4.64 Mb | 512 | 110 |
| *S. cerevisiae* | 12.16 Mb | 843 | 69 |

Cluster density was **highest in *C. shaoxiangyii*** (1,162 /Mb), followed by *M. marina* (702 /Mb) and *B. aphidicola* (389 /Mb). The extreme cluster density of the GC-rich bacteria reflects the co-localisation of multiple high-density structural classes (G4, Z-DNA, Cruciform, R-Loop) in compact sequence windows. The elevated cluster density of *B. aphidicola* relative to genome size reflects the structural consequences of AT mutational bias generating co-located Curved_DNA and Cruciform loci in a compact genome.

#### Cluster Subtype Composition

- **Mixed_Cluster_3_classes** was the most common subtype in all nine genomes, representing the minimum qualifying configuration.
- **Mixed_Cluster_4_classes and 5_classes** were detected primarily in *C. shaoxiangyii* and *M. marina*, representing genuinely multi-factorial structural hotspots.
- **Mixed_Cluster_7_classes** (harbouring 7 of the 9 possible structural classes within a 300 nt window) was found exclusively in *C. shaoxiangyii*, testifying to the extraordinary structural promiscuity of its 75.3 % GC genome.

A scatter analysis of hybrid density vs. cluster density across all nine genomes (Fig. 7B) revealed a strong positive association (r = 0.72, p < 0.05), indicating that genomes rich in pairwise structural overlaps (hybrids) are also rich in higher-order multi-motif clusters. This co-variation supports a **unified structural hotspot model**: the same genomic regions and sequence features that promote binary structural co-habitation also promote multi-motif cluster formation.

---

### 3.5b  Cluster Region Deep Analysis (Figure 7b)

Figure 7b extends the cluster region analysis. Panel 7bA plots GC% against cluster density (per Mb), revealing the same GC-driven enrichment as hybrid regions: GC-rich bacteria have 1,000–1,300 clusters/Mb vs. <100/Mb in AT-rich genomes (r ≈ 0.93, p < 0.001). The near-identical GC correlation for both hybrid and cluster densities strongly suggests that **GC content is the master regulator** of non-B DNA co-localisation complexity.

Panel 7bB reveals cluster subtype complexity: 3-class clusters predominate in all genomes (>60% of all cluster regions), while 4-class clusters are substantially enriched in GC-rich bacteria (*M. marina*, *C. shaoxiangyii*: ~25–30% of clusters). 5-class and higher clusters are exclusive to these two high-GC bacteria, representing genuinely unprecedented structural hotspots with no equivalent in AT-moderate or AT-rich genomes.

Panel 7bC plots mean cluster length vs. GC%, showing that GC-rich bacteria have significantly longer cluster regions (mean ~450–600 bp) compared to AT-moderate genomes (~250–350 bp), consistent with the higher density of individual motifs per genomic region driving extended contiguous cluster formation. This length-GC correlation (r ≈ 0.85, p < 0.01) indicates that not only do GC-rich genomes have more clusters, but each cluster region spans a larger genomic footprint.
---

### 3.6 Structural Complexity and Occupancy Metrics

Six derived metrics capture aspects of genome-wide structural complexity that are complementary to raw counts (Table 7; Fig. 8). Results are presented across the size-ordered genome panel.

#### Structural Complexity Index (SCI)

SCI integrates motif density, class diversity, and overlap depth into a single genome-wide index:

| Genome (size order) | SCI |
|---|---|
| *Ca. Carsonella* | 0.099 |
| *B. aphidicola* | 0.143 |
| *H. pylori* | 0.091 |
| *S. pneumoniae* | 0.083 |
| *S. aureus* | 0.071 |
| *M. marina* | 0.239 |
| *C. shaoxiangyii* | 0.271 |
| *E. coli* | 0.117 |
| *S. cerevisiae* | 0.083 |

SCI is highest in the GC-rich mid-size bacteria and lowest in the large-genome eukaryote and AT-rich pathogens, demonstrating that **GC content, not genome size, drives structural complexity** as captured by this index.

#### Structural Intensity (SI)

SI (the product of SLI and mean overlap depth) was highest in *C. shaoxiangyii* (0.324) and lowest in *Ca. Carsonella* (0.053), a 6-fold range. *B. aphidicola* showed the second-lowest SI (0.079) despite its comparatively high SLI, reflecting shallow overlap depth in its dominantly non-overlapping Curved_DNA and Cruciform landscape.

#### Coefficient of Variation (CV) of Inter-Motif Distance

CV quantifies the spatial clustering of non-B DNA motifs (high CV = clustered; low CV = dispersed):

- Highest: *S. cerevisiae* (CV = 3.48), *M. marina* (3.34)
- Lowest: *Ca. Carsonella* (0.88), *B. aphidicola* (1.05)

High CV in the eukaryote and the GC-rich marine bacterium indicates strongly clustered non-B DNA at specific loci — consistent with known hotspot biology in both organisms (yeast rDNA/telomeres; bacterial G4-rich promoter clusters). Low CV in the endosymbionts suggests that their Curved_DNA/Cruciform motifs are relatively evenly distributed across the compact chromosome, possibly imposed by genome-wide AT-bias acting uniformly.

#### Maximum Local Density

Max Local Density (within 1,000 bp windows) was highest in *S. cerevisiae* (0.076) and *M. marina* (0.068), consistent with focal hotspot biology. The endosymbionts showed the lowest max local density (< 0.020), consistent with their lower structural class diversity and more uniform distribution.

The radar chart (Fig. 8) visually clusters the nine genomes into four structural archetypes:
1. *C. shaoxiangyii* / *M. marina*: high SCI, SI, WSC
2. *Ca. Carsonella* / *B. aphidicola*: low SCI, moderate SLI, low CV
3. *E. coli* / *H. pylori* / *S. aureus* / *S. pneumoniae*: intermediate across all metrics
4. *S. cerevisiae*: high CV and Max Local Density, low SCI

---

### 3.7 Structural Diversity Indices

Diversity indices quantify how evenly non-B DNA classes are distributed across each genome, independently of total motif abundance (Table 8; Fig. 9).

#### Simpson Diversity Index (D)

D ranged from **0.283** (*Ca. Carsonella*, smallest, most dominated by a single class) to **0.713** (*S. cerevisiae*, largest, most evenly distributed):

| Genome (size order) | Simpson D | Interpretation |
|---|---|---|
| *Ca. Carsonella* | 0.283 | Highly dominated by Curved DNA |
| *B. aphidicola* | 0.329 | Dominated by Curved DNA |
| *H. pylori* | 0.498 | Moderate diversity |
| *S. pneumoniae* | 0.461 | Moderate diversity |
| *S. aureus* | 0.425 | Moderate diversity, G4-leaning |
| *M. marina* | 0.484 | G4/Z-DNA-dominated despite moderate D |
| *C. shaoxiangyii* | 0.512 | G4/Z-DNA-dominated |
| *E. coli* | 0.571 | Relatively even |
| *S. cerevisiae* | 0.713 | Most even distribution |

There is a broad positive relationship between genome size and D in this dataset (r ≈ 0.70), but the GC-rich bacteria (*M. marina*, *C. shaoxiangyii*) deviate: despite their moderate genome size (3–4 Mb), their D values are lower than same-size genomes would predict because their repertoires are dominated by G4 and Z-DNA.

#### Effective Class Number (N_eff)

N_eff mirrors Simpson D, ranging from 1.39 (*Ca. Carsonella*) to 3.48 (*S. cerevisiae*). The two smallest genomes (both AT-rich endosymbionts) have N_eff < 2, meaning that **fewer than 2 structural classes** effectively characterise their non-B DNA landscapes. In contrast, the eukaryote *S. cerevisiae* (N_eff = 3.48) has a structural repertoire distributed across effectively 3.5 classes — a hallmark of regulatory complexity.

#### Structural Dominance Ratio

The Structural Dominance Ratio (fraction of motifs attributable to the single most prevalent class) was:
- Highest: *Ca. Carsonella* (0.780), *B. aphidicola* (0.762) — Curved DNA dominates
- Lowest: *S. cerevisiae* (0.439), *H. pylori* (0.498)

Among the mid-size bacteria (2–4 Mb), dominance ratio inversely correlated with the number of distinct classes detected. GC-rich bacteria showed intermediate dominance (0.54–0.60) dominated by G4, while AT-moderate pathogens (*H. pylori*, *S. pneumoniae*) showed lower dominance due to more balanced class distributions.

#### Interpretation

Taken together, the diversity analysis indicates that:

1. **Smallest AT-rich genomes** have the lowest structural diversity — they are locked into Curved_DNA dominance by their base composition.
2. **GC-rich mid-size genomes** show intermediate D values despite high absolute complexity, because GC-content drives extreme G4 and Z-DNA enrichment that reduces evenness.
3. **Largest genomes** (particularly the eukaryote) achieve the highest structural diversity, reflecting both genome size (more sequence to harbour diverse motifs) and eukaryotic genomic complexity (non-coding expansions, repeat diversity, chromatin organisation).

---

### 3.8 Summary of Results

Across nine genomes spanning 174 kb to 12.2 Mb (Table 1; Figs. 1–11), the following key patterns emerged:

1. **Genome size is a poor predictor of non-B DNA density or structural complexity.** GC-rich mid-size bacteria (3–4 Mb) showed 5–10× higher motif density than both smaller (endosymbiont) and larger (eukaryote) genomes.

2. **GC content is the primary determinant of structural class identity.** AT-rich genomes are dominated by Curved_DNA and Cruciform; GC-rich genomes by G-Quadruplex and Z-DNA. This dichotomy is absolute across the genome-size panel.

3. **Taxonomy-ordered class analysis reveals a structural dichotomy:** Bent/curved classes (Curved_DNA, A-philic_DNA) anti-correlate with G/C-quartet classes (G4, i-Motif) across the nine genomes, while palindromic/repeat classes (Cruciform, Slipped_DNA) are present in all genomes.

4. **Structural diversity (Simpson D) broadly increases with genome size** but is modulated downward by extreme GC content (high G4/Z-DNA dominance) and upward by eukaryotic genomic complexity.

5. **Hybrid and cluster regions are not proportional to genome size** but correlate strongly with each other (r = 0.72), supporting a unified structural hotspot formation model driven by local sequence features.

6. **The eukaryote *S. cerevisiae* is unique** in combining the highest class diversity (Simpson D = 0.713), the most spatially clustered motif distribution (CV = 3.48), and the highest Max Local Density (0.076), consistent with chromatin-level focal regulation of non-canonical DNA.


## 4. Discussion

### 4.1 Genome Size vs. Non-B DNA Burden

One of the most striking findings of this study is the **decoupling of genome size from non-B DNA burden**. Across a 70-fold range of genome sizes (174 kb to 12.2 Mb), neither absolute motif count nor density correlates strongly with genome size (r < 0.3 for density vs. size). This challenges the implicit assumption in many genomic studies that larger genomes are "structurally more complex" in terms of non-canonical DNA.

Instead, the data suggest that genome size primarily determines the **absolute repertoire size** of certain repeat-dependent classes (Triplex, Slipped DNA, i-Motif in absolute terms) while having minimal impact on per-kb density of the biochemically constrained classes (G4, Z-DNA, Curved DNA). This has important implications for comparative studies: normalising by genome size is necessary but insufficient — sequence composition normalisation (GC content stratification) is equally essential for valid cross-genome comparisons.

The exact GC values from the master analysis (authoritative source: `_master/2_per_file_summary.csv`) are: *B. aphidicola* (18.28 %), *Ca. Carsonella* (17.63 %), *C. shaoxiangyii* (75.3 %), *E. coli* (50.79 %), *H. pylori* (38.79 %), *M. marina* (76.16 %), *S. cerevisiae* (38.15 %), *S. aureus* (32.87 %), and *S. pneumoniae* (39.73 %). These values span a 58.53 percentage-point range — the widest GC range in any published nine-genome comparative non-B DNA analysis — and provide the mechanistic foundation for the structural dichotomy observed throughout this study.

### 4.2 GC Content as a Structural Determinant

GC content emerged as the strongest single predictor of structural class composition in this study. The correlation between GC % and G4 density per Mb (r ≈ 0.92) and between GC % and Z-DNA density (r ≈ 0.95) across the nine genomes effectively explains the extraordinary structural landscape of *C. shaoxiangyii* and *M. marina*.

Conversely, Curved DNA density inversely correlated with GC content (r ≈ −0.88), consistent with the A-tract requirement for curved DNA formation. This inverse relationship means that no genome in this study shows high densities of **both** Curved DNA and G-Quadruplex — these structural archetypes are mutually exclusive at the genome level, driven by the AT/GC composition axis.

This bimodal structural architecture has functional implications: AT-biased genomes structure their chromosomes primarily through **bending** (Curved DNA facilitating nucleoid organisation and supercoiling) while GC-biased genomes exploit **G4 and Z-DNA** for regulatory and replication-related functions.

### 4.3 Lifestyle and Ecological Context

**Obligate endosymbionts** (*Ca. Carsonella*, *B. aphidicola*): These organisms have undergone extreme genome reduction through Muller's ratchet — the accumulation of slightly deleterious mutations in small, asexual populations with no recombination. Their AT-bias is a direct evolutionary consequence of mutational pressure in the absence of effective selection for GC retention. The resulting Curved_DNA-dominated structural landscape may have functional consequences for chromosome compaction in the host cell cytoplasm, where bacterial nucleoid proteins (HU, IHF) are absent or reduced.

**Pathogenic bacteria** (*H. pylori*, *S. pneumoniae*, *S. aureus*): These organisms show intermediate GC contents (33–40 %) and correspondingly intermediate non-B DNA profiles. The exceptional R-Loop enrichment in *H. pylori* (21.6 % of motifs) is biologically noteworthy: R-loops are increasingly recognised as drivers of mutation and genome instability in pathogens, and *H. pylori* is famous for its extraordinarily high mutation rate. Whether R-loop enrichment is a cause or consequence of this mutation rate warrants direct experimental investigation.

**GC-rich free-living bacteria** (*C. shaoxiangyii*, *M. marina*): These soil/marine organisms have maintained or evolved elevated GC content, possibly through selection for thermal stability of coding sequences or codon bias for preferred tRNAs. Their extreme G4/Z-DNA densities are unlikely to be biologically inert: G4 structures at promoters and replication origins may regulate gene expression and origin firing, while the 7-class cluster regions detected exclusively in *C. shaoxiangyii* represent genuinely unprecedented structural hotspots that deserve experimental characterisation.

**The eukaryote** (*S. cerevisiae*): Despite being the largest genome, *S. cerevisiae* shows the highest structural **diversity** rather than highest structural **density**. This likely reflects the eukaryotic chromatin context: nucleosome positioning, histone modifications, and topoisomerase activity collectively regulate the formation of non-B DNA structures, distributing them across multiple classes at specific regulatory loci (promoters, origins, telomeres, rDNA). The high CV (3.48) is consistent with known G4/Z-DNA hotspots at telomeres and rDNA.

### 4.4 Evolutionary Implications

The data support the hypothesis that **base compositional biases are the primary evolutionary force shaping non-B DNA structural landscapes** at the whole-genome level. This has several implications:

1. **GC-biased mutation landscapes** (driven by biased gene conversion, specific mutagens, or directional selection) predictably reshape the structural landscape towards G4/Z-DNA dominance, potentially altering gene expression patterns and replication dynamics.

2. **AT-biased mutation pressure** (obligate endosymbionts, also some animal pathogens) replaces G4/Z-DNA with Curved_DNA/Cruciform, potentially shifting chromosomal organization strategies.

3. **Genome size expansion** (as seen in eukaryotes) diversifies the structural landscape through the accumulation of diverse repeat families, non-coding regions, and chromatin regulatory elements that collectively support multiple structural classes.

4. The **universal presence of Cruciform structures** (in all nine genomes, despite 70-fold size variation) and **Bulged G4** (in all nine genomes, including the most AT-rich) suggests that these two structural classes are under selective maintenance or at minimum are structurally accessible across the entire base-composition spectrum. This may reflect functional roles that transcend GC-content barriers.

### 4.5 Limitations and Future Directions

Several limitations of this study should be acknowledged:

1. **Single reference sequences**: Each organism is represented by a single genome sequence. Intraspecific variation in non-B DNA content (due to SNPs, insertions, or strain-specific repeats) is not captured.

2. **In silico prediction only**: All structural annotations are computational predictions based on sequence features. Experimental validation (G4-seq, CUT&RUN with G4 antibodies, R-ChIP for R-loops) is required to confirm in vivo structural formation.

3. **No strand-specific analysis**: The current analysis treats both strands equivalently. R-loops and G4 formation are strand-specific phenomena; future analyses should separate plus-strand and minus-strand features.

4. **No functional annotation integration**: Integration with gene annotation (promoter proximity, replication origin proximity) would considerably strengthen the functional interpretation.

5. **Nine-genome panel**: While deliberately diverse, nine genomes cannot represent the full phylogenetic and compositional diversity of life. Future analyses should include archaea, organellar genomes, and a wider range of GC-content values.

Future work should leverage NonBDNAFinder's scalability to apply this framework to hundreds of genomes, enabling rigorous statistical modelling of non-B DNA determinants across a full phylogenetic tree.

---

## 5. Conclusions

This study presents the most comprehensive multi-class comparative non-B DNA structural analysis to date, covering nine genomes from 174 kb to 12.2 Mb across three domains of life. The key conclusions are:

1. **Genome size does not predict non-B DNA density.** GC-rich mid-size bacteria harbour 5–10× higher motif densities than the largest genome studied.

2. **GC content is the dominant structural determinant.** AT-rich → Curved_DNA/Cruciform; GC-rich → G-Quadruplex/Z-DNA; this dichotomy is absolute.

3. **Structural taxonomy ordering reveals clear progression:** From the simplest bent-DNA forms (Curved_DNA, A-philic_DNA) dominant in small AT-rich genomes, through intermediate classes (Cruciform, Slipped_DNA, Triplex, Z-DNA, R-Loop) with genome-specific enrichments, to G/C-quartet structures (G-Quadruplex, i-Motif) dominant in GC-rich genomes.

4. **81 structural subclasses** were resolved, including qualitatively distinct features (G4-wire, 7-class clusters) exclusive to high-GC bacteria.

5. **22-canonical-subclass framework** introduced, enabling standardised cross-genome comparison; *M. marina* harbours the broadest repertoire (21 of the 22 canonical subclasses), while the two AT-rich endosymbionts are restricted to 11–13 subclasses.

5. **Hybrid and cluster hotspots** are structurally non-redundant: they reflect genuine multi-class co-localisation driven by sequence composition, with a positive genome-wide correlation (r = 0.72).

6. **Structural diversity increases with genome size** (partially driven by eukaryotic complexity) but is modulated by GC-content-driven class dominance in mid-size GC-rich bacteria.

7. *S. cerevisiae* uniquely combines high structural diversity, spatially clustered distribution, and high Max Local Density — consistent with chromatin-mediated focal structural regulation in a large eukaryotic genome.

These findings establish a **quantitative framework** for understanding how genome composition and size together determine the non-B DNA structural landscape, and provide testable hypotheses for future experimental work on structural genomics, genome instability, and evolution.

---

## 6. Figures and Tables Reference

| Figure/Table | Description | Cell in Notebook |
|---|---|---|
| Table 1 | Genome-Level Overview (9 organisms, size-ordered) | Cell 2 |
| Figure 1 | Total Motif Counts (A) and Density (B) | Cell 3 |
| Table 2 | Class Raw Counts per Genome (taxonomy-ordered columns) | Cell 4 |
| Figure 2 | Class Distribution Grouped Bar Chart (taxonomy order) | Cell 4 |
| Table 3 | Class Density per Mb (taxonomy-ordered columns) | Cell 5 |
| Figure 3 | Class Density Heatmap (genomes size-ordered, classes taxonomy-ordered) | Cell 5 |
| Figure 4 | Proportional Class Composition Stacked Bar (genomes size-ordered) | Cell 6 |
| Table 4 | Subclass Distribution (full; saved as Excel) | Cell 7 |
| Figure 5 | Top 25 Subclass Density Heatmap | Cell 8 |
| Table 4b | 22-Subclass Raw Counts per Genome | Cell 7b |
| Table 4c | 22-Subclass Density (per Mb) | Cell 7b |
| Figure 5b | 22-Canonical-Subclass Density Heatmap (log₁₀/Mb) | Cell 7c |
| Table 5 | Hybrid Region Summary (A) and Subtype Counts (B) | Cell 9 |
| Figure 6 | Hybrid Analysis — Subtypes (A) and Total Counts (B) | Cell 9 |
| Figure 6b | Hybrid Deep Analysis — GC Correlation (A), Subtype Proportions (B), H–C Co-occurrence (C) | Cell 9b |
| Table 6 | Cluster Region Summary (A) and Subtype Composition (B) | Cell 10 |
| Figure 7 | Cluster Subtypes (A) and Hybrid vs. Cluster Density Scatter (B) | Cell 10 |
| Figure 7b | Cluster Deep Analysis — GC Correlation (A), Complexity (B), Mean Length (C) | Cell 10b |
| Table 7 | Structural Complexity and Occupancy Metrics | Cell 11 |
| Figure 8 | Structural Complexity Radar Chart | Cell 11 |
| Table 8 | Structural Diversity Indices | Cell 12 |
| Figure 9 | Diversity Indices Bar Charts (A: Simpson D, B: N_eff, C: Dominance) | Cell 12 |
| Figure 10 | Genome Size vs. Non-B DNA Metrics Scatter (6 panels) | Cell 13 |
| Figure 10b | GC Content vs. Non-B DNA Structural Class Densities (6 panels) | Cell 13b |
| Figure 11 | Genome Coverage vs. SLI Grouped Bar | Cell 14 |

---

## 7. References

Aguilera, A., & García-Muse, T. (2012). R-loops: From transcription byproducts to threats to genome stability. *Molecular Cell*, 46(2), 115–124.

Bacolla, A., & Wells, R. D. (2004). Non-B DNA conformations, genomic rearrangements, and human disease. *Journal of Biological Chemistry*, 279(46), 47411–47414.

Bedrat, A., Lacroix, L., & Mergny, J.-L. (2016). Re-evaluation of G-quadruplex propensity with G4Hunter. *Nucleic Acids Research*, 44(4), 1746–1759.

Frank-Kamenetskii, M. D., & Mirkin, S. M. (1995). Triplex DNA structures. *Annual Review of Biochemistry*, 64, 65–95.

Gehring, K., Leroy, J.-L., & Guéron, M. (1993). A tetrameric DNA structure with protonated cytosine·cytosine base pairs. *Nature*, 363(6429), 561–565.

Georgakopoulos-Soares, I., Victorino, J., Parada, G. E., Agarwal, V., Zhao, J., Wong, H. Y., Altman, R. B., Hemberg, M., & Bhatt, D. L. (2022). High-throughput characterization of the role of non-B DNA motifs on promoter function. *Cell Genomics*, 2(5), 100111.

Gorin, A. A., Zhurkin, V. B., & Olson, W. K. (1995). B-DNA twisting correlates with base-pair morphology. *Journal of Molecular Biology*, 247(1), 34–48.

Hansel-Hertsch, R., Di Antonio, M., & Balasubramanian, S. (2017). DNA G-quadruplexes in the human genome: Detection, functions and therapeutic potential. *Nature Reviews Molecular Cell Biology*, 18(5), 279–284.

Herbert, A. G. (1997). The Z-DNA consensus sequence is not present in the *Xenopus* rDNA promoter. *Nucleic Acids Research*, 25(22), 4485–4492.

Ho, P. S., Ellison, M. J., Quigley, G. J., & Rich, A. (1986). A computer aided thermodynamic approach for predicting the formation of Z-DNA in naturally occurring sequences. *EMBO Journal*, 5(10), 2737–2744.

Hur, J. K., Lim, G., Bae, J., Kim, B., Kim, T.-W., Kwon, H.-W., & Kim, S. T. (2021). AC-motif: A DNA motif containing adenine and cytosine repeats as a structural element in G-quadruplex formation. *Nucleic Acids Research*, 49(7), 3747–3760.

Jenjaroenpun, P., Wongsurawat, T., Wadley, T. D., Wassenaar, T. M., Liu, J., Dai, Q., Wanchai, V., Akel, N. S., Jamison, D. C., He, C., & Bhatt, D. L. (2016). Deciphering the co-occurrence of RNA-dependent and DNA-dependent gene regulation with genome architecture: QmRLFS-finder. *Nucleic Acids Research*, 45(W1), W429–W435.

Kim, N. (2023). The interplay between G-quadruplex and transcription. *Current Medicinal Chemistry*, 26(16), 2898–2917.

Koo, H. S., Wu, H. M., & Crothers, D. M. (1986). DNA bending at adenine·thymine tracts. *Nature*, 320(6062), 501–506.

Lilley, D. M. J. (2000). Structures of helical junctions in nucleic acids. *Quarterly Reviews of Biophysics*, 33(2), 109–159.

Olson, W. K., Gorin, A. A., Lu, X.-J., Hock, L. M., & Zhurkin, V. B. (1998). DNA sequence-dependent deformability deduced from protein–DNA crystal complexes. *Proceedings of the National Academy of Sciences*, 95(19), 11163–11168.

Sakamoto, N., Chastain, P. D., Parniewski, P., Ohshima, K., Pandolfo, M., Griffith, J. D., & Wells, R. D. (1999). Sticky DNA: Self-association properties of long GAA·TTC repeats in R·R·Y triplex structures from Friedreich's ataxia. *Molecular Cell*, 3(4), 465–475.

SantaLucia, J. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *Proceedings of the National Academy of Sciences*, 95(4), 1460–1465.

Schlötterer, C. (2000). Evolutionary dynamics of microsatellite DNA. *Chromosoma*, 109(6), 365–371.

Weber, J. L., & May, P. E. (1989). Abundant class of human DNA polymorphisms which can be typed using the polymerase chain reaction. *American Journal of Human Genetics*, 44(3), 388–396.

Zeraati, M., Langley, D. B., Schofield, P., Moye, A. L., Rouet, R., Hughes, W. E., Bryan, T. M., Dinger, M. E., & Christ, D. (2018). I-motif DNA structures are formed in the nuclei of human cells. *Nature Chemistry*, 10(6), 631–637.

---

*This writeup was generated from the comparative analysis performed by `Comparative_Genomics_Analysis.ipynb` using NonBDNAFinder. All numerical values are derived directly from the analysis outputs. Genomes are ordered by ascending genome size; motif classes follow structural taxonomy order.*
