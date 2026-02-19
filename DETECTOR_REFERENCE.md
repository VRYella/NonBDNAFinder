# NonBDNAFinder Detector Reference Guide

**Version:** 2025.1  
**Last Updated:** February 19, 2026  
**Purpose:** Comprehensive technical reference for all Non-B DNA detector classes, subclasses, detection methods, scoring systems, parameters, and outputs.

---

## Table of Contents

1. [A-philic DNA](#1-a-philic-dna)
2. [Curved DNA](#2-curved-dna)
3. [G-Quadruplex](#3-g-quadruplex)
4. [i-Motif](#4-i-motif)
5. [Z-DNA](#5-z-dna)
6. [R-loops](#6-r-loops)
7. [Cruciform](#7-cruciform)
8. [Triplex DNA](#8-triplex-dna)
9. [Slipped DNA](#9-slipped-dna)
10. [Hybrid Motifs](#10-hybrid-motifs)
11. [Non-B DNA Clusters](#11-non-b-dna-clusters)
12. [Cross-Detector Comparison](#cross-detector-comparison)
13. [Performance & Acceleration](#performance--acceleration)

---

## 1. A-philic DNA

### Overview
**Detection Method:** 10-mer propensity table scoring (Vinogradov 2003)  
**Biological Basis:** DNA sequences with high propensity for adopting the A-form helix conformation  
**File:** `Detectors/aphilic/detector.py`

### Subclasses/Variants

| Subclass | Pattern Type | Description |
|----------|-------------|-------------|
| A-philic DNA | 10-mer table | Single class based on sequence propensity scores |

### Detection Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MIN_SUM_LOG2` | 0.5 | Minimum sum of log2 scores for region acceptance |
| **10-mer table** | 208 patterns | Vinogradov 2003 propensity scores |

### Scoring System

| Component | Formula/Method | Range |
|-----------|----------------|-------|
| **Base scoring** | Sum of log2 propensity scores from 10-mer table | Variable |
| **Per-base contribution** | Each base receives score/10 from overlapping 10-mers | Distributed |
| **Region score** | Sum of per-base contributions across merged regions | ≥0.5 |
| **Normalization** | Linear sum; per-base averaging for distribution | [0, ∞) |

**Scoring Formula:**
```
region_score = Σ(per_base_contribution[start:end])
per_base_contribution[i] = Σ(10mer_log2_score / 10) for all 10-mers covering position i
```

### Normalization Method

| Method | Description | Implementation |
|--------|-------------|----------------|
| **Per-base distribution** | Each 10-mer score distributed equally across 10 positions | `_build_per_base_contrib()` |
| **Region merging** | Adjacent/overlapping 10-mers merged into continuous regions | Interval merging |
| **Acceleration** | NumPy vectorization for sequences >1000bp | Optional (2-3x speedup) |

### Output Fields

| Field | Type | Description |
|-------|------|-------------|
| `sum_log2` | float | Total sum of log2 scores for the region |
| `mean_log2_per_10mer` | float | Average log2 score per 10-mer in region |
| `n_10mers` | int | Number of contributing 10-mers |
| `contributing_10mers` | list | List of 10-mers with positions and scores |
| `GC_Content` | float | Percentage of G+C bases (%) |
| `AT_Content` | float | Percentage of A+T bases (%) |
| `Type_Of_Repeat` | string | Classification (High AT, Moderate AT) |
| `disease_relevance` | string | Disease associations (AT-rich expansion disorders) |

### Acceleration Methods

| Technology | Speedup | Availability |
|------------|---------|--------------|
| **Hyperscan** | 2-5x | Optional (pattern matching) |
| **NumPy vectorization** | 2-3x | Optional (large sequences >1KB) |
| **Pure Python fallback** | 1x | Always available |

### References
- Vinogradov AE (2003). *DNA helix: the importance of being GC-rich.* Nucleic Acids Research 31(7):1838-1844
- Bolshoy A et al. (1991). *Curved DNA without A-A.* PNAS 88(6):2312-2316
- Rohs R et al. (2009). *Origins of specificity in protein-DNA recognition.* Annual Review of Biochemistry

---

## 2. Curved DNA

### Overview
**Detection Method:** A-tract and T-tract phasing analysis relative to 10.5bp helical periodicity  
**Biological Basis:** DNA bending induced by A/T homopolymers with phased spacing  
**File:** `Detectors/curved/detector.py`

### Subclasses/Variants

| Subclass | Pattern Type | Detection Method | Reference |
|----------|-------------|------------------|-----------|
| **Global Curvature** | A-phased Repeats (APR) | ≥3 A/T tracts with ~10.5bp spacing | Koo 1986 |
| **Local Curvature** | Long A/T-tracts | Single A≥7 or T≥7 homopolymer tracts | Olson 1998 |

### Detection Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MIN_AT_TRACT` | 3 | Minimum A/T tract length for APR detection |
| `PHASING_CENTER_SPACING` | 11.0 bp | Ideal spacing for A-phased repeats |
| `PHASING_TOL_LOW` | 9.9 bp | Lower tolerance for phasing deviation |
| `PHASING_TOL_HIGH` | 11.1 bp | Upper tolerance for phasing deviation |
| `MIN_APR_TRACTS` | 3 | Minimum number of phased tracts required |
| `LOCAL_LONG_TRACT` | 7 | Minimum length for local curvature tracts |
| `SCORE_THRESHOLD` | 0.1 | Minimum score for motif acceptance |

### Scoring System

| Component | Formula | Range |
|-----------|---------|-------|
| **APR phasing score** | `1.0 - (mean_deviation / max_allowed_dev)` | [0, 1] |
| **Local tract score** | `tract_length / (tract_length + 6)` | [0, 1] |
| **Phasing deviation** | Distance from 10.5bp helical periodicity | Variable |
| **Normalization** | Linear scaling to [0, 1] range | [0, 1] |

**APR Scoring Formula:**
```
mean_spacing = average(spacing_between_tract_centers)
deviation = |mean_spacing - PHASING_CENTER_SPACING|
max_allowed_dev = max(PHASING_CENTER_SPACING - PHASING_TOL_LOW, 
                      PHASING_TOL_HIGH - PHASING_CENTER_SPACING)
phasing_score = 1.0 - (deviation / max_allowed_dev)
```

**Local Scoring Formula:**
```
local_score = tract_length / (tract_length + 6)  # Saturation model
```

### Normalization Method

| Method | Type | Description |
|--------|------|-------------|
| **APR phasing** | Gaussian-like | Deviation from ideal 10.5bp periodicity normalized |
| **Local tract** | Saturation curve | Length-based sigmoid-like normalization |
| **Output range** | [0, 1] | Both scoring methods produce normalized scores |

### Output Fields

| Field | Type | Description |
|-------|------|-------------|
| `A_Tracts` | list | List of A-tract sequences |
| `T_Tracts` | list | List of T-tract sequences |
| `Num_A_Tracts` | int | Count of A-tracts in region |
| `Num_T_Tracts` | int | Count of T-tracts in region |
| `A_Tract_Lengths` | list | Lengths of individual A-tracts |
| `T_Tract_Lengths` | list | Lengths of individual T-tracts |
| `Center_Positions` | list | Center positions of phased tracts |
| `spacing_info` | string | Spacing statistics for APR tracts |
| `GC_Content` | float | Percentage of G+C bases (%) |
| `AT_Content` | float | Percentage of A+T bases (%) |
| `Type_Of_Repeat` | string | "A/T-tract phased repeat (global APR)" or "Long A-tract/T-tract" |
| `Regions_Involved` | string | Summary of tracts and spacing |

### References
- Koo HS et al. (1986). *DNA bending at adenine-thymine tracts.* Nature 320:501-506
- Olson WK et al. (1998). *DNA sequence-dependent deformability.* Philosophical Transactions A
- Crothers DM et al. (1992). *Intrinsically bent DNA.* Journal of Biological Chemistry

---

## 3. G-Quadruplex

### Overview
**Detection Method:** G4Hunter sliding window scoring + pattern-based classification  
**Biological Basis:** Four-stranded structures formed by guanine-rich sequences  
**File:** `Detectors/gquad/detector.py`

### Subclasses/Variants

| Subclass | Pattern | G-tracts | Loop Length | Priority | Reference |
|----------|---------|----------|-------------|----------|-----------|
| **Telomeric G4** | `(TTAGGG)≥4` | 4 | Fixed 3bp (TTA) | 1 | Blackburn 1991 |
| **Stacked canonical G4s** | `((G3+[ACGT]1-7)3G3+)2+` | 8+ | 1-7 bp | 2 | Huppert 2005 |
| **Stacked G4s with linker** | Multiple G4 + 0-20bp linker | 8+ | 1-7 bp | 3 | Bedrat 2016 |
| **Canonical G4** | `G3+[ACGT]1-7` × 4 tracts | 4 | 1-7 bp | 4 | Huppert 2005 |
| **Extended-loop G4** | `G3+[ACGT]1-12` × 4 tracts | 4 | 1-12 bp | 5 | Bedrat 2016 |
| **Higher-order G4/G-wire** | `(G3+[ACGT]1-7)7+` | 7+ | 1-7 bp | 6 | Varizhuk 2019 |
| **G-triplex** | `G3+[ACGT]1-7` × 3 tracts | 3 | 1-7 bp | 7 | Webba da Silva 2007 |
| **Weak PQS** | `G2+[ACGT]1-7` × 4 tracts | 4 | 1-7 bp | 8 | Huppert 2005 |

### Detection Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `WINDOW_SIZE_DEFAULT` | 25 bp | Sliding window size for G4Hunter scoring |
| `MIN_REGION_LEN` | 8 bp | Minimum region length for acceptance |
| `CLASS_PRIORITY` | 8-tier list | Priority order for overlap resolution |

### Scoring System

| Component | Formula/Method | Description |
|-----------|----------------|-------------|
| **G4Hunter score** | G-only density in sliding window | C bases ignored, G bases counted |
| **Sliding window max** | `max_G_count / window_size` | Maximum G-density across all windows |
| **Region scaling** | `normalized_score × (region_length / window_size)` | Length-adjusted scoring |
| **Final score** | G-density × length factor | [0, 1+] range |

**G4Hunter Scoring Formula:**
```
per_position_G_score[i] = 1 if seq[i] == 'G' else 0
max_window_sum = max(Σ(G_score[i:i+window_size]))
normalized_score = max_window_sum / window_size
region_score = normalized_score × (region_length / window_size)
```

### Normalization Method

| Method | Description | Range |
|--------|-------------|-------|
| **Window normalization** | G-count divided by window size | [0, 1] |
| **Length scaling** | Multiply by (region_length / window_size) | [0, ∞) |
| **JIT compilation** | Numba-accelerated sliding window computation | 2-5x speedup |

### Output Fields

| Field | Type | Description |
|-------|------|-------------|
| `Num_Tracts` | int | Number of G-tracts (runs of Gs) |
| `G_Tract_Lengths` | list | Lengths of individual G-tracts |
| `Min_Tract_Length` | int | Shortest G-tract in motif |
| `Max_Tract_Length` | int | Longest G-tract in motif |
| `Avg_Tract_Length` | float | Mean G-tract length |
| `Num_Loops` | int | Number of loops between G-tracts |
| `Loop_Lengths` | list | Lengths of individual loops |
| `Min_Loop_Length` | int | Shortest loop in motif |
| `Max_Loop_Length` | int | Longest loop in motif |
| `Avg_Loop_Length` | float | Mean loop length |
| `GC_Content` | float | Percentage of G+C bases (%) |
| `Priority` | int | Classification priority (1-8) |

### Overlap Resolution

| Method | Description |
|--------|-------------|
| **Priority-based** | CLASS_PRIORITY list determines precedence (telomeric > stacked > canonical > weak) |
| **Greedy selection** | Higher priority variants selected first, lower priority removed if overlapping |

### References
- Huppert JL (2005). *Prevalence of quadruplexes in the human genome.* Nucleic Acids Research
- Bedrat A et al. (2016). *Re-evaluation of G-quadruplex propensity.* Nucleic Acids Research
- Parkinson GN et al. (2002). *Crystal structure of parallel quadruplexes.* Nature
- Mergny JL & Sen D (2019). *DNA quadruple helices in nanotechnology.* Chemical Reviews

---

## 4. i-Motif

### Overview
**Detection Method:** Pattern-based C-tract detection + thermodynamic scoring  
**Biological Basis:** Four-stranded C-rich structures formed at acidic pH  
**File:** `Detectors/imotif/detector.py`

### Subclasses/Variants

| Subclass | Pattern | C-tracts | Loop Length | Reference |
|----------|---------|----------|-------------|-----------|
| **Canonical i-motif** | `C3+[ATGC]1-7` × 4 | 4 | 1-7 bp | Gehring 1993 |
| **HUR AC-motif** (6 variants) | `A3-[linker]-C3-[linker]-C3-[linker]-C3` | 3 C-tracts + 1 A-tract | 4, 5, or 6 bp | Hur 2021 |

### Detection Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MIN_REGION_LEN` | 10 bp | Minimum region length for acceptance |
| `CLASS_PRIORITIES` | canonical: 1, HUR: 2 | Priority order for overlap resolution |
| **Validated sequences** | Reference table | Known i-motif sequences (IM_VAL_001, IM_VAL_002) |

### Scoring System

| Subclass | Formula | Components | Range |
|----------|---------|------------|-------|
| **Canonical i-motif** | `C_density + tract_bonus` | C-density + 0.12×(n_tracts - 2), max 0.4 | [0, 1] |
| **HUR AC-motif** | `base + tract_score + linker_boost` | AC_fraction×0.8 + tract bonus + linker weight | [0, 1] |

**Canonical Scoring Formula:**
```
base_score = C_count / region_length
tract_bonus = 0.12 × min((n_tracts - 2), (0.4 / 0.12))
final_score = min(1.0, base_score + tract_bonus)
```

**HUR AC-motif Scoring Formula:**
```
base = min(0.6, (A_count + C_count) / length × 0.8)
tract_score = min(0.2, 0.12 × (n_C_tracts - 1))
linker_boost = 0.25 (for 4bp/5bp linkers) or 0.12 (for 6bp linkers)
final_score = min(1.0, base + tract_score + linker_boost)
```

### Normalization Method

| Method | Description | Range |
|--------|-------------|-------|
| **Density normalization** | C-count divided by region length | [0, 1] |
| **Tract bonuses** | Fixed increments for structural features | Capped |
| **Clamping** | Final score clamped to [0, 1] | [0, 1] |

### Output Fields

| Field | Type | Description |
|-------|------|-------------|
| `Stems` | list | C-rich stem sequences |
| `Loops` | list | Loop sequences between stems |
| `Num_Stems` | int | Number of C-rich stems (3-4) |
| `Num_Loops` | int | Number of loops (2-3) |
| `Stem_Lengths` | list | Lengths of individual stems |
| `Loop_Lengths` | list | Lengths of individual loops |
| `GC_Content` | float | Total G+C percentage (%) |
| `GC_Content_Stems` | float | G+C percentage in stems only (%) |
| `AC_Content` | float | A+C percentage (for HUR motifs) (%) |
| `Type_Of_Repeat` | string | "Canonical i-motif" or "HUR AC-motif" |

### Overlap Resolution

| Method | Description |
|--------|-------------|
| **Class priority** | Canonical i-motif (priority 1) > HUR AC-motif (priority 2) |
| **Score-based** | Within same class, higher score wins |
| **Greedy selection** | Process by priority, remove overlapping lower-priority motifs |

### References
- Gehring K et al. (1993). *A tetrameric DNA structure with protonated cytosine-cytosine base pairs.* Nature
- Zeraati M et al. (2018). *I-motif DNA structures are formed in the nuclei of human cells.* Nature Chemistry
- Hur JH et al. (2021). *AC-motif: A DNA motif containing adenine and cytosine repeat plays a role in gene regulation.* Nucleic Acids Research

---

## 5. Z-DNA

### Overview
**Detection Method:** 10-mer propensity table (Ho 1986) + extruded-G repeat patterns  
**Biological Basis:** Left-handed helix conformation favored by alternating purines/pyrimidines  
**File:** `Detectors/zdna/detector.py`

### Subclasses/Variants

| Subclass | Pattern Type | Detection Method | Reference |
|----------|-------------|------------------|-----------|
| **Z-DNA** | 10-mer table | Ho 1986 propensity scoring | Ho 1986 |
| **eGZ-DNA** (Extruded-G) | CGG/GGC/CCG/GCC repeats | Trinucleotide repeat detection (≥3 copies) | Herbert 1997 |

### Detection Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MIN_Z_SCORE` | 50.0 | Minimum cumulative Z-score for acceptance |
| `MIN_EGZ_REPEATS` | 3 | Minimum copies of CGG/GGC/CCG/GCC for eGZ |
| `EGZ_BASE_SCORE` | 0.85 | Base score threshold for eGZ motifs |
| `EGZ_MIN_SCORE_THRESHOLD` | 0.80 | Minimum score for eGZ acceptance |
| **10-mer table** | 48 patterns | Ho 1986 Z-DNA propensity scores |

### Scoring System

| Subclass | Formula | Description | Range |
|----------|---------|-------------|-------|
| **Z-DNA** | Sum of 10-mer scores | Per-base contribution summed across region | ≥50.0 |
| **eGZ-DNA** | `threshold × (repeat_count / MIN_EGZ_REPEATS)` | Repeat count normalized scoring | [0.8, ∞) |

**Z-DNA Scoring Formula:**
```
per_base_contribution[i] = Σ(10mer_score / 10) for all 10-mers covering position i
region_score = Σ(per_base_contribution[start:end])
acceptance: region_score ≥ MIN_Z_SCORE (50.0)
```

**eGZ Scoring Formula:**
```
base_threshold = EGZ_BASE_SCORE  # 0.85
repeat_factor = repeat_count / MIN_EGZ_REPEATS
score = base_threshold × repeat_factor
acceptance: score ≥ EGZ_MIN_SCORE_THRESHOLD (0.80)
```

### Normalization Method

| Method | Description | Implementation |
|--------|-------------|----------------|
| **Per-base distribution** | 10-mer scores distributed across 10 positions | Same as A-philic DNA |
| **Region summation** | Per-base contributions summed | Cumulative scoring |
| **Repeat normalization** | eGZ scores normalized by minimum repeat count | Linear scaling |

### Output Fields

| Field | Type | Description |
|-------|------|-------------|
| `Contributing_10mers` | list | Z-DNA: 10-mers with scores and positions |
| `Mean_10mer_Score` | float | Z-DNA: Average 10-mer score |
| `CG_Dinucleotides` | int | Count of CG dinucleotides |
| `AT_Dinucleotides` | int | Count of AT dinucleotides |
| `Alternating_CG_Regions` | int | Count of alternating CG regions |
| `Alternating_AT_Regions` | int | Count of alternating AT regions |
| `Repeat_Unit` | string | eGZ: CGG/GGC/CCG/GCC |
| `Repeat_Count` | int | eGZ: Number of repeat copies |
| `GC_Content` | float | Percentage of G+C bases (%) |

### Acceleration Methods

| Technology | Speedup | Availability |
|------------|---------|--------------|
| **Hyperscan** | 3-5x | Optional (pattern matching) |
| **NumPy vectorization** | 2-3x | Optional (per-base calculation) |
| **Pure Python fallback** | 1x | Always available |

### References
- Ho PS et al. (1986). *Crystal structure of a B-DNA to Z-DNA junction.* EMBO Journal
- Wang AH et al. (2010). *Left-handed Z-DNA: Structure, topology, biological roles.* Journal of Biomolecular Structure & Dynamics
- Herbert A et al. (1997, 2019). *Z-DNA and Z-RNA in human disease.* Communications Biology

---

## 6. R-loops

### Overview
**Detection Method:** QmRLFS algorithm (RIZ + REZ G-rich zone detection)  
**Biological Basis:** Three-stranded RNA-DNA hybrid structures with displaced ssDNA  
**File:** `Detectors/rloop/detector.py`

### Subclasses/Variants

| Component | Pattern | Description | Reference |
|-----------|---------|-------------|-----------|
| **RIZ (RNA Invasion Zone)** | Model 1: G3+ runs | ≥50% G content, G3+ tracts | Jenjaroenpun 2016 |
| **RIZ (RNA Invasion Zone)** | Model 2: G4+ runs | ≥50% G content, G4+ tracts | Jenjaroenpun 2016 |
| **REZ (RNA Exit Zone)** | G-rich window | ≥40% G content, downstream of RIZ | Jenjaroenpun 2016 |

### Detection Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MIN_PERC_G_RIZ` | 50% | Minimum G percentage in RIZ |
| `NUM_LINKER` | 50 bp | Maximum distance between RIZ and REZ |
| `WINDOW_STEP` | 100 bp | Step size for REZ scanning window |
| `MAX_LENGTH_REZ` | 2000 bp | Maximum REZ length |
| `MIN_PERC_G_REZ` | 40% | Minimum G percentage in REZ |
| `QUALITY_THRESHOLD` | 0.4 | Minimum combined score for acceptance |

### Scoring System

| Component | Formula | Range |
|-----------|---------|-------|
| **RIZ scoring** | `RIZ_%G / 100` | [0.5, 1.0] |
| **REZ scoring** | `REZ_%G / 100` | [0.4, 1.0] |
| **Combined score** | `min(1.0, (RIZ_%G + REZ_%G) / 200)` | [0, 1] |
| **Quality filter** | Combined score ≥ 0.4 | Accept/Reject |

**Combined Scoring Formula:**
```
riz_contribution = RIZ_perc_G / 100
rez_contribution = REZ_perc_G / 100
combined_score = min(1.0, riz_contribution + rez_contribution)
acceptance: combined_score ≥ QUALITY_THRESHOLD (0.4)
```

### Normalization Method

| Method | Description | Implementation |
|--------|-------------|----------------|
| **Percentage-based** | G-content normalized to [0, 1] by dividing by 100 | Linear scaling |
| **Seed acceleration** | Prefix-sum G-count array for O(1) window queries | NumPy optimization |
| **Threshold filtering** | RIZ ≥50% G, REZ ≥40% G | Hard thresholds |

### Output Fields

| Field | Type | Description |
|-------|------|-------------|
| `RIZ_Length` | int | Length of RNA invasion zone (bp) |
| `RIZ_Perc_G` | float | G percentage in RIZ (%) |
| `REZ_Length` | int | Length of RNA exit zone (bp) |
| `REZ_Perc_G` | float | G percentage in REZ (%) |
| `Linker_Length` | int | Distance between RIZ and REZ (bp) |
| `GC_Content` | float | Total G+C percentage (%) |
| `GC_Skew` | float | (G-C)/(G+C) ratio |
| `Model` | string | "M1" (G3+) or "M2" (G4+) |

### Acceleration Methods

| Technology | Speedup | Availability |
|------------|---------|--------------|
| **Hyperscan** | 2-4x | Optional (RIZ pattern matching) |
| **NumPy prefix-sum** | 10-20x | Optional (REZ G-window queries) |
| **Pure Python fallback** | 1x | Always available |

### References
- Jenjaroenpun P et al. (2016). *QmRLFS-finder: A model, web server and stand-alone tool for prediction and analysis of R-loop forming sequences.* Nucleic Acids Research
- Aguilera A & García-Muse T (2012). *R loops: From transcription byproducts to threats to genome stability.* Molecular Cell

---

## 7. Cruciform

### Overview
**Detection Method:** Seed-and-extend inverted repeat detection with thermodynamic validation  
**Biological Basis:** Four-way junction structures formed by palindromic sequences  
**File:** `Detectors/cruciform/detector.py`

### Subclasses/Variants

| Subclass | Description | Requirements |
|----------|-------------|--------------|
| **Cruciform** | Thermodynamic palindromes | ΔG ≤ -5.0 kcal/mol, perfect mirror matching |

### Detection Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MIN_ARM` | 8 bp | Minimum arm length for cruciform formation |
| `MAX_ARM` | 50 bp | Maximum arm length to scan |
| `MAX_LOOP` | 12 bp | Maximum loop size at hairpin tip |
| `MAX_MISMATCHES` | 0 | Perfect complementarity required |
| `SEED_SIZE` | 6 bp | Initial seed size for arm matching |
| `DELTA_G_THRESHOLD` | -5.0 kcal/mol | Minimum thermodynamic stability |

### Scoring System

| Component | Formula | Reference |
|-----------|---------|-----------|
| **Nearest-neighbor ΔG** | SantaLucia 1998 parameters | 16 dinucleotide terms |
| **Loop penalty** | `1.75 + 0.6 × ln(loop_length)` | Entropic cost |
| **Total ΔG** | `ΔG_stem + ΔG_loop` | kcal/mol |
| **Normalized score** | `min(1, max(0, -ΔG / 20.0))` | [0, 1] |

**Thermodynamic Scoring Formula:**
```
ΔG_stem = Σ(NN_ENERGY[dinucleotide_i]) for all base pairs in arms
ΔG_loop = 1.75 + 0.6 × ln(loop_length)
ΔG_total = ΔG_stem + ΔG_loop
normalized_score = min(1.0, max(0.0, -ΔG_total / 20.0))
acceptance: ΔG_total ≤ -5.0 kcal/mol
```

**Nearest-Neighbor Energy Parameters (kcal/mol at 37°C):**
| Dinucleotide | ΔG° | Dinucleotide | ΔG° |
|--------------|-----|--------------|-----|
| AA/TT | -1.00 | AG/CT | -1.29 |
| AC/GT | -1.45 | AT/AT | -0.88 |
| CA/TG | -1.45 | CC/GG | -1.84 |
| CG/CG | -2.17 | GA/TC | -1.30 |
| GC/GC | -2.36 | TA/TA | -0.58 |

### Normalization Method

| Method | Description | Range |
|--------|-------------|-------|
| **ΔG normalization** | Free energy divided by -20.0 kcal/mol reference | [0, 1] |
| **Clamping** | Scores clamped to [0, 1] bounds | [0, 1] |
| **Seed-and-extend** | Exact mirror matching algorithm | Binary (match/no match) |

### Output Fields

| Field | Type | Description |
|-------|------|-------------|
| `Left_Arm` | string | Sequence of left arm |
| `Right_Arm` | string | Sequence of right arm (reverse complement of left) |
| `Loop_Seq` | string | Loop sequence at hairpin tip |
| `Arm_Length` | int | Length of stem arms (bp) |
| `Loop_Length` | int | Length of loop (bp) |
| `Mismatches` | int | Number of mismatches in arms (always 0) |
| `Match_Fraction` | float | Fraction of matching base pairs (always 1.0) |
| `DeltaG` | float | Free energy of formation (kcal/mol) |
| `GC_Content` | float | Total G+C percentage (%) |
| `GC_Content_Arms` | float | G+C percentage in arms (%) |
| `GC_Content_Loop` | float | G+C percentage in loop (%) |

### References
- SantaLucia J Jr (1998). *A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics.* PNAS
- Lilley DM (2000). *Structures of helical junctions in nucleic acids.* Quarterly Reviews of Biophysics
- Mizuuchi K et al. (1982). *Cruciform structures in palindromic DNA are favored by DNA supercoiling.* Journal of Molecular Biology

---

## 8. Triplex DNA

### Overview
**Detection Method:** Seed-and-extend mirror repeat detection (H-DNA) + GAA/TTC pattern matching (Sticky DNA)  
**Biological Basis:** Three-stranded structures formed by mirror repeats or disease-associated trinucleotide repeats  
**File:** `Detectors/triplex/detector.py`

### Subclasses/Variants

| Subclass | Pattern | Requirements | Disease Association |
|----------|---------|--------------|---------------------|
| **Triplex (H-DNA)** | Mirror repeats with Pu/Py bias | Arms: 10-100bp, Loop: ≤8bp, Purity ≥90% | General H-DNA |
| **Sticky DNA** | GAA/TTC trinucleotide repeats | ≥4 copies (12bp minimum) | Friedreich Ataxia (FRDA) |

### Detection Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MIN_ARM` | 10 bp | Minimum arm length for H-DNA |
| `MAX_ARM` | 100 bp | Maximum arm length for H-DNA |
| `MAX_LOOP` | 8 bp | Maximum loop size |
| `PURITY_THRESHOLD` | 0.90 | Minimum purine/pyrimidine bias |
| `SEED_SIZE` | 6 bp | Initial seed for mirror matching |
| `MIN_STICKY_COPIES` | 4 | Minimum GAA/TTC copies |
| `SCORE_THRESHOLD` | 0.25 | Minimum score for acceptance |

**Sticky DNA Disease Thresholds:**
| Threshold | Copies | Clinical Significance |
|-----------|--------|----------------------|
| `STICKY_REPLICATION_MIN` | 20 | Replication blockage onset |
| `STICKY_STABLE_MIN` | 40 | Stable expanded allele |
| `STICKY_PATHOGENIC_MIN` | 60 | Pathogenic FRDA threshold |

### Scoring System

| Subclass | Formula | Components | Range |
|----------|---------|------------|-------|
| **H-DNA** | `1.0 + 2.0 × min(1, mechanistic_score)` | 4 weighted factors | [1.0, 3.0] |
| **Sticky DNA** | Piecewise linear by copy number | 4 threshold-based ranges | [1.0, 3.0+] |

**H-DNA Mechanistic Scoring (Frank-Kamenetskii 1995):**
```
L = min(1.0, log(arm_length) / log(35)) × 0.35        # Arm stability
H = exp(-0.4 × loop_length) × 0.20                     # Loop entropy
P = max(0, (purity - 0.8) / 0.2) × 0.30               # Pu/Py bias
I = 1 / (1 + interruptions) × 0.15                     # Structural integrity
raw_score = L + H + P + I
final_score = 1.0 + 2.0 × min(1.0, raw_score)
```

**Sticky DNA Piecewise Scoring (Sakamoto 1999):**
```
Copies 4-19:   score = 1.0 + 0.015 × n
Copies 20-39:  score = 1.3 + 0.03 × (n - 20)    # Replication blockage
Copies 40-59:  score = 2.0 + 0.02 × (n - 40)    # Sticky DNA threshold
Copies 60+:    score = 2.6 + 0.01 × (n - 60)    # Pathogenic FRDA
```

### Normalization Method

| Method | Description | Range |
|--------|-------------|-------|
| **H-DNA log-scaling** | Arm length normalized via log transform | [0, 1] |
| **H-DNA exponential** | Loop penalty exponential decay | [0, 1] |
| **Sticky piecewise** | Linear segments with threshold breakpoints | [1.0, 3.0+] |
| **JIT compilation** | Numba-accelerated scoring functions | 2-3x speedup |

### Output Fields

#### H-DNA (Mirror Triplex)
| Field | Type | Description |
|-------|------|-------------|
| `Arm_Length` | int | Length of mirror repeat arms (bp) |
| `Loop_Length` | int | Length of central loop (bp) |
| `Purity` | float | Purine or pyrimidine bias (fraction) |
| `Interruptions` | int | Number of bases breaking mirror symmetry |
| `GC_Content` | float | Percentage of G+C bases (%) |

#### Sticky DNA (GAA/TTC)
| Field | Type | Description |
|-------|------|-------------|
| `Repeat_Unit` | string | "GAA" or "TTC" |
| `Copy_Number` | int | Number of trinucleotide repeat copies |
| `Replication_Blockage` | boolean | Copy number ≥20 |
| `Sticky_DNA` | boolean | Copy number ≥40 |
| `Pathogenic_FRDA` | boolean | Copy number ≥60 |
| `GC_Content` | float | Percentage of G+C bases (%) |

### Overlap Resolution

| Method | Description |
|--------|-------------|
| **Non-overlapping** | Sticky DNA and H-DNA reported independently |
| **Greedy by score** | Within each subclass, highest score wins overlaps |

### References
- Frank-Kamenetskii MD & Mirkin SM (1995). *Triplex DNA structures.* Annual Review of Biochemistry
- Sakamoto N et al. (1999). *Sticky DNA: Self-association properties of long GAA/TTC repeats in R·R·Y triplex structures.* Journal of Molecular Biology
- Campuzano V et al. (1996). *Friedreich's ataxia: Autosomal recessive disease caused by an intronic GAA triplet repeat expansion.* Science

---

## 9. Slipped DNA

### Overview
**Detection Method:** Unified STR and Direct Repeat detection with mechanistic slippage scoring  
**Biological Basis:** Tandem repeats prone to strand slippage during replication  
**File:** `Detectors/slipped/detector.py`

### Subclasses/Variants

| Subclass | Unit Size (k) | Repeat Types | Clinical Relevance |
|----------|--------------|--------------|-------------------|
| **STR (Short Tandem Repeat)** | 1-6 bp | Mono-, di-, tri-, tetra-, penta-, hexanucleotide | Highly disease-associated |
| **Direct Repeat** | ≥7 bp | Variable unit size repeats | Genomic instability |

### Detection Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MIN_TRACT_LENGTH` | 20 bp | Minimum total tract length |
| `MIN_PURITY` | 0.90 | Minimum repeat purity (fraction of perfect repeat) |
| `MAX_UNIT_SIZE` | 100 bp | Maximum unit size to detect |
| `STR_DIRECT_THRESHOLD` | 7 bp | Unit size boundary (STR <7, Direct Repeat ≥7) |
| `ENTROPY_MIN` | 0.5 bits | Minimum Shannon entropy |
| **CORE mode** | k=1-4, ≥6 copies | Stringent STR detection |
| **LENIENT mode** | k=1-6, ≥4 copies | Relaxed STR detection |
| `MIN_COPIES_DR` | 2 | Minimum copies for direct repeats |

**Disease-Associated Motif Thresholds:**
| Motif | Disease | Pathogenic Threshold |
|-------|---------|---------------------|
| CAG | Huntington Disease | >36 copies |
| CTG | Myotonic Dystrophy | >50 copies |
| CGG | Fragile X Syndrome | >200 copies |
| GAA | Friedreich Ataxia | >66 copies |
| GGGGCC | C9orf72 ALS/FTD | >30 copies |

### Scoring System

| Component | Weight | Formula | Description |
|-----------|--------|---------|-------------|
| **Tract length (L)** | 0.30 | `min(1, log(tract_len) / log(25) / 2)` | Log-scaled length factor |
| **Copy number (C)** | 0.30 | `min(1, log(copy_num) / log(4) / 2)` | Log-scaled repeat count |
| **Unit size (U)** | 0.15 | Unit-dependent weighting | 1.0 (k=2-4), 0.75 (k=1,5-6), 0.5 (k≤20), 0.3 (k>20) |
| **Purity (P)** | 0.15 | `purity²` | Quadratic penalty for imperfect repeats |
| **GC stability (G)** | 0.10 | `0.6 + 0.4 × gc_fraction` | GC-content weighting |
| **Disease bonus** | ×1.15 | Applied to CAG/CTG/CGG/GAA/TTC | Emphasis on clinical relevance |

**Mechanistic Slippage Scoring Formula (Sinden 1994):**
```
L = min(1.0, log(tract_length) / log(25) / 2.0) × 0.30
C = min(1.0, log(copy_number) / log(4) / 2.0) × 0.30
U = unit_size_weight(k) × 0.15
P = purity² × 0.15
G = (0.6 + 0.4 × gc_fraction) × 0.10
raw_score = L + C + U + P + G
base_score = 1.0 + 2.0 × min(1.0, raw_score)
final_score = base_score × disease_motif_bonus(1.15 for CAG/CTG/CGG/GAA/TTC)
```

### Normalization Method

| Method | Description | Range |
|--------|-------------|-------|
| **Log-scaling** | Length and copy number normalized via log transform | [0, 1] |
| **Primitivity** | Longest irreducible repeat unit determined | N/A |
| **Purity calculation** | Alignment-based perfect repeat fraction | [0, 1] |
| **Entropy** | Shannon entropy in bits | [0, 2] |
| **JIT compilation** | Numba-accelerated scoring and entropy | 2-3x speedup |

### Output Fields

| Field | Type | Description |
|-------|------|-------------|
| `Repeat_Unit` | string | Primitive (irreducible) repeat unit |
| `Unit_Size` | int | Length of repeat unit (bp) |
| `Copy_Number` | float | Number of repeat copies (can be fractional) |
| `Purity` | float | Fraction of bases matching perfect repeat |
| `Entropy` | float | Shannon entropy (bits) |
| `GC_Content` | float | Percentage of G+C bases (%) |
| `Slippage_Score` | float | Mechanistic slippage propensity score |
| `Type_Of_Repeat` | string | "mononucleotide", "dinucleotide", "trinucleotide", "tetranucleotide", "pentanucleotide", "hexanucleotide", or "Direct Repeat" |
| `Disease_Relevance` | string | Known disease associations (HD, DM1, FXS, FRDA, C9orf72) |
| **Disease flags** | boolean | Pathogenic threshold indicators per motif type |

### Overlap Resolution

| Method | Description |
|--------|-------------|
| **Max-k dominance** | Longest primitive unit wins overlapping regions |
| **Primitivity check** | Ensures shortest irreducible unit reported |
| **No nested repeats** | Child repeats absorbed into parent primitive unit |

### References
- Sinden RR (1994). *DNA Structure and Function.* Academic Press
- Pearson CE et al. (2005). *Repeat instability: Mechanisms of dynamic mutations.* Nature Reviews Genetics
- Mirkin SM (2007). *Expandable DNA repeats and human disease.* Nature Reviews Genetics

---

## 10. Hybrid Motifs

### Overview
**Detection Method:** Post-processing analysis of overlapping primary motifs  
**Biological Basis:** Regions where multiple Non-B DNA structures co-occur, potentially indicating structural hotspots  
**File:** `Utilities/nonbscanner.py` (method: `_detect_hybrid_motifs`)

### Subclasses/Variants

| Subclass | Pattern | Description |
|----------|---------|-------------|
| **Class1_Class2_Overlap** | Dynamic (any two classes) | Named by constituent classes (e.g., "G-Quadruplex_Z-DNA_Overlap") |

### Detection Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `HYBRID_MIN_OVERLAP` | 0.50 (50%) | Minimum overlap fraction for hybrid detection |
| `HYBRID_MAX_OVERLAP` | 0.99 (99%) | Maximum overlap fraction (excludes complete overlaps) |
| **Minimum motifs** | 2 | At least 2 primary motifs required in sequence |

### Detection Algorithm

```
For each pair of motifs (i, j) where i < j:
    1. Check: different classes (Class[i] ≠ Class[j])
    2. Calculate: overlap_length = min(End[i], End[j]) - max(Start[i], Start[j])
    3. Calculate: overlap_fraction = overlap_length / min(Length[i], Length[j])
    4. Accept if: HYBRID_MIN_OVERLAP < overlap_fraction < HYBRID_MAX_OVERLAP
    5. Create: hybrid motif spanning union of both motifs
```

### Scoring System

| Component | Formula | Description |
|-----------|---------|-------------|
| **Hybrid score** | `(Score[motif1] + Score[motif2]) / 2` | Average of constituent motif scores |
| **No independent scoring** | Derived from primary detectors | Inherits scoring from constituent motifs |

### Normalization Method

| Method | Description |
|--------|-------------|
| **Average scoring** | Simple arithmetic mean of two constituent scores |
| **Range preservation** | Maintains score interpretation from primary detectors |

### Output Fields

| Field | Type | Description |
|-------|------|-------------|
| `Class` | string | Always "Hybrid" |
| `Subclass` | string | Format: "{Class1}_{Class2}_Overlap" |
| `Component_Classes` | list | List of constituent motif classes [Class1, Class2] |
| `Score` | float | Average of constituent motif scores |
| `Start` | int | Start of hybrid region (minimum of constituent starts) |
| `End` | int | End of hybrid region (maximum of constituent ends) |
| `Length` | int | Total length spanning both motifs |

### Biological Significance

Hybrid motifs may represent:
- **Structural competition zones:** Regions where multiple conformations compete
- **Regulatory hotspots:** Areas with enhanced regulatory potential
- **Genomic instability:** Sites prone to multiple structural transitions
- **Cooperative structures:** Adjacent structures that may influence each other

### Overlap Resolution

| Method | Description |
|--------|-------------|
| **De-duplication** | Same hybrid (start, end, class pair) reported only once |
| **Independent from primary** | Does not affect primary motif reporting |
| **All valid overlaps** | All qualifying overlaps reported (not greedy) |

---

## 11. Non-B DNA Clusters

### Overview
**Detection Method:** Sliding window analysis of primary motif density and diversity  
**Biological Basis:** Regions with high concentration of multiple Non-B DNA structure types  
**File:** `Utilities/nonbscanner.py` (method: `_detect_clusters`)

### Subclasses/Variants

| Subclass | Pattern | Description |
|----------|---------|-------------|
| **Mixed_Cluster_N_classes** | Variable | Dynamic naming by number of classes (N = 3, 4, 5, etc.) |

### Detection Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `CLUSTER_WINDOW_SIZE` | 300 bp | Sliding window size for cluster detection |
| `CLUSTER_MIN_MOTIFS` | 4 | Minimum number of motifs in window |
| `CLUSTER_MIN_CLASSES` | 3 | Minimum number of different motif classes |
| **Minimum motifs total** | 3 | At least 3 primary motifs required in sequence |

### Detection Algorithm

```
For each position i in sequence:
    1. Define window: [i, i + CLUSTER_WINDOW_SIZE]
    2. Count motifs: N_motifs in window
    3. Count classes: N_classes (unique motif classes)
    4. Accept if: N_motifs ≥ CLUSTER_MIN_MOTIFS AND N_classes ≥ CLUSTER_MIN_CLASSES
    5. Create: cluster motif spanning [min(motif_starts), max(motif_ends)]
    6. De-duplicate: Skip if window already processed
```

### Scoring System

| Component | Formula | Description |
|-----------|---------|-------------|
| **Cluster score** | `Σ(motif_scores) / N_motifs` | Average score of constituent motifs |
| **No independent scoring** | Derived from primary detectors | Aggregate measure of cluster quality |

### Normalization Method

| Method | Description |
|--------|-------------|
| **Average scoring** | Arithmetic mean of all motifs in cluster window |
| **Complexity indicator** | Higher diversity (more classes) suggests structural complexity |

### Output Fields

| Field | Type | Description |
|-------|------|-------------|
| `Class` | string | Always "Non-B_DNA_Clusters" |
| `Subclass` | string | Format: "Mixed_Cluster_{N}_classes" where N = number of classes |
| `Motif_Count` | int | Total number of motifs in cluster window |
| `Class_Diversity` | int | Number of different motif classes in cluster |
| `Component_Classes` | list | List of motif classes present in cluster |
| `Score` | float | Average score of constituent motifs |
| `Start` | int | Start of cluster region (minimum motif start) |
| `End` | int | End of cluster region (maximum motif end) |
| `Length` | int | Total length of cluster region |

### Biological Significance

Clusters may indicate:
- **Genomic hotspots:** Regions prone to multiple structural transitions
- **Regulatory hubs:** Areas with enhanced regulatory complexity
- **Instability zones:** Sites with elevated recombination/mutation risk
- **Chromatin remodeling sites:** Regions requiring complex chromatin modifications
- **Transcriptional complexity:** Areas with multiple regulatory mechanisms

### Overlap Resolution

| Method | Description |
|--------|-------------|
| **Window-based de-duplication** | Each window start position processed only once |
| **Binary search optimization** | O(n log n) complexity via sorted motif list |
| **No greedy selection** | All qualifying windows reported independently |

---

## Cross-Detector Comparison

### Summary Table

| Class | Subclasses | Score Range | Primary Method | Acceleration | Disease Relevance |
|-------|-----------|-------------|----------------|--------------|-------------------|
| **A-philic DNA** | 1 | [0, ∞) | 10-mer log2 table | Hyperscan, NumPy | AT-rich expansions |
| **Curved DNA** | 2 (APR, Local) | [0, 1] | A/T-tract phasing | None | Regulatory regions |
| **G-Quadruplex** | 8 variants | [0, ∞) | G4Hunter + patterns | Numba JIT | Telomeres, oncogenes |
| **i-Motif** | 2 (Canonical, HUR) | [0, 1] | C-density + tract bonus | None | Gene regulation |
| **Z-DNA** | 2 (Z-DNA, eGZ) | [50, ∞) / [0.8, ∞) | 10-mer table + patterns | Hyperscan, NumPy | Transcription |
| **R-loops** | 2 models (M1, M2) | [0, 1] | QmRLFS (RIZ + REZ) | Hyperscan, NumPy | Genomic instability |
| **Cruciform** | 1 | [0, 1] (ΔG: [-∞, -5]) | Thermodynamic ΔG | None | Recombination hotspots |
| **Triplex** | 2 (H-DNA, Sticky) | [1, 3] | Mechanistic + piecewise | Numba JIT | FRDA (GAA/TTC) |
| **Slipped DNA** | 2 (STR, DR) | [1, 3+] | Mechanistic slippage | Numba JIT | HD, DM1, FXS, FRDA, ALS |
| **Hybrid Motifs** | Dynamic (N² pairs) | Average of pair | Overlap detection | None | Structural hotspots |
| **Non-B DNA Clusters** | Dynamic (by diversity) | Average of cluster | Density + diversity | Binary search | Instability zones |

### Normalization Methods Comparison

| Class | Normalization Type | Method | Output Range |
|-------|-------------------|--------|--------------|
| **A-philic DNA** | Per-base contribution | Linear sum | [0, ∞) |
| **Curved DNA** | Gaussian + saturation | Phasing deviation + length curve | [0, 1] |
| **G-Quadruplex** | Window + length scaling | G-density × region factor | [0, ∞) |
| **i-Motif** | Density + bonuses | C-density with tract/linker bonuses | [0, 1] |
| **Z-DNA** | Per-base contribution | Same as A-philic | [50, ∞) / [0.8, ∞) |
| **R-loops** | Percentage-based | (RIZ_%G + REZ_%G) / 200 | [0, 1] |
| **Cruciform** | Thermodynamic | -ΔG / 20.0 kcal/mol | [0, 1] |
| **Triplex H-DNA** | Mechanistic weighted | 4-factor composite (log, exp, linear) | [1, 3] |
| **Sticky DNA** | Piecewise linear | 4-segment threshold-based | [1, 3+] |
| **Slipped DNA** | Mechanistic weighted | 5-factor composite (2× log-scaled) | [1, 3+] |
| **Hybrid Motifs** | Average scoring | Arithmetic mean of constituent scores | Inherited from primaries |
| **Non-B DNA Clusters** | Average scoring | Arithmetic mean of cluster motifs | Inherited from primaries |

### Overlap Resolution Strategies

| Class | Strategy | Description |
|-------|----------|-------------|
| **A-philic DNA** | Merging | Adjacent/overlapping 10-mers merged into continuous regions |
| **Curved DNA** | Independent | APR and Local reported separately |
| **G-Quadruplex** | Priority-based | 8-tier priority (telomeric > stacked > canonical > weak) |
| **i-Motif** | Priority + score | Canonical (priority 1) > HUR (priority 2), then by score |
| **Z-DNA** | Subclass separation | Z-DNA and eGZ reported independently |
| **R-loops** | Non-overlapping | RIZ-REZ pairs non-overlapping by construction |
| **Cruciform** | Greedy max-score | Highest scoring palindrome wins overlaps |
| **Triplex** | Subclass separation | H-DNA and Sticky DNA independent |
| **Slipped DNA** | Max-k dominance | Longest primitive unit wins (absorbs nested repeats) |
| **Hybrid Motifs** | De-duplication | Same hybrid pair reported once; independent from primaries |
| **Non-B DNA Clusters** | Window-based | Each window start processed once; no greedy selection |

---

## Performance & Acceleration

### Computational Complexity

| Class | Complexity | Bottleneck | Optimization |
|-------|-----------|------------|--------------|
| **A-philic DNA** | O(n) | 10-mer scanning | Hyperscan, NumPy vectorization |
| **Curved DNA** | O(n) | Pattern matching + phasing | None (fast) |
| **G-Quadruplex** | O(n × w) | Sliding window scoring | Numba JIT (2-5x) |
| **i-Motif** | O(n) | Pattern matching | None (fast) |
| **Z-DNA** | O(n) | 10-mer scanning | Hyperscan, NumPy vectorization |
| **R-loops** | O(n) | RIZ pattern + REZ scanning | Hyperscan, NumPy prefix-sum |
| **Cruciform** | O(n²) worst | Seed-and-extend palindrome | Seed pruning (→ O(n)) |
| **Triplex** | O(n²) worst | Seed-and-extend mirror | Seed pruning (→ O(n)), Numba JIT |
| **Slipped DNA** | O(n × k) | Repeat unit extraction | Numba JIT (entropy, scoring) |
| **Hybrid Motifs** | O(m²) | Pairwise overlap detection | Early exit on no overlap |
| **Non-B DNA Clusters** | O(m log m) | Sliding window + sorting | Binary search optimization |

*Note: n = sequence length, w = window size, k = max unit size, m = number of primary motifs*

### Acceleration Technologies

| Technology | Detectors Using | Speedup | Description |
|------------|----------------|---------|-------------|
| **Numba JIT** | G4, Triplex, Slipped | 2-5x | Just-in-time compilation for scoring functions |
| **Hyperscan** | A-philic, Z-DNA, R-loops | 2-5x | High-performance regex engine (optional) |
| **NumPy vectorization** | A-philic, Z-DNA, R-loops | 2-3x | Vector operations for large sequences (>1KB) |
| **Seed-and-extend** | Cruciform, Triplex | 10-100x | Reduces O(n²) to O(n) via early pruning |
| **Prefix-sum arrays** | R-loops (REZ) | 10-20x | O(1) window queries vs O(n) scanning |

### Memory Efficiency

| Feature | Implementation | Benefit |
|---------|----------------|---------|
| **Disk-based storage** | Sequences saved to temp files | Constant ~70MB RAM for any genome size |
| **Streaming FASTA** | Chunk-based parsing | 50-90% memory reduction for large files |
| **Result streaming** | JSONL format on disk | No in-memory accumulation |
| **Chunk analysis** | 50KB chunks with 2KB overlap | Large genome processing (>200MB) |

### Parallel Processing

| Strategy | Implementation | Speedup | Threshold |
|----------|----------------|---------|-----------|
| **Multi-sequence** | ThreadPoolExecutor | 2-8x | ≥2 sequences in FASTA |
| **Chunk parallelization** | ProcessPoolExecutor | 4-12x | Sequences >1MB |
| **Detector parallelization** | ThreadPoolExecutor | 1.5-2x | Sequences >50KB |

---

## Detection Method Details

### Pattern-Based Detection

**Used by:** G-Quadruplex, i-Motif, Z-DNA (eGZ), Triplex (Sticky), Slipped DNA

| Class | Pattern Type | Regex Examples |
|-------|-------------|----------------|
| **G-Quadruplex** | G-tract repeats | `G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}` |
| **i-Motif** | C-tract repeats | `C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}` |
| **Z-DNA (eGZ)** | Trinucleotide | `(CGG){3,}`, `(GGC){3,}`, `(CCG){3,}`, `(GCC){3,}` |
| **Triplex Sticky** | GAA/TTC | `(GAA){4,}`, `(TTC){4,}` |
| **Slipped DNA** | Tandem repeats | Dynamic unit extraction (k=1-100) |

### Table-Based Scoring

**Used by:** A-philic DNA, Z-DNA

| Class | Table Size | Reference | Score Type |
|-------|-----------|-----------|------------|
| **A-philic DNA** | 208 10-mers | Vinogradov 2003 | log2 propensity |
| **Z-DNA** | 48 10-mers | Ho 1986 | Z-DNA propensity |

### Thermodynamic Scoring

**Used by:** Cruciform, R-loops, Triplex (H-DNA)

| Class | Method | Parameters | Energy Type |
|-------|--------|------------|-------------|
| **Cruciform** | SantaLucia ΔG | 16 NN terms + loop penalty | kcal/mol |
| **R-loops** | G-content threshold | RIZ ≥50% G, REZ ≥40% G | Percentage-based |
| **Triplex H-DNA** | Mechanistic model | 4 weighted factors | Composite score |

### Mechanistic Scoring

**Used by:** Curved DNA (phasing), Triplex (H-DNA, Sticky), Slipped DNA

| Class | Model Basis | Key Factors | Reference |
|-------|------------|-------------|-----------|
| **Curved DNA** | Helical phasing | A/T-tract spacing deviation from 10.5bp | Koo 1986 |
| **Triplex H-DNA** | H-DNA stability | Arm length, loop size, purity, interruptions | Frank-Kamenetskii 1995 |
| **Triplex Sticky** | Disease thresholds | Piecewise linear by copy number (4 ranges) | Sakamoto 1999 |
| **Slipped DNA** | Slippage propensity | Length, copies, unit size, purity, GC content | Sinden 1994 |

---

## Output Schema Comparison

### Core Fields (All Detectors)

| Field | Type | Description |
|-------|------|-------------|
| `Sequence_Name` | string | Name of the input sequence |
| `Class` | string | Motif class (e.g., "G-Quadruplex") |
| `Subclass` | string | Motif subclass (e.g., "Canonical G4") |
| `Start` | int | Start position (1-indexed) |
| `End` | int | End position (inclusive) |
| `Length` | int | Motif length in base pairs |
| `Strand` | string | "+" or "-" (or "both") |
| `Score` | float | Detector-specific score |
| `Method` | string | Detection method identifier |
| `Pattern_ID` | string | Unique pattern identifier |

### Detector-Specific Fields

| Detector | Unique Fields | Description |
|----------|--------------|-------------|
| **A-philic** | `sum_log2`, `mean_log2_per_10mer`, `n_10mers`, `contributing_10mers` | 10-mer scoring details |
| **Curved** | `A_Tracts`, `T_Tracts`, `Center_Positions`, `spacing_info` | Tract and phasing information |
| **G-Quadruplex** | `Num_Tracts`, `Loop_Lengths`, `G_Tract_Lengths`, `Priority` | G4 structural features |
| **i-Motif** | `Stems`, `Loops`, `Num_Stems`, `GC_Content_Stems`, `AC_Content` | C-rich structural details |
| **Z-DNA** | `Contributing_10mers`, `CG_Dinucleotides`, `Alternating_CG_Regions`, `Repeat_Unit` (eGZ) | Z-forming sequence features |
| **R-loops** | `RIZ_Length`, `RIZ_Perc_G`, `REZ_Length`, `REZ_Perc_G`, `Linker_Length`, `GC_Skew`, `Model` | Three-zone structure |
| **Cruciform** | `Left_Arm`, `Right_Arm`, `Loop_Seq`, `DeltaG`, `Match_Fraction` | Palindrome thermodynamics |
| **Triplex H-DNA** | `Arm_Length`, `Loop_Length`, `Purity`, `Interruptions` | Mirror repeat structure |
| **Triplex Sticky** | `Repeat_Unit`, `Copy_Number`, `Replication_Blockage`, `Sticky_DNA`, `Pathogenic_FRDA` | Disease thresholds |
| **Slipped DNA** | `Repeat_Unit`, `Unit_Size`, `Copy_Number`, `Purity`, `Entropy`, `Slippage_Score` | Tandem repeat features |

### Confidence/Quality Tiers

| Detector | Tiers | Threshold Criteria |
|----------|-------|-------------------|
| **A-philic DNA** | 2 tiers | High AT (>80%), Moderate AT (≤80%) |
| **Curved DNA** | 2 types | Global (APR ≥3 tracts), Local (single tract ≥7bp) |
| **G-Quadruplex** | 8 priorities | Telomeric (1) → Weak PQS (8) |
| **i-Motif** | 2 classes | Canonical (priority 1), HUR AC-motif (priority 2) |
| **Z-DNA** | Score-based | Z-DNA ≥50.0, eGZ ≥0.80 |
| **R-loops** | Quality filter | Combined score ≥0.4 |
| **Cruciform** | ΔG threshold | ΔG ≤ -5.0 kcal/mol |
| **Triplex H-DNA** | Score-based | Score ≥0.25 (1-3 scale) |
| **Triplex Sticky** | 4 clinical ranges | Weak (4-19), Replication (20-39), Sticky (40-59), Pathogenic (60+) |
| **Slipped DNA** | Disease thresholds | Motif-specific pathogenic copy numbers |

---

## Parameter Summary Tables

### Length Constraints

| Detector | Min Length | Max Length | Typical Length | Notes |
|----------|-----------|------------|----------------|-------|
| **A-philic** | 10 bp (1 10-mer) | No limit | 50-200 bp | Merged regions |
| **Curved APR** | ~25 bp (3×3bp + spacing) | No limit | 30-60 bp | Phased tracts |
| **Curved Local** | 7 bp | No limit | 7-15 bp | Single tract |
| **G-Quadruplex** | 15 bp (G3N1G3N1G3N1G3) | ~100 bp | 20-40 bp | Loop-dependent |
| **i-Motif** | 10 bp | ~50 bp | 15-30 bp | 4 C-tracts + loops |
| **Z-DNA** | 10 bp (1 10-mer) | No limit | 20-50 bp | Alternating Pu/Py |
| **Z-DNA eGZ** | 9 bp (3 copies) | No limit | 9-30 bp | CGG/GGC repeats |
| **R-loops** | ~20 bp (RIZ+REZ min) | ~2050 bp (REZ max 2000) | 100-500 bp | RIZ+linker+REZ |
| **Cruciform** | 16 bp (2×8bp arms) | 112 bp (2×50bp + 12bp loop) | 20-40 bp | Palindrome |
| **Triplex H-DNA** | 20 bp (2×10bp arms) | 208 bp (2×100bp + 8bp loop) | 25-60 bp | Mirror repeat |
| **Triplex Sticky** | 12 bp (4 copies) | No limit | 12-200 bp | GAA/TTC |
| **Slipped STR** | 20 bp | Variable | 20-100 bp | k=1-6 |
| **Slipped DR** | 20 bp | Variable | 20-200 bp | k≥7 |

### Purity/Quality Thresholds

| Detector | Purity Requirement | Quality Metric | Threshold |
|----------|-------------------|----------------|-----------|
| **A-philic** | N/A | sum_log2 | ≥0.5 |
| **Curved** | N/A | Phasing/tract score | ≥0.1 |
| **G-Quadruplex** | N/A | G4Hunter score | ≥0.5 typical |
| **i-Motif** | N/A | C-density + bonuses | Variable |
| **Z-DNA** | N/A | Cumulative score | ≥50.0 |
| **R-loops** | RIZ ≥50% G, REZ ≥40% G | Combined score | ≥0.4 |
| **Cruciform** | Perfect match (0 mismatches) | ΔG | ≤-5.0 kcal/mol |
| **Triplex H-DNA** | ≥90% Pu or Py | Purity | ≥0.90 |
| **Triplex Sticky** | Implicit (perfect GAA/TTC) | Copy number | ≥4 |
| **Slipped DNA** | ≥90% | Purity | ≥0.90 |

### GC Content Ranges

| Detector | Typical GC% | GC Preference | Notes |
|----------|-------------|--------------|-------|
| **A-philic** | 20-40% | AT-rich | Low GC by definition |
| **Curved** | 20-45% | AT-rich | A/T-tract based |
| **G-Quadruplex** | 60-100% | Very GC-rich | G-rich by definition |
| **i-Motif** | 60-100% | Very GC-rich | C-rich by definition |
| **Z-DNA** | 70-100% | GC-rich | Alternating CG favored |
| **R-loops** | 60-90% | GC-rich | G-rich RIZ/REZ |
| **Cruciform** | Variable | Any | Palindrome-dependent |
| **Triplex** | Variable | Any | Mirror/repeat-dependent |
| **Slipped DNA** | Variable | Any | Repeat-unit dependent |

---

## Scoring System Summary

### Scoring Philosophies

| Approach | Detectors | Description |
|----------|-----------|-------------|
| **Table-based propensity** | A-philic, Z-DNA | Pre-computed sequence-specific scores from literature |
| **Density scoring** | G4, i-Motif, R-loops | Base composition within windows/regions |
| **Structural geometry** | Curved (phasing) | Geometric relationships (spacing, periodicity) |
| **Thermodynamic** | Cruciform (ΔG) | Free energy calculations |
| **Mechanistic modeling** | Triplex, Slipped | Multi-factor composite scoring based on biological mechanisms |
| **Pattern matching** | All primary detectors | Initial seed/candidate identification |
| **Derived metrics** | Hybrid, Clusters | Post-processing aggregate metrics |

### Score Interpretation

| Range | Interpretation | Detectors |
|-------|---------------|-----------|
| **[0, 1]** | Normalized probability/confidence | Curved, i-Motif, R-loops, Cruciform (normalized) |
| **[0, ∞)** | Unbounded cumulative score | A-philic, G4, Z-DNA |
| **[1, 3]** | Mechanistic risk scale | Triplex (H-DNA, Sticky), Slipped DNA |
| **Negative** | Free energy (kcal/mol) | Cruciform (ΔG, more negative = more stable) |
| **Inherited** | Derived from constituents | Hybrid Motifs, Clusters (average of components) |

---

## Validation & Quality Assurance

### Reference Data Sources

| Detector | Reference Database | Validation Status |
|----------|-------------------|-------------------|
| **A-philic** | Vinogradov 2003 table | ✓ Published propensities |
| **Curved** | Koo 1986, Olson 1998 | ✓ Experimental validation |
| **G-Quadruplex** | Huppert 2005, Bedrat 2016 | ✓ G4-seq experimental |
| **i-Motif** | Zeraati 2018, Hur 2021 | ✓ In vivo detection |
| **Z-DNA** | Ho 1986, Herbert 2019 | ✓ Crystal structures |
| **R-loops** | DRIP/DRIPc-seq data | ✓ Genome-wide maps |
| **Cruciform** | SantaLucia 1998 ΔG | ✓ Thermodynamic tables |
| **Triplex** | Frank-Kamenetskii 1995 | ✓ Biochemical studies |
| **Slipped** | Sinden 1994, Pearson 2005 | ✓ Clinical repeat databases |

### Testing Coverage

| Detector | Unit Tests | Integration Tests | Validation Data |
|----------|-----------|-------------------|-----------------|
| **All detectors** | ✓ Basic functionality | ✓ End-to-end workflows | ✓ Known motif examples |
| **Slipped DNA** | ✓ Disease motifs | ✓ Pathogenic thresholds | ✓ Clinical databases |
| **G-Quadruplex** | ✓ 8 variant types | ✓ Priority resolution | ✓ Telomeric sequences |
| **Cruciform** | ✓ Thermodynamics | ✓ ΔG calculations | ✓ Palindrome library |

---

## Clinical & Disease Relevance

### Disease-Associated Motifs

| Motif | Detector | Disease | Pathogenic Threshold | Gene Examples |
|-------|----------|---------|---------------------|---------------|
| **CAG** | Slipped DNA | Huntington Disease | >36 copies | HTT |
| **CTG** | Slipped DNA | Myotonic Dystrophy | >50 copies | DMPK |
| **CGG** | Slipped DNA, Z-DNA eGZ | Fragile X Syndrome | >200 copies | FMR1 |
| **GAA** | Triplex Sticky, Slipped | Friedreich Ataxia | >66 copies | FXN |
| **TTC** | Triplex Sticky, Slipped | Friedreich Ataxia | >66 copies | FXN (complement) |
| **GGGGCC** | Slipped DNA | C9orf72 ALS/FTD | >30 copies | C9orf72 |
| **A-rich** | A-philic DNA | AT expansion disorders | Context-dependent | Various |

### Clinical Scoring Features

| Feature | Implementation | Benefit |
|---------|----------------|---------|
| **Disease motif bonus** | 1.15× score multiplier | Emphasizes clinical relevance |
| **Pathogenic flags** | Boolean indicators | Direct clinical interpretation |
| **Copy number thresholds** | Literature-based cutoffs | Aligns with diagnostic criteria |
| **Repeat unit annotation** | Explicit disease motif identification | Facilitates genetic counseling |

---

## Software Architecture

### Base Detector Interface

All detectors inherit from `BaseMotifDetector` with the following interface:

| Method | Purpose | Return Type |
|--------|---------|-------------|
| `get_motif_class_name()` | Return class name | string |
| `get_patterns()` | Return detection patterns | Dict[str, List[Tuple]] |
| `detect_motifs()` | Main detection method | List[Dict[str, Any]] |
| `calculate_score()` | Compute sequence score | float |
| `annotate_sequence()` | Return structural annotations | List[Dict[str, Any]] |

### Detection Pipeline

```
Input Sequence
    ↓
Pattern Matching / Table Lookup
    ↓
Candidate Region Identification
    ↓
Scoring & Annotation
    ↓
Quality Filtering (thresholds)
    ↓
Overlap Resolution
    ↓
Output Formatting
    ↓
Motif List (JSON/CSV)
```

---

## Usage Guidelines

### Choosing Detection Parameters

| Goal | Recommendation | Detectors Affected |
|------|---------------|-------------------|
| **High sensitivity** | Lower score thresholds | All (adjust MIN/THRESHOLD params) |
| **High specificity** | Raise score thresholds | All (stricter filtering) |
| **Disease screening** | Use pathogenic thresholds | Slipped, Triplex Sticky |
| **Regulatory elements** | Focus on G4, i-Motif, Curved | Promoter/enhancer analysis |
| **Genome stability** | All detectors | Comprehensive profiling |

### Computational Considerations

| Sequence Size | Recommended Approach | Expected Time |
|--------------|---------------------|---------------|
| **< 50 KB** | Sequential detection | Seconds |
| **50 KB - 1 MB** | Parallel detectors | Seconds to minutes |
| **1 MB - 100 MB** | Chunking + parallel | Minutes |
| **> 100 MB** | Disk-based + chunking | 5-10 minutes |

---

## Quality Standards

### Publication Requirements

All detectors meet the following quality standards:

- ✅ **Peer-reviewed foundations**: Every algorithm based on published research
- ✅ **Thermodynamic grounding**: Scoring systems aligned with biophysical principles
- ✅ **Reproducibility**: All parameters documented and accessible
- ✅ **Validation**: Tested against experimental data and reference databases
- ✅ **Output quality**: Meets Nature/Science/Cell journal standards

### Code Quality Features

- ✅ **Type hints**: Full type annotation for all public methods
- ✅ **Documentation**: Comprehensive docstrings with references
- ✅ **Testing**: Unit tests and integration tests
- ✅ **Performance**: JIT compilation and vectorization where applicable
- ✅ **Maintainability**: Modular architecture with clear separation of concerns

---

## Quick Reference

### One-Line Summaries

| Detector | One-Line Description | Key Output |
|----------|---------------------|-----------|
| **A-philic** | 10-mer propensity table (Vinogradov 2003), log2 scoring | sum_log2, n_10mers |
| **Curved** | A/T-tract phasing (10.5bp periodicity) + long tracts (≥7bp) | Phasing score, tract counts |
| **G-Quadruplex** | G4Hunter sliding window + 8 pattern variants | G-tract structure, priority |
| **i-Motif** | C-density + tract bonus; canonical & HUR AC-motifs | Stem/loop structure |
| **Z-DNA** | 10-mer table (Ho 1986) + eGZ trinucleotide repeats | 10-mer scores, repeat unit |
| **R-loops** | QmRLFS: RIZ (G-rich invasion) + REZ (G-rich exit) | RIZ/REZ lengths, %G |
| **Cruciform** | Seed-and-extend palindromes, SantaLucia ΔG ≤ -5.0 kcal/mol | ΔG, arm/loop structure |
| **Triplex H-DNA** | Mirror repeats, mechanistic 4-factor scoring | Purity, interruptions |
| **Triplex Sticky** | GAA/TTC repeats, piecewise clinical thresholds | Copy number, FRDA flags |
| **Slipped DNA** | STR (k=1-6) + DR (k≥7), mechanistic slippage scoring | Repeat unit, copy number, disease flags |
| **Hybrid Motifs** | Overlapping primary motifs (50-99% overlap); averaged scoring | Component classes, overlap fraction |
| **Non-B DNA Clusters** | High-density regions (≥4 motifs, ≥3 classes in 300bp window) | Motif count, class diversity |

---

## Appendix: Complete Parameter List

### A-philic DNA Parameters
```python
MIN_SUM_LOG2 = 0.5  # Minimum sum of log2 scores
# 10-mer table: 208 patterns from Vinogradov 2003
```

### Curved DNA Parameters
```python
MIN_AT_TRACT = 3
PHASING_CENTER_SPACING = 11.0  # bp
PHASING_TOL_LOW = 9.9          # bp
PHASING_TOL_HIGH = 11.1        # bp
MIN_APR_TRACTS = 3
LOCAL_LONG_TRACT = 7           # bp
SCORE_THRESHOLD = 0.1
```

### G-Quadruplex Parameters
```python
WINDOW_SIZE_DEFAULT = 25  # bp
MIN_REGION_LEN = 8        # bp
CLASS_PRIORITY = [
    "telomeric_g4",           # Priority 1
    "stacked_canonical_g4s",  # Priority 2
    "stacked_g4s_linker",     # Priority 3
    "canonical_g4",           # Priority 4
    "extended_loop_g4",       # Priority 5
    "higher_order_g4",        # Priority 6
    "g_triplex",              # Priority 7
    "weak_pqs"                # Priority 8
]
```

### i-Motif Parameters
```python
MIN_REGION_LEN = 10  # bp
CLASS_PRIORITIES = {
    'canonical_imotif': 1,
    'hur_ac_motif': 2
}
# HUR linker lengths: 4, 5, 6 bp
```

### Z-DNA Parameters
```python
MIN_Z_SCORE = 50.0             # Cumulative 10-mer score threshold
MIN_EGZ_REPEATS = 3            # Minimum CGG/GGC/CCG/GCC copies
EGZ_BASE_SCORE = 0.85
EGZ_MIN_SCORE_THRESHOLD = 0.80
# 10-mer table: 48 patterns from Ho 1986
```

### R-loops Parameters
```python
MIN_PERC_G_RIZ = 50      # %
NUM_LINKER = 50          # bp
WINDOW_STEP = 100        # bp
MAX_LENGTH_REZ = 2000    # bp
MIN_PERC_G_REZ = 40      # %
QUALITY_THRESHOLD = 0.4
```

### Cruciform Parameters
```python
MIN_ARM = 8              # bp
MAX_ARM = 50             # bp
MAX_LOOP = 12            # bp
MAX_MISMATCHES = 0
SEED_SIZE = 6            # bp
DELTA_G_THRESHOLD = -5.0 # kcal/mol
```

### Triplex Parameters
```python
# H-DNA
MIN_ARM = 10              # bp
MAX_ARM = 100             # bp
MAX_LOOP = 8              # bp
PURITY_THRESHOLD = 0.90
SEED_SIZE = 6             # bp
H_REF_ARM = 35            # Reference arm length
H_LOOP_ALPHA = 0.4        # Loop penalty exponent
H_WEIGHT_L = 0.35         # Arm weight
H_WEIGHT_H = 0.20         # Loop weight
H_WEIGHT_P = 0.30         # Purity weight
H_WEIGHT_I = 0.15         # Interruption weight

# Sticky DNA
MIN_STICKY_COPIES = 4
STICKY_REPLICATION_MIN = 20   # copies
STICKY_STABLE_MIN = 40        # copies
STICKY_PATHOGENIC_MIN = 60    # copies (FRDA)
```

### Slipped DNA Parameters
```python
MIN_TRACT_LENGTH = 20      # bp
MIN_PURITY = 0.90
MAX_UNIT_SIZE = 100        # bp
STR_DIRECT_THRESHOLD = 7   # bp (STR <7, DR ≥7)
MIN_COPIES_STR_CORE = 6    # CORE mode
MIN_COPIES_STR_RELAXED = 4 # LENIENT mode
MIN_COPIES_DR = 2
ENTROPY_MIN = 0.5          # bits

# Disease thresholds (copies)
CAG_HD_PATHOGENIC = 36
CTG_DM1_PATHOGENIC = 50
CGG_FXS_PATHOGENIC = 200
GAA_FRDA_PATHOGENIC = 66
GGGGCC_C9_PATHOGENIC = 30
```

### Hybrid Motifs Parameters
```python
HYBRID_MIN_OVERLAP = 0.50  # 50% minimum overlap fraction
HYBRID_MAX_OVERLAP = 0.99  # 99% maximum (excludes complete overlaps)
# Minimum 2 primary motifs required for hybrid detection
```

### Non-B DNA Clusters Parameters
```python
CLUSTER_WINDOW_SIZE = 300      # bp
CLUSTER_MIN_MOTIFS = 4         # Minimum motifs in window
CLUSTER_MIN_CLASSES = 3        # Minimum different classes
# Minimum 3 primary motifs total required for cluster detection
```

---

## Version Information

| Component | Version | Last Updated |
|-----------|---------|--------------|
| **NonBDNAFinder** | 2025.1 | December 2025 |
| **Documentation** | 1.0 | February 2026 |
| **Detector API** | Stable | Backward compatible |
| **Output Schema** | Minimal v1 | Production |

---

## Additional Resources

- **[README.md](./README.md)**: General overview and quick start
- **[OUTPUT_SCHEMA.md](./OUTPUT_SCHEMA.md)**: Detailed output field documentation
- **[IMPLEMENTATION_SUMMARY.md](./IMPLEMENTATION_SUMMARY.md)**: Technical implementation notes
- **[NOTEBOOK_GUIDE.md](./NOTEBOOK_GUIDE.md)**: Jupyter notebook usage guide
- **[PARALLEL_PROCESSING.md](./docs/PARALLEL_PROCESSING.md)**: Performance optimization details

---

**Document Status:** ✅ Complete and Production-Ready  
**Maintained by:** Dr. Venkata Rajesh Yella  
**Last Review:** February 2026
