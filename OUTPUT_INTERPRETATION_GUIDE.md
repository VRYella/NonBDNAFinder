# Output Interpretation Guide
## Understanding NonBDNAFinder Results - Detailed Motif Information

**Version:** 2025.1  
**Last Updated:** December 29, 2025

---

## Quick Reference: Score Interpretation

### Confidence Tiers (1-3 Scale)

| Score Range | Confidence | Interpretation | Action |
|-------------|-----------|----------------|--------|
| **2.5 - 3.0** | 🟢 HIGH | Strong evidence for motif formation | High priority for experimental validation |
| **1.5 - 2.5** | 🟡 MODERATE | Likely motif with moderate confidence | Consider for validation based on context |
| **1.0 - 1.5** | 🟠 LOW | Possible motif, lower confidence | Secondary priority, may need additional evidence |
| **< 1.0** | 🔴 VERY LOW | Weak evidence | Use with caution, likely false positive |

---

## Understanding Motif Detail Fields

### Core Information (Always Present)

#### 1. **Class**
The primary motif category (11 classes total):
- `G-Quadruplex`: Four-stranded G-rich structures
- `Slipped_DNA`: Direct repeats and STRs
- `Z-DNA`: Left-handed double helix
- `Cruciform`: Inverted palindromic repeats
- `i-Motif`: C-rich four-stranded structures
- `R-Loop`: RNA-DNA hybrids
- `Triplex`: Three-stranded structures
- `Curved_DNA`: Intrinsically bent DNA
- `A-philic DNA`: A-tract regions with high protein affinity
- `Hybrid`: Overlapping different motif classes
- `Non-B_DNA_Clusters`: High-density regions

#### 2. **Subclass**
More specific variant within the class:
- Example (G4): `Canonical intramolecular G4`, `Telomeric G4`, `Extended-loop canonical`
- Example (Slipped): `STR`, `Direct Repeat`
- Example (Z-DNA): `Z-DNA`, `eGZ (Extruded-G) DNA`

#### 3. **Position**
- `Start`: 1-based start position in the sequence
- `End`: End position (inclusive)
- `Length`: Motif length in base pairs

#### 4. **Score**
Universal confidence score (1-3 scale):
- Based on thermodynamic properties, sequence characteristics, and pattern strength
- Normalized across all motif classes for comparability
- Higher score = stronger evidence for motif formation

#### 5. **Method**
Detection algorithm used:
- `G4Hunter_detection`: G4-quadruplex scoring algorithm
- `Slipped_DNA_detection`: Repeat unit identification
- `Z-DNA_detection`: 10-mer scoring for Z-DNA potential
- `i-Motif_detection`: C-run pattern matching
- `Hybrid_Detection`: Multi-class overlap analysis
- `Cluster_Detection`: Density-based clustering

#### 6. **Pattern_ID**
Internal pattern identifier from registry:
- Format: `{CLASS}_{NUMBER}` (e.g., `G4_CAN`, `SLIPPED_1`, `ZDNA_1`)
- Links to specific detection pattern in registry

---

## Class-Specific Output Fields

### G-Quadruplex Output

**Structural Components:**
- `Stems`: G-rich tracts (e.g., `GGG, GGG, GGG, GGG`)
- `Loops`: Connecting sequences between stems
- `Num_Stems`: Number of G-tracts (usually 4)
- `Num_Loops`: Number of connecting loops (usually 3)
- `Stem_Lengths`: Length of each G-tract
- `Loop_Lengths`: Length of each loop

**Composition:**
- `GC_Total`: Overall GC content (%)
- `GC_Stems`: GC content in G-tracts (usually 100%)
- `Stem_Length`: Average stem length
- `Loop_Length`: Average loop length

**Example Output:**
```
Class: G-Quadruplex
Subclass: Canonical intramolecular G4
Score: 2.029 (HIGH CONFIDENCE)
Stems: GGG, GGG, GGG, GGG
Loops: T, TT, T
Num_Stems: 4
GC_Stems: 100.0%
```

---

### Slipped DNA (STR) Output

**Repeat Characterization:**
- `Repeat_Unit`: The repeating sequence (e.g., `CAG`, `CGG`, `CTG`)
- `Unit_Size`: Length of repeat unit in bp
- `Copy_Number`: Number of repeat copies (float, allows partial repeats)
- `Purity`: Repeat purity (0-1), where 1.0 = perfect repeat

**Scoring:**
- `Slippage_Energy_Score`: Energy-based score for slippage potential
- Higher purity + more copies = higher score

**Medical Relevance:**
| Repeat | Disease | Normal Range | Pathogenic Range |
|--------|---------|--------------|------------------|
| CAG | Huntington's | 6-35 | >36 |
| CGG | Fragile X | 5-44 | >200 |
| CTG | Myotonic Dystrophy | 5-37 | >50 |
| GAA | Friedreich's Ataxia | 5-33 | >66 |

**Example Output:**
```
Class: Slipped_DNA
Subclass: STR
Score: 2.581 (HIGH CONFIDENCE)
Repeat_Unit: CAG
Copy_Number: 20.0
Purity: 1.000 (perfect repeat)
Slippage_Energy_Score: 2.581
```

**Interpretation:** 20 CAG repeats = normal Huntingtin range (healthy)

---

### Z-DNA Output

**10-mer Analysis:**
- `Contributing_10mers`: Number of 10-mer windows analyzed
- `Mean_10mer_Score`: Average Z-DNA propensity score per 10-mer
- `Raw_Score`: Sum of all 10-mer scores

**Sequence Characteristics:**
- `CG_Dinucleotides`: Count of CG dinucleotides
- `AT_Dinucleotides`: Count of AT dinucleotides
- `Alternating_CG_Regions`: Number of (CG)n regions
- `Alternating_AT_Regions`: Number of (AT)n regions
- `GC_Content`: Overall GC percentage

**Z-DNA Potential Indicators:**
- High alternating CG/CA content
- High mean 10-mer score (>50)
- High overall raw score

**Example Output:**
```
Class: Z-DNA
Subclass: Z-DNA
Score: 2.420 (HIGH CONFIDENCE)
Contributing_10mers: 11
Mean_10mer_Score: 63.0
CG_Dinucleotides: 19
Alternating_CG_Regions: 2
GC_Content: 100.0%
Raw_Score: 693.0
```

---

### i-Motif Output

**Similar to G4, but for C-rich sequences:**
- `Stems`: C-rich tracts (e.g., `CCC, CCC, CCC, CCC`)
- `Loops`: Connecting sequences
- `Num_Stems`: Number of C-tracts
- `GC_Stems`: GC content in C-tracts (usually 100%)
- `Raw_Score`: i-Motif formation score

**pH Sensitivity:**
i-Motifs are pH-sensitive structures that form preferentially at slightly acidic pH (~6.5)

---

### Cruciform Output

**Structural Features:**
- `Arm_Length`: Length of palindromic arms
- `Loop_Length`: Size of hairpin loops
- `Num_Stems`: Number of stem regions (usually 2)
- `AT_Content`: A+T percentage (affects stability)

**Stability Factors:**
- Longer arms = more stable
- High AT content = easier to unwind (more cruciform-prone)

---

### Triplex DNA Output

**Mirror Repeat Properties:**
- `Arm_Length`: Length of mirror repeat arms
- `Loop_Length`: Spacer between arms
- `Purine_Content`: Percentage of purines (A+G)
- `Pyrimidine_Content`: Percentage of pyrimidines (C+T)

**Triplex Formation Requirements:**
- High purity in purine or pyrimidine content (>90%)
- Sufficient arm length (>10 bp)

---

### Curved DNA Output

**A-tract Analysis:**
- `Tract_Type`: Type of tract (`A-tract` or `T-tract`)
- `Tract_Length`: Length of A/T tract
- `AT_Content`: Percentage of A+T
- `Raw_Score`: Curvature propensity score

**Phasing:**
- `Phasing_Quality`: Whether A-tracts are in-phase with helix repeat (~10-11 bp)
- In-phase A-tracts cause cumulative bending

---

### R-Loop Output

**Thermodynamic Properties:**
- `GC_Skew`: (G-C)/(G+C) on template strand
- `RIZ_Length`: RNA-DNA hybrid initiation zone length
- `REZ_Length`: RNA-DNA hybrid extension zone length

**Formation Potential:**
- Negative GC skew favors R-loop formation
- Longer RIZ/REZ = more stable R-loop

---

## Interpreting Hybrid and Cluster Motifs

### Hybrid Motifs
**Definition:** Regions where different motif classes overlap (50-99% overlap)

**Output Fields:**
- `Component_Classes`: List of overlapping classes (e.g., `[Z-DNA, Slipped_DNA]`)
- `Score`: Average of component scores

**Biological Significance:**
Hybrid regions may represent:
- Competition between alternative structures
- Context-dependent conformational switches
- Regulatory elements with multiple functional states

---

### Non-B DNA Clusters
**Definition:** High-density regions with ≥4 motifs from ≥3 different classes within 300 bp window

**Output Fields:**
- `Motif_Count`: Total motifs in cluster
- `Class_Diversity`: Number of different motif classes
- `Component_Classes`: List of classes present

**Biological Significance:**
- Hotspots of structural instability
- Potential sites for chromosomal rearrangements
- Regulatory hubs with multiple structural switches

---

## Using the Enhanced Motif Detail Viewer

### Accessing Detailed Information

1. **Navigate to Results Tab**
   - Run analysis first in "Upload & Analyze" tab
   - Switch to "Results" tab

2. **Select Motif from Dropdown**
   - Use dropdown showing: `Class | Subclass | Position | Score`
   - Example: `G-Quadruplex | Canonical intramolecular G4 | Pos: 1-21 | Score: 2.03`

3. **View Organized Information**
   - **Metrics**: Class, Subclass, Position, Score, Length
   - **Confidence Tier**: HIGH/MODERATE/LOW classification
   - **Sequence**: DNA sequence with smart truncation for long motifs
   - **Detection Method**: Algorithm and pattern ID used

4. **Enable "Show all fields"**
   - Reveals categorized detailed fields:
     - **Scoring Components**: All score-related calculations
     - **Structural Properties**: Stems, loops, tracts, repeats
     - **Additional Properties**: Class-specific metadata

5. **Read Biological Context**
   - Class-specific description of biological roles
   - Disease associations where relevant
   - Functional implications

---

## Tips for Interpretation

### 1. Prioritize by Confidence
Start with HIGH confidence motifs (score ≥2.5) for experimental validation

### 2. Consider Biological Context
- CAG/CGG/CTG repeats: Check copy number against disease thresholds
- G4s in promoters: May regulate gene expression
- Z-DNA near transcription start sites: Transcriptional regulation

### 3. Check Structural Features
- Perfect STRs (purity = 1.0) are most likely to cause disease
- G4s with 4 stems and short loops are most stable
- Longer cruciform arms indicate more stable structures

### 4. Review Detection Method
- Different methods have different strengths
- G4Hunter: Thermodynamically-grounded G4 scoring
- K-mer indexing: Fast, comprehensive repeat detection
- 10-mer scoring: Sensitive to local sequence properties

### 5. Cross-reference with Visualizations
- Use Manhattan plots to identify density hotspots
- Check co-occurrence matrix for motif clustering
- Examine GC correlation plots for sequence bias

---

## Common Questions

**Q: Why do some motifs have scores >3?**  
A: Raw scores before normalization may exceed 3. The displayed "Score" field is normalized to 1-3 scale.

**Q: Can a sequence contain multiple overlapping motifs?**  
A: Yes. Overlap resolution removes duplicates *within* the same class, but different classes can overlap (shown as Hybrid motifs).

**Q: What does "Purity: 1.000" mean for STRs?**  
A: Perfect repeat with no interruptions (e.g., `(CAG)20` = 100% pure CAG repeat)

**Q: Why are some expected motifs not detected?**  
A: Possible reasons:
- Below minimum length threshold for that class
- Score below confidence cutoff
- Another motif type scored higher in overlap region

**Q: How reliable are LOW confidence motifs?**  
A: Use with caution. They may represent:
- Weak/partial motifs
- Edge cases that barely meet criteria
- Context-dependent structures

---

## Export Options for Detailed Data

### CSV Export
- All core columns + selected additional columns
- Suitable for spreadsheet analysis

### Excel Export
- Multiple sheets: Consolidated + per-class detailed sheets
- All fields preserved in class-specific sheets
- Recommended for comprehensive analysis

### JSON Export
- All fields preserved with metadata
- Suitable for programmatic analysis

### BED Format
- Genome browser compatible
- Position and score information only

---

## Further Resources

- **README.md**: Overview and quick start
- **OUTPUT_SCHEMA.md**: Detailed column descriptions
- **VALIDATION_REPORT.md**: Accuracy validation results
- **Consolidated_registry.json**: Pattern database

---

**End of Output Interpretation Guide**
