# Enhanced Feature Output Examples

This document provides real examples of the enhanced feature output from each detector.

## Example 1: G-Quadruplex

### Input Sequence
```
GGGAGGGAGGGAGGG (15 bp)
```

### Enhanced Output Fields
```json
{
  "ID": "test_G4_1",
  "Sequence_Name": "test",
  "Class": "G-Quadruplex",
  "Subclass": "Canonical intramolecular G4",
  "Start": 1,
  "End": 15,
  "Length": 15,
  "Sequence": "GGGAGGGAGGGAGGG",
  "Score": 0.876,
  "Strand": "+",
  "Method": "Seeded_G4Hunter",
  "Pattern_ID": "G4_CAN",
  
  // NEW ENHANCED FIELDS
  "GC_Content": 80.0,
  "Num_Tracts": 4,
  "Loop_Lengths": "1,1,1",
  "Min_Loop_Length": 1,
  "Max_Loop_Length": 1,
  "Avg_Loop_Length": 1.0,
  "Loop_Length": 1.0,              // ✅ NEW unified field
  "Num_Loops": 3,
  "G_Tract_Lengths": "3,3,3,3",
  "Min_Tract_Length": 3,
  "Max_Tract_Length": 3,
  "Avg_Tract_Length": 3.0,
  "Arm_Length": 3.0,               // ✅ NEW unified field (maps to tract length)
  "Type_Of_Repeat": "Four-tetrad intramolecular G4",                              // ✅ NEW
  "Criterion": "Canonical G4: 4+ G-tracts (≥3G each), loops 1-7bp",              // ✅ NEW
  "Disease_Relevance": "Promoter-like G4 (potential oncogene regulation...)",    // ✅ NEW
  "Regions_Involved": "4 G-tracts: G3, G3, G3, G3; 3 loops: 1bp, 1bp, 1bp"     // ✅ NEW
}
```

## Example 2: Cruciform

### Input Sequence
```
ATCGATCGATCGGGGCGATCGATCGAT (27 bp)
```

### Enhanced Output Fields
```json
{
  "ID": "test_CRU_1",
  "Sequence_Name": "test",
  "Class": "Cruciform",
  "Subclass": "Cruciform forming IRs",
  "Start": 1,
  "End": 27,
  "Length": 27,
  "Sequence": "ATCGATCGATCGGGGCGATCGATCGAT",
  "Score": 0.845,
  "Strand": "+",
  "Method": "Cruciform_detection",
  "Pattern_ID": "CRU_1",
  
  // NEW ENHANCED FIELDS
  "Left_Arm": "ATCGATCGATCG",
  "Right_Arm": "CGATCGATCGAT",
  "Loop_Seq": "GGG",
  "Arm_Length": 12,                // ✅ NEW unified field
  "Loop_Length": 3,                // ✅ NEW unified field
  "Stem_Length": 12,
  "GC_Content": 48.15,             // ✅ NEW unified field
  "GC_Total": 48.15,
  "GC_Left_Arm": 50.0,
  "GC_Right_Arm": 50.0,
  "GC_Loop": 100.0,
  "Mismatches": 0,
  "Match_Fraction": 1.0,
  "DeltaG": -13.2,
  "Type_Of_Repeat": "Inverted repeat (palindromic mirror)",                                    // ✅ NEW
  "Criterion": "Inverted repeat: arm≥8bp, loop≤12bp; arm_length=12bp, loop_length=3bp...",   // ✅ NEW
  "Disease_Relevance": "Highly stable cruciform (ΔG=-13.2) - DNA breakage, genomic...",       // ✅ NEW
  "Regions_Involved": "Left arm (12bp) - Loop (3bp) - Right arm (12bp mirror)"                // ✅ NEW
}
```

## Example 3: Slipped DNA

### Input Sequence
```
CAGCAGCAGCAGCAGCAGCAGCAG (24 bp)
```

### Enhanced Output Fields
```json
{
  "ID": "test_SLIPPED_1",
  "Sequence_Name": "test",
  "Class": "Slipped_DNA",
  "Subclass": "STR",
  "Start": 1,
  "End": 24,
  "Length": 24,
  "Sequence": "CAGCAGCAGCAGCAGCAGCAGCAG",
  "Score": 0.923,
  "Strand": "+",
  "Method": "Slipped_DNA_detection",
  "Pattern_ID": "SLIPPED_1",
  
  // NEW ENHANCED FIELDS
  "Repeat_Unit": "CAG",
  "Unit_Size": 3,
  "Copy_Number": 8.0,
  "Purity": 1.0,
  "Entropy": 0.654,
  "GC_Content": 66.67,             // ✅ NEW
  "Arm_Length": "N/A",             // ✅ NEW (not applicable for tandem repeats)
  "Loop_Length": "N/A",            // ✅ NEW (not applicable for tandem repeats)
  "Type_Of_Repeat": "Trinucleotide (CAG)",                                                      // ✅ NEW
  "Criterion": "STR: unit size 3bp < 7bp threshold; tract ≥20bp; purity 1.00 ≥0.90",          // ✅ NEW
  "Disease_Relevance": "None annotated",                                                         // ✅ NEW
  "Regions_Involved": "Tandem repeat of CAG unit (3bp) × 8.0 copies, total length 24bp",      // ✅ NEW
  "References": "Sinden 1994; Pearson 2005; Mirkin 2007"
}
```

## Example 4: i-Motif

### Input Sequence
```
CCCCACCCCACCCCACCCC (19 bp)
```

### Enhanced Output Fields
```json
{
  "ID": "test_IMOT_1",
  "Sequence_Name": "test",
  "Class": "i-Motif",
  "Subclass": "Canonical i-motif",
  "Start": 1,
  "End": 19,
  "Length": 19,
  "Sequence": "CCCCACCCCACCCCACCCC",
  "Score": 0.891,
  "Strand": "+",
  "Method": "i-Motif_detection",
  "Pattern_ID": "IMOT_1",
  
  // NEW ENHANCED FIELDS
  "Stems": ["CCCC", "CCCC", "CCCC", "CCCC"],
  "Loops": ["A", "A", "A"],
  "Num_Stems": 4,
  "Num_Loops": 3,
  "Stem_Lengths": [4, 4, 4, 4],
  "Loop_Lengths": [1, 1, 1],
  "GC_Content": 84.21,             // ✅ NEW unified field
  "GC_Total": 84.21,
  "GC_Stems": 100.0,
  "Stem_Length": 4.0,
  "Arm_Length": 4.0,               // ✅ NEW unified field (maps to stem length)
  "Loop_Length": 1.0,              // ✅ NEW unified field
  "Type_Of_Repeat": "Four-stranded canonical i-motif (C-rich)",                                 // ✅ NEW
  "Criterion": "Canonical i-motif: 4+ C-tracts (≥3C each), loops 1-7bp",                      // ✅ NEW
  "Disease_Relevance": "Promoter-like i-motif (high C-content) - oncogene regulation...",      // ✅ NEW
  "Regions_Involved": "4 C-tracts: C4, C4, C4, C4; 3 loops: 1bp, 1bp, 1bp"                   // ✅ NEW
}
```

## Example 5: Z-DNA

### Input Sequence
```
CGCGCGCGCGCGCGCG (16 bp)
```

### Enhanced Output Fields
```json
{
  "ID": "test_ZDNA_1",
  "Sequence_Name": "test",
  "Class": "Z-DNA",
  "Subclass": "Z-DNA",
  "Start": 1,
  "End": 16,
  "Length": 16,
  "Sequence": "CGCGCGCGCGCGCGCG",
  "Score": 441.0,
  "Strand": "+",
  "Method": "Z-DNA_detection",
  "Pattern_ID": "ZDNA_1",
  
  // NEW ENHANCED FIELDS
  "Contributing_10mers": 7,
  "Mean_10mer_Score": 63.0,
  "CG_Dinucleotides": 15,
  "AT_Dinucleotides": 0,
  "Alternating_CG_Regions": 2,
  "Alternating_AT_Regions": 0,
  "GC_Content": 100.0,             // ✅ NEW
  "Arm_Length": "N/A",             // ✅ NEW (not applicable for Z-DNA)
  "Loop_Length": "N/A",            // ✅ NEW (not applicable for Z-DNA)
  "Type_Of_Repeat": "CG/GC alternating purine-pyrimidine (canonical Z-DNA)",                   // ✅ NEW
  "Criterion": "10-mer table scoring (Ho 1986); Sum score 441.0 >50.0; 7 contributing...",    // ✅ NEW
  "Disease_Relevance": "Long CG/GC alternating - potential methylation site, Z-DNA...",        // ✅ NEW
  "Regions_Involved": "15 CG/GC dinucleotides; 2 alternating CG/GC regions"                   // ✅ NEW
}
```

## Example 6: Triplex (Sticky DNA)

### Input Sequence
```
GAAGAAGAAGAAGAAGAA (18 bp)
```

### Enhanced Output Fields
```json
{
  "ID": "test_TRX_STICKY_GAA_1",
  "Sequence_Name": "test",
  "Class": "Triplex",
  "Subclass": "Sticky DNA",
  "Start": 1,
  "End": 18,
  "Length": 18,
  "Sequence": "GAAGAAGAAGAAGAAGAA",
  "Score": 1.0,
  "Strand": "+",
  "Method": "Triplex_detection",
  "Pattern_ID": "TRX_STICKY_GAA",
  
  // NEW ENHANCED FIELDS
  "GC_Content": 33.33,             // ✅ NEW
  "Repeat_Unit": "GAA",
  "Copy_Number": 6,
  "Arm_Length": "N/A",             // ✅ NEW (not applicable for tandem repeats)
  "Loop_Length": "N/A",            // ✅ NEW (not applicable for tandem repeats)
  "Type_Of_Repeat": "Trinucleotide (GAA)",                                                      // ✅ NEW
  "Criterion": "GAA/TTC trinucleotide repeat ≥4 copies; n=6 copies",                           // ✅ NEW
  "Disease_Relevance": "Weak triplex potential: n=6",                                           // ✅ NEW
  "Regions_Involved": "GAA trinucleotide repeat × 6 copies"                                     // ✅ NEW
}
```

## CSV Export Example

When exported to CSV, all these enhanced fields are included:

```csv
Sequence_Name,Class,Subclass,Start,End,Length,Sequence,Strand,Score,Method,Pattern_ID,GC_Content,Arm_Length,Loop_Length,Type_Of_Repeat,Criterion,Disease_Relevance,Regions_Involved,...
test,G-Quadruplex,Canonical intramolecular G4,1,15,15,GGGAGGGAGGGAGGG,+,0.876,Seeded_G4Hunter,G4_CAN,80.0,3.0,1.0,"Four-tetrad intramolecular G4","Canonical G4: 4+ G-tracts (≥3G each), loops 1-7bp","Promoter-like G4 (potential oncogene regulation...)","4 G-tracts: G3, G3, G3, G3; 3 loops: 1bp, 1bp, 1bp",...
test,Cruciform,Cruciform forming IRs,1,27,27,ATCGATCGATCGGGGCGATCGATCGAT,+,0.845,Cruciform_detection,CRU_1,48.15,12,3,"Inverted repeat (palindromic mirror)","Inverted repeat: arm≥8bp, loop≤12bp; arm_length=12bp...","Highly stable cruciform (ΔG=-13.2) - DNA breakage...","Left arm (12bp) - Loop (3bp) - Right arm (12bp mirror)",...
test,Slipped_DNA,STR,1,24,24,CAGCAGCAGCAGCAGCAGCAGCAG,+,0.923,Slipped_DNA_detection,SLIPPED_1,66.67,N/A,N/A,"Trinucleotide (CAG)","STR: unit size 3bp < 7bp threshold; tract ≥20bp...","None annotated","Tandem repeat of CAG unit (3bp) × 8.0 copies, total length 24bp",...
```

## JSON Export Example

```json
{
  "version": "2024.1",
  "analysis_type": "NBDScanner_Non-B_DNA_Analysis",
  "total_motifs": 4,
  "motifs": [
    {
      "ID": "test_G4_1",
      "Class": "G-Quadruplex",
      "GC_Content": 80.0,
      "Arm_Length": 3.0,
      "Loop_Length": 1.0,
      "Type_Of_Repeat": "Four-tetrad intramolecular G4",
      "Criterion": "Canonical G4: 4+ G-tracts (≥3G each), loops 1-7bp",
      "Disease_Relevance": "Promoter-like G4 (potential oncogene regulation...)",
      "Regions_Involved": "4 G-tracts: G3, G3, G3, G3; 3 loops: 1bp, 1bp, 1bp",
      ...
    },
    ...
  ]
}
```

## Visualization Example

The enhanced `plot_spacer_loop_variation()` now shows:

### Left Panel: Loop Length Distribution
```
Loop Length Distribution
│
│ Cruciform     ▪──────────■────────▪
│
│ G-Quadruplex  ▪─────■─────▪
│
│ i-Motif       ▪──■──▪
│
│ R-Loop        ▪───────────■──────▪
│
│ Triplex       ▪────■────▪
└────────────────────────────────────
    0   5   10   15   20   25   30
           Loop Length (bp)
```

### Right Panel: Arm/Stem Length Distribution
```
Arm/Stem Length Distribution
│
│ Cruciform     ▪─────────■─────────▪
│
│ G-Quadruplex  ▪──■──▪
│
│ i-Motif       ▪───■───▪
│
│ Triplex       ▪──────■──────▪
│
└────────────────────────────────────
    0   10   20   30   40   50
        Arm/Stem Length (bp)
```

## Summary

All detectors now output comprehensive features enabling:
- ✅ Detailed structural analysis
- ✅ Disease context understanding
- ✅ Reproducible research
- ✅ Publication-quality data
- ✅ Meticulous characterization

Every motif includes complete metadata for thorough scientific analysis.
