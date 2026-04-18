# NBDFinder vs NBST: Concordance and Differences

This comparison is based on the provided `comparison_summary.tsv` and `subclass_comparison_summary.tsv` outputs in `examples/data`.

## Overall concordance
- Total matched calls (TP across mapped classes): **27**
- Total NBDFinder-unique calls (not matched to NBST): **220**
- Total NBST-unique calls (not matched to NBDFinder): **69**
- Micro-precision: **0.109**
- Micro-recall: **0.281**
- Micro-F1: **0.157**

## Where tools agree most (class-level)
- Slipped_DNA (DR): F1=0.556, Precision=0.500, Recall=0.625, TP=5
- Slipped_DNA (STR): F1=0.300, Precision=0.900, Recall=0.180, TP=9
- Z-DNA (Z): F1=0.222, Precision=0.333, Recall=0.167, TP=1

## Where tools differ most (class-level)
- Curved_DNA (curved): F1=0.000, Precision=0.000, Recall=0.000, FP=46, FN=1
- Cruciform (MR): F1=0.069, Precision=0.050, Recall=0.111, FP=19, FN=8
- G-Quadruplex (GQ): F1=0.122, Precision=0.070, Recall=0.500, FP=147, FN=11

## Key differences observed
- NBDFinder reports substantially more G-Quadruplex and Curved_DNA candidates than NBST, indicating broader sensitivity or different thresholding/definitions.
- STR subclass shows strong positional overlap for shared calls, but NBST reports many additional STR calls not captured by NBDFinder in this dataset.
- Cruciform and Curved_DNA classes show low class-level concordance under current matching criteria.
- Z-DNA and Direct Repeat show partial overlap, suggesting differences are concentrated in boundary choices and/or motif model specificity.

## Subclass-level highlights
- STR (mapped to STR): F1=0.305, TP=9, Mean Jaccard=0.948
- Direct Repeat (mapped to DR): F1=0.222, TP=1, Mean Jaccard=1.000
- Z-DNA (mapped to Z): F1=0.222, TP=1, Mean Jaccard=0.867
- Extended-loop canonical (mapped to GQ): F1=0.200, TP=4, Mean Jaccard=0.682
- Canonical intramolecular G4 (mapped to GQ): F1=0.167, TP=2, Mean Jaccard=1.000
