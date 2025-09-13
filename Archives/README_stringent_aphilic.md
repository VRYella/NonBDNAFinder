# Stringent A-philic DNA Detection

This directory contains `check_a_philic_stringent.py`, a tool for stringent A-philic 10-mer detection using tri/tetra propensity tables.

## Usage

### Basic usage with a sequence:
```bash
python check_a_philic_stringent.py \
  --tetra example_tetra_propensities.csv \
  --tri example_tri_propensities.csv \
  --seq "GATCGAAAATTTTATTTAAATTTAAATTTGGGTTA"
```

### Process FASTA file:
```bash
python check_a_philic_stringent.py \
  --tetra example_tetra_propensities.csv \
  --tri example_tri_propensities.csv \
  --fasta input.fasta \
  --out results.tsv
```

### Adjust stringency parameters:
```bash
python check_a_philic_stringent.py \
  --tetra example_tetra_propensities.csv \
  --tri example_tri_propensities.csv \
  --seq "AAAATTTTAAAATTTT" \
  --min_step_value 0.5 \
  --min_sum_tetra 10.0 \
  --min_sum_tri 8.0 \
  --nucleation_peak 2.0
```

## Algorithm Steps

1. **Load propensity tables** from CSV files
2. **Score sequence** with tri/tetra propensities  
3. **Apply stringent validation** to each 10-mer:
   - All tri/tetra steps must exceed min_step_value
   - Sum of tetra/tri scores must exceed thresholds
   - Must have nucleation peak (high tri score)
4. **Extend regions** rightward while maintaining criteria
5. **Report non-overlapping** maximal regions

## Parameters

- `--min_step_value`: Minimum score for every tri/tetra step (default: 0.0)
- `--min_sum_tetra`: Minimum sum of tetra scores across 10-mer (default: 5.0)
- `--min_sum_tri`: Minimum sum of tri scores across 10-mer (default: 4.0)
- `--nucleation_peak`: Minimum single tri peak for nucleation (default: 1.5)
- `--min_region_len`: Minimum reported region length (default: 10)

## Input Format

CSV files with columns: `Step,Log2_Odds`
- Comma or tab/space separated
- Missing steps treated as 0.0
- Comments start with `#`

## Output Format

TSV with columns: `chrom,start,end,length,sum_tet,sum_tri,tet_min,tri_min,nucleation_max_tri,sequence`