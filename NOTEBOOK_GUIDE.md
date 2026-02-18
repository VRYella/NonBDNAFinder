# NonBDNAFinder Jupyter Notebooks Guide

## Overview

NonBDNAFinder provides two Jupyter notebooks for different use cases:

1. **NonBDNAFinder_SingleCell.ipynb** - Production-ready single-cell batch processor
2. **NonBDNAFinder_Demo.ipynb** - Interactive step-by-step demo

---

## NonBDNAFinder_SingleCell.ipynb

### Purpose
A minimal, efficient single-cell notebook for batch processing multiple DNA sequence files on local machines.

### Key Features
- ✅ **Single cell execution** - Edit config, run once, get results
- ✅ **Batch processing** - Process multiple files automatically
- ✅ **Flexible input** - Single file, list of files, or glob patterns
- ✅ **Automatic export** - CSV results with summary statistics and PDF visualizations
- ✅ **Production ready** - Minimal code, maximum efficiency

### Usage

1. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Open the notebook:**
   ```bash
   jupyter notebook NonBDNAFinder_SingleCell.ipynb
   ```

3. **Edit the configuration** (lines 36-40):
   ```python
   # Single file
   file_patterns = 'genome.fasta'
   
   # Multiple files
   file_patterns = ['file1.fasta', 'file2.fasta']
   
   # Glob pattern (all FASTA files in directory)
   file_patterns = 'data/*.fasta'
   
   # Recursive glob (all .fa files in subdirectories)
   file_patterns = 'data/**/*.fa'
   
   # Output prefix
   output_prefix = 'nonbdna'  # Creates nonbdna_results.csv, nonbdna_summary.csv, and nonbdna_visualization.pdf
   ```

4. **Run the cell** (Shift+Enter)

5. **Check outputs:**
   - `{output_prefix}_results.csv` - All detected motifs with full details
   - `{output_prefix}_summary.csv` - Summary statistics by motif class
   - `{output_prefix}_visualization.pdf` - Visual plots and charts (multi-page PDF)

### Input Format
- Accepts standard FASTA files (.fasta, .fa, .fna)
- Can process multiple sequences per file
- Handles both single and multi-sequence FASTA files

### Output Format

#### Detailed Results (`*_results.csv`)
Contains all detected motifs with columns:
- Sequence_Name, Class, Subclass, Start, End, Length, Strand, Score, Method, Pattern_ID
- Source_File (automatically added)
- Additional motif-specific columns (e.g., Num_Tracts for G-Quadruplex)

#### Summary Statistics (`*_summary.csv`)
Aggregated statistics by motif class:
- Count - Total motifs detected
- Avg_Length, Min_Length, Max_Length
- Avg_Score, Min_Score, Max_Score

### Example Output
```
📁 Found 2 file(s) to process
   • example.fasta
   • genome.fasta

================================================================================
Processing file 1/2: example.fasta
================================================================================
📊 Found 2 sequence(s) in file

  🔬 Analyzing: test_sequence_1 (362 bp)
  ✓ Detected 29 motifs

  🔬 Analyzing: test_sequence_2 (183 bp)
  ✓ Detected 12 motifs

================================================================================
📊 ANALYSIS COMPLETE
================================================================================

✓ Total motifs detected: 41
✓ Motif classes found: 11

Motifs by class:
  • Non-B_DNA_Clusters: 15
  • Triplex: 6
  • G-Quadruplex: 5
  ...

💾 Detailed results saved to: nonbdna_results.csv
💾 Summary statistics saved to: nonbdna_summary.csv

📊 Generating visualizations...
📊 Visualizations saved to: nonbdna_visualization.pdf

✨ Analysis complete! Check the output files for results.
```

### Error Handling
The notebook includes comprehensive error handling:
- File not found errors with helpful messages
- Invalid FASTA format detection
- Per-file error handling (continues processing other files)
- Clear error messages with troubleshooting steps

### Performance
- Processes ~13,000-20,000 bp/second
- Handles genomes up to 200MB+
- Memory-efficient streaming for large files
- Parallel processing for sequences >50KB

---

## NonBDNAFinder_Demo.ipynb

### Purpose
An interactive step-by-step demonstration of NonBDNAFinder with educational content.

### Key Features
- 📤 Interactive file upload widget
- 📚 Detailed explanations and documentation
- 🔬 Cell-by-cell execution
- 📊 Results visualization
- 🎓 Educational comments

### Usage
Perfect for:
- Learning how NonBDNAFinder works
- Exploring different features
- Interactive analysis sessions
- Teaching and demonstrations

---

## Comparison

| Feature | SingleCell | Demo |
|---------|-----------|------|
| **Cells** | 1 (single cell) | 3 (multiple cells) |
| **Use Case** | Production batch processing | Learning & exploration |
| **Input Method** | Edit config variable | File upload widget |
| **Batch Processing** | ✅ Yes | ❌ No (single file) |
| **Glob Patterns** | ✅ Yes | ❌ No |
| **Auto Export** | ✅ Yes (CSV + PDF) | ✅ Yes (CSV) |
| **Documentation** | Minimal (succinct) | Extensive (educational) |
| **Best For** | Automated workflows | Interactive analysis |

---

## Troubleshooting

### Common Issues

**Issue: "No files found matching pattern"**
- Check file path is correct
- Verify file exists in specified location
- Ensure glob pattern syntax is correct (use `*.fasta` not `*.fasta/`)

**Issue: "ModuleNotFoundError"**
- Install dependencies: `pip install -r requirements.txt`
- Ensure you're in the repository root directory

**Issue: "Invalid DNA sequence"**
- FASTA file must contain valid DNA characters (ATCG)
- Check for invalid characters or formatting issues
- Ensure file is in proper FASTA format

**Issue: "No motifs detected"**
- Sequences may be too short (need sufficient length for motifs)
- Sequences may lack Non-B DNA forming patterns
- This is normal for some sequences

### Getting Help

1. Check the error message - they include troubleshooting steps
2. Review the README.md for usage examples
3. Consult the full documentation in the repository
4. Report issues on GitHub

---

## Tips & Best Practices

### Performance Optimization
1. Use glob patterns for batch processing multiple files
2. Process multiple sequences in one execution
3. For very large files (>100MB), the tool automatically uses chunking
4. Consider parallelization for multi-file processing

### File Organization
```
project/
├── data/
│   ├── genome1.fasta
│   ├── genome2.fasta
│   └── ...
├── NonBDNAFinder_SingleCell.ipynb
└── results/
    ├── analysis1_results.csv
    └── analysis1_summary.csv
```

### Workflow Example
```python
# 1. Process all FASTA files in a directory
file_patterns = 'data/*.fasta'
output_prefix = 'batch_analysis'

# 2. Run cell - generates:
#    - batch_analysis_results.csv
#    - batch_analysis_summary.csv

# 3. Further analysis in pandas
import pandas as pd
df = pd.read_csv('batch_analysis_results.csv')

# Filter for specific motif classes
g4_motifs = df[df['Class'] == 'G-Quadruplex']

# Analyze by source file
by_file = df.groupby('Source_File')['Class'].value_counts()
```

---

## Advanced Usage

### Custom Output Analysis
```python
# After running the single-cell notebook
import pandas as pd
import matplotlib.pyplot as plt

# Load results
df = pd.read_csv('nonbdna_results.csv')

# Visualize motif distribution
df['Class'].value_counts().plot(kind='bar')
plt.title('Motif Class Distribution')
plt.xlabel('Motif Class')
plt.ylabel('Count')
plt.tight_layout()
plt.savefig('motif_distribution.png')

# Analyze by chromosome/sequence
by_seq = df.groupby('Sequence_Name')['Class'].value_counts()
print(by_seq)
```

### Integration with Workflows
```bash
# 1. Run notebook programmatically
jupyter nbconvert --execute NonBDNAFinder_SingleCell.ipynb \
    --to notebook --inplace

# 2. Process outputs
python analyze_results.py

# 3. Generate report
python generate_report.py
```

---

## License

MIT License - See LICENSE file for details

## Citation

If you use NonBDNAFinder in your research, please cite:
```bibtex
@software{nonbdnafinder2025,
  author = {Yella, Venkata Rajesh},
  title = {NonBDNAFinder 2025.1: Nobel-Level Quality DNA Motif Detection},
  year = {2025},
  url = {https://github.com/VRYella/NonBDNAFinder}
}
```
