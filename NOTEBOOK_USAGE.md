# NonBDNAFinder Jupyter Notebook

## Quick Start

The `NonBDNAFinder_Demo.ipynb` notebook provides a standalone demonstration of the NonBDNAFinder tool with 3 independent cells.

### Prerequisites

1. **Install Jupyter** (if not already installed):
   ```bash
   pip install jupyter
   ```

2. **Install Dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

3. **Prepare a FASTA file**: You'll need a DNA sequence in FASTA format to analyze.

### Running the Notebook

1. **Launch Jupyter**:
   ```bash
   jupyter notebook NonBDNAFinder_Demo.ipynb
   ```

2. **Run Cell 1**: Upload your FASTA file using the file upload widget
   - Click the "Upload FASTA" button
   - Select your FASTA file (.fasta, .fa, .fna, or .txt)
   - Wait for the upload confirmation and sequence preview

3. **Run Cell 2**: Analyze the uploaded sequence

4. **Run Cell 3**: View detailed results

### Notebook Structure

The notebook contains **3 standalone code cells** that should be run in sequence:

#### Cell 1: File Upload
- Displays a file upload widget
- Accepts FASTA format files (.fasta, .fa, .fna, .txt)
- Parses the uploaded file using Biopython
- Shows sequence information and preview
- **Note**: If the file contains multiple sequences, only the first one will be used

#### Cell 2: Analysis
- Runs NonBDNAFinder analysis on the uploaded sequence
- Detects all Non-B DNA motifs (11 classes)
- Displays motif counts by class

#### Cell 3: Results Display
- Shows detailed results in a formatted table
- Provides summary statistics
- Exports results to CSV file

### Creating a FASTA File

If you don't have a FASTA file, you can create one with this format:

```
>My_DNA_Sequence
GGGTTAGGGTTAGGGTTAGGG
AAAAATAAAAAAAAAAA
CCCCCCCCCCCCCCCC
CGCGCGCGCGCGCGCG
```

Save this text in a file with a `.fasta` or `.fa` extension.

### Expected Output

The analysis results will vary based on your input sequence. The output includes:
- Number of motifs detected
- Motif classes found (from 11 possible classes)
- Classes may include: G-Quadruplex, Curved DNA, i-Motif, Z-DNA, R-Loop, Slipped DNA, A-philic DNA, Hybrid motifs, and DNA Clusters
- Results exported to `nonbdna_results.csv`

### Notes

- **File Upload Only**: This notebook uses a file upload widget as the only input method
- Each cell is designed to work in sequence (run Cell 1 first, then Cell 2, then Cell 3)
- The cells use standard Python libraries (pandas, ipywidgets, biopython, nonbscanner)
- Requires `ipywidgets` for the file upload functionality
- Maximum file size depends on your Jupyter configuration and available memory

### Troubleshooting

**Widget not displaying?**
- Make sure ipywidgets is installed: `pip install ipywidgets>=8.0.0`
- Enable the widget extension: `jupyter nbextension enable --py widgetsnbextension`

**Upload button not working?**
- Try restarting the Jupyter kernel
- Check that your file is in valid FASTA format

**Parse error?**
- Ensure your file uses the standard FASTA format with `>` prefix for headers
- Check that the file encoding is UTF-8

For more information, see the main [README.md](README.md).
