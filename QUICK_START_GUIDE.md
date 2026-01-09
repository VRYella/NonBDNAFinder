# Quick Start Guide: NonBDNAFinder

## 🎯 What's This Tool?

NonBDNAFinder is a comprehensive tool for detecting and analyzing Non-B DNA structures in genomic sequences. It identifies 11 motif classes with high accuracy and provides publication-ready visualizations.

## 🚀 How to Use

### Step 1: Install Dependencies
```bash
pip install -r requirements.txt
```

### Step 2: Run Analysis
1. Open NonBDNAFinder app: `streamlit run app.py`
2. Upload your sequence (FASTA format) or use example data
3. Click "Run Complete Motif Analysis"
4. Wait for analysis to complete (~13,000 bp/s typical speed)

### Step 3: View Results
1. Go to **"Results"** tab
2. Explore multiple visualization panels:
   - **Figure 1**: Global landscape of detected motifs
   - **Figure 2**: Clustering and co-occurrence patterns
   - **Figure 3**: Structural constraints (optional)

### Step 4: Download Data
1. Go to **"Download"** tab
2. Choose your export format:
   - **CSV**: Simple tabular format
   - **BED**: Genome browser format
   - **Excel**: Multi-sheet workbook with statistics
   - **JSON**: Structured data format

## 📊 What You'll See

### Motif Detection
- **11 motif classes**: A-philic, Curved DNA, G-Quadruplex, i-Motif, Z-DNA, R-loops, Cruciform, Triplex, Slipped DNA, Hybrid, and Clusters
- **Confidence scores**: 1-3 scale for each detection
- **Detailed annotations**: Position, length, strand, and motif-specific features

### Visualizations
- **Distribution plots**: See where motifs are located
- **Density heatmaps**: Identify hotspots
- **Co-occurrence matrices**: Understand motif interactions
- **Class/subclass breakdowns**: Detailed classification

## 💡 Quick Tips

1. **File Formats**:
   - Supports FASTA, multi-FASTA
   - Handles large files (100+ MB) automatically
   - Accepts compressed files (.gz)

2. **Performance**:
   - Typical speed: ~13,000 bp/s
   - Large sequences processed in chunks
   - All visualizations are publication-ready (300 DPI)

3. **Understanding Results**:
   - Score 3.0: High confidence detection
   - Score 2.0-2.99: Medium confidence
   - Score 1.0-1.99: Low confidence but detectable

## 🔍 Example Interpretation

**If you see**:
- 50 G-Quadruplex motifs → Your sequence has strong G-quadruplex forming potential
- High density in promoter region → Potential regulatory significance
- Multiple hybrid zones → Complex structural interactions

## 📚 Learn More

- **README.md**: Full documentation and features
- **MODERNIZATION_SUMMARY.md**: Architecture details
- **Pattern registry**: pattern_registry2.xlsx

## ✅ Verification

Run this test to verify installation:
```bash
python3 -c "
from nonbscanner import analyze_sequence

test_seq = 'GGGTTAGGGTTAGGGTTAGGG' * 20
motifs = analyze_sequence(test_seq, 'test')

print('✅ NonBDNAFinder working correctly!')
print(f'Detected {len(motifs)} motifs')
print(f'Classes found: {set(m[\"Class\"] for m in motifs)}')
"
```

Expected output: No errors, motifs detected successfully

---

**Enjoy your NonBDNAFinder analysis!** 🎉
