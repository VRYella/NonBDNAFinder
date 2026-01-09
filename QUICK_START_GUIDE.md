# Quick Start Guide: New Statistical Features

## 🎯 What's New?

NonBDNAFinder now includes advanced statistical enrichment and structural pattern analysis!

## 🚀 How to Use

### Step 1: Run Analysis
1. Open NonBDNAFinder app: `streamlit run app.py`
2. Upload your sequence (or use example data)
3. Click "Run Complete Motif Analysis"
4. **New**: Wait ~5-10 seconds for statistical analysis to complete

### Step 2: View Statistical Results
1. Go to **"Results"** tab
2. Click on **"Figure 4: Statistical Enrichment & Structural Analysis"**
3. Explore two panels:
   - **Panel G**: Statistical enrichment (p-values, z-scores)
   - **Panel H**: Structural patterns (blocks, hybrids, clusters)

### Step 3: Download Enhanced Data
1. Go to **"Download"** tab
2. Scroll to **"Download Enrichment & Structural Analysis Data"**
3. Choose your export:
   - **Enrichment Analysis (CSV)**: Statistical metrics
   - **Structural Analysis (Excel)**: Blocks, hybrids, clusters
   - **Complete Package (Excel)**: Everything combined

## 📊 What You'll See

### Enrichment Analysis
- **100× shuffle-based null model** for statistical rigor
- **P-values** showing if patterns are statistically significant
- **Z-scores** quantifying enrichment strength
- **Violin plots** comparing observed vs. expected distributions

### Structural Patterns
- **Pattern-rich blocks**: High-density genomic regions
- **Hybrid zones**: Areas where multiple motif classes overlap
- **Clusters**: Spatial groupings of motifs
- **Interaction matrices**: Which classes co-occur most often

## 💡 Quick Tips

1. **Statistical Significance**:
   - P-value < 0.05: Significant enrichment
   - P-value < 0.01: Highly significant enrichment

2. **Structural Insights**:
   - Blocks show where motifs concentrate
   - Hybrids reveal structural interactions
   - Clusters indicate regulatory hotspots

3. **Performance**:
   - Analysis adds ~5-10 seconds per sequence
   - All visualizations are publication-ready (300 DPI)
   - Exports include all detailed metrics

## 🔍 Example Interpretation

**If you see**:
- P-value = 0.02 for pattern count → Your sequence has **significantly more** motifs than expected by chance
- 5 hybrid zones detected → Multiple motif classes **frequently overlap**
- 3 clusters with stability > 0.8 → **Stable regulatory hotspots** identified

## 📚 Learn More

- **INTEGRATION_COMPLETE.md**: Full technical documentation
- **INTEGRATION_GUIDE.md**: Developer integration guide
- **Module documentation**: enrichment.py, structural_analysis.py, visualizations_enhanced.py

## ✅ Verification

Run this test to verify installation:
```bash
python3 -c "
from enrichment import run_enrichment_analysis
from structural_analysis import run_structural_analysis
from nonbscanner import analyze_sequence

test_seq = 'GGGTTAGGGTTAGGGTTAGGG' * 20
motifs = analyze_sequence(test_seq, 'test')
enr = run_enrichment_analysis(test_seq, motifs, n_shuffles=10)
struct = run_structural_analysis(test_seq, motifs)

print('✅ All modules working correctly!')
print(f'Detected {len(motifs)} motifs')
print(f'Enrichment p-value: {enr[\"enrichment_scores\"][\"pattern_count\"][\"pvalue\"]:.4f}')
print(f'Structural features: {struct[\"summary\"][\"num_blocks\"]} blocks')
"
```

Expected output: No errors, motifs detected, p-value calculated

---

**Enjoy your enhanced NonBDNAFinder with statistical rigor!** 🎉
