# Nature-Level Publication Figures - Quick Start

## New Feature: Multi-Panel Publication Figures 🎨

NBDScanner now includes comprehensive multi-panel figure generation for Nature-level publications!

### What's New

✨ **4 Publication-Ready Multi-Panel Figures**
- Figure 1: Comprehensive Overview (2×2 layout)
- Figure 2: Statistical Characterization (2×2 layout)  
- Figure 3: Genomic Landscape (3×1 layout)
- Figure 4: Detailed Subclass Analysis (2×2 layout)

✅ **Publication Standards Met**
- 300 DPI resolution
- Vector formats (PDF/SVG)
- Colorblind-friendly (Wong 2011 palette)
- Nature column widths (89mm, 120mm, 183mm)
- Professional panel labeling (A, B, C, D)

### Quick Start

```bash
# Install dependencies (if needed)
pip install matplotlib seaborn pandas numpy scipy

# Generate all publication figures
python generate_publication_figures.py --output-dir ./figures --format pdf

# Output: figure1_overview.pdf, figure2_statistics.pdf,
#         figure3_landscape.pdf, figure4_subclass.pdf
#         + figure_captions.md
```

### Python Usage

```python
from publication_figures import (
    create_figure_1_overview,
    create_figure_2_statistics,
    create_figure_3_landscape,
    create_figure_4_subclass_analysis
)
import nonbscanner as nbs

# Analyze your sequence
sequence = "ATCG..." # Your DNA sequence
motifs = nbs.analyze_sequence(sequence, "my_sequence")

# Generate figures
fig1 = create_figure_1_overview(motifs, len(sequence))
fig1.savefig('figure1.pdf', dpi=300, bbox_inches='tight')

fig2 = create_figure_2_statistics(motifs)
fig2.savefig('figure2.pdf', dpi=300, bbox_inches='tight')

fig3 = create_figure_3_landscape(motifs, len(sequence))
fig3.savefig('figure3.pdf', dpi=300, bbox_inches='tight')

fig4 = create_figure_4_subclass_analysis(motifs, len(sequence))
fig4.savefig('figure4.pdf', dpi=300, bbox_inches='tight')
```

### Figure Gallery

#### Figure 1: Comprehensive Overview
Multi-panel layout showing:
- (A) Motif class distribution
- (B) Class-subclass hierarchy (nested donut)
- (C) Genomic coverage map

#### Figure 2: Statistical Characterization
Statistical analysis with:
- (A) Score distributions (violin plots)
- (B) Length distributions (box plots)
- (C) Score-length correlation (scatter)
- (D) Summary statistics table

#### Figure 3: Genomic Landscape
Genome-wide visualization:
- (A) Density heatmap across genome
- (B) Manhattan plot showing hotspots
- (C) Cumulative distribution curves

#### Figure 4: Subclass Analysis
Detailed subclass-level analysis:
- (A) Top 15 subclasses (horizontal bars)
- (B) Subclass density heatmap
- (C) Length distributions (KDE curves)
- (D) Motif co-occurrence matrix

### Documentation

📖 **Full Guide**: See `docs/PUBLICATION_FIGURES_GUIDE.md`

Topics covered:
- Figure standards and requirements
- Multi-panel layouts and customization
- Export options (PDF, SVG, PNG)
- Submission checklist for journals
- Troubleshooting common issues
- Advanced usage examples

### Key Features

✅ **Automatic Panel Labeling**: Bold uppercase (A, B, C, D) per Nature style  
✅ **Consistent Styling**: Uniform fonts, colors, and spacing  
✅ **Colorblind-Safe**: Wong (2011) palette verified for accessibility  
✅ **High Resolution**: 300 DPI for publication quality  
✅ **Vector Format**: PDF/SVG for editability and scalability  
✅ **Comprehensive Captions**: Detailed figure descriptions included  
✅ **Automated Generation**: Single command creates all figures  

### Suitable For

Journals including:
- Nature Methods ✓
- Nature Genetics ✓
- Cell ✓
- Science ✓
- PLOS journals ✓
- BMC journals ✓
- Nucleic Acids Research ✓

### Requirements

```python
matplotlib >= 3.5.0
seaborn >= 0.11.0
pandas >= 1.3.0
numpy >= 1.21.0
scipy >= 1.7.0  # Optional, for KDE plots
```

### File Structure

```
NonBDNAFinder/
├── publication_figures.py           # Multi-panel figure functions (850 lines)
├── generate_publication_figures.py  # Automated generation script
├── docs/
│   └── PUBLICATION_FIGURES_GUIDE.md # Comprehensive documentation
└── visualizations.py                # Enhanced with helper functions
```

### Examples

See generated examples:
- `figure1_overview.pdf` - Comprehensive overview
- `figure2_statistics.pdf` - Statistical characterization
- `figure3_landscape.pdf` - Genomic landscape
- `figure4_subclass.pdf` - Subclass analysis
- `figure_captions.md` - Pre-written captions

### Citation

When using these figures in publications, please cite:

```bibtex
@software{nbdscanner2024,
  title = {NBDScanner: Comprehensive Non-B DNA Motif Detection System},
  author = {Yella, Venkata Rajesh},
  year = {2024},
  version = {2024.1.5},
  url = {https://github.com/VRYella/NonBDNAFinder}
}
```

And reference the colorblind-safe palette:
> Wong, B. (2011). Points of view: Color blindness. Nature Methods, 8(441).

### Support

- **Issues**: https://github.com/VRYella/NonBDNAFinder/issues
- **Documentation**: See `docs/` directory
- **Examples**: Run `generate_publication_figures.py`

---

**Version**: 2024.1.5  
**License**: MIT  
**Author**: Dr. Venkata Rajesh Yella
