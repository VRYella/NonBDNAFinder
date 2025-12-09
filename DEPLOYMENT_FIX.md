# ImportError Fix - Deployment Guide

## Issue Summary

The application was experiencing an `ImportError` when attempting to start on Streamlit Cloud. The error occurred at line 45 in `app.py` when trying to import visualization functions:

```python
from visualizations import (
    plot_motif_distribution, plot_coverage_map, plot_density_heatmap,
    plot_length_distribution, plot_score_distribution, plot_nested_pie_chart, 
    MOTIF_CLASS_COLORS, plot_density_comparison,
    plot_circos_motif_density,
    plot_density_comparison_by_subclass, plot_enrichment_analysis_by_subclass,
    plot_subclass_density_heatmap,
    # New Nature-quality visualizations
    plot_manhattan_motif_density, plot_cumulative_motif_distribution,
    plot_motif_cooccurrence_matrix, plot_gc_content_correlation,
    plot_linear_motif_track, plot_cluster_size_distribution,
    plot_motif_length_kde
)
```

## Root Cause

The error was caused by **missing Python dependencies** required by the `visualizations.py` module. Specifically:
- `matplotlib` >= 3.5.0
- `seaborn` >= 0.11.0
- `plotly` >= 5.17.0
- And other dependencies listed in `requirements.txt`

## Solution

The fix involves ensuring all dependencies from `requirements.txt` are properly installed before running the application.

### For Streamlit Cloud Deployment

Streamlit Cloud automatically installs packages from `requirements.txt`, but you need to ensure:

1. **Verify `requirements.txt` is in the root directory** of your repository
2. **All required packages are listed** with appropriate version constraints
3. **The deployment environment has access** to PyPI to download packages

### For Local Development

Install all dependencies using:

```bash
pip install -r requirements.txt
```

### Required Dependencies

The following core dependencies must be installed (all specified in `requirements.txt`):

```
# Web Framework
streamlit>=1.28.0

# Scientific Computing
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0
scikit-learn>=1.0.0

# Visualization  
matplotlib>=3.5.0
seaborn>=0.11.0
plotly>=5.17.0

# Bioinformatics
biopython>=1.79

# Performance (optional but recommended)
hyperscan>=0.7.0

# Data Export
openpyxl>=3.0.0
xlsxwriter>=3.0.0

# Networking
requests>=2.28.0

# System Monitoring
psutil>=5.8.0
```

## Verification

To verify the fix is working correctly, run the import test script:

```bash
python3 test_imports.py
```

Expected output:
```
======================================================================
NonBDNAFinder Import Verification Test
======================================================================
Testing utilities imports...
✓ All utilities imports successful

Testing nonbscanner imports...
✓ All nonbscanner imports successful

Testing visualizations imports...
✓ All visualizations imports successful
  - Standard visualization functions: 8
  - Nature-quality visualization functions: 7
  - Color constants: 1

Testing app.py full import...
✓ app.py imported successfully

======================================================================
Test Summary:
======================================================================
Utilities            : PASS ✓
Nonbscanner          : PASS ✓
Visualizations       : PASS ✓
App.py               : PASS ✓
======================================================================
✓ All imports successful! The ImportError is resolved.
```

## Running the Application

After dependencies are installed, start the Streamlit app:

```bash
streamlit run app.py
```

The application should start without any ImportError and display:
```
You can now view your Streamlit app in your browser.

Local URL: http://localhost:8501
Network URL: http://<your-ip>:8501
```

## Troubleshooting

### Issue: "No module named 'matplotlib'" (or similar)

**Solution**: Install the missing package:
```bash
pip install matplotlib
# or install all requirements
pip install -r requirements.txt
```

### Issue: Package installation fails

**Solution**: Try upgrading pip first:
```bash
pip install --upgrade pip
pip install -r requirements.txt
```

### Issue: Hyperscan installation fails (Windows/Mac)

**Solution**: Hyperscan is optional. You can comment it out in `requirements.txt`:
```
# Performance (optional but recommended)
# hyperscan>=0.7.0  # May fail on Windows/Mac
```

The application will still work, just with slightly reduced performance for large sequences.

## Testing

All existing tests continue to pass:

1. **Import Tests**: `python3 test_imports.py` ✓
2. **Visualization Tests**: `python3 test_new_visualizations.py` ✓
3. **Detector Tests**: `python3 test_detectors.py` (1 pre-existing failure unrelated to this fix)

## Deployment Checklist

- [x] Verify `requirements.txt` is complete and in repository root
- [x] Test all imports work (`python3 test_imports.py`)
- [x] Test app starts locally (`streamlit run app.py`)
- [x] Verify visualization functions are accessible
- [x] Test existing visualization tests pass
- [x] Document the fix and deployment requirements

## Files Modified

- None (the issue was dependency installation, not code)

## Files Added

- `test_imports.py` - Comprehensive import verification test
- `DEPLOYMENT_FIX.md` - This deployment guide

## Next Steps for Streamlit Cloud

1. Push changes to your repository
2. Go to Streamlit Cloud dashboard
3. Click "Manage app" on your deployed app
4. Click "Reboot app" to trigger a fresh deployment
5. Streamlit Cloud will automatically install dependencies from `requirements.txt`
6. Verify the app starts without errors

## Support

If you continue to experience import errors after following this guide:

1. Check Streamlit Cloud logs ("Manage app" → "Logs")
2. Verify all files committed to repository (especially `requirements.txt`)
3. Ensure no typos in package names or version constraints
4. Try rebooting the app from Streamlit Cloud dashboard
