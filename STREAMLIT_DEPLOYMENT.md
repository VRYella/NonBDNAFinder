# Streamlit Cloud Deployment Guide

## Overview

This document explains how to deploy NonBScanner to Streamlit Cloud and addresses common deployment issues.

## Deployment Files

### Required Files

1. **requirements.txt** - Python package dependencies
   - Contains all Python packages needed by the application
   - Streamlit Cloud automatically installs these during deployment

2. **packages.txt** - System-level dependencies (APT packages)
   - Required for compiling native extensions like hyperscan
   - Contains: build-essential, cmake, libboost-all-dev, ragel
   - Streamlit Cloud installs these before installing Python packages

3. **app.py** - Main Streamlit application entry point

4. **.streamlit/config.toml** - Streamlit configuration
   - Configures upload size limits, CORS, and other settings

## Common Deployment Issues

### ImportError at Line 37 in app.py

**Symptom:** 
```
ImportError: This app has encountered an error
Traceback: File "/mount/src/nonbdnafinder/app.py", line 37
```

**Cause:** 
Missing system-level dependencies required for hyperscan compilation.

**Solution:**
The `packages.txt` file provides the necessary system dependencies. Ensure it's present in the repository root.

### Hyperscan Installation Failures

**Note:** Hyperscan is optional for NonBScanner. The application will work without it, though with reduced performance for some operations.

If hyperscan installation continues to fail on Streamlit Cloud:

1. The code already handles hyperscan as optional (using try-except blocks)
2. The app will fall back to regex-based scanning
3. No code changes are needed - the fallback is automatic

### Testing Deployment Locally

Before deploying to Streamlit Cloud, test locally:

```bash
# Install dependencies
pip install -r requirements.txt

# Run the test script
python test_app_imports.py

# Run the app locally
streamlit run app.py
```

## Streamlit Cloud Configuration

When deploying to Streamlit Cloud:

1. **Connect your GitHub repository**
   - Streamlit Cloud will automatically detect app.py

2. **Advanced Settings (if needed):**
   - Python version: 3.9+ (tested with 3.12)
   - Main file path: `app.py`

3. **Files checked by Streamlit Cloud:**
   - `requirements.txt` - Python dependencies
   - `packages.txt` - System dependencies
   - `.streamlit/config.toml` - App configuration

## Verification

After deployment, verify:

1. ✅ App loads without import errors
2. ✅ All modules import successfully
3. ✅ File upload works (1GB limit configured)
4. ✅ Analysis functions work correctly

## Troubleshooting

### Check Logs

If deployment fails:
1. Click "Manage app" in lower right
2. View logs for detailed error messages
3. Look for import errors or package installation failures

### Common Solutions

1. **Clear cache:** Use "Clear cache" in Streamlit Cloud settings
2. **Rebuild:** Trigger a rebuild by pushing a small change
3. **Check requirements:** Ensure all packages have compatible versions

## Testing

Use the provided test script to verify all imports:

```bash
python test_app_imports.py
```

This tests:
- All Python dependencies are installed
- utilities module imports work
- nonbscanner module imports work
- visualizations module imports work

## Support

For issues related to:
- **Streamlit Cloud:** https://docs.streamlit.io/streamlit-community-cloud
- **NonBScanner:** Open an issue in the GitHub repository
