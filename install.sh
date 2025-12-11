#!/bin/bash
# Installation script for NonBDNAFinder
# This script installs core requirements and attempts to install optional dependencies

set -e  # Exit on error for core requirements

echo "=============================================="
echo "NonBDNAFinder Installation"
echo "=============================================="
echo ""

# Install core requirements
echo "📦 Installing core requirements..."
pip install -r requirements.txt
echo "✅ Core requirements installed successfully"
echo ""

# Try to install optional requirements (don't fail if this doesn't work)
echo "📦 Attempting to install optional performance dependencies..."
echo "(Hyperscan requires system libraries: build-essential, cmake, libboost-all-dev, ragel)"
echo ""

if pip install -r requirements-optional.txt 2>/dev/null; then
    echo "✅ Optional dependencies installed successfully"
    echo "   Hyperscan is available for high-performance pattern matching"
else
    echo "⚠️  Optional dependencies could not be installed"
    echo "   This is OK - the application will use fallback implementations"
    echo "   Performance may be slightly reduced but all features will work"
fi

echo ""
echo "=============================================="
echo "Installation Complete!"
echo "=============================================="
echo ""
echo "To run the application:"
echo "  streamlit run app.py"
echo ""
echo "To test the installation:"
echo "  python test_app_imports.py     # Quick import test"
echo "  python test_deployment.py      # Comprehensive deployment test (recommended)"
echo ""
