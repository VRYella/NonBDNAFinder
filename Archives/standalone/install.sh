#!/bin/bash
# Installation script for NonBDNA Finder Standalone

echo "🧬 NonBDNA Finder Standalone - Installation Script"
echo "=================================================="

# Check Python version
echo "🐍 Checking Python version..."
python_version=$(python3 --version 2>&1 | cut -d' ' -f2 | cut -d'.' -f1,2)
required_version="3.7"

if python3 -c "import sys; exit(0 if sys.version_info >= (3,7) else 1)"; then
    echo "   ✅ Python $python_version detected (compatible)"
else
    echo "   ❌ Python 3.7+ required, found $python_version"
    exit 1
fi

# Install core dependencies
echo ""
echo "📦 Installing Python dependencies..."
pip3 install -r requirements.txt

if [ $? -eq 0 ]; then
    echo "   ✅ Dependencies installed successfully"
else
    echo "   ❌ Failed to install dependencies"
    exit 1
fi

# Try to install Jupyter if not present
echo ""
echo "📓 Checking Jupyter installation..."
if command -v jupyter &> /dev/null; then
    echo "   ✅ Jupyter already installed"
else
    echo "   📦 Installing Jupyter..."
    pip3 install jupyter
    if [ $? -eq 0 ]; then
        echo "   ✅ Jupyter installed successfully"
    else
        echo "   ⚠️ Failed to install Jupyter, but core functionality should work"
    fi
fi

# Try to install optional hyperscan
echo ""
echo "⚡ Installing optional performance enhancement..."
pip3 install hyperscan &> /dev/null
if [ $? -eq 0 ]; then
    echo "   ✅ Hyperscan installed for optimized performance"
else
    echo "   ⚠️ Hyperscan not available (normal for Windows), using standard regex"
fi

# Run tests
echo ""
echo "🧪 Running tests..."
python3 test_standalone.py
if [ $? -eq 0 ]; then
    echo ""
    echo "🎉 Installation completed successfully!"
    echo ""
    echo "Next steps:"
    echo "1. Launch Jupyter: jupyter notebook"
    echo "2. Open: NonBDNA_Finder_Notebook.ipynb"
    echo "3. Upload your FASTA file and analyze!"
    echo ""
    echo "📚 See README.md for detailed usage instructions"
else
    echo ""
    echo "❌ Installation tests failed. Please check error messages above."
    exit 1
fi