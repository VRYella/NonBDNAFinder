#!/usr/bin/env python3
"""
Quick verification script for NonBDNA Finder Standalone
Run this to verify the installation is working correctly
"""

import sys
import os

def check_installation():
    """Quick check of standalone installation"""
    
    print("🧬 NonBDNA Finder Standalone - Installation Check")
    print("=" * 50)
    
    # Check Python version
    print("🐍 Python version check...")
    if sys.version_info >= (3, 7):
        print(f"   ✅ Python {sys.version_info.major}.{sys.version_info.minor} (compatible)")
    else:
        print(f"   ❌ Python {sys.version_info.major}.{sys.version_info.minor} (requires 3.7+)")
        return False
    
    # Check core module
    print("\n📦 Core module check...")
    try:
        from nbdfinder_core import StandaloneNonBDNAFinder, analyze_fasta_file
        print("   ✅ nbdfinder_core module loaded")
    except ImportError as e:
        print(f"   ❌ Failed to import core module: {e}")
        return False
    
    # Check dependencies
    print("\n📚 Dependencies check...")
    deps = {
        'pandas': 'Data manipulation',
        'numpy': 'Numerical computing',
        'Bio': 'Bioinformatics (biopython)',
        'matplotlib': 'Basic plotting',
        'openpyxl': 'Excel support',
        'plotly': 'Interactive plots',
    }
    
    failed_deps = []
    for dep, desc in deps.items():
        try:
            if dep == 'Bio':
                from Bio import SeqIO
            else:
                __import__(dep)
            print(f"   ✅ {dep} - {desc}")
        except ImportError:
            print(f"   ❌ {dep} - {desc} (MISSING)")
            failed_deps.append(dep)
    
    if failed_deps:
        print(f"\n   ⚠️ Missing dependencies: {', '.join(failed_deps)}")
        print("   💡 Run: pip install -r requirements.txt")
        return False
    
    # Check optional dependencies
    print("\n⚡ Optional enhancements...")
    try:
        import hyperscan
        print("   ✅ hyperscan - Performance optimization")
    except ImportError:
        print("   ⚠️ hyperscan - Not available (normal on Windows)")
    
    try:
        import ipywidgets
        print("   ✅ ipywidgets - Jupyter notebook widgets")
    except ImportError:
        print("   ⚠️ ipywidgets - Missing (needed for notebook interface)")
    
    # Quick functionality test
    print("\n🧪 Functionality test...")
    try:
        test_seq = ">test\nGGGTTTTGGGTTTTGGGTTTTGGG"
        df, zip_data = analyze_fasta_file(test_seq)
        if len(df) > 0:
            print(f"   ✅ Analysis working (detected {len(df)} motifs)")
        else:
            print("   ⚠️ Analysis working but no motifs detected")
    except Exception as e:
        print(f"   ❌ Analysis failed: {e}")
        return False
    
    # Check file structure
    print("\n📁 File structure check...")
    files = {
        'nbdfinder_core.py': 'Core analysis engine',
        'NonBDNA_Finder_Notebook.ipynb': 'Jupyter notebook interface',
        'launcher.py': 'Command-line interface',
        'requirements.txt': 'Python dependencies',
        'README.md': 'Documentation',
        'test_standalone.py': 'Test suite'
    }
    
    missing_files = []
    for filename, desc in files.items():
        if os.path.exists(filename):
            print(f"   ✅ {filename} - {desc}")
        else:
            print(f"   ❌ {filename} - {desc} (MISSING)")
            missing_files.append(filename)
    
    if missing_files:
        print(f"\n   ⚠️ Missing files: {', '.join(missing_files)}")
        return False
    
    print("\n🎉 Installation verification completed!")
    print("\n📋 SUMMARY:")
    print("   ✅ All core dependencies available")
    print("   ✅ Analysis engine functional")  
    print("   ✅ All required files present")
    print("\n🚀 Ready to use! Try:")
    print("   • Jupyter notebook: jupyter notebook NonBDNA_Finder_Notebook.ipynb")
    print("   • Command line: python launcher.py your_file.fasta")
    print("   • Full test: python test_standalone.py")
    
    return True

if __name__ == "__main__":
    success = check_installation()
    if not success:
        print("\n❌ Installation issues detected. Please address the problems above.")
        sys.exit(1)
    else:
        print("\n✅ Installation verified successfully!")
        sys.exit(0)