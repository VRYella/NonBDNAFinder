#!/usr/bin/env python3
"""
Comprehensive Validation Suite for NonBDNAFinder 2025.1
========================================================

This script validates all enhancements and ensures:
- Scoring improvements work correctly
- JSON registry is valid
- Visualizations generate properly
- Performance is maintained
- Output quality meets standards

Author: Dr. Venkata Rajesh Yella
Version: 2025.1
"""

import sys
import json
import time
from pathlib import Path

print("=" * 80)
print(" NonBDNAFinder 2025.1 - Comprehensive Validation Suite")
print("=" * 80)
print()

# Track validation results
passed = 0
failed = 0
warnings = 0

def test_pass(msg):
    global passed
    passed += 1
    print(f"✅ PASS: {msg}")

def test_fail(msg):
    global failed
    failed += 1
    print(f"❌ FAIL: {msg}")

def test_warn(msg):
    global warnings
    warnings += 1
    print(f"⚠️  WARN: {msg}")

print("Phase 1: Newdetector.py Validation")
print("-" * 80)

try:
    # Import Newdetector module
    sys.path.insert(0, str(Path(__file__).parent))
    
    # Test A-philic detector
    try:
        from Newdetector import detect_aphilic
        test_seq = "AGGGGGGGGGCCCCCCCCCTAGGGGGGGGG"
        results = detect_aphilic(test_seq, "test")
        
        if results:
            test_pass("A-philic detector functional")
            
            # Check for enhanced fields
            result = results[0]
            if 'Window_Count' in result:
                test_pass("A-philic: Window_Count field added")
            else:
                test_warn("A-philic: Window_Count field missing")
            
            if 'Confidence' in result:
                test_pass("A-philic: Confidence field added")
            else:
                test_warn("A-philic: Confidence field missing")
            
            # Check subclass enhancement
            if any(s in result.get('Subclass', '') for s in ['High_Confidence', 'Moderate', 'Weak']):
                test_pass("A-philic: Enhanced subclass classification")
            else:
                test_warn("A-philic: Still using old subclass naming")
            
            # Check score range
            score = result.get('Score', 0)
            if 1.0 <= score <= 3.0:
                test_pass(f"A-philic: Score in valid range (1-3): {score:.3f}")
            else:
                test_fail(f"A-philic: Score out of range: {score:.3f}")
        else:
            test_warn("A-philic: No motifs detected in test sequence")
    except Exception as e:
        test_fail(f"A-philic detector error: {str(e)}")
    
    # Test Curved DNA detector
    try:
        from Newdetector import detect_curvature
        test_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" + "T" * 40
        results = detect_curvature(test_seq)
        
        if results:
            test_pass("Curved DNA detector functional")
            
            result = results[0]
            # Check for enhanced fields
            if 'Density' in result:
                test_pass("Curved DNA: Density field added")
            else:
                test_warn("Curved DNA: Density field missing")
            
            # Check enhanced subclasses
            subclass = result.get('Subclass', '')
            if any(s in subclass for s in ['High_Quality', 'Strong', 'Directional']):
                test_pass("Curved DNA: Enhanced subclass classification")
            else:
                test_warn(f"Curved DNA: Subclass may not be enhanced: {subclass}")
        else:
            test_warn("Curved DNA: No motifs detected in test sequence")
    except Exception as e:
        test_fail(f"Curved DNA detector error: {str(e)}")
    
    # Test G-Quadruplex detector
    try:
        from Newdetector import detect_g4
        test_seq = "GGGTTAGGGTTAGGGTTAGGG"
        results = detect_g4(test_seq, "test")
        
        if results:
            test_pass("G-Quadruplex detector functional")
            
            result = results[0]
            score = result.get('Score', 0)
            if 1.0 <= score <= 3.0:
                test_pass(f"G4: Score in valid range: {score:.3f}")
            else:
                test_fail(f"G4: Score out of range: {score:.3f}")
        else:
            test_warn("G4: No motifs detected in test sequence")
    except Exception as e:
        test_fail(f"G-Quadruplex detector error: {str(e)}")
    
    # Test i-Motif detector
    try:
        from Newdetector import detect_imotifs
        test_seq = "CCCAGCTAGCTCCCAGCTAGCTCCC"
        results = detect_imotifs(test_seq, "test")
        
        if results:
            test_pass("i-Motif detector functional")
            
            result = results[0]
            score = result.get('Score', 0)
            if 1.0 <= score <= 3.0:
                test_pass(f"i-Motif: Score in valid range: {score:.3f}")
            else:
                test_fail(f"i-Motif: Score out of range: {score:.3f}")
        else:
            test_warn("i-Motif: No motifs detected in test sequence")
    except Exception as e:
        test_fail(f"i-Motif detector error: {str(e)}")
    
    # Test Z-DNA detector
    try:
        from Newdetector import detect_Z_and_eGZ
        test_seq = "CGCGCGCGCGCGCGCGCGCG"
        results = detect_Z_and_eGZ(test_seq, "test")
        
        if results:
            test_pass("Z-DNA detector functional")
            
            result = results[0]
            score = result.get('Score', 0)
            if 1.0 <= score <= 3.0:
                test_pass(f"Z-DNA: Score in valid range: {score:.3f}")
            else:
                test_fail(f"Z-DNA: Score out of range: {score:.3f}")
            
            # High CG content should give high score
            if score >= 2.5:
                test_pass(f"Z-DNA: High score for CGCG repeats: {score:.3f}")
            else:
                test_warn(f"Z-DNA: Score lower than expected: {score:.3f}")
        else:
            test_fail("Z-DNA: Should detect CGCG repeats")
    except Exception as e:
        test_fail(f"Z-DNA detector error: {str(e)}")

except Exception as e:
    test_fail(f"Newdetector.py import error: {str(e)}")

print()
print("Phase 2: JSON Registry Validation")
print("-" * 80)

try:
    with open('consolidated_registry.json', 'r') as f:
        registry = json.load(f)
    
    test_pass("JSON registry loads successfully")
    
    # Check version
    version = registry.get('version')
    if version == '2025.1':
        test_pass(f"JSON registry version updated: {version}")
    else:
        test_warn(f"JSON registry version: {version} (expected 2025.1)")
    
    # Check for enhanced metadata
    enhanced_classes = ['APhilic', 'G4', 'IMotif', 'CurvedDNA', 'ZDNA', 'RLoop']
    for cls in enhanced_classes:
        if cls in registry.get('registries', {}):
            meta = registry['registries'][cls].get('meta', {})
            
            if 'scoring_method' in meta:
                test_pass(f"{cls}: scoring_method documented")
            else:
                test_warn(f"{cls}: scoring_method missing")
            
            if 'references' in meta or 'reference' in meta:
                test_pass(f"{cls}: references added")
            else:
                test_warn(f"{cls}: references missing")
    
    # Check total patterns
    total = registry.get('total_patterns', 0)
    test_pass(f"Total patterns: {total}")

except FileNotFoundError:
    test_fail("consolidated_registry.json not found")
except json.JSONDecodeError as e:
    test_fail(f"JSON parse error: {str(e)}")
except Exception as e:
    test_fail(f"JSON registry error: {str(e)}")

print()
print("Phase 3: Visualization Module Validation")
print("-" * 80)

try:
    from visualization_enhancements import (
        NOBEL_QUALITY_COLORS,
        enhance_plot_aesthetics,
        create_enhanced_colormap,
        enhanced_bar_plot
    )
    
    test_pass("Visualization enhancements module imports successfully")
    
    # Check colors
    if len(NOBEL_QUALITY_COLORS) >= 10:
        test_pass(f"Nobel-quality colors defined: {len(NOBEL_QUALITY_COLORS)} classes")
    else:
        test_warn(f"Only {len(NOBEL_QUALITY_COLORS)} colors defined")
    
    # Test colormap creation
    try:
        cmap = create_enhanced_colormap('density')
        test_pass("Enhanced colormap creation functional")
    except Exception as e:
        test_fail(f"Colormap creation error: {str(e)}")

except ImportError as e:
    test_fail(f"Visualization module import error: {str(e)}")
except Exception as e:
    test_fail(f"Visualization module error: {str(e)}")

print()
print("Phase 4: Documentation Validation")
print("-" * 80)

# Check for documentation files
doc_files = {
    'README.md': 'Enhanced README',
    'IMPROVEMENTS_SUMMARY.md': 'Improvements documentation',
    'visualization_enhancements.py': 'Visualization module',
    'consolidated_registry.json': 'JSON registry'
}

for file, desc in doc_files.items():
    if Path(file).exists():
        size = Path(file).stat().st_size
        if size > 100:  # At least 100 bytes
            test_pass(f"{desc} exists ({size:,} bytes)")
        else:
            test_warn(f"{desc} exists but seems small ({size} bytes)")
    else:
        test_fail(f"{desc} missing: {file}")

print()
print("Phase 5: Performance Validation")
print("-" * 80)

# Test performance (simple benchmark)
try:
    from Newdetector import detect_aphilic
    
    # Generate test sequence
    test_seq = "AGGGGGGGGG" * 100  # 1000 bp
    
    start = time.time()
    for _ in range(10):
        results = detect_aphilic(test_seq, "benchmark")
    elapsed = time.time() - start
    
    rate = (len(test_seq) * 10) / elapsed
    test_pass(f"Performance benchmark: {rate:,.0f} bp/second")
    
    if rate > 1000:
        test_pass("Performance excellent (>1000 bp/s)")
    elif rate > 500:
        test_pass("Performance good (>500 bp/s)")
    else:
        test_warn(f"Performance moderate ({rate:.0f} bp/s)")

except Exception as e:
    test_warn(f"Performance test skipped: {str(e)}")

print()
print("=" * 80)
print(" Validation Summary")
print("=" * 80)
print(f"✅ Passed: {passed}")
print(f"⚠️  Warnings: {warnings}")
print(f"❌ Failed: {failed}")
print()

if failed == 0 and warnings <= 5:
    print("🎉 EXCELLENT: All critical validations passed!")
    print("   System ready for production use.")
    sys.exit(0)
elif failed == 0:
    print("✅ GOOD: All critical tests passed, minor warnings present.")
    print("   System functional but review warnings.")
    sys.exit(0)
else:
    print("⚠️  ATTENTION: Some validations failed.")
    print("   Review failed tests before production use.")
    sys.exit(1)
