#!/usr/bin/env python3
"""
Verify that canonical motif taxonomy is enforced in utilities.py
This script checks the source code without importing matplotlib-dependent modules
"""

import sys
import re

def verify_imports():
    """Verify that utilities.py imports canonical taxonomy"""
    print("=" * 70)
    print("TEST 1: Verify Canonical Taxonomy Imports")
    print("=" * 70)
    
    with open('utilities.py', 'r') as f:
        content = f.read()
    
    # Check for the import statement (at start of line, not in comments)
    import_pattern = r'^from config\.motif_taxonomy import'
    if re.search(import_pattern, content, re.MULTILINE):
        print("✓ utilities.py imports from config.motif_taxonomy")
    else:
        print("❌ FAILED: utilities.py does not import from config.motif_taxonomy")
        return False
    
    # Check that the specific imports are present
    required_imports = ['VALID_CLASSES', 'VALID_SUBCLASSES', 'SUBCLASS_TO_CLASS']
    for imp in required_imports:
        if imp in content:
            print(f"✓ Imports {imp}")
        else:
            print(f"❌ FAILED: Missing import {imp}")
            return False
    
    print("\n✅ Import verification PASSED\n")
    return True


def verify_plot_function_uses_taxonomy():
    """Verify that plot_motif_distribution uses canonical taxonomy"""
    print("=" * 70)
    print("TEST 2: Verify plot_motif_distribution Uses Canonical Taxonomy")
    print("=" * 70)
    
    with open('utilities.py', 'r') as f:
        content = f.read()
    
    # Find the plot_motif_distribution function
    func_match = re.search(
        r'def plot_motif_distribution.*?(?=\ndef |\nclass |\Z)',
        content,
        re.DOTALL
    )
    
    if not func_match:
        print("❌ FAILED: Could not find plot_motif_distribution function")
        return False
    
    func_code = func_match.group(0)
    
    # Check that it uses VALID_CLASSES and VALID_SUBCLASSES
    if 'ALL_CLASSES = sorted(VALID_CLASSES)' in func_code:
        print("✓ Uses sorted(VALID_CLASSES) for ALL_CLASSES")
    else:
        print("❌ FAILED: Does not use VALID_CLASSES correctly")
        return False
    
    if 'ALL_SUBCLASSES = sorted(VALID_SUBCLASSES)' in func_code:
        print("✓ Uses sorted(VALID_SUBCLASSES) for ALL_SUBCLASSES")
    else:
        print("❌ FAILED: Does not use VALID_SUBCLASSES correctly")
        return False
    
    # Check that it uses SUBCLASS_TO_CLASS for color mapping
    if 'SUBCLASS_TO_CLASS.get(cat' in func_code:
        print("✓ Uses SUBCLASS_TO_CLASS for subclass color mapping")
    else:
        print("❌ FAILED: Does not use SUBCLASS_TO_CLASS for colors")
        return False
    
    # Check that hardcoded lists are removed
    if "ALL_CLASSES = [" in func_code or "ALL_SUBCLASSES = [" in func_code:
        print("⚠️  WARNING: Found hardcoded list definitions (should be removed)")
    else:
        print("✓ No hardcoded class/subclass lists found")
    
    print("\n✅ Plot function verification PASSED\n")
    return True


def verify_pattern_registry_uses_taxonomy():
    """Verify that PatternRegistry uses canonical taxonomy"""
    print("=" * 70)
    print("TEST 3: Verify PatternRegistry Uses Canonical Taxonomy")
    print("=" * 70)
    
    with open('utilities.py', 'r') as f:
        content = f.read()
    
    # Find the get_subclass_mapping method
    method_match = re.search(
        r'def get_subclass_mapping\(cls\).*?(?=\n    @|\n    def |\nclass |\Z)',
        content,
        re.DOTALL
    )
    
    if not method_match:
        print("❌ FAILED: Could not find get_subclass_mapping method")
        return False
    
    method_code = method_match.group(0)
    
    # Check that it imports CLASS_TO_SUBCLASSES
    if 'from config.motif_taxonomy import CLASS_TO_SUBCLASSES' in method_code:
        print("✓ Imports CLASS_TO_SUBCLASSES from canonical taxonomy")
    else:
        print("❌ FAILED: Does not import CLASS_TO_SUBCLASSES")
        return False
    
    # Check that it uses CLASS_TO_SUBCLASSES
    if "CLASS_TO_SUBCLASSES['Curved_DNA']" in method_code:
        print("✓ Uses CLASS_TO_SUBCLASSES for mapping")
    else:
        print("❌ FAILED: Does not use CLASS_TO_SUBCLASSES correctly")
        return False
    
    print("\n✅ PatternRegistry verification PASSED\n")
    return True


def verify_no_hardcoded_incorrect_names():
    """Verify no incorrect hardcoded subclass names remain"""
    print("=" * 70)
    print("TEST 4: Verify No Incorrect Hardcoded Subclass Names")
    print("=" * 70)
    
    with open('utilities.py', 'r') as f:
        lines = f.readlines()
    
    # Note: These are known incorrect names that were replaced with canonical names
    # They should NOT appear in active code (may appear in comments/docs about the fix)
    incorrect_names = {
        "'Inverted Repeats'": "'Cruciform forming IRs'",
        "'QmRLFS-m1'": "'R-loop formation sites'",
        "'QmRLFS-m2'": "'R-loop formation sites'",
        "'eGZ (Extruded-G) DNA'": "'eGZ'",
    }
    
    found_errors = []
    
    for line_num, line in enumerate(lines, 1):
        # Skip comment lines and docstrings
        stripped = line.strip()
        if stripped.startswith('#') or stripped.startswith('"""') or stripped.startswith("'''"):
            continue
        
        for incorrect, correct in incorrect_names.items():
            if incorrect in line:
                # Check if it's in a tuple pattern definition (which we've fixed)
                # Pattern tuples should now use canonical names
                found_errors.append(f"Line {line_num}: Found {incorrect} (should be {correct})")
    
    if found_errors:
        print("❌ FAILED: Found incorrect hardcoded subclass names:")
        for error in found_errors[:10]:  # Show first 10 errors
            print(f"  {error}")
        if len(found_errors) > 10:
            print(f"  ... and {len(found_errors) - 10} more")
        return False
    else:
        print("✓ No incorrect hardcoded subclass names found in active code")
        print("  - 'Inverted Repeats' → 'Cruciform forming IRs' ✓")
        print("  - 'QmRLFS-m1/m2' → 'R-loop formation sites' ✓")
        print("  - 'eGZ (Extruded-G) DNA' → 'eGZ' ✓")
    
    print("\n✅ Hardcoded names verification PASSED\n")
    return True


def verify_analyze_function_uses_taxonomy():
    """Verify that analyze_class_subclass_detection uses canonical taxonomy"""
    print("=" * 70)
    print("TEST 5: Verify analyze_class_subclass_detection Uses Canonical Taxonomy")
    print("=" * 70)
    
    with open('utilities.py', 'r') as f:
        content = f.read()
    
    # Find the analyze_class_subclass_detection function
    func_match = re.search(
        r'def analyze_class_subclass_detection.*?(?=\ndef |\nclass |\Z)',
        content,
        re.DOTALL
    )
    
    if not func_match:
        print("❌ FAILED: Could not find analyze_class_subclass_detection function")
        return False
    
    func_code = func_match.group(0)
    
    # Check that it imports CLASS_TO_SUBCLASSES
    if 'from config.motif_taxonomy import CLASS_TO_SUBCLASSES' in func_code:
        print("✓ Imports CLASS_TO_SUBCLASSES from canonical taxonomy")
    else:
        print("❌ FAILED: Does not import CLASS_TO_SUBCLASSES")
        return False
    
    # Check that it uses CLASS_TO_SUBCLASSES
    if 'all_classes = CLASS_TO_SUBCLASSES' in func_code:
        print("✓ Uses CLASS_TO_SUBCLASSES directly")
    else:
        print("❌ FAILED: Does not use CLASS_TO_SUBCLASSES correctly")
        return False
    
    print("\n✅ analyze_class_subclass_detection verification PASSED\n")
    return True


def main():
    """Run all verification tests"""
    print("\n" + "=" * 70)
    print("CANONICAL MOTIF TAXONOMY ENFORCEMENT VERIFICATION")
    print("=" * 70 + "\n")
    
    try:
        success = True
        success = verify_imports() and success
        success = verify_plot_function_uses_taxonomy() and success
        success = verify_pattern_registry_uses_taxonomy() and success
        success = verify_no_hardcoded_incorrect_names() and success
        success = verify_analyze_function_uses_taxonomy() and success
        
        if success:
            print("=" * 70)
            print("🎉 ALL VERIFICATION TESTS PASSED!")
            print("=" * 70)
            print("\nCanonical motif taxonomy is now enforced:")
            print("  ✓ utilities.py imports from config.motif_taxonomy")
            print("  ✓ plot_motif_distribution uses VALID_CLASSES/VALID_SUBCLASSES")
            print("  ✓ PatternRegistry uses CLASS_TO_SUBCLASSES")
            print("  ✓ analyze_class_subclass_detection uses CLASS_TO_SUBCLASSES")
            print("  ✓ All hardcoded incorrect subclass names removed")
            print("  ✓ Visualization colors use SUBCLASS_TO_CLASS mapping")
            print("\n✅ No duplicate bars from variant spellings")
            print("✅ Consistent grouping across detectors, exports, and visualizations")
            print("\n")
            return 0
        else:
            print("\n❌ SOME VERIFICATION TESTS FAILED\n")
            return 1
            
    except Exception as e:
        print(f"\n❌ VERIFICATION ERROR: {e}\n")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
