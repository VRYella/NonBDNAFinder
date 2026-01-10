#!/usr/bin/env python3
"""
Modular Architecture Verification Script
=========================================

Verifies the modular architecture implementation by:
1. Testing example modules work correctly
2. Checking module structure and conventions
3. Verifying documentation completeness
4. Reporting on migration progress

Usage:
    python verify_modules.py [--verbose]
"""

import os
import sys
import importlib.util
import ast
from pathlib import Path
from typing import Tuple, List, Dict

# ANSI color codes for terminal output
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
BLUE = '\033[94m'
RESET = '\033[0m'


def check_mark(passed: bool) -> str:
    """Return colored checkmark or X"""
    return f"{GREEN}✓{RESET}" if passed else f"{RED}✗{RESET}"


def verify_directory_structure() -> Tuple[bool, List[str]]:
    """Verify all required directories exist"""
    print(f"\n{BLUE}=== Verifying Directory Structure ==={RESET}")
    
    required_dirs = [
        'ui',
        'engine',
        'engine/detectors',
        'utils',
        'utils/plotting',
    ]
    
    results = []
    all_passed = True
    
    for dir_path in required_dirs:
        exists = os.path.isdir(dir_path)
        has_init = os.path.exists(os.path.join(dir_path, '__init__.py'))
        
        if exists and has_init:
            results.append(f"{check_mark(True)} {dir_path}/ exists with __init__.py")
        elif exists:
            results.append(f"{check_mark(False)} {dir_path}/ exists but missing __init__.py")
            all_passed = False
        else:
            results.append(f"{check_mark(False)} {dir_path}/ does not exist")
            all_passed = False
    
    for result in results:
        print(f"  {result}")
    
    return all_passed, results


def verify_module_file(filepath: str) -> Tuple[bool, Dict[str, any]]:
    """Verify a Python module file meets quality standards"""
    if not os.path.exists(filepath):
        return False, {'error': 'File does not exist'}
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Parse the module
    try:
        tree = ast.parse(content)
    except SyntaxError as e:
        return False, {'error': f'Syntax error: {e}'}
    
    # Check module docstring
    has_module_docstring = (
        isinstance(tree.body[0], ast.Expr) and
        isinstance(tree.body[0].value, ast.Constant) and
        isinstance(tree.body[0].value.value, str)
    )
    
    # Count functions and classes
    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    classes = [node for node in tree.body if isinstance(node, ast.ClassDef)]
    
    # Check function docstrings
    functions_with_docstrings = sum(
        1 for func in functions
        if (func.body and 
            isinstance(func.body[0], ast.Expr) and
            isinstance(func.body[0].value, ast.Constant) and
            isinstance(func.body[0].value.value, str))
    )
    
    # Check type hints
    functions_with_hints = sum(
        1 for func in functions
        if func.returns is not None or any(arg.annotation for arg in func.args.args)
    )
    
    # Count lines
    lines = content.split('\n')
    total_lines = len(lines)
    code_lines = len([line for line in lines if line.strip() and not line.strip().startswith('#')])
    
    results = {
        'has_module_docstring': has_module_docstring,
        'total_lines': total_lines,
        'code_lines': code_lines,
        'num_functions': len(functions),
        'num_classes': len(classes),
        'functions_with_docstrings': functions_with_docstrings,
        'functions_with_hints': functions_with_hints,
    }
    
    # Module passes if:
    # - Has module docstring
    # - Under 400 lines (target ~200, max ~300, but flexible)
    # - Functions have docstrings
    # - Functions have type hints
    passed = (
        has_module_docstring and
        total_lines < 400 and
        (len(functions) == 0 or functions_with_docstrings == len(functions)) and
        (len(functions) == 0 or functions_with_hints >= len(functions) * 0.8)
    )
    
    return passed, results


def verify_example_modules() -> Tuple[bool, List[str]]:
    """Verify example module implementations"""
    print(f"\n{BLUE}=== Verifying Example Modules ==={RESET}")
    
    modules = [
        ('engine/sequence_ops.py', 'Sequence operations'),
        ('utils/fasta.py', 'FASTA parsing'),
        ('utils/validation.py', 'Validation functions'),
    ]
    
    results = []
    all_passed = True
    
    for filepath, description in modules:
        passed, info = verify_module_file(filepath)
        
        if passed:
            status = f"{check_mark(True)} {filepath}"
            details = f"    {info['total_lines']} lines, {info['num_functions']} functions"
        else:
            status = f"{check_mark(False)} {filepath}"
            if 'error' in info:
                details = f"    Error: {info['error']}"
            else:
                issues = []
                if not info['has_module_docstring']:
                    issues.append("missing module docstring")
                if info['total_lines'] >= 400:
                    issues.append(f"too long ({info['total_lines']} lines)")
                if info['num_functions'] > 0 and info['functions_with_docstrings'] < info['num_functions']:
                    issues.append(f"missing function docstrings")
                if info['num_functions'] > 0 and info['functions_with_hints'] < info['num_functions'] * 0.8:
                    issues.append(f"missing type hints")
                details = f"    Issues: {', '.join(issues)}"
            all_passed = False
        
        results.append(status)
        print(f"  {status}")
        print(f"  {details}")
    
    return all_passed, results


def verify_documentation() -> Tuple[bool, List[str]]:
    """Verify documentation files exist and are complete"""
    print(f"\n{BLUE}=== Verifying Documentation ==={RESET}")
    
    docs = [
        ('MODULAR_ARCHITECTURE_GUIDE.md', 'Architecture specification'),
        ('IMPLEMENTATION_SUMMARY.md', 'Implementation summary'),
        ('MODULAR_QUICKSTART.md', 'Quick start guide'),
        ('migrate_to_modules.py', 'Migration script'),
    ]
    
    results = []
    all_passed = True
    
    for filepath, description in docs:
        exists = os.path.exists(filepath)
        
        if exists:
            # Check file size (should not be empty)
            size = os.path.getsize(filepath)
            if size > 100:  # At least 100 bytes
                results.append(f"{check_mark(True)} {filepath} ({size:,} bytes)")
            else:
                results.append(f"{check_mark(False)} {filepath} (too small: {size} bytes)")
                all_passed = False
        else:
            results.append(f"{check_mark(False)} {filepath} (not found)")
            all_passed = False
    
    for result in results:
        print(f"  {result}")
    
    return all_passed, results


def test_module_imports() -> Tuple[bool, List[str]]:
    """Test that example modules can be imported"""
    print(f"\n{BLUE}=== Testing Module Imports ==={RESET}")
    
    modules_to_test = [
        'engine.sequence_ops',
        'utils.fasta',
        'utils.validation',
    ]
    
    results = []
    all_passed = True
    
    for module_name in modules_to_test:
        try:
            # Import the module
            module = importlib.import_module(module_name)
            
            # Check for expected functions
            functions = [name for name in dir(module) if not name.startswith('_')]
            
            if functions:
                results.append(
                    f"{check_mark(True)} {module_name} imports successfully "
                    f"({len(functions)} public functions)"
                )
            else:
                results.append(
                    f"{check_mark(False)} {module_name} imports but has no public functions"
                )
                all_passed = False
                
        except ImportError as e:
            results.append(f"{check_mark(False)} {module_name} failed to import: {e}")
            all_passed = False
        except Exception as e:
            results.append(f"{check_mark(False)} {module_name} error: {e}")
            all_passed = False
    
    for result in results:
        print(f"  {result}")
    
    return all_passed, results


def calculate_progress() -> Tuple[int, Dict[str, any]]:
    """Calculate overall migration progress"""
    print(f"\n{BLUE}=== Migration Progress ==={RESET}")
    
    # Count what's done
    done = {
        'directories': 5,  # ui, engine, engine/detectors, utils, utils/plotting
        'example_modules': 3,  # sequence_ops, fasta, validation
        'documentation': 3,  # 3 MD files
        'tooling': 1,  # migration script
    }
    
    # Estimate what's remaining (from IMPLEMENTATION_SUMMARY.md)
    remaining = {
        'engine_modules': 5,  # detection, scoring, merging, chunking, patterns
        'detector_classes': 10,  # 9 detectors + base class
        'utility_modules': 10,  # registry, caching, state, export + 6 plotting
        'ui_modules': 6,  # layout, formatting, metrics, progress, inputs, downloads
    }
    
    total_items = sum(done.values()) + sum(remaining.values())
    completed_items = sum(done.values())
    progress_percent = (completed_items / total_items) * 100
    
    print(f"  {GREEN}Completed:{RESET}")
    print(f"    • Directory structure: {done['directories']} packages")
    print(f"    • Example modules: {done['example_modules']} modules")
    print(f"    • Documentation files: {done['documentation']} documents")
    print(f"    • Migration tooling: {done['tooling']} script")
    print(f"    {GREEN}Total completed: {completed_items} items{RESET}")
    
    print(f"\n  {YELLOW}Remaining:{RESET}")
    print(f"    • Engine modules: {remaining['engine_modules']} modules")
    print(f"    • Detector classes: {remaining['detector_classes']} classes")
    print(f"    • Utility modules: {remaining['utility_modules']} modules")
    print(f"    • UI modules: {remaining['ui_modules']} modules")
    print(f"    {YELLOW}Total remaining: {sum(remaining.values())} items{RESET}")
    
    print(f"\n  {BLUE}Overall Progress: {progress_percent:.1f}% ({completed_items}/{total_items}){RESET}")
    
    return int(progress_percent), {'done': done, 'remaining': remaining}


def main():
    """Run all verification checks"""
    print(f"{BLUE}{'='*60}{RESET}")
    print(f"{BLUE}Modular Architecture Verification{RESET}")
    print(f"{BLUE}{'='*60}{RESET}")
    
    results = {}
    
    # Run verification checks
    results['structure'] = verify_directory_structure()
    results['modules'] = verify_example_modules()
    results['docs'] = verify_documentation()
    results['imports'] = test_module_imports()
    
    # Calculate progress
    progress, details = calculate_progress()
    
    # Summary
    print(f"\n{BLUE}{'='*60}{RESET}")
    print(f"{BLUE}Verification Summary{RESET}")
    print(f"{BLUE}{'='*60}{RESET}")
    
    all_checks = [
        ('Directory Structure', results['structure'][0]),
        ('Example Modules', results['modules'][0]),
        ('Documentation', results['docs'][0]),
        ('Module Imports', results['imports'][0]),
    ]
    
    for check_name, passed in all_checks:
        print(f"  {check_mark(passed)} {check_name}")
    
    all_passed = all(passed for _, passed in all_checks)
    
    print(f"\n{BLUE}Overall Status:{RESET}")
    if all_passed:
        print(f"  {GREEN}✓ All verification checks passed!{RESET}")
        print(f"  {GREEN}✓ Foundation is complete and ready for full migration{RESET}")
        exit_code = 0
    else:
        print(f"  {YELLOW}⚠ Some checks failed{RESET}")
        print(f"  {YELLOW}⚠ Review output above for details{RESET}")
        exit_code = 1
    
    print(f"\n{BLUE}Next Steps:{RESET}")
    print(f"  1. Review IMPLEMENTATION_SUMMARY.md for detailed migration plan")
    print(f"  2. Use migrate_to_modules.py to extract remaining modules")
    print(f"  3. Test continuously after each extraction")
    print(f"  4. Update imports in main files")
    
    print(f"\n{BLUE}{'='*60}{RESET}\n")
    
    return exit_code


if __name__ == '__main__':
    sys.exit(main())
