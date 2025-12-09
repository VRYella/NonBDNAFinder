"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    CODE QUALITY IMPROVEMENTS SCRIPT                           ║
║     Automated improvements for consistency, clarity, and documentation       ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: code_quality_improvements.py
AUTHOR: Dr. Venkata Rajesh Yella  
VERSION: 2024.1
LICENSE: MIT

DESCRIPTION:
    Automated script to improve code quality across the NonBDNAFinder codebase:
    - Add missing docstrings
    - Improve type hints
    - Standardize formatting
    - Add comprehensive comments
    - Optimize imports
    
USAGE:
    python3 code_quality_improvements.py --check    # Check issues only
    python3 code_quality_improvements.py --fix      # Apply fixes
"""

import os
import re
import ast
from typing import List, Dict, Tuple
from pathlib import Path


class CodeQualityChecker:
    """Check and improve code quality across Python files"""
    
    def __init__(self, project_root: str = "."):
        self.project_root = Path(project_root)
        self.issues = []
        
    def check_docstring_coverage(self, filename: str) -> List[str]:
        """Check which functions/classes are missing docstrings"""
        issues = []
        
        with open(filename, 'r') as f:
            try:
                tree = ast.parse(f.read(), filename=filename)
            except SyntaxError:
                return [f"Syntax error in {filename}"]
        
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                # Skip private methods
                if node.name.startswith('_') and not node.name.startswith('__'):
                    continue
                    
                # Check if has docstring
                docstring = ast.get_docstring(node)
                if not docstring:
                    issues.append(f"Missing docstring: function {node.name} at line {node.lineno}")
                elif len(docstring) < 10:
                    issues.append(f"Minimal docstring: function {node.name} at line {node.lineno}")
                    
            elif isinstance(node, ast.ClassDef):
                docstring = ast.get_docstring(node)
                if not docstring:
                    issues.append(f"Missing docstring: class {node.name} at line {node.lineno}")
                    
        return issues
    
    def check_type_hints(self, filename: str) -> List[str]:
        """Check for missing type hints in function signatures"""
        issues = []
        
        with open(filename, 'r') as f:
            try:
                tree = ast.parse(f.read(), filename=filename)
            except SyntaxError:
                return []
        
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                # Skip private methods and __init__
                if node.name.startswith('_'):
                    continue
                
                # Check return type
                if node.returns is None:
                    issues.append(f"Missing return type hint: {node.name} at line {node.lineno}")
                
                # Check parameter types
                for arg in node.args.args:
                    if arg.annotation is None and arg.arg != 'self' and arg.arg != 'cls':
                        issues.append(f"Missing type hint for parameter '{arg.arg}' in {node.name} at line {node.lineno}")
        
        return issues
    
    def check_code_complexity(self, filename: str) -> List[str]:
        """Check for overly complex functions (high cyclomatic complexity)"""
        issues = []
        
        with open(filename, 'r') as f:
            try:
                tree = ast.parse(f.read(), filename=filename)
            except SyntaxError:
                return []
        
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                # Count decision points (if, for, while, except, etc.)
                complexity = self._calculate_complexity(node)
                if complexity > 15:
                    issues.append(f"High complexity ({complexity}): {node.name} at line {node.lineno} - consider refactoring")
                    
        return issues
    
    def _calculate_complexity(self, node: ast.AST) -> int:
        """Calculate cyclomatic complexity of a function"""
        complexity = 1  # Base complexity
        
        for child in ast.walk(node):
            # Each branch point adds 1 to complexity
            if isinstance(child, (ast.If, ast.For, ast.While, ast.ExceptHandler)):
                complexity += 1
            elif isinstance(child, ast.BoolOp):
                # And/Or operators
                complexity += len(child.values) - 1
                
        return complexity
    
    def check_file(self, filename: str) -> Dict[str, List[str]]:
        """Run all checks on a file"""
        print(f"Checking {filename}...")
        
        results = {
            'docstrings': self.check_docstring_coverage(filename),
            'type_hints': self.check_type_hints(filename),
            'complexity': self.check_code_complexity(filename)
        }
        
        return results
    
    def run_all_checks(self, files: List[str]) -> Dict[str, Dict[str, List[str]]]:
        """Run all checks on multiple files"""
        all_results = {}
        
        for filename in files:
            if os.path.exists(filename):
                all_results[filename] = self.check_file(filename)
        
        return all_results
    
    def print_report(self, results: Dict[str, Dict[str, List[str]]]):
        """Print formatted report of all issues"""
        print("\n" + "="*80)
        print("CODE QUALITY REPORT")
        print("="*80)
        
        total_issues = 0
        
        for filename, file_results in results.items():
            file_issues = sum(len(issues) for issues in file_results.values())
            if file_issues > 0:
                print(f"\n📄 {filename}: {file_issues} issues")
                
                for category, issues in file_results.items():
                    if issues:
                        print(f"\n  {category.upper()}:")
                        for issue in issues[:10]:  # Show first 10
                            print(f"    - {issue}")
                        if len(issues) > 10:
                            print(f"    ... and {len(issues) - 10} more")
                
                total_issues += file_issues
        
        print("\n" + "="*80)
        print(f"TOTAL ISSUES: {total_issues}")
        print("="*80)
        
        return total_issues


def check_visualization_clarity() -> List[str]:
    """Check visualizations for common clarity issues"""
    issues = []
    
    with open('visualizations.py', 'r') as f:
        content = f.read()
    
    # Check for overlapping text
    if 'adjust_text' not in content:
        issues.append("Consider using 'adjustText' library to prevent label overlap")
    
    # Check for consistent figure sizes
    figsize_pattern = re.findall(r'figsize=\((\d+),\s*(\d+)\)', content)
    unique_sizes = set(figsize_pattern)
    if len(unique_sizes) > 8:
        issues.append(f"Too many different figure sizes ({len(unique_sizes)}) - consider standardizing")
    
    # Check for colorblind accessibility
    if 'colorblind' not in content.lower():
        issues.append("No mention of colorblind-friendly palettes - consider adding")
    
    # Check for font size consistency
    fontsize_pattern = re.findall(r'fontsize=(\d+)', content)
    if len(set(fontsize_pattern)) > 6:
        issues.append("Consider standardizing font sizes for consistency")
    
    return issues


def main():
    """Main function to run code quality checks"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Check and improve code quality')
    parser.add_argument('--check', action='store_true', help='Check for issues only')
    parser.add_argument('--fix', action='store_true', help='Apply fixes (not implemented yet)')
    args = parser.parse_args()
    
    # Files to check
    core_files = [
        'nonbscanner.py',
        'detectors.py',
        'utilities.py',
        'visualizations.py',
        'app.py'
    ]
    
    checker = CodeQualityChecker()
    results = checker.run_all_checks(core_files)
    
    # Print report
    total_issues = checker.print_report(results)
    
    # Check visualization clarity
    print("\n" + "="*80)
    print("VISUALIZATION CLARITY CHECK")
    print("="*80)
    viz_issues = check_visualization_clarity()
    if viz_issues:
        for issue in viz_issues:
            print(f"  - {issue}")
    else:
        print("  ✓ No visualization clarity issues found")
    
    print("\n" + "="*80)
    print(f"SUMMARY: {total_issues + len(viz_issues)} total issues found")
    print("="*80)
    
    if args.fix:
        print("\n⚠  --fix flag is not yet implemented")
        print("Please review issues and fix manually")
    
    return 0 if total_issues == 0 else 1


if __name__ == '__main__':
    exit(main())
