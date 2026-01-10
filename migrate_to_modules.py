#!/usr/bin/env python3
"""
Modular Architecture Migration Script
======================================

This script automates the extraction of functionality from monolithic files
into focused, modular components following the architecture design.

Usage:
    python migrate_to_modules.py [--dry-run] [--module MODULE_NAME]

Options:
    --dry-run: Show what would be done without making changes
    --module: Extract only the specified module (e.g., 'validation')
"""

import os
import re
import ast
import argparse
from pathlib import Path
from typing import List, Dict, Tuple, Set
import shutil


class ModuleExtractor:
    """Extracts functions and classes into separate modules."""
    
    def __init__(self, source_file: str, target_dir: str, dry_run: bool = False):
        self.source_file = source_file
        self.target_dir = target_dir
        self.dry_run = dry_run
        self.extracted_items = set()
    
    def extract_functions(self, function_names: List[str], output_file: str, 
                         module_docstring: str) -> bool:
        """
        Extract specified functions from source file to new module.
        
        Args:
            function_names: List of function names to extract
            output_file: Target module filename
            module_docstring: Docstring for the new module
            
        Returns:
            True if successful
        """
        if self.dry_run:
            print(f"[DRY RUN] Would extract {len(function_names)} functions to {output_file}")
            print(f"  Functions: {', '.join(function_names)}")
            return True
        
        # Parse source file
        with open(self.source_file, 'r') as f:
            source_code = f.read()
        
        try:
            tree = ast.parse(source_code)
        except SyntaxError as e:
            print(f"Error parsing {self.source_file}: {e}")
            return False
        
        # Extract function definitions
        extracted_code = [f'"""{module_docstring}"""\n\n']
        imports = set()
        
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef) and node.name in function_names:
                # Get the function source code
                func_code = ast.get_source_segment(source_code, node)
                if func_code:
                    extracted_code.append(func_code)
                    extracted_code.append('\n\n')
                    self.extracted_items.add(node.name)
            
            # Collect imports
            elif isinstance(node, (ast.Import, ast.ImportFrom)):
                import_code = ast.get_source_segment(source_code, node)
                if import_code:
                    imports.add(import_code)
        
        # Write to target file
        target_path = os.path.join(self.target_dir, output_file)
        os.makedirs(os.path.dirname(target_path), exist_ok=True)
        
        with open(target_path, 'w') as f:
            # Write imports first
            if imports:
                f.write('\n'.join(sorted(imports)))
                f.write('\n\n')
            
            # Write extracted functions
            f.write(''.join(extracted_code))
        
        print(f"✓ Created {output_file} with {len(self.extracted_items)} functions")
        return True


def create_validation_module():
    """Extract validation functions to utils/validation.py"""
    functions = [
        'validate_sequence',
        'validate_motif',
        'validate_score',
        'validate_coordinates',
    ]
    
    docstring = """Sequence Validation Module

Functions for validating DNA sequences, motifs, scores, and coordinates.
Extracted from utilities.py for focused validation logic.
"""
    
    return functions, 'utils/validation.py', docstring


def create_export_module():
    """Extract export functions to utils/export.py"""
    functions = [
        'export_to_csv',
        'export_to_bed',
        'export_to_json',
        'export_to_excel',
        'export_to_pdf',
    ]
    
    docstring = """Data Export Module

Functions for exporting motif data to various formats (CSV, BED, JSON, Excel, PDF).
Extracted from utilities.py for focused export functionality.
"""
    
    return functions, 'utils/export.py', docstring


def create_registry_module():
    """Extract registry functions to utils/registry.py"""
    functions = [
        'load_db_for_class',
        'load_registry_for_class',
        'get_cached_registry',
        'scan_with_registry',
        '_load_consolidated_registry',
        'clear_pattern_registry_cache',
    ]
    
    docstring = """Pattern Registry Module

Functions for loading and managing pattern registries and Hyperscan databases.
Extracted from utilities.py for focused registry management.
"""
    
    return functions, 'utils/registry.py', docstring


def create_caching_module():
    """Extract caching functions to utils/caching.py"""
    functions = [
        'get_cached_registry',
        'clear_pattern_registry_cache',
        '_get_cached_scanner',
    ]
    
    docstring = """Caching Module

Functions for caching frequently used objects and data structures.
Provides performance optimization through intelligent caching.
"""
    
    return functions, 'utils/caching.py', docstring


def main():
    parser = argparse.ArgumentParser(description='Migrate to modular architecture')
    parser.add_argument('--dry-run', action='store_true', 
                       help='Show what would be done without making changes')
    parser.add_argument('--module', type=str,
                       help='Extract only the specified module')
    parser.add_argument('--backup', action='store_true',
                       help='Create backup of original files')
    
    args = parser.parse_args()
    
    # Define module extraction tasks
    extraction_tasks = [
        create_validation_module(),
        create_export_module(),
        create_registry_module(),
        create_caching_module(),
    ]
    
    # Filter by module name if specified
    if args.module:
        extraction_tasks = [
            task for task in extraction_tasks 
            if args.module.lower() in task[1].lower()
        ]
        
        if not extraction_tasks:
            print(f"No module matching '{args.module}' found")
            return 1
    
    # Create backup if requested
    if args.backup and not args.dry_run:
        backup_dir = 'backup_pre_modular'
        os.makedirs(backup_dir, exist_ok=True)
        
        for file in ['utilities.py', 'detectors.py', 'nonbscanner.py', 'app.py']:
            if os.path.exists(file):
                shutil.copy2(file, os.path.join(backup_dir, file))
                print(f"✓ Backed up {file}")
    
    # Execute extractions
    print(f"\n{'='*60}")
    print(f"Module Extraction {'(DRY RUN)' if args.dry_run else ''}")
    print(f"{'='*60}\n")
    
    for functions, output_file, docstring in extraction_tasks:
        print(f"Processing {output_file}...")
        
        # Determine source file
        if 'detector' in output_file.lower():
            source_file = 'detectors.py'
        elif output_file.startswith('ui/'):
            source_file = 'app.py'
        else:
            source_file = 'utilities.py'
        
        if not os.path.exists(source_file):
            print(f"  ⚠ Source file {source_file} not found, skipping")
            continue
        
        extractor = ModuleExtractor(source_file, '.', dry_run=args.dry_run)
        success = extractor.extract_functions(functions, output_file, docstring)
        
        if success:
            print(f"  ✓ Extracted {len(functions)} functions")
        else:
            print(f"  ✗ Extraction failed")
        
        print()
    
    print(f"{'='*60}")
    print(f"Migration {'simulation' if args.dry_run else 'complete'}!")
    print(f"{'='*60}\n")
    
    if args.dry_run:
        print("Run without --dry-run to apply changes")
    else:
        print("✓ Modules created successfully")
        print("✓ Update imports in main files")
        print("✓ Run tests to verify functionality")
    
    return 0


if __name__ == '__main__':
    exit(main())
