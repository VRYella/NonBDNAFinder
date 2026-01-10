#!/usr/bin/env python3
"""
Automated Utility Module Extraction Script
==========================================

Extracts utility functions from utilities.py and app.py into modular components.
"""

import os
import re


def extract_functions_by_line_range(source_file: str, start_line: int, end_line: int, output_file: str, module_doc: str):
    """
    Extract functions from a line range and create a new module.
    
    Args:
        source_file: Path to source file
        start_line: Starting line number (1-based, inclusive)
        end_line: Ending line number (1-based, exclusive)
        output_file: Target module path
        module_doc: Module docstring
    """
    with open(source_file, 'r') as f:
        lines = f.readlines()
    
    # Extract the content
    content_lines = lines[start_line-1:end_line-1]
    content = ''.join(content_lines)
    
    # Collect imports - look for imports at the beginning of source file
    imports = []
    for line in lines[:100]:  # Check first 100 lines for imports
        if line.strip().startswith('import ') or line.strip().startswith('from '):
            # Include common imports
            if any(pkg in line for pkg in ['typing', 'Dict', 'List', 'Any', 'pandas', 'numpy', 
                                           'plotly', 'matplotlib', 'streamlit', 'pathlib', 'os']):
                imports.append(line.rstrip())
    
    # Build final content
    final_content = f'''"""{module_doc}"""

'''
    
    if imports:
        final_content += '\n'.join(sorted(set(imports))) + '\n\n'
    
    final_content += content
    
    # Write to output file
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        f.write(final_content)
    
    print(f"✓ Created {output_file}")
    return True


def extract_registry_module():
    """Extract registry functions from utilities.py"""
    # Registry functions are scattered, so we'll extract by function names
    with open('utilities.py', 'r') as f:
        content = f.read()
        lines = content.split('\n')
    
    # Find the functions we need
    function_names = [
        '_load_consolidated_registry_from_excel',
        '_load_consolidated_registry',
        '_load_registry',
        'clear_pattern_registry_cache',
        'load_db_for_class',
        'load_registry_for_class',
        'get_cached_registry',
        'scan_with_registry',
        'get_pattern_registry',
        'get_hs_db_for_class',
    ]
    
    # Extract functions (lines 649-1755 contain registry-related functions)
    extract_functions_by_line_range(
        'utilities.py',
        649, 
        1755,
        'utils/registry.py',
        """Pattern Registry Module
=======================

Functions for loading and managing pattern registries and Hyperscan databases.
Provides centralized registry management for all detector classes.

Extracted from utilities.py for modular architecture.
"""
    )
    return True


def extract_plotting_modules():
    """Extract plotting functions into categorized modules"""
    
    # Distribution plots (plot_motif_distribution, plot_score_distribution, plot_length_distribution)
    # Lines 3975-4710
    extract_functions_by_line_range(
        'utilities.py',
        3975,
        4710,
        'utils/plotting/distributions.py',
        """Distribution Plots Module
=========================

Visualization functions for motif distributions, scores, and lengths.
Includes bar charts, histograms, and comparative visualizations.

Extracted from utilities.py for modular architecture.
"""
    )
    
    # Coverage and density plots (plot_coverage_map, plot_density_heatmap)
    # Lines 4335-4525
    extract_functions_by_line_range(
        'utilities.py',
        4335,
        4525,
        'utils/plotting/coverage.py',
        """Coverage Maps Module
===================

Visualization functions for motif coverage and spatial distribution.
Creates coverage maps and position-based visualizations.

Extracted from utilities.py for modular architecture.
"""
    )
    
    # Density plots (plot_density_heatmap, plot_density_comparison, plot_circos_motif_density)
    # Lines 4424-6335
    extract_functions_by_line_range(
        'utilities.py',
        5456,
        6335,
        'utils/plotting/density.py',
        """Density Plots Module
===================

Visualization functions for motif density analysis.
Includes heatmaps, circos plots, and density comparisons.

Extracted from utilities.py for modular architecture.
"""
    )
    
    # Statistical plots (plot_class_analysis_comprehensive, plot_score_statistics_by_class)
    # Lines 4988-5456
    extract_functions_by_line_range(
        'utilities.py',
        4988,
        5456,
        'utils/plotting/statistical.py',
        """Statistical Plots Module
========================

Statistical analysis and visualization functions.
Includes comprehensive class analysis and statistical comparisons.

Extracted from utilities.py for modular architecture.
"""
    )
    
    # Genomic plots (plot_manhattan_motif_density, plot_genome_landscape_track, etc.)
    # Lines 6335-7400
    extract_functions_by_line_range(
        'utilities.py',
        6335,
        7400,
        'utils/plotting/genomic.py',
        """Genomic Plots Module
====================

Genome-wide visualization functions.
Includes Manhattan plots, landscape tracks, and genomic overview visualizations.

Extracted from utilities.py for modular architecture.
"""
    )
    
    return True


def extract_ui_modules():
    """Extract UI components from app.py"""
    
    # First, let's check what functions exist in app.py
    with open('app.py', 'r') as f:
        lines = f.readlines()
    
    # Formatting functions (format_time, format_time_scientific, etc.)
    # Lines 1118-1166
    extract_functions_by_line_range(
        'app.py',
        1118,
        1365,
        'ui/formatting.py',
        """UI Formatting Module
====================

Text formatting and display helper functions.
Provides consistent formatting across the UI.

Extracted from app.py for modular architecture.
"""
    )
    
    # Download functions (generate_excel_bytes, etc.)
    # Lines 1367-1410
    extract_functions_by_line_range(
        'app.py',
        1367,
        1540,
        'ui/downloads.py',
        """UI Downloads Module
===================

Download and export functionality for UI.
Handles file generation and download buttons.

Extracted from app.py for modular architecture.
"""
    )
    
    return True


def main():
    """Extract all utility modules."""
    
    print("="*60)
    print("Extracting Utility Modules")
    print("="*60)
    print()
    
    print("Phase 4: Utility Modules")
    print("-" * 60)
    
    # Registry module
    print("\n[1/6] Extracting registry module...")
    extract_registry_module()
    
    # Plotting modules
    print("\n[2/6] Extracting plotting modules...")
    extract_plotting_modules()
    
    # UI modules
    print("\n[3/6] Extracting UI modules...")
    extract_ui_modules()
    
    # Create __init__.py for plotting subpackage
    os.makedirs('utils/plotting', exist_ok=True)
    with open('utils/plotting/__init__.py', 'w') as f:
        f.write('''"""
Plotting Module
==============

Visualization functions for Non-B DNA motif analysis.

Submodules:
- distributions: Distribution plots
- coverage: Coverage maps
- density: Density analysis
- statistical: Statistical plots
- genomic: Genome-wide visualizations
"""

__version__ = "2025.1"
''')
    
    print()
    print("="*60)
    print("Extraction Complete!")
    print("="*60)
    
    return 0


if __name__ == '__main__':
    exit(main())
