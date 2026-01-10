#!/usr/bin/env python3
"""
Automated Detector Extraction Script
=====================================

Extracts individual detector classes from detectors.py into separate module files.
"""

import os
import re
from typing import List, Tuple

def extract_detector_class(source_file: str, class_name: str, start_line: int, end_line: int, output_file: str):
    """
    Extract a detector class and its dependencies to a new module file.
    
    Args:
        source_file: Path to detectors.py
        class_name: Name of the detector class
        start_line: Starting line number (1-based)
        end_line: Ending line number (1-based, exclusive)
        output_file: Target module path
    """
    with open(source_file, 'r') as f:
        lines = f.readlines()
    
    # Extract the class definition
    class_lines = lines[start_line-1:end_line-1]
    class_content = ''.join(class_lines)
    
    # Extract helper functions needed
    helpers = []
    if 'revcomp' in class_content:
        helpers.append(('revcomp', 316, 319))
    if '_generate_phased_repeat_patterns' in class_content:
        helpers.append(('_generate_phased_repeat_patterns', 321, 377))
    
    # Create module docstring
    module_doc = f'''"""
{class_name} Module
{'=' * (len(class_name) + 7)}

Detector class for {class_name.replace('Detector', '').replace('DNA', ' DNA')} motifs.
Extracted from detectors.py for modular architecture.

Provides specialized detection algorithms and pattern matching for this motif type.
"""

'''
    
    # Collect imports
    imports = []
    imports.append('import re')
    imports.append('from abc import ABC, abstractmethod')
    imports.append('from typing import List, Dict, Any, Tuple, Optional')
    imports.append('')
    
    # Import base detector
    imports.append('from .base import BaseMotifDetector')
    imports.append('')
    
    # Check for specific imports needed
    if 'hyperscan' in class_content.lower():
        imports.append('import logging')
        imports.append('logger = logging.getLogger(__name__)')
        imports.append('')
        imports.append('_HYPERSCAN_AVAILABLE = False')
        imports.append('try:')
        imports.append('    import hyperscan')
        imports.append('    _HYPERSCAN_AVAILABLE = True')
        imports.append('except ImportError:')
        imports.append('    pass')
        imports.append('')
    
    # Build final content
    content_parts = [module_doc]
    content_parts.append('\n'.join(imports))
    content_parts.append('\n')
    
    # Add helper functions if needed
    for helper_name, h_start, h_end in helpers:
        helper_lines = lines[h_start-1:h_end]
        content_parts.append(''.join(helper_lines))
        content_parts.append('\n\n')
    
    # Add the class
    content_parts.append(class_content)
    
    final_content = ''.join(content_parts)
    
    # Write to output file
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        f.write(final_content)
    
    print(f"✓ Extracted {class_name} to {output_file}")
    return True


def main():
    """Extract all detector classes."""
    
    source_file = 'detectors.py'
    
    # Define detector classes to extract (name, start_line, end_line, output_file)
    detectors = [
        ('CurvedDNADetector', 379, 938, 'engine/detectors/curved_dna.py'),
        ('ZDNADetector', 938, 1564, 'engine/detectors/z_dna.py'),
        ('APhilicDetector', 1564, 2100, 'engine/detectors/a_philic.py'),
        ('SlippedDNADetector', 2100, 2644, 'engine/detectors/slipped_dna.py'),
        ('CruciformDetector', 2644, 3155, 'engine/detectors/cruciform.py'),
        ('RLoopDetector', 3155, 3671, 'engine/detectors/r_loop.py'),
        ('TriplexDetector', 3671, 4156, 'engine/detectors/triplex.py'),
        ('GQuadruplexDetector', 4156, 4577, 'engine/detectors/g_quadruplex.py'),
        ('IMotifDetector', 4577, 4857, 'engine/detectors/i_motif.py'),
    ]
    
    print("="*60)
    print("Extracting Detector Classes")
    print("="*60)
    print()
    
    for class_name, start, end, output in detectors:
        extract_detector_class(source_file, class_name, start, end, output)
    
    print()
    print("="*60)
    print("Extraction Complete!")
    print("="*60)
    
    return 0


if __name__ == '__main__':
    exit(main())
