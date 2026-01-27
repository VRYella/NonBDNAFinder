"""
Curved DNA Pattern Generation
Dr. Venkata Rajesh Yella | 2024.1 | MIT License

Helper functions for generating phased repeat patterns for curved DNA detection.
"""

from typing import List, Tuple


def _generate_phased_repeat_patterns(base: str, num_tracts: int, tract_sizes: range, 
                                      id_start: int) -> List[Tuple[str, str, str, str, int, str, float, str, str]]:
    """
    Generate phased repeat patterns programmatically.
    
    # Output Pattern Structure:
    
    Args:
        base: Base nucleotide ('A' or 'T')
        num_tracts: Number of tracts (3, 4, or 5)
        tract_sizes: Range of tract sizes (e.g., range(3, 10))
        id_start: Starting ID number for pattern naming
    
    Returns:
        List of 9-field pattern tuples
    """
    patterns = []
    label = 'APR' if base == 'A' else 'TPR'
    scores = {3: 0.90, 4: 0.92, 5: 0.95}
    min_lens = {3: 20, 4: 25, 5: 30}
    
    pattern_id = id_start
    for size in tract_sizes:
        spacing_min = max(0, 11 - size)
        spacing_max = min(8, 11 - size + 2)
        
        repeat_unit = f'{base}{{{size}}}[ACGT]{{{spacing_min},{spacing_max}}}'
        pattern = f'(?:{repeat_unit * (num_tracts - 1)}{base}{{{size}}})'
        
        name = f'{base}{size}-{label}' + (f'-{num_tracts}' if num_tracts > 3 else '')
        patterns.append((
            pattern,
            f'CRV_{pattern_id:03d}',
            name,
            'Global Curvature',
            min_lens[num_tracts],
            'phasing_score',
            scores[num_tracts],
            f'{num_tracts}-tract {label}',
            'Koo 1986'
        ))
        pattern_id += 1
    
    return patterns
