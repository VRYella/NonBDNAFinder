"""
Export configuration for NBDScanner.

This module contains export settings:
- Available export formats
- Default export options
- Core columns for exports
"""

# ==================== EXPORT FORMATS ====================
# Control data export options and default settings
EXPORT_CONFIG = {
    # Available export formats
    # Add or remove formats as needed
    'available_formats': ['CSV', 'Excel', 'JSON', 'BED'],
    
    # Default export options (can be overridden by user)
    'include_sequences': True,     # Include full motif sequences in exports
    'excel_multi_sheet': True,     # Create multi-sheet Excel workbooks (one per motif class)
    'json_pretty': True,           # Pretty-print JSON exports (more readable)
    
    # Column selection for exports
    # These are the core columns that appear in all export formats
    'core_columns': [
        'Sequence_Name',  # Name or accession of the sequence
        'Source',         # Source database or experiment
        'Class',          # Motif class (G4, Z-DNA, etc.)
        'Subclass',       # Motif subclass or subtype
        'Start',          # Start position (bp)
        'End',            # End position (bp)
        'Length',         # Motif length (bp)
        'Sequence',       # DNA sequence of motif
        'Score',          # Confidence score (0-3 scale)
    ],
}
