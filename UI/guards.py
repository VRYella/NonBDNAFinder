"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Guards Module - Validation and guard utilities for NBDScanner                │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
│ Functions for motif subclass validation, Excel file generation               │
└──────────────────────────────────────────────────────────────────────────────┘
"""
# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
import os
import tempfile
from Utilities.utilities import export_to_excel, export_multifasta_to_excel
from UI.storage_helpers import get_results, get_sequence_length

# ═══════════════════════════════════════════════════════════════════════════════
# TUNABLE PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════════
DEFAULT_SUBCLASS = 'Other'; EXCEL_SUFFIX = '.xlsx'
# ═══════════════════════════════════════════════════════════════════════════════


def ensure_subclass(motif):
    """Guarantee every motif has a string 'Subclass'"""
    if isinstance(motif, dict):
        if 'Subclass' not in motif or motif['Subclass'] is None:
            motif['Subclass'] = motif.get('Subtype', DEFAULT_SUBCLASS)
        return motif
    else:
        return {'Subclass': DEFAULT_SUBCLASS, 'Motif': motif}


def generate_excel_bytes(motifs, simple_format=True):
    """
    Generate Excel file as bytes for Streamlit download button.
    
    Args:
        motifs: List of motif dictionaries
        simple_format: If True, use 2-tab format (NonOverlappingConsolidated, OverlappingAll)
        
    Returns:
        bytes: Excel file data as bytes
    """
    # Create a temporary file
    with tempfile.NamedTemporaryFile(mode='wb', suffix='.xlsx', delete=False) as tmp:
        tmp_path = tmp.name
    
    try:
        # Use the existing export_to_excel function to create the file
        export_to_excel(motifs, tmp_path, simple_format=simple_format)
        
        # Read the file as bytes
        with open(tmp_path, 'rb') as f:
            excel_bytes = f.read()
        
        return excel_bytes
    finally:
        # Clean up the temporary file
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def generate_multifasta_excel_bytes(names, lengths, seq_count):
    """
    Generate MultiFASTA Excel file as bytes for Streamlit download button.
    
    Args:
        names: List of sequence names
        lengths: List of sequence lengths
        seq_count: Number of sequences
        
    Returns:
        bytes: Excel file data as bytes
    """
    # Collect all motifs by sequence
    annotations_by_sequence = {}
    sequence_lengths = {}
    
    for i in range(seq_count):
        fasta_id = names[i]
        motifs = get_results(i)
        annotations_by_sequence[fasta_id] = motifs
        sequence_lengths[fasta_id] = lengths[i]
    
    # Check if all sequences have equal length
    equal_length = len(set(lengths)) == 1
    seq_length = lengths[0] if equal_length else None
    
    # Create a temporary file
    with tempfile.NamedTemporaryFile(mode='wb', suffix='.xlsx', delete=False) as tmp:
        tmp_path = tmp.name
    
    try:
        # Use the export_multifasta_to_excel function to create the file
        export_multifasta_to_excel(
            annotations_by_sequence,
            sequence_lengths,
            tmp_path,
            equal_length,
            seq_length
        )
        
        # Read the file as bytes
        with open(tmp_path, 'rb') as f:
            excel_bytes = f.read()
        
        return excel_bytes
    finally:
        # Clean up the temporary file
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

