"""
Analysis configuration for NBDScanner.

This module contains analysis parameters:
- Sequence processing thresholds
- Performance and display settings
- Motif filtering options
- File upload limits
- GC balance thresholds
- Results display settings

PARAMETERS USED IN THIS MODULE
--------------------------------
CHUNK_THRESHOLD = 50_000 bp
DEFAULT_CHUNK_SIZE = 50_000 bp
DEFAULT_CHUNK_OVERLAP = 5_000 bp
GC_BALANCE_MIN = 30%
GC_BALANCE_MAX = 70%
MAX_OVERLAP_DISPLAY = 10
"""

# ==================== ANALYSIS PARAMETERS ====================
# Control sequence processing and analysis behavior
ANALYSIS_CONFIG = {
    # Sequence processing thresholds
    'chunk_threshold': 1_000_000,     # Sequences > this size use chunking (bp) - Changed from 50KB to 1MB
    'default_chunk_size': 50_000,  # Default chunk size for large sequences (bp)
    'default_chunk_overlap': 5_000,  # Overlap between chunks to catch motifs at boundaries (bp)
    
    # Performance and display settings
    'max_sequences_preview': 3,    # Number of sequences to show in file preview
    'rows_per_page': 100,          # Pagination size for large result tables
    'update_interval': 5,          # Progress update frequency (in sequences)
    
    # Motif filtering and display
    'min_score_threshold': 0.0,    # Minimum score to display (0 = show all)
    
    # *** IMPORTANT: Control which motifs appear in distributions ***
    # Set to True to include Hybrid and Cluster motifs in all visualizations
    # Set to False to exclude them from distribution plots (they'll still appear in dedicated tab)
    'include_hybrid_in_distribution': True,   # Include Hybrid motifs in plots
    'include_clusters_in_distribution': True, # Include Cluster motifs in plots
    
    # File upload limits
    'max_file_size_mb': 1024,      # Maximum file size in MB (1 GB default)
}

# ==================== TRIPLE ADAPTIVE CHUNKING CONFIG ====================
# Three-tier hierarchical chunking for genome-scale analysis
# Automatically selects optimal strategy based on sequence size
CHUNKING_CONFIG = {
    # Micro-tier (base analysis level)
    'micro_chunk_size': 50_000,       # 50KB chunks for fast analysis
    'micro_overlap': 2_000,            # 2KB overlap to catch boundary motifs
    
    # Meso-tier (memory management level)
    'meso_chunk_size': 5_000_000,     # 5MB chunks for memory efficiency
    'meso_overlap': 5_000,             # 5KB overlap between meso-chunks
    
    # Macro-tier (parallelization level)
    'macro_chunk_size': 50_000_000,   # 50MB chunks for parallel processing
    'macro_overlap': 10_000,           # 10KB overlap between macro-chunks
    
    # Adaptive thresholds for automatic strategy selection
    'direct_threshold': 50_000,              # <50KB: direct analysis (no chunking) - CORRECTED
    'single_tier_threshold': 1_000_000,      # 50KB-1MB: micro-tier only - CORRECTED
    'double_tier_threshold': 100_000_000,    # >=100MB: triple-tier begins (macro+meso+micro) - CORRECTED
    # 1MB-100MB: double-tier (meso+micro tiers)
    
    # Parallel detector execution
    'enable_parallel_detectors': True,       # Run detectors in parallel
    'detectors_per_chunk_parallel': True,    # Parallel per-chunk execution
    
    # Performance tuning
    'enable_adaptive': False,         # Disable adaptive strategy selection - use simple chunking from 1MB
    'max_workers': None,              # None = auto-detect CPU count
}

# ==================== SEQUENCE COMPOSITION THRESHOLDS ====================
# GC balance thresholds for genomic DNA analysis
# Typical balanced genomic DNA falls within 30-70% GC content
GC_BALANCE_MIN = 30   # Minimum GC% for balanced genome
GC_BALANCE_MAX = 70   # Maximum GC% for balanced genome

# ==================== RESULTS DISPLAY SETTINGS ====================
# Maximum number of overlap pairs to display in overlap matrices
MAX_OVERLAP_DISPLAY = 10
