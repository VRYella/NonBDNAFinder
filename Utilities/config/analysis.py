"""
Analysis configuration for NBDScanner.

This module contains analysis parameters:
- Sequence processing thresholds
- Performance and display settings
- Motif filtering options
- File upload limits
- GC balance thresholds
- Results display settings

IMPORTANT THRESHOLDS DISTINCTION
--------------------------------
TWO DIFFERENT thresholds control TWO DIFFERENT features:

1. DETECTOR PARALLELIZATION THRESHOLD (in nonbscanner.py):
   - CHUNK_THRESHOLD = 50,000 bp (50KB)
   - Triggers parallel detector execution for sequences > 50KB
   - Runs 9 detectors in parallel using ThreadPoolExecutor
   - Added in PR#5 for 1.5-2x speedup

2. SEQUENCE CHUNKING THRESHOLD (in this config):
   - chunk_threshold = 1,000,000 bp (1MB)  
   - Triggers sequence splitting for sequences > 1MB
   - Splits into 50KB chunks using ProcessPoolExecutor
   - Updated to 1MB in PR#7 (was 50KB)

PERFORMANCE BEHAVIOR
--------------------
- < 50KB: Sequential detectors, no chunking
- 50KB - 1MB: Parallel detectors (1.5-2x), no chunking
- > 1MB: Parallel detectors + chunking (4-12x speedup)

OTHER PARAMETERS
----------------
GC_BALANCE_MIN = 30%
GC_BALANCE_MAX = 70%
MAX_OVERLAP_DISPLAY = 10
"""

# ==================== ANALYSIS PARAMETERS ====================
# Control sequence processing and analysis behavior
ANALYSIS_CONFIG = {
    # Sequence processing thresholds
    'chunk_threshold': 1_000_000,     # Sequences > this size use chunking (bp) - Changed from 50KB to 1MB
    'default_chunk_size': 50_000,     # Default chunk size for large sequences (bp) - ALWAYS use 50Kbp chunks
    'default_chunk_overlap': 2_000,   # 2Kbp (2,000bp) overlap ensures motifs at boundaries are captured (balanced for performance/accuracy)
    
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
    'max_file_size_mb': 100,       # Maximum file size in MB (100 MB ceiling for optimal performance)
}

# ==================== TRIPLE ADAPTIVE CHUNKING CONFIG ====================
# Three-tier hierarchical chunking for genome-scale analysis
# Automatically selects optimal strategy based on sequence size
CHUNKING_CONFIG = {
    # Micro-tier (base analysis level) - ALWAYS USE 50Kbp CHUNKS
    'micro_chunk_size': 50_000,       # 50Kbp chunks for fast analysis
    'micro_overlap': 2_000,           # 2Kbp overlap to catch boundary motifs (optimized)
    
    # Meso-tier (memory management level) - ALWAYS USE 50Kbp CHUNKS
    'meso_chunk_size': 50_000,        # 50Kbp chunks for consistency
    'meso_overlap': 2_000,            # 2Kbp overlap between meso-chunks (optimized)
    
    # Macro-tier (parallelization level) - ALWAYS USE 50Kbp CHUNKS
    'macro_chunk_size': 50_000,       # 50Kbp chunks for consistency
    'macro_overlap': 2_000,           # 2Kbp overlap between macro-chunks (optimized)
    
    # Adaptive thresholds for automatic strategy selection
    'direct_threshold': 50_000,              # <=50Kbp: direct analysis (no chunking) - CORRECTED
    'single_tier_threshold': 1_000_000,      # 50Kbp-1MB: micro-tier only - CORRECTED
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
