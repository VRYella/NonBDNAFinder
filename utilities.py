"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    CONSOLIDATED UTILITIES MODULE                              ║
║     Utility Functions for Sequence Processing, Export, and Visualization     ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: utilities.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.1 - Consolidated with Visualizations
LICENSE: MIT

DESCRIPTION:
    Consolidated module containing all utility functions for sequence processing,
    pattern loading, data export, visualization, and validation. Combines functionality
    from utils, load_hsdb, load_regex_registry, motif_patterns, canonicalize_motif,
    and visualizations modules into a single cohesive unit.

MAIN FUNCTIONS:
    Sequence Processing:
    - parse_fasta(): Parse FASTA format sequences
    - gc_content(): Calculate GC content
    - reverse_complement(): Generate reverse complement
    - validate_sequence(): Validate DNA sequence

    Data Export:
    - export_to_csv(): Export to CSV format
    - export_to_bed(): Export to BED format
    - export_to_json(): Export to JSON format
    - export_to_excel(): Export to Excel format

    Pattern Loading:
    - load_hyperscan_db(): Load pre-compiled Hyperscan databases
    - load_regex_registry(): Load regex pattern registries
    - get_motif_patterns(): Get patterns for specific motif classes

    Statistics:
    - get_basic_stats(): Calculate sequence statistics
    - quality_check_motifs(): Validate motif quality
    
    Visualizations:
    - plot_motif_distribution(): Distribution plots
    - plot_coverage_map(): Coverage visualization
    - plot_density_heatmap(): Density heatmaps
    - plot_nested_pie_chart(): Hierarchical charts
    - plot_manhattan_motif_density(): Manhattan plots
    - and 20+ more visualization functions
"""

from typing import Dict, Any, List, Optional, Tuple, Union
import json
import re
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import pandas as pd
import numpy as np
from collections import Counter, defaultdict
import warnings
warnings.filterwarnings("ignore")

# Try to import plotly for interactive plots
try:
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

# =============================================================================
# CONSTANTS FOR DATA EXPORT
# =============================================================================

# Core columns for output tables (as per requirements)
# These columns appear in display tables, CSV exports, and are the primary columns
# Additional detailed columns only appear in Excel download per-motif sheets
CORE_OUTPUT_COLUMNS = [
    'Sequence_Name',  # Sequence Name (or Accession)
    'Source',         # Source (e.g., genome, experiment, study)
    'Class',          # Motif Class
    'Subclass',       # Motif Subclass
    'Start',          # Start Position
    'End',            # End Position
    'Length',         # Length (bp)
    'Sequence',       # Sequence
    'Score',          # Motif Score
]

def canonicalize_motif(m: Dict[str, Any]) -> Dict[str, Any]:
    """
    Canonicalize motif dictionary to standard format with consistent field names.
    
    This function ensures all motif dictionaries follow the same structure regardless
    of their source, making them compatible with export and visualization functions.
    Maintains backward compatibility with legacy field names.
    
    Args:
        m: Motif dictionary with potentially inconsistent field names
        
    Returns:
        Standardized motif dictionary with canonical field names:
        - Class: Motif class (e.g., 'G-Quadruplex', 'Z-DNA')
        - Subclass: Motif subclass/subtype
        - Start: Start position (int)
        - End: End position (int)
        - Length: Motif length in bp (int)
        - Score: Normalized score (float)
        - Actual_Score: Raw detector score (float)
        - Normalized_Score: Legacy field, set to 0 if not present
        - Motif: Actual sequence string
        - Sequence_Name: Name of parent sequence
        
    Example:
        >>> motif = {'Type': 'G-Quadruplex', 'Start': 100, 'End': 120}
        >>> canonical = canonicalize_motif(motif)
        >>> canonical['Class']
        'G-Quadruplex'
        >>> canonical['Length']
        20
    """
    mapping = {
        'Actual Score': 'Actual_Score',
        'ActualScore': 'Actual_Score',
        'Score': 'Score',
        'Normalized Score': 'Normalized_Score',
        'Normalized_Score': 'Normalized_Score',
        'Class': 'Class',
        'Type': 'Class',
        'Subclass': 'Subclass',
        'Subtype': 'Subclass',
        'Start': 'Start',
        'End': 'End',
        'Length': 'Length',
        'Sequence_Name': 'Sequence_Name',
        'Motif': 'Motif'
    }
    out = {}
    for new_key in mapping.values():
        aliases = [k for k, v in mapping.items() if v == new_key]
        val = None
        for alias in aliases:
            if alias in m:
                val = m[alias]
                break
        out[new_key] = val
    if out.get('Start') and out.get('End') and not out.get('Length'):
        out['Length'] = int(out['End']) - int(out['Start'])
    out['Class'] = str(out.get('Class') or 'Unknown')
    out['Subclass'] = str(out.get('Subclass') or 'Other')
    out['Actual_Score'] = float(out.get('Actual_Score') or m.get('Actual_Score') or m.get('Score') or 0.0)
    out['Score'] = float(out.get('Score') or out['Actual_Score'])
    # Keep Normalized_Score for backward compatibility if it exists, but set to 0 if not present
    out['Normalized_Score'] = float(out.get('Normalized_Score') or m.get('Normalized_Score') or 0.0)
    out['Motif'] = m.get('Motif') or m.get('matched_seq') or ''
    out['Sequence_Name'] = out.get('Sequence_Name') or m.get('sequence_name') or ''
    return out


# =============================================================================
# SCORE NORMALIZATION (1-3 SCALE)
# =============================================================================
"""
Scientific Score Normalization System
=====================================

Normalizes raw detector scores to a 1-3 scale where:
  - 1 = Minimal (weak motif formation potential)
  - 2 = Moderate (typical motif formation potential)
  - 3 = Highest (strong motif formation potential)

Each motif class has specific normalization parameters based on:
  - Scientific literature on formation potential
  - Biological significance thresholds
  - Experimental validation data
"""

# Scientific normalization parameters per motif class
# Format: {'class_name': {'min_raw': min_threshold, 'max_raw': max_threshold, 'method': str}}
SCORE_NORMALIZATION_PARAMS = {
    'Curved_DNA': {
        'min_raw': 0.1,   # Minimal curvature threshold
        'max_raw': 0.95,  # Strong curvature (APR phasing)
        'method': 'linear',
        'ref': 'Koo 1986, Olson 1998'
    },
    'Slipped_DNA': {
        'min_raw': 0.3,   # Minimal repeat instability
        'max_raw': 0.98,  # High repeat instability (STRs)
        'method': 'linear',
        'ref': 'Schlötterer 2000, Weber 1989'
    },
    'Cruciform': {
        'min_raw': 0.5,   # Minimum palindrome stability
        'max_raw': 0.95,  # Strong inverted repeat
        'method': 'linear',
        'ref': 'Lilley 2000, Sinden 1994'
    },
    'R-Loop': {
        'min_raw': 0.4,   # Minimal R-loop potential
        'max_raw': 0.95,  # High R-loop formation
        'method': 'linear',
        'ref': 'Aguilera 2012, Jenjaroenpun 2016'
    },
    'Triplex': {
        'min_raw': 0.5,   # Minimal triplex potential
        'max_raw': 0.95,  # Strong H-DNA formation
        'method': 'linear',
        'ref': 'Frank-Kamenetskii 1995'
    },
    'G-Quadruplex': {
        'min_raw': 0.5,   # G4Hunter threshold ~1.2 normalized
        'max_raw': 1.0,   # Perfect G4 structure
        'method': 'g4hunter',  # Special G4Hunter-based normalization
        'ref': 'Bedrat 2016, Burge 2006'
    },
    'i-Motif': {
        'min_raw': 0.4,   # Minimal i-motif stability
        'max_raw': 0.98,  # High C-tract density
        'method': 'linear',
        'ref': 'Gehring 1993, Zeraati 2018'
    },
    # NOTE: Z-DNA uses cumulative 10-mer scoring (Ho et al. 1986)
    # Individual 10-mer scores range from ~50-63, but regions with multiple
    # overlapping 10-mers accumulate scores, hence the larger max_raw value.
    # This is intentionally different from other 0-1 normalized scores.
    'Z-DNA': {
        'min_raw': 50.0,    # Single Z-DNA 10-mer score minimum (Ho 1986)
        'max_raw': 2000.0,  # Cumulative max for long Z-DNA regions (~30+ bp)
        'method': 'zdna_cumulative',  # Uses log-linear scaling for cumulative scores
        'ref': 'Ho 1986, Rich 1984'
    },
    'A-philic_DNA': {
        'min_raw': 0.5,   # Minimal A-philic propensity
        'max_raw': 0.95,  # High A-tract density
        'method': 'linear',
        'ref': 'Gorin 1995, Vinogradov 2003'
    },
    'Hybrid': {
        'min_raw': 0.3,   # Minimal hybrid score
        'max_raw': 0.9,   # Strong hybrid overlap
        'method': 'linear',
        'ref': 'This work'
    },
    'Non-B_DNA_Clusters': {
        'min_raw': 0.3,   # Minimal cluster score
        'max_raw': 0.9,   # High-density cluster
        'method': 'linear',
        'ref': 'This work'
    }
}


def normalize_score_to_1_3(raw_score: float, motif_class: str) -> float:
    """
    Normalize raw detector score to 1-3 scale.
    
    Scientific basis:
    - 1.0 = Minimal (below typical formation threshold)
    - 1.5 = Low (marginal formation potential)
    - 2.0 = Moderate (typical formation potential)
    - 2.5 = High (strong formation potential)
    - 3.0 = Highest (exceptional formation potential)
    
    Args:
        raw_score: Raw score from detector (typically 0-1 or class-specific range)
        motif_class: Motif class name for class-specific normalization
        
    Returns:
        Normalized score between 1.0 and 3.0
        
    Example:
        >>> normalize_score_to_1_3(0.8, 'G-Quadruplex')
        2.6
    """
    # Get normalization parameters for this class
    params = SCORE_NORMALIZATION_PARAMS.get(motif_class)
    
    if params is None:
        # Default linear normalization for unknown classes
        params = {'min_raw': 0.0, 'max_raw': 1.0, 'method': 'linear'}
    
    min_raw = params['min_raw']
    max_raw = params['max_raw']
    method = params.get('method', 'linear')
    
    # Handle special normalization methods
    if method == 'g4hunter':
        # G4Hunter scores: abs value indicates formation potential
        # Typical range is 0-2+, with 1.2+ being significant
        # Always use absolute value for G4Hunter scores
        effective_score = abs(raw_score)
        if effective_score < 0.5:
            return 1.0
        elif effective_score >= 1.0:
            return min(3.0, 2.0 + (effective_score - 0.5) * 2.0)
        else:
            return 1.0 + (effective_score - 0.5) * 2.0
    
    elif method == 'zdna_cumulative':
        # Z-DNA cumulative scores (sum of multiple 10-mer matches)
        # Use log-linear scaling for better distribution
        import math
        if raw_score < min_raw:
            return 1.0
        elif raw_score >= max_raw:
            return 3.0
        else:
            # Log-scale normalization for cumulative scores
            log_score = math.log10(raw_score + 1)
            log_min = math.log10(min_raw + 1)
            log_max = math.log10(max_raw + 1)
            normalized = (log_score - log_min) / (log_max - log_min)
            return 1.0 + normalized * 2.0
    
    elif method == 'zdna_10mer':
        # Z-DNA 10-mer scores range from ~50 to 63 (deprecated, use zdna_cumulative)
        if raw_score < min_raw:
            return 1.0
        elif raw_score >= max_raw:
            return 3.0
        else:
            normalized = (raw_score - min_raw) / (max_raw - min_raw)
            return 1.0 + normalized * 2.0
    
    else:  # Linear normalization (default)
        if raw_score <= min_raw:
            return 1.0
        elif raw_score >= max_raw:
            return 3.0
        else:
            # Linear interpolation between 1 and 3
            normalized = (raw_score - min_raw) / (max_raw - min_raw)
            return 1.0 + normalized * 2.0


def normalize_motif_scores(motifs: list) -> list:
    """
    Normalize all motif scores in a list to 1-3 scale.
    
    Args:
        motifs: List of motif dictionaries with 'Score' and 'Class' fields
        
    Returns:
        List of motifs with updated 'Score' field (1-3 scale)
        
    Note:
        Original raw score is preserved in 'Raw_Score' field
    """
    normalized_motifs = []
    
    for motif in motifs:
        m = motif.copy()
        raw_score = m.get('Score', 0.0)
        motif_class = m.get('Class', 'Unknown')
        
        # Store raw score for reference
        m['Raw_Score'] = raw_score
        
        # Normalize to 1-3 scale
        m['Score'] = round(normalize_score_to_1_3(raw_score, motif_class), 2)
        
        normalized_motifs.append(m)
    
    return normalized_motifs


"""
Helper to load Hyperscan DB or compile from patterns when needed.

API:
    load_db_for_class(class_name: str, registry_dir: str) -> (db, id_to_pattern, id_to_score)
      - db: hyperscan.Database instance (or None if hyperscan not available)
      - id_to_pattern: dict mapping id -> pattern string (tenmer for 10-mer, regex for others)
      - id_to_score: dict mapping id -> score (float)
      
    Supports both:
      - 10-mer patterns (ZDNA, APhilic): patterns have 'tenmer' key
      - Regex patterns (G4, IMotif, etc.): patterns have 'pattern' key
"""

import os
import logging
import json

logger = logging.getLogger(__name__)

try:
    import hyperscan
    _HYPERSCAN_AVAILABLE = True
except Exception:
    hyperscan = None
    _HYPERSCAN_AVAILABLE = False

# Cache for consolidated registry
_CONSOLIDATED_REGISTRY = None

# Try to import pandas for Excel support
_PANDAS_AVAILABLE = False
_PANDAS_WARNING_LOGGED = False
try:
    import pandas as pd
    _PANDAS_AVAILABLE = True
except ImportError:
    pass  # Will log warning on first use attempt


# Pattern registry version
PATTERN_REGISTRY_VERSION = "2024.2"

def _load_consolidated_registry_from_excel():
    """Load the consolidated registry from Excel file and cache it
    
    Returns:
        dict: Registry data in the same format as consolidated_registry.json
    """
    global _PANDAS_WARNING_LOGGED
    
    if not _PANDAS_AVAILABLE:
        if not _PANDAS_WARNING_LOGGED:
            logger.warning("pandas not available - Excel pattern loading disabled. Install with: pip install pandas openpyxl")
            _PANDAS_WARNING_LOGGED = True
        return None
    
    # Try to load pattern_registry2.xlsx from current directory (updated with normalized scores)
    excel_path = "pattern_registry2.xlsx"
    if not os.path.isfile(excel_path):
        return None
    
    try:
        # Read all sheets from Excel
        excel_data = pd.read_excel(excel_path, sheet_name=None)  # Load all sheets
        
        # Build consolidated registry structure
        registries = {}
        for sheet_name, df in excel_data.items():
            if sheet_name == 'Summary':
                continue  # Skip summary sheet
            
            # Convert DataFrame to list of dictionaries
            patterns = df.to_dict('records')
            
            # Clean up NaN values
            for pattern in patterns:
                for key, value in list(pattern.items()):
                    if pd.isna(value):
                        pattern[key] = None
            
            # Create registry structure
            registries[sheet_name] = {
                'class': sheet_name,
                'n_patterns': len(patterns),
                'patterns': patterns
            }
        
        registry_data = {
            'version': f'{PATTERN_REGISTRY_VERSION}-excel',
            'source': excel_path,
            'registries': registries,
            'total_patterns': sum(r['n_patterns'] for r in registries.values()),
            'total_classes': len(registries)
        }
        
        logger.info(f"Loaded consolidated registry from {excel_path} ({registry_data['total_patterns']} patterns)")
        return registry_data
        
    except Exception as e:
        logger.error(f"Failed to load Excel registry from {excel_path}: {e}")
        return None


def _load_consolidated_registry_from_json():
    """Load the consolidated registry from JSON file (fallback)
    
    Returns:
        dict: Registry data from consolidated_registry.json
    """
    # Try to load consolidated_registry.json from current directory
    consolidated_path = "consolidated_registry.json"
    if not os.path.isfile(consolidated_path):
        return None
    
    try:
        with open(consolidated_path, "r") as fh:
            registry = json.load(fh)
            logger.info("Loaded consolidated registry from consolidated_registry.json")
            return registry
    except Exception as e:
        logger.error(f"Failed to load JSON registry from {consolidated_path}: {e}")
        return None


def _load_consolidated_registry():
    """Load the consolidated registry file once and cache it
    
    Tries to load from Excel first (pattern_registry2.xlsx with updated normalized scores), then falls back to JSON.
    
    Returns:
        dict: Consolidated registry data
    """
    global _CONSOLIDATED_REGISTRY
    if _CONSOLIDATED_REGISTRY is not None:
        return _CONSOLIDATED_REGISTRY
    
    # Try Excel first (preferred)
    _CONSOLIDATED_REGISTRY = _load_consolidated_registry_from_excel()
    if _CONSOLIDATED_REGISTRY is not None:
        return _CONSOLIDATED_REGISTRY
    
    # Fall back to JSON
    _CONSOLIDATED_REGISTRY = _load_consolidated_registry_from_json()
    if _CONSOLIDATED_REGISTRY is not None:
        return _CONSOLIDATED_REGISTRY
    
    return None


def _load_registry(registry_dir: str, class_name: str):
    """Load registry from consolidated file (Excel or JSON)"""
    # Load from consolidated registry
    consolidated = _load_consolidated_registry()
    if consolidated and "registries" in consolidated:
        if class_name in consolidated["registries"]:
            logger.debug(f"Loading {class_name} from consolidated registry")
            return consolidated["registries"][class_name]
    
    raise FileNotFoundError(f"No registry found for {class_name} in pattern registry (tried Excel and JSON)")


def clear_pattern_registry_cache():
    """
    Clear the cached pattern registry data.
    
    This is mainly useful for testing or when you need to reload
    the pattern data after modifying the Excel or JSON files.
    
    Example:
        >>> from utilities import clear_pattern_registry_cache, _load_consolidated_registry
        >>> clear_pattern_registry_cache()
        >>> registry = _load_consolidated_registry()  # Reloads from file
    """
    global _CONSOLIDATED_REGISTRY
    _CONSOLIDATED_REGISTRY = None
    logger.debug("Pattern registry cache cleared")


def load_db_for_class(class_name: str, registry_dir: str = "registry"):
    """
    Returns (db, id_to_pattern, id_to_score).
    db is a hyperscan.Database instance (or None if hyperscan not available).
    id_to_pattern and id_to_score are dicts mapping integer id -> pattern/tenmer / score.
    
    Handles both 10-mer patterns (with 'tenmer' key) and regex patterns (with 'pattern' key).
    """
    reg = _load_registry(registry_dir, class_name)
    patterns = reg.get("patterns", [])
    
    # Handle both 10-mer registries (tenmer) and regex registries (pattern)
    id_to_pattern = {}
    for p in patterns:
        pattern_id = int(p["id"])
        # Try 'tenmer' first (for 10-mer patterns like ZDNA, APhilic)
        if "tenmer" in p:
            id_to_pattern[pattern_id] = p["tenmer"]
        # Fall back to 'pattern' (for regex patterns like G4, IMotif, etc.)
        elif "pattern" in p:
            id_to_pattern[pattern_id] = p["pattern"]
        else:
            logger.warning(f"Pattern {pattern_id} in {class_name} has neither 'tenmer' nor 'pattern' key")
    
    id_to_score = {int(p["id"]): float(p.get("score", 0.0)) for p in patterns}

    db = None
    hsdb_path = os.path.join(registry_dir, f"{class_name}.hsdb")
    if _HYPERSCAN_AVAILABLE:
        try:
            db = hyperscan.Database()
            if os.path.isfile(hsdb_path):
                # Many hyperscan Python bindings provide deserialize()
                with open(hsdb_path, "rb") as fh:
                    raw = fh.read()
                try:
                    db.deserialize(raw)
                    logger.info(f"Loaded serialized DB for {class_name} from {hsdb_path}")
                except Exception:
                    # Fallback: compile from patterns
                    expressions = [id_to_pattern[i].encode("ascii") for i in sorted(id_to_pattern.keys())]
                    ids = sorted(id_to_pattern.keys())
                    db.compile(expressions=expressions, ids=ids, elements=len(expressions))
                    logger.warning(f"deserialize() failed; compiled DB for {class_name} from patterns")
            else:
                # No serialized DB; compile from patterns
                expressions = [id_to_pattern[i].encode("ascii") for i in sorted(id_to_pattern.keys())]
                ids = sorted(id_to_pattern.keys())
                db.compile(expressions=expressions, ids=ids, elements=len(expressions))
                logger.info(f"Compiled DB for {class_name} from patterns in {registry_dir}")
        except Exception as e:
            logger.error(f"Hyperscan operations failed: {e}")
            db = None
    else:
        logger.debug("Hyperscan not installed; returning None DB (use pure-Python matcher).")

    return db, id_to_pattern, id_to_score
"""
Helper to load regex pattern registries and compile them with Hyperscan.

This module provides utilities for loading and using regex-based pattern registries
(CurvedDNA, G4, IMotif) similar to how the 10-mer registries work for A-philic and Z-DNA.

API:
    load_registry_for_class(class_name: str, registry_dir: str) 
        -> (db, id_to_pattern, id_to_subclass, id_to_score)
      - db: hyperscan.Database instance (or None if hyperscan not available)
      - id_to_pattern: dict mapping id -> regex pattern (string)
      - id_to_subclass: dict mapping id -> subclass name (string)
      - id_to_score: dict mapping id -> score (float)
"""


def load_registry_for_class(class_name: str, registry_dir: str = "registry"):
    """
    Load regex pattern registry and compile with Hyperscan.
    
    Returns (db, id_to_pattern, id_to_subclass, id_to_score).
    - db is a hyperscan.Database instance (or None if hyperscan not available).
    - id_to_pattern: dict mapping integer id -> regex pattern string
    - id_to_subclass: dict mapping integer id -> subclass name
    - id_to_score: dict mapping integer id -> score value
    """
    # Load registry data
    reg = _load_registry(registry_dir, class_name)
    patterns = reg.get("patterns", [])
    
    # Extract mappings
    id_to_pattern = {int(p["id"]): p["pattern"] for p in patterns}
    id_to_subclass = {int(p["id"]): p.get("subclass", "unknown") for p in patterns}
    id_to_score = {int(p["id"]): float(p.get("score", 0.0)) for p in patterns}
    
    # Compile Hyperscan database if available
    db = None
    if _HYPERSCAN_AVAILABLE:
        try:
            # Try to load pre-compiled DB first
            hsdb_path = os.path.join(registry_dir, f"{class_name}.hsdb")
            if os.path.isfile(hsdb_path):
                try:
                    db = hyperscan.Database()
                    with open(hsdb_path, "rb") as fh:
                        raw = fh.read()
                    db.deserialize(raw)
                    logger.info(f"Loaded serialized Hyperscan DB for {class_name}")
                except Exception as e:
                    logger.warning(f"Failed to deserialize {class_name}.hsdb: {e}")
                    db = None
            
            # If no pre-compiled DB, compile from patterns
            if db is None:
                expressions = []
                ids = []
                flags = []
                
                for i in sorted(id_to_pattern.keys()):
                    expressions.append(id_to_pattern[i].encode("ascii"))
                    ids.append(i)
                    # Use CASELESS and DOTALL flags for DNA matching
                    flags.append(hyperscan.HS_FLAG_CASELESS | hyperscan.HS_FLAG_DOTALL)
                
                db = hyperscan.Database()
                db.compile(
                    expressions=expressions,
                    ids=ids,
                    elements=len(expressions),
                    flags=flags
                )
                logger.info(f"Compiled Hyperscan DB for {class_name} from {len(expressions)} patterns")
                
        except Exception as e:
            logger.error(f"Hyperscan compilation failed for {class_name}: {e}")
            db = None
    else:
        logger.debug(f"Hyperscan not available for {class_name}; use pure-Python matching")
    
    return db, id_to_pattern, id_to_subclass, id_to_score


# Cache compiled databases to avoid recompilation
_REGISTRY_CACHE = {}


def get_cached_registry(class_name: str, registry_dir: str = "registry"):
    """
    Get cached registry or load and compile it.
    Returns (db, id_to_pattern, id_to_subclass, id_to_score).
    """
    cache_key = f"{registry_dir}/{class_name}"
    
    if cache_key not in _REGISTRY_CACHE:
        _REGISTRY_CACHE[cache_key] = load_registry_for_class(class_name, registry_dir)
    
    return _REGISTRY_CACHE[cache_key]


def scan_with_registry(class_name: str, sequence: str, registry_dir: str = "registry"):
    """
    Scan a sequence using a registry's compiled Hyperscan DB.
    
    Returns list of (start, end, pattern_id, subclass) matches.
    Falls back to pure-Python regex matching if Hyperscan not available.
    """
    db, id_to_pattern, id_to_subclass, id_to_score = get_cached_registry(class_name, registry_dir)
    matches = []
    
    if db is not None and _HYPERSCAN_AVAILABLE:
        # Use Hyperscan for fast matching
        def on_match(pattern_id, start, end, flags, context):
            subclass = id_to_subclass.get(pattern_id, "unknown")
            matches.append((start, end, pattern_id, subclass))
        
        try:
            db.scan(sequence.encode(), match_event_handler=on_match)
        except Exception as e:
            logger.error(f"Hyperscan scan failed for {class_name}: {e}")
            return []
    else:
        # Fallback to pure-Python regex matching
        import re
        for pattern_id in sorted(id_to_pattern.keys()):
            pattern = id_to_pattern[pattern_id]
            subclass = id_to_subclass[pattern_id]
            
            try:
                compiled_pattern = re.compile(pattern, re.IGNORECASE)
                for match in compiled_pattern.finditer(sequence):
                    matches.append((match.start(), match.end(), pattern_id, subclass))
            except re.error as e:
                logger.warning(f"Invalid regex pattern {pattern_id} in {class_name}: {e}")
                continue
    
    # Sort by start position
    matches.sort(key=lambda x: x[0])
    return matches
"""
Non-B DNA Motif Patterns Registry
=================================

Comprehensive pattern library with scoring algorithms for all 11 motif classes.
Optimized for scientific accuracy and Hyperscan compatibility.

PATTERN CLASSIFICATION TABLE:
=============================
Class | Patterns | Subclasses | Scoring Method            | References
------|----------|------------|---------------------------|------------------
  1   |    15    |     2      | A-tract curvature         | Olson 1998
  2   |    37    |     2      | Repeat instability        | Wells 2005  
  3   |    89    |     1      | Palindrome stability      | Lilley 2000
  4   |    11    |     1      | R-loop formation          | Skourti 2019
  5   |    14    |     2      | Triplex potential         | Frank-Kamenetskii 1995
  6   |    16    |     7      | G4Hunter algorithm        | Bedrat 2016
  7   |    12    |     3      | i-motif pH stability      | Zeraati 2018
  8   |    13    |     2      | Z-DNA transition          | Ho 1986
  9   |     8    |     1      | A-philic propensity       | Gorin 1995
 10   |  Dynamic |  Dynamic   | Overlap analysis          | This work
 11   |  Dynamic |  Dynamic   | Clustering algorithm      | This work

Total: 207+ patterns across 22+ subclasses
"""

import os
import re
import numpy as np
from typing import Dict, List, Tuple, Any, Optional
from collections import defaultdict

# Try to import Hyperscan for high-performance pattern matching
try:
    import hyperscan
    HYPERSCAN_AVAILABLE = True
except (ImportError, Exception):
    HYPERSCAN_AVAILABLE = False

# =============================================================================
# COMPREHENSIVE PATTERN REGISTRY
# =============================================================================

class PatternRegistry:
    """Complete registry of all Non-B DNA patterns with metadata"""
    
    # Class 1: Curved DNA - A-tract mediated DNA bending (optimized)
    CURVED_DNA_PATTERNS = {
        'a_tracts': [
            # (pattern, pattern_id, name, subclass, min_len, scoring_method, confidence, biological_significance, reference)
            (r'A{4,15}', 'CRV_1_1', 'A-tract', 'Local Curvature', 4, 'curvature_score', 0.95, 'DNA bending', 'Crothers 1992'),
            (r'T{4,15}', 'CRV_1_2', 'T-tract', 'Local Curvature', 4, 'curvature_score', 0.95, 'DNA bending', 'Crothers 1992'),
            (r'(?:A{3,8}T{1,5}){2,}', 'CRV_1_3', 'AT-rich tract', 'Local Curvature', 8, 'curvature_score', 0.85, 'Sequence-dependent bending', 'Hagerman 1986'),
        ],
        'phased_a_tracts': [
            (r'(?:A{3,8}.{8,12}){3,}A{3,8}', 'CRV_1_4', 'Phased A-tracts', 'Global Curvature', 20, 'phasing_score', 0.90, 'Macroscopic curvature', 'Koo 1986'),
            (r'(?:T{3,8}.{8,12}){3,}T{3,8}', 'CRV_1_5', 'Phased T-tracts', 'Global Curvature', 20, 'phasing_score', 0.90, 'Macroscopic curvature', 'Koo 1986'),
        ]
    }
    
    # Class 2: Slipped DNA - Tandem repeat-induced slippage
    SLIPPED_DNA_PATTERNS = {
        'short_tandem_repeats': [
            (r'([ATGC])\1{9,}', 'SLP_2_1', 'Mononucleotide repeat', 'STR', 10, 'instability_score', 0.98, 'Replication slippage', 'Schlötterer 2000'),
            (r'([ATGC]{2})\1{4,}', 'SLP_2_2', 'Dinucleotide repeat', 'STR', 10, 'instability_score', 0.95, 'Microsatellite instability', 'Weber 1989'),
            (r'([ATGC]{3})\1{3,}', 'SLP_2_3', 'Trinucleotide repeat', 'STR', 12, 'instability_score', 0.92, 'Expansion diseases', 'Ashley 1993'),
            (r'([ATGC]{4})\1{2,}', 'SLP_2_4', 'Tetranucleotide repeat', 'STR', 12, 'instability_score', 0.85, 'Genetic polymorphisms', 'Edwards 1991'),
            (r'(CA)\1{4,}', 'SLP_2_5', 'CA repeat', 'STR', 10, 'instability_score', 0.95, 'Common microsatellite', 'Weber 1989'),
            (r'(CGG)\1{3,}', 'SLP_2_6', 'CGG repeat', 'STR', 12, 'instability_score', 0.90, 'Fragile X syndrome', 'Verkerk 1991'),
        ],
        'direct_repeats': [
            (r'([ATGC]{5,20})(?:[ATGC]{0,100})\1', 'SLP_2_7', 'Direct repeat', 'Direct Repeat', 10, 'repeat_score', 0.80, 'Recombination hotspots', 'Jeffreys 1985'),
            (r'([ATGC]{10,50})(?:[ATGC]{0,200})\1', 'SLP_2_8', 'Long direct repeat', 'Direct Repeat', 20, 'repeat_score', 0.75, 'Genomic instability', 'Lupski 1998'),
        ]
    }
    
    # Class 3: Cruciform DNA - Inverted repeat-induced four-way junctions  
    CRUCIFORM_PATTERNS = {
        'inverted_repeats': [
            (r'([ATGC]{6,20})[ATGC]{0,50}', 'CRU_3_1', 'Potential palindrome', 'Inverted Repeats', 12, 'cruciform_stability', 0.95, 'DNA secondary structure', 'Lilley 2000'),
            (r'([ATGC]{8,15})[ATGC]{2,20}([ATGC]{8,15})', 'CRU_3_2', 'Inverted repeat candidate', 'Inverted Repeats', 16, 'cruciform_stability', 0.80, 'Secondary structure prone', 'Pearson 1996'),
            (r'([ATGC]{4,10})[ATGC]{0,10}([ATGC]{4,10})', 'CRU_3_3', 'Short inverted repeat', 'Inverted Repeats', 8, 'cruciform_stability', 0.70, 'Local secondary structure', 'Sinden 1994'),
        ]
    }
    
    # Class 4: R-loop - RNA-DNA hybrid structures
    R_LOOP_PATTERNS = {
        'r_loop_formation_sites': [
            (r'[GC]{10,}[AT]{2,10}[GC]{10,}', 'RLP_4_1', 'GC-rich R-loop site', 'R-loop formation sites', 20, 'r_loop_potential', 0.85, 'Transcription-replication conflicts', 'Aguilera 2012'),
            (r'G{5,}[ATGC]{10,100}C{5,}', 'RLP_4_2', 'G-C rich region', 'R-loop formation sites', 20, 'r_loop_potential', 0.80, 'R-loop prone regions', 'Ginno 2012'),
            (r'[GC]{6,}[AT]{1,5}[GC]{6,}', 'RLP_4_3', 'GC-AT pattern', 'R-loop formation sites', 15, 'r_loop_potential', 0.75, 'Transcriptional pausing', 'Skourti-Stathaki 2011'),
        ],
        'qmrlfs_model_1': [
            (r'G{3,}[ATCGU]{1,10}?G{3,}(?:[ATCGU]{1,10}?G{3,}){1,}?', 'QmRLFS_4_1', 'QmRLFS Model 1', 'QmRLFS-m1', 25, 'qmrlfs_score', 0.90, 'RIZ detection with 3+ G tracts', 'Jenjaroenpun 2016'),
        ],
        'qmrlfs_model_2': [
            (r'G{4,}(?:[ATCGU]{1,10}?G{4,}){1,}?', 'QmRLFS_4_2', 'QmRLFS Model 2', 'QmRLFS-m2', 30, 'qmrlfs_score', 0.95, 'RIZ detection with 4+ G tracts', 'Jenjaroenpun 2016'),
        ]
    }
    
    # Class 5: Triplex DNA - Three-stranded structures (optimized)
    TRIPLEX_PATTERNS = {
        'triplex_forming_sequences': [
            (r'[GA]{15,}', 'TRX_5_1', 'Homopurine tract', 'Triplex', 15, 'triplex_potential', 0.90, 'H-DNA formation', 'Frank-Kamenetskii 1995'),
            (r'[CT]{15,}', 'TRX_5_2', 'Homopyrimidine tract', 'Triplex', 15, 'triplex_potential', 0.90, 'H-DNA formation', 'Frank-Kamenetskii 1995'),
            (r'(?:GA){6,}[GA]*(?:TC){6,}', 'TRX_5_3', 'Mirror repeat', 'Triplex', 24, 'triplex_potential', 0.85, 'Intermolecular triplex', 'Beal 1996'),
            (r'(?:GAA){4,}', 'TRX_5_4', 'GAA repeat', 'Sticky DNA', 12, 'sticky_dna_score', 0.95, 'Disease-associated repeats', 'Sakamoto 1999'),
            (r'(?:TTC){4,}', 'TRX_5_5', 'TTC repeat', 'Sticky DNA', 12, 'sticky_dna_score', 0.95, 'Disease-associated repeats', 'Sakamoto 1999'),
        ]
    }
    
    # Class 6: G-Quadruplex Family - Four-stranded G-rich structures (7 subclasses)
    G_QUADRUPLEX_PATTERNS = {
        'canonical_g4': [
            (r'G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'G4_6_1', 'Canonical G4', 'Canonical G4', 15, 'g4hunter_score', 0.95, 'Stable G4 structures', 'Burge 2006'),
            (r'G{4,}[ATGC]{1,5}G{4,}[ATGC]{1,5}G{4,}[ATGC]{1,5}G{4,}', 'G4_6_2', 'High-density G4', 'Canonical G4', 16, 'g4hunter_score', 0.98, 'Very stable G4', 'Todd 2005'),
        ],
        'relaxed_g4': [
            (r'G{2,}[ATGC]{1,12}G{2,}[ATGC]{1,12}G{2,}[ATGC]{1,12}G{2,}', 'G4_6_3', 'Relaxed G4', 'Relaxed G4', 12, 'g4hunter_score', 0.80, 'Potential G4 structures', 'Huppert 2005'),
            (r'G{3,}[ATGC]{8,15}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'G4_6_4', 'Long-loop G4', 'Relaxed G4', 18, 'g4hunter_score', 0.75, 'Alternative G4 topology', 'Phan 2006'),
        ],
        'bulged_g4': [
            (r'G{3,}[ATGC]{8,25}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'G4_6_5', 'Bulged G4', 'Bulged G4', 20, 'g4hunter_score', 0.85, 'G4 with bulge loops', 'Lim 2009'),
            (r'G{2,}[ATGC]{15,40}G{2,}[ATGC]{1,7}G{2,}[ATGC]{1,7}G{2,}', 'G4_6_6', 'Large bulge G4', 'Bulged G4', 25, 'g4hunter_score', 0.70, 'Extended bulge G4', 'Adrian 2014'),
        ],
        'bipartite_g4': [
            (r'G{2,}[ATGC]{15,70}G{2,}[ATGC]{1,7}G{2,}[ATGC]{1,7}G{2,}', 'G4_6_7', 'Bipartite G4', 'Bipartite G4', 30, 'g4hunter_score', 0.75, 'Two-block G4 structures', 'Guédin 2010'),
        ],
        'multimeric_g4': [
            (r'(?:G{3,}[ATGC]{1,7}){4,}G{3,}', 'G4_6_8', 'Multimeric G4', 'Multimeric G4', 25, 'g4hunter_score', 0.90, 'Multiple G4 units', 'Phan 2007'),
            (r'(?:G{2,}[ATGC]{1,10}){5,}G{2,}', 'G4_6_9', 'Extended multimeric G4', 'Multimeric G4', 30, 'g4hunter_score', 0.85, 'Long G4 arrays', 'Maizels 2006'),
        ],
        'imperfect_g4': [
            (r'G{2,}[ATGC]{1,10}[AG]G{1,3}[ATGC]{1,10}G{2,}[ATGC]{1,10}G{2,}', 'G4_6_10', 'Imperfect G4', 'Imperfect G4', 15, 'g4hunter_score', 0.65, 'G4-like with interruptions', 'Kuryavyi 2010'),
            (r'G{3,}[ATGC]{1,7}[AG]{1,2}G{2,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'G4_6_11', 'G-rich imperfect', 'Imperfect G4', 18, 'g4hunter_score', 0.70, 'Interrupted G-tracts', 'Webba da Silva 2007'),
        ],
        'g_triplex': [
            (r'G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}', 'G4_6_12', 'G-Triplex', 'G-Triplex intermediate', 12, 'g_triplex_score', 0.80, 'Three G-tract structures', 'Sen 1988'),
            (r'G{4,}[ATGC]{1,5}G{4,}[ATGC]{1,5}G{4,}', 'G4_6_13', 'High-density G-triplex', 'G-Triplex intermediate', 14, 'g_triplex_score', 0.85, 'Stable three-tract G-structure', 'Williamson 1989'),
        ]
    }
    
    # Class 7: i-Motif Family - C-rich structures (3 subclasses)
    I_MOTIF_PATTERNS = {
        'canonical_imotif': [
            (r'C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}', 'IM_7_1', 'Canonical i-motif', 'Canonical i-motif', 15, 'imotif_score', 0.95, 'pH-dependent C-rich structure', 'Gehring 1993'),
            (r'C{4,}[ATGC]{1,5}C{4,}[ATGC]{1,5}C{4,}[ATGC]{1,5}C{4,}', 'IM_7_2', 'High-density i-motif', 'Canonical i-motif', 16, 'imotif_score', 0.98, 'Stable i-motif', 'Leroy 1993'),
        ],
        'relaxed_imotif': [
            (r'C{2,}[ATGC]{1,12}C{2,}[ATGC]{1,12}C{2,}[ATGC]{1,12}C{2,}', 'IM_7_3', 'Relaxed i-motif', 'Relaxed i-motif', 12, 'imotif_score', 0.80, 'Potential i-motif structures', 'Mergny 1995'),
            (r'C{3,}[ATGC]{8,15}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}', 'IM_7_4', 'Long-loop i-motif', 'Relaxed i-motif', 18, 'imotif_score', 0.75, 'Alternative i-motif topology', 'Phan 2002'),
        ],
        'ac_motif': [
            (r'(?:AC){4,}|(?:CA){4,}', 'IM_7_5', 'AC-motif', 'AC-motif', 8, 'ac_motif_score', 0.85, 'AC alternating motif', 'Kang 2009'),
            (r'(?:A{2,3}C{2,3}){3,}', 'IM_7_6', 'Extended AC-motif', 'AC-motif', 12, 'ac_motif_score', 0.80, 'Variable AC motif', 'Zhou 2010'),
            (r'(?:ACG){3,}|(?:GCA){3,}', 'IM_7_7', 'ACG-motif', 'AC-motif', 9, 'ac_motif_score', 0.75, 'Trinucleotide AC motif', 'Liu 2012'),
        ]
    }
    
    # Class 8: Z-DNA - Left-handed double helix (2 subclasses)
    Z_DNA_PATTERNS = {
        'z_dna_canonical': [
            (r'(?:CG){4,}|(?:GC){4,}', 'ZDN_8_1', 'CG alternating', 'Z-DNA', 8, 'z_dna_score', 0.95, 'Classical Z-DNA sequence', 'Rich 1984'),
            (r'(?:CA){4,}(?:TG){4,}|(?:TG){4,}(?:CA){4,}', 'ZDN_8_2', 'CA-TG alternating', 'Z-DNA', 16, 'z_dna_score', 0.85, 'Alternative Z-DNA', 'Nordheim 1981'),
            (r'(?:AT){6,}|(?:TA){6,}', 'ZDN_8_3', 'AT alternating', 'Z-DNA', 12, 'z_dna_score', 0.75, 'AT-rich Z-DNA', 'Ellison 1986'),
        ],
        'egz_dna': [
            (r'[CG]{8,}', 'ZDN_8_4', 'CG-rich region', 'eGZ (Extruded-G) DNA', 8, 'egz_score', 0.90, 'Extruded-G Z-DNA', 'Herbert 1997'),
            (r'G{4,}C{4,}|C{4,}G{4,}', 'ZDN_8_5', 'GC clusters', 'eGZ (Extruded-G) DNA', 8, 'egz_score', 0.85, 'Clustered GC Z-DNA', 'Liu 1999'),
        ]
    }
    
    # Class 9: A-philic DNA - A-rich structural motifs
    A_PHILIC_PATTERNS = {
        'a_philic_tracts': [
            (r'A{6,}', 'APH_9_1', 'Poly-A tract', 'A-philic DNA', 6, 'a_philic_score', 0.95, 'A-rich structural element', 'Gorin 1995'),
            (r'(?:A{3,}[AT]){3,}', 'APH_9_2', 'A-rich region', 'A-philic DNA', 12, 'a_philic_score', 0.85, 'Mixed A-rich motif', 'Nelson 1987'),
            (r'(?:AA[AT]){4,}', 'APH_9_3', 'AAT/AAA motif', 'A-philic DNA', 12, 'a_philic_score', 0.80, 'Structured A-rich element', 'Crothers 1990'),
        ]
    }

    @classmethod
    def get_all_patterns(cls) -> Dict[str, Dict[str, List[Tuple]]]:
        """Return complete pattern registry"""
        return {
            'curved_dna': cls.CURVED_DNA_PATTERNS,
            'slipped_dna': cls.SLIPPED_DNA_PATTERNS,  
            'cruciform': cls.CRUCIFORM_PATTERNS,
            'r_loop': cls.R_LOOP_PATTERNS,
            'triplex': cls.TRIPLEX_PATTERNS,
            'g_quadruplex': cls.G_QUADRUPLEX_PATTERNS,
            'i_motif': cls.I_MOTIF_PATTERNS,
            'z_dna': cls.Z_DNA_PATTERNS,
            'a_philic': cls.A_PHILIC_PATTERNS
        }
    
    @classmethod
    def get_pattern_count(cls) -> Dict[str, int]:
        """Get pattern count statistics"""
        all_patterns = cls.get_all_patterns()
        counts = {}
        total = 0
        
        for motif_class, pattern_groups in all_patterns.items():
            class_count = sum(len(patterns) for patterns in pattern_groups.values())
            counts[motif_class] = class_count
            total += class_count
        
        counts['total'] = total
        return counts
    
    @classmethod 
    def get_subclass_mapping(cls) -> Dict[str, List[str]]:
        """Get mapping of classes to subclasses"""
        mapping = {
            'curved_dna': ['Global curvature', 'Local Curvature'],
            'slipped_dna': ['Direct Repeat', 'STR'],
            'cruciform': ['Inverted Repeats'],
            'r_loop': ['R-loop formation sites', 'QmRLFS-m1', 'QmRLFS-m2'],
            'triplex': ['Triplex', 'Sticky DNA'],
            'g_quadruplex': ['Multimeric G4', 'Canonical G4', 'Relaxed G4', 'Bulged G4', 
                           'Bipartite G4', 'Imperfect G4', 'G-Triplex intermediate'],
            'i_motif': ['Canonical i-motif', 'Relaxed i-motif', 'AC-motif'],
            'z_dna': ['Z-DNA', 'eGZ (Extruded-G) DNA'],
            'a_philic': ['A-philic DNA'],
            'hybrid': ['Dynamic overlaps'],
            'cluster': ['Dynamic clusters']
        }
        return mapping

# =============================================================================
# SCIENTIFIC SCORING ALGORITHMS
# =============================================================================

class MotifScoring:
    """Comprehensive scoring algorithms for all motif classes"""
    
    @staticmethod
    def g4hunter_score(sequence: str, window_size: int = 25) -> float:
        """
        G4Hunter algorithm for G-quadruplex scoring (Bedrat et al., 2016)
        
        Args:
            sequence: DNA sequence
            window_size: Sliding window size
            
        Returns:
            G4Hunter score (normalized)
        """
        if len(sequence) < window_size:
            window_size = len(sequence)
        
        scores = []
        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i + window_size]
            score = 0
            
            for base in window:
                if base == 'G':
                    score += 1
                elif base == 'C':
                    score -= 1
            
            scores.append(score)
        
        if not scores:
            return 0.0
        
        max_score = max(abs(s) for s in scores)
        return max_score / window_size
    
    @staticmethod
    def imotif_score(sequence: str) -> float:
        """
        i-motif scoring based on C-tract analysis (Zeraati et al., 2018)
        
        Args:
            sequence: DNA sequence
            
        Returns:
            i-motif formation score
        """
        if len(sequence) < 12:
            return 0.0
        
        # Count C-tracts
        c_tracts = re.findall(r'C{2,}', sequence)
        if len(c_tracts) < 3:
            return 0.0
        
        # Calculate score based on C-tract density and length
        total_c_length = sum(len(tract) for tract in c_tracts)
        c_density = total_c_length / len(sequence)
        tract_bonus = len(c_tracts) / 4  # Bonus for multiple tracts
        
        return min(c_density + tract_bonus * 0.2, 1.0)
    
    @staticmethod
    def z_dna_score(sequence: str) -> float:
        """
        Z-DNA scoring based on alternating purine-pyrimidine content (Ho et al., 1986)
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Z-DNA formation probability
        """
        if len(sequence) < 6:
            return 0.0
        
        alternating_score = 0
        total_pairs = 0
        
        for i in range(len(sequence) - 1):
            curr_base = sequence[i]
            next_base = sequence[i + 1]
            total_pairs += 1
            
            # Score alternating purine-pyrimidine pattern
            if ((curr_base in 'AG' and next_base in 'CT') or 
                (curr_base in 'CT' and next_base in 'AG')):
                alternating_score += 1
                
                # Bonus for CG steps (classical Z-DNA)
                if (curr_base == 'C' and next_base == 'G') or (curr_base == 'G' and next_base == 'C'):
                    alternating_score += 0.5
        
        return alternating_score / total_pairs if total_pairs > 0 else 0.0
    
    @staticmethod
    def curvature_score(sequence: str) -> float:
        """
        DNA curvature scoring based on A-tract analysis (Olson et al., 1998)
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Intrinsic curvature score
        """
        if len(sequence) < 4:
            return 0.0
        
        # Find A/T tracts
        a_tracts = re.findall(r'A{3,}', sequence)
        t_tracts = re.findall(r'T{3,}', sequence)
        
        # Calculate curvature based on tract length and frequency
        curvature = 0
        for tract in a_tracts + t_tracts:
            curvature += len(tract) ** 1.5  # Non-linear length dependency
        
        return min(curvature / len(sequence), 1.0)
    
    @staticmethod
    def triplex_potential(sequence: str) -> float:
        """
        Triplex formation potential (Frank-Kamenetskii & Mirkin, 1995)
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Triplex formation score
        """
        if len(sequence) < 15:
            return 0.0
        
        # Calculate homopurine and homopyrimidine content
        purine_runs = re.findall(r'[AG]{5,}', sequence)
        pyrimidine_runs = re.findall(r'[CT]{5,}', sequence)
        
        # Score based on run lengths and purity
        purine_score = sum(len(run) ** 1.2 for run in purine_runs)
        pyrimidine_score = sum(len(run) ** 1.2 for run in pyrimidine_runs)
        
        max_score = max(purine_score, pyrimidine_score)
        return min(max_score / (len(sequence) ** 1.2), 1.0)
    
    @staticmethod
    def r_loop_potential(sequence: str) -> float:
        """
        R-loop formation potential (Aguilera & García-Muse, 2012)
        
        Args:
            sequence: DNA sequence
            
        Returns:
            R-loop formation score
        """
        if len(sequence) < 20:
            return 0.0
        
        # GC skew calculation
        gc_skew = 0
        for base in sequence:
            if base == 'G':
                gc_skew += 1
            elif base == 'C':
                gc_skew -= 1
        
        # Normalize by sequence length
        gc_skew = abs(gc_skew) / len(sequence)
        
        # GC content (R-loops prefer GC-rich regions)
        gc_content = len(re.findall(r'[GC]', sequence)) / len(sequence)
        
        return min(gc_skew * 0.6 + gc_content * 0.4, 1.0)
    
    @staticmethod
    def qmrlfs_score(sequence: str) -> float:
        """
        QmRLFS-based R-loop formation scoring (Jenjaroenpun & Wongsurawat, 2016)
        
        Implements the QmRLFS algorithm for R-loop forming sequence detection
        based on RIZ (R-loop Initiating Zone) and REZ (R-loop Extending Zone)
        
        Args:
            sequence: DNA sequence
            
        Returns:
            QmRLFS score (0.0-1.0)
        """
        try:
            from qmrlfs_finder import QmRLFSDetector
            
            if len(sequence) < 20:
                return 0.0
            
            # Use quick mode for scoring to avoid performance issues
            detector = QmRLFSDetector(quick_mode=True)
            results = detector.analyze_sequence(sequence, analyze_both_strands=False)
            
            if not results:
                return 0.0
            
            # Return the highest scoring RLFS
            max_score = max(result["qmrlfs_score"] for result in results)
            return max_score
            
        except ImportError:
            # Fallback to simple scoring if QmRLFS module not available
            return MotifScoring.r_loop_potential(sequence)
    
    @staticmethod
    def instability_score(sequence: str) -> float:
        """
        Repeat instability scoring for slipped DNA structures
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Instability score based on repeat characteristics
        """
        if len(sequence) < 6:
            return 0.0
        
        # Find repeating units
        max_instability = 0
        
        for unit_length in range(1, min(7, len(sequence) // 3)):
            for i in range(len(sequence) - unit_length * 2):
                unit = sequence[i:i + unit_length]
                count = 1
                
                # Count consecutive repeats
                pos = i + unit_length
                while pos + unit_length <= len(sequence) and sequence[pos:pos + unit_length] == unit:
                    count += 1
                    pos += unit_length
                
                if count >= 3:  # At least 3 repeats
                    instability = count * (unit_length ** 0.5)
                    max_instability = max(max_instability, instability)
        
        return min(max_instability / 10, 1.0)
    
    @staticmethod
    def cruciform_stability(sequence: str) -> float:
        """
        Cruciform stability based on palindrome characteristics
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Stability score for cruciform formation
        """
        if len(sequence) < 8:
            return 0.0
        
        def reverse_complement(seq):
            complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            return ''.join(complement.get(base, base) for base in reversed(seq))
        
        max_palindrome = 0
        
        # Look for palindromic regions
        for i in range(len(sequence)):
            for j in range(i + 8, len(sequence) + 1):
                subseq = sequence[i:j]
                if subseq == reverse_complement(subseq):
                    palindrome_length = len(subseq)
                    stability = palindrome_length ** 1.5 / len(sequence)
                    max_palindrome = max(max_palindrome, stability)
        
        return min(max_palindrome, 1.0)
    
    @staticmethod
    def a_philic_score(sequence: str) -> float:
        """
        A-philic DNA scoring based on A-tract characteristics
        
        Args:
            sequence: DNA sequence
            
        Returns:
            A-philic propensity score
        """
        if len(sequence) < 6:
            return 0.0
        
        a_content = len(re.findall(r'A', sequence)) / len(sequence)
        a_tracts = re.findall(r'A{3,}', sequence)
        
        # Score based on A content and tract formation
        tract_bonus = sum(len(tract) ** 1.2 for tract in a_tracts) / len(sequence)
        
        return min(a_content * 0.7 + tract_bonus * 0.3, 1.0)

# =============================================================================
# HYPERSCAN INTEGRATION & OPTIMIZATION
# =============================================================================

class HyperscanManager:
    """Hyperscan database management for high-performance pattern matching"""
    
    def __init__(self):
        self.compiled_db = None
        self.pattern_info = {}
        self.hyperscan_available = HYPERSCAN_AVAILABLE
    
    def compile_database(self, patterns: List[Tuple[str, str]]) -> bool:
        """
        Compile Hyperscan database from patterns
        
        Args:
            patterns: List of (pattern, identifier) tuples
            
        Returns:
            True if compilation successful
        """
        if not self.hyperscan_available:
            return False
        
        try:
            # Prepare patterns for Hyperscan
            hyperscan_patterns = []
            for i, (pattern, pattern_id) in enumerate(patterns):
                hyperscan_patterns.append((pattern.encode(), i, hyperscan.HS_FLAG_CASELESS))
                self.pattern_info[i] = pattern_id
            
            # Compile database
            self.compiled_db = hyperscan.hs_compile_multi(hyperscan_patterns)
            return True
            
        except Exception as e:
            print(f"Hyperscan compilation failed: {e}")
            return False
    
    def scan_sequence(self, sequence: str) -> List[Tuple[int, int, str]]:
        """
        Scan sequence with compiled Hyperscan database
        
        Args:
            sequence: DNA sequence to scan
            
        Returns:
            List of (start, end, pattern_id) matches
        """
        if not self.compiled_db:
            return []
        
        matches = []
        
        def match_handler(pattern_id: int, start: int, end: int, flags: int, context=None):
            pattern_info = self.pattern_info.get(pattern_id, f'pattern_{pattern_id}')
            matches.append((start, end, pattern_info))
        
        try:
            hyperscan.hs_scan(self.compiled_db, sequence.encode(), match_handler, None)
        except Exception as e:
            print(f"Hyperscan scanning failed: {e}")
        
        return matches

# =============================================================================
# PATTERN VALIDATION & TESTING
# =============================================================================

def validate_all_patterns() -> Dict[str, Any]:
    """
    Validate all patterns for correctness and Hyperscan compatibility
    
    Returns:
        Validation results dictionary
    """
    results = {
        'total_patterns': 0,
        'valid_patterns': 0,
        'invalid_patterns': [],
        'hyperscan_compatible': 0,
        'pattern_counts': {},
        'validation_passed': True
    }
    
    all_patterns = PatternRegistry.get_all_patterns()
    
    for motif_class, pattern_groups in all_patterns.items():
        class_count = 0
        for pattern_group, patterns in pattern_groups.items():
            for pattern_tuple in patterns:
                if len(pattern_tuple) >= 9:  # Full tuple
                    pattern = pattern_tuple[0]
                    pattern_id = pattern_tuple[1]
                    
                    results['total_patterns'] += 1
                    class_count += 1
                    
                    # Test regex compilation
                    try:
                        re.compile(pattern, re.IGNORECASE)
                        results['valid_patterns'] += 1
                        
                        # Test Hyperscan compatibility
                        if not any(incompatible in pattern for incompatible in ['\\b', '\\B', '(?=', '(?!', '(?<=', '(?<!', '\\1', '\\2']):
                            results['hyperscan_compatible'] += 1
                        
                    except re.error as e:
                        results['invalid_patterns'].append((motif_class, pattern_id, str(e)))
                        results['validation_passed'] = False
        
        results['pattern_counts'][motif_class] = class_count
    
    return results

def run_pattern_tests() -> bool:
    """Run comprehensive pattern tests"""
    print("Running Non-B DNA Pattern Validation...")
    
    # Test sequences for each class
    test_sequences = {
        'curved_dna': 'AAAAAAAATTTTTTTAAAAATTTT',  # A/T tracts
        'slipped_dna': 'CACACACACACACACACACA',     # CA repeats
        'cruciform': 'ATGCATGCATGCATGC',          # Palindrome
        'r_loop': 'GGGCCCGGGATGCCCGGG',           # GC-rich R-loop site
        'triplex': 'GAAAGAAAGAAAGAAAGAAA',        # Homopurine tract
        'g_quadruplex': 'GGGTTAGGGTTAGGGTTAGGG',   # Canonical G4
        'i_motif': 'CCCTAACCCTAACCCTAACCC',       # Canonical i-motif
        'z_dna': 'CGCGCGCGCGCGCGCGCG',            # CG alternating
        'a_philic': 'AAAAAAAAAAAAAAAAA'           # Poly-A
    }
    
    validation_results = validate_all_patterns()
    print(f"Pattern validation: {validation_results['valid_patterns']}/{validation_results['total_patterns']} patterns valid")
    print(f"Hyperscan compatible: {validation_results['hyperscan_compatible']}/{validation_results['total_patterns']} patterns")
    
    # Test pattern matching
    all_patterns = PatternRegistry.get_all_patterns()
    # Instantiate scoring system to validate initialization
    scoring = MotifScoring()
    
    test_results = {}
    for motif_class, test_seq in test_sequences.items():
        if motif_class in all_patterns:
            matches_found = 0
            
            for pattern_group, patterns in all_patterns[motif_class].items():
                for pattern_tuple in patterns:
                    if len(pattern_tuple) >= 3:
                        pattern = pattern_tuple[0]
                        try:
                            compiled_pattern = re.compile(pattern, re.IGNORECASE)
                            if compiled_pattern.search(test_seq):
                                matches_found += 1
                                break
                        except:
                            continue
                            
                if matches_found > 0:
                    break
            
            test_results[motif_class] = matches_found > 0
            print(f"[OK] {motif_class}: {'PASS' if matches_found > 0 else 'FAIL'}")
    
    all_passed = all(test_results.values()) and validation_results['validation_passed']
    print(f"\nOverall validation: {'PASSED' if all_passed else 'FAILED'}")
    
    return all_passed

# =============================================================================
# HYPERSCAN REGISTRY LOADER INTEGRATION
# =============================================================================

try:
    # from .load_hsdb import load_db_for_class
    _LOAD_HSDB_AVAILABLE = True
except ImportError:
    try:
        # Fallback: direct import
        import sys
        import os as _os
        sys.path.insert(0, _os.path.dirname(__file__))
        from load_hsdb import load_db_for_class
        _LOAD_HSDB_AVAILABLE = True
    except ImportError:
        _LOAD_HSDB_AVAILABLE = False
        load_db_for_class = None

# In-memory cache to avoid repeated compiles/deserializes
_HS_DB_CACHE = {}
_REGISTRY_CACHE = {}


def get_pattern_registry(class_name: str, registry_dir: str = "registry"):
    """
    Returns parsed registry dict (as saved by generator): contains 'patterns' list etc.
    Caches the result.
    """
    if not _LOAD_HSDB_AVAILABLE:
        raise ImportError("load_hsdb module not available")
    
    key = f"{registry_dir}/{class_name}"
    if key in _REGISTRY_CACHE:
        return _REGISTRY_CACHE[key]
    
    # Read registry directly from file
    import pickle
    import json
    pkl_path = os.path.join(registry_dir, f"{class_name}_registry.pkl")
    json_path = os.path.join(registry_dir, f"{class_name}_registry.json")
    
    if os.path.isfile(pkl_path):
        with open(pkl_path, "rb") as fh:
            full = pickle.load(fh)
    elif os.path.isfile(json_path):
        with open(json_path, "r") as fh:
            full = json.load(fh)
    else:
        raise FileNotFoundError(f"No registry found for {class_name} in {registry_dir}")
    
    _REGISTRY_CACHE[key] = full
    return full


def get_hs_db_for_class(class_name: str, registry_dir: str = "registry"):
    """
    Return (db, id_to_pattern, id_to_score). Caches DB in process memory.
    db is a hyperscan.Database instance (or None if hyperscan not available).
    id_to_pattern: dict mapping id -> pattern string (tenmer for 10-mer, regex for others)
    """
    if not _LOAD_HSDB_AVAILABLE:
        raise ImportError("load_hsdb module not available")
    
    key = f"{registry_dir}/{class_name}"
    if key in _HS_DB_CACHE:
        return _HS_DB_CACHE[key]
    
    db, id_to_pattern, id_to_score = load_db_for_class(class_name, registry_dir)
    _HS_DB_CACHE[key] = (db, id_to_pattern, id_to_score)
    return db, id_to_pattern, id_to_score


# =============================================================================
# PATTERN STATISTICS & INFORMATION
# =============================================================================

def get_pattern_statistics() -> Dict[str, Any]:
    """Get comprehensive pattern statistics"""
    counts = PatternRegistry.get_pattern_count()
    subclasses = PatternRegistry.get_subclass_mapping()
    
    stats = {
        'total_patterns': counts['total'],
        'total_classes': 11,
        'total_subclasses': sum(len(subs) for subs in subclasses.values()),
        'class_breakdown': {k: v for k, v in counts.items() if k != 'total'},
        'subclass_breakdown': subclasses,
        'scoring_methods': [
            'g4hunter_score', 'imotif_score', 'z_dna_score', 'curvature_score',
            'triplex_potential', 'r_loop_potential', 'qmrlfs_score', 'instability_score', 
            'cruciform_stability', 'a_philic_score'
        ]
    }
    
    return stats

if __name__ == "__main__":
    # Run validation and display statistics
    print("="*60)
    print("NON-B DNA MOTIF PATTERNS REGISTRY")
    print("="*60)
    
    success = run_pattern_tests()
    
    if success:
        stats = get_pattern_statistics()
        print(f"\nPattern Registry Statistics:")
        print(f"Total Patterns: {stats['total_patterns']}")
        print(f"Total Classes: {stats['total_classes']}")
        print(f"Total Subclasses: {stats['total_subclasses']}")
        
        print("\nClass Breakdown:")
        for class_name, count in stats['class_breakdown'].items():
            print(f"  {class_name:<15}: {count:>3} patterns")
        
        print("\nScoring Methods Available:")
        for method in stats['scoring_methods']:
            print(f"  • {method}")
            
    print("="*60)


# ============================================================================
# NBDScanner Utilities - Helper Functions & I/O Operations
# ============================================================================
"""
NBDScanner Utilities - Helper Functions & I/O Operations
========================================================

Consolidated utility functions for sequence processing, I/O operations,
statistical calculations, and data formatting for the NBDScanner system.

UTILITY FUNCTIONS TABLE:
========================
Category          | Functions                    | Description
------------------|------------------------------|----------------------------------
Sequence I/O      | parse_fasta, write_fasta    | FASTA file handling
Sequence Utils    | reverse_complement, gc_content | Basic sequence operations  
Statistics        | get_basic_stats, motif_stats | Sequence and motif statistics
Formatting        | wrap, format_results        | Output formatting
Validation        | validate_sequence, quality_check | Input validation
Export            | to_bed, to_csv, to_json     | Data export formats

Author: Dr. Venkata Rajesh Yella
License: MIT
Version: 2024.1
"""

import re
import os
import json
import csv
import random
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Tuple
from collections import Counter, defaultdict
from io import StringIO
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# SEQUENCE I/O OPERATIONS
# =============================================================================

def parse_fasta(fasta_content: str) -> Dict[str, str]:
    """
    Parse FASTA format content into sequences dictionary
    
    Args:
        fasta_content: FASTA format string content
        
    Returns:
        Dictionary of {sequence_name: sequence}
    """
    sequences = {}
    current_name = None
    current_seq = []
    
    lines = fasta_content.strip().split('\n')
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        if line.startswith('>'):
            # Save previous sequence
            if current_name and current_seq:
                sequences[current_name] = ''.join(current_seq)
            
            # Start new sequence
            current_name = line[1:].strip()
            if not current_name:
                current_name = f"sequence_{len(sequences) + 1}"
            current_seq = []
        else:
            # Add to current sequence
            current_seq.append(line.upper())
    
    # Save last sequence
    if current_name and current_seq:
        sequences[current_name] = ''.join(current_seq)
    
    return sequences

def parse_fasta_chunked(file_object, chunk_size_mb: int = 2):
    """
    Memory-efficient FASTA parser using chunked reading for large files.
    Yields (name, sequence) tuples one at a time to avoid loading entire file in memory.
    
    This version is optimized to use less memory than loading the entire file at once
    by processing sequences one at a time and yielding them immediately.
    
    Args:
        file_object: File-like object (from st.file_uploader or open())
        chunk_size_mb: Size of read chunks in MB (default: 2MB - optimized for memory efficiency)
        
    Yields:
        Tuple of (sequence_name, sequence_string)
    """
    chunk_size = chunk_size_mb * 1024 * 1024  # Convert to bytes
    current_name = None
    current_seq_parts = []  # Use list to accumulate, join only when yielding
    
    # Read file in smaller chunks to avoid memory overflow
    buffer = ""
    while True:
        chunk = file_object.read(chunk_size)
        if not chunk:
            break
        
        # Decode if bytes
        if isinstance(chunk, bytes):
            chunk = chunk.decode('utf-8', errors='ignore')
        
        buffer += chunk
        lines = buffer.split('\n')
        
        # Keep the last incomplete line in buffer
        buffer = lines[-1]
        lines = lines[:-1]
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                # Yield previous sequence
                if current_name and current_seq_parts:
                    # Join and yield, then clear parts to free memory immediately
                    full_seq = ''.join(current_seq_parts)
                    current_seq_parts.clear()  # Free memory immediately
                    yield (current_name, full_seq)
                    del full_seq  # Help garbage collector
                
                # Start new sequence
                current_name = line[1:].strip()
                if not current_name:
                    current_name = f"sequence_{id(line)}"
                current_seq_parts = []
            else:
                # Add to current sequence
                current_seq_parts.append(line.upper())
    
    # Process remaining buffer
    if buffer.strip() and not buffer.startswith('>'):
        current_seq_parts.append(buffer.strip().upper())
    
    # Yield last sequence
    if current_name and current_seq_parts:
        full_seq = ''.join(current_seq_parts)
        current_seq_parts.clear()
        yield (current_name, full_seq)
        del full_seq

def get_file_preview(file_object, max_sequences: int = 3, max_preview_chars: int = 400):
    """
    Get a preview of FASTA file content without loading entire file.
    
    Args:
        file_object: File-like object from st.file_uploader
        max_sequences: Maximum number of sequences to preview (default: 3)
        max_preview_chars: Maximum characters to show per sequence (default: 400)
        
    Returns:
        Dict with keys: 'num_sequences', 'total_bp', 'previews' (list of dicts)
    """
    file_object.seek(0)  # Reset to beginning
    
    previews = []
    total_bp = 0
    num_sequences = 0
    
    for name, seq in parse_fasta_chunked(file_object):
        num_sequences += 1
        seq_len = len(seq)
        total_bp += seq_len
        
        if len(previews) < max_sequences:
            preview_seq = seq[:max_preview_chars]
            if len(seq) > max_preview_chars:
                preview_seq += "..."
            
            previews.append({
                'name': name,
                'length': seq_len,
                'preview': preview_seq
            })
    
    file_object.seek(0)  # Reset for subsequent use
    
    return {
        'num_sequences': num_sequences,
        'total_bp': total_bp,
        'previews': previews
    }

def write_fasta(sequences: Dict[str, str], filename: str) -> bool:
    """
    Write sequences to FASTA format file
    
    Args:
        sequences: Dictionary of {name: sequence}
        filename: Output filename
        
    Returns:
        True if successful
    """
    try:
        with open(filename, 'w') as f:
            for name, seq in sequences.items():
                f.write(f">{name}\n")
                # Write sequence in 80-character lines
                for i in range(0, len(seq), 80):
                    f.write(f"{seq[i:i+80]}\n")
        return True
    except Exception as e:
        print(f"Error writing FASTA file {filename}: {e}")
        return False

def read_fasta_file(filename: str) -> Dict[str, str]:
    """
    Read FASTA file and return sequences dictionary
    
    Args:
        filename: Path to FASTA file
        
    Returns:
        Dictionary of {sequence_name: sequence}
    """
    try:
        with open(filename, 'r') as f:
            content = f.read()
        return parse_fasta(content)
    except Exception as e:
        print(f"Error reading FASTA file {filename}: {e}")
        return {}

# =============================================================================
# SEQUENCE MANIPULATION & VALIDATION
# =============================================================================

def validate_sequence(sequence: str) -> Tuple[bool, str]:
    """
    Validate DNA sequence format and content
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        (is_valid, error_message)
    """
    if not sequence:
        return False, "Empty sequence"
    
    if not isinstance(sequence, str):
        return False, "Sequence must be string"
    
    # Check for valid DNA characters
    valid_chars = set('ATGCRYSWKMBDHVN-')  # Include ambiguous bases
    invalid_chars = set(sequence.upper()) - valid_chars
    
    if invalid_chars:
        return False, f"Invalid characters found: {invalid_chars}"
    
    if len(sequence) < 10:
        return False, "Sequence too short (minimum 10 bp)"
    
    return True, "Valid sequence"

def reverse_complement(sequence: str) -> str:
    """
    Generate reverse complement of DNA sequence
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        Reverse complement sequence
    """
    complement_map = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'R': 'Y', 'Y': 'R', 'S': 'W', 'W': 'S',
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
        'D': 'H', 'H': 'D', 'N': 'N', '-': '-'
    }
    
    return ''.join(complement_map.get(base.upper(), base) for base in reversed(sequence))

def gc_content(sequence: str) -> float:
    """
    Calculate GC content of sequence
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        GC content as percentage (0-100)
    """
    if not sequence:
        return 0.0
    
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100

def at_content(sequence: str) -> float:
    """
    Calculate AT content of sequence
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        AT content as percentage (0-100)
    """
    return 100.0 - gc_content(sequence)

def is_palindrome(sequence: str) -> bool:
    """
    Check if sequence is palindromic (same as reverse complement)
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        True if palindromic
    """
    return sequence.upper() == reverse_complement(sequence).upper()

def calculate_tm(sequence: str) -> float:
    """
    Calculate melting temperature using nearest neighbor method (simplified)
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        Melting temperature in Celsius
    """
    if len(sequence) < 2:
        return 0.0
    
    # Simplified Tm calculation
    gc = gc_content(sequence)
    length = len(sequence)
    
    if length <= 13:
        # For short sequences
        tm = (sequence.upper().count('A') + sequence.upper().count('T')) * 2 + \
             (sequence.upper().count('G') + sequence.upper().count('C')) * 4
    else:
        # For longer sequences
        tm = 64.9 + 41 * (gc / 100) - 650 / length
    
    return tm

def shuffle_sequence(sequence: str, seed: Optional[int] = None) -> str:
    """
    Generate shuffled version of sequence (preserving composition)
    
    Args:
        sequence: DNA sequence string
        seed: Random seed for reproducibility
        
    Returns:
        Shuffled sequence
    """
    import random
    if seed is not None:
        random.seed(seed)
    
    seq_list = list(sequence.upper())
    random.shuffle(seq_list)
    return ''.join(seq_list)

# =============================================================================
# STATISTICS & ANALYSIS FUNCTIONS
# =============================================================================

def get_basic_stats(sequence: str, motifs: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
    """
    Calculate basic sequence statistics
    
    Args:
        sequence: DNA sequence string
        motifs: Optional list of detected motifs
        
    Returns:
        Dictionary of statistics
    """
    if not sequence:
        return {}
    
    seq = sequence.upper()
    length = len(seq)
    
    # Base composition
    base_counts = Counter(seq)
    
    stats = {
        'Length': length,
        'A': base_counts.get('A', 0),
        'T': base_counts.get('T', 0),
        'G': base_counts.get('G', 0),
        'C': base_counts.get('C', 0),
        'N': base_counts.get('N', 0),
        'GC%': round(gc_content(seq), 2),
        'AT%': round(at_content(seq), 2),
        'Tm': round(calculate_tm(seq), 1)
    }
    
    # Motif statistics if provided
    if motifs:
        stats.update(calculate_motif_statistics(motifs, length))
    
    return stats

def calculate_motif_statistics(motifs: List[Dict[str, Any]], sequence_length: int) -> Dict[str, Any]:
    """
    Calculate comprehensive motif statistics
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of analyzed sequence
        
    Returns:
        Dictionary of motif statistics
    """
    if not motifs:
        return {
            'Total_Motifs': 0,
            'Coverage%': 0.0,
            'Density': 0.0,
            'Classes_Detected': 0,
            'Subclasses_Detected': 0
        }
    
    # Count by class and subclass
    class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    subclass_counts = Counter(m.get('Subclass', 'Unknown') for m in motifs)
    
    # Calculate coverage
    # COORDINATE SYSTEM: Motifs use 1-based INCLUSIVE coordinates
    # Example: Start=1, End=20 means positions 1-20 inclusive (20 bases)
    # Convert to 0-based half-open for Python range(): range(0, 20)
    covered_positions = set()
    for motif in motifs:
        start = motif.get('Start', 0) - 1  # Convert 1-based to 0-based
        end = motif.get('End', 0)  # Inclusive end becomes exclusive in range()
        covered_positions.update(range(start, end))
    
    coverage_percent = (len(covered_positions) / sequence_length * 100) if sequence_length > 0 else 0
    density = len(motifs) / (sequence_length / 1000) if sequence_length > 0 else 0  # Motifs per kb
    
    stats = {
        'Total_Motifs': len(motifs),
        'Coverage%': round(coverage_percent, 2),
        'Density': round(density, 2),
        'Classes_Detected': len(class_counts),
        'Subclasses_Detected': len(subclass_counts),
        'Class_Distribution': dict(class_counts),
        'Subclass_Distribution': dict(subclass_counts)
    }
    
    # Score statistics
    scores = [m.get('Score', 0) for m in motifs if isinstance(m.get('Score'), (int, float))]
    if scores:
        stats.update({
            'Score_Mean': round(np.mean(scores), 3),
            'Score_Std': round(np.std(scores), 3),
            'Score_Min': round(min(scores), 3),
            'Score_Max': round(max(scores), 3)
        })
    
    # Length statistics
    lengths = [m.get('Length', 0) for m in motifs if isinstance(m.get('Length'), int)]
    if lengths:
        stats.update({
            'Length_Mean': round(np.mean(lengths), 1),
            'Length_Std': round(np.std(lengths), 1),
            'Length_Min': min(lengths),
            'Length_Max': max(lengths)
        })
    
    return stats


def analyze_class_subclass_detection(motifs: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Analyze which classes and subclasses were detected and which were not.
    Provides comprehensive report on all 11 Non-B DNA classes.
    
    Args:
        motifs: List of detected motif dictionaries
        
    Returns:
        Dictionary with detailed class/subclass detection analysis
    """
    # Define all expected Non-B DNA classes and their subclasses
    all_classes = {
        'Curved DNA': ['Global Curvature', 'Local Curvature'],
        'Slipped DNA': ['Direct Repeat', 'STR'],
        'Cruciform DNA': ['Inverted Repeats'],
        'R-loop': ['R-loop formation sites', 'QmRLFS-m1', 'QmRLFS-m2'],
        'Triplex': ['Triplex', 'Sticky DNA'],
        'G-Quadruplex': [
            'Multimeric G4', 'Canonical G4', 'Relaxed G4', 'Bulged G4',
            'Bipartite G4', 'Imperfect G4', 'G-Triplex intermediate', 'Long-loop G4'
        ],
        'i-Motif': ['Canonical i-motif', 'Relaxed i-motif', 'AC-motif'],
        'Z-DNA': ['Z-DNA', 'eGZ'],
        'A-philic DNA': ['A-philic DNA'],
        'Hybrid': ['Dynamic overlaps'],
        'Non-B_DNA_Clusters': ['Dynamic clusters']
    }
    
    # Count detected classes and subclasses
    detected_classes = defaultdict(int)
    detected_subclasses = defaultdict(lambda: defaultdict(int))
    
    for motif in motifs:
        cls = motif.get('Class', 'Unknown')
        subcls = motif.get('Subclass', 'Unknown')
        detected_classes[cls] += 1
        detected_subclasses[cls][subcls] += 1
    
    # Analyze detection status
    detection_report = {
        'total_classes': len(all_classes),
        'detected_classes': len(detected_classes),
        'not_detected_classes': [],
        'class_details': {},
        'summary': {}
    }
    
    for cls, expected_subclasses in all_classes.items():
        if cls in detected_classes:
            # Class was detected
            detected_subs = list(detected_subclasses[cls].keys())
            not_detected_subs = [sub for sub in expected_subclasses 
                                if sub not in detected_subs and sub not in ['Dynamic overlaps', 'Dynamic clusters']]
            
            detection_report['class_details'][cls] = {
                'status': 'DETECTED',
                'count': detected_classes[cls],
                'expected_subclasses': expected_subclasses,
                'detected_subclasses': detected_subs,
                'not_detected_subclasses': not_detected_subs,
                'subclass_counts': dict(detected_subclasses[cls])
            }
        else:
            # Class was not detected
            detection_report['not_detected_classes'].append(cls)
            detection_report['class_details'][cls] = {
                'status': 'NOT_DETECTED',
                'count': 0,
                'expected_subclasses': expected_subclasses,
                'detected_subclasses': [],
                'not_detected_subclasses': expected_subclasses,
                'subclass_counts': {}
            }
    
    # Create summary
    detection_report['summary'] = {
        'Total Classes': len(all_classes),
        'Detected Classes': len(detected_classes),
        'Not Detected Classes': len(detection_report['not_detected_classes']),
        'Total Motifs': len(motifs),
        'Detection Rate': f"{len(detected_classes) / len(all_classes) * 100:.1f}%"
    }
    
    return detection_report


def print_detection_report(detection_report: Dict[str, Any]) -> str:
    """
    Format detection report as readable text
    
    Args:
        detection_report: Report from analyze_class_subclass_detection()
        
    Returns:
        Formatted text report
    """
    lines = []
    lines.append("="*80)
    lines.append("NON-B DNA MOTIF DETECTION ANALYSIS REPORT")
    lines.append("="*80)
    lines.append("")
    
    # Summary
    lines.append("SUMMARY:")
    for key, value in detection_report['summary'].items():
        lines.append(f"  {key}: {value}")
    lines.append("")
    
    # Detected classes
    lines.append("DETECTED CLASSES:")
    lines.append("-"*80)
    for cls, details in sorted(detection_report['class_details'].items()):
        if details['status'] == 'DETECTED':
            lines.append(f"\n{cls} ({details['count']} motifs):")
            lines.append(f"  Expected subclasses: {len(details['expected_subclasses'])}")
            lines.append(f"  Detected subclasses: {len(details['detected_subclasses'])}")
            
            if details['detected_subclasses']:
                lines.append("  [OK] Detected:")
                for subcls in details['detected_subclasses']:
                    count = details['subclass_counts'].get(subcls, 0)
                    lines.append(f"    - {subcls} ({count} motifs)")
            
            if details['not_detected_subclasses']:
                lines.append("  ✗ Not Detected:")
                for subcls in details['not_detected_subclasses']:
                    lines.append(f"    - {subcls}")
    
    lines.append("")
    
    # Not detected classes
    if detection_report['not_detected_classes']:
        lines.append("NOT DETECTED CLASSES:")
        lines.append("-"*80)
        for cls in detection_report['not_detected_classes']:
            details = detection_report['class_details'][cls]
            lines.append(f"\n✗ {cls}:")
            lines.append(f"  Expected subclasses: {', '.join(details['expected_subclasses'])}")
            lines.append(f"  Reason: No motifs of this class were found in the sequence")
    
    lines.append("")
    lines.append("="*80)
    
    return '\n'.join(lines)


# =============================================================================
# FORMATTING & OUTPUT UTILITIES
# =============================================================================

def wrap(sequence: str, width: int = 80) -> str:
    """
    Wrap sequence to specified width
    
    Args:
        sequence: DNA sequence string
        width: Line width for wrapping
        
    Returns:
        Wrapped sequence string
    """
    if not sequence or width <= 0:
        return sequence
    
    return '\n'.join(sequence[i:i+width] for i in range(0, len(sequence), width))

def format_motif_rows(motifs: List[Dict[str, Any]]) -> List[List[str]]:
    """
    Format motifs for tabular display
    
    Args:
        motifs: List of motif dictionaries
        
    Returns:
        List of formatted rows
    """
    if not motifs:
        return []
    
    rows = []
    headers = ['Class', 'Subclass', 'Start', 'End', 'Length', 'Sequence', 'Score', 'Strand']
    
    for motif in motifs:
        row = []
        for header in headers:
            value = motif.get(header, 'N/A')
            if header == 'Sequence' and isinstance(value, str) and len(value) > 50:
                value = value[:47] + '...'
            elif header == 'Score' and isinstance(value, (int, float)):
                value = f"{value:.3f}"
            row.append(str(value))
        rows.append(row)
    
    return rows

def create_summary_table(sequences: Dict[str, str], results: Dict[str, List[Dict[str, Any]]]) -> pd.DataFrame:
    """
    Create summary table for multiple sequence analysis
    
    Args:
        sequences: Dictionary of {name: sequence}
        results: Dictionary of {name: motifs_list}
        
    Returns:
        Summary DataFrame
    """
    summary_data = []
    
    for name, sequence in sequences.items():
        motifs = results.get(name, [])
        stats = get_basic_stats(sequence, motifs)
        
        summary_data.append({
            'Sequence_Name': name,
            'Length_bp': stats.get('Length', 0),
            'GC_Content': stats.get('GC%', 0),
            'Total_Motifs': stats.get('Total_Motifs', 0),
            'Coverage_Percent': stats.get('Coverage%', 0),
            'Motif_Density': stats.get('Density', 0),
            'Classes_Detected': stats.get('Classes_Detected', 0),
            'Subclasses_Detected': stats.get('Subclasses_Detected', 0)
        })
    
    return pd.DataFrame(summary_data)

# =============================================================================
# DATA EXPORT FUNCTIONS
# =============================================================================

def export_to_bed(motifs: List[Dict[str, Any]], sequence_name: str = "sequence", 
                  filename: Optional[str] = None) -> str:
    """
    Export motifs to BED format
    
    Args:
        motifs: List of motif dictionaries
        sequence_name: Name of the sequence
        filename: Optional output filename
        
    Returns:
        BED format string
    """
    bed_lines = []
    bed_lines.append("track name=NBDScanner_motifs description=\"Non-B DNA motifs\" itemRgb=On")
    
    # Color mapping for different classes
    class_colors = {
        'Curved_DNA': '255,182,193',      # Light pink
        'Slipped_DNA': '255,218,185',     # Peach
        'Cruciform': '173,216,230',       # Light blue
        'R-Loop': '144,238,144',          # Light green
        'Triplex': '221,160,221',         # Plum
        'G-Quadruplex': '255,215,0',      # Gold
        'i-Motif': '255,165,0',           # Orange
        'Z-DNA': '138,43,226',            # Blue violet
        'A-philic_DNA': '230,230,250',    # Lavender
        'Hybrid': '192,192,192',          # Silver
        'Non-B_DNA_Clusters': '128,128,128'  # Gray
    }
    
    for motif in motifs:
        chrom = sequence_name
        start = max(0, motif.get('Start', 1) - 1)  # Convert to 0-based
        end = motif.get('End', start + 1)
        name = f"{motif.get('Class', 'Unknown')}_{motif.get('Subclass', 'Unknown')}"
        score = int(min(1000, max(0, motif.get('Score', 0) * 1000)))  # Scale to 0-1000
        strand = motif.get('Strand', '+')
        color = class_colors.get(motif.get('Class'), '128,128,128')
        
        bed_line = f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{start}\t{end}\t{color}"
        bed_lines.append(bed_line)
    
    bed_content = '\n'.join(bed_lines)
    
    if filename:
        try:
            with open(filename, 'w') as f:
                f.write(bed_content)
        except Exception as e:
            print(f"Error writing BED file {filename}: {e}")
    
    return bed_content

def export_to_csv(motifs: List[Dict[str, Any]], filename: Optional[str] = None) -> str:
    """
    Export motifs to CSV format with CORE fields only (as per requirements).
    
    Output tables should only have: Sequence_Name, Source, Class, Subclass, 
    Start, End, Length, Sequence, Score
    
    Additional detailed columns are only shown in Excel download per-motif sheets.
    
    Args:
        motifs: List of motif dictionaries
        filename: Optional output filename
        
    Returns:
        CSV format string
    """
    if not motifs:
        return "No motifs to export"
    
    # Use core columns constant
    core_columns = CORE_OUTPUT_COLUMNS
    
    output = StringIO()
    writer = csv.DictWriter(output, fieldnames=core_columns)
    writer.writeheader()
    
    for motif in motifs:
        # Create row with only core fields
        row = {}
        for col in core_columns:
            value = motif.get(col, 'NA')
            
            # If empty, set to NA
            if value == '' or value is None:
                value = 'NA'
            
            row[col] = value
        
        writer.writerow(row)
    
    csv_content = output.getvalue()
    output.close()
    
    if filename:
        try:
            with open(filename, 'w', newline='') as f:
                f.write(csv_content)
        except Exception as e:
            print(f"Error writing CSV file {filename}: {e}")
    
    return csv_content

def export_to_json(motifs: List[Dict[str, Any]], filename: Optional[str] = None, 
                   pretty: bool = True) -> str:
    """
    Export motifs to JSON format
    
    Args:
        motifs: List of motif dictionaries
        filename: Optional output filename
        pretty: Whether to format JSON prettily
        
    Returns:
        JSON format string
    """
    json_data = {
        'version': '2024.1',
        'analysis_type': 'NBDScanner_Non-B_DNA_Analysis',
        'total_motifs': len(motifs),
        'motifs': motifs
    }
    
    if pretty:
        json_content = json.dumps(json_data, indent=2, ensure_ascii=False)
    else:
        json_content = json.dumps(json_data, ensure_ascii=False)
    
    if filename:
        try:
            with open(filename, 'w') as f:
                f.write(json_content)
        except Exception as e:
            print(f"Error writing JSON file {filename}: {e}")
    
    return json_content


def export_to_excel(motifs: List[Dict[str, Any]], filename: str = "nonbscanner_results.xlsx") -> str:
    """
    Export motifs to Excel format with multiple sheets:
    - First sheet: Consolidated non-overlapping motifs
    - Subsequent sheets: Individual motif classes/subclasses
    
    Args:
        motifs: List of motif dictionaries
        filename: Output Excel filename (default: "nonbscanner_results.xlsx")
        
    Returns:
        Success message string
    """
    try:
        import openpyxl
    except ImportError:
        raise ImportError("openpyxl is required for Excel export. Install with: pip install openpyxl")
    
    if not motifs:
        return "No motifs to export"
    
    # Comprehensive column order
    comprehensive_columns = [
        'ID', 'Sequence_Name', 'Source', 'Class', 'Subclass', 'Pattern_ID',
        'Start', 'End', 'Length', 'Sequence', 'Method', 'Score',
        'Repeat_Type', 'Left_Arm', 'Right_Arm', 'Loop_Seq',
        'Arm_Length', 'Loop_Length', 'Stem_Length', 'Unit_Length',
        'Number_Of_Copies', 'Spacer_Length', 'Spacer_Sequence',
        'GC_Content', 'Structural_Features', 'Strand'
    ]
    
    # Get all unique keys from motifs
    all_keys = set()
    for motif in motifs:
        all_keys.update(motif.keys())
    
    # Add any additional fields found
    columns = comprehensive_columns.copy()
    for key in sorted(all_keys):
        if key not in columns:
            columns.append(key)
    
    # Prepare data rows
    def prepare_row(motif):
        row = {}
        for col in columns:
            value = motif.get(col, 'NA')
            
            # Map alternative field names
            if value == 'NA' or value == '' or value is None:
                if col == 'Number_Of_Copies' and 'Repeat_Units' in motif:
                    value = motif['Repeat_Units']
                elif col == 'Repeat_Type' and 'Tract_Type' in motif:
                    value = motif['Tract_Type']
                elif col == 'GC_Content' and 'GC_Total' in motif:
                    value = motif['GC_Total']
                elif col == 'Structural_Features':
                    features = []
                    if 'Tract_Type' in motif and motif['Tract_Type'] not in ['', 'NA', None]:
                        features.append(f"Tract:{motif['Tract_Type']}")
                    if 'Curvature_Score' in motif and motif['Curvature_Score'] not in ['', 'NA', None]:
                        features.append(f"Curvature:{motif['Curvature_Score']}")
                    if 'Z_Score' in motif and motif['Z_Score'] not in ['', 'NA', None]:
                        features.append(f"Z-Score:{motif['Z_Score']}")
                    value = '; '.join(features) if features else 'NA'
                
                if value == '' or value is None:
                    value = 'NA'
            
            row[col] = value
        return row
    
    # Create Excel writer
    with pd.ExcelWriter(filename, engine='openpyxl') as writer:
        # Sheet 1: Consolidated non-overlapping motifs (exclude Hybrid and Cluster)
        consolidated_motifs = [m for m in motifs 
                             if m.get('Class') not in ['Hybrid', 'Non-B_DNA_Clusters']]
        
        if consolidated_motifs:
            consolidated_data = [prepare_row(m) for m in consolidated_motifs]
            df_consolidated = pd.DataFrame(consolidated_data, columns=columns)
            df_consolidated.to_excel(writer, sheet_name='Consolidated_NonOverlapping', index=False)
        
        # Group motifs by class
        class_groups = defaultdict(list)
        for motif in motifs:
            cls = motif.get('Class', 'Unknown')
            class_groups[cls].append(motif)
        
        # Create separate sheets for each class
        for cls, class_motifs in sorted(class_groups.items()):
            # Sanitize sheet name (Excel has 31 character limit)
            sheet_name = cls.replace('/', '_').replace(' ', '_').replace('-', '_')[:31]
            
            class_data = [prepare_row(m) for m in class_motifs]
            df_class = pd.DataFrame(class_data, columns=columns)
            df_class.to_excel(writer, sheet_name=sheet_name, index=False)
            
            # Also create sheets by subclass if there are multiple subclasses
            subclass_groups = defaultdict(list)
            for motif in class_motifs:
                subclass = motif.get('Subclass', 'Other')
                subclass_groups[subclass].append(motif)
            
            # Only create subclass sheets if there are multiple subclasses
            if len(subclass_groups) > 1:
                for subclass, subclass_motifs in sorted(subclass_groups.items()):
                    # Sanitize subclass sheet name
                    subclass_name = f"{cls}_{subclass}".replace('/', '_').replace(' ', '_').replace('-', '_')[:31]
                    
                    subclass_data = [prepare_row(m) for m in subclass_motifs]
                    df_subclass = pd.DataFrame(subclass_data, columns=columns)
                    
                    # Try to add sheet, skip if name collision
                    try:
                        df_subclass.to_excel(writer, sheet_name=subclass_name, index=False)
                    except:
                        pass  # Skip if sheet name collision
    
    return f"Excel file exported successfully to {filename}"


def export_to_gff3(motifs: List[Dict[str, Any]], sequence_name: str = "sequence", 
                   filename: Optional[str] = None) -> str:
    """
    Export motifs to GFF3 format
    
    Args:
        motifs: List of motif dictionaries
        sequence_name: Name of the sequence
        filename: Optional output filename
        
    Returns:
        GFF3 format string
    """
    gff_lines = []
    gff_lines.append("##gff-version 3")
    gff_lines.append(f"##sequence-region {sequence_name} 1 {len(sequence_name)}")
    
    for i, motif in enumerate(motifs, 1):
        seqid = sequence_name
        source = "NBDScanner"
        feature_type = "Non_B_DNA_motif"
        start = motif.get('Start', 1)
        end = motif.get('End', start)
        score = motif.get('Score', '.')
        strand = motif.get('Strand', '+')
        phase = "."
        
        # Attributes
        attributes = [
            f"ID=motif_{i}",
            f"Name={motif.get('Class', 'Unknown')}_{motif.get('Subclass', 'Unknown')}",
            f"motif_class={motif.get('Class', 'Unknown')}",
            f"motif_subclass={motif.get('Subclass', 'Unknown')}",
            f"length={motif.get('Length', 0)}",
            f"method={motif.get('Method', 'NBDScanner')}"
        ]
        
        attributes_str = ';'.join(attributes)
        
        gff_line = f"{seqid}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes_str}"
        gff_lines.append(gff_line)
    
    gff_content = '\n'.join(gff_lines)
    
    if filename:
        try:
            with open(filename, 'w') as f:
                f.write(gff_content)
        except Exception as e:
            print(f"Error writing GFF3 file {filename}: {e}")
    
    return gff_content

# =============================================================================
# QUALITY CONTROL & FILTERING
# =============================================================================

def quality_check_motifs(motifs: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Perform quality checks on detected motifs
    
    Args:
        motifs: List of motif dictionaries
        
    Returns:
        Quality check report
    """
    if not motifs:
        return {'status': 'No motifs to check', 'passed': True, 'issues': []}
    
    issues = []
    
    # Check for required fields
    required_fields = ['Class', 'Subclass', 'Start', 'End', 'Sequence']
    for i, motif in enumerate(motifs):
        missing_fields = [field for field in required_fields if field not in motif]
        if missing_fields:
            issues.append(f"Motif {i+1}: Missing fields {missing_fields}")
    
    # Check coordinate consistency
    for i, motif in enumerate(motifs):
        start = motif.get('Start')
        end = motif.get('End')
        length = motif.get('Length')
        sequence = motif.get('Sequence', '')
        
        if start and end and start >= end:
            issues.append(f"Motif {i+1}: Invalid coordinates (start >= end)")
        
        if start and end and length and (end - start + 1) != length:
            issues.append(f"Motif {i+1}: Length inconsistent with coordinates")
        
        if sequence and length and len(sequence) != length:
            issues.append(f"Motif {i+1}: Sequence length doesn't match reported length")
    
    # Check for overlaps within same class
    class_groups = defaultdict(list)
    for motif in motifs:
        class_groups[motif.get('Class')].append(motif)
    
    for class_name, class_motifs in class_groups.items():
        for i, motif1 in enumerate(class_motifs):
            for motif2 in class_motifs[i+1:]:
                if (motif1.get('Start', 0) < motif2.get('End', 0) and 
                    motif2.get('Start', 0) < motif1.get('End', 0)):
                    issues.append(f"Overlapping motifs in class {class_name}")
                    break
    
    report = {
        'total_motifs': len(motifs),
        'issues_found': len(issues),
        'passed': len(issues) == 0,
        'issues': issues[:10],  # Limit to first 10 issues
        'status': 'PASSED' if len(issues) == 0 else f'FAILED ({len(issues)} issues)'
    }
    
    return report

def filter_motifs_by_score(motifs: List[Dict[str, Any]], min_score: float = 0.0) -> List[Dict[str, Any]]:
    """
    Filter motifs by minimum score threshold
    
    Args:
        motifs: List of motif dictionaries
        min_score: Minimum score threshold
        
    Returns:
        Filtered motifs list
    """
    return [m for m in motifs if m.get('Score', 0) >= min_score]

def filter_motifs_by_length(motifs: List[Dict[str, Any]], 
                           min_length: int = 0, max_length: int = float('inf')) -> List[Dict[str, Any]]:
    """
    Filter motifs by length range
    
    Args:
        motifs: List of motif dictionaries
        min_length: Minimum length
        max_length: Maximum length
        
    Returns:
        Filtered motifs list
    """
    return [m for m in motifs if min_length <= m.get('Length', 0) <= max_length]

def filter_motifs_by_class(motifs: List[Dict[str, Any]], 
                          allowed_classes: List[str]) -> List[Dict[str, Any]]:
    """
    Filter motifs by allowed classes
    
    Args:
        motifs: List of motif dictionaries
        allowed_classes: List of allowed class names
        
    Returns:
        Filtered motifs list
    """
    return [m for m in motifs if m.get('Class') in allowed_classes]

# =============================================================================
# CROSS-DETECTOR OVERLAP RESOLUTION (Hyperscan Integration Pattern)
# =============================================================================

def resolve_cross_class_overlaps(motifs: List[Dict[str, Any]], 
                                 mode: str = 'strict') -> List[Dict[str, Any]]:
    """
    Resolve overlaps across different motif classes using deterministic selection.
    
    Implements the overlap resolution strategy from the Hyperscan integration plan:
    - Option A (strict): Select highest normalized score among overlapping regions
    - Option B (hybrid): Keep both but mark as Hybrid (handled by detector pipeline)
    
    Algorithm:
    1. Sort candidates by score (descending), then by length (descending)
    2. Greedily select highest-scoring non-overlapping motifs
    3. This ensures deterministic, reproducible output
    
    Args:
        motifs: List of motif dictionaries from multiple detectors
        mode: 'strict' for non-overlapping selection, 'hybrid' for hybrid marking
        
    Returns:
        List of resolved motifs (non-overlapping if mode='strict')
        
    References:
        - Hyperscan integration pattern for Non-B DNA detection
        - Score-aware greedy selection (similar to Cruciform._remove_overlaps)
    """
    if not motifs:
        return []
    
    if mode == 'hybrid':
        # Return all motifs - hybrid detection happens in the scanner pipeline
        return motifs
    
    # Mode 'strict': Select highest-scoring non-overlapping set
    # Sort by score (descending), then by length (descending) for tie-breaking
    sorted_motifs = sorted(motifs, 
                          key=lambda x: (-x.get('Score', 0), 
                                       -(x.get('End', 0) - x.get('Start', 0))))
    
    selected = []
    
    for candidate in sorted_motifs:
        # Check if this candidate overlaps with any already selected motif
        overlaps = False
        cand_start = candidate.get('Start', 0)
        cand_end = candidate.get('End', 0)
        
        for selected_motif in selected:
            sel_start = selected_motif.get('Start', 0)
            sel_end = selected_motif.get('End', 0)
            
            # Check for overlap: two regions overlap if neither is completely before the other
            if not (cand_end <= sel_start or cand_start >= sel_end):
                overlaps = True
                break
        
        if not overlaps:
            selected.append(candidate)
    
    # Sort by start position for final output
    selected.sort(key=lambda x: x.get('Start', 0))
    
    return selected

def merge_detector_results(detector_results: Dict[str, List[Dict[str, Any]]],
                          overlap_mode: str = 'strict') -> List[Dict[str, Any]]:
    """
    Merge results from multiple detectors with overlap resolution.
    
    This function implements the complete Hyperscan integration pattern:
    1. Collect results from all detectors (already scored and filtered)
    2. Resolve cross-class overlaps according to specified mode
    3. Return unified, non-redundant motif list
    
    Args:
        detector_results: Dictionary mapping detector_name -> list of motifs
        overlap_mode: 'strict' for non-overlapping, 'hybrid' for overlap annotation
        
    Returns:
        Merged list of motifs with overlaps resolved
        
    Example:
        >>> results = {
        ...     'a_philic': [motif1, motif2],
        ...     'z_dna': [motif3, motif4],
        ...     'g_quadruplex': [motif5]
        ... }
        >>> merged = merge_detector_results(results, overlap_mode='strict')
    """
    # Flatten all detector results into single list
    all_motifs = []
    for detector_name, motifs in detector_results.items():
        # Add detector source to each motif for tracking
        for motif in motifs:
            if 'Detector' not in motif:
                motif['Detector'] = detector_name
            all_motifs.append(motif)
    
    # Resolve overlaps according to mode
    resolved = resolve_cross_class_overlaps(all_motifs, mode=overlap_mode)
    
    return resolved

def export_results_to_dataframe(motifs: List[Dict[str, Any]]) -> pd.DataFrame:
    """Convert motif results to pandas DataFrame with CORE fields only for display tables.
    
    Per requirements, output tables should only show:
    - Sequence_Name, Source, Class, Subclass, Start, End, Length, Sequence, Score
    
    Additional detailed columns (Repeat_Type, Left_Arm, Right_Arm, Loop_Seq, etc.)
    are only shown in Excel download per-motif sheets, not in display tables.
    """
    if not motifs:
        return pd.DataFrame()
    
    df = pd.DataFrame(motifs)
    
    # Use core columns constant
    core_columns = CORE_OUTPUT_COLUMNS
    
    # Ensure all core columns are present, fill missing with 'NA'
    for col in core_columns:
        if col not in df.columns:
            df[col] = 'NA'
    
    # Fill all NaN/None values with 'NA' string
    result_df = df[core_columns].fillna('NA')
    
    return result_df


# =============================================================================
# TESTING & EXAMPLES
# =============================================================================

def test_utilities():
    """Test utility functions with example data"""
    print("Testing NBDScanner utilities...")
    
    # Test sequence validation
    test_sequences = [
        "ATGCATGCATGC",     # Valid
        "ATGCXYZ",          # Invalid characters
        "ATG",              # Too short
        "",                 # Empty
    ]
    
    print("\nSequence validation tests:")
    for seq in test_sequences:
        valid, msg = validate_sequence(seq)
        print(f"  '{seq[:20]}': {'VALID' if valid else 'INVALID'} - {msg}")
    
    # Test basic statistics
    test_seq = "GGGTTAGGGTTAGGGTTAGGGAAAAATTTTCGCGCGCGCG"
    stats = get_basic_stats(test_seq)
    print(f"\nBasic stats for test sequence:")
    for key, value in stats.items():
        print(f"  {key}: {value}")
    
    # Test FASTA parsing
    fasta_content = """>test_sequence
GGGTTAGGGTTAGGGTTAGGG
>another_sequence
AAAAATTTTCCCCGGGG"""
    
    sequences = parse_fasta(fasta_content)
    print(f"\nFASTA parsing test: {len(sequences)} sequences parsed")
    for name, seq in sequences.items():
        print(f"  {name}: {seq}")
    
    # Test motif formatting
    example_motifs = [
        {
            'Class': 'G-Quadruplex',
            'Subclass': 'Canonical G4',
            'Start': 1,
            'End': 21,
            'Length': 21,
            'Sequence': 'GGGTTAGGGTTAGGGTTAGGG',
            'Score': 0.857,
            'Strand': '+'
        }
    ]
    
    quality_report = quality_check_motifs(example_motifs)
    print(f"\nQuality check: {quality_report['status']}")
    
    print("[OK] All utility tests completed")

# =============================================================================
# ENHANCED STATISTICS: DENSITY AND ENRICHMENT ANALYSIS
# =============================================================================

# Constants for enrichment analysis
DEFAULT_FOLD_ENRICHMENT_WHEN_ZERO_BACKGROUND = 1.0  # When background is zero


def shuffle_sequence(sequence: str, preserve_composition: bool = True) -> str:
    """
    Shuffle a DNA sequence while preserving nucleotide composition.
    
    Args:
        sequence: DNA sequence to shuffle
        preserve_composition: If True, maintains exact nucleotide counts
        
    Returns:
        Shuffled sequence
    """
    seq_list = list(sequence.upper())
    random.shuffle(seq_list)
    return ''.join(seq_list)


def calculate_genomic_density(motifs: List[Dict[str, Any]], 
                               sequence_length: int,
                               by_class: bool = True,
                               by_subclass: bool = False) -> Dict[str, float]:
    """
    Calculate genomic density (coverage) for motifs.
    
    Genomic Density (σ_G) = (Total unique bp covered by motifs / 
                             Total length in bp of analyzed region) × 100
    
    IMPORTANT: Uses set-based overlap handling to ensure coverage never exceeds 100%.
    If motifs overlap, only unique positions are counted.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total length of analyzed sequence
        by_class: If True, calculate density per motif class
        by_subclass: If True, calculate density per motif subclass (takes precedence over by_class)
        
    Returns:
        Dictionary with density metrics (percentage, capped at 100%)
        - If by_subclass=True: keys are 'Class:Subclass' format
        - If by_class=True: keys are class names
        - Otherwise: key is 'Overall'
    """
    if not motifs or sequence_length == 0:
        return {'Overall': 0.0}
    
    if not by_class and not by_subclass:
        # Overall density using set-based coverage (handles overlaps correctly)
        covered_positions = set()
        for motif in motifs:
            start = motif.get('Start', 0) - 1  # Convert to 0-based
            end = motif.get('End', 0)
            covered_positions.update(range(start, end))
        
        overall_density = min((len(covered_positions) / sequence_length) * 100, 100.0)
        return {'Overall': round(overall_density, 4)}
    
    # Density per subclass using set-based coverage
    if by_subclass:
        density_by_subclass = {}
        subclass_groups = defaultdict(list)
        
        for motif in motifs:
            class_name = motif.get('Class', 'Unknown')
            subclass_name = motif.get('Subclass', 'Unknown')
            key = f"{class_name}:{subclass_name}"
            subclass_groups[key].append(motif)
        
        # Calculate per-subclass density with overlap handling
        for subclass_key, subclass_motifs in subclass_groups.items():
            covered_positions = set()
            for motif in subclass_motifs:
                start = motif.get('Start', 0) - 1  # Convert to 0-based
                end = motif.get('End', 0)
                covered_positions.update(range(start, end))
            
            subclass_density = min((len(covered_positions) / sequence_length) * 100, 100.0)
            density_by_subclass[subclass_key] = round(subclass_density, 4)
        
        # Calculate overall density (all motifs combined)
        all_covered_positions = set()
        for motif in motifs:
            start = motif.get('Start', 0) - 1  # Convert to 0-based
            end = motif.get('End', 0)
            all_covered_positions.update(range(start, end))
        
        overall_density = min((len(all_covered_positions) / sequence_length) * 100, 100.0)
        density_by_subclass['Overall'] = round(overall_density, 4)
        
        return density_by_subclass
    
    # Density per class using set-based coverage
    density_by_class = {}
    class_groups = defaultdict(list)
    
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        class_groups[class_name].append(motif)
    
    # Calculate per-class density with overlap handling
    for class_name, class_motifs in class_groups.items():
        covered_positions = set()
        for motif in class_motifs:
            start = motif.get('Start', 0) - 1  # Convert to 0-based
            end = motif.get('End', 0)
            covered_positions.update(range(start, end))
        
        class_density = min((len(covered_positions) / sequence_length) * 100, 100.0)
        density_by_class[class_name] = round(class_density, 4)
    
    # Calculate overall density (all motifs combined)
    all_covered_positions = set()
    for motif in motifs:
        start = motif.get('Start', 0) - 1  # Convert to 0-based
        end = motif.get('End', 0)
        all_covered_positions.update(range(start, end))
    
    overall_density = min((len(all_covered_positions) / sequence_length) * 100, 100.0)
    density_by_class['Overall'] = round(overall_density, 4)
    
    return density_by_class


def calculate_positional_density(motifs: List[Dict[str, Any]], 
                                  sequence_length: int,
                                  unit: str = 'Mbp',
                                  by_class: bool = True,
                                  by_subclass: bool = False) -> Dict[str, float]:
    """
    Calculate positional density (frequency) for motifs.
    
    Positional Density (λ) = Total count of predicted motifs / 
                             Total length (in kbp or Mbp) of analyzed region
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total length of analyzed sequence in bp
        unit: 'kbp' or 'Mbp' for reporting units
        by_class: If True, calculate density per motif class
        by_subclass: If True, calculate density per motif subclass (takes precedence over by_class)
        
    Returns:
        Dictionary with positional density (motifs per unit)
        - If by_subclass=True: keys are 'Class:Subclass' format
        - If by_class=True: keys are class names
        - Otherwise: key is 'Overall'
    """
    if not motifs or sequence_length == 0:
        return {'Overall': 0.0}
    
    # Convert to appropriate unit
    if unit == 'kbp':
        sequence_length_unit = sequence_length / 1000
    elif unit == 'Mbp':
        sequence_length_unit = sequence_length / 1000000
    else:
        sequence_length_unit = sequence_length
    
    if not by_class and not by_subclass:
        # Overall positional density
        overall_density = len(motifs) / sequence_length_unit
        return {'Overall': round(overall_density, 2)}
    
    # Positional density per subclass
    if by_subclass:
        density_by_subclass = {}
        subclass_counts = Counter()
        
        for motif in motifs:
            class_name = motif.get('Class', 'Unknown')
            subclass_name = motif.get('Subclass', 'Unknown')
            key = f"{class_name}:{subclass_name}"
            subclass_counts[key] += 1
        
        for subclass_key, count in subclass_counts.items():
            subclass_density = count / sequence_length_unit
            density_by_subclass[subclass_key] = round(subclass_density, 2)
        
        # Also add overall
        density_by_subclass['Overall'] = round(len(motifs) / sequence_length_unit, 2)
        
        return density_by_subclass
    
    # Positional density per class
    density_by_class = {}
    class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    
    for class_name, count in class_counts.items():
        class_density = count / sequence_length_unit
        density_by_class[class_name] = round(class_density, 2)
    
    # Also add overall
    density_by_class['Overall'] = round(len(motifs) / sequence_length_unit, 2)
    
    return density_by_class


def calculate_enrichment_with_shuffling(motifs: List[Dict[str, Any]], 
                                       sequence: str,
                                       n_shuffles: int = 100,
                                       by_class: bool = True,
                                       by_subclass: bool = False,
                                       progress_callback=None) -> Dict[str, Any]:
    """
    Calculate fold enrichment and statistical significance using sequence shuffling.
    
    Compares observed motif density against background (shuffled sequences).
    
    Fold Enrichment = D_Observed / D_Background
    where D = Motif Density = Total bp of motif / Total bp of region
    
    Args:
        motifs: List of detected motifs in original sequence
        sequence: Original DNA sequence
        n_shuffles: Number of shuffling iterations (default: 100)
        by_class: If True, calculate enrichment per motif class
        by_subclass: If True, calculate enrichment per motif subclass (takes precedence over by_class)
        progress_callback: Optional callback function for progress updates
        
    Returns:
        Dictionary with enrichment metrics including:
        - fold_enrichment: Observed/Background ratio
        - p_value: Statistical significance
        - observed_density: Density in original sequence
        - background_mean: Mean density in shuffled sequences
        - background_std: Std deviation of background
    """
    # Import analyze_sequence locally to avoid circular dependency
    # (nonbscanner imports from utilities, so we can't import at module level)
    try:
        from nonbscanner import analyze_sequence
    except ImportError:
        logger.error("Failed to import analyze_sequence from nonbscanner")
        return {}
    
    if not motifs or not sequence:
        return {}
    
    sequence_length = len(sequence)
    
    # Calculate observed densities (genomic density is used for enrichment)
    observed_genomic_density = calculate_genomic_density(motifs, sequence_length, 
                                                         by_class=by_class, 
                                                         by_subclass=by_subclass)
    
    # Initialize background storage
    if by_subclass:
        # Create keys for all unique subclass combinations
        subclass_keys = set()
        for m in motifs:
            class_name = m.get('Class', 'Unknown')
            subclass_name = m.get('Subclass', 'Unknown')
            key = f"{class_name}:{subclass_name}"
            subclass_keys.add(key)
        background_densities = {key: [] for key in subclass_keys}
        background_densities['Overall'] = []
    elif by_class:
        class_names = set(m.get('Class', 'Unknown') for m in motifs)
        background_densities = {cls: [] for cls in class_names}
        background_densities['Overall'] = []
    else:
        background_densities = {'Overall': []}
    
    # Perform shuffling and detection
    for i in range(n_shuffles):
        if progress_callback:
            progress_callback(i + 1, n_shuffles)
        
        # Shuffle sequence
        shuffled_seq = shuffle_sequence(sequence)
        
        # Detect motifs in shuffled sequence
        try:
            shuffled_motifs = analyze_sequence(shuffled_seq, f"shuffled_{i}")
            
            # Calculate density for shuffled sequence
            shuffled_density = calculate_genomic_density(shuffled_motifs, sequence_length, 
                                                        by_class=by_class,
                                                        by_subclass=by_subclass)
            
            # Store background densities
            for key in background_densities.keys():
                background_densities[key].append(shuffled_density.get(key, 0.0))
        
        except Exception as e:
            # If shuffled analysis fails, record zero density
            logger.warning(f"Shuffled analysis {i} failed: {e}")
            for key in background_densities.keys():
                background_densities[key].append(0.0)
    
    # Calculate enrichment statistics
    enrichment_results = {}
    
    for group_name, bg_densities in background_densities.items():
        bg_densities = [d for d in bg_densities if d is not None]
        
        if not bg_densities:
            continue
        
        obs_density = observed_genomic_density.get(group_name, 0.0)
        bg_mean = np.mean(bg_densities)
        bg_std = np.std(bg_densities)
        
        # Calculate fold enrichment
        if bg_mean > 0:
            fold_enrichment = obs_density / bg_mean
        else:
            # When background is zero, use infinity if observed > 0, else default value
            fold_enrichment = float('inf') if obs_density > 0 else DEFAULT_FOLD_ENRICHMENT_WHEN_ZERO_BACKGROUND
        
        # Calculate p-value (proportion of shuffled >= observed)
        if bg_densities:
            p_value = sum(1 for bg in bg_densities if bg >= obs_density) / len(bg_densities)
        else:
            p_value = 1.0
        
        # Count motifs for this group
        if by_subclass and group_name != 'Overall':
            # Split the key to get class and subclass
            class_part, subclass_part = group_name.split(':', 1)
            observed_count = len([m for m in motifs 
                                if m.get('Class', 'Unknown') == class_part 
                                and m.get('Subclass', 'Unknown') == subclass_part])
        elif by_class and group_name != 'Overall':
            observed_count = len([m for m in motifs if m.get('Class', 'Unknown') == group_name])
        else:
            observed_count = len(motifs)
        
        enrichment_results[group_name] = {
            'observed_density': round(obs_density, 4),
            'background_mean': round(bg_mean, 4),
            'background_std': round(bg_std, 4),
            'fold_enrichment': round(fold_enrichment, 2) if not np.isinf(fold_enrichment) else 'Inf',
            'p_value': round(p_value, 4),
            'n_shuffles': n_shuffles,
            'observed_count': observed_count
        }
    
    return enrichment_results


def calculate_enhanced_statistics(motifs: List[Dict[str, Any]], 
                                  sequence: str,
                                  include_enrichment: bool = True,
                                  n_shuffles: int = 100,
                                  progress_callback=None) -> Dict[str, Any]:
    """
    Calculate comprehensive statistics including density and enrichment analysis.
    
    Args:
        motifs: List of detected motifs
        sequence: Original DNA sequence
        include_enrichment: Whether to perform enrichment analysis (time-consuming)
        n_shuffles: Number of shuffles for enrichment analysis
        progress_callback: Optional callback for progress updates
        
    Returns:
        Dictionary with comprehensive statistics
    """
    sequence_length = len(sequence)
    
    # Basic statistics
    basic_stats = calculate_motif_statistics(motifs, sequence_length)
    
    # Density calculations
    genomic_density = calculate_genomic_density(motifs, sequence_length, by_class=True)
    positional_density_kbp = calculate_positional_density(motifs, sequence_length, unit='kbp', by_class=True)
    positional_density_mbp = calculate_positional_density(motifs, sequence_length, unit='Mbp', by_class=True)
    
    enhanced_stats = {
        'basic': basic_stats,
        'genomic_density': genomic_density,
        'positional_density_per_kbp': positional_density_kbp,
        'positional_density_per_mbp': positional_density_mbp
    }
    
    # Enrichment analysis (optional, time-consuming)
    if include_enrichment and motifs:
        enrichment_results = calculate_enrichment_with_shuffling(
            motifs, sequence, n_shuffles=n_shuffles, 
            by_class=True, progress_callback=progress_callback
        )
        enhanced_stats['enrichment'] = enrichment_results
    
    return enhanced_stats


if __name__ == "__main__":
    test_utilities()

# =============================================================================
# VISUALIZATION SUITE (merged from visualizations.py)
# =============================================================================
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    NBDSCANNER VISUALIZATION SUITE                             ║
║        Comprehensive Plotting and Visualization Functions                    ║
╚══════════════════════════════════════════════════════════════════════════════╝

VISUALIZATION FUNCTIONS (integrated into utilities.py)
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.1
LICENSE: MIT

DESCRIPTION:
    Comprehensive plotting and visualization functions for Non-B DNA motif analysis.
    Includes both static and interactive visualizations with scientific styling.

VISUALIZATION TABLE:
┌──────────────────────────────────────────────────────────────────────────────┐
│ Category     │ Functions                 │ Description                       │
├──────────────┼───────────────────────────┼───────────────────────────────────┤
│ Distribution │ plot_motif_distribution   │ Class/subclass distribution plots │
│ Coverage     │ plot_coverage_map         │ Sequence coverage visualization   │
│ Statistics   │ plot_length_distribution  │ Length distributions              │
│ Hierarchy    │ plot_nested_pie_chart     │ Class-subclass hierarchy          │
│ Export       │ save_all_plots            │ Batch plot export                 │
└──────────────────────────────────────────────────────────────────────────────┘

ADVANCED VISUALIZATIONS (Publication-Quality):
┌──────────────────────────────────────────────────────────────────────────────┐
│ - plot_genome_landscape_track      │ Horizontal genomic ruler with glyphs │
│ - plot_sliding_window_heat_ribbon  │ 1D density heatmap with score        │
│ - plot_ridge_plots_length_by_class │ Stacked density ridges (joyplots)    │
│ - plot_sunburst_treemap            │ Hierarchical composition             │
│ - plot_hexbin_start_vs_score       │ 2D hexbin with marginals             │
│ - plot_upset_intersection          │ UpSet plot for motif overlaps        │
│ - plot_score_violin_beeswarm       │ Score distributions with points      │
│ - plot_cluster_hotspot_map         │ Cluster regions with annotations     │
└──────────────────────────────────────────────────────────────────────────────┘

FEATURES:
    - Publication-quality static plots
    - Interactive visualizations (Plotly)
    - Colorblind-friendly palettes
    - SVG/PNG export at 300 DPI
    - Customizable styling
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Tuple, Union
from collections import Counter, defaultdict
import warnings
warnings.filterwarnings("ignore")

# Try to import plotly for interactive plots
try:
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

# =============================================================================
# NATURE-LEVEL PUBLICATION STYLING & CONFIGURATION
# =============================================================================
# Following Nature Methods/Nature Genetics style guidelines:
# - Sans-serif fonts (Arial/Helvetica preferred)
# - High DPI (300+ for print)
# - Clean, minimal design with proper spacing
# - Colorblind-friendly palettes

# Default DPI for publication quality (Nature requires 300 DPI minimum)
PUBLICATION_DPI = 300

# Nature-level color palette for motif classes (colorblind-friendly)
# Based on Wong, B. (2011) Nature Methods colorblind-safe palette
# Each motif class has a unique color for distinguishability
# Reference: Wong, B. (2011) Points of view: Color blindness. Nat Methods 8, 441
MOTIF_CLASS_COLORS = {
    'Curved_DNA': '#CC79A7',          # Reddish Purple
    'Slipped_DNA': '#E69F00',         # Orange
    'Cruciform': '#56B4E9',           # Sky Blue
    'R-Loop': '#009E73',              # Bluish Green
    'Triplex': '#F0E442',             # Yellow
    'G-Quadruplex': '#0072B2',        # Blue
    'i-Motif': '#D55E00',             # Vermillion
    'Z-DNA': '#882255',               # Wine (distinct from Curved_DNA)
    'A-philic_DNA': '#44AA99',        # Teal (distinct from Cruciform)
    'Hybrid': '#999999',              # Gray - neutral
    'Non-B_DNA_Clusters': '#666666'   # Dark Gray - professional
}

# Helper function to format display names
def _format_display_name(name: str) -> str:
    """Convert internal names to publication-ready display format.
    
    Replaces underscores with spaces for better readability.
    
    Args:
        name: Internal name (e.g., 'Curved_DNA')
        
    Returns:
        Display name (e.g., 'Curved DNA')
    """
    return name.replace('_', ' ')

# Nature-level scientific styling configuration for publication-quality plots
# Reference: Nature author guidelines for figure preparation
_NATURE_STYLE_PARAMS = {
    # Typography - Arial/Helvetica as per Nature guidelines
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,  # Nature recommends 5-7pt, we use 8pt for readability
    
    # Title and labels
    'axes.titlesize': 9,
    'axes.titleweight': 'bold',
    'axes.labelsize': 8,
    'axes.labelweight': 'normal',
    
    # Tick labels
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    
    # Legend
    'legend.fontsize': 7,
    'legend.frameon': False,
    'legend.borderpad': 0.4,
    
    # Figure
    'figure.titlesize': 10,
    'figure.titleweight': 'bold',
    'figure.dpi': PUBLICATION_DPI,
    'figure.facecolor': 'white',
    'savefig.dpi': PUBLICATION_DPI,
    'savefig.format': 'pdf',  # Vector format for publication
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    
    # Axes - minimal, clean design
    'axes.grid': False,  # Nature typically uses minimal grids
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.spines.left': True,
    'axes.spines.bottom': True,
    'axes.linewidth': 0.8,
    'axes.edgecolor': 'black',
    'axes.facecolor': 'white',
    
    # Lines
    'lines.linewidth': 1.0,
    'lines.markersize': 4,
    
    # PDF/SVG output
    'pdf.fonttype': 42,  # TrueType for editability
    'ps.fonttype': 42,
}

# Apply Nature-level styling at module load
plt.rcParams.update(_NATURE_STYLE_PARAMS)

# Constants for enrichment analysis visualization
INFINITE_FOLD_ENRICHMENT_CAP = 100  # Cap for infinite fold enrichment values in plots

# Figure size constants (in inches) following Nature guidelines
# Nature column widths: single=89mm(3.5in), 1.5=120mm(4.7in), double=183mm(7.2in)
FIGURE_SIZES = {
    'single_column': (3.5, 2.8),      # Single column
    'one_and_half': (4.7, 3.5),       # 1.5 column
    'double_column': (7.2, 4.5),      # Double column
    'square': (3.5, 3.5),             # Square figure
    'wide': (7.2, 3.0),               # Wide panoramic
}


def set_scientific_style(style: str = 'nature'):
    """
    Apply scientific publication-ready styling.
    
    Args:
        style: Style preset ('nature', 'default', 'presentation')
               - 'nature': Nature/Science journal style (default)
               - 'default': Standard seaborn whitegrid
               - 'presentation': Larger fonts for presentations
    """
    if style == 'nature':
        plt.rcParams.update(_NATURE_STYLE_PARAMS)
        sns.set_style("white")
        sns.set_context("paper")
    elif style == 'presentation':
        plt.rcParams.update(_NATURE_STYLE_PARAMS)
        plt.rcParams.update({
            'font.size': 12,
            'axes.titlesize': 14,
            'axes.labelsize': 12,
            'xtick.labelsize': 11,
            'ytick.labelsize': 11,
            'legend.fontsize': 10,
        })
        sns.set_style("whitegrid")
        sns.set_context("talk")
    else:
        sns.set_style("whitegrid")
        sns.set_palette("husl")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def _apply_nature_style(ax):
    """Apply Nature journal style to matplotlib axis (remove top/right spines).
    
    | Parameter | Type | Description            |
    |-----------|------|------------------------|
    | ax        | Axes | Matplotlib axis object |
    """
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# =============================================================================
# DISTRIBUTION PLOTS (Nature-Level Quality)
# =============================================================================

def plot_motif_distribution(motifs: List[Dict[str, Any]], 
                           by: str = 'Class',
                           title: Optional[str] = None,
                           figsize: Tuple[float, float] = None,
                           style: str = 'nature') -> plt.Figure:
    """
    Plot distribution of motifs by class or subclass.
    
    Publication-quality bar chart following Nature Methods guidelines.
    Shows all classes/subclasses even when count is 0, ensuring comprehensive
    visualization of both detected and undetected motif types.
    
    Args:
        motifs: List of motif dictionaries.
        by: Group by 'Class' or 'Subclass'.
        title: Custom plot title.
        figsize: Figure size (width, height) in inches. Uses Nature standard if None.
        style: Style preset ('nature', 'default', 'presentation').
        
    Returns:
        Matplotlib figure object (publication-ready at 300 DPI)
    """
    set_scientific_style(style)
    
    # Use Nature-appropriate figure size
    if figsize is None:
        figsize = FIGURE_SIZES['double_column'] if by == 'Subclass' else FIGURE_SIZES['one_and_half']
    
    # Define ALL expected classes and subclasses (always show these)
    ALL_CLASSES = [
        'Curved_DNA', 'Slipped_DNA', 'Cruciform', 'R-Loop', 'Triplex',
        'G-Quadruplex', 'i-Motif', 'Z-DNA', 'A-philic_DNA', 'Hybrid', 'Non-B_DNA_Clusters'
    ]
    
    ALL_SUBCLASSES = [
        'Global Curvature', 'Local Curvature',  # Curved DNA
        'Direct Repeat', 'STR',  # Slipped DNA
        'Inverted Repeats',  # Cruciform
        'R-loop formation sites', 'QmRLFS-m1', 'QmRLFS-m2',  # R-Loop
        'Triplex', 'Sticky DNA',  # Triplex
        'Canonical G4', 'Relaxed G4', 'Bulged G4', 'Bipartite G4',  # G-Quadruplex
        'Multimeric G4', 'Imperfect G4', 'G-Triplex intermediate',
        'Canonical i-motif', 'Relaxed i-motif', 'AC-motif',  # i-Motif
        'Z-DNA', 'eGZ (Extruded-G) DNA',  # Z-DNA
        'A-philic DNA',  # A-philic
    ]
    
    # Count motifs by specified grouping
    counts = Counter(m.get(by, 'Unknown') for m in motifs) if motifs else Counter()
    
    # Prepare data with all categories
    if by == 'Class':
        categories = ALL_CLASSES
    else:
        categories = ALL_SUBCLASSES
    
    # Get counts (0 if not present)
    values = [counts.get(cat, 0) for cat in categories]
    
    # Get colors (colorblind-friendly)
    if by == 'Class':
        colors = [MOTIF_CLASS_COLORS.get(cat, '#808080') for cat in categories]
    else:
        # Map subclasses to their parent class colors
        subclass_to_class = {
            'Global Curvature': 'Curved_DNA', 'Local Curvature': 'Curved_DNA',
            'Direct Repeat': 'Slipped_DNA', 'STR': 'Slipped_DNA',
            'Inverted Repeats': 'Cruciform',
            'R-loop formation sites': 'R-Loop', 'QmRLFS-m1': 'R-Loop', 'QmRLFS-m2': 'R-Loop',
            'Triplex': 'Triplex', 'Sticky DNA': 'Triplex',
            'Canonical G4': 'G-Quadruplex', 'Relaxed G4': 'G-Quadruplex',
            'Bulged G4': 'G-Quadruplex', 'Bipartite G4': 'G-Quadruplex',
            'Multimeric G4': 'G-Quadruplex', 'Imperfect G4': 'G-Quadruplex',
            'G-Triplex intermediate': 'G-Quadruplex',
            'Canonical i-motif': 'i-Motif', 'Relaxed i-motif': 'i-Motif', 'AC-motif': 'i-Motif',
            'Z-DNA': 'Z-DNA', 'eGZ (Extruded-G) DNA': 'Z-DNA',
            'A-philic DNA': 'A-philic_DNA',
        }
        colors = [MOTIF_CLASS_COLORS.get(subclass_to_class.get(cat, ''), '#808080') for cat in categories]
    
    # Create figure with high DPI
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    # Create bars with Nature-style aesthetics
    bars = ax.bar(range(len(categories)), values, color=colors, 
                  edgecolor='black', linewidth=0.5, width=0.8)
    
    # Customize axes (Nature style - minimal, clean)
    ax.set_xlabel(f'Motif {by}', fontweight='normal')
    ax.set_ylabel('Count', fontweight='normal')
    if title:
        # Replace underscores with spaces in title
        display_title = title.replace('_', ' ')
        ax.set_title(display_title, fontweight='bold', pad=10)
    ax.set_xticks(range(len(categories)))
    
    # Replace underscores with spaces in category labels
    # Apply 45° rotation for all categories (Nature style - consistent readability)
    display_categories = [cat.replace('_', ' ') for cat in categories]
    ax.set_xticklabels(display_categories, rotation=45, ha='right')
    
    # Add count labels on bars (small font, positioned above)
    if len(categories) <= 20:
        max_val = max(values) if max(values) > 0 else 1
        for bar, count in zip(bars, values):
            if count > 0:  # Only label non-zero bars
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + max_val * 0.02,
                        str(count), ha='center', va='bottom', fontsize=6)
    
    # Apply Nature journal style
    _apply_nature_style(ax)
    
    plt.tight_layout()
    return fig


def plot_class_subclass_sunburst(motifs: List[Dict[str, Any]], 
                                 title: str = "Motif Class-Subclass Distribution",
                                 figsize: Tuple[float, float] = None) -> Union[plt.Figure, Any]:
    """
    Create sunburst plot showing class-subclass hierarchy.
    
    Publication-quality hierarchical visualization.
    
    Args:
        motifs: List of motif dictionaries
        title: Plot title
        figsize: Figure size (width, height) in inches
        
    Returns:
        Plotly figure if available, otherwise matplotlib figure
    """
    if figsize is None:
        figsize = FIGURE_SIZES['square']
        
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center', 
                transform=ax.transAxes)
        ax.axis('off')
        if title:
            ax.set_title(title)
        return fig
    
    if not PLOTLY_AVAILABLE:
        # Fallback to matplotlib nested pie chart
        return plot_nested_pie_chart(motifs, title)
    
    # Build hierarchical data for sunburst
    class_subclass_counts = defaultdict(lambda: defaultdict(int))
    
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        subclass_name = motif.get('Subclass', 'Unknown')
        class_subclass_counts[class_name][subclass_name] += 1
    
    # Prepare data for plotly
    ids = []
    labels = []
    parents = []
    values = []
    colors = []
    
    # Add classes (inner ring)
    for class_name, subclasses in class_subclass_counts.items():
        total_class_count = sum(subclasses.values())
        ids.append(class_name)
        labels.append(f"{class_name}<br>({total_class_count})")
        parents.append("")
        values.append(total_class_count)
        colors.append(MOTIF_CLASS_COLORS.get(class_name, '#808080'))
    
    # Add subclasses (outer ring)
    for class_name, subclasses in class_subclass_counts.items():
        for subclass_name, count in subclasses.items():
            ids.append(f"{class_name}_{subclass_name}")
            labels.append(f"{subclass_name}<br>({count})")
            parents.append(class_name)
            values.append(count)
            # Lighter shade of class color for subclasses
            base_color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
            colors.append(base_color + '80')  # Add transparency
    
    fig = go.Figure(go.Sunburst(
        ids=ids,
        labels=labels,
        parents=parents,
        values=values,
        branchvalues="total",
        marker=dict(colors=colors, line=dict(color="#FFFFFF", width=1)),
        hovertemplate='<b>%{label}</b><br>Count: %{value}<extra></extra>',
    ))
    
    # Publication-quality layout
    fig.update_layout(
        title=dict(text=title, font=dict(size=10, family='Arial')),
        font=dict(size=8, family='Arial'),
        width=int(figsize[0] * 100),
        height=int(figsize[1] * 100),
        margin=dict(t=30, l=10, r=10, b=10)
    )
    
    return fig


def plot_nested_pie_chart(motifs: List[Dict[str, Any]], 
                         title: str = "Motif Distribution",
                         figsize: Tuple[float, float] = None) -> plt.Figure:
    """
    Create nested donut chart with improved text placement to avoid overlapping labels.
    
    Publication-quality hierarchical pie chart following Nature guidelines.
    
    Args:
        motifs: List of motif dictionaries
        title: Plot title
        figsize: Figure size (width, height) in inches
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style('nature')
    
    if figsize is None:
        figsize = FIGURE_SIZES['square']
    
    # Count by class and subclass
    class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    class_subclass_counts = defaultdict(lambda: defaultdict(int))
    
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        subclass_name = motif.get('Subclass', 'Unknown')
        class_subclass_counts[class_name][subclass_name] += 1
    
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    # Inner donut (classes)
    class_names = list(class_counts.keys())
    class_values = list(class_counts.values())
    class_colors = [MOTIF_CLASS_COLORS.get(name, '#808080') for name in class_names]
    
    # Create inner donut with Nature-style clean design
    # Use labels=None to manually place labels later for better control
    wedges1, texts1, autotexts1 = ax.pie(
        class_values, 
        labels=None,  # We'll add labels manually
        colors=class_colors,
        radius=0.65,
        autopct=lambda pct: f'{pct:.1f}%' if pct > 3 else '',  # Show more percentage labels with 1 decimal
        pctdistance=0.80,
        startangle=90,
        wedgeprops=dict(width=0.35, edgecolor='white', linewidth=2)  # Thicker edge for better clarity
    )
    
    # Manually add class labels with better positioning (replace underscores with spaces)
    for i, (wedge, class_name) in enumerate(zip(wedges1, class_names)):
        angle = (wedge.theta2 + wedge.theta1) / 2
        x = 0.5 * np.cos(np.radians(angle))
        y = 0.5 * np.sin(np.radians(angle))
        # Replace underscores with spaces in labels
        display_name = class_name.replace('_', ' ')
        # Add white background box for better readability
        ax.text(x, y, display_name, ha='center', va='center', fontsize=8, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.8))
    
    # Outer donut (subclasses)
    all_subclass_counts = []
    all_subclass_colors = []
    all_subclass_labels = []
    
    for class_name in class_names:
        subclass_dict = class_subclass_counts[class_name]
        base_color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        
        for subclass_name, count in subclass_dict.items():
            all_subclass_counts.append(count)
            # Truncate long names for clean appearance and replace underscores with spaces
            # Use consistent truncation length (15 chars max, including ellipsis)
            display_name = subclass_name.replace('_', ' ')
            MAX_LABEL_LENGTH = 15
            label = display_name if len(display_name) <= MAX_LABEL_LENGTH else display_name[:MAX_LABEL_LENGTH-1] + '…'
            all_subclass_labels.append(label)
            all_subclass_colors.append(base_color)
    
    # Use smarter labeling strategy to avoid overlap
    # For many subclasses, hide labels and rely on legend instead
    if len(all_subclass_labels) > 10:
        # Hide outer ring labels when there are too many
        wedges2, texts2 = ax.pie(
            all_subclass_counts,
            labels=None,  # No labels for cleaner appearance
            colors=all_subclass_colors,
            radius=1.0,
            startangle=90,
            wedgeprops=dict(width=0.35, edgecolor='white', linewidth=1.5)  # Thicker edge for clarity
        )
        
        # Add a legend for subclasses instead
        # Note: All subclasses of the same parent class share the same color by design.
        # This ensures visual grouping in the nested donut chart.
        # Track first occurrence of each unique label for the legend.
        seen_labels = {}
        legend_handles = []
        legend_labels = []
        for i, (label, color) in enumerate(zip(all_subclass_labels, all_subclass_colors)):
            if label not in seen_labels:
                seen_labels[label] = color
                legend_handles.append(plt.Rectangle((0,0),1,1, fc=color, ec='white', lw=1))
                legend_labels.append(label)
                if len(legend_labels) >= 12:  # Limit to 12 for better display
                    break
        
        ax.legend(legend_handles, legend_labels, loc='center left', bbox_to_anchor=(1, 0.5),
                 fontsize=7, frameon=True, title='Top Subclasses', title_fontsize=8,
                 framealpha=0.95, edgecolor='lightgray')
    else:
        # For fewer subclasses, show labels with improved spacing
        wedges2, texts2 = ax.pie(
            all_subclass_counts,
            labels=all_subclass_labels,
            colors=all_subclass_colors,
            radius=1.0,
            labeldistance=1.18,  # Push labels further out to avoid overlap
            startangle=90,
            wedgeprops=dict(width=0.35, edgecolor='white', linewidth=1.5),  # Thicker edge
            textprops={'fontsize': 7, 'weight': 'medium'}  # Larger, bolder text
        )
        
        # Adjust label positions to avoid overlap
        for text in texts2:
            text.set_fontsize(7)
            text.set_weight('medium')
            # Add slight rotation for better readability
            angle = text.get_rotation()
            if 90 < angle < 270:
                text.set_rotation(angle - 180)
    
    if title:
        # Replace underscores with spaces in title
        display_title = title.replace('_', ' ')
        ax.set_title(display_title, fontweight='bold', fontsize=12, pad=15)
    
    # Style percentage labels - larger and bolder
    for autotext in autotexts1:
        autotext.set_fontsize(7)
        autotext.set_weight('bold')
        autotext.set_color('white')
    
    return fig


# =============================================================================
# COVERAGE & POSITIONAL PLOTS (Nature-Level Quality)
# =============================================================================

def plot_coverage_map(motifs: List[Dict[str, Any]], 
                     sequence_length: int,
                     title: Optional[str] = None,
                     figsize: Tuple[float, float] = None,
                     style: str = 'nature') -> plt.Figure:
    """
    Plot motif coverage map showing positions along sequence.
    
    Publication-quality genomic track visualization.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of the analyzed sequence
        title: Custom plot title
        figsize: Figure size (width, height) in inches
        style: Style preset ('nature', 'default', 'presentation')
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style(style)
    
    if figsize is None:
        figsize = FIGURE_SIZES['wide']
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center', 
                transform=ax.transAxes)
        ax.axis('off')
        if title:
            ax.set_title(title)
        return fig
    
    # Group motifs by class
    class_motifs = defaultdict(list)
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        class_motifs[class_name].append(motif)
    
    # Create figure with high DPI
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    y_pos = 0
    class_positions = {}
    
    for class_name, class_motif_list in class_motifs.items():
        class_positions[class_name] = y_pos
        color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        
        for motif in class_motif_list:
            start = motif.get('Start', 0) - 1  # Convert to 0-based
            end = motif.get('End', start + 1)
            length = end - start
            
            # Draw motif as rectangle (Nature style - clean edges)
            rect = patches.Rectangle(
                (start, y_pos - 0.35), length, 0.7,
                facecolor=color, edgecolor='black', linewidth=0.3
            )
            ax.add_patch(rect)
        
        y_pos += 1
    
    # Customize axes (Nature style)
    ax.set_xlim(0, sequence_length)
    ax.set_ylim(-0.5, len(class_motifs) - 0.5)
    ax.set_xlabel('Position (bp)')
    ax.set_ylabel('Motif Class')
    if title:
        # Replace underscores with spaces in title
        display_title = title.replace('_', ' ')
        ax.set_title(display_title, fontweight='bold', pad=10)
    
    # Set y-axis labels with underscores replaced by spaces
    ax.set_yticks(list(class_positions.values()))
    display_labels = [label.replace('_', ' ') for label in class_positions.keys()]
    ax.set_yticklabels(display_labels)
    
    # Clean x-axis ticks
    ax.ticklabel_format(style='sci', axis='x', scilimits=(3, 3))
    
    # Apply Nature journal style
    _apply_nature_style(ax)
    
    plt.tight_layout()
    return fig


def plot_density_heatmap(motifs: List[Dict[str, Any]], 
                        sequence_length: int,
                        window_size: int = 1000,
                        title: Optional[str] = None,
                        figsize: Tuple[float, float] = None,
                        style: str = 'nature') -> plt.Figure:
    """
    Plot motif density heatmap along sequence.
    
    Publication-quality density visualization.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of the analyzed sequence
        window_size: Window size for density calculation
        title: Custom plot title
        figsize: Figure size (width, height) in inches
        style: Style preset
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style(style)
    
    if figsize is None:
        figsize = FIGURE_SIZES['wide']
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center', 
                transform=ax.transAxes)
        ax.axis('off')
        if title:
            ax.set_title(title)
        return fig
    
    # Calculate windows
    num_windows = max(1, sequence_length // window_size)
    windows = np.linspace(0, sequence_length, num_windows + 1)
    
    # Get unique classes
    classes = sorted(set(m.get('Class', 'Unknown') for m in motifs))
    
    # Calculate density matrix
    density_matrix = np.zeros((len(classes), num_windows))
    
    for i, class_name in enumerate(classes):
        class_motifs = [m for m in motifs if m.get('Class') == class_name]
        
        for j in range(num_windows):
            window_start = windows[j]
            window_end = windows[j + 1]
            
            # Count motifs in window
            count = 0
            for motif in class_motifs:
                motif_start = motif.get('Start', 0)
                motif_end = motif.get('End', 0)
                
                # Check if motif overlaps with window
                if not (motif_end <= window_start or motif_start >= window_end):
                    count += 1
            
            density_matrix[i, j] = count
    
    # Create heatmap with publication quality
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    # Use colorblind-friendly colormap
    im = ax.imshow(density_matrix, cmap='viridis', aspect='auto', interpolation='nearest')
    
    # Customize axes (Nature style)
    ax.set_xlabel(f'Position (kb)')
    ax.set_ylabel('Motif Class')
    if title:
        # Replace underscores with spaces in title
        display_title = title.replace('_', ' ')
        ax.set_title(display_title, fontweight='bold', pad=10)
    
    # Set ticks and labels with underscores replaced by spaces
    ax.set_yticks(range(len(classes)))
    display_classes = [cls.replace('_', ' ') for cls in classes]
    ax.set_yticklabels(display_classes)
    
    # Clean x-axis with kb units
    x_ticks = np.arange(0, num_windows, max(1, num_windows // 5))
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f'{int(windows[i]/1000)}' for i in x_ticks])
    
    # Add colorbar with proper label
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Count', fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    
    plt.tight_layout()
    return fig


# =============================================================================
# STATISTICAL PLOTS (Nature-Level Quality)
# =============================================================================

def plot_score_distribution(motifs: List[Dict[str, Any]], 
                           by_class: bool = True,
                           title: Optional[str] = None,
                           figsize: Tuple[float, float] = None,
                           style: str = 'nature') -> plt.Figure:
    """
    Plot distribution of motif scores.
    
    Publication-quality box plot visualization.
    
    Args:
        motifs: List of motif dictionaries
        by_class: Whether to separate by motif class
        title: Custom plot title
        figsize: Figure size (width, height) in inches
        style: Style preset
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style(style)
    
    if figsize is None:
        figsize = FIGURE_SIZES['one_and_half']
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center', 
                transform=ax.transAxes)
        ax.axis('off')
        if title:
            ax.set_title(title)
        return fig
    
    # Extract scores
    scores_data = []
    for motif in motifs:
        score = motif.get('Score', motif.get('Normalized_Score'))
        if isinstance(score, (int, float)):
            if by_class:
                scores_data.append({
                    'Score': score,
                    'Class': motif.get('Class', 'Unknown')
                })
            else:
                scores_data.append(score)
    
    if not scores_data:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No score data available', ha='center', va='center', 
                transform=ax.transAxes)
        ax.axis('off')
        if title:
            display_title = title.replace('_', ' ')
            ax.set_title(display_title)
        return fig
    
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    if by_class and isinstance(scores_data[0], dict):
        # Create DataFrame for seaborn
        df = pd.DataFrame(scores_data)
        
        # Box plot by class (Nature style - clean, minimal)
        colors = [MOTIF_CLASS_COLORS.get(cls, '#808080') for cls in df['Class'].unique()]
        sns.boxplot(data=df, x='Class', y='Score', ax=ax, palette=colors,
                   linewidth=0.8, fliersize=2)
        # Replace underscores with spaces in x-tick labels
        ax.set_xticklabels([label.get_text().replace('_', ' ') for label in ax.get_xticklabels()], 
                          rotation=45, ha='right')
        ax.set_ylabel('Score')
        ax.set_xlabel('Motif Class')
    else:
        # Simple histogram
        ax.hist(scores_data, bins=20, edgecolor='black', linewidth=0.5, color='#56B4E9')
        ax.set_xlabel('Score')
        ax.set_ylabel('Frequency')
    
    if title:
        # Replace underscores with spaces in title
        display_title = title.replace('_', ' ')
        ax.set_title(display_title, fontweight='bold', pad=10)
    
    # Apply Nature journal style
    _apply_nature_style(ax)
    
    plt.tight_layout()
    return fig


def plot_length_distribution(motifs: List[Dict[str, Any]], 
                           by_class: bool = True,
                           title: Optional[str] = None,
                           figsize: Tuple[float, float] = None,
                           style: str = 'nature') -> plt.Figure:
    """
    Plot distribution of motif lengths.
    
    Publication-quality violin plot visualization.
    
    Args:
        motifs: List of motif dictionaries
        by_class: Whether to separate by motif class
        title: Custom plot title
        figsize: Figure size (width, height) in inches
        style: Style preset
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style(style)
    
    if figsize is None:
        figsize = FIGURE_SIZES['one_and_half']
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center', 
                transform=ax.transAxes)
        ax.axis('off')
        if title:
            ax.set_title(title)
        return fig
    
    # Extract lengths
    length_data = []
    for motif in motifs:
        length = motif.get('Length')
        if isinstance(length, int) and length > 0:
            if by_class:
                length_data.append({
                    'Length': length,
                    'Class': motif.get('Class', 'Unknown')
                })
            else:
                length_data.append(length)
    
    if not length_data:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No length data available', ha='center', va='center', 
                transform=ax.transAxes)
        ax.axis('off')
        if title:
            display_title = title.replace('_', ' ')
            ax.set_title(display_title)
        return fig
    
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    if by_class and isinstance(length_data[0], dict):
        # Create DataFrame for seaborn
        df = pd.DataFrame(length_data)
        
        # Violin plot by class (Nature style - clean)
        colors = [MOTIF_CLASS_COLORS.get(cls, '#808080') for cls in df['Class'].unique()]
        sns.violinplot(data=df, x='Class', y='Length', ax=ax, palette=colors,
                      linewidth=0.8, inner='box')
        # Replace underscores with spaces in x-tick labels
        ax.set_xticklabels([label.get_text().replace('_', ' ') for label in ax.get_xticklabels()], 
                          rotation=45, ha='right')
        ax.set_ylabel('Length (bp)')
        ax.set_xlabel('Motif Class')
    else:
        # Simple histogram
        ax.hist(length_data, bins=20, edgecolor='black', linewidth=0.5, color='#009E73')
        ax.set_xlabel('Length (bp)')
        ax.set_ylabel('Frequency')
    
    if title:
        # Replace underscores with spaces in title
        display_title = title.replace('_', ' ')
        ax.set_title(display_title, fontweight='bold', pad=10)
    
    # Apply Nature journal style
    _apply_nature_style(ax)
    
    plt.tight_layout()
    return fig


# =============================================================================
# COMPARISON PLOTS (Nature-Level Quality)
# =============================================================================

def plot_class_comparison(results: Dict[str, List[Dict[str, Any]]], 
                         metric: str = 'count',
                         title: Optional[str] = None,
                         figsize: Tuple[float, float] = None,
                         style: str = 'nature') -> plt.Figure:
    """
    Compare motif classes across multiple sequences/samples.
    
    Publication-quality heatmap comparison.
    
    Args:
        results: Dictionary of {sample_name: motifs_list}
        metric: Comparison metric ('count', 'coverage', 'density')
        title: Custom plot title
        figsize: Figure size (width, height) in inches
        style: Style preset
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style(style)
    
    if figsize is None:
        figsize = FIGURE_SIZES['double_column']
    
    if not results:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No data to display', ha='center', va='center', 
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title or 'Class Comparison')
        return fig
    
    # Collect all classes across samples
    all_classes = set()
    for motifs in results.values():
        all_classes.update(m.get('Class', 'Unknown') for m in motifs)
    all_classes = sorted(all_classes)
    
    # Calculate metrics for each sample
    comparison_data = []
    
    for sample_name, motifs in results.items():
        class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
        
        for class_name in all_classes:
            count = class_counts.get(class_name, 0)
            
            if metric == 'count':
                value = count
            elif metric == 'coverage':
                # Calculate coverage percentage (simplified)
                class_motifs = [m for m in motifs if m.get('Class') == class_name]
                covered_length = sum(m.get('Length', 0) for m in class_motifs)
                # Assume 10kb sequence for percentage calculation
                value = (covered_length / 10000) * 100
            elif metric == 'density':
                # Motifs per kb (assume 10kb sequence)
                value = count / 10
            else:
                value = count
            
            comparison_data.append({
                'Sample': sample_name,
                'Class': class_name,
                'Value': value
            })
    
    # Create DataFrame and pivot for heatmap
    df = pd.DataFrame(comparison_data)
    pivot_df = df.pivot(index='Sample', columns='Class', values='Value').fillna(0)
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=figsize)
    
    sns.heatmap(pivot_df, annot=True, fmt='.1f', cmap='YlOrRd', 
                ax=ax, cbar_kws={'label': f'Motif {metric.title()}'})
    
    ax.set_title(title or f'Motif Class Comparison by {metric.title()}')
    ax.set_xlabel('Motif Class')
    ax.set_ylabel('Sample')
    
    plt.tight_layout()
    return fig

# =============================================================================
# INTERACTIVE PLOTS (PLOTLY)
# =============================================================================

def create_interactive_coverage_plot(motifs: List[Dict[str, Any]], 
                                   sequence_length: int,
                                   title: str = "Interactive Motif Coverage") -> Union[Any, plt.Figure]:
    """
    Create interactive coverage plot using Plotly
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of the analyzed sequence
        title: Plot title
        
    Returns:
        Plotly figure if available, otherwise matplotlib figure
    """
    if not PLOTLY_AVAILABLE:
        return plot_coverage_map(motifs, sequence_length, title)
    
    if not motifs:
        fig = go.Figure()
        fig.add_annotation(text="No motifs to display", 
                          xref="paper", yref="paper",
                          x=0.5, y=0.5, showarrow=False)
        fig.update_layout(title=title)
        return fig
    
    fig = go.Figure()
    
    # Group motifs by class
    class_motifs = defaultdict(list)
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        class_motifs[class_name].append(motif)
    
    y_pos = 0
    for class_name, class_motif_list in class_motifs.items():
        color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        
        x_starts = []
        x_ends = []
        y_positions = []
        hover_texts = []
        
        for motif in class_motif_list:
            start = motif.get('Start', 0)
            end = motif.get('End', start + 1)
            
            x_starts.append(start)
            x_ends.append(end)
            y_positions.extend([y_pos - 0.4, y_pos - 0.4, y_pos + 0.4, y_pos + 0.4, None])
            
            hover_text = f"{class_name}<br>{motif.get('Subclass', '')}<br>" + \
                        f"Position: {start}-{end}<br>Length: {motif.get('Length', 0)} bp<br>" + \
                        f"Score: {motif.get('Score', 'N/A')}"
            hover_texts.append(hover_text)
        
        # Add rectangles for motifs
        for i, (start, end) in enumerate(zip(x_starts, x_ends)):
            fig.add_shape(
                type="rect",
                x0=start, y0=y_pos - 0.4,
                x1=end, y1=y_pos + 0.4,
                fillcolor=color,
                opacity=0.7,
                line=dict(color="black", width=1)
            )
            
            # Add invisible scatter points for hover
            fig.add_trace(go.Scatter(
                x=[(start + end) / 2],
                y=[y_pos],
                mode='markers',
                marker=dict(size=0.1, opacity=0),
                hoverinfo='text',
                hovertext=hover_texts[i],
                showlegend=False
            ))
        
        # Add class label
        fig.add_trace(go.Scatter(
            x=[sequence_length * 1.02],
            y=[y_pos],
            mode='text',
            text=[class_name],
            textposition="middle right",
            showlegend=False,
            hoverinfo='skip'
        ))
        
        y_pos += 1
    
    fig.update_layout(
        title=title,
        xaxis_title="Sequence Position (bp)",
        yaxis_title="Motif Class",
        xaxis=dict(range=[0, sequence_length * 1.15]),
        yaxis=dict(range=[-0.5, len(class_motifs) - 0.5], showticklabels=False),
        showlegend=False,
        height=max(400, len(class_motifs) * 60)
    )
    
    return fig

# =============================================================================
# BATCH EXPORT FUNCTIONS (Publication-Ready)
# =============================================================================

def save_all_plots(motifs: List[Dict[str, Any]], 
                   sequence_length: int,
                   output_dir: str = "plots",
                   file_format: str = "pdf",
                   dpi: int = 300,
                   style: str = 'nature') -> Dict[str, str]:
    """
    Generate and save all standard plots in publication quality.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of analyzed sequence
        output_dir: Output directory for plots
        file_format: File format ('pdf', 'png', 'svg') - pdf recommended for publication
        dpi: Resolution for raster formats (default 300 DPI for publication)
        style: Style preset ('nature', 'default')
        
    Returns:
        Dictionary of {plot_name: file_path}
    """
    import os
    
    # Apply style
    set_scientific_style(style)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    saved_files = {}
    
    # Use publication DPI if not specified
    if dpi < PUBLICATION_DPI:
        dpi = PUBLICATION_DPI
    
    # List of plots to generate (diverse, non-repetitive visualization types)
    plots_to_generate = [
        ("motif_distribution_class", lambda: plot_motif_distribution(motifs, by='Class', style=style)),
        ("coverage_map", lambda: plot_coverage_map(motifs, sequence_length, style=style)),
        ("density_heatmap", lambda: plot_density_heatmap(motifs, sequence_length, style=style)),
        ("score_distribution", lambda: plot_score_distribution(motifs, by_class=True, style=style)),
        ("length_distribution", lambda: plot_length_distribution(motifs, by_class=True, style=style)),
        ("nested_donut_chart", lambda: plot_nested_pie_chart(motifs))
    ]
    
    for plot_name, plot_func in plots_to_generate:
        try:
            fig = plot_func()
            filename = f"{plot_name}.{file_format}"
            filepath = os.path.join(output_dir, filename)
            
            if hasattr(fig, 'savefig'):  # Matplotlib figure
                # Publication-quality save options
                fig.savefig(filepath, format=file_format, dpi=dpi, 
                           bbox_inches='tight', pad_inches=0.05,
                           facecolor='white', edgecolor='none')
                plt.close(fig)
            else:  # Plotly figure
                if file_format.lower() == 'png':
                    fig.write_image(filepath, scale=2)  # Higher resolution
                elif file_format.lower() in ('pdf', 'svg'):
                    fig.write_image(filepath)
                elif file_format.lower() == 'html':
                    fig.write_html(filepath)
                else:
                    # Unsupported format - try default write_image
                    print(f"⚠ Unsupported format '{file_format}' for Plotly, attempting save anyway")
                    try:
                        fig.write_image(filepath)
                    except Exception as fmt_err:
                        raise ValueError(f"Unsupported file format: {file_format}. Use 'png', 'pdf', 'svg', or 'html'.") from fmt_err
            
            saved_files[plot_name] = filepath
            print(f"[OK] Saved {plot_name} to {filepath}")
            
        except Exception as e:
            print(f"✗ Failed to generate {plot_name}: {e}")
    
    return saved_files


# =============================================================================
# ENHANCED SCIENTIFIC VISUALIZATION FUNCTIONS
# =============================================================================

def plot_class_analysis_comprehensive(motifs: List[Dict[str, Any]], 
                                     figsize: Tuple[int, int] = (16, 12)) -> plt.Figure:
    """
    Comprehensive class-level analysis with multiple subplots.
    Shows distribution, statistics, and comparison of all 11 Non-B DNA classes.
    Highlights which classes were detected and which were not.
    
    Args:
        motifs: List of motif dictionaries
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object with multiple subplots
    """
    set_scientific_style()
    
    # Define all 11 Non-B DNA classes
    all_classes = [
        'Curved_DNA', 'Slipped_DNA', 'Cruciform', 'R-Loop', 'Triplex',
        'G-Quadruplex', 'i-Motif', 'Z-DNA', 'A-philic_DNA', 
        'Hybrid', 'Non-B_DNA_Clusters'
    ]
    
    # Count motifs by class
    detected_classes = Counter(m.get('Class', 'Unknown') for m in motifs)
    
    # Identify detected vs not detected
    detected = [cls for cls in all_classes if detected_classes.get(cls, 0) > 0]
    not_detected = [cls for cls in all_classes if detected_classes.get(cls, 0) == 0]
    
    # Create figure with subplots
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)
    
    # 1. Main distribution bar chart
    ax1 = fig.add_subplot(gs[0, :])
    counts = [detected_classes.get(cls, 0) for cls in all_classes]
    colors = [MOTIF_CLASS_COLORS.get(cls, '#808080') for cls in all_classes]
    bars = ax1.bar(range(len(all_classes)), counts, color=colors, alpha=0.8, 
                   edgecolor='black', linewidth=1.5)
    
    ax1.set_xlabel('Non-B DNA Class', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax1.set_title('Distribution of All 11 Non-B DNA Classes', fontsize=14, fontweight='bold')
    ax1.set_xticks(range(len(all_classes)))
    ax1.set_xticklabels(all_classes, rotation=45, ha='right', fontsize=10)
    ax1.grid(axis='y', alpha=0.3)
    
    # Add count labels on bars
    for bar, count in zip(bars, counts):
        if count > 0:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{count}', ha='center', va='bottom', fontweight='bold', fontsize=9)
    
    # 2. Detected vs Not Detected pie chart
    ax2 = fig.add_subplot(gs[1, 0])
    detection_counts = [len(detected), len(not_detected)]
    detection_labels = [f'Detected\n({len(detected)} classes)', 
                       f'Not Detected\n({len(not_detected)} classes)']
    colors_pie = ['#4CAF50', '#FF5722']
    ax2.pie(detection_counts, labels=detection_labels, autopct='%1.1f%%',
            colors=colors_pie, startangle=90, textprops={'fontsize': 11, 'fontweight': 'bold'})
    ax2.set_title('Class Detection Status', fontsize=12, fontweight='bold')
    
    # 3. Statistics table for detected classes
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axis('off')
    
    if detected:
        # Calculate statistics for detected classes
        class_stats = []
        for cls in detected[:5]:  # Show top 5
            cls_motifs = [m for m in motifs if m.get('Class') == cls]
            count = len(cls_motifs)
            avg_length = np.mean([m.get('Length', 0) for m in cls_motifs]) if cls_motifs else 0
            avg_score = np.mean([m.get('Score', 0) for m in cls_motifs]) if cls_motifs else 0
            class_stats.append([cls[:15], count, f'{avg_length:.1f}', f'{avg_score:.3f}'])
        
        table = ax3.table(cellText=class_stats,
                         colLabels=['Class', 'Count', 'Avg Len', 'Avg Score'],
                         cellLoc='left', loc='center',
                         colWidths=[0.4, 0.2, 0.2, 0.2])
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 2)
        
        # Style header
        for i in range(4):
            table[(0, i)].set_facecolor('#E0E0E0')
            table[(0, i)].set_text_props(weight='bold')
        
        ax3.set_title('Top Detected Classes (Statistics)', fontsize=12, fontweight='bold', pad=20)
    else:
        ax3.text(0.5, 0.5, 'No classes detected', ha='center', va='center',
                transform=ax3.transAxes, fontsize=12)
    
    # 4. List of not detected classes
    ax4 = fig.add_subplot(gs[2, :])
    ax4.axis('off')
    
    if not_detected:
        not_detected_text = 'Classes NOT Detected:\n' + ', '.join(not_detected)
        ax4.text(0.5, 0.5, not_detected_text, ha='center', va='center',
                transform=ax4.transAxes, fontsize=11, 
                bbox=dict(boxstyle='round', facecolor='#FFEBEE', alpha=0.8),
                wrap=True)
    else:
        ax4.text(0.5, 0.5, 'All 11 Non-B DNA classes detected! [OK]', 
                ha='center', va='center', transform=ax4.transAxes, fontsize=12,
                fontweight='bold', color='green')
    
    plt.suptitle(f'Comprehensive Class Analysis ({len(motifs)} total motifs)', 
                fontsize=16, fontweight='bold', y=0.98)
    
    return fig


def plot_subclass_analysis_comprehensive(motifs: List[Dict[str, Any]], 
                                        figsize: Tuple[int, int] = (18, 14)) -> plt.Figure:
    """
    Comprehensive subclass-level analysis showing all detected subclasses
    organized by their parent class.
    
    Args:
        motifs: List of motif dictionaries
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object with subclass analysis
    """
    set_scientific_style()
    
    # Group motifs by class and subclass
    class_subclass_counts = defaultdict(lambda: defaultdict(int))
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        subclass_name = motif.get('Subclass', 'Unknown')
        class_subclass_counts[class_name][subclass_name] += 1
    
    # Prepare data for visualization
    all_subclasses = []
    all_counts = []
    all_classes = []
    
    for class_name in sorted(class_subclass_counts.keys()):
        for subclass_name, count in sorted(class_subclass_counts[class_name].items()):
            all_subclasses.append(f"{class_name}:{subclass_name}")
            all_counts.append(count)
            all_classes.append(class_name)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, 
                                   gridspec_kw={'height_ratios': [2, 1]})
    
    # 1. Subclass distribution bar chart
    colors = [MOTIF_CLASS_COLORS.get(cls, '#808080') for cls in all_classes]
    x_pos = range(len(all_subclasses))
    bars = ax1.barh(x_pos, all_counts, color=colors, alpha=0.8, 
                    edgecolor='black', linewidth=0.5)
    
    ax1.set_yticks(x_pos)
    ax1.set_yticklabels(all_subclasses, fontsize=9)
    ax1.set_xlabel('Count', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Class:Subclass', fontsize=12, fontweight='bold')
    ax1.set_title('Distribution of All Detected Subclasses', fontsize=14, fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)
    
    # Add count labels on bars
    for bar, count in zip(bars, all_counts):
        width = bar.get_width()
        ax1.text(width, bar.get_y() + bar.get_height()/2.,
                f' {count}', ha='left', va='center', fontweight='bold', fontsize=8)
    
    # 2. Subclass summary by class
    ax2.axis('off')
    
    # Create summary text
    summary_lines = ['Subclass Summary by Class:\n']
    for class_name in sorted(class_subclass_counts.keys()):
        subclasses = class_subclass_counts[class_name]
        n_subclasses = len(subclasses)
        total_count = sum(subclasses.values())
        summary_lines.append(f'{class_name}: {n_subclasses} subclass(es), {total_count} motifs')
    
    summary_text = '\n'.join(summary_lines)
    ax2.text(0.1, 0.5, summary_text, ha='left', va='center',
            transform=ax2.transAxes, fontsize=10, family='monospace',
            bbox=dict(boxstyle='round', facecolor='#F5F5F5', alpha=0.8))
    
    plt.suptitle(f'Comprehensive Subclass Analysis ({len(motifs)} total motifs)', 
                fontsize=16, fontweight='bold', y=0.99)
    plt.tight_layout()
    
    return fig


def plot_score_statistics_by_class(motifs: List[Dict[str, Any]], 
                                   figsize: Tuple[int, int] = (14, 8)) -> plt.Figure:
    """
    Advanced statistical visualization of scores by class.
    Shows box plots, violin plots, and statistical annotations.
    
    Args:
        motifs: List of motif dictionaries
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object with score statistics
    """
    set_scientific_style()
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No motifs to analyze', ha='center', va='center',
               transform=ax.transAxes, fontsize=14)
        return fig
    
    # Prepare data
    df_data = []
    for motif in motifs:
        df_data.append({
            'Class': motif.get('Class', 'Unknown'),
            'Score': motif.get('Score', 0),
            'Length': motif.get('Length', 0)
        })
    df = pd.DataFrame(df_data)
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, 
                                   gridspec_kw={'height_ratios': [2, 1]})
    
    # 1. Violin plot with box plot overlay
    classes = sorted(df['Class'].unique())
    positions = range(len(classes))
    
    # Create violin plot
    parts = ax1.violinplot([df[df['Class'] == cls]['Score'].values for cls in classes],
                          positions=positions, widths=0.7, showmeans=True, showmedians=True)
    
    # Color violins by class
    for i, pc in enumerate(parts['bodies']):
        cls = classes[i]
        color = MOTIF_CLASS_COLORS.get(cls, '#808080')
        pc.set_facecolor(color)
        pc.set_alpha(0.6)
    
    # Overlay box plots
    bp = ax1.boxplot([df[df['Class'] == cls]['Score'].values for cls in classes],
                     positions=positions, widths=0.3, patch_artist=True,
                     boxprops=dict(facecolor='white', alpha=0.7),
                     medianprops=dict(color='red', linewidth=2))
    
    ax1.set_xticks(positions)
    ax1.set_xticklabels(classes, rotation=45, ha='right', fontsize=10)
    ax1.set_ylabel('Score', fontsize=12, fontweight='bold')
    ax1.set_title('Score Distribution by Class (Violin + Box Plot)', fontsize=14, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3)
    
    # Add statistical annotations
    for i, cls in enumerate(classes):
        scores = df[df['Class'] == cls]['Score'].values
        if len(scores) > 0:
            mean_score = np.mean(scores)
            std_score = np.std(scores)
            ax1.text(i, ax1.get_ylim()[1] * 0.95, 
                    f'μ={mean_score:.2f}\nσ={std_score:.2f}',
                    ha='center', va='top', fontsize=8,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    # 2. Statistical summary table
    ax2.axis('off')
    
    # Calculate statistics
    stats_data = []
    for cls in classes:
        cls_scores = df[df['Class'] == cls]['Score'].values
        if len(cls_scores) > 0:
            stats_data.append([
                cls[:15],
                len(cls_scores),
                f'{np.mean(cls_scores):.3f}',
                f'{np.median(cls_scores):.3f}',
                f'{np.std(cls_scores):.3f}',
                f'{np.min(cls_scores):.3f}',
                f'{np.max(cls_scores):.3f}'
            ])
    
    if stats_data:
        table = ax2.table(cellText=stats_data,
                         colLabels=['Class', 'N', 'Mean', 'Median', 'Std', 'Min', 'Max'],
                         cellLoc='center', loc='center',
                         colWidths=[0.25, 0.1, 0.13, 0.13, 0.13, 0.13, 0.13])
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 1.8)
        
        # Style header
        for i in range(7):
            table[(0, i)].set_facecolor('#E0E0E0')
            table[(0, i)].set_text_props(weight='bold')
    
    plt.suptitle('Score Statistics Analysis', fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    return fig


def plot_length_statistics_by_class(motifs: List[Dict[str, Any]], 
                                   figsize: Tuple[int, int] = (14, 10)) -> plt.Figure:
    """
    Advanced visualization of motif length distributions by class.
    
    Args:
        motifs: List of motif dictionaries
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object with length statistics
    """
    set_scientific_style()
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No motifs to analyze', ha='center', va='center',
               transform=ax.transAxes, fontsize=14)
        return fig
    
    # Prepare data
    df_data = []
    for motif in motifs:
        df_data.append({
            'Class': motif.get('Class', 'Unknown'),
            'Length': motif.get('Length', 0)
        })
    df = pd.DataFrame(df_data)
    
    # Create figure with subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=figsize, 
                                       gridspec_kw={'height_ratios': [2, 1.5, 1]})
    
    # 1. Histogram with KDE overlay for each class
    classes = sorted(df['Class'].unique())
    
    for cls in classes:
        cls_lengths = df[df['Class'] == cls]['Length'].values
        if len(cls_lengths) > 1:
            color = MOTIF_CLASS_COLORS.get(cls, '#808080')
            ax1.hist(cls_lengths, bins=20, alpha=0.4, label=cls, color=color, edgecolor='black')
    
    ax1.set_xlabel('Length (bp)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax1.set_title('Length Distribution by Class', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=9, ncol=2)
    ax1.grid(axis='y', alpha=0.3)
    
    # 2. Box plot comparison
    bp = ax2.boxplot([df[df['Class'] == cls]['Length'].values for cls in classes],
                     labels=classes, patch_artist=True, vert=True)
    
    # Color boxes by class
    for i, (patch, cls) in enumerate(zip(bp['boxes'], classes)):
        color = MOTIF_CLASS_COLORS.get(cls, '#808080')
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    
    ax2.set_xticklabels(classes, rotation=45, ha='right', fontsize=10)
    ax2.set_ylabel('Length (bp)', fontsize=12, fontweight='bold')
    ax2.set_title('Length Comparison (Box Plot)', fontsize=12, fontweight='bold')
    ax2.grid(axis='y', alpha=0.3)
    
    # 3. Statistical summary table
    ax3.axis('off')
    
    stats_data = []
    for cls in classes:
        cls_lengths = df[df['Class'] == cls]['Length'].values
        if len(cls_lengths) > 0:
            stats_data.append([
                cls[:15],
                len(cls_lengths),
                f'{np.mean(cls_lengths):.1f}',
                f'{np.median(cls_lengths):.1f}',
                f'{np.std(cls_lengths):.1f}',
                f'{np.min(cls_lengths):.0f}',
                f'{np.max(cls_lengths):.0f}'
            ])
    
    if stats_data:
        table = ax3.table(cellText=stats_data,
                         colLabels=['Class', 'N', 'Mean', 'Median', 'Std', 'Min', 'Max'],
                         cellLoc='center', loc='center',
                         colWidths=[0.25, 0.1, 0.13, 0.13, 0.13, 0.13, 0.13])
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 1.8)
        
        # Style header
        for i in range(7):
            table[(0, i)].set_facecolor('#E0E0E0')
            table[(0, i)].set_text_props(weight='bold')
    
    plt.suptitle('Length Statistics Analysis', fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    return fig


# =============================================================================
# TESTING & EXAMPLES
# =============================================================================

def test_visualizations():
    """Test visualization functions with example data"""
    print("Testing NBDScanner visualizations...")
    
    # Create example motif data
    example_motifs = [
        {'Class': 'G-Quadruplex', 'Subclass': 'Canonical G4', 'Start': 1, 'End': 21, 'Length': 21, 'Score': 0.85},
        {'Class': 'G-Quadruplex', 'Subclass': 'Relaxed G4', 'Start': 45, 'End': 60, 'Length': 16, 'Score': 0.72},
        {'Class': 'Curved_DNA', 'Subclass': 'A-tract', 'Start': 80, 'End': 95, 'Length': 16, 'Score': 0.65},
        {'Class': 'Z-DNA', 'Subclass': 'CG alternating', 'Start': 120, 'End': 135, 'Length': 16, 'Score': 0.90},
        {'Class': 'i-Motif', 'Subclass': 'Canonical i-motif', 'Start': 160, 'End': 180, 'Length': 21, 'Score': 0.78}
    ]
    
    sequence_length = 200
    
    print(f"\nTesting with {len(example_motifs)} example motifs:")
    for motif in example_motifs:
        print(f"  {motif['Class']} at {motif['Start']}-{motif['End']}")
    
    # Test basic plots
    try:
        fig1 = plot_motif_distribution(example_motifs, by='Class')
        plt.close(fig1)
        print("[OK] Motif distribution plot: PASS")
    except Exception as e:
        print(f"✗ Motif distribution plot: FAIL - {e}")
    
    try:
        fig2 = plot_coverage_map(example_motifs, sequence_length)
        plt.close(fig2)
        print("[OK] Coverage map plot: PASS")
    except Exception as e:
        print(f"✗ Coverage map plot: FAIL - {e}")
    
    try:
        fig3 = plot_score_distribution(example_motifs)
        plt.close(fig3)
        print("[OK] Score distribution plot: PASS")
    except Exception as e:
        print(f"✗ Score distribution plot: FAIL - {e}")
    
    try:
        fig4 = plot_nested_pie_chart(example_motifs)
        plt.close(fig4)
        print("[OK] Nested pie chart: PASS")
    except Exception as e:
        print(f"✗ Nested pie chart: FAIL - {e}")
    
    print(f"\n[OK] Visualization testing completed")
    print(f"Plotly available: {'Yes' if PLOTLY_AVAILABLE else 'No'}")


# =============================================================================
# ENHANCED STATISTICS VISUALIZATIONS: DENSITY AND ENRICHMENT
# =============================================================================

def plot_density_comparison(genomic_density: Dict[str, float],
                            positional_density: Dict[str, float],
                            title: str = "Motif Density Analysis",
                            figsize: Tuple[int, int] = (14, 6)) -> plt.Figure:
    """
    Plot comparison of genomic density (coverage %) and positional density (motifs/kbp).
    
    Args:
        genomic_density: Dictionary of class -> genomic density (%)
        positional_density: Dictionary of class -> positional density (motifs/unit)
        title: Plot title
        figsize: Figure size
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    # Remove 'Overall' for class-specific comparison
    classes = [k for k in genomic_density.keys() if k != 'Overall']
    if not classes:
        classes = list(genomic_density.keys())
    
    # Sort classes alphabetically
    classes = sorted(classes)
    
    genomic_vals = [genomic_density.get(c, 0) for c in classes]
    positional_vals = [positional_density.get(c, 0) for c in classes]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, dpi=PUBLICATION_DPI)
    
    # Genomic Density (Coverage %)
    colors1 = [MOTIF_CLASS_COLORS.get(c, '#808080') for c in classes]
    
    # Replace underscores with spaces in labels
    display_classes = [c.replace('_', ' ') for c in classes]
    
    bars1 = ax1.barh(display_classes, genomic_vals, color=colors1, alpha=0.8, 
                     edgecolor='black', linewidth=0.5)
    ax1.set_xlabel('Genomic Density (Coverage %)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Motif Class', fontsize=11, fontweight='bold')
    ax1.set_title('A. Genomic Density (σ_G)', fontsize=12, fontweight='bold', pad=10)
    
    # Add value labels with safe max calculation
    max_genomic = max(genomic_vals) if genomic_vals and max(genomic_vals) > 0 else 1
    for i, (bar, val) in enumerate(zip(bars1, genomic_vals)):
        if val > 0:
            ax1.text(val + max_genomic * 0.01, i, f'{val:.3f}%', 
                    va='center', fontsize=9, fontweight='bold')
    
    # Positional Density (Frequency)
    colors2 = [MOTIF_CLASS_COLORS.get(c, '#808080') for c in classes]
    bars2 = ax2.barh(display_classes, positional_vals, color=colors2, alpha=0.8, 
                     edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('Positional Density (motifs/kbp)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Motif Class', fontsize=11, fontweight='bold')
    ax2.set_title('B. Positional Density (λ)', fontsize=12, fontweight='bold', pad=10)
    
    # Add value labels with safe max calculation
    max_positional = max(positional_vals) if positional_vals and max(positional_vals) > 0 else 1
    for i, (bar, val) in enumerate(zip(bars2, positional_vals)):
        if val > 0:
            ax2.text(val + max_positional * 0.01, i, f'{val:.2f}', 
                    va='center', fontsize=9, fontweight='bold')
    
    # Apply Nature journal style
    _apply_nature_style(ax1)
    _apply_nature_style(ax2)
    
    plt.suptitle(title, fontsize=14, fontweight='bold', y=1.00)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    
    return fig


def plot_enrichment_analysis(enrichment_results: Dict[str, Dict[str, Any]],
                             title: str = "Motif Enrichment Analysis",
                             figsize: Tuple[int, int] = (14, 8)) -> plt.Figure:
    """
    Plot enrichment analysis results with fold enrichment and p-values.
    
    Args:
        enrichment_results: Dictionary from calculate_enrichment_with_shuffling
        title: Plot title
        figsize: Figure size
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    # Extract data (exclude 'Overall' for class-specific view)
    classes = [k for k in enrichment_results.keys() if k != 'Overall']
    if not classes:
        classes = list(enrichment_results.keys())
    
    fold_enrichments = []
    p_values = []
    observed_densities = []
    background_means = []
    
    for cls in classes:
        result = enrichment_results[cls]
        fe = result.get('fold_enrichment', 0)
        # Handle both string 'Inf' and float infinity values robustly
        if fe == 'Inf' or (isinstance(fe, float) and np.isinf(fe)):
            fe = INFINITE_FOLD_ENRICHMENT_CAP  # Cap infinite values for visualization
        fold_enrichments.append(fe)
        p_values.append(result.get('p_value', 1.0))
        observed_densities.append(result.get('observed_density', 0))
        background_means.append(result.get('background_mean', 0))
    
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)
    
    # 1. Fold Enrichment
    ax1 = fig.add_subplot(gs[0, 0])
    colors = [MOTIF_CLASS_COLORS.get(c, '#808080') for c in classes]
    bars1 = ax1.barh(classes, fold_enrichments, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax1.axvline(x=1.0, color='red', linestyle='--', linewidth=2, label='No enrichment (FE=1)')
    ax1.set_xlabel('Fold Enrichment', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Motif Class', fontsize=11, fontweight='bold')
    ax1.set_title('A. Fold Enrichment', fontsize=12, fontweight='bold')
    ax1.legend(loc='best', fontsize=9)
    
    # Add value labels
    for i, (bar, val) in enumerate(zip(bars1, fold_enrichments)):
        label_text = f'{val:.2f}' if val < INFINITE_FOLD_ENRICHMENT_CAP else 'Inf'
        ax1.text(val + max(fold_enrichments) * 0.01, i, label_text, 
                va='center', fontsize=9, fontweight='bold')
    
    # 2. P-values
    ax2 = fig.add_subplot(gs[0, 1])
    # Color code by significance
    p_colors = ['green' if p < 0.05 else 'orange' if p < 0.1 else 'red' for p in p_values]
    bars2 = ax2.barh(classes, p_values, color=p_colors, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax2.axvline(x=0.05, color='green', linestyle='--', linewidth=2, label='p=0.05')
    ax2.axvline(x=0.1, color='orange', linestyle='--', linewidth=1.5, label='p=0.1')
    ax2.set_xlabel('P-value', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Motif Class', fontsize=11, fontweight='bold')
    ax2.set_title('B. Statistical Significance', fontsize=12, fontweight='bold')
    ax2.set_xlim(0, max(1.0, max(p_values) * 1.1))
    ax2.legend(loc='best', fontsize=9)
    
    # Add value labels
    for i, (bar, val) in enumerate(zip(bars2, p_values)):
        ax2.text(val + 0.02, i, f'{val:.3f}', va='center', fontsize=9, fontweight='bold')
    
    # 3. Observed vs Background Density
    ax3 = fig.add_subplot(gs[1, :])
    x = np.arange(len(classes))
    width = 0.35
    
    bars3a = ax3.bar(x - width/2, observed_densities, width, label='Observed', 
                    color='steelblue', alpha=0.8, edgecolor='black', linewidth=0.5)
    bars3b = ax3.bar(x + width/2, background_means, width, label='Background (Mean)', 
                    color='coral', alpha=0.8, edgecolor='black', linewidth=0.5)
    
    ax3.set_xlabel('Motif Class', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Density (%)', fontsize=11, fontweight='bold')
    ax3.set_title('C. Observed vs. Background Density Comparison', fontsize=12, fontweight='bold')
    ax3.set_xticks(x)
    ax3.set_xticklabels(classes, rotation=45, ha='right', fontsize=9)
    ax3.legend(loc='best', fontsize=10)
    ax3.grid(axis='y', alpha=0.3)
    
    plt.suptitle(title, fontsize=14, fontweight='bold', y=0.99)
    
    return fig


def plot_enrichment_summary_table(enrichment_results: Dict[str, Dict[str, Any]],
                                  title: str = "Enrichment Summary Statistics") -> plt.Figure:
    """
    Create a summary table visualization for enrichment results.
    
    Args:
        enrichment_results: Dictionary from calculate_enrichment_with_shuffling
        title: Plot title
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    # Prepare data for table
    classes = [k for k in enrichment_results.keys() if k != 'Overall']
    if not classes:
        return None
    
    table_data = []
    for cls in classes:
        result = enrichment_results[cls]
        fe = result.get('fold_enrichment', 0)
        fe_str = f"{fe:.2f}" if fe != 'Inf' else 'Inf'
        
        row = [
            cls,
            result.get('observed_count', 0),
            f"{result.get('observed_density', 0):.4f}%",
            f"{result.get('background_mean', 0):.4f}%",
            fe_str,
            f"{result.get('p_value', 1.0):.4f}",
            '***' if result.get('p_value', 1.0) < 0.001 else 
            '**' if result.get('p_value', 1.0) < 0.01 else 
            '*' if result.get('p_value', 1.0) < 0.05 else 'ns'
        ]
        table_data.append(row)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, max(6, len(classes) * 0.4)))
    ax.axis('tight')
    ax.axis('off')
    
    # Create table
    headers = ['Class', 'Count', 'Observed\nDensity', 'Background\nMean', 
               'Fold\nEnrichment', 'P-value', 'Sig.']
    
    table = ax.table(cellText=table_data, colLabels=headers, 
                    cellLoc='center', loc='center',
                    colWidths=[0.20, 0.10, 0.15, 0.15, 0.15, 0.12, 0.08])
    
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    # Style header
    for i, header in enumerate(headers):
        cell = table[(0, i)]
        cell.set_facecolor('#2196F3')
        cell.set_text_props(weight='bold', color='white', fontsize=11)
    
    # Style rows with alternating colors
    for i in range(1, len(table_data) + 1):
        for j in range(len(headers)):
            cell = table[(i, j)]
            if i % 2 == 0:
                cell.set_facecolor('#F5F5F5')
            else:
                cell.set_facecolor('white')
            
            # Highlight significant results
            if j == 6:  # Significance column
                if table_data[i-1][j] in ['***', '**', '*']:
                    cell.set_facecolor('#C8E6C9')
                    cell.set_text_props(weight='bold', color='green')
    
    plt.title(title, fontsize=14, fontweight='bold', pad=20)
    plt.tight_layout()
    
    return fig


# =============================================================================
# CIRCOS PLOT FOR NON-B DNA MOTIF DENSITY
# =============================================================================

# Motif classes to exclude from Circos visualization (dynamic classes)
CIRCOS_EXCLUDED_CLASSES = ['Hybrid', 'Non-B_DNA_Clusters']


def plot_circos_motif_density(motifs: List[Dict[str, Any]], 
                               sequence_length: int,
                               title: str = "Non-B DNA Motif Density Circos Plot",
                               figsize: Tuple[int, int] = (12, 12),
                               window_size: int = None) -> plt.Figure:
    """
    Create a circular Circos-style plot showing non-B DNA motif class density.
    
    The plot shows:
    - Outer ring: Sequence position ruler
    - Inner rings: One ring per motif class showing density
    - Center: Summary statistics
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total length of analyzed sequence in bp
        title: Plot title
        figsize: Figure size (width, height)
        window_size: Window size for density calculation (auto-calculated if None)
        
    Returns:
        Matplotlib figure object with Circos-style visualization
    """
    set_scientific_style()
    
    if not motifs or sequence_length == 0:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    # Auto-calculate window size if not provided
    if window_size is None:
        window_size = max(100, sequence_length // 50)
    
    # Calculate number of windows
    num_windows = max(1, sequence_length // window_size)
    
    # Get unique classes (excluding Hybrid and Clusters for cleaner visualization)
    classes = sorted(set(m.get('Class', 'Unknown') for m in motifs 
                        if m.get('Class') not in CIRCOS_EXCLUDED_CLASSES))
    
    if not classes:
        classes = sorted(set(m.get('Class', 'Unknown') for m in motifs))
    
    # Calculate density per window per class
    class_densities = {}
    for class_name in classes:
        densities = []
        class_motifs = [m for m in motifs if m.get('Class') == class_name]
        
        for i in range(num_windows):
            window_start = i * window_size
            window_end = (i + 1) * window_size
            
            # Count motifs in this window
            count = 0
            for motif in class_motifs:
                motif_start = motif.get('Start', 0) - 1  # 0-based
                motif_end = motif.get('End', 0)
                # Check overlap with window
                if not (motif_end <= window_start or motif_start >= window_end):
                    count += 1
            
            # Convert to density (motifs per kb)
            window_kb = window_size / 1000
            densities.append(count / window_kb if window_kb > 0 else 0)
        
        class_densities[class_name] = densities
    
    # Create figure
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='polar')
    
    # Calculate angles for each window
    theta = np.linspace(0, 2 * np.pi, num_windows, endpoint=False)
    
    # Width of each bar
    width = 2 * np.pi / num_windows * 0.8
    
    # Ring configuration
    ring_width = 0.12
    inner_radius = 0.3
    
    # Plot each class as a ring
    for i, class_name in enumerate(classes):
        densities = class_densities[class_name]
        
        # Normalize densities for this class (0-1 scale for bar height)
        max_density = max(densities) if max(densities) > 0 else 1
        normalized = [d / max_density for d in densities]
        
        # Calculate ring position
        ring_bottom = inner_radius + i * ring_width
        
        # Get color for this class
        color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        
        # Plot bars for this ring
        heights = [n * ring_width * 0.9 for n in normalized]
        # Replace underscores with spaces in legend label
        display_name = class_name.replace('_', ' ')
        bars = ax.bar(theta, heights, width=width, bottom=ring_bottom,
                     color=color, alpha=0.7, edgecolor='white', linewidth=0.5,
                     label=f'{display_name} (max: {max_density:.1f}/kb)')
    
    # Add outer position ruler
    outer_radius = inner_radius + len(classes) * ring_width + 0.05
    ruler_theta = np.linspace(0, 2 * np.pi, 12, endpoint=False)
    ruler_labels = [f'{int(i * sequence_length / 12 / 1000)}kb' for i in range(12)]
    
    ax.set_xticks(ruler_theta)
    ax.set_xticklabels(ruler_labels, fontsize=9, fontweight='bold')
    
    # Remove radial labels
    ax.set_yticklabels([])
    
    # Set limits
    ax.set_ylim(0, outer_radius + 0.1)
    
    # Add legend
    ax.legend(loc='center', bbox_to_anchor=(0.5, 0.5), fontsize=8, 
             framealpha=0.9, ncol=1)
    
    # Add title (replace underscores with spaces)
    display_title = title.replace('_', ' ')
    fig.suptitle(display_title, fontsize=14, fontweight='bold', y=0.98)
    
    # Add center statistics
    total_motifs = len([m for m in motifs if m.get('Class') not in CIRCOS_EXCLUDED_CLASSES])
    center_text = f"Total: {total_motifs}\n{len(classes)} classes\n{sequence_length/1000:.1f} kb"
    ax.text(0, 0, center_text, ha='center', va='center', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    return fig


def plot_radial_class_density(motifs: List[Dict[str, Any]], 
                               sequence_length: int,
                               title: str = "Radial Motif Class Density",
                               figsize: Tuple[int, int] = (10, 10)) -> plt.Figure:
    """
    Create a radial bar chart showing motif density per class.
    
    A simpler alternative to full Circos plot, showing aggregate density
    per motif class in a radial/polar layout.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total length of analyzed sequence in bp
        title: Plot title
        figsize: Figure size
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    if not motifs or sequence_length == 0:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    # Calculate density per class (motifs per kb)
    sequence_kb = sequence_length / 1000
    class_counts = Counter(m.get('Class', 'Unknown') for m in motifs
                          if m.get('Class') not in CIRCOS_EXCLUDED_CLASSES)
    
    if not class_counts:
        class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    
    classes = list(class_counts.keys())
    densities = [class_counts[c] / sequence_kb for c in classes]
    
    # Create polar plot
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='polar')
    
    # Calculate angles
    num_classes = len(classes)
    theta = np.linspace(0, 2 * np.pi, num_classes, endpoint=False)
    width = 2 * np.pi / num_classes * 0.7
    
    # Get colors
    colors = [MOTIF_CLASS_COLORS.get(c, '#808080') for c in classes]
    
    # Plot bars
    bars = ax.bar(theta, densities, width=width, color=colors, alpha=0.8,
                 edgecolor='white', linewidth=2)
    
    # Add class labels (replace underscores with spaces)
    ax.set_xticks(theta)
    display_classes = [cls.replace('_', ' ') for cls in classes]
    ax.set_xticklabels(display_classes, fontsize=10, fontweight='bold')
    
    # Add value labels on bars
    for angle, density, bar in zip(theta, densities, bars):
        ax.text(angle, density + max(densities) * 0.05, f'{density:.1f}',
               ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # Style (replace underscores with spaces in title)
    ax.set_ylabel('Density (motifs/kb)', labelpad=30, fontsize=11, fontweight='bold')
    display_title = title.replace('_', ' ')
    ax.set_title(display_title, fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    return fig


def plot_stacked_density_track(motifs: List[Dict[str, Any]], 
                                sequence_length: int,
                                title: str = "Stacked Motif Density Track",
                                figsize: Tuple[int, int] = (14, 6),
                                window_size: int = None) -> plt.Figure:
    """
    Create a stacked area chart showing motif density along the sequence.
    
    A linear (non-circular) alternative to Circos plot showing how different
    motif classes are distributed along the sequence.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total length of analyzed sequence in bp
        title: Plot title
        figsize: Figure size
        window_size: Window size for density calculation (auto if None)
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    if not motifs or sequence_length == 0:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    # Auto-calculate window size
    if window_size is None:
        window_size = max(100, sequence_length // 100)
    
    num_windows = max(1, sequence_length // window_size)
    
    # Get unique classes
    classes = sorted(set(m.get('Class', 'Unknown') for m in motifs
                        if m.get('Class') not in CIRCOS_EXCLUDED_CLASSES))
    
    if not classes:
        classes = sorted(set(m.get('Class', 'Unknown') for m in motifs))
    
    # Calculate density per window per class
    positions = np.arange(num_windows) * window_size / 1000  # In kb
    class_densities = {}
    
    for class_name in classes:
        densities = []
        class_motifs = [m for m in motifs if m.get('Class') == class_name]
        
        for i in range(num_windows):
            window_start = i * window_size
            window_end = (i + 1) * window_size
            
            count = 0
            for motif in class_motifs:
                motif_start = motif.get('Start', 0) - 1
                motif_end = motif.get('End', 0)
                if not (motif_end <= window_start or motif_start >= window_end):
                    count += 1
            
            densities.append(count)
        
        class_densities[class_name] = densities
    
    # Create stacked area plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Stack the densities
    colors = [MOTIF_CLASS_COLORS.get(c, '#808080') for c in classes]
    
    # Create arrays for stacking (replace underscores with spaces in labels)
    density_arrays = [np.array(class_densities[c]) for c in classes]
    display_classes = [cls.replace('_', ' ') for cls in classes]
    
    ax.stackplot(positions, *density_arrays, labels=display_classes, colors=colors, alpha=0.8)
    
    # Styling (replace underscores with spaces in title)
    ax.set_xlabel('Position (kb)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Motif Count per Window', fontsize=12, fontweight='bold')
    display_title = title.replace('_', ' ')
    ax.set_title(display_title, fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9, framealpha=0.9)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    return fig


# =============================================================================
# SUBCLASS-LEVEL DENSITY AND ENRICHMENT VISUALIZATIONS
# =============================================================================

def plot_density_comparison_by_subclass(genomic_density: Dict[str, float],
                                        positional_density: Dict[str, float],
                                        title: str = "Motif Density Analysis (by Subclass)",
                                        figsize: Tuple[int, int] = (16, 10)) -> plt.Figure:
    """
    Plot comparison of genomic density (coverage %) and positional density (motifs/kbp)
    at the subclass level.
    
    Args:
        genomic_density: Dictionary of 'Class:Subclass' -> genomic density (%)
        positional_density: Dictionary of 'Class:Subclass' -> positional density (motifs/unit)
        title: Plot title
        figsize: Figure size
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    # Remove 'Overall' for subclass-specific comparison
    subclasses = [k for k in genomic_density.keys() if k != 'Overall' and ':' in k]
    if not subclasses:
        # Fallback to regular class-level if no subclass data
        subclasses = [k for k in genomic_density.keys() if k != 'Overall']
    
    # Sort by class, then by subclass
    subclasses.sort()
    
    genomic_vals = [genomic_density.get(c, 0) for c in subclasses]
    positional_vals = [positional_density.get(c, 0) for c in subclasses]
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, dpi=PUBLICATION_DPI)
    
    # Genomic Density (Coverage %)
    y_pos = np.arange(len(subclasses))
    
    # Get colors based on parent class
    colors1 = []
    for subclass in subclasses:
        if ':' in subclass:
            parent_class = subclass.split(':')[0]
            colors1.append(MOTIF_CLASS_COLORS.get(parent_class, '#808080'))
        else:
            colors1.append(MOTIF_CLASS_COLORS.get(subclass, '#808080'))
    
    bars1 = ax1.barh(y_pos, genomic_vals, color=colors1, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax1.set_xlabel('Genomic Density (Coverage %)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Motif Subclass', fontsize=11, fontweight='bold')
    ax1.set_title('A. Genomic Density (σ_G) by Subclass', fontsize=12, fontweight='bold', pad=10)
    ax1.set_yticks(y_pos)
    
    # Replace underscores and format labels
    display_labels = [label.replace('_', ' ').replace(':', ': ') for label in subclasses]
    ax1.set_yticklabels(display_labels, fontsize=8)
    
    # Add value labels
    max_val = max(genomic_vals) if genomic_vals else 1
    for i, (bar, val) in enumerate(zip(bars1, genomic_vals)):
        if val > 0:
            ax1.text(val + max_val * 0.01, i, f'{val:.3f}%', 
                    va='center', fontsize=7, fontweight='bold')
    
    # Positional Density (Frequency)
    bars2 = ax2.barh(y_pos, positional_vals, color=colors1, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('Positional Density (motifs/kbp)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Motif Subclass', fontsize=11, fontweight='bold')
    ax2.set_title('B. Positional Density (λ) by Subclass', fontsize=12, fontweight='bold', pad=10)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(display_labels, fontsize=8)
    
    # Add value labels
    max_val2 = max(positional_vals) if positional_vals else 1
    for i, (bar, val) in enumerate(zip(bars2, positional_vals)):
        if val > 0:
            ax2.text(val + max_val2 * 0.01, i, f'{val:.2f}', 
                    va='center', fontsize=7, fontweight='bold')
    
    plt.suptitle(title, fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    
    return fig


def plot_enrichment_analysis_by_subclass(enrichment_results: Dict[str, Dict[str, Any]],
                                         title: str = "Motif Enrichment Analysis (by Subclass)",
                                         figsize: Tuple[int, int] = (16, 12)) -> plt.Figure:
    """
    Plot enrichment analysis results at subclass level with fold enrichment and p-values.
    
    Args:
        enrichment_results: Dictionary from calculate_enrichment_with_shuffling with by_subclass=True
        title: Plot title
        figsize: Figure size
        
    Returns:
        Matplotlib figure object
    """
    set_scientific_style()
    
    # Extract data (exclude 'Overall' for subclass-specific view)
    subclasses = [k for k in enrichment_results.keys() if k != 'Overall' and ':' in k]
    if not subclasses:
        # Fallback to class level
        subclasses = [k for k in enrichment_results.keys() if k != 'Overall']
    
    # Sort alphabetically
    subclasses.sort()
    
    fold_enrichments = []
    p_values = []
    observed_densities = []
    background_means = []
    
    for subclass in subclasses:
        result = enrichment_results[subclass]
        fe = result.get('fold_enrichment', 0)
        # Handle both string 'Inf' and float infinity values robustly
        if fe == 'Inf' or (isinstance(fe, float) and np.isinf(fe)):
            fe = INFINITE_FOLD_ENRICHMENT_CAP  # Cap infinite values for visualization
        fold_enrichments.append(fe)
        p_values.append(result.get('p_value', 1.0))
        observed_densities.append(result.get('observed_density', 0))
        background_means.append(result.get('background_mean', 0))
    
    # Create figure with subplots
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 1, hspace=0.35, height_ratios=[2, 2, 1.5])
    
    y_pos = np.arange(len(subclasses))
    
    # Get colors based on parent class
    colors = []
    for subclass in subclasses:
        if ':' in subclass:
            parent_class = subclass.split(':')[0]
            colors.append(MOTIF_CLASS_COLORS.get(parent_class, '#808080'))
        else:
            colors.append(MOTIF_CLASS_COLORS.get(subclass, '#808080'))
    
    # 1. Fold Enrichment
    ax1 = fig.add_subplot(gs[0])
    bars1 = ax1.barh(y_pos, fold_enrichments, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    ax1.axvline(x=1.0, color='red', linestyle='--', linewidth=2, label='No enrichment (FE=1)')
    ax1.set_xlabel('Fold Enrichment', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Motif Subclass', fontsize=11, fontweight='bold')
    ax1.set_title('A. Fold Enrichment by Subclass', fontsize=12, fontweight='bold')
    ax1.set_yticks(y_pos)
    
    # Format labels
    display_labels = [label.replace('_', ' ').replace(':', ': ') for label in subclasses]
    ax1.set_yticklabels(display_labels, fontsize=8)
    ax1.legend(loc='best', fontsize=9)
    
    # Add value labels
    max_fe = max(fold_enrichments) if fold_enrichments else 1
    for i, (bar, val) in enumerate(zip(bars1, fold_enrichments)):
        label_text = f'{val:.2f}' if val < INFINITE_FOLD_ENRICHMENT_CAP else 'Inf'
        ax1.text(val + max_fe * 0.01, i, label_text, 
                va='center', fontsize=7, fontweight='bold')
    
    # 2. P-values
    ax2 = fig.add_subplot(gs[1])
    # Color code by significance
    p_colors = ['green' if p < 0.05 else 'orange' if p < 0.1 else 'red' for p in p_values]
    bars2 = ax2.barh(y_pos, p_values, color=p_colors, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax2.axvline(x=0.05, color='green', linestyle='--', linewidth=2, label='p=0.05')
    ax2.axvline(x=0.1, color='orange', linestyle='--', linewidth=1.5, label='p=0.1')
    ax2.set_xlabel('P-value', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Motif Subclass', fontsize=11, fontweight='bold')
    ax2.set_title('B. Statistical Significance by Subclass', fontsize=12, fontweight='bold')
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(display_labels, fontsize=8)
    ax2.set_xlim(0, max(1.0, max(p_values) * 1.1))
    ax2.legend(loc='best', fontsize=9)
    
    # Add value labels
    for i, (bar, val) in enumerate(zip(bars2, p_values)):
        ax2.text(val + 0.02, i, f'{val:.3f}', va='center', fontsize=7, fontweight='bold')
    
    # 3. Observed vs Background Density (grouped bar chart)
    ax3 = fig.add_subplot(gs[2])
    x = np.arange(len(subclasses))
    width = 0.35
    
    bars3a = ax3.bar(x - width/2, observed_densities, width, label='Observed', 
                    color='steelblue', alpha=0.8, edgecolor='black', linewidth=0.5)
    bars3b = ax3.bar(x + width/2, background_means, width, label='Background (Mean)', 
                    color='coral', alpha=0.8, edgecolor='black', linewidth=0.5)
    
    ax3.set_xlabel('Motif Subclass', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Density (%)', fontsize=11, fontweight='bold')
    ax3.set_title('C. Observed vs. Background Density Comparison', fontsize=12, fontweight='bold')
    ax3.set_xticks(x)
    ax3.set_xticklabels(display_labels, rotation=45, ha='right', fontsize=7)
    ax3.legend(loc='best', fontsize=10)
    ax3.grid(axis='y', alpha=0.3)
    
    plt.suptitle(title, fontsize=14, fontweight='bold', y=0.995)
    
    return fig


def plot_subclass_density_heatmap(motifs: List[Dict[str, Any]], 
                                  sequence_length: int,
                                  window_size: int = 1000,
                                  title: str = "Subclass Density Heatmap",
                                  figsize: Tuple[int, int] = (16, 12)) -> plt.Figure:
    """
    Create a heatmap showing density of each subclass across the sequence.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Length of the analyzed sequence
        window_size: Window size for density calculation
        title: Plot title
        figsize: Figure size
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style()
    
    if not motifs or sequence_length == 0:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    # Calculate windows
    num_windows = max(1, sequence_length // window_size)
    windows = np.linspace(0, sequence_length, num_windows + 1)
    
    # Get unique subclasses
    subclass_groups = defaultdict(list)
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        subclass_name = motif.get('Subclass', 'Unknown')
        key = f"{class_name}:{subclass_name}"
        subclass_groups[key].append(motif)
    
    subclasses = sorted(subclass_groups.keys())
    
    # Calculate density matrix
    density_matrix = np.zeros((len(subclasses), num_windows))
    
    for i, subclass_key in enumerate(subclasses):
        subclass_motifs = subclass_groups[subclass_key]
        
        for j in range(num_windows):
            window_start = windows[j]
            window_end = windows[j + 1]
            
            # Count motifs in window
            count = 0
            for motif in subclass_motifs:
                motif_start = motif.get('Start', 0) - 1  # 0-based
                motif_end = motif.get('End', 0)
                
                # Check if motif overlaps with window
                if not (motif_end <= window_start or motif_start >= window_end):
                    count += 1
            
            density_matrix[i, j] = count
    
    # Create heatmap with publication quality
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    # Use colorblind-friendly colormap
    im = ax.imshow(density_matrix, cmap='viridis', aspect='auto', interpolation='nearest')
    
    # Customize axes (Nature style)
    ax.set_xlabel(f'Position (kb)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Motif Subclass', fontsize=12, fontweight='bold')
    
    # Replace underscores with spaces in title
    display_title = title.replace('_', ' ')
    ax.set_title(display_title, fontsize=14, fontweight='bold', pad=10)
    
    # Set ticks and labels with underscores replaced by spaces
    ax.set_yticks(range(len(subclasses)))
    display_subclasses = [sc.replace('_', ' ').replace(':', ': ') for sc in subclasses]
    ax.set_yticklabels(display_subclasses, fontsize=8)
    
    # Clean x-axis with kb units
    x_ticks = np.arange(0, num_windows, max(1, num_windows // 10))
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f'{int(windows[i]/1000)}' for i in x_ticks])
    
    # Add colorbar with proper label
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Motif Count', fontsize=10, fontweight='bold')
    cbar.ax.tick_params(labelsize=8)
    
    plt.tight_layout()
    return fig


# =============================================================================
# NATURE-LEVEL PUBLICATION VISUALIZATIONS (GENOME-WIDE)
# =============================================================================

def plot_manhattan_motif_density(motifs: List[Dict[str, Any]], 
                                  sequence_length: int,
                                  window_size: int = None,
                                  score_type: str = 'density',
                                  title: str = "Manhattan Plot - Motif Distribution",
                                  figsize: Tuple[int, int] = None) -> plt.Figure:
    """
    Create Manhattan plot showing motif density or score across genomic coordinates.
    
    Ideal for highlighting hotspots, clusters, and hybrid zones in large genomes.
    Publication-quality visualization following Nature Methods guidelines.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total length of analyzed sequence in bp
        window_size: Window size for density calculation (auto if None)
        score_type: 'density' for motif count or 'score' for average score
        title: Plot title
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object (publication-ready at 300 DPI)
    """
    set_scientific_style()
    
    if figsize is None:
        figsize = FIGURE_SIZES['wide']
    
    if not motifs or sequence_length == 0:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    # Auto-calculate window size (1% of sequence or minimum 1kb)
    if window_size is None:
        window_size = max(1000, sequence_length // 100)
    
    num_windows = max(1, sequence_length // window_size)
    
    # Get unique classes (exclude synthetic classes for cleaner visualization)
    classes = sorted(set(m.get('Class', 'Unknown') for m in motifs
                        if m.get('Class') not in CIRCOS_EXCLUDED_CLASSES))
    
    if not classes:
        classes = sorted(set(m.get('Class', 'Unknown') for m in motifs))
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    # Calculate metric for each window and class
    positions_kb = []
    values = []
    colors = []
    
    for class_name in classes:
        class_motifs = [m for m in motifs if m.get('Class') == class_name]
        color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        
        for i in range(num_windows):
            window_start = i * window_size
            window_end = (i + 1) * window_size
            window_center_kb = (window_start + window_end) / 2 / 1000
            
            # Count motifs in this window
            window_motifs = []
            for motif in class_motifs:
                motif_start = motif.get('Start', 0) - 1
                motif_end = motif.get('End', 0)
                if not (motif_end <= window_start or motif_start >= window_end):
                    window_motifs.append(motif)
            
            if window_motifs:
                if score_type == 'density':
                    # Density: motifs per kb
                    value = len(window_motifs) / (window_size / 1000)
                else:  # score
                    # Average score in window
                    scores = [m.get('Score', 0) for m in window_motifs if isinstance(m.get('Score'), (int, float))]
                    value = np.mean(scores) if scores else 0
                
                positions_kb.append(window_center_kb)
                values.append(value)
                colors.append(color)
    
    # Plot points with class-specific colors
    ax.scatter(positions_kb, values, c=colors, s=20, alpha=0.6, edgecolors='black', linewidth=0.3)
    
    # Styling
    ax.set_xlabel('Genomic Position (kb)', fontsize=12, fontweight='bold')
    y_label = 'Motif Density (motifs/kb)' if score_type == 'density' else 'Average Motif Score'
    ax.set_ylabel(y_label, fontsize=12, fontweight='bold')
    
    # Replace underscores with spaces in title
    display_title = title.replace('_', ' ')
    ax.set_title(display_title, fontsize=14, fontweight='bold', pad=10)
    
    # Add horizontal grid for readability
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    
    # Create legend with class names
    legend_elements = []
    display_classes = []
    for class_name in classes:
        if any(c == MOTIF_CLASS_COLORS.get(class_name, '#808080') for c in colors):
            legend_elements.append(plt.scatter([], [], c=MOTIF_CLASS_COLORS.get(class_name, '#808080'), 
                                             s=50, alpha=0.6, edgecolors='black', linewidth=0.5))
            display_classes.append(class_name.replace('_', ' '))
    
    if legend_elements:
        ax.legend(legend_elements, display_classes, loc='upper right', 
                 fontsize=8, framealpha=0.9, ncol=min(3, len(legend_elements)))
    
    # Apply Nature journal style
    _apply_nature_style(ax)
    
    plt.tight_layout()
    return fig


def plot_cumulative_motif_distribution(motifs: List[Dict[str, Any]], 
                                       sequence_length: int,
                                       title: str = "Cumulative Motif Distribution",
                                       figsize: Tuple[int, int] = None,
                                       by_class: bool = True) -> plt.Figure:
    """
    Create cumulative distribution plot showing running sum of motifs over genome.
    
    Useful for comparing motifs or samples, showing how motif accumulation
    varies across the sequence. Publication-quality following Nature guidelines.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total length of analyzed sequence in bp
        title: Plot title
        figsize: Figure size (width, height)
        by_class: Whether to separate by motif class
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style()
    
    if figsize is None:
        figsize = FIGURE_SIZES['wide']
    
    if not motifs or sequence_length == 0:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    if by_class:
        # Group by class
        classes = sorted(set(m.get('Class', 'Unknown') for m in motifs
                            if m.get('Class') not in CIRCOS_EXCLUDED_CLASSES))
        
        if not classes:
            classes = sorted(set(m.get('Class', 'Unknown') for m in motifs))
        
        for class_name in classes:
            class_motifs = [m for m in motifs if m.get('Class') == class_name]
            
            # Sort by start position
            class_motifs_sorted = sorted(class_motifs, key=lambda m: m.get('Start', 0))
            
            # Calculate cumulative count
            positions = [0]
            cumulative = [0]
            
            for i, motif in enumerate(class_motifs_sorted, 1):
                positions.append(motif.get('Start', 0) / 1000)  # Convert to kb
                cumulative.append(i)
            
            # Add final point at sequence end
            positions.append(sequence_length / 1000)
            cumulative.append(len(class_motifs))
            
            # Plot with class color
            color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
            display_name = class_name.replace('_', ' ')
            ax.plot(positions, cumulative, color=color, linewidth=1.5, 
                   label=display_name, alpha=0.8)
    else:
        # Overall cumulative
        motifs_sorted = sorted(motifs, key=lambda m: m.get('Start', 0))
        
        positions = [0]
        cumulative = [0]
        
        for i, motif in enumerate(motifs_sorted, 1):
            positions.append(motif.get('Start', 0) / 1000)
            cumulative.append(i)
        
        positions.append(sequence_length / 1000)
        cumulative.append(len(motifs))
        
        ax.plot(positions, cumulative, color='#0072B2', linewidth=2, alpha=0.8)
    
    # Styling
    ax.set_xlabel('Genomic Position (kb)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Cumulative Motif Count', fontsize=12, fontweight='bold')
    
    # Replace underscores in title
    display_title = title.replace('_', ' ')
    ax.set_title(display_title, fontsize=14, fontweight='bold', pad=10)
    
    # Add grid
    ax.grid(alpha=0.3, linestyle='--', linewidth=0.5)
    
    if by_class:
        ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
    
    # Apply Nature journal style
    _apply_nature_style(ax)
    
    plt.tight_layout()
    return fig


def plot_motif_cooccurrence_matrix(motifs: List[Dict[str, Any]], 
                                   title: str = "Motif Co-occurrence Matrix",
                                   figsize: Tuple[int, int] = None,
                                   overlap_threshold: int = 1) -> plt.Figure:
    """
    Create heatmap showing co-occurrence frequency between motif classes.
    
    Shows which motif classes tend to appear together (within overlap_threshold bp).
    Excellent for publication figures showing motif relationships.
    
    Args:
        motifs: List of motif dictionaries
        title: Plot title
        figsize: Figure size (width, height)
        overlap_threshold: Maximum distance (bp) to consider as co-occurrence
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style()
    
    if figsize is None:
        figsize = FIGURE_SIZES['square']
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    # Get unique classes (exclude synthetic ones)
    classes = sorted(set(m.get('Class', 'Unknown') for m in motifs
                        if m.get('Class') not in CIRCOS_EXCLUDED_CLASSES))
    
    if not classes:
        classes = sorted(set(m.get('Class', 'Unknown') for m in motifs))
    
    # Initialize co-occurrence matrix
    n_classes = len(classes)
    cooccurrence_matrix = np.zeros((n_classes, n_classes))
    
    # Calculate co-occurrences
    for i, class_i in enumerate(classes):
        motifs_i = [m for m in motifs if m.get('Class') == class_i]
        
        for j, class_j in enumerate(classes):
            motifs_j = [m for m in motifs if m.get('Class') == class_j]
            
            # Count overlaps
            count = 0
            for mi in motifs_i:
                start_i = mi.get('Start', 0)
                end_i = mi.get('End', 0)
                
                for mj in motifs_j:
                    start_j = mj.get('Start', 0)
                    end_j = mj.get('End', 0)
                    
                    # Check if they overlap or are within threshold
                    distance = max(0, max(start_i, start_j) - min(end_i, end_j))
                    
                    if distance <= overlap_threshold:
                        count += 1
            
            cooccurrence_matrix[i, j] = count
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    # Use colorblind-friendly colormap
    im = ax.imshow(cooccurrence_matrix, cmap='YlOrRd', aspect='auto', interpolation='nearest')
    
    # Set ticks and labels
    ax.set_xticks(range(n_classes))
    ax.set_yticks(range(n_classes))
    
    # Replace underscores with spaces in labels
    display_classes = [c.replace('_', ' ') for c in classes]
    ax.set_xticklabels(display_classes, rotation=45, ha='right', fontsize=9)
    ax.set_yticklabels(display_classes, fontsize=9)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Co-occurrence Count', fontsize=10, fontweight='bold')
    cbar.ax.tick_params(labelsize=8)
    
    # Add values to cells (if not too many)
    if n_classes <= 10:
        for i in range(n_classes):
            for j in range(n_classes):
                value = int(cooccurrence_matrix[i, j])
                if value > 0:
                    text_color = 'white' if value > cooccurrence_matrix.max() / 2 else 'black'
                    ax.text(j, i, str(value), ha='center', va='center', 
                           color=text_color, fontsize=8, fontweight='bold')
    
    # Title and labels
    display_title = title.replace('_', ' ')
    ax.set_title(display_title, fontsize=14, fontweight='bold', pad=10)
    ax.set_xlabel('Motif Class', fontsize=11, fontweight='bold')
    ax.set_ylabel('Motif Class', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    return fig


def plot_gc_content_correlation(motifs: List[Dict[str, Any]], 
                                sequence: str,
                                window_size: int = 1000,
                                title: str = "Motif Density vs GC Content",
                                figsize: Tuple[int, int] = None) -> plt.Figure:
    """
    Create scatter plot showing correlation between GC content and motif density.
    
    Shows GC-driven motif enrichment patterns. One dot per genomic window.
    Publication-quality visualization with regression line.
    
    Args:
        motifs: List of motif dictionaries
        sequence: DNA sequence string
        window_size: Window size for GC and density calculation
        title: Plot title
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style()
    
    if figsize is None:
        figsize = FIGURE_SIZES['one_and_half']
    
    sequence_length = len(sequence)
    
    if not motifs or sequence_length == 0:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    num_windows = max(1, sequence_length // window_size)
    
    # Calculate GC% and motif density for each window
    gc_percentages = []
    motif_densities = []
    
    for i in range(num_windows):
        window_start = i * window_size
        window_end = min((i + 1) * window_size, sequence_length)
        
        # Calculate GC%
        window_seq = sequence[window_start:window_end].upper()
        gc_count = window_seq.count('G') + window_seq.count('C')
        window_length = len(window_seq)
        gc_pct = (gc_count / window_length * 100) if window_length > 0 else 0
        
        # Count motifs in window
        motif_count = 0
        for motif in motifs:
            motif_start = motif.get('Start', 0) - 1
            motif_end = motif.get('End', 0)
            if not (motif_end <= window_start or motif_start >= window_end):
                motif_count += 1
        
        # Density: motifs per kb
        density = motif_count / (window_size / 1000)
        
        gc_percentages.append(gc_pct)
        motif_densities.append(density)
    
    # Create scatter plot
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    ax.scatter(gc_percentages, motif_densities, alpha=0.5, s=30, 
              color='#0072B2', edgecolors='black', linewidth=0.3)
    
    # Add regression line
    if len(gc_percentages) > 1:
        z = np.polyfit(gc_percentages, motif_densities, 1)
        p = np.poly1d(z)
        x_line = np.linspace(min(gc_percentages), max(gc_percentages), 100)
        ax.plot(x_line, p(x_line), 'r--', linewidth=1.5, alpha=0.8, 
               label=f'Linear fit: y={z[0]:.2f}x+{z[1]:.2f}')
        
        # Calculate correlation coefficient
        correlation = np.corrcoef(gc_percentages, motif_densities)[0, 1]
        ax.text(0.05, 0.95, f'R = {correlation:.3f}', transform=ax.transAxes,
               fontsize=10, fontweight='bold', verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Styling
    ax.set_xlabel('GC Content (%)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Motif Density (motifs/kb)', fontsize=12, fontweight='bold')
    
    display_title = title.replace('_', ' ')
    ax.set_title(display_title, fontsize=14, fontweight='bold', pad=10)
    
    ax.grid(alpha=0.3, linestyle='--', linewidth=0.5)
    
    if len(gc_percentages) > 1:
        ax.legend(loc='upper right', fontsize=9, framealpha=0.9)
    
    # Apply Nature journal style
    _apply_nature_style(ax)
    
    plt.tight_layout()
    return fig


def plot_linear_motif_track(motifs: List[Dict[str, Any]], 
                            sequence_length: int,
                            region_start: int = 0,
                            region_end: int = None,
                            title: str = "Linear Motif Track",
                            figsize: Tuple[int, int] = None,
                            show_labels: bool = True) -> plt.Figure:
    """
    Create horizontal graphical track with colored blocks for motifs.
    
    Best for visualizing <10kb regions. Colored blocks show motif positions
    with class-specific colors. Publication-quality linear genome browser view.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total length of analyzed sequence in bp
        region_start: Start position of region to display (0-based)
        region_end: End position of region to display (None = full sequence)
        title: Plot title
        figsize: Figure size (width, height)
        show_labels: Whether to show motif class labels
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style()
    
    if figsize is None:
        figsize = FIGURE_SIZES['wide']
    
    if region_end is None:
        region_end = sequence_length
    
    # Filter motifs in region
    region_motifs = [m for m in motifs 
                     if m.get('End', 0) > region_start and m.get('Start', 0) < region_end]
    
    if not region_motifs:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, f'No motifs in region {region_start}-{region_end}', 
               ha='center', va='center', transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    # Group by class
    class_motifs = defaultdict(list)
    for motif in region_motifs:
        class_name = motif.get('Class', 'Unknown')
        class_motifs[class_name].append(motif)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    # Plot each class on a separate track
    classes = sorted(class_motifs.keys())
    track_height = 0.6
    track_spacing = 1.0
    
    for i, class_name in enumerate(classes):
        y_pos = i * track_spacing
        color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        
        for motif in class_motifs[class_name]:
            start = max(region_start, motif.get('Start', 0))
            end = min(region_end, motif.get('End', 0))
            length = end - start
            
            # Draw motif as rectangle
            rect = patches.Rectangle(
                (start, y_pos - track_height/2), length, track_height,
                facecolor=color, edgecolor='black', linewidth=0.5, alpha=0.8
            )
            ax.add_patch(rect)
            
            # Add score if available (small label inside block)
            if length > (region_end - region_start) * 0.02:  # Only if block is large enough
                score = motif.get('Score')
                if isinstance(score, (int, float)):
                    ax.text(start + length/2, y_pos, f'{score:.2f}', 
                           ha='center', va='center', fontsize=6, color='white', fontweight='bold')
        
        # Add class label on the left
        if show_labels:
            display_name = class_name.replace('_', ' ')
            ax.text(region_start - (region_end - region_start) * 0.02, y_pos, 
                   display_name, ha='right', va='center', fontsize=9, fontweight='bold')
    
    # Styling
    ax.set_xlim(region_start, region_end)
    ax.set_ylim(-0.5, len(classes) * track_spacing - 0.5)
    
    ax.set_xlabel('Position (bp)', fontsize=12, fontweight='bold')
    ax.set_yticks([])
    
    display_title = title.replace('_', ' ')
    ax.set_title(display_title, fontsize=14, fontweight='bold', pad=10)
    
    # Add position ruler at top
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    
    # Format x-axis
    ax.ticklabel_format(style='plain', axis='x')
    
    # Apply Nature journal style (minimal spines)
    for spine in ['left', 'right', 'bottom']:
        ax.spines[spine].set_visible(False)
    ax.spines['top'].set_visible(True)
    
    plt.tight_layout()
    return fig


def plot_cluster_size_distribution(motifs: List[Dict[str, Any]], 
                                   title: str = "Cluster Size Distribution",
                                   figsize: Tuple[int, int] = None) -> plt.Figure:
    """
    Plot distribution of cluster sizes (number of motifs per cluster).
    
    Shows histogram and statistics of cluster composition.
    Publication-quality visualization.
    
    Args:
        motifs: List of motif dictionaries (should include cluster motifs)
        title: Plot title
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style()
    
    if figsize is None:
        figsize = FIGURE_SIZES['one_and_half']
    
    # Extract cluster motifs
    cluster_motifs = [m for m in motifs if m.get('Class') == 'Non-B_DNA_Clusters']
    
    if not cluster_motifs:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No cluster motifs found', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    # Extract cluster sizes
    cluster_sizes = []
    cluster_diversities = []
    
    for motif in cluster_motifs:
        size = motif.get('Motif_Count', 0)
        diversity = motif.get('Class_Diversity', 0)
        if size > 0:
            cluster_sizes.append(size)
            cluster_diversities.append(diversity)
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, dpi=PUBLICATION_DPI)
    
    # 1. Cluster size histogram
    ax1.hist(cluster_sizes, bins=min(20, max(cluster_sizes) if cluster_sizes else 1), 
            edgecolor='black', linewidth=0.5, color='#4ecdc4', alpha=0.7)
    
    ax1.set_xlabel('Motifs per Cluster', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Frequency', fontsize=11, fontweight='bold')
    ax1.set_title('Cluster Size Distribution', fontsize=12, fontweight='bold')
    
    # Add statistics
    if cluster_sizes:
        mean_size = np.mean(cluster_sizes)
        median_size = np.median(cluster_sizes)
        ax1.axvline(mean_size, color='red', linestyle='--', linewidth=1.5, 
                   label=f'Mean: {mean_size:.1f}')
        ax1.axvline(median_size, color='orange', linestyle='--', linewidth=1.5, 
                   label=f'Median: {median_size:.1f}')
        ax1.legend(fontsize=9, framealpha=0.9)
    
    ax1.grid(axis='y', alpha=0.3)
    _apply_nature_style(ax1)
    
    # 2. Class diversity histogram
    if cluster_diversities:
        ax2.hist(cluster_diversities, bins=min(10, max(cluster_diversities)), 
                edgecolor='black', linewidth=0.5, color='#95e1d3', alpha=0.7)
        
        ax2.set_xlabel('Class Diversity per Cluster', fontsize=11, fontweight='bold')
        ax2.set_ylabel('Frequency', fontsize=11, fontweight='bold')
        ax2.set_title('Cluster Diversity Distribution', fontsize=12, fontweight='bold')
        
        mean_div = np.mean(cluster_diversities)
        ax2.axvline(mean_div, color='red', linestyle='--', linewidth=1.5, 
                   label=f'Mean: {mean_div:.1f}')
        ax2.legend(fontsize=9, framealpha=0.9)
        
        ax2.grid(axis='y', alpha=0.3)
        _apply_nature_style(ax2)
    else:
        ax2.text(0.5, 0.5, 'No diversity data', ha='center', va='center',
                transform=ax2.transAxes, fontsize=12)
        ax2.axis('off')
    
    display_title = title.replace('_', ' ')
    plt.suptitle(display_title, fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    return fig


def plot_motif_length_kde(motifs: List[Dict[str, Any]], 
                          by_class: bool = True,
                          title: str = "Motif Length Distribution (KDE)",
                          figsize: Tuple[int, int] = None) -> plt.Figure:
    """
    Plot kernel density estimation of motif length distributions.
    
    Shows smooth probability density curves for motif lengths,
    useful for comparing length patterns across classes.
    Publication-quality visualization.
    
    Args:
        motifs: List of motif dictionaries
        by_class: Whether to separate by motif class
        title: Plot title
        figsize: Figure size (width, height)
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style()
    
    if figsize is None:
        figsize = FIGURE_SIZES['one_and_half']
    
    if not motifs:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
    
    if by_class:
        # Group by class
        classes = sorted(set(m.get('Class', 'Unknown') for m in motifs
                            if m.get('Class') not in CIRCOS_EXCLUDED_CLASSES))
        
        if not classes:
            classes = sorted(set(m.get('Class', 'Unknown') for m in motifs))
        
        for class_name in classes:
            class_motifs = [m for m in motifs if m.get('Class') == class_name]
            lengths = [m.get('Length', 0) for m in class_motifs if m.get('Length', 0) > 0]
            
            if len(lengths) > 1:
                color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
                display_name = class_name.replace('_', ' ')
                
                # Plot KDE
                try:
                    from scipy import stats
                    kde = stats.gaussian_kde(lengths)
                    x_range = np.linspace(min(lengths), max(lengths), 200)
                    density = kde(x_range)
                    ax.plot(x_range, density, color=color, linewidth=2, 
                           label=display_name, alpha=0.8)
                    ax.fill_between(x_range, density, alpha=0.2, color=color)
                except:
                    # Fallback to histogram if KDE fails
                    ax.hist(lengths, bins=20, alpha=0.3, color=color, 
                           label=display_name, density=True, edgecolor='black', linewidth=0.5)
    else:
        # Overall distribution
        lengths = [m.get('Length', 0) for m in motifs if m.get('Length', 0) > 0]
        
        if len(lengths) > 1:
            try:
                from scipy import stats
                kde = stats.gaussian_kde(lengths)
                x_range = np.linspace(min(lengths), max(lengths), 200)
                density = kde(x_range)
                ax.plot(x_range, density, color='#0072B2', linewidth=2.5, alpha=0.8)
                ax.fill_between(x_range, density, alpha=0.3, color='#0072B2')
            except:
                ax.hist(lengths, bins=30, alpha=0.7, color='#0072B2', 
                       density=True, edgecolor='black', linewidth=0.5)
    
    # Styling
    ax.set_xlabel('Motif Length (bp)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Probability Density', fontsize=12, fontweight='bold')
    
    display_title = title.replace('_', ' ')
    ax.set_title(display_title, fontsize=14, fontweight='bold', pad=10)
    
    ax.grid(alpha=0.3, linestyle='--', linewidth=0.5)
    
    if by_class:
        ax.legend(loc='upper right', fontsize=8, framealpha=0.9, ncol=2)
    
    # Apply Nature journal style
    _apply_nature_style(ax)
    
    plt.tight_layout()
    return fig


# =============================================================================
# FUNCTION ALIASES FOR BACKWARD COMPATIBILITY
# =============================================================================

# Aliases for functions with different naming conventions
plot_comprehensive_class_analysis = plot_class_analysis_comprehensive
plot_comprehensive_subclass_analysis = plot_subclass_analysis_comprehensive


# =============================================================================
# ADDITIONAL VISUALIZATION FUNCTIONS FOR NOTEBOOK COMPATIBILITY
# =============================================================================

def plot_genome_landscape_track(motifs: List[Dict[str, Any]], 
                               sequence_length: int,
                               title: str = "Genome Landscape Track",
                               figsize: Tuple[int, int] = None,
                               window_size: int = None) -> plt.Figure:
    """
    Create genome landscape track visualization showing motif distribution along sequence.
    
    This is a simplified horizontal track view showing motif positions and density.
    Similar to plot_linear_motif_track but with additional density information.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total length of analyzed sequence in bp
        title: Plot title
        figsize: Figure size (width, height)
        window_size: Window size for density calculation (auto if None)
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style()
    
    if figsize is None:
        figsize = FIGURE_SIZES['wide']
    
    if not motifs or sequence_length == 0:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    # Auto-calculate window size
    if window_size is None:
        window_size = max(1000, sequence_length // 50)
    
    # Create figure with two subplots: density track + motif positions
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, dpi=PUBLICATION_DPI,
                                    gridspec_kw={'height_ratios': [1, 2]})
    
    # Top panel: Density line plot
    num_windows = max(1, sequence_length // window_size)
    positions = []
    densities = []
    
    for i in range(num_windows):
        window_start = i * window_size
        window_end = (i + 1) * window_size
        window_center = (window_start + window_end) / 2 / 1000  # In kb
        
        # Count motifs in window
        count = 0
        for motif in motifs:
            motif_start = motif.get('Start', 0) - 1
            motif_end = motif.get('End', 0)
            if not (motif_end <= window_start or motif_start >= window_end):
                count += 1
        
        density = count / (window_size / 1000)  # motifs per kb
        positions.append(window_center)
        densities.append(density)
    
    # Plot density line
    ax1.fill_between(positions, densities, alpha=0.3, color='#0072B2')
    ax1.plot(positions, densities, color='#0072B2', linewidth=1.5)
    ax1.set_ylabel('Density\n(motifs/kb)', fontsize=10, fontweight='bold')
    ax1.set_xlim(0, sequence_length / 1000)
    ax1.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    _apply_nature_style(ax1)
    
    # Bottom panel: Motif track by class
    class_motifs = defaultdict(list)
    for motif in motifs:
        class_name = motif.get('Class', 'Unknown')
        if class_name not in CIRCOS_EXCLUDED_CLASSES:
            class_motifs[class_name].append(motif)
    
    if not class_motifs:
        # Include all classes if none pass the filter
        for motif in motifs:
            class_name = motif.get('Class', 'Unknown')
            class_motifs[class_name].append(motif)
    
    classes = sorted(class_motifs.keys())
    track_height = 0.5
    track_spacing = 1.0
    
    for i, class_name in enumerate(classes):
        y_pos = i * track_spacing
        color = MOTIF_CLASS_COLORS.get(class_name, '#808080')
        
        for motif in class_motifs[class_name]:
            start_kb = motif.get('Start', 0) / 1000
            end_kb = motif.get('End', 0) / 1000
            length_kb = end_kb - start_kb
            
            # Draw motif as rectangle
            rect = patches.Rectangle(
                (start_kb, y_pos - track_height/2), length_kb, track_height,
                facecolor=color, edgecolor='black', linewidth=0.3, alpha=0.8
            )
            ax2.add_patch(rect)
    
    # Customize bottom panel
    ax2.set_xlim(0, sequence_length / 1000)
    ax2.set_ylim(-0.5, len(classes) * track_spacing - 0.5)
    ax2.set_xlabel('Position (kb)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Motif Class', fontsize=11, fontweight='bold')
    ax2.set_yticks([i * track_spacing for i in range(len(classes))])
    display_classes = [c.replace('_', ' ') for c in classes]
    ax2.set_yticklabels(display_classes, fontsize=9)
    _apply_nature_style(ax2)
    
    # Overall title
    display_title = title.replace('_', ' ')
    fig.suptitle(display_title, fontsize=14, fontweight='bold', y=0.98)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


def plot_sliding_window_heat_ribbon(motifs: List[Dict[str, Any]], 
                                    sequence_length: int,
                                    title: str = "Sliding Window Heat Ribbon",
                                    figsize: Tuple[int, int] = None,
                                    window_size: int = None) -> plt.Figure:
    """
    Create a 1D heatmap ribbon showing motif density along the sequence.
    
    This creates a horizontal ribbon colored by density values, with an
    accompanying line plot showing the density profile.
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total length of analyzed sequence in bp
        title: Plot title
        figsize: Figure size (width, height)
        window_size: Window size for density calculation (auto if None)
        
    Returns:
        Matplotlib figure object (publication-ready)
    """
    set_scientific_style()
    
    if figsize is None:
        figsize = FIGURE_SIZES['wide']
    
    if not motifs or sequence_length == 0:
        fig, ax = plt.subplots(figsize=figsize, dpi=PUBLICATION_DPI)
        ax.text(0.5, 0.5, 'No motifs to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        ax.set_title(title)
        return fig
    
    # Auto-calculate window size
    if window_size is None:
        window_size = max(1000, sequence_length // 100)
    
    num_windows = max(1, sequence_length // window_size)
    
    # Calculate density for each window
    density_values = []
    positions = []
    
    for i in range(num_windows):
        window_start = i * window_size
        window_end = (i + 1) * window_size
        
        # Count motifs in window
        count = 0
        for motif in motifs:
            motif_start = motif.get('Start', 0) - 1
            motif_end = motif.get('End', 0)
            if not (motif_end <= window_start or motif_start >= window_end):
                count += 1
        
        density = count / (window_size / 1000)  # motifs per kb
        density_values.append(density)
        positions.append((window_start + window_end) / 2 / 1000)  # Center position in kb
    
    # Create figure with two subplots
    fig = plt.figure(figsize=figsize, dpi=PUBLICATION_DPI)
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 3], hspace=0.3)
    
    # Top: Heat ribbon (1D heatmap)
    ax1 = fig.add_subplot(gs[0])
    
    # Create 2D array for heatmap (1 row)
    density_array = np.array(density_values).reshape(1, -1)
    
    # Plot heatmap
    im = ax1.imshow(density_array, cmap='YlOrRd', aspect='auto', 
                    extent=[0, sequence_length / 1000, 0, 1],
                    interpolation='bilinear')
    
    ax1.set_yticks([])
    ax1.set_ylabel('Density\nHeatmap', fontsize=10, fontweight='bold')
    ax1.set_xticks([])
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax1, orientation='vertical', pad=0.02, shrink=0.8)
    cbar.set_label('Motifs/kb', fontsize=9, fontweight='bold')
    cbar.ax.tick_params(labelsize=8)
    
    # Bottom: Line plot
    ax2 = fig.add_subplot(gs[1])
    
    ax2.plot(positions, density_values, color='#0072B2', linewidth=2, alpha=0.8)
    ax2.fill_between(positions, density_values, alpha=0.3, color='#0072B2')
    
    # Mark peaks (top 10% density)
    if density_values:
        threshold = np.percentile(density_values, 90)
        peak_indices = [i for i, d in enumerate(density_values) if d >= threshold]
        if peak_indices:
            peak_positions = [positions[i] for i in peak_indices]
            peak_densities = [density_values[i] for i in peak_indices]
            ax2.scatter(peak_positions, peak_densities, color='red', s=50, 
                       zorder=5, alpha=0.7, edgecolors='black', linewidth=0.5,
                       label=f'High density (top 10%)')
    
    # Styling
    ax2.set_xlabel('Position (kb)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Motif Density (motifs/kb)', fontsize=11, fontweight='bold')
    ax2.set_xlim(0, sequence_length / 1000)
    ax2.grid(alpha=0.3, linestyle='--', linewidth=0.5)
    
    if peak_indices:
        ax2.legend(loc='upper right', fontsize=9, framealpha=0.9)
    
    _apply_nature_style(ax2)
    
    # Overall title
    display_title = title.replace('_', ' ')
    fig.suptitle(display_title, fontsize=14, fontweight='bold', y=0.98)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


if __name__ == "__main__":
    test_visualizations()