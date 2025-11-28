"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                    CONSOLIDATED UTILITIES MODULE                              ║
║        Utility Functions for Sequence Processing and Data Export             ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: utilities.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.1 - Consolidated
LICENSE: MIT

DESCRIPTION:
    Consolidated module containing all utility functions for sequence processing,
    pattern loading, data export, and validation. Combines functionality from
    utils, load_hsdb, load_regex_registry, motif_patterns, and canonicalize_motif.

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

    Pattern Loading:
    - load_hyperscan_db(): Load pre-compiled Hyperscan databases
    - load_regex_registry(): Load regex pattern registries
    - get_motif_patterns(): Get patterns for specific motif classes

    Statistics:
    - get_basic_stats(): Calculate sequence statistics
    - quality_check_motifs(): Validate motif quality
"""

def canonicalize_motif(m):
    """
    Canonicalize motif dictionary to standard format.
    Maintains backward compatibility with Normalized_Score but doesn't require it.
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
import pickle
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


def _load_consolidated_registry():
    """Load the consolidated registry file once and cache it"""
    global _CONSOLIDATED_REGISTRY
    if _CONSOLIDATED_REGISTRY is not None:
        return _CONSOLIDATED_REGISTRY
    
    # Try to load consolidated_registry.json from current directory
    consolidated_path = "consolidated_registry.json"
    if os.path.isfile(consolidated_path):
        with open(consolidated_path, "r") as fh:
            _CONSOLIDATED_REGISTRY = json.load(fh)
            logger.info("Loaded consolidated registry from consolidated_registry.json")
            return _CONSOLIDATED_REGISTRY
    
    return None


def _load_registry(registry_dir: str, class_name: str):
    """Load registry from consolidated file"""
    # Load from consolidated registry
    consolidated = _load_consolidated_registry()
    if consolidated and "registries" in consolidated:
        if class_name in consolidated["registries"]:
            logger.debug(f"Loading {class_name} from consolidated registry")
            return consolidated["registries"][class_name]
    
    raise FileNotFoundError(f"No registry found for {class_name} in consolidated_registry.json")


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
            print(f"✓ {motif_class}: {'PASS' if matches_found > 0 else 'FAIL'}")
    
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
from typing import Dict, List, Any, Optional, Union, Tuple
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
                lines.append("  ✓ Detected:")
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
    Export motifs to CSV format with comprehensive fields
    
    Args:
        motifs: List of motif dictionaries
        filename: Optional output filename
        
    Returns:
        CSV format string
    """
    if not motifs:
        return "No motifs to export"
    
    # Comprehensive column order matching user requirements
    comprehensive_columns = [
        'ID',
        'Sequence_Name',  # Sequence Name (or Accession)
        'Source',  # Source (e.g., genome, experiment, study)
        'Class',  # Motif Class
        'Subclass',  # Motif Subclass
        'Pattern_ID',  # Pattern/Annotation ID
        'Start',  # Start Position
        'End',  # End Position
        'Length',  # Length (bp)
        'Sequence',  # Sequence
        'Method',  # Detection Method
        'Score',  # Motif Score
        'Repeat_Type',  # Repeat/Tract Type
        'Left_Arm',  # Left Arm Sequence
        'Right_Arm',  # Right Arm Sequence
        'Loop_Seq',  # Loop Sequence
        'Arm_Length',  # Arm Length
        'Loop_Length',  # Loop Length
        'Stem_Length',  # Stem Length(s)
        'Unit_Length',  # Unit/Repeat Length
        'Number_Of_Copies',  # Number of Copies/Repeats
        'Spacer_Length',  # Spacer Length
        'Spacer_Sequence',  # Spacer Sequence
        'GC_Content',  # GC Content (%)
        'Structural_Features',  # Structural Features (e.g., Tract Type, Curvature Score)
        'Strand'  # Strand information
    ]
    
    # Get all unique keys from motifs to include additional fields
    all_keys = set()
    for motif in motifs:
        all_keys.update(motif.keys())
    
    # Start with comprehensive columns, then add any additional fields found
    columns = comprehensive_columns.copy()
    for key in sorted(all_keys):
        if key not in columns:
            columns.append(key)
    
    output = StringIO()
    writer = csv.DictWriter(output, fieldnames=columns)
    writer.writeheader()
    
    for motif in motifs:
        # Create row with comprehensive field mappings
        row = {}
        for col in columns:
            value = motif.get(col, 'NA')
            
            # Map alternative field names to comprehensive columns
            if value == 'NA' or value == '' or value is None:
                # Try alternative mappings
                if col == 'Number_Of_Copies' and 'Repeat_Units' in motif:
                    value = motif['Repeat_Units']
                elif col == 'Repeat_Type' and 'Tract_Type' in motif:
                    value = motif['Tract_Type']
                elif col == 'GC_Content' and 'GC_Total' in motif:
                    value = motif['GC_Total']
                elif col == 'Structural_Features':
                    # Combine relevant structural features
                    features = []
                    if 'Tract_Type' in motif and motif['Tract_Type'] not in ['', 'NA', None]:
                        features.append(f"Tract:{motif['Tract_Type']}")
                    if 'Curvature_Score' in motif and motif['Curvature_Score'] not in ['', 'NA', None]:
                        features.append(f"Curvature:{motif['Curvature_Score']}")
                    if 'Z_Score' in motif and motif['Z_Score'] not in ['', 'NA', None]:
                        features.append(f"Z-Score:{motif['Z_Score']}")
                    value = '; '.join(features) if features else 'NA'
                
                # If still empty, set to NA
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
    """Convert motif results to pandas DataFrame with comprehensive fields"""
    if not motifs:
        return pd.DataFrame()
    
    df = pd.DataFrame(motifs)
    
    # Comprehensive column list based on user requirements
    comprehensive_columns = [
        'ID',
        'Sequence_Name',  # Sequence Name (or Accession)
        'Source',  # Source (e.g., genome, experiment, study)
        'Class',  # Motif Class
        'Subclass',  # Motif Subclass
        'Pattern_ID',  # Pattern/Annotation ID
        'Start',  # Start Position
        'End',  # End Position
        'Length',  # Length (bp)
        'Sequence',  # Sequence
        'Method',  # Detection Method
        'Score',  # Motif Score
        'Repeat_Type',  # Repeat/Tract Type
        'Left_Arm',  # Left Arm Sequence
        'Right_Arm',  # Right Arm Sequence
        'Loop_Seq',  # Loop Sequence
        'Arm_Length',  # Arm Length
        'Loop_Length',  # Loop Length
        'Stem_Length',  # Stem Length(s)
        'Unit_Length',  # Unit/Repeat Length
        'Number_Of_Copies',  # Number of Copies/Repeats
        'Spacer_Length',  # Spacer Length
        'Spacer_Sequence',  # Spacer Sequence
        'GC_Content',  # GC Content (%)
        'Structural_Features',  # Structural Features (e.g., Tract Type, Curvature Score)
        'Strand'  # Strand information
    ]
    
    # Ensure all comprehensive columns are present, fill missing with 'NA'
    for col in comprehensive_columns:
        if col not in df.columns:
            df[col] = 'NA'
    
    # Map existing fields to comprehensive column names if they differ
    column_mappings = {
        'Repeat_Units': 'Number_Of_Copies',
        'Tract_Type': 'Repeat_Type',
        'GC_Total': 'GC_Content',
        'Gc_Total': 'GC_Content',
        'Curvature_Score': 'Structural_Features',
        'Spacer': 'Spacer_Length',
        'Spacer_Seq': 'Spacer_Sequence'
    }
    
    for old_col, new_col in column_mappings.items():
        if old_col in df.columns and new_col in comprehensive_columns:
            # Only map if the new column is 'NA' (empty)
            df[new_col] = df.apply(
                lambda row: row[old_col] if pd.isna(row[new_col]) or row[new_col] == 'NA' else row[new_col],
                axis=1
            )
    
    # Fill all NaN/None values with 'NA' string
    result_df = df[comprehensive_columns].fillna('NA')
    
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
    
    print("✓ All utility tests completed")

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
                               by_class: bool = True) -> Dict[str, float]:
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
        
    Returns:
        Dictionary with density metrics (percentage, capped at 100%)
    """
    if not motifs or sequence_length == 0:
        return {'Overall': 0.0}
    
    if not by_class:
        # Overall density using set-based coverage (handles overlaps correctly)
        covered_positions = set()
        for motif in motifs:
            start = motif.get('Start', 0) - 1  # Convert to 0-based
            end = motif.get('End', 0)
            covered_positions.update(range(start, end))
        
        overall_density = min((len(covered_positions) / sequence_length) * 100, 100.0)
        return {'Overall': round(overall_density, 4)}
    
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
                                  by_class: bool = True) -> Dict[str, float]:
    """
    Calculate positional density (frequency) for motifs.
    
    Positional Density (λ) = Total count of predicted motifs / 
                             Total length (in kbp or Mbp) of analyzed region
    
    Args:
        motifs: List of motif dictionaries
        sequence_length: Total length of analyzed sequence in bp
        unit: 'kbp' or 'Mbp' for reporting units
        by_class: If True, calculate density per motif class
        
    Returns:
        Dictionary with positional density (motifs per unit)
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
    
    if not by_class:
        # Overall positional density
        overall_density = len(motifs) / sequence_length_unit
        return {'Overall': round(overall_density, 2)}
    
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
    observed_genomic_density = calculate_genomic_density(motifs, sequence_length, by_class=by_class)
    
    # Initialize background storage
    if by_class:
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
            shuffled_density = calculate_genomic_density(shuffled_motifs, sequence_length, by_class=by_class)
            
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
    
    for class_name, bg_densities in background_densities.items():
        bg_densities = [d for d in bg_densities if d is not None]
        
        if not bg_densities:
            continue
        
        obs_density = observed_genomic_density.get(class_name, 0.0)
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
        
        enrichment_results[class_name] = {
            'observed_density': round(obs_density, 4),
            'background_mean': round(bg_mean, 4),
            'background_std': round(bg_std, 4),
            'fold_enrichment': round(fold_enrichment, 2) if not np.isinf(fold_enrichment) else 'Inf',
            'p_value': round(p_value, 4),
            'n_shuffles': n_shuffles,
            'observed_count': len([m for m in motifs if m.get('Class', 'Unknown') == class_name]) if class_name != 'Overall' else len(motifs)
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