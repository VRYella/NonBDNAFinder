#!/usr/bin/env python3
"""
enhanced_aphilic_hs.py

Enhanced A-philic region detection using Hyperscan and Kadane's algorithm.

Integrates the high-performance k-mer detection approach from detect_Aphilic_hyperscan.py
with the optimal subarray detection using Kadane's algorithm from zdna_hs.py.

Key improvements:
- Uses Kadane's algorithm for optimal region detection instead of fixed 10-mer windows
- Leverages Hyperscan for high-performance tri/tetra k-mer matching
- Creates scoring arrays from k-mer propensities
- Applies nucleation filtering and combined scoring

Author: Enhanced integration for NonBDNAFinder
"""

from __future__ import annotations
import argparse
from pathlib import Path
import csv
import sys
import numpy as np
from typing import List, Tuple, Dict, Optional
import logging

# Try to import hyperscan; otherwise fallback
try:
    import hyperscan as hs
    HYPERSCAN_AVAILABLE = True
except Exception:
    HYPERSCAN_AVAILABLE = False
    hs = None

# Import from existing zdna_hs infrastructure
try:
    from zdna_hs import read_fasta, Params
except ImportError:
    # Fallback implementations
    def read_fasta(path: Path):
        """Simple FASTA reader"""
        with open(path, 'rt', encoding='utf8') as fh:
            header = None
            seq_parts = []
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if header:
                        yield header, ''.join(seq_parts).upper()
                    header = line[1:].split()[0]
                    seq_parts = []
                else:
                    seq_parts.append(line)
            if header:
                yield header, ''.join(seq_parts).upper()
    
    class Params:
        """Basic parameters class"""
        def __init__(self):
            self.threshold = 5.0
            self.drop_threshold = 10.0

class AphilicParams:
    """Parameters for A-philic detection"""
    def __init__(self):
        self.min_tetra_run = 7
        self.tri_nucleation_min = 1
        self.tri_weight = 1.0
        self.tetra_weight = 1.0
        self.threshold = 5.0  # Minimum score for Kadane's algorithm
        self.drop_threshold = 10.0  # Drop threshold for region termination
        self.min_region_length = 10  # Minimum nucleotide length of regions

def load_propensity_table(path: Path, k_expected: int) -> Dict[str, float]:
    """
    Loads CSV table with at least columns: Step, log2_odds (or Log2_Odds_laplace)
    Returns dict step -> log2_odds for steps of length k_expected.
    """
    d = {}
    with open(path, newline='', encoding='utf8') as fh:
        rdr = csv.DictReader(fh)
        for row in rdr:
            # Try multiple column name variations
            step = (row.get('Step') or row.get('step') or 
                   row.get('kmer') or row.get('Kmer'))
            if step is None:
                continue
            step = step.strip().upper()
            if len(step) != k_expected:
                continue
            
            # Try multiple log2_odds column name variations
            try:
                log2 = float(row.get('Log2_Odds_laplace') or 
                           row.get('log2_odds') or 
                           row.get('Log2_Odds') or 
                           row.get('log2') or '0')
            except (ValueError, TypeError):
                log2 = 0.0
            d[step] = log2
    return d

class HSMatcher:
    """Hyperscan matcher for k-mer patterns"""
    def __init__(self, patterns: Dict[str, float]):
        self.patterns = patterns
        self._id_to_pat = []
        self._pat_to_id = {}
        self.db = None
        
        if not patterns:
            return

        if not HYPERSCAN_AVAILABLE:
            # Will use fallback
            return

        expressions = []
        flags = []
        ids = []
        for i, pat in enumerate(patterns.keys()):
            expressions.append(pat.encode('ascii'))
            flags.append(hs.HS_FLAG_CASELESS)
            ids.append(i)
            self._id_to_pat.append(pat)
            self._pat_to_id[pat] = i

        self.db = hs.Database()
        try:
            self.db.compile(expressions=expressions, ids=ids, 
                          elements=len(expressions), flags=flags)
        except Exception as e:
            logging.warning(f"Hyperscan compile failed: {e}")
            self.db = None

    def scan(self, sequence: str):
        """Returns list of (pattern, start, end) matches"""
        matches = []
        
        if not HYPERSCAN_AVAILABLE or self.db is None:
            # Fallback to naive search
            return self._naive_scan(sequence)

        def on_match(id, start, end, flags, context):
            matches.append((self._id_to_pat[id], start, end))
            return 0

        try:
            scratch = hs.Scratch(self.db)
            self.db.scan(sequence.encode('ascii'), 
                        match_event_handler=on_match, scratch=scratch)
        except Exception as e:
            logging.warning(f"Hyperscan scan failed: {e}")
            return self._naive_scan(sequence)
        
        return matches
    
    def _naive_scan(self, sequence: str):
        """Fallback naive pattern matching"""
        matches = []
        seq = sequence.upper()
        for pattern in self.patterns.keys():
            for i in range(len(seq) - len(pattern) + 1):
                if seq[i:i+len(pattern)] == pattern:
                    matches.append((pattern, i, i+len(pattern)))
        return matches

def create_scoring_array_from_kmers(sequence: str, 
                                   tetra_matches: List[Tuple[str, int, int]],
                                   tri_matches: List[Tuple[str, int, int]],
                                   positive_tetra: Dict[str, float],
                                   positive_tri: Dict[str, float],
                                   params: AphilicParams) -> np.ndarray:
    """
    Create a scoring array where each position gets score based on overlapping k-mers.
    Uses a nucleotide-level scoring approach compatible with Kadane's algorithm.
    """
    n = len(sequence)
    if n < 4:
        return np.zeros(max(0, n-1), dtype=float)
    
    # Create transition-based scoring array (like Z-DNA)
    scoring_array = np.zeros(n - 1, dtype=float)
    
    # Map k-mer matches to positions and score them
    tetra_positions = {}
    for pattern, start, end in tetra_matches:
        score = positive_tetra.get(pattern, 0.0)
        if score > 0:
            # Score the transitions covered by this tetramer
            for pos in range(start, min(start + 3, len(scoring_array))):
                if pos >= 0:
                    scoring_array[pos] += score * params.tetra_weight
            tetra_positions[start] = score
    
    tri_positions = {}
    for pattern, start, end in tri_matches:
        score = positive_tri.get(pattern, 0.0)
        if score > 0:
            # Score the transitions covered by this trimer
            for pos in range(start, min(start + 2, len(scoring_array))):
                if pos >= 0:
                    scoring_array[pos] += score * params.tri_weight
            tri_positions[start] = score
    
    return scoring_array

def subarrays_above_threshold_aphilic(scoring_array: np.ndarray, 
                                     sequence: str, 
                                     params: AphilicParams,
                                     tetra_matches: List[Tuple[str, int, int]],
                                     tri_matches: List[Tuple[str, int, int]],
                                     positive_tetra: Dict[str, float],
                                     positive_tri: Dict[str, float]) -> List[Tuple[int, int, float, str, dict]]:
    """
    Modified Kadane's algorithm for A-philic regions with nucleation filtering.
    Returns list of (start, end, score, sequence, details) tuples.
    """
    if scoring_array.size == 0:
        return []
    
    subarrays = []
    max_ending_here = scoring_array[0]
    start_idx = end_idx = 0
    current_max = 0
    candidate_array = None
    
    for i in range(1, scoring_array.size):
        num = scoring_array[i]
        if num >= max_ending_here + num:
            start_idx = i
            end_idx = i + 1
            max_ending_here = num
        else:
            max_ending_here += num
            end_idx = i + 1
        
        if max_ending_here >= params.threshold and current_max < max_ending_here:
            # Convert transition indices to nucleotide indices
            start_nt = start_idx
            end_nt = end_idx  # end_idx is already the next position
            
            # Apply nucleation filter and get region details
            is_valid, details = apply_nucleation_filter(
                start_nt, end_nt, sequence, tetra_matches, tri_matches,
                positive_tetra, positive_tri, params)
            
            if is_valid and end_nt - start_nt >= params.min_region_length:
                substring = sequence[start_nt:end_nt + 1]
                candidate_array = (start_nt, end_nt, max_ending_here, substring, details)
                current_max = max_ending_here
        
        if (candidate_array and 
            (max_ending_here < 0 or current_max - max_ending_here >= params.drop_threshold)):
            subarrays.append(candidate_array)
            candidate_array = None
            max_ending_here = current_max = 0
    
    if candidate_array:
        subarrays.append(candidate_array)
    
    return subarrays

def apply_nucleation_filter(start_nt: int, end_nt: int, sequence: str,
                           tetra_matches: List[Tuple[str, int, int]],
                           tri_matches: List[Tuple[str, int, int]],
                           positive_tetra: Dict[str, float],
                           positive_tri: Dict[str, float],
                           params: AphilicParams) -> Tuple[bool, dict]:
    """
    Apply nucleation filtering similar to the original stringent approach.
    Returns (is_valid, details_dict).
    """
    # Count k-mers inside the region
    tetra_count = 0
    tetra_sum = 0.0
    tri_count = 0
    tri_sum = 0.0
    
    for pattern, start, end in tetra_matches:
        if start_nt <= start < end_nt:
            tetra_count += 1
            tetra_sum += positive_tetra.get(pattern, 0.0)
    
    for pattern, start, end in tri_matches:
        if start_nt <= start < end_nt:
            tri_count += 1
            tri_sum += positive_tri.get(pattern, 0.0)
    
    details = {
        "tet_pos_count": tetra_count,
        "tri_pos_count": tri_count,
        "tetra_sum": tetra_sum,
        "tri_sum": tri_sum,
        "combined_score": params.tetra_weight * tetra_sum + params.tri_weight * tri_sum,
        "length": end_nt - start_nt
    }
    
    # Apply nucleation filter
    is_valid = tri_count >= params.tri_nucleation_min
    
    return is_valid, details

def detect_aphilic_regions(sequence: str, rec_id: str,
                          tetra_table: Dict[str, float],
                          tri_table: Dict[str, float],
                          params: AphilicParams) -> List[dict]:
    """
    Main detection function using enhanced Hyperscan + Kadane approach.
    """
    # Filter positive k-mers
    positive_tetra = {k: v for k, v in tetra_table.items() if v > 0}
    positive_tri = {k: v for k, v in tri_table.items() if v > 0}
    
    if not positive_tetra and not positive_tri:
        return []
    
    # Create matchers
    tetra_matcher = HSMatcher(positive_tetra) if positive_tetra else None
    tri_matcher = HSMatcher(positive_tri) if positive_tri else None
    
    # Find k-mer matches
    tetra_matches = tetra_matcher.scan(sequence) if tetra_matcher else []
    tri_matches = tri_matcher.scan(sequence) if tri_matcher else []
    
    if not tetra_matches and not tri_matches:
        return []
    
    # Create scoring array
    scoring_array = create_scoring_array_from_kmers(
        sequence, tetra_matches, tri_matches, 
        positive_tetra, positive_tri, params)
    
    # Apply Kadane's algorithm with nucleation filtering
    regions = subarrays_above_threshold_aphilic(
        scoring_array, sequence, params, tetra_matches, tri_matches,
        positive_tetra, positive_tri)
    
    # Convert to result format
    results = []
    for start, end, score, seq_str, details in regions:
        results.append({
            "Chromosome": rec_id,
            "Start": start,
            "End": end,
            "length": details["length"],
            "tet_pos_count": details["tet_pos_count"],
            "tri_pos_count": details["tri_pos_count"],
            "tetra_sum": f"{details['tetra_sum']:.6f}",
            "tri_sum": f"{details['tri_sum']:.6f}",
            "combined_score": f"{details['combined_score']:.6f}",
            "sequence": seq_str
        })
    
    return results

def main():
    parser = argparse.ArgumentParser(
        description="Enhanced A-philic detection using Hyperscan and Kadane's algorithm")
    parser.add_argument("--fasta", required=True, help="FASTA file with sequences")
    parser.add_argument("--tetra_table", required=True, help="CSV table of 4-mers with log2_odds column")
    parser.add_argument("--tri_table", required=True, help="CSV table of 3-mers with log2_odds column")
    parser.add_argument("--outdir", default="enhanced_aphilic_out", help="Output directory")
    parser.add_argument("--min_tetra_run", type=int, default=7,
                        help="Minimum consecutive positive tetramer starts (legacy parameter)")
    parser.add_argument("--tri_nucleation_min", type=int, default=1,
                        help="Minimum number of positive trimers for nucleation")
    parser.add_argument("--tri_weight", type=float, default=1.0, help="Weight for tri-sum")
    parser.add_argument("--tetra_weight", type=float, default=1.0, help="Weight for tetra-sum")
    parser.add_argument("--threshold", type=float, default=5.0, help="Minimum score threshold for Kadane's algorithm")
    parser.add_argument("--drop_threshold", type=float, default=10.0, help="Drop threshold for region termination")
    parser.add_argument("--min_region_length", type=int, default=10, help="Minimum region length in nucleotides")
    
    args = parser.parse_args()
    
    # Setup parameters
    params = AphilicParams()
    params.min_tetra_run = args.min_tetra_run
    params.tri_nucleation_min = args.tri_nucleation_min
    params.tri_weight = args.tri_weight
    params.tetra_weight = args.tetra_weight
    params.threshold = args.threshold
    params.drop_threshold = args.drop_threshold
    params.min_region_length = args.min_region_length
    
    # Load propensity tables
    try:
        tetra_table = load_propensity_table(Path(args.tetra_table), 4)
        tri_table = load_propensity_table(Path(args.tri_table), 3)
        print(f"Loaded {len(tetra_table)} tetramers, {len(tri_table)} trimers")
        
        positive_tetra = {k: v for k, v in tetra_table.items() if v > 0}
        positive_tri = {k: v for k, v in tri_table.items() if v > 0}
        print(f"Positive: {len(positive_tetra)} tetramers, {len(positive_tri)} trimers")
        
    except Exception as e:
        print(f"Error loading propensity tables: {e}")
        return 1
    
    # Setup output
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    fasta_path = Path(args.fasta)
    out_csv = outdir / f"{fasta_path.stem}_enhanced_aphilic_regions.csv"
    
    # Process sequences
    fieldnames = ["Chromosome", "Start", "End", "length", "tet_pos_count", "tri_pos_count",
                  "tetra_sum", "tri_sum", "combined_score", "sequence"]
    
    with open(out_csv, "w", newline='', encoding='utf8') as outfh:
        writer = csv.DictWriter(outfh, fieldnames=fieldnames)
        writer.writeheader()
        
        total_regions = 0
        for rec_id, sequence in read_fasta(fasta_path):
            if len(sequence) < params.min_region_length:
                continue
                
            regions = detect_aphilic_regions(sequence, rec_id, tetra_table, tri_table, params)
            
            for region in regions:
                writer.writerow(region)
                total_regions += 1
    
    print(f"Detection complete. Found {total_regions} regions.")
    print(f"Results written to: {out_csv}")
    
    if not HYPERSCAN_AVAILABLE:
        print("Note: Hyperscan not available — used slower fallback matcher.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())