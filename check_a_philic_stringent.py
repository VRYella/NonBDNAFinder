#!/usr/bin/env python3
"""
check_a_philic_stringent.py

Stringent A-philic 10-mer detection using tri/tetra propensity tables.

Usage example:
  python check_a_philic_stringent.py --tetra tetra.csv --tri tri.csv \
    --seq "GATCGAAAATTTTATTTAAATTTAAATTTGGGTTA" --out results.tsv

The script's default behavior is stricter than 'all tetrasteps > 0':
- requires all tetrasteps AND tristeps inside a 10-mer to be > min_step_value (default 0.0)
- requires summed tetra and tri propensities across the 10-mer to exceed min_sum_tetra/min_sum_tri
- requires at least one trinucleotide nucleation peak >= nucleation_peak

CSV input (tri/tetra) accepted as comma or tab separated with columns: Step,Log2_Odds
Missing steps are treated as 0.0 by default.

Clean Steps:
1. Load propensity tables from CSV files
2. Score sequence with tri/tetra propensities  
3. Apply stringent 10-mer validation checks
4. Extend valid regions rightward while maintaining criteria
5. Report non-overlapping maximal regions
"""

from pathlib import Path
import argparse
import csv
from typing import Dict, List, Tuple
import sys

# ---------- Step 1: Utilities for loading and scoring ----------
def load_propensity_table(path: Path) -> Dict[str, float]:
    """
    Step 1a: Load tri/tetra propensity table from CSV file.
    
    Supports comma-separated or whitespace-separated formats.
    Expected columns: Step, Log2_Odds
    Missing steps default to 0.0 score.
    """
    prop = {}
    with path.open() as fh:
        # attempt to sniff delimiter: treat whitespace or comma/tab
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            # split by comma first, otherwise whitespace
            if "," in line:
                parts = [p.strip() for p in line.split(",") if p.strip()]
            else:
                parts = line.split()
            if len(parts) < 2:
                continue
            step = parts[0].upper()
            try:
                val = float(parts[1])
            except ValueError:
                # skip header or invalid lines
                continue
            prop[step] = val
    return prop

def score_sequence_steps(seq: str, tetra: Dict[str,float], tri: Dict[str,float]) -> Tuple[List[float], List[float]]:
    """
    Step 1b: Score all overlapping tetra/tri steps in sequence.
    
    Returns:
        tetra_scores: List of tetranuc-step scores (len L-3)
        tri_scores: List of trinuc-step scores (len L-2)
    
    Missing steps in propensity tables are scored as 0.0.
    """
    seq = seq.upper()
    L = len(seq)
    tetra_scores = [0.0] * max(0, L - 3)
    tri_scores = [0.0] * max(0, L - 2)
    for i in range(0, L - 3):
        step = seq[i:i+4]
        tetra_scores[i] = tetra.get(step, 0.0)
    for i in range(0, L - 2):
        step = seq[i:i+3]
        tri_scores[i] = tri.get(step, 0.0)
    return tetra_scores, tri_scores

# ---------- Step 2: Stringent 10-mer validation ----------
def tenmer_checks(start: int,
                  tetra_scores: List[float],
                  tri_scores: List[float],
                  seq_len: int,
                  min_step_value: float,
                  min_sum_tetra: float,
                  min_sum_tri: float,
                  nucleation_peak: float) -> Tuple[bool, dict]:
    """
    Step 2: Apply stringent validation checks to a 10-mer starting at `start`.
    
    Validation criteria:
    - All tetrasteps AND tristeps inside 10-mer must be > min_step_value
    - Sum of tetra propensities across 10-mer must exceed min_sum_tetra
    - Sum of tri propensities across 10-mer must exceed min_sum_tri  
    - At least one trinucleotide peak must be >= nucleation_peak
    
    Returns:
        pass_bool: True if all criteria met
        details_dict: Diagnostic values for each test
    """
    L = seq_len
    if start < 0 or start + 10 > L:
        return False, {}

    # Indices of tet steps inside 10-mer: start .. start+6 (7 values)
    tet_window = tetra_scores[start : start + 7] if start + 7 <= len(tetra_scores) else []
    # Indices of tri steps inside 10-mer: start .. start+7 (8 values)
    tri_window = tri_scores[start : start + 8] if start + 8 <= len(tri_scores) else []

    details = {
        "tet_count": len(tet_window),
        "tri_count": len(tri_window),
        "tet_min": min(tet_window) if tet_window else None,
        "tri_min": min(tri_window) if tri_window else None,
        "sum_tet": sum(tet_window) if tet_window else 0.0,
        "sum_tri": sum(tri_window) if tri_window else 0.0,
        "nucleation_max_tri": max(tri_window) if tri_window else 0.0
    }

    # Requirement 1: step-wise positivity threshold for every tet and tri
    if not tet_window or not tri_window:
        return False, details
    if any(x <= min_step_value for x in tet_window):
        return False, details
    if any(x <= min_step_value for x in tri_window):
        return False, details

    # Requirement 2: summed thresholds
    if details["sum_tet"] < min_sum_tetra:
        return False, details
    if details["sum_tri"] < min_sum_tri:
        return False, details

    # Requirement 3: nucleation peak (at least one tri >= nucleation_peak)
    if details["nucleation_max_tri"] < nucleation_peak:
        return False, details

    return True, details

# ---------- Step 3: Region expansion and Step 4: Sequence scanning ----------
def extend_region_right(start: int, end: int, seq_len: int,
                        tetra_scores: List[float],
                        tri_scores: List[float],
                        **check_kwargs) -> Tuple[int, int]:
    """
    Step 3: Extend region rightward while maintaining stringent criteria.
    
    Given an initial region [start,end) where end = start+10, attempt to extend 
    rightwards one base at a time while ensuring EVERY 10-mer fully inside the 
    extended region meets the tenmer_checks.
    
    Returns: Updated (start, end) coordinates
    """
    # current region spans positions [start, end)
    cur_end = end
    while cur_end < seq_len:
        # when extending by 1 base, the new 10-mer to check is the one that ends at cur_end (start index = cur_end-9)
        new_t10_start = cur_end - 9
        # we must ensure new_t10_start >= start (i.e., the new 10-mer is inside region after extension)
        # but for our rule we require EVERY 10-mer inside the region passes; checking the newly added 10-mer suffices because previous 10-mers already passed
        pass_new, _ = tenmer_checks(new_t10_start, tetra_scores, tri_scores, seq_len, **check_kwargs)
        if pass_new:
            cur_end += 1
        else:
            break
    return start, cur_end

def scan_sequence_for_regions(seq_id: str,
                              seq: str,
                              tetra: Dict[str,float],
                              tri: Dict[str,float],
                              min_step_value: float,
                              min_sum_tetra: float,
                              min_sum_tri: float,
                              nucleation_peak: float,
                              min_region_length: int = 10) -> List[dict]:
    """
    Step 4: Scan sequence for stringent A-philic regions.
    
    Process:
    1. Score sequence with tri/tetra propensities
    2. Check each potential 10-mer start position
    3. Apply stringent validation criteria
    4. Extend valid regions rightward maximally
    5. Report non-overlapping regions
    
    Returns:
        List of region dicts with keys: chrom,start,end,length,sequence,
        sum_tet,sum_tri,tet_min,tri_min,nucleation_max_tri
    """
    seq = seq.upper()
    seq_len = len(seq)
    tetra_scores, tri_scores = score_sequence_steps(seq, tetra, tri)
    results = []
    i = 0
    check_kwargs = {
        "min_step_value": min_step_value,
        "min_sum_tetra": min_sum_tetra,
        "min_sum_tri": min_sum_tri,
        "nucleation_peak": nucleation_peak
    }
    while i <= seq_len - min_region_length:
        # check 10-mer starting at i
        ok, details = tenmer_checks(i, tetra_scores, tri_scores, seq_len, **check_kwargs)
        if ok:
            # extend rightwards as much as possible (keeps non-overlap)
            start = i
            end = i + 10
            start, end = extend_region_right(start, end, seq_len, tetra_scores, tri_scores, **check_kwargs)
            region_seq = seq[start:end]
            region_tri_sum = sum(tri_scores[start:start + max(0, (end - start) - 2)])  # sum tri steps fully inside region
            region_tet_sum = sum(tetra_scores[start:start + max(0, (end - start) - 3)])  # sum tet steps fully inside region
            # recompute region-level details for reporting
            region_details = {
                "chrom": seq_id,
                "start": start,
                "end": end,
                "length": end - start,
                "sequence": region_seq,
                "sum_tet": round(region_tet_sum, 6),
                "sum_tri": round(region_tri_sum, 6),
                "tet_min": details.get("tet_min"),
                "tri_min": details.get("tri_min"),
                "nucleation_max_tri": details.get("nucleation_max_tri")
            }
            results.append(region_details)
            i = end  # enforce non-overlap
        else:
            i += 1
    return results

# ---------- Step 5: Command-line interface and main execution ----------
def parse_args():
    """Step 5a: Parse command-line arguments for stringent A-philic detection."""
    p = argparse.ArgumentParser(description="Stringent A-philic region detection (contiguous-10 rule; adjustable thresholds).")
    p.add_argument('--tetra', required=True, type=Path, help='Tetranucleotide propensity CSV (Step,Log2_Odds)')
    p.add_argument('--tri', required=True, type=Path, help='Trinucleotide propensity CSV (Step,Log2_Odds)')
    p.add_argument('--seq', required=False, help='Single sequence to check (letters ACGT).')
    p.add_argument('--fasta', required=False, type=Path, help='FASTA file of sequences to scan (multi-FASTA allowed).')
    p.add_argument('--out', required=False, type=Path, help='Output TSV (if scanning FASTA).')
    # Stringency params (tweakable)
    p.add_argument('--min_step_value', type=float, default=0.0, help='Minimum value for every tri/tetra step inside 10-mer (default 0.0).')
    p.add_argument('--min_sum_tetra', type=float, default=5.0, help='Minimum sum of tetra propensities across 10-mer (default 5.0).')
    p.add_argument('--min_sum_tri', type=float, default=4.0, help='Minimum sum of tri propensities across 10-mer (default 4.0).')
    p.add_argument('--nucleation_peak', type=float, default=1.5, help='Minimum single tri-step peak inside 10-mer to call nucleation (default 1.5).')
    p.add_argument('--min_region_len', type=int, default=10, help='Minimum reported region length (default 10).')
    return p.parse_args()

def main():
    """
    Step 5b: Main execution function.
    
    Execution flow:
    1. Parse command-line arguments
    2. Load propensity tables
    3. Process input sequence(s) 
    4. Apply stringent detection criteria
    5. Output results in TSV format
    """
    args = parse_args()
    tetra = load_propensity_table(args.tetra)
    tri = load_propensity_table(args.tri)

    if not args.seq and not args.fasta:
        print("Error: provide either --seq or --fasta", file=sys.stderr)
        sys.exit(1)

    all_results = []
    if args.seq:
        s = args.seq.strip().upper()
        if any(c not in "ACGT" for c in s):
            print("Warning: seq contains non-ACGT characters (they will be left as-is and steps missing in tables treated as 0.0).", file=sys.stderr)
        regs = scan_sequence_for_regions("query_seq", s, tetra, tri,
                                         min_step_value=args.min_step_value,
                                         min_sum_tetra=args.min_sum_tetra,
                                         min_sum_tri=args.min_sum_tri,
                                         nucleation_peak=args.nucleation_peak,
                                         min_region_length=args.min_region_len)
        all_results.extend(regs)

    if args.fasta:
        if not args.out:
            print("Error: when scanning a FASTA, provide --out for TSV results.", file=sys.stderr)
            sys.exit(1)
        out_fh = args.out.open("w")
        header = "\t".join(["chrom","start","end","length","sum_tet","sum_tri","tet_min","tri_min","nucleation_max_tri","sequence"])
        out_fh.write(header + "\n")
        # streaming FASTA
        try:
            from Bio import SeqIO
        except Exception:
            print("Biopython required to parse FASTA. Install with `pip install biopython`.", file=sys.stderr)
            sys.exit(1)
        for rec in SeqIO.parse(str(args.fasta), "fasta"):
            seqid = rec.id
            seqstr = str(rec.seq).upper()
            regs = scan_sequence_for_regions(seqid, seqstr, tetra, tri,
                                             min_step_value=args.min_step_value,
                                             min_sum_tetra=args.min_sum_tetra,
                                             min_sum_tri=args.min_sum_tri,
                                             nucleation_peak=args.nucleation_peak,
                                             min_region_length=args.min_region_len)
            for r in regs:
                out_fh.write("\t".join([
                    r["chrom"], str(r["start"]), str(r["end"]), str(r["length"]),
                    f"{r['sum_tet']:.4f}", f"{r['sum_tri']:.4f}",
                    str(r['tet_min']), str(r['tri_min']), str(r['nucleation_max_tri']),
                    r['sequence']
                ]) + "\n")
        out_fh.close()
        print(f"Wrote results to {args.out}")

    else:
        # print results to stdout if single seq
        if not all_results:
            print("No A-philic regions detected under current stringent parameters.")
        else:
            print("chrom\tstart\tend\tlength\tsum_tet\tsum_tri\ttet_min\ttri_min\tnucleation_max_tri\tsequence")
            for r in all_results:
                print("\t".join([
                    r["chrom"], str(r["start"]), str(r["end"]), str(r["length"]),
                    f"{r['sum_tet']:.4f}", f"{r['sum_tri']:.4f}",
                    str(r['tet_min']), str(r['tri_min']), str(r['nucleation_max_tri']),
                    r['sequence']
                ]))

if __name__ == "__main__":
    main()