#!/usr/bin/env python3
"""
zdna_hs.py

Single-file Z-DNA extractor & annotator using Intel Hyperscan for transition detection.

Author: Adapted for user request (based on prior code snippet).
Date: 2025-09-09

Requirements:
  - Python 3.8+
  - numpy, pandas, biopython, hyperscan, termcolor
  - Hyperscan C library installed
"""

from __future__ import annotations
import argparse
import gzip
import os
from pathlib import Path
import time
import sys
from typing import Iterable, Iterator, Optional, Tuple, List, Dict, Any
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import pandas as pd
import concurrent.futures
import multiprocessing as mp
import logging
import bisect
from termcolor import colored

# -------------------------
# Hyperscan availability
# -------------------------
try:
    import hyperscan
    HYPERSCAN_AVAILABLE = True
except Exception:
    hyperscan = None
    HYPERSCAN_AVAILABLE = False

# Globals for per-process Hyperscan DB
_hs_db = None
_hs_id_to_key = None

_HS_TRANSITIONS = {"GC": 1, "CG": 2, "GT": 3, "TG": 4, "AC": 5, "CA": 6, "AT": 7, "TA": 8}

def _build_hyperscan_db():
    """Compile Hyperscan DB for dinucleotide transitions in the current process."""
    global _hs_db, _hs_id_to_key
    if not HYPERSCAN_AVAILABLE:
        _hs_db = None
        _hs_id_to_key = None
        return
    if _hs_db is not None:
        return
    db = hyperscan.Database()
    expressions = []
    ids = []
    flags = []
    for key, pid in _HS_TRANSITIONS.items():
        expressions.append(key.encode("ascii"))
        ids.append(pid)
        flags.append(hyperscan.HS_FLAG_UTF8)
    try:
        db.compile(expressions=expressions, ids=ids, flags=flags)
        _hs_db = db
        _hs_id_to_key = {v: k for k, v in _HS_TRANSITIONS.items()}
    except Exception as e:
        print("Hyperscan compile failed: ", e)
        _hs_db = None
        _hs_id_to_key = None

# -------------------------
# Params dataclass
# -------------------------
class Params:
    # default numerical parameters (inferred / sensible defaults)
    GC_weight: float = 3.0
    AT_weight: float = 1.0
    GT_weight: float = 2.0
    AC_weight: float = 2.0
    cadence_reward: float = 0.2
    mismatch_penalty_starting_value: int = 3
    mismatch_penalty_linear_delta: int = 1
    mismatch_penalty_type: str = "linear"  # or "exponential"
    mismatch_penalty_choices = ("linear", "exponential")
    n_jobs: int = -1
    threshold: float = 5.0
    consecutive_AT_scoring: Tuple[float, ...] = (3.0, 1.5, 0.7)
    display_sequence_score: int = 0
    output_dir: str = "zdna_extractions"
    method: str = "transitions"  # "transitions", "coverage", "layered"
    drop_threshold: float = 50.0
    total_sequence_scoring: bool = False

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            if hasattr(self, k):
                setattr(self, k, v)
        # Derived dict for CLI print
        self.__new_dict__ = {k: getattr(self, k) for k in dir(self) if not k.startswith("_") and not callable(getattr(self, k))}

# -------------------------
# FASTA reader
# -------------------------
def read_fasta(fasta: Path | str) -> Iterator[Tuple[str, str]]:
    p = Path(fasta)
    if p.name.endswith(".gz"):
        fh = gzip.open(p, 'rt')
    else:
        fh = open(p, encoding='utf-8', mode='r')
    with fh:
        for name, seq in SimpleFastaParser(fh):
            yield name.split(" ")[0], seq
    # file closed by context

# -------------------------
# Utility / timing decorator
# -------------------------
def timeit(func):
    def wrapper(*args, **kwargs):
        print(colored("Z-DNA extraction initialized.", "magenta"))
        t0 = time.perf_counter()
        res = func(*args, **kwargs)
        t1 = time.perf_counter()
        print(colored(f"Process finished within {t1-t0:.2f} second(s).", "magenta"))
        return res
    return wrapper

# -------------------------
# Original Python per-transition scorer (fallback)
# -------------------------
def _zdna_transitions_python(s: str, params: Params) -> np.ndarray:
    """Pure-Python scoring array identical to original semantics."""
    s = s.upper()
    n = len(s)
    if n < 2:
        return np.empty(0, dtype=float)
    scoring_array = np.empty(n - 1, dtype=float)
    mismatches_counter = 0
    consecutive_AT_counter = 0
    for i in range(n - 1):
        transition = s[i] + s[i+1]
        if transition in ("GC", "CG"):
            scoring_array[i] = params.GC_weight
            mismatches_counter = 0
            consecutive_AT_counter = 0
        elif transition in ("GT", "TG"):
            scoring_array[i] = params.GT_weight
            mismatches_counter = 0
            consecutive_AT_counter = 0
        elif transition in ("AC", "CA"):
            scoring_array[i] = params.AC_weight
            mismatches_counter = 0
            consecutive_AT_counter = 0
        elif transition in ("AT", "TA"):
            adjusted_weight = params.AT_weight
            if consecutive_AT_counter < len(params.consecutive_AT_scoring):
                adjusted_weight += params.consecutive_AT_scoring[consecutive_AT_counter]
            else:
                adjusted_weight += params.consecutive_AT_scoring[-1]
            scoring_array[i] = adjusted_weight
            consecutive_AT_counter += 1
            mismatches_counter = 0
        else:
            mismatches_counter += 1
            consecutive_AT_counter = 0
            if params.mismatch_penalty_type == "exponential":
                scoring_array[i] = - (params.mismatch_penalty_starting_value ** mismatches_counter) if mismatches_counter < 15 else -32000
            elif params.mismatch_penalty_type == "linear":
                scoring_array[i] = - (params.mismatch_penalty_starting_value + params.mismatch_penalty_linear_delta * (mismatches_counter - 1))
            else:
                raise ValueError(f"Unknown mismatch_penalty_type: {params.mismatch_penalty_type}")
        # cadence reward
        if transition in ("GC","CG","GT","TG","AC","CA","AT","TA"):
            scoring_array[i] += params.cadence_reward
    return scoring_array

# -------------------------
# Hyperscan-based transitions scorer
# -------------------------
def zdna_scoring_with_hyperscan(sequence: str, params: Params) -> np.ndarray:
    """
    Use Hyperscan to detect the 2-base transitions; then apply identical
    post-processing (consecutive AT scoring, mismatch penalties, cadence reward).
    """
    if not HYPERSCAN_AVAILABLE:
        return _zdna_transitions_python(sequence, params)

    # ensure DB present
    global _hs_db, _hs_id_to_key
    if _hs_db is None:
        _build_hyperscan_db()
    if _hs_db is None:
        # failed build -> fallback
        return _zdna_transitions_python(sequence, params)

    s = sequence.upper()
    n = len(s)
    if n < 2:
        return np.empty(0, dtype=float)

    scoring_array = np.zeros(n - 1, dtype=float)
    matched_mask = np.zeros(n - 1, dtype=bool)
    matches: List[Tuple[int, int]] = []  # (start_index, pid)

    def on_match(hs_id, start, end, flags, context):
        # For 2-mer patterns, start index corresponds to transition index
        matches.append((int(start), int(hs_id)))
        return 0

    _hs_db.scan(s.encode("ascii"), match_event_handler=on_match)

    # Fill matches into scoring_array (last assignment to index wins)
    for start_idx, pid in matches:
        key = _hs_id_to_key.get(pid)
        if key is None:
            continue
        if key in ("GC","CG"):
            weight = params.GC_weight
        elif key in ("GT","TG"):
            weight = params.GT_weight
        elif key in ("AC","CA"):
            weight = params.AC_weight
        elif key in ("AT","TA"):
            weight = params.AT_weight
        else:
            weight = 0.0
        scoring_array[start_idx] = weight
        matched_mask[start_idx] = True

    # cadence reward
    if getattr(params, "cadence_reward", 0):
        scoring_array[matched_mask] += params.cadence_reward

    # consecutive AT scoring
    at_positions = [idx for (idx, pid) in matches if _hs_id_to_key.get(pid) in ("AT","TA")]
    if at_positions:
        at_positions = sorted(set(at_positions))
        groups = []
        cur = [at_positions[0]]
        for idx in at_positions[1:]:
            if idx == cur[-1] + 1:
                cur.append(idx)
            else:
                groups.append(cur)
                cur = [idx]
        groups.append(cur)
        for grp in groups:
            for j, pos in enumerate(grp):
                add = params.consecutive_AT_scoring[j] if j < len(params.consecutive_AT_scoring) else params.consecutive_AT_scoring[-1]
                scoring_array[pos] += add

    # mismatch penalties for unmatched regions
    unmatched = ~matched_mask
    if unmatched.any():
        u = unmatched.astype(int)
        diffs = np.diff(np.concatenate(([0], u, [0])))
        starts = np.where(diffs == 1)[0]
        ends = np.where(diffs == -1)[0]
        for s_idx, e_idx in zip(starts, ends):
            run_len = e_idx - s_idx
            if run_len <= 0:
                continue
            if params.mismatch_penalty_type == "exponential":
                kidx = np.arange(1, run_len + 1, dtype=float)
                penalties = -np.minimum(params.mismatch_penalty_starting_value ** kidx, 32000.0)
            else:
                kidx = np.arange(0, run_len, dtype=float)
                penalties = -(params.mismatch_penalty_starting_value + params.mismatch_penalty_linear_delta * kidx)
            scoring_array[s_idx:e_idx] = penalties

    return scoring_array

# -------------------------
# Layered & coverage wrappers
# -------------------------
def zdna_scoring_layered(sequence: str, params: Params) -> np.ndarray:
    """
    Layered: use either hyperscan transitions then apply layered adjustments (same cadence).
    This implementation uses the transitions scoring as base (HS or python) then returns it.
    If additional layered adjustments required, they can be added here.
    """
    base = zdna_scoring_with_hyperscan(sequence, params) if HYPERSCAN_AVAILABLE else _zdna_transitions_python(sequence, params)
    # (original code added cadence and other small adjustments already)
    return base

# -------------------------
# Subarray detection (Kadane-like with drop_threshold)
# -------------------------
def subarrays_above_threshold_from_scoring(scoring_array: np.ndarray, s: str, params: Params) -> List[Tuple[int,int,float,str]]:
    """
    Given scoring_array (len = n-1) and original sequence s (len n),
    return list of tuples (start, end, score, substring) matching original semantics:
      - start and end are indices in sequence (start..end inclusive?)
    In original code candidate_array was (start_idx, end_idx, score, substring) where substring = self[start_idx: end_idx + 1]
    We will interpret start=end indexes consistent with that approach.
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
            # in original they built substring as self[start_idx: end_idx + 1]
            # map to sequence substring: transitions index start_idx corresponds to bases start_idx .. end_idx+1
            substring = s[start_idx: end_idx + 1]
            candidate_array = (start_idx, end_idx, max_ending_here, substring)
            current_max = max_ending_here

        if candidate_array and (max_ending_here < 0 or current_max - max_ending_here >= params.drop_threshold):
            subarrays.append(candidate_array)
            candidate_array = None
            max_ending_here = current_max = 0

    if candidate_array:
        subarrays.append(candidate_array)
    # Convert start/end from transition indices to base indices for consistency with earlier code if needed.
    # Earlier code seems to use start,end directly as indices into sequence in returned DataFrame.
    # We'll return them as base indices: start (base) = start_idx, end (base) = end_idx + 1
    converted = []
    for st, ed, sc, substr in subarrays:
        start_base = st
        end_base = ed + 1  # because substring used end_idx+1
        converted.append((start_base, end_base, sc, substr))
    return converted

# -------------------------
# ZDNACalculatorSeq-like behavior
# -------------------------
class ZDNACalculatorSeq:
    def __init__(self, data: str, params: Optional[Params] = None):
        self.seq = data.upper()
        self.params = Params() if params is None else params
        self.scoring_array: np.ndarray = np.empty(0, dtype=float)

    def zdna_calculator_transitions(self) -> np.ndarray:
        if HYPERSCAN_AVAILABLE:
            try:
                # ensure DB in this worker
                if _hs_db is None:
                    _build_hyperscan_db()
                arr = zdna_scoring_with_hyperscan(self.seq, self.params)
                self.scoring_array = arr
                return arr
            except Exception as e:
                logging.warning("Hyperscan scoring failed; falling back to python scorer: %s", e)
                arr = _zdna_transitions_python(self.seq, self.params)
                self.scoring_array = arr
                return arr
        else:
            arr = _zdna_transitions_python(self.seq, self.params)
            self.scoring_array = arr
            return arr

    def zdna_calculator_layered(self) -> np.ndarray:
        arr = zdna_scoring_layered(self.seq, self.params)
        self.scoring_array = arr
        return arr

    def zdna_calculator_coverage(self) -> np.ndarray:
        # For completeness include a simple placeholder similar to original coverage logic,
        # but we will fallback to layered if not implemented fully.
        return self.zdna_calculator_transitions()

    def subarrays_above_threshold(self) -> List[Tuple[int,int,float,str]]:
        method = self.params.method
        if method == "transitions":
            self.scoring_array = self.zdna_calculator_transitions()
        elif method == "coverage":
            self.scoring_array = self.zdna_calculator_coverage()
        elif method == "layered":
            self.scoring_array = self.zdna_calculator_layered()
        else:
            raise ValueError(f"Method {method} not supported")
        return subarrays_above_threshold_from_scoring(self.scoring_array, self.seq, self.params)

# -------------------------
# Extraction helpers for file-level operations (similar to original)
# -------------------------
def _extract(ID: str, sequence: str, params: Params) -> pd.DataFrame:
    zcalc = ZDNACalculatorSeq(data=sequence, params=params)
    subarrays = zcalc.subarrays_above_threshold()
    scoring_array = zcalc.scoring_array
    headers = ["Chromosome", "Start", "End", "Z-DNA Score", "Sequence", "totalSequenceScore"]
    if not subarrays:
        # produce a NaN row with totalSequenceScore
        total = float(np.sum(scoring_array)) if scoring_array.size>0 else 0.0
        df = pd.DataFrame([[ID, np.nan, np.nan, np.nan, np.nan, total]], columns=headers)
        return df

    rows = []
    total = float(np.sum(scoring_array))
    for (st, ed, sc, substr) in subarrays:
        rows.append([ID, int(st), int(ed), float(sc), substr, total])
    df = pd.DataFrame(rows, columns=headers)
    return df

def _subextraction(submission_forms: List[Dict[str,Any]], params: Params) -> pd.DataFrame:
    outputs = []
    for form in submission_forms:
        try:
            df = _extract(ID=form["recordID"], sequence=form["recordSeq"], params=params)
            outputs.append(df)
        except Exception as e:
            logging.error("Failed processing %s: %s", form.get("recordID","?"), e)
            continue
    if outputs:
        return pd.concat(outputs, axis=0).reset_index(drop=True)
    else:
        return pd.DataFrame(columns=["Chromosome","Start","End","Z-DNA Score","Sequence","totalSequenceScore"])

def assign_tasks(tasks: List[Any], total_buckets: int) -> List[List[Any]]:
    total = len(tasks)
    if total_buckets <= 0:
        total_buckets = 1
    if total_buckets > total:
        total_buckets = total
    step = total // total_buckets
    remainder = total % total_buckets
    assigned = []
    inf = 0
    for b in range(total_buckets):
        add = step + (1 if remainder>0 else 0)
        assigned.append(tasks[inf:inf+add])
        inf += add
        if remainder>0:
            remainder -= 1
    return assigned

def extract_zdna_v2(fasta: str | Path, params: Params) -> pd.DataFrame:
    t0 = time.perf_counter()
    forms = []
    for rec_id, rec_seq in read_fasta(fasta):
        forms.append({"recordID": rec_id, "recordSeq": rec_seq.upper(), "input_fasta": str(fasta)})
    if params.total_sequence_scoring:
        rows = []
        for f in forms:
            calc = ZDNACalculatorSeq(f["recordSeq"], params)
            scoring_arr = calc.zdna_calculator_transitions()
            total = float(np.sum(scoring_arr))
            rows.append({"Chromosome": f["recordID"], "Start": 0, "End": len(f["recordSeq"])-1, "Z-DNA Score": total, "Sequence": f["recordSeq"]})
        df = pd.DataFrame(rows, columns=["Chromosome","Start","End","Z-DNA Score","Sequence"])
        print(colored(f"Process finished within {time.perf_counter()-t0:.2f} second(s).", "magenta"))
        return df

    n_jobs = params.n_jobs if params.n_jobs>0 else mp.cpu_count()
    assigned = assign_tasks(forms, n_jobs)
    outputs = []
    # ensure Hyperscan DB built per worker
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_jobs, mp_context=mp.get_context("spawn"), initializer=_build_hyperscan_db) as executor:
        futures = [executor.submit(_subextraction, chunk, params) for chunk in assigned]
        for fut in concurrent.futures.as_completed(futures):
            try:
                res = fut.result()
                if res is not None and not res.empty:
                    outputs.append(res)
            except Exception as e:
                logging.error("Worker failed: %s", e)
    if outputs:
        df = pd.concat(outputs, axis=0).reset_index(drop=True)
    else:
        df = pd.DataFrame(columns=["Chromosome","Start","End","Z-DNA Score","Sequence","totalSequenceScore"])
    print(colored(f"Process finished within {time.perf_counter()-t0:.2f} second(s).", "magenta"))
    return df

# -------------------------
# CLI & main
# -------------------------
@timeit
def transform(path: Path | str, params: Params, gff_file: Optional[Path] = None) -> pd.DataFrame:
    path = Path(path)
    path = path.resolve()
    outdir = Path(params.output_dir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    zdna_df = extract_zdna_v2(path, params)
    # Save CSV
    out_file = outdir.joinpath(f"{path.stem}_zdna_score.csv")
    zdna_df.to_csv(out_file, index=False)
    logging.info("Saved results to %s", out_file)
    return zdna_df

def parse_args():
    p = argparse.ArgumentParser(description="Z-DNA extraction using Hyperscan-accelerated transitions")
    p.add_argument("--fasta", type=str, required=True)
    p.add_argument("--GC_weight", type=float, default=Params().GC_weight)
    p.add_argument("--AT_weight", type=float, default=Params().AT_weight)
    p.add_argument("--GT_weight", type=float, default=Params().GT_weight)
    p.add_argument("--AC_weight", type=float, default=Params().AC_weight)
    p.add_argument("--cadence_reward", type=float, default=Params().cadence_reward)
    p.add_argument("--mismatch_penalty_starting_value", type=int, default=Params().mismatch_penalty_starting_value)
    p.add_argument("--mismatch_penalty_linear_delta", type=int, default=Params().mismatch_penalty_linear_delta)
    p.add_argument("--mismatch_penalty_type", choices=Params().mismatch_penalty_choices, default=Params().mismatch_penalty_type)
    p.add_argument("--n_jobs", type=int, default=Params().n_jobs)
    p.add_argument("--threshold", type=float, default=Params().threshold)
    p.add_argument("--consecutive_AT_scoring", type=str, default="3.0,1.5,0.7",
                   help="Comma-separated floats for consecutive AT scoring, e.g. '3.0,1.5,0.7'")
    p.add_argument("--display_sequence_score", type=int, choices=[0,1], default=Params().display_sequence_score)
    p.add_argument("--output_dir", type=str, default=Params().output_dir)
    p.add_argument("--method", choices=["transitions","coverage","layered"], default=Params().method)
    p.add_argument("--drop_threshold", type=float, default=Params().drop_threshold)
    p.add_argument("--total_sequence_scoring", action="store_true")
    return p.parse_args()

def parse_consec_at(s: str) -> Tuple[float, ...]:
    return tuple(float(x.strip()) for x in s.split(","))

def main():
    args = parse_args()
    consecutive = parse_consec_at(args.consecutive_AT_scoring)
    params = Params(
        GC_weight=args.GC_weight,
        AT_weight=args.AT_weight,
        GT_weight=args.GT_weight,
        AC_weight=args.AC_weight,
        cadence_reward=args.cadence_reward,
        mismatch_penalty_starting_value=args.mismatch_penalty_starting_value,
        mismatch_penalty_linear_delta=args.mismatch_penalty_linear_delta,
        mismatch_penalty_type=args.mismatch_penalty_type,
        n_jobs=args.n_jobs,
        threshold=args.threshold,
        consecutive_AT_scoring=consecutive,
        display_sequence_score=args.display_sequence_score,
        output_dir=args.output_dir,
        method=args.method,
        drop_threshold=args.drop_threshold,
        total_sequence_scoring=args.total_sequence_scoring
    )

    fasta = Path(args.fasta).expanduser().resolve()
    assert fasta.is_file(), f"No file {fasta} found"

    print(colored("Process parameters", "magenta"))
    for k, v in params.__new_dict__.items():
        print(colored(f"{k}: {v}", "magenta"))

    logging.basicConfig(level=logging.WARNING, format="%(levelname)s:%(message)s")
    # Build HS DB in main (child processes will be initialized by initializer)
    if HYPERSCAN_AVAILABLE:
        _build_hyperscan_db()
        if _hs_db is None:
            print(colored("Warning: Hyperscan requested but DB failed to compile. Falling back to Python scoring.", "yellow"))

    transform(path=fasta, params=params)

if __name__ == "__main__":
    main()