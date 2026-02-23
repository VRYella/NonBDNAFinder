"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ NonBScanner - Non-B DNA Motif Detection Suite                                │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘
"""
# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
import os
import warnings
import time
import threading
import logging
import bisect
import multiprocessing
import pandas as pd
from typing import List, Dict, Any, Optional, Union, Tuple, Callable, overload, Literal
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)

# Detector imports
from Detectors import CurvedDNADetector, SlippedDNADetector, CruciformDetector, RLoopDetector, TriplexDetector, GQuadruplexDetector, IMotifDetector, ZDNADetector, APhilicDetector
from Utilities.utilities import parse_fasta, read_fasta_file, validate_sequence, export_to_csv, export_to_bed, export_to_json, export_to_excel, calculate_motif_statistics, normalize_motif_scores

# Optional progress tracking support (for Streamlit UI integration)
try:
    from scientific_progress import StreamlitProgressPanel, create_streamlit_progress_callback, create_compact_progress_callback
    STREAMLIT_PROGRESS_AVAILABLE = True
except ImportError:
    STREAMLIT_PROGRESS_AVAILABLE = False

# ═══════════════════════════════════════════════════════════════════════════════
# TUNABLE PARAMETERS - All configuration values at the top for easy modification
# ═══════════════════════════════════════════════════════════════════════════════
__version__ = "2024.2"; __author__ = "Dr. Venkata Rajesh Yella"

# === DETECTOR PARALLELIZATION THRESHOLD (PR#5) ===
# CHUNK_THRESHOLD = 50KB triggers parallel detector execution (NOT sequence chunking!)
# For sequences > 50KB, runs 9 detectors in parallel using ThreadPoolExecutor (1.5-2x speedup)
# This is SEPARATE from sequence chunking threshold (1MB in config/UI)
CHUNK_THRESHOLD = 50000; DEFAULT_CHUNK_SIZE = 50000; DEFAULT_CHUNK_OVERLAP = 2000

# === SEQUENCE CHUNKING THRESHOLD ===
# SEQUENCE_CHUNKING_THRESHOLD = 0 so ALL sequences use 50KB chunks with 2KB overlap
# regardless of sequence size (RAM-only processing)
SEQUENCE_CHUNKING_THRESHOLD = 0  # Always chunk

# Parallel detector execution for maximum performance (enabled by default for sequences >50KB)
# MAX_DETECTOR_WORKERS limited to 9 because there are exactly 9 detector types in the system
USE_PARALLEL_DETECTORS = True; MAX_DETECTOR_WORKERS = min(9, os.cpu_count() or 4)  # Up to 9 detectors (one per detector type)

# === MOTIF CLUSTERING & OVERLAP PARAMETERS ===
HYBRID_MIN_OVERLAP = 0.50; HYBRID_MAX_OVERLAP = 0.99
CLUSTER_WINDOW_SIZE = 300; CLUSTER_MIN_MOTIFS = 4; CLUSTER_MIN_CLASSES = 3

# === DETECTOR NAME MAPPINGS ===
DETECTOR_DISPLAY_NAMES = {'curved_dna': 'Curved DNA', 'slipped_dna': 'Slipped DNA', 'cruciform': 'Cruciform', 'r_loop': 'R-Loop', 'triplex': 'Triplex', 'g_quadruplex': 'G-Quadruplex', 'i_motif': 'i-Motif', 'z_dna': 'Z-DNA', 'a_philic': 'A-philic DNA'}
CLASS_TO_DETECTOR = {'Curved_DNA': 'curved_dna', 'Slipped_DNA': 'slipped_dna', 'Cruciform': 'cruciform', 'R-Loop': 'r_loop', 'Triplex': 'triplex', 'G-Quadruplex': 'g_quadruplex', 'i-Motif': 'i_motif', 'Z-DNA': 'z_dna', 'A-philic_DNA': 'a_philic'}
DETECTOR_TO_CLASS = {v: k for k, v in CLASS_TO_DETECTOR.items()}
# ═══════════════════════════════════════════════════════════════════════════════

_DETECTOR_TIMINGS = {}; _TIMINGS_LOCK = threading.Lock(); _CACHED_SCANNER = None; _SCANNER_LOCK = threading.Lock()

class AnalysisProgress:
    STAGES = {'INPUT': {'name': 'Input Parsing', 'description': 'Parse FASTA and validate sequences', 'icon': '📥'}, 'DETECTION': {'name': 'Motif Detection', 'description': 'Run 9 specialized detectors', 'icon': '[DETECT]'}, 'PROCESSING': {'name': 'Post-Processing', 'description': 'Overlap removal, hybrid/cluster detection', 'icon': '[PROCESS]'}, 'OUTPUT': {'name': 'Output Generation', 'description': 'Score normalization and formatting', 'icon': '[OUTPUT]'}}
    
    def __init__(self, sequence_length: int = 0, sequence_name: str = "sequence"):
        self.sequence_length = sequence_length; self.sequence_name = sequence_name; self.start_time = time.time(); self.current_stage = None; self.stage_times = {}; self.total_motifs = 0; self.total_bp_processed = 0; self.errors = []
        self.detector_status = {name: {'status': 'pending', 'elapsed': 0.0, 'motifs': 0, 'speed': 0.0} for name in DETECTOR_DISPLAY_NAMES}
    
    def start_stage(self, stage: str) -> None:
        if stage in self.STAGES: self.current_stage = stage; self.stage_times[stage] = {'start': time.time(), 'end': None}
    
    def end_stage(self, stage: str) -> None:
        if stage in self.stage_times: self.stage_times[stage]['end'] = time.time()
    
    def update_detector(self, detector_name: str, completed: bool = False, motifs: int = 0, elapsed: float = 0.0, error: Optional[str] = None) -> None:
        if detector_name not in self.detector_status: return
        if error: self.detector_status[detector_name]['status'] = 'error'; self.errors.append(f"{detector_name}: {error}")
        elif completed: self.detector_status[detector_name]['status'] = 'complete'; self.detector_status[detector_name]['elapsed'] = elapsed; self.detector_status[detector_name]['motifs'] = motifs; self.detector_status[detector_name]['speed'] = self.sequence_length / elapsed if elapsed > 0 else 0; self.total_motifs += motifs
        else: self.detector_status[detector_name]['status'] = 'running'
    
    def get_progress_percentage(self) -> float: return (sum(1 for d in self.detector_status.values() if d['status'] in ('complete', 'error')) / len(self.detector_status)) * 100
    def get_elapsed_time(self) -> float: return time.time() - self.start_time
    def get_throughput(self) -> float: elapsed = self.get_elapsed_time(); return self.sequence_length / elapsed if elapsed > 0 else 0.0
    def get_summary(self) -> Dict[str, Any]: return {'sequence_name': self.sequence_name, 'sequence_length': self.sequence_length, 'percentage': self.get_progress_percentage(), 'elapsed': self.get_elapsed_time(), 'throughput': self.get_throughput(), 'current_stage': self.current_stage, 'detectors': self.detector_status.copy(), 'total_motifs': self.total_motifs, 'errors': self.errors.copy()}
    def format_progress_bar(self, width: int = 50) -> str: pct = self.get_progress_percentage(); filled = int((pct / 100) * width); return f"[{'█' * filled}{'░' * (width - filled)}] {pct:.1f}%"
    
    def format_detector_table(self) -> str:
        lines = ["┌─────────────────┬───────────┬─────────┬────────┬─────────────┐", "│ Detector        │ Status    │ Time    │ Motifs │ Speed (bp/s)│", "├─────────────────┼───────────┼─────────┼────────┼─────────────┤"]
        for name, display_name in DETECTOR_DISPLAY_NAMES.items():
            s = self.detector_status[name]; icon, status_text = ('[OK]', 'Complete') if s['status'] == 'complete' else ('►', 'Running') if s['status'] == 'running' else ('✗', 'Error') if s['status'] == 'error' else ('○', 'Pending')
            time_str = f"{s['elapsed']:.3f}s" if s['elapsed'] > 0 else '-'; motifs_str = str(s['motifs']) if s['status'] == 'complete' else '-'; speed_str = f"{s['speed']:,.0f}" if s['speed'] > 0 else '-'
            lines.append(f"│ {icon} {display_name:<14}│ {status_text:<9} │ {time_str:>7} │ {motifs_str:>6} │ {speed_str:>11} │")
        lines.append("└─────────────────┴───────────┴─────────┴────────┴─────────────┘"); return '\n'.join(lines)
    
    def format_pipeline_status(self) -> str:
        lines = ["Pipeline Status:", "════════════════"]
        for stage_id, stage_info in self.STAGES.items():
            if stage_id in self.stage_times and self.stage_times[stage_id]['end']: status = f"Complete ({self.stage_times[stage_id]['end'] - self.stage_times[stage_id]['start']:.2f}s)"
            elif self.current_stage == stage_id or (stage_id in self.stage_times and not self.stage_times[stage_id]['end']): status = "► Running..."
            else: status = "○ Pending"
            lines.append(f"  {stage_info['icon']} {stage_info['name']}: {status}")
        return '\n'.join(lines)

def create_progress_callback(progress: AnalysisProgress) -> Callable:
    def callback(detector_name: str, completed: int, total: int, elapsed: float, motif_count: int): progress.update_detector(detector_name, completed=True, motifs=motif_count, elapsed=elapsed)
    return callback

def create_enhanced_progress_callback(progress: AnalysisProgress, print_updates: bool = False, on_detector_complete: Optional[Callable] = None) -> Callable:
    def callback(detector_name: str, completed: int, total: int, elapsed: float, motif_count: int):
        progress.update_detector(detector_name, completed=True, motifs=motif_count, elapsed=elapsed)
        if print_updates:
            print(f"\r  {progress.format_progress_bar(40)} | {DETECTOR_DISPLAY_NAMES.get(detector_name, detector_name)}: {elapsed:.3f}s ({motif_count} motifs)", end='')
            if completed == total: print()
        if on_detector_complete: on_detector_complete(detector_name, completed, total, elapsed, motif_count)
    return callback

def get_last_detector_timings() -> Dict[str, float]:
    with _TIMINGS_LOCK: return dict(_DETECTOR_TIMINGS)

def get_detector_display_names() -> Dict[str, str]: return dict(DETECTOR_DISPLAY_NAMES)

def _update_detector_timing(detector_name: str, elapsed: float) -> None:
    global _DETECTOR_TIMINGS
    with _TIMINGS_LOCK: _DETECTOR_TIMINGS[detector_name] = elapsed

def _reset_detector_timings() -> None:
    global _DETECTOR_TIMINGS
    with _TIMINGS_LOCK: _DETECTOR_TIMINGS = {}

def _get_cached_scanner() -> 'NonBScanner':
    global _CACHED_SCANNER
    if _CACHED_SCANNER is None:
        with _SCANNER_LOCK:
            if _CACHED_SCANNER is None: _CACHED_SCANNER = NonBScanner(enable_all_detectors=True)
    return _CACHED_SCANNER

class NonBScanner:
    def __init__(self, enable_all_detectors: bool = True):
        self.detectors = {'curved_dna': CurvedDNADetector(), 'slipped_dna': SlippedDNADetector(), 'cruciform': CruciformDetector(), 'r_loop': RLoopDetector(), 'triplex': TriplexDetector(), 'g_quadruplex': GQuadruplexDetector(), 'i_motif': IMotifDetector(), 'z_dna': ZDNADetector(), 'a_philic': APhilicDetector()} if enable_all_detectors else {}
    
    def analyze_sequence(self, sequence: str, sequence_name: str = "sequence", progress_callback: Optional[Callable[[str, int, int, float, int], None]] = None, enabled_classes: Optional[List[str]] = None, use_parallel_detectors: bool = None, use_preprocessing: bool = False) -> List[Dict[str, Any]]:
        """Optimized sequence analysis with optional parallel detector execution.
        
        Args:
            sequence: DNA sequence to analyze
            sequence_name: Name identifier for the sequence
            progress_callback: Optional callback for progress updates
            enabled_classes: Optional list of motif classes to detect (None = all)
            use_parallel_detectors: Enable parallel detector execution (None = auto based on sequence length)
            use_preprocessing: Enable comprehensive preprocessing with validation (default: False for backward compatibility)
        
        Returns:
            List of detected motif dictionaries
        """
        # Optional: Use comprehensive preprocessing pipeline
        if use_preprocessing:
            from Utilities.sequence_preprocessor import preprocess_sequence
            
            preprocessing_result = preprocess_sequence(sequence)
            
            # Check for errors
            if preprocessing_result.validation_status == "error":
                error_msg = "; ".join(preprocessing_result.errors)
                raise ValueError(f"Sequence preprocessing failed: {error_msg}")
            
            # Log warnings (e.g., N's present)
            if preprocessing_result.warnings:
                for warning in preprocessing_result.warnings:
                    logger.warning(f"Sequence '{sequence_name}': {warning}")
            
            # Use cleaned sequence
            sequence = preprocessing_result.sequence
        else:
            # Original validation (backward compatibility)
            sequence = sequence.upper().strip()
            is_valid, msg = validate_sequence(sequence)
            if not is_valid:
                if "too short" in msg.lower():
                    logger.warning(f"Skipping sequence '{sequence_name}': {msg}")
                    return []
                raise ValueError(f"Invalid sequence: {msg}")
        
        # Auto-enable parallel detectors for large sequences if not specified
        if use_parallel_detectors is None:
            use_parallel_detectors = len(sequence) >= CHUNK_THRESHOLD and USE_PARALLEL_DETECTORS
        
        all_motifs = []
        if enabled_classes:
            enabled_detectors = {CLASS_TO_DETECTOR.get(c) for c in enabled_classes if CLASS_TO_DETECTOR.get(c) in self.detectors}
            detectors_to_run = {k: v for k, v in self.detectors.items() if k in enabled_detectors}
        else: detectors_to_run = self.detectors
        total_detectors = len(detectors_to_run); _reset_detector_timings()
        
        if use_parallel_detectors and total_detectors > 1:
            # Parallel detector execution using ThreadPoolExecutor
            all_motifs = self._analyze_parallel_detectors(sequence, sequence_name, detectors_to_run, progress_callback)
        else:
            # Sequential detector execution (original implementation)
            for idx, (detector_name, detector) in enumerate(detectors_to_run.items()):
                try:
                    start_time = time.time(); motifs = detector.detect_motifs(sequence, sequence_name); elapsed = time.time() - start_time; motif_count = len(motifs)
                    _update_detector_timing(detector_name, elapsed); all_motifs.extend(motifs)
                    if progress_callback is not None: progress_callback(detector_name, idx + 1, total_detectors, elapsed, motif_count)
                except Exception as e:
                    warnings.warn(f"Error in {detector_name} detector: {e}")
                    if progress_callback is not None: progress_callback(detector_name, idx + 1, total_detectors, 0.0, 0)
        
        # Consolidated filtering - do overlap removal once on all motifs
        filtered_motifs = self._remove_overlaps(all_motifs)
        hybrid_motifs = self._detect_hybrid_motifs(filtered_motifs, sequence)
        cluster_motifs = self._detect_clusters(filtered_motifs, sequence)
        # Combine and sort (normalization now handled by detectors)
        # NOTE: normalize_motif_scores() deprecated - detectors self-normalize scores
        final_motifs = filtered_motifs + hybrid_motifs + cluster_motifs
        final_motifs.sort(key=lambda x: x.get('Start', 0))
        return final_motifs
    
    def _analyze_parallel_detectors(self, sequence: str, sequence_name: str, detectors_to_run: Dict, progress_callback: Optional[Callable] = None) -> List[Dict[str, Any]]:
        """Execute detectors in parallel using ThreadPoolExecutor.
        
        Uses thread-based parallelism since detector operations are mostly computational
        and NumPy operations release the GIL, allowing true parallel execution.
        
        Args:
            sequence: DNA sequence to analyze
            sequence_name: Name identifier for the sequence
            detectors_to_run: Dictionary of detector_name -> detector_instance
            progress_callback: Optional callback for progress updates
        
        Returns:
            List of all detected motifs from all detectors
        """
        all_motifs = []
        completed_count = 0  # Must be modified only within results_lock to prevent race conditions
        total_detectors = len(detectors_to_run)
        results_lock = threading.Lock()
        
        def run_detector(detector_name: str, detector):
            """Worker function to run a single detector."""
            nonlocal completed_count  # All modifications of completed_count must occur within results_lock
            try:
                start_time = time.time()
                motifs = detector.detect_motifs(sequence, sequence_name)
                elapsed = time.time() - start_time
                motif_count = len(motifs)
                
                _update_detector_timing(detector_name, elapsed)
                
                with results_lock:
                    completed_count += 1  # Thread-safe increment
                    if progress_callback is not None:
                        progress_callback(detector_name, completed_count, total_detectors, elapsed, motif_count)
                
                return motifs
            except Exception as e:
                warnings.warn(f"Error in {detector_name} detector: {e}")
                with results_lock:
                    completed_count += 1  # Thread-safe increment
                    if progress_callback is not None:
                        progress_callback(detector_name, completed_count, total_detectors, 0.0, 0)
                return []
        
        # Use ThreadPoolExecutor for parallel detector execution
        max_workers = min(MAX_DETECTOR_WORKERS, total_detectors)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all detector jobs
            future_to_detector = {
                executor.submit(run_detector, detector_name, detector): detector_name
                for detector_name, detector in detectors_to_run.items()
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_detector):
                motifs = future.result()
                all_motifs.extend(motifs)
        
        return all_motifs
    
    def _remove_overlaps(self, motifs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """O(n log n) overlap removal using binary search on sorted non-overlapping intervals.

        The overlap check is O(log n) per candidate via bisect; list.insert() shifts
        are O(n) but are C-level memcpy — the same cost as the previous bisect.insort()
        call on the single-tuple list.  The former O(n) Python comparison loop is gone.
        """
        if not motifs: return motifs
        groups = defaultdict(list)
        for motif in motifs: groups[f"{motif.get('Class', '')}-{motif.get('Subclass', '')}"].append(motif)
        filtered_motifs = []
        for group_motifs in groups.values():
            if len(group_motifs) <= 1: filtered_motifs.extend(group_motifs); continue
            # Sort by score (desc), then length (desc) to keep best motifs
            group_motifs.sort(key=lambda x: (-x.get('Score', 0), -x.get('Length', 0)))
            non_overlapping = []
            # Maintain parallel sorted lists of accepted starts/ends for O(log n) overlap check.
            # Accepted intervals are non-overlapping, so both lists are monotonically sorted.
            acc_starts: List[int] = []
            acc_ends: List[int] = []
            for motif in group_motifs:
                start, end = motif.get('Start', 0), motif.get('End', 0)
                # Bisect by start (not end) to find the correct neighbor interval.
                # Use >= / <= for 1-based inclusive coordinates.
                idx = bisect.bisect_left(acc_starts, start)
                overlaps = (idx > 0 and acc_ends[idx - 1] >= start) or \
                           (idx < len(acc_starts) and acc_starts[idx] <= end)
                if not overlaps:
                    non_overlapping.append(motif)
                    ins = bisect.bisect_left(acc_starts, start)
                    acc_starts.insert(ins, start)
                    acc_ends.insert(ins, end)
            filtered_motifs.extend(non_overlapping)
        return filtered_motifs
    
    def _calculate_overlap(self, motif1: Dict[str, Any], motif2: Dict[str, Any]) -> float:
        start1, end1 = motif1.get('Start', 0), motif1.get('End', 0); start2, end2 = motif2.get('Start', 0), motif2.get('End', 0)
        if end1 < start2 or end2 < start1: return 0.0
        overlap_length = min(end1, end2) - max(start1, start2) + 1; min_length = min(end1 - start1 + 1, end2 - start2 + 1)
        return overlap_length / min_length if min_length > 0 else 0.0
    
    def _detect_hybrid_motifs(self, motifs: List[Dict[str, Any]], sequence: str) -> List[Dict[str, Any]]:
        """Optimized hybrid motif detection with cached lookups and early exit."""
        if len(motifs) < 2: return []
        hybrid_motifs = []; seen_hybrids = set(); sorted_motifs = sorted(motifs, key=lambda x: (x.get('Start', 0), x.get('End', 0))); n = len(sorted_motifs)
        # Cache dict lookups for hot path
        motif_data = [(m.get('Start', 0), m.get('End', 0), m.get('Class', ''), m.get('Score', 0), m.get('Sequence_Name', 'seq')) for m in sorted_motifs]
        seq_len = len(sequence)
        for i in range(n):
            start1, end1, class1, score1, seq_name = motif_data[i]
            # Binary search for potential overlaps instead of linear scan
            for j in range(i + 1, n):
                start2, end2, class2, score2, _ = motif_data[j]
                if start2 > end1: break  # No more overlaps possible (1-based inclusive)
                if class1 == class2: continue
                # Inline overlap calculation for 1-based inclusive coordinates
                if end1 < start2 or end2 < start1: continue
                overlap_len = min(end1, end2) - max(start1, start2) + 1; min_len = min(end1 - start1 + 1, end2 - start2 + 1)
                overlap = overlap_len / min_len if min_len > 0 else 0.0
                if HYBRID_MIN_OVERLAP < overlap < HYBRID_MAX_OVERLAP:
                    start, end = min(start1, start2), max(end1, end2); hybrid_key = (start, end, frozenset([class1, class2]))
                    if hybrid_key in seen_hybrids: continue
                    seen_hybrids.add(hybrid_key); avg_score = (score1 + score2) / 2
                    seq_text = sequence[start-1:end] if 0 <= start - 1 < seq_len and 0 < end <= seq_len else 'HYBRID_REGION'
                    hybrid_motifs.append({'ID': f"{seq_name}_HYBRID_{start}", 'Sequence_Name': seq_name, 'Class': 'Hybrid', 'Subclass': f"{class1}_{class2}_Overlap", 'Start': start, 'End': end, 'Length': end - start + 1, 'Sequence': seq_text, 'Score': round(avg_score, 3), 'Strand': '+', 'Method': 'Hybrid_Detection', 'Component_Classes': [class1, class2]})
        return hybrid_motifs
    
    def _detect_clusters(self, motifs: List[Dict[str, Any]], sequence: str) -> List[Dict[str, Any]]:
        """Optimized cluster detection using binary search for O(n log n) complexity."""
        if len(motifs) < 3: return []
        cluster_motifs = []; sorted_motifs = sorted(motifs, key=lambda x: x.get('Start', 0)); n = len(sorted_motifs); detected_window_starts = set()
        # Pre-cache motif positions and create separate start positions list for binary search
        motif_positions = [(m.get('Start', 0), m.get('End', 0), m.get('Class'), m.get('Score', 0), m.get('Sequence_Name', 'seq')) for m in sorted_motifs]
        start_positions = [p[0] for p in motif_positions]  # Created once - O(n)
        seq_len = len(sequence)
        for i in range(n):
            window_start, _, _, _, seq_name = motif_positions[i]; window_end = window_start + CLUSTER_WINDOW_SIZE
            if window_start in detected_window_starts: continue
            # Binary search to find end of window - O(log n)
            end_idx = bisect.bisect_right(start_positions, window_end, lo=i)
            window_motifs_data = motif_positions[i:end_idx]
            if len(window_motifs_data) >= CLUSTER_MIN_MOTIFS:
                classes = set(p[2] for p in window_motifs_data)
                if len(classes) >= CLUSTER_MIN_CLASSES:
                    actual_start = min(p[0] for p in window_motifs_data); actual_end = max(p[1] for p in window_motifs_data)
                    detected_window_starts.add(window_start); avg_score = sum(p[3] for p in window_motifs_data) / len(window_motifs_data)
                    seq_text = sequence[actual_start-1:actual_end] if 0 <= actual_start - 1 < seq_len and 0 < actual_end <= seq_len else 'CLUSTER_REGION'
                    cluster_motifs.append({'ID': f"{seq_name}_CLUSTER_{actual_start}", 'Sequence_Name': seq_name, 'Class': 'Non-B_DNA_Clusters', 'Subclass': f'Mixed_Cluster_{len(classes)}_classes', 'Start': actual_start, 'End': actual_end, 'Length': actual_end - actual_start + 1, 'Sequence': seq_text, 'Score': round(avg_score, 3), 'Strand': '+', 'Method': 'Cluster_Detection', 'Motif_Count': len(window_motifs_data), 'Class_Diversity': len(classes), 'Component_Classes': list(classes)})
        return cluster_motifs
    
    def get_detector_info(self) -> Dict[str, Any]:
        info = {'total_detectors': len(self.detectors), 'detectors': {}, 'total_patterns': 0}
        for name, detector in self.detectors.items(): stats = detector.get_statistics(); info['detectors'][name] = stats; info['total_patterns'] += stats['total_patterns']
        return info

def analyze_sequence(sequence: str, sequence_name: str = "sequence", use_fast_mode: bool = True, use_chunking: bool = None, chunk_size: int = None, chunk_overlap: int = None, progress_callback: Optional[Callable[[int, int, int, float, float], None]] = None, use_parallel_chunks: bool = True, use_parallel_detectors: bool = None, enabled_classes: Optional[List[str]] = None) -> List[Dict[str, Any]]:
    """
    Analyze DNA sequence for non-B DNA motifs with robust error handling.
    
    Args:
        sequence: DNA sequence string to analyze
        sequence_name: Name/identifier for the sequence
        use_fast_mode: Whether to use parallel scanner if available
        use_chunking: Whether to use chunked analysis (auto-enabled for large sequences if None)
        chunk_size: Size of chunks for chunked analysis (default: DEFAULT_CHUNK_SIZE)
        chunk_overlap: Overlap between chunks (default: DEFAULT_CHUNK_OVERLAP)
        progress_callback: Callback function for progress updates (chunk_num, total_chunks, bp_processed, elapsed, throughput)
        use_parallel_chunks: Whether to process chunks in parallel
        use_parallel_detectors: Whether to run detectors in parallel (auto-enabled for sequences >50KB if None)
        enabled_classes: List of motif classes to detect (None = all classes)
    
    Returns:
        List of detected motif dictionaries
        Returns empty list if sequence is None, empty, or invalid
    
    Raises:
        TypeError: If sequence is not a string or None
    """
    # Validate input - check type first before attempting len()
    if sequence is None or not sequence or len(sequence) == 0:
        logger.warning(f"Empty or None sequence provided for {sequence_name}")
        return []
    if not isinstance(sequence, str):
        raise TypeError(f"Sequence must be string, got {type(sequence)}")
    
    seq_len = len(sequence); chunk_size = chunk_size or DEFAULT_CHUNK_SIZE; chunk_overlap = chunk_overlap or DEFAULT_CHUNK_OVERLAP
    if use_chunking is None: use_chunking = seq_len > SEQUENCE_CHUNKING_THRESHOLD
    if not use_chunking or seq_len <= chunk_size:
        if use_fast_mode:
            try: from parallel_scanner import analyze_sequence_parallel; return analyze_sequence_parallel(sequence, sequence_name, use_parallel=True, enabled_classes=enabled_classes)
            except ImportError: warnings.warn("Fast mode not available, falling back to standard mode")
        return _get_cached_scanner().analyze_sequence(sequence, sequence_name, enabled_classes=enabled_classes, use_parallel_detectors=use_parallel_detectors)
    return _analyze_sequence_chunked(sequence, sequence_name, chunk_size, chunk_overlap, progress_callback, use_parallel_chunks, enabled_classes, use_parallel_detectors)

def _process_chunk_worker(chunk_info: Tuple[int, Tuple[int, int]], sequence: str, sequence_name: str, enabled_classes: Optional[List[str]], use_parallel_detectors: bool = True) -> Tuple[int, int, List[Dict[str, Any]]]:
    """
    Worker function for parallel chunk processing.
    Must be at module level to be picklable by ProcessPoolExecutor.
    
    Args:
        chunk_info: Tuple of (chunk_idx, (chunk_start, chunk_end))
        sequence: Full DNA sequence string
        sequence_name: Name/identifier for the sequence
        enabled_classes: List of motif classes to detect
        use_parallel_detectors: Whether to run detectors in parallel within the chunk
    
    Returns:
        Tuple of (chunk_idx, chunk_length, chunk_motifs)
    """
    chunk_idx, (chunk_start, chunk_end) = chunk_info
    chunk_seq = sequence[chunk_start:chunk_end]
    scanner = _get_cached_scanner()
    chunk_motifs = scanner.analyze_sequence(chunk_seq, sequence_name, enabled_classes=enabled_classes, use_parallel_detectors=use_parallel_detectors)
    # Adjust positions relative to full sequence
    for motif in chunk_motifs:
        motif['Start'] += chunk_start
        motif['End'] += chunk_start
    return chunk_idx, chunk_end - chunk_start, chunk_motifs

def _analyze_sequence_chunked(sequence: str, sequence_name: str, chunk_size: int, chunk_overlap: int, progress_callback: Optional[Callable[[int, int, int, float, float], None]] = None, use_parallel_chunks: bool = True, enabled_classes: Optional[List[str]] = None, use_parallel_detectors: bool = None) -> List[Dict[str, Any]]:
    """Optimized chunked analysis with ProcessPoolExecutor for true CPU parallelism."""
    # Validate input - check for None and empty sequences
    if sequence is None or not sequence or len(sequence) == 0:
        logger.warning(f"Empty or None sequence provided for chunked analysis: {sequence_name}")
        return []
    
    def _throughput(bp, elapsed): return bp / elapsed if elapsed > 0 else 0
    seq_len = len(sequence); chunks = []; start = 0
    while start < seq_len:
        end = min(start + chunk_size, seq_len); chunks.append((start, end))
        if end >= seq_len: break
        start = end - chunk_overlap
    total_chunks = len(chunks); all_motifs = []; start_time = time.time(); bp_processed = 0
    if use_parallel_chunks and total_chunks > 1:
        max_workers = min(total_chunks, os.cpu_count() or 4)
        
        try:
            # Try true parallelism with ProcessPoolExecutor
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                futures = {
                    executor.submit(_process_chunk_worker, (i, chunk), sequence, sequence_name, enabled_classes, use_parallel_detectors): i 
                    for i, chunk in enumerate(chunks)
                }
                results_by_idx = {}
                for future in as_completed(futures):
                    chunk_idx, chunk_len, chunk_motifs = future.result()
                    results_by_idx[chunk_idx] = chunk_motifs
                    bp_processed += chunk_len
                    if progress_callback:
                        elapsed = time.time() - start_time
                        progress_callback(chunk_idx + 1, total_chunks, bp_processed, elapsed, _throughput(bp_processed, elapsed))
                for i in range(total_chunks):
                    all_motifs.extend(results_by_idx.get(i, []))
        except (RuntimeError, OSError, AttributeError, BrokenProcessPool) as e:
            # Fallback to sequential if multiprocessing fails (e.g., restricted environments or pickle errors)
            logger.warning(f"ProcessPoolExecutor failed ({e}), falling back to sequential processing")
            scanner = _get_cached_scanner()
            for chunk_idx, (chunk_start, chunk_end) in enumerate(chunks):
                chunk_seq = sequence[chunk_start:chunk_end]
                chunk_motifs = scanner.analyze_sequence(chunk_seq, sequence_name, enabled_classes=enabled_classes, use_parallel_detectors=use_parallel_detectors)
                for motif in chunk_motifs: motif['Start'] += chunk_start; motif['End'] += chunk_start
                all_motifs.extend(chunk_motifs); bp_processed += chunk_end - chunk_start
                if progress_callback: elapsed = time.time() - start_time; progress_callback(chunk_idx + 1, total_chunks, bp_processed, elapsed, _throughput(bp_processed, elapsed))
    else:
        # Sequential processing
        scanner = _get_cached_scanner()
        for chunk_idx, (chunk_start, chunk_end) in enumerate(chunks):
            chunk_seq = sequence[chunk_start:chunk_end]
            chunk_motifs = scanner.analyze_sequence(chunk_seq, sequence_name, enabled_classes=enabled_classes, use_parallel_detectors=use_parallel_detectors)
            for motif in chunk_motifs: motif['Start'] += chunk_start; motif['End'] += chunk_start
            all_motifs.extend(chunk_motifs); bp_processed += chunk_end - chunk_start
            if progress_callback: elapsed = time.time() - start_time; progress_callback(chunk_idx + 1, total_chunks, bp_processed, elapsed, _throughput(bp_processed, elapsed))
    deduplicated_motifs = _deduplicate_motifs(all_motifs); deduplicated_motifs.sort(key=lambda x: x.get('Start', 0))
    return deduplicated_motifs

def _deduplicate_motifs(motifs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Remove duplicate motifs using overlap-based matching for boundary motifs.
    
    Key Changes:
    - Sort by POSITION first (not class)
    - Use 50% overlap threshold to detect duplicates
    - Keep motif with highest score when duplicates found
    
    Scientific Basis:
    - Motifs in chunk overlaps (2KB) can differ by up to overlap size
    - 50% threshold balances precision (no false merges) vs recall (catch all duplicates)
    - Position-first sorting ensures adjacent duplicates are compared
    
    Args:
        motifs: List of all motifs from all chunks (may contain duplicates)
        
    Returns:
        Deduplicated list sorted by genomic position
    """
    if not motifs:
        return motifs
    
    # Sort by POSITION first (critical for boundary detection)
    sorted_motifs = sorted(motifs, key=lambda x: (
        x.get('Start', 0),        # Position FIRST
        x.get('End', 0),
        x.get('Class', '')        # Then class
    ))
    
    deduplicated = []
    
    for motif in sorted_motifs:
        is_duplicate = False
        best_match_idx = -1
        best_match_score = -1
        motif_start = motif.get('Start', 0)
        
        # Check overlap with already-added motifs.
        # Duplicates can only arise within the chunk-overlap window (DEFAULT_CHUNK_OVERLAP).
        # Iterate in reverse so we check the most recent (closest) motifs first and
        # break early once we've moved farther back than the overlap window.
        for i in range(len(deduplicated) - 1, -1, -1):
            existing = deduplicated[i]
            # Early-exit: existing motif starts more than overlap-window before current
            if existing.get('Start', 0) < motif_start - DEFAULT_CHUNK_OVERLAP:
                break
            # Must be same motif type
            if (motif.get('Class') != existing.get('Class') or
                motif.get('Subclass') != existing.get('Subclass')):
                continue
            
            # Calculate overlap using 1-based inclusive coordinates
            start1, end1 = motif_start, motif.get('End', 0)
            start2, end2 = existing.get('Start', 0), existing.get('End', 0)
            
            overlap_start = max(start1, start2)
            overlap_end = min(end1, end2)
            overlap_len = max(0, overlap_end - overlap_start + 1)
            
            # Calculate overlap percentage relative to shorter motif (1-based inclusive)
            len1 = end1 - start1 + 1
            len2 = end2 - start2 + 1
            min_len = min(len1, len2)
            
            if min_len > 0:
                overlap_pct = overlap_len / min_len
                
                # 50% overlap threshold
                if overlap_pct >= 0.5:
                    is_duplicate = True
                    # Track best match (highest score among duplicates)
                    if existing.get('Score', 0) > best_match_score:
                        best_match_idx = i
                        best_match_score = existing.get('Score', 0)
        
        # Handle duplicate: replace with higher scoring motif if needed
        if is_duplicate and best_match_idx >= 0:
            if motif.get('Score', 0) > best_match_score:
                deduplicated[best_match_idx] = motif
        elif not is_duplicate:
            deduplicated.append(motif)
    
    # Final sort by position
    deduplicated.sort(key=lambda x: x.get('Start', 0))
    
    return deduplicated

def analyze_fasta(fasta_content: str) -> Dict[str, List[Dict[str, Any]]]:
    sequences = parse_fasta(fasta_content); scanner = _get_cached_scanner()
    return {name: scanner.analyze_sequence(seq, name) for name, seq in sequences.items()}

def analyze_file(filename: str) -> Dict[str, List[Dict[str, Any]]]:
    if not os.path.exists(filename): raise FileNotFoundError(f"File not found: {filename}")
    sequences = read_fasta_file(filename); scanner = _get_cached_scanner()
    return {name: scanner.analyze_sequence(seq, name) for name, seq in sequences.items()}

def _analyze_sequence_worker(args: Tuple[str, str]) -> Tuple[str, List[Dict[str, Any]]]:
    name, seq = args; import nonbscanner; return (name, nonbscanner.analyze_sequence(seq, name))

def analyze_multiple_sequences_parallel(sequences: Dict[str, str], num_processes: Optional[int] = None, preserve_order: bool = True) -> Dict[str, List[Dict[str, Any]]]:
    if not sequences: return {}
    # Handle single sequence case
    if len(sequences) == 1:
        name, seq = next(iter(sequences.items()))
        return {name: _get_cached_scanner().analyze_sequence(seq, name)}
    num_processes = num_processes or min(multiprocessing.cpu_count(), len(sequences))
    try:
        sequence_items = list(sequences.items())
        with multiprocessing.Pool(processes=num_processes) as pool: results_list = pool.map(_analyze_sequence_worker, sequence_items)
        results = dict(results_list)
        if preserve_order: results = {name: results[name] for name, _ in sequence_items}
        logger.info(f"Parallel analysis completed: {len(sequences)} sequences using {num_processes} processes"); return results
    except Exception as e:
        logger.warning(f"Parallel processing failed, falling back to sequential: {e}"); scanner = _get_cached_scanner()
        return {name: scanner.analyze_sequence(seq, name) for name, seq in sequences.items()}

def analyze_fasta_parallel(fasta_content: str, num_processes: Optional[int] = None, preserve_order: bool = True) -> Dict[str, List[Dict[str, Any]]]:
    return analyze_multiple_sequences_parallel(parse_fasta(fasta_content), num_processes, preserve_order)

def analyze_file_parallel(filename: str, num_processes: Optional[int] = None, preserve_order: bool = True) -> Dict[str, List[Dict[str, Any]]]:
    if not os.path.exists(filename): raise FileNotFoundError(f"File not found: {filename}")
    return analyze_multiple_sequences_parallel(read_fasta_file(filename), num_processes, preserve_order)

def get_motif_info() -> Dict[str, Any]:
    from Utilities.config.motif_taxonomy import MOTIF_CLASSIFICATION, VALID_SUBCLASSES
    classification_legacy = {motif_id: {'name': motif_data['class'], 'subclasses': motif_data['subclasses']} for motif_id, motif_data in MOTIF_CLASSIFICATION.items()}
    return {'version': __version__, 'author': __author__, 'total_classes': 11, 'total_subclasses': len(VALID_SUBCLASSES), 'classification': classification_legacy}

def export_results(motifs: List[Dict[str, Any]], format: str = 'csv', filename: Optional[str] = None, **kwargs) -> str:
    if format.lower() == 'csv': return export_to_csv(motifs, filename)
    elif format.lower() == 'bed': return export_to_bed(motifs, kwargs.get('sequence_name', 'sequence'), filename)
    elif format.lower() == 'json': return export_to_json(motifs, filename, kwargs.get('pretty', True))
    elif format.lower() in ['excel', 'xlsx']: return export_to_excel(motifs, filename or 'nonbscanner_results.xlsx')
    else: raise ValueError(f"Unsupported format: {format}. Use 'csv', 'bed', 'json', or 'excel'")

@overload
def analyze_with_progress(sequence: str, sequence_name: str = ..., print_progress: bool = ..., return_progress: Literal[False] = ...) -> List[Dict[str, Any]]: ...
@overload
def analyze_with_progress(sequence: str, sequence_name: str = ..., print_progress: bool = ..., return_progress: Literal[True] = ...) -> Tuple[List[Dict[str, Any]], AnalysisProgress]: ...

def analyze_with_progress(sequence: str, sequence_name: str = "sequence", print_progress: bool = True, return_progress: bool = False) -> Union[List[Dict[str, Any]], Tuple[List[Dict[str, Any]], AnalysisProgress]]:
    progress = AnalysisProgress(sequence_length=len(sequence), sequence_name=sequence_name)
    if print_progress: print(f"\n{'='*70}\nNonBScanner - Non-B DNA Motif Detection\n{'='*70}\nSequence: {sequence_name} ({len(sequence):,} bp)\n\nAnalysis Progress:")
    callback = create_enhanced_progress_callback(progress, print_updates=print_progress); progress.start_stage('DETECTION')
    motifs = _get_cached_scanner().analyze_sequence(sequence, sequence_name, progress_callback=callback); progress.end_stage('DETECTION')
    if print_progress: print(f"\n{progress.format_detector_table()}\n\nSummary:\n  Total motifs: {len(motifs)}\n  Elapsed time: {progress.get_elapsed_time():.2f}s\n  Throughput: {progress.get_throughput():,.0f} bp/s\n{'='*70}\n")
    return (motifs, progress) if return_progress else motifs

def get_pipeline_info() -> Dict[str, Any]:
    return {'stages': [{'id': 'INPUT', 'name': 'Input Parsing', 'icon': '📥', 'description': 'Parse FASTA format and validate DNA sequences', 'complexity': 'O(n)', 'operations': ['Parse FASTA headers', 'Extract sequences', 'Validate ATGC characters']}, {'id': 'DETECTION', 'name': 'Motif Detection', 'icon': '[DETECT]', 'description': 'Run 9 specialized detectors to find Non-B DNA motifs', 'complexity': 'O(n) per detector', 'operations': ['Run parallel detectors', 'Pattern matching', 'Score calculation']}, {'id': 'PROCESSING', 'name': 'Post-Processing', 'icon': '[PROCESS]', 'description': 'Process and refine detected motifs', 'complexity': 'O(m log m) to O(m²)', 'operations': ['Remove overlapping motifs', 'Detect hybrid regions', 'Identify clusters']}, {'id': 'OUTPUT', 'name': 'Output Generation', 'icon': '[OUTPUT]', 'description': 'Normalize scores and format results', 'complexity': 'O(m)', 'operations': ['Normalize scores to 1-3 scale', 'Sort by position', 'Generate output']}], 'detectors': [{'name': 'Curved DNA', 'method': 'A-tract phasing', 'subclasses': 2}, {'name': 'Slipped DNA', 'method': 'K-mer indexing', 'subclasses': 2}, {'name': 'Cruciform', 'method': 'Inverted repeat detection', 'subclasses': 1}, {'name': 'R-Loop', 'method': 'QmRLFS algorithm', 'subclasses': 3}, {'name': 'Triplex', 'method': 'Mirror repeat + purine runs', 'subclasses': 2}, {'name': 'G-Quadruplex', 'method': 'G4Hunter + patterns', 'subclasses': 7}, {'name': 'i-Motif', 'method': 'C-run patterns', 'subclasses': 3}, {'name': 'Z-DNA', 'method': 'CG/CA repeat scoring', 'subclasses': 2}, {'name': 'A-philic DNA', 'method': '10-mer patterns', 'subclasses': 1}], 'total_detectors': 9, 'total_subclasses': 23}

def analyze_multiple_sequences(sequences: Dict[str, str], use_multiprocessing: bool = False) -> Dict[str, List[Dict[str, Any]]]:
    scanner = NonBScanner()
    if use_multiprocessing:
        from concurrent.futures import ProcessPoolExecutor, as_completed; import multiprocessing as mp
        with ProcessPoolExecutor(max_workers=min(len(sequences), mp.cpu_count())) as executor:
            future_to_name = {executor.submit(scanner.analyze_sequence, seq, name): name for name, seq in sequences.items()}
            return {future_to_name[f]: f.result() if not f.exception() else [] for f in as_completed(future_to_name)}
    return {name: scanner.analyze_sequence(seq, name) for name, seq in sequences.items()}

def get_summary_statistics(results: Dict[str, List[Dict[str, Any]]]) -> pd.DataFrame:
    summary_data = []
    for name, motifs in results.items():
        stats = calculate_motif_statistics(motifs, 0)
        summary_data.append({'Sequence_Name': name, 'Total_Motifs': stats.get('Total_Motifs', 0), 'Classes_Detected': stats.get('Classes_Detected', 0), 'Subclasses_Detected': stats.get('Subclasses_Detected', 0), 'Score_Mean': stats.get('Score_Mean', 0), 'Length_Mean': stats.get('Length_Mean', 0)})
    return pd.DataFrame(summary_data)

def main():
    print(f"{'='*70}\nNonBScanner - Non-B DNA Motif Detection Suite\nVersion {__version__} by {__author__}\n{'='*70}")
    info = get_motif_info(); print(f"\nDetection Capabilities:\n  Total Classes: {info['total_classes']}\n  Total Subclasses: {info['total_subclasses']}")
    demo_seq = "GGGTTAGGGTTAGGGTTAGGGAAAAATTTTCGCGCGCGCGATATATATATCCCCTAACCCTAACCCTAACCC"
    print(f"\nDemo Analysis:\n  Sequence: {demo_seq[:50]}... ({len(demo_seq)} bp)")
    motifs = analyze_sequence(demo_seq, "demo"); print(f"\nResults: {len(motifs)} motifs detected")
    for motif in motifs[:5]: print(f"  {motif['Class']:<15} {motif['Subclass']:<20} Pos:{motif['Start']:>3}-{motif['End']:<3} Score:{motif['Score']:.2f}")
    if len(motifs) > 5: print(f"  ... and {len(motifs) - 5} more")
    print(f"\nNonBScanner ready for analysis!\n{'='*70}")

if __name__ == "__main__": main()

# ═══════════════════════════════════════════════════════════════════════════════
# STREAMLIT CLOUD OPTIMIZATION
# ═══════════════════════════════════════════════════════════════════════════════

# Check if running on Streamlit Cloud or optimization explicitly requested
USE_OPTIMIZED = os.environ.get('NONBDNA_OPTIMIZED', 'true').lower() == 'true'

if USE_OPTIMIZED:
    try:
        from Utilities.nonbscanner_optimized import NonBScannerOptimized
        
        # Replace cached scanner with optimized version
        def _get_cached_scanner_optimized() -> 'NonBScannerOptimized':
            """Get cached optimized scanner instance."""
            global _CACHED_SCANNER
            if _CACHED_SCANNER is None:
                with _SCANNER_LOCK:
                    if _CACHED_SCANNER is None:
                        _CACHED_SCANNER = NonBScannerOptimized(enable_all_detectors=True)
                        logger.info("✓ Using OPTIMIZED NonBScanner (50-200x faster)")
            return _CACHED_SCANNER
        
        # Replace the standard cached scanner function with optimized version
        _get_cached_scanner = _get_cached_scanner_optimized
        
        logger.info("NonBDNAFinder: Aho-Corasick optimization ENABLED")
    
    except ImportError as e:
        logger.warning(f"Optimization not available: {e}. Using standard implementation.")
else:
    logger.info("NonBDNAFinder: Using standard implementation (set NONBDNA_OPTIMIZED=true to enable)")
