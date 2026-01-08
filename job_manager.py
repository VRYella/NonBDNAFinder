"""
Job Management Module for NBDScanner
=====================================

Provides job ID generation, result persistence, and job lookup functionality.

Features:
- Short, URL-safe unique job IDs (10-character hex)
- Disk-based result persistence under results/<job_id>/
- Job metadata storage (timestamp, sequence info)
- Job lookup and retrieval
"""

import os
import json
import uuid
from datetime import datetime
from typing import Dict, List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)

# Base directory for all job results
RESULTS_BASE_DIR = "results"


def generate_job_id() -> str:
    """
    Generate a unique, URL-safe job ID.
    
    Uses the first 10 characters of a UUID hex string, providing
    2^40 (≈1 trillion) possible values for collision resistance.
    
    Returns:
        str: 10-character hexadecimal job ID
    
    Example:
        >>> job_id = generate_job_id()
        >>> len(job_id)
        10
        >>> all(c in '0123456789abcdef' for c in job_id)
        True
    """
    return uuid.uuid4().hex[:10]


def get_job_directory(job_id: str) -> str:
    """
    Get the file system path for a job's results directory.
    
    Args:
        job_id: The job identifier (must be 10-char hex string)
        
    Returns:
        str: Absolute path to job directory
        
    Raises:
        ValueError: If job_id format is invalid
    """
    # Validate job_id to prevent path traversal attacks
    if not job_id or len(job_id) != 10:
        raise ValueError(f"Invalid job_id length: expected 10, got {len(job_id) if job_id else 0}")
    
    if not all(c in '0123456789abcdef' for c in job_id):
        raise ValueError(f"Invalid job_id format: must contain only hexadecimal characters")
    
    return os.path.join(RESULTS_BASE_DIR, job_id)


def ensure_job_directory(job_id: str) -> str:
    """
    Create job directory if it doesn't exist.
    
    Args:
        job_id: The job identifier
        
    Returns:
        str: Path to job directory
    """
    job_dir = get_job_directory(job_id)
    os.makedirs(job_dir, exist_ok=True)
    return job_dir


def save_job_results(
    job_id: str,
    results: List[List[Dict]],
    sequences: List[str],
    sequence_names: List[str],
    metadata: Optional[Dict] = None
) -> bool:
    """
    Save analysis results to disk under the job ID.
    
    Creates the following files in results/<job_id>/:
    - results.json: All motif detection results
    - metadata.json: Job metadata (timestamp, sequence info)
    - sequences.json: Original sequences (for re-analysis if needed)
    
    Args:
        job_id: The job identifier
        results: List of motif lists (one per sequence)
        sequences: List of DNA sequences analyzed
        sequence_names: List of sequence names
        metadata: Optional additional metadata (analysis params, etc.)
        
    Returns:
        bool: True if save succeeded, False otherwise
    """
    try:
        job_dir = ensure_job_directory(job_id)
        
        # Save results
        results_path = os.path.join(job_dir, "results.json")
        with open(results_path, 'w') as f:
            json.dump(results, f, indent=2)
        
        # Save sequences
        sequences_path = os.path.join(job_dir, "sequences.json")
        sequences_data = {
            "sequences": sequences,
            "names": sequence_names
        }
        with open(sequences_path, 'w') as f:
            json.dump(sequences_data, f, indent=2)
        
        # Save metadata
        metadata_path = os.path.join(job_dir, "metadata.json")
        job_metadata = {
            "job_id": job_id,
            "timestamp": datetime.now().isoformat(),
            "num_sequences": len(sequences),
            "sequence_names": sequence_names,
            "total_bp": sum(len(seq) for seq in sequences),
            "total_motifs": sum(len(motif_list) for motif_list in results)
        }
        
        # Add any additional metadata provided (user metadata in separate namespace)
        if metadata:
            # Prevent overwriting critical system fields
            reserved_keys = {'job_id', 'timestamp', 'num_sequences', 'sequence_names', 'total_bp', 'total_motifs'}
            for key, value in metadata.items():
                if key not in reserved_keys:
                    job_metadata[key] = value
                else:
                    logger.warning(f"Ignoring metadata key '{key}' - reserved for system use")
        
        with open(metadata_path, 'w') as f:
            json.dump(job_metadata, f, indent=2)
        
        logger.info(f"Job {job_id} saved successfully to {job_dir}")
        return True
        
    except Exception as e:
        logger.error(f"Failed to save job {job_id}: {e}")
        return False


def load_job_results(job_id: str) -> Optional[Tuple[List[List[Dict]], List[str], List[str], Dict]]:
    """
    Load analysis results from disk for a given job ID.
    
    Args:
        job_id: The job identifier
        
    Returns:
        Tuple of (results, sequences, sequence_names, metadata) if job exists,
        None if job not found or error occurred
    """
    try:
        job_dir = get_job_directory(job_id)
        
        if not os.path.exists(job_dir):
            logger.warning(f"Job {job_id} not found at {job_dir}")
            return None
        
        # Load results
        results_path = os.path.join(job_dir, "results.json")
        if not os.path.exists(results_path):
            logger.error(f"Results file missing for job {job_id}")
            return None
        
        with open(results_path, 'r') as f:
            results = json.load(f)
        
        # Load sequences
        sequences_path = os.path.join(job_dir, "sequences.json")
        if not os.path.exists(sequences_path):
            logger.error(f"Sequences file missing for job {job_id}")
            return None
        
        with open(sequences_path, 'r') as f:
            sequences_data = json.load(f)
        
        sequences = sequences_data.get("sequences", [])
        sequence_names = sequences_data.get("names", [])
        
        # Load metadata
        metadata_path = os.path.join(job_dir, "metadata.json")
        metadata = {}
        if os.path.exists(metadata_path):
            with open(metadata_path, 'r') as f:
                metadata = json.load(f)
        
        logger.info(f"Job {job_id} loaded successfully from {job_dir}")
        return results, sequences, sequence_names, metadata
        
    except Exception as e:
        logger.error(f"Failed to load job {job_id}: {e}")
        return None


def job_exists(job_id: str) -> bool:
    """
    Check if a job exists on disk.
    
    Args:
        job_id: The job identifier
        
    Returns:
        bool: True if job directory and results exist
    """
    job_dir = get_job_directory(job_id)
    results_path = os.path.join(job_dir, "results.json")
    return os.path.exists(job_dir) and os.path.exists(results_path)


def list_all_jobs() -> List[Dict]:
    """
    List all available jobs with their metadata.
    
    Returns:
        List of job metadata dictionaries, sorted by timestamp (newest first)
    """
    jobs = []
    
    if not os.path.exists(RESULTS_BASE_DIR):
        return jobs
    
    try:
        for job_id in os.listdir(RESULTS_BASE_DIR):
            job_dir = os.path.join(RESULTS_BASE_DIR, job_id)
            
            if not os.path.isdir(job_dir):
                continue
            
            metadata_path = os.path.join(job_dir, "metadata.json")
            if os.path.exists(metadata_path):
                try:
                    with open(metadata_path, 'r') as f:
                        metadata = json.load(f)
                        jobs.append(metadata)
                except Exception as e:
                    logger.warning(f"Failed to load metadata for job {job_id}: {e}")
        
        # Sort by timestamp, newest first
        jobs.sort(key=lambda x: x.get('timestamp', ''), reverse=True)
        
    except Exception as e:
        logger.error(f"Failed to list jobs: {e}")
    
    return jobs


def get_job_summary(job_id: str) -> Optional[Dict]:
    """
    Get a summary of a job without loading full results.
    
    Args:
        job_id: The job identifier
        
    Returns:
        Dict with job metadata, or None if job not found
    """
    try:
        job_dir = get_job_directory(job_id)
        metadata_path = os.path.join(job_dir, "metadata.json")
        
        if not os.path.exists(metadata_path):
            return None
        
        with open(metadata_path, 'r') as f:
            return json.load(f)
            
    except Exception as e:
        logger.error(f"Failed to get summary for job {job_id}: {e}")
        return None
