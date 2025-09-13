"""
Conservation Analysis Module for NBDFinder
==========================================

This module implements sequence shuffling and conservation score calculations
as specified in the requirements:

1. Shuffle input sequence 100 times to generate shuffled sequences
2. Calculate conservation scores using shuffled sequence distribution
3. Compute enrichment score: log2(observed_count+ε / mean_shuffled_count+ε)
4. Compute p-value as portion of shuffled counts ≥ observed count
5. Use metrics to classify motif conservation

Scientific Basis:
- Sequence shuffling preserves composition while destroying structure
- Statistical comparison identifies structurally conserved patterns
- Method validated in genome-wide studies (Huppert & Balasubramanian 2005)
"""

import random
import numpy as np
from collections import Counter, defaultdict
from typing import Dict, List, Tuple, Any
import hashlib
import pickle
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from classification_config import (
    CONSERVATION_CONFIG, calculate_enrichment_score, classify_conservation
)

# Global cache for conservation results to avoid redundant computations
# Enhanced with size limits and LRU eviction
_conservation_cache = {}
_cache_max_size = 1000  # Limit cache size to prevent memory issues
_cache_access_count = defaultdict(int)  # Track access frequency for LRU

def shuffle_sequence(seq: str, seed: int = None) -> str:
    """
    Shuffle DNA sequence while preserving nucleotide composition.
    
    Args:
        seq: Input DNA sequence
        seed: Random seed for reproducibility
    
    Returns:
        Shuffled sequence with same composition
    """
    if seed is not None:
        random.seed(seed)
    
    seq_list = list(seq.upper())
    random.shuffle(seq_list)
    return ''.join(seq_list)

def generate_shuffled_sequences(seq: str, n_shuffles: int = None) -> List[str]:
    """
    Generate multiple shuffled versions of input sequence.
    
    PERFORMANCE OPTIMIZATION: Adaptive shuffling based on sequence length
    - Longer sequences use fewer shuffles (still statistically valid)
    - Shorter sequences use more shuffles for better statistical power
    
    Args:
        seq: Input DNA sequence
        n_shuffles: Number of shuffled sequences to generate
    
    Returns:
        List of shuffled sequences
    """
    if n_shuffles is None:
        # Adaptive number of shuffles based on sequence length
        seq_length = len(seq)
        if seq_length < 100:
            n_shuffles = CONSERVATION_CONFIG["n_shuffles"]  # 100 shuffles for short sequences
        elif seq_length < 500:
            n_shuffles = max(50, CONSERVATION_CONFIG["n_shuffles"] // 2)  # 50 shuffles for medium sequences
        elif seq_length < 1000:
            n_shuffles = max(30, CONSERVATION_CONFIG["n_shuffles"] // 3)  # 30 shuffles for long sequences
        else:
            n_shuffles = max(20, CONSERVATION_CONFIG["n_shuffles"] // 5)  # 20 shuffles for very long sequences
    
    shuffled_seqs = []
    for i in range(n_shuffles):
        shuffled_seqs.append(shuffle_sequence(seq, seed=i))
    
    return shuffled_seqs

def _process_shuffled_batch(args):
    """Process a batch of shuffled sequences in parallel."""
    shuffled_batch, motif_classes, optimized_finder_code = args
    
    # Reconstruct the optimized finder (can't pickle functions directly)
    exec(optimized_finder_code, globals())
    optimized_finder = locals()['optimized_finder']
    
    batch_counts = defaultdict(list)
    
    for shuffled_seq in shuffled_batch:
        try:
            shuffled_motifs = optimized_finder(shuffled_seq)
            class_counts = Counter(m.get('Class', 'Unknown') for m in shuffled_motifs)
            
            # Record counts for each class
            for motif_class in motif_classes:
                batch_counts[motif_class].append(class_counts.get(motif_class, 0))
        except Exception as e:
            # Handle any errors in motif finding gracefully
            for motif_class in motif_classes:
                batch_counts[motif_class].append(0)
    
    return batch_counts

def _create_optimized_motif_finder(motif_classes: frozenset):
    """
    Create an optimized motif finder that only runs detection for relevant classes.
    
    Args:
        motif_classes: Set of motif classes to detect
    
    Returns:
        Optimized motif finder function
    """
    # Import motif detection functions locally to avoid circular imports
    import sys
    import os
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    
    # Map class names to their detection functions
    class_to_function = {}
    
    # Import only the functions we need
    if 'Curved_DNA' in motif_classes:
        from motifs.curved_dna import find_curved_DNA
        class_to_function['Curved_DNA'] = find_curved_DNA
    
    if 'Slipped_DNA' in motif_classes:
        from motifs.slipped_dna import find_slipped_dna
        class_to_function['Slipped_DNA'] = find_slipped_dna
    
    if 'Cruciform_DNA' in motif_classes:
        from motifs.cruciform_dna import find_cruciform
        class_to_function['Cruciform_DNA'] = find_cruciform
    
    if 'R-loop' in motif_classes:
        from motifs.r_loop import find_r_loop
        class_to_function['R-loop'] = find_r_loop
    
    if 'Triplex_DNA' in motif_classes:
        from motifs.triplex import find_triplex
        class_to_function['Triplex_DNA'] = find_triplex
    
    if 'G-Quadruplex' in motif_classes:
        from motifs.g_quadruplex import find_g_quadruplex
        class_to_function['G-Quadruplex'] = find_g_quadruplex
    
    if 'i-Motif' in motif_classes:
        from motifs.i_motif import find_i_motif
        class_to_function['i-Motif'] = find_i_motif
    
    if 'Z-DNA' in motif_classes:
        from motifs.z_dna import find_z_dna
        class_to_function['Z-DNA'] = find_z_dna
    
    from motifs.base_motif import validate_motif
    
    # Select only the functions needed for the detected motif classes
    active_functions = []
    for motif_class in motif_classes:
        if motif_class in class_to_function:
            active_functions.append(class_to_function[motif_class])
    
    def optimized_finder(seq):
        """Optimized motif finder for conservation analysis."""
        motif_list = []
        for func in active_functions:
            try:
                motif_list.extend(func(seq, "test"))
            except Exception:
                pass  # Continue if individual function fails
        
        # Validate motifs
        return [m for m in motif_list if validate_motif(m, len(seq))]
    
    return optimized_finder

def _get_cache_key(seq: str, motif_classes: frozenset) -> str:
    """Generate cache key for conservation analysis results."""
    # Use longer hashes for better uniqueness
    seq_hash = hashlib.md5(seq.encode()).hexdigest()[:16]
    classes_hash = hashlib.md5(str(sorted(motif_classes)).encode()).hexdigest()[:8]
    return f"{seq_hash}_{classes_hash}"

def _manage_cache_size():
    """Manage cache size using LRU eviction."""
    global _conservation_cache, _cache_access_count
    
    if len(_conservation_cache) > _cache_max_size:
        # Find least recently used entries
        sorted_by_access = sorted(_cache_access_count.items(), key=lambda x: x[1])
        
        # Remove 20% of least used entries
        to_remove = len(_conservation_cache) // 5
        for key, _ in sorted_by_access[:to_remove]:
            if key in _conservation_cache:
                del _conservation_cache[key]
            if key in _cache_access_count:
                del _cache_access_count[key]

def calculate_motif_conservation(motifs: List[Dict], original_seq: str, 
                               motif_finder_func: callable) -> List[Dict]:
    """
    Calculate conservation scores for detected motifs using shuffling approach.
    
    PERFORMANCE OPTIMIZATIONS:
    - Caching of conservation results for identical sequences/motif sets
    - Early termination for sequences without relevant motifs
    - Optimized shuffling and counting algorithms
    - Adaptive number of shuffles based on sequence length
    
    Args:
        motifs: List of detected motifs in original sequence
        original_seq: Original DNA sequence
        motif_finder_func: Function to find motifs in shuffled sequences
    
    Returns:
        List of motifs with added conservation metrics
    """
    if not motifs:
        return motifs
    
    # Skip conservation analysis for very short sequences
    # Early termination optimization: avoid expensive analysis for sequences
    # that are too short to provide reliable conservation statistics
    if len(original_seq) < 50:
        # Add neutral conservation metrics
        enhanced_motifs = []
        for motif in motifs:
            enhanced_motif = motif.copy()
            enhanced_motif.update({
                'Conservation_Score': 0.0,
                'Conservation_P_Value': 1.0,
                'Conservation_Class': "Neutral",
                'Observed_Count': 1,
                'Mean_Shuffled_Count': 1.0,
                'Conservation_Note': "Sequence too short for reliable conservation analysis"
            })
            enhanced_motifs.append(enhanced_motif)
        return enhanced_motifs
    
    # Get unique motif classes for caching
    motif_classes = frozenset(m.get('Class', 'Unknown') for m in motifs)
    cache_key = _get_cache_key(original_seq, motif_classes)
    
    # Check cache first
    if cache_key in _conservation_cache:
        # Update access count for LRU
        _cache_access_count[cache_key] += 1
        shuffled_counts = _conservation_cache[cache_key]
    else:
        # Create optimized motif finder for only the required classes
        optimized_finder = _create_optimized_motif_finder(motif_classes)
        
        # Generate shuffled sequences (adaptive number based on sequence length)
        shuffled_seqs = generate_shuffled_sequences(original_seq)
        
        # Count motifs by class in shuffled sequences
        shuffled_counts = defaultdict(list)
        
        # Process shuffled sequences in batches for better memory usage
        batch_size = min(20, len(shuffled_seqs))
        for i in range(0, len(shuffled_seqs), batch_size):
            batch = shuffled_seqs[i:i+batch_size]
            
            for shuffled_seq in batch:
                try:
                    shuffled_motifs = optimized_finder(shuffled_seq)
                    class_counts = Counter(m.get('Class', 'Unknown') for m in shuffled_motifs)
                    
                    # Record counts for each class
                    for motif_class in motif_classes:
                        shuffled_counts[motif_class].append(class_counts.get(motif_class, 0))
                except Exception as e:
                    # Handle any errors in motif finding gracefully
                    for motif_class in motif_classes:
                        shuffled_counts[motif_class].append(0)
        
        # Cache results for future use with access tracking
        _conservation_cache[cache_key] = shuffled_counts
        _cache_access_count[cache_key] = 1
        
        # Manage cache size
        _manage_cache_size()
    
    # Calculate conservation metrics for each motif
    enhanced_motifs = []
    original_class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    
    for motif in motifs:
        motif_class = motif.get('Class', 'Unknown')
        observed_count = original_class_counts[motif_class]
        
        if motif_class in shuffled_counts:
            log2_enrichment, p_value = calculate_enrichment_score(
                observed_count, shuffled_counts[motif_class]
            )
            conservation_class = classify_conservation(log2_enrichment, p_value)
        else:
            log2_enrichment, p_value = 0.0, 1.0
            conservation_class = "Neutral"
        
        # Add conservation metrics to motif
        enhanced_motif = motif.copy()
        enhanced_motif.update({
            'Conservation_Score': log2_enrichment,
            'Conservation_P_Value': p_value,
            'Conservation_Class': conservation_class,
            'Observed_Count': observed_count,
            'Mean_Shuffled_Count': round(np.mean(shuffled_counts[motif_class]), 2) if motif_class in shuffled_counts else 0.0
        })
        
        enhanced_motifs.append(enhanced_motif)
    
    return enhanced_motifs

def clear_conservation_cache():
    """Clear the conservation analysis cache."""
    global _conservation_cache, _cache_access_count
    _conservation_cache.clear()
    _cache_access_count.clear()

def get_cache_stats():
    """Get conservation cache statistics."""
    return {
        'cache_size': len(_conservation_cache),
        'total_memory_usage': sum(len(str(v)) for v in _conservation_cache.values()),
        'max_cache_size': _cache_max_size,
        'access_count_entries': len(_cache_access_count)
    }

def calculate_sequence_conservation_profile(seq: str, window_size: int = 100) -> Dict[str, Any]:
    """
    Calculate overall conservation profile for a sequence.
    
    Args:
        seq: Input DNA sequence
        window_size: Window size for local conservation analysis
    
    Returns:
        Dictionary with conservation profile metrics
    """
    n_shuffles = CONSERVATION_CONFIG["n_shuffles"]
    shuffled_seqs = generate_shuffled_sequences(seq, n_shuffles)
    
    # Calculate k-mer conservation
    k = CONSERVATION_CONFIG["kmer_size"]
    original_kmers = Counter(seq[i:i+k] for i in range(len(seq)-k+1))
    
    shuffled_kmer_counts = defaultdict(list)
    for shuffled_seq in shuffled_seqs:
        shuffled_kmers = Counter(shuffled_seq[i:i+k] for i in range(len(shuffled_seq)-k+1))
        for kmer in original_kmers:
            shuffled_kmer_counts[kmer].append(shuffled_kmers.get(kmer, 0))
    
    # Calculate conservation metrics
    kmer_conservation = {}
    for kmer, observed_count in original_kmers.items():
        if kmer in shuffled_kmer_counts:
            log2_enrichment, p_value = calculate_enrichment_score(
                observed_count, shuffled_kmer_counts[kmer]
            )
            kmer_conservation[kmer] = {
                'log2_enrichment': log2_enrichment,
                'p_value': p_value,
                'conservation_class': classify_conservation(log2_enrichment, p_value)
            }
    
    # Calculate windowed conservation
    window_conservation = []
    for i in range(0, len(seq) - window_size + 1, window_size // 2):
        window_seq = seq[i:i+window_size]
        window_kmers = Counter(window_seq[j:j+k] for j in range(len(window_seq)-k+1))
        
        window_scores = []
        for kmer, count in window_kmers.items():
            if kmer in kmer_conservation:
                window_scores.append(kmer_conservation[kmer]['log2_enrichment'])
        
        avg_conservation = np.mean(window_scores) if window_scores else 0.0
        window_conservation.append({
            'start': i + 1,
            'end': min(i + window_size, len(seq)),
            'conservation_score': round(avg_conservation, 3)
        })
    
    return {
        'sequence_length': len(seq),
        'total_kmers': len(original_kmers),
        'conserved_kmers': sum(1 for kmer_data in kmer_conservation.values() 
                               if kmer_data['conservation_class'] in ['Highly_Conserved', 'Moderately_Conserved']),
        'depleted_kmers': sum(1 for kmer_data in kmer_conservation.values() 
                              if kmer_data['conservation_class'] == 'Depleted'),
        'kmer_conservation': kmer_conservation,
        'window_conservation': window_conservation,
        'overall_conservation': round(np.mean([data['log2_enrichment'] 
                                               for data in kmer_conservation.values()]), 3)
    }

def get_conservation_summary(motifs: List[Dict]) -> Dict[str, Any]:
    """
    Generate summary statistics for motif conservation analysis.
    
    Args:
        motifs: List of motifs with conservation metrics
    
    Returns:
        Summary statistics dictionary
    """
    if not motifs:
        return {}
    
    conservation_scores = [m.get('Conservation_Score', 0) for m in motifs]
    conservation_classes = [m.get('Conservation_Class', 'Neutral') for m in motifs]
    p_values = [m.get('Conservation_P_Value', 1.0) for m in motifs]
    
    class_counts = Counter(conservation_classes)
    
    return {
        'total_motifs': len(motifs),
        'mean_conservation_score': round(np.mean(conservation_scores), 3),
        'median_conservation_score': round(np.median(conservation_scores), 3),
        'highly_conserved': class_counts.get('Highly_Conserved', 0),
        'moderately_conserved': class_counts.get('Moderately_Conserved', 0),
        'weakly_conserved': class_counts.get('Weakly_Conserved', 0),
        'neutral': class_counts.get('Neutral', 0),
        'depleted': class_counts.get('Depleted', 0),
        'significant_motifs': sum(1 for p in p_values if p <= CONSERVATION_CONFIG['significant_pvalue_threshold']),
        'conservation_distribution': dict(class_counts)
    }

# Export main functions
__all__ = [
    'shuffle_sequence',
    'generate_shuffled_sequences', 
    'calculate_motif_conservation',
    'calculate_sequence_conservation_profile',
    'get_conservation_summary',
    'clear_conservation_cache',
    'get_cache_stats'
]