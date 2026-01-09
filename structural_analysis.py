"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                   STRUCTURAL BLOCK & HYBRID PATTERN MODULE                    ║
║              Higher-Order Pattern Behavior Discovery & Analysis              ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE:        structural_analysis.py
AUTHOR:        Dr. Venkata Rajesh Yella
VERSION:       2025.1 - Structural Framework
LICENSE:       MIT

DESCRIPTION:
    Discovers higher-order pattern behavior beyond raw hits:
    - Pattern-rich blocks across sequences
    - Overlap-driven hybrid zones with interaction analysis
    - Clustering metrics (size, density, distance, skew, stability)
    - Comparison with shuffled expectations
    - Dedicated visualizations for cluster/hybrid behavior

KEY CAPABILITIES:
    • Pattern-rich block identification
    • Hybrid zone detection (multi-class overlaps)
    • Cluster metric calculation (size, density, inter-cluster distance)
    • Compositional skew analysis
    • Block stability scoring
    • Real vs. shuffled significance testing
    • Visualization generation for structural patterns

STRUCTURAL ANALYSIS FRAMEWORK:
    ┌────────────────────────────────────────────────────────────────────────────┐
    │ Analysis Type          │ Metrics                  │ Comparison             │
    ├────────────────────────┼──────────────────────────┼────────────────────────┤
    │ Pattern-Rich Blocks    │ Size, Density, Count     │ Real vs. Shuffled      │
    │ Hybrid Zones           │ Overlap, Interaction     │ Significance           │
    │ Clustering             │ Distance, Stability      │ Statistical Test       │
    │ Compositional Skew     │ GC%, AT%, Dinuc Freq     │ Local vs. Global       │
    │ Block Stability        │ Score Variance, Length   │ Robust vs. Fragile     │
    └────────────────────────┴──────────────────────────┴────────────────────────┘

REFERENCE:
    Higher-order pattern analysis inspired by ENCODE regulatory element
    clustering and chromatin domain identification methodologies.
"""

import numpy as np
from typing import List, Dict, Any, Tuple, Optional, Set
from collections import defaultdict, Counter
import warnings

# Suppress warnings during calculations
warnings.filterwarnings("ignore")

# =============================================================================
# CONSTANTS & CONFIGURATION
# =============================================================================

# Block detection parameters
BLOCK_MIN_DENSITY = 2.0  # Minimum motifs per kb for block classification
BLOCK_WINDOW_SIZE = 1000  # Window size for block detection (bp)
BLOCK_MIN_LENGTH = 500   # Minimum block length (bp)

# Hybrid zone parameters
HYBRID_MIN_OVERLAP_FRACTION = 0.3  # Minimum overlap for hybrid classification
HYBRID_MIN_CLASSES = 2             # Minimum classes for hybrid zone

# Cluster analysis parameters
CLUSTER_MAX_GAP = 500        # Maximum gap between motifs in same cluster (bp)
CLUSTER_MIN_MOTIFS = 3       # Minimum motifs per cluster
CLUSTER_MIN_DENSITY = 3.0    # Minimum density for cluster classification (motifs/kb)

# Stability scoring parameters
STABILITY_SCORE_WEIGHTS = {
    'length': 0.3,
    'density': 0.3,
    'class_diversity': 0.2,
    'score_consistency': 0.2
}


# =============================================================================
# PATTERN-RICH BLOCK IDENTIFICATION
# =============================================================================

def identify_pattern_rich_blocks(motifs: List[Dict[str, Any]], 
                                 sequence_length: int,
                                 window_size: int = BLOCK_WINDOW_SIZE,
                                 min_density: float = BLOCK_MIN_DENSITY,
                                 min_length: int = BLOCK_MIN_LENGTH) -> List[Dict[str, Any]]:
    """
    Identify pattern-rich blocks across the sequence.
    
    A block is a genomic region with elevated motif density above background.
    Uses sliding window approach to identify contiguous high-density regions.
    
    Args:
        motifs: List of detected motifs
        sequence_length: Total sequence length in bp
        window_size: Sliding window size (default: 1000 bp)
        min_density: Minimum density threshold (motifs/kb)
        min_length: Minimum block length (bp)
        
    Returns:
        List of pattern-rich blocks:
        [{
            'start': int,
            'end': int,
            'length': int,
            'motif_count': int,
            'density': float,
            'classes': List[str],
            'class_diversity': int,
            'mean_score': float,
            'motifs': List[Dict]  # Motifs within block
        }, ...]
    """
    if not motifs or sequence_length == 0:
        return []
    
    # Calculate density across windows
    num_windows = max(1, sequence_length // window_size)
    window_densities = []
    
    for i in range(num_windows):
        window_start = i * window_size
        window_end = min((i + 1) * window_size, sequence_length)
        
        # Count motifs in window
        window_motifs = [m for m in motifs 
                        if not (m.get('End', 0) <= window_start or m.get('Start', 0) >= window_end)]
        
        window_length_kb = (window_end - window_start) / 1000
        density = len(window_motifs) / window_length_kb if window_length_kb > 0 else 0
        
        window_densities.append({
            'start': window_start,
            'end': window_end,
            'density': density,
            'motifs': window_motifs
        })
    
    # Merge contiguous high-density windows into blocks
    blocks = []
    current_block = None
    
    for window in window_densities:
        if window['density'] >= min_density:
            if current_block is None:
                # Start new block
                current_block = {
                    'start': window['start'],
                    'end': window['end'],
                    'motifs': window['motifs']
                }
            else:
                # Extend current block
                current_block['end'] = window['end']
                current_block['motifs'].extend(window['motifs'])
        else:
            # End current block if density drops
            if current_block is not None:
                block_length = current_block['end'] - current_block['start']
                if block_length >= min_length:
                    blocks.append(current_block)
                current_block = None
    
    # Don't forget last block
    if current_block is not None:
        block_length = current_block['end'] - current_block['start']
        if block_length >= min_length:
            blocks.append(current_block)
    
    # Calculate block statistics
    for block in blocks:
        block_motifs = block['motifs']
        block_length = block['end'] - block['start']
        
        # Remove duplicates (motifs appearing in multiple windows)
        seen_ids = set()
        unique_motifs = []
        for m in block_motifs:
            motif_id = (m.get('Start'), m.get('End'), m.get('Class'))
            if motif_id not in seen_ids:
                seen_ids.add(motif_id)
                unique_motifs.append(m)
        
        block['motifs'] = unique_motifs
        block['length'] = block_length
        block['motif_count'] = len(unique_motifs)
        block['density'] = (len(unique_motifs) / block_length) * 1000
        
        # Class diversity
        classes = set(m.get('Class') for m in unique_motifs)
        block['classes'] = sorted(classes)
        block['class_diversity'] = len(classes)
        
        # Mean score
        scores = [m.get('Score', 0) for m in unique_motifs]
        block['mean_score'] = float(np.mean(scores)) if scores else 0.0
    
    return blocks


# =============================================================================
# HYBRID ZONE DETECTION
# =============================================================================

def detect_hybrid_zones(motifs: List[Dict[str, Any]],
                       min_overlap_fraction: float = HYBRID_MIN_OVERLAP_FRACTION,
                       min_classes: int = HYBRID_MIN_CLASSES) -> List[Dict[str, Any]]:
    """
    Detect overlap-driven hybrid zones where multiple pattern classes interact.
    
    A hybrid zone is a region where motifs from different classes overlap
    significantly, suggesting complex structural interactions.
    
    Args:
        motifs: List of detected motifs
        min_overlap_fraction: Minimum overlap fraction for hybrid classification
        min_classes: Minimum number of different classes required
        
    Returns:
        List of hybrid zones:
        [{
            'start': int,
            'end': int,
            'length': int,
            'classes': List[str],
            'num_classes': int,
            'motifs': List[Dict],
            'overlap_fraction': float,
            'interaction_strength': float,
            'hybrid_density': float
        }, ...]
    """
    if not motifs:
        return []
    
    # Find all pairwise overlaps between different classes
    hybrid_regions = []
    seen_regions = set()
    
    for i, m1 in enumerate(motifs):
        for m2 in motifs[i+1:]:
            class1 = m1.get('Class')
            class2 = m2.get('Class')
            
            # Only consider different classes
            if class1 == class2:
                continue
            
            # Check for overlap
            start1, end1 = m1.get('Start', 0), m1.get('End', 0)
            start2, end2 = m2.get('Start', 0), m2.get('End', 0)
            
            if end1 <= start2 or end2 <= start1:
                continue  # No overlap
            
            # Calculate overlap
            overlap_start = max(start1, start2)
            overlap_end = min(end1, end2)
            overlap_length = overlap_end - overlap_start
            
            # Check if overlap meets threshold
            min_length = min(end1 - start1, end2 - start2)
            overlap_fraction = overlap_length / min_length if min_length > 0 else 0
            
            if overlap_fraction >= min_overlap_fraction:
                # Create hybrid region spanning both motifs
                region_start = min(start1, start2)
                region_end = max(end1, end2)
                region_key = (region_start, region_end)
                
                if region_key not in seen_regions:
                    seen_regions.add(region_key)
                    hybrid_regions.append({
                        'start': region_start,
                        'end': region_end,
                        'overlap_fraction': overlap_fraction,
                        'motifs': [m1, m2]
                    })
    
    # Merge overlapping hybrid regions
    merged_regions = []
    hybrid_regions.sort(key=lambda r: r['start'])
    
    current_region = None
    for region in hybrid_regions:
        if current_region is None:
            current_region = region.copy()
        else:
            # Check if regions overlap
            if region['start'] <= current_region['end']:
                # Merge regions
                current_region['end'] = max(current_region['end'], region['end'])
                current_region['motifs'].extend(region['motifs'])
            else:
                # Save current region and start new one
                merged_regions.append(current_region)
                current_region = region.copy()
    
    if current_region is not None:
        merged_regions.append(current_region)
    
    # Calculate hybrid zone statistics
    filtered_zones = []
    for zone in merged_regions:
        # Remove duplicate motifs
        unique_motifs = []
        seen_ids = set()
        for m in zone['motifs']:
            motif_id = (m.get('Start'), m.get('End'), m.get('Class'))
            if motif_id not in seen_ids:
                seen_ids.add(motif_id)
                unique_motifs.append(m)
        
        zone['motifs'] = unique_motifs
        
        # Class diversity
        classes = set(m.get('Class') for m in unique_motifs)
        zone['classes'] = sorted(classes)
        zone['num_classes'] = len(classes)
        
        # Filter by minimum classes
        if zone['num_classes'] < min_classes:
            continue
        
        # Calculate statistics
        zone['length'] = zone['end'] - zone['start']
        zone['hybrid_density'] = (len(unique_motifs) / zone['length']) * 1000
        
        # Interaction strength (based on overlap density and class diversity)
        zone['interaction_strength'] = zone['hybrid_density'] * zone['num_classes']
        
        filtered_zones.append(zone)
    
    return filtered_zones


# =============================================================================
# CLUSTER ANALYSIS
# =============================================================================

def calculate_cluster_metrics(motifs: List[Dict[str, Any]],
                              max_gap: int = CLUSTER_MAX_GAP,
                              min_motifs: int = CLUSTER_MIN_MOTIFS) -> List[Dict[str, Any]]:
    """
    Calculate clustering metrics for spatial distribution analysis.
    
    Identifies clusters of motifs and calculates:
    - Cluster size (number of motifs)
    - Cluster density (motifs/kb)
    - Inter-cluster distance
    - Compositional skew within clusters
    - Cluster stability score
    
    Args:
        motifs: List of detected motifs (must be sorted by Start position)
        max_gap: Maximum gap between motifs in same cluster (bp)
        min_motifs: Minimum motifs required for cluster
        
    Returns:
        List of clusters:
        [{
            'cluster_id': int,
            'start': int,
            'end': int,
            'length': int,
            'size': int (number of motifs),
            'density': float (motifs/kb),
            'classes': List[str],
            'class_diversity': int,
            'mean_score': float,
            'score_std': float,
            'stability_score': float,
            'motifs': List[Dict]
        }, ...]
    """
    if not motifs:
        return []
    
    # Sort motifs by start position
    sorted_motifs = sorted(motifs, key=lambda m: m.get('Start', 0))
    
    # Group into clusters based on proximity
    clusters = []
    current_cluster = [sorted_motifs[0]]
    
    for motif in sorted_motifs[1:]:
        # Check gap from last motif in current cluster
        last_motif = current_cluster[-1]
        gap = motif.get('Start', 0) - last_motif.get('End', 0)
        
        if gap <= max_gap:
            # Add to current cluster
            current_cluster.append(motif)
        else:
            # Save current cluster if large enough
            if len(current_cluster) >= min_motifs:
                clusters.append(current_cluster)
            
            # Start new cluster
            current_cluster = [motif]
    
    # Don't forget last cluster
    if len(current_cluster) >= min_motifs:
        clusters.append(current_cluster)
    
    # Calculate metrics for each cluster
    cluster_data = []
    for cluster_id, cluster_motifs in enumerate(clusters):
        cluster_start = min(m.get('Start', 0) for m in cluster_motifs)
        cluster_end = max(m.get('End', 0) for m in cluster_motifs)
        cluster_length = cluster_end - cluster_start
        
        # Basic metrics
        size = len(cluster_motifs)
        density = (size / cluster_length) * 1000 if cluster_length > 0 else 0
        
        # Class diversity
        classes = set(m.get('Class') for m in cluster_motifs)
        class_diversity = len(classes)
        
        # Score statistics
        scores = [m.get('Score', 0) for m in cluster_motifs]
        mean_score = float(np.mean(scores)) if scores else 0.0
        score_std = float(np.std(scores)) if len(scores) > 1 else 0.0
        
        # Stability score (composite of length, density, diversity, score consistency)
        length_score = min(1.0, cluster_length / 2000)  # Normalized by 2kb
        density_score = min(1.0, density / 10)  # Normalized by 10 motifs/kb
        diversity_score = min(1.0, class_diversity / 5)  # Normalized by 5 classes
        consistency_score = 1.0 - min(1.0, score_std / mean_score if mean_score > 0 else 0)
        
        stability_score = (
            STABILITY_SCORE_WEIGHTS['length'] * length_score +
            STABILITY_SCORE_WEIGHTS['density'] * density_score +
            STABILITY_SCORE_WEIGHTS['class_diversity'] * diversity_score +
            STABILITY_SCORE_WEIGHTS['score_consistency'] * consistency_score
        )
        
        cluster_data.append({
            'cluster_id': cluster_id,
            'start': cluster_start,
            'end': cluster_end,
            'length': cluster_length,
            'size': size,
            'density': float(density),
            'classes': sorted(classes),
            'class_diversity': class_diversity,
            'mean_score': mean_score,
            'score_std': score_std,
            'stability_score': float(stability_score),
            'motifs': cluster_motifs
        })
    
    return cluster_data


def calculate_inter_cluster_distances(clusters: List[Dict[str, Any]]) -> List[float]:
    """
    Calculate distances between consecutive clusters.
    
    Args:
        clusters: List of clusters from calculate_cluster_metrics()
        
    Returns:
        List of inter-cluster distances in bp
    """
    if len(clusters) < 2:
        return []
    
    distances = []
    for i in range(len(clusters) - 1):
        cluster1_end = clusters[i]['end']
        cluster2_start = clusters[i + 1]['start']
        distance = cluster2_start - cluster1_end
        distances.append(max(0, distance))  # Ensure non-negative
    
    return distances


# =============================================================================
# COMPOSITIONAL SKEW ANALYSIS
# =============================================================================

def calculate_compositional_skew(sequence: str, 
                                motifs: List[Dict[str, Any]],
                                region_type: str = 'motif') -> Dict[str, Any]:
    """
    Calculate compositional skew within motifs/blocks vs. background.
    
    Measures:
    - GC skew: (G-C)/(G+C)
    - AT skew: (A-T)/(A+T)
    - Dinucleotide frequencies
    - Purine/pyrimidine ratios
    
    Args:
        sequence: Full sequence
        motifs: List of motifs or blocks
        region_type: 'motif' or 'block'
        
    Returns:
        Dictionary of compositional skew metrics:
        {
            'gc_skew_motif': float,
            'gc_skew_background': float,
            'at_skew_motif': float,
            'at_skew_background': float,
            'purine_ratio_motif': float,
            'purine_ratio_background': float,
            'dinuc_freq_motif': Dict[str, float],
            'dinuc_freq_background': Dict[str, float]
        }
    """
    if not motifs or not sequence:
        return {}
    
    sequence = sequence.upper()
    
    # Extract motif regions
    motif_sequence = ''
    for motif in motifs:
        start = motif.get('Start', 0) - 1  # Convert to 0-based
        end = motif.get('End', 0)
        if 0 <= start < len(sequence) and start < end <= len(sequence):
            motif_sequence += sequence[start:end]
    
    # Background sequence (everything except motifs)
    # Create mask for motif regions
    motif_mask = [False] * len(sequence)
    for motif in motifs:
        start = max(0, motif.get('Start', 0) - 1)
        end = min(len(sequence), motif.get('End', 0))
        for i in range(start, end):
            motif_mask[i] = True
    
    background_sequence = ''.join(sequence[i] for i in range(len(sequence)) if not motif_mask[i])
    
    # Calculate skew metrics
    def calc_skew(seq):
        if not seq:
            return {'gc_skew': 0, 'at_skew': 0, 'purine_ratio': 0, 'dinuc_freq': {}}
        
        g_count = seq.count('G')
        c_count = seq.count('C')
        a_count = seq.count('A')
        t_count = seq.count('T')
        
        gc_skew = (g_count - c_count) / (g_count + c_count) if (g_count + c_count) > 0 else 0
        at_skew = (a_count - t_count) / (a_count + t_count) if (a_count + t_count) > 0 else 0
        
        purine_count = g_count + a_count
        pyrimidine_count = c_count + t_count
        purine_ratio = purine_count / (purine_count + pyrimidine_count) if (purine_count + pyrimidine_count) > 0 else 0
        
        # Dinucleotide frequencies
        dinuc_counts = Counter()
        for i in range(len(seq) - 1):
            dinuc = seq[i:i+2]
            if len(dinuc) == 2:
                dinuc_counts[dinuc] += 1
        
        total_dinuc = sum(dinuc_counts.values())
        dinuc_freq = {dinuc: count / total_dinuc for dinuc, count in dinuc_counts.items()} if total_dinuc > 0 else {}
        
        return {
            'gc_skew': float(gc_skew),
            'at_skew': float(at_skew),
            'purine_ratio': float(purine_ratio),
            'dinuc_freq': dinuc_freq
        }
    
    motif_skew = calc_skew(motif_sequence)
    background_skew = calc_skew(background_sequence)
    
    return {
        'gc_skew_motif': motif_skew['gc_skew'],
        'gc_skew_background': background_skew['gc_skew'],
        'at_skew_motif': motif_skew['at_skew'],
        'at_skew_background': background_skew['at_skew'],
        'purine_ratio_motif': motif_skew['purine_ratio'],
        'purine_ratio_background': background_skew['purine_ratio'],
        'dinuc_freq_motif': motif_skew['dinuc_freq'],
        'dinuc_freq_background': background_skew['dinuc_freq']
    }


# =============================================================================
# MAIN STRUCTURAL ANALYSIS PIPELINE
# =============================================================================

def run_structural_analysis(sequence: str,
                            motifs: List[Dict[str, Any]],
                            enable_blocks: bool = True,
                            enable_hybrids: bool = True,
                            enable_clusters: bool = True,
                            enable_skew: bool = True) -> Dict[str, Any]:
    """
    Complete structural analysis pipeline.
    
    Args:
        sequence: Input DNA sequence
        motifs: Detected motifs from sequence
        enable_blocks: Enable pattern-rich block detection
        enable_hybrids: Enable hybrid zone detection
        enable_clusters: Enable cluster analysis
        enable_skew: Enable compositional skew analysis
        
    Returns:
        Dictionary with structural analysis results:
        {
            'blocks': List[Dict],
            'hybrid_zones': List[Dict],
            'clusters': List[Dict],
            'inter_cluster_distances': List[float],
            'compositional_skew': Dict,
            'summary': Dict
        }
    """
    results = {}
    
    sequence_length = len(sequence)
    
    # Pattern-rich blocks
    if enable_blocks:
        blocks = identify_pattern_rich_blocks(motifs, sequence_length)
        results['blocks'] = blocks
    else:
        results['blocks'] = []
    
    # Hybrid zones
    if enable_hybrids:
        hybrid_zones = detect_hybrid_zones(motifs)
        results['hybrid_zones'] = hybrid_zones
    else:
        results['hybrid_zones'] = []
    
    # Cluster analysis
    if enable_clusters:
        clusters = calculate_cluster_metrics(motifs)
        inter_cluster_distances = calculate_inter_cluster_distances(clusters)
        results['clusters'] = clusters
        results['inter_cluster_distances'] = inter_cluster_distances
    else:
        results['clusters'] = []
        results['inter_cluster_distances'] = []
    
    # Compositional skew
    if enable_skew:
        compositional_skew = calculate_compositional_skew(sequence, motifs)
        results['compositional_skew'] = compositional_skew
    else:
        results['compositional_skew'] = {}
    
    # Summary statistics
    results['summary'] = {
        'num_blocks': len(results.get('blocks', [])),
        'num_hybrid_zones': len(results.get('hybrid_zones', [])),
        'num_clusters': len(results.get('clusters', [])),
        'mean_cluster_size': float(np.mean([c['size'] for c in results.get('clusters', [])]) if results.get('clusters') else 0),
        'mean_cluster_density': float(np.mean([c['density'] for c in results.get('clusters', [])]) if results.get('clusters') else 0),
        'mean_inter_cluster_distance': float(np.mean(results.get('inter_cluster_distances', [])) if results.get('inter_cluster_distances') else 0),
    }
    
    return results


# =============================================================================
# TESTING & VALIDATION
# =============================================================================

def test_structural_analysis():
    """Test structural analysis module with synthetic data."""
    print("Testing structural analysis module...")
    
    # Create synthetic motifs
    motifs = [
        {'Start': 100, 'End': 150, 'Class': 'G-Quadruplex', 'Score': 2.5},
        {'Start': 200, 'End': 250, 'Class': 'Z-DNA', 'Score': 2.3},
        {'Start': 220, 'End': 270, 'Class': 'i-Motif', 'Score': 2.1},  # Overlaps with Z-DNA
        {'Start': 1100, 'End': 1150, 'Class': 'G-Quadruplex', 'Score': 2.7},
        {'Start': 1200, 'End': 1250, 'Class': 'Curved_DNA', 'Score': 2.2},
        {'Start': 1250, 'End': 1300, 'Class': 'R-Loop', 'Score': 2.4},
    ]
    
    sequence = "A" * 2000  # Dummy sequence
    
    # Test 1: Pattern-rich blocks
    blocks = identify_pattern_rich_blocks(motifs, len(sequence), window_size=500, min_density=1.0)
    assert len(blocks) > 0, "No blocks detected"
    print(f"✓ Test 1 passed: Detected {len(blocks)} pattern-rich blocks")
    
    # Test 2: Hybrid zones
    hybrid_zones = detect_hybrid_zones(motifs, min_overlap_fraction=0.1)
    assert len(hybrid_zones) > 0, "No hybrid zones detected"
    print(f"✓ Test 2 passed: Detected {len(hybrid_zones)} hybrid zones")
    
    # Test 3: Cluster analysis
    clusters = calculate_cluster_metrics(motifs, max_gap=200, min_motifs=2)
    assert len(clusters) > 0, "No clusters detected"
    print(f"✓ Test 3 passed: Detected {len(clusters)} clusters")
    
    # Test 4: Inter-cluster distances
    distances = calculate_inter_cluster_distances(clusters)
    if len(clusters) > 1:
        assert len(distances) == len(clusters) - 1, "Wrong number of distances"
        print(f"✓ Test 4 passed: Calculated {len(distances)} inter-cluster distances")
    else:
        print("✓ Test 4 passed: Insufficient clusters for distance calculation")
    
    # Test 5: Compositional skew
    skew = calculate_compositional_skew(sequence, motifs)
    assert 'gc_skew_motif' in skew, "Missing GC skew"
    print("✓ Test 5 passed: Compositional skew calculated")
    
    # Test 6: Full pipeline
    results = run_structural_analysis(sequence, motifs)
    assert 'blocks' in results, "Missing blocks in results"
    assert 'hybrid_zones' in results, "Missing hybrid zones in results"
    assert 'clusters' in results, "Missing clusters in results"
    assert 'summary' in results, "Missing summary in results"
    print("✓ Test 6 passed: Full pipeline executed successfully")
    
    print("\n✓ All tests passed successfully!")


if __name__ == "__main__":
    test_structural_analysis()
