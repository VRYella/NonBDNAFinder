"""SMID Clustering: O(n) 1D density clustering and window extension/merge.

No third-party dependencies beyond NumPy.

Public interface
----------------
    clusters = linear_1d_clustering(positions, epsilon, min_samples)
    windows  = extend_and_merge(clusters, max_extension, seq_len)
"""

from __future__ import annotations

from typing import List, Tuple

import numpy as np


def linear_1d_clustering(
    positions: np.ndarray,
    epsilon: int,
    min_samples: int,
) -> List[Tuple[int, int]]:
    """Cluster 1-D positions by proximity in O(n) time.

    Algorithm
    ---------
    1. Sort positions.
    2. Walk through sorted positions; start a new cluster whenever the gap
       to the next position exceeds *epsilon*.
    3. Discard clusters with fewer than *min_samples* hits.

    Parameters
    ----------
    positions : array-like of int
        0-based seed-hit start positions (unsorted).
    epsilon : int
        Maximum gap (bp) between consecutive hits to remain in same cluster.
    min_samples : int
        Minimum number of hits in a cluster to keep it.

    Returns
    -------
    list of (start, end) tuples
        Both coordinates are inclusive and in the same coordinate space as
        *positions*.  Empty list if no qualifying cluster is found.
    """
    if len(positions) == 0:
        return []

    sorted_pos = np.sort(np.asarray(positions, dtype=np.int64))

    clusters: List[Tuple[int, int]] = []
    cluster_start = int(sorted_pos[0])
    cluster_end = int(sorted_pos[0])
    count = 1

    for pos in sorted_pos[1:]:
        pos = int(pos)
        if pos - cluster_end <= epsilon:
            cluster_end = pos
            count += 1
        else:
            if count >= min_samples:
                clusters.append((cluster_start, cluster_end))
            cluster_start = pos
            cluster_end = pos
            count = 1

    # Don't forget the last open cluster
    if count >= min_samples:
        clusters.append((cluster_start, cluster_end))

    return clusters


def extend_and_merge(
    clusters: List[Tuple[int, int]],
    max_extension: int,
    seq_len: int,
) -> List[Tuple[int, int]]:
    """Extend cluster boundaries by *max_extension* bp and merge overlaps.

    Parameters
    ----------
    clusters : list of (start, end)
        Raw cluster windows (0-based, end-inclusive).
    max_extension : int
        Number of bp to add to each side before merging.
    seq_len : int
        Total sequence length; used to clip coordinates.

    Returns
    -------
    list of (start, end)
        Merged, extended windows (0-based, end-inclusive), clipped to
        [0, seq_len - 1].
    """
    if not clusters:
        return []

    # Extend each cluster, clipping to sequence bounds
    extended = [
        (max(0, s - max_extension), min(seq_len - 1, e + max_extension))
        for s, e in clusters
    ]

    # Sort by start, then merge overlapping / adjacent intervals
    extended.sort()
    merged: List[Tuple[int, int]] = [extended[0]]
    for start, end in extended[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))

    return merged
