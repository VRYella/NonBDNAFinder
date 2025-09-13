#!/usr/bin/env python3
"""
Enhanced overlap resolution system for Non-B DNA motif detection.

This module implements sophisticated algorithms for resolving overlapping motifs
within the same subclass based on scoring criteria and scientific logic.

Key Features:
- Score-based selection for overlapping motifs within same subclass
- Dynamic merging strategies based on motif types
- Integration with existing Hyperscan detection pipeline
- Configurable overlap thresholds and resolution strategies

Author: Dr. Venkata Rajesh Yella
"""

from typing import List, Dict, Set, Tuple, Optional, Union
import logging
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum

# Add parent directory to path to import from other modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from intervaltree import IntervalTree, Interval
    INTERVALTREE_AVAILABLE = True
except ImportError:
    INTERVALTREE_AVAILABLE = False
    logging.warning("intervaltree not available - using simple overlap detection")

from REGISTRIES.motifs.base import Candidate

logger = logging.getLogger(__name__)


class OverlapStrategy(Enum):
    """Strategies for resolving overlapping motifs"""
    HIGHEST_SCORE = "highest_score"  # Keep motif with highest score
    LONGEST_MOTIF = "longest_motif"  # Keep longest motif
    MERGE_COMPATIBLE = "merge_compatible"  # Merge compatible overlapping motifs
    SCIENTIFIC_PRIORITY = "scientific_priority"  # Use scientific hierarchy
    KEEP_ALL = "keep_all"  # Keep all overlaps (mark as hybrid)


@dataclass
class OverlapConfig:
    """Configuration for overlap resolution"""
    strategy: OverlapStrategy = OverlapStrategy.HIGHEST_SCORE
    min_overlap_percent: float = 0.1  # Minimum overlap percentage to consider
    merge_threshold: float = 0.8  # Overlap threshold for merging
    same_class_only: bool = True  # Only resolve overlaps within same class
    preserve_hybrid: bool = True  # Preserve hybrid motifs


class EnhancedOverlapResolver:
    """
    Advanced overlap resolution system for motif detection.
    
    This class provides comprehensive overlap resolution capabilities including:
    - Detection of overlapping motifs within and across classes
    - Score-based selection strategies
    - Dynamic merging of compatible motifs
    - Scientific priority-based resolution
    """
    
    def __init__(self, config: Optional[OverlapConfig] = None):
        """
        Initialize the overlap resolver.
        
        Args:
            config: Configuration for overlap resolution strategy
        """
        self.config = config or OverlapConfig()
        self.logger = logging.getLogger(self.__class__.__name__)
        
    def resolve_overlaps(self, candidates: List[Candidate]) -> List[Candidate]:
        """
        Resolve overlapping motifs based on configured strategy.
        
        Args:
            candidates: List of motif candidates to process
            
        Returns:
            List of candidates with overlaps resolved
        """
        if not candidates:
            return candidates
            
        self.logger.info(f"Resolving overlaps for {len(candidates)} candidates")
        
        # Group candidates by sequence and class/subclass
        grouped_candidates = self._group_candidates(candidates)
        
        resolved_candidates = []
        for group_key, group_candidates in grouped_candidates.items():
            if len(group_candidates) <= 1:
                resolved_candidates.extend(group_candidates)
                continue
                
            # Detect overlaps within group
            overlaps = self._detect_overlaps(group_candidates)
            
            # Resolve overlaps based on strategy
            resolved_group = self._resolve_group_overlaps(group_candidates, overlaps)
            resolved_candidates.extend(resolved_group)
        
        self.logger.info(f"Resolved to {len(resolved_candidates)} candidates")
        return resolved_candidates
    
    def _group_candidates(self, candidates: List[Candidate]) -> Dict[tuple, List[Candidate]]:
        """
        Group candidates by sequence, class, and subclass for overlap resolution.
        
        Args:
            candidates: List of candidates to group
            
        Returns:
            Dictionary mapping group keys to candidate lists
        """
        groups = defaultdict(list)
        
        for candidate in candidates:
            if self.config.same_class_only:
                # Group by sequence and class only (not subclass to allow intra-class resolution)
                group_key = (
                    candidate.sequence_name,
                    candidate.class_name
                )
            else:
                # Group only by sequence for cross-class resolution
                group_key = (candidate.sequence_name,)
            
            groups[group_key].append(candidate)
        
        return dict(groups)
    
    def _detect_overlaps(self, candidates: List[Candidate]) -> List[Tuple[int, int, float]]:
        """
        Detect overlapping candidates within a group.
        
        Args:
            candidates: List of candidates to check for overlaps
            
        Returns:
            List of (index1, index2, overlap_percent) tuples
        """
        overlaps = []
        
        if INTERVALTREE_AVAILABLE:
            overlaps = self._detect_overlaps_intervaltree(candidates)
        else:
            overlaps = self._detect_overlaps_simple(candidates)
        
        # Filter by minimum overlap percentage
        filtered_overlaps = [
            (i, j, percent) for i, j, percent in overlaps
            if percent >= self.config.min_overlap_percent
        ]
        
        return filtered_overlaps
    
    def _detect_overlaps_intervaltree(self, candidates: List[Candidate]) -> List[Tuple[int, int, float]]:
        """
        Detect overlaps using IntervalTree for efficiency.
        
        Args:
            candidates: List of candidates to check
            
        Returns:
            List of overlap tuples
        """
        overlaps = []
        intervals = IntervalTree()
        
        # Build interval tree
        for i, candidate in enumerate(candidates):
            intervals[candidate.start:candidate.end + 1] = i
        
        # Find overlaps
        for i, candidate in enumerate(candidates):
            overlapping_intervals = intervals[candidate.start:candidate.end + 1]
            
            for interval in overlapping_intervals:
                j = interval.data
                if i < j:  # Avoid duplicate pairs
                    overlap_percent = self._calculate_overlap_percent(
                        candidates[i], candidates[j]
                    )
                    overlaps.append((i, j, overlap_percent))
        
        return overlaps
    
    def _detect_overlaps_simple(self, candidates: List[Candidate]) -> List[Tuple[int, int, float]]:
        """
        Simple O(n²) overlap detection when IntervalTree is not available.
        
        Args:
            candidates: List of candidates to check
            
        Returns:
            List of overlap tuples
        """
        overlaps = []
        
        for i in range(len(candidates)):
            for j in range(i + 1, len(candidates)):
                overlap_percent = self._calculate_overlap_percent(
                    candidates[i], candidates[j]
                )
                if overlap_percent > 0:
                    overlaps.append((i, j, overlap_percent))
        
        return overlaps
    
    def _calculate_overlap_percent(self, candidate1: Candidate, candidate2: Candidate) -> float:
        """
        Calculate the percentage overlap between two candidates.
        
        Args:
            candidate1: First candidate
            candidate2: Second candidate
            
        Returns:
            Overlap percentage (0.0 to 1.0)
        """
        # Calculate overlap region
        overlap_start = max(candidate1.start, candidate2.start)
        overlap_end = min(candidate1.end, candidate2.end)
        
        if overlap_start > overlap_end:
            return 0.0
        
        overlap_length = overlap_end - overlap_start + 1
        
        # Calculate as percentage of smaller motif
        min_length = min(candidate1.length, candidate2.length)
        return overlap_length / min_length if min_length > 0 else 0.0
    
    def _resolve_group_overlaps(self, candidates: List[Candidate], 
                              overlaps: List[Tuple[int, int, float]]) -> List[Candidate]:
        """
        Resolve overlaps within a group based on configured strategy.
        
        Args:
            candidates: List of candidates in the group
            overlaps: List of detected overlaps
            
        Returns:
            List of resolved candidates
        """
        if not overlaps:
            return candidates
        
        if self.config.strategy == OverlapStrategy.HIGHEST_SCORE:
            return self._resolve_by_highest_score(candidates, overlaps)
        elif self.config.strategy == OverlapStrategy.LONGEST_MOTIF:
            return self._resolve_by_longest_motif(candidates, overlaps)
        elif self.config.strategy == OverlapStrategy.MERGE_COMPATIBLE:
            return self._resolve_by_merging(candidates, overlaps)
        elif self.config.strategy == OverlapStrategy.SCIENTIFIC_PRIORITY:
            return self._resolve_by_scientific_priority(candidates, overlaps)
        elif self.config.strategy == OverlapStrategy.KEEP_ALL:
            return self._resolve_keep_all(candidates, overlaps)
        else:
            self.logger.warning(f"Unknown strategy {self.config.strategy}, using HIGHEST_SCORE")
            return self._resolve_by_highest_score(candidates, overlaps)
    
    def _resolve_by_highest_score(self, candidates: List[Candidate], 
                                overlaps: List[Tuple[int, int, float]]) -> List[Candidate]:
        """
        Resolve overlaps by keeping candidates with highest scores.
        
        Args:
            candidates: List of candidates
            overlaps: List of overlap information
            
        Returns:
            List of resolved candidates
        """
        to_remove = set()
        
        for i, j, overlap_percent in overlaps:
            if i in to_remove or j in to_remove:
                continue
            
            candidate_i = candidates[i]
            candidate_j = candidates[j]
            
            # Compare scores (handle None scores)
            score_i = candidate_i.raw_score or 0.0
            score_j = candidate_j.raw_score or 0.0
            
            if score_i > score_j:
                to_remove.add(j)
                # Update overlap classes for the kept candidate
                if candidate_j.class_name not in (candidate_i.overlap_classes or []):
                    if candidate_i.overlap_classes is None:
                        candidate_i.overlap_classes = []
                    candidate_i.overlap_classes.append(candidate_j.class_name)
            elif score_j > score_i:
                to_remove.add(i)
                # Update overlap classes for the kept candidate
                if candidate_i.class_name not in (candidate_j.overlap_classes or []):
                    if candidate_j.overlap_classes is None:
                        candidate_j.overlap_classes = []
                    candidate_j.overlap_classes.append(candidate_i.class_name)
            else:
                # Equal scores, keep longer motif
                if candidate_i.length >= candidate_j.length:
                    to_remove.add(j)
                else:
                    to_remove.add(i)
        
        # Return candidates not marked for removal
        resolved = [candidates[i] for i in range(len(candidates)) if i not in to_remove]
        
        self.logger.debug(f"Score-based resolution: {len(candidates)} -> {len(resolved)}")
        return resolved
    
    def _resolve_by_longest_motif(self, candidates: List[Candidate], 
                                overlaps: List[Tuple[int, int, float]]) -> List[Candidate]:
        """
        Resolve overlaps by keeping longest motifs.
        
        Args:
            candidates: List of candidates
            overlaps: List of overlap information
            
        Returns:
            List of resolved candidates
        """
        to_remove = set()
        
        for i, j, overlap_percent in overlaps:
            if i in to_remove or j in to_remove:
                continue
            
            candidate_i = candidates[i]
            candidate_j = candidates[j]
            
            if candidate_i.length > candidate_j.length:
                to_remove.add(j)
            elif candidate_j.length > candidate_i.length:
                to_remove.add(i)
            else:
                # Equal length, use score as tiebreaker
                score_i = candidate_i.raw_score or 0.0
                score_j = candidate_j.raw_score or 0.0
                if score_i >= score_j:
                    to_remove.add(j)
                else:
                    to_remove.add(i)
        
        resolved = [candidates[i] for i in range(len(candidates)) if i not in to_remove]
        
        self.logger.debug(f"Length-based resolution: {len(candidates)} -> {len(resolved)}")
        return resolved
    
    def _resolve_by_merging(self, candidates: List[Candidate], 
                          overlaps: List[Tuple[int, int, float]]) -> List[Candidate]:
        """
        Resolve overlaps by merging compatible motifs.
        
        Args:
            candidates: List of candidates
            overlaps: List of overlap information
            
        Returns:
            List of resolved candidates
        """
        # This is a simplified implementation
        # In practice, merging would require sophisticated biological logic
        
        merged_groups = []
        processed = set()
        
        for i, j, overlap_percent in overlaps:
            if i in processed or j in processed:
                continue
            
            if overlap_percent >= self.config.merge_threshold:
                # Merge the candidates
                candidate_i = candidates[i]
                candidate_j = candidates[j]
                
                merged = self._merge_candidates(candidate_i, candidate_j)
                merged_groups.append(merged)
                processed.add(i)
                processed.add(j)
        
        # Add unprocessed candidates
        for i, candidate in enumerate(candidates):
            if i not in processed:
                merged_groups.append(candidate)
        
        self.logger.debug(f"Merge-based resolution: {len(candidates)} -> {len(merged_groups)}")
        return merged_groups
    
    def _merge_candidates(self, candidate1: Candidate, candidate2: Candidate) -> Candidate:
        """
        Merge two compatible candidates into one.
        
        Args:
            candidate1: First candidate
            candidate2: Second candidate
            
        Returns:
            Merged candidate
        """
        # Use boundaries of the larger motif
        start = min(candidate1.start, candidate2.start)
        end = max(candidate1.end, candidate2.end)
        
        # Use higher score
        score1 = candidate1.raw_score or 0.0
        score2 = candidate2.raw_score or 0.0
        higher_score_candidate = candidate1 if score1 >= score2 else candidate2
        
        # Create merged candidate
        merged = Candidate(
            sequence_name=candidate1.sequence_name,
            contig=candidate1.contig,
            class_id=higher_score_candidate.class_id,
            class_name=higher_score_candidate.class_name,
            subclass=f"merged_{candidate1.subclass}_{candidate2.subclass}",
            motif_id=max(candidate1.motif_id, candidate2.motif_id),
            start=start,
            end=end,
            length=end - start + 1,
            matched_seq=higher_score_candidate.matched_seq,  # Use sequence from higher-scoring candidate
            pattern_name=f"merged_{candidate1.pattern_name}_{candidate2.pattern_name}",
            raw_score=max(score1, score2),
            scoring_method=f"merged_{candidate1.scoring_method}_{candidate2.scoring_method}",
            overlap_classes=[candidate1.class_name, candidate2.class_name]
        )
        
        return merged
    
    def _resolve_by_scientific_priority(self, candidates: List[Candidate], 
                                      overlaps: List[Tuple[int, int, float]]) -> List[Candidate]:
        """
        Resolve overlaps using scientific priority hierarchy.
        
        The priority order is based on biological significance:
        1. G-quadruplex (most studied)
        2. i-motif (C-rich quadruplex)
        3. Z-DNA (alternate conformation)
        4. Triplex DNA (triple helix)
        5. Cruciform (DNA crossover)
        6. R-loops (RNA-DNA hybrid)
        7. Curved DNA (bent structure)
        8. Slipped DNA (repeat expansion)
        
        Args:
            candidates: List of candidates
            overlaps: List of overlap information
            
        Returns:
            List of resolved candidates
        """
        priority_order = {
            'g_quadruplex': 1,
            'i_motif': 2,
            'z_dna': 3,
            'triplex': 4,
            'cruciform': 5,
            'r_loop': 6,
            'curved_dna': 7,
            'slipped_dna': 8,
            'hybrid': 9,
            'cluster': 10
        }
        
        to_remove = set()
        
        for i, j, overlap_percent in overlaps:
            if i in to_remove or j in to_remove:
                continue
            
            candidate_i = candidates[i]
            candidate_j = candidates[j]
            
            priority_i = priority_order.get(candidate_i.class_name, 99)
            priority_j = priority_order.get(candidate_j.class_name, 99)
            
            if priority_i < priority_j:
                to_remove.add(j)
            elif priority_j < priority_i:
                to_remove.add(i)
            else:
                # Same priority, use score
                score_i = candidate_i.raw_score or 0.0
                score_j = candidate_j.raw_score or 0.0
                if score_i >= score_j:
                    to_remove.add(j)
                else:
                    to_remove.add(i)
        
        resolved = [candidates[i] for i in range(len(candidates)) if i not in to_remove]
        
        self.logger.debug(f"Priority-based resolution: {len(candidates)} -> {len(resolved)}")
        return resolved
    
    def _resolve_keep_all(self, candidates: List[Candidate], 
                         overlaps: List[Tuple[int, int, float]]) -> List[Candidate]:
        """
        Keep all overlapping motifs but mark them appropriately.
        
        Args:
            candidates: List of candidates
            overlaps: List of overlap information
            
        Returns:
            List of candidates with overlap information updated
        """
        # Create overlap mapping
        overlap_map = defaultdict(set)
        for i, j, overlap_percent in overlaps:
            overlap_map[i].add(j)
            overlap_map[j].add(i)
        
        # Update overlap classes for all candidates
        for i, candidate in enumerate(candidates):
            if i in overlap_map:
                overlapping_classes = []
                for j in overlap_map[i]:
                    overlapping_classes.append(candidates[j].class_name)
                
                if candidate.overlap_classes is None:
                    candidate.overlap_classes = []
                candidate.overlap_classes.extend(overlapping_classes)
        
        self.logger.debug(f"Keep-all resolution: maintained {len(candidates)} candidates")
        return candidates


def resolve_motif_overlaps(candidates: List[Candidate], 
                          strategy: OverlapStrategy = OverlapStrategy.HIGHEST_SCORE,
                          config: Optional[OverlapConfig] = None) -> List[Candidate]:
    """
    Convenience function for resolving motif overlaps.
    
    Args:
        candidates: List of motif candidates
        strategy: Overlap resolution strategy
        config: Optional configuration override
        
    Returns:
        List of candidates with overlaps resolved
    """
    if config is None:
        config = OverlapConfig(strategy=strategy)
    else:
        config.strategy = strategy
    
    resolver = EnhancedOverlapResolver(config)
    return resolver.resolve_overlaps(candidates)