"""
Z-DNA calculator implementing Kadane algorithm for maximum subarray detection.

This module provides the core Z-DNA detection functionality using dinucleotide
transition scoring and a modified Kadane algorithm to find high-scoring
contiguous regions that are likely to form Z-DNA structures.
"""

from Bio.Seq import Seq
import numpy as np
from typing import Optional, List, Tuple
from constants import Params


class ZDNACalculatorSeq(Seq):
    """
    Z-DNA calculator that extends Bio.Seq with scoring and subarray detection capabilities.
    
    This class implements multiple scoring methods and uses a Kadane-like algorithm
    to find maximum scoring subarrays representing potential Z-DNA forming regions.
    """

    def __init__(self, data: str, params: Optional[Params] = None) -> None:
        """
        Initialize Z-DNA calculator.
        
        Args:
            data: DNA sequence string
            params: Parameters object containing scoring configuration
        """
        super().__init__(data.upper())
        self.scoring_array: np.ndarray = np.array([])
        self.params = params if params is not None else Params()

    def zdna_calculator_layered(self) -> np.ndarray:
        """
        Calculate scoring array using layered approach with separate components.
        
        This method calculates scores for dinucleotide transitions, applies
        mismatch penalties, and adds cadence rewards and consecutive AT scoring.
        
        Returns:
            np.ndarray: Scoring array of length (sequence_length - 1)
        """
        if len(self) < 2:
            return np.array([])
            
        scoring_array = np.empty(len(self) - 1, dtype=float)
        mismatches_counter = 0
        consecutive_AT_counter = 0

        for i in range(len(self) - 1):
            transition = self[i] + self[i + 1]
            
            # Score the transition
            if transition in ("GC", "CG"):
                scoring_array[i] = self.params.GC_weight
                mismatches_counter = 0
                consecutive_AT_counter = 0
            elif transition in ("GT", "TG"):
                scoring_array[i] = self.params.GT_weight
                mismatches_counter = 0
                consecutive_AT_counter = 0
            elif transition in ("AC", "CA"):
                scoring_array[i] = self.params.AC_weight
                mismatches_counter = 0
                consecutive_AT_counter = 0
            elif transition in ("AT", "TA"):
                adjusted_weight = self.params.AT_weight
                # Add consecutive AT bonus/penalty
                if consecutive_AT_counter < len(self.params.consecutive_AT_scoring):
                    adjusted_weight += self.params.consecutive_AT_scoring[consecutive_AT_counter]
                else:
                    adjusted_weight += self.params.consecutive_AT_scoring[-1]
                scoring_array[i] = adjusted_weight
                consecutive_AT_counter += 1
                mismatches_counter = 0
            else:
                # Apply mismatch penalty
                mismatches_counter += 1
                consecutive_AT_counter = 0
                
                if self.params.mismatch_penalty_type == "exponential":
                    penalty = self.params.mismatch_penalty_starting_value ** mismatches_counter
                    scoring_array[i] = -min(penalty, 32000)  # Cap at 32000 to prevent overflow
                elif self.params.mismatch_penalty_type == "linear":
                    penalty = (self.params.mismatch_penalty_starting_value + 
                             self.params.mismatch_penalty_linear_delta * (mismatches_counter - 1))
                    scoring_array[i] = -penalty
                else:
                    raise ValueError(f"Unknown mismatch penalty type: {self.params.mismatch_penalty_type}")

            # Add cadence reward for valid transitions
            if transition in ("GC", "CG", "GT", "TG", "AC", "CA", "AT", "TA"):
                scoring_array[i] += self.params.cadence_reward

        return scoring_array

    def zdna_calculator_transitions(self) -> np.ndarray:
        """
        Calculate scoring array treating each transition individually.
        
        This is the primary scoring method that evaluates each dinucleotide
        transition and applies appropriate weights and penalties.
        
        Returns:
            np.ndarray: Scoring array of length (sequence_length - 1)
        """
        if len(self) < 2:
            return np.array([])
            
        scoring_array = np.empty(len(self) - 1, dtype=float)
        mismatches_counter = 0
        consecutive_AT_counter = 0

        for i in range(len(self) - 1):
            transition = self[i] + self[i + 1]
            
            if transition in ("GC", "CG"):
                scoring_array[i] = self.params.GC_weight
                mismatches_counter = 0
                consecutive_AT_counter = 0
            elif transition in ("GT", "TG"):
                scoring_array[i] = self.params.GT_weight
                mismatches_counter = 0
                consecutive_AT_counter = 0
            elif transition in ("AC", "CA"):
                scoring_array[i] = self.params.AC_weight
                mismatches_counter = 0
                consecutive_AT_counter = 0
            elif transition in ("AT", "TA"):
                adjusted_weight = self.params.AT_weight
                # Apply consecutive AT scoring
                if consecutive_AT_counter < len(self.params.consecutive_AT_scoring):
                    adjusted_weight += self.params.consecutive_AT_scoring[consecutive_AT_counter]
                else:
                    adjusted_weight += self.params.consecutive_AT_scoring[-1]
                scoring_array[i] = adjusted_weight
                consecutive_AT_counter += 1
                mismatches_counter = 0
            else:
                # Handle mismatch
                mismatches_counter += 1
                consecutive_AT_counter = 0
                
                if self.params.mismatch_penalty_type == "exponential":
                    penalty = self.params.mismatch_penalty_starting_value ** mismatches_counter
                    scoring_array[i] = -min(penalty, 32000)
                elif self.params.mismatch_penalty_type == "linear":
                    penalty = (self.params.mismatch_penalty_starting_value + 
                             self.params.mismatch_penalty_linear_delta * (mismatches_counter - 1))
                    scoring_array[i] = -penalty
                else:
                    raise ValueError(f"Unknown mismatch penalty type: {self.params.mismatch_penalty_type}")

        return scoring_array

    def subarrays_above_threshold(self) -> List[Tuple[int, int, float, str]]:
        """
        Find all subarrays above threshold using modified Kadane algorithm.
        
        This method implements a modified version of Kadane's maximum subarray
        algorithm to find all contiguous subsequences with scores above the
        specified threshold.
        
        Returns:
            List of tuples (start, end, score, substring) where:
            - start: starting position in sequence (0-based)
            - end: ending position in sequence (0-based, inclusive)
            - score: total score of the subarray
            - substring: the actual DNA sequence substring
        """
        # First, calculate the scoring array using transitions method
        self.scoring_array = self.zdna_calculator_transitions()
        
        if self.scoring_array.size == 0:
            return []

        subarrays = []
        n = len(self.scoring_array)
        
        # Initialize Kadane algorithm variables
        max_ending_here = 0
        current_start = 0
        candidate_start = 0
        candidate_score = 0
        max_score_seen = 0
        in_valid_subarray = False
        
        for i in range(n):
            # Extend current subarray or start new one
            if max_ending_here + self.scoring_array[i] < self.scoring_array[i]:
                max_ending_here = self.scoring_array[i]
                current_start = i
            else:
                max_ending_here += self.scoring_array[i]
            
            # Check if we enter a valid subarray
            if not in_valid_subarray and max_ending_here >= self.params.threshold:
                in_valid_subarray = True
                candidate_start = current_start
                candidate_score = max_ending_here
                max_score_seen = max_ending_here
            elif in_valid_subarray:
                # Update the maximum score seen in this subarray
                if max_ending_here > max_score_seen:
                    max_score_seen = max_ending_here
                    candidate_score = max_ending_here
            
            # Check if we should end current subarray
            if in_valid_subarray and (max_ending_here < 0 or 
                                     (max_score_seen - max_ending_here >= self.params.drop_threshold)):
                # Map from scoring array indices to sequence indices
                # scoring_array[i] represents transition from base i to base i+1
                # So subarray from candidate_start to i represents bases candidate_start to i+1
                start_base = candidate_start
                end_base = i
                substring = str(self[start_base:end_base + 1])
                subarrays.append((start_base, end_base, candidate_score, substring))
                
                # Reset for next subarray
                in_valid_subarray = False
                max_ending_here = 0
                current_start = i + 1
                max_score_seen = 0

        # Add final candidate if it exists and is valid
        if in_valid_subarray:
            start_base = candidate_start
            end_base = n  # Last transition index maps to last base
            substring = str(self[start_base:end_base + 1])
            subarrays.append((start_base, end_base, candidate_score, substring))

        return subarrays

    def zdna_calculator_coverage(self) -> np.ndarray:
        """
        Calculate coverage-based scoring (placeholder for future implementation).
        
        Currently falls back to transitions method.
        
        Returns:
            np.ndarray: Scoring array
        """
        return self.zdna_calculator_transitions()