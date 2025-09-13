"""
Constants and parameters for Z-DNA detection using Kadane algorithm approach.

This module defines the configuration parameters used for scoring dinucleotide
transitions and finding Z-DNA forming regions using a modified Kadane algorithm.
"""

from typing import Tuple
from attrs import define, field
import multiprocessing as mp


@define(kw_only=True, slots=True, frozen=True)
class Params:
    """
    Parameters for Z-DNA detection algorithm.
    
    This class contains all the configuration parameters needed for:
    - Scoring dinucleotide transitions
    - Applying mismatch penalties
    - Finding subarrays using Kadane-like algorithm
    - Parallel processing and output formatting
    """
    
    # Transition weights
    GC_weight: float = field(default=3.0, converter=float)
    AT_weight: float = field(default=1.0, converter=float) 
    GT_weight: float = field(default=2.0, converter=float)
    AC_weight: float = field(default=2.0, converter=float)
    
    # Additional scoring parameters
    cadence_reward: float = field(default=0.2, converter=float)
    
    # Mismatch penalty parameters
    mismatch_penalty_starting_value: int = field(default=3, converter=int)
    mismatch_penalty_linear_delta: int = field(default=1, converter=int)
    mismatch_penalty_type: str = field(default="linear")
    mismatch_penalty_choices: Tuple[str, ...] = field(default=("linear", "exponential"), init=False)
    
    # Processing parameters
    n_jobs: int = field(default=-1, converter=int)
    threshold: float = field(default=5.0, converter=float)
    drop_threshold: float = field(default=50.0, converter=float)
    
    # Consecutive AT scoring - penalties/bonuses for consecutive AT/TA transitions
    consecutive_AT_scoring: Tuple[float, ...] = field(default=(3.0, 1.5, 0.7))
    
    # Output parameters
    display_sequence_score: int = field(default=0, converter=int)
    output_dir: str = field(default="zdna_extractions", converter=str)
    total_sequence_scoring: bool = field(default=False)
    
    # Column headers for output DataFrame
    headers: Tuple[str, ...] = field(
        default=("Chromosome", "Start", "End", "Z-DNA Score", "Sequence", "totalSequenceScore"),
        init=False
    )
    
    def __attrs_post_init__(self):
        """Validate parameters after initialization."""
        if self.mismatch_penalty_type not in self.mismatch_penalty_choices:
            raise ValueError(
                f"mismatch_penalty_type must be one of {self.mismatch_penalty_choices}, "
                f"got {self.mismatch_penalty_type}"
            )
        
        # Handle n_jobs = -1 (use all cores)
        if self.n_jobs == -1:
            object.__setattr__(self, 'n_jobs', mp.cpu_count())
        elif self.n_jobs <= 0:
            object.__setattr__(self, 'n_jobs', 1)
    
    @property 
    def __new_dict__(self):
        """Dictionary representation for CLI display."""
        return {
            field.name: getattr(self, field.name) 
            for field in self.__class__.__attrs_attrs__ 
            if field.init  # Only include fields that can be set during init
        }