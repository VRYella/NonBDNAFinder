# Enhanced A-philic Detection with Kadane's Algorithm

## Overview

This implementation enhances A-philic DNA detection by integrating the high-performance k-mer detection approach from the provided `detect_Aphilic_hyperscan.py` script with Kadane's algorithm for optimal subarray detection, as used in the existing Z-DNA implementation.

## Key Improvements

### 1. Dynamic Region Detection
- **Before**: Fixed 10-mer sliding windows with manual merging
- **After**: Kadane's algorithm finds optimal subarrays that maximize combined tri/tetra propensity scores
- **Benefit**: More accurate region boundaries and better detection of variable-length A-philic regions

### 2. Enhanced Scoring Integration
- **Scoring Array Creation**: Converts tri/tetra propensity scores into a transition-based scoring array
- **Combined Scoring**: Integrates both tetramer and trimer scores with configurable weights
- **Nucleation Filtering**: Applies nucleation requirements during optimization rather than as post-processing

### 3. Framework Integration
- **Existing Infrastructure**: Enhances the existing `APhilicDetector` class in `motif_detectors.py`
- **Backward Compatibility**: Maintains compatibility with the existing `Candidate` structure
- **Enhanced Metadata**: Provides detailed scoring information for downstream analysis

## Implementation Files

### Core Implementation
- **`motif_detectors.py`**: Enhanced `APhilicDetector` class with Kadane's algorithm
- **`enhanced_aphilic_hs.py`**: Standalone script demonstrating the approach

### Test and Validation Files
- **`test_kadane_aphilic.py`**: Unit test for the enhanced detector
- **`comprehensive_comparison.py`**: Comparison of all three approaches
- **`demo_enhanced_integration.py`**: Framework integration demonstration

### Test Data
- **`test_tetra_propensities.csv`**: Tetramer propensity data from problem statement
- **`test_tri_propensities.csv`**: Trimer propensity data from problem statement
- **`test_sequence.fasta`**: Test sequence for validation

## Algorithm Details

### Scoring Array Creation
```python
def _create_scoring_array_from_propensities(self, seq: str):
    """Create scoring array from tri/tetra propensities for Kadane's algorithm"""
    scoring_array = np.zeros(n - 1, dtype=float)
    
    # Score tetranucleotides
    for i in range(n - 3):
        tetra = seq[i:i+4]
        tetra_score = self.TETRA_LOG2.get(tetra, 0.0)
        if tetra_score > 0:
            # Distribute score across transitions
            for pos in range(i, min(i + 3, len(scoring_array))):
                scoring_array[pos] += tetra_score
    
    # Score trinucleotides with reduced weight
    for i in range(n - 2):
        tri = seq[i:i+3]
        tri_score = self.TRI_LOG2.get(tri, 0.0)
        if tri_score > 0:
            for pos in range(i, min(i + 2, len(scoring_array))):
                scoring_array[pos] += tri_score * 0.5
```

### Kadane's Algorithm Application
```python
def _kadane_aphilic_regions(self, scoring_array, seq: str, 
                           threshold: float = 5.0, drop_threshold: float = 10.0):
    """Apply Kadane's algorithm to find optimal A-philic regions"""
    max_ending_here = scoring_array[0]
    start_idx = end_idx = 0
    current_max = 0
    candidate_region = None
    
    for i in range(1, scoring_array.size):
        num = scoring_array[i]
        if num >= max_ending_here + num:
            start_idx = i
            end_idx = i + 1
            max_ending_here = num
        else:
            max_ending_here += num
            end_idx = i + 1
        
        # Apply threshold and nucleation filtering
        if max_ending_here >= threshold and current_max < max_ending_here:
            # Create candidate region with detailed scoring
            details = self._get_region_details(seq, start_idx, end_idx)
            if self._passes_nucleation_filter(region_seq, details):
                candidate_region = (start_idx, end_idx, max_ending_here, region_seq, details)
                current_max = max_ending_here
```

## Performance Characteristics

### Computational Complexity
- **Time Complexity**: O(n) for Kadane's algorithm + O(n) for k-mer scoring = O(n)
- **Space Complexity**: O(n) for scoring array
- **Improvement**: More efficient than O(n²) sliding window approach

### Detection Quality
- **Better Boundaries**: Dynamic region detection finds optimal start/end points
- **Variable Length**: Supports regions from minimum length to full sequence
- **Combined Optimization**: Maximizes both tetramer and trimer contributions simultaneously

## Usage Examples

### Basic Detection
```python
from motif_detectors import APhilicDetector

detector = APhilicDetector()
candidates = detector.detect(sequence, "seq_name", "contig", 0)
scored_candidates = detector.score(candidates)
```

### Standalone Script
```bash
python3 enhanced_aphilic_hs.py \
    --fasta input.fasta \
    --tetra_table tetra_propensities.csv \
    --tri_table tri_propensities.csv \
    --outdir results \
    --threshold 10.0 \
    --tri_nucleation_min 2
```

## Validation Results

The enhanced implementation successfully:

1. **Identifies optimal regions**: Uses Kadane's algorithm to find regions that maximize combined scores
2. **Maintains nucleation requirements**: Applies trinucleotide nucleation filtering
3. **Integrates with framework**: Works seamlessly with existing motif detection pipeline
4. **Supports both pattern types**: Handles both GC-rich and AT-rich A-philic patterns
5. **Provides enhanced metadata**: Includes detailed scoring information for analysis

## Comparison with Original Methods

| Feature | Original Stringent | Enhanced Standalone | Enhanced Detector |
|---------|-------------------|-------------------|------------------|
| Region Detection | Fixed 10-mer windows | Kadane's algorithm | Kadane's algorithm |
| K-mer Matching | Direct scanning | Hyperscan + fallback | Direct scanning |
| Region Boundaries | Fixed length | Dynamic optimal | Dynamic optimal |
| Framework Integration | None | Standalone | Full integration |
| Nucleation Filtering | Post-processing | During optimization | During optimization |
| Performance | O(n²) sliding | O(n) optimal | O(n) optimal |

The enhanced implementation provides significant improvements in detection accuracy and computational efficiency while maintaining compatibility with the existing framework.