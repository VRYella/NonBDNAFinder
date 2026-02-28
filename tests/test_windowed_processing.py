"""Tests for memory-safe windowed processing of large sequences.

Validates that:
1. Sequences <= 10 Mb are processed directly (no windowing).
2. Sequences > 10 Mb are split into 5 Mb windows with 2 kb overlap.
3. Genomic coordinates are correctly offset per window.
4. Overlap-induced duplicate motifs are removed.
5. Multi-chromosome FASTA files process each sequence independently.
6. The genome_worker parameters match the required specification.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest

# ---------------------------------------------------------------------------
# Constants matching the problem-statement specification
# ---------------------------------------------------------------------------
THRESHOLD_BP    = 10_000_000   # 10 Mb
WINDOW_SIZE_BP  = 5_000_000    # 5 Mb
OVERLAP_BP      = 2_000        # 2 kb
STEP_BP         = WINDOW_SIZE_BP - OVERLAP_BP


# ---------------------------------------------------------------------------
# Helpers (inline, no import required)
# ---------------------------------------------------------------------------

def _windowed_scan(seq, detect_fn, threshold=THRESHOLD_BP,
                   window_size=WINDOW_SIZE_BP, overlap=OVERLAP_BP):
    """Windowed scan matching the logic added to the Analysis notebook."""
    if len(seq) <= threshold:
        return detect_fn(seq, 0)

    step = window_size - overlap
    all_motifs = []
    for wstart in range(0, len(seq), step):
        wend = min(wstart + window_size, len(seq))
        wseq = seq[wstart:wend]
        wmotifs = detect_fn(wseq, wstart)
        for m in wmotifs:
            m['Start'] = m.get('Start', 0) + wstart
            m['End']   = m.get('End',   0) + wstart
        all_motifs.extend(wmotifs)

    # Dedup
    seen = {}
    for m in all_motifs:
        k = (m.get('Class', ''), m.get('Subclass', ''), m.get('Start', 0), m.get('End', 0))
        if k not in seen or m.get('Score', 0) > seen[k].get('Score', 0):
            seen[k] = m
    return sorted(seen.values(), key=lambda x: x.get('Start', 0))


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestWindowParameters:
    """Step size and window boundaries match the specification."""

    def test_step_size(self):
        assert STEP_BP == 4_998_000

    def test_window_count_for_exact_multiple(self):
        seq_len = 10_000_000  # exactly at threshold — processed directly (not windowed)
        calls = []

        def mock_detect(subseq, offset):
            calls.append((len(subseq), offset))
            return []

        _windowed_scan('A' * seq_len, mock_detect)
        assert calls == [(seq_len, 0)], "Sequence at threshold must be processed as a single call"

    def test_window_count_just_above_threshold(self):
        seq_len = 10_000_001
        starts = list(range(0, seq_len, STEP_BP))
        assert len(starts) == 3   # windows starting at 0, 4998000, 9996000

    def test_last_window_clipped_to_sequence_length(self):
        seq_len = 12_000_000
        starts = list(range(0, seq_len, STEP_BP))
        ends = [min(s + WINDOW_SIZE_BP, seq_len) for s in starts]
        assert ends[-1] == seq_len


class TestSmallSequenceDirectProcessing:
    """Sequences at or below 10 Mb are not windowed."""

    def test_10mb_processed_directly(self):
        calls = []

        def mock_detect(subseq, offset):
            calls.append((len(subseq), offset))
            return []

        seq = 'A' * THRESHOLD_BP
        _windowed_scan(seq, mock_detect)
        # Single call with the whole sequence and offset 0
        assert calls == [(THRESHOLD_BP, 0)]

    def test_small_seq_processed_directly(self):
        calls = []

        def mock_detect(subseq, offset):
            calls.append(len(subseq))
            return []

        seq = 'A' * 1_000_000   # 1 Mb
        _windowed_scan(seq, mock_detect)
        assert calls == [1_000_000]


class TestLargeSequenceWindowing:
    """Sequences > 10 Mb are split with the correct window / overlap parameters."""

    def test_multiple_windows_created(self):
        seq_len = 15_000_000
        calls = []

        def mock_detect(subseq, offset):
            calls.append((len(subseq), offset))
            return []

        _windowed_scan('A' * seq_len, mock_detect)
        # Expected windows: 0, 4998000, 9996000, 14994000 → 4 calls
        assert len(calls) == len(list(range(0, seq_len, STEP_BP)))

    def test_window_offsets_match_step(self):
        seq_len = 12_000_000
        offsets = []

        def mock_detect(subseq, offset):
            offsets.append(offset)
            return []

        _windowed_scan('A' * seq_len, mock_detect)
        assert offsets[0] == 0
        assert offsets[1] == STEP_BP

    def test_coordinate_offset_applied(self):
        seq_len = 12_000_000

        def mock_detect(subseq, offset):
            # Return one fake motif at local position 100
            return [{'Class': 'Z-DNA', 'Subclass': 'Z-DNA', 'Start': 100, 'End': 200, 'Score': 1.0}]

        motifs = _windowed_scan('A' * seq_len, mock_detect)
        starts = [m['Start'] for m in motifs]
        # First window: offset=0, start should be 100
        assert 100 in starts
        # Second window: offset=STEP_BP, start should be STEP_BP + 100
        assert STEP_BP + 100 in starts

    def test_overlap_duplicates_removed(self):
        """A motif returned by two windows (due to overlap) appears only once."""
        seq_len = 12_000_000
        # Motif that falls in the overlap zone between window 0 and window 1
        overlap_start = STEP_BP - 500   # just before window-1 starts
        overlap_end   = overlap_start + 100

        def mock_detect(subseq, offset):
            local_start = overlap_start - offset
            local_end   = overlap_end   - offset
            if 0 <= local_start < len(subseq):
                return [{'Class': 'Z-DNA', 'Subclass': 'Z-DNA',
                         'Start': local_start, 'End': local_end, 'Score': 1.0}]
            return []

        motifs = _windowed_scan('A' * seq_len, mock_detect)
        # After dedup, the motif with absolute coords should appear exactly once
        matching = [m for m in motifs
                    if m['Start'] == overlap_start and m['End'] == overlap_end]
        assert len(matching) == 1


class TestMultiChromosomeIndependence:
    """Each chromosome is processed independently; large ones are windowed."""

    def _process_fasta(self, seqs):
        """Simulate per-chromosome processing with the windowed strategy."""
        results = {}
        for name, seq in seqs.items():
            results[name] = _windowed_scan(seq, lambda s, off: [])
        return results

    def test_small_chromosomes_processed_directly(self):
        seqs = {'chr1': 'A' * 1_000_000, 'chr2': 'C' * 5_000_000}
        results = self._process_fasta(seqs)
        assert set(results.keys()) == {'chr1', 'chr2'}

    def test_large_and_small_chromosomes_mixed(self):
        seqs = {
            'small': 'A' * 5_000_000,
            'large': 'G' * 15_000_000,
        }
        results = self._process_fasta(seqs)
        assert set(results.keys()) == {'small', 'large'}

    def test_chromosomes_independent(self):
        """Processing one chromosome does not affect another."""
        processed = []

        def mock_detect(subseq, offset):
            processed.append(len(subseq))
            return []

        for name, seq in {'chr1': 'A' * 1_000_000, 'chr2': 'T' * 1_000_000}.items():
            _windowed_scan(seq, mock_detect)

        assert len(processed) == 2


class TestGenomeWorkerParameters:
    """genome_worker module constants match the specification."""

    def test_process_chromosome_importable(self):
        from Utilities.genome_worker import process_chromosome
        assert callable(process_chromosome)

    def test_genome_worker_chunk_function(self):
        """_chunk_sequence with window/overlap params produces correct tiles."""
        from Utilities.genome_worker import _chunk_sequence
        tiles = _chunk_sequence('A' * 12_000_000, WINDOW_SIZE_BP, OVERLAP_BP)
        assert tiles[0] == (0, WINDOW_SIZE_BP)
        assert tiles[1][0] == STEP_BP


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
