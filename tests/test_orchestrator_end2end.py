"""
End-to-end tests for the orchestrator pipeline.
"""

import pytest
import tempfile
import os
from pathlib import Path
import pandas as pd

from orchestrator import run_pipeline, load_fasta_sequences, candidates_to_dataframe
from motifs.base import Candidate


class TestOrchestratorEndToEnd:
    """Test complete pipeline execution"""
    
    def create_test_fasta(self, sequences):
        """Create temporary FASTA file for testing"""
        tmp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
        
        for i, seq in enumerate(sequences):
            tmp_file.write(f">test_sequence_{i}\n{seq}\n")
        
        tmp_file.close()
        return tmp_file.name
    
    def test_fasta_loading(self):
        """Test FASTA file loading"""
        test_sequences = [
            "GGGTTTGGGTTTGGGTTTGGG",
            "ATCGATCGATCGATCGATCG",
            "AAAAAAAAATGCGTAAAAAAAAAA"
        ]
        
        fasta_file = self.create_test_fasta(test_sequences)
        
        try:
            sequences = load_fasta_sequences(fasta_file)
            
            assert len(sequences) == 3
            for i, (seq_id, sequence) in enumerate(sequences):
                assert seq_id == f"test_sequence_{i}"
                assert sequence == test_sequences[i]
        
        finally:
            os.unlink(fasta_file)
    
    def test_small_pipeline_run(self):
        """Test complete pipeline on small test data"""
        # Create test FASTA with known motif-containing sequences
        test_sequences = [
            "GGGTTTGGGTTTGGGTTTGGGAAATTTCCCAAATTTCCCAAATTTCCC",  # G4 + i-motif
            "AAAAAAAAATGCGTAAAAAAAAAATGCGT",  # Curved DNA
            "AGAGAGAGAGAGAGAGAGAGAGAGAGAGA"   # Potential triplex
        ]
        
        fasta_file = self.create_test_fasta(test_sequences)
        output_prefix = tempfile.mktemp()
        
        try:
            # Run pipeline
            output_files = run_pipeline(
                fasta_path=fasta_file,
                output_prefix=output_prefix,
                max_workers=2,
                chunk_size=1000,
                detector_classes=['g_quadruplex', 'curved_dna', 'triplex']
            )
            
            # Check that output files were created
            assert 'csv' in output_files
            assert os.path.exists(output_files['csv'])
            
            # Load and validate results
            df = pd.read_csv(output_files['csv'])
            
            # Basic validation
            assert len(df) > 0
            assert 'Class' in df.columns
            assert 'Start' in df.columns
            assert 'End' in df.columns
            assert 'Normalized_Score' in df.columns
            
            # Check coordinate validity
            for _, row in df.iterrows():
                assert row['Start'] > 0  # 1-based coordinates
                assert row['End'] >= row['Start']
                assert row['Length'] == row['End'] - row['Start'] + 1
            
            # Check that we found some G4s (most likely)
            classes = df['Class'].unique()
            assert 'g_quadruplex' in classes
            
        finally:
            # Cleanup
            os.unlink(fasta_file)
            for file_path in output_files.values():
                if os.path.exists(file_path):
                    os.unlink(file_path)
    
    def test_candidates_to_dataframe(self):
        """Test conversion of candidates to DataFrame"""
        # Create test candidates
        candidates = [
            Candidate(
                sequence_name="test_seq",
                contig="test_contig",
                class_id=6,
                class_name="g_quadruplex",
                subclass="canonical",
                motif_id=1,
                start=10,
                end=30,
                length=21,
                matched_seq=b"GGGTTTGGGTTTGGGTTTGGG",
                pattern_name="test_pattern",
                raw_score=0.5,
                normalized_score=0.7,
                scoring_method="G4Hunter",
                gc_content=0.6
            ),
            Candidate(
                sequence_name="test_seq",
                contig="test_contig", 
                class_id=1,
                class_name="curved_dna",
                subclass="phased",
                motif_id=2,
                start=50,
                end=70,
                length=21,
                matched_seq=b"AAAAAAAAATGCGTAAAAAAA",
                pattern_name="a_tract",
                raw_score=0.3,
                normalized_score=0.4,
                scoring_method="Curvature",
                gc_content=0.2
            )
        ]
        
        df = candidates_to_dataframe(candidates)
        
        # Check DataFrame structure
        assert len(df) == 2
        assert 'S.No' in df.columns
        assert 'Sequence_Name' in df.columns
        assert 'Class' in df.columns
        assert 'Subclass' in df.columns
        
        # Check data integrity
        assert df.iloc[0]['Class'] == 'g_quadruplex'
        assert df.iloc[1]['Class'] == 'curved_dna'
        assert df.iloc[0]['Start'] == 10
        assert df.iloc[1]['Start'] == 50
    
    def test_empty_results(self):
        """Test handling of sequences with no motifs"""
        # Sequence unlikely to contain any motifs
        test_sequences = ["ATATATATATATATATATATA"]
        
        fasta_file = self.create_test_fasta(test_sequences)
        output_prefix = tempfile.mktemp()
        
        try:
            output_files = run_pipeline(
                fasta_path=fasta_file,
                output_prefix=output_prefix,
                max_workers=1,
                detector_classes=['triplex']  # Very specific class
            )
            
            # Should still create output files
            assert 'csv' in output_files
            assert os.path.exists(output_files['csv'])
            
            # Results may be empty or minimal
            df = pd.read_csv(output_files['csv'])
            # Just ensure it's a valid DataFrame
            assert isinstance(df, pd.DataFrame)
        
        finally:
            os.unlink(fasta_file)
            for file_path in output_files.values():
                if os.path.exists(file_path):
                    os.unlink(file_path)
    
    def test_multiple_sequences(self):
        """Test pipeline with multiple sequences"""
        test_sequences = [
            "GGGTTTGGGTTTGGGTTTGGG",  # G4 sequence
            "ATCGATCGATCGATCGATCG",   # Control sequence
            "CCCAAACCCAAACCCAAACCC",  # i-motif sequence
            "AAAAAAAAATGCGTAAAAAAA"   # Curved DNA
        ]
        
        fasta_file = self.create_test_fasta(test_sequences)
        output_prefix = tempfile.mktemp()
        
        try:
            output_files = run_pipeline(
                fasta_path=fasta_file,
                output_prefix=output_prefix,
                max_workers=2,
                detector_classes=['g_quadruplex', 'i_motif', 'curved_dna']
            )
            
            df = pd.read_csv(output_files['csv'])
            
            # Should have results from multiple sequences
            unique_sequences = df['Sequence_Name'].unique()
            assert len(unique_sequences) >= 2  # At least some sequences with motifs
            
            # Check that sequence names are correct
            for seq_name in unique_sequences:
                assert seq_name.startswith('test_sequence_')
        
        finally:
            os.unlink(fasta_file)
            for file_path in output_files.values():
                if os.path.exists(file_path):
                    os.unlink(file_path)
    
    def test_chunked_sequence_processing(self):
        """Test processing of long sequences with chunking"""
        # Create a longer sequence with embedded motifs
        long_sequence = ("ATCGATCG" * 1000 +  # 8kb of random
                        "GGGTTTGGGTTTGGGTTTGGG" +  # G4 motif
                        "ATCGATCG" * 1000)  # 8kb more
        
        fasta_file = self.create_test_fasta([long_sequence])
        output_prefix = tempfile.mktemp()
        
        try:
            output_files = run_pipeline(
                fasta_path=fasta_file,
                output_prefix=output_prefix,
                max_workers=2,
                chunk_size=5000,  # Force chunking
                detector_classes=['g_quadruplex']
            )
            
            df = pd.read_csv(output_files['csv'])
            
            # Should find the embedded G4
            g4_results = df[df['Class'] == 'g_quadruplex']
            assert len(g4_results) > 0
            
            # Check that coordinates are reasonable
            for _, row in g4_results.iterrows():
                assert row['Start'] > 0
                assert row['End'] <= len(long_sequence)
        
        finally:
            os.unlink(fasta_file)
            for file_path in output_files.values():
                if os.path.exists(file_path):
                    os.unlink(file_path)