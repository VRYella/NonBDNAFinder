"""
Test configuration and shared fixtures.
"""

import pytest
import tempfile
import shutil
import os


@pytest.fixture(scope="session")
def temp_dir():
    """Create temporary directory for test files"""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)


@pytest.fixture
def sample_g4_sequence():
    """Sample G4-forming sequence"""
    return "GGGTTTGGGTTTGGGTTTGGG"


@pytest.fixture
def sample_curved_sequence():
    """Sample curved DNA sequence"""
    return "AAAAAAAAATGCGTAAAAAAAAAATGCGTAAAAAAAAAA"


@pytest.fixture  
def sample_triplex_sequence():
    """Sample triplex-forming sequence"""
    return "AAAAAGGGGGGAAAAGGGGGGAAAAGGGGGG"


@pytest.fixture
def mixed_motif_sequence():
    """Sequence containing multiple motif types"""
    return ("GGGTTTGGGTTTGGGTTTGGG" +  # G4
            "AAAAAAAAAA" +              # A-tract  
            "CCCAAACCCAAACCCAAACCC")   # i-motif


@pytest.fixture
def create_test_fasta():
    """Factory fixture for creating test FASTA files"""
    created_files = []
    
    def _create_fasta(sequences, seq_names=None):
        tmp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
        
        if seq_names is None:
            seq_names = [f"seq_{i}" for i in range(len(sequences))]
        
        for name, seq in zip(seq_names, sequences):
            tmp_file.write(f">{name}\n{seq}\n")
        
        tmp_file.close()
        created_files.append(tmp_file.name)
        return tmp_file.name
    
    yield _create_fasta
    
    # Cleanup
    for file_path in created_files:
        if os.path.exists(file_path):
            os.unlink(file_path)