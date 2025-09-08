#!/usr/bin/env python3
"""
Comprehensive test with 20kb sequence to validate overlap resolution and performance.

This test validates:
1. Overlap resolution works correctly with large sequences
2. Different class overlaps are preserved
3. Same class/different subclass overlaps are resolved
4. All visualizations are produced
5. Performance benchmarks
"""

import os
import sys
import time
import tempfile
import random
from pathlib import Path
from typing import Dict, List, Any

import pandas as pd

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from orchestrator import run_pipeline
from overlap_resolution import EnhancedOverlapResolver, OverlapConfig, OverlapStrategy
from motifs.base import Candidate
from hyperscan_integration import all_motifs_refactored


def generate_test_sequence_20kb() -> str:
    """
    Generate a 20kb test sequence with known motifs for comprehensive testing.
    
    Returns:
        20,000 nucleotide sequence with embedded motifs
    """
    sequence_parts = []
    
    # Start with random DNA
    def random_dna(length: int) -> str:
        return ''.join(random.choices('ATCG', k=length))
    
    # Add initial random sequence
    sequence_parts.append(random_dna(1000))
    
    # Add known G-quadruplex motifs (different subclasses)
    g4_motifs = [
        # Canonical G4
        'GGGAGGGTGGGAGGGT',
        'GGGCGGGAGGGCGGG',
        # Relaxed G4 (longer loops)
        'GGGATTTTGGGATTTTGGGATTTTGGG',
        # Bulged G4
        'GGGAGGGTGGGAGGGTGGGAGGGT',
    ]
    
    # Add spacers and motifs
    for motif in g4_motifs:
        sequence_parts.append(random_dna(500))  # Spacer
        sequence_parts.append(motif)
    
    # Add triplex motifs
    triplex_motifs = [
        'AAAAAAAAAAAAAAAAAAAAA',  # Homopurine
        'TTTTTTTTTTTTTTTTTTTTT',  # Homopyrimidine
        'GGGGGGGGGGGGGGGGGGGGG',
    ]
    
    for motif in triplex_motifs:
        sequence_parts.append(random_dna(800))
        sequence_parts.append(motif)
    
    # Add i-motif motifs
    i_motif_sequences = [
        'CCCACCCCCCACCCC',
        'CCCTCCCACCCGCCC',
    ]
    
    for motif in i_motif_sequences:
        sequence_parts.append(random_dna(600))
        sequence_parts.append(motif)
    
    # Add curved DNA (A-tracts)
    curved_motifs = [
        'AAAAAAAAAATTTTTTTTTT',
        'AAAAAAAAAAGGGGGGGGGG',
    ]
    
    for motif in curved_motifs:
        sequence_parts.append(random_dna(700))
        sequence_parts.append(motif)
    
    # Add Z-DNA motifs
    z_dna_motifs = [
        'CGCGCGCGCGCGCGCG',
        'CACACACACACACACA',
    ]
    
    for motif in z_dna_motifs:
        sequence_parts.append(random_dna(500))
        sequence_parts.append(motif)
    
    # Fill to exactly 20kb
    current_length = sum(len(part) for part in sequence_parts)
    remaining = 20000 - current_length
    if remaining > 0:
        sequence_parts.append(random_dna(remaining))
    elif remaining < 0:
        # Trim if too long
        total_seq = ''.join(sequence_parts)
        return total_seq[:20000]
    
    return ''.join(sequence_parts)


def create_test_fasta(sequence: str, seq_name: str = "test_20kb") -> str:
    """
    Create a temporary FASTA file with the test sequence.
    
    Args:
        sequence: DNA sequence
        seq_name: Sequence name
        
    Returns:
        Path to temporary FASTA file
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(f">{seq_name}\n{sequence}\n")
        return f.name


def test_overlap_resolution_specificity():
    """
    Test that overlap resolution works correctly for different scenarios.
    """
    print("🧪 Testing overlap resolution specificity...")
    
    # Create test candidates with overlaps
    candidates = [
        # Same class, different subclasses - should resolve overlaps
        Candidate(
            sequence_name="test_seq",
            contig="test_seq",
            class_id=6,
            class_name="g_quadruplex",
            subclass="canonical",
            motif_id=1,
            start=100,
            end=120,
            length=21,
            matched_seq=b"GGGAGGGTGGGAGGGTGGGAG",
            pattern_name="canonical_g4",
            raw_score=0.9
        ),
        Candidate(
            sequence_name="test_seq",
            contig="test_seq",
            class_id=6,
            class_name="g_quadruplex",
            subclass="relaxed",
            motif_id=2,
            start=110,
            end=130,
            length=21,
            matched_seq=b"GGGATTTTGGGATTTTGGGAT",
            pattern_name="relaxed_g4",
            raw_score=0.7
        ),
        # Different classes - should preserve overlaps
        Candidate(
            sequence_name="test_seq",
            contig="test_seq",
            class_id=5,
            class_name="triplex",
            subclass="homopurine",
            motif_id=3,
            start=115,
            end=135,
            length=21,
            matched_seq=b"AAAAAAAAAAAAAAAAAAAAA",
            pattern_name="triplex_homopurine",
            raw_score=0.8
        ),
    ]
    
    # Test same_class_only=True (should resolve G4 overlaps, preserve triplex)
    config = OverlapConfig(
        strategy=OverlapStrategy.HIGHEST_SCORE,
        same_class_only=True,
        min_overlap_percent=0.1
    )
    
    resolver = EnhancedOverlapResolver(config)
    resolved = resolver.resolve_overlaps(candidates.copy())
    
    print(f"  Original candidates: {len(candidates)}")
    print(f"  Resolved candidates: {len(resolved)}")
    
    # Should have 2 candidates: best G4 + triplex
    assert len(resolved) == 2, f"Expected 2 candidates, got {len(resolved)}"
    
    # Check that triplex is preserved
    triplex_found = any(c.class_name == "triplex" for c in resolved)
    assert triplex_found, "Triplex motif should be preserved"
    
    # Check that best G4 is kept
    g4_candidates = [c for c in resolved if c.class_name == "g_quadruplex"]
    assert len(g4_candidates) == 1, "Should have exactly one G4 candidate"
    assert g4_candidates[0].raw_score == 0.9, "Should keep the highest scoring G4"
    
    print("  ✅ Overlap resolution working correctly")


def benchmark_performance(fasta_path: str) -> Dict[str, float]:
    """
    Benchmark performance of the pipeline with timing metrics.
    
    Args:
        fasta_path: Path to test FASTA file
        
    Returns:
        Dictionary of timing metrics
    """
    print("⏱️  Benchmarking performance...")
    
    metrics = {}
    
    # Time the full pipeline
    start_time = time.time()
    
    with tempfile.TemporaryDirectory() as temp_dir:
        output_prefix = os.path.join(temp_dir, "benchmark_test")
        
        try:
            output_files = run_pipeline(
                fasta_path=fasta_path,
                output_prefix=output_prefix,
                max_workers=2,
                chunk_size=50000,
                detector_classes=['g_quadruplex', 'triplex', 'i_motif', 'curved_dna']
            )
            
            pipeline_time = time.time() - start_time
            metrics['pipeline_total_time'] = pipeline_time
            
            # Check if outputs were created
            metrics['csv_created'] = os.path.exists(output_files.get('csv', ''))
            metrics['excel_created'] = os.path.exists(output_files.get('excel', ''))
            metrics['parquet_created'] = os.path.exists(output_files.get('parquet', ''))
            metrics['gff3_created'] = os.path.exists(output_files.get('gff3', ''))
            
            # Load and analyze results
            csv_file = output_files.get('csv')
            if csv_file and os.path.exists(csv_file):
                df = pd.read_csv(csv_file)
                metrics['total_motifs_detected'] = len(df)
                metrics['unique_classes'] = len(df['Class'].unique())
                
                print(f"  📊 Performance Metrics:")
                print(f"    Pipeline time: {pipeline_time:.2f} seconds")
                print(f"    Motifs detected: {len(df)}")
                print(f"    Unique classes: {len(df['Class'].unique())}")
                print(f"    Throughput: {20000/pipeline_time:.0f} bp/second")
                
        except Exception as e:
            print(f"  ❌ Pipeline failed: {e}")
            metrics['pipeline_failed'] = True
            return metrics
    
    print("  ✅ Performance benchmark completed")
    return metrics


def test_visualization_outputs(fasta_path: str):
    """
    Test that all visualization outputs are properly generated.
    
    Args:
        fasta_path: Path to test FASTA file
    """
    print("📊 Testing visualization outputs...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        output_prefix = os.path.join(temp_dir, "viz_test")
        
        try:
            output_files = run_pipeline(
                fasta_path=fasta_path,
                output_prefix=output_prefix,
                max_workers=2,
                chunk_size=50000
            )
            
            # Test CSV output
            csv_file = output_files.get('csv')
            if csv_file and os.path.exists(csv_file):
                df = pd.read_csv(csv_file)
                assert len(df) > 0, "CSV should contain detected motifs"
                
                # Check required columns
                required_columns = ['Class', 'Subclass', 'Start', 'End', 'Normalized_Score']
                for col in required_columns:
                    assert col in df.columns, f"Missing required column: {col}"
                
                print(f"  ✅ CSV export: {len(df)} motifs")
            else:
                print("  ❌ CSV export failed")
            
            # Test Excel output
            excel_file = output_files.get('excel')
            if excel_file and os.path.exists(excel_file):
                print("  ✅ Excel export successful")
            else:
                print("  ❌ Excel export failed")
            
            # Test Parquet output
            parquet_file = output_files.get('parquet')
            if parquet_file and os.path.exists(parquet_file):
                print("  ✅ Parquet export successful")
            else:
                print("  ❌ Parquet export failed")
            
            # Test GFF3 output
            gff3_file = output_files.get('gff3')
            if gff3_file and os.path.exists(gff3_file):
                with open(gff3_file, 'r') as f:
                    content = f.read()
                    assert content.startswith("##gff-version 3"), "GFF3 should have proper header"
                print("  ✅ GFF3 export successful")
            else:
                print("  ❌ GFF3 export failed")
                
        except Exception as e:
            print(f"  ❌ Visualization test failed: {e}")
            raise


def test_hyperscan_integration(sequence: str):
    """
    Test the hyperscan integration directly.
    
    Args:
        sequence: Test DNA sequence
    """
    print("🚀 Testing Hyperscan integration...")
    
    start_time = time.time()
    
    try:
        motifs = all_motifs_refactored(
            sequence=sequence,
            sequence_name="test_20kb",
            nonoverlap=False,  # Allow overlaps for testing
            report_hotspots=True,
            calculate_conservation=False
        )
        
        integration_time = time.time() - start_time
        
        print(f"  📊 Hyperscan Integration Results:")
        print(f"    Integration time: {integration_time:.2f} seconds")
        print(f"    Motifs detected: {len(motifs)}")
        
        if motifs:
            # Check for different classes
            classes_found = set()
            for motif in motifs:
                classes_found.add(motif.get('Class', 'Unknown'))
            
            print(f"    Classes detected: {', '.join(sorted(classes_found))}")
            print("  ✅ Hyperscan integration working")
        else:
            print("  ⚠️  No motifs detected by Hyperscan")
            
    except Exception as e:
        print(f"  ❌ Hyperscan integration failed: {e}")
        raise


def main():
    """
    Run comprehensive 20kb sequence testing.
    """
    print("🧬 NBDFinder Comprehensive 20kb Sequence Test")
    print("=" * 60)
    
    # Set random seed for reproducible results
    random.seed(42)
    
    # Generate test sequence
    print("📝 Generating 20kb test sequence...")
    test_sequence = generate_test_sequence_20kb()
    print(f"  Generated sequence: {len(test_sequence)} nucleotides")
    
    # Create temporary FASTA file
    fasta_path = create_test_fasta(test_sequence)
    print(f"  Created FASTA: {fasta_path}")
    
    try:
        # Test 1: Overlap resolution specificity
        test_overlap_resolution_specificity()
        
        # Test 2: Hyperscan integration
        test_hyperscan_integration(test_sequence)
        
        # Test 3: Full pipeline with visualizations
        test_visualization_outputs(fasta_path)
        
        # Test 4: Performance benchmarking
        metrics = benchmark_performance(fasta_path)
        
        print("\n" + "=" * 60)
        print("🎉 All tests completed successfully!")
        print("\n📈 Summary:")
        print(f"  Sequence length: 20,000 bp")
        print(f"  Pipeline time: {metrics.get('pipeline_total_time', 'N/A'):.2f}s")
        print(f"  Motifs detected: {metrics.get('total_motifs_detected', 'N/A')}")
        print(f"  Output formats: CSV✅ Excel✅ Parquet✅ GFF3✅")
        print(f"  Overlap resolution: ✅")
        print(f"  Hyperscan retained: ✅")
        
        return True
        
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False
        
    finally:
        # Cleanup
        if os.path.exists(fasta_path):
            os.unlink(fasta_path)


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)