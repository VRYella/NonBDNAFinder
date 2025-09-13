#!/usr/bin/env python3
"""
Test the integration of enhanced A-philic detector with the orchestrator and app
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import tempfile
import json

def test_orchestrator_integration():
    """Test that the enhanced A-philic detector works with the orchestrator"""
    print("Testing orchestrator integration...")
    
    try:
        from orchestrator import run_pipeline
        
        # Create a test FASTA file
        test_sequence = (
            "ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC"
            "ATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT"
            "GGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGTAGGAACCCCTAACCCCTAA"
            "GAAAGAAGAAGAAGAAGAAGAAAGGAAGGAAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGG"
            "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"
            "GAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA"
            "CTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT"
        )
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(">test_sequence\n")
            f.write(test_sequence + "\n")
            fasta_path = f.name
        
        try:
            # Run pipeline with A-philic detector enabled
            output_files = run_pipeline(
                fasta_path=fasta_path,
                output_prefix="/tmp/test_aphilic",
                max_workers=1,
                chunk_size=10000,
                detector_classes=['a_philic']  # Only test A-philic
            )
            
            print(f"Pipeline completed. Output files: {output_files}")
            
            # Check if results were generated
            if 'csv' in output_files and os.path.exists(output_files['csv']):
                import pandas as pd
                df = pd.read_csv(output_files['csv'])
                print(f"Found {len(df)} motifs in CSV output")
                
                # Check for A-philic motifs
                aphilic_motifs = df[df['Class'] == 'a_philic']
                print(f"Found {len(aphilic_motifs)} A-philic motifs")
                
                if len(aphilic_motifs) > 0:
                    print("Sample A-philic motif:")
                    for col in ['Start', 'End', 'Length', 'Normalized_Score']:
                        if col in aphilic_motifs.columns:
                            print(f"  {col}: {aphilic_motifs.iloc[0][col]}")
                
                return len(aphilic_motifs) > 0
            else:
                print("No CSV output file generated")
                return False
                
        finally:
            # Cleanup
            try:
                os.unlink(fasta_path)
                for file_path in output_files.values():
                    if os.path.exists(file_path):
                        os.unlink(file_path)
            except:
                pass
                
    except Exception as e:
        print(f"Error during orchestrator integration test: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_hyperscan_integration():
    """Test that hyperscan integration works"""
    print("\nTesting hyperscan integration...")
    
    try:
        from hyperscan_integration import all_motifs_refactored
        
        test_sequence = (
            "ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC"
            "ATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT"
        )
        
        # Test the hyperscan integration function
        motifs = all_motifs_refactored(
            sequence=test_sequence,
            sequence_name="test_hyperscan",
            nonoverlap=False,
            report_hotspots=True,
            calculate_conservation=False
        )
        
        print(f"Hyperscan integration returned {len(motifs)} motifs")
        
        # Look for A-philic motifs
        aphilic_motifs = [m for m in motifs if m.get('Class') == 'A-philic DNA']
        print(f"Found {len(aphilic_motifs)} A-philic motifs via hyperscan integration")
        
        if aphilic_motifs:
            print("Sample A-philic motif from hyperscan:")
            sample = aphilic_motifs[0]
            for key in ['Start', 'End', 'Length', 'Normalized_Score', 'Class', 'Subclass']:
                if key in sample:
                    print(f"  {key}: {sample[key]}")
                    
        return True
        
    except Exception as e:
        print(f"Error during hyperscan integration test: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("Enhanced A-philic Integration Tests")
    print("=" * 50)
    
    success1 = test_orchestrator_integration()
    success2 = test_hyperscan_integration()
    
    if success1 and success2:
        print("\n✅ All integration tests passed!")
    else:
        print("\n❌ Some integration tests failed!")
        if not success1:
            print("  - Orchestrator integration failed")
        if not success2:
            print("  - Hyperscan integration failed")
        sys.exit(1)