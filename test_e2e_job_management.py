#!/usr/bin/env python3
"""
End-to-End Integration Test
============================

Tests the complete job management flow:
1. Generate job ID
2. Simulate analysis results
3. Save to disk
4. Load from disk
5. Verify data integrity
"""

import sys
import os
import shutil

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from job_manager import (
    generate_job_id, save_job_results, load_job_results,
    job_exists, get_job_summary
)

def test_end_to_end_workflow():
    """Test complete workflow as it would happen in the app"""
    print("=" * 60)
    print("End-to-End Job Management Test")
    print("=" * 60)
    
    # Step 1: Generate job ID (happens when analysis starts)
    print("\n1. Generating Job ID...")
    job_id = generate_job_id()
    print(f"   ✓ Job ID: {job_id}")
    
    # Step 2: Simulate analysis results
    print("\n2. Simulating analysis...")
    sequences = [
        "GGGTTAGGGTTAGGGTTAGGG" * 10,  # G4-rich sequence
        "ATCGATCGATCGATCGATCG" * 15    # Simple repeat
    ]
    sequence_names = ["Test_G4_Sequence", "Test_Simple_Repeat"]
    
    results = [
        [
            {
                "Class": "G4",
                "Subclass": "G4_Intramolecular",
                "Start": 10,
                "End": 30,
                "Length": 20,
                "Score": 0.95,
                "Sequence_Name": "Test_G4_Sequence"
            },
            {
                "Class": "G4",
                "Subclass": "G4_Parallel",
                "Start": 50,
                "End": 72,
                "Length": 22,
                "Score": 0.88,
                "Sequence_Name": "Test_G4_Sequence"
            }
        ],
        [
            {
                "Class": "STR",
                "Subclass": "Dinucleotide_Repeat",
                "Start": 20,
                "End": 60,
                "Length": 40,
                "Score": 0.75,
                "Sequence_Name": "Test_Simple_Repeat"
            }
        ]
    ]
    
    metadata = {
        "analysis_time": 2.5,
        "speed_bp_per_sec": 1200,
        "test_mode": True
    }
    
    print(f"   ✓ Generated {len(sequences)} sequences")
    print(f"   ✓ Generated {sum(len(r) for r in results)} motifs")
    
    # Step 3: Save results (happens after analysis completes)
    print("\n3. Saving results to disk...")
    save_success = save_job_results(job_id, results, sequences, sequence_names, metadata)
    assert save_success, "Save should succeed"
    print(f"   ✓ Results saved to results/{job_id}/")
    
    # Step 4: Verify job exists
    print("\n4. Verifying job exists...")
    exists = job_exists(job_id)
    assert exists, "Job should exist after save"
    print(f"   ✓ Job {job_id} exists")
    
    # Step 5: Get job summary
    print("\n5. Retrieving job summary...")
    summary = get_job_summary(job_id)
    assert summary is not None, "Summary should be available"
    print(f"   ✓ Summary retrieved:")
    print(f"      - Sequences: {summary['num_sequences']}")
    print(f"      - Motifs: {summary['total_motifs']}")
    print(f"      - Total BP: {summary['total_bp']:,}")
    print(f"      - Timestamp: {summary['timestamp']}")
    
    # Step 6: Load complete results (happens on job lookup)
    print("\n6. Loading complete results...")
    loaded = load_job_results(job_id)
    assert loaded is not None, "Load should return data"
    
    loaded_results, loaded_seqs, loaded_names, loaded_metadata = loaded
    
    # Step 7: Verify data integrity
    print("\n7. Verifying data integrity...")
    assert loaded_results == results, "Results should match"
    assert loaded_seqs == sequences, "Sequences should match"
    assert loaded_names == sequence_names, "Names should match"
    assert loaded_metadata['test_mode'] == True, "Metadata should match"
    print("   ✓ All data verified - perfect match!")
    
    # Cleanup
    print("\n8. Cleaning up test data...")
    job_dir = f"results/{job_id}"
    if os.path.exists(job_dir):
        shutil.rmtree(job_dir)
    print("   ✓ Test data removed")
    
    print("\n" + "=" * 60)
    print("✅ END-TO-END TEST PASSED")
    print("=" * 60)
    return True

if __name__ == "__main__":
    try:
        success = test_end_to_end_workflow()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
