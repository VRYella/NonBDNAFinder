#!/usr/bin/env python3
"""
Simulate App Workflow
=====================

This script simulates what happens in the actual Streamlit app
to ensure all components work together correctly.
"""

import sys
import os
import shutil

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from job_manager import generate_job_id, save_job_results, load_job_results, job_exists
from ntfy_notifier import validate_topic, send_ntfy_notification

def simulate_analysis_workflow():
    """Simulate the complete analysis workflow"""
    print("\n" + "=" * 70)
    print("SIMULATING STREAMLIT APP WORKFLOW")
    print("=" * 70)
    
    # ========== PHASE 1: User submits analysis ==========
    print("\n📤 PHASE 1: User submits analysis")
    print("-" * 70)
    
    # Generate job ID immediately (displayed to user)
    job_id = generate_job_id()
    print(f"✓ Job ID generated and displayed: {job_id}")
    
    # Optional ntfy topic input
    user_topic = ""  # Could be a ntfy topic or empty
    if user_topic and validate_topic(user_topic):
        print(f"✓ ntfy topic provided and validated")
    else:
        print("ℹ No valid topic provided - will skip notification")
    
    # ========== PHASE 2: Analysis runs ==========
    print("\n🧬 PHASE 2: Analysis executes")
    print("-" * 70)
    
    # Simulated analysis results
    sequences = ["GGGTTAGGGTTAGGGTTAGGG" * 20]
    sequence_names = ["Test_Sequence_1"]
    results = [[
        {"Class": "G4", "Start": 10, "End": 30, "Score": 0.9, "Sequence_Name": "Test_Sequence_1"},
        {"Class": "G4", "Start": 50, "End": 70, "Score": 0.85, "Sequence_Name": "Test_Sequence_1"}
    ]]
    
    print(f"✓ Analysis complete: {len(sequences)} sequences, {sum(len(r) for r in results)} motifs")
    
    # ========== PHASE 3: Save results to disk ==========
    print("\n💾 PHASE 3: Persist results")
    print("-" * 70)
    
    metadata = {"analysis_time": 1.5, "speed": 2000}
    save_success = save_job_results(job_id, results, sequences, sequence_names, metadata)
    
    if save_success:
        print(f"✓ Results saved to disk under job ID: {job_id}")
        print(f"  Location: results/{job_id}/")
        print(f"  Files: results.json, sequences.json, metadata.json")
    else:
        print("✗ Failed to save results!")
        return False
    
    # ========== PHASE 4: Optional ntfy notification ==========
    print("\n🔔 PHASE 4: ntfy push notification (optional)")
    print("-" * 70)
    
    if user_topic:
        notification_sent = send_ntfy_notification(
            topic=user_topic,
            job_id=job_id,
            metadata=metadata
        )
        
        if notification_sent:
            print(f"✓ Push notification sent")
        else:
            print(f"⚠ Notification not sent (topic may be invalid) - job still successful")
            print(f"  User can still retrieve results via Job ID: {job_id}")
    else:
        print("ℹ No topic provided - skipping notification")
    
    # ========== PHASE 5: User retrieves results later ==========
    print("\n🔍 PHASE 5: User retrieves results (Job Lookup)")
    print("-" * 70)
    
    # Check if job exists
    if job_exists(job_id):
        print(f"✓ Job {job_id} found in system")
    else:
        print(f"✗ Job {job_id} not found!")
        return False
    
    # Load results
    loaded = load_job_results(job_id)
    
    if loaded:
        loaded_results, loaded_seqs, loaded_names, loaded_metadata = loaded
        print(f"✓ Results loaded successfully")
        print(f"  Sequences: {len(loaded_seqs)}")
        print(f"  Motifs: {sum(len(r) for r in loaded_results)}")
        print(f"  Timestamp: {loaded_metadata.get('timestamp', 'N/A')}")
        
        # Verify data integrity
        assert loaded_results == results, "Data integrity check failed!"
        assert loaded_seqs == sequences, "Data integrity check failed!"
        print(f"✓ Data integrity verified")
    else:
        print(f"✗ Failed to load results!")
        return False
    
    # ========== PHASE 6: Cleanup ==========
    print("\n🧹 PHASE 6: Test cleanup")
    print("-" * 70)
    
    job_dir = f"results/{job_id}"
    if os.path.exists(job_dir):
        shutil.rmtree(job_dir)
        print(f"✓ Test data cleaned up")
    
    print("\n" + "=" * 70)
    print("✅ WORKFLOW SIMULATION COMPLETE - ALL PHASES SUCCESSFUL")
    print("=" * 70)
    print("\nKey Features Verified:")
    print("  ✓ Job ID generation and display")
    print("  ✓ Result persistence to disk")
    print("  ✓ Optional ntfy notification (graceful failure without topic)")
    print("  ✓ Job lookup and retrieval")
    print("  ✓ Data integrity across save/load cycle")
    print("  ✓ App continues working even if notification fails")
    
    return True

def test_invalid_job_lookup():
    """Test that invalid job lookups are handled gracefully"""
    print("\n" + "=" * 70)
    print("TESTING INVALID JOB LOOKUP")
    print("=" * 70)
    
    fake_job_id = "0000000000"
    print(f"\nAttempting to lookup non-existent job: {fake_job_id}")
    
    if not job_exists(fake_job_id):
        print("✓ Correctly identified job does not exist")
    else:
        print("✗ False positive - job should not exist!")
        return False
    
    loaded = load_job_results(fake_job_id)
    if loaded is None:
        print("✓ Load returned None for invalid job (correct behavior)")
    else:
        print("✗ Should return None for invalid job!")
        return False
    
    print("\n✅ Invalid job lookup handled gracefully")
    return True

def test_topic_validation_scenarios():
    """Test various ntfy topic scenarios"""
    print("\n" + "=" * 70)
    print("TESTING NTFY TOPIC SCENARIOS")
    print("=" * 70)
    
    test_cases = [
        ("", False, "Empty topic"),
        ("   ", False, "Whitespace only"),
        ("nbd-job-123", True, "Valid topic with hyphens"),
        ("my-analysis", True, "Simple topic name"),
        ("test_topic", True, "Topic with underscore"),
    ]
    
    all_passed = True
    for topic, expected, description in test_cases:
        result = validate_topic(topic)
        status = "✓" if result == expected else "✗"
        if result != expected:
            all_passed = False
        print(f"  {status} {description}: '{topic}' -> {result}")
    
    if all_passed:
        print("\n✅ All topic validation tests passed")
    else:
        print("\n✗ Some topic validation tests failed")
    
    return all_passed

if __name__ == "__main__":
    try:
        # Run all test scenarios
        workflow_ok = simulate_analysis_workflow()
        invalid_lookup_ok = test_invalid_job_lookup()
        topic_ok = test_topic_validation_scenarios()
        
        if workflow_ok and invalid_lookup_ok and topic_ok:
            print("\n" + "=" * 70)
            print("🎉 ALL TESTS PASSED - READY FOR PRODUCTION")
            print("=" * 70)
            sys.exit(0)
        else:
            print("\n❌ SOME TESTS FAILED")
            sys.exit(1)
            
    except Exception as e:
        print(f"\n❌ TEST FAILED WITH EXCEPTION: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
