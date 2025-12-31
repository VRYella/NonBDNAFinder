#!/usr/bin/env python3
"""
Tests for Job Management and Email Notification
================================================

Tests the job ID generation, result persistence, and optional email notification.
"""

import os
import sys
import json
import shutil

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import job_manager
import email_notifier


def test_job_id_generation():
    """Test job ID generation produces valid, unique IDs"""
    print("\n=== Testing Job ID Generation ===")
    
    # Define valid hex characters
    HEX_CHARS = '0123456789abcdef'
    
    # Generate multiple IDs
    job_ids = [job_manager.generate_job_id() for _ in range(100)]
    
    # Check format
    for job_id in job_ids:
        assert len(job_id) == 10, f"Job ID should be 10 chars, got {len(job_id)}"
        assert all(c in HEX_CHARS for c in job_id), f"Job ID should be hex, got {job_id}"
    
    # Check uniqueness
    assert len(set(job_ids)) == len(job_ids), "All job IDs should be unique"
    
    print(f"✓ Generated {len(job_ids)} unique job IDs")
    print(f"  Example IDs: {job_ids[:3]}")
    return True


def test_result_persistence():
    """Test saving and loading job results"""
    print("\n=== Testing Result Persistence ===")
    
    # Create test data
    job_id = job_manager.generate_job_id()
    test_results = [
        [
            {"Class": "G4", "Start": 10, "End": 30, "Score": 0.95},
            {"Class": "Z_DNA", "Start": 50, "End": 70, "Score": 0.85}
        ],
        [
            {"Class": "Cruciform", "Start": 100, "End": 150, "Score": 0.75}
        ]
    ]
    test_sequences = ["ATCGATCGATCG" * 10, "GCTAGCTAGCTA" * 8]
    test_names = ["Seq1", "Seq2"]
    test_metadata = {"test_param": "test_value"}
    
    # Save results
    print(f"  Saving job {job_id}...")
    success = job_manager.save_job_results(
        job_id, test_results, test_sequences, test_names, test_metadata
    )
    assert success, "Save should succeed"
    
    # Check job exists
    assert job_manager.job_exists(job_id), "Job should exist after save"
    
    # Load results
    print(f"  Loading job {job_id}...")
    loaded = job_manager.load_job_results(job_id)
    assert loaded is not None, "Load should return data"
    
    loaded_results, loaded_seqs, loaded_names, loaded_metadata = loaded
    
    # Verify data
    assert loaded_results == test_results, "Results should match"
    assert loaded_seqs == test_sequences, "Sequences should match"
    assert loaded_names == test_names, "Names should match"
    assert loaded_metadata['test_param'] == test_metadata['test_param'], "Metadata should match"
    assert 'job_id' in loaded_metadata, "Metadata should include job_id"
    assert 'timestamp' in loaded_metadata, "Metadata should include timestamp"
    
    # Get summary
    summary = job_manager.get_job_summary(job_id)
    assert summary is not None, "Summary should be available"
    assert summary['num_sequences'] == 2, "Summary should show 2 sequences"
    assert summary['total_motifs'] == 3, "Summary should show 3 motifs"
    
    # Cleanup
    job_dir = job_manager.get_job_directory(job_id)
    if os.path.exists(job_dir):
        shutil.rmtree(job_dir)
    
    print(f"✓ Job persistence working correctly")
    print(f"  Saved and loaded {len(test_results)} sequence results")
    print(f"  Metadata: {summary['num_sequences']} sequences, {summary['total_motifs']} motifs")
    return True


def test_job_listing():
    """Test listing all available jobs"""
    print("\n=== Testing Job Listing ===")
    
    # Create multiple test jobs
    job_ids = []
    for i in range(3):
        job_id = job_manager.generate_job_id()
        job_ids.append(job_id)
        
        test_results = [[{"Class": f"Test{i}", "Start": i*10, "End": i*10+5}]]
        test_sequences = [f"ATCG" * (i+1)]
        test_names = [f"TestSeq{i}"]
        
        job_manager.save_job_results(job_id, test_results, test_sequences, test_names)
    
    # List all jobs
    all_jobs = job_manager.list_all_jobs()
    
    assert len(all_jobs) >= 3, f"Should have at least 3 jobs, got {len(all_jobs)}"
    
    # Verify all our jobs are in the list
    listed_ids = [job['job_id'] for job in all_jobs]
    for job_id in job_ids:
        assert job_id in listed_ids, f"Job {job_id} should be in list"
    
    # Cleanup
    for job_id in job_ids:
        job_dir = job_manager.get_job_directory(job_id)
        if os.path.exists(job_dir):
            shutil.rmtree(job_dir)
    
    print(f"✓ Job listing working correctly")
    print(f"  Found {len(all_jobs)} total jobs")
    return True


def test_invalid_job_lookup():
    """Test behavior with non-existent job ID"""
    print("\n=== Testing Invalid Job Lookup ===")
    
    fake_job_id = "0000000000"
    
    # Check existence
    assert not job_manager.job_exists(fake_job_id), "Fake job should not exist"
    
    # Try to load
    result = job_manager.load_job_results(fake_job_id)
    assert result is None, "Loading fake job should return None"
    
    # Try to get summary
    summary = job_manager.get_job_summary(fake_job_id)
    assert summary is None, "Summary of fake job should return None"
    
    print(f"✓ Invalid job lookup handled gracefully")
    return True


def test_job_id_validation():
    """Test job ID validation prevents path traversal attacks"""
    print("\n=== Testing Job ID Validation ===")
    
    # Test path traversal attempts
    malicious_ids = [
        "../../../etc/passwd",
        "..%2F..%2F..%2Fetc%2Fpasswd",
        "....//....//....//etc//passwd",
        "/absolute/path",
        "valid12345/../../bad",
        "abc",  # Too short
        "abcdefghijklmnop",  # Too long
        "gggggggggg",  # Invalid chars
    ]
    
    for malicious_id in malicious_ids:
        try:
            job_manager.get_job_directory(malicious_id)
            print(f"✗ Should have rejected: {malicious_id}")
            return False
        except ValueError:
            # Expected - malicious ID rejected
            pass
    
    # Test valid ID works
    valid_id = "a1b2c3d4e5"
    try:
        job_dir = job_manager.get_job_directory(valid_id)
        assert valid_id in job_dir, "Valid ID should be accepted"
    except ValueError:
        print(f"✗ Valid ID should be accepted: {valid_id}")
        return False
    
    print(f"✓ Job ID validation prevents path traversal attacks")
    return True


def test_email_validation():
    """Test email address format validation"""
    print("\n=== Testing Email Validation ===")
    
    valid_emails = [
        "user@example.com",
        "test.user@domain.org",
        "admin@mail.university.edu",
    ]
    
    invalid_emails = [
        "",
        "notanemail",
        "@example.com",
        "user@",
        "user@@example.com",
        "user@noextension",
    ]
    
    for email in valid_emails:
        assert email_notifier.validate_email_format(email), f"{email} should be valid"
    
    for email in invalid_emails:
        assert not email_notifier.validate_email_format(email), f"{email} should be invalid"
    
    print(f"✓ Email validation working correctly")
    print(f"  Validated {len(valid_emails)} valid and {len(invalid_emails)} invalid emails")
    return True


def test_email_notification_without_config():
    """Test that email notification fails gracefully without config"""
    print("\n=== Testing Email Notification (No Config) ===")
    
    # Create a mock secrets object without email config
    class MockSecrets:
        pass
    
    secrets = MockSecrets()
    
    # Try to send notification
    result = email_notifier.send_job_notification(
        to_email="test@example.com",
        job_id="test123456",
        secrets=secrets
    )
    
    # Should fail gracefully (return False, not raise exception)
    assert result is False, "Should return False when no config available"
    
    print(f"✓ Email notification fails gracefully without config")
    return True


def run_all_tests():
    """Run all tests and report results"""
    print("=" * 60)
    print("NBDScanner Job Management & Email Notification Tests")
    print("=" * 60)
    
    tests = [
        ("Job ID Generation", test_job_id_generation),
        ("Result Persistence", test_result_persistence),
        ("Job Listing", test_job_listing),
        ("Invalid Job Lookup", test_invalid_job_lookup),
        ("Job ID Validation (Security)", test_job_id_validation),
        ("Email Validation", test_email_validation),
        ("Email Without Config", test_email_notification_without_config),
    ]
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        try:
            if test_func():
                passed += 1
            else:
                failed += 1
                print(f"✗ {test_name} FAILED")
        except Exception as e:
            failed += 1
            print(f"✗ {test_name} FAILED with exception: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60)
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
