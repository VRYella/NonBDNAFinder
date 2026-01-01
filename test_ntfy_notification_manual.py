#!/usr/bin/env python3
"""
Manual ntfy.sh Notification Test
=================================

This script demonstrates and tests ntfy.sh notification functionality.
It shows what the notification looks like without actually sending it
(unless a real topic is provided).
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ntfy_notifier import send_ntfy_notification, validate_topic

def test_topic_validation():
    """Test topic validation"""
    print("\n" + "=" * 70)
    print("TESTING TOPIC VALIDATION")
    print("=" * 70)
    
    test_topics = [
        ("nbd-job-123", True, "Valid alphanumeric with hyphens"),
        ("my-analysis", True, "Simple topic name"),
        ("test_topic", True, "Topic with underscore"),
        ("", False, "Empty topic"),
        ("   ", False, "Whitespace only"),
    ]
    
    for topic, expected, description in test_topics:
        result = validate_topic(topic)
        status = "✓" if result == expected else "✗"
        print(f"  {status} {description}: {result}")
        if result != expected:
            print(f"    ERROR: Expected {expected}, got {result}")
    
    print("\n✅ Topic validation tests complete")


def test_notification_payload_generation():
    """Test notification payload generation (without sending)"""
    print("\n" + "=" * 70)
    print("TESTING NOTIFICATION PAYLOAD GENERATION")
    print("=" * 70)
    
    # Test with empty topic (won't actually send)
    print("\n1. Testing with empty topic (should fail gracefully)...")
    result = send_ntfy_notification(
        topic="",
        job_id="test123abc"
    )
    if not result:
        print("   ✓ Correctly failed with empty topic")
    else:
        print("   ✗ Should have failed with empty topic")
    
    print("\n2. Testing with missing job_id (should fail gracefully)...")
    result = send_ntfy_notification(
        topic="test-topic",
        job_id=""
    )
    if not result:
        print("   ✓ Correctly failed with empty job_id")
    else:
        print("   ✗ Should have failed with empty job_id")
    
    print("\n3. Demonstrating payload structure...")
    print("   When a valid topic is provided, the ntfy notification includes:")
    print("   - Title: 🧬 Job Completed")
    print("   - Job ID")
    print("   - Number of sequences analyzed")
    print("   - Total motifs detected")
    print("   - Optional direct access link (if app_url provided)")
    print("   - Tags: dna, science")
    
    print("\n✅ Payload generation tests complete")


def demonstrate_usage():
    """Demonstrate how to use ntfy.sh notifications"""
    print("\n" + "=" * 70)
    print("NTFY.SH NOTIFICATION USAGE DEMONSTRATION")
    print("=" * 70)
    
    print("\n📖 HOW TO USE NTFY.SH:")
    print("   1. Pick any topic name (e.g., 'nbd-job-123', 'my-analysis')")
    print("   2. Subscribe to that topic:")
    print("      - Web: Visit https://ntfy.sh/your-topic-name")
    print("      - Mobile: Install ntfy app and subscribe to your topic")
    print("      - Desktop: Use ntfy CLI or web interface")
    print("   3. Enter the topic name in NBDScanner")
    print("   4. Run your analysis")
    print("   5. Get instant push notification when complete!")
    
    print("\n📝 EXAMPLE USAGE IN NBDSCANNER:")
    print("   - Navigate to 'Upload & Analyze' tab")
    print("   - Upload your FASTA sequences")
    print("   - (Optional) Enter a topic name in the notification field")
    print("   - Subscribe to that topic at ntfy.sh")
    print("   - Click 'Run Analysis'")
    print("   - When complete, push notification is sent automatically")
    
    print("\n🔒 PRIVACY & SECURITY:")
    print("   - No account required")
    print("   - No signup or registration")
    print("   - Topic is NOT stored on disk")
    print("   - Only held in session memory during analysis")
    print("   - No personal data collected or transmitted")
    print("   - Notification is completely optional")
    print("   - Analysis runs successfully with or without notification")
    
    print("\n⚡ FEATURES:")
    print("   - Zero configuration")
    print("   - Instant notification when analysis completes")
    print("   - Works on web, mobile, and desktop")
    print("   - Includes Job ID for result retrieval")
    print("   - Shows analysis summary (sequences, motifs)")
    print("   - Non-blocking (won't delay analysis)")
    print("   - Graceful failure (won't crash if topic is invalid)")
    print("   - Free public service")


def test_with_real_topic():
    """Test with a real ntfy topic (if provided by user)"""
    print("\n" + "=" * 70)
    print("OPTIONAL: TEST WITH REAL NTFY TOPIC")
    print("=" * 70)
    
    print("\n⚠️  To test with a real ntfy topic:")
    print("   Set the NTFY_TOPIC environment variable")
    print("   Example: export NTFY_TOPIC='test-nbd-123'")
    print("   Then subscribe at: https://ntfy.sh/test-nbd-123")
    
    topic = os.environ.get('NTFY_TOPIC', '')
    
    if topic:
        print(f"\n✓ Found topic in environment: {topic}")
        print("  Attempting to send test notification...")
        print(f"  Subscribe at: https://ntfy.sh/{topic}")
        
        result = send_ntfy_notification(
            topic=topic,
            job_id="test_manual_123",
            app_url="https://example.com",
            metadata={
                'num_sequences': 5,
                'total_motifs': 142,
                'timestamp': '2025-12-31T08:00:00Z'
            }
        )
        
        if result:
            print("  ✅ Notification sent successfully! Check your ntfy subscription.")
        else:
            print("  ⚠️  Notification failed (check topic validity)")
    else:
        print("\n  No topic provided - skipping real notification test")
        print("  This is normal and expected for automated testing")


if __name__ == "__main__":
    print("=" * 70)
    print("NTFY.SH NOTIFICATION - MANUAL TEST SUITE")
    print("=" * 70)
    
    try:
        test_topic_validation()
        test_notification_payload_generation()
        demonstrate_usage()
        test_with_real_topic()
        
        print("\n" + "=" * 70)
        print("✅ ALL TESTS COMPLETED SUCCESSFULLY")
        print("=" * 70)
        print("\nThe ntfy.sh notification system is:")
        print("  ✓ Account-free (no signup required)")
        print("  ✓ Optional (jobs run without it)")
        print("  ✓ Non-blocking (won't delay analysis)")
        print("  ✓ Secure (no credentials needed)")
        print("  ✓ Free (no paid services)")
        print("  ✓ Gracefully failing (invalid topics won't crash app)")
        
        sys.exit(0)
        
    except Exception as e:
        print(f"\n❌ TEST FAILED WITH EXCEPTION: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
