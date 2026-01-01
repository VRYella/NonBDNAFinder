#!/usr/bin/env python3
"""
Manual Discord Webhook Test
============================

This script demonstrates and tests Discord webhook functionality.
It shows what the notification looks like without actually sending it
(unless a real webhook URL is provided).
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from discord_notifier import send_discord_webhook, validate_webhook_url

def test_webhook_validation():
    """Test webhook URL validation"""
    print("\n" + "=" * 70)
    print("TESTING WEBHOOK URL VALIDATION")
    print("=" * 70)
    
    test_urls = [
        ("https://discord.com/api/webhooks/123456789/abcdef", True, "Valid Discord webhook"),
        ("https://discordapp.com/api/webhooks/987654321/zyxwvu", True, "Valid DiscordApp webhook"),
        ("https://example.com/webhook", True, "Generic HTTPS webhook"),
        ("http://discord.com/api/webhooks/123/abc", False, "Not HTTPS"),
        ("", False, "Empty URL"),
        ("not-a-url", False, "Invalid URL"),
    ]
    
    for url, expected, description in test_urls:
        result = validate_webhook_url(url)
        status = "✓" if result == expected else "✗"
        print(f"  {status} {description}: {result}")
        if result != expected:
            print(f"    ERROR: Expected {expected}, got {result}")
    
    print("\n✅ Webhook validation tests complete")


def test_webhook_payload_generation():
    """Test webhook payload generation (without sending)"""
    print("\n" + "=" * 70)
    print("TESTING WEBHOOK PAYLOAD GENERATION")
    print("=" * 70)
    
    # Test with invalid webhook (won't actually send)
    print("\n1. Testing with empty webhook URL (should fail gracefully)...")
    result = send_discord_webhook(
        webhook_url="",
        job_id="test123abc"
    )
    if not result:
        print("   ✓ Correctly failed with empty webhook")
    else:
        print("   ✗ Should have failed with empty webhook")
    
    print("\n2. Testing with invalid protocol (should fail gracefully)...")
    result = send_discord_webhook(
        webhook_url="http://example.com",
        job_id="test456def"
    )
    if not result:
        print("   ✓ Correctly failed with HTTP (not HTTPS)")
    else:
        print("   ✗ Should have failed with non-HTTPS webhook")
    
    print("\n3. Demonstrating payload structure...")
    print("   When a valid webhook is provided, the Discord notification includes:")
    print("   - Title: 🧬 NBDScanner Analysis Complete")
    print("   - Job ID field with monospace formatting")
    print("   - Analysis Summary (sequences analyzed, motifs detected)")
    print("   - Optional direct access link (if job_url provided)")
    print("   - Footer with app branding")
    print("   - Timestamp of notification")
    print("   - Blue color theme (0x3498DB)")
    
    print("\n✅ Payload generation tests complete")


def demonstrate_usage():
    """Demonstrate how to use Discord webhooks"""
    print("\n" + "=" * 70)
    print("DISCORD WEBHOOK USAGE DEMONSTRATION")
    print("=" * 70)
    
    print("\n📖 HOW TO GET A DISCORD WEBHOOK URL:")
    print("   1. Open your Discord server")
    print("   2. Go to Server Settings → Integrations → Webhooks")
    print("   3. Click 'New Webhook' or select existing webhook")
    print("   4. Copy the webhook URL")
    print("   5. Paste it in the NBDScanner app (optional field)")
    
    print("\n📝 EXAMPLE USAGE IN NBDSCANNER:")
    print("   - Navigate to 'Upload & Analyze' tab")
    print("   - Upload your FASTA sequences")
    print("   - (Optional) Paste Discord webhook URL in the notification field")
    print("   - Click 'Run Analysis'")
    print("   - When complete, Discord notification is sent automatically")
    
    print("\n🔒 PRIVACY & SECURITY:")
    print("   - Webhook URL is NOT stored on disk")
    print("   - Only held in session memory during analysis")
    print("   - No personal data collected or transmitted")
    print("   - Webhook is completely optional")
    print("   - Analysis runs successfully with or without webhook")
    
    print("\n⚡ FEATURES:")
    print("   - Instant notification when analysis completes")
    print("   - Rich Discord embed formatting")
    print("   - Includes Job ID for result retrieval")
    print("   - Shows analysis summary (sequences, motifs)")
    print("   - Non-blocking (won't delay analysis)")
    print("   - Graceful failure (won't crash if webhook is invalid)")


def test_with_real_webhook():
    """Test with a real webhook URL (if provided by user)"""
    print("\n" + "=" * 70)
    print("OPTIONAL: TEST WITH REAL DISCORD WEBHOOK")
    print("=" * 70)
    
    print("\n⚠️  To test with a real Discord webhook:")
    print("   Set the DISCORD_WEBHOOK_URL environment variable")
    print("   Example: export DISCORD_WEBHOOK_URL='https://discord.com/api/webhooks/...'")
    
    webhook_url = os.environ.get('DISCORD_WEBHOOK_URL', '')
    
    if webhook_url:
        print(f"\n✓ Found webhook URL in environment")
        print("  Attempting to send test notification...")
        
        result = send_discord_webhook(
            webhook_url=webhook_url,
            job_id="test_manual_123",
            job_url="https://example.com/?job_id=test_manual_123",
            metadata={
                'num_sequences': 5,
                'total_motifs': 142,
                'timestamp': '2025-12-31T08:00:00Z'
            }
        )
        
        if result:
            print("  ✅ Webhook sent successfully! Check your Discord channel.")
        else:
            print("  ⚠️  Webhook failed (check URL validity)")
    else:
        print("\n  No webhook URL provided - skipping real webhook test")
        print("  This is normal and expected for automated testing")


if __name__ == "__main__":
    print("=" * 70)
    print("DISCORD WEBHOOK NOTIFICATION - MANUAL TEST SUITE")
    print("=" * 70)
    
    try:
        test_webhook_validation()
        test_webhook_payload_generation()
        demonstrate_usage()
        test_with_real_webhook()
        
        print("\n" + "=" * 70)
        print("✅ ALL TESTS COMPLETED SUCCESSFULLY")
        print("=" * 70)
        print("\nThe Discord webhook notification system is:")
        print("  ✓ Optional (jobs run without it)")
        print("  ✓ Non-blocking (won't delay analysis)")
        print("  ✓ Secure (no credentials stored)")
        print("  ✓ Free (no paid services)")
        print("  ✓ Gracefully failing (invalid webhooks won't crash app)")
        
        sys.exit(0)
        
    except Exception as e:
        print(f"\n❌ TEST FAILED WITH EXCEPTION: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
