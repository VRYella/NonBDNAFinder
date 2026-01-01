"""
ntfy.sh Push Notification Module for NBDScanner
================================================

Provides account-free, zero-configuration push notifications via ntfy.sh.

Features:
- No signup or credentials required
- Simple topic-based messaging
- HTTP POST to public ntfy.sh service
- Graceful failure (non-blocking, no exceptions raised)
- Timeout protection (≤ 5 seconds)
- Free public usage

Usage:
    from ntfy_notifier import send_ntfy_notification
    
    success = send_ntfy_notification(
        topic="nbd-job-123",
        job_id="abc123def4",
        app_url="https://app-url",
        metadata={"num_sequences": 5, "total_motifs": 142}
    )

Topic Format:
    - Simple string identifier (alphanumeric + hyphens)
    - Example: "nbd-job-abc123", "my-analysis-results"
    - No validation or registration needed
    - Anyone with the topic name can subscribe
"""

import requests
import logging
from typing import Optional, Dict

logger = logging.getLogger(__name__)

# ntfy.sh public server
NTFY_SERVER = "https://ntfy.sh"

# Request timeout (in seconds)
NTFY_TIMEOUT = 5


def send_ntfy_notification(
    topic: str,
    job_id: str,
    app_url: Optional[str] = None,
    metadata: Optional[Dict] = None
) -> bool:
    """
    Send ntfy.sh push notification about completed job.
    
    This function is best-effort only and never raises exceptions.
    All errors are caught, logged, and result in returning False.
    
    Args:
        topic: ntfy topic name (user-provided string identifier)
        job_id: Job identifier
        app_url: Optional base URL of the application
        metadata: Optional job metadata (num_sequences, total_motifs, etc.)
        
    Returns:
        bool: True if notification sent successfully, False otherwise
    """
    # Validate inputs
    if not topic or not topic.strip():
        logger.debug("Empty topic, skipping notification")
        return False
    
    topic = topic.strip()
    
    if not job_id:
        logger.warning("Cannot send notification without job ID")
        return False
    
    try:
        # Extract metadata for display
        num_sequences = metadata.get('num_sequences', 'N/A') if metadata else 'N/A'
        total_motifs = metadata.get('total_motifs', 'N/A') if metadata else 'N/A'
        
        # Build notification message
        # Simple plain-text format for universal compatibility
        title = "🧬 Job Completed"
        
        message = f"Job ID: {job_id}\n"
        message += f"Sequences: {num_sequences}\n"
        message += f"Motifs: {total_motifs}\n"
        
        # Add results URL if app_url provided
        if app_url:
            # Construct direct link to results
            results_url = f"{app_url}?job_id={job_id}"
            message += f"\nResults: {results_url}"
        else:
            message += f"\nRetrieve results with Job ID in app"
        
        # ntfy.sh uses simple HTTP POST with headers
        # Title and message are sent as headers for better mobile display
        url = f"{NTFY_SERVER}/{topic}"
        
        headers = {
            "Title": title,
            "Priority": "default",
            "Tags": "dna,science"
        }
        
        # Send notification with timeout protection
        response = requests.post(
            url,
            data=message.encode('utf-8'),
            headers=headers,
            timeout=NTFY_TIMEOUT
        )
        
        # Check response status
        if response.status_code == 200:
            logger.info(f"ntfy notification sent successfully for job {job_id} to topic {topic}")
            return True
        else:
            logger.warning(f"ntfy notification failed with status {response.status_code} for job {job_id}")
            return False
            
    except requests.exceptions.Timeout:
        logger.warning(f"ntfy notification timed out after {NTFY_TIMEOUT}s for job {job_id}")
        return False
    except requests.exceptions.RequestException as e:
        logger.warning(f"ntfy notification network error for job {job_id}: {e}")
        return False
    except Exception as e:
        logger.warning(f"ntfy notification unexpected error for job {job_id}: {e}")
        return False


def validate_topic(topic: str) -> bool:
    """
    Basic ntfy topic validation.
    
    Performs minimal validation only (no aggressive checks as per requirements).
    
    Args:
        topic: Topic name to validate
        
    Returns:
        bool: True if topic appears valid
    """
    if not topic or not topic.strip():
        return False
    
    # ntfy.sh accepts most reasonable topic names
    # We don't enforce strict validation per requirements
    return True
