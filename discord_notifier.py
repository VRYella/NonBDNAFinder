"""
Discord Webhook Notification Module for NBDScanner
===================================================

Provides optional Discord webhook notification functionality for job completion.

Features:
- Simple HTTP POST to Discord webhook URL
- Graceful failure (non-blocking, no exceptions raised)
- No credentials required (webhook URL is the only secret)
- Timeout protection (≤ 5 seconds)
- Rich Discord embed formatting

Usage:
    from discord_notifier import send_discord_webhook
    
    success = send_discord_webhook(
        webhook_url="https://discord.com/api/webhooks/...",
        job_id="abc123def4",
        job_url="https://app-url/?job_id=abc123def4",
        metadata={"num_sequences": 5, "total_motifs": 142}
    )
"""

import requests
import logging
from typing import Optional, Dict
from datetime import datetime

logger = logging.getLogger(__name__)

# Discord webhook timeout (in seconds)
WEBHOOK_TIMEOUT = 5


def send_discord_webhook(
    webhook_url: str,
    job_id: str,
    job_url: Optional[str] = None,
    metadata: Optional[Dict] = None
) -> bool:
    """
    Send Discord webhook notification about completed job.
    
    This function is best-effort only and never raises exceptions.
    All errors are caught, logged, and result in returning False.
    
    Args:
        webhook_url: Discord webhook URL (must start with https://discord.com/api/webhooks/)
        job_id: Job identifier
        job_url: Optional public URL for job lookup
        metadata: Optional job metadata (num_sequences, total_motifs, etc.)
        
    Returns:
        bool: True if webhook sent successfully, False otherwise
    """
    # Validate inputs
    if not webhook_url or not webhook_url.strip():
        logger.debug("Empty webhook URL, skipping notification")
        return False
    
    webhook_url = webhook_url.strip()
    
    # Basic URL validation (simple presence check only, as per requirements)
    if not webhook_url.startswith('https://'):
        logger.debug("Webhook URL must use HTTPS")
        return False
    
    if not job_id:
        logger.warning("Cannot send notification without job ID")
        return False
    
    try:
        # Extract metadata for display
        num_sequences = metadata.get('num_sequences', 'N/A') if metadata else 'N/A'
        total_motifs = metadata.get('total_motifs', 'N/A') if metadata else 'N/A'
        timestamp = metadata.get('timestamp', datetime.now().isoformat()) if metadata else datetime.now().isoformat()
        
        # Format timestamp for better readability
        try:
            dt = datetime.fromisoformat(timestamp.replace('Z', '+00:00'))
            formatted_time = dt.strftime('%Y-%m-%d %H:%M:%S UTC')
        except:
            formatted_time = timestamp
        
        # Build Discord embed payload
        # Discord embeds provide rich formatting with colors, fields, and structure
        embed = {
            "title": "🧬 NBDScanner Analysis Complete",
            "description": f"Your Non-B DNA motif analysis has completed successfully.",
            "color": 3447003,  # Blue color (0x3498DB)
            "fields": [
                {
                    "name": "Job ID",
                    "value": f"`{job_id}`",
                    "inline": False
                },
                {
                    "name": "📊 Analysis Summary",
                    "value": f"**Sequences Analyzed:** {num_sequences}\n**Total Motifs Detected:** {total_motifs}",
                    "inline": False
                }
            ],
            "footer": {
                "text": "NonBDNAFinder - Non-B DNA Motif Detection System"
            },
            "timestamp": datetime.now().isoformat()
        }
        
        # Add job URL field if provided
        if job_url:
            embed["fields"].insert(1, {
                "name": "🔗 Direct Access",
                "value": f"[Click here to view results]({job_url})",
                "inline": False
            })
        else:
            embed["fields"].insert(1, {
                "name": "📥 Retrieve Results",
                "value": "Return to the NBDScanner app and enter your Job ID to download results",
                "inline": False
            })
        
        # Discord webhook payload
        payload = {
            "embeds": [embed]
        }
        
        # Send webhook with timeout protection
        response = requests.post(
            webhook_url,
            json=payload,
            timeout=WEBHOOK_TIMEOUT,
            headers={'Content-Type': 'application/json'}
        )
        
        # Check response status
        if response.status_code in [200, 204]:
            logger.info(f"Discord webhook sent successfully for job {job_id}")
            return True
        else:
            logger.warning(f"Discord webhook failed with status {response.status_code} for job {job_id}")
            return False
            
    except requests.exceptions.Timeout:
        logger.warning(f"Discord webhook timed out after {WEBHOOK_TIMEOUT}s for job {job_id}")
        return False
    except requests.exceptions.RequestException as e:
        logger.warning(f"Discord webhook network error for job {job_id}: {e}")
        return False
    except Exception as e:
        logger.warning(f"Discord webhook unexpected error for job {job_id}: {e}")
        return False


def validate_webhook_url(url: str) -> bool:
    """
    Basic Discord webhook URL validation.
    
    Performs simple presence check only (no aggressive validation as per requirements).
    
    Args:
        url: Webhook URL to validate
        
    Returns:
        bool: True if URL format appears valid
    """
    if not url or not url.strip():
        return False
    
    url = url.strip()
    
    # Must be HTTPS
    if not url.startswith('https://'):
        return False
    
    # Should contain discord.com domain (optional check, not enforced)
    # This is just a helpful hint, not a hard requirement
    if 'discord.com/api/webhooks/' in url or 'discordapp.com/api/webhooks/' in url:
        return True
    
    # Accept any HTTPS URL (allows for Discord-compatible webhook services)
    return True
