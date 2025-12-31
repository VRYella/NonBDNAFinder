"""
Email Notification Module for NBDScanner
=========================================

Provides optional email notification functionality for job completion.

Features:
- SMTP configuration from Streamlit secrets (optional)
- Graceful failure if email not configured
- Simple text-based notification with job ID and lookup link
- Non-blocking operation (no UI delays)

Configuration (optional, via .streamlit/secrets.toml):
    [email]
    smtp_host = "smtp.gmail.com"
    smtp_port = 587
    smtp_user = "your-email@gmail.com"
    smtp_password = "your-app-password"
    from_address = "noreply@nbdscanner.org"
    support_email = "raazbiochem@gmail.com"  # Support contact
"""

import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from typing import Optional
import logging

logger = logging.getLogger(__name__)

# Default support email (can be overridden in secrets)
DEFAULT_SUPPORT_EMAIL = "raazbiochem@gmail.com"


def get_email_config(secrets) -> Optional[dict]:
    """
    Extract email configuration from Streamlit secrets.
    
    Args:
        secrets: Streamlit secrets object
        
    Returns:
        dict with email config if available, None otherwise
    """
    try:
        if hasattr(secrets, 'email'):
            config = {
                'smtp_host': secrets.email.get('smtp_host'),
                'smtp_port': secrets.email.get('smtp_port', 587),
                'smtp_user': secrets.email.get('smtp_user'),
                'smtp_password': secrets.email.get('smtp_password'),
                'from_address': secrets.email.get('from_address', secrets.email.get('smtp_user')),
                'support_email': secrets.email.get('support_email', DEFAULT_SUPPORT_EMAIL)
            }
            
            # Validate required fields
            if config['smtp_host'] and config['smtp_user'] and config['smtp_password']:
                return config
        
        return None
        
    except Exception as e:
        logger.debug(f"Email configuration not available: {e}")
        return None


def send_job_notification(
    to_email: str,
    job_id: str,
    job_url: Optional[str] = None,
    secrets=None,
    metadata: Optional[dict] = None
) -> bool:
    """
    Send email notification about completed job.
    
    This function fails gracefully - returns False on any error without raising exceptions.
    
    Args:
        to_email: Recipient email address
        job_id: Job identifier
        job_url: Optional public URL for job lookup
        secrets: Streamlit secrets object (for SMTP config)
        metadata: Optional job metadata (sequence count, motif count, etc.)
        
    Returns:
        bool: True if email sent successfully, False otherwise
    """
    # Validate input
    if not to_email or not to_email.strip():
        logger.debug("Empty email address, skipping notification")
        return False
    
    if not job_id:
        logger.warning("Cannot send notification without job ID")
        return False
    
    # Get email configuration
    if not secrets:
        logger.debug("No secrets provided, cannot send email")
        return False
    
    email_config = get_email_config(secrets)
    if not email_config:
        logger.debug("Email not configured in secrets, skipping notification")
        return False
    
    try:
        # Create message
        msg = MIMEMultipart('alternative')
        msg['Subject'] = f'NBDScanner Analysis Complete - Job {job_id}'
        msg['From'] = email_config['from_address']
        msg['To'] = to_email
        
        # Build email body
        num_sequences = metadata.get('num_sequences', 'N/A') if metadata else 'N/A'
        total_motifs = metadata.get('total_motifs', 'N/A') if metadata else 'N/A'
        timestamp = metadata.get('timestamp', 'N/A') if metadata else 'N/A'
        support_email = email_config.get('support_email', DEFAULT_SUPPORT_EMAIL)
        
        text_body = f"""
NBDScanner Analysis Complete
============================

Your Non-B DNA motif analysis has completed successfully.

Job ID: {job_id}

Analysis Summary:
- Sequences Analyzed: {num_sequences}
- Total Motifs Detected: {total_motifs}
- Completion Time: {timestamp}

To retrieve your results:
1. Return to the NBDScanner application
2. Navigate to the "Home" tab
3. Enter your Job ID: {job_id}
4. Download your results in your preferred format (CSV, Excel, JSON, BED, PDF)

"""
        
        if job_url:
            text_body += f"\nDirect Link: {job_url}\n"
        
        text_body += f"""

Important Notes:
- Results are stored securely and accessible only via your Job ID
- No personal data is retained beyond job results
- Keep your Job ID for future reference

Questions or issues? Contact: {support_email}

---
NBDScanner - Non-B DNA Motif Detection System
Developed by Dr. Venkata Rajesh Yella
"""
        
        # Attach text part
        text_part = MIMEText(text_body, 'plain')
        msg.attach(text_part)
        
        # Send email
        with smtplib.SMTP(email_config['smtp_host'], email_config['smtp_port']) as server:
            server.starttls()
            server.login(email_config['smtp_user'], email_config['smtp_password'])
            server.send_message(msg)
        
        logger.info(f"Email notification sent successfully to {to_email} for job {job_id}")
        return True
        
    except smtplib.SMTPAuthenticationError:
        logger.warning(f"SMTP authentication failed for {to_email} - check credentials")
        return False
    except smtplib.SMTPException as e:
        logger.warning(f"SMTP error sending notification to {to_email}: {e}")
        return False
    except Exception as e:
        logger.warning(f"Failed to send email notification to {to_email}: {e}")
        return False


def validate_email_format(email: str) -> bool:
    """
    Basic email format validation.
    
    Args:
        email: Email address to validate
        
    Returns:
        bool: True if format appears valid
    """
    if not email or not email.strip():
        return False
    
    email = email.strip()
    
    # Basic checks
    if '@' not in email:
        return False
    
    if email.count('@') != 1:
        return False
    
    local, domain = email.split('@')
    
    if not local or not domain:
        return False
    
    if '.' not in domain:
        return False
    
    return True
