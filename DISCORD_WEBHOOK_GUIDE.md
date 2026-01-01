# Discord Webhook Notifications for NBDScanner

## Overview

NBDScanner now supports **optional Discord webhook notifications** for job completion. This provides instant, free notifications without requiring email configuration or SMTP credentials.

## Features

✅ **Free** - No paid services required  
✅ **Instant** - Real-time notifications via Discord  
✅ **Optional** - Jobs run successfully with or without webhook  
✅ **Non-blocking** - Won't delay analysis  
✅ **Secure** - No credentials stored, webhook URL not persisted  
✅ **Graceful** - Invalid webhooks won't crash the app  

## How to Use

### Step 1: Get a Discord Webhook URL

1. Open your Discord server
2. Go to **Server Settings** → **Integrations** → **Webhooks**
3. Click **New Webhook** (or select an existing one)
4. Copy the **Webhook URL**

### Step 2: Use in NBDScanner

1. Navigate to the **"Upload & Analyze"** tab
2. Upload your FASTA sequences
3. **Optional:** Paste your Discord webhook URL in the notification field
4. Click **"Run Analysis"**
5. When complete, a Discord notification is sent automatically

## Notification Content

The Discord notification includes:

- 🧬 **Job ID** - For retrieving results
- 📊 **Analysis Summary** - Sequences analyzed, motifs detected
- 🔗 **Direct Link** - Quick access to results (if configured)
- ⏰ **Timestamp** - When analysis completed

## Example Notification

```
🧬 NBDScanner Analysis Complete

Job ID: `abc123def4`

📊 Analysis Summary
Sequences Analyzed: 5
Total Motifs Detected: 142

📥 Retrieve Results
Return to the NBDScanner app and enter your Job ID to download results
```

## Privacy & Security

- ✅ Webhook URL is **NOT stored on disk**
- ✅ Only held in session memory during analysis
- ✅ No personal data collected or transmitted
- ✅ Webhook is completely **optional**
- ✅ Analysis runs successfully with or without webhook

## Technical Details

### Implementation

The Discord webhook notification system:

- Uses HTTPS POST with JSON payload
- Timeout protection (≤ 5 seconds)
- Rich Discord embed formatting
- Graceful error handling (all exceptions caught)

### Testing

Run the test suite to verify functionality:

```bash
# Test basic webhook functionality
python test_job_management.py

# Test complete workflow
python test_app_workflow.py

# Test end-to-end integration
python test_e2e_job_management.py

# Manual demonstration and testing
python test_discord_webhook_manual.py
```

### Optional: Test with Real Webhook

To test with a real Discord webhook:

```bash
export DISCORD_WEBHOOK_URL='https://discord.com/api/webhooks/...'
python test_discord_webhook_manual.py
```

## Comparison: Discord vs Email

| Feature | Discord Webhook | Email (Previous) |
|---------|----------------|------------------|
| **Cost** | Free | Requires SMTP service |
| **Setup** | Copy/paste URL | Configure SMTP credentials |
| **Speed** | Instant | Variable (SMTP delays) |
| **Reliability** | High | Depends on SMTP config |
| **Security** | No credentials | Requires stored secrets |
| **Privacy** | No personal data | Email address required |
| **Maintenance** | Zero | SMTP config updates |

## Migration from Email

The email notification system has been replaced with Discord webhooks:

- **Old:** Email input field (`job_email`)
- **New:** Discord webhook URL input (`discord_webhook_url`)

Existing functionality preserved:
- ✅ Job ID generation and display
- ✅ Result persistence to disk
- ✅ Job lookup and retrieval
- ✅ Optional notifications (now via Discord)

## Troubleshooting

### Webhook not sending?

1. **Check URL format** - Must start with `https://`
2. **Verify webhook is active** - Test in Discord settings
3. **Check permissions** - Ensure webhook has post permissions
4. **Network connectivity** - Webhook requires internet access

### Invalid webhook error?

- This is **normal and expected** if webhook URL is incorrect
- Analysis will **still complete successfully**
- Results are **always accessible via Job ID**

## Support

For issues or questions:
- 📧 Email: yvrajesh_bt@kluniversity.in
- 🐙 GitHub: [VRYella/NonBDNAFinder](https://github.com/VRYella/NonBDNAFinder)

---

**Note:** Discord webhook notifications are optional. All analysis features work perfectly without providing a webhook URL.
