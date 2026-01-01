# ntfy.sh Push Notification Guide

## Overview

NonBDNAFinder supports **optional** push notifications via [ntfy.sh](https://ntfy.sh) - a free, account-free notification service. Get instant alerts when your analysis completes, without any signup or configuration.

## Why ntfy.sh?

- ✅ **Zero Setup**: No account or credentials required
- ✅ **Free Forever**: Public service, no paid tiers
- ✅ **Cross-Platform**: Works on web, mobile (iOS/Android), and desktop
- ✅ **Privacy-Friendly**: No personal data collected
- ✅ **Simple**: Just pick a topic name and subscribe

## How to Use

### Step 1: Choose a Topic Name

Pick any unique identifier for your notification topic:
- Examples: `nbd-job-123`, `my-analysis`, `dna-results-2026`
- No registration needed
- Anyone with the topic name can subscribe (keep it unique!)

### Step 2: Subscribe to Your Topic

**Web Browser:**
```
Visit: https://ntfy.sh/your-topic-name
```

**Mobile App (iOS/Android):**
1. Install the ntfy app from App Store or Google Play
2. Add subscription → Enter your topic name
3. Done!

**Desktop/CLI:**
```bash
# Linux/Mac
ntfy subscribe your-topic-name

# Or via curl
curl -s ntfy.sh/your-topic-name/json
```

### Step 3: Enter Topic in NBDScanner

1. Navigate to the **"Upload & Analyze"** tab
2. Upload your sequences
3. In the **"🔔 Optional: Push Notification"** section:
   - Enter your topic name (e.g., `nbd-job-123`)
4. Run your analysis

### Step 4: Get Notified!

When your analysis completes, you'll receive a push notification with:
- Job ID
- Number of sequences analyzed
- Total motifs detected
- Link to retrieve results (if deployed)

## Example Workflow

```bash
# 1. Subscribe to a topic (web browser)
Open: https://ntfy.sh/my-dna-analysis-2026

# 2. In NBDScanner:
#    - Upload: my_genome.fasta
#    - Topic: my-dna-analysis-2026
#    - Click: "Run Analysis"

# 3. Receive notification:
#    🧬 Job Completed
#    Job ID: abc123def4
#    Sequences: 5
#    Motifs: 142
#    Results: https://app-url?job_id=abc123def4
```

## Privacy & Security

- **No Account**: ntfy.sh doesn't require signup
- **Not Stored**: Topic name is only held in session memory
- **No Disk**: Topics are never saved to disk
- **Ephemeral**: Topic is forgotten when session ends
- **Optional**: Analysis works perfectly without notifications

## Important Notes

1. **Topic Security**: Choose unique topic names to prevent others from seeing your notifications
2. **Public Service**: ntfy.sh is a free public service - treat topics as public URLs
3. **Non-Blocking**: Notifications never delay or block analysis
4. **Graceful Failure**: If notification fails, your results are still saved and accessible

## Advanced: Self-Hosted ntfy

For private deployments, you can host your own ntfy server:

```bash
# Docker
docker run -d --name ntfy -p 80:80 binwiederhier/ntfy serve

# Then use your server URL
# Topic: https://your-server.com/your-topic
```

See [ntfy.sh documentation](https://docs.ntfy.sh) for more details.

## FAQ

**Q: Do I need an account?**  
A: No! ntfy.sh is completely account-free.

**Q: Is my topic stored anywhere?**  
A: No. The topic is only kept in your browser session and is never saved to disk.

**Q: What if notification fails?**  
A: Your analysis continues normally. Results are always saved and accessible via Job ID.

**Q: Can I use special characters in topic names?**  
A: Use alphanumeric characters, hyphens, and underscores for best compatibility.

**Q: Is this free?**  
A: Yes! ntfy.sh is free and open source. No paid tiers or subscriptions.

**Q: How long are notifications stored?**  
A: ntfy.sh keeps notifications for 12 hours by default. After that, they're deleted.

## Alternative: No Notifications

You can always skip notifications entirely! Just:
1. Don't enter a topic name
2. Save your Job ID when analysis completes
3. Return to the app anytime to retrieve results

---

**Learn More**: https://ntfy.sh  
**Source Code**: https://github.com/binwiederhier/ntfy
