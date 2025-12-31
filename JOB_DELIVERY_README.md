# Job-ID Based Result Delivery

This document describes the new job-ID based result delivery system with optional email notification.

## Overview

The NBDScanner application now provides **persistent result storage** with unique Job IDs, allowing users to retrieve analysis results at any time without requiring authentication or accounts.

## Key Features

### 🆔 Unique Job IDs
- **Format**: 10-character hexadecimal strings (e.g., `a1b2c3d4e5`)
- **Generation**: Automatically created when analysis starts
- **Collision resistance**: 2^40 (≈1 trillion) possible values
- **Non-guessable**: Cryptographically random using UUID4

### 💾 Result Persistence
- **Storage location**: `results/<job_id>/`
- **Files created**:
  - `results.json`: Complete motif detection results
  - `sequences.json`: Original DNA sequences
  - `metadata.json`: Job metadata (timestamp, statistics)
- **Isolation**: Each job stored in separate directory
- **Retention**: Persists across app restarts

### 🔍 Job Lookup
- **Location**: Home tab → "Retrieve Previous Results" section
- **Input**: Enter 10-character Job ID
- **Output**: Complete job summary and downloadable results
- **Validation**: Graceful error handling for invalid IDs

### 📧 Optional Email Notification
- **Configuration**: Via `.streamlit/secrets.toml` (optional)
- **Behavior**: Gracefully fails if not configured
- **Content**: Job ID, summary statistics, retrieval instructions
- **Privacy**: Email address not stored on disk

## User Experience Flow

### 1. Submit Analysis
```
User uploads sequence → Enters optional email → Clicks "Run Analysis"
```

### 2. Job ID Display
```
✅ Job ID generated immediately: a1b2c3d4e5
📋 User saves Job ID for later retrieval
```

### 3. Analysis Execution
```
🧬 Analysis runs → Results generated → Progress updates
```

### 4. Result Persistence
```
💾 Results saved to disk under Job ID
📧 Optional email sent (if configured and provided)
```

### 5. Download Results
```
📥 Download tab shows Job ID
📊 Download buttons for CSV, Excel, JSON, BED, PDF
```

### 6. Later Retrieval
```
🔍 User returns to app → Enters Job ID on Home tab
✅ Previous results loaded → Can view and download
```

## Configuration

### Email Notification (Optional)

Create `.streamlit/secrets.toml`:

```toml
[email]
smtp_host = "smtp.gmail.com"
smtp_port = 587
smtp_user = "your-email@gmail.com"
smtp_password = "your-app-password"
from_address = "noreply@nbdscanner.org"
support_email = "raazbiochem@gmail.com"  # Optional, defaults to raazbiochem@gmail.com
```

**Important**: The app works fully without this configuration. Email is completely optional.

### Gmail App Password Setup

If using Gmail for notifications:

1. Enable 2-Factor Authentication on your Google account
2. Go to https://myaccount.google.com/apppasswords
3. Generate an app password for "Mail"
4. Use this password in `smtp_password` field

## API Reference

### Job Manager Module (`job_manager.py`)

#### `generate_job_id() -> str`
Generate a unique 10-character job ID.

#### `save_job_results(job_id, results, sequences, names, metadata) -> bool`
Save analysis results to disk under job ID.

**Parameters**:
- `job_id`: Job identifier
- `results`: List of motif detection results
- `sequences`: List of DNA sequences
- `names`: List of sequence names
- `metadata`: Optional additional metadata

**Returns**: True if successful

#### `load_job_results(job_id) -> Optional[Tuple]`
Load results for a given job ID.

**Returns**: Tuple of (results, sequences, names, metadata) or None

#### `job_exists(job_id) -> bool`
Check if a job exists on disk.

#### `get_job_summary(job_id) -> Optional[Dict]`
Get job metadata without loading full results.

### Email Notifier Module (`email_notifier.py`)

#### `send_job_notification(to_email, job_id, job_url, secrets, metadata) -> bool`
Send email notification about completed job.

**Parameters**:
- `to_email`: Recipient email address
- `job_id`: Job identifier
- `job_url`: Optional public URL for job lookup
- `secrets`: Streamlit secrets object
- `metadata`: Optional job metadata

**Returns**: True if email sent, False otherwise (graceful failure)

#### `validate_email_format(email) -> bool`
Basic email format validation.

## Security & Privacy

### Job ID Security
- **Non-sequential**: UUIDs prevent enumeration attacks
- **Large space**: 2^40 possible IDs (1,099,511,627,776)
- **Random**: Cryptographically secure random generation
- **No session binding**: Jobs accessible only via ID

### Data Privacy
- **No PII storage**: Email addresses not saved to disk
- **Job isolation**: Results stored in separate directories
- **No authentication**: Access control via Job ID only
- **Minimal metadata**: Only analysis statistics stored

### Recommendations
1. Treat Job IDs like passwords - don't share publicly
2. Consider implementing job expiry for production deployments
3. Monitor `results/` directory size in production
4. Add rate limiting if deploying publicly

## Testing

### Automated Tests

Run all tests:
```bash
python test_job_management.py       # Unit tests (6 tests)
python test_e2e_job_management.py   # Integration test
python test_app_workflow.py         # Workflow simulation
```

All tests should pass with "✅" indicators.

### Manual Testing

See `TESTING_CHECKLIST.md` for comprehensive manual testing procedures.

Quick smoke test:
1. Load example sequence
2. Run analysis (note Job ID)
3. Check Download tab (Job ID visible)
4. Go to Home tab
5. Enter Job ID → Load results
6. Verify results match

## Troubleshooting

### Job Not Found
- **Cause**: Invalid Job ID or results not saved
- **Solution**: Check Job ID spelling (case-sensitive, 10 chars)
- **Check**: Verify `results/<job_id>/` directory exists

### Email Not Sent
- **Cause**: Email configuration missing or invalid
- **Solution**: Check `.streamlit/secrets.toml` exists and is valid
- **Note**: App still works - results saved and accessible via Job ID

### Results Directory Growing
- **Cause**: No automatic cleanup (by design)
- **Solution**: Implement periodic cleanup script (future enhancement)
- **Temporary**: Manually delete old `results/<job_id>/` directories

## Future Enhancements

Potential improvements (out of scope for this PR):

1. **Job Expiry**: Auto-delete jobs after N days
2. **Job History UI**: List all user's previous jobs
3. **ZIP Downloads**: Bundle all formats in single download
4. **HTML Email**: Rich-formatted email templates
5. **Progress Tracking**: Real-time progress via WebSocket
6. **Rate Limiting**: Prevent abuse in public deployments
7. **Cloud Storage**: S3/Azure Blob integration for scalability
8. **Job Sharing**: Generate shareable links with access codes

## Architecture Decisions

### Why Disk Storage?
- **Simplicity**: No database required
- **Portability**: Works on any file system
- **Scalability**: Easy to migrate to cloud storage
- **Backup**: Standard file backup tools work

### Why Optional Email?
- **Public Tool**: Many users won't provide email
- **Privacy**: Don't force email collection
- **Reliability**: App works even if email fails
- **Configuration**: Avoid SMTP complexity

### Why 10-Character IDs?
- **Balance**: Short enough to type, long enough for security
- **URL-safe**: Can be used in query parameters
- **Collision-resistant**: 2^40 space sufficient for millions of jobs
- **User-friendly**: Easier to communicate than longer IDs

## Contributing

When extending this feature:

1. **Maintain backwards compatibility**: Old job files must remain loadable
2. **Test graceful failures**: Email, storage, etc. should fail gracefully
3. **Preserve privacy**: Don't store unnecessary personal data
4. **Document changes**: Update this README for API changes
5. **Add tests**: All new functionality needs unit tests

## License

Same as main NBDScanner project (MIT License).

## Contact

For questions or issues:
- **Developer**: Dr. Venkata Rajesh Yella
- **Email**: raazbiochem@gmail.com
- **GitHub**: https://github.com/VRYella/NonBDNAFinder

---

**Version**: 1.0
**Last Updated**: 2025-12-31
**Status**: Production Ready ✅
