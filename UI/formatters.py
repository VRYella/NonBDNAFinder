"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Formatters Module - Formatting utilities for NBDScanner                      │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
│ Functions for formatting time values, sequence limits                        │
└──────────────────────────────────────────────────────────────────────────────┘
"""
# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
# No external imports required

# ═══════════════════════════════════════════════════════════════════════════════
# TUNABLE PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════════
TIME_FORMAT_PRECISION = 1  # Decimal places for seconds display
# ═══════════════════════════════════════════════════════════════════════════════


def format_time_scientific(seconds: float) -> str:
    """
    Format elapsed time in simple MM:SS format.
    
    This format provides:
    - Human-readable minutes and seconds
    - No hours or microseconds (simplified display)
    - Consistent display across all workflows
    
    Args:
        seconds: Elapsed time in seconds (float)
        
    Returns:
        Formatted time string (e.g., "02:15" or "125:32")
        
    Examples:
        >>> format_time_scientific(0.234)
        "00:00"
        >>> format_time_scientific(135.678)
        "02:15"
        >>> format_time_scientific(5432.123)
        "90:32"
    """
    minutes = int(seconds // 60)
    secs = int(seconds % 60)
    
    return f"{minutes:02d}:{secs:02d}"


def format_time_compact(seconds: float) -> str:
    """
    Format elapsed time in MM:SS format for compact displays.
    
    Simple minutes:seconds format for all durations.
    
    Args:
        seconds: Elapsed time in seconds (float)
        
    Returns:
        Formatted time string (e.g., "02:15" or "125:32")
    """
    minutes = int(seconds // 60)
    secs = int(seconds % 60)
    return f"{minutes:02d}:{secs:02d}"


def format_time(seconds):
    """Format time in seconds to a human-readable string (compact format).
    
    This format is more compact than format_time_human() and is suitable for
    inline display in UI messages and progress updates.
    
    Args:
        seconds: Time in seconds (float or int)
        
    Returns:
        Formatted string (e.g., "45.3s", "12m 30s", "2h 15m")
    
    Examples:
        >>> format_time(45.3)
        '45.3s'
        >>> format_time(750)
        '12m 30s'
        >>> format_time(7800)
        '2h 10m'
    
    See Also:
        format_time_human() - More detailed format for performance reports
    """
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        mins = int(seconds // 60)
        secs = int(seconds % 60)
        return f"{mins}m {secs}s"
    else:
        hours = int(seconds // 3600)
        mins = int((seconds % 3600) // 60)
        return f"{hours}h {mins}m"


def format_time_human(seconds):
    """Format time in seconds to a detailed human-readable string.
    
    Provides more detailed breakdown than format_time() for performance reports
    and summary displays where space is not a constraint.
    
    Args:
        seconds: Time in seconds (float or int)
        
    Returns:
        Formatted string with full detail (e.g., "2 hours, 15 minutes, 30 seconds")
    
    Examples:
        >>> format_time_human(45.3)
        '45.3 seconds'
        >>> format_time_human(750)
        '12 minutes, 30 seconds'
        >>> format_time_human(7800)
        '2 hours, 10 minutes'
    
    See Also:
        format_time() - Compact format for inline UI messages
    """
    if seconds < 1:
        return f"{seconds*1000:.0f} milliseconds"
    elif seconds < 60:
        return f"{seconds:.1f} seconds"
    elif seconds < 3600:
        mins = int(seconds // 60)
        secs = seconds % 60
        if secs >= 1:
            return f"{mins} minute{'s' if mins != 1 else ''}, {secs:.1f} seconds"
        else:
            return f"{mins} minute{'s' if mins != 1 else ''}"
    else:
        hours = int(seconds // 3600)
        mins = int((seconds % 3600) // 60)
        secs = seconds % 60
        result = f"{hours} hour{'s' if hours != 1 else ''}"
        if mins > 0:
            result += f", {mins} minute{'s' if mins != 1 else ''}"
        if secs >= 1:
            result += f", {secs:.0f} seconds"
        return result


def format_sequence_limit():
    """Format the sequence limit for display - now shows 'unlimited' since limit is removed"""
    return "unlimited (chunked processing enabled)"
