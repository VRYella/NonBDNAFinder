"""
Animation configuration for NBDScanner.

This module contains animation and transition settings:
- Transition durations
- Easing functions
- Animation names
- Motion preferences
"""

# ==================== ANIMATION & TRANSITION SETTINGS ====================
# Control animation timing and effects for smooth UI interactions
ANIMATION_CONFIG = {
    # Transition durations (in seconds)
    'transition_fast': '0.15s',
    'transition_normal': '0.2s',
    'transition_slow': '0.3s',
    
    # Easing functions
    'easing_smooth': 'ease',
    'easing_in': 'ease-in',
    'easing_out': 'ease-out',
    'easing_in_out': 'ease-in-out',
    
    # Animation names (CSS keyframes defined in styles.css)
    'fade_in': 'fade-in',
    'pulse': 'pulse-dot',
    'shimmer': 'progress-shimmer',
    
    # Enable/disable animations globally
    'enable_animations': True,
    'reduce_motion': False,  # Respect user's prefers-reduced-motion setting
}
