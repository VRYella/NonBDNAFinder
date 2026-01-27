"""
Typography configuration for NBDScanner.

This module contains font settings for the application including:
- Font families (primary and monospace)
- Font sizes for different heading levels and body text
- Font weights from light to extrabold
"""

# ==================== TYPOGRAPHY & FONTS ====================
# Control all font settings for the application - MODERN & READABLE
# Optimized for modern high-resolution displays with excellent readability
FONT_CONFIG = {
    # Primary font families (in order of preference)
    # The browser will use the first available font in the list
    'primary_font': "'Inter', 'IBM Plex Sans', 'Segoe UI', system-ui, -apple-system, sans-serif",
    'monospace_font': "'JetBrains Mono', 'Fira Code', 'Consolas', monospace",
    
    # Font sizes (in rem units, where 1rem ≈ 16px in most browsers)
    # Enhanced sizes for modern, bold, research-quality appearance
    'h1_size': '2.75rem',     # Main page headers - bold, impactful
    'h2_size': '2.0rem',      # Section headers - clear hierarchy
    'h3_size': '1.5rem',      # Subsection headers - organized structure
    'h4_size': '1.25rem',     # Small headers - subtle distinction
    'body_size': '1.0rem',    # Body text, paragraphs - optimal readability
    'small_size': '0.9rem',   # Small text, notes - clear but compact
    'caption_size': '0.8rem', # Captions, footnotes - supporting information
    
    # Font weights (100-900, where 400 is normal and 700 is bold)
    'light_weight': 300,
    'normal_weight': 400,
    'medium_weight': 500,
    'semibold_weight': 600,
    'bold_weight': 700,
    'extrabold_weight': 800,
}
