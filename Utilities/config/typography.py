"""
Typography configuration for NBDScanner.

Refined edition: Elegant, research-grade, publication-quality typography.
Optimized for clarity, hierarchy, and long-form scientific readability.
"""

# ==================== TYPOGRAPHY & FONTS ====================
FONT_CONFIG = {
    # ==================== FONT FAMILIES ====================
    # Clean, modern, highly legible for scientific interfaces
    'primary_font': "'Inter', 'IBM Plex Sans', 'Segoe UI', system-ui, -apple-system, sans-serif",
    'monospace_font': "'JetBrains Mono', 'IBM Plex Mono', 'Fira Code', 'Consolas', monospace",
    
    # ==================== FONT SIZES ====================
    # Carefully balanced scale for analytical interfaces
    # Designed for calm hierarchy — not marketing emphasis
    
    'h1_size': '2.25rem',      # Page titles (~36px)
    'h2_size': '1.875rem',     # Major section headers (~30px)
    'h3_size': '1.5rem',       # Subsection headers (~24px)
    'h4_size': '1.25rem',      # Minor headers (~20px)
    'body_size': '1.0625rem',  # Body text (~17px) – optimal reading size
    'small_size': '0.95rem',   # Secondary UI text
    'caption_size': '0.875rem',# Figure captions, metadata
    
    # ==================== LINE HEIGHTS ====================
    # Essential for scientific readability and density control
    
    'heading_line_height': 1.2,   # Tight but breathable
    'body_line_height': 1.65,     # Ideal for long-form reading
    'compact_line_height': 1.45,  # Tables / dense UI areas
    
    # ==================== LETTER SPACING ====================
    # Subtle tightening improves authority and polish
    
    'heading_letter_spacing': '-0.01em',
    'body_letter_spacing': '0em',
    'caption_letter_spacing': '0.01em',
    
    # ==================== FONT WEIGHTS ====================
    # True typographic hierarchy (not inflated)
    
    'light_weight': 300,       # Rarely used (large display only)
    'normal_weight': 400,      # Standard body text
    'medium_weight': 500,      # UI emphasis
    'semibold_weight': 600,    # Section headers
    'bold_weight': 700,        # Strong emphasis
    'extrabold_weight': 800,   # Reserved for h1 only
    
    # ==================== MONOSPACE SETTINGS ====================
    # For code, genomic sequences, and data blocks
    
    'mono_size': '0.95rem',
    'mono_line_height': 1.6,
    'mono_letter_spacing': '0.01em',
}
