"""
Layout configuration for NBDScanner.

This module contains layout and spacing settings:
- Layout mode (wide/centered)
- Border radius values
- Padding and margins
- Column gaps
- Card styling
- Z-index layers
"""

# ==================== LAYOUT & SPACING ====================
# Control spacing, padding, and visual structure
LAYOUT_CONFIG = {
    # Streamlit layout mode
    # 'wide': Uses full browser width (recommended for data apps)
    # 'centered': Centers content with maximum width
    'layout_mode': 'wide',
    
    # Border radius values (in pixels) for rounded corners
    'border_radius': {
        'small': '8px',    # Small elements (buttons, tags)
        'medium': '12px',  # Medium elements (cards, inputs)
        'large': '16px',   # Large elements (panels, sections)
        'pill': '50px',    # Fully rounded (pills, progress bars)
    },
    
    # Padding and margins (in rem units)
    'padding': {
        'small': '0.5rem',   # Tight spacing
        'medium': '1rem',    # Standard spacing
        'large': '1.5rem',   # Generous spacing
        'xlarge': '2rem',    # Extra spacing for sections
    },
    
    # Margins (in rem units)
    'margin': {
        'none': '0',
        'small': '0.5rem',
        'medium': '1rem',
        'large': '1.5rem',
        'xlarge': '2rem',
    },
    
    # Column gaps for multi-column layouts
    # Options: 'small', 'medium', 'large'
    'column_gap': 'large',
    
    # Card styling for panels and containers
    'card_shadow': '0 4px 20px rgba(0,0,0,0.08)',  # Soft shadow for depth
    'card_border': '1px solid #e5e7eb',             # Subtle border
    'card_shadow_hover': '0 8px 30px rgba(0,0,0,0.12)',  # Elevated shadow on hover
    
    # Section spacing
    'section_spacing': '2rem',      # Space between major sections
    'subsection_spacing': '1.5rem', # Space between subsections
    'element_spacing': '1rem',      # Space between UI elements
    
    # Container widths
    'container_max_width': '1400px',  # Maximum width for centered content
    'sidebar_width': '300px',          # Sidebar width
    
    # Z-index layers (for stacking context)
    'z_index': {
        'base': 1,
        'dropdown': 100,
        'sticky': 200,
        'modal': 1000,
        'tooltip': 2000,
    },
}
