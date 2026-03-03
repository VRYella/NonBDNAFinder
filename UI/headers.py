"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Headers Module - Uniform section heading component                           │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
│ Provides consistent heading style across all pages                           │
└──────────────────────────────────────────────────────────────────────────────┘
"""
# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
import streamlit as st
import html
from Utilities.config.colors import HOME_COLORS, GLOBAL_COLORS

# ═══════════════════════════════════════════════════════════════════════════════
# TUNABLE PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════════
PAGE_HEADING_COLORS = {
    'Home': {'primary': '#2C3E50', 'pattern': 'linear-gradient(135deg, #2C3E50 0%, #34495E 50%, #5D8AA8 100%)'},
    'Upload & Analyze': {'primary': '#4F46E5', 'pattern': 'linear-gradient(135deg, #4F46E5 0%, #6366F1 50%, #818CF8 100%)'},
    'Results': {'primary': '#0EA5A4', 'pattern': 'linear-gradient(135deg, #0891B2 0%, #0EA5A4 50%, #06B6D4 100%)'},
    'Downloads': {'primary': '#16A34A', 'pattern': 'linear-gradient(135deg, #15803D 0%, #16A34A 50%, #22C55E 100%)'},
    'Documentation': {'primary': '#B91C1C', 'pattern': 'linear-gradient(135deg, #991B1B 0%, #B91C1C 50%, #DC2626 100%)'},
}
HEADING_PADDING = "0.15rem 1.2rem"; HEADING_FONT_SIZE = "1.5rem"; HEADING_BORDER_RADIUS = "8px"
# ═══════════════════════════════════════════════════════════════════════════════


def render_section_heading(title: str, page: str = None):
    """Render a uniform section heading with page-specific background colors.
    
    This is the canonical heading component for NonBDNAFinder. All section
    headings should use this function to maintain visual consistency.
    
    Features:
    - Thin box with page-specific gradient background color patterns
    - Increased font size for better readability
    - White glowing text with text-shadow effect
    - Compact, polished appearance
    
    Args:
        title: The heading text to display
        page: Optional page name for page-specific colors (Home, Upload & Analyze, Results, Downloads, Documentation)
    """
    # Escape HTML to prevent XSS
    safe_title = html.escape(title)
    
    # Determine colors based on page or use default
    if page and page in PAGE_HEADING_COLORS:
        colors = PAGE_HEADING_COLORS[page]
        background = colors['pattern']
        border_color = colors['primary']
    else:
        # Default to blue gradient
        primary_color = HOME_COLORS['primary']  # #0091FF
        secondary_color = HOME_COLORS['secondary']  # #00B4FF
        background = f"linear-gradient(135deg, {primary_color} 0%, {secondary_color} 100%)"
        border_color = primary_color
    
    white = GLOBAL_COLORS['white']
    
    st.markdown(f"""
    <div class="section-heading-box" style="
        background: {background};
        padding: {HEADING_PADDING};
        border-radius: {HEADING_BORDER_RADIUS};
        margin: 0.4rem 0 0.6rem 0;
        border: 2px solid {border_color};
        text-align: center;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.25);
        width: 100%;
    ">
        <div title="{safe_title}" style="
            margin: 0;
            font-size: {HEADING_FONT_SIZE};
            font-weight: 800;
            color: {white};
            text-shadow: 0 0 10px rgba(255,255,255,0.8), 0 0 20px rgba(255,255,255,0.6);
            letter-spacing: -0.01em;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        ">
            {safe_title}
        </div>
    </div>
    """, unsafe_allow_html=True)
