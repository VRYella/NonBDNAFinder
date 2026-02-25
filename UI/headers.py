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
    'Home': {'primary': '#1e3a8a', 'pattern': 'linear-gradient(135deg, #1e3a8a 0%, #1e40af 100%)'},
    'Upload & Analyze': {'primary': '#065f46', 'pattern': 'linear-gradient(135deg, #065f46 0%, #047857 100%)'},
    'Results': {'primary': '#3730a3', 'pattern': 'linear-gradient(135deg, #3730a3 0%, #4338ca 100%)'},
    'Downloads': {'primary': '#0c4a6e', 'pattern': 'linear-gradient(135deg, #0c4a6e 0%, #0369a1 100%)'},
    'Documentation': {'primary': '#78350f', 'pattern': 'linear-gradient(135deg, #78350f 0%, #92400e 100%)'},
}
HEADING_PADDING = "0.25rem 1.2rem"; HEADING_FONT_SIZE = "1.15rem"; HEADING_BORDER_RADIUS = "6px"
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
    <div style="
        background: {background};
        padding: 0.25rem 1.2rem;
        border-radius: 6px;
        margin: 0.4rem 0 0.6rem 0;
        border: 1px solid {border_color};
        text-align: center;
        box-shadow: 0 0 12px rgba(0, 145, 255, 0.3);
        width: 100%;
    ">
        <h2 style="
            margin: 0;
            font-size: 1.15rem;
            font-weight: 700;
            color: {white};
            text-shadow: 0 0 10px rgba(255,255,255,0.8), 0 0 20px rgba(255,255,255,0.6);
            letter-spacing: -0.01em;
        ">
            {safe_title}
        </h2>
    </div>
    """, unsafe_allow_html=True)
