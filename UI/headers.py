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
    'Home': {'primary': '#DC2626', 'border': '#fca5a5', 'pattern': 'linear-gradient(135deg, #7f1d1d 0%, #B91C1C 40%, #DC2626 75%, #EF4444 100%)'},
    'Upload & Analyze': {'primary': '#4F46E5', 'border': '#C7D2FE', 'pattern': 'linear-gradient(135deg, #1e1b4b 0%, #3730A3 40%, #4F46E5 75%, #818CF8 100%)'},
    'Results': {'primary': '#0EA5A4', 'border': '#67e8f9', 'pattern': 'linear-gradient(135deg, #164e63 0%, #0891B2 40%, #0EA5A4 75%, #22D3EE 100%)'},
    'Downloads': {'primary': '#16A34A', 'border': '#86efac', 'pattern': 'linear-gradient(135deg, #14532d 0%, #16A34A 40%, #22C55E 75%, #4ADE80 100%)'},
    'Documentation': {'primary': '#EA580C', 'border': '#FDBA74', 'pattern': 'linear-gradient(135deg, #7c2d12 0%, #EA580C 40%, #FB923C 75%, #FBBF24 100%)'},
}
HEADING_PADDING = "0.25rem 1.4rem"; HEADING_FONT_SIZE = "1.65rem"; HEADING_BORDER_RADIUS = "12px"
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
        border_color = colors['border']
    else:
        # Fallback to HOME_COLORS gradient when page-specific colors are not defined
        primary_color = HOME_COLORS['primary']
        secondary_color = HOME_COLORS['secondary']
        background = f"linear-gradient(135deg, #9A3412 0%, {primary_color} 50%, {secondary_color} 100%)"
        border_color = '#FDBA74'
    
    white = GLOBAL_COLORS['white']
    
    st.markdown(f"""
    <div class="section-heading-box" style="
        background: {background};
        padding: {HEADING_PADDING};
        border-radius: {HEADING_BORDER_RADIUS};
        margin: 0.4rem 0 0.6rem 0;
        border: 3px solid #1e293b;
        text-align: center;
        box-shadow: 0 3px 10px rgba(0, 0, 0, 0.15);
        width: 100%;
        position: relative;
        overflow: hidden;
    ">
        <div title="{safe_title}" style="
            margin: 0;
            font-size: {HEADING_FONT_SIZE};
            font-weight: 900;
            color: {white};
            text-shadow: 0 1px 3px rgba(0,0,0,0.35);
            letter-spacing: 0.03em;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            position: relative;
            z-index: 1;
        ">
            {safe_title}
        </div>
    </div>
    """, unsafe_allow_html=True)
