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
    'Home': {'primary': '#EA580C', 'border': '#FDBA74', 'pattern': 'linear-gradient(135deg, #9A3412 0%, #EA580C 40%, #FB923C 75%, #FBBF24 100%)'},
    'Upload & Analyze': {'primary': '#CA8A04', 'border': '#FEF08A', 'pattern': 'linear-gradient(135deg, #713F12 0%, #CA8A04 40%, #EAB308 75%, #FACC15 100%)'},
    'Results': {'primary': '#0EA5A4', 'border': '#67e8f9', 'pattern': 'linear-gradient(135deg, #0891B2 0%, #0EA5A4 40%, #06B6D4 75%, #22D3EE 100%)'},
    'Downloads': {'primary': '#16A34A', 'border': '#86efac', 'pattern': 'linear-gradient(135deg, #16A34A 0%, #22C55E 40%, #4ADE80 75%, #86EFAC 100%)'},
    'Documentation': {'primary': '#DC2626', 'border': '#fca5a5', 'pattern': 'linear-gradient(135deg, #B91C1C 0%, #DC2626 40%, #EF4444 75%, #FCA5A5 100%)'},
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
        border: 4px solid {border_color};
        text-align: center;
        box-shadow: 0 6px 28px rgba(0, 0, 0, 0.35), 0 0 20px rgba(255,255,255,0.08) inset, 0 0 0 1px rgba(255,255,255,0.1);
        width: 100%;
        position: relative;
        overflow: hidden;
    ">
        <div title="{safe_title}" style="
            margin: 0;
            font-size: {HEADING_FONT_SIZE};
            font-weight: 900;
            color: {white};
            text-shadow: 0 0 18px rgba(255,255,255,0.95), 0 0 35px rgba(255,255,255,0.6), 0 2px 4px rgba(0,0,0,0.4);
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
