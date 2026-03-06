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
    'Home': {'primary': '#2C3E50', 'border': '#E5E7EB', 'pattern': 'transparent'},
    'Upload & Analyze': {'primary': '#2C3E50', 'border': '#E5E7EB', 'pattern': 'transparent'},
    'Results': {'primary': '#2C3E50', 'border': '#E5E7EB', 'pattern': 'transparent'},
    'Downloads': {'primary': '#2C3E50', 'border': '#E5E7EB', 'pattern': 'transparent'},
    'Documentation': {'primary': '#2C3E50', 'border': '#E5E7EB', 'pattern': 'transparent'},
}
HEADING_PADDING = "0.25rem 0 0.75rem 0"; HEADING_FONT_SIZE = "2.0rem"; HEADING_BORDER_RADIUS = "0"
# ═══════════════════════════════════════════════════════════════════════════════


def render_section_heading(title: str, page: str = None):
    """Render a uniform section heading with a clean scientific style.
    
    This is the canonical heading component for NonBDNAFinder. All section
    headings should use this function to maintain visual consistency.
    
    Features:
    - Clean bottom border separator (no heavy gradient banner)
    - Consistent typography aligned with the scientific theme
    - Primary color (#2C3E50) text for strong visual hierarchy
    - Subtitle text supported via optional description
    
    Args:
        title: The heading text to display
        page: Optional page name for page-specific colors (Home, Upload & Analyze, Results, Downloads, Documentation)
    """
    # Escape HTML to prevent XSS
    safe_title = html.escape(title)
    
    # Scientific theme colors
    primary_color = '#2C3E50'
    border_color = '#E5E7EB'
    
    st.markdown(f"""
    <div class="section-heading-box" style="
        background: transparent;
        padding: {HEADING_PADDING};
        border-radius: {HEADING_BORDER_RADIUS};
        margin: 0.4rem 0 0.8rem 0;
        border: none;
        border-bottom: 2px solid {border_color};
        text-align: left;
        box-shadow: none;
        width: 100%;
        position: relative;
    ">
        <div title="{safe_title}" style="
            margin: 0;
            font-size: {HEADING_FONT_SIZE};
            font-weight: 600;
            color: {primary_color};
            text-shadow: none;
            letter-spacing: -0.01em;
            white-space: normal;
            overflow: hidden;
            text-overflow: ellipsis;
            position: relative;
            z-index: 1;
        ">
            {safe_title}
        </div>
    </div>
    """, unsafe_allow_html=True)
