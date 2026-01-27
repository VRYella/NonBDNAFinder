"""
CSS and styling utilities for NBDScanner.

This module contains functions for:
- Loading CSS with theme variables
- Color conversion utilities
- SVG pattern generation
- Page color retrieval
"""

import os
import streamlit as st

from config.colors import (
    GLOBAL_COLORS,
    HOME_COLORS,
    INPUT_COLORS,
    ANALYSIS_COLORS,
    RESULTS_COLORS,
    VISUALIZATION_COLORS,
    DOWNLOAD_COLORS,
    DOCUMENTATION_COLORS,
    SEMANTIC_COLORS,
)
from config.themes import COLOR_THEMES
from config.typography import FONT_CONFIG
from config.layout import LAYOUT_CONFIG
from config.animation import ANIMATION_CONFIG


def hex_to_rgb(hex_color: str) -> tuple:
    """Convert hex color to RGB tuple for CSS rgba() usage."""
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))


def get_dna_pattern_svg(stroke_color: str) -> str:
    """Generate subtle DNA helix SVG pattern for background."""
    # Compact and URL-encoded SVG for background pattern
    svg = (
        "%3Csvg xmlns='http://www.w3.org/2000/svg' width='60' height='60' viewBox='0 0 60 60'%3E"
        f"%3Cg fill='none' stroke='%23{stroke_color}' stroke-width='0.6' opacity='0.12'%3E"
        "%3Cpath d='M10 30 C 18 12, 42 12, 50 30'/%3E"
        "%3Cpath d='M10 30 C 18 48, 42 48, 50 30'/%3E"
        "%3C/g%3E%3C/svg%3E"
    )
    return f"url(\"data:image/svg+xml,{svg}\")"


def load_css(theme_name=None):
    """
    Load external CSS file and inject dynamic theme variables.
    This makes the app.py file more succinct by separating styling concerns.
    If theme_name provided, use per-page theme, otherwise use session color_theme.
    """
    theme_to_use = COLOR_THEMES.get(theme_name, COLOR_THEMES.get(st.session_state.color_theme, COLOR_THEMES['scientific_blue']))
    is_dark = st.session_state.theme_mode == 'dark'
    
    # Read the external CSS file if present
    css_file_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'styles.css')
    try:
        with open(css_file_path, 'r') as f:
            css_content = f.read()
    except Exception:
        css_content = ""
    
    # Compute some derived values
    p_rgb = hex_to_rgb(theme_to_use['primary'])
    s_rgb = hex_to_rgb(theme_to_use['secondary'])
    tab_bg_color_local = theme_to_use.get('tab_bg', theme_to_use['bg_card'])
    tab_active_color_local = theme_to_use.get('tab_active', theme_to_use['primary'])
    dna_svg_local = get_dna_pattern_svg('1e3a5f' if is_dark else 'bbdefb')
    
    # Inject centralized configuration into CSS variables
    theme_vars = f"""
    <style>
    :root {{
        /* Theme Colors */
        --primary-color: {theme_to_use['primary']};
        --secondary-color: {theme_to_use['secondary']};
        --accent-color: {theme_to_use['accent']};
        --bg-light: {theme_to_use['bg_light']};
        --bg-card: {theme_to_use['bg_card']};
        --text-color: {theme_to_use['text']};
        --tab-bg: {tab_bg_color_local};
        --tab-active: {tab_active_color_local};
        --shadow-color: {theme_to_use['shadow']};
        --primary-rgb: {p_rgb[0]}, {p_rgb[1]}, {p_rgb[2]};
        --secondary-rgb: {s_rgb[0]}, {s_rgb[1]}, {s_rgb[2]};
        
        /* Typography from FONT_CONFIG */
        --font-primary: {FONT_CONFIG['primary_font']};
        --font-monospace: {FONT_CONFIG['monospace_font']};
        --font-h1: {FONT_CONFIG['h1_size']};
        --font-h2: {FONT_CONFIG['h2_size']};
        --font-h3: {FONT_CONFIG['h3_size']};
        --font-h4: {FONT_CONFIG['h4_size']};
        --font-body: {FONT_CONFIG['body_size']};
        --font-small: {FONT_CONFIG['small_size']};
        --font-caption: {FONT_CONFIG['caption_size']};
        --font-weight-light: {FONT_CONFIG['light_weight']};
        --font-weight-normal: {FONT_CONFIG['normal_weight']};
        --font-weight-medium: {FONT_CONFIG['medium_weight']};
        --font-weight-semibold: {FONT_CONFIG['semibold_weight']};
        --font-weight-bold: {FONT_CONFIG['bold_weight']};
        --font-weight-extrabold: {FONT_CONFIG['extrabold_weight']};
        
        /* Layout from LAYOUT_CONFIG */
        --border-radius-sm: {LAYOUT_CONFIG['border_radius']['small']};
        --border-radius-md: {LAYOUT_CONFIG['border_radius']['medium']};
        --border-radius-lg: {LAYOUT_CONFIG['border_radius']['large']};
        --border-radius-pill: {LAYOUT_CONFIG['border_radius']['pill']};
        --spacing-xs: 0.25rem;
        --spacing-sm: {LAYOUT_CONFIG['padding']['small']};
        --spacing-md: {LAYOUT_CONFIG['padding']['medium']};
        --spacing-lg: {LAYOUT_CONFIG['padding']['large']};
        --spacing-xl: {LAYOUT_CONFIG['padding']['xlarge']};
        --spacing-2xl: 2.5rem;
        
        /* Transitions from ANIMATION_CONFIG */
        --transition-fast: {ANIMATION_CONFIG['transition_fast']} {ANIMATION_CONFIG['easing_smooth']};
        --transition-normal: {ANIMATION_CONFIG['transition_normal']} {ANIMATION_CONFIG['easing_smooth']};
        --transition-slow: {ANIMATION_CONFIG['transition_slow']} {ANIMATION_CONFIG['easing_smooth']};
        
        /* Theme State */
        --dark-mode: {1 if is_dark else 0};
        --dna-pattern: {dna_svg_local};
    }}
    {css_content}
    </style>
    """
    
    st.markdown(theme_vars, unsafe_allow_html=True)


def get_page_colors(page_name='Home'):
    """
    Get color dictionary for inline HTML styles based on page context.
    
    NOTE: Inline HTML styles require literal color values and cannot use CSS variables.
    This function returns the current page's colors from the centralized token system
    for use in inline styles. All values are derived from the token block at the top.
    
    Args:
        page_name: Name of the page ('Home', 'Upload & Analyze', 'Results', etc.)
    
    Returns:
        Dictionary with page-specific colors from centralized tokens
    """
    # Map page names to their color palettes from centralized tokens
    page_color_map = {
        'Home': HOME_COLORS,
        'Upload & Analyze': INPUT_COLORS,
        'Analysis': ANALYSIS_COLORS,
        'Results': RESULTS_COLORS,
        'Visualization': VISUALIZATION_COLORS,
        'Download': DOWNLOAD_COLORS,
        'Documentation': DOCUMENTATION_COLORS,
    }
    
    # Get the page-specific palette
    page_palette = page_color_map.get(page_name, HOME_COLORS)
    
    # Return comprehensive color set combining page colors, global colors, and semantic colors
    return {
        # Page-specific colors
        **page_palette,
        # Global colors for consistent elements
        'white': GLOBAL_COLORS['white'],
        'neutral_50': GLOBAL_COLORS['neutral_50'],
        'neutral_100': GLOBAL_COLORS['neutral_100'],
        'neutral_200': GLOBAL_COLORS['neutral_200'],
        'neutral_500': GLOBAL_COLORS['neutral_500'],
        'neutral_600': GLOBAL_COLORS['neutral_600'],
        'neutral_700': GLOBAL_COLORS['neutral_700'],
        # Semantic colors for status indicators
        'success': SEMANTIC_COLORS['success'],
        'warning': SEMANTIC_COLORS['warning'],
        'error': SEMANTIC_COLORS['error'],
        'info': SEMANTIC_COLORS['info'],
    }
