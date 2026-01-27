"""
Theme configuration for NBDScanner.

This module defines color themes and tab themes:
- COLOR_THEMES: Different color themes (scientific_blue, nature_green, etc.)
- TAB_THEMES: Maps each tab to a specific theme
"""

from config.colors import (
    GLOBAL_COLORS,
    HOME_COLORS,
    INPUT_COLORS,
    RESULTS_COLORS,
    VISUALIZATION_COLORS,
    DOCUMENTATION_COLORS,
)

# ==================== COLOR THEMES ====================
# Define color themes for different moods and contexts - VIBRANT EDITION
# Each theme has primary, secondary, accent, backgrounds, text, and shadow colors
# All values now reference the centralized VIBRANT color tokens defined above
# IMPORTANT: No hardcoded hex values - all colors come from token system
COLOR_THEMES = {
    'scientific_blue': {
        'primary': HOME_COLORS['primary'],
        'secondary': HOME_COLORS['secondary'],
        'accent': HOME_COLORS['accent'],
        'bg_light': HOME_COLORS['lighter'],
        'bg_card': HOME_COLORS['light'],
        'text': HOME_COLORS['text'],
        'tab_bg': GLOBAL_COLORS['neutral_100'],
        'tab_active': HOME_COLORS['primary'],
        'shadow': 'rgba(0, 145, 255, 0.35)'  # VIBRANT blue shadow
    },
    'nature_green': {
        'primary': INPUT_COLORS['primary'],
        'secondary': INPUT_COLORS['secondary'],
        'accent': INPUT_COLORS['accent'],
        'bg_light': INPUT_COLORS['lighter'],
        'bg_card': INPUT_COLORS['light'],
        'text': INPUT_COLORS['text'],
        'tab_bg': GLOBAL_COLORS['neutral_100'],
        'tab_active': INPUT_COLORS['primary'],
        'shadow': 'rgba(0, 230, 118, 0.35)'  # VIBRANT green shadow
    },
    'genomic_purple': {
        'primary': RESULTS_COLORS['primary'],
        'secondary': RESULTS_COLORS['secondary'],
        'accent': RESULTS_COLORS['accent'],
        'bg_light': RESULTS_COLORS['lighter'],
        'bg_card': RESULTS_COLORS['light'],
        'text': RESULTS_COLORS['text'],
        'tab_bg': GLOBAL_COLORS['neutral_100'],
        'tab_active': RESULTS_COLORS['primary'],
        'shadow': 'rgba(213, 0, 249, 0.35)'  # VIBRANT purple shadow
    },
    'clinical_teal': {
        'primary': VISUALIZATION_COLORS['primary'],
        'secondary': VISUALIZATION_COLORS['secondary'],
        'accent': VISUALIZATION_COLORS['accent'],
        'bg_light': VISUALIZATION_COLORS['lighter'],
        'bg_card': VISUALIZATION_COLORS['light'],
        'text': VISUALIZATION_COLORS['text'],
        'tab_bg': GLOBAL_COLORS['neutral_100'],
        'tab_active': VISUALIZATION_COLORS['primary'],
        'shadow': 'rgba(0, 229, 255, 0.4)'  # VIBRANT cyan shadow
    },
    'midnight': {
        'primary': DOCUMENTATION_COLORS['primary'],
        'secondary': DOCUMENTATION_COLORS['secondary'],
        'accent': DOCUMENTATION_COLORS['accent'],
        'bg_light': DOCUMENTATION_COLORS['lighter'],
        'bg_card': DOCUMENTATION_COLORS['light'],
        'text': DOCUMENTATION_COLORS['text'],
        'tab_bg': GLOBAL_COLORS['dark_100'],
        'tab_active': DOCUMENTATION_COLORS['primary'],
        'shadow': 'rgba(0, 0, 0, 0.5)'  # Strong shadow for dark mode depth
    }
}

# ==================== PAGE THEMES PER TAB ====================
# Assign a specific theme to each tab for visual distinction
# This creates a unique color experience for each section of the app
# Change values to any theme name defined in COLOR_THEMES above
TAB_THEMES = {
    'Home': 'scientific_blue',          # Homepage uses professional blue
    'Upload & Analyze': 'nature_green',  # Upload tab uses natural green
    'Results': 'genomic_purple',        # Results tab uses genomic purple
    'Download': 'clinical_teal',        # Download tab uses clinical teal
    'Documentation': 'midnight'         # Documentation uses dark midnight theme
}
