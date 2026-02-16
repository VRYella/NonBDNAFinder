"""
Color palette definitions for NBDScanner.

This module contains all color constants used throughout the application:
- GLOBAL_COLORS: Foundation colors for neutral elements
- HOME_COLORS, INPUT_COLORS, etc.: Page-specific accent palettes
- SEMANTIC_COLORS: Status colors (success, warning, error, info, progress)
- VISUALIZATION_PALETTE: Chart and plot colors
"""

# ==================== GLOBAL BASE COLORS ====================
# Foundation colors used throughout the application
# These provide the core visual structure and hierarchy
GLOBAL_COLORS = {
    # Neutral backgrounds - Light mode foundation
    'neutral_50': '#FAFAFA',      # Lightest background for page base
    'neutral_100': '#F5F5F5',     # Very light background for sections
    'neutral_200': '#E5E5E5',     # Light borders and dividers
    'neutral_300': '#D4D4D4',     # Medium borders
    'neutral_400': '#A3A3A3',     # Disabled/muted elements
    'neutral_500': '#737373',     # Secondary text
    'neutral_600': '#525252',     # Body text
    'neutral_700': '#404040',     # Primary text
    'neutral_800': '#262626',     # Headers and emphasis
    'neutral_900': '#171717',     # Maximum contrast text
    
    # Dark mode backgrounds
    'dark_50': '#18181B',         # Darkest background for dark mode
    'dark_100': '#1F1F23',        # Very dark background
    'dark_200': '#27272A',        # Dark card background
    'dark_300': '#3F3F46',        # Dark borders
    'dark_400': '#52525B',        # Dark muted elements
    
    # White and black
    'white': '#FFFFFF',           # Pure white for cards and highlights
    'black': '#000000',           # Pure black for maximum contrast
}

# ==================== PAGE-SPECIFIC ACCENT PALETTES ====================
# Each page/workflow has a distinct color identity while maintaining consistency
# These colors provide subtle visual navigation cues

# HOME / OVERVIEW PAGE - Confident Scientific Blue (Trust, Clarity, Science)
# HIGHLY VIBRANT blue palette for eye-catching professional impact
HOME_COLORS = {
    'primary': '#0091FF',         # ELECTRIC BLUE - ultra vibrant primary actions
    'secondary': '#00B4FF',       # CYAN ELECTRIC - super bright hover states
    'accent': '#66D9FF',          # BRILLIANT SKY - vivid emphasis
    'light': '#CCF2FF',           # BRIGHT CYAN TINT - energetic backgrounds
    'lighter': '#E5F9FF',         # ULTRA LIGHT CYAN - vibrant page base
    'border': '#80E5FF',          # VIVID CYAN borders
    'text': '#003D82',            # DEEP VIBRANT BLUE - strong readability
    'shadow': 'rgba(0, 145, 255, 0.35)',  # STRONGER blue shadow for depth
}

# INPUT / UPLOAD PAGE - Fresh Natural Green (Growth, Initiation, Start)
# HIGHLY VIBRANT green palette for explosive fresh energy
INPUT_COLORS = {
    'primary': '#00E676',         # NEON GREEN - ultra vivid primary actions
    'secondary': '#1DE9B6',       # ELECTRIC MINT - brilliant hover states
    'accent': '#69F0AE',          # BRIGHT LIME - energetic emphasis
    'light': '#B9F6CA',           # VIVID MINT TINT - fresh backgrounds
    'lighter': '#E0FFF4',         # ULTRA LIGHT MINT - vibrant page base
    'border': '#7FFF9F',          # BRILLIANT GREEN borders
    'text': '#00612E',            # DEEP FOREST - strong readability
    'shadow': 'rgba(0, 230, 118, 0.35)',  # STRONGER green shadow for depth
}

# ANALYSIS / COMPUTATION PAGE - Energetic Orange (Energy, Processing, Activity)
# EXPLOSIVE orange palette for maximum energy and dynamism
ANALYSIS_COLORS = {
    'primary': '#FF6D00',         # BLAZING ORANGE - ultra bold primary actions
    'secondary': '#FF9100',       # ELECTRIC GOLD - brilliant hover states
    'accent': '#FFAB00',          # VIVID AMBER - striking emphasis
    'light': '#FFE57F',           # BRIGHT YELLOW TINT - energetic backgrounds
    'lighter': '#FFF9E6',         # ULTRA LIGHT GOLD - vibrant page base
    'border': '#FFCA28',          # BRILLIANT GOLD borders
    'text': '#BF360C',            # DEEP FLAME - strong readability
    'shadow': 'rgba(255, 109, 0, 0.4)',  # INTENSE orange shadow for depth
}

# RESULTS / TABLES PAGE - Refined Vibrant Purple (Insight, Data, Discovery)
# ELECTRIC purple palette for maximum visual impact and insight
RESULTS_COLORS = {
    'primary': '#D500F9',         # ELECTRIC PURPLE - ultra vivid primary actions
    'secondary': '#E040FB',       # NEON MAGENTA - brilliant hover states
    'accent': '#EA80FC',          # BRIGHT VIOLET - striking emphasis
    'light': '#F3E5F5',           # ELEGANT LAVENDER TINT - refined backgrounds
    'lighter': '#FCF2FF',         # ULTRA LIGHT VIOLET - vibrant page base
    'border': '#E1BEE7',          # VIVID LILAC borders
    'text': '#4A148C',            # DEEP ROYAL PURPLE - strong readability
    'shadow': 'rgba(213, 0, 249, 0.35)',  # INTENSE purple shadow for depth
}

# VISUALIZATION / PLOTS PAGE - Clinical Vibrant Teal (Precision, Clarity, Analysis)
# BRILLIANT teal palette for maximum clarity and visual precision
VISUALIZATION_COLORS = {
    'primary': '#00E5FF',         # NEON CYAN - ultra brilliant primary actions
    'secondary': '#18FFFF',       # ELECTRIC AQUA - super bright hover states
    'accent': '#84FFFF',          # VIVID SKY CYAN - striking emphasis
    'light': '#B2FFFF',           # BRIGHT AQUA TINT - energetic backgrounds
    'lighter': '#E0FFFF',         # ULTRA LIGHT CYAN - vibrant page base
    'border': '#76FFE7',          # BRILLIANT CYAN borders
    'text': '#004D5A',            # DEEP TEAL - strong readability
    'shadow': 'rgba(0, 229, 255, 0.4)',  # INTENSE cyan shadow for depth
}

# DOWNLOAD / EXPORT PAGE - Professional Vibrant Indigo (Completion, Authority, Final)
# BOLD indigo palette for commanding authority and completion
DOWNLOAD_COLORS = {
    'primary': '#536DFE',         # ELECTRIC INDIGO - ultra bold primary actions
    'secondary': '#5E72FF',       # NEON PERIWINKLE - brilliant hover states
    'accent': '#8C9EFF',          # BRIGHT LAVENDER BLUE - vivid emphasis
    'light': '#C5CAE9',           # VIVID INDIGO TINT - strong backgrounds
    'lighter': '#E8EAFF',         # ULTRA LIGHT INDIGO - vibrant page base
    'border': '#9FA8DA',          # BRILLIANT IRIS borders
    'text': '#1A237E',            # DEEP NAVY - strong readability
    'shadow': 'rgba(83, 109, 254, 0.4)',  # INTENSE indigo shadow for depth
}

# DOCUMENTATION PAGE - Vibrant Deep Purple Theme (Depth, Reference, Technical)
# Rich purple theme for technical documentation with high contrast
DOCUMENTATION_COLORS = {
    'primary': '#7C4DFF',         # Vivid electric purple - primary actions (Material Deep Purple A200)
    'secondary': '#B388FF',       # Bright lavender - hover states (Material Deep Purple A100)
    'accent': '#D1C4E9',          # Soft periwinkle - emphasis (Material Deep Purple 100)
    'light': '#0E1726',           # Dark background for documentation
    'lighter': '#0B1220',         # Darkest background for page base
    'border': '#1F2937',          # Dark borders
    'text': '#E5E7EB',            # Light text for dark background
    'shadow': 'rgba(0, 0, 0, 0.5)',  # Strong shadow for dark mode depth
}

# ==================== SEMANTIC STATUS COLORS ====================
# Consistent meaning across all pages - Universal visual language (ULTRA VIBRANT)
# EXPLOSIVE vibrant colors for immediate, unmistakable visual feedback
SEMANTIC_COLORS = {
    # Success states - Positive outcomes, completion, validation
    'success': '#00E676',         # NEON SUCCESS GREEN - ultra vivid
    'success_light': '#B9F6CA',   # BRIGHT success background
    'success_dark': '#00612E',    # DEEP success text
    'success_border': '#69F0AE',  # BRILLIANT success border
    
    # Warning states - Caution, important notices, attention needed
    'warning': '#FF9100',         # BLAZING WARNING ORANGE - maximum attention
    'warning_light': '#FFE57F',   # VIVID warning background
    'warning_dark': '#BF360C',    # DEEP warning text
    'warning_border': '#FFAB00',  # ELECTRIC warning border
    
    # Error states - Problems, failures, invalid inputs
    'error': '#FF1744',           # NEON ERROR RED - ultra striking
    'error_light': '#FFCDD2',     # BRIGHT error background
    'error_dark': '#B71C1C',      # DEEP error text
    'error_border': '#FF5252',    # BRILLIANT error border
    
    # Info states - Neutral information, tips, explanations
    'info': '#00B4FF',            # ELECTRIC INFO BLUE - ultra vivid
    'info_light': '#CCF2FF',      # BRIGHT info background
    'info_dark': '#003D82',       # DEEP info text
    'info_border': '#66D9FF',     # BRILLIANT info border
    
    # Progress states - Processing, loading, intermediate states
    'progress': '#E040FB',        # NEON PROGRESS PURPLE - ultra striking
    'progress_light': '#F3E5F5',  # Light progress background
    'progress_dark': '#4A148C',   # DEEP progress text
    'progress_border': '#EA80FC', # BRILLIANT progress border
}

# ==================== VISUALIZATION COLOR PALETTE ====================
# Scientific color scheme for charts and plots - ULTRA VIBRANT & ACCESSIBLE
# EXPLOSIVE colorblind-friendly palette with maximum saturation and visual impact
VISUALIZATION_PALETTE = {
    'chart_1': '#FF6D00',         # BLAZING ORANGE - Ultra high contrast
    'chart_2': '#0091FF',         # ELECTRIC BLUE - Ultra distinct
    'chart_3': '#00E676',         # NEON GREEN - Ultra clear
    'chart_4': '#FFEA00',         # BRILLIANT YELLOW - Ultra striking
    'chart_5': '#0043A8',         # DEEP ELECTRIC BLUE - Ultra professional
    'chart_6': '#FF1744',         # NEON RED - Ultra energetic
    'chart_7': '#FF00AA',         # ELECTRIC PINK - Ultra unique
    'chart_8': '#76FF03',         # VIVID LIME - Ultra natural
    'chart_9': '#D500F9',         # ELECTRIC PURPLE - Ultra elegant
    'chart_10': '#00E5FF',        # NEON CYAN - Ultra fresh
    'chart_11': '#FFC400',        # GOLD FLASH - Ultra warm
    'chart_12': '#546E7A',        # STEEL BLUE - Balanced contrast
}

# ==================== UNIFIED MOTIF CLASS COLORS ====================
# SINGLE SOURCE OF TRUTH for all motif class colors across the entire application.
# 
# This unified palette is used for:
# - Scientific visualizations (plots, charts, tracks)
# - UI elements (cards, checkboxes, headers)
# - Documentation and home page displays
#
# Design principles:
# 1. Each of the 11 classes has a UNIQUE, distinct color
# 2. Colorblind-friendly (tested with protanopia, deuteranopia, tritanopia)
# 3. High contrast for both screen and print
# 4. Inspired by Wong, B. (2011) Nature Methods colorblind-safe palette
#    Reference: Wong, B. (2011) Points of view: Color blindness. Nat Methods 8, 441
#
# Non-B DNA Finder detects 11 major classes and 24 subclasses:
# 1. Curved DNA (Global Curvature, Local Curvature)
# 2. Slipped DNA (Direct Repeats, Short Tandem Repeats)
# 3. Cruciform (Cruciform forming IRs)
# 4. R-Loop (R-loop formation sites)
# 5. Triplex/H-DNA (Triplex, Sticky DNA)
# 6. G-Quadruplex (8 subclasses including Telomeric G4, Canonical G4, etc.)
# 7. i-Motif (Canonical i-motif, Relaxed i-motif, AC-motif)
# 8. Z-DNA (Z-DNA, eGZ motif)
# 9. A-philic DNA (A-philic DNA)
# 10. Hybrid (Dynamic overlaps)
# 11. Non-B_DNA_Clusters (Dynamic clusters)

UNIFIED_MOTIF_COLORS = {
    # ---- Primary Structural Classes ----
    'Curved_DNA': '#06b6d4',          # Cyan - DNA bending/curvature
    'Slipped_DNA': '#f59e0b',         # Amber - Direct repeats & STRs
    'Cruciform': '#ef4444',           # Red - Inverted repeats/hairpins
    'R-Loop': '#8b5cf6',              # Violet - RNA:DNA hybrids
    'Triplex': '#ec4899',             # Pink/Magenta - Triple-stranded H-DNA
    
    # ---- G-rich & C-rich Structures ----
    'G-Quadruplex': '#10b981',        # Emerald - G4 structures (8 subclasses)
    'i-Motif': '#22c55e',             # Green - C-rich intercalated structures
    
    # ---- Alternative Helices ----
    'Z-DNA': '#6366f1',               # Indigo - Left-handed Z-form helix
    'A-philic_DNA': '#f97316',        # Orange - A-form propensity regions
    
    # ---- Meta-classes (Overlaps & Clusters) ----
    'Hybrid': '#64748b',              # Slate - Multi-class overlapping regions
    'Non-B_DNA_Clusters': '#334155',  # Dark Slate - High-density Non-B hotspots
}

# ==================== VISUALIZATION MOTIF COLORS (PUBLICATION-GRADE) ====================
# Re-exports UNIFIED_MOTIF_COLORS for backward compatibility.
# Used by utilities.py, visualization/standards.py, and all plotting functions.
#
# NOTE: This is now an alias to UNIFIED_MOTIF_COLORS to ensure consistency.
VISUALIZATION_MOTIF_COLORS = UNIFIED_MOTIF_COLORS.copy()

# ==================== MOTIF CLASS COLOR ENCODING (UI) ====================
# Re-exports UNIFIED_MOTIF_COLORS for backward compatibility.
# Used by UI components (upload.py checkboxes, home.py cards, etc.)
#
# NOTE: This is now an alias to UNIFIED_MOTIF_COLORS to ensure consistency.
CLASS_COLORS = UNIFIED_MOTIF_COLORS.copy()

# Streamlit color names for :color[text] syntax (matches CLASS_COLORS)
# Used for colorizing submotif checkbox labels in the UI
CLASS_COLOR_NAMES = {
    'Curved_DNA': 'blue',             # Cyan → blue (closest Streamlit color)
    'Slipped_DNA': 'orange',          # Amber → orange
    'Cruciform': 'red',               # Red → red
    'R-Loop': 'violet',               # Violet → violet
    'Triplex': 'violet',              # Pink → violet (closest)
    'G-Quadruplex': 'green',          # Emerald → green
    'i-Motif': 'green',               # Green → green
    'Z-DNA': 'blue',                  # Indigo → blue
    'A-philic_DNA': 'orange',         # Orange → orange
    'Hybrid': 'gray',                 # Slate → gray
    'Non-B_DNA_Clusters': 'gray',     # Dark slate → gray
}

# ==================== MOTIF CLASS CARD COLORS ====================
# Pre-computed light/dark color variations for UI cards.
# Each class has:
#   - 'primary': The main UNIFIED_MOTIF_COLORS color (for border)
#   - 'light': Light background gradient start
#   - 'lighter': Light background gradient end  
#   - 'dark': Dark text color
# These ensure uniform styling across home.py and documentation.py
MOTIF_CARD_COLORS = {
    'Curved_DNA': {
        'primary': '#06b6d4',   # Cyan
        'light': '#cffafe',     # cyan-100
        'lighter': '#a5f3fc',   # cyan-200
        'dark': '#164e63',      # cyan-900
    },
    'Slipped_DNA': {
        'primary': '#f59e0b',   # Amber
        'light': '#fef3c7',     # amber-100
        'lighter': '#fde68a',   # amber-200
        'dark': '#78350f',      # amber-900
    },
    'Cruciform': {
        'primary': '#ef4444',   # Red
        'light': '#fee2e2',     # red-100
        'lighter': '#fecaca',   # red-200
        'dark': '#7f1d1d',      # red-900
    },
    'R-Loop': {
        'primary': '#8b5cf6',   # Violet
        'light': '#ede9fe',     # violet-100
        'lighter': '#ddd6fe',   # violet-200
        'dark': '#4c1d95',      # violet-900
    },
    'Triplex': {
        'primary': '#ec4899',   # Pink
        'light': '#fce7f3',     # pink-100
        'lighter': '#fbcfe8',   # pink-200
        'dark': '#831843',      # pink-900
    },
    'G-Quadruplex': {
        'primary': '#10b981',   # Emerald
        'light': '#d1fae5',     # emerald-100
        'lighter': '#a7f3d0',   # emerald-200
        'dark': '#064e3b',      # emerald-900
    },
    'i-Motif': {
        'primary': '#22c55e',   # Green
        'light': '#dcfce7',     # green-100
        'lighter': '#bbf7d0',   # green-200
        'dark': '#14532d',      # green-900
    },
    'Z-DNA': {
        'primary': '#6366f1',   # Indigo
        'light': '#e0e7ff',     # indigo-100
        'lighter': '#c7d2fe',   # indigo-200
        'dark': '#312e81',      # indigo-900
    },
    'A-philic_DNA': {
        'primary': '#f97316',   # Orange
        'light': '#ffedd5',     # orange-100
        'lighter': '#fed7aa',   # orange-200
        'dark': '#7c2d12',      # orange-900
    },
    'Hybrid': {
        'primary': '#64748b',   # Slate
        'light': '#f1f5f9',     # slate-100
        'lighter': '#e2e8f0',   # slate-200
        'dark': '#1e293b',      # slate-800
    },
    'Non-B_DNA_Clusters': {
        'primary': '#334155',   # Dark Slate
        'light': '#e2e8f0',     # slate-200
        'lighter': '#cbd5e1',   # slate-300
        'dark': '#0f172a',      # slate-900
    },
}

# ==================== MOTIF CLASS DISPLAY INFO ====================
# Complete display information for all 11 motif classes.
# Used by home.py and documentation.py for consistent display.
MOTIF_CLASS_INFO = [
    {'key': 'Curved_DNA', 'name': 'Curved DNA', 'subtitle': 'A-tract curvature', 'num': 1},
    {'key': 'Slipped_DNA', 'name': 'Slipped DNA', 'subtitle': 'Direct repeats, STRs', 'num': 2},
    {'key': 'Cruciform', 'name': 'Cruciform', 'subtitle': 'Inverted repeats', 'num': 3},
    {'key': 'R-Loop', 'name': 'R-Loop', 'subtitle': 'RNA-DNA hybrids', 'num': 4},
    {'key': 'Triplex', 'name': 'Triplex/H-DNA', 'subtitle': 'Mirror repeats, Sticky DNA', 'num': 5},
    {'key': 'G-Quadruplex', 'name': 'G-Quadruplex', 'subtitle': '8 subclasses', 'num': 6},
    {'key': 'i-Motif', 'name': 'i-Motif', 'subtitle': 'C-rich structures', 'num': 7},
    {'key': 'Z-DNA', 'name': 'Z-DNA', 'subtitle': 'Left-handed helix', 'num': 8},
    {'key': 'A-philic_DNA', 'name': 'A-philic DNA', 'subtitle': 'A/T-rich regions', 'num': 9},
    {'key': 'Hybrid', 'name': 'Hybrid', 'subtitle': 'Multi-class overlap', 'num': 10},
    {'key': 'Non-B_DNA_Clusters', 'name': 'Clusters', 'subtitle': 'Motif hotspots', 'num': 11},
]


def get_motif_card_style(class_key: str) -> dict:
    """
    Get CSS style properties for a motif class card.
    
    Args:
        class_key: The motif class key (e.g., 'Curved_DNA', 'G-Quadruplex')
        
    Returns:
        Dictionary with 'background', 'border', 'text' CSS color values
    """
    colors = MOTIF_CARD_COLORS.get(class_key, MOTIF_CARD_COLORS['Hybrid'])
    return {
        'background': f"linear-gradient(135deg, {colors['light']} 0%, {colors['lighter']} 100%)",
        'border': colors['primary'],
        'text': colors['dark'],
    }


def get_motif_color(class_key: str) -> str:
    """
    Get the primary color for a motif class.
    
    Args:
        class_key: The motif class key (e.g., 'Curved_DNA', 'G-Quadruplex')
        
    Returns:
        Hex color string
    """
    return UNIFIED_MOTIF_COLORS.get(class_key, '#808080')
