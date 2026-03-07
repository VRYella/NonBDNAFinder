"""
Color palette definitions for NBDScanner.

Refined edition: Elegant, scientifically disciplined color scheme.
All variable names and structure preserved.
"""

# ==================== GLOBAL BASE COLORS ====================
GLOBAL_COLORS = {
    # Neutral backgrounds - Light mode foundation
    'neutral_50': '#F6F7FB',
    'neutral_100': '#EEF2F7',
    'neutral_200': '#E4E4E7',
    'neutral_300': '#D4D4D8',
    'neutral_400': '#A1A1AA',
    'neutral_500': '#64748B',
    'neutral_600': '#475569',
    'neutral_700': '#0F172A',
    'neutral_800': '#27272A',
    'neutral_900': '#18181B',

    # Dark mode backgrounds (more atmospheric, less blue-tinted)
    'dark_50': '#18181B',
    'dark_100': '#1F1F23',
    'dark_200': '#27272A',
    'dark_300': '#3F3F46',
    'dark_400': '#52525B',

    'white': '#FFFFFF',
    'black': '#000000',
}

# ==================== GLOBAL NEUTRAL SYSTEM ====================

NEUTRAL_SYSTEM = {
    # Surfaces
    'bg_main': '#F6F7FB',
    'bg_surface': '#FFFFFF',
    'bg_muted': '#EEF2F7',

    # Borders
    'border_subtle': 'rgba(0,0,0,0.08)',
    'border_strong': 'rgba(0,0,0,0.16)',

    # Text
    'text_primary': '#0F172A',
    'text_secondary': '#475569',
    'text_muted': '#64748B',

    # Shadows
    'shadow_soft': '0 6px 16px rgba(0,0,0,0.08)',
}

# ==================== METRIC / SEMANTIC INDICATOR COLORS ====================

METRIC_COLORS = {
    'good': '#059669',
    'warning': '#D97706',
    'bad': '#DC2626',
    'info': '#0284C7',
}

# ==================== PAGE-SPECIFIC ACCENT PALETTES ====================

HOME_COLORS = {
    'primary': '#0284C7',      # Sky Blue
    'primary_light': '#38BDF8',
    'primary_dark': '#0369A1',
    'primary_soft': '#E0F2FE',
    'secondary': '#38BDF8',    # Sky Blue Light
    'accent': '#7DD3FC',       # Sky Blue Accent
    'light': '#F0F9FF',        # Sky-50 light background
    'lighter': '#F7FCFF',      # Near-white sky background
    'medium': '#E0F2FE',       # Sky-100 medium background
    'border': '#BAE6FD',       # Sky border
    'text': '#0C4A6E',         # Deep sky text
    'shadow': 'rgba(2, 132, 199, 0.18)',
}

INPUT_COLORS = {
    'primary': '#4F46E5',      # Royal Indigo
    'primary_light': '#6366F1',
    'primary_dark': '#3730A3',
    'primary_soft': '#E0E7FF',
    'secondary': '#6366F1',    # Bright Indigo
    'accent': '#818CF8',       # Light Indigo accent
    'light': '#EEF2FF',        # Indigo-50 light background
    'lighter': '#F5F7FF',      # Near-white indigo background
    'medium': '#E0E7FF',       # Indigo-100 medium background
    'border': '#C7D2FE',       # Indigo border
    'text': '#1e1b4b',         # Deep navy text
    'shadow': 'rgba(79, 70, 229, 0.18)',
}

ANALYSIS_COLORS = {
    'primary': '#D97706',
    'primary_light': '#F59E0B',
    'primary_dark': '#92400E',
    'primary_soft': '#FEF3C7',
    'secondary': '#F59E0B',
    'accent': '#FBBF24',
    'light': '#FFFBEB',        # Amber-50 light background
    'lighter': '#FFFDF5',      # Near-white amber background
    'medium': '#FEF3C7',       # Amber-100 medium background
    'border': '#FDE68A',
    'text': '#7C2D12',
    'shadow': 'rgba(217, 119, 6, 0.12)',
}

RESULTS_COLORS = {
    'primary': '#0EA5A4',      # Scientific Teal
    'primary_light': '#2DD4BF',
    'primary_dark': '#0F766E',
    'primary_soft': '#CCFBF1',
    'secondary': '#06B6D4',
    'accent': '#22D3EE',
    'light': '#ECFEFF',        # Cyan-50 light background
    'lighter': '#F5FEFF',      # Near-white cyan background
    'medium': '#CFFAFE',       # Cyan-100 medium background
    'border': '#A5F3FC',
    'text': '#164E63',
    'shadow': 'rgba(14, 165, 164, 0.12)',
}

VISUALIZATION_COLORS = {
    'primary': '#0F766E',
    'primary_light': '#14B8A6',
    'primary_dark': '#134E4A',
    'primary_soft': '#CCFBF1',
    'secondary': '#14B8A6',
    'accent': '#5EEAD4',
    'light': '#F0FDFA',        # Teal-50 light background
    'lighter': '#F7FFFD',      # Near-white teal background
    'medium': '#CCFBF1',       # Teal-100 medium background
    'border': '#99F6E4',
    'text': '#134E4A',
    'shadow': 'rgba(15, 118, 110, 0.12)',
}

DOWNLOAD_COLORS = {
    'primary': '#16A34A',      # Emerald Green
    'primary_light': '#22C55E',
    'primary_dark': '#166534',
    'primary_soft': '#DCFCE7',
    'secondary': '#22C55E',
    'accent': '#4ADE80',
    'light': '#ECFDF5',        # Green-50 light background
    'lighter': '#F5FEF8',      # Near-white green background
    'medium': '#DCFCE7',       # Green-100 medium background
    'border': '#BBF7D0',
    'text': '#14532D',
    'shadow': 'rgba(22, 163, 74, 0.12)',
}

DOCUMENTATION_COLORS = {
    'primary': '#7C3AED',      # Violet/Purple
    'primary_light': '#A78BFA',
    'primary_dark': '#6D28D9',
    'primary_soft': '#EDE9FE',
    'secondary': '#A78BFA',    # Bright Violet
    'accent': '#C4B5FD',       # Violet accent
    'light': '#F5F3FF',        # Violet-50 light background
    'lighter': '#FAF9FF',      # Near-white violet background
    'medium': '#EDE9FE',       # Violet-100 medium background
    'border': '#DDD6FE',       # Violet border
    'text': '#4C1D95',         # Deep violet text
    'shadow': 'rgba(124, 58, 237, 0.18)',
}

# ==================== SEMANTIC STATUS COLORS ====================

SEMANTIC_COLORS = {
    'success': '#15803D',
    'success_light': '#ECFDF5',
    'success_dark': '#14532D',
    'success_border': '#86EFAC',
    
    'warning': '#B45309',
    'warning_light': '#FFFBEB',
    'warning_dark': '#78350F',
    'warning_border': '#FCD34D',
    
    'error': '#B91C1C',
    'error_light': '#FEF2F2',
    'error_dark': '#7F1D1D',
    'error_border': '#FCA5A5',
    
    'info': '#1D4ED8',
    'info_light': '#EFF6FF',
    'info_dark': '#1E3A8A',
    'info_border': '#93C5FD',
    
    'progress': '#6D28D9',
    'progress_light': '#F5F3FF',
    'progress_dark': '#4C1D95',
    'progress_border': '#C4B5FD',
}

# ==================== VISUALIZATION COLOR PALETTE ====================

VISUALIZATION_PALETTE = {
    'chart_1': '#D97706',
    'chart_2': '#2563EB',
    'chart_3': '#16A34A',
    'chart_4': '#CA8A04',
    'chart_5': '#1E40AF',
    'chart_6': '#B91C1C',
    'chart_7': '#BE185D',
    'chart_8': '#15803D',
    'chart_9': '#7C3AED',
    'chart_10': '#0F766E',
    'chart_11': '#B45309',
    'chart_12': '#475569',
}

# ==================== UNIFIED MOTIF CLASS COLORS ====================

UNIFIED_MOTIF_COLORS = {
    'Curved_DNA': '#0891B2',
    'Slipped_DNA': '#D97706',
    'Cruciform': '#DC2626',
    'R-Loop': '#7C3AED',
    'Triplex': '#BE185D',
    'G-Quadruplex': '#059669',
    'i-Motif': '#16A34A',
    'Z-DNA': '#4F46E5',
    'A-philic_DNA': '#EA580C',
    'Hybrid': '#64748B',
    'Non-B_DNA_Clusters': '#334155',
}

VISUALIZATION_MOTIF_COLORS = UNIFIED_MOTIF_COLORS.copy()
CLASS_COLORS = UNIFIED_MOTIF_COLORS.copy()

CLASS_COLOR_NAMES = {
    'Curved_DNA': 'blue',
    'Slipped_DNA': 'orange',
    'Cruciform': 'red',
    'R-Loop': 'violet',
    'Triplex': 'violet',
    'G-Quadruplex': 'green',
    'i-Motif': 'green',
    'Z-DNA': 'blue',
    'A-philic_DNA': 'orange',
    'Hybrid': 'gray',
    'Non-B_DNA_Clusters': 'gray',
}

# ==================== MOTIF CLASS CARD COLORS ====================

MOTIF_CARD_COLORS = {
    'Curved_DNA': {'primary': '#0891B2','light': '#ECFEFF','lighter': '#CFFAFE','dark': '#164E63'},
    'Slipped_DNA': {'primary': '#D97706','light': '#FFFBEB','lighter': '#FEF3C7','dark': '#78350F'},
    'Cruciform': {'primary': '#DC2626','light': '#FEF2F2','lighter': '#FEE2E2','dark': '#7F1D1D'},
    'R-Loop': {'primary': '#7C3AED','light': '#F5F3FF','lighter': '#EDE9FE','dark': '#4C1D95'},
    'Triplex': {'primary': '#BE185D','light': '#FDF2F8','lighter': '#FCE7F3','dark': '#831843'},
    'G-Quadruplex': {'primary': '#059669','light': '#ECFDF5','lighter': '#D1FAE5','dark': '#064E3B'},
    'i-Motif': {'primary': '#16A34A','light': '#F0FDF4','lighter': '#DCFCE7','dark': '#14532D'},
    'Z-DNA': {'primary': '#4F46E5','light': '#EEF2FF','lighter': '#E0E7FF','dark': '#312E81'},
    'A-philic_DNA': {'primary': '#EA580C','light': '#FFF7ED','lighter': '#FFEDD5','dark': '#7C2D12'},
    'Hybrid': {'primary': '#F97316','light': '#FFF7ED','lighter': '#FFEDD5','dark': '#7C2D12'},
    'Non-B_DNA_Clusters': {'primary': '#EC4899','light': '#FDF2F8','lighter': '#FCE7F3','dark': '#831843'},
}

# ==================== MOTIF CLASS DISPLAY INFO ====================

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

# ==================== TOP NAVIGATION TAB COLORS ====================

TAB_NAV_COLORS = {
    'Home': {'bg': '#0284C7', 'hover': '#0369A1', 'text': '#E0F2FE', 'glow': 'rgba(2, 132, 199, 0.70)'},
    'Upload & Analyze': {'bg': '#4F46E5', 'hover': '#3730A3', 'text': '#EEF2FF', 'glow': 'rgba(99, 102, 241, 0.65)'},
    'Results': {'bg': '#0EA5A4', 'hover': '#0D9488', 'text': '#FFFFFF', 'glow': 'rgba(34, 211, 238, 0.65)'},
    'Download': {'bg': '#16A34A', 'hover': '#15803D', 'text': '#FFFFFF', 'glow': 'rgba(74, 222, 128, 0.65)'},
    'Documentation': {'bg': '#7C3AED', 'hover': '#6D28D9', 'text': '#F5F3FF', 'glow': 'rgba(167, 139, 250, 0.65)'},
}

# ==================== DOCUMENTATION SUB-TAB COLORS ====================

DOC_SUBTAB_COLORS = {
    'Overview & Architecture': {'bg': '#4F46E5', 'hover': '#4338CA', 'text': '#EEF2FF'},
    'Motif Library & Algorithms': {'bg': '#4F46E5', 'hover': '#4338CA', 'text': '#EEF2FF'},
    'Scoring & Analysis': {'bg': '#0F766E', 'hover': '#14B8A6', 'text': '#F0FDFA'},
    'Statistics Guide': {'bg': '#10B981', 'hover': '#059669', 'text': '#ECFDF5'},
    'References & Citation': {'bg': '#0284C7', 'hover': '#0369A1', 'text': '#E0F2FE'},
}

def get_motif_card_style(class_key: str) -> dict:
    colors = MOTIF_CARD_COLORS.get(class_key, MOTIF_CARD_COLORS['Hybrid'])
    return {
        'background': f"linear-gradient(135deg, {colors['light']} 0%, {colors['lighter']} 100%)",
        'border': colors['primary'],
        'text': colors['dark'],
    }

def get_motif_color(class_key: str) -> str:
    return UNIFIED_MOTIF_COLORS.get(class_key, '#808080')
