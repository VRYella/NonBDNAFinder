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
    'primary': '#AC0808',      # Deep Crimson Red
    'primary_light': '#D32F2F',
    'primary_dark': '#7B0000',
    'primary_soft': '#FFEBEE',
    'secondary': '#C62828',    # Dark Red
    'accent': '#EF9A9A',       # Light Red accent
    'light': '#FFF5F5',        # Red-50 light background
    'lighter': '#FFF9F9',      # Near-white red background
    'medium': '#FFEBEE',       # Red-100 medium background
    'border': '#FFCDD2',       # Red border
    'text': '#7B0000',         # Deep red text
    'shadow': 'rgba(172, 8, 8, 0.18)',
}

INPUT_COLORS = {
    'primary': '#6B21A8',      # Deep Scientific Purple
    'primary_light': '#9333EA',
    'primary_dark': '#581C87',
    'primary_soft': '#F3E8FF',
    'secondary': '#9333EA',    # Bright Purple
    'accent': '#C084FC',       # Light Purple accent
    'light': '#FAF5FF',        # Purple-50 light background
    'lighter': '#FDFAFF',      # Near-white purple background
    'medium': '#F3E8FF',       # Purple-100 medium background
    'border': '#E9D5FF',       # Purple border
    'text': '#3B0764',         # Deep violet text
    'shadow': 'rgba(107, 33, 168, 0.18)',
}

ANALYSIS_COLORS = {
    'primary': '#7E22CE',      # Medium Scientific Purple (analysis variant)
    'primary_light': '#A855F7',
    'primary_dark': '#6B21A8',
    'primary_soft': '#F3E8FF',
    'secondary': '#A855F7',
    'accent': '#D8B4FE',
    'light': '#FAF5FF',        # Purple-50 light background
    'lighter': '#FDFAFF',      # Near-white purple background
    'medium': '#F3E8FF',       # Purple-100 medium background
    'border': '#E9D5FF',
    'text': '#4C1D95',
    'shadow': 'rgba(126, 34, 206, 0.12)',
}

RESULTS_COLORS = {
    'primary': '#15803D',      # Forest Emerald Green
    'primary_light': '#22C55E',
    'primary_dark': '#14532D',
    'primary_soft': '#DCFCE7',
    'secondary': '#22C55E',    # Bright Green
    'accent': '#4ADE80',       # Light Green accent
    'light': '#F0FDF4',        # Green-50 light background
    'lighter': '#F7FEFA',      # Near-white green background
    'medium': '#DCFCE7',       # Green-100 medium background
    'border': '#BBF7D0',       # Green border
    'text': '#052E16',         # Deep forest text
    'shadow': 'rgba(21, 128, 61, 0.18)',
}

VISUALIZATION_COLORS = {
    'primary': '#166534',      # Deep Forest Green (visualization variant)
    'primary_light': '#15803D',
    'primary_dark': '#14532D',
    'primary_soft': '#DCFCE7',
    'secondary': '#15803D',
    'accent': '#4ADE80',
    'light': '#F0FDF4',        # Green-50 light background
    'lighter': '#F7FEFA',      # Near-white green background
    'medium': '#DCFCE7',       # Green-100 medium background
    'border': '#BBF7D0',
    'text': '#052E16',
    'shadow': 'rgba(22, 101, 52, 0.12)',
}

DOWNLOAD_COLORS = {
    'primary': '#C68642',      # Sandy Brown
    'primary_light': '#D4935E',
    'primary_dark': '#A0632A',
    'primary_soft': '#FEF9F0',
    'secondary': '#D4935E',    # Medium Sandy Brown
    'accent': '#F4A460',       # Sandy Brown accent
    'light': '#FEF6E8',        # Sandy-50 light background
    'lighter': '#FFFDF8',      # Near-white sandy background
    'medium': '#FDF0DC',       # Sandy-100 medium background
    'border': '#F7D9B0',       # Sandy border
    'text': '#5C3A1E',         # Deep sandy text
    'shadow': 'rgba(198, 134, 66, 0.18)',
}

DOCUMENTATION_COLORS = {
    'primary': '#0EA5E9',      # Sky Blue
    'primary_light': '#38BDF8',
    'primary_dark': '#0369A1',
    'primary_soft': '#F0F9FF',
    'secondary': '#38BDF8',    # Bright Sky Blue
    'accent': '#7DD3FC',       # Light Sky Blue accent
    'light': '#F0F9FF',        # Sky-50 light background
    'lighter': '#F8FBFF',      # Near-white sky background
    'medium': '#E0F2FE',       # Sky-100 medium background
    'border': '#BAE6FD',       # Sky border
    'text': '#0C4A6E',         # Deep sky text
    'shadow': 'rgba(14, 165, 233, 0.18)',
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
    'Home': {'bg': '#AC0808', 'hover': '#7B0000', 'text': '#FFEBEE', 'glow': 'rgba(172, 8, 8, 0.70)'},
    'Upload & Analyze': {'bg': '#6B21A8', 'hover': '#581C87', 'text': '#F3E8FF', 'glow': 'rgba(107, 33, 168, 0.65)'},
    'Results': {'bg': '#15803D', 'hover': '#14532D', 'text': '#FFFFFF', 'glow': 'rgba(21, 128, 61, 0.65)'},
    'Download': {'bg': '#C68642', 'hover': '#A0632A', 'text': '#FEF9F0', 'glow': 'rgba(198, 134, 66, 0.65)'},
    'Documentation': {'bg': '#0EA5E9', 'hover': '#0369A1', 'text': '#F0F9FF', 'glow': 'rgba(14, 165, 233, 0.65)'},
}

# ==================== DOCUMENTATION SUB-TAB COLORS ====================

DOC_SUBTAB_COLORS = {
    'Overview & Architecture': {'bg': '#AC0808', 'hover': '#7B0000', 'text': '#FFEBEE'},
    'Motif Library & Algorithms': {'bg': '#6B21A8', 'hover': '#581C87', 'text': '#F3E8FF'},
    'Scoring & Analysis': {'bg': '#15803D', 'hover': '#14532D', 'text': '#F0FDF4'},
    'Statistics Guide': {'bg': '#B45309', 'hover': '#92400E', 'text': '#FFFBEB'},
    'References & Citation': {'bg': '#0EA5E9', 'hover': '#0369A1', 'text': '#F0F9FF'},
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
