"""
Color palette definitions for NBDScanner.

Refined edition: Elegant, publication-grade, scientifically disciplined.
All variable names and structure preserved.
"""

# ==================== GLOBAL BASE COLORS ====================
GLOBAL_COLORS = {
    # Neutral backgrounds - Light mode foundation
    'neutral_50': '#FAFAF9',
    'neutral_100': '#F4F4F5',
    'neutral_200': '#E4E4E7',
    'neutral_300': '#D4D4D8',
    'neutral_400': '#A1A1AA',
    'neutral_500': '#71717A',
    'neutral_600': '#52525B',
    'neutral_700': '#3F3F46',
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

# ==================== PAGE-SPECIFIC ACCENT PALETTES ====================

HOME_COLORS = {
    'primary': '#2563EB',      # Refined scientific blue
    'secondary': '#3B82F6',
    'accent': '#60A5FA',
    'light': '#EFF6FF',
    'lighter': '#F8FAFF',
    'border': '#BFDBFE',
    'text': '#1E3A8A',
    'shadow': 'rgba(37, 99, 235, 0.12)',
}

INPUT_COLORS = {
    'primary': '#16A34A',
    'secondary': '#22C55E',
    'accent': '#4ADE80',
    'light': '#ECFDF5',
    'lighter': '#F6FFFA',
    'border': '#BBF7D0',
    'text': '#14532D',
    'shadow': 'rgba(22, 163, 74, 0.12)',
}

ANALYSIS_COLORS = {
    'primary': '#D97706',
    'secondary': '#F59E0B',
    'accent': '#FBBF24',
    'light': '#FFFBEB',
    'lighter': '#FFFEF7',
    'border': '#FDE68A',
    'text': '#7C2D12',
    'shadow': 'rgba(217, 119, 6, 0.12)',
}

RESULTS_COLORS = {
    'primary': '#7C3AED',
    'secondary': '#8B5CF6',
    'accent': '#A78BFA',
    'light': '#F5F3FF',
    'lighter': '#FAF8FF',
    'border': '#DDD6FE',
    'text': '#4C1D95',
    'shadow': 'rgba(124, 58, 237, 0.12)',
}

VISUALIZATION_COLORS = {
    'primary': '#0F766E',
    'secondary': '#14B8A6',
    'accent': '#5EEAD4',
    'light': '#F0FDFA',
    'lighter': '#F7FEFC',
    'border': '#99F6E4',
    'text': '#134E4A',
    'shadow': 'rgba(15, 118, 110, 0.12)',
}

DOWNLOAD_COLORS = {
    'primary': '#4338CA',
    'secondary': '#6366F1',
    'accent': '#818CF8',
    'light': '#EEF2FF',
    'lighter': '#F6F8FF',
    'border': '#C7D2FE',
    'text': '#1E1B4B',
    'shadow': 'rgba(67, 56, 202, 0.12)',
}

DOCUMENTATION_COLORS = {
    'primary': '#6D28D9',
    'secondary': '#8B5CF6',
    'accent': '#C4B5FD',
    'light': '#111827',
    'lighter': '#0B1120',
    'border': '#1F2937',
    'text': '#E5E7EB',
    'shadow': 'rgba(0, 0, 0, 0.45)',
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
    'Hybrid': {'primary': '#64748B','light': '#F8FAFC','lighter': '#F1F5F9','dark': '#1E293B'},
    'Non-B_DNA_Clusters': {'primary': '#334155','light': '#F1F5F9','lighter': '#E2E8F0','dark': '#0F172A'},
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

def get_motif_card_style(class_key: str) -> dict:
    colors = MOTIF_CARD_COLORS.get(class_key, MOTIF_CARD_COLORS['Hybrid'])
    return {
        'background': f"linear-gradient(135deg, {colors['light']} 0%, {colors['lighter']} 100%)",
        'border': colors['primary'],
        'text': colors['dark'],
    }

def get_motif_color(class_key: str) -> str:
    return UNIFIED_MOTIF_COLORS.get(class_key, '#808080')
