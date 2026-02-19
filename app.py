"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ NBDScanner - Non-B DNA Motif Finder                                          │
│ Main Application Entry Point                                                 │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.4            │
└──────────────────────────────────────────────────────────────────────────────┘
"""
# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
import streamlit as st
import pandas as pd
import os
import sys
import logging

_current_dir = os.path.dirname(os.path.abspath(__file__))
if _current_dir not in sys.path:
    sys.path.insert(0, _current_dir)

from Utilities.config.text import UI_TEXT
from Utilities.config.layout import LAYOUT_CONFIG
from Utilities.config.themes import TAB_THEMES
from UI.css import load_css
from UI import home, upload, results, download, documentation
from Utilities.nonbscanner import get_motif_info as get_motif_classification_info

# BioPython availability check
try:
    from Bio import Entrez, SeqIO
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False

# Hyperscan availability check (optional high-performance pattern matching)
logger = logging.getLogger(__name__)
HYPERSCAN_AVAILABLE = False
HYPERSCAN_VERSION = None
HYPERSCAN_ERROR = None
try:
    import hyperscan
    HYPERSCAN_AVAILABLE = True
    HYPERSCAN_VERSION = getattr(hyperscan, '__version__', 'unknown')
    logger.info(f"Hyperscan loaded (v{HYPERSCAN_VERSION})")
except ImportError as e:
    HYPERSCAN_ERROR = f"Hyperscan not installed: {e}"
    logger.info(f"Using pure Python fallback. {HYPERSCAN_ERROR}")
except Exception as e:
    HYPERSCAN_ERROR = f"Hyperscan init failed: {e}"
    logger.warning(f"Using pure Python fallback. {HYPERSCAN_ERROR}")

# ═══════════════════════════════════════════════════════════════════════════════
# TUNABLE PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════════
ENTREZ_EMAIL = "raazbiochem@gmail.com"
ENTREZ_API_KEY = None
DEFAULT_THEME_MODE = 'light'
DEFAULT_TABLE_DENSITY = 'relaxed'
DEFAULT_COLOR_THEME = 'scientific_blue'
PAGES = {"Home": "Overview", "Upload & Analyze": "Sequence Upload and Motif Analysis", "Results": "Analysis Results and Visualization", "Download": "Export Data", "Documentation": "Scientific Documentation & References"}
SESSION_DEFAULTS = {
    'seqs': [],  # Legacy: kept for backward compatibility
    'names': [],  # Legacy: kept for backward compatibility
    'results': [],  # Legacy: kept for backward compatibility
    'summary_df': pd.DataFrame(),
    'analysis_status': "Ready",
    'selected_classes': [],
    'selected_subclasses': [],
    'selected_classes_used': [],
    'selected_subclasses_used': [],
    'current_job_id': None,
    # Disk-based storage disabled: all processing uses RAM with 50KB/2KB chunking
    'use_disk_storage': False,
    'seq_storage': None,
    'seq_ids': [],
    'results_storage': {}
}
# ═══════════════════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════════════════
# APPLICATION INITIALIZATION
# ═══════════════════════════════════════════════════════════════════════════════
st.set_page_config(page_title=f"{UI_TEXT['app_title']} - Non-B DNA Motif Finder", layout=LAYOUT_CONFIG['layout_mode'], page_icon=None, menu_items={'About': f"NBDScanner | Developed by {UI_TEXT['author']}"})
CLASSIFICATION_INFO = get_motif_classification_info()

# Configure Entrez for NCBI API access (if BioPython is available)
if BIO_AVAILABLE:
    Entrez.email = ENTREZ_EMAIL
    Entrez.api_key = ENTREZ_API_KEY

# Initialize session state defaults
_THEME_DEFAULTS = {'theme_mode': DEFAULT_THEME_MODE, 'table_density': DEFAULT_TABLE_DENSITY, 'color_theme': DEFAULT_COLOR_THEME}
for k in _THEME_DEFAULTS:
    if k not in st.session_state:
        st.session_state[k] = _THEME_DEFAULTS[k]
for k, v in SESSION_DEFAULTS.items():
    if k not in st.session_state:
        st.session_state[k] = v

# ═══════════════════════════════════════════════════════════════════════════════
# MAIN APPLICATION
# ═══════════════════════════════════════════════════════════════════════════════
tabs = st.tabs(list(PAGES.keys()))
tab_pages = dict(zip(PAGES.keys(), tabs))
with tab_pages["Home"]: home.render()
with tab_pages["Upload & Analyze"]: upload.render()
with tab_pages["Results"]: results.render()
with tab_pages["Download"]: download.render()
with tab_pages["Documentation"]: documentation.render()
