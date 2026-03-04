"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Home Page - NonBDNAFinder Landing Page                                       │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘
"""
# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
import streamlit as st
import os
from Utilities.config.text import UI_TEXT
from Utilities.config.themes import TAB_THEMES
from Utilities.config.colors import MOTIF_CLASS_INFO, get_motif_card_style
from UI.css import load_css, get_page_colors
from UI.headers import render_section_heading

# ═══════════════════════════════════════════════════════════════════════════════
# TUNABLE PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════════
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IMAGE_PATHS = [os.path.join(_PROJECT_ROOT, "nbdcircle.JPG"), "nbdcircle.JPG", "archive/nbdcircle.JPG", "./nbdcircle.JPG"]
# ═══════════════════════════════════════════════════════════════════════════════

def _build_motif_class_card(info):
    s = get_motif_card_style(info['key'])
    return f"<div style='padding:0.3rem 0.6rem;background:{s['background']};border-radius:10px;border:1px solid #1e293b;border-left:4px solid {s['border']};box-shadow:0 1px 4px rgba(0,0,0,0.07);transition:all 0.2s ease;'><div style='font-weight:700;color:{s['text']};font-size:0.87rem;line-height:1.3;'>{info['num']}. {info['name']}</div><div style='color:{s['text']};font-size:0.79rem;margin-top:0.05rem;line-height:1.2;opacity:0.85;'>{info['subtitle']}</div></div>"

def render():
    load_css(TAB_THEMES.get('Home', 'scientific_blue')); c = get_page_colors('Home')
    render_section_heading("Non-B DNA Finder", page="Home")
    col_nbd, col_science = st.columns([1, 1], gap="large")
    with col_nbd:
        img_found = False
        for p in IMAGE_PATHS:
            if os.path.exists(p): st.image(p, caption=UI_TEXT['home_image_caption'], use_container_width=True); img_found = True; break
        if not img_found: st.markdown(f"<div style='background:linear-gradient(135deg,{c['primary']} 0%,{c['secondary']} 100%);border-radius:15px;padding:40px;text-align:center;color:{c['white']};margin-bottom:1rem;'><h2 style='margin:0;color:{c['white']};font-size:2.0rem;'>{UI_TEXT['home_image_fallback_title']}</h2><h3 style='margin:8px 0 0 0;color:{c['white']};font-size:1.3rem;'>{UI_TEXT['home_image_fallback_subtitle']}</h3><p style='margin:5px 0 0 0;color:rgba(255,255,255,0.9);font-size:1.05rem;'>{UI_TEXT['home_image_fallback_caption']}</p></div>", unsafe_allow_html=True)
    with col_science:
        st.markdown(f"<div style='background:{c['white']};padding:0.5rem 1rem 0.7rem 1rem;border-radius:14px;box-shadow:0 2px 8px rgba(0,0,0,0.08);border:1.5px solid #1e293b;height:100%;'><h2 style='color:{c['text']};font-size:1.35rem;margin:0 0 0.45rem 0;font-weight:700;border-left:4px solid {c['accent']};padding-left:0.6rem;'>Scientific Foundation</h2><p style='color:{c['neutral_700']};font-size:1.0rem;line-height:1.6;margin-bottom:0.65rem;'><b style='color:{c['primary']};'>Non-canonical DNA structures</b> are critical regulatory elements implicated in genome stability, transcriptional regulation, replication, and disease mechanisms.</p><ul style='color:{c['neutral_700']};font-size:0.95rem;line-height:1.55;padding-left:1.5rem;margin-bottom:0.7rem;'><li><b>Genome Instability:</b> Hotspots for mutations and chromosomal rearrangements</li><li><b>Gene Regulation:</b> Promoter and enhancer activity modulation</li><li><b>DNA Replication:</b> Origins of replication and fork progression</li><li><b>Disease Association:</b> Cancer, neurological disorders, and aging</li></ul><div style='background:linear-gradient(135deg,{c['primary']} 0%,{c['secondary']} 60%,#FBBF24 100%);padding:0.4rem 0.9rem;border-radius:10px;text-align:center;box-shadow:0 2px 8px rgba(0,0,0,0.12);border:1.5px solid #1e293b;'><h3 style='color:{c['white']};margin:0 0 0.15rem 0;font-size:1.05rem;text-shadow:0 1px 2px rgba(0,0,0,0.3);'>{UI_TEXT['home_call_to_action_title']}</h3><p style='color:rgba(255,255,255,0.95);margin:0 0 0.35rem 0;font-size:0.9rem;'>{UI_TEXT['home_call_to_action_text']}</p><div style='background:{c['white']};color:{c['primary']};padding:0.35rem 1rem;border-radius:8px;display:inline-block;font-weight:700;font-size:0.9rem;box-shadow:0 1px 4px rgba(0,0,0,0.12);'>{UI_TEXT['home_call_to_action_button']}</div></div></div>", unsafe_allow_html=True)
    cards_html = ''.join(_build_motif_class_card(info) for info in MOTIF_CLASS_INFO)
    st.markdown(f"<div style='background:{c['white']};padding:0.5rem 1rem 0.7rem 1rem;border-radius:14px;box-shadow:0 2px 8px rgba(0,0,0,0.08);border:1.5px solid #1e293b;margin-top:0.75rem;'><h2 style='color:{c['text']};font-size:1.35rem;margin:0 0 0.45rem 0;font-weight:700;border-left:4px solid {c['accent']};padding-left:0.6rem;'>Detected Motif Classes</h2><div style='display:grid;grid-template-columns:repeat(auto-fit,minmax(185px,1fr));gap:0.45rem;margin-top:0.4rem;'>{cards_html}</div></div>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown(f"<div style='background:{c['white']};padding:0.5rem 1rem 0.7rem 1rem;border-radius:14px;box-shadow:0 2px 8px rgba(0,0,0,0.08);border:1.5px solid #1e293b;margin-top:0.75rem;'><h2 style='color:{c['primary']};font-size:1.35rem;margin:0 0 0.45rem 0;font-weight:700;border-left:4px solid {c['primary']};padding-left:0.6rem;'>How to Cite</h2><div style='background:{c['neutral_50']};padding:0.4rem 0.9rem;border-radius:10px;border-left:4px solid {c['primary']};font-family:\"Courier New\",monospace;font-size:0.9rem;line-height:1.55;color:{c['neutral_700']};box-shadow:0 1px 4px rgba(0,0,0,0.06);'><b>NonBFinder: Comprehensive Detection and Analysis of Non-B DNA Motifs</b><br>Dr. Venkata Rajesh Yella<br>GitHub: <a href=\"https://github.com/VRYella/NonBFinder\" style=\"color:{c['primary']};font-weight:700;\">https://github.com/VRYella/NonBFinder</a><br>Email: yvrajesh_bt@kluniversity.in</div><p style='color:{c['neutral_600']};font-size:0.9rem;margin-top:0.45rem;line-height:1.5;'>If you use NonBFinder in your research, please cite this resource. For methodology references, see the <b>Documentation</b> tab.</p></div>", unsafe_allow_html=True)

