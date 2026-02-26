"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Download Page - Export Data and Results                                      │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘
"""
# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
import streamlit as st
import pandas as pd
import numpy as np
import re
import io
import time
import logging
from collections import Counter
from Utilities.config.text import UI_TEXT
from Utilities.config.themes import TAB_THEMES
from UI.css import load_css
from UI.headers import render_section_heading
from UI.guards import generate_excel_bytes, generate_multifasta_excel_bytes
from UI.storage_helpers import has_results, get_sequences_info, get_results
from Utilities.utilities import export_to_csv, export_to_json, export_to_excel, export_to_pdf, export_to_bed
from Utilities.export.export_validator import validate_export_data

# ═══════════════════════════════════════════════════════════════════════════════
# TUNABLE PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════════
FILENAME_MAX_LENGTH = 50
_logger = logging.getLogger(__name__)
# ═══════════════════════════════════════════════════════════════════════════════

@st.cache_data(show_spinner=False)
def generate_statistics(all_motifs, names, lengths, seq_count):
    """
    Compute class- and subclass-level distribution statistics plus the combined
    stats Excel workbook.  Decorated with @st.cache_data so that download-button
    reruns never trigger a recomputation – the result is locked until the caller
    passes different arguments (new upload or reset both call st.cache_data.clear).
    """
    dist_df = pd.DataFrame()
    sub_df = pd.DataFrame()
    stats_excel = b""

    if all_motifs:
        # Group motifs by sequence name once (O(n+m)) to avoid repeated O(n*m) filters
        name_to_length = dict(zip(names, lengths))
        motifs_by_name: dict = {name: [] for name in names}
        for m in all_motifs:
            seq_name = m.get('Sequence_Name')
            if seq_name in motifs_by_name:
                motifs_by_name[seq_name].append(m)

        try:
            dist_data = []
            for name in names:
                slen = name_to_length.get(name, 0)
                motifs = motifs_by_name[name]
                cls_cnt = Counter(m.get('Class', 'Unknown') for m in motifs)
                for cname, cnt in cls_cnt.items():
                    cls_pos = set()
                    for m in motifs:
                        if m.get('Class') == cname:
                            s = m.get('Start', 0) - 1
                            e = m.get('End', 0)
                            if e > s:
                                cls_pos.update(range(s, e))
                    cls_covered = len(cls_pos)
                    gd = (cls_covered / slen * 100) if slen > 0 else 0
                    mkb = (cnt / slen * 1000) if slen > 0 else 0
                    avl = np.mean([m.get('Length', 0) for m in motifs if m.get('Class') == cname])
                    dist_data.append({'Sequence Name': name, 'Motif Class': cname.replace('_', ' '), 'Count': cnt, 'Coverage (%)': f"{gd:.4f}", 'Motifs per kbp': f"{mkb:.2f}", 'Average Length (bp)': f"{avl:.1f}", 'Covered Bases (bp)': cls_covered})
            dist_df = pd.DataFrame(dist_data)

            sub_data = []
            for name in names:
                slen = name_to_length.get(name, 0)
                motifs = motifs_by_name[name]
                sub_cnt = Counter(m.get('Subclass', 'Unknown') for m in motifs)
                for sname, cnt in sub_cnt.items():
                    pcls = next((m.get('Class') for m in motifs if m.get('Subclass') == sname), 'Unknown')
                    sub_pos = set()
                    for m in motifs:
                        if m.get('Subclass') == sname:
                            s = m.get('Start', 0) - 1
                            e = m.get('End', 0)
                            if e > s:
                                sub_pos.update(range(s, e))
                    sub_covered = len(sub_pos)
                    gd = (sub_covered / slen * 100) if slen > 0 else 0
                    mkb = (cnt / slen * 1000) if slen > 0 else 0
                    avl = np.mean([m.get('Length', 0) for m in motifs if m.get('Subclass') == sname])
                    sub_data.append({'Sequence Name': name, 'Motif Class': pcls.replace('_', ' '), 'Motif Subclass': sname.replace('_', ' '), 'Count': cnt, 'Coverage (%)': f"{gd:.4f}", 'Motifs per kbp': f"{mkb:.2f}", 'Average Length (bp)': f"{avl:.1f}", 'Covered Bases (bp)': sub_covered})
            sub_df = pd.DataFrame(sub_data)
        except Exception as exc:
            _logger.error("Statistics computation failed: %s", exc)

        if not dist_df.empty and not sub_df.empty:
            try:
                out = io.BytesIO()
                with pd.ExcelWriter(out, engine='openpyxl') as w:
                    dist_df.to_excel(w, sheet_name='Class Statistics', index=False)
                    sub_df.to_excel(w, sheet_name='Subclass Statistics', index=False)
                out.seek(0)
                stats_excel = out.getvalue()
            except Exception as exc:
                _logger.error("Stats Excel generation failed: %s", exc)

    return dist_df, sub_df, stats_excel


@st.cache_data(show_spinner=False)
def generate_all_exports(all_motifs, names, lengths, seq_count):
    export_times = {}
    csv_data = None
    excel_data = None
    excel_label = "📊 Excel"
    excel_fname = None
    json_data = None
    bed_data = None
    pdf_data = None
    excel_error = None
    pdf_error = None

    if not all_motifs:
        return {}

    try:
        _t0 = time.time()
        csv_data = export_to_csv(all_motifs, non_overlapping_only=False).encode('utf-8')
        export_times['csv'] = time.time() - _t0
    except Exception:
        pass

    try:
        _t0 = time.time()
        if seq_count > 1:
            excel_data = generate_multifasta_excel_bytes(names, lengths, seq_count)
            excel_label = "📊 Excel (MultiFASTA)"
            excel_fname = "multifasta_results.xlsx"
        else:
            excel_data = generate_excel_bytes(all_motifs, simple_format=True)
            excel_fname = "results.xlsx"
        export_times['excel'] = time.time() - _t0
    except Exception as e:
        excel_error = str(e)

    try:
        _t0 = time.time()
        json_data = export_to_json(all_motifs, pretty=True).encode('utf-8')
        export_times['json'] = time.time() - _t0
    except Exception:
        pass

    try:
        if names:
            _t0 = time.time()
            bed_data = export_to_bed(all_motifs, names[0]).encode('utf-8')
            export_times['bed'] = time.time() - _t0
    except Exception:
        pass

    try:
        if lengths and lengths[0] > 0:
            _t0 = time.time()
            pdf_data = export_to_pdf(all_motifs, lengths[0], names[0])
            export_times['pdf'] = time.time() - _t0
    except Exception as e:
        pdf_error = str(e)

    return {
        'csv': csv_data,
        'excel': excel_data,
        'excel_label': excel_label,
        'excel_fname': excel_fname,
        'json': json_data,
        'bed': bed_data,
        'pdf': pdf_data,
        'excel_error': excel_error,
        'pdf_error': pdf_error,
        'export_times': export_times,
    }

def render():
    load_css(TAB_THEMES.get('Download', 'clinical_teal')); render_section_heading("Download & Export Results", page="Downloads")
    if not has_results(): st.info(UI_TEXT['download_no_results']); st.markdown("<div style='background:linear-gradient(135deg,#f0f9ff 0%,#e0f2fe 100%);padding:1.5rem;border-radius:12px;margin-top:0.8rem;border:1px solid #bae6fd;text-align:center;'><h3 style='color:#0284c7;margin:0 0 0.7rem 0;font-size:1.4rem;'>Export Formats Available</h3><p style='color:#6b7280;margin:0 0 0.7rem 0;font-size:1.0rem;'>Once you analyze a sequence, export results in:</p><div style='display:flex;justify-content:center;gap:0.6rem;flex-wrap:wrap;'><span style='background:#0ea5e9;color:white;padding:0.45rem 0.9rem;border-radius:8px;font-weight:600;font-size:1.0rem;'>CSV</span><span style='background:#0ea5e9;color:white;padding:0.45rem 0.9rem;border-radius:8px;font-weight:600;font-size:1.0rem;'>Excel</span><span style='background:#0ea5e9;color:white;padding:0.45rem 0.9rem;border-radius:8px;font-weight:600;font-size:1.0rem;'>JSON</span><span style='background:#0ea5e9;color:white;padding:0.45rem 0.9rem;border-radius:8px;font-weight:600;font-size:1.0rem;'>BED</span><span style='background:#0ea5e9;color:white;padding:0.45rem 0.9rem;border-radius:8px;font-weight:600;font-size:1.0rem;'>PDF</span></div></div>", unsafe_allow_html=True); return
    
    # Get sequence information
    names, lengths, seq_count = get_sequences_info()
    seq_name = names[0] if names else "Unknown"; safe_fn = re.sub(r'[^\w\-]', '_', seq_name)[:FILENAME_MAX_LENGTH].strip('_')
    cls_used = st.session_state.get('selected_classes_used', []); sub_used = st.session_state.get('selected_subclasses_used', []); mode_used = st.session_state.get('analysis_mode_used', 'Motif Level')
    if cls_used: st.markdown(f"<div style='background:linear-gradient(135deg,#fefce8 0%,#fef9c3 100%);padding:0.6rem;border-radius:6px;margin-bottom:0.5rem;border-left:3px solid #ca8a04;'><p style='color:#713f12;margin:0;font-size:0.95rem;'><strong>Data Filter:</strong> Exports contain selected classes ({len(cls_used)}) and subclasses ({len(sub_used)}) from {mode_used} analysis</p></div>", unsafe_allow_html=True)
    
    # Collect all motifs from all sequences
    all_motifs = []
    for i in range(seq_count):
        motifs = get_results(i)  # Get all results for export
        for m in motifs: 
            em = m.copy()
            em.setdefault('Sequence_Name', names[i])
            all_motifs.append(em)
    
    try: all_motifs = validate_export_data(all_motifs, auto_normalize=True, strict=False)
    except Exception as e: st.warning(f"Some motifs normalized: {e}")
    
    # ═══════════════════════════════════════════════════════════════════════════════
    # STANDARD EXPORT FORMATS
    # Use @st.cache_data-backed generator so that reruns triggered by download
    # button clicks never recompute exports and never cause buttons to disappear.
    # ═══════════════════════════════════════════════════════════════════════════════
    exports = generate_all_exports(tuple(all_motifs), tuple(names), tuple(lengths), seq_count)
    csv_data = exports.get('csv')
    excel_data = exports.get('excel')
    excel_label = exports.get('excel_label', "📊 Excel")
    excel_fname = exports.get('excel_fname', f"{safe_fn}_results.xlsx")
    json_data = exports.get('json')
    bed_data = exports.get('bed')
    pdf_data = exports.get('pdf')
    excel_error = exports.get('excel_error')
    pdf_error = exports.get('pdf_error')
    export_times = exports.get('export_times', {})

    st.markdown("### Export Options"); st.markdown("<div style='background:linear-gradient(135deg,#f0f9ff 0%,#e0f2fe 100%);padding:0.8rem;border-radius:10px;margin-bottom:0.8rem;border-left:4px solid #0ea5e9;'><p style='color:#0c4a6e;margin:0;font-size:0.95rem;'><strong>Quick Export:</strong> Choose your preferred format for data and visualizations.</p></div>", unsafe_allow_html=True)
    st.markdown("#### Data Formats"); c1, c2, c3, c4, c5 = st.columns(5)
    with c1:
        st.download_button("📄 CSV", data=csv_data or b"", file_name=f"{safe_fn}_all_motifs.csv", mime="text/csv", use_container_width=True, type="primary", help="CSV with all motifs", disabled=(csv_data is None), key="dl_csv")
    with c2:
        st.download_button(excel_label or "📊 Excel", data=excel_data or b"", file_name=excel_fname or f"{safe_fn}_results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", use_container_width=True, type="primary", help="Excel workbook", disabled=(excel_data is None), key="dl_excel")
        if excel_error: st.error(f"Excel error: {excel_error}")
    with c3:
        st.download_button("🔗 JSON", data=json_data or b"", file_name=f"{safe_fn}_results.json", mime="application/json", use_container_width=True, type="primary", disabled=(json_data is None), key="dl_json")
    with c4:
        st.download_button("🧬 BED", data=bed_data or b"", file_name=f"{safe_fn}_results.bed", mime="text/plain", use_container_width=True, type="primary", disabled=(bed_data is None), key="dl_bed")
    with c5:
        st.download_button("📑 PDF", data=pdf_data or b"", file_name=f"{safe_fn}_viz.pdf", mime="application/pdf", use_container_width=True, type="primary", disabled=(pdf_data is None), key="dl_pdf")
        if pdf_error: st.error(f"PDF error: {pdf_error}")
        elif all_motifs and not (lengths and lengths[0] > 0): st.warning("No sequence for PDF")
    
    # ═══════════════════════════════════════════════════════════════════════════════
    # EXPORT TIMING SUMMARY
    # ═══════════════════════════════════════════════════════════════════════════════
    if export_times:
        with st.expander("⏱️ Export Preparation Times", expanded=False):
            _fmt_labels = {'csv': '📄 CSV', 'excel': '📊 Excel', 'json': '🔗 JSON', 'bed': '🧬 BED', 'pdf': '📑 PDF'}
            _cols = st.columns(len(export_times))
            for _col, (_fmt, _elapsed) in zip(_cols, export_times.items()):
                _ms = _elapsed * 1000
                _label = _fmt_labels.get(_fmt, _fmt.upper())
                _col.metric(_label, f"{_ms:.0f} ms")
    
    # ═══════════════════════════════════════════════════════════════════════════════
    # STATISTICAL ANALYSIS TABLES
    # Use @st.cache_data-backed generator so that reruns triggered by download
    # button clicks never recompute statistics and never cause buttons to flicker.
    # ═══════════════════════════════════════════════════════════════════════════════
    st.markdown("---"); st.markdown("### Statistical Analysis Tables"); st.markdown("<div style='background:linear-gradient(135deg,#f0fdf4 0%,#dcfce7 100%);padding:0.8rem;border-radius:10px;margin-bottom:0.8rem;border-left:4px solid #22c55e;'><p style='color:#14532d;margin:0;font-size:0.95rem;'><strong>Advanced Analytics:</strong> Detailed distribution and density statistics</p></div>", unsafe_allow_html=True)
    dist_df, sub_df, _stats_excel = generate_statistics(tuple(all_motifs), tuple(names), tuple(lengths), seq_count)
    st.markdown("#### Class-Level Distribution Statistics")
    if not dist_df.empty:
        st.dataframe(dist_df.head(10), use_container_width=True, height=300); st.caption(f"Showing first 10 of {len(dist_df)} records")
    else:
        st.info("No class-level statistics available.")
    st.markdown("#### Subclass-Level Distribution Statistics")
    if not sub_df.empty:
        st.dataframe(sub_df.head(10), use_container_width=True, height=300); st.caption(f"Showing first 10 of {len(sub_df)} records")
    else:
        st.info("No subclass-level statistics available.")
    _dist_csv = dist_df.to_csv(index=False).encode('utf-8') if not dist_df.empty else b""
    _sub_csv = sub_df.to_csv(index=False).encode('utf-8') if not sub_df.empty else b""
    s1, s2, s3 = st.columns(3)
    with s1:
        st.download_button("📈 Class Statistics (CSV)", data=_dist_csv, file_name=f"{safe_fn}_class_statistics.csv", mime="text/csv", use_container_width=True, disabled=dist_df.empty, key="dl_class_stats_csv")
    with s2:
        st.download_button("📉 Subclass Statistics (CSV)", data=_sub_csv, file_name=f"{safe_fn}_subclass_statistics.csv", mime="text/csv", use_container_width=True, disabled=sub_df.empty, key="dl_subclass_stats_csv")
    with s3:
        st.download_button("📊 All Statistics (Excel)", data=_stats_excel, file_name=f"{safe_fn}_all_statistics.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", use_container_width=True, disabled=(not _stats_excel), key="dl_all_stats_xlsx")
    
    # ═══════════════════════════════════════════════════════════════════════════════
    # VISUALIZATION NOTES
    # ═══════════════════════════════════════════════════════════════════════════════
    st.markdown("---")
    st.markdown("<div style='background:linear-gradient(135deg,#faf5ff 0%,#f3e8ff 100%);padding:1.0rem;border-radius:10px;border-left:4px solid #a855f7;'><p style='color:#6b21a8;margin:0;font-size:0.95rem;'><strong>📊 Visualization Note:</strong> All visualizations in the Results tab follow Nature publication standards with uniform axis labels, clean aesthetics, and colorblind-friendly palettes. The PDF export includes all primary visualizations for publication use.</p></div>", unsafe_allow_html=True)
