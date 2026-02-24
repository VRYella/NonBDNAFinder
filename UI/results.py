"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Results Page - Analysis Results and Visualization                            │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘
"""
# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import logging
from collections import Counter
from Utilities.config.text import UI_TEXT
from Utilities.config.themes import TAB_THEMES
from Utilities.config.analysis import MAX_OVERLAP_DISPLAY
from UI.css import load_css
from UI.headers import render_section_heading
from UI.storage_helpers import (
    has_results, get_sequences_info, get_sequence_length,
    get_results, get_results_summary, get_results_dataframe
)
from Utilities.utilities import (
    get_basic_stats, export_results_to_dataframe, optimize_dataframe_memory, 
    calculate_genomic_density, calculate_positional_density, calculate_motif_statistics,
    plot_motif_distribution, plot_nested_pie_chart, plot_linear_motif_track, 
    plot_linear_subclass_track, plot_density_comparison, plot_cluster_size_distribution, 
    plot_motif_cooccurrence_matrix, plot_motif_length_kde, plot_score_distribution,
    plot_score_violin, plot_structural_heatmap, plot_motif_network,
    plot_chromosome_density, plot_spacer_loop_variation, plot_motif_clustering_distance,
    plot_structural_competition_upset, compute_comprehensive_genome_stats
)
from Utilities.visualization import NATURE_MOTIF_COLORS
from Utilities.multifasta_engine import MultiFastaEngine

logger = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════════
# TUNABLE PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════════
CLUSTER_CLASSES = ['Hybrid', 'Non-B_DNA_Clusters']
# ═══════════════════════════════════════════════════════════════════════════════

def _render_section_divider(label): st.markdown(f"<div style='display:flex;align-items:center;gap:8px;padding:2px 0;margin-top:6px;'><span style='font-size:0.8rem;color:#64748b;font-weight:600;'>{label}</span><div style='flex:1;height:1px;background:linear-gradient(90deg,#a855f7 0%,transparent 100%);'></div></div>", unsafe_allow_html=True)

def _render_metric_panel(header: str, rows: list, note: str = "") -> None:
    """Render a structured metric panel as an HTML table with Metric | Value | Definition columns.

    Args:
        header: Section header string (e.g. "🟦 I. Genome Overview").
        rows: List of (metric_name, value_str, definition_str) tuples.
        note: Optional italicised footnote shown below the table.
    """
    th_style = "background:#f1f5f9;color:#334155;font-size:0.75rem;font-weight:700;padding:5px 10px;text-align:left;border-bottom:2px solid #e2e8f0;"
    td_style = "font-size:0.78rem;color:#1e293b;padding:5px 10px;border-bottom:1px solid #f1f5f9;vertical-align:top;"
    val_style = "font-size:0.78rem;font-weight:600;color:#7c3aed;padding:5px 10px;border-bottom:1px solid #f1f5f9;white-space:nowrap;"
    def_style = "font-size:0.75rem;color:#64748b;padding:5px 10px;border-bottom:1px solid #f1f5f9;font-style:italic;vertical-align:top;"
    rows_html = "".join(
        f"<tr><td style='{td_style}'>{m}</td><td style='{val_style}'>{v}</td><td style='{def_style}'>{d}</td></tr>"
        for m, v, d in rows
    )
    note_html = f"<p style='font-size:0.73rem;color:#64748b;margin:4px 0 0 2px;font-style:italic;'>{note}</p>" if note else ""
    st.markdown(
        f"<div style='margin-top:10px;'>"
        f"<div style='font-size:0.82rem;font-weight:700;color:#475569;margin-bottom:4px;'>{header}</div>"
        f"<table style='width:100%;border-collapse:collapse;border:1px solid #e2e8f0;border-radius:6px;overflow:hidden;'>"
        f"<thead><tr><th style='{th_style}'>Metric</th><th style='{th_style}'>Value</th><th style='{th_style}'>Definition</th></tr></thead>"
        f"<tbody>{rows_html}</tbody></table>"
        f"{note_html}</div>",
        unsafe_allow_html=True,
    )

def _render_analysis_summary_box(cov, den, cnt, slen): st.markdown(f"<div style='display:flex;flex-wrap:wrap;gap:4px;padding:5px 10px;background:linear-gradient(135deg,#faf5ff 0%,#f3e8ff 100%);border-radius:6px;border:1px solid #e9d5ff;margin-bottom:8px;justify-content:space-around;align-items:center;'><div style='display:flex;flex-direction:column;align-items:center;padding:1px 8px;'><span style='font-size:0.95rem;font-weight:800;background:linear-gradient(135deg,#a855f7,#8b5cf6);-webkit-background-clip:text;-webkit-text-fill-color:transparent;'>{cov:.2f}%</span><span style='font-size:0.6rem;color:#64748b;text-transform:uppercase;'>Coverage</span></div><div style='display:flex;flex-direction:column;align-items:center;padding:1px 8px;'><span style='font-size:0.95rem;font-weight:800;background:linear-gradient(135deg,#a855f7,#8b5cf6);-webkit-background-clip:text;-webkit-text-fill-color:transparent;'>{den:.2f}</span><span style='font-size:0.6rem;color:#64748b;text-transform:uppercase;'>Motifs/kb</span></div><div style='display:flex;flex-direction:column;align-items:center;padding:1px 8px;'><span style='font-size:0.95rem;font-weight:800;background:linear-gradient(135deg,#a855f7,#8b5cf6);-webkit-background-clip:text;-webkit-text-fill-color:transparent;'>{cnt:,}</span><span style='font-size:0.6rem;color:#64748b;text-transform:uppercase;'>Motifs</span></div><div style='display:flex;flex-direction:column;align-items:center;padding:1px 8px;'><span style='font-size:0.95rem;font-weight:800;background:linear-gradient(135deg,#a855f7,#8b5cf6);-webkit-background-clip:text;-webkit-text-fill-color:transparent;'>{slen:,}</span><span style='font-size:0.6rem;color:#64748b;text-transform:uppercase;'>bp</span></div></div>", unsafe_allow_html=True)

def _calculate_overlaps(motifs, by='Class'):
    overlaps = {}; sorted_m = sorted(motifs, key=lambda m: m.get('Start', 0))
    for i, m1 in enumerate(sorted_m):
        for m2 in sorted_m[i+1:]:
            if m2.get('Start', 0) >= m1.get('End', 0): break
            k1, k2 = m1.get(by, 'Unknown'), m2.get(by, 'Unknown')
            if k1 != k2: pair = tuple(sorted([k1, k2])); overlaps[pair] = overlaps.get(pair, 0) + 1
    return overlaps

def _render_overlap_matrix(overlaps, title):
    if not overlaps: return
    sorted_ovl = sorted(overlaps.items(), key=lambda x: x[1], reverse=True)[:MAX_OVERLAP_DISPLAY]
    html = f"<div style='font-size:0.8rem;font-weight:600;margin-bottom:4px;color:#334155;'>{title}</div><div style='font-size:0.75rem;color:#64748b;'>"
    for (k1, k2), c in sorted_ovl: html += f'<div style="padding:2px 0;">{k1.replace("_"," ")} ↔ {k2.replace("_"," ")}: <strong>{c}</strong></div>'
    html += "</div>"; st.markdown(html, unsafe_allow_html=True)

def render():
    load_css(TAB_THEMES.get('Results', 'genomic_purple')); render_section_heading("Analysis Results & Visualization", page="Results")
    if not has_results(): st.info(UI_TEXT['status_no_results']); st.info("Run analysis first in 'Upload & Analyze' tab."); st.markdown("<div style='background:linear-gradient(135deg,#faf5ff 0%,#f3e8ff 100%);padding:1.2rem;border-radius:12px;margin-top:0.8rem;border:1px solid #e9d5ff;text-align:center;'><h3 style='color:#7c3aed;margin:0 0 0.6rem 0;'>Visualization Preview</h3><p style='color:#6b7280;margin:0;'>Upload and analyze a sequence to see motif track visualizations, distributions, density plots, and more.</p></div>", unsafe_allow_html=True); return
    with st.expander(f"Summary ({len(st.session_state.summary_df)} rows)", expanded=False): st.dataframe(st.session_state.summary_df, use_container_width=True)
    
    # Get sequence information
    names, lengths, seq_count = get_sequences_info()
    
    # ═══════════════════════════════════════════════════════════════════════════════
    # AGGREGATED MODE IF >20 SEQUENCES
    # ═══════════════════════════════════════════════════════════════════════════════
    
    if seq_count > 20:
        st.info("📊 Aggregated visualization mode enabled (>20 sequences detected)")
        
        # Collect annotations by sequence
        annotations_by_sequence = {}
        sequence_lengths = {}
        
        for i in range(seq_count):
            fasta_id = names[i]
            motifs_i = get_results(i)
            annotations_by_sequence[fasta_id] = motifs_i
            sequence_lengths[fasta_id] = lengths[i]
        
        engine = MultiFastaEngine(
            annotations_by_sequence,
            sequence_lengths
        )
        
        # ---------------------------------------------------------
        # 1. GLOBAL CLASS DISTRIBUTION
        # ---------------------------------------------------------
        st.subheader("Global Motif Class Distribution")
        class_dist = engine.global_class_distribution()
        if class_dist:
            st.bar_chart(class_dist)
        else:
            st.info("No motifs detected across sequences.")
        
        # ---------------------------------------------------------
        # 2. CLASS COVERAGE (% sequences)
        # ---------------------------------------------------------
        st.subheader("% Sequences Containing Each Motif Class")
        coverage = engine.class_sequence_coverage()
        if coverage:
            st.bar_chart(coverage)
        else:
            st.info("No class coverage data available.")
        
        # ---------------------------------------------------------
        # 3. MOTIF DENSITY HISTOGRAM
        # ---------------------------------------------------------
        st.subheader("Motif Density Distribution (motifs per kb)")
        densities = engine.motif_density_per_sequence()
        if densities:
            density_df = pd.DataFrame({"Density": densities})
            st.bar_chart(density_df["Density"])
        else:
            st.info("No density data available.")
        
        # ---------------------------------------------------------
        # 4. POSITIONAL CONSERVATION (if equal length)
        # ---------------------------------------------------------
        equal_length = len(set(lengths)) == 1
        
        if equal_length:
            st.subheader("Positional Conservation")
            seq_length = lengths[0]
            conservation = engine.positional_conservation(seq_length)
            if conservation:
                st.line_chart(conservation)
            else:
                st.info("No conservation data available.")
        else:
            st.info("⚠️ Positional conservation analysis requires equal-length sequences.")
        
        return
    
    # ═══════════════════════════════════════════════════════════════════════════════
    # MULTIFASTA UNIFIED VISUALIZATIONS (≤20 SEQUENCES)
    # ═══════════════════════════════════════════════════════════════════════════════
    if seq_count > 1:
        st.markdown("---")
        st.markdown("### 📊 MultiFASTA Unified Analysis")
        st.markdown("<div style='background:linear-gradient(135deg,#fef3c7 0%,#fde68a 100%);padding:0.6rem;border-radius:10px;margin-bottom:0.8rem;border-left:4px solid #f59e0b;'><p style='color:#92400e;margin:0;font-size:0.8rem;'><strong>Unified View:</strong> Comparing motif patterns across all sequences in your MultiFASTA input.</p></div>", unsafe_allow_html=True)
        
        # Check if all sequences have equal length
        equal_length = len(set(lengths)) == 1
        
        # Show equal-length indicator
        if equal_length:
            st.info(f"✨ All sequences have equal length ({lengths[0]:,} bp). Positional occurrence analysis is available below.")
        
        # Collect annotations by sequence
        from Utilities.multifasta_visualizer import MultiFastaVisualizer
        annotations_by_sequence = {}
        sequence_lengths = {}
        
        for i in range(seq_count):
            fasta_id = names[i]
            motifs_i = get_results(i)
            annotations_by_sequence[fasta_id] = motifs_i
            sequence_lengths[fasta_id] = lengths[i]
        
        visualizer = MultiFastaVisualizer(annotations_by_sequence)
        
        # Generate unified summary statistics
        summary_stats = visualizer.generate_unified_summary()
        
        # Display sequence-level statistics
        with st.expander("📋 Sequence Statistics", expanded=True):
            seq_stats = []
            for stat in summary_stats['sequence_stats']:
                fasta_id = stat['FASTA_ID']
                total = stat['Total_Motifs']
                seq_len = sequence_lengths[fasta_id]
                density = (total / seq_len * 1000) if seq_len > 0 else 0
                seq_motif_stats = calculate_motif_statistics(annotations_by_sequence[fasta_id], seq_len)
                coverage_pct = seq_motif_stats.get('Coverage%', 0)
                seq_stats.append({
                    'Sequence': fasta_id,
                    'Length (bp)': f"{seq_len:,}",
                    'Total Motifs': total,
                    'Motifs/kb': f"{density:.2f}",
                    'Coverage%': f"{coverage_pct:.2f}"
                })
            st.dataframe(pd.DataFrame(seq_stats), use_container_width=True)
        
        # Unified visualizations
        unified_tabs = st.tabs(["Class Distribution", "Density Heatmap", "Positional Analysis"])
        
        with unified_tabs[0]:
            st.markdown("#### Motif Class Distribution per Sequence")
            try:
                fig = visualizer.generate_class_distribution_plot(figsize=(12, 6))
                st.pyplot(fig)
                plt.close(fig)
            except Exception as e:
                st.error(f"Class distribution plot error: {e}")
        
        with unified_tabs[1]:
            st.markdown("#### Motif Density Heatmap")
            try:
                fig = visualizer.generate_density_heatmap(sequence_lengths, figsize=(12, 6))
                st.pyplot(fig)
                plt.close(fig)
            except Exception as e:
                st.error(f"Density heatmap error: {e}")
        
        with unified_tabs[2]:
            if equal_length:
                st.markdown("#### Positional Occurrence Analysis")
                st.markdown("<div style='background:linear-gradient(135deg,#dbeafe 0%,#bfdbfe 100%);padding:0.6rem;border-radius:10px;margin-bottom:0.8rem;'><p style='color:#1e3a8a;margin:0;font-size:0.75rem;'><strong>Note:</strong> Shows how many sequences have a motif at each position. Useful for identifying conserved regions.</p></div>", unsafe_allow_html=True)
                try:
                    seq_length = lengths[0]
                    figures = visualizer.generate_positional_panels(seq_length, smooth=True, figsize=(12, 6))
                    
                    if 'overall' in figures:
                        st.markdown("**Overall Positional Occurrence**")
                        st.pyplot(figures['overall'])
                        plt.close(figures['overall'])
                    
                    if 'by_class' in figures:
                        st.markdown("**Positional Occurrence by Class**")
                        st.pyplot(figures['by_class'])
                        plt.close(figures['by_class'])
                except Exception as e:
                    st.error(f"Positional analysis error: {e}")
            else:
                st.info("⚠️ Positional occurrence analysis is only available when all sequences have equal length. Your sequences have different lengths.")
        
        st.markdown("---")
        st.markdown("### 📄 Individual Sequence Analysis")
        st.markdown("<div style='background:linear-gradient(135deg,#f3e8ff 0%,#e9d5ff 100%);padding:0.6rem;border-radius:10px;margin-bottom:0.8rem;'><p style='color:#581c87;margin:0;font-size:0.8rem;'>Select a sequence below to view detailed per-sequence visualizations.</p></div>", unsafe_allow_html=True)
    
    seq_idx = 0
    if seq_count > 1:
        try: sel = st.pills("Sequence:", options=list(range(seq_count)), format_func=lambda i: names[i], selection_mode="single", default=0); seq_idx = sel or 0
        except: seq_idx = 0
    
    # Get results for selected sequence with adaptive loading strategy
    # First, get the total count to decide on strategy
    summary = get_results_summary(seq_idx)
    total_count = summary.get('total_count', 0)
    
    # Adaptive loading: full load for <100K motifs, limited load for >100K
    if total_count < 100_000:
        # Full load for accurate visualization
        motifs = get_results(seq_idx, limit=None)
    else:
        # Load first 100,000 motifs for performance (sequential, not random sampling)
        motifs = get_results(seq_idx, limit=100_000)
        # Inform user about limited loading
        st.info(f"📊 Displaying limited data: showing first 100,000 of {total_count:,} motifs for visualization performance. Statistics are computed from complete dataset.")
    
    slen = get_sequence_length(seq_idx)
    
    if not motifs: st.warning("No motifs detected."); return
    df = pd.DataFrame(motifs) if motifs else pd.DataFrame()
    if len(df) > 1000: df = optimize_dataframe_memory(df)
    
    # Get sequence for stats calculation
    # For disk storage mode, we need to load sequence for stats
    if st.session_state.get('use_disk_storage') and st.session_state.get('seq_ids'):
        seq_id = st.session_state.seq_ids[seq_idx]
        # For large sequences, calculate stats from metadata instead of loading full sequence
        if slen > 10_000_000:  # 10MB threshold
            # Use results summary for stats
            summary = get_results_summary(seq_idx)
            coverage_pct = (summary['coverage_bp'] / slen * 100) if slen > 0 else 0
            density = (summary['total_count'] / (slen / 1000)) if slen > 0 else 0
            _render_analysis_summary_box(coverage_pct, density, summary['total_count'], slen)
        else:
            # Load sequence for accurate stats
            seq = st.session_state.seq_storage.get_sequence_chunk(seq_id, 0, slen)
            stats = get_basic_stats(seq, motifs)
            _render_analysis_summary_box(stats.get("Coverage%", 0), stats.get("Density", 0), len(motifs), slen)
    else:
        # Legacy mode - sequence is in memory
        stats = get_basic_stats(st.session_state.seqs[seq_idx], motifs)
        _render_analysis_summary_box(stats.get("Coverage%", 0), stats.get("Density", 0), len(motifs), slen)
    
    has_clusters = any(m.get('Class') == 'Non-B_DNA_Clusters' for m in motifs); has_hybrids = any(m.get('Class') == 'Hybrid' for m in motifs)
    
    # ═══════════════════════════════════════════════════════════════════════════════
    # COMPREHENSIVE GENOME STATISTICS
    # ═══════════════════════════════════════════════════════════════════════════════
    try:
        gstats = compute_comprehensive_genome_stats(motifs, slen)
        with st.expander("🧬 GENOME STRUCTURAL LANDSCAPE REPORT", expanded=False):
            st.markdown(
                "<div style='color:#64748b;font-size:0.78rem;margin-bottom:8px;'>"
                "Overall genome coverage excludes Hybrid and Cluster regions (reported individually below)."
                "</div>",
                unsafe_allow_html=True,
            )
            # ── Section I: Genome Overview ────────────────────────────────────
            _render_metric_panel(
                "🟦 I. Genome Overview",
                [
                    ("Genome Length", f"{gstats['genome_length']:,} bp", "Total sequence length analyzed"),
                    ("Motifs (excl. Hybrid/Cluster)", f"{gstats['n_motifs']:,}", "Individual structural motifs detected"),
                    ("Motifs (incl. Hybrid/Cluster)", f"{gstats['n_motifs_all']:,}", "Total including merged structural regions"),
                    ("Motif Classes (C)", f"{gstats['n_classes']}", "Distinct structural motif types"),
                    ("Motif Density", f"{gstats['density_per_kb']:.4f} / kb", "Motifs per kilobase"),
                ],
            )
            # ── Section II: Structural Coverage ──────────────────────────────
            _render_metric_panel(
                "🟩 II. Structural Coverage",
                [
                    ("Total Covered Bases", f"{gstats['total_covered_bases']:,} bp", "Unique genomic bases overlapping motifs"),
                    ("Coverage Fraction", f"{gstats['coverage_fraction']:.6f}", "Covered bases / genome length"),
                    ("Coverage (%)", f"{gstats['coverage_pct']:.4f}%", "% of genome structurally annotated"),
                ],
                note="🔎 Coverage reflects structural footprint, excluding Hybrid/Cluster merged regions.",
            )
            # ── Section III: Occupancy Metrics ────────────────────────────────
            _render_metric_panel(
                "🟨 III. Occupancy Metrics",
                [
                    ("Raw Occupancy", f"{gstats['raw_occupancy_bp']:,} bp", "Total motif bases including overlaps"),
                    ("Normalized Occupancy (SLI)", f"{gstats['normalized_occupancy']:.6f}", "Structural Load Index (overlap-adjusted)"),
                    ("Mean Overlap Depth", f"{gstats['mean_overlap_depth']:.4f}", "Average motif stacking depth"),
                ],
                note="📌 SLI measures genome-wide structural burden.",
            )
            # ── Section IV: Class-Specific Coverage ──────────────────────────
            if gstats['class_covered_bases']:
                td_s = "font-size:0.76rem;color:#1e293b;padding:5px 8px;border-bottom:1px solid #f1f5f9;"
                th_s = "background:#f1f5f9;color:#334155;font-size:0.74rem;font-weight:700;padding:5px 8px;text-align:left;border-bottom:2px solid #e2e8f0;"
                cls_html = "".join(
                    f"<tr>"
                    f"<td style='{td_s}'>{cls}</td>"
                    f"<td style='{td_s}'>{gstats['class_covered_bases'].get(cls, 0):,}</td>"
                    f"<td style='{td_s}'>{gstats['class_coverage_pct'].get(cls, 0):.4f}%</td>"
                    f"<td style='{td_s}'>{gstats['class_contribution'].get(cls, 0):.2f}%</td>"
                    f"</tr>"
                    for cls in sorted(
                        gstats['class_covered_bases'],
                        key=lambda c: gstats['class_covered_bases'].get(c, 0),
                        reverse=True,
                    )
                )
                st.markdown(
                    f"<div style='margin-top:10px;'>"
                    f"<div style='font-size:0.82rem;font-weight:700;color:#475569;margin-bottom:4px;'>🟧 IV. Class-Specific Coverage</div>"
                    f"<table style='width:100%;border-collapse:collapse;border:1px solid #e2e8f0;border-radius:6px;overflow:hidden;'>"
                    f"<thead><tr>"
                    f"<th style='{th_s}'>Class</th>"
                    f"<th style='{th_s}'>Covered (bp)</th>"
                    f"<th style='{th_s}'>Coverage (%)</th>"
                    f"<th style='{th_s}'>Contribution</th>"
                    f"</tr></thead>"
                    f"<tbody>{cls_html}</tbody></table>"
                    f"<p style='font-size:0.73rem;color:#64748b;margin:4px 0 0 2px;font-style:italic;'>"
                    f"📊 Contribution indicates proportional share of total structural coverage.</p>"
                    f"</div>",
                    unsafe_allow_html=True,
                )
            # ── Section V: Structural Load Metrics ───────────────────────────
            _render_metric_panel(
                "🟥 V. Structural Load Metrics",
                [
                    ("SLI", f"{gstats['sli']:.6f}", "Normalized genome structural burden"),
                    ("Structural Intensity", f"{gstats['structural_intensity']:.6f}", "Weighted motif stacking strength"),
                    ("Weighted Structural Coverage", f"{gstats['weighted_structural_coverage']:.6f}", "Overlap-weighted coverage fraction"),
                ],
            )
            # ── Section VI: Spatial Distribution ─────────────────────────────
            _render_metric_panel(
                "🟪 VI. Spatial Distribution",
                [
                    ("Mean Inter-Motif Distance", f"{gstats['mean_inter_motif_distance']:.2f} bp", "Average spacing between motifs"),
                    ("CV (Clustering Coefficient)", f"{gstats['cv_spatial_clustering']:.4f}", "Spatial clustering variability"),
                ],
                note="CV > 1 indicates clustered distribution.",
            )
            # ── Section VII: Hotspot / Cluster Metrics ────────────────────────
            _render_metric_panel(
                "🟫 VII. Hotspot / Cluster Metrics",
                [
                    (f"Max Local Density (W={gstats['window_size']:,} bp)", f"{gstats['max_local_density']:.6f}", "Highest motif density window"),
                    ("Max Class Diversity", f"{gstats['max_class_diversity_window']}", "Distinct classes in single window"),
                    ("Max Cluster Score", f"{gstats['max_cluster_score']:.6f}", "Peak structural aggregation score"),
                ],
            )
            # ── Section VIII: Hybrid & Cluster (Reported Separately) ─────────
            hc_td = "font-size:0.76rem;color:#1e293b;padding:5px 8px;border-bottom:1px solid #f1f5f9;"
            hc_th = "background:#f1f5f9;color:#334155;font-size:0.74rem;font-weight:700;padding:5px 8px;text-align:left;border-bottom:2px solid #e2e8f0;"
            hc_html = (
                f"<tr><td style='{hc_td}'>Hybrid Regions</td><td style='{hc_td}'>{gstats['hybrid_count']:,}</td></tr>"
                f"<tr><td style='{hc_td}'>Hybrid Coverage</td><td style='{hc_td}'>{gstats['hybrid_coverage_pct']:.4f}%</td></tr>"
                f"<tr><td style='{hc_td}'>Cluster Regions</td><td style='{hc_td}'>{gstats['cluster_count']:,}</td></tr>"
                f"<tr><td style='{hc_td}'>Cluster Coverage</td><td style='{hc_td}'>{gstats['cluster_coverage_pct']:.4f}%</td></tr>"
                f"<tr><td style='{hc_td}'>Mean Overlap Fraction</td><td style='{hc_td}'>{gstats['mean_overlap_fraction']:.4f}</td></tr>"
            )
            st.markdown(
                f"<div style='margin-top:10px;'>"
                f"<div style='font-size:0.82rem;font-weight:700;color:#475569;margin-bottom:4px;'>🔶 VIII. Hybrid &amp; Cluster Regions (Reported Separately)</div>"
                f"<table style='width:100%;border-collapse:collapse;border:1px solid #e2e8f0;border-radius:6px;overflow:hidden;'>"
                f"<thead><tr><th style='{hc_th}'>Metric</th><th style='{hc_th}'>Value</th></tr></thead>"
                f"<tbody>{hc_html}</tbody></table>"
                f"<p style='font-size:0.73rem;color:#64748b;margin:4px 0 0 2px;font-style:italic;'>"
                f"⚠ Hybrid/Cluster regions are excluded from primary genome coverage metrics.</p>"
                f"</div>",
                unsafe_allow_html=True,
            )
            # ── Section IX: Structural Diversity ─────────────────────────────
            _render_metric_panel(
                "🟦 IX. Structural Diversity",
                [
                    ("Simpson Diversity Index (D)", f"{gstats['simpson_diversity_index']:.4f}", "Probability two motifs differ in class"),
                    ("Effective Class Number (Neff)", f"{gstats['effective_class_number']:.4f}", "Functional structural diversity"),
                ],
            )
            # ── Section X: Genome-Scale Comparative Indices ───────────────────
            _render_metric_panel(
                "🟩 X. Genome-Scale Comparative Indices",
                [
                    ("SCI (Structural Complexity Index)", f"{gstats['sci']:.4f}", "Composite structural diversity × density"),
                    ("Structural Dominance Ratio", f"{gstats['dominance_ratio']:.4f}", "Dominance of most abundant class"),
                ],
            )
    except Exception as _gse:
        logger.warning(f"Comprehensive genome stats error: {_gse}")

    # ═══════════════════════════════════════════════════════════════════════════════
    # TWO-TAB VISUALIZATION LAYOUT (Nature Publication Standard)
    # ═══════════════════════════════════════════════════════════════════════════════
    viz_tabs = st.tabs(["Primary Visualizations", "Advanced Analysis"])

    # Pre-generated figure bytes (cached during analysis phase for instant display)
    _cf = st.session_state.get('cached_visualizations', {}).get(f"seq_{seq_idx}", {}).get('figures', {})

    def _show_fig(cache_key, fig_fn, error_msg):
        """Display a pre-generated PNG or fall back to on-demand matplotlib render."""
        if cache_key in _cf:
            st.image(_cf[cache_key], use_container_width=True)
        else:
            try:
                fig = fig_fn()
                st.pyplot(fig)
                plt.close(fig)
            except Exception as e:
                st.error(f"{error_msg}: {e}")
    
    # ═══════════════════════════════════════════════════════════════════════════════
    # TAB 1: PRIMARY VISUALIZATIONS
    # ═══════════════════════════════════════════════════════════════════════════════
    with viz_tabs[0]:
        _render_section_divider("Genome Track")
        _show_fig('class_track', lambda: plot_linear_motif_track(motifs, slen, title="Class Track"), "Track error")
        
        _render_section_divider("Subclass Track"); sm = [m for m in motifs if m.get('Class') not in CLUSTER_CLASSES]
        if sm:
            _show_fig('subclass_track', lambda: plot_linear_subclass_track(sm, slen, title="Subclass Track"), "Subclass track error")
        else: st.info("No non-cluster motifs.")
        
        _render_section_divider("Distribution"); c1, c2 = st.columns(2)
        with c1:
            _show_fig('class_dist', lambda: plot_motif_distribution(motifs, by='Class', title="Class Distribution"), "Class dist error")
        with c2:
            _show_fig('subclass_dist', lambda: plot_motif_distribution(motifs, by='Subclass', title="Subclass Distribution"), "Subclass dist error")
        
        _render_section_divider("Density Analysis")
        try:
            gd = None
            pd_kbp = None
            vk = f"seq_{seq_idx}"; cv = st.session_state.get('cached_visualizations', {}).get(vk, {}); cd = cv.get('densities', {})
            if cd.keys() >= {'class_genomic', 'class_positional'}: 
                gd, pd_kbp = cd['class_genomic'], cd['class_positional']
            else: 
                gd = calculate_genomic_density(motifs, slen, by_class=True)
                pd_kbp = calculate_positional_density(motifs, slen, unit='kbp', by_class=True)
            if gd is not None and pd_kbp is not None:
                _show_fig('density_comparison', lambda: plot_density_comparison(gd, pd_kbp, title="Density Analysis"), "Density error")
        except Exception as e: st.error(f"Density error: {e}")
        
        _render_section_divider("Length Distribution")
        _show_fig('length_kde', lambda: plot_motif_length_kde(motifs, by_class=True, title="Length Distribution"), "Length dist error")
        
        _render_section_divider("Score Distribution")
        _show_fig('score_violin', lambda: plot_score_violin(motifs, by_class=True, title="Score Distribution"), "Score error")
        
        _render_section_divider("Composition")
        _show_fig('nested_pie', lambda: plot_nested_pie_chart(motifs, title="Class → Subclass"), "Pie error")
    
    # ═══════════════════════════════════════════════════════════════════════════════
    # TAB 2: ADVANCED ANALYSIS
    # ═══════════════════════════════════════════════════════════════════════════════
    with viz_tabs[1]:
        _render_section_divider("Structural Potential Heatmap")
        _show_fig('struct_heatmap', lambda: plot_structural_heatmap(motifs, slen, title="Structural Potential Heatmap"), "Heatmap error")
        
        _render_section_divider("Co-occurrence Network")
        _show_fig('motif_network', lambda: plot_motif_network(motifs, title="Motif Co-occurrence Network"), "Network error")
        
        _render_section_divider("Co-occurrence Matrix")
        _show_fig('cooccurrence_mat', lambda: plot_motif_cooccurrence_matrix(motifs, title="Co-occurrence Matrix"), "Co-occurrence error")
        
        _render_section_divider("Chromosome Density")
        _show_fig('chrom_density', lambda: plot_chromosome_density(motifs, title="Motif Density by Class"), "Chromosome density error")
        
        _render_section_divider("Clustering Distance")
        _show_fig('cluster_dist_fig', lambda: plot_motif_clustering_distance(motifs, title="Inter-Motif Distance"), "Clustering distance error")
        
        _render_section_divider("Structural Features (Loop/Arm Lengths)")
        _show_fig('spacer_loop', lambda: plot_spacer_loop_variation(motifs, title="Structural Features Distribution"), "Structural features error")
        
        if has_clusters or has_hybrids:
            _render_section_divider("Clusters & Hybrids"); chm = [m for m in motifs if m.get('Class') in CLUSTER_CLASSES]
            _show_fig('cluster_track', lambda: plot_linear_motif_track(chm, slen, title="Hybrid & Cluster Track"), "Cluster track error")
            if has_clusters:
                _render_section_divider("Cluster Statistics")
                _show_fig('cluster_size', lambda: plot_cluster_size_distribution(motifs, title="Cluster Statistics"), "Cluster size error")
        
        _render_section_divider("Structural Competition (UpSet)")
        _show_fig('upset_fig', lambda: plot_structural_competition_upset(motifs, title="Structural Competition"), "UpSet error")
        
        _render_section_divider("Overlaps Summary"); co, so = _calculate_overlaps(motifs, by='Class'), _calculate_overlaps(motifs, by='Subclass'); c1, c2 = st.columns(2)
        with c1:
            if co: _render_overlap_matrix(co, "Class Overlaps")
            else: st.markdown('<div style="color:#64748b;font-size:0.8rem;">No class overlaps</div>', unsafe_allow_html=True)
        with c2:
            if so: _render_overlap_matrix(so, "Subclass Overlaps")
            else: st.markdown('<div style="color:#64748b;font-size:0.8rem;">No subclass overlaps</div>', unsafe_allow_html=True)
