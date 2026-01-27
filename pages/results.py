"""
Results page for NonBDNAFinder application.
Displays analysis results with publication-quality visualizations.
"""

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import logging
from collections import Counter

from config.text import UI_TEXT
from config.themes import TAB_THEMES
from ui.css import load_css
from utilities import (
    get_basic_stats, 
    export_results_to_dataframe, 
    CORE_OUTPUT_COLUMNS,
    optimize_dataframe_memory,
    calculate_genomic_density, 
    calculate_positional_density,
    # Visualization functions
    plot_motif_distribution, 
    plot_nested_pie_chart,
    plot_manhattan_motif_density,
    plot_linear_motif_track,
    plot_density_comparison,
    plot_cluster_size_distribution,
    plot_motif_cooccurrence_matrix,
    plot_motif_length_kde,
    plot_score_distribution
)
from nonbscanner import get_motif_info as get_motif_classification_info
from visualization_standards import (
    NATURE_MOTIF_COLORS, 
    PlotDominance, 
    FigurePanel, 
    MetricFilter,
    LabelPolicy, 
    UILayout, 
    TRANSPARENCY_NOTE, 
    SUPPLEMENTARY_NOTE,
    should_show_plot, 
    get_nature_style_params
)

logger = logging.getLogger(__name__)


def render():
    """Render the Results tab content."""
    # Apply Results tab theme based on configuration
    load_css(TAB_THEMES.get('Results', 'genomic_purple'))
    st.markdown(f'<h2>{UI_TEXT["heading_analysis_results"]}</h2>', unsafe_allow_html=True)
    
    # Deterministic Results Page: Only render, never compute
    # If results are missing, show info and stop
    if not st.session_state.results:
        st.info(UI_TEXT['status_no_results'])
        st.info("💡 Run analysis first in the 'Upload & Analyze' tab")
        st.stop()  # Explicit stop to prevent any further execution
    
    # Performance metrics display if available
    if st.session_state.get('performance_metrics'):
        metrics = st.session_state.performance_metrics
        st.markdown(f"""
        <div class='progress-panel progress-panel--metrics'>
            <h3 class='progress-panel__title'>Performance Metrics</h3>
            <div class='stats-grid stats-grid--wide'>
                <div class='stat-card'>
                    <h2 class='stat-card__value'>{metrics['total_time']:.2f}s</h2>
                    <p class='stat-card__label'>Processing Time</p>
                </div>
                <div class='stat-card'>
                    <h2 class='stat-card__value'>{metrics['total_bp']:,}</h2>
                    <p class='stat-card__label'>Base Pairs</p>
                </div>
                <div class='stat-card'>
                    <h2 class='stat-card__value'>{metrics['speed']:,.0f}</h2>
                    <p class='stat-card__label'>bp/second</p>
                </div>
                <div class='stat-card'>
                    <h2 class='stat-card__value'>{metrics.get('detector_count', 9)}</h2>
                    <p class='stat-card__label'>Detector Processes</p>
                </div>
                <div class='stat-card'>
                    <h2 class='stat-card__value'>{metrics['sequences']}</h2>
                    <p class='stat-card__label'>Sequences</p>
                </div>
                <div class='stat-card'>
                    <h2 class='stat-card__value'>{metrics['total_motifs']}</h2>
                    <p class='stat-card__label'>Total Motifs</p>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    # Enhanced summary display
    st.markdown(f"### {UI_TEXT['heading_analysis_summary']}")
    st.dataframe(st.session_state.summary_df, use_container_width=True)
    
    # Sequence selection for detailed analysis using pills for better UX
    seq_idx = 0
    if len(st.session_state.seqs) > 1:
        # Use pills for sequence selection - a more modern and visual alternative to dropdown
        try:
            selected_seq = st.pills(
                "Choose Sequence for Details:",
                options=list(range(len(st.session_state.seqs))),
                format_func=lambda i: st.session_state.names[i],
                selection_mode="single",
                default=0,
                help="Select a sequence to view detailed analysis results"
            )
            seq_idx = selected_seq or 0
        except Exception:
            seq_idx = 0
    
    motifs = st.session_state.results[seq_idx]
    sequence_length = len(st.session_state.seqs[seq_idx])
    sequence_name = st.session_state.names[seq_idx]
    
    if not motifs:
        st.warning("No motifs detected for this sequence.")
    else:
        # Show all motifs including hybrid/cluster motifs
        # No filtering is applied - all results are displayed
        filtered_motifs = motifs
        hybrid_cluster_motifs = [m for m in motifs if m.get('Class') in ['Hybrid', 'Non-B_DNA_Clusters']]
        
        # Create enhanced motifs DataFrame
        df = pd.DataFrame(filtered_motifs) if filtered_motifs else pd.DataFrame()
        
        # Memory optimization: Optimize DataFrame for large result sets
        if len(df) > 1000:
            df = optimize_dataframe_memory(df)
            logger.debug(f"Optimized DataFrame memory for {len(df)} motifs")
        
        # Calculate and display enhanced coverage statistics (using filtered motifs)
        stats = get_basic_stats(st.session_state.seqs[seq_idx], filtered_motifs)
        
        motif_count = len(filtered_motifs)
        hybrid_cluster_count = len(hybrid_cluster_motifs)
        coverage_pct = stats.get("Motif Coverage %", 0)
        non_b_density = (motif_count / sequence_length * 1000) if sequence_length > 0 else 0
        
        # Enhanced summary card with modern research-quality styling
        st.markdown(f"""
        <div class='progress-panel progress-panel--results'>
            <h3 class='progress-panel__title progress-panel__title--large'>
                NBDScanner Analysis Results
            </h3>
            <div class='stats-grid stats-grid--extra-wide'>
                <div class='stat-card stat-card--large'>
                    <h2 class='stat-card__value stat_card_value_large'>
                        {stats.get("Coverage%", 0):.2f}%
                    </h2>
                    <p class='stat-card__label stat-card__label--large'>
                        Sequence Coverage
                    </p>
                </div>
                <div class='stat-card stat-card--large'>
                    <h2 class='stat-card__value stat-card__value--large'>
                        {stats.get("Density", 0):.2f}
                    </h2>
                    <p class='stat-card__label stat-card__label--large'>
                        Motif Density<br>(motifs/kb)
                    </p>
                </div>
                <div class='stat-card stat-card--large'>
                    <h2 class='stat-card__value stat-card__value--large'>
                        {motif_count}
                    </h2>
                    <p class='stat-card__label stat-card__label--large'>
                        Total Motifs
                    </p>
                </div>
                <div class='stat-card stat-card--large'>
                    <h2 class='stat-card__value stat-card__value--large'>
                        {sequence_length:,}
                    </h2>
                    <p class='stat-card__label stat-card__label--large'>
                        Sequence Length (bp)
                    </p>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        # Add info about hybrid/cluster motifs being shown separately
        if hybrid_cluster_count > 0:
            st.info(f"ℹ️ {hybrid_cluster_count} Hybrid/Cluster motifs detected. View them in the 'Cluster/Hybrid' tab below.")
        
        # Show cached visualization summary if available
        viz_cache_key = f"seq_{seq_idx}"
        cached_viz = st.session_state.get('cached_visualizations', {}).get(viz_cache_key, {})
        if cached_viz.get('summary'):
            viz_summary = cached_viz['summary']
            st.success(f"""✅ **Pre-generated Analysis Ready:** 
            {viz_summary['unique_classes']} unique classes, 
            {viz_summary['unique_subclasses']} unique subclasses analyzed
            """)
        
        # Enhanced motif table with canonical core reporting schema
        # Enhanced motif table with new columns and pagination for large datasets
        st.markdown("### 📋 All Detected Motifs")
        st.caption("**Core Reporting Schema** -- Publication-grade fields aligned with Nature/NAR/Genome Biology standards")
        
        # Info box explaining export-only advanced features
        st.info("""
        ℹ️ **Advanced Features:** All motif-specific details (ΔG components, dinucleotide counts, structural features, etc.) 
        are retained internally and available via exports (CSV/Excel/JSON). This view shows only publication-relevant fields 
        for clarity and interpretability.
        """)
        
        # Display only CORE_OUTPUT_COLUMNS (no column selection UI)
        # All advanced features remain available in exports
        display_columns = [col for col in CORE_OUTPUT_COLUMNS if col in df.columns]
        
        # Pagination for large datasets (improves performance)
        ROWS_PER_PAGE = 100
        total_rows = len(df)
        
        if total_rows > ROWS_PER_PAGE:
            # Add pagination controls
            col_pag1, col_pag2, col_pag3 = st.columns([2, 3, 2])
            with col_pag2:
                max_pages = (total_rows + ROWS_PER_PAGE - 1) // ROWS_PER_PAGE
                page_num = st.number_input(
                    f"Page (showing {ROWS_PER_PAGE} rows per page)",
                    min_value=1,
                    max_value=max_pages,
                    value=1,
                    step=1,
                    help=f"Total: {total_rows} motifs across {max_pages} pages"
                )
                start_idx = (page_num - 1) * ROWS_PER_PAGE
                end_idx = min(start_idx + ROWS_PER_PAGE, total_rows)
                
                st.caption(f"Showing motifs {start_idx + 1} to {end_idx} of {total_rows}")
            
            # Slice dataframe for current page
            df_page = df.iloc[start_idx:end_idx]
        else:
            df_page = df
        
        if display_columns:
            # Display core schema columns only
            filtered_df = df_page[display_columns].copy()
            # Replace underscores with spaces in column names for display
            filtered_df.columns = [col.replace('_', ' ') for col in filtered_df.columns]
            st.dataframe(filtered_df, use_container_width=True, height=360)
        else:
            # Fallback: display full dataframe if no core columns available
            display_df = df_page.copy()
            display_df.columns = [col.replace('_', ' ') for col in display_df.columns]
            st.dataframe(display_df, use_container_width=True, height=360)
        
        # NATURE-READY VISUALIZATION SUITE
        st.markdown(f'<h3>{UI_TEXT["heading_results_viz"]}</h3>', unsafe_allow_html=True)
        
        # Scientific Transparency Badge
        st.info(TRANSPARENCY_NOTE)
        
        # Create simplified visualization tabs based on Figure Panel layout
        viz_tabs = st.tabs([
            "📊 Figure 1: Global Landscape", 
            "🔗 Figure 2: Clustering & Co-occurrence",
            "📏 Figure 3: Structural Constraints (Optional)"
        ])
        
        # Check if clusters exist
        has_clusters = any(m.get('Class') == 'Non-B_DNA_Clusters' for m in filtered_motifs)
        
        # =================================================================
        # FIGURE 1: Global Non-B DNA Landscape
        # =================================================================
        with viz_tabs[0]:
            st.markdown("#### Figure 1: Global Non-B DNA Landscape")
            st.caption("*Purpose: What structures exist, and where?*")
            
            # Panel A: Motif Composition (nested donut + bar plots)
            st.markdown("##### Panel A: Motif Composition (Class → Subclass)")
            
            # Nested pie chart for hierarchical view
            try:
                fig_composition = plot_nested_pie_chart(
                    filtered_motifs, 
                    title=f"Motif Composition - {sequence_name}"
                )
                st.pyplot(fig_composition)
                plt.close(fig_composition)
            except Exception as e:
                st.error(f"Error generating composition plot: {e}")
            
            # Bar plots for class and subclass distributions
            col1, col2 = st.columns(2)
            with col1:
                try:
                    fig_class_bar = plot_motif_distribution(
                        filtered_motifs,
                        by='Class',
                        title=f"Motif Class Distribution - {sequence_name}"
                    )
                    st.pyplot(fig_class_bar)
                    plt.close(fig_class_bar)
                except Exception as e:
                    st.error(f"Error generating class bar plot: {e}")
            
            with col2:
                try:
                    fig_subclass_bar = plot_motif_distribution(
                        filtered_motifs,
                        by='Subclass',
                        title=f"Motif Subclass Distribution - {sequence_name}"
                    )
                    st.pyplot(fig_subclass_bar)
                    plt.close(fig_subclass_bar)
                except Exception as e:
                    st.error(f"Error generating subclass bar plot: {e}")
            
            # Panel B: Genome-Scale Localization (size-dependent: Manhattan OR Linear)
            st.markdown("##### Panel B: Genome-Scale Localization")
            try:
                if sequence_length > 50000:
                    # Large sequences: Manhattan plot
                    st.caption("*Manhattan plot for large genome (>50kb)*")
                    fig_position = plot_manhattan_motif_density(
                        filtered_motifs, sequence_length,
                        title=f"Motif Density Hotspots - {sequence_name}"
                    )
                    st.pyplot(fig_position)
                    plt.close(fig_position)
                else:
                    # Small sequences: Linear track
                    st.caption("*Linear track for short sequence (≤50kb)*")
                    fig_position = plot_linear_motif_track(
                        filtered_motifs, sequence_length,
                        title=f"Motif Track - {sequence_name}"
                    )
                    st.pyplot(fig_position)
                    plt.close(fig_position)
            except Exception as e:
                st.error(f"Error generating localization plot: {e}")
            
            # Panel C: Genome Coverage (compact bar chart)
            st.markdown("##### Panel C: Genome Coverage (% per motif class)")
            try:
                # Calculate or retrieve density metrics
                viz_cache_key = f"seq_{seq_idx}"
                cached_viz = st.session_state.get('cached_visualizations', {}).get(viz_cache_key, {})
                cached_densities = cached_viz.get('densities', {})
                
                if cached_densities:
                    genomic_density = cached_densities['class_genomic']
                    positional_density_kbp = cached_densities['class_positional']
                else:
                    genomic_density = calculate_genomic_density(filtered_motifs, sequence_length, by_class=True)
                    positional_density_kbp = calculate_positional_density(filtered_motifs, sequence_length, unit='kbp', by_class=True)
                
                # Compact visualization: density comparison
                fig_density = plot_density_comparison(
                    genomic_density, positional_density_kbp,
                    title="Motif Density Analysis"
                )
                st.pyplot(fig_density)
                plt.close(fig_density)
                
                # Compact density table
                density_data = []
                for class_name in sorted([k for k in genomic_density.keys() if k != 'Overall']):
                    density_data.append({
                        'Motif Class': class_name.replace('_', ' '),
                        'Coverage (%)': f"{genomic_density.get(class_name, 0):.3f}",
                        'Motifs/kb': f"{positional_density_kbp.get(class_name, 0):.2f}"
                    })
                
                if density_data:
                    density_df = pd.DataFrame(density_data)
                    st.dataframe(density_df, use_container_width=True, height=200)
                    
            except Exception as e:
                st.error(f"Error calculating density metrics: {e}")
        
        # =================================================================
        # FIGURE 2: Structural Clustering & Co-Occurrence
        # =================================================================
        with viz_tabs[1]:
            st.markdown("#### Figure 2: Structural Clustering & Co-Occurrence")
            st.caption("*Purpose: Which structures co-localize and form regulatory hotspots?*")
            
            # Panel D: Cluster Size Distribution (conditional on cluster existence)
            if has_clusters:
                st.markdown("##### Panel D: Cluster Size Distribution")
                try:
                    fig_cluster = plot_cluster_size_distribution(
                        filtered_motifs,
                        title=f"Cluster Statistics - {sequence_name}"
                    )
                    st.pyplot(fig_cluster)
                    plt.close(fig_cluster)
                except Exception as e:
                    st.error(f"Error generating cluster plot: {e}")
            else:
                st.info("**Panel D**: No clusters detected (requires multiple motifs in close proximity)")
            
            # Panel E: Motif Co-occurrence Matrix (always shown)
            st.markdown("##### Panel E: Motif Co-occurrence Matrix")
            st.caption("*Shows which motif classes tend to appear together (overlapping or within 1bp)*")
            try:
                fig_cooccur = plot_motif_cooccurrence_matrix(
                    filtered_motifs,
                    title=f"Co-occurrence Matrix - {sequence_name}"
                )
                st.pyplot(fig_cooccur)
                plt.close(fig_cooccur)
            except Exception as e:
                st.error(f"Error generating co-occurrence matrix: {e}")
        
        # =================================================================
        # FIGURE 3: Structural Constraints (Optional/Toggle)
        # =================================================================
        with viz_tabs[2]:
            st.markdown("#### Figure 3: Structural Constraints")
            st.caption("*Purpose: What are the physical constraints shaping each motif class?*")
            st.info("📌 **Optional Figure**: Toggle to include in main report or move to supplementary materials.")
            
            # User toggle for including this figure
            show_fig3 = st.checkbox(
                "Include Figure 3 in main report", 
                value=True,
                help="Enable to show structural constraint analysis in main figures"
            )
            
            if show_fig3:
                # Panel F: Length Distribution by Class (KDE only - no histogram redundancy)
                st.markdown("##### Panel F: Length Distributions by Class")
                try:
                    fig_length = plot_motif_length_kde(
                        filtered_motifs,
                        by_class=True,
                        title=f"Length Distribution (KDE) - {sequence_name}"
                    )
                    st.pyplot(fig_length)
                    plt.close(fig_length)
                except Exception as e:
                    st.error(f"Error generating length KDE: {e}")
                
                # Optional: Score distribution
                st.markdown("**Additional: Score Distribution by Class**")
                try:
                    fig_score = plot_score_distribution(
                        filtered_motifs, by_class=True,
                        title="Score Distribution (1-3 Scale)"
                    )
                    st.pyplot(fig_score)
                    plt.close(fig_score)
                except Exception as e:
                    st.error(f"Error generating score distribution: {e}")
            else:
                st.info("✓ Figure 3 hidden from main report (available in supplementary exports)")
                st.caption(SUPPLEMENTARY_NOTE)
