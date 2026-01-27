"""Download tab content for NonBDNAFinder"""

import streamlit as st
import pandas as pd
import numpy as np
import re
import io
import traceback
from collections import Counter
from openpyxl import Workbook

from config.text import UI_TEXT
from config.themes import TAB_THEMES
from ui.css import load_css
from ui.guards import generate_excel_bytes
from utilities import (
    export_to_csv,
    export_to_json,
    export_to_excel,
    export_to_pdf,
    export_to_bed
)
from export.export_validator import validate_export_data


def render():
    """Render the Download tab content"""
    # Apply Download tab theme based on configuration
    load_css(TAB_THEMES.get('Download', 'clinical_teal'))
    
    if not st.session_state.results:
        st.info(UI_TEXT['download_no_results'])
    else:
        # Display current Job ID prominently
        current_job_id = st.session_state.get('current_job_id')
        if current_job_id:
            st.markdown(f"""
            <div style='background: linear-gradient(135deg, #0ea5e9 0%, #0284c7 100%); 
                        color: white; padding: 1rem; border-radius: 12px; margin-bottom: 1.5rem;
                        box-shadow: 0 4px 12px rgba(14, 165, 233, 0.3);'>
                <div style='display: flex; justify-content: space-between; align-items: center;'>
                    <div>
                        <div style='font-size: 0.85rem; opacity: 0.9;'>Job ID for this analysis:</div>
                        <div style='font-size: 1.5rem; font-weight: 700; letter-spacing: 0.1em; 
                                   font-family: "Courier New", monospace; margin-top: 0.3rem;'>
                            {current_job_id}
                        </div>
                    </div>
                    <div style='text-align: right;'>
                        <div style='font-size: 0.85rem; opacity: 0.9;'>
                            💾 Results saved<br/>
                            🔍 Retrieve anytime with this ID
                        </div>
                    </div>
                </div>
            </div>
            """, unsafe_allow_html=True)
        
        primary_sequence_name = st.session_state.names[0] if st.session_state.names else "Unknown Sequence"
        analysis_time = st.session_state.get('analysis_time', 0)
        
        # Sanitize sequence name for safe filenames
        safe_filename = re.sub(r'[^\w\-]', '_', primary_sequence_name)[:50].strip('_')
        

        
        # Prepare motif data and validate
        all_motifs = []
        for i, motifs in enumerate(st.session_state.results):
            for m in motifs:
                export_motif = m.copy()
                if 'Sequence_Name' not in export_motif:
                    export_motif['Sequence_Name'] = st.session_state.names[i]
                all_motifs.append(export_motif)
        
        # Validate motifs before export (auto-normalize if needed)
        try:
            all_motifs = validate_export_data(all_motifs, auto_normalize=True, strict=False)
        except Exception as e:
            st.warning(f"⚠️ Some motifs had invalid class/subclass data and were normalized: {e}")
        
        # Individual file downloads as main option
        st.markdown("### 📥 Download Results")
        st.markdown("""
        <div style='background: #f0f9ff; padding: 1rem; border-radius: 8px; margin-bottom: 1rem;
                    border-left: 4px solid #0ea5e9;'>
            <p style='color: #0c4a6e; margin: 0;'>
                💡 <strong>Download your results</strong> in different file formats
            </p>
        </div>
        """, unsafe_allow_html=True)
        
        col1, col2, col3, col4, col5 = st.columns(5)
        
        with col1:
            # CSV Export (All Motifs)
            if all_motifs:
                csv_data = export_to_csv(all_motifs, non_overlapping_only=False)
                st.download_button(
                    "📋 CSV (All Motifs)", 
                    data=csv_data.encode('utf-8'), 
                    file_name=f"{safe_filename}_all_motifs.csv", 
                    mime="text/csv",
                    use_container_width=True,
                    type="primary",
                    help="Download CSV with all detected motifs including Hybrid and Clusters"
                )
        
        with col2:
            # Excel Export (2-tab format)
            if all_motifs:
                try:
                    excel_bytes = generate_excel_bytes(all_motifs, simple_format=True)
                    st.download_button(
                        "📊 Excel (2 tabs)", 
                        data=excel_bytes, 
                        file_name=f"{safe_filename}_results.xlsx", 
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        use_container_width=True,
                        type="primary",
                        help="Download Excel with 2 tabs: NonOverlappingConsolidated, OverlappingAll"
                    )
                except Exception as e:
                    st.error(f"Excel export error: {str(e)}")
        
        with col3:
            # JSON Export  
            if all_motifs:
                json_data = export_to_json(all_motifs, pretty=True)
                st.download_button(
                    "📄 JSON", 
                    data=json_data.encode('utf-8'), 
                    file_name=f"{safe_filename}_results.json", 
                    mime="application/json",
                    use_container_width=True,
                    type="primary",
                    help="Download results in JSON format"
                )
        
        with col4:
            # BED Export
            if all_motifs and st.session_state.names:
                bed_data = export_to_bed(all_motifs, st.session_state.names[0])
                st.download_button(
                    "🧬 BED", 
                    data=bed_data.encode('utf-8'), 
                    file_name=f"{safe_filename}_results.bed", 
                    mime="text/plain",
                    use_container_width=True,
                    type="primary",
                    help="Download results in BED format for genome browsers"
                )
        
        with col5:
            # PDF Export for visual summaries
            if all_motifs and st.session_state.seqs:
                try:
                    # Get sequence length from the first sequence (with safe None check)
                    sequence_length = len(st.session_state.seqs[0]) if st.session_state.seqs and st.session_state.seqs[0] else 0
                    if sequence_length > 0:
                        pdf_data = export_to_pdf(all_motifs, sequence_length, primary_sequence_name)
                        st.download_button(
                            "📄 PDF (Visualizations)",
                            data=pdf_data,
                            file_name=f"{safe_filename}_visualizations.pdf",
                            mime="application/pdf",
                            use_container_width=True,
                            type="primary",
                            help="Download PDF containing graphical summarizations of sequence analyses"
                        )
                    else:
                        st.warning("No sequence data available for PDF generation")
                except Exception as e:
                    st.error(f"PDF export error: {str(e)}")
        
        # Add Distribution & Statistics Tables Download Section
        st.markdown("---")
        st.markdown("### 📊 Download Distribution & Statistics Tables")
        st.markdown("""
        <div style='background: #f0fdf4; padding: 1rem; border-radius: 8px; margin-bottom: 1rem;
                    border-left: 4px solid #22c55e;'>
            <p style='color: #14532d; margin: 0;'>
                📈 <strong>Download detailed distribution and density statistics</strong> for publication and analysis
            </p>
        </div>
        """, unsafe_allow_html=True)
        
        # Calculate distribution and density statistics for all sequences
        if all_motifs and st.session_state.seqs:
            try:
                # Prepare distribution statistics
                distribution_data = []
                for seq_idx, (seq, name, motifs) in enumerate(zip(st.session_state.seqs, st.session_state.names, st.session_state.results)):
                    sequence_length = len(seq)
                    
                    # Show all motifs including hybrid/cluster motifs
                    # No filtering is applied - all results are included in statistics
                    filtered_motifs = motifs
                    
                    # Calculate class-level statistics
                    class_counts = Counter(m.get('Class', 'Unknown') for m in filtered_motifs)
                    for class_name, count in class_counts.items():
                        genomic_density = (sum(m.get('Length', 0) for m in filtered_motifs if m.get('Class') == class_name) / sequence_length * 100) if sequence_length > 0 else 0
                        motifs_per_kbp = (count / sequence_length * 1000) if sequence_length > 0 else 0
                        avg_length = np.mean([m.get('Length', 0) for m in filtered_motifs if m.get('Class') == class_name])
                        
                        distribution_data.append({
                            'Sequence Name': name,
                            'Motif Class': class_name.replace('_', ' '),
                            'Count': count,
                            'Genomic Density (%)': f"{genomic_density:.4f}",
                            'Motifs per kbp': f"{motifs_per_kbp:.2f}",
                            'Average Length (bp)': f"{avg_length:.1f}",
                            'Total Coverage (bp)': sum(m.get('Length', 0) for m in filtered_motifs if m.get('Class') == class_name)
                        })
                
                distribution_df = pd.DataFrame(distribution_data)
                
                # Prepare subclass-level statistics
                subclass_data = []
                for seq_idx, (seq, name, motifs) in enumerate(zip(st.session_state.seqs, st.session_state.names, st.session_state.results)):
                    sequence_length = len(seq)
                    
                    # Show all motifs including hybrid/cluster motifs
                    # No filtering is applied - all results are included in statistics
                    filtered_motifs = motifs
                    
                    # Calculate subclass-level statistics
                    subclass_counts = Counter(m.get('Subclass', 'Unknown') for m in filtered_motifs)
                    for subclass_name, count in subclass_counts.items():
                        parent_class = next((m.get('Class') for m in filtered_motifs if m.get('Subclass') == subclass_name), 'Unknown')
                        genomic_density = (sum(m.get('Length', 0) for m in filtered_motifs if m.get('Subclass') == subclass_name) / sequence_length * 100) if sequence_length > 0 else 0
                        motifs_per_kbp = (count / sequence_length * 1000) if sequence_length > 0 else 0
                        avg_length = np.mean([m.get('Length', 0) for m in filtered_motifs if m.get('Subclass') == subclass_name])
                        
                        subclass_data.append({
                            'Sequence Name': name,
                            'Motif Class': parent_class.replace('_', ' '),
                            'Motif Subclass': subclass_name.replace('_', ' '),
                            'Count': count,
                            'Genomic Density (%)': f"{genomic_density:.4f}",
                            'Motifs per kbp': f"{motifs_per_kbp:.2f}",
                            'Average Length (bp)': f"{avg_length:.1f}",
                            'Total Coverage (bp)': sum(m.get('Length', 0) for m in filtered_motifs if m.get('Subclass') == subclass_name)
                        })
                
                subclass_df = pd.DataFrame(subclass_data)
                
                # Display preview of tables
                st.markdown("#### Class-Level Distribution Statistics")
                st.dataframe(distribution_df.head(10), use_container_width=True, height=300)
                st.caption(f"Showing first 10 of {len(distribution_df)} total records")
                
                st.markdown("#### Subclass-Level Distribution Statistics")
                st.dataframe(subclass_df.head(10), use_container_width=True, height=300)
                st.caption(f"Showing first 10 of {len(subclass_df)} total records")
                
                # Download buttons for statistics tables
                col_stat1, col_stat2, col_stat3 = st.columns(3)
                
                with col_stat1:
                    # Class-level CSV
                    class_csv = distribution_df.to_csv(index=False)
                    st.download_button(
                        "📊 Class Statistics (CSV)",
                        data=class_csv.encode('utf-8'),
                        file_name=f"{safe_filename}_class_statistics.csv",
                        mime="text/csv",
                        use_container_width=True,
                        help="Download class-level distribution statistics"
                    )
                
                with col_stat2:
                    # Subclass-level CSV
                    subclass_csv = subclass_df.to_csv(index=False)
                    st.download_button(
                        "📊 Subclass Statistics (CSV)",
                        data=subclass_csv.encode('utf-8'),
                        file_name=f"{safe_filename}_subclass_statistics.csv",
                        mime="text/csv",
                        use_container_width=True,
                        help="Download subclass-level distribution statistics"
                    )
                
                with col_stat3:
                    # Combined Excel with both sheets
                    try:
                        output = io.BytesIO()
                        with pd.ExcelWriter(output, engine='openpyxl') as writer:
                            distribution_df.to_excel(writer, sheet_name='Class Statistics', index=False)
                            subclass_df.to_excel(writer, sheet_name='Subclass Statistics', index=False)
                        
                        output.seek(0)
                        st.download_button(
                            "📊 All Statistics (Excel)",
                            data=output.getvalue(),
                            file_name=f"{safe_filename}_all_statistics.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                            use_container_width=True,
                            help="Download all statistics in Excel format with separate sheets"
                        )
                    except Exception as e:
                        st.error(f"Excel generation error: {str(e)}")
                
            except Exception as e:
                st.error(f"Error generating distribution statistics: {str(e)}")
                st.code(traceback.format_exc(), language="python")
