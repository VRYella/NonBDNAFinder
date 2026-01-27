"""
UI text content for NBDScanner.

This module contains all text strings used in the application UI:
- Application metadata (title, author, version, etc.)
- Page titles and section headers
- Content for each page
- Status messages and tooltips
- Column names and labels
"""

# ==================== UI TEXT CONTENT ====================
# All text strings for the application UI
# Edit these to customize labels, messages, and help text
UI_TEXT = {
    # ===== Application Metadata =====
    'app_title': 'NonBDNAFinder: A Comprehensive Non B-DNA forming motif detection system',
    'author': 'Dr. Venkata Rajesh Yella',
    'author_email': 'yvrajesh_bt@kluniversity.in',
    'github_profile': 'VRYella',
    'github_repo': 'NonBFinder',
    'github_url': 'https://github.com/VRYella/NonBFinder',
    'version': '2024.1',
    
    # ===== Page Titles =====
    'home_title': 'NonBDNA Motif Detection System',
    'upload_title': 'Sequence Upload and Motif Analysis',
    'results_title': 'Analysis Results and Visualization',
    'download_title': 'Export Data',
    'documentation_title': 'Scientific Documentation & References',
    
    # ===== Section Headers =====
    'section_scientific_foundation': 'Scientific Foundation',
    'section_motif_classes': 'Detected Motif Classes',
    'section_key_features': 'Key Features & Capabilities',
    'section_how_to_cite': 'How to Cite',
    'section_sequence_upload': 'Sequence Upload',
    'section_analysis_run': 'Analysis & Run',
    'section_quick_options': 'Quick Options',
    'section_analysis_summary': 'Analysis Summary',
    'section_visualizations': 'Visualizations',
    'section_export_options': 'Export Options',
    'section_export_preview': 'Export Preview',
    'section_download_files': 'Download Files',
    'section_additional_exports': 'Additional Exports',
    
    # ===== Home Page Content =====
    'home_system_status_hyperscan': 'Performance Mode: Hyperscan acceleration active for high-speed pattern matching',
    'home_system_status_standard': 'Standard Mode: Using regex-based pattern matching (all features fully functional)',
    'home_publication_ready': 'Publication Ready formats with 300 DPI resolution',
    'home_call_to_action_title': 'Ready to Analyze?',
    'home_call_to_action_text': 'Upload your FASTA sequences to begin comprehensive Non-B DNA motif detection',
    'home_call_to_action_button': '→ Go to "Upload & Analyze" tab',
    'home_image_caption': 'Non-B DNA Structural Diversity',
    'home_image_fallback_title': 'DNA',
    'home_image_fallback_subtitle': 'Non-B DNA Structures',
    'home_image_fallback_caption': 'Structural Diversity Database',
    
    # ===== Upload & Analyze Page =====
    'upload_input_method_prompt': 'Choose your input method:',
    'upload_method_file': 'Upload FASTA File',
    'upload_method_paste': 'Paste Sequence',
    'upload_method_example': 'Example Data',
    'upload_method_ncbi': 'NCBI Fetch',
    'upload_file_prompt': 'Drag and drop FASTA/multi-FASTA file here',
    'upload_file_help': 'Upload a FASTA file containing DNA sequences',
    'upload_processing': 'Processing',
    'upload_file_valid': 'Valid',
    'upload_preview_button': 'Preview Sequences',
    'upload_no_sequences': 'No sequences found in file.',
    'upload_paste_prompt': 'Paste single or multi-FASTA here:',
    'upload_paste_placeholder': 'Paste your DNA sequence(s) here...',
    'upload_paste_help': 'Paste DNA sequences in FASTA format',
    'upload_example_type_prompt': 'Example Type:',
    'upload_example_single': 'Single Example',
    'upload_example_multi': 'Multi-FASTA Example',
    'upload_example_help': 'Load example sequences for testing',
    'upload_load_single_button': 'Load Single Example',
    'upload_load_multi_button': 'Load Multi-FASTA Example',
    'upload_example_single_success': 'Single example sequence loaded.',
    'upload_example_multi_success': 'Multi-FASTA example loaded with {count} sequences.',
    'upload_ncbi_db_prompt': 'NCBI Database',
    'upload_ncbi_help': 'Only nucleotide and gene databases are applicable for DNA motif analysis',
    'upload_ncbi_query_prompt': 'Enter query (accession, gene, etc.):',
    'upload_ncbi_max_records': 'Max Records',
    'upload_ncbi_fetch_button': 'Fetch from NCBI',
    'upload_ncbi_success': 'Fetched {count} sequences.',
    'upload_ncbi_error': 'NCBI fetch failed: {error}',
    'upload_ncbi_empty_warning': 'Enter a query before fetching.',
    'upload_quick_options_note': 'All 11 motif classes with 22+ subclasses are detected automatically',
    'upload_parallel_note': 'Parallel scanner works best on sequences >100kb with multiple CPU cores',
    'upload_run_analysis_button': 'Run Complete Motif Analysis',
    'upload_no_sequences_error': 'Please upload or input sequences before running analysis.',
    
    # ===== Analysis Progress Messages =====
    'progress_validating': 'Validating results for consistency and quality...',
    'progress_generating_viz': 'Generating comprehensive visualizations for all classes and subclasses...',
    'progress_all_detectors': 'All detectors process in parallel | Followed by overlap resolution & clustering',
    'progress_parallel_fallback': 'Parallel scanner failed, falling back to standard: {error}',
    
    # ===== Status Messages =====
    'status_no_results': 'No analysis results. Please run motif analysis first.',
    'status_analysis_ready': 'Ready to analyze',
    'status_analysis_running': 'Analysis in progress...',
    'status_analysis_complete': 'Analysis complete!',
    'status_validation_passed': 'Validation passed: No consistency issues found',
    'status_validation_issues': 'Validation found {count} potential issues:',
    'status_viz_prepared': 'All visualizations prepared: {count} components in {time:.2f}s',
    'status_pre_generated_ready': 'Pre-generated Analysis Ready: {classes} unique classes, {subclasses} unique subclasses analyzed',
    
    # ===== Results Page =====
    'results_performance_title': 'Performance Metrics',
    'results_no_motifs': 'No motifs detected for this sequence.',
    'results_sequence_selector': 'Choose Sequence for Details:',
    'results_hybrid_cluster_info': '{count} Hybrid/Cluster motifs detected. View them in the \'Cluster/Hybrid\' tab below.',
    'results_columns_display': 'Columns to display',
    'results_page_label': 'Page (showing {rows} rows per page)',
    'results_page_info': 'Showing motifs {start} to {end} of {total}',
    'results_viz_title': 'Visualizations',
    
    # ===== Visualization Tab Names =====
    'viz_tab_distribution': 'Distribution & Statistics',
    'viz_tab_coverage': 'Coverage & Density',
    'viz_tab_genome_wide': 'Genome-Wide Analysis',
    'viz_tab_hybrid_cluster': 'Cluster/Hybrid',
    
    # ===== Visualization Section Headers =====
    'viz_motif_distribution': 'Motif Distribution',
    'viz_statistical_analysis': 'Statistical Analysis',
    'viz_density_metrics': 'Density Metrics',
    'viz_distributions': 'Distributions',
    'viz_sequence_coverage': 'Sequence Coverage Analysis',
    'viz_circos_density': 'Circos Density Plot',
    'viz_genome_wide_title': 'Genome-Wide Motif Analysis',
    'viz_genome_wide_desc': 'Publication-quality genome-scale visualizations showing motif distribution patterns across the entire sequence.',
    'viz_manhattan_plot': 'Manhattan Plot - Motif Density Hotspots',
    'viz_cumulative_dist': 'Cumulative Motif Distribution',
    'viz_linear_track': 'Linear Motif Track Viewer',
    'viz_linear_track_info': 'Linear motif track is shown for sequences < 50kb. For large sequences, use Manhattan plot or Coverage map.',
    'viz_hybrid_cluster_title': 'Hybrid & Cluster Motif Analysis',
    'viz_advanced_title': 'Advanced Statistical Visualizations',
    'viz_advanced_desc': 'Advanced publication-quality visualizations for in-depth analysis and manuscript figures.',
    'viz_cooccurrence': 'Motif Co-occurrence Matrix',
    'viz_cooccurrence_desc': 'Shows which motif classes tend to appear together (overlapping or within 1bp)',
    'viz_gc_correlation': 'GC Content vs Motif Density Correlation',
    'viz_gc_correlation_desc': 'Scatter plot showing relationship between GC content and motif density',
    'viz_length_kde': 'Motif Length Distribution (Kernel Density)',
    'viz_length_kde_desc': 'Smooth probability density curves showing length patterns by class',
    'viz_cluster_size': 'Cluster Size & Diversity Distribution',
    'viz_cluster_size_desc': 'Distribution of motif counts and class diversity within clusters',
    
    # ===== Analysis Section Headers =====
    'analysis_section_title': 'Analysis & Run',
    'analysis_quick_options_title': 'Quick Options',
    'analysis_run_button': 'Run NBDScanner Analysis',
    'analysis_run_button_disabled': 'Run NBDScanner Analysis (Disabled)',
    'analysis_run_button_disabled_note': 'Please upload or paste a valid sequence first',
    'analysis_pipeline_title': 'Analysis Pipeline',
    'analysis_progress_title': 'NonBFinder Analysis',
    
    # ===== Hybrid/Cluster Specific =====
    'hybrid_all_tab': 'All',
    'hybrid_motifs_tab': 'Hybrid Motifs',
    'hybrid_clusters_tab': 'Cluster Motifs',
    'hybrid_motifs_title': 'Hybrid Motifs (Overlapping Different Classes)',
    'hybrid_motifs_info': 'Hybrid motifs are regions where different Non-B DNA motif classes overlap. Found {count} hybrid regions.',
    'hybrid_clusters_title': 'Non-B DNA Clusters (High-Density Regions)',
    'hybrid_clusters_info': 'DNA Clusters are high-density regions with multiple Non-B DNA motif classes. Found {count} cluster regions.',
    'hybrid_no_motifs': 'No hybrid or cluster motifs detected in this sequence. Hybrid motifs occur when different Non-B DNA classes overlap, and clusters form when multiple motifs are found in close proximity.',
    'hybrid_no_hybrid': 'No hybrid motifs detected in this sequence.',
    'hybrid_no_clusters': 'No DNA clusters detected in this sequence.',
    
    # ===== Download Page =====
    'download_no_results': 'No results available to download.',
    'download_config_used': 'Analysis Configuration Used:',
    'download_config_overlap': 'Overlap Handling: {option}',
    'download_config_classes': 'Motif Classes: {count} classes selected',
    'download_config_total': 'Total Motifs Found: {count}',
    'download_export_config': 'Export Configuration',
    'download_include_sequences': 'Include Full Sequences',
    'download_include_help': 'Include full motif sequences in export',
    'download_excluded_info': '{count} Hybrid/Cluster motifs are excluded from downloads based on configuration. These are shown only in the Cluster/Hybrid visualization tab.',
    'download_included_info': 'All {count} Hybrid/Cluster motifs are included in downloads.',
    'download_preview_caption': 'Showing first 10 of {count} total records{excluded}',
    'download_excel_info_title': 'Excel Format',
    'download_excel_info': 'Downloads a multi-sheet workbook with:\n• Consolidated Sheet: All non-overlapping motifs with core columns\n• Class Sheets: Separate sheets for each motif class with ALL detailed columns\n• Subclass Sheets: Detailed breakdown by subclass with ALL detailed columns',
    'download_excel_note': 'Additional detailed columns are only shown in per-motif-class/subclass sheets, not in display tables.',
    'download_button_csv': 'Download CSV',
    'download_button_excel': 'Download Excel',
    'download_button_json': 'Download JSON',
    'download_button_bed': 'Download BED',
    'download_button_config': 'Download Config',
    'download_help_csv': 'Comma-separated values format',
    'download_help_excel': 'Excel format with multiple sheets (Consolidated + per-class)',
    'download_help_json': 'JSON format with metadata',
    'download_help_bed': 'BED format for genome browsers',
    'download_help_config': 'Analysis configuration and metadata',
    'download_error_excel': 'Excel export error: {error}',
    
    # ===== Documentation Page =====
    'doc_motif_classes_title': 'Motif Classes Detected:',
    'doc_references_title': 'References:',
    'doc_config_details': 'Scoring Configuration Details',
    'doc_length_constraints': 'Motif Length Constraints',
    'doc_scoring_methods': 'Scoring Methods',
    'doc_developed_by': 'Developed by',
    
    # ===== Tooltips and Help Text =====
    'help_detailed_analysis': 'Include comprehensive motif metadata in results',
    'help_quality_validation': 'Validate detected motifs for consistency',
    'help_parallel_scanner': 'Enable experimental parallel chunk-based scanner (>100kb sequences)',
    'help_chunk_progress': 'Display detailed progress for each processing chunk',
    'help_sequence_selector': 'Select a sequence to view detailed analysis results',
    
    # ===== Metric Labels =====
    'metric_sequence_coverage': 'Sequence Coverage',
    'metric_motif_density': 'Motif Density<br>(motifs/kb)',
    'metric_total_motifs': 'Total Motifs',
    'metric_sequence_length': 'Sequence Length (bp)',
    'metric_processing_time': 'Processing Time',
    'metric_base_pairs': 'Base Pairs',
    'metric_speed': 'bp/second',
    'metric_detector_processes': 'Detector Processes',
    'metric_sequences': 'Sequences',
    'metric_total_motifs_found': 'Total Motifs',
    'metric_hybrid_motifs': 'Hybrid Motifs',
    'metric_dna_clusters': 'DNA Clusters',
    'metric_avg_length': 'Avg Length (bp)',
    'metric_total': 'Total',
    
    # ===== Column Names (for display) =====
    'col_motif_class': 'Motif Class',
    'col_motif_subclass': 'Motif Subclass',
    'col_genomic_density': 'Genomic Density (%)',
    'col_motifs_per_kbp': 'Motifs/kbp',
    'col_sequence_name': 'Sequence Name',
    'col_source': 'Source',
    'col_class': 'Class',
    'col_subclass': 'Subclass',
    'col_start': 'Start',
    'col_end': 'End',
    'col_length': 'Length',
    'col_sequence': 'Sequence',
    'col_score': 'Score',
    
    # ===== Generic UI Elements =====
    'button_load': 'Load',
    'button_fetch': 'Fetch',
    'button_analyze': 'Analyze',
    'button_download': 'Download',
    'button_export': 'Export',
    'button_preview': 'Preview',
    'label_success': 'Success',
    'label_error': 'Error',
    'label_warning': 'Warning',
    'label_info': 'Info',
    'label_tip': 'Tip',
    'label_note': 'Note',
    'label_progress': 'Progress',
    'label_complete': 'Complete',
    'label_ready': 'Ready',
    'label_loading': 'Loading',
    'label_processing': 'Processing',
    'label_validating': 'Validating',
    'label_valid': 'Valid',
    
    # ===== Analysis Messages =====
    'analysis_no_sequences_warning': 'No sequences found.',
    'analysis_all_detectors_parallel': 'All detectors process in parallel | Followed by overlap resolution & clustering',
    
    # ===== Tooltips =====
    'tooltip_detailed_analysis': 'Include comprehensive motif metadata',
    'tooltip_quality_validation': 'Validate detected motifs',
    'tooltip_chunk_progress': 'Display detailed progress for each processing chunk',
    'tooltip_parallel_scanner': 'Enable experimental parallel chunk-based scanner (>100kb sequences)',
    
    # ===== Headings for Results Page =====
    'heading_results_viz': 'Visualizations',
    'heading_analysis_results': 'Analysis Results and Visualization',
    'heading_analysis_summary': 'Analysis Summary',
    'heading_motif_distribution': 'Motif Distribution',
    'heading_statistical_analysis': 'Statistical Analysis',
    'heading_density_metrics': 'Density Metrics',
    'heading_class_level': 'Class Level Analysis',
    'heading_subclass_level': 'Subclass Level Analysis',
    'heading_distributions': 'Distributions',
    'heading_sequence_coverage': 'Sequence Coverage Analysis',
    'heading_circos_density': 'Circos Density Plot',
    'heading_genome_wide': 'Genome-Wide Motif Analysis',
    'heading_hybrid_cluster': 'Hybrid & Cluster Motif Analysis',
    'heading_advanced_viz': 'Advanced Statistical Visualizations',
    
    # ===== Headings for Analysis Section =====
    'heading_analysis_run': 'Analysis & Run',
    
    # ===== Export Page Headings =====
    'heading_export_data': 'Export Data',
    'heading_export_options': 'Export Options',
    'heading_export_preview': 'Export Preview',
    'heading_download_files': 'Download Files',
    'heading_export_config': 'Export Configuration',
    'heading_additional_exports': 'Additional Exports',
    
    # ===== Documentation Page Headings =====
    'heading_documentation': 'Scientific Documentation & References',
    'heading_motif_classes': 'Motif Classes Detected:',
    'heading_references': 'References:',
    'heading_config_details': 'Scoring Configuration Details',
    'heading_length_constraints': 'Motif Length Constraints',
    'heading_scoring_methods': 'Scoring Methods',
}
