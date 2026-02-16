UI_TEXT = {
    # ===== Application Metadata =====
    'app_title': 'NonBDNAFinder: Comprehensive Detection of Non-B DNA Forming Motifs',
    'author': 'Dr. Venkata Rajesh Yella',
    'author_email': 'yvrajesh_bt@kluniversity.in',
    'github_profile': 'VRYella',
    'github_repo': 'NonBFinder',
    'github_url': 'https://github.com/VRYella/NonBFinder',
    'version': '2024.1',

    # ===== Page Titles =====
    'home_title': 'Non-B DNA Motif Detection System',
    'upload_title': 'Sequence Upload and Motif Analysis',
    'results_title': 'Analysis Results and Visualization',
    'download_title': 'Export Results',
    'documentation_title': 'Scientific Documentation & References',

    # ===== Section Headers =====
    'section_scientific_foundation': 'Scientific Background',
    'section_motif_classes': 'Detected Motif Classes',
    'section_key_features': 'Key Capabilities',
    'section_how_to_cite': 'How to Cite',
    'section_sequence_upload': 'Sequence Input',
    'section_analysis_run': 'Analysis & Execution',
    'section_quick_options': 'Execution Options',
    'section_analysis_summary': 'Analysis Summary',
    'section_visualizations': 'Visualizations',
    'section_export_options': 'Export Options',
    'section_export_preview': 'Export Preview',
    'section_download_files': 'Download Files',
    'section_additional_exports': 'Additional Exports',

    # ===== Home Page Content =====
    'home_system_status_hyperscan':
        'Performance mode enabled: Hyperscan acceleration used for high-speed pattern matching',

    'home_system_status_standard':
        'Standard mode: Regex-based pattern matching (full functionality enabled)',

    'home_publication_ready':
        'Outputs suitable for publication-quality figures and tables',

    'home_call_to_action_title': 'Start an Analysis',

    'home_call_to_action_text':
        'Upload DNA sequences in FASTA format to identify and analyze Non-B DNA motifs',

    'home_call_to_action_button': 'Go to Upload & Analyze',

    'home_image_caption': 'Structural diversity of Non-B DNA conformations',
    'home_image_fallback_title': 'DNA',
    'home_image_fallback_subtitle': 'Non-B DNA Structures',
    'home_image_fallback_caption': 'Structural Motif Diversity',

    # ===== Upload & Analyze Page =====
    'upload_input_method_prompt': 'Select input source:',
    'upload_method_file': 'Upload FASTA File',
    'upload_method_paste': 'Paste FASTA Sequence',
    'upload_method_example': 'Example Data',
    'upload_method_ncbi': 'NCBI Fetch',

    'upload_file_prompt': 'Upload FASTA or multi-FASTA file',
    'upload_file_help': 'FASTA-formatted DNA sequences (single or multiple records)',
    'upload_processing': 'Parsing input sequences',
    'upload_file_valid': 'Valid FASTA',
    'upload_preview_button': 'Preview parsed sequences',
    'upload_no_sequences': 'No valid sequences detected in the input.',

    'upload_paste_prompt': 'Paste FASTA-formatted sequence(s):',
    'upload_paste_placeholder': 'Paste DNA sequences in FASTA format',
    'upload_paste_help': 'Supports single-sequence or multi-FASTA input',

    'upload_example_type_prompt': 'Example dataset:',
    'upload_example_single': 'Single sequence example',
    'upload_example_multi': 'Multi-FASTA example',
    'upload_example_help': 'Load example datasets for demonstration or testing',
    'upload_load_single_button': 'Load single example',
    'upload_load_multi_button': 'Load multi-FASTA example',
    'upload_example_single_success': 'Single example sequence loaded successfully.',
    'upload_example_multi_success': 'Multi-FASTA example loaded ({count} sequences).',

    'upload_ncbi_db_prompt': 'NCBI database',
    'upload_ncbi_help':
        'Only nucleotide and gene databases are suitable for DNA motif analysis',
    'upload_ncbi_query_prompt': 'Enter accession number or gene identifier:',
    'upload_ncbi_max_records': 'Maximum records',
    'upload_ncbi_fetch_button': 'Fetch sequences',
    'upload_ncbi_success': 'Successfully retrieved {count} sequences.',
    'upload_ncbi_error': 'NCBI retrieval failed: {error}',
    'upload_ncbi_empty_warning': 'Please provide a valid query before fetching.',

    'upload_quick_options_note':
        'All supported motif classes and subclasses are scanned by default',

    'upload_parallel_note':
        'Parallel scanning is recommended for large sequences (>100 kb)',

    'upload_run_analysis_button': 'Run Motif Detection',
    'upload_no_sequences_error':
        'Provide at least one valid DNA sequence before running analysis.',

    # ===== Analysis Progress Messages =====
    'progress_validating':
        'Validating detected motifs for consistency and redundancy',

    'progress_generating_viz':
        'Preparing class-level and subclass-level visualizations',

    'progress_all_detectors':
        'Running all motif detectors in parallel, followed by overlap resolution and clustering',

    'progress_parallel_fallback':
        'Parallel scanning failed; reverting to standard scanning mode ({error})',

    # ===== Status Messages =====
    'status_no_results':
        'No analysis results available. Run motif detection to proceed.',

    'status_analysis_ready': 'Ready for analysis',
    'status_analysis_running': 'Analysis in progress',
    'status_analysis_complete': 'Analysis completed successfully',

    'status_validation_passed':
        'Validation completed: no inconsistencies detected',

    'status_validation_issues':
        'Validation completed with {count} potential issues',

    'status_viz_prepared':
        'Visualization data prepared ({count} components, {time:.2f}s)',

    'status_pre_generated_ready':
        'Pre-generated results available: {classes} classes, {subclasses} subclasses',

    # ===== Results Page =====
    'results_performance_title': 'Performance Metrics',
    'results_no_motifs': 'No motifs detected for this sequence.',
    'results_sequence_selector': 'Select sequence for detailed results:',
    'results_hybrid_cluster_info':
        '{count} hybrid or cluster regions detected (see Cluster/Hybrid tab)',

    'results_columns_display': 'Displayed columns',
    'results_page_label': 'Page ({rows} rows per page)',
    'results_page_info': 'Displaying motifs {start}–{end} of {total}',
    'results_viz_title': 'Visualizations',

    # ===== Visualization Tabs =====
    'viz_tab_distribution': 'Distribution & Statistics',
    'viz_tab_coverage': 'Coverage & Density',
    'viz_tab_genome_wide': 'Genome-Wide Analysis',
    'viz_tab_hybrid_cluster': 'Cluster / Hybrid',

    # ===== Visualization Sections =====
    'viz_motif_distribution': 'Motif Distribution',
    'viz_statistical_analysis': 'Statistical Analysis',
    'viz_density_metrics': 'Density Metrics',
    'viz_distributions': 'Distributions',
    'viz_sequence_coverage': 'Sequence Coverage Analysis',
    'viz_circos_density': 'Circos Density Plot',

    'viz_genome_wide_title': 'Genome-Wide Motif Analysis',
    'viz_genome_wide_desc':
        'Genome-scale visualization of motif distribution across the full sequence',

    'viz_manhattan_plot': 'Manhattan Plot: Motif Density Hotspots',
    'viz_cumulative_dist': 'Cumulative Motif Distribution',

    'viz_linear_track': 'Linear Motif Track',
    'viz_linear_track_info':
        'Linear motif tracks are shown for sequences shorter than 50 kb',

    'viz_hybrid_cluster_title': 'Hybrid and Cluster Motif Analysis',

    'viz_advanced_title': 'Advanced Statistical Visualizations',
    'viz_advanced_desc':
        'Advanced visualizations suitable for in-depth analysis and manuscript figures',

    'viz_cooccurrence': 'Motif Co-occurrence Matrix',
    'viz_cooccurrence_desc':
        'Co-occurrence of motif classes within overlapping or adjacent regions',

    'viz_gc_correlation': 'GC Content vs Motif Density',
    'viz_gc_correlation_desc':
        'Relationship between GC content and motif density',

    'viz_length_kde': 'Motif Length Distribution',
    'viz_length_kde_desc':
        'Kernel density estimation of motif length distributions by class',

    'viz_cluster_size': 'Cluster Size and Diversity',
    'viz_cluster_size_desc':
        'Distribution of motif counts and class diversity within clusters',

    # ===== Analysis Section =====
    'analysis_section_title': 'Analysis & Execution',
    'analysis_quick_options_title': 'Execution Options',
    'analysis_run_button': 'Run Motif Detection',
    'analysis_run_button_disabled': 'Run Motif Detection (disabled)',
    'analysis_run_button_disabled_note':
        'Upload or paste a valid DNA sequence to continue',

    'analysis_pipeline_title': 'Analysis Pipeline',
    'analysis_progress_title': 'NonBDNAFinder Execution',

    # ===== Hybrid / Cluster =====
    'hybrid_all_tab': 'All',
    'hybrid_motifs_tab': 'Hybrid Motifs',
    'hybrid_clusters_tab': 'Cluster Motifs',

    'hybrid_motifs_title':
        'Hybrid Motifs (overlapping different motif classes)',

    'hybrid_motifs_info':
        'Hybrid motifs represent regions where different Non-B DNA classes overlap ({count} regions)',

    'hybrid_clusters_title':
        'Non-B DNA Clusters (high-density regions)',

    'hybrid_clusters_info':
        'Clusters are regions with a high density of multiple Non-B DNA motifs ({count} regions)',

    'hybrid_no_motifs':
        'No hybrid or cluster regions detected for this sequence.',

    'hybrid_no_hybrid': 'No hybrid motifs detected.',
    'hybrid_no_clusters': 'No cluster motifs detected.',

    # ===== Download Page =====
    'download_no_results': 'No results available for download.',
    'download_config_used': 'Analysis parameters applied:',
    'download_config_overlap': 'Overlap handling: {option}',
    'download_config_classes': 'Motif classes selected: {count}',
    'download_config_total': 'Total motifs detected: {count}',

    'download_export_config': 'Export Configuration',
    'download_include_sequences': 'Include full motif sequences',
    'download_include_help': 'Include nucleotide sequences in exported files',

    'download_excluded_info':
        '{count} hybrid or cluster motifs excluded based on configuration',

    'download_included_info':
        'All {count} hybrid or cluster motifs included',

    'download_preview_caption':
        'Showing first 10 of {count} records{excluded}',

    'download_excel_info_title': 'Excel Format',
    'download_excel_info':
        'Exports a multi-sheet workbook:\n'
        '• Consolidated sheet with non-overlapping motifs\n'
        '• Separate sheets for each motif class\n'
        '• Subclass-level detailed sheets',

    'download_excel_note':
        'Detailed columns appear only in class- and subclass-specific sheets',

    'download_button_csv': 'Download CSV',
    'download_button_excel': 'Download Excel',
    'download_button_json': 'Download JSON',
    'download_button_bed': 'Download BED',
    'download_button_config': 'Download Configuration',

    'download_help_csv': 'Comma-separated values',
    'download_help_excel': 'Excel workbook (multi-sheet)',
    'download_help_json': 'JSON with metadata',
    'download_help_bed': 'BED format for genome browsers (UCSC, IGV)',
    'download_help_config': 'Analysis parameters and metadata',

    'download_error_excel': 'Excel export failed: {error}',

    # ===== Documentation Page =====
    'doc_motif_classes_title': 'Detected Motif Classes',
    'doc_references_title': 'References',
    'doc_config_details': 'Scoring Configuration',
    'doc_length_constraints': 'Motif Length Constraints',
    'doc_scoring_methods': 'Scoring Methods',
    'doc_developed_by': 'Developed by',

    # ===== Tooltips =====
    'help_detailed_analysis': 'Include full motif metadata in output',
    'help_quality_validation': 'Check detected motifs for internal consistency',
    'help_parallel_scanner': 'Enable chunk-based parallel scanning for large sequences',
    'help_chunk_progress': 'Display per-chunk progress information',
    'help_sequence_selector': 'Select a sequence to view detailed results',

    # ===== Metrics =====
    'metric_sequence_coverage': 'Sequence Coverage',
    'metric_motif_density': 'Motif Density (motifs/kb)',
    'metric_total_motifs': 'Total Motifs',
    'metric_sequence_length': 'Sequence Length (bp)',
    'metric_processing_time': 'Processing Time',
    'metric_base_pairs': 'Base Pairs',
    'metric_speed': 'bp / second',
    'metric_detector_processes': 'Detector Processes',
    'metric_sequences': 'Sequences',
    'metric_total_motifs_found': 'Total Motifs',
    'metric_hybrid_motifs': 'Hybrid Motifs',
    'metric_dna_clusters': 'DNA Clusters',
    'metric_avg_length': 'Average Length (bp)',
    'metric_total': 'Total',

    # ===== Column Labels =====
    'col_motif_class': 'Motif Class',
    'col_motif_subclass': 'Motif Subclass',
    'col_genomic_density': 'Genomic Density (%)',
    'col_motifs_per_kbp': 'Motifs / kbp',
    'col_sequence_name': 'Sequence Name',
    'col_source': 'Source',
    'col_class': 'Class',
    'col_subclass': 'Subclass',
    'col_start': 'Start',
    'col_end': 'End',
    'col_length': 'Length',
    'col_sequence': 'Sequence',
    'col_score': 'Score',

    # ===== Generic Labels =====
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
    'analysis_no_sequences_warning': 'No valid sequences detected.',
    'analysis_all_detectors_parallel':
        'All motif detectors execute in parallel, followed by overlap resolution and clustering',

    # ===== Headings =====
    'heading_results_viz': 'Visualizations',
    'heading_analysis_results': 'Analysis Results',
    'heading_analysis_summary': 'Analysis Summary',
    'heading_motif_distribution': 'Motif Distribution',
    'heading_statistical_analysis': 'Statistical Analysis',
    'heading_density_metrics': 'Density Metrics',
    'heading_class_level': 'Class-Level Analysis',
    'heading_subclass_level': 'Subclass-Level Analysis',
    'heading_distributions': 'Distributions',
    'heading_sequence_coverage': 'Sequence Coverage',
    'heading_circos_density': 'Circos Density',
    'heading_genome_wide': 'Genome-Wide Analysis',
    'heading_hybrid_cluster': 'Hybrid and Cluster Analysis',
    'heading_advanced_viz': 'Advanced Visualizations',

    'heading_analysis_run': 'Analysis & Execution',

    'heading_export_data': 'Export Data',
    'heading_export_options': 'Export Options',
    'heading_export_preview': 'Export Preview',
    'heading_download_files': 'Download Files',
    'heading_export_config': 'Export Configuration',
    'heading_additional_exports': 'Additional Exports',

    'heading_documentation': 'Scientific Documentation & References',
    'heading_motif_classes': 'Detected Motif Classes',
    'heading_references': 'References',
    'heading_config_details': 'Scoring Configuration',
    'heading_length_constraints': 'Motif Length Constraints',
    'heading_scoring_methods': 'Scoring Methods',
}
