"""Documentation Page - Scientific Documentation & References."""

from typing import Optional

import pandas as pd
import streamlit as st

from Utilities.config.analysis import ANALYSIS_CONFIG
from Utilities.config.motif_taxonomy import MOTIF_CLASSIFICATION
from Utilities.config.text import UI_TEXT
from Utilities.config.themes import TAB_THEMES
from UI.css import get_page_colors, load_css

GITHUB_URL = UI_TEXT["github_url"]
MAX_SEQUENCES = ANALYSIS_CONFIG.get("max_sequences_per_upload", 25)
MAX_FILE_SIZE_MB = ANALYSIS_CONFIG.get("max_file_size_mb", 10)
CITATION_TEXT = (
    "Yella VR (2025). NonBDNAFinder: Comprehensive Detection and Analysis of "
    "Non-B DNA Forming Motifs. GitHub repository: https://github.com/VRYella/NonBDNAFinder"
)

MOTIF_LIBRARY = [
    {
        "class": "Curved DNA",
        "principle": "A/T-tract phasing and local tract detection",
        "scoring": "Tract length and phasing-based normalized confidence",
        "subclasses": "Global Curvature; Local Curvature",
        "details": "Detects phased A/T-tract organization associated with intrinsic curvature.",
    },
    {
        "class": "Slipped DNA",
        "principle": "k-mer repeat scanning for STRs and direct repeats",
        "scoring": "Repeat length, copy number, purity, entropy, and GC-aware scoring",
        "subclasses": "Direct Repeat; STR",
        "details": "Captures repeat architectures associated with strand slippage.",
    },
    {
        "class": "Cruciform",
        "principle": "Seed-and-extend inverted repeat detection",
        "scoring": "Arm length, symmetry, GC contribution, and ΔG filtering",
        "subclasses": "Cruciform forming IRs",
        "details": "Reports palindromic inverted repeats compatible with cruciform extrusion.",
    },
    {
        "class": "R-Loop",
        "principle": "QmRLFS-related initiation and elongation zone detection",
        "scoring": "GC content and G-density normalized confidence",
        "subclasses": "R-loop formation sites",
        "details": "Models RNA-DNA hybrid-forming sequence patterns using RIZ/REZ logic.",
    },
    {
        "class": "Triplex",
        "principle": "Mirror-repeat detection with Sticky DNA repeat detection",
        "scoring": "Arm length, loop tolerance, purity, and interruption-aware scoring",
        "subclasses": "H-DNA; Sticky DNA",
        "details": "Separates mirror-repeat triplex candidates from GAA/TTC expansion motifs.",
    },
    {
        "class": "G-Quadruplex",
        "principle": "Seeded G4Hunter-derived scanning plus subclass pattern matching",
        "scoring": "G4Hunter-derived normalized confidence",
        "subclasses": "Telomeric G4; Stacked G4; Canonical intramolecular G4; Extended-loop canonical; Higher-order G4 array/G4-wire; Intramolecular G-triplex; Two-tetrad weak PQS; Bulged G4",
        "details": "Supports canonical and higher-order G-rich structural subclasses.",
    },
    {
        "class": "i-Motif",
        "principle": "C-tract pattern matching plus AC-motif detection",
        "scoring": "C-density, tract regularity, and loop-aware normalized confidence",
        "subclasses": "Canonical i-motif; Relaxed i-motif; AC-motif",
        "details": "Captures canonical and relaxed i-motif patterns alongside AC-motif signatures.",
    },
    {
        "class": "Z-DNA",
        "principle": "10-mer propensity scoring with eGZ repeat detection",
        "scoring": "Cumulative Z-DNA propensity and eGZ repeat-based normalized confidence",
        "subclasses": "Z-DNA; eGZ",
        "details": "Combines adapted Z-DNA propensity scoring with extruded-guanine Z-DNA detection.",
    },
    {
        "class": "A-philic DNA",
        "principle": "10-mer A-form propensity scoring",
        "scoring": "Merged cumulative A-form propensity scoring",
        "subclasses": "A-philic DNA",
        "details": "Reports sequence-derived propensity toward A-form geometry rather than experimental confirmation.",
    },
    {
        "class": "Hybrid",
        "principle": "Inter-class overlap annotation",
        "scoring": "Composite overlap summary score",
        "subclasses": "Dynamic overlaps",
        "details": "Represents overlapping prediction criteria from distinct motif classes.",
    },
    {
        "class": "Non-B DNA Clusters",
        "principle": "Density-based hotspot detection",
        "scoring": "Cluster-level aggregation summary score",
        "subclasses": "Dynamic clusters",
        "details": "Flags dense windows containing multiple motifs from multiple classes.",
    },
]

ALGORITHM_SUMMARIES = {
    "Curved DNA": "Curved DNA is identified from phased A/T-tract organization and long local tracts. The implementation reports global curvature and local curvature subclasses.",
    "Slipped DNA": "Slipped DNA detection scans repeat architectures for STRs and direct repeats and assigns normalized confidence from repeat composition, repeat count, and tract characteristics.",
    "Cruciform": "Cruciform detection uses inverted-repeat search with thermodynamic filtering to retain structurally plausible candidates.",
    "R-Loop": "R-loop prediction follows QmRLFS-related sequence logic by combining guanine-rich initiation zones with downstream elongation support.",
    "Triplex": "Triplex analysis distinguishes H-DNA mirror-repeat candidates from Sticky DNA repeat expansions using subclass-specific rules already implemented in the repository.",
    "G-Quadruplex": "G-quadruplex analysis uses G4Hunter-derived scoring together with subclass resolution for telomeric, canonical, bulged, higher-order, stacked, extended-loop, triplex-like, and weak PQS patterns.",
    "i-Motif": "i-motif analysis reports canonical i-motif, relaxed i-motif, and AC-motif subclasses using C-tract and alternating-pattern logic.",
    "Z-DNA": "Z-DNA analysis applies cumulative 10-mer propensity scoring and separately reports eGZ repeat candidates.",
    "A-philic DNA": "A-philic DNA analysis estimates sequence-derived propensity toward A-form geometry from positive 10-mer propensity regions.",
}

REFERENCE_LIST = [
    "Bedrat A, Lacroix L, Mergny JL. Re-evaluation of G-quadruplex propensity with G4Hunter. Nucleic Acids Research (2016).",
    "Jenjaroenpun P, Wongsurawat T, Sutheeworapong S, et al. QmRLFS-finder: a model, web server and stand-alone tool. Nucleic Acids Research (2015/2016).",
    "Frank-Kamenetskii MD, Mirkin SM. Triplex DNA structures. Annual Review of Biochemistry (1995).",
    "Ho PS, Frederick CA, Saal D, et al. Z-DNA propensity framework underlying the implemented 10-mer approach. Journal of Biomolecular Structure and Dynamics (1986).",
    "Gehring K, Leroy JL, Gueron M and later i-motif studies including Zeraati M, Langley DB, et al. Nature Chemistry (2018).",
    "Crothers DM, Koo HS, Olson WK and related curvature studies used for A/T-tract phasing interpretation.",
    "Vinogradov AE and related conformational propensity studies underlying the implemented A-philic scoring approach.",
]

FAQS = [
    (
        "What are the public web application limits?",
        f"The public web application currently supports a maximum of {MAX_SEQUENCES} sequences per upload and a maximum genome/input size of {MAX_FILE_SIZE_MB} MB.",
    ),
    (
        "Can I analyze larger genomes or datasets?",
        f"Yes. For larger genomes or datasets exceeding the web application limits, use the local Streamlit application, the included Jupyter notebook, or the Python API described in {GITHUB_URL}.",
    ),
    (
        "Are confidence scores directly comparable across motif classes?",
        "No. Normalized confidence values are primarily intended for relative ranking within the relevant structural class or subclass and should not be interpreted as directly equivalent measures of structural stability across different motif classes.",
    ),
    (
        "What does a hybrid overlap mean?",
        "Hybrid overlap represents overlapping sequence-level prediction criteria and does not by itself demonstrate simultaneous physical formation of multiple DNA structures at the same locus. G-quadruplex/i-motif overlap should be interpreted as co-localized structural potential on complementary strands.",
    ),
    (
        "How should A-philic DNA predictions be interpreted?",
        "A-philic DNA predictions represent sequence-derived conformational propensity and should not be interpreted as direct experimental confirmation of A-DNA formation.",
    ),
    (
        "Are the predictions experimental confirmation?",
        "No. NonBDNAFinder provides sequence-based computational predictions. Biological validation still depends on the relevant experimental context.",
    ),
]


def _colors():
    return get_page_colors("Documentation")


def _render_header():
    c = _colors()
    st.markdown(
        f"""
        <div style='background:linear-gradient(135deg,#f0f9ff 0%,#e0f2fe 55%,#ecfeff 100%);border:1px solid {c['border']};border-radius:18px;padding:1.15rem 1.25rem 1rem 1.25rem;margin-bottom:1rem;'>
            <div style='color:{c['primary']};font-size:0.8rem;font-weight:700;letter-spacing:0.08em;text-transform:uppercase;margin-bottom:0.35rem;'>Documentation</div>
            <div style='color:#0f172a;font-size:1.9rem;font-weight:800;line-height:1.15;margin-bottom:0.25rem;'>Scientific Documentation</div>
            <div style='color:#475569;font-size:1rem;line-height:1.55;'>Methods, algorithms, analysis, usage, and references</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def _section_title(title: str, subtitle: Optional[str] = None):
    c = _colors()
    subtitle_html = (
        f"<div style='color:#64748b;font-size:0.92rem;margin-top:0.15rem;'>{subtitle}</div>"
        if subtitle
        else ""
    )
    st.markdown(
        f"""
        <div style='margin:0.15rem 0 0.9rem 0;'>
            <div style='color:#0f172a;font-size:1.22rem;font-weight:800;line-height:1.2;'>{title}</div>
            <div style='width:4.2rem;height:2px;background:{c['primary']};margin-top:0.35rem;border-radius:999px;'></div>
            {subtitle_html}
        </div>
        """,
        unsafe_allow_html=True,
    )


def _callout(title: str, body: str, tone: str = "info"):
    tones = {
        "info": ("#eff6ff", "#bfdbfe", "#0369a1"),
        "success": ("#ecfdf5", "#bbf7d0", "#047857"),
        "warning": ("#fff7ed", "#fed7aa", "#c2410c"),
    }
    bg, border, text = tones[tone]
    st.markdown(
        f"""
        <div style='background:{bg};border:1px solid {border};border-left:4px solid {text};border-radius:12px;padding:0.85rem 0.95rem;margin:0.7rem 0;'>
            <div style='color:{text};font-weight:700;font-size:0.95rem;margin-bottom:0.25rem;'>{title}</div>
            <div style='color:#334155;font-size:0.92rem;line-height:1.6;'>{body}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def _workflow_diagram():
    stages = [
        "Input FASTA / NCBI Fetch",
        "Sequence Processing",
        "Motif Detection",
        "Motif-specific Algorithms",
        "Scoring / Classification",
        "Hybrid / Cluster Analysis",
        "Results & Visualization",
        "Download / Export",
    ]
    items = []
    for idx, stage in enumerate(stages):
        items.append(
            "<div style='display:flex;align-items:center;justify-content:center;min-height:54px;padding:0.7rem 0.9rem;background:#ffffff;border:1px solid #dbeafe;border-radius:14px;color:#0f172a;font-weight:700;text-align:center;'>"
            + stage
            + "</div>"
        )
        if idx < len(stages) - 1:
            items.append(
                "<div style='text-align:center;color:#0ea5e9;font-size:1.1rem;font-weight:700;'>↓</div>"
            )
    st.markdown(
        "<div style='background:#f8fbff;border:1px solid #dbeafe;border-radius:18px;padding:1rem;display:grid;gap:0.45rem;'>"
        + "".join(items)
        + "</div>",
        unsafe_allow_html=True,
    )


def _overview_section():
    _section_title(
        "Overview",
        "Repository-aligned summary of the implemented workflow and scientific scope.",
    )
    st.markdown(
        "NonBDNAFinder integrates nine primary sequence-based detector families and two derived summary classes to identify non-B DNA-forming sequence motifs, annotate subclasses, normalize class-specific confidence scores, summarize hybrid overlaps, identify cluster hotspots, and provide publication-oriented visualizations and exports."
    )
    _workflow_diagram()
    col1, col2, col3 = st.columns(3)
    col1.metric("Primary detector classes", "9")
    col2.metric("Derived output classes", "2")
    col3.metric(
        "Canonical subclasses",
        str(sum(len(entry["subclasses"]) for entry in MOTIF_CLASSIFICATION.values())),
    )
    _callout(
        "Interpretation and limitations",
        "These predictions are sequence-based computational calls. They do not automatically capture DNA supercoiling, transcriptional state, chromatin state, DNA methylation, or histone modifications unless those contextual variables are explicitly incorporated by the current implementation. Computational prediction is not experimental confirmation.",
        "warning",
    )


def _motif_library_section():
    _section_title(
        "Motif Library",
        "Motif classes, detection principles, scoring summaries, and subclasses used by the implementation.",
    )
    df = pd.DataFrame(MOTIF_LIBRARY)
    st.dataframe(
        df[["class", "principle", "scoring", "subclasses"]],
        use_container_width=True,
        hide_index=True,
    )
    for item in MOTIF_LIBRARY:
        with st.expander(item["class"], expanded=False):
            st.write(item["details"])


def _algorithms_section():
    _section_title(
        "Algorithms",
        "Concise explanations of the implemented detector logic without changing any underlying methods.",
    )
    for name, summary in ALGORITHM_SUMMARIES.items():
        with st.expander(name, expanded=False):
            st.write(summary)
    _callout(
        "Implemented workflow note",
        "The application uses class-specific detectors first, then performs score normalization, overlap-aware hybrid annotation, and density-based cluster detection. This documentation reflects the current repository implementation and does not describe future or hypothetical workflows.",
    )


def _scoring_section():
    _section_title(
        "Scoring & Analysis",
        "How normalized confidence, hybrid overlap, and A-philic terminology should be interpreted.",
    )
    st.markdown("**Universal normalized confidence scale:** 1.0–3.0")
    st.code("Score = 1.0 + 2.0 × (raw − floor) / (ceiling − floor)", language="text")
    st.markdown(
        "Normalized confidence values are primarily intended for relative ranking within the relevant structural class or subclass. They should not be interpreted as directly equivalent measures of structural stability across different motif classes."
    )
    _callout(
        "Hybrid prediction interpretation",
        "Hybrid overlap represents overlapping sequence-level prediction criteria and does not by itself demonstrate simultaneous physical formation of multiple DNA structures at the same locus. For G-quadruplex/i-motif overlap, interpret the result as co-localized structural potential on complementary strands.",
        "info",
    )
    _callout(
        "A-philic DNA terminology",
        "A-philic DNA predictions represent sequence-derived conformational propensity and should not be interpreted as direct experimental confirmation of A-DNA formation. Use 'A-form DNA propensity' or 'propensity toward A-form geometry' when discussing these calls scientifically.",
        "info",
    )


def _statistics_section():
    _section_title("Statistics", "What the reported summary metrics and plots represent.")
    cols = st.columns(4)
    metric_cards = [
        ("Coverage (%)", "Fraction of unique sequence bases covered by predicted motifs."),
        ("Motif density", "Number of motif calls per kilobase of analyzed sequence."),
        (
            "Class contribution",
            "Relative contribution of each motif class to the motif set or coverage footprint.",
        ),
        (
            "Hybrid / cluster counts",
            "Counts of derived overlap and hotspot annotations reported separately from primary classes.",
        ),
    ]
    for col, (title, desc) in zip(cols, metric_cards):
        col.markdown(
            f"<div style='background:#ffffff;border:1px solid #e2e8f0;border-radius:14px;padding:0.85rem;height:100%;'><div style='font-weight:700;color:#0f172a;margin-bottom:0.25rem;'>{title}</div><div style='color:#475569;font-size:0.9rem;line-height:1.55;'>{desc}</div></div>",
            unsafe_allow_html=True,
        )
    st.markdown("")
    stats_df = pd.DataFrame(
        [
            {
                "Metric": "Coverage (%)",
                "Represents": "Unique sequence fraction covered by predicted motifs",
            },
            {
                "Metric": "Total motifs",
                "Represents": "Absolute number of reported motif intervals",
            },
            {
                "Metric": "Simpson diversity / effective class number",
                "Represents": "How broadly motif counts are distributed across classes",
            },
            {
                "Metric": "Length and score distributions",
                "Represents": "How motif span and normalized confidence vary across classes",
            },
            {
                "Metric": "Density and track plots",
                "Represents": "Where motifs accumulate across the analyzed sequence",
            },
        ]
    )
    st.dataframe(stats_df, use_container_width=True, hide_index=True)
    _callout(
        "Statistics guide",
        "The Statistics section improves presentation only. It does not recalculate metrics or alter existing numerical results.",
        "success",
    )


def _user_guide_section():
    _section_title(
        "User Guide",
        "Input requirements, workflow steps, application limits, exports, and local processing guidance.",
    )
    st.markdown(
        """
1. Upload FASTA, paste FASTA, load example data, or fetch a sequence from NCBI.
2. Select the motif classes and subclasses to include.
3. Run the existing analysis options.
4. Review results, visualizations, and summary statistics.
5. Download the generated exports.
        """
    )
    st.markdown("**Input requirements**")
    st.markdown(
        "- FASTA headers identify each sequence record.\n"
        "- Sequence names are taken from the FASTA header.\n"
        "- Nucleotide sequences accept IUPAC DNA characters; pasted RNA is normalized by converting U to T.\n"
        "- Results use 1-based inclusive coordinates."
    )
    limits_df = pd.DataFrame(
        [
            {"Limit": "Sequences per upload", "Maximum": str(MAX_SEQUENCES)},
            {"Limit": "Genome/input size", "Maximum": f"{MAX_FILE_SIZE_MB} MB"},
        ]
    )
    st.markdown("**Application Limits**")
    st.table(limits_df)
    _callout(
        "Larger datasets",
        f"Need to analyze larger genomes or datasets? Please visit the Non-B DNA Finder GitHub repository for local installation and processing instructions: <a href='{GITHUB_URL}' target='_blank' style='color:#0369a1;font-weight:700;'>{GITHUB_URL}</a>. The repository includes a local Streamlit application, a Jupyter notebook, and Python API entry points.",
        "warning",
    )
    exports_df = pd.DataFrame(
        [
            {"Export": "CSV", "Availability": "Web application download"},
            {"Export": "Excel (XLSX)", "Availability": "Web application download"},
            {"Export": "JSON", "Availability": "Web application download"},
            {"Export": "BED", "Availability": "Web application download"},
            {"Export": "GFF3", "Availability": "Web application download"},
            {
                "Export": "PDF",
                "Availability": "Web application download for visualization package",
            },
        ]
    )
    st.markdown("**Output formats**")
    st.dataframe(exports_df, use_container_width=True, hide_index=True)


def _faq_section():
    _section_title("FAQ", "Common interpretation and usage questions.")
    for question, answer in FAQS:
        with st.expander(question, expanded=False):
            st.write(answer)


def _references_section():
    _section_title(
        "References & Citation",
        "Citation text, key algorithm references, and repository link.",
    )
    st.markdown("**Cite Non-B DNA Finder**")
    st.code(CITATION_TEXT, language="text")
    st.markdown("**Repository**")
    st.markdown(f"- [{GITHUB_URL}]({GITHUB_URL})")
    st.markdown("**Algorithm and background references**")
    for ref in REFERENCE_LIST:
        st.markdown(f"- {ref}")
    _callout(
        "Citation note",
        "No DOI is documented in this repository at present, so this documentation does not invent one.",
        "info",
    )


def render():
    load_css(TAB_THEMES.get("Documentation", "orchid_docs"))
    _render_header()
    tabs = st.tabs(
        [
            "Overview",
            "Motif Library",
            "Algorithms",
            "Scoring & Analysis",
            "Statistics",
            "User Guide",
            "FAQ",
            "References & Citation",
        ]
    )
    with tabs[0]:
        _overview_section()
    with tabs[1]:
        _motif_library_section()
    with tabs[2]:
        _algorithms_section()
    with tabs[3]:
        _scoring_section()
    with tabs[4]:
        _statistics_section()
    with tabs[5]:
        _user_guide_section()
    with tabs[6]:
        _faq_section()
    with tabs[7]:
        _references_section()
