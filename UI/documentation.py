"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Documentation Page - Scientific Documentation & References                    │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘
"""
# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
import os
import streamlit as st
import pandas as pd
from Utilities.config.text import UI_TEXT
from Utilities.config.themes import TAB_THEMES
from Utilities.config.colors import UNIFIED_MOTIF_COLORS, MOTIF_CLASS_INFO
from UI.css import load_css
from UI.headers import render_section_heading

# ═══════════════════════════════════════════════════════════════════════════════
# TUNABLE PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════════
MOTIF_PARAMETERS = {"Curved DNA": {"min_length": "8 bp (local) / ~31 bp (global)", "max_length": "50 bp (local) / 120 bp (global)", "algorithm": "A-tract phasing detection", "scoring": "Tract length × phasing score", "refs": ["Crothers et al., 1990", "Koo et al., 1986"]}, "G-Quadruplex": {"min_length": "15 bp", "max_length": "100 bp", "algorithm": "G4Hunter + Regex pattern matching", "scoring": "G4Hunter score (threshold ≥0.5)", "refs": ["Bedrat et al., 2016", "Huppert & Balasubramanian, 2005"]}, "Z-DNA": {"min_length": "8 bp", "max_length": "300 bp", "algorithm": "10-mer propensity table (Ho 1986) + eGZ regex", "scoring": "Cumulative 10-mer score (threshold ≥50.0); linear normalization", "refs": ["Ho et al., 1986", "Wang et al., 1979"]}, "Cruciform": {"min_length": "8 bp arm", "max_length": "50 bp arm", "algorithm": "Inverted repeat detection", "scoring": "Arm length × stem stability (ΔG ≤ −5.0 kcal/mol)", "refs": ["Lilley, 1980", "Panayotatos & Wells, 1981"]}, "R-Loop": {"min_length": "15 bp", "max_length": "2000 bp", "algorithm": "QmRLFS model", "scoring": "GC content + G-density (quality ≥0.4)", "refs": ["Jenjaroenpun et al., 2015", "Aguilera & García-Muse, 2012"]}, "Triplex": {"min_length": "10 bp", "max_length": "100 bp", "algorithm": "Mirror repeat detection", "scoring": "Purine/pyrimidine purity ≥90%", "refs": ["Mirkin et al., 1987", "Frank-Kamenetskii & Mirkin, 1995"]}, "i-Motif": {"min_length": "12 bp", "max_length": "60 bp", "algorithm": "C-tract pattern matching", "scoring": "C-density + C-tract count bonus", "refs": ["Zeraati et al., 2018", "Day et al., 2014"]}, "Slipped DNA": {"min_length": "1 bp unit", "max_length": "100 bp unit", "algorithm": "k-mer STR + direct repeat scan", "scoring": "Length × copy number × unit size × purity × GC", "refs": ["Pearson et al., 2005", "Wells, 2007"]}, "A-philic DNA": {"min_length": "10 bp", "max_length": "300 bp", "algorithm": "10-mer log₂ A-form propensity scoring", "scoring": "Cumulative log₂ sum (threshold ≥0.5)", "refs": ["Vinogradov, 2003", "Rohs et al., 2009"]}}
REFERENCES = [{"authors": "Bedrat A, Lacroix L, Mergny JL", "year": 2016, "title": "Re-evaluation of G-quadruplex propensity with G4Hunter", "journal": "Nucleic Acids Res", "volume": "44(4):1746-59", "doi": "10.1093/nar/gkw006"}, {"authors": "Huppert JL, Balasubramanian S", "year": 2005, "title": "Prevalence of quadruplexes in the human genome", "journal": "Nucleic Acids Res", "volume": "33(9):2908-16", "doi": "10.1093/nar/gki609"}, {"authors": "Zeraati M, Langley DB, et al.", "year": 2018, "title": "I-motif DNA structures are formed in the nuclei of human cells", "journal": "Nat Chem", "volume": "10(6):631-637", "doi": "10.1038/s41557-018-0046-3"}, {"authors": "Ho PS, Frederick CA, Saal D, et al.", "year": 1986, "title": "The interactions of ruthenium hexaammine with Z-DNA", "journal": "J Biomol Struct Dyn", "volume": "4(3):521-34", "doi": "10.1080/07391102.1986.10506363"}, {"authors": "Jenjaroenpun P, Wongsurawat T, et al.", "year": 2015, "title": "QmRLFS-finder: a model, web server and stand-alone tool", "journal": "Nucleic Acids Res", "volume": "43(W1):W527-34", "doi": "10.1093/nar/gkv344"}, {"authors": "Aguilera A, García-Muse T", "year": 2012, "title": "R loops: from transcription byproducts to threats to genome stability", "journal": "Mol Cell", "volume": "46(2):115-24", "doi": "10.1016/j.molcel.2012.04.009"}, {"authors": "Frank-Kamenetskii MD, Mirkin SM", "year": 1995, "title": "Triplex DNA structures", "journal": "Annu Rev Biochem", "volume": "64:65-95", "doi": "10.1146/annurev.bi.64.070195.000433"}, {"authors": "Crothers DM, Drak J, et al.", "year": 1990, "title": "DNA bending, flexibility, and helical repeat", "journal": "Methods Enzymol", "volume": "212:3-29", "doi": "10.1016/0076-6879(92)12003-9"}, {"authors": "Vinogradov AE", "year": 2003, "title": "DNA helix: the importance of being GC-rich", "journal": "Nucleic Acids Res", "volume": "31(7):1838-44", "doi": "10.1093/nar/gkg296"}, {"authors": "Rohs R, West SM, et al.", "year": 2009, "title": "The role of DNA shape in protein-DNA recognition", "journal": "Nature", "volume": "461(7268):1248-53", "doi": "10.1038/nature08473"}, {"authors": "Pearson CE, Nichol Edamura K, Cleary JD", "year": 2005, "title": "Repeat instability: mechanisms of dynamic mutations", "journal": "Nat Rev Genet", "volume": "6(10):729-42", "doi": "10.1038/nrg1689"}, {"authors": "Bacolla A, Wells RD", "year": 2004, "title": "Non-B DNA conformations, genomic rearrangements, and human disease", "journal": "J Biol Chem", "volume": "279(46):47411-4", "doi": "10.1074/jbc.R400028200"}]
MOTIF_DESCRIPTIONS = {'Curved_DNA': "Intrinsic DNA curvature from phased A-tracts", 'Slipped_DNA': "Slippage-mediated repeat expansions", 'Cruciform': "Hairpin structures from palindromic sequences", 'R-Loop': "Co-transcriptional R-loop formation sites", 'Triplex': "Triple-stranded DNA from mirror repeats", 'G-Quadruplex': "Four-stranded G-rich secondary structures", 'i-Motif': "Intercalated cytosine structures", 'Z-DNA': "Alternating purine-pyrimidine sequences", 'A-philic_DNA': "A-form DNA propensity regions", 'Hybrid': "Regions with multiple motif types", 'Non-B_DNA_Clusters': "High-density Non-B DNA regions"}
# ═══════════════════════════════════════════════════════════════════════════════

_PIPELINE_IMAGE = os.path.join(os.path.dirname(__file__), "pipelines.png")

_PROSE_STYLE = "font-family:Georgia,serif;line-height:1.8;color:#1e293b;font-size:1.0rem;"
_H3_STYLE = "color:#1e40af;font-size:1.2rem;margin-top:1.4rem;margin-bottom:0.4rem;"
_H4_STYLE = "color:#334155;font-size:1.05rem;margin-top:1.1rem;margin-bottom:0.3rem;"

# ─── Flow diagram builder ──────────────────────────────────────────────────────

_FLOW_COLORS = {
    'input':   ('#dbeafe', '#2563eb', '#1e40af'),
    'process': ('#f0f9ff', '#0284c7', '#0c4a6e'),
    'filter':  ('#fefce8', '#ca8a04', '#78350f'),
    'score':   ('#f0fdf4', '#16a34a', '#14532d'),
    'output':  ('#faf5ff', '#7c3aed', '#4c1d95'),
}


def _build_flow_diagram(stages, caption=""):
    """Build a horizontal scrollable pipeline flow diagram.

    stages: list of dicts with keys:
        'title'  – short label shown in the box
        'desc'   – one-line description (optional)
        'color'  – key in _FLOW_COLORS ('input', 'process', 'filter', 'score', 'output')
    caption: optional italic caption rendered below the diagram
    """
    items = []
    for i, s in enumerate(stages):
        bg, bdr, txt = _FLOW_COLORS.get(s.get('color', 'process'), _FLOW_COLORS['process'])
        desc_html = (
            f"<div style='color:#475569;font-size:0.78rem;margin-top:0.15rem;"
            f"line-height:1.3;'>{s['desc']}</div>"
            if s.get('desc') else ""
        )
        items.append(
            f"<div style='flex-shrink:0;background:{bg};border:2px solid {bdr};"
            f"border-radius:8px;padding:0.6rem 0.9rem;min-width:140px;max-width:190px;"
            f"text-align:center;'>"
            f"<div style='color:{bdr};font-size:0.72rem;font-weight:700;'>STEP {i + 1}</div>"
            f"<div style='color:{txt};font-size:0.85rem;font-weight:600;line-height:1.3;'>"
            f"{s['title']}</div>"
            f"{desc_html}</div>"
        )
        if i < len(stages) - 1:
            items.append(
                "<div style='flex-shrink:0;display:flex;align-items:center;"
                "color:#94a3b8;font-size:1.2rem;padding:0 0.3rem;'>&#8594;</div>"
            )
    caption_html = (
        f"<p style='text-align:center;color:#64748b;font-size:0.72rem;"
        f"margin-top:0.3rem;font-style:italic;'>{caption}</p>"
        if caption else ""
    )
    return (
        "<div style='overflow-x:auto;margin-bottom:1rem;'>"
        "<div style='display:inline-flex;align-items:stretch;gap:0;"
        "padding:0.8rem 1rem;background:#f8fafc;border-radius:10px;"
        "border:1px solid #e2e8f0;min-width:max-content;'>"
        + "".join(items)
        + "</div>"
        + caption_html
        + "</div>"
    )


def _build_motif_card(n, sub, col, desc):
    return (
        f"<div style='background:white;padding:0.9rem;border-radius:8px;"
        f"box-shadow:0 1px 4px rgba(0,0,0,0.07);border-left:4px solid {col};'>"
        f"<strong style='color:#1e293b;font-size:1.0rem;'>{n}</strong>"
        f"<div style='color:{col};font-size:0.85rem;font-weight:600;margin:0.15rem 0;'>{sub}</div>"
        f"<div style='color:#64748b;font-size:0.85rem;line-height:1.4;'>{desc}</div></div>"
    )


def _build_reference_card(r):
    return (
        f"<div style='padding:0.65rem 0;border-bottom:1px solid #e2e8f0;'>"
        f"<div style='font-weight:600;color:#1e293b;font-size:0.95rem;margin-bottom:0.1rem;'>"
        f"{r['authors']} ({r['year']})</div>"
        f"<div style='color:#334155;font-size:0.9rem;font-style:italic;margin-bottom:0.1rem;'>{r['title']}</div>"
        f"<div style='color:#64748b;font-size:0.85rem;'><strong>{r['journal']}</strong> {r['volume']} · "
        f"<a href=\"https://doi.org/{r['doi']}\" target=\"_blank\" style=\"color:#2563eb;\">DOI: {r['doi']}</a>"
        f"</div></div>"
    )


def _section_heading(text):
    st.markdown(
        f"<h3 style='{_H3_STYLE}'>{text}</h3>",
        unsafe_allow_html=True,
    )


def _prose(html):
    st.markdown(
        f"<div style='{_PROSE_STYLE}'>{html}</div>",
        unsafe_allow_html=True,
    )


# ─── Tab renderers ────────────────────────────────────────────────────────────

def _tab_overview():
    st.markdown(
        "<p style='color:#334155;font-size:1.05rem;line-height:1.7;margin-bottom:1rem;'>"
        "Comprehensive platform for <strong>genome-wide detection</strong> of Non-B DNA structures. "
        "Implements <strong>11 motif classes</strong> with <strong>24 subclasses</strong>, "
        "validated against G4Hunter, QmRLFS, and Z-Seeker.</p>",
        unsafe_allow_html=True,
    )
    st.markdown(
        "<div style='display:flex;gap:0.5rem;flex-wrap:wrap;margin-bottom:1.2rem;'>"
        "<span style='background:#dbeafe;color:#1d4ed8;padding:0.35rem 0.8rem;border-radius:16px;font-size:0.88rem;font-weight:600;'>24,674 bp/s</span>"
        "<span style='background:#dbeafe;color:#1d4ed8;padding:0.35rem 0.8rem;border-radius:16px;font-size:0.88rem;font-weight:600;'>200MB+ sequences</span>"
        "<span style='background:#dbeafe;color:#1d4ed8;padding:0.35rem 0.8rem;border-radius:16px;font-size:0.88rem;font-weight:600;'>25+ visualizations</span>"
        "<span style='background:#dbeafe;color:#1d4ed8;padding:0.35rem 0.8rem;border-radius:16px;font-size:0.88rem;font-weight:600;'>Nature-ready output</span>"
        "</div>",
        unsafe_allow_html=True,
    )
    _section_heading("Detected Non-B DNA Motif Classes")
    cards_html = '<div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(240px,1fr));gap:0.5rem;margin-bottom:1rem;">'
    for info in MOTIF_CLASS_INFO:
        cards_html += _build_motif_card(
            info['name'], info['subtitle'],
            UNIFIED_MOTIF_COLORS[info['key']],
            MOTIF_DESCRIPTIONS.get(info['key'], info['subtitle']),
        )
    cards_html += '</div>'
    st.markdown(cards_html, unsafe_allow_html=True)


def _tab_architecture():
    if os.path.isfile(_PIPELINE_IMAGE):
        st.image(_PIPELINE_IMAGE, use_container_width=True)
        st.markdown(
            "<p style='text-align:center;color:#64748b;font-size:0.88rem;margin-top:0.25rem;margin-bottom:1.2rem;'>"
            "<em>Figure 1. NonBDNAFinder analytical pipeline — from sequence ingestion through "
            "motif detection, scoring, hybrid annotation, and structured export.</em></p>",
            unsafe_allow_html=True,
        )
    _section_heading("2.1 System Architecture and Workflow")
    _prose(
        "<p>NBDFinder is a two-layer computational framework designed for genome-scale detection and "
        "integration of non-B DNA–forming sequence motifs. The backend consists of nine modular "
        "Python-based detector classes — covering Curved DNA, Slipped DNA, Cruciform DNA, Triplex DNA, "
        "R-loops, Z-DNA, G-Quadruplex, i-Motif, and A-philic DNA — each implementing biologically "
        "informed regular expressions, compositional heuristics, and propensity-based scoring models. "
        "The frontend is implemented using the Streamlit framework and supports sequence submission, "
        "parameter selection, execution, visualization, and export. Two integrative post-processing "
        "stages — hybrid-region annotation and density-based cluster identification — operate across "
        "the outputs of all nine primary detectors.</p>"
        "<p>Input sequences may be provided as FASTA files, direct nucleotide strings, or retrieved "
        "programmatically from NCBI using Biopython (Entrez and SeqIO modules). Upon ingestion, "
        "sequences are normalized to uppercase, header lines are removed, and non-ATGC characters are "
        "filtered; all downstream composition calculations are performed exclusively on the resulting "
        "canonical base set. GC percentage is accordingly computed as (G + C) / (A + T + G + C) × 100, "
        "matching the NCBI and Ensembl standard by excluding ambiguous bases from both the numerator "
        "and the denominator. All coordinates are reported using 1-based inclusive indexing.</p>"
        "<p>The analytical workflow consists of: (i) sequence preprocessing and validation; "
        "(ii) primary detection of canonical motifs; (iii) detection of relaxed or variant subclasses; "
        "(iv) class-specific scoring and normalization; (v) intra-class overlap resolution; "
        "(vi) inter-class hybrid detection; (vii) density-based hotspot identification; and "
        "(viii) structured output generation. For sequences shorter than 100,000 bp, analysis proceeds "
        "directly using all detectors in a single pass. Sequences between 100,000 bp and 5,000,000 bp "
        "are processed using a two-worker ProcessPoolExecutor with 50 kb chunks and 2 kb overlaps. "
        "Sequences exceeding 5,000,000 bp are handled by a disk-streaming strategy that maintains "
        "approximately constant RAM usage (~70 MB) regardless of genome size.</p>"
    )


def _tab_motif_library():
    _section_heading("2.2 Motif Library and Structural Classification")
    _prose(
        "<p>The motif library was curated through systematic literature review and organized into nine "
        "principal structural detector classes, each with biologically supported subclasses: Curved DNA "
        "(global and local curvature subtypes), Slipped DNA (short tandem repeats and direct repeats), "
        "Cruciform DNA (palindromic inverted repeats), Triplex DNA (H-DNA mirror repeats and Sticky DNA), "
        "R-loops (promoter-proximal RNA–DNA hybrid formation sites), Z-DNA (left-handed helix and "
        "extruded-G Z-DNA motifs), G-Quadruplex (eight subclasses ranging from telomeric repeats to "
        "G-triplexes), i-Motif (canonical cytosine-intercalated structures and AC-motifs), and A-philic "
        "DNA (A-form propensity regions). Two integrative classes — hybrid motifs (overlapping multi-class "
        "intervals) and non-B DNA clusters (high-density structural hotspots) — are derived from the "
        "primary detector outputs. Disease-associated repeat expansions were incorporated where "
        "experimentally validated structural consequences have been reported.</p>"
    )
    _section_heading("Detection Parameters & Algorithms")
    td = "padding:0.45rem 0.6rem;border-bottom:1px solid #e2e8f0;"
    rows_html = "".join(
        f"<tr>"
        f"<td style='{td}font-weight:600;color:#1e293b;'>{motif}</td>"
        f"<td style='{td}'>{p['min_length']}</td>"
        f"<td style='{td}'>{p['max_length']}</td>"
        f"<td style='{td}'>{p['algorithm']}</td>"
        f"<td style='{td}'>{p['scoring']}</td>"
        f"</tr>"
        for motif, p in MOTIF_PARAMETERS.items()
    )
    th = "text-align:left;padding:0.5rem 0.6rem;border-bottom:2px solid #cbd5e1;background:#f1f5f9;"
    table_html = (
        "<table style='width:100%;border-collapse:collapse;font-size:0.8rem;font-family:Georgia,serif;'>"
        "<thead><tr>"
        f"<th style='{th}'>Motif Class</th>"
        f"<th style='{th}'>Min Length</th>"
        f"<th style='{th}'>Max Length</th>"
        f"<th style='{th}'>Algorithm</th>"
        f"<th style='{th}'>Scoring Method</th>"
        "</tr></thead>"
        f"<tbody>{rows_html}</tbody></table>"
    )
    st.markdown(table_html, unsafe_allow_html=True)


def _tab_detection_algorithms():
    subtabs = st.tabs([
        "Curved DNA", "Slipped DNA", "Cruciform", "Triplex",
        "R-Loop", "Z-DNA", "G-Quadruplex", "i-Motif", "A-philic DNA",
    ])

    with subtabs[0]:
        _section_heading("2.3.1 Curved DNA")
        st.markdown(
            _build_flow_diagram([
                {'title': 'Input Sequence', 'desc': 'Canonical ATGC bases', 'color': 'input'},
                {'title': 'A/T-tract Enumeration', 'desc': 'Runs ≥3 bp, both strands', 'color': 'process'},
                {'title': 'Center Detection', 'desc': 'A-runs exceed T-runs by ≥3 nt', 'color': 'process'},
                {'title': 'Phasing Analysis', 'desc': 'Helical spacing 9.9–11.1 bp', 'color': 'filter'},
                {'title': 'Local Tract Detection', 'desc': 'Isolated A/T-tracts ≥8 bp', 'color': 'process'},
                {'title': 'Score Normalization', 'desc': '1.0–3.0 confidence scale', 'color': 'score'},
                {'title': 'Curved DNA Hits', 'desc': 'APRs + local curvature motifs', 'color': 'output'},
            ], caption="Figure: Curved DNA detection pipeline — two-stage A-tract phasing analysis."),
            unsafe_allow_html=True,
        )
        _prose(
            "<p>Curvature-prone regions were identified by a two-stage A/T-tract analysis. In the first "
            "stage, all contiguous runs of A and T bases meeting a minimum length of 3 bp are enumerated. "
            "Each candidate window is analyzed on both strands to determine the dominant A-tract center "
            "position using a sliding algorithm that tracks A-run length, T-run length, and AT/TA "
            "transition states. Only windows where the longest A-rich run exceeds the longest T-only run "
            "by at least 3 nt are retained as called tract centers, and a fractional center position is "
            "computed for each.</p>"
            "<p>In the second stage, called centers are grouped into A-phased repeats (APRs) by requiring "
            "at least three consecutive centers with helical-phasing spacing between 9.9 and 11.1 bp "
            "(reference spacing 11.0 bp, approximating the B-DNA helical repeat). A phasing score is "
            "computed as 1 − (mean deviation from ideal spacing / maximum allowed deviation), rewarding "
            "tighter helical periodicity. Separately, long isolated A-tracts or T-tracts of at least 8 bp "
            "are identified as local curvature motifs. Their scores are computed as L/(L + 7), saturating "
            "asymptotically toward 1.0 with increasing tract length. All raw scores are normalized linearly "
            "to a 1.0–3.0 confidence scale bounded by empirically defined floor (0.1) and ceiling (0.95) "
            "values.</p>"
        )

    with subtabs[1]:
        _section_heading("2.3.2 Slipped DNA and Tandem Repeats")
        st.markdown(
            _build_flow_diagram([
                {'title': 'Input Sequence', 'desc': 'Canonical ATGC bases', 'color': 'input'},
                {'title': 'STR Pattern Scan', 'desc': 'Units 1–6 nt, ≥4–6 copies', 'color': 'process'},
                {'title': 'Direct Repeat Scan', 'desc': 'Units 10–50 nt, ≥2 copies', 'color': 'process'},
                {'title': 'Low-complexity Filter', 'desc': 'Shannon entropy ≥0.5 bits', 'color': 'filter'},
                {'title': 'Piecewise Scoring', 'desc': 'Length × copy × purity × GC', 'color': 'score'},
                {'title': 'Score Normalization', 'desc': '1.0–3.0 confidence scale', 'color': 'score'},
                {'title': 'Slipped DNA Hits', 'desc': 'STRs + direct repeats', 'color': 'output'},
            ], caption="Figure: Slipped DNA detection pipeline — tandem repeat scan with entropy filtering."),
            unsafe_allow_html=True,
        )
        _prose(
            "<p>Direct repeats and short tandem repeats (STRs) were identified as tandemly arrayed sequence "
            "units using Python regular expressions compiled with a per-unit-size pattern cache. STRs were "
            "defined using core parameters as repeat units of 1–4 nt occurring at least 6 times "
            "consecutively, and under relaxed parameters as units of 1–9 nt occurring at least 4 times. "
            "Direct repeats were defined as repeat units of 10–50 nt occurring at least 2 times, with a "
            "minimum total tract length of 20 bp. The boundary between STR and direct repeat classification "
            "is set at a unit length threshold of 10 nt. Low-complexity artifacts were filtered using "
            "Shannon entropy computed over A, T, G, and C base frequencies, with a minimum threshold of "
            "0.5 bits. Scores are computed by a piecewise mechanistic model integrating JIT-compiled "
            "functions for tract length, copy number, unit size, base purity, and GC fraction, normalized "
            "linearly to the 1.0–3.0 confidence scale.</p>"
        )

    with subtabs[2]:
        _section_heading("2.3.3 Cruciform-Forming Inverted Repeats")
        st.markdown(
            _build_flow_diagram([
                {'title': 'Input Sequence', 'desc': 'Canonical ATGC bases', 'color': 'input'},
                {'title': '6-mer Seed Index', 'desc': 'All hexamer positions indexed', 'color': 'process'},
                {'title': 'Rev. Complement Search', 'desc': 'Locate matching arm partner', 'color': 'process'},
                {'title': 'Bidirectional Extension', 'desc': 'Arms 8–50 bp, spacer ≤12 bp', 'color': 'process'},
                {'title': 'ΔG Stability Filter', 'desc': 'ΔG ≤ −5.0 kcal/mol (NN model)', 'color': 'filter'},
                {'title': 'Score Normalization', 'desc': 'Arm × symmetry × GC → 1.0–3.0', 'color': 'score'},
                {'title': 'Cruciform Hits', 'desc': 'Inverted repeats with ΔG filter', 'color': 'output'},
            ], caption="Figure: Cruciform detection pipeline — seed-and-extend with thermodynamic validation."),
            unsafe_allow_html=True,
        )
        _prose(
            "<p>Cruciform motifs were detected via a seed-and-extend algorithm applied to inverted repeat "
            "pairs. An initial 6-mer seed index is constructed over the sequence; for each seed hit, the "
            "reverse complement arm is located and extended bidirectionally until a mismatch is encountered "
            "or arm length limits are reached. Arm lengths between 8 and 50 bp with spacer loops up to "
            "12 bp and zero mismatches are retained. Thermodynamic stability is assessed using a "
            "nearest-neighbor energy model (SantaLucia 1998) with parameters covering all 16 dinucleotide "
            "steps, and only inverted repeats with predicted ΔG ≤ −5.0 kcal/mol are reported. Scoring "
            "integrates arm length, symmetry, and GC content, and is normalized to the 1.0–3.0 confidence "
            "scale using linear interpolation between raw score floor 0.5 and ceiling 0.95.</p>"
        )

    with subtabs[3]:
        _section_heading("2.3.4 Triplex and Sticky DNA")
        st.markdown(
            _build_flow_diagram([
                {'title': 'Input Sequence', 'desc': 'Canonical ATGC bases', 'color': 'input'},
                {'title': 'Mirror Repeat Scan', 'desc': 'Homopurine/homopyrimidine arms', 'color': 'process'},
                {'title': 'Purity Filter', 'desc': 'Purine/pyrimidine fraction ≥90%', 'color': 'filter'},
                {'title': 'H-DNA Scoring', 'desc': 'Arm + loop + purity weights', 'color': 'score'},
                {'title': 'Sticky DNA Scan', 'desc': '(GAA)n/(TTC)n ≥20 copies', 'color': 'process'},
                {'title': 'Score Normalization', 'desc': '1.0–3.0 confidence scale', 'color': 'score'},
                {'title': 'Triplex Hits', 'desc': 'H-DNA + Sticky DNA candidates', 'color': 'output'},
            ], caption="Figure: Triplex/Sticky DNA detection pipeline — mirror repeat scan with purity filtering."),
            unsafe_allow_html=True,
        )
        _prose(
            "<p>Triplex-forming mirror repeats were detected as homopurine or homopyrimidine mirror "
            "sequences with arm lengths of 10–100 bp and spacer tolerance up to 8 bp, requiring a purine "
            "or pyrimidine fraction of at least 0.90. A JIT-compiled purity calculation determines whether "
            "each candidate arm meets the strand-purity threshold. The H-DNA scoring model integrates arm "
            "length relative to a 35 bp reference, a loop penalty term modulated by a 0.4 exponent on "
            "spacer length, a purity weight (0.30), an arm-length weight (0.35), a loop weight (0.20), and "
            "an interruption penalty weight (0.15). Raw scores are normalized linearly to the 1.0–3.0 "
            "confidence scale using floor 0.5 and ceiling 0.95.</p>"
            "<p>Sticky DNA was operationalized as uninterrupted (GAA)n or (TTC)n tracts using a piecewise "
            "linear scoring model calibrated to three biological thresholds: a replication-blockage "
            "threshold (≥20 copies, base score 1.3), a stable-formation threshold (≥40 copies, base score "
            "2.0), and a pathogenic-expansion threshold (≥60 copies, base score 2.6). These cutoffs are "
            "consistent with clinical Friedreich's ataxia expansion thresholds reported in the "
            "literature.</p>"
        )

    with subtabs[4]:
        _section_heading("2.3.5 R-Loop–Forming Sequences")
        st.markdown(
            _build_flow_diagram([
                {'title': 'Input Sequence', 'desc': 'Canonical ATGC bases', 'color': 'input'},
                {'title': 'RIZ Detection', 'desc': 'G-cluster models, G-content ≥50%', 'color': 'process'},
                {'title': 'Model 1: G{3+} Clusters', 'desc': 'Separated by ≤10 intervening bases', 'color': 'process'},
                {'title': 'Model 2: G{4+} Tracts', 'desc': 'Dense guanine tracts', 'color': 'process'},
                {'title': 'REZ Extension', 'desc': 'GC-rich zone ≤2000 bp, GC ≥40%', 'color': 'filter'},
                {'title': 'Score Normalization', 'desc': 'GC × G-density → 1.0–3.0', 'color': 'score'},
                {'title': 'R-Loop Hits', 'desc': 'RIZ + REZ combined intervals', 'color': 'output'},
            ], caption="Figure: R-Loop detection pipeline — two-stage QmRLFS-adapted initiation/elongation model."),
            unsafe_allow_html=True,
        )
        _prose(
            "<p>R-loop prediction followed a two-stage framework adapted from the QmRLFS model "
            "(Jenjaroenpun et al. 2015, 2016). In the first stage, initiation zones (R-loop initiation "
            "zones, RIZ) are identified as regions where G-run–enriched patterns match either of two "
            "models: Model 1 captures G{3+} clusters separated by up to 10 intervening bases, while "
            "Model 2 captures denser G{4+} tracts. A guanine content threshold of ≥50% is required for "
            "RIZ candidacy. In the second stage, each valid RIZ is extended downstream into a GC-rich "
            "elongation zone (REZ) of up to 2,000 bp, requiring a GC content of at least 40% across the "
            "combined region. Final R-loop scores integrate GC content and G-run density across both zones "
            "and are normalized linearly to the 1.0–3.0 confidence scale.</p>"
        )

    with subtabs[5]:
        _section_heading("2.3.6 Z-DNA and eGZ Motifs")
        st.markdown(
            _build_flow_diagram([
                {'title': 'Input Sequence', 'desc': 'Canonical ATGC bases', 'color': 'input'},
                {'title': '10-mer Table Matching', 'desc': 'Hyperscan / Python fallback', 'color': 'process'},
                {'title': 'Per-base Score Accumulation', 'desc': 'Uniform contribution over 10 bp', 'color': 'process'},
                {'title': 'Region Merging', 'desc': 'Adjacent/overlapping 10-mers merged', 'color': 'process'},
                {'title': 'Threshold Filter', 'desc': 'Cumulative score ≥50.0', 'color': 'filter'},
                {'title': 'eGZ Detection', 'desc': 'CGG/GGC/CCG/GCC repeats ≥4 copies', 'color': 'process'},
                {'title': 'Z-DNA Hits', 'desc': 'Z-DNA + eGZ candidates', 'color': 'output'},
            ], caption="Figure: Z-DNA/eGZ detection pipeline — 10-mer propensity table with linear normalization."),
            unsafe_allow_html=True,
        )
        _prose(
            "<p>Z-DNA propensity was estimated using a precomputed 10-mer scoring table derived from the "
            "dinucleotide propensity matrix of Ho et al. (1986), which assigns higher scores to CG, GC, "
            "AC, and GT dinucleotide steps. All exact 10-mer matches in the query sequence are located "
            "using either the Hyperscan high-performance regex engine (when available) or a pure-Python "
            "sliding-window fallback. Matched 10-mers are merged if they are adjacent or overlapping, and "
            "per-base contributions are aggregated by distributing each 10-mer score uniformly across its "
            "ten positions, optionally using NumPy vectorization for sequences longer than 1,000 bp. Merged "
            "regions must exceed a cumulative sum score threshold (default 50.0) and contain at least one "
            "contributing 10-mer to be reported. Scores are normalized linearly between raw "
            "floor 50.0 and ceiling 2,000.0.</p>"
            "<p>Extruded-G Z-DNA (eGZ) motifs were independently detected as trinucleotide repeats of CGG, "
            "GGC, CCG, or GCC with at least 4 consecutive copies, consistent with the structural model of "
            "Herbert (1997). Each eGZ repeat is scored proportionally to repeat count relative to the "
            "minimum threshold and annotated for association with disease-relevant CGG/CCG expansions "
            "(Fragile X syndrome) when copy numbers exceed established clinical cutoffs.</p>"
        )

    with subtabs[6]:
        _section_heading("2.3.7 G-Quadruplex Family")
        st.markdown(
            _build_flow_diagram([
                {'title': 'Input Sequence', 'desc': 'Canonical ATGC bases', 'color': 'input'},
                {'title': 'G-seed Location', 'desc': 'All G{3+} positions indexed', 'color': 'process'},
                {'title': 'Scan Region Merging', 'desc': 'Contiguous G-rich windows', 'color': 'process'},
                {'title': '8-Subclass Matching', 'desc': 'Canonical, telomeric, bulged, etc.', 'color': 'process'},
                {'title': 'G4Hunter Scoring', 'desc': 'Sliding window 25 nt, JIT-compiled', 'color': 'score'},
                {'title': 'Overlap Resolution', 'desc': 'Priority: telomeric > higher-order…', 'color': 'filter'},
                {'title': 'G4 Hits', 'desc': '8 subclass G-Quadruplex candidates', 'color': 'output'},
            ], caption="Figure: G-Quadruplex detection pipeline — seeded scan with 8-subclass pattern matching and G4Hunter scoring."),
            unsafe_allow_html=True,
        )
        _prose(
            "<p>Canonical G4 motifs were identified using the pattern G{3+}N{1–7}G{3+}N{1–7}G{3+}N{1–7}"
            "G{3+}. Additional detectors captured telomeric repeats ((TTAGGG){4+}), extended-loop G4s "
            "(loops 1–12 nt), bulged G4s (G-tracts containing a single-nucleotide internal bulge), "
            "higher-order G4 arrays and G-wires (seven or more consecutive G-tract repeats), stacked "
            "multi-quadruplex assemblies (two or more canonical G4 units separated by up to 20 bp linkers), "
            "G-triplex intermediates (three G-runs), and weak PQS motifs (four G-runs of at least 2 G "
            "each).</p>"
            "<p>Detection proceeds via a seeded scanning strategy: all G{3+} positions are first located, "
            "adjacent seed windows are merged into contiguous scan regions to eliminate redundant processing "
            "in GC-rich stretches, and each pattern is matched only within those regions. Each candidate is "
            "scored using a G-only sliding-window G4Hunter algorithm (window size 25 nt by default), "
            "computing the maximum normalized G density and scaling it by the region length relative to the "
            "window size. Numba JIT compilation is applied to the sliding-window maximum subarray function "
            "for sequences exceeding 50 bp. Overlap resolution retains the highest-scoring representative "
            "from each genomic interval, with priority assigned in the order: telomeric > higher-order > "
            "stacked > canonical > bulged > extended-loop > G-triplex > weak PQS.</p>"
        )

    with subtabs[7]:
        _section_heading("2.3.8 i-Motif and AC-Motif Detection")
        st.markdown(
            _build_flow_diagram([
                {'title': 'Input Sequence', 'desc': 'Canonical ATGC bases', 'color': 'input'},
                {'title': 'C-tract Detection', 'desc': 'C{3+}N{1–7} × 4, loops 1–7 nt', 'color': 'process'},
                {'title': 'AC-motif Scan', 'desc': 'Phased A/C tracts (Hur et al. 2021)', 'color': 'process'},
                {'title': 'Subclass Assignment', 'desc': 'Short-loop vs long-loop variants', 'color': 'process'},
                {'title': 'Overlap Resolution', 'desc': 'Canonical i-motif > AC-motif', 'color': 'filter'},
                {'title': 'Score Normalization', 'desc': 'C-density × cytosine% × loop → 1.0–3.0', 'color': 'score'},
                {'title': 'i-Motif Hits', 'desc': 'Canonical + AC-motif candidates', 'color': 'output'},
            ], caption="Figure: i-Motif detection pipeline — C-tract and AC-motif scan with subclass assignment."),
            unsafe_allow_html=True,
        )
        _prose(
            "<p>i-motifs were defined as sequences containing four or more cytosine tracts, each consisting "
            "of at least 3 consecutive cytosines, separated by 1–7 nt loops, matching the pattern "
            "C{3+}N{1–7}C{3+}N{1–7}C{3+}N{1–7}C{3+}. Scores integrate C-run density, cytosine fraction, "
            "and loop compactness. Canonical subclasses with short loops (1–7 nt) are distinguished from "
            "long-loop variants during subclass assignment. AC-motifs were identified via phased alternating "
            "A-tract and C-tract consensus patterns as described by Hur et al. (2021), covering spacer "
            "lengths of 4, 5, and 6 bp between tract anchors, and reported qualitatively with separate "
            "pattern identifiers. Overlap resolution retains canonical i-motifs over AC-motifs when both "
            "patterns share the same genomic interval.</p>"
        )

    with subtabs[8]:
        _section_heading("2.3.9 A-Philic DNA Propensity Modeling")
        st.markdown(
            _build_flow_diagram([
                {'title': 'Input Sequence', 'desc': 'Canonical ATGC bases', 'color': 'input'},
                {'title': '10-mer Table Matching', 'desc': 'A-form log₂ odds propensity table', 'color': 'process'},
                {'title': 'Per-base Score Accumulation', 'desc': 'NumPy vectorization ≥1000 bp', 'color': 'process'},
                {'title': 'Segment Merging', 'desc': 'Adjacent high-scoring regions joined', 'color': 'process'},
                {'title': 'Threshold Filter', 'desc': 'Cumulative log₂ sum ≥0.5', 'color': 'filter'},
                {'title': 'Score Normalization', 'desc': '1.0–3.0 confidence scale', 'color': 'score'},
                {'title': 'A-philic DNA Hits', 'desc': 'A-form propensity regions', 'color': 'output'},
            ], caption="Figure: A-philic DNA detection pipeline — 10-mer A-form propensity scoring with segment merging."),
            unsafe_allow_html=True,
        )
        _prose(
            "<p>A-philic DNA propensity was estimated using a precomputed 10-mer scoring table of log₂ "
            "enrichment values derived from A-form and B-form DNA crystal structures in the Nucleic Acid "
            "Knowledgebase. Only canonical ATGC sequences were retained; tetranucleotide step frequencies "
            "were stratified by conformational class, and signed propensity scores were computed as log₂ "
            "odds ratios with pseudocount 0.5. The scoring table covers all possible 10-mer sequences with "
            "positive A-form propensity. Detection proceeds by locating all exact 10-mer matches using the "
            "same Hyperscan or pure-Python fallback used for Z-DNA, accumulating per-base contributions "
            "from each matching 10-mer into a per-position score array (optionally vectorized with NumPy), "
            "and merging adjacent high-scoring segments into contiguous A-philic regions. A minimum "
            "cumulative log₂ sum threshold of 0.5 is applied. Raw scores are normalized linearly to the "
            "1.0–3.0 confidence scale.</p>"
        )


def _tab_scoring():
    _section_heading("Universal Confidence Scale")
    _prose(
        "<p>All nine primary detector classes report scores normalized to a universal 1.0–3.0 confidence "
        "scale using class-specific raw score bounds, enabling direct cross-class comparison. "
        "Normalization uses linear interpolation between empirically defined floor and ceiling values for "
        "each class. Scores below the floor are clipped to 1.0; scores above the ceiling are clipped to "
        "3.0. A length-aware component scales the score proportionally to motif length relative to the "
        "structural cap for each class.</p>"
    )
    _section_heading("Per-Class Scoring Methods")
    rows = []
    for motif, p in MOTIF_PARAMETERS.items():
        rows.append(
            f"<tr style='border-bottom:1px solid #e2e8f0;'>"
            f"<td style='padding:0.45rem 0.6rem;font-weight:600;color:#1e293b;'>{motif}</td>"
            f"<td style='padding:0.45rem 0.6rem;color:#334155;'>{p['scoring']}</td>"
            f"</tr>"
        )
    st.markdown(
        "<table style='width:100%;border-collapse:collapse;font-size:0.82rem;font-family:Georgia,serif;'>"
        "<thead><tr style='background:#f1f5f9;'>"
        "<th style='text-align:left;padding:0.5rem 0.6rem;border-bottom:2px solid #cbd5e1;'>Motif Class</th>"
        "<th style='text-align:left;padding:0.5rem 0.6rem;border-bottom:2px solid #cbd5e1;'>Scoring Method</th>"
        "</tr></thead><tbody>" + "".join(rows) + "</tbody></table>",
        unsafe_allow_html=True,
    )


def _tab_hybrid_clustering():
    _section_heading("2.4 Hybrid Motifs and Cluster Detection")
    _prose(
        "<p>Hybrid regions were defined as genomic intervals where two or more distinct motif classes "
        "overlapped. Two motifs from different classes are considered overlapping when their shared bases "
        "constitute at least 50% of the length of the shorter interval (HYBRID_MIN_OVERLAP = 0.50). "
        "Near-identical overlaps exceeding 99% are excluded (HYBRID_MAX_OVERLAP = 0.99) to avoid "
        "reporting nearly identical motifs as hybrids. Overlapping pairs are merged into a single hybrid "
        "annotation that records all contributing class labels and is assigned a composite score computed "
        "as the mean score of the two contributing motifs.</p>"
        "<p>Structural hotspots (clusters) were identified by scanning the sequence using a 300 bp window "
        "anchored at each detected motif start position and counting overlapping motif calls from any "
        "class within each window. Windows containing at least 4 total motifs from at least 3 unique "
        "structural classes were recorded as cluster candidates. Each reported cluster is annotated "
        "with its genomic span, total motif count, class diversity, and a score reflecting the average "
        "score of the contributing motifs.</p>"
    )


def _tab_optimization():
    _section_heading("2.5 Algorithmic Optimization and Statistical Controls")
    _prose(
        "<p>Motif scanning for Z-DNA and A-philic DNA was accelerated using Hyperscan (Intel's "
        "high-performance regular expression engine), which provides linear-time exact-string matching "
        "over large 10-mer pattern sets via compiled automata. When Hyperscan is unavailable, a "
        "pure-Python sliding-window fallback is invoked transparently. For G-Quadruplex, Triplex, and "
        "Slipped DNA detection, computationally intensive scoring kernels (sliding-window maximum subarray, "
        "purity calculation, slippage base score) are compiled using Numba just-in-time (JIT) compilation "
        "with caching, providing 2–5× speedup over interpreted Python on the hot path. NumPy vectorization "
        "is applied to per-base contribution array construction for Z-DNA and A-philic DNA for sequences "
        "longer than 1,000 bp.</p>"
        "<p>For multi-detector execution, sequences exceeding 50,000 bp are analyzed using a "
        "ThreadPoolExecutor with up to nine concurrent worker threads (one per detector type), yielding a "
        "1.5–2× throughput improvement on multi-core hardware. All motif calls are validated for minimum "
        "length, length normalization, and GC composition adjustment. Shannon entropy (minimum 0.5 bits "
        "for Slipped DNA) is applied to filter low-complexity artifacts. Scores are normalized to a "
        "universal 1.0–3.0 confidence scale using class-specific raw score bounds, enabling cross-class "
        "comparison.</p>"
    )


def _tab_validation():
    _section_heading("2.6 Output Structure and Implementation")
    _prose(
        "<p>For each sequence, NBDFinder reports a tabular record for every detected motif containing: "
        "motif coordinates (Start, End, Length), subclass label, class-specific score, pattern identifier, "
        "strand orientation, detection method, and class-specific supplementary fields such as G-tract "
        "counts and loop lengths for G-Quadruplex, arm length and thermodynamic stability for Cruciform, "
        "repeat unit and copy count for Slipped DNA, and GC skew metrics for R-loops. Per-class counts "
        "and motif coverage fractions are summarized in a companion statistics table. Outputs are available "
        "in CSV and XLSX formats and include summary statistics and visualization-ready annotations "
        "compatible with the integrated Streamlit visualization pipeline (Manhattan plots, density "
        "heatmaps, circos-style diagrams).</p>"
        "<p>The platform is implemented in Python 3 using pandas, NumPy, matplotlib, seaborn, and the "
        "regex library. Optional acceleration dependencies include Numba (≥0.56.0) for JIT compilation, "
        "Hyperscan for high-performance pattern matching, and PyFastx for accelerated FASTA parsing. "
        "The application is open-source and available at "
        "<a href='https://github.com/VRYella/NonBDNAFinder' target='_blank' style='color:#1e40af;'>"
        "https://github.com/VRYella/NonBDNAFinder</a>.</p>"
    )


def _tab_references():
    _prose(
        "<p style='margin-bottom:0.8rem;'>NonBDNAFinder implements algorithms validated in peer-reviewed "
        "publications:</p>"
    )
    mid = len(REFERENCES) // 2
    col_a, col_b = st.columns(2)
    with col_a:
        st.markdown(
            '<div style="display:flex;flex-direction:column;">'
            + "".join(_build_reference_card(r) for r in REFERENCES[:mid])
            + "</div>",
            unsafe_allow_html=True,
        )
    with col_b:
        st.markdown(
            '<div style="display:flex;flex-direction:column;">'
            + "".join(_build_reference_card(r) for r in REFERENCES[mid:])
            + "</div>",
            unsafe_allow_html=True,
        )


def _tab_citation():
    _section_heading("How to Cite")
    st.markdown(
        "<div style='font-family:\"Courier New\",monospace;font-size:0.78rem;line-height:1.5;"
        "color:#334155;background:#f8fafc;padding:0.8rem 1rem;border-left:4px solid #2563eb;"
        "border-radius:4px;margin-bottom:1rem;'>"
        "<strong>Yella VR</strong> (2025). NonBDNAFinder: Comprehensive Detection and Analysis of "
        "Non-B DNA Forming Motifs.<br>"
        "GitHub: <a href='https://github.com/VRYella/NonBDNAFinder' style='color:#2563eb;'>"
        "https://github.com/VRYella/NonBDNAFinder</a></div>",
        unsafe_allow_html=True,
    )
    _section_heading("Developed by")
    st.markdown(
        f"<p style='font-size:0.85rem;color:#334155;line-height:1.7;'>"
        f"<strong>{UI_TEXT['author']}</strong><br>"
        f"Email: <a href='mailto:{UI_TEXT['author_email']}' style='color:#2563eb;'>"
        f"{UI_TEXT['author_email']}</a><br>"
        f"<a href='https://github.com/VRYella' target='_blank' style='color:#2563eb;'>"
        f"GitHub: VRYella</a></p>",
        unsafe_allow_html=True,
    )


# ─── Public entry point ───────────────────────────────────────────────────────

def render():
    load_css(TAB_THEMES.get('Documentation', 'scientific_blue'))
    render_section_heading("Scientific Documentation & References", page="Documentation")

    tabs = st.tabs([
        "Overview",
        "Architecture",
        "Motif Library",
        "Detection Algorithms",
        "Scoring & Normalization",
        "Hybrid & Clustering",
        "Optimization",
        "Validation",
        "References",
        "Citation",
    ])

    with tabs[0]:
        _tab_overview()
    with tabs[1]:
        _tab_architecture()
    with tabs[2]:
        _tab_motif_library()
    with tabs[3]:
        _tab_detection_algorithms()
    with tabs[4]:
        _tab_scoring()
    with tabs[5]:
        _tab_hybrid_clustering()
    with tabs[6]:
        _tab_optimization()
    with tabs[7]:
        _tab_validation()
    with tabs[8]:
        _tab_references()
    with tabs[9]:
        _tab_citation()
