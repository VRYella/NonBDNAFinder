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
            f"border-radius:8px;padding:0.3rem 0.9rem;min-width:140px;max-width:190px;"
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
        "padding:0.4rem 1rem;background:#f8fafc;border-radius:10px;"
        "border:1px solid #e2e8f0;min-width:max-content;'>"
        + "".join(items)
        + "</div>"
        + caption_html
        + "</div>"
    )


def _build_motif_card(n, sub, col, desc):
    return (
        f"<div style='background:white;padding:0.25rem 0.7rem;border-radius:8px;"
        f"box-shadow:0 1px 4px rgba(0,0,0,0.07);border-left:4px solid {col};'>"
        f"<strong style='color:#1e293b;font-size:0.95rem;'>{n}</strong>"
        f"<div style='color:{col};font-size:0.85rem;font-weight:600;margin:0.1rem 0;'>{sub}</div>"
        f"<div style='color:#64748b;font-size:0.85rem;line-height:1.3;'>{desc}</div></div>"
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


def _tab_statistical_guide():
    _section_heading("How to Interpret Your Results")
    _prose(
        "<p>This guide provides precise definitions, mathematical formulas, valid value ranges, and "
        "biological significance for every statistical metric reported by NonBDNAFinder. Understanding "
        "these metrics allows you to assess the structural complexity of your sequence, compare results "
        "across genomes, and identify biologically meaningful patterns.</p>"
    )

    # ── 1. Core output metrics ────────────────────────────────────────────────
    _section_heading("1. Core Output Metrics")
    th = "text-align:left;padding:0.5rem 0.75rem;border-bottom:2px solid #cbd5e1;background:#f1f5f9;font-size:0.8rem;color:#334155;"
    td = "padding:0.45rem 0.75rem;border-bottom:1px solid #e2e8f0;font-size:0.8rem;vertical-align:top;"
    def _metric_table(rows):
        header = (
            "<table style='width:100%;border-collapse:collapse;font-family:Georgia,serif;margin-bottom:1.2rem;'>"
            "<thead><tr>"
            f"<th style='{th}'>Metric</th>"
            f"<th style='{th}'>Formula</th>"
            f"<th style='{th}'>Range</th>"
            f"<th style='{th}'>Significance</th>"
            "</tr></thead><tbody>"
        )
        body = "".join(
            f"<tr>"
            f"<td style='{td}font-weight:600;color:#1e40af;'>{m}</td>"
            f"<td style='{td}font-family:\"Courier New\",monospace;color:#0f766e;'>{f}</td>"
            f"<td style='{td}color:#7c3aed;'>{r}</td>"
            f"<td style='{td}color:#334155;'>{s}</td>"
            "</tr>"
            for m, f, r, s in rows
        )
        return header + body + "</tbody></table>"

    st.markdown(_metric_table([
        (
            "Total Motifs (N)",
            "count of all detected motif intervals",
            "0 – ∞ (genome-size dependent)",
            "Absolute number of non-B DNA–forming sites in the analysed region. "
            "Higher counts suggest a structurally complex or repeat-rich genome.",
        ),
        (
            "Motif Density (D)",
            "D = N / (G / 1 000)  [motifs per kb]",
            "Typical prokaryotes: 1–10 /kb; "
            "eukaryotes: 5–40 /kb",
            "Normalizes motif abundance to sequence length, enabling fair comparison across "
            "genomes of different sizes. Values >20/kb indicate a structurally dense sequence.",
        ),
        (
            "Coverage %",
            "Cov% = (∑ unique covered bases / G) × 100",
            "0 – 100 %",
            "Fraction of the genome occupied by at least one non-B motif (after merging "
            "overlapping intervals). Values >10% flag highly non-B–enriched sequences.",
        ),
        (
            "Total Covered Bases",
            "∑ union of all motif intervals (bp)",
            "0 – G (genome length)",
            "Absolute base-pair footprint of non-B DNA. Useful for estimating the proportion "
            "of the sequence that may adopt alternative secondary structures in vivo.",
        ),
        (
            "Raw Occupancy (bp)",
            "∑ Li  (sum of all motif lengths, counting overlaps)",
            "0 – ∞",
            "Total base-pair length of all motifs including overlap regions. Raw occupancy "
            "exceeds covered bases whenever motifs overlap.",
        ),
        (
            "Normalized Occupancy",
            "NO = (∑ Li) / G",
            "0 – ∞ (values >1 indicate average overlap depth >1×)",
            "Mean motif 'pile-up' depth per base. NO = 0.5 means, on average, every other "
            "base is covered by a motif.",
        ),
        (
            "Mean Overlap Depth",
            "Ratio of Raw Occupancy to Covered Bases",
            "1.0 – ∞",
            "Average number of motifs stacked at each covered base. Values near 1.0 indicate "
            "non-overlapping motifs; high values (>3) suggest hotspot regions.",
        ),
        (
            "GC Content",
            "GC% = (G + C) / (A + T + G + C) × 100",
            "0 – 100 %; balanced genomes: 30–70 %",
            "Only canonical ATGC bases are counted in both numerator and denominator. "
            "Deviations outside 30–70% correlate with higher densities of Z-DNA, G4, "
            "and i-Motif structures.",
        ),
    ]), unsafe_allow_html=True)

    # ── 2. Confidence Score (1.0–3.0 scale) ──────────────────────────────────
    _section_heading("2. Confidence Score (Universal 1.0–3.0 Scale)")
    _prose(
        "<p>Every reported motif carries a normalized confidence score on a universal "
        "<strong>1.0–3.0</strong> scale, enabling direct comparison across all nine structural classes. "
        "Normalization uses linear interpolation between empirically defined class-specific floor and "
        "ceiling values:</p>"
        "<p style='font-family:\"Courier New\",monospace;background:#f8fafc;padding:0.5rem 1rem;"
        "border-left:4px solid #2563eb;border-radius:4px;'>"
        "Score = 1.0 + 2.0 × (raw − floor) / (ceiling − floor)"
        "</p>"
        "<p>Scores below the floor are clipped to 1.0; scores above the ceiling are clipped to 3.0.</p>"
    )
    st.markdown(_metric_table([
        ("Score 1.0 – 1.5", "Low confidence", "Any value ≥1.0",
         "Weak pattern match; sequence satisfies minimum motif criteria but may be a "
         "low-probability structural form. Treat with caution in downstream analyses."),
        ("Score 1.5 – 2.0", "Moderate confidence", "—",
         "Canonical motif with typical structural parameters. Represents the majority of "
         "detected sites in a typical genomic sequence."),
        ("Score 2.0 – 2.5", "High confidence", "—",
         "Strong structural candidate with favourable length, composition, and thermodynamic "
         "parameters. Prioritise these sites for experimental validation."),
        ("Score 2.5 – 3.0", "Very high confidence", "Maximum = 3.0",
         "Optimal or extended motif with multiple reinforcing features. Highest-priority sites "
         "for ChIP-seq, DMS-MaPseq, or biophysical studies."),
    ]), unsafe_allow_html=True)

    # ── 3. Structural Load Metrics ────────────────────────────────────────────
    _section_heading("3. Structural Load Metrics")
    _prose(
        "<p>Structural load metrics quantify the overall burden of non-B DNA structures relative "
        "to genome size, weighting contributions by both motif length and motif quality.</p>"
    )
    st.markdown(_metric_table([
        (
            "Structural Load Index (SLI)",
            "SLI = ∑ Li / G",
            "0 – ∞ (typical: 0.01 – 2.0)",
            "Identical to Normalized Occupancy. SLI = 0.1 means 10% of the genome is covered "
            "by non-B DNA motif bases on average. Values >0.5 indicate a heavily loaded sequence.",
        ),
        (
            "Structural Intensity (SI)",
            "SI = ∑(Score_i × Li) / G",
            "0 – 3.0 (bounded by score ceiling × SLI)",
            "Score-weighted coverage: rewards long, high-confidence motifs. "
            "More informative than SLI alone for assessing biological impact.",
        ),
        (
            "Weighted Structural Coverage (WSC)",
            "WSC = ∑(Score_i / Score_max × Li) / G",
            "0 – 1.0",
            "Normalizes SI by the maximum observed score, producing a dimensionless "
            "quality-weighted coverage fraction. WSC = 1.0 would mean the entire genome is "
            "covered by maximum-confidence motifs.",
        ),
    ]), unsafe_allow_html=True)

    # ── 4. Diversity Metrics ──────────────────────────────────────────────────
    _section_heading("4. Structural Diversity Metrics")
    _prose(
        "<p>Diversity metrics describe how motif counts are distributed across structural classes, "
        "analogous to species diversity measures in ecology.</p>"
    )
    st.markdown(_metric_table([
        (
            "Simpson Diversity Index (D)",
            "D = 1 − ∑(n_c / N)²",
            "0 – 1  (higher = more diverse)",
            "D = 0: all motifs belong to one class (no diversity). "
            "D → 1: motifs are evenly distributed across many classes. "
            "D > 0.7 indicates a highly diverse structural landscape. "
            "Directly analogous to the ecological Simpson index.",
        ),
        (
            "Effective Class Number (N_eff)",
            "N_eff = 1 / ∑(n_c / N)²",
            "1 – C  (where C = number of classes detected)",
            "Equivalent number of equally-abundant classes that would produce the observed "
            "diversity. N_eff = 3.5 means the distribution is as diverse as if 3.5 classes "
            "were equally represented.",
        ),
        (
            "Classes Detected (C)",
            "Count of distinct motif classes with ≥ 1 hit",
            "0 – 11  (0 = no classes detected; 11 = all 9 primary + Hybrid + Clusters)",
            "Raw class richness. C = 9 indicates all primary structural classes are present, "
            "suggesting a genomically complex sequence.",
        ),
        (
            "Subclasses Detected",
            "Count of distinct subclass labels with ≥ 1 hit",
            "0 – 24  (0 = none detected; 24 = all canonical subclasses present)",
            "Subclass richness across all 24 canonical subclasses. Higher values indicate "
            "structural variety even within individual motif families.",
        ),
    ]), unsafe_allow_html=True)

    # ── 5. Comparative and Complexity Metrics ─────────────────────────────────
    _section_heading("5. Comparative & Complexity Metrics")
    st.markdown(_metric_table([
        (
            "Structural Complexity Index (SCI)",
            "SCI = Coverage_Fraction × N_eff",
            "0 – C (typical: 0 – 5)",
            "Composite metric integrating coverage breadth with structural diversity. "
            "SCI = 2.0 means the sequence covers 50% of the genome with an effective "
            "class richness of 4, or 100% coverage with an effective class richness of 2.",
        ),
        (
            "Dominance Ratio (DR)",
            "DR = max(n_c) / N",
            "0 – 1  (lower = more even distribution)",
            "Fraction of all motifs attributable to the single most abundant class. "
            "DR > 0.8 indicates one structural class dominates. Low DR (< 0.3) suggests "
            "a balanced multi-class distribution.",
        ),
        (
            "Class Contribution (%)",
            "(n_c / N) × 100 per class",
            "0 – 100 % per class; all classes sum to 100 %",
            "Relative proportion of each structural class in the total motif pool, useful "
            "for identifying the dominant structural theme of a sequence.",
        ),
        (
            "Class Coverage (%)",
            "(unique bases covered by class c / G) × 100",
            "0 – 100 % per class",
            "Individual genomic footprint of each class independent of other classes. "
            "Useful for assessing which structural type has the greatest spatial spread.",
        ),
    ]), unsafe_allow_html=True)

    # ── 6. Spatial distribution metrics ──────────────────────────────────────
    _section_heading("6. Spatial Distribution Metrics")
    st.markdown(_metric_table([
        (
            "Mean Inter-Motif Distance",
            "mean(Start_{i+1} − End_i)  across consecutive motifs",
            "0 – G bp  (0 = adjacent motifs)",
            "Average gap (in bp) between consecutive motifs. Short mean distances indicate "
            "clustering; large distances suggest sparse, widely-spaced structural sites.",
        ),
        (
            "CV Spatial Clustering",
            "CV = σ(inter-motif gaps) / μ(inter-motif gaps)",
            "0 – ∞  (CV = 0: perfectly regular; CV >> 1: highly clustered)",
            "Coefficient of variation of inter-motif distances. CV < 0.5 indicates a nearly "
            "regular, periodic distribution. CV > 2 indicates strong spatial clustering with "
            "hotspot regions interspersed by large motif-free stretches.",
        ),
    ]), unsafe_allow_html=True)

    # ── 7. Hybrid & Cluster metrics ───────────────────────────────────────────
    _section_heading("7. Hybrid Motif & Cluster Metrics")
    _prose(
        "<p>Hybrid motifs are derived intervals where two or more distinct structural classes "
        "overlap by at least 50% of the shorter interval. Clusters are high-density windows "
        "containing ≥4 motifs from ≥3 classes within a 300 bp window.</p>"
    )
    st.markdown(_metric_table([
        (
            "Hybrid Count",
            "count of inter-class overlapping intervals",
            "0 – N",
            "Number of genomic sites where two or more non-B structures co-occupy the same "
            "region. High hybrid counts may indicate regulatory hotspots or fragile sites.",
        ),
        (
            "Hybrid Coverage %",
            "(unique bases in hybrid intervals / G) × 100",
            "0 – 100 %",
            "Fraction of the genome exhibiting simultaneous structural multiplicity. "
            "Values >1% are noteworthy in compact prokaryotic genomes.",
        ),
        (
            "Hybrid Density",
            "Hybrid Count / G  [hybrids per bp]",
            "0 – ∞",
            "Per-base frequency of hybrid structural sites.",
        ),
        (
            "Cluster Count",
            "count of 300 bp windows with ≥4 motifs from ≥3 classes",
            "0 – N",
            "Number of identified structural hotspot windows. Clusters coincide with "
            "replication origins, promoters, and recombination hotspots in model genomes.",
        ),
        (
            "Cluster Coverage %",
            "(unique bases in cluster windows / G) × 100",
            "0 – 100 %",
            "Spatial extent of structural hotspot regions.",
        ),
    ]), unsafe_allow_html=True)

    # ── 8. Class-specific scoring parameters ─────────────────────────────────
    _section_heading("8. Class-Specific Scoring Parameters")
    _prose(
        "<p>Each detector employs a biologically informed raw scoring model before normalization "
        "to the universal 1.0–3.0 scale. The table below summarises the key formula components "
        "and threshold values for each class.</p>"
    )
    th2 = "text-align:left;padding:0.5rem 0.75rem;border-bottom:2px solid #cbd5e1;background:#f1f5f9;font-size:0.8rem;color:#334155;"
    td2 = "padding:0.45rem 0.75rem;border-bottom:1px solid #e2e8f0;font-size:0.8rem;vertical-align:top;"
    class_rows = [
        ("Curved DNA",
         "Score = 1 − (|spacing − 11.0| / 1.1) for APRs; L/(L+7) for local tracts",
         "1.0 – 3.0", "9.9 – 11.1 bp phasing window; ≥3 A-tract centers required for APR",
         "Crothers et al. 1990"),
        ("Slipped DNA",
         "Score = f(L, copies, unit_size, purity, GC); Shannon entropy ≥ 0.5 bits",
         "1.0 – 3.0", "STR: units 1–9 nt, ≥4 copies; DR: units 10–50 nt, ≥2 copies",
         "Pearson et al. 2005"),
        ("Cruciform",
         "Score = arm_length × symmetry × GC_factor; ΔG ≤ −5.0 kcal/mol required",
         "1.0 – 3.0", "Arm 8–50 bp; spacer ≤12 bp; nearest-neighbour ΔG filter",
         "Lilley 1980"),
        ("R-Loop",
         "Score = GC_content × G_density across RIZ + REZ",
         "1.0 – 3.0", "G-content ≥50% in RIZ; GC ≥40% in REZ (≤2000 bp)",
         "Jenjaroenpun et al. 2015"),
        ("Triplex",
         "Score = 0.35×arm + 0.20×loop_penalty + 0.30×purity + 0.15×interruptions",
         "1.0 – 3.0", "Arm 10–100 bp; purine/pyrimidine purity ≥90%; spacer ≤8 bp",
         "Frank-Kamenetskii & Mirkin 1995"),
        ("G-Quadruplex",
         "G4Hunter: sliding-window 25 nt; score = max_G_density × length_factor",
         "1.0 – 3.0 (raw G4Hunter ≥0.5)", "8 subclasses; Telomeric >(TTAGGG){4}; canonical G{3+}N{1–7}×4",
         "Bedrat et al. 2016"),
        ("i-Motif",
         "Score = C_density + cytosine_fraction + loop_compactness",
         "1.0 – 3.0", "C{3+}N{1–7}×4; subclasses: canonical (loops 1–7 nt) vs relaxed (1–12 nt)",
         "Zeraati et al. 2018"),
        ("Z-DNA",
         "Cumulative 10-mer propensity score ≥50.0; eGZ: CGG/GGC/CCG/GCC ≥4 repeats",
         "1.0 – 3.0 (raw 50–2000)", "10-mer propensity table from Ho et al. 1986",
         "Ho et al. 1986"),
        ("A-philic DNA",
         "Cumulative log₂ A-form propensity sum ≥0.5 per merged region",
         "1.0 – 3.0", "log₂ odds from A-form vs B-form crystal structures; pseudocount 0.5",
         "Vinogradov 2003"),
    ]
    rows_html = "".join(
        f"<tr>"
        f"<td style='{td2}font-weight:600;color:#1e40af;'>{c}</td>"
        f"<td style='{td2}font-family:\"Courier New\",monospace;font-size:0.75rem;color:#0f766e;'>{fo}</td>"
        f"<td style='{td2}color:#7c3aed;'>{sc}</td>"
        f"<td style='{td2}color:#334155;'>{th_val}</td>"
        f"<td style='{td2}color:#64748b;font-style:italic;'>{ref}</td>"
        "</tr>"
        for c, fo, sc, th_val, ref in class_rows
    )
    st.markdown(
        "<table style='width:100%;border-collapse:collapse;font-family:Georgia,serif;margin-bottom:1.2rem;'>"
        "<thead><tr>"
        f"<th style='{th2}'>Class</th>"
        f"<th style='{th2}'>Scoring Formula</th>"
        f"<th style='{th2}'>Score Range</th>"
        f"<th style='{th2}'>Key Threshold(s)</th>"
        f"<th style='{th2}'>Reference</th>"
        "</tr></thead>"
        f"<tbody>{rows_html}</tbody></table>",
        unsafe_allow_html=True,
    )

    # ── 9. Quick-reference interpretation guide ───────────────────────────────
    _section_heading("9. Quick-Reference Interpretation Guide")
    cards = [
        ("#dbeafe", "#1e40af",
         "Low non-B DNA burden",
         "Coverage% < 5%, Density < 2 /kb, SLI < 0.05",
         "Typical for compact genomes with few repeats (e.g., some bacteria). "
         "Limited propensity for structural mutagenesis."),
        ("#dcfce7", "#15803d",
         "Moderate non-B DNA burden",
         "Coverage% 5–20%, Density 2–15 /kb, SLI 0.05–0.20",
         "Common in most prokaryotes and simple eukaryotes. Balance between "
         "regulatory flexibility and genomic stability."),
        ("#fef9c3", "#a16207",
         "High non-B DNA burden",
         "Coverage% 20–40%, Density 15–40 /kb, SLI 0.20–0.50",
         "Typical for repeat-rich organisms (e.g., S. cerevisiae). "
         "Elevated risk of replication stress and mutagenesis."),
        ("#fce7f3", "#9d174d",
         "Very high non-B DNA burden",
         "Coverage% > 40%, Density > 40 /kb, SLI > 0.50",
         "Indicates extreme repeat content or structural enrichment. "
         "Found in pathological expansions or highly GC-biased genomes."),
    ]
    cards_html = '<div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(260px,1fr));gap:0.6rem;margin-bottom:1.2rem;">'
    for bg, col, title, metrics, desc in cards:
        cards_html += (
            f"<div style='background:{bg};border:1.5px solid {col};border-radius:8px;"
            f"padding:0.35rem 0.9rem;'>"
            f"<div style='color:{col};font-weight:700;font-size:0.95rem;margin-bottom:0.3rem;'>{title}</div>"
            f"<div style='color:#334155;font-size:0.78rem;font-weight:600;margin-bottom:0.3rem;'>{metrics}</div>"
            f"<div style='color:#475569;font-size:0.78rem;line-height:1.4;'>{desc}</div>"
            "</div>"
        )
    cards_html += "</div>"
    st.markdown(cards_html, unsafe_allow_html=True)

    # ── 10. Diversity interpretation ──────────────────────────────────────────
    _section_heading("10. Structural Diversity Interpretation")
    st.markdown(_metric_table([
        ("D = 0.0 – 0.3", "Low diversity", "1–2 dominant classes",
         "The structural landscape is dominated by one or two motif types. "
         "Common in AT-rich genomes (Curved DNA, Slipped DNA dominate) or GC-rich genomes "
         "(G4, Z-DNA, i-Motif dominate)."),
        ("D = 0.3 – 0.6", "Moderate diversity", "3–5 roughly equal classes",
         "Multiple structural classes contribute substantially. "
         "Typical for balanced eukaryotic genomes."),
        ("D = 0.6 – 0.8", "High diversity", "6–8 effectively equal classes",
         "Rich multi-class structural landscape. Often associated with complex regulatory "
         "regions such as promoters, origins of replication, and enhancers."),
        ("D = 0.8 – 1.0", "Very high diversity", "≥9 effectively equal classes",
         "Exceptionally balanced distribution across all structural classes. "
         "Rare in natural genomic DNA; may indicate a synthetic or shuffled sequence."),
    ]), unsafe_allow_html=True)


def _tab_novelty_scope():
    _section_heading("Novelty, Prior-Work Comparison & Future Scope")
    _prose(
        "<p>This page documents how NonBDNAFinder's comparative nine-genome analysis advances the field "
        "beyond existing non-B DNA tools and studies, and outlines the most productive research directions "
        "identified from the results.</p>"
    )

    # ── 1. What Makes This Analysis Novel ────────────────────────────────────
    _section_heading("1. What Makes This Analysis Novel")
    _prose(
        "<p>The present comparative study breaks new ground across five dimensions:</p>"
    )
    _prose(
        "<ul style='line-height:1.9;'>"
        "<li><strong>Simultaneous nine-class detection.</strong> All nine primary structural classes "
        "(Curved DNA, A-philic DNA, Cruciform, Slipped DNA, Triplex, Z-DNA, R-Loop, G-Quadruplex, "
        "i-Motif) are detected in a single pipeline run, enabling discovery of hybrid overlaps and "
        "cluster hotspots that are invisible when classes are analysed independently.</li>"
        "<li><strong>22-canonical-subclass framework.</strong> A fixed vocabulary of 22 biologically "
        "meaningful subclass labels allows standardised cross-genome benchmarking — the first published "
        "comparative subclass-level non-B DNA landscape across nine organisms.</li>"
        "<li><strong>Seven derived structural metrics.</strong> Structural Complexity Index (SCI), "
        "Weighted Structural Coverage (WSC), Coefficient of Variation of inter-motif distance (CV), "
        "Structural Load Index (SLI), Structural Intensity (SI), Simpson Diversity Index (D), and "
        "Effective Class Number (N_eff) characterise each genome beyond simple motif counts.</li>"
        "<li><strong>Widest GC range in any comparative non-B DNA study.</strong> 17.6 % "
        "(<em>Ca. Carsonella</em>) to 76.2 % (<em>M. marina</em>) — a 58.6 percentage-point span — "
        "reveals composition–structure relationships invisible in studies restricted to organisms "
        "within a 30–60 % GC band.</li>"
        "<li><strong>Hybrid and cluster co-occurrence analysis.</strong> Hybrid regions and "
        "Non-B DNA Cluster regions are formally detected and characterised as analytical layers rather "
        "than filtered artefacts. The strong hybrid–cluster co-variation (r = 0.72) is the first "
        "quantitative evidence that binary co-habituation and higher-order clustering are mechanistically "
        "coupled genome-wide.</li>"
        "</ul>"
    )

    # ── 2. Comparison with Prior Non-B DNA Studies ────────────────────────────
    _section_heading("2. Comparison with Prior Non-B DNA Studies")

    th_c = "text-align:left;padding:0.5rem 0.75rem;border-bottom:2px solid #cbd5e1;background:#f1f5f9;font-size:0.8rem;color:#334155;"
    td_c = "padding:0.45rem 0.75rem;border-bottom:1px solid #e2e8f0;font-size:0.8rem;vertical-align:top;"

    comparison_rows = [
        ("Huppert &amp; Balasubramanian (2005)",
         "1 (human)", "1 (G4)", "None", "No", "~41 %",
         "Regex-based G4 census; no thermodynamic scoring or cross-class analysis."),
        ("Cer et al. — NonB DB (2013)",
         "1 (human)", "7", "None", "No", "~41 %",
         "Multi-class human-only catalogue; no composite metrics, hybrids, or clusters."),
        ("Georgakopoulos-Soares et al. (2022)",
         "1 (human)", "4", "None", "No", "~41 %",
         "Promoter-focused MPRA study; genome-wide baseline and prokaryotes absent."),
        ("Wolfl et al. (1995)",
         "3 prokaryotes", "1 (Z-DNA)", "None", "No", "40–65 %",
         "Z-DNA propensity survey; no multi-class or diversity metrics."),
        ("<strong>This study (NonBDNAFinder)</strong>",
         "<strong>9 (3 domains)</strong>",
         "<strong>9 + hybrid + cluster</strong>",
         "<strong>7 metrics</strong>",
         "<strong>Yes</strong>",
         "<strong>18–76 %</strong>",
         "<strong>First complete nine-class comparative landscape across a 58.6 pp GC range.</strong>"),
    ]
    rows_html = "".join(
        f"<tr>"
        f"<td style='{td_c}color:#1e40af;'>{study}</td>"
        f"<td style='{td_c}'>{orgs}</td>"
        f"<td style='{td_c}'>{cls}</td>"
        f"<td style='{td_c}'>{metrics}</td>"
        f"<td style='{td_c}'>{hybrid}</td>"
        f"<td style='{td_c}color:#7c3aed;'>{gc}</td>"
        f"<td style='{td_c}color:#334155;font-style:italic;'>{note}</td>"
        "</tr>"
        for study, orgs, cls, metrics, hybrid, gc, note in comparison_rows
    )
    st.markdown(
        "<table style='width:100%;border-collapse:collapse;font-family:Georgia,serif;margin-bottom:1.2rem;'>"
        "<thead><tr>"
        f"<th style='{th_c}'>Study</th>"
        f"<th style='{th_c}'>Organisms</th>"
        f"<th style='{th_c}'>Classes</th>"
        f"<th style='{th_c}'>Composite Metrics</th>"
        f"<th style='{th_c}'>Hybrid / Cluster</th>"
        f"<th style='{th_c}'>GC Range</th>"
        f"<th style='{th_c}'>Key Gap Addressed</th>"
        "</tr></thead>"
        f"<tbody>{rows_html}</tbody></table>",
        unsafe_allow_html=True,
    )

    _prose(
        "<p><strong>G4-only studies</strong> (Huppert &amp; Balasubramanian, 2005; "
        "Hansel-Hertsch <em>et al.</em>, 2017) demonstrate that G4 abundance correlates with cancer "
        "transcriptional dysregulation in human cells. NonBDNAFinder's nine-genome results add the "
        "prokaryotic perspective: G4 density in GC-rich bacteria (<em>C. shaoxiangyii</em>, "
        "<em>M. marina</em>) exceeds published human G4 densities by &gt;3-fold per Mb, implying that "
        "bacterial structural genomics is an underexplored source of functional G4 biology.</p>"
        "<p><strong>Z-DNA surveys</strong> (Wolfl <em>et al.</em>, 1995; Zhang <em>et al.</em>, 2006) "
        "placed typical Z-DNA density at 0.5–5 % of prokaryotic genomes. Our GC-rich bacteria show "
        "densities 10–100× higher — establishing a new empirical upper bound and demonstrating that "
        "Z-DNA content is primarily a compositional rather than a functional selection signature.</p>"
        "<p><strong>NonB DB</strong> (Cer <em>et al.</em>, 2013) is restricted to the human genome and "
        "reports raw class counts without composite structural metrics. The present study applies an "
        "equivalent multi-class detection to nine organisms spanning three domains of life, showing that "
        "the structural landscapes of GC-rich bacteria rival or exceed the human genome in certain class "
        "densities.</p>"
    )

    # ── 3. Future Research Directions ────────────────────────────────────────
    _section_heading("3. Future Research Directions")

    future_items = [
        ("Near-term", "#dbeafe", "#1e40af", [
            ("Hundred-genome cohort",
             "Apply NonBDNAFinder to 100–500 bacterial genomes stratified by GC content (10 % bins, 20–80 %) "
             "to derive robust regression models of non-B DNA density vs. base composition."),
            ("Functional annotation overlay",
             "Integrate RefSeq gene annotation (CDS, promoters, replication origins) with non-B DNA coordinates "
             "to test whether G4/Z-DNA are enriched at promoters in GC-rich bacteria."),
            ("Archaea and organellar genomes",
             "Extend kingdom representation to archaea and mitochondria/plastids to test whether "
             "endosymbiotic reductive evolution universally drives Curved DNA dominance."),
        ]),
        ("Medium-term", "#dcfce7", "#15803d", [
            ("Strand-specific analysis",
             "Separate plus- and minus-strand detection for R-loops, G4, and i-Motif — inherently "
             "strand-specific structures — and compute strand-asymmetry indices."),
            ("Multi-strain comparison",
             "Compare non-B DNA landscapes across strains of the same organism (e.g., H. pylori reference "
             "vs. clinical isolates) to identify strain-specific structural hotspots linked to virulence."),
            ("Transcriptomics / epigenomics integration",
             "Cross-reference non-B DNA coordinates with RNA-seq, ATAC-seq, and BS-seq to determine whether "
             "high-structural-burden regions are transcriptionally active, accessible, or methylated."),
        ]),
        ("Long-term", "#fce7f3", "#9d174d", [
            ("Machine learning for structural prediction",
             "Train classification models on the nine-genome dataset to predict structural class composition "
             "from sequence features alone, enabling rapid screening of newly sequenced genomes."),
            ("Experimental validation of predicted hotspots",
             "Prioritise 7-class cluster regions in C. shaoxiangyii and R-loop-rich regions in H. pylori "
             "for DMS-MaPseq, G4-ChIP-seq, or R-ChIP experimental validation."),
            ("Clinical antimicrobial application",
             "Develop a pipeline that flags non-B DNA hotspots in pathogen genomes as candidate antimicrobial "
             "target sites, leveraging the observation that H. pylori's R-loop burden is exceptionally high."),
        ]),
    ]

    for phase, bg, col, items in future_items:
        st.markdown(
            f"<div style='background:{bg};border-left:4px solid {col};border-radius:6px;"
            f"padding:0.3rem 1rem;margin-bottom:0.6rem;'>"
            f"<div style='color:{col};font-weight:700;font-size:0.95rem;margin-bottom:0.4rem;'>"
            f"{phase}</div>"
            + "".join(
                f"<div style='margin-bottom:0.5rem;'>"
                f"<span style='color:{col};font-weight:600;font-size:0.85rem;'>▸ {title}:</span> "
                f"<span style='color:#334155;font-size:0.82rem;'>{desc}</span></div>"
                for title, desc in items
            )
            + "</div>",
            unsafe_allow_html=True,
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
        "color:#334155;background:#f8fafc;padding:0.4rem 1rem;border-left:4px solid #2563eb;"
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
        "🔬 Overview & Architecture",
        "🧬 Motif Library & Algorithms",
        "📊 Scoring & Analysis",
        "📈 Statistics Guide",
        "📚 References & Citation",
    ])

    with tabs[0]:
        _tab_overview()
        st.markdown("<hr style='border:1px solid #e2e8f0;margin:1.5rem 0;'>", unsafe_allow_html=True)
        _tab_architecture()

    with tabs[1]:
        _tab_motif_library()
        st.markdown("<hr style='border:1px solid #e2e8f0;margin:1.5rem 0;'>", unsafe_allow_html=True)
        _tab_detection_algorithms()

    with tabs[2]:
        _tab_scoring()
        st.markdown("<hr style='border:1px solid #e2e8f0;margin:1.5rem 0;'>", unsafe_allow_html=True)
        _tab_hybrid_clustering()
        st.markdown("<hr style='border:1px solid #e2e8f0;margin:1.5rem 0;'>", unsafe_allow_html=True)
        _tab_optimization()
        st.markdown("<hr style='border:1px solid #e2e8f0;margin:1.5rem 0;'>", unsafe_allow_html=True)
        _tab_validation()

    with tabs[3]:
        _tab_statistical_guide()

    with tabs[4]:
        _tab_references()
        st.markdown("<hr style='border:1px solid #e2e8f0;margin:1.5rem 0;'>", unsafe_allow_html=True)
        _tab_citation()
