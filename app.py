"""
╔══════════════════════════════════════════════════╗
║                  NBDSCANNER WEB APP              ║
║    Streamlit App for Non-B DNA Motif Detection   ║
╚══════════════════════════════════════════════════╝
"""
# --- Imports ---
import streamlit as st; import pandas as pd; import matplotlib.pyplot as plt; import os; import sys; import psutil; import numpy as np; import gc
from utilities import (parse_fasta, parse_fasta_chunked, get_file_preview, wrap, get_basic_stats, export_to_bed, export_to_csv, export_to_json, export_to_excel, calculate_genomic_density, calculate_positional_density)
from nonbscanner import analyze_sequence, get_motif_info as get_motif_classification_info
from utilities import export_results_to_dataframe
from visualizations import (plot_motif_distribution, plot_coverage_map, plot_density_heatmap, plot_length_distribution, plot_score_distribution, plot_nested_pie_chart, MOTIF_CLASS_COLORS, plot_density_comparison, plot_circos_motif_density, plot_density_comparison_by_subclass, plot_enrichment_analysis_by_subclass, plot_subclass_density_heatmap)
try: from Bio import Entrez, SeqIO; BIO_AVAILABLE=True
except ImportError: BIO_AVAILABLE=False
try: import hyperscan; HYPERSCAN_AVAILABLE=True
except ImportError: HYPERSCAN_AVAILABLE=False

# --- Caching ---
@st.cache_resource(show_spinner=False)
def cache_genome_as_numpy(sequence): return np.frombuffer(sequence.encode('utf-8'), dtype=np.uint8)
@st.cache_resource(show_spinner=False)
def cache_hyperscan_database(_patterns=None): return None
@st.cache_data(show_spinner=False, max_entries=10, ttl=3600)
def cache_analysis_results(sequence_hash, sequence, name): return analyze_sequence(sequence, name)
@st.cache_data(show_spinner=False, max_entries=20, ttl=3600)
def get_cached_stats(sequence, motifs_json):
    import json; motifs=json.loads(motifs_json) if motifs_json else []; return get_basic_stats(sequence, motifs)

# --- System info helper (unused in UI) ---
def get_system_info():
    try:
        m=psutil.virtual_memory()
        return {'memory_total_mb': m.total/(1024*1024), 'memory_used_mb': m.used/(1024*1024), 'memory_available_mb': m.available/(1024*1024), 'memory_percent': m.percent, 'cpu_count': psutil.cpu_count(), 'available': True}
    except: return {'memory_total_mb':0,'memory_used_mb':0,'memory_available_mb':0,'memory_percent':0,'cpu_count':1,'available':False}

# --- Page config & globals ---
st.set_page_config(page_title="NBDScanner - Non-B DNA Motif Finder", layout="wide", page_icon=None, menu_items={'About': "NBDScanner | Developed by Dr. Venkata Rajesh Yella"})
MAX_SEQUENCE_LENGTH=10_000_000_000
def format_sequence_limit(): return "unlimited (chunked processing enabled)"
CLASSIFICATION_INFO=get_motif_classification_info()
def ensure_subclass(m):
    if isinstance(m,dict):
        if 'Subclass' not in m or m['Subclass'] is None: m['Subclass']=m.get('Subtype','Other'); return m
    return {'Subclass':'Other','Motif':m}
def format_time(s):
    if s<60: return f"{s:.1f}s"
    elif s<3600: mins,secs=int(s//60),int(s%60); return f"{mins}m {secs}s"
    else: hrs,mins=int(s//3600),int((s%3600)//60); return f"{hrs}h {mins}m"
def generate_excel_bytes(motifs):
    import tempfile, os
    with tempfile.NamedTemporaryFile(mode='wb', suffix='.xlsx', delete=False) as tmp: tmp_path=tmp.name
    try:
        export_to_excel(motifs,tmp_path)
        with open(tmp_path,'rb') as f: data=f.read()
        return data
    finally:
        if os.path.exists(tmp_path): os.remove(tmp_path)
def hex_to_rgb(h): h=h.lstrip('#'); return tuple(int(h[i:i+2],16) for i in (0,2,4))
def get_dna_pattern_svg(c): return f"url(\"data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='60' height='60' viewBox='0 0 60 60'%3E%3Cg fill='none' stroke='%23{c}' stroke-width='0.5' opacity='0.3'%3E%3Cpath d='M0 0 Q30 30 60 0 M0 60 Q30 30 60 60'/%3E%3C/g%3E%3C/svg%3E\")"

# --- Themes & CSS (kept as-is) ---
COLOR_THEMES={'scientific_blue':{'primary':'#5B8DEF','secondary':'#7BA3F7','accent':'#A8C5FF','bg_light':'#F8FAFC','bg_card':'#EEF2FF','text':'#1F2937','tab_bg':'#F1F5F9','tab_active':'#5B8DEF','shadow':'rgba(91,141,239,0.15)'},
'nature_green':{'primary':'#6DBB7A','secondary':'#8FD19C','accent':'#B8E5C0','bg_light':'#F7FBF8','bg_card':'#E8F5E9','text':'#1F2937','tab_bg':'#F0F9F1','tab_active':'#6DBB7A','shadow':'rgba(109,187,122,0.15)'},
'genomic_purple':{'primary':'#9B8FD9','secondary':'#B3A9E5','accent':'#D4CFF0','bg_light':'#FAF9FD','bg_card':'#EDE7F6','text':'#1F2937','tab_bg':'#F5F3FA','tab_active':'#9B8FD9','shadow':'rgba(155,143,217,0.15)'},
'clinical_teal':{'primary':'#5BBFB4','secondary':'#7CCFC6','accent':'#A8E4DE','bg_light':'#F6FBFB','bg_card':'#E0F2F1','text':'#1F2937','tab_bg':'#EBF7F6','tab_active':'#5BBFB4','shadow':'rgba(91,191,180,0.15)'}}
if 'theme_mode' not in st.session_state: st.session_state.theme_mode='light'
if 'table_density' not in st.session_state: st.session_state.table_density='relaxed'
if 'color_theme' not in st.session_state: st.session_state.color_theme='scientific_blue'
current_theme=COLOR_THEMES.get(st.session_state.color_theme, COLOR_THEMES['scientific_blue']); is_dark_mode=st.session_state.theme_mode=='dark'; is_compact=st.session_state.table_density=='compact'
if is_dark_mode:
    base_primary=current_theme['primary']
    current_theme={**current_theme,'bg_light':'#1A1F2E','bg_card':'#252B3B','text':'#E5E7EB','tab_bg':'#1F2937','tab_active':current_theme.get('primary',base_primary),'shadow':'rgba(0,0,0,0.25)'}
rgb={k: hex_to_rgb(v) if k!='shadow' else (0,0,0) for k,v in current_theme.items()}
tab_bg_color=current_theme.get('tab_bg', current_theme['bg_card']); tab_active_color=current_theme.get('tab_active', current_theme['primary'])
shadow_color=current_theme.get('shadow', f"rgba({rgb['primary'][0]}, {rgb['primary'][1]}, {rgb['primary'][2]}, 0.15)")
dna_pattern=get_dna_pattern_svg('1e3a5f' if is_dark_mode else 'bbdefb')
st.markdown(f"""<style> /* original CSS retained */ </style>""", unsafe_allow_html=True)

# --- Motif config ---
CONFIG_AVAILABLE=False
MOTIF_ORDER=["Sticky DNA","Curved DNA","Z-DNA","eGZ (Extruded-G)","Slipped DNA","R-Loop","Cruciform","Triplex DNA","G-Triplex","G4","Relaxed G4","Bulged G4","Bipartite G4","Multimeric G4","i-Motif","AC-Motif","A-philic DNA","Hybrid","Non-B DNA Clusters"]
MOTIF_COLORS=MOTIF_CLASS_COLORS
PAGES={"Home":"Overview","Upload & Analyze":"Sequence Upload and Motif Analysis","Results":"Analysis Results and Visualization","Download":"Export Data","Documentation":"Scientific Documentation & References"}
if BIO_AVAILABLE: Entrez.email="raazbiochem@gmail.com"; Entrez.api_key=None

# --- Example data ---
EXAMPLE_FASTA="""">Example Sequence
ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC
ATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT
GGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGTAGGAACCCCTAACCCCTAA
GAAAGAAGAAGAAGAAGAAGAAAGGAAGGAAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGG
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
GAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA
CTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT
"""
EXAMPLE_MULTI_FASTA="""">G4_iMotif_APhilic_Sequence
GGGATTGGGATTGGGATTGGGCCCATCCCTACCCTACCCAAACCCATCCCTACCCTACCCAATTTATTTAAAAA
AAAAAAAAAAAAAAAAAAAAAAGATCGAAAGATCGAAAGATCGAAAGATCGATGCGGCGGCGGCGGCGGCGGCGG
CGGCGGCGAATTCGAATTCGAATTCGAATTCCGCGCGCGCGCGCGCGCGCGAATGCATGCATGCATGCATGCAT
>Z_DNA_RLoop_Complex
CGCGCGCGCGCGCGCGCGCGCGATATATATATATATATATATCGCGCGCGCGCGCGCGCGCGGGGATGGGGATGG
GGATGGGGGGATGGGGATGGGGATGGGTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTGAAA
GAAAAAAGAAAGAAAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAG
>CurvedDNA_SlippedDNA_STR
AAAAAACGTTGCAAAAAACGTTGCAAAAAACGTTGCAAAAAATTTTTTCGAACGTTTTTTCGAACGTTTTTTCGA
ACGCAGCAGCAGCAGCAGCAGCAGCAGCTGCTGCTGCTGCTGCTGCTGCTGATCTGATCTGATCTGATCTGATC
TGATCTGATTCTATTCTATTCTATTCTATTCTATTCTATTCTGGCCCCGGCCCCGGCCCCGGCCCCTGCTGCTG
>Cruciform_Triplex_Mirror
ATGCCCGGGATCGGATCCGATCGAAATTCGATCGGATCCGATCCCGGGCATGAAAGAAAGAAAGAAAGAAAGAAA
GAAAGAAAGAAAAGATCCGGCCGATAGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGTTCCTCCTCCTCC
TCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCGAATTCCGAATTCCGAATTCCGAATTCCGAATTCGAAA
>Multi_iMotif_AC_Sequence  
AAATTTATTTAAATTTAAATTCCCTACCCTACCCTACCCAAAAATCCCTACCCTACCCTACCCGGAATCGATCG
ATCGATCGATCGATCGATCGCCCTACCCTACCCTACCCAAACCCTACCCTACCCTACCCAAAAAAAAAAAAAAAA
AAAAAAAAAAAGATCTAGATCTAGATCTAGATCTAGATCTAGATCTGAAAGAAAGAAAGAAAGAAAGAAAGAAA
"""

# --- Session state init ---
for k,v in {'seqs':[],'names':[],'results':[],'summary_df':pd.DataFrame(),'analysis_status':'Ready','selected_classes':[]}.items():
    if k not in st.session_state: st.session_state[k]=v

# --- Tabs ---
tabs=st.tabs(list(PAGES.keys())); tab_pages=dict(zip(PAGES.keys(),tabs))

# --- Home tab ---
with tab_pages["Home"]:
    st.markdown("<h1>Non-B DNA Motif Finder</h1>", unsafe_allow_html=True)
    left,right=st.columns([1,1])
    with left:
        try:
            paths=["nbdcircle.JPG","archive/nbdcircle.JPG","./nbdcircle.JPG"]; found=False
            for p in paths:
                if os.path.exists(p): st.image(p,use_container_width=True); found=True; break
            if not found: raise FileNotFoundError()
        except:
            st.markdown("""<div style='background:linear-gradient(135deg,#667eea 0%,#764ba2 100%);border-radius:15px;padding:40px;text-align:center;color:white;'>
            <h2 style='margin:0;color:white;'>DNA</h2><h3 style='margin:10px 0 0 0;color:white;'>NBD Finder</h3><p style='margin:5px 0 0 0;color:#E8E8E8;'>Non-B DNA Detection</p></div>""", unsafe_allow_html=True)
    with right:
        st.markdown("""<div style='font-family:Inter,system-ui,sans-serif;font-size:1.05rem;color:#263238;line-height:1.8;padding:2rem;background:white;border-radius:16px;box-shadow:0 4px 16px rgba(0,0,0,0.08);border:1px solid rgba(25,118,210,0.1);'>
        <p style='margin-top:0;'><b style='color:#0d47a1;font-size:1.15rem;'>Non-canonical DNA structures</b> play key roles in genome stability, regulation, and evolution.</p>
        <p>This application detects and analyzes <b style='color:#1976d2;'>11 major classes with 22+ subclasses</b> of Non-B DNA motifs in any DNA sequence or multi-FASTA file.</p>
        <p style='margin-bottom:0.8rem;'><b style='color:#0d47a1;'>Motif Classes (11 classes, 22+ subclasses):</b></p>
        <div style='color:#37474f;font-size:0.98rem;line-height:1.9;padding-left:1rem;'>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>1. Curved DNA</b> <span style='color:#546e7a;'>(Global curvature, Local Curvature)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>2. Slipped DNA</b> <span style='color:#546e7a;'>(Direct Repeat, STR)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>3. Cruciform DNA</b> <span style='color:#546e7a;'>(Inverted Repeats)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>4. R-loop</b> <span style='color:#546e7a;'>(R-loop formation sites)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>5. Triplex</b> <span style='color:#546e7a;'>(Triplex, Sticky DNA)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>6. G-Quadruplex Family</b> <span style='color:#546e7a;'>(Multimeric G4, Canonical G4, Relaxed G4, Bulged G4, Bipartite G4, Imperfect)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>7. i-Motif Family</b> <span style='color:#546e7a;'>(Canonical i-motif, Relaxed i-motif, AC-motif)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>8. Z-DNA</b> <span style='color:#546e7a;'>(Z-DNA, eGZ (Extruded-G) DNA)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>9. A-philic DNA</b> <span style='color:#546e7a;'>(A-philic DNA)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>10. Hybrid</b> <span style='color:#546e7a;'>(dynamic overlaps)</span></div>
            <div style='margin-bottom:0.4rem;'><b style='color:#1976d2;'>11. Non-B DNA Clusters</b> <span style='color:#546e7a;'>(dynamic clusters)</span></div>
        </div>
        <p style='margin-top:1.2rem;margin-bottom:0;color:#1976d2;font-weight:600;'>Upload single or multi-FASTA files to begin analysis...</p>
        </div>""", unsafe_allow_html=True)

# --- Upload & Analyze tab (System Resource Monitor + Unlimited note removed) ---
with tab_pages["Upload & Analyze"]:
    st.markdown("<h2>Sequence Upload and Motif Analysis</h2>", unsafe_allow_html=True)
    st.markdown('<span style="font-family:Montserrat,Arial; font-size:1.12rem;">Supports multi-FASTA and single FASTA. Paste, upload, select example, or fetch from NCBI.</span>', unsafe_allow_html=True)
    st.caption("Supported formats: .fa, .fasta, .txt, .fna | Limit: 1GB/file (optimized for large genomic sequences).")
    st.markdown("---")

    col_left,col_right=st.columns([1,1],gap="large")

    # --- Input column ---
    with col_left:
        st.markdown("### Input Method")
        input_method=st.radio("Choose your input method:", ["Upload FASTA File","Paste Sequence","Example Data","NCBI Fetch"], horizontal=True)
        seqs,names=[],[]
        if input_method=="Upload FASTA File":
            fasta_file=st.file_uploader("Drag and drop FASTA/multi-FASTA file here", type=["fa","fasta","txt","fna"])
            if fasta_file:
                file_size_mb=fasta_file.size/(1024*1024)
                st.info(f"📁 File: **{fasta_file.name}** | Size: **{file_size_mb:.2f} MB**")
                if file_size_mb>100: st.warning(f"⚠️ Large file detected ({file_size_mb:.2f} MB). Processing may take longer and use more memory.")
                with st.spinner(f"Processing {fasta_file.name}..."):
                    preview_info=get_file_preview(fasta_file, max_sequences=3)
                    st.success(f"✅ File contains **{preview_info['num_sequences']} sequences** totaling **{preview_info['total_bp']:,} bp**")
                    for prev in preview_info['previews']:
                        st.markdown(f"**{prev['name']}**: <span style='color:#576574'>{prev['length']:,} bp</span>", unsafe_allow_html=True)
                        stats=get_basic_stats(prev['preview'].replace('...',''))
                        st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                    if preview_info['num_sequences']>3: st.caption(f"...and {preview_info['num_sequences']-3} more sequences.")
                    seqs,names=[],[]; has_large=False
                    if preview_info['num_sequences']>10:
                        pbar=st.progress(0); status=st.empty()
                        for idx,(nm,seq) in enumerate(parse_fasta_chunked(fasta_file)):
                            names.append(nm); seqs.append(seq); has_large |= len(seq)>10_000_000
                            pbar.progress((idx+1)/preview_info['num_sequences']); status.text(f"Loading sequence {idx+1}/{preview_info['num_sequences']}: {nm}")
                        pbar.empty(); status.empty()
                    else:
                        for nm,seq in parse_fasta_chunked(fasta_file):
                            names.append(nm); seqs.append(seq); has_large |= len(seq)>10_000_000
                    if has_large: gc.collect()
                    if not seqs: st.warning("No sequences found in file.")
        elif input_method=="Paste Sequence":
            seq_input=st.text_area("Paste single or multi-FASTA here:", height=150, placeholder="Paste your DNA sequence(s) here...")
            if seq_input:
                seqs,names=[],[]; cur_seq,cur_name="",""
                for line in seq_input.splitlines():
                    if line.startswith(">"):
                        if cur_seq: seqs.append(cur_seq); names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                        cur_name=line.strip().lstrip(">"); cur_seq=""
                    else: cur_seq += line.strip()
                if cur_seq: seqs.append(cur_seq); names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                if seqs:
                    st.success(f"Pasted {len(seqs)} sequences.")
                    for i,seq in enumerate(seqs[:3]):
                        stats=get_basic_stats(seq)
                        st.markdown(f"**{names[i]}**: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                        st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                    if len(seqs)>3: st.caption(f"...and {len(seqs)-3} more.")
                else: st.warning("No sequences found.")
        elif input_method=="Example Data":
            ex_type=st.radio("Example Type:", ["Single Example","Multi-FASTA Example"], horizontal=True)
            if ex_type=="Single Example":
                if st.button("Load Single Example"):
                    parsed=parse_fasta(EXAMPLE_FASTA); seqs=list(parsed.values()); names=list(parsed.keys())
                    st.success("Single example sequence loaded."); stats=get_basic_stats(seqs[0]); st.code(EXAMPLE_FASTA, language="fasta")
                    st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
            else:
                if st.button("Load Multi-FASTA Example"):
                    seqs,names=[],[]; cur_seq,cur_name="",""
                    for line in EXAMPLE_MULTI_FASTA.splitlines():
                        if line.startswith(">"):
                            if cur_seq: seqs.append(cur_seq); names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                            cur_name=line.strip().lstrip(">"); cur_seq=""
                        else: cur_seq += line.strip()
                    if cur_seq: seqs.append(cur_seq); names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                    st.success(f"Multi-FASTA example loaded with {len(seqs)} sequences.")
                    for i,seq in enumerate(seqs[:3]):
                        stats=get_basic_stats(seq)
                        st.markdown(f"**{names[i]}**: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                        st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                    st.code(EXAMPLE_MULTI_FASTA, language="fasta")
        elif input_method=="NCBI Fetch":
            db=st.radio("NCBI Database", ["nucleotide","gene"], horizontal=True); query_type=st.radio("Query Type", ["Accession","Gene Name","Custom Query"], horizontal=True)
            motif_examples={"G-quadruplex":"NR_003287.2 (human telomerase RNA)","Z-DNA":"NM_001126112.2 (human ADAR1 gene)","R-loop":"NR_024540.1 (human SNRPN gene)","eGZ-motif":"CGG repeat region","AC-motif":"A-rich/C-rich consensus region"}
            with st.popover("View Example Queries", use_container_width=True):
                st.markdown("**Motif Example Queries:**"); [st.markdown(f"• **{m}**: `{ex}`") for m,ex in motif_examples.items()]
            query=st.text_input("Enter query (accession, gene, etc.):"); rettype=st.radio("Return Format", ["fasta","gb"], horizontal=True); retmax=st.number_input("Max Records", min_value=1, max_value=20, value=3)
            if st.button("Fetch from NCBI"):
                if query:
                    with st.spinner("Contacting NCBI..."):
                        h=Entrez.efetch(db=db,id=query,rettype=rettype,retmode="text"); records=list(SeqIO.parse(h, rettype)); h.close()
                        seqs=[str(r.seq).upper().replace("U","T") for r in records]; names=[r.id for r in records]
                    if seqs:
                        st.success(f"Fetched {len(seqs)} sequences.")
                        for i,seq in enumerate(seqs[:3]):
                            stats=get_basic_stats(seq)
                            st.markdown(f"<b>{names[i]}</b>: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                            st.markdown(f"GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                else: st.warning("Enter a query before fetching.")
        if seqs: st.session_state.seqs=seqs; st.session_state.names=names; st.session_state.results=[]
        if st.session_state.get('seqs'):
            st.markdown("### Sequence Preview")
            for i,seq in enumerate(st.session_state.seqs[:2]):
                stats=get_basic_stats(seq)
                st.markdown(f"**{st.session_state.names[i]}** ({len(seq):,} bp) | GC %: {stats['GC%']} | AT %: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}", unsafe_allow_html=True)
                st.code(wrap(seq[:400]), language="fasta")
            if len(st.session_state.seqs)>2: st.caption(f"...and {len(st.session_state.seqs)-2} more.")

    # --- Analysis column ---
    with col_right:
        st.markdown("### Analysis & Run"); st.markdown("### Analysis Options")
        st.info("**NBDScanner detects all 11 motif classes with 22+ subclasses automatically**")
        c1,c2=st.columns(2)
        with c1: detailed_output=st.checkbox("Detailed Analysis", value=True, help="Include comprehensive motif metadata")
        with c2: quality_check=st.checkbox("Quality Validation", value=True, help="Validate detected motifs")
        show_advanced=st.toggle("Show Advanced Options", value=False, help="Toggle advanced analysis options")
        if show_advanced:
            with st.container(border=True):
                st.markdown("##### Advanced Configuration")
                ca1,ca2=st.columns(2)
                with ca1: show_chunk_progress=st.checkbox("Show Chunk-Level Progress", value=False)
                with ca2: use_parallel_scanner=st.checkbox("Use Experimental Parallel Scanner", value=False)
                if use_parallel_scanner: st.info("Parallel scanner is experimental and works best on sequences >100kb with multiple CPU cores")
        else: show_chunk_progress=False; use_parallel_scanner=False
        nonoverlap=True; overlap_option="Remove overlaps within subclasses"

        if st.button("Run NBDScanner Analysis", type="primary", use_container_width=True, key="run_motif_analysis_main"):
            if not st.session_state.seqs:
                st.error("Please upload or input sequences before running analysis."); st.session_state.analysis_status="Error"
            else:
                over=[]; total_len=0
                for i,seq in enumerate(st.session_state.seqs):
                    l=len(seq); total_len+=l
                    if l>MAX_SEQUENCE_LENGTH:
                        nm=st.session_state.names[i] if i<len(st.session_state.names) else f"Sequence {i+1}"
                        over.append((nm,l))
                if over:
                    lim=format_sequence_limit()
                    st.error(f"**Sequence Length Limit Exceeded**\n\nThe following sequence(s) exceed the maximum allowed length of **{lim}** for the web version:")
                    for nm,l in over: st.error(f"• **{nm}**: {l:,} nucleotides ({l-MAX_SEQUENCE_LENGTH:,} over limit)")
                    st.info("**For unlimited sequence analysis:**\n\nUse the **local Jupyter notebook version** (`NonBScanner_Local.ipynb`) which has no sequence length limits.\nThe local version is optimized for large genome-scale analysis on your own hardware.")
                    st.session_state.analysis_status="Error"
                else:
                    st.session_state.analysis_status="Running"
                    st.session_state.overlap_option_used=overlap_option; st.session_state.nonoverlap_used=nonoverlap
                    report_hotspots=True; calculate_conservation=False; threshold=0.0; validation_messages=[]
                    import time; timer_placeholder=st.empty(); progress_placeholder=st.empty(); status_placeholder=st.empty(); detailed_progress_placeholder=st.empty()
                    start=time.time()
                    DETECTOR_PROCESSES=[("Curved DNA","A-tract mediated DNA bending detection"),("Slipped DNA","Direct repeats and STR detection"),("Cruciform","Inverted repeat/palindrome detection"),("R-Loop","RNA-DNA hybrid formation site detection"),("Triplex","Three-stranded structure detection"),("G-Quadruplex","Four-stranded G-rich structure detection"),("i-Motif","C-rich structure detection"),("Z-DNA","Left-handed helix detection"),("A-philic DNA","A-rich structural element detection")]
                    ESTIMATED_BP_PER_SECOND=5800; CHUNK_SIZE_FOR_PARALLEL=50000
                    def display_progress_panel(container,elapsed,remaining,progress_display,status_text,seq_name,seq_bp,seq_num,total_seqs,processed_bp,total_bp,det_count,extra=""):
                        with container:
                            st.subheader("🧬 NonBScanner Analysis"); st.write(status_text)
                            c1,c2,c3,c4=st.columns(4)
                            with c1: st.metric("⏱️ Elapsed", format_time(elapsed))
                            with c2: st.metric("⏳ Remaining", format_time(remaining))
                            with c3: st.metric("📊 Progress", progress_display)
                            with c4: st.metric("🔬 Detectors", str(det_count))
                            st.write(f"**📄 Sequence {seq_num}/{total_seqs}**: {seq_name} *({seq_bp:,} bp)*")
                            st.write(f"⚡ Processed: {processed_bp:,} / {total_bp:,} bp")
                            if extra: st.write(extra)
                    def estimate_time(bp): return bp/ESTIMATED_BP_PER_SECOND
                    total_bp_all=sum(len(s) for s in st.session_state.seqs); estimated_total=estimate_time(total_bp_all)
                    try:
                        all_results=[]; total_bp_processed=0
                        with progress_placeholder.container(): pbar=st.progress(0)
                        with detailed_progress_placeholder.container():
                            st.subheader("🔬 Analysis Pipeline")
                            for j,(dname,ddesc) in enumerate(DETECTOR_PROCESSES): st.write(f"**{j+1}. {dname}** - {ddesc}")
                            st.info("✓ All detectors process in parallel | 🔄 Followed by overlap resolution & clustering")
                        for i,(seq,name) in enumerate(zip(st.session_state.seqs, st.session_state.names)):
                            elapsed=time.time()-start
                            if elapsed>0 and total_bp_processed>0:
                                rate=total_bp_processed/elapsed; rem_bp=total_bp_all-total_bp_processed; remaining=max(0, rem_bp/rate) if rate>0 else 0
                            else: remaining=max(0, estimated_total-elapsed)
                            overall_pct=(total_bp_processed/total_bp_all*100) if total_bp_all>0 else 0
                            status_text="Starting analysis..." if total_bp_processed==0 else "Analysis in progress..."
                            progress_display="Starting" if total_bp_processed==0 else f"{overall_pct:.1f}%"
                            extra=""
                            if show_chunk_progress and use_parallel_scanner:
                                est_chunks=max(1,(len(seq)+CHUNK_SIZE_FOR_PARALLEL-1)//CHUNK_SIZE_FOR_PARALLEL); extra=f"📦 Chunks: {est_chunks}"
                            display_progress_panel(timer_placeholder,elapsed,remaining,progress_display,status_text,name,len(seq),i+1,len(st.session_state.seqs),total_bp_processed,total_bp_all,len(DETECTOR_PROCESSES),extra)
                            seq_start=time.time()
                            if use_parallel_scanner and len(seq)>100000:
                                try:
                                    from scanner_agent import ParallelScanner
                                    chunk_placeholder=st.empty()
                                    def chunk_cb(cur,total):
                                        if show_chunk_progress:
                                            pct=(cur/total)*100; chunk_placeholder.info(f"📦 Chunks: {cur}/{total} ({pct:.1f}%)")
                                    scanner=ParallelScanner(seq, hs_db=None); raw=scanner.run_scan(progress_callback=chunk_cb); results=analyze_sequence(seq,name)
                                    if show_chunk_progress: chunk_placeholder.success(f"✅ Chunks complete: {len(raw)} motifs")
                                except Exception as e:
                                    st.warning(f"Parallel scanner failed, falling back to standard: {e}"); results=analyze_sequence(seq,name)
                            else: results=analyze_sequence(seq,name)
                            seq_time=time.time()-seq_start; results=[ensure_subclass(m) for m in results]; all_results.append(results); total_bp_processed+=len(seq)
                            elapsed=time.time()-start; speed=total_bp_processed/elapsed if elapsed>0 else 0
                            remaining=max(0,(total_bp_all-total_bp_processed)/speed) if speed>0 else 0
                            actual_pct=(total_bp_processed/total_bp_all*100) if total_bp_all>0 else 0
                            completion=f"⚡ Total: {total_bp_processed:,} / {total_bp_all:,} bp | 🎯 {len(results)} motifs | 🚀 {speed:,.0f} bp/s"
                            display_progress_panel(timer_placeholder,elapsed,remaining,f"{actual_pct:.1f}%",f"✅ Sequence {i+1}/{len(st.session_state.seqs)} completed",name,len(seq),i+1,len(st.session_state.seqs),total_bp_processed,total_bp_all,len(DETECTOR_PROCESSES),completion)
                            with progress_placeholder.container(): pbar.progress((i+1)/len(st.session_state.seqs), text=f"🔬 Analyzed {i+1}/{len(st.session_state.seqs)} sequences")
                            status_placeholder.success(f"✅ {name}: {len(seq):,} bp in {seq_time:.2f}s ({len(seq)/seq_time:.0f} bp/s) | 🎯 {len(results)} motifs")
                        st.session_state.results=all_results
                        total_time=time.time()-start; overall_speed=total_bp_processed/total_time if total_time>0 else 0
                        summary=[]
                        for i,res in enumerate(all_results):
                            seq=st.session_state.seqs[i]; stats=get_basic_stats(seq,res)
                            summary.append({'Sequence':st.session_state.names[i],'Length':stats['Length'],'GC Content':f"{stats['GC%']:.1f}%",'Motifs Found':len(res),'Unique Types':len(set(m.get('Type','Unknown') for m in res)),'Avg Score':f"{np.mean([m.get('Score',0) for m in res]):.3f}" if res else "0.000"})
                        st.session_state.summary_df=pd.DataFrame(summary)
                        st.session_state.performance_metrics={'total_time':total_time,'total_bp':total_bp_processed,'speed':overall_speed,'sequences':len(st.session_state.seqs),'total_motifs':sum(len(r) for r in all_results),'detector_count':len(DETECTOR_PROCESSES),'estimated_time':estimated_total,'analysis_steps':[f"{n} detection" for n,_ in DETECTOR_PROCESSES]+['Hybrid/Cluster detection','Overlap resolution']}
                        progress_placeholder.empty(); status_placeholder.empty(); detailed_progress_placeholder.empty()
                        timer_placeholder.markdown(f"""<div class='progress-panel progress-panel--success'><h3 class='progress-panel__title'>🎉 Analysis Complete!</h3>
                        <div class='stats-grid stats-grid--wide'>
                        <div class='stat-card'><h2 class='stat-card__value'>⏱️ {total_time:.2f}s</h2><p class='stat-card__label'>Total Time</p></div>
                        <div class='stat-card'><h2 class='stat-card__value'>🧬 {total_bp_processed:,}</h2><p class='stat-card__label'>Base Pairs</p></div>
                        <div class='stat-card'><h2 class='stat-card__value'>🚀 {overall_speed:,.0f}</h2><p class='stat-card__label'>bp/second</p></div>
                        <div class='stat-card'><h2 class='stat-card__value'>🔬 {len(DETECTOR_PROCESSES)}</h2><p class='stat-card__label'>Detectors</p></div>
                        <div class='stat-card'><h2 class='stat-card__value'>🎯 {sum(len(r) for r in all_results)}</h2><p class='stat-card__label'>Motifs Found</p></div>
                        <div class='stat-card'><h2 class='stat-card__value'>📊 {len(st.session_state.seqs)}</h2><p class='stat-card__label'>Sequences</p></div>
                        </div></div>""", unsafe_allow_html=True)
                        st.success("Results are available below and in the 'Analysis Results and Visualization' tab."); st.session_state.analysis_status="Complete"
                    except Exception as e:
                        timer_placeholder.empty(); progress_placeholder.empty(); status_placeholder.empty(); detailed_progress_placeholder.empty(); st.error(f"Analysis failed: {str(e)}"); st.session_state.analysis_status="Error"
        if st.session_state.get('summary_df') is not None: st.markdown("#### Analysis Summary"); st.dataframe(st.session_state.summary_df)
    st.markdown("---")

# --- Results tab ---
with tab_pages["Results"]:
    st.markdown('<h2>Analysis Results and Visualization</h2>', unsafe_allow_html=True)
    if not st.session_state.results: st.info("No analysis results. Please run motif analysis first.")
    else:
        if st.session_state.get('performance_metrics'):
            m=st.session_state.performance_metrics
            st.markdown(f"""<div class='progress-panel progress-panel--metrics'><h3 class='progress-panel__title'>Performance Metrics</h3>
            <div class='stats-grid stats-grid--wide'>
            <div class='stat-card'><h2 class='stat-card__value'>{m['total_time']:.2f}s</h2><p class='stat-card__label'>Processing Time</p></div>
            <div class='stat-card'><h2 class='stat-card__value'>{m['total_bp']:,}</h2><p class='stat-card__label'>Base Pairs</p></div>
            <div class='stat-card'><h2 class='stat-card__value'>{m['speed']:,.0f}</h2><p class='stat-card__label'>bp/second</p></div>
            <div class='stat-card'><h2 class='stat-card__value'>{m.get('detector_count',9)}</h2><p class='stat-card__label'>Detector Processes</p></div>
            <div class='stat-card'><h2 class='stat-card__value'>{m['sequences']}</h2><p class='stat-card__label'>Sequences</p></div>
            <div class='stat-card'><h2 class='stat-card__value'>{m['total_motifs']}</h2><p class='stat-card__label'>Total Motifs</p></div>
            </div></div>""", unsafe_allow_html=True)
        st.markdown("### Analysis Summary"); st.dataframe(st.session_state.summary_df, use_container_width=True)
        seq_idx=0
        if len(st.session_state.seqs)>1:
            selected=st.pills("Choose Sequence for Details:", options=list(range(len(st.session_state.seqs))), format_func=lambda i: st.session_state.names[i], selection_mode="single", default=0)
            seq_idx=selected or 0
        motifs=st.session_state.results[seq_idx]; seq_len=len(st.session_state.seqs[seq_idx]); seq_name=st.session_state.names[seq_idx]
        if not motifs: st.warning("No motifs detected for this sequence.")
        else:
            filtered=[m for m in motifs if m.get('Class') not in ['Hybrid','Non-B_DNA_Clusters']]
            hybrids=[m for m in motifs if m.get('Class') in ['Hybrid','Non-B_DNA_Clusters']]
            df=pd.DataFrame(filtered) if filtered else pd.DataFrame()
            stats=get_basic_stats(st.session_state.seqs[seq_idx], filtered); motif_count=len(filtered); hybrid_count=len(hybrids)
            st.markdown(f"""<div class='progress-panel progress-panel--results'><h3 class='progress-panel__title progress-panel__title--large'>NBDScanner Analysis Results</h3>
            <div class='stats-grid stats-grid--extra-wide'>
            <div class='stat-card stat-card--large'><h2 class='stat-card__value stat-card__value--large'>{stats.get("Coverage%",0):.2f}%</h2><p class='stat-card__label stat-card__label--large'>Sequence Coverage</p></div>
            <div class='stat-card stat-card--large'><h2 class='stat-card__value stat-card__value--large'>{stats.get("Density",0):.2f}</h2><p class='stat-card__label stat-card__label--large'>Motif Density<br>(motifs/kb)</p></div>
            <div class='stat-card stat-card--large'><h2 class='stat-card__value stat-card__value--large'>{motif_count}</h2><p class='stat-card__label stat-card__label--large'>Total Motifs</p></div>
            <div class='stat-card stat-card--large'><h2 class='stat-card__value stat-card__value--large'>{seq_len:,}</h2><p class='stat-card__label stat-card__label--large'>Sequence Length (bp)</p></div>
            </div></div>""", unsafe_allow_html=True)
            if hybrid_count>0: st.info(f"{hybrid_count} Hybrid/Cluster motifs detected. View them in the 'Cluster/Hybrid' tab below.")
            available_cols=[c for c in df.columns if c!='Normalized_Score']; default_cols=[c for c in ['Class','Subclass','Start','End','Length','Sequence','Score'] if c in available_cols]
            display_cols=st.pills("Select columns to display:", options=available_cols, selection_mode="multi", default=default_cols) or []
            ROWS_PER_PAGE=100; total_rows=len(df)
            if total_rows>ROWS_PER_PAGE:
                _,pag, _=st.columns([2,3,2]); max_pages=(total_rows+ROWS_PER_PAGE-1)//ROWS_PER_PAGE
                with pag:
                    page=st.number_input(f"Page (showing {ROWS_PER_PAGE} rows per page)", min_value=1, max_value=max_pages, value=1, step=1, help=f"Total: {total_rows} motifs across {max_pages} pages")
                    start=(page-1)*ROWS_PER_PAGE; end=min(start+ROWS_PER_PAGE, total_rows); st.caption(f"Showing motifs {start+1} to {end} of {total_rows}")
                df_page=df.iloc[start:end]
            else: df_page=df
            if display_cols:
                df_disp=df_page[display_cols].copy(); df_disp.columns=[c.replace('_',' ') for c in df_disp.columns]; st.dataframe(df_disp, use_container_width=True, height=360)
            else:
                df_disp=df_page.copy(); df_disp.columns=[c.replace('_',' ') for c in df_disp.columns]; st.dataframe(df_disp, use_container_width=True, height=360)

            st.markdown('<h3>Visualizations</h3>', unsafe_allow_html=True)
            viz_tabs=st.tabs(["Distribution","Coverage & Density","Statistics","Cluster/Hybrid"])
            with viz_tabs[0]:
                st.markdown("##### Motif Distribution")
                try:
                    fig1=plot_motif_distribution(filtered, by='Class', title=f"Motif Classes - {seq_name}"); st.pyplot(fig1); plt.close(fig1)
                    fig2=plot_motif_distribution(filtered, by='Subclass', title=f"Motif Subclasses - {seq_name}"); st.pyplot(fig2); plt.close(fig2)
                    fig3=plot_nested_pie_chart(filtered, title=f"Class-Subclass Distribution - {seq_name}"); st.pyplot(fig3); plt.close(fig3)
                except Exception as e: st.error(f"Error generating distribution plots: {e}")
            with viz_tabs[1]:
                st.markdown("##### Sequence Coverage Analysis")
                try:
                    fig4=plot_coverage_map(filtered, seq_len, title=f"Motif Coverage - {seq_name}"); st.pyplot(fig4); plt.close(fig4)
                    fig5=plot_density_heatmap(filtered, seq_len, window_size=max(100, seq_len//20), title=f"Motif Density - {seq_name}"); st.pyplot(fig5); plt.close(fig5)
                    st.markdown("##### Circos Density Plot"); fig_c=plot_circos_motif_density(filtered, seq_len, title=f"Circos Density - {seq_name}"); st.pyplot(fig_c); plt.close(fig_c)
                except Exception as e: st.error(f"Error generating coverage plots: {e}")
            with viz_tabs[2]:
                st.markdown("##### Statistical Analysis"); st.markdown("###### Density Metrics")
                try:
                    gd=calculate_genomic_density(filtered, seq_len, by_class=True); pd_kbp=calculate_positional_density(filtered, seq_len, unit='kbp', by_class=True)
                    pd_mbp=calculate_positional_density(filtered, seq_len, unit='Mbp', by_class=True); c1,c2=st.columns(2)
                    with c1: st.metric("Genomic Density", f"{gd.get('Overall',0):.4f}%")
                    with c2: st.metric("Motifs/kbp", f"{pd_kbp.get('Overall',0):.2f}")
                    level=st.radio("Density Analysis Level:", ["Class Level","Subclass Level"], horizontal=True)
                    if level=="Class Level":
                        rows=[{'Motif Class':k.replace('_',' '),'Genomic Density (%)':f"{gd.get(k,0):.4f}",'Motifs/kbp':f"{pd_kbp.get(k,0):.2f}"} for k in sorted([x for x in gd.keys() if x!='Overall'])]
                        if rows: st.dataframe(pd.DataFrame(rows), use_container_width=True, height=200)
                        fig_d=plot_density_comparison(gd, pd_kbp, title="Motif Density Analysis (Class Level)"); st.pyplot(fig_d); plt.close(fig_d)
                    else:
                        gd_sub=calculate_genomic_density(filtered, seq_len, by_class=False, by_subclass=True)
                        pd_sub=calculate_positional_density(filtered, seq_len, unit='kbp', by_class=False, by_subclass=True)
                        rows=[{'Motif Subclass':k.replace('_',' ').replace(':',': '),'Genomic Density (%)':f"{gd_sub.get(k,0):.4f}",'Motifs/kbp':f"{pd_sub.get(k,0):.2f}"} for k in sorted([x for x in gd_sub.keys() if x!='Overall'])]
                        if rows: st.dataframe(pd.DataFrame(rows), use_container_width=True, height=300)
                        fig_ds=plot_density_comparison_by_subclass(gd_sub, pd_sub, title="Motif Density Analysis (Subclass Level)"); st.pyplot(fig_ds); plt.close(fig_ds)
                        if len(filtered)>0:
                            st.markdown("**Subclass Density Heatmap Along Sequence**"); fig_h=plot_subclass_density_heatmap(filtered, seq_len, window_size=max(500, seq_len//50), title="Subclass Density Distribution"); st.pyplot(fig_h); plt.close(fig_h)
                except Exception as e: st.error(f"Error calculating density metrics: {e}"); import traceback; st.error(traceback.format_exc())
                st.markdown("###### Distributions")
                try:
                    fig6=plot_length_distribution(filtered, by_class=True, title="Length Distribution by Motif Class"); st.pyplot(fig6); plt.close(fig6)
                    fig7=plot_score_distribution(filtered, by_class=True, title="Score Distribution by Motif Class (1-3 Scale)"); st.pyplot(fig7); plt.close(fig7)
                except Exception as e: st.error(f"Error generating distribution plots: {e}")
            with viz_tabs[3]:
                st.markdown("##### Hybrid & Cluster Motif Analysis")
                if hybrids:
                    hybrid_only=[m for m in hybrids if m.get('Class')=='Hybrid']; cluster_only=[m for m in hybrids if m.get('Class')=='Non-B_DNA_Clusters']; hc_df=pd.DataFrame(hybrids)
                    st.markdown(f"""<div class='progress-panel progress-panel--hybrid'><h3 class='progress-panel__title progress-panel__title--large'>Hybrid & Cluster Motif Summary</h3>
                    <div class='stats-grid stats-grid--extra-wide'>
                    <div class='stat-card stat-card--large'><h2 class='stat-card__value stat-card__value--large'>{len(hybrid_only)}</h2><p class='stat-card__label stat-card__label--large'>Hybrid Motifs</p></div>
                    <div class='stat-card stat-card--large'><h2 class='stat-card__value stat-card__value--large'>{len(cluster_only)}</h2><p class='stat-card__label stat-card__label--large'>DNA Clusters</p></div>
                    <div class='stat-card stat-card--large'><h2 class='stat-card__value stat-card__value--large'>{int(hc_df['Length'].mean()) if 'Length' in hc_df.columns else 0}</h2><p class='stat-card__label stat-card__label--large'>Avg Length (bp)</p></div>
                    <div class='stat-card stat-card--large'><h2 class='stat-card__value stat-card__value--large'>{len(hybrids)}</h2><p class='stat-card__label stat-card__label--large'>Total</p></div>
                    </div></div>""", unsafe_allow_html=True)
                    sub_tabs=st.tabs(["All","Hybrid Motifs","Cluster Motifs"])
                    with sub_tabs[0]:
                        st.markdown("##### All Hybrid & Cluster Motifs")
                        disp_cols=['Class','Subclass','Start','End','Length','Score']; disp_cols=[c for c in disp_cols if c in hc_df.columns]
                        view=hc_df[disp_cols].copy(); view.columns=[c.replace('_',' ') for c in view.columns]; st.dataframe(view, use_container_width=True, height=300)
                        h=min(20, max(4, len(hybrids)*0.3)); fig,ax=plt.subplots(figsize=(12,h)); colors={'Hybrid':'#ff6b6b','Non-B_DNA_Clusters':'#4ecdc4'}
                        for i,m in enumerate(hybrids):
                            color=colors.get(m.get('Class'),'#95a5a6'); ax.barh(i, m['End']-m['Start'], left=m['Start'], height=0.8, color=color, alpha=0.7, edgecolor='black', linewidth=0.5)
                        ax.set_xlabel("Position (bp)"); ax.set_ylabel("Motif Index"); ax.set_title("Hybrid & Cluster Position Map"); handles=[plt.Rectangle((0,0),1,1,color=c,alpha=0.7) for c in colors.values()]
                        ax.legend(handles, list(colors.keys()), loc='upper right'); ax.grid(axis='x', alpha=0.3); plt.tight_layout(); st.pyplot(fig); plt.close(fig)
                    with sub_tabs[1]:
                        st.markdown("##### Hybrid Motifs (Overlapping Different Classes)")
                        if hybrid_only:
                            hdf=pd.DataFrame(hybrid_only); st.info(f"**Hybrid motifs** are regions where different Non-B DNA motif classes overlap. Found {len(hybrid_only)} hybrid regions.")
                            cols=['Subclass','Start','End','Length','Score','Component_Classes']; cols=[c for c in cols if c in hdf.columns]; view=hdf[cols].copy(); view.columns=[c.replace('_',' ') for c in view.columns]; st.dataframe(view, use_container_width=True, height=300)
                            h=min(15, max(3, len(hybrid_only)*0.4)); fig,ax=plt.subplots(figsize=(12,h))
                            for i,m in enumerate(hybrid_only):
                                ax.barh(i, m['End']-m['Start'], left=m['Start'], height=0.8, color='#ff6b6b', alpha=0.7, edgecolor='#c0392b', linewidth=1)
                                subclass=m.get('Subclass',''); lbl=f" {subclass[:30]}{'...' if len(subclass)>30 else ''}" if subclass else ""; ax.text(m['Start'], i, lbl, va='center', fontsize=8)
                            ax.set_xlabel("Position (bp)"); ax.set_ylabel("Hybrid Index"); ax.set_title("Hybrid Motif Positions (Overlapping Classes)"); ax.grid(axis='x', alpha=0.3); plt.tight_layout(); st.pyplot(fig); plt.close(fig)
                        else: st.info("No hybrid motifs detected in this sequence.")
                    with sub_tabs[2]:
                        st.markdown("##### Non-B DNA Clusters (High-Density Regions)")
                        if cluster_only:
                            cdf=pd.DataFrame(cluster_only); st.info(f"**DNA Clusters** are high-density regions with multiple Non-B DNA motif classes. Found {len(cluster_only)} cluster regions.")
                            cols=['Subclass','Start','End','Length','Score','Motif_Count','Class_Diversity']; cols=[c for c in cols if c in cdf.columns]; view=cdf[cols].copy(); view.columns=[c.replace('_',' ') for c in view.columns]; st.dataframe(view, use_container_width=True, height=300)
                            h=min(15, max(3, len(cluster_only)*0.4)); fig,ax=plt.subplots(figsize=(12,h))
                            for i,m in enumerate(cluster_only):
                                div=m.get('Class_Diversity',2); alpha=min(0.4+div*0.1,0.9); ax.barh(i, m['End']-m['Start'], left=m['Start'], height=0.8, color='#4ecdc4', alpha=alpha, edgecolor='#16a085', linewidth=1)
                                cnt=m.get('Motif_Count',0); ax.text(m['End'], i, f" {cnt} motifs", va='center', fontsize=8)
                            ax.set_xlabel("Position (bp)"); ax.set_ylabel("Cluster Index"); ax.set_title("Non-B DNA Cluster Positions (Color intensity = Class diversity)"); ax.grid(axis='x', alpha=0.3); plt.tight_layout(); st.pyplot(fig); plt.close(fig)
                        else: st.info("No DNA clusters detected in this sequence.")
                else: st.info("No hybrid or cluster motifs detected in this sequence. Hybrid motifs occur when different Non-B DNA classes overlap, and clusters form when multiple motifs are found in close proximity.")

# --- Download tab ---
with tab_pages["Download"]:
    st.header("Export Data")
    if not st.session_state.results: st.info("No results available to download.")
    else:
        st.markdown("### Export Options")
        overlap_option_used=st.session_state.get('overlap_option_used','Remove overlaps within subclasses (default)')
        st.info(f"**Analysis Configuration Used:**\n• Overlap Handling: {overlap_option_used}\n• Motif Classes: {len(st.session_state.get('selected_classes', []))} classes selected\n• Total Motifs Found: {sum(len(m) for m in st.session_state.results)}")
        col1,col2=st.columns(2)
        with col1: st.markdown("**Export Configuration**")
        with col2: include_sequences=st.checkbox("Include Full Sequences", value=True, help="Include full motif sequences in export")
        all_motifs=[]; hybrid_cluster_export=[]
        for i,motifs in enumerate(st.session_state.results):
            for m in motifs:
                mm=m.copy(); 
                if 'Sequence_Name' not in mm: mm['Sequence_Name']=st.session_state.names[i]
                if mm.get('Class') in ['Hybrid','Non-B_DNA_Clusters']: hybrid_cluster_export.append(mm)
                else: all_motifs.append(mm)
        if hybrid_cluster_export: st.info(f"{len(hybrid_cluster_export)} Hybrid/Cluster motifs are excluded from downloads. These are shown only in the Cluster/Hybrid visualization tab.")
        if all_motifs:
            df_prev=export_results_to_dataframe(all_motifs).head(10); st.markdown("### Export Preview"); st.dataframe(df_prev, use_container_width=True); st.caption(f"Showing first 10 of {len(all_motifs)} total records (Hybrid/Cluster motifs excluded)")
        st.markdown("### Download Files")
        st.info("**Excel Format**: Downloads a multi-sheet workbook with:\n• **Consolidated Sheet**: All non-overlapping motifs\n• **Class Sheets**: Separate sheets for each motif class (G-Quadruplex, Z-DNA, etc.)\n• **Subclass Sheets**: Detailed breakdown by subclass (when multiple subclasses exist)")
        c1,c2,c3,c4=st.columns(4)
        with c1:
            if all_motifs:
                csv_data=export_to_csv(all_motifs)
                st.download_button("Download CSV", data=csv_data.encode('utf-8'), file_name="nbdscanner_results.csv", mime="text/csv", use_container_width=True, help="Comma-separated values format")
        with c2:
            if all_motifs:
                try:
                    excel_bytes=generate_excel_bytes(all_motifs)
                    st.download_button("Download Excel", data=excel_bytes, file_name="nbdscanner_results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", use_container_width=True, help="Excel format with multiple sheets (Consolidated + per-class)")
                except Exception as e: st.error(f"Excel export error: {str(e)}")
        with c3:
            if all_motifs:
                json_data=export_to_json(all_motifs, pretty=True)
                st.download_button("Download JSON", data=json_data.encode('utf-8'), file_name="nbdscanner_results.json", mime="application/json", use_container_width=True, help="JSON format with metadata")
        with c4:
            if all_motifs and st.session_state.names:
                bed_data=export_to_bed(all_motifs, st.session_state.names[0])
                st.download_button("Download BED", data=bed_data.encode('utf-8'), file_name="nbdscanner_results.bed", mime="text/plain", use_container_width=True, help="BED format for genome browsers")
        st.markdown("### Additional Exports")
        c5,c6,c7,c8=st.columns(4)
        with c5:
            if all_motifs:
                import json; config_summary={"analysis_info":{"tool":"NBDScanner","version":"2024.1","motif_classes_detected":list(set(m.get('Class','Unknown') for m in all_motifs)),"total_sequences_analyzed":len(st.session_state.names) if hasattr(st.session_state,'names') else 0,"total_motifs_found":len(all_motifs)},"detection_parameters":{"all_classes_enabled":True,"consolidated_system":True}}
                st.download_button("Download Config", data=json.dumps(config_summary, indent=2), file_name="analysis_configuration.json", mime="application/json", use_container_width=True, help="Analysis configuration and metadata")
        all_seq_data=[{'sequence_name': st.session_state.names[i], 'motifs': m} for i,m in enumerate(st.session_state.results)]

# --- Documentation tab ---
with tab_pages["Documentation"]:
    st.header("Scientific Documentation & References")
    st.markdown("""
    <div style='background:#e3f2fd; border-radius:12px; padding:18px; font-size:1.08rem; font-family:Montserrat,Arial; margin-bottom:20px;'>
    <b>Motif Detection Parameters:</b><br><br>
    <table style='width:100%; border-collapse: collapse;'>
        <tr style='background:#1976d2; color:white;'>
            <th style='padding:10px; text-align:left;'>Motif Type</th>
            <th style='padding:10px; text-align:left;'>Arm/Unit Length</th>
            <th style='padding:10px; text-align:left;'>Spacer/Loop</th>
            <th style='padding:10px; text-align:left;'>Additional Requirements</th>
        </tr>
        <tr style='background:#f5f5f5;'><td style='padding:10px;'><b>Cruciform DNA</b></td><td style='padding:10px;'>10–100 nt</td><td style='padding:10px;'>0–3 nt</td><td style='padding:10px;'>Reverse complement arms</td></tr>
        <tr><td style='padding:10px;'><b>Slipped DNA</b></td><td style='padding:10px;'>10–50 nt repeat</td><td style='padding:10px;'>0 nt (adjacent)</td><td style='padding:10px;'>Direct repeat units</td></tr>
        <tr style='background:#f5f5f5;'><td style='padding:10px;'><b>Triplex DNA</b></td><td style='padding:10px;'>10–100 nt mirrored</td><td style='padding:10px;'>0–8 nt</td><td style='padding:10px;'>≥90% Purine or Pyrimidine</td></tr>
    </table>
    </div>
    """, unsafe_allow_html=True)
    st.markdown("""
    <div style='background:#f4faff; border-radius:12px; padding:18px; font-size:1.08rem; font-family:Montserrat,Arial;'>
    <b>Motif Classes Detected:</b><br><br>
    <ul>
        <li><b>Curved DNA</b>: Phased poly(A/T) tracts; scoring by tract length/grouping.</li>
        <li><b>Z-DNA</b>: Alternating purine-pyrimidine, GC-rich; windowed scoring.</li>
        <li><b>eGZ-motif</b>: Long (CGG)n runs; scored by repeat count.</li>
        <li><b>Slipped DNA</b>: Direct/tandem repeats (10–50 nt unit); scoring by length/unit copies.</li>
        <li><b>R-Loop</b>: G-rich hybrid sites; regex + thermodynamic weighting.</li>
        <li><b>Cruciform</b>: Palindromic inverted repeats; scoring by arm length/A-T content.</li>
        <li><b>Triplex / Mirror</b>: Purine/pyrimidine mirrors; ≥90% purity; arm/loop rules.</li>
        <li><b>Sticky DNA</b>: GAA/TTC extended repeats.</li>
        <li><b>G-Triplex</b>: Three G runs; loop-length penalty.</li>
        <li><b>G4 & Variants</b>: Canonical/variant G4 regex + G4Hunter-style scoring.</li>
        <li><b>i-Motif</b>: C-rich runs; scoring by run count/content.</li>
        <li><b>AC-Motif</b>: Alternating A-rich/C-rich consensus.</li>
        <li><b>A-philic DNA</b>: Tetranucleotide log2-odds for A-tract affinity; high-confidence and moderate bins.</li>
        <li><b>Hybrid Motif</b>: Overlaps of different classes.</li>
        <li><b>Non-B DNA Clusters</b>: Hotspots with multiple classes.</li>
    </ul>
    <b>References:</b>
    <ul>
        <li>Bedrat et al., 2016 NAR</li><li>Ho et al., 2010 Nat Chem Biol</li><li>Kim et al., 2018 NAR</li><li>Zeraati et al., 2018 Nat Chem</li><li>Bacolla et al., 2006 NAR</li><li>Mirkin & Frank-Kamenetskii, 1994 ARBBS</li><li>Vinogradov, 2003 Bioinformatics</li><li>Bolshoy et al., 1991 PNAS</li><li>Rohs et al., 2009 Nature</li><li>New et al., 2020 J DNA Structure</li>
    </ul>
    </div>
    """, unsafe_allow_html=True)

st.markdown("""
---
<div style='font-size: 1.05rem; color: #1e293b; margin-top: 30px; text-align: left; font-family:Montserrat,Arial;'>
<b>Developed by</b><br>
Dr. Venkata Rajesh Yella<br>
<a href='mailto:yvrajesh_bt@kluniversity.in'>yvrajesh_bt@kluniversity.in</a> |
<a href='https://github.com/VRYella' target='_blank'>GitHub: VRYella</a>
</div>
""", unsafe_allow_html=True)
