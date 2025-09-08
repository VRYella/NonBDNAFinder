---
title: Non-B DNA Motif Finder
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: streamlit
sdk_version: "1.38.0"
app_file: app.py
pinned: false
---

# Non-B DNA Motif Finder

This tool identifies and visualizes diverse non-B DNA structures such as:
- A-Phased Repeats
- Z-DNA
- Slipped DNA
- R-Loops
- Cruciforms
- H-DNA
- Sticky DNA
- G-Triplexes
- G-Quadruplex variants (canonical, relaxed, bulged, bipartite, multimeric)
- i-Motifs
- Hybrids
- Motif hotspots

## Features
- Sequence input (FASTA or plain text)
- Regex-based motif detection
- Scoring models like NBST, G4Hunter, and Z-Seeker
- Visualization of motif locations
- Hotspot identification
- Export in BED, GFF, or FASTA format

## How to Run Locally
```bash
pip install -r requirements.txt
streamlit run app.py
