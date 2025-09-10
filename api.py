#!/usr/bin/env python3
"""
NBDFinder REST API
==================

Production-ready REST API endpoints for programmatic access to NBDFinder motif detection.
Provides JSON responses with comprehensive motif analysis results.

Usage:
    python api.py  # Start FastAPI server on localhost:8000
    
Endpoints:
    GET /api/v1/health - Health check
    POST /api/v1/analyze - Analyze DNA sequence for motifs
    GET /api/v1/stats - Get analysis statistics
    GET /api/v1/motif-info - Get supported motif class information
"""

from fastapi import FastAPI, HTTPException, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from typing import List, Dict, Any, Optional
import json
import time
import hashlib
from datetime import datetime
import uvicorn

# Import NBDFinder modules
from all_motifs_refactored import all_motifs_refactored
from utils import get_basic_stats, parse_fasta

# Import configuration if available
try:
    from classification_config import MOTIF_LENGTH_LIMITS, SCORING_METHODS
except ImportError:
    # Fallback configuration based on NBDFinder specifications
    MOTIF_LENGTH_LIMITS = {
        "Curved_DNA": [10, 200],
        "A_Phiic_DNA": [15, 500],
        "Z_DNA": [6, 100], 
        "Slipped_DNA": [6, 500],
        "R_Loop": [12, 500],
        "Cruciform_DNA": [10, 200],
        "Triplex_DNA": [10, 200],
        "G_Quadruplex": [15, 50],
        "i_Motif": [12, 40],
        "AC_Motif": [10, 100],
        "Hybrid_Motif": [10, 500],
        "Non_B_DNA_Cluster": [50, 2000]
    }
    
    SCORING_METHODS = {
        "Curved_DNA": "Curvature-based geometric scoring",
        "A_Phiic_DNA": "Tetranucleotide propensity scoring", 
        "Z_DNA": "GC-content and alternating pattern scoring",
        "Slipped_DNA": "Repeat length and copy number scoring",
        "R_Loop": "Thermodynamic hybrid stability scoring",
        "Cruciform_DNA": "Palindrome arm length and AT-content scoring",
        "Triplex_DNA": "Purine/pyrimidine purity scoring",
        "G_Quadruplex": "G4Hunter thermodynamic scoring",
        "i_Motif": "C-run count and pH-dependent scoring",
        "AC_Motif": "Alternating pattern recognition scoring",
        "Hybrid_Motif": "Multi-motif overlap scoring",
        "Non_B_DNA_Cluster": "Motif density and diversity scoring"
    }

# API Models
class SequenceRequest(BaseModel):
    """Request model for sequence analysis"""
    sequence: str = Field(..., description="DNA sequence to analyze (FASTA or plain text)")
    sequence_name: Optional[str] = Field("API_Sequence", description="Name for the sequence")
    nonoverlap: Optional[bool] = Field(False, description="Remove overlapping motifs")
    report_hotspots: Optional[bool] = Field(True, description="Detect motif hotspots")
    calculate_conservation: Optional[bool] = Field(False, description="Calculate conservation analysis")
    score_threshold: Optional[float] = Field(0.0, description="Minimum normalized score threshold")

class MotifResult(BaseModel):
    """Response model for individual motif"""
    s_no: int
    sequence_name: str
    chromosome: str
    motif_class: str
    subclass: str
    motif_id: str
    start: int
    end: int
    length: int
    normalized_score: float
    actual_score: float
    scoring_method: str
    gc_content: float
    sequence: str
    overlap_classes: str

class AnalysisResponse(BaseModel):
    """Response model for sequence analysis"""
    success: bool
    analysis_id: str
    timestamp: str
    sequence_info: Dict[str, Any]
    motifs: List[Dict[str, Any]]
    statistics: Dict[str, Any]
    processing_time: float
    
class HealthResponse(BaseModel):
    """Health check response"""
    status: str
    version: str
    timestamp: str
    motif_classes_available: int

class StatsResponse(BaseModel):
    """Statistics response"""
    total_analyses: int
    motif_classes_supported: List[str]
    scoring_methods: Dict[str, str]
    length_limits: Dict[str, List[int]]

# Initialize FastAPI app
app = FastAPI(
    title="NBDFinder API",
    description="REST API for Non-B DNA motif detection and analysis",
    version="1.0.0",
    contact={
        "name": "Dr. Venkata Rajesh Yella",
        "email": "yvrajesh_bt@kluniversity.in",
        "url": "https://github.com/VRYella/NBDFinder"
    }
)

# Add CORS middleware for web browser access
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Global statistics tracking
analysis_counter = 0
analysis_cache = {}

@app.get("/api/v1/health", response_model=HealthResponse)
async def health_check():
    """Health check endpoint"""
    return HealthResponse(
        status="healthy",
        version="1.0.0",
        timestamp=datetime.now().isoformat(),
        motif_classes_available=len(MOTIF_LENGTH_LIMITS)
    )

@app.get("/api/v1/stats", response_model=StatsResponse)
async def get_stats():
    """Get API usage statistics and configuration"""
    return StatsResponse(
        total_analyses=analysis_counter,
        motif_classes_supported=list(MOTIF_LENGTH_LIMITS.keys()),
        scoring_methods=SCORING_METHODS,
        length_limits=MOTIF_LENGTH_LIMITS
    )

@app.get("/api/v1/motif-info")
async def get_motif_info():
    """Get detailed information about supported motif classes with scientific descriptions"""
    motif_info = {
        "classes": len(MOTIF_LENGTH_LIMITS),
        "details": {},
        "api_version": "1.0.0",
        "last_updated": "2024",
        "documentation": "Comprehensive motif class information with scientific literature references"
    }
    
    for motif_class, limits in MOTIF_LENGTH_LIMITS.items():
        description_data = get_motif_description(motif_class)
        motif_info["details"][motif_class] = {
            "min_length": limits[0],
            "max_length": limits[1],
            "scoring_method": SCORING_METHODS.get(motif_class, "Pattern-based"),
            "description": description_data.get("description", "Non-B DNA motif"),
            "mechanism": description_data.get("mechanism", "Alternative DNA structure formation"),
            "biological_significance": description_data.get("biological_significance", "Potential regulatory role"),
            "detection_method": description_data.get("detection_method", "Computational pattern matching"),
            "literature_references": description_data.get("literature_references", [])
        }
    
    return motif_info

def get_motif_description(motif_class: str) -> Dict[str, Any]:
    """Get comprehensive description for motif class with scientific details"""
    descriptions = {
        "Curved_DNA": {
            "description": "DNA sequences with intrinsic curvature formed by phased A-tracts or T-tracts. These sequences exhibit significant bending (up to 18° per A-tract) due to altered backbone geometry and groove widths.",
            "mechanism": "Phased A/T-tracts spaced at 10-11 bp intervals create cumulative DNA bending. The narrow minor groove and altered sugar pucker contribute to the curvature.",
            "biological_significance": "Involved in nucleosome positioning, transcriptional regulation, and protein-DNA recognition. Common in promoter regions and origins of replication.",
            "detection_method": "Regex-based identification of phased poly(A) or poly(T) tracts with specific spacing rules. Scoring based on tract length and periodicity.",
            "literature_references": [
                "Crothers et al. (1990) J Biol Chem 265:7093-7096",
                "Koo et al. (1986) Nature 320:501-506", 
                "Trifonov & Sussman (1980) Proc Natl Acad Sci USA 77:3816-3820"
            ]
        },
        "A_Phiic_DNA": {
            "description": "DNA regions with high affinity for A-tract formation, characterized by specific tetranucleotide and trinucleotide propensity patterns that favor A-philic conformations.",
            "mechanism": "Tetranucleotide and trinucleotide log2-odds scoring identifies regions with A-tract propensity. Window merging creates maximal non-overlapping A-philic regions.",
            "biological_significance": "Associated with nucleosome positioning signals and chromatin organization. Influences DNA accessibility and protein binding.",
            "detection_method": "Statistical scoring using tetranucleotide/trinucleotide propensity matrices with window-based detection and merging algorithms.",
            "literature_references": [
                "Segal et al. (2006) Nature 442:772-778",
                "Kaplan et al. (2009) Nature 458:362-366"
            ]
        },
        "Z_DNA": {
            "description": "Left-handed double helix formed by alternating purine-pyrimidine sequences, particularly (CG)n and (CA/TG)n repeats. Z-DNA has a zigzag backbone and altered major/minor groove geometry.",
            "mechanism": "High GC content and alternating purine-pyrimidine sequences promote left-handed helix formation. Negative supercoiling and methylation can stabilize Z-DNA.",
            "biological_significance": "Implicated in transcriptional regulation, chromatin structure, and genomic instability. Z-DNA binding proteins regulate gene expression.",
            "detection_method": "Combined approach using regex for dinucleotide repeats and windowed GC-content scoring. Includes detection of eGZ (extended G-rich Z-DNA) motifs.",
            "literature_references": [
                "Rich & Zhang (2003) Nat Rev Genet 4:566-572",
                "Wang et al. (1979) Nature 282:680-686",
                "Ho et al. (2010) Nat Chem Biol 6:261-267"
            ]
        },
        "Slipped_DNA": {
            "description": "Secondary structures formed by direct tandem repeats that can loop out during replication, creating slipped-strand intermediates with potential for expansion/contraction.",
            "mechanism": "Repetitive sequences form hairpin loops during DNA replication or repair, facilitated by complementary base pairing between repeat units.",
            "biological_significance": "Associated with microsatellite instability, trinucleotide repeat expansion diseases, and evolutionary genome dynamics.",
            "detection_method": "Pattern matching for direct and tandem repeats using repeat-unit identification. Scoring based on repeat length and copy number.",
            "literature_references": [
                "Pearson & Sinden (1996) Biochemistry 35:5041-5053",
                "Wells et al. (2005) J Biol Chem 280:2743-2746",
                "Bacolla & Wells (2004) J Biol Chem 279:47411-47414"
            ]
        },
        "R_Loop": {
            "description": "Three-stranded nucleic acid structures where RNA hybridizes to complementary DNA, displacing the non-template DNA strand as a single-stranded loop.",
            "mechanism": "G-rich DNA regions facilitate stable RNA-DNA hybrid formation through enhanced thermodynamic stability compared to DNA-DNA duplexes.",
            "biological_significance": "Play roles in transcription termination, DNA repair, genome instability, and regulation of gene expression.",
            "detection_method": "Detection of G-rich regions using RLFS (R-Loop Forming Sequence) models and regex patterns. Thermodynamic scoring for hybrid stability assessment.",
            "literature_references": [
                "Ginno et al. (2012) Mol Cell 45:814-825",
                "Skourti-Stathaki et al. (2011) Mol Cell 42:794-805",
                "Santos-Pereira & Aguilera (2015) Nat Rev Genet 16:583-597"
            ]
        },
        "Cruciform_DNA": {
            "description": "Four-way junction structures formed by palindromic inverted repeat sequences. These structures resemble Holliday junctions and can adopt various conformations.",
            "mechanism": "Palindromic sequences with appropriate spacer regions can extrude to form cruciform structures, particularly under negative supercoiling conditions.",
            "biological_significance": "Potential sites for recombination, DNA damage, and genomic instability. May play roles in replication fork stalling and chromosome rearrangements.",
            "detection_method": "Identification of palindromic inverted repeats with spacers using regex and reverse complement matching. Scoring based on arm length and A/T content.",
            "literature_references": [
                "Lilley & Kemper (1984) Cell 36:413-422",
                "Pearson et al. (1996) Nucleic Acids Res 24:1399-1405",
                "Bacolla et al. (2006) Nucleic Acids Res 34:2898-2905"
            ]
        },
        "Triplex_DNA": {
            "description": "Three-stranded DNA structures formed when a pyrimidine or purine-rich single strand binds in the major groove of a duplex DNA through Hoogsteen base pairing.",
            "mechanism": "Triplex formation requires purine/pyrimidine mirror repeats or specific triplex motifs. Hoogsteen hydrogen bonding creates the third strand interaction.",
            "biological_significance": "Involved in gene regulation, DNA damage, and potential therapeutic targeting. Can cause replication fork stalling and genomic instability.",
            "detection_method": "Regex-based identification of purine/pyrimidine mirror repeats and triplex motifs. Scoring by base composition and sequence purity.",
            "literature_references": [
                "Frank-Kamenetskii & Mirkin (1995) Annu Rev Biochem 64:65-95",
                "Vasquez & Wilson (1998) Trends Biochem Sci 23:4-9",
                "Bacolla et al. (2001) J Biol Chem 276:18597-18605"
            ]
        },
        "G_Quadruplex": {
            "description": "Four-stranded DNA structures formed by G-rich sequences through Hoogsteen hydrogen bonding between guanine bases, creating stacked G-quartets stabilized by monovalent cations.",
            "mechanism": "Multiple G-runs connected by loops fold into G-quartets through cyclic Hoogsteen bonding. K+ or Na+ ions coordinate in the central channel providing stability.",
            "biological_significance": "Present in telomeres, promoters, and origins of replication. Involved in telomerase regulation, transcriptional control, and genome stability.",
            "detection_method": "Regex-based detection of canonical and variant G4 motifs using G-run and loop length criteria. G4Hunter algorithm provides thermodynamic scoring.",
            "literature_references": [
                "Burge et al. (2006) Nucleic Acids Res 34:5402-5415",
                "Todd et al. (2005) Nucleic Acids Res 33:2901-2907",
                "Bedrat et al. (2016) Nucleic Acids Res 44:1746-1759"
            ]
        },
        "i_Motif": {
            "description": "Four-stranded DNA structures formed by C-rich sequences through intercalated cytosine-cytosine+ base pairs, creating a compact intercalated structure stable at acidic pH.",
            "mechanism": "Cytosine-rich strands form intercalated parallel-stranded structures through protonated cytosine-cytosine base pairs, requiring acidic conditions for stability.",
            "biological_significance": "Complementary to G-quadruplexes, potentially involved in pH-sensitive gene regulation and chromatin dynamics. May play roles in transcriptional switching.",
            "detection_method": "Regex-based identification of C-rich sequences capable of i-motif formation. Scoring based on C-run count, length, and sequence composition.",
            "literature_references": [
                "Zeraati et al. (2018) Nat Chem 10:631-637",
                "Wright et al. (2018) Nat Commun 9:4709",
                "Abou Assi et al. (2018) Nucleic Acids Res 46:8038-8056"
            ]
        },
        "AC_Motif": {
            "description": "Alternating purine-pyrimidine sequences with specific A-rich and C-rich consensus patterns that can adopt non-B DNA conformations under certain conditions.",
            "mechanism": "Alternating A-rich/C-rich regions create specific sequence contexts that may facilitate alternative DNA structures through altered stacking and pairing.",
            "biological_significance": "May contribute to chromatin organization and DNA accessibility. Potentially involved in nucleosome positioning and protein recognition.",
            "detection_method": "Regex-based identification of alternating A-rich/C-rich consensus regions. Scoring based on pattern presence and sequence characteristics.",
            "literature_references": [
                "Travers & Klug (1987) Nature 327:280-288",
                "Drew & Travers (1985) Nucleic Acids Res 13:4445-4467"
            ]
        },
        "Hybrid_Motif": {
            "description": "Genomic regions where multiple non-B DNA motif classes overlap, creating complex structural environments with potential for multiple alternative conformations.",
            "mechanism": "Sequence regions capable of adopting multiple non-B DNA structures depending on cellular conditions, supercoiling, and protein interactions.",
            "biological_significance": "Hotspots for genomic instability, complex regulation, and structure-function relationships. May serve as integration points for multiple regulatory pathways.",
            "detection_method": "Computational identification through interval intersection of different motif classes. Scoring based on motif diversity and overlap complexity.",
            "literature_references": [
                "Wells et al. (2009) DNA Repair 8:1118-1126",
                "Bacolla & Wells (2009) J Biol Chem 284:23412-23416"
            ]
        },
        "Non_B_DNA_Cluster": {
            "description": "Genomic hotspot regions containing multiple non-B DNA motifs within close proximity, representing areas of high structural diversity and potential instability.",
            "mechanism": "Clustering of different non-B DNA motifs within defined genomic windows, creating regions of enhanced structural complexity and regulatory potential.",
            "biological_significance": "Associated with fragile sites, recombination hotspots, and genomic instability. May serve as regulatory hubs and evolutionary breakpoints.",
            "detection_method": "Sliding window algorithm to identify regions with multiple motif types. Scoring based on motif count, diversity, and spatial distribution.",
            "literature_references": [
                "Wang & Vasquez (2014) DNA Repair 19:143-151",
                "Zhao et al. (2010) Nucleic Acids Res 38:6527-6536"
            ]
        }
    }
    
    # Return comprehensive description or default minimal info
    return descriptions.get(motif_class, {
        "description": f"Non-B DNA structural motif: {motif_class}",
        "mechanism": "Alternative DNA structure formation",
        "biological_significance": "Potential regulatory or structural role",
        "detection_method": "Computational pattern matching",
        "literature_references": ["Multiple references in non-B DNA literature"]
    })

@app.post("/api/v1/analyze", response_model=AnalysisResponse)
async def analyze_sequence(request: SequenceRequest, background_tasks: BackgroundTasks):
    """Analyze DNA sequence for non-B DNA motifs"""
    global analysis_counter
    
    start_time = time.time()
    
    try:
        # Parse and validate sequence
        sequence = parse_fasta(request.sequence)
        if len(sequence) < 10:
            raise HTTPException(status_code=400, detail="Sequence too short (minimum 10 bp)")
        if len(sequence) > 100000:
            raise HTTPException(status_code=400, detail="Sequence too long (maximum 100,000 bp)")
        
        # Generate analysis ID
        analysis_id = hashlib.md5(f"{request.sequence}{time.time()}".encode()).hexdigest()[:12]
        
        # Run motif analysis
        motifs_result = all_motifs_refactored(
            sequence,
            nonoverlap=request.nonoverlap,
            report_hotspots=request.report_hotspots,
            sequence_name=request.sequence_name,
            calculate_conservation=request.calculate_conservation
        )
        
        # Filter by score threshold if specified
        if request.score_threshold > 0:
            motifs_result = [
                m for m in motifs_result 
                if float(m.get('Normalized_Score', 0)) >= request.score_threshold
            ]
        
        # Get sequence statistics
        stats = get_basic_stats(sequence, motifs_result)
        
        processing_time = time.time() - start_time
        analysis_counter += 1
        
        # Prepare response
        response = AnalysisResponse(
            success=True,
            analysis_id=analysis_id,
            timestamp=datetime.now().isoformat(),
            sequence_info={
                "name": request.sequence_name,
                "length": len(sequence),
                "gc_content": stats.get("GC_Content", 0),
                "composition": stats.get("composition", {})
            },
            motifs=motifs_result,
            statistics={
                "total_motifs": len(motifs_result),
                "motif_classes_found": len(set(m.get('Class', '') for m in motifs_result)),
                "sequence_coverage": stats.get("coverage_percentage", 0),
                "motif_density": len(motifs_result) / len(sequence) * 1000,  # per kb
                "class_distribution": stats.get("class_counts", {}),
                "scoring_summary": {
                    "min_score": min([float(m.get('Normalized_Score', 0)) for m in motifs_result] + [0]),
                    "max_score": max([float(m.get('Normalized_Score', 0)) for m in motifs_result] + [0]),
                    "mean_score": sum([float(m.get('Normalized_Score', 0)) for m in motifs_result]) / max(len(motifs_result), 1)
                }
            },
            processing_time=processing_time
        )
        
        # Cache result for potential reuse
        analysis_cache[analysis_id] = response
        
        return response
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Analysis failed: {str(e)}")

@app.get("/api/v1/analyze/{analysis_id}")
async def get_cached_analysis(analysis_id: str):
    """Retrieve cached analysis result"""
    if analysis_id not in analysis_cache:
        raise HTTPException(status_code=404, detail="Analysis not found or expired")
    return analysis_cache[analysis_id]

@app.get("/")
async def root():
    """API root with documentation link"""
    return {
        "message": "NBDFinder REST API",
        "documentation": "/docs",
        "health": "/api/v1/health",
        "version": "1.0.0"
    }

# Cleanup old cache entries periodically
async def cleanup_cache():
    """Remove old cache entries to prevent memory buildup"""
    global analysis_cache
    # Keep only last 100 analyses
    if len(analysis_cache) > 100:
        sorted_cache = sorted(analysis_cache.items(), key=lambda x: x[1].timestamp)
        analysis_cache = dict(sorted_cache[-50:])  # Keep most recent 50

if __name__ == "__main__":
    print("🧬 Starting NBDFinder REST API Server...")
    print("📚 API Documentation: http://localhost:8000/docs")
    print("🔍 Health Check: http://localhost:8000/api/v1/health")
    
    uvicorn.run(
        "api:app",
        host="0.0.0.0",
        port=8000,
        reload=True,
        log_level="info"
    )