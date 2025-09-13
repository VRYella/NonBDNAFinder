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
from classification_config import MOTIF_LENGTH_LIMITS, SCORING_METHODS

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
    """Get detailed information about supported motif classes"""
    motif_info = {
        "classes": len(MOTIF_LENGTH_LIMITS),
        "details": {}
    }
    
    for motif_class, limits in MOTIF_LENGTH_LIMITS.items():
        motif_info["details"][motif_class] = {
            "min_length": limits[0],
            "max_length": limits[1],
            "scoring_method": SCORING_METHODS.get(motif_class, "Pattern-based"),
            "description": get_motif_description(motif_class)
        }
    
    return motif_info

def get_motif_description(motif_class: str) -> str:
    """Get description for motif class"""
    descriptions = {
        "Curved_DNA": "DNA sequences with intrinsic curvature from phased A-tracts",
        "Slipped_DNA": "Direct and tandem repeats that can form slipped structures",
        "Cruciform_DNA": "Palindromic sequences forming cruciform structures",
        "R-Loop": "G-rich sequences prone to RNA-DNA hybrid formation",
        "Triplex_DNA": "Purine/pyrimidine-rich sequences forming triplex structures",
        "G-Quadruplex": "G-rich sequences forming four-stranded G4 structures",
        "i-Motif": "C-rich sequences forming i-motif structures",
        "Z-DNA": "Alternating purine-pyrimidine sequences adopting Z-form",
        "Hybrid": "Overlapping regions between different motif classes",
        "Cluster": "Hotspot regions with multiple motif types"
    }
    return descriptions.get(motif_class, "Non-B DNA structural motif")

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