"""
Utilities package for NonBDNAFinder.

Contains utility modules for:
- Configuration (config/)
- Core functionality (core/)
- Export utilities (export/)
- Visualization (visualization/)
- Main scanner API (nonbscanner.py)
- Detector utilities (detectors_utils.py)
- Job management (job_manager.py)
- Scanner agent (scanner_agent.py)
- General utilities (utilities.py)

Production-grade adaptive chunking architecture:
- SystemResourceInspector  – RAM / CPU / disk availability
- AdaptiveChunkPlanner     – Dynamic chunk-size and worker selection
- ChunkGenerator           – Overlapping genome segment iterator (2 000 bp overlap)
- OverlapDeduplicator      – Remove duplicate motifs at chunk boundaries
- DetectorRunner           – Run detectors on a chunk, write results to disk
- ParallelChunkExecutor    – Parallel / sequential chunk execution
- FinalExporter            – Stream-merge chunk results for full-table export
"""
