# NonBDNAFinder - Complete Documentation

## Overview
NonBDNAFinder is a comprehensive tool for detecting and analyzing 11 distinct classes of Non-B DNA structures in genomic sequences. The tool has been completely reorganized for better maintainability and performance.

## New Directory Structure

### Main Files
- **app.py**: Main Streamlit application - the only entry point
- **requirements.txt**: Python dependencies
- **documentation.md**: This comprehensive documentation

### Core Directories

#### REGISTRIES/
Contains motif-specific registry scripts and pattern definitions:
- `regex_registry.py`: Central pattern registry
- `motifs/`: Motif definitions and configurations
  - `base.py`: Core data structures and candidate definitions
  - `registry.py`: Pattern management for hyperscan integration
  - `visualization.py`: Motif visualization utilities
  - `enhanced_visualization.py`: Advanced plotting functions
  - `__init__.py`: Package initialization

#### HYPERSCAN/
High-performance hyperscan implementation:
- `hyperscan_integration.py`: Main integration layer for hyperscan
- `hs_registry_manager.py`: Hyperscan registry management
- `production_hyperscan_streamlit.py`: Production hyperscan implementation

#### NONHYPERSCAN/
Alternative motif detection without hyperscan dependency:
- `orchestrator.py`: Main pipeline coordinator
- `motif_detectors.py`: Core detection algorithms for all 11 motif classes
- `overlap_resolution.py`: Advanced overlap resolution algorithms
- `motif_classification.py`: Motif classification utilities
- `all_motifs_refactored.py`: Refactored motif detection wrapper

#### SCORINGS/
Motif-specific scoring systems:
- `zdna_calculator.py`: Z-DNA scoring algorithms
- `zdna_hs.py`: Z-DNA hyperscan-based detection
- `conservation_analysis.py`: Conservation scoring across sequences

#### BASE_CODES/
Utility and base functionality:
- `utils.py`: Core utility functions (FASTA parsing, sequence manipulation)
- `constants.py`: Application constants and configurations
- `enhanced_cache.py`: Performance caching system
- `export_utils.py`: Data export utilities
- `viz_tools.py`: Enhanced visualization tools

#### Archives/
Deprecated and archived files:
- Historical test files
- Old implementation versions
- Documentation archives
- Example data files
- Standalone implementations

## Supported Non-B DNA Classes

### Official 11 Classes
1. **Curved DNA**: Intrinsically curved DNA sequences with phased A/T tracts
2. **Slipped DNA**: Direct repeats and short tandem repeats
3. **Cruciform DNA**: Palindromic inverted repeats
4. **R-loop**: RNA-DNA hybrid structures
5. **Triplex**: Triple-stranded DNA structures
6. **G-Quadruplex Family**: Four-stranded G-rich structures
7. **i-motif family**: C-rich structures forming tetraplex
8. **Z-DNA**: Left-handed DNA helices
9. **A-philic DNA**: A-tract containing regions
10. **Hybrid**: Overlapping motif regions
11. **Non-B DNA cluster regions**: Hotspots with multiple motifs

## Key Features

### High Performance
- **Hyperscan Integration**: Ultra-fast pattern matching using Intel Hyperscan
- **Parallel Processing**: Multi-threaded detection algorithms
- **Smart Caching**: Enhanced caching system for repeated analyses
- **Memory Optimization**: Efficient memory usage for large sequences

### Scientific Accuracy
- **Normalized Scoring**: Length-constraint based normalization (0-1 scale)
- **Evidence-based Methods**: Scoring based on published literature
- **Validation**: Comprehensive validation of detected motifs
- **Conservation Analysis**: Cross-sequence conservation scoring

### Advanced Analytics
- **Overlap Resolution**: Sophisticated algorithms for handling overlapping motifs
- **Hotspot Detection**: Identification of Non-B DNA clustering regions
- **Statistical Analysis**: Comprehensive statistical metrics
- **Comparative Analysis**: Multi-sequence comparison tools

### Visualization Suite
- **21+ Chart Types**: Comprehensive visualization options
- **Interactive Plots**: Plotly-powered interactive visualizations
- **Information-based Organization**: Charts organized by information type
- **Export Capabilities**: Multiple export formats (CSV, Excel, images)

## Usage

### Basic Usage
```bash
streamlit run app.py
```

### Input Methods
1. **File Upload**: Support for FASTA/multi-FASTA files
2. **Text Input**: Direct sequence pasting
3. **Example Data**: Built-in example sequences
4. **NCBI Fetch**: Direct download from NCBI databases

### Analysis Parameters
- **Comprehensive Analysis**: All 11 motif classes analyzed by default
- **Normalized Scoring**: Automatic score normalization for fair comparison
- **Overlap Handling**: Complete motif landscape with overlap detection
- **Hotspot Detection**: Automatic cluster region identification

## Technical Specifications

### Dependencies
- Python 3.7+
- Streamlit for web interface
- Hyperscan for high-performance pattern matching
- BioPython for sequence handling
- NumPy/Pandas for data processing
- Matplotlib/Plotly for visualization

### Performance
- **Speed**: 10-100x faster than regex-only approaches
- **Memory**: Optimized for sequences up to 1MB+
- **Scalability**: Handles multi-FASTA files efficiently
- **Caching**: Smart caching reduces repeated computation

### Output Formats
- **Interactive Tables**: Sortable, filterable results
- **Visualizations**: 21+ chart types
- **Export Options**: CSV, Excel, JSON formats
- **Genome Browser**: Compatible export formats

## Development and Maintenance

### Code Organization
- **Modular Design**: Clear separation of concerns
- **Package Structure**: Proper Python package organization
- **Import Management**: Systematic import handling
- **Documentation**: Comprehensive inline documentation

### Testing
- **Unit Tests**: Individual component testing
- **Integration Tests**: End-to-end functionality testing
- **Performance Tests**: Speed and memory validation
- **Scientific Validation**: Accuracy verification

### Extensibility
- **Plugin Architecture**: Easy addition of new motif types
- **Configurable Scoring**: Customizable scoring algorithms
- **Flexible Visualization**: Extensible chart types
- **API Ready**: Prepared for API integration

## Scientific References

### Core Algorithms
- Bedrat et al. (2016) Nucleic Acids Research - G-quadruplex detection
- Ho et al. (2010) Nature Chemical Biology - Z-DNA characterization
- Kim et al. (2018) Nucleic Acids Research - R-loop prediction
- Zeraati et al. (2018) Nature Chemistry - i-motif structures
- Bacolla et al. (2006) Nucleic Acids Research - Cruciform DNA
- Mirkin & Frank-Kamenetskii (1994) Annual Review of Biophysics - Slipped structures

### Methodology
- Length-constraint normalization for fair motif comparison
- Evidence-based scoring thresholds
- Conservative overlap resolution strategies
- Comprehensive validation protocols

## Support and Contact

**Developer**: Dr. Venkata Rajesh Yella  
**Email**: yvrajesh_bt@kluniversity.in  
**GitHub**: https://github.com/VRYella/NonBDNAFinder  

## License
This project is licensed under the terms specified in the LICENSE file.

---

*Last Updated: September 2024*  
*Version: 2.0 - Reorganized Architecture*