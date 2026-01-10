# Modular Architecture Migration Guide

## Overview

This document describes the modular architecture refactoring of NonBDNAFinder to improve code organization, maintainability, and scalability.

## Architecture Structure

```
NonBDNAFinder/
├── app.py                   # Main Streamlit UI (minimal orchestration)
├── nonbscanner.py           # Main API (minimal orchestration)
├── ui/                      # UI Components Module
│   ├── __init__.py
│   ├── layout.py           # Page layout components
│   ├── formatting.py       # Text formatting helpers
│   ├── metrics.py          # Metrics display components
│   ├── progress.py         # Progress indicators
│   ├── inputs.py           # Input widgets
│   └── downloads.py        # Download buttons
├── engine/                  # Core Detection Engine
│   ├── __init__.py
│   ├── detection.py        # Main detection orchestration
│   ├── scoring.py          # Scoring algorithms
│   ├── patterns.py         # Pattern definitions
│   ├── merging.py          # Overlap merging logic
│   ├── chunking.py         # Sequence chunking
│   ├── sequence_ops.py     # Sequence operations ✓ Created
│   └── detectors/          # Individual Detector Classes
│       ├── __init__.py
│       ├── base.py         # Base detector class
│       ├── curved_dna.py   # Curved DNA detector
│       ├── z_dna.py        # Z-DNA detector
│       ├── a_philic.py     # A-philic detector
│       ├── slipped_dna.py  # Slipped DNA detector
│       ├── cruciform.py    # Cruciform detector
│       ├── r_loop.py       # R-loop detector
│       ├── triplex.py      # Triplex detector
│       ├── g_quadruplex.py # G-quadruplex detector
│       └── i_motif.py      # i-Motif detector
└── utils/                   # Shared Utilities
    ├── __init__.py
    ├── registry.py         # Pattern registry loading
    ├── caching.py          # Caching utilities
    ├── state.py            # State management
    ├── export.py           # Data export (CSV, BED, JSON, Excel)
    ├── fasta.py            # FASTA parsing ✓ Created
    ├── validation.py       # Sequence validation
    └── plotting/           # Visualization Functions
        ├── __init__.py
        ├── distributions.py # Distribution plots
        ├── coverage.py     # Coverage maps
        ├── density.py      # Density heatmaps
        ├── statistical.py  # Statistical plots
        └── genomic.py      # Genome-wide visualizations
```

## Module Guidelines

### Size Limits
- Target: ~200 lines per module maximum
- Actual: Flexible based on cohesion
- Split modules that exceed ~300 lines

### Naming Conventions
- Use snake_case for filenames and functions
- Use PascalCase for class names
- Clear, descriptive names indicating purpose

### Documentation
- Each module has a docstring describing its purpose
- Functions have docstrings with Args, Returns, Examples
- Lightweight comments for complex logic only

### Imports
- Absolute imports from root: `from engine.sequence_ops import validate_sequence`
- Relative imports within module: `from .base import BaseDetector`
- Group imports: stdlib, third-party, local

## Migration Steps

### Phase 1: Core Infrastructure (Priority)
1. ✅ Create directory structure
2. ✅ Create __init__.py files
3. ✅ Extract sequence_ops module
4. ✅ Extract fasta module
5. Extract validation module
6. Extract export module

### Phase 2: Engine Modules
7. Extract base detector class
8. Split individual detector classes
9. Extract detection orchestration
10. Extract scoring logic
11. Extract merging logic
12. Extract chunking logic

### Phase 3: UI Modules
13. Extract layout components
14. Extract formatting helpers
15. Extract metrics displays
16. Extract progress indicators
17. Extract input widgets
18. Extract download buttons

### Phase 4: Utility Modules
19. Extract registry loading
20. Extract caching utilities
21. Extract state management
22. Split visualization functions

### Phase 5: Integration
23. Update app.py imports
24. Update nonbscanner.py imports
25. Update tests
26. Verify backwards compatibility

### Phase 6: Optimization
27. Deduplicate repeated code
28. Preload reusable structures
29. Performance benchmarking
30. Documentation updates

## Benefits

### Maintainability
- Easier to locate and update specific functionality
- Clear separation of concerns
- Reduced cognitive load

### Testability
- Individual modules can be tested in isolation
- Easier to mock dependencies
- Faster test execution

### Scalability
- New detectors can be added without modifying existing code
- Parallel processing easier to implement
- Async upgrades more straightforward

### Performance
- Lazy loading of unused modules
- Better caching opportunities
- Reduced import overhead

## Backwards Compatibility

All existing public APIs remain unchanged:
- `nonbscanner.analyze_sequence()`
- `nonbscanner.analyze_file()`
- `nonbscanner.analyze_fasta()`
- Detector class interfaces
- Export functions

Internal refactoring does not affect external usage.

## Testing Strategy

1. **Unit Tests**: Test each module independently
2. **Integration Tests**: Test module interactions
3. **Regression Tests**: Ensure existing behavior preserved
4. **Performance Tests**: Verify no performance degradation

## Migration Checklist for Developers

When extracting a new module:

- [ ] Create module file with descriptive docstring
- [ ] Move related functions/classes to module
- [ ] Add type hints for all function signatures
- [ ] Write docstrings with Args, Returns, Examples
- [ ] Update imports in dependent modules
- [ ] Add unit tests for the module
- [ ] Update this migration guide
- [ ] Verify module is under ~200 lines

## Performance Considerations

### Preloading
- Pattern registries cached on first use
- Detector instances reused across analyses
- Compiled regex patterns stored

### Chunking
- Large sequences split efficiently
- Overlap handling optimized
- Parallel chunk processing where beneficial

### Memory
- Explicit garbage collection at key points
- DataFrame memory optimization
- Streaming for large exports

## Future Enhancements

### Ready for Implementation
- Async detector execution
- Parallel chunk processing
- Distributed computing support
- Plugin architecture for custom detectors

### Architecture Support
- Clear module boundaries enable parallel execution
- State isolation supports async operations
- Caching layer ready for distributed cache
- API stable for external integrations

## References

- Original codebase: ~19K lines in 4 files
- Target: ~80+ modules under 200 lines each
- Estimated migration effort: 2-3 weeks full-time
- Benefits: Immediate for maintainability, long-term for scalability
