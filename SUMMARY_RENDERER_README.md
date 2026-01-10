# Summary Renderer Module

## Overview

The `summary_renderer.py` module provides robust HTML rendering utilities for the NonBDNAFinder Streamlit application. It ensures proper HTML display without escaping issues and provides guaranteed-visible progress bars.

## Problem Solved

This module addresses two critical UI issues:

1. **HTML Escaping**: Previously, HTML blocks in summary sections were being escaped by Streamlit, displaying `<div>` tags instead of rendering them properly.
2. **Progress Bar Visibility**: Progress bars were failing to appear due to unsafe HTML rendering inside containers and CSS conflicts.

## Key Features

### 1. Safe HTML Rendering
- No HTML escaping (renders actual HTML tags)
- Automatic whitespace stripping
- Basic HTML validation
- Compatible with CSS overrides

### 2. Guaranteed-Visible Progress Bars
- Native Streamlit widget (not HTML-based)
- Isolated containers prevent CSS conflicts
- Works inside columns, tabs, and expanders
- Visible against all themes

## API Reference

### `render_summary_block(html, container=None)`

Renders HTML summary blocks with proper formatting.

**Parameters:**
- `html` (str): HTML string to render
- `container` (optional): Streamlit container to render into

**Example:**
```python
from summary_renderer import render_summary_block

render_summary_block("""
    <div style="background: #f0f9ff; padding: 1rem;">
        <b>Processing time:</b> 00:03:42
        <br>
        <b>Motifs detected:</b> 33
    </div>
""")
```

### `render_progress(value, label=None, container=None)`

Renders a guaranteed-visible progress bar.

**Parameters:**
- `value` (float): Progress value between 0.0 and 1.0
- `label` (str, optional): Text label above progress bar
- `container` (optional): Streamlit container to render into

**Example:**
```python
from summary_renderer import render_progress

render_progress(0.75, label="Analysis Progress")
```

### `render_metric_card(label, value, delta=None, color="#10b981", container=None)`

Renders a styled metric card.

**Parameters:**
- `label` (str): Metric label
- `value` (str): Metric value
- `delta` (str, optional): Delta indicator
- `color` (str): Primary color (hex format)
- `container` (optional): Streamlit container

**Example:**
```python
from summary_renderer import render_metric_card

render_metric_card("Total Motifs", "1,234", delta="+15%")
```

### `render_info_box(title, content, box_type="info", container=None)`

Renders an information box with icon.

**Parameters:**
- `title` (str): Box title
- `content` (str): Box content
- `box_type` (str): Type - "info", "warning", "error", "success"
- `container` (optional): Streamlit container

**Example:**
```python
from summary_renderer import render_info_box

render_info_box(
    "Important", 
    "Please validate your results", 
    box_type="warning"
)
```

## Implementation Details

### HTML Validation
The module performs basic HTML validation by checking for balanced `<div>` tags. If unmatched tags are detected, a warning is displayed.

### Whitespace Handling
Multiline HTML strings are automatically cleaned:
1. Strip leading/trailing whitespace
2. Detect minimum indentation
3. Remove common leading whitespace from all lines

### Progress Value Clamping
Progress values are automatically clamped to the valid [0.0, 1.0] range to prevent errors.

## Integration with NonBDNAFinder

The module is integrated into `app.py` at the following locations:

1. **Performance Metrics Section** (line ~2872): Uses `render_summary_block()`
2. **Results Summary Section** (line ~2950): Uses `render_summary_block()`
3. **Chunk Progress Display** (line ~2528): Uses `render_summary_block()` + native `st.progress()`
4. **Download Page Info Boxes** (lines ~3372, 3465, 3597): Use `render_info_box()`

## Benefits

✅ **Correct HTML Rendering**: No more escaped HTML tags
✅ **Progress Bar Visibility**: Always visible across all themes
✅ **Consistent Styling**: Centralized rendering logic
✅ **Better Maintainability**: Single source of truth for UI components
✅ **Type Safety**: Clear function signatures
✅ **Error Prevention**: Automatic validation and clamping

## Testing

All rendering functions have been validated with unit tests:
- Import validation
- HTML tag balancing
- Whitespace stripping
- Progress value clamping

## Version History

- **v2025.1**: Initial release
  - Created `summary_renderer.py` module
  - Implemented 4 core rendering functions
  - Integrated into NonBDNAFinder application
  - Fixed HTML escaping and progress bar visibility issues
