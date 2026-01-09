# Badge Tabs Feature - Quick Reference

## Overview
The badge-tabs feature provides professional, reusable status badges for displaying project information throughout the NonBDNAFinder application.

## Quick Start

### Basic Usage
```html
<div class="badge-tabs-container badge-tabs-container--centered">
    <div class="badge-tab badge-tab--version badge-tab--split">
        <span class="badge-tab__label">Version</span>
        <span class="badge-tab__value">2025.1</span>
    </div>
</div>
```

### In Streamlit
```python
st.markdown("""
    <div class="badge-tabs-container badge-tabs-container--centered">
        <div class="badge-tab badge-tab--version badge-tab--split">
            <span class="badge-tab__label">Version</span>
            <span class="badge-tab__value">2025.1</span>
        </div>
    </div>
""", unsafe_allow_html=True)
```

## Available Badge Types
- **version** (Blue) - Version information
- **quality** (Gold) - Quality ratings
- **license** (Green) - License information
- **status** (Teal) - Status indicators
- **build** (Purple) - Build status
- **coverage** (Indigo) - Code coverage
- **downloads** (Orange) - Download counts
- **rating** (Pink) - Ratings/reviews
- **success** (Bright Green) - Success messages
- **warning** (Amber) - Warnings
- **error** (Red) - Error messages
- **info** (Light Blue) - Information
- **neutral** (Gray) - Neutral/default

## Visual Examples
See `badge_tabs_demo.html` for a complete interactive demonstration, or view the screenshots in `badge_tabs_screenshot.png`.

## Documentation
- **Complete Guide**: `BADGE_TABS_GUIDE.md`
- **Implementation Details**: `BADGE_TABS_IMPLEMENTATION.md`
- **Live Demo**: `badge_tabs_demo.html`

## Where It's Used
Currently integrated in:
- Home tab header (app.py) - displaying Version, Quality, License, and Status

## Key Features
✅ 14 color variants
✅ Split-style (label+value) and simple formats
✅ Three size options (small, default, large)
✅ Fully responsive
✅ Hover effects and animations
✅ Icon support
✅ Print-friendly

## Browser Compatibility
- Chrome/Edge 90+
- Firefox 88+
- Safari 14+
- All modern mobile browsers
