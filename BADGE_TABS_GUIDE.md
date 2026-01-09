# Badge Tabs Styling Guide

## Overview

This guide explains how to use the reusable badge-tab CSS classes implemented in `styles.css`. These badge-tabs provide a professional way to display status information like version numbers, quality ratings, license information, and more.

## Features

- **14 Color Variants**: version, quality, license, status, build, coverage, downloads, rating, success, warning, error, info, neutral
- **Two Display Styles**: Simple single-value badges or split label+value badges
- **Responsive Design**: Adapts to different screen sizes
- **Interactive States**: Hover effects, click animations
- **Flexible Layout**: Multiple container alignment options

## Basic Usage

### Simple Badge

```html
<div class="badge-tabs-container">
    <div class="badge-tab badge-tab--version">
        v2025.1
    </div>
</div>
```

### Split Badge (Label + Value)

```html
<div class="badge-tabs-container">
    <div class="badge-tab badge-tab--version badge-tab--split">
        <span class="badge-tab__label">Version</span>
        <span class="badge-tab__value">2025.1</span>
    </div>
</div>
```

### Badge with Icon

```html
<div class="badge-tabs-container">
    <div class="badge-tab badge-tab--quality">
        <span class="badge-tab__icon">⭐</span>
        Nobel-Level
    </div>
</div>
```

## Color Variants

### Version Badge (Blue)
```html
<div class="badge-tab badge-tab--version badge-tab--split">
    <span class="badge-tab__label">Version</span>
    <span class="badge-tab__value">2025.1</span>
</div>
```

### Quality Badge (Gold)
```html
<div class="badge-tab badge-tab--quality badge-tab--split">
    <span class="badge-tab__label">Quality</span>
    <span class="badge-tab__value">Nobel-Level</span>
</div>
```

### License Badge (Green)
```html
<div class="badge-tab badge-tab--license badge-tab--split">
    <span class="badge-tab__label">License</span>
    <span class="badge-tab__value">MIT</span>
</div>
```

### Status Badge (Teal)
```html
<div class="badge-tab badge-tab--status badge-tab--split">
    <span class="badge-tab__label">Status</span>
    <span class="badge-tab__value">Production</span>
</div>
```

### Build Badge (Purple)
```html
<div class="badge-tab badge-tab--build badge-tab--split">
    <span class="badge-tab__label">Build</span>
    <span class="badge-tab__value">Passing</span>
</div>
```

### Coverage Badge (Indigo)
```html
<div class="badge-tab badge-tab--coverage badge-tab--split">
    <span class="badge-tab__label">Coverage</span>
    <span class="badge-tab__value">95%</span>
</div>
```

### Downloads Badge (Orange)
```html
<div class="badge-tab badge-tab--downloads badge-tab--split">
    <span class="badge-tab__label">Downloads</span>
    <span class="badge-tab__value">10K+</span>
</div>
```

### Rating Badge (Pink)
```html
<div class="badge-tab badge-tab--rating badge-tab--split">
    <span class="badge-tab__label">Rating</span>
    <span class="badge-tab__value">⭐⭐⭐⭐⭐</span>
</div>
```

### Success Badge (Bright Green)
```html
<div class="badge-tab badge-tab--success">
    ✓ Validated
</div>
```

### Warning Badge (Amber)
```html
<div class="badge-tab badge-tab--warning">
    ⚠ Beta
</div>
```

### Error Badge (Red)
```html
<div class="badge-tab badge-tab--error">
    ✗ Failed
</div>
```

### Info Badge (Light Blue)
```html
<div class="badge-tab badge-tab--info">
    ℹ Info
</div>
```

### Neutral Badge (Gray)
```html
<div class="badge-tab badge-tab--neutral">
    Unknown
</div>
```

## Container Layouts

### Default (Left-aligned)
```html
<div class="badge-tabs-container">
    <!-- badges -->
</div>
```

### Centered
```html
<div class="badge-tabs-container badge-tabs-container--centered">
    <!-- badges -->
</div>
```

### Space Between
```html
<div class="badge-tabs-container badge-tabs-container--spaced">
    <!-- badges -->
</div>
```

## Size Variants

### Small Badge
```html
<div class="badge-tab badge-tab--version badge-tab--small">
    v2025.1
</div>
```

### Default Badge
```html
<div class="badge-tab badge-tab--version">
    v2025.1
</div>
```

### Large Badge
```html
<div class="badge-tab badge-tab--version badge-tab--large">
    v2025.1
</div>
```

## Clickable/Interactive Badges

For badges that should be clickable (e.g., links):

```html
<a href="https://example.com" class="badge-tab badge-tab--version badge-tab--clickable">
    v2025.1
</a>
```

## Complete Example

Here's a complete example showing multiple badges as they might appear in a project header:

```html
<div class="badge-tabs-container badge-tabs-container--centered">
    <!-- Version Badge -->
    <div class="badge-tab badge-tab--version badge-tab--split">
        <span class="badge-tab__label">Version</span>
        <span class="badge-tab__value">2025.1</span>
    </div>
    
    <!-- Quality Badge -->
    <div class="badge-tab badge-tab--quality badge-tab--split">
        <span class="badge-tab__label">Quality</span>
        <span class="badge-tab__value">Nobel-Level</span>
    </div>
    
    <!-- License Badge -->
    <div class="badge-tab badge-tab--license badge-tab--split">
        <span class="badge-tab__label">License</span>
        <span class="badge-tab__value">MIT</span>
    </div>
    
    <!-- Status Badge -->
    <div class="badge-tab badge-tab--status badge-tab--split">
        <span class="badge-tab__label">Status</span>
        <span class="badge-tab__value">Production</span>
    </div>
    
    <!-- Build Badge -->
    <div class="badge-tab badge-tab--build badge-tab--split">
        <span class="badge-tab__label">Build</span>
        <span class="badge-tab__value">Passing</span>
    </div>
</div>
```

## Usage in Streamlit (app.py)

To use these badges in a Streamlit application, use `st.markdown` with `unsafe_allow_html=True`:

```python
import streamlit as st

# Single badge
st.markdown("""
    <div class="badge-tabs-container badge-tabs-container--centered">
        <div class="badge-tab badge-tab--version badge-tab--split">
            <span class="badge-tab__label">Version</span>
            <span class="badge-tab__value">2025.1</span>
        </div>
    </div>
""", unsafe_allow_html=True)

# Multiple badges
st.markdown("""
    <div class="badge-tabs-container badge-tabs-container--centered">
        <div class="badge-tab badge-tab--version badge-tab--split">
            <span class="badge-tab__label">Version</span>
            <span class="badge-tab__value">2025.1</span>
        </div>
        <div class="badge-tab badge-tab--quality badge-tab--split">
            <span class="badge-tab__label">Quality</span>
            <span class="badge-tab__value">Nobel-Level</span>
        </div>
        <div class="badge-tab badge-tab--license badge-tab--split">
            <span class="badge-tab__label">License</span>
            <span class="badge-tab__value">MIT</span>
        </div>
    </div>
""", unsafe_allow_html=True)
```

## Styling Customization

All badge-tabs use CSS custom properties (CSS variables) defined in `styles.css`:

- `--primary-color`: Primary theme color
- `--secondary-color`: Secondary theme color
- `--border-radius-pill`: Border radius for pill-shaped badges
- `--font-primary`: Primary font family
- `--transition-normal`: Transition timing

You can override these variables to customize the appearance globally.

## Responsive Behavior

The badge-tabs automatically adapt to smaller screens:

- **Desktop**: Full size with normal padding
- **Mobile (≤768px)**: Reduced font size and padding for better fit
- Container wraps badges to multiple lines as needed

## Accessibility

- All badges use semantic HTML
- Color is not the only means of conveying information (text labels are present)
- Interactive badges have proper cursor indicators
- Sufficient color contrast for readability

## Browser Compatibility

The badge-tabs use modern CSS features but are compatible with:
- Chrome/Edge 90+
- Firefox 88+
- Safari 14+
- All modern mobile browsers

## Notes

- Badge-tabs are display elements only and don't have built-in functionality
- Use appropriate semantic HTML (`<a>`, `<button>`, etc.) for interactive badges
- The split badge style (label + value) is recommended for clarity
- Icons can be emoji or font icons (Font Awesome, Material Icons, etc.)
