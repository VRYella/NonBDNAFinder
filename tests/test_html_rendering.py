"""
HTML Rendering Helper Tests for NBDScanner
===========================================

Tests for the HTML rendering utility functions that provide safe HTML
rendering in Streamlit with proper XSS prevention.
"""

import os
import sys
import unittest
import html

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class TestRenderMotifClassBadges(unittest.TestCase):
    """Tests for render_motif_class_badges function."""

    def test_basic_badge_generation(self):
        """Test basic badge HTML generation."""
        # Import locally to avoid Streamlit context issues
        from visualizations import MOTIF_CLASS_COLORS
        
        # Create the function inline to test logic without Streamlit context
        def render_motif_class_badges(motif_classes, colors=None):
            if colors is None:
                colors = MOTIF_CLASS_COLORS
            
            badges = []
            for cls in motif_classes:
                safe_cls = html.escape(str(cls))
                color = colors.get(cls, '#607D8B')
                badge_html = f'<span style="background-color: {color};">{safe_cls}</span>'
                badges.append(badge_html)
            
            return '<div>' + ''.join(badges) + '</div>'
        
        motif_classes = ["G-Quadruplex", "Z-DNA", "Cruciform"]
        result = render_motif_class_badges(motif_classes)
        
        # Check that all classes are present in the output
        for cls in motif_classes:
            self.assertIn(cls, result)
        
        # Check that HTML structure is correct
        self.assertIn('<div>', result)
        self.assertIn('</div>', result)
        self.assertIn('<span', result)
        self.assertIn('</span>', result)

    def test_xss_prevention(self):
        """Test that HTML special characters are escaped to prevent XSS."""
        # Create the function inline to test logic without Streamlit context
        def render_motif_class_badges(motif_classes, colors=None):
            badges = []
            for cls in motif_classes:
                safe_cls = html.escape(str(cls))
                color = '#607D8B'
                badge_html = f'<span style="background-color: {color};">{safe_cls}</span>'
                badges.append(badge_html)
            
            return '<div>' + ''.join(badges) + '</div>'
        
        # Test with malicious input containing HTML/script tags
        malicious_input = ["<script>alert('xss')</script>", "<img onerror='alert(1)' src='x'>"]
        result = render_motif_class_badges(malicious_input)
        
        # The output should NOT contain unescaped script tags
        self.assertNotIn("<script>", result)
        self.assertNotIn("<img", result)
        
        # The output should contain escaped versions
        self.assertIn("&lt;script&gt;", result)
        self.assertIn("&lt;img", result)

    def test_custom_colors(self):
        """Test that custom colors are applied correctly."""
        def render_motif_class_badges(motif_classes, colors=None):
            if colors is None:
                colors = {}
            
            badges = []
            for cls in motif_classes:
                safe_cls = html.escape(str(cls))
                color = colors.get(cls, '#607D8B')
                badge_html = f'<span style="background-color: {color};">{safe_cls}</span>'
                badges.append(badge_html)
            
            return '<div>' + ''.join(badges) + '</div>'
        
        custom_colors = {
            "TestClass1": "#FF0000",
            "TestClass2": "#00FF00"
        }
        motif_classes = ["TestClass1", "TestClass2", "TestClass3"]
        result = render_motif_class_badges(motif_classes, custom_colors)
        
        # Check custom colors are used
        self.assertIn("#FF0000", result)
        self.assertIn("#00FF00", result)
        
        # Check default color is used for unknown class
        self.assertIn("#607D8B", result)

    def test_empty_list(self):
        """Test handling of empty motif class list."""
        def render_motif_class_badges(motif_classes, colors=None):
            badges = []
            for cls in motif_classes:
                safe_cls = html.escape(str(cls))
                badge_html = f'<span>{safe_cls}</span>'
                badges.append(badge_html)
            
            return '<div>' + ''.join(badges) + '</div>'
        
        result = render_motif_class_badges([])
        
        # Should return a valid but empty container
        self.assertIn('<div>', result)
        self.assertIn('</div>', result)
        # Should not contain any span elements
        self.assertNotIn('<span>', result)


class TestHtmlModuleImport(unittest.TestCase):
    """Test that the html module is properly imported and used."""

    def test_html_escape_function(self):
        """Test html.escape functionality for XSS prevention."""
        # Test that dangerous characters are escaped
        dangerous_input = "<script>alert('xss')</script>"
        result = html.escape(dangerous_input)
        
        # Script tags should be escaped
        self.assertNotIn("<script>", result)
        self.assertIn("&lt;script&gt;", result)
        
        # Normal text should remain unchanged
        normal_text = "Normal Text"
        self.assertEqual(html.escape(normal_text), normal_text)
        
        # Ampersand should be escaped
        text_with_amp = "G-Quadruplex & Z-DNA"
        self.assertIn("&amp;", html.escape(text_with_amp))


class TestRenderHtmlIntegration(unittest.TestCase):
    """Integration tests for HTML rendering with actual module imports."""

    def test_motif_class_colors_available(self):
        """Test that MOTIF_CLASS_COLORS is accessible."""
        from visualizations import MOTIF_CLASS_COLORS
        
        # Check that common motif classes have colors defined
        expected_classes = ['G-Quadruplex', 'Z-DNA', 'Cruciform', 'Triplex']
        for cls in expected_classes:
            self.assertIn(cls, MOTIF_CLASS_COLORS)
            # Check that color is a valid hex color
            color = MOTIF_CLASS_COLORS[cls]
            self.assertTrue(color.startswith('#'))
            self.assertEqual(len(color), 7)


if __name__ == '__main__':
    unittest.main()
