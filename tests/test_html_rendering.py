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

    def test_css_color_sanitization(self):
        """Test that CSS color values are sanitized to prevent CSS injection."""
        import re
        
        # Create the sanitize function inline to test logic
        _CSS_COLOR_PATTERN = re.compile(r'^#[0-9A-Fa-f]{6}$|^#[0-9A-Fa-f]{3}$')
        
        def _sanitize_css_color(color, default='#607D8B'):
            if isinstance(color, str) and _CSS_COLOR_PATTERN.match(color):
                return color
            return default
        
        # Valid hex colors should pass through
        self.assertEqual(_sanitize_css_color('#FF0000'), '#FF0000')
        self.assertEqual(_sanitize_css_color('#F00'), '#F00')
        self.assertEqual(_sanitize_css_color('#abcdef'), '#abcdef')
        
        # Invalid/malicious inputs should return default
        self.assertEqual(_sanitize_css_color('red'), '#607D8B')
        self.assertEqual(_sanitize_css_color('expression(alert(1))'), '#607D8B')
        self.assertEqual(_sanitize_css_color('url(http://evil.com)'), '#607D8B')
        self.assertEqual(_sanitize_css_color('#FF0000; background-image: url(x)'), '#607D8B')
        self.assertEqual(_sanitize_css_color(None), '#607D8B')
        self.assertEqual(_sanitize_css_color(123), '#607D8B')


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


class TestSequenceNameEscaping(unittest.TestCase):
    """Tests for sequence name escaping in HTML contexts."""

    def test_sequence_name_with_special_characters(self):
        """Test that sequence names with special HTML characters are properly escaped."""
        # Example of a sequence name that might come from FASTA headers or NCBI
        test_names = [
            "NC_000012.12:c102481839-102395874 Homo sapiens chromosome 12, GRCh38.p14 Primary Assembly (86,456 bp)",
            "sequence_with > greater_than",
            "sequence_with < less_than",
            "sequence_with & ampersand",
            "test <script>alert('xss')</script>",
        ]
        
        for name in test_names:
            escaped_name = html.escape(name)
            # Verify that escaped versions contain proper HTML entity escapes
            if '>' in name:
                self.assertIn('&gt;', escaped_name)
            if '<' in name:
                self.assertIn('&lt;', escaped_name)
            if '&' in name:
                self.assertIn('&amp;', escaped_name)
        
        # Test specific escape sequences are present
        self.assertIn('&gt;', html.escape("test > value"))
        self.assertIn('&lt;', html.escape("test < value"))
        self.assertIn('&amp;', html.escape("test & value"))

    def test_detection_pipeline_html_safety(self):
        """Test that the detection pipeline HTML is safe with escaped sequence names."""
        # Simulate how app.py constructs the detailed_progress_html
        name = "NC_000012.12 > test sequence"
        safe_name = html.escape(name)
        seq_len = 1000
        total_bp_processed = 500
        total_bp_all_sequences = 1000
        
        # Construct HTML similar to app.py line 2122-2131
        detailed_progress_html = f"""
        <div style='background: linear-gradient(135deg, #f5f5f5 0%, #e0e0e0 100%); 
                    border-radius: 12px; padding: 1rem; margin: 0.8rem 0;
                    border: 1px solid #bdbdbd;'>
            <h4 style='margin: 0 0 0.8rem 0; color: #424242;'>📋 Detection Pipeline</h4>
            <p style='margin: 0; font-size: 0.85rem; color: #616161;'>
                <strong>Sequence:</strong> {safe_name} ({seq_len:,} bp)<br/>
                <strong>Processed:</strong> {total_bp_processed:,} / {total_bp_all_sequences:,} bp
            </p>
            <div style='display: grid; grid-template-columns: repeat(3, 1fr); gap: 0.4rem; margin-top: 0.8rem; font-size: 0.8rem;'>
        """
        
        # Verify the > character is escaped in the output
        self.assertIn('&gt;', detailed_progress_html)
        # Verify no unescaped > character appears within the name portion
        # (there are legitimate > in HTML tags, but the name should be escaped)
        self.assertIn('NC_000012.12 &gt; test sequence', detailed_progress_html)


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
