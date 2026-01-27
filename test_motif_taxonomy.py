"""
Test suite for motif taxonomy and normalization
Tests the canonical motif classification system and enforcement
"""

import pytest
from config.motif_taxonomy import (
    MOTIF_CLASSIFICATION,
    VALID_CLASSES,
    VALID_SUBCLASSES,
    CLASS_TO_SUBCLASSES,
    SUBCLASS_TO_CLASS,
    get_all_classes,
    get_subclasses_for_class,
    get_class_for_subclass,
    is_valid_class,
    is_valid_subclass,
    is_valid_pairing,
)
from core.motif_normalizer import (
    MotifNormalizationError,
    normalize_class_name,
    normalize_subclass_name,
    normalize_class_subclass,
    normalize_motif_dict,
    validate_motif_dict,
)
from export.export_validator import (
    ExportValidationError,
    validate_single_motif,
    validate_export_data,
    get_export_summary,
    check_class_completeness,
)


class TestMotifTaxonomy:
    """Test canonical motif taxonomy structure"""
    
    def test_taxonomy_has_11_classes(self):
        """Verify we have exactly 11 motif classes"""
        assert len(MOTIF_CLASSIFICATION) == 11
        assert len(VALID_CLASSES) == 11
    
    def test_all_classes_unique(self):
        """Verify all class names are unique"""
        class_names = [entry['class'] for entry in MOTIF_CLASSIFICATION.values()]
        assert len(class_names) == len(set(class_names))
    
    def test_all_subclasses_unique(self):
        """Verify all subclass names are unique across taxonomy"""
        all_subclasses = []
        for entry in MOTIF_CLASSIFICATION.values():
            all_subclasses.extend(entry['subclasses'])
        assert len(all_subclasses) == len(set(all_subclasses))
    
    def test_required_classes_present(self):
        """Verify all required classes are in taxonomy"""
        required_classes = [
            'Curved_DNA',
            'Slipped_DNA',
            'Cruciform',
            'R-Loop',
            'Triplex',
            'G-Quadruplex',
            'i-Motif',
            'Z-DNA',
            'A-philic_DNA',
            'Hybrid',
            'Non-B_DNA_Clusters'
        ]
        for class_name in required_classes:
            assert class_name in VALID_CLASSES, f"Missing required class: {class_name}"
    
    def test_gquad_has_8_subclasses(self):
        """Verify G-Quadruplex has exactly 8 subclasses"""
        gquad_subclasses = get_subclasses_for_class('G-Quadruplex')
        assert len(gquad_subclasses) == 8
    
    def test_imotif_has_3_subclasses(self):
        """Verify i-Motif has exactly 3 subclasses"""
        imotif_subclasses = get_subclasses_for_class('i-Motif')
        assert len(imotif_subclasses) == 3
        assert 'Canonical i-motif' in imotif_subclasses
        assert 'Relaxed i-motif' in imotif_subclasses
        assert 'AC-motif' in imotif_subclasses
    
    def test_subclass_to_class_mapping(self):
        """Verify reverse mapping from subclass to class"""
        assert get_class_for_subclass('Telomeric G4') == 'G-Quadruplex'
        assert get_class_for_subclass('Z-DNA') == 'Z-DNA'
        assert get_class_for_subclass('Cruciform forming IRs') == 'Cruciform'
    
    def test_class_to_subclasses_mapping(self):
        """Verify forward mapping from class to subclasses"""
        curved_subs = CLASS_TO_SUBCLASSES['Curved_DNA']
        assert 'Global Curvature' in curved_subs
        assert 'Local Curvature' in curved_subs


class TestClassNormalization:
    """Test class name normalization"""
    
    def test_canonical_class_unchanged(self):
        """Canonical class names should pass through unchanged"""
        assert normalize_class_name('G-Quadruplex') == 'G-Quadruplex'
        assert normalize_class_name('Curved_DNA') == 'Curved_DNA'
    
    def test_lowercase_normalization(self):
        """Lowercase variants should normalize to canonical"""
        assert normalize_class_name('g-quadruplex') == 'G-Quadruplex'
        assert normalize_class_name('curved dna') == 'Curved_DNA'
    
    def test_alias_normalization(self):
        """Alias variants should normalize to canonical"""
        assert normalize_class_name('gquadruplex') == 'G-Quadruplex'
        assert normalize_class_name('z dna') == 'Z-DNA'
        assert normalize_class_name('i motif') == 'i-Motif'
    
    def test_invalid_class_strict_mode(self):
        """Invalid class should raise error in strict mode"""
        with pytest.raises(MotifNormalizationError):
            normalize_class_name('InvalidClass', strict=True)
    
    def test_invalid_class_non_strict(self):
        """Invalid class should return as-is in non-strict mode"""
        result = normalize_class_name('InvalidClass', strict=False)
        assert result == 'InvalidClass'


class TestSubclassNormalization:
    """Test subclass name normalization"""
    
    def test_canonical_subclass_unchanged(self):
        """Canonical subclass names should pass through unchanged"""
        assert normalize_subclass_name('Telomeric G4') == 'Telomeric G4'
        assert normalize_subclass_name('Global Curvature') == 'Global Curvature'
    
    def test_lowercase_normalization(self):
        """Lowercase variants should normalize to canonical"""
        assert normalize_subclass_name('telomeric g4') == 'Telomeric G4'
        assert normalize_subclass_name('global curvature') == 'Global Curvature'
    
    def test_alias_normalization(self):
        """Alias variants should normalize to canonical"""
        assert normalize_subclass_name('canonical_imotif') == 'Canonical i-motif'
        assert normalize_subclass_name('egz') == 'eGZ'
        assert normalize_subclass_name('inverted repeat') == 'Cruciform forming IRs'
    
    def test_invalid_subclass_strict_mode(self):
        """Invalid subclass should raise error in strict mode"""
        with pytest.raises(MotifNormalizationError):
            normalize_subclass_name('InvalidSubclass', strict=True)


class TestPairingValidation:
    """Test class/subclass pairing validation"""
    
    def test_valid_pairings(self):
        """Valid class/subclass pairings should be accepted"""
        assert is_valid_pairing('G-Quadruplex', 'Telomeric G4') is True
        assert is_valid_pairing('Cruciform', 'Cruciform forming IRs') is True
        assert is_valid_pairing('Z-DNA', 'eGZ') is True
    
    def test_invalid_pairings(self):
        """Invalid class/subclass pairings should be rejected"""
        # Telomeric G4 belongs to G-Quadruplex, not Z-DNA
        assert is_valid_pairing('Z-DNA', 'Telomeric G4') is False
        # Global Curvature belongs to Curved_DNA, not Triplex
        assert is_valid_pairing('Triplex', 'Global Curvature') is False
    
    def test_normalize_with_autocorrect(self):
        """Auto-correct should fix mismatched class"""
        # Telomeric G4 belongs to G-Quadruplex
        canonical_class, canonical_subclass = normalize_class_subclass(
            'Z-DNA',  # Wrong class
            'Telomeric G4',  # Correct subclass
            auto_correct=True
        )
        assert canonical_class == 'G-Quadruplex'  # Should be corrected
        assert canonical_subclass == 'Telomeric G4'
    
    def test_normalize_without_autocorrect_strict(self):
        """Strict mode without auto-correct should raise error on mismatch"""
        with pytest.raises(MotifNormalizationError):
            normalize_class_subclass(
                'Z-DNA',
                'Telomeric G4',
                strict=True,
                auto_correct=False
            )


class TestMotifDictNormalization:
    """Test normalization of motif dictionaries"""
    
    def test_normalize_valid_motif(self):
        """Valid motif dict should normalize successfully"""
        motif = {
            'Class': 'g-quadruplex',
            'Subclass': 'telomeric g4',
            'Start': 100,
            'End': 125
        }
        normalized = normalize_motif_dict(motif)
        assert normalized['Class'] == 'G-Quadruplex'
        assert normalized['Subclass'] == 'Telomeric G4'
    
    def test_normalize_with_autocorrect(self):
        """Mismatched class should be auto-corrected"""
        motif = {
            'Class': 'Triplex',  # Wrong
            'Subclass': 'Telomeric G4',  # Belongs to G-Quadruplex
        }
        normalized = normalize_motif_dict(motif, auto_correct=True)
        assert normalized['Class'] == 'G-Quadruplex'
    
    def test_validate_valid_motif(self):
        """Valid motif should pass validation"""
        motif = {
            'Class': 'G-Quadruplex',
            'Subclass': 'Telomeric G4',
        }
        assert validate_motif_dict(motif, raise_on_error=False) is True
    
    def test_validate_invalid_motif(self):
        """Invalid motif should fail validation"""
        motif = {
            'Class': 'InvalidClass',
            'Subclass': 'InvalidSubclass',
        }
        assert validate_motif_dict(motif, raise_on_error=False) is False


class TestExportValidation:
    """Test export validation"""
    
    def test_validate_empty_list(self):
        """Empty motif list should return empty"""
        result = validate_export_data([])
        assert result == []
    
    def test_validate_valid_motifs(self):
        """Valid motifs should pass through"""
        motifs = [
            {'Class': 'G-Quadruplex', 'Subclass': 'Telomeric G4', 'Start': 0, 'End': 25},
            {'Class': 'Z-DNA', 'Subclass': 'Z-DNA', 'Start': 50, 'End': 75},
        ]
        result = validate_export_data(motifs, strict=False)
        assert len(result) == 2
    
    def test_validate_with_normalization(self):
        """Invalid motifs should be normalized when auto_normalize=True"""
        motifs = [
            {'Class': 'g-quadruplex', 'Subclass': 'telomeric g4', 'Start': 0, 'End': 25},
        ]
        result = validate_export_data(motifs, auto_normalize=True, strict=False)
        assert len(result) == 1
        assert result[0]['Class'] == 'G-Quadruplex'
        assert result[0]['Subclass'] == 'Telomeric G4'
    
    def test_export_summary(self):
        """Export summary should provide correct statistics"""
        motifs = [
            {'Class': 'G-Quadruplex', 'Subclass': 'Telomeric G4'},
            {'Class': 'G-Quadruplex', 'Subclass': 'Canonical intramolecular G4'},
            {'Class': 'Z-DNA', 'Subclass': 'Z-DNA'},
        ]
        summary = get_export_summary(motifs)
        assert summary['total_motifs'] == 3
        assert summary['unique_classes'] == 2
        assert 'G-Quadruplex' in summary['classes']
        assert 'Z-DNA' in summary['classes']
    
    def test_class_completeness(self):
        """Class completeness check should identify missing classes"""
        motifs = [
            {'Class': 'G-Quadruplex', 'Subclass': 'Telomeric G4'},
            {'Class': 'Z-DNA', 'Subclass': 'Z-DNA'},
        ]
        completeness = check_class_completeness(motifs)
        assert len(completeness['present_classes']) == 2
        assert len(completeness['missing_classes']) == 9  # 11 total - 2 present
        assert 'Curved_DNA' in completeness['missing_classes']


class TestIntegration:
    """Integration tests for end-to-end workflow"""
    
    def test_detector_to_export_workflow(self):
        """Simulate detector output → normalization → export workflow"""
        # Simulate detector output (may have variant naming)
        detector_output = [
            {'Class': 'g-quadruplex', 'Subclass': 'telomeric g4', 'Start': 0, 'End': 25},
            {'Class': 'Curved_DNA', 'Subclass': 'Global Curvature', 'Start': 30, 'End': 55},
        ]
        
        # Normalize
        normalized = [normalize_motif_dict(m, auto_correct=True) for m in detector_output]
        
        # Validate for export
        validated = validate_export_data(normalized, auto_normalize=True, strict=False)
        
        # Check results
        assert len(validated) == 2
        assert validated[0]['Class'] == 'G-Quadruplex'
        assert validated[0]['Subclass'] == 'Telomeric G4'
        assert validated[1]['Class'] == 'Curved_DNA'
        assert validated[1]['Subclass'] == 'Global Curvature'


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
