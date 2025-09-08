"""
Official Non-B DNA Classification Structure for NBDFinder
========================================================

Scientific Classification Structure:
1. Curved DNA (2 subclasses)
2. Slipped DNA (2 subclasses) 
3. Cruciform DNA (1 subclass)
4. R-loop (1 subclass)
5. Triplex (2 subclasses)
6. G-Quadruplex Family (7 subclasses)
7. i-motif family (3 subclasses)
8. Z-DNA (2 subclasses)
9. Hybrid (variable subclasses based on overlaps)
10. Non-B DNA cluster regions (variable subclasses)

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

# Official 10 Non-B DNA Classes and 22 Subclasses Configuration
OFFICIAL_CLASSIFICATION = {
    1: {
        "class_name": "Curved DNA",
        "subclasses": [
            "Global curvature",
            "Local Curvature"
        ]
    },
    2: {
        "class_name": "Slipped DNA", 
        "subclasses": [
            "Slipped DNA [Direct Repeat]",
            "Slipped DNA [STR]"
        ]
    },
    3: {
        "class_name": "Cruciform DNA",
        "subclasses": [
            "Cruciform DNA [IR]/HairPin [IR]"
        ]
    },
    4: {
        "class_name": "R-loop",
        "subclasses": [
            "R-loop"
        ]
    },
    5: {
        "class_name": "Triplex",
        "subclasses": [
            "Triplex",
            "sticky DNA"
        ]
    },
    6: {
        "class_name": "G-Quadruplex Family",
        "subclasses": [
            "Multimeric G4",
            "Canonical G4", 
            "Relaxed G4",
            "Bulged G4",
            "Bipartite G4",
            "Imperfect G4",
            "G-Triplex intermediate"
        ],
        "priority_order": [
            "Multimeric G4",
            "Canonical G4",
            "Relaxed G4", 
            "Bulged G4",
            "Bipartite G4",
            "Imperfect G4",
            "G-Triplex intermediate"
        ]
    },
    7: {
        "class_name": "i-motif family",
        "subclasses": [
            "Canonical i-motif",
            "Relaxed i-motif", 
            "AC-motif"
        ]
    },
    8: {
        "class_name": "Z-DNA",
        "subclasses": [
            "Z-DNA",
            "eGZ (Extruded-G) DNA"
        ]
    },
    9: {
        "class_name": "Hybrid",
        "subclasses": []  # Dynamic based on overlaps between any two classes
    },
    10: {
        "class_name": "Non-B DNA cluster regions",
        "subclasses": []  # Dynamic based on occurring motifs: any three classes occurring 3+ times in 100 nt
    }
}

# Motif ID mapping for each subclass
MOTIF_IDS = {
    # Class 1: Curved DNA
    "Global curvature": "1.1",
    "Local Curvature": "1.2",
    
    # Class 2: Slipped DNA
    "Slipped DNA [Direct Repeat]": "2.1",
    "Slipped DNA [STR]": "2.2",
    
    # Class 3: Cruciform DNA
    "Cruciform DNA [IR]/HairPin [IR]": "3.1",
    
    # Class 4: R-loop
    "R-loop": "4.1",
    
    # Class 5: Triplex
    "Triplex": "5.1",
    "sticky DNA": "5.2",
    
    # Class 6: G-Quadruplex Family
    "Multimeric G4": "6.1",
    "Canonical G4": "6.2", 
    "Relaxed G4": "6.3",
    "Bulged G4": "6.4",
    "Bipartite G4": "6.5",
    "Imperfect G4": "6.6",
    "G-Triplex intermediate": "6.7",
    
    # Class 7: i-motif family
    "Canonical i-motif": "7.1",
    "Relaxed i-motif": "7.2", 
    "AC-motif": "7.3",
    
    # Class 8: Z-DNA
    "Z-DNA": "8.1",
    "eGZ (Extruded-G) DNA": "8.2",
    
    # Class 9 & 10 will be dynamic
}

# Current implementation mapping (for compatibility)
CURRENT_TO_OFFICIAL = {
    # Curved DNA mappings
    "Global_Array": "Global curvature",
    "Local_Tract": "Local Curvature",
    
    # Slipped DNA mappings
    "Direct_Repeat": "Slipped DNA [Direct Repeat]",
    "STR": "Slipped DNA [STR]",
    
    # Cruciform DNA mappings
    "Cruciform_DNA": "Cruciform DNA [IR]/HairPin [IR]",
    
    # R-loop mappings
    "R-loop": "R-loop",
    "RLFS_m1": "R-loop",
    "RLFS_m2": "R-loop",
    
    # Triplex mappings
    "Triplex": "Triplex",
    "Sticky_DNA": "sticky DNA",
    "Mirror_Repeat": "Triplex",
    
    # G-Quadruplex mappings
    "Multimeric_G4": "Multimeric G4",
    "Canonical_G4": "Canonical G4",
    "Relaxed_G4": "Relaxed G4",
    "Bulged_G4": "Bulged G4",
    "Bipartite_G4": "Bipartite G4",
    "Imperfect_G4": "Imperfect G4",
    "G-Triplex_intermediate": "G-Triplex intermediate",
    
    # i-motif mappings
    "Canonical_iMotif": "Canonical i-motif",
    "Relaxed_iMotif": "Relaxed i-motif",
    "AC-motif": "AC-motif",
    "Other_iMotif": "Relaxed i-motif",  # Map unclassified i-motifs to relaxed category
    
    # Z-DNA mappings
    "Z-DNA": "Z-DNA",
    "eGZ": "eGZ (Extruded-G) DNA",
}

def get_motif_id(subclass_name):
    """Get the official motif ID for a subclass"""
    # First try to map current implementation names to official names
    official_name = CURRENT_TO_OFFICIAL.get(subclass_name, subclass_name)
    return MOTIF_IDS.get(official_name, "0.0")

def get_official_subclass_name(current_name):
    """Get the official subclass name from current implementation name"""
    return CURRENT_TO_OFFICIAL.get(current_name, current_name)

def update_motif_with_ids(motif):
    """Update a motif dictionary with official classification and motif IDs"""
    current_subclass = motif.get("Subclass", "")
    official_subclass = get_official_subclass_name(current_subclass)
    motif_id = get_motif_id(current_subclass)
    
    # Update the motif with official classification
    motif["Subclass"] = official_subclass
    motif["Motif_ID"] = f"{motif.get('Class', 'Unknown')}_{motif_id}_{motif.get('Start', '')}-{motif.get('End', '')}"
    
    return motif

# Export main functions
__all__ = [
    'OFFICIAL_CLASSIFICATION',
    'MOTIF_IDS',
    'CURRENT_TO_OFFICIAL', 
    'get_motif_id',
    'get_official_subclass_name',
    'update_motif_with_ids'
]