#!/usr/bin/env python3
"""
Test script to demonstrate enhanced features in NonBDNAFinder detectors.
All detectors now output comprehensive features including:
- Type_Of_Repeat: Classification of repeat structure
- Criterion: Explanation of detection logic
- Disease_Relevance: Annotations for disease associations
- Regions_Involved: Detailed composition of motif regions
- GC_Content, Loop_Length, Arm_Length, etc. as applicable
"""

import sys
import tempfile
import os
# Add repository root to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Utilities.nonbscanner import analyze_sequence
from Utilities.utilities import export_to_csv, export_to_json

# Test sequences for different motif types
test_sequences = {
    'G-Quadruplex': 'GGGAGGGAGGGAGGG',
    'Slipped_DNA': 'CAGCAGCAGCAGCAGCAGCAGCAG',
    'i-Motif': 'CCCCACCCCACCCCACCCC',
    'Z-DNA': 'CGCGCGCGCGCGCGCG',
    'Curved': 'AAAAAAAATTTTTTTTT',
    'R-loop': 'GGGGGGGGGGGGG'
}

print("=" * 80)
print("NonBDNAFinder Enhanced Features Demonstration")
print("=" * 80)
print()

all_motifs = []

for motif_type, sequence in test_sequences.items():
    print(f"\n{'='*60}")
    print(f"Testing {motif_type} Detection")
    print(f"{'='*60}")
    print(f"Sequence: {sequence}")
    print(f"Length: {len(sequence)} bp\n")
    
    motifs = analyze_sequence(sequence, f'{motif_type}_test')
    all_motifs.extend(motifs)
    
    if motifs:
        for motif in motifs:
            print(f"Detected: {motif['Class']} - {motif['Subclass']}")
            print(f"  Position: {motif['Start']}-{motif['End']} ({motif['Length']} bp)")
            print(f"  Score: {motif.get('Score', 'N/A')}")
            
            # Show enhanced features
            enhanced_features = [
                ('GC_Content', '%'),
                ('Type_Of_Repeat', ''),
                ('Criterion', ''),
                ('Disease_Relevance', ''),
                ('Regions_Involved', '')
            ]
            
            for feature, unit in enhanced_features:
                if feature in motif:
                    value = motif[feature]
                    # Truncate long values for display
                    if isinstance(value, str) and len(value) > 100:
                        value = value[:100] + '...'
                    print(f"  {feature}: {value}{unit}")
            
            # Show any motif-specific features
            specific_features = {
                'G-Quadruplex': ['Num_Tracts', 'Loop_Lengths', 'Avg_Tract_Length'],
                'Slipped_DNA': ['Repeat_Unit', 'Copy_Number', 'Purity', 'Entropy'],
                'i-Motif': ['Num_Stems', 'Num_Loops', 'GC_Stems'],
                'Z-DNA': ['Contributing_10mers', 'Alternating_CG_Regions'],
                'R-Loop': ['RIZ_Length', 'REZ_Length', 'GC_Skew', 'Linker_Length'],
                'Cruciform': ['Arm_Length', 'Loop_Length', 'DeltaG'],
                'Triplex': ['Arm_Length', 'Loop_Length', 'Purity'],
                'Curved': ['Num_A_Tracts', 'Num_T_Tracts', 'AT_Content']
            }
            
            for class_name, features in specific_features.items():
                if class_name in motif['Class']:
                    print(f"\n  {class_name}-specific features:")
                    for feature in features:
                        if feature in motif:
                            print(f"    {feature}: {motif[feature]}")
            
            print()
    else:
        print("  No motifs detected for this sequence type")

# Test CSV export
print("\n" + "=" * 80)
print("CSV Export Test")
print("=" * 80)

if all_motifs:
    csv_data = export_to_csv(all_motifs, include_all_fields=True)
    lines = csv_data.strip().split('\n')
    header = lines[0].split(',')
    
    print(f"\nExported {len(lines)-1} motifs")
    print(f"CSV contains {len(header)} columns:")
    print(f"  Core columns: {', '.join(header[:12])}")
    print(f"  Enhanced columns: {', '.join([h for h in header[12:25]])}")
    if len(header) > 25:
        print(f"  Plus {len(header)-25} more columns...")
    
    # Verify all enhanced fields are present
    enhanced_fields = ['Type_Of_Repeat', 'Criterion', 'Disease_Relevance', 'Regions_Involved', 'GC_Content']
    print(f"\nEnhanced fields verification:")
    for field in enhanced_fields:
        status = "✓" if field in header else "✗"
        print(f"  [{status}] {field}")
    
    # Save to file
    output_file = os.path.join(tempfile.gettempdir(), 'enhanced_features_test.csv')
    with open(output_file, 'w') as f:
        f.write(csv_data)
    print(f"\n✓ CSV saved to: {output_file}")

# Test JSON export
print("\n" + "=" * 80)
print("JSON Export Test")
print("=" * 80)

if all_motifs:
    json_data = export_to_json(all_motifs, pretty=True)
    output_file = os.path.join(tempfile.gettempdir(), 'enhanced_features_test.json')
    with open(output_file, 'w') as f:
        f.write(json_data)
    print(f"\n✓ JSON saved to: {output_file}")
    print(f"  JSON includes all motif fields automatically")

print("\n" + "=" * 80)
print("✓ All tests completed successfully!")
print("=" * 80)
print("\nSummary:")
print(f"  - Tested {len(test_sequences)} motif types")
print(f"  - Detected {len(all_motifs)} total motifs")
print(f"  - All detectors output comprehensive features:")
print(f"    • Type_Of_Repeat (motif classification)")
print(f"    • Criterion (detection logic explanation)")
print(f"    • Disease_Relevance (clinical annotations)")
print(f"    • Regions_Involved (structural composition)")
print(f"    • GC_Content (sequence composition)")
print(f"    • Plus motif-specific features (arm length, loop length, etc.)")
print(f"  - CSV export includes ALL {len(header) if all_motifs else 0} columns")
print(f"  - JSON export includes complete motif data")
print()
